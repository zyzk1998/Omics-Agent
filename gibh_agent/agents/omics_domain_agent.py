"""
三大组学（基因组 / 蛋白组 / 表观）共用领域智能体基类：SOP 规划、文件注册、系统提示注册。
"""
from __future__ import annotations

import logging
import os
from typing import Any, AsyncIterator, Dict, List, Optional

from .base_agent import BaseAgent
from ..core.llm_client import LLMClient
from ..core.prompt_manager import PromptManager
from ..core.tool_retriever import ToolRetriever
from ..core.planner import SOPPlanner
from ..core.workflows import WorkflowRegistry
from ..tools.general.file_inspector import inspect_file

logger = logging.getLogger(__name__)


class OmicsDomainAgent(BaseAgent):
    """基于 WorkflowRegistry + SOPPlanner 的组学智能体基类。"""

    def __init__(
        self,
        llm_client: Optional[LLMClient],
        prompt_manager: PromptManager,
        domain_name: str,
        expert_role: str,
        system_prompt: str,
    ) -> None:
        super().__init__(llm_client, prompt_manager, expert_role)
        self._domain_name = domain_name
        self._system_prompt_text = system_prompt
        if self.prompt_manager:
            self.prompt_manager.register_template(f"{expert_role}_system", system_prompt)

    def _make_sop_planner(self) -> SOPPlanner:
        llm = self.llm_client
        if llm is None:
            from ..core.llm_client import LLMClientFactory
            from ..core.llm_cloud_providers import get_default_chat_model

            llm = LLMClientFactory.create_for_model(get_default_chat_model())
        return SOPPlanner(ToolRetriever(), llm)

    def _omics_route_in_query(self, query: str) -> bool:
        return "[Omics_Route:" in (query or "")

    def _is_workflow_request(self, query: str, file_paths: List[str]) -> bool:
        q = (query or "").strip()
        ql = q.lower()
        if self._omics_route_in_query(q):
            return True
        keys = (
            "workflow",
            "pipeline",
            "全流程",
            "分析",
            "规划",
            "执行",
            "跑一下",
            "跑流",
            "标准流程",
            "质控",
            "比对",
            "变异",
            "peak",
            "搜库",
        )
        if any(k.lower() in ql for k in keys):
            return True
        if file_paths and (not q or len(q) < 6):
            return True
        return False

    async def process_query(
        self,
        query: str,
        history: List[Dict[str, str]] = None,
        uploaded_files: List[Dict[str, str]] = None,
        **kwargs,
    ) -> Dict[str, Any]:
        file_paths = self.get_file_paths(uploaded_files or [])
        query = query or ""

        if uploaded_files:
            self._refresh_context_for_new_files(uploaded_files)
            for fi in uploaded_files:
                if isinstance(fi, dict):
                    fn = fi.get("name") or fi.get("path") or "unknown"
                    fp = fi.get("path") or fi.get("file_id") or fn
                else:
                    fn = getattr(fi, "name", None) or "unknown"
                    fp = getattr(fi, "path", None) or fn
                abs_p = fp
                for p in file_paths:
                    if fn in p or str(fp).split("/")[-1] in p:
                        abs_p = p
                        break
                if file_paths and abs_p == fp:
                    abs_p = file_paths[-1]
                self.register_file(fn, abs_p, None)
            names = [
                (f.get("name") if isinstance(f, dict) else getattr(f, "name", ""))
                for f in uploaded_files
            ]
            if names:
                self.set_active_file(names[-1])

        if not file_paths:
            af = self.get_active_file_info()
            if af and af.get("path"):
                file_paths = [af["path"]]

        query_lower = query.lower().strip()

        intent = "chat"
        try:
            intent_result = await self._detect_intent_with_llm(query, file_paths, uploaded_files or [])
            intent = intent_result.get("intent", "chat")
        except Exception as e:
            logger.warning("意图检测失败，使用关键词回退: %s", e)
            intent = None

        if intent == "explain_file":
            if not file_paths:
                return {"type": "chat", "response": self._stream_string_response("请先上传数据文件后再提问。")}
            inp = file_paths[-1]
            try:
                inspection = inspect_file(inp)
                if inspection.get("error") or inspection.get("status") == "error":
                    return {
                        "type": "chat",
                        "response": self._stream_string_response(
                            f"文件检查失败: {inspection.get('message') or inspection.get('error')}"
                        ),
                    }
                explanation = await self._explain_file_with_llm(query, inspection, inp)
                return {"type": "chat", "response": self._stream_string_response(explanation)}
            except Exception as e:
                logger.error("文件解释失败: %s", e, exc_info=True)
                return {"type": "chat", "response": self._stream_string_response(f"文件处理出错: {e}")}

        if intent is None or intent == "chat":
            if self._is_workflow_request(query_lower, file_paths) or self._omics_route_in_query(query):
                return await self._generate_workflow_config(query, file_paths)
            return {"type": "chat", "response": self._stream_chat_response(query, file_paths)}

        if intent == "run_workflow":
            return await self._generate_workflow_config(query, file_paths)

        return {"type": "chat", "response": self._stream_chat_response(query, file_paths)}

    def _stream_string_response(self, text: str) -> AsyncIterator[str]:
        async def _gen():
            yield text

        return _gen()

    async def _stream_chat_response(self, query: str, file_paths: List[str]) -> AsyncIterator[str]:
        ctx = {"file_paths": file_paths}
        async for chunk in self.chat(query, context=ctx, stream=True):
            yield chunk

    async def _explain_file_with_llm(
        self, query: str, inspection_result: Dict[str, Any], file_path: str
    ) -> str:
        import json

        summary = json.dumps(inspection_result, ensure_ascii=False, indent=2)[:12000]
        messages = [
            {"role": "system", "content": self._system_prompt_text[:8000]},
            {
                "role": "user",
                "content": f"用户问题: {query}\n路径: {file_path}\ninspect:\n{summary}\n请用中文简要解释。",
            },
        ]
        try:
            _llm = self.llm_for_request()
            completion = await _llm.achat(messages, temperature=0.2, max_tokens=1200)
            _, response = _llm.extract_think_and_content(completion)
            return response or ""
        except Exception as e:
            return f"解释生成失败: {e}"

    async def _detect_intent_with_llm(
        self,
        query: str,
        file_paths: List[str],
        uploaded_files: Optional[List[Dict[str, str]]],
    ) -> Dict[str, Any]:
        import json

        names: List[str] = []
        if uploaded_files:
            for f in uploaded_files:
                if isinstance(f, dict):
                    names.append(f.get("name") or f.get("path") or "")
        elif file_paths:
            names = [os.path.basename(p) for p in file_paths]
        prompt = f"""User: {query}
Files: {', '.join(names) or 'None'}
Classify intent JSON only:
{{"intent":"explain_file"|"run_workflow"|"chat","reasoning":"..."}}"""
        messages = [
            {"role": "system", "content": "Return JSON only."},
            {"role": "user", "content": prompt},
        ]
        _llm = self.llm_for_request()
        completion = await _llm.achat(messages, temperature=0.1, max_tokens=128)
        _, response = _llm.extract_think_and_content(completion)
        text = response.strip()
        if "```" in text:
            text = text.split("```")[1].split("```")[0].strip()
        try:
            return json.loads(text)
        except Exception:
            return {"intent": "chat", "reasoning": "parse_failed"}

    async def _generate_workflow_config(self, query: str, file_paths: List[str]) -> Dict[str, Any]:
        planner = self._make_sop_planner()
        wf_registry = WorkflowRegistry()
        workflow = wf_registry.get_workflow(self._domain_name)
        file_metadata = None
        if file_paths:
            try:
                from ..core.file_inspector import FileInspector

                upload_dir = os.getenv("UPLOAD_DIR", "/app/uploads")
                inspector = FileInspector(upload_dir)
                file_metadata = inspector.inspect_file(file_paths[-1])
                if file_metadata.get("status") != "success" and not file_metadata.get("file_path"):
                    fp = file_paths[-1]
                    file_metadata = {"file_path": fp, "success": True, "status": "success"}
            except Exception as e:
                logger.warning("FileInspector 失败，使用裸路径: %s", e)
                file_metadata = {"file_path": file_paths[-1], "success": True, "status": "success"}
        try:
            plan_result: Optional[Dict[str, Any]] = None
            async for ev, data in planner.generate_plan(
                user_query=query,
                file_metadata=file_metadata,
                domain_name=self._domain_name,
            ):
                if ev == "workflow":
                    plan_result = data
            if plan_result and plan_result.get("type") != "error":
                if file_metadata and workflow:
                    try:
                        step_ids = []
                        for s in (plan_result.get("workflow_data") or {}).get("steps") or []:
                            if isinstance(s, dict):
                                sid = s.get("step_id") or s.get("id")
                                if sid:
                                    step_ids.append(sid)
                        diag = await self._perform_data_diagnosis(
                            file_metadata=file_metadata,
                            omics_type=self._domain_name,
                            system_instruction=self._system_prompt_text,
                            workflow_for_whitelist=workflow,
                            target_step_ids_for_whitelist=step_ids or None,
                        )
                        if diag:
                            plan_result["diagnosis_report"] = diag
                    except Exception as e:
                        logger.warning("诊断附加失败: %s", e)
                return plan_result
        except Exception as e:
            logger.error("SOPPlanner 失败: %s", e, exc_info=True)

        if not workflow:
            return {"type": "error", "error": "no workflow", "message": self._domain_name}
        tmpl = workflow.generate_template(
            target_steps=list(workflow.get_steps_dag().keys()),
            file_metadata=file_metadata,
        )
        return tmpl
