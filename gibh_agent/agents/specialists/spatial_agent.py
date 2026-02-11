"""空间组学智能体（Spatial Omics Agent）"""
import json
import logging
import os
from typing import Dict, Any, List, Optional, AsyncIterator

from ..base_agent import BaseAgent
from ...core.llm_client import LLMClient
from ...core.prompt_manager import PromptManager
from ...core.workflows.spatial_workflow import SpatialWorkflow

logger = logging.getLogger(__name__)


# 领域系统指令（与 RNA/Metabolomics 一致的模式）
SPATIAL_INSTRUCTION = """You are a Senior Bioinformatician specializing in Spatial Transcriptomics (e.g. 10x Visium).

**CAPABILITIES:**
- Analyze 10x Visium (Space Ranger) data: load matrices and spatial coordinates, quality control, and visualization.
- Identify spatial domains and gene expression patterns across tissue (e.g. Moran's I for spatially variable genes).
- Interpret spot-level data, histological images alignment, and spatial autocorrelation results.

**TERMINOLOGY:**
- Spot, Spots, Visium, Space Ranger
- obsm['spatial'], spatial_connectivities
- Spatially Variable Genes (SVGs), Moran's I, spatial autocorrelation
- Spatial domain, tissue region

**CONTEXT:**
Data are spatial transcriptomics (gene expression per spot with 2D coordinates). Generate diagnosis and recommendations in Simplified Chinese (简体中文). Focus on spatial-specific metrics: spot count, genes, spatial graph, and SVG detection."""


class SpatialAgent(BaseAgent):
    """空间组学智能体：处理 Visium 等空间转录组分析。"""

    def __init__(
        self,
        llm_client: LLMClient,
        prompt_manager: PromptManager,
        spatial_config: Optional[Dict[str, Any]] = None,
    ):
        super().__init__(llm_client, prompt_manager, "spatial_expert")
        self.name = "Spatial"
        self.workflow = SpatialWorkflow()
        self.spatial_config = spatial_config or {}

    async def process_query(
        self,
        query: str,
        history: List[Dict[str, str]] = None,
        uploaded_files: List[Dict[str, str]] = None,
        **kwargs,
    ) -> Dict[str, Any]:
        """
        处理用户查询：意图识别 -> 解释文件 / 生成工作流 / 聊天。
        """
        query_lower = (query or "").lower().strip()
        file_paths = self.get_file_paths(uploaded_files or [])

        if uploaded_files:
            for file_info in uploaded_files:
                if isinstance(file_info, dict):
                    fn = file_info.get("name") or file_info.get("path") or file_info.get("file_id", "unknown")
                    fp = file_info.get("path") or file_info.get("file_id", fn)
                else:
                    fn = getattr(file_info, "name", None) or getattr(file_info, "path", None) or "unknown"
                    fp = getattr(file_info, "path", None) or fn
                abs_path = fp
                if file_paths and fp not in file_paths:
                    for p in file_paths:
                        if p.endswith(os.path.basename(fp)) or fp in p:
                            abs_path = p
                            break
                self.register_file(fn, abs_path, None)
            if uploaded_files:
                latest = uploaded_files[-1]
                name = latest.get("name") or latest.get("path", "unknown") if isinstance(latest, dict) else getattr(latest, "name", "unknown")
                self.set_active_file(name)
                if file_paths:
                    file_paths = [file_paths[-1]]

        if not file_paths:
            active = self.get_active_file_info()
            if active:
                file_paths = [active["path"]]

        intent = "chat"
        try:
            intent_result = await self._detect_intent(query, query_lower, file_paths, uploaded_files)
            intent = intent_result.get("intent", "chat")
        except Exception as e:
            logger.warning("SpatialAgent 意图检测失败，使用回退: %s", e)

        if intent == "explain_file":
            if not file_paths:
                return {"type": "chat", "response": self._stream_string_response("请先上传文件后再询问。")}
            return await self._handle_explain_file(query, file_paths[-1])

        if intent == "run_workflow" or (file_paths and self._is_workflow_request(query_lower, file_paths)):
            return await self._generate_workflow_config(query, file_paths)

        return {
            "type": "chat",
            "response": self._stream_chat_response(query, file_paths),
        }

    async def _detect_intent(
        self,
        query: str,
        query_lower: str,
        file_paths: List[str],
        uploaded_files: Optional[List[Dict[str, str]]],
    ) -> Dict[str, Any]:
        """LLM 意图分类：explain_file | run_workflow | chat。"""
        file_names = []
        if uploaded_files:
            for f in uploaded_files:
                file_names.append((f.get("name") or f.get("path") or ""))
        elif file_paths:
            file_names = [os.path.basename(p) for p in file_paths]
        files_str = ", ".join(file_names) if file_names else "None"

        prompt = f"""分析用户输入，判断意图。

User Input: {query}
Uploaded Files: {files_str}

意图分类为以下之一：
1. "explain_file" - 用户想了解文件内容/结构（如：这是什么文件、文件里有什么）
2. "run_workflow" - 用户想执行分析流程（如：分析一下、跑流程、做空间转录组分析）
3. "chat" - 普通对话或咨询

只返回 JSON：
{{"intent": "explain_file"|"run_workflow"|"chat", "reasoning": "简短理由"}}"""

        messages = [
            {"role": "system", "content": "You are an intent classifier. Return JSON only."},
            {"role": "user", "content": prompt},
        ]
        try:
            completion = await self.llm_client.achat(messages, temperature=0.1, max_tokens=128)
            _, content = self.llm_client.extract_think_and_content(completion)
            raw = content.strip()
            if "```json" in raw:
                raw = raw.split("```json")[1].split("```")[0].strip()
            elif "```" in raw:
                raw = raw.split("```")[1].split("```")[0].strip()
            result = json.loads(raw)
            if result.get("intent") not in ("explain_file", "run_workflow", "chat"):
                result["intent"] = "chat"
            return result
        except Exception as e:
            logger.warning("SpatialAgent 意图解析失败: %s", e)
            return {"intent": "chat", "reasoning": str(e)}

    async def _handle_explain_file(self, query: str, input_path: str) -> Dict[str, Any]:
        """解释文件：检查 + LLM 总结。"""
        try:
            from ...core.file_inspector import FileInspector
            upload_dir = os.getenv("UPLOAD_DIR", "/app/uploads")
            inspector = FileInspector(upload_dir)
            inspection = inspector.inspect_file(input_path)
            if inspection.get("status") != "success":
                return {"type": "chat", "response": self._stream_string_response(inspection.get("error", "检查失败"))}
            summary = (
                f"文件路径: {input_path}\n"
                f"类型/形状等信息: {inspection.get('file_type', 'N/A')} / {inspection.get('shape', {})}\n"
            )
            return {"type": "chat", "response": self._stream_string_response(summary)}
        except Exception as e:
            logger.exception("SpatialAgent explain_file 失败: %s", e)
            return {"type": "chat", "response": self._stream_string_response(f"解释文件时出错: {e}")}

    def _is_workflow_request(self, query_lower: str, file_paths: List[str]) -> bool:
        """是否为工作流类请求。"""
        workflow_kw = [
            "规划", "流程", "workflow", "pipeline", "分析", "run", "执行",
            "plan", "做一下", "跑一下", "分析一下", "全流程", "visium", "空间转录",
        ]
        if query_lower and any(kw in query_lower for kw in workflow_kw):
            return True
        if file_paths and (not query_lower or len(query_lower.strip()) < 5):
            non_w = ["你好", "hello", "hi", "介绍", "你是谁", "who are you"]
            if not query_lower or query_lower.strip().lower() not in [q.lower() for q in non_w]:
                return True
        return False

    async def _generate_workflow_config(self, query: str, file_paths: List[str]) -> Dict[str, Any]:
        """使用 SpatialWorkflow 生成工作流配置。"""
        file_metadata = None
        if file_paths:
            try:
                from ...core.file_inspector import FileInspector
                upload_dir = os.getenv("UPLOAD_DIR", "/app/uploads")
                inspector = FileInspector(upload_dir)
                file_metadata = inspector.inspect_file(file_paths[-1])
                if file_metadata.get("status") != "success":
                    file_metadata = {"file_path": file_paths[-1]}
            except Exception as e:
                logger.warning("SpatialAgent 文件检查失败，使用路径: %s", e)
                file_metadata = {"file_path": file_paths[-1]}
        if not file_metadata:
            file_metadata = {}
        workflow_result = self.workflow.generate_template(
            target_steps=None,
            file_metadata=file_metadata,
        )
        return workflow_result

    def _stream_string_response(self, text: str) -> AsyncIterator[str]:
        """将字符串包装为异步生成器。"""
        async def _gen():
            yield text
        return _gen()

    async def _stream_chat_response(
        self,
        query: str,
        file_paths: List[str],
    ) -> AsyncIterator[str]:
        """流式聊天，使用 SPATIAL_INSTRUCTION 作为系统提示。"""
        system_prompt = SPATIAL_INSTRUCTION
        messages = [
            {"role": "system", "content": system_prompt},
            {"role": "user", "content": query},
        ]
        try:
            async for chunk in self.llm_client.astream(messages):
                if chunk.choices and chunk.choices[0].delta.content:
                    c = chunk.choices[0].delta.content
                    if c:
                        yield c
        except Exception as e:
            logger.exception("SpatialAgent 流式聊天失败: %s", e)
            yield f"\n\n错误: {str(e)}"
