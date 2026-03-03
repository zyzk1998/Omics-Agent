"""影像组学智能体（Radiomics Agent）"""
import json
import logging
import os
from typing import Dict, Any, List, Optional, AsyncIterator

from ..base_agent import BaseAgent
from ...core.llm_client import LLMClient
from ...core.prompt_manager import PromptManager
from ...core.workflows.radiomics_workflow import RadiomicsWorkflow

logger = logging.getLogger(__name__)


RADIOMICS_INSTRUCTION = """You are a Senior Radiomics Expert specializing in Medical Imaging, CT/MRI, Texture Analysis, Biomarker Discovery, and Rad-Score Modeling.

**CAPABILITIES:**
- Full clinical pipeline: Load image -> Preprocessing (resampling, normalization) -> Feature extraction -> Feature selection -> Rad-Score calculation -> Risk probability (Sigmoid) -> Score visualization.
- Analyze CT, MRI, and other medical imaging (NIfTI, DICOM).
- Interpret radiomics features: shape, first-order statistics, and texture (GLCM, GLRLM, GLDM, NGTDM, GLSZM).
- Preprocessing: isotropic resampling (e.g. 1×1×1 mm) and intensity normalization before extraction.
- Rad-Score: pre-defined signature (weighted sum of key features, e.g. Entropy, Sphericity) and risk probability.
- Support tumor heterogeneity assessment and imaging biomarker extraction.
- Explain mid-slice previews, feature CSV, Rad-Score, and risk in plain language.

**TERMINOLOGY:**
- NIfTI (.nii, .nii.gz), DICOM (.dcm)
- ROI, mask, segmentation
- PyRadiomics, SimpleITK
- Preprocessing, resampling, normalization
- Rad-Score, signature, risk probability, Sigmoid
- Shape features, first-order, texture, GLCM, Haralick

**CONTEXT:**
Data are 3D medical images and optional segmentation masks. Generate diagnosis and recommendations in Simplified Chinese (简体中文). Focus on image dimensions, spacing, preprocessing, radiomics features, and Rad-Score interpretation."""


# 前端/编排器展示用名称，避免被错误显示为「转录组」
MODALITY_DISPLAY = "影像组 (Radiomics)"


class RadiomicsAgent(BaseAgent):
    """影像组学智能体：处理医学影像与影像组学分析。"""

    def __init__(
        self,
        llm_client: LLMClient,
        prompt_manager: PromptManager,
        radiomics_config: Optional[Dict[str, Any]] = None,
    ):
        super().__init__(llm_client, prompt_manager, "radiomics_expert")
        self.name = "Radiomics"
        self.modality_display = MODALITY_DISPLAY
        self.workflow = RadiomicsWorkflow()
        self.radiomics_config = radiomics_config or {}

    async def process_query(
        self,
        query: str,
        history: List[Dict[str, str]] = None,
        uploaded_files: List[Dict[str, str]] = None,
        **kwargs,
    ) -> Dict[str, Any]:
        """处理用户查询：意图识别 -> 解释文件 / 生成工作流 / 聊天。"""
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
            logger.warning("RadiomicsAgent 意图检测失败，使用回退: %s", e)

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
1. "explain_file" - 用户想了解文件内容/结构（如：这是什么影像、文件里有什么）
2. "run_workflow" - 用户想执行分析流程（如：分析一下、提取影像组学特征、做 Radiomics）
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
            logger.warning("RadiomicsAgent 意图解析失败: %s", e)
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
                f"类型: {inspection.get('file_type', 'N/A')} / 域: {inspection.get('domain', 'N/A')}\n"
                f"形状/元数据: {inspection.get('shape', {})}\n"
            )
            return {"type": "chat", "response": self._stream_string_response(summary)}
        except Exception as e:
            logger.exception("RadiomicsAgent explain_file 失败: %s", e)
            return {"type": "chat", "response": self._stream_string_response(f"解释文件时出错: {e}")}

    def _is_workflow_request(self, query_lower: str, file_paths: List[str]) -> bool:
        """是否为工作流类请求。"""
        workflow_kw = [
            "规划", "流程", "workflow", "pipeline", "分析", "run", "执行",
            "plan", "做一下", "跑一下", "分析一下", "全流程", "radiomics", "影像组学",
            "提取特征", "extract feature", "nifti", "dicom", "ct", "mri",
        ]
        if query_lower and any(kw in query_lower for kw in workflow_kw):
            return True
        if file_paths and (not query_lower or len(query_lower.strip()) < 5):
            non_w = ["你好", "hello", "hi", "介绍", "你是谁", "who are you"]
            if not query_lower or query_lower.strip().lower() not in [q.lower() for q in non_w]:
                return True
        return False

    async def _generate_workflow_config(self, query: str, file_paths: List[str]) -> Dict[str, Any]:
        """使用 RadiomicsWorkflow 生成工作流配置。"""
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
                logger.warning("RadiomicsAgent 文件检查失败，使用路径: %s", e)
                file_metadata = {"file_path": file_paths[-1]}
        if not file_metadata:
            file_metadata = {}
        return self.workflow.generate_template(target_steps=None, file_metadata=file_metadata)

    def _stream_string_response(self, text: str) -> AsyncIterator[str]:
        async def _gen():
            yield text
        return _gen()

    async def _generate_analysis_summary(
        self,
        steps_results: List[Dict[str, Any]],
        omics_type: str = "Radiomics",
        workflow_name: str = "影像组学流程",
        summary_context: Optional[Dict[str, Any]] = None,
        output_dir: Optional[str] = None,
    ) -> Optional[str]:
        """Hard override: do NOT use BaseAgent. Generate summary with strictly Radiomics prompt via LLM. FORBIDDEN: Metabolomics, Genes, Cells, LC-MS."""
        import json
        system_prompt = """You are a Medical Imaging Expert. Summarize the radiomics analysis of the NIfTI/medical image.

You MUST report:
1. Image dimensions and mask status (alignment, ROI).
2. Preprocessing (resampling, normalization).
3. Extracted texture features (shape, first-order, GLCM etc.) and Rad-Score.
4. Risk level or probability if available.

Use ONLY these terms: NIfTI, DICOM, ROI, mask, PyRadiomics, texture features, Rad-Score, risk probability, image size, preprocessing.

FORBIDDEN WORDS (do not use): Metabolomics, metabolites, Genes, gene expression, Cells, single cell, LC-MS, GC-MS, PCA (unless for imaging), pathway enrichment, differential metabolites, VIP, volcano plot.

Output in Simplified Chinese (简体中文), Markdown, with clear sections. Be concise (300–600 words)."""

        steps_text = []
        for i, sr in enumerate(steps_results or []):
            name = sr.get("step_name", sr.get("name", f"Step {i+1}"))
            status = sr.get("status", "")
            data = sr.get("data", {})
            steps_text.append(f"- **{name}** ({status}): {json.dumps(data, ensure_ascii=False)[:500]}")
        user_content = f"Workflow: {workflow_name}\n\nSteps:\n" + "\n".join(steps_text) + "\n\nSummarize the above radiomics run in Chinese. Report image size, mask status, Rad-Score, and risk. Do not mention metabolomics, genes, or LC-MS."
        # 🔥 User-Facing Error Translation: 若有失败步骤，追加执行日志并要求输出「数据诊断与优化建议」章节
        failed_steps_full = (summary_context or {}).get("failed_steps", [])
        if failed_steps_full:
            _max_tb = 1000
            execution_log = ""
            for fd in failed_steps_full:
                name = fd.get("name", fd.get("step_id", "Unknown"))
                err = fd.get("error", fd.get("message", ""))
                tb = fd.get("traceback", fd.get("debug_info", ""))
                if isinstance(tb, str) and len(tb) > _max_tb:
                    tb = tb[:_max_tb] + "\n... (已截断)"
                execution_log += f"\n--- {name} ---\n错误: {err}\n" + (f"Traceback:\n{tb}\n" if tb else "")
            user_content += f"\n\n【执行日志与异常诊断】\n以下是本次工作流的执行日志（含报错）：\n{execution_log}\n请将上述报错翻译为生物学小白能听懂的建议，并在报告中增加章节「⚠️ 数据诊断与优化建议」；区分数据类问题与系统类问题，委婉专业地给出建议。即使部分步骤失败，也基于已成功步骤提供有价值分析。"

        messages = [
            {"role": "system", "content": system_prompt},
            {"role": "user", "content": user_content},
        ]
        try:
            completion = await self.llm_client.achat(messages, temperature=0.2, max_tokens=1500)
            if not completion or not getattr(completion, "choices", None):
                return None
            if hasattr(self.llm_client, "extract_think_and_content"):
                _, text = self.llm_client.extract_think_and_content(completion)
            else:
                block = completion.choices[0]
                msg = getattr(block, "message", None)
                text = getattr(msg, "content", None) if msg else getattr(block, "content", "")
            return (text or "").strip() or None
        except Exception as e:
            logger.exception("RadiomicsAgent _generate_analysis_summary failed: %s", e)
            return None

    async def _stream_chat_response(
        self,
        query: str,
        file_paths: List[str],
    ) -> AsyncIterator[str]:
        """流式聊天，使用 RADIOMICS_INSTRUCTION 作为系统提示。"""
        messages = [
            {"role": "system", "content": RADIOMICS_INSTRUCTION},
            {"role": "user", "content": query},
        ]
        try:
            async for chunk in self.llm_client.astream(messages):
                if chunk.choices and chunk.choices[0].delta.content:
                    c = chunk.choices[0].delta.content
                    if c:
                        yield c
        except Exception as e:
            logger.exception("RadiomicsAgent 流式聊天失败: %s", e)
            yield f"\n\n错误: {str(e)}"
