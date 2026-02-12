"""
基础智能体抽象类
所有领域智能体都继承此类
"""
from abc import ABC, abstractmethod
from typing import Dict, Any, Optional, List, AsyncIterator
import logging
from openai import AuthenticationError, APIError
from ..core.llm_client import LLMClient
from ..core.prompt_manager import PromptManager, DATA_DIAGNOSIS_PROMPT
from ..core.data_diagnostician import DataDiagnostician
from ..core.stream_utils import strip_suggestions_from_text

logger = logging.getLogger(__name__)


class BaseAgent(ABC):
    """
    基础智能体抽象类
    
    所有领域智能体都应该继承此类并实现：
    - process_query: 处理用户查询
    - generate_workflow: 生成工作流（如果需要）
    """
    
    def __init__(
        self,
        llm_client: LLMClient,
        prompt_manager: PromptManager,
        expert_role: str
    ):
        """
        初始化基础智能体
        
        Args:
            llm_client: LLM 客户端
            prompt_manager: 提示管理器
            expert_role: 专家角色名称（如 "rna_expert"）
        """
        self.llm_client = llm_client
        self.prompt_manager = prompt_manager
        self.expert_role = expert_role
        self.diagnostician = DataDiagnostician()
        # 🔥 架构重构：会话级文件注册表
        self.context: Dict[str, Any] = {
            "file_registry": {},  # Key: filename, Value: {path, metadata, timestamp}
            "active_file": None   # 当前活动的文件名
        }
    
    def register_file(
        self,
        filename: str,
        file_path: str,
        file_metadata: Optional[Dict[str, Any]] = None
    ) -> None:
        """
        注册文件到会话注册表
        
        🔥 架构重构：维护文件历史，而不是清除
        
        Args:
            filename: 文件名（用作注册表的 key）
            file_path: 文件的绝对路径
            file_metadata: 文件元数据（可选，稍后可以更新）
        """
        import time
        if "file_registry" not in self.context:
            self.context["file_registry"] = {}
        
        self.context["file_registry"][filename] = {
            "path": file_path,
            "metadata": file_metadata,
            "timestamp": time.time(),
            "registered_at": time.strftime("%Y-%m-%d %H:%M:%S")
        }
        logger.info(f"📝 [FileRegistry] Registered file: {filename} (Total: {len(self.context['file_registry'])} files)")
    
    def set_active_file(self, filename: str) -> None:
        """
        设置当前活动的文件
        
        Args:
            filename: 文件名（必须在注册表中存在）
        """
        if filename not in self.context.get("file_registry", {}):
            logger.warning(f"⚠️ [FileRegistry] File {filename} not in registry. Registering...")
            # 如果文件不在注册表中，尝试注册（使用路径作为文件名）
            self.register_file(filename, filename)
        
        old_active = self.context.get("active_file")
        self.context["active_file"] = filename
        
        if old_active != filename:
            logger.info(f"🔄 [FileRegistry] Active file changed: {old_active} -> {filename}")
        else:
            logger.debug(f"✅ [FileRegistry] Active file unchanged: {filename}")
    
    def get_active_file_info(self) -> Optional[Dict[str, Any]]:
        """
        获取当前活动文件的信息
        
        🔥 架构重构：统一接口获取活动文件信息
        
        Returns:
            包含 path 和 metadata 的字典，如果没有活动文件返回 None
        """
        active_file = self.context.get("active_file")
        if not active_file:
            logger.debug("⚠️ [FileRegistry] No active file set")
            return None
        
        registry = self.context.get("file_registry", {})
        if active_file not in registry:
            logger.warning(f"⚠️ [FileRegistry] Active file '{active_file}' not found in registry")
            return None
        
        file_info = registry[active_file]
        logger.debug(f"✅ [FileRegistry] Retrieved active file info: {active_file}")
        return {
            "filename": active_file,
            "path": file_info.get("path"),
            "metadata": file_info.get("metadata"),
            "timestamp": file_info.get("timestamp")
        }
    
    def _refresh_context_for_new_files(self, uploaded_files: List[Dict[str, str]]) -> None:
        """
        刷新上下文以处理新文件
        
        🔥 修复：当新文件上传时，清除旧的上下文，确保使用新文件作为单一数据源
        
        Args:
            uploaded_files: 当前请求中的文件列表
        """
        if uploaded_files and len(uploaded_files) > 0:
            # 提取文件名用于日志
            file_names = []
            for file_info in uploaded_files:
                if isinstance(file_info, dict):
                    name = file_info.get("name") or file_info.get("path") or file_info.get("file_id", "unknown")
                else:
                    name = getattr(file_info, "name", None) or getattr(file_info, "path", None) or "unknown"
                file_names.append(name)
            
            # 清除旧的上下文
            old_file_paths = self.context.get("file_paths", [])
            old_file_metadata = self.context.get("file_metadata")
            
            if old_file_paths or old_file_metadata:
                logger.info(f"🔄 [System] Context refreshed. Clearing old context:")
                logger.info(f"   Old files: {old_file_paths}")
                logger.info(f"   New active files: {file_names}")
            
            # 清除文件相关的上下文
            self.context.pop("file_paths", None)
            self.context.pop("file_metadata", None)
            self.context.pop("diagnosis_report", None)
            self.context.pop("diagnosis_stats", None)
            
            logger.info(f"✅ [System] Context refreshed. New active file: {file_names[0] if file_names else 'None'}")
    
    @abstractmethod
    async def process_query(
        self,
        query: str,
        history: List[Dict[str, str]] = None,
        uploaded_files: List[Dict[str, str]] = None,
        **kwargs
    ) -> Dict[str, Any]:
        """
        处理用户查询（抽象方法）
        
        Args:
            query: 用户查询文本
            history: 对话历史
            uploaded_files: 上传的文件列表
            **kwargs: 其他参数
        
        Returns:
            处理结果字典
        """
        pass
    
    async def chat(
        self,
        query: str,
        context: Dict[str, Any] = None,
        stream: bool = False
    ) -> AsyncIterator[str]:
        """
        通用聊天方法
        
        Args:
            query: 用户查询
            context: 上下文信息
            stream: 是否流式输出
        
        Yields:
            响应文本块
        """
        context = context or {}
        
        # 获取系统提示词
        system_prompt = self.prompt_manager.get_system_prompt(
            self.expert_role,
            context
        )
        
        messages = [
            {"role": "system", "content": system_prompt},
            {"role": "user", "content": query}
        ]
        
        try:
            if stream:
                # 流式输出：直接传递内容，让前端处理 think 标签
                # DeepSeek-R1 的 think 过程会以 <think>...</think> 标签形式返回
                # 也支持旧协议的 <think>...</think> 标签
                has_yielded = False
                try:
                    async for chunk in self.llm_client.astream(messages):
                        if chunk.choices and chunk.choices[0].delta.content:
                            content = chunk.choices[0].delta.content
                            if content:
                                # 直接传递内容，前端会检测和处理 think 标签（<think> 或 <think>）
                                yield content
                                has_yielded = True
                except Exception as stream_error:
                    logger.error(f"❌ 流式响应错误: {stream_error}", exc_info=True)
                    if not has_yielded:
                        yield f"\n\n❌ 错误: {str(stream_error)}\n\n请检查服务器日志获取详细信息。"
                    else:
                        # 如果已经有一些输出，只记录错误，不重复输出错误信息
                        logger.warning(f"⚠️ 流式响应中断，但已有部分输出")
            else:
                completion = await self.llm_client.achat(messages)
                # 提取 think 过程和实际内容
                think_content, actual_content = self.llm_client.extract_think_and_content(completion)
                
                # 如果有 think 内容，包装在标签中
                if think_content:
                    yield f"<think>{think_content}</think>\n\n{actual_content}"
                else:
                    yield actual_content
        except AuthenticationError as e:
            error_msg = (
                f"\n\n❌ 认证错误 (Error code: 401 - Invalid token)\n"
                f"请检查 API 密钥是否正确设置。\n"
                f"设置方法: export SILICONFLOW_API_KEY='your_api_key_here'\n"
                f"详细错误: {str(e)}"
            )
            logger.error(f"API 认证失败: {e}")
            yield error_msg
        except APIError as e:
            error_msg = (
                f"\n\n❌ API 错误 (Error code: {getattr(e, 'status_code', 'unknown')})\n"
                f"详细错误: {str(e)}"
            )
            logger.error(f"API 调用失败: {e}")
            yield error_msg
        except Exception as e:
            error_msg = f"\n\n❌ 错误: {str(e)}"
            logger.error(f"聊天处理失败: {e}", exc_info=True)
            yield error_msg
    
    def get_file_paths(self, uploaded_files: List[Dict[str, str]]) -> List[str]:
        """
        从上传文件列表中提取文件路径，并转换为绝对路径
        
        核心原则：智能体只处理文件路径（字符串），不处理二进制数据
        **关键修复**：确保返回绝对路径，避免 "File Not Found" 错误
        
        Args:
            uploaded_files: 文件列表（可能包含相对路径或 file_id）
        
        Returns:
            绝对文件路径列表
        """
        import os
        from pathlib import Path
        
        # 获取上传目录（与 server.py 保持一致）
        upload_dir = Path(os.getenv("UPLOAD_DIR", "/app/uploads"))
        
        paths = []
        for file_info in uploaded_files:
            if isinstance(file_info, dict):
                path = file_info.get("path") or file_info.get("name") or file_info.get("file_path") or file_info.get("file_id")
            else:
                path = getattr(file_info, "path", None) or getattr(file_info, "name", None) or getattr(file_info, "file_path", None) or getattr(file_info, "file_id", None)
            
            if not path:
                continue
            
            # 🔥 修复：转换为绝对路径
            path_obj = Path(path)
            
            # 如果已经是绝对路径，直接使用
            if path_obj.is_absolute():
                absolute_path = str(path_obj.resolve())
            else:
                # 如果是相对路径，拼接 UPLOAD_DIR
                absolute_path = str((upload_dir / path_obj).resolve())
            
            # 验证路径是否存在（如果不存在，记录警告但继续处理，让调用方处理错误）
            if not os.path.exists(absolute_path):
                logger.warning(f"⚠️ 文件路径不存在: {absolute_path} (原始路径: {path})")
                # 仍然添加到列表，让调用方处理（可能文件稍后会被创建）
            
            paths.append(absolute_path)
        
        return paths
    
    def detect_file_type(self, file_path: str) -> str:
        """
        检测文件类型
        
        Args:
            file_path: 文件路径或目录路径
        
        Returns:
            文件类型（如 "fastq", "bam", "h5ad", "10x_mtx"）
        """
        import os
        
        # 如果是目录，检查是否是 FASTQ 目录或 10x MTX 目录
        if os.path.isdir(file_path):
            # 检查是否是 FASTQ 目录（包含 .fastq 或 .fq 文件）
            fastq_files = [f for f in os.listdir(file_path) if f.endswith(('.fastq', '.fq', '.fastq.gz', '.fq.gz'))]
            if fastq_files:
                return "fastq"
            
            # 检查是否是 10x MTX 目录（包含 matrix.mtx 或 matrix.mtx.gz）
            mtx_files = [f for f in os.listdir(file_path) if 'matrix.mtx' in f.lower()]
            if mtx_files:
                return "10x_mtx"
            
            # 检查是否是 Cell Ranger 输出目录（包含 filtered_feature_bc_matrix）
            if 'filtered_feature_bc_matrix' in os.listdir(file_path) or \
               any('filtered_feature_bc_matrix' in subdir for subdir in os.listdir(file_path) if os.path.isdir(os.path.join(file_path, subdir))):
                return "10x_mtx"
            
            return "directory"
        
        # 如果是文件，检查扩展名
        ext = file_path.split('.')[-1].lower()
        
        type_mapping = {
            "fastq": ["fastq", "fq"],
            "bam": ["bam"],
            "h5ad": ["h5ad"],
            "mtx": ["mtx"],
            "vcf": ["vcf"],
            "bed": ["bed"],
            "bw": ["bw", "bigwig"],
            "sam": ["sam"],
            "csv": ["csv"]  # 代谢组学数据通常使用 CSV 格式
        }
        
        for file_type, extensions in type_mapping.items():
            if ext in extensions:
                return file_type
        
        return "unknown"
    
    async def _perform_data_diagnosis(
        self,
        file_metadata: Dict[str, Any],
        omics_type: str,
        dataframe: Optional[Any] = None,
        system_instruction: Optional[str] = None
    ) -> Optional[str]:
        """
        执行数据诊断并生成 Markdown 报告
        
        🔥 架构重构：使用策略模式，接受 domain-specific system_instruction
        
        这是统一的数据诊断入口，所有 Agent 都应该调用此方法。
        
        🔥 TASK 3 & 4: 集成诊断缓存机制
        - 首先检查是否有已保存的诊断结果
        - 如果有，直接返回缓存的诊断报告
        - 如果没有，执行诊断并保存结果
        
        Args:
            file_metadata: FileInspector 返回的文件元数据
            omics_type: 组学类型（"scRNA", "Metabolomics", "BulkRNA", "default"）
            dataframe: 可选的数据预览（DataFrame 或 AnnData）
            system_instruction: 领域特定的系统指令（由各个 Agent 提供）
        
        Returns:
            Markdown 格式的诊断报告，如果失败返回 None
        """
        try:
            # 🔥 TASK 4: 检查诊断缓存
            from ..core.diagnosis_cache import get_diagnosis_cache
            cache = get_diagnosis_cache()
            
            file_path = file_metadata.get("file_path", "")
            if file_path:
                cached_diagnosis = cache.load_diagnosis(file_path)
                if cached_diagnosis:
                    logger.info(f"✅ [DataDiagnostician] 使用缓存的诊断结果: {file_path}")
                    # 返回缓存的诊断报告
                    return cached_diagnosis.get("diagnosis_report")
            
            logger.info(f"🔍 [DataDiagnostician] 开始数据诊断 - 组学类型: {omics_type}")
            
            # Step 1: 使用 DataDiagnostician 计算统计事实
            diagnosis_result = self.diagnostician.analyze(
                file_metadata=file_metadata,
                omics_type=omics_type,
                dataframe=dataframe
            )
            
            if diagnosis_result.get("status") != "success":
                logger.warning(f"⚠️ 数据诊断失败: {diagnosis_result.get('error')}")
                return None
            
            stats = diagnosis_result.get("stats", {})
            omics_type = diagnosis_result.get("omics_type") or ""
            logger.info(f"✅ [DataDiagnostician] 统计计算完成: {len(stats)} 个指标, omics_type={omics_type}")

            # Step 2: 构建 LLM Prompt（按领域选择模板，Radiomics 使用影像专用模板）
            # 将统计事实格式化为 JSON 字符串
            import json
            try:
                stats_json = json.dumps(stats, ensure_ascii=False, indent=2)
                logger.debug(f"📝 [DEBUG] Stats JSON length: {len(stats_json)}")
            except Exception as json_err:
                logger.error(f"❌ [DataDiagnostician] JSON 序列化失败: {json_err}")
                stats_json = json.dumps({"error": "无法序列化统计信息"}, ensure_ascii=False)
            
            # 🔥 修复：安全地截断 JSON 字符串（而不是字典）
            # 如果 JSON 太长，截断它（但保留完整的结构）
            max_json_length = 2000  # 限制 JSON 长度
            if len(stats_json) > max_json_length:
                logger.warning(f"⚠️ Stats JSON 太长 ({len(stats_json)} 字符)，截断到 {max_json_length} 字符")
                # 截断字符串，但确保 JSON 结构完整
                truncated_json = stats_json[:max_json_length]
                # 尝试找到最后一个完整的 JSON 对象/数组边界
                last_brace = truncated_json.rfind('}')
                last_bracket = truncated_json.rfind(']')
                last_comma = max(truncated_json.rfind(','), truncated_json.rfind('\n'))
                # 选择最接近末尾的边界
                cut_point = max(last_brace, last_bracket, last_comma)
                if cut_point > max_json_length * 0.8:  # 如果截断点不太早
                    stats_json = truncated_json[:cut_point + 1] + "\n  ... (truncated)"
                else:
                    stats_json = truncated_json + "\n  ... (truncated)"
            
            # 🔥 安全地提取文件预览信息（如果可用）
            # 注意：file_metadata 是字典，不能直接切片
            head_preview = ""
            try:
                head_data = file_metadata.get("head", {})
                if isinstance(head_data, dict):
                    # head_data 是字典，包含 "markdown" 或 "json" 键
                    if "markdown" in head_data:
                        head_preview = head_data["markdown"]
                    elif "json" in head_data:
                        # 如果是 JSON 格式，转换为字符串
                        head_preview = json.dumps(head_data["json"], ensure_ascii=False, indent=2)
                    else:
                        head_preview = str(head_data)
                elif isinstance(head_data, str):
                    # 如果已经是字符串，直接使用
                    head_preview = head_data
                else:
                    head_preview = str(head_data)
                
                # 🔥 安全地截断字符串预览（不是字典）
                if len(head_preview) > 1000:
                    head_preview = head_preview[:1000] + "\n... (truncated)"
            except Exception as head_err:
                logger.warning(f"⚠️ 提取文件预览失败: {head_err}")
                head_preview = "无法提取数据预览"
            
            # 使用 PromptManager 获取诊断模板（按 omics_type 选择领域专用模板）
            try:
                if (omics_type or "").lower() in ("radiomics", "medical_image", "imaging") or stats.get("_imaging_only"):
                    from ..core.prompts.radiomics_prompts import RADIOMICS_DIAGNOSIS_TEMPLATE
                    dimensions = stats.get("dimensions_str", "N/A")
                    spacing = stats.get("spacing_str", "N/A")
                    mask_status = "已提供" if stats.get("mask_present") else "未提供"
                    origin = stats.get("origin")
                    origin_str = str(origin) if origin is not None else "N/A"
                    prompt = self.prompt_manager.get_prompt(
                        "data_diagnosis_radiomics",
                        {
                            "dimensions": dimensions,
                            "spacing": spacing,
                            "mask_status": mask_status,
                            "origin": origin_str,
                            "inspection_data": stats_json,
                        },
                        fallback=RADIOMICS_DIAGNOSIS_TEMPLATE
                    )
                else:
                    prompt = self.prompt_manager.get_prompt(
                        "data_diagnosis",
                        {
                            "inspection_data": stats_json,
                            "head_preview": head_preview[:500] if head_preview else "",
                        },
                        fallback=DATA_DIAGNOSIS_PROMPT.format(inspection_data=stats_json)
                    )
                logger.debug(f"📝 [DEBUG] Prompt length: {len(prompt)}")
            except Exception as prompt_err:
                logger.warning(f"⚠️ 获取诊断模板失败，使用默认模板: {prompt_err}")
                try:
                    if not isinstance(stats_json, str):
                        stats_json = json.dumps(stats_json, ensure_ascii=False)
                    if (omics_type or "").lower() in ("radiomics", "medical_image", "imaging") or stats.get("_imaging_only"):
                        from ..core.prompts.radiomics_prompts import RADIOMICS_DIAGNOSIS_TEMPLATE
                        dimensions = stats.get("dimensions_str", "N/A")
                        spacing = stats.get("spacing_str", "N/A")
                        mask_status = "已提供" if stats.get("mask_present") else "未提供"
                        origin_str = str(stats.get("origin", "N/A"))
                        prompt = RADIOMICS_DIAGNOSIS_TEMPLATE.replace("{{ dimensions }}", dimensions).replace("{{ spacing }}", spacing).replace("{{ mask_status }}", mask_status).replace("{{ origin }}", origin_str).replace("{{ inspection_data }}", stats_json)
                    else:
                        prompt = DATA_DIAGNOSIS_PROMPT.format(inspection_data=stats_json)
                except Exception as format_err:
                    logger.error(f"❌ [DataDiagnostician] Prompt 格式化失败: {format_err}")
                    if not isinstance(stats_json, str):
                        stats_json = json.dumps(stats_json, ensure_ascii=False)
                    prompt = f"""You are a Senior Bioinformatician specializing in {omics_type}.

Based on the following data statistics:
{stats_json}

Please generate a data diagnosis and parameter recommendation report in Simplified Chinese (简体中文).

Format:
### 🔍 数据体检报告
- **数据规模**: [样本数、代谢物数]
- **数据特征**: [缺失值率、数据范围等]
- **数据质量**: [质量评估]

### 💡 参数推荐
Create a Markdown table with parameter recommendations.

Use Simplified Chinese for all content."""
            
            # Step 3: 调用 LLM 生成 Markdown 报告
            # 🔥 CRITICAL FIX: 强制注入统计数据到系统提示，防止 LLM 产生幻觉（按领域区分）
            stats_facts = []
            if (omics_type or "").lower() in ("radiomics", "medical_image", "imaging") or stats.get("_imaging_only"):
                stats_facts.append(f"影像尺寸: {stats.get('dimensions_str', 'N/A')}；层厚/间距: {stats.get('spacing_str', 'N/A')}；掩膜: {'已提供' if stats.get('mask_present') else '未提供'}。")
            elif omics_type.lower() in ["metabolomics", "metabolomic", "metabonomics"]:
                n_samples = stats.get("n_samples", 0)
                n_metabolites = stats.get("n_metabolites", 0)
                missing_rate = stats.get("missing_rate", 0)
                stats_facts.append(f"数据集包含 {n_samples} 个样本和 {n_metabolites} 个代谢物。")
                if missing_rate > 0:
                    stats_facts.append(f"缺失值率为 {missing_rate:.2f}%。")
            elif omics_type.lower() in ["scrna", "scrna-seq", "single_cell", "single-cell"]:
                n_cells = stats.get("n_cells", 0)
                n_genes = stats.get("n_genes", 0)
                stats_facts.append(f"数据集包含 {n_cells} 个细胞和 {n_genes} 个基因。")
            else:
                n_rows = stats.get("n_rows", stats.get("n_samples", 0))
                n_cols = stats.get("n_cols", stats.get("n_features", 0))
                stats_facts.append(f"数据集包含 {n_rows} 行和 {n_cols} 列。")
            
            # 构建强制事实字符串
            facts_str = " ".join(stats_facts) if stats_facts else "统计数据已提供在用户提示中。"
            
            # 🔥 CRITICAL DEBUGGING: 检查数据是否缺失，如果缺失则注入调试跟踪
            n_samples_value = stats.get("n_samples", stats.get("n_cells", stats.get("n_rows", 0)))
            debug_trace = file_metadata.get("debug_trace")
            
            # 如果数据缺失（0 samples），强制注入调试跟踪到系统提示
            debug_section = ""
            if n_samples_value == 0 and debug_trace:
                debug_section = f"""

**🔍 CRITICAL FAILURE: 数据检查返回 0 个样本**

这是一个严重错误。数据检查器无法正确读取数据文件。请在诊断报告末尾添加一个名为 "🔍 调试日志 (Debug Log)" 的章节，并将以下执行跟踪完整复制到该章节中：

```
{debug_trace}
```

这个调试日志将帮助诊断问题所在。"""
                logger.warning(f"⚠️ [DataDiagnostician] Detected 0 samples, injecting debug trace into prompt")
            
            # 🔥 架构重构：使用策略模式，从 Agent 传入 system_instruction
            if system_instruction:
                # 使用 Agent 提供的领域特定指令，并强制注入统计数据
                system_prompt = f"""{system_instruction}

**CRITICAL: 数据事实（必须严格遵循，不得产生幻觉）**
{facts_str}
请确保诊断报告中的数字与上述事实完全一致。不要猜测或编造不同的数字。{debug_section}"""
                logger.debug(f"✅ [DataDiagnostician] Using domain-specific system instruction with facts (length: {len(system_prompt)})")
            else:
                # 回退到通用指令（向后兼容），但也注入统计数据
                system_prompt = f"""You are a Senior Bioinformatician. Generate data diagnosis and parameter recommendations in Simplified Chinese.

**CRITICAL: 数据事实（必须严格遵循，不得产生幻觉）**
{facts_str}
请确保诊断报告中的数字与上述事实完全一致。不要猜测或编造不同的数字。{debug_section}"""
                logger.warning(f"⚠️ [DataDiagnostician] No system_instruction provided, using generic prompt with facts")
            
            # 🔥 架构重构：将 system_instruction 前置到用户 prompt（确保上下文隔离）
            if system_instruction:
                # 在用户 prompt 前添加系统指令，确保 LLM 理解领域约束
                prompt = f"""{system_instruction}

**数据事实（必须严格遵循）：**
{facts_str}

{prompt}"""
            
            messages = [
                {"role": "system", "content": system_prompt},
                {"role": "user", "content": prompt}
            ]
            
            # 🔥 Step 3: 调用 LLM 生成 Markdown 报告
            # 🔥 CRITICAL DEBUGGING: 包装在详细的 try-except 中
            # 🔥 TASK 2: 统一LLM客户端获取逻辑（与规划阶段一致）
            try:
                # 🔥 TASK 2: 如果 self.llm_client 不可用，使用工厂方法创建（与规划阶段一致）
                llm_client_to_use = self.llm_client
                if not llm_client_to_use:
                    logger.warning("⚠️ [DataDiagnostician] self.llm_client 不可用，使用 LLMClientFactory.create_default()")
                    from gibh_agent.core.llm_client import LLMClientFactory
                    llm_client_to_use = LLMClientFactory.create_default()
                    logger.info(f"✅ [DataDiagnostician] 已创建默认LLM客户端: {llm_client_to_use.base_url}")
                
                logger.info(f"📞 [DataDiagnostician] 调用 LLM 生成报告...")
                logger.info(f"📊 [DataDiagnostician] 统计数据摘要: n_samples={stats.get('n_samples', stats.get('n_cells', stats.get('n_rows', 0)))}, n_features={stats.get('n_features', stats.get('n_genes', stats.get('n_metabolites', stats.get('n_cols', 0))))}")
                logger.info(f"📊 [DataDiagnostician] Stats JSON 长度: {len(stats_json)} 字符")
                if len(stats_json) < 100:
                    logger.warning(f"⚠️ [DataDiagnostician] 警告：Stats JSON 过短，可能缺少关键数据")
                logger.debug(f"📝 [DEBUG] LLM Client type: {type(llm_client_to_use)}")
                logger.debug(f"📝 [DEBUG] LLM Client methods: {dir(llm_client_to_use)}")
                
                completion = await llm_client_to_use.achat(messages, temperature=0.3, max_tokens=1500)
                logger.info(f"✅ [DataDiagnostician] LLM调用完成，开始解析响应...")
                
                logger.debug(f"📝 [DEBUG] LLM completion type: {type(completion)}")
                logger.debug(f"📝 [DEBUG] LLM completion: {completion}")
                
                think_content, response = llm_client_to_use.extract_think_and_content(completion)
                # 🔥 DRY: Strip <<<SUGGESTIONS>>> block so it never reaches the frontend
                if response:
                    response, suggestions = strip_suggestions_from_text(response)
                    self.context["diagnosis_suggestions"] = suggestions
                # 🔥 DEBUG: 打印诊断报告信息
                if response:
                    logger.info(f"✅ [DataDiagnostician] 诊断报告生成成功，长度: {len(response)}")
                    logger.debug(f"📝 [DEBUG] Diagnosis report preview: {response[:200]}...")
                else:
                    logger.warning(f"⚠️ [DataDiagnostician] 诊断报告为空")
                    logger.warning(f"⚠️ [DEBUG] Think content: {think_content[:200] if think_content else 'None'}")
                
                # Step 4: 从诊断报告中提取参数推荐
                # 🔥 TASK 5: 解析诊断报告中的参数推荐表格
                recommendation = self._extract_parameter_recommendations(response, omics_type, stats)
                
                # Step 5: 保存到上下文（供 UI 和后续步骤使用）
                self.context["diagnosis_report"] = response
                self.context["diagnosis_stats"] = stats
                if recommendation:
                    self.context["parameter_recommendation"] = recommendation
                
                # 🔥 TASK 4: 保存诊断结果到缓存
                if file_path and response:
                    cache_data = {
                        "diagnosis_report": response,
                        "stats": stats,
                        "recommendation": recommendation,
                        "omics_type": omics_type
                    }
                    cache.save_diagnosis(file_path, cache_data, file_metadata)
                    logger.info(f"✅ [DataDiagnostician] 诊断结果已保存到缓存")
                
                return response
                
            except AttributeError as attr_err:
                # LLM 客户端方法不存在
                import traceback
                error_msg = (
                    f"❌ [DataDiagnostician] LLM 客户端方法调用失败\n"
                    f"━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"
                    f"错误类型: AttributeError\n"
                    f"错误信息: {str(attr_err)}\n"
                    f"━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"
                    f"调用上下文:\n"
                    f"  - LLM Client type: {type(llm_client_to_use)}\n"
                    f"  - LLM Client base_url: {llm_client_to_use.base_url if hasattr(llm_client_to_use, 'base_url') else 'N/A'}\n"
                    f"  - LLM Client model: {llm_client_to_use.model if hasattr(llm_client_to_use, 'model') else 'N/A'}\n"
                    f"  - Available methods: {[m for m in dir(llm_client_to_use) if not m.startswith('_')]}\n"
                    f"━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"
                    f"调用参数:\n"
                    f"  - messages数量: {len(messages)}\n"
                    f"  - temperature: 0.3\n"
                    f"  - max_tokens: 1500\n"
                    f"━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"
                    f"完整堆栈:\n{traceback.format_exc()}\n"
                    f"━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
                )
                logger.error(error_msg)
                return f"⚠️ **诊断报告生成失败**\n\nLLM 客户端错误: {str(attr_err)}\n\n请检查服务器日志获取详细信息。"
                
            except Exception as llm_err:
                # LLM 调用失败 - 🔥 TASK 2: 增强错误日志输出
                import traceback
                error_type = type(llm_err).__name__
                error_msg = (
                    f"❌ [DataDiagnostician] LLM 调用失败\n"
                    f"━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"
                    f"错误类型: {error_type}\n"
                    f"错误信息: {str(llm_err)}\n"
                    f"━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"
                    f"调用上下文:\n"
                    f"  - LLM Client base_url: {llm_client_to_use.base_url if hasattr(llm_client_to_use, 'base_url') else 'N/A'}\n"
                    f"  - LLM Client model: {llm_client_to_use.model if hasattr(llm_client_to_use, 'model') else 'N/A'}\n"
                    f"  - API Key: {'已设置' if hasattr(llm_client_to_use, 'api_key') and llm_client_to_use.api_key else '未设置'}\n"
                    f"━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"
                    f"调用参数:\n"
                    f"  - messages数量: {len(messages)}\n"
                    f"  - system message长度: {len(messages[0]['content']) if messages else 0} 字符\n"
                    f"  - user message长度: {len(messages[1]['content']) if len(messages) > 1 else 0} 字符\n"
                    f"  - temperature: 0.3\n"
                    f"  - max_tokens: 1500\n"
                    f"━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"
                    f"可能原因:\n"
                    f"  - API密钥无效或过期\n"
                    f"  - 网络连接问题\n"
                    f"  - API服务暂时不可用\n"
                    f"  - 请求超时\n"
                    f"━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"
                    f"完整堆栈:\n{traceback.format_exc()}\n"
                    f"━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
                )
                logger.error(error_msg)
                
                # 🔥 TASK: 将详细错误信息存储到context，供orchestrator通过SSE发送到前端
                self.context["last_llm_error"] = {
                    "error_type": error_type,
                    "error_message": str(llm_err),
                    "error_details": error_msg,
                    "context": {
                        "llm_client_base_url": llm_client_to_use.base_url if hasattr(llm_client_to_use, 'base_url') else 'N/A',
                        "llm_client_model": llm_client_to_use.model if hasattr(llm_client_to_use, 'model') else 'N/A',
                        "api_key_set": bool(hasattr(llm_client_to_use, 'api_key') and llm_client_to_use.api_key),
                        "messages_count": len(messages),
                        "temperature": 0.3,
                        "max_tokens": 1500
                    },
                    "possible_causes": [
                        "API密钥无效或过期",
                        "网络连接问题",
                        "API服务暂时不可用",
                        "请求超时"
                    ]
                }
                
                return f"⚠️ **诊断报告生成失败**\n\n错误: {str(llm_err)}\n\n请检查服务器日志获取详细信息。"
            
        except Exception as e:
            # 整体异常处理
            import traceback
            error_msg = (
                f"数据诊断过程失败: {str(e)}\n"
                f"Error type: {type(e).__name__}\n"
                f"Stack trace:\n{traceback.format_exc()}"
            )
            logger.error(f"❌ [DataDiagnostician] {error_msg}")
            # 🔥 返回详细的错误信息，而不是 None，这样用户可以在 UI 中看到
            return f"⚠️ **诊断报告生成失败**\n\n错误: {str(e)}\n\n请检查服务器日志获取详细信息。"
    
    def _read_execution_results(self, output_dir: Optional[str] = None) -> Dict[str, Any]:
        """
        🔥 TASK 3: Read actual execution results from generated files
        
        Args:
            output_dir: Output directory path (from executor)
        
        Returns:
            Dictionary containing:
            - csv_files: List of CSV files with head rows and describe()
            - image_files: List of generated PNG images
            - summary: Text summary of findings
        """
        import pandas as pd
        from pathlib import Path
        import os
        
        results_context = {
            "csv_files": [],
            "image_files": [],
            "summary": ""
        }
        
        if not output_dir or not os.path.exists(output_dir):
            logger.warning(f"⚠️ [ResultReader] 输出目录不存在: {output_dir}")
            return results_context
        
        try:
            output_path = Path(output_dir)
            
            # 🔥 TASK 3: Read CSV files (e.g., deg_results.csv, pathway_results.csv)
            csv_files = list(output_path.glob("*.csv"))
            for csv_file in csv_files[:5]:  # Limit to first 5 CSV files
                try:
                    df = pd.read_csv(csv_file, nrows=100)  # Read first 100 rows
                    head_rows = df.head(10).to_dict(orient='records')
                    
                    # Get describe() statistics
                    numeric_cols = df.select_dtypes(include=['number']).columns
                    describe_stats = {}
                    if len(numeric_cols) > 0:
                        describe_stats = df[numeric_cols].describe().to_dict()
                    
                    results_context["csv_files"].append({
                        "filename": csv_file.name,
                        "path": str(csv_file),
                        "shape": f"{len(df)} rows × {len(df.columns)} columns",
                        "columns": list(df.columns),
                        "head_rows": head_rows[:5],  # First 5 rows
                        "statistics": describe_stats,
                        "summary": f"Found {len(df)} rows with {len(df.columns)} columns. Top columns: {', '.join(df.columns[:5])}"
                    })
                    logger.info(f"✅ [ResultReader] 读取CSV文件: {csv_file.name} ({len(df)} rows)")
                except Exception as e:
                    logger.warning(f"⚠️ [ResultReader] 无法读取CSV文件 {csv_file.name}: {e}")
            
            # 🔥 TASK 3: List image files (PNG, JPG, etc.)
            image_extensions = ['.png', '.jpg', '.jpeg', '.pdf', '.svg']
            for ext in image_extensions:
                image_files = list(output_path.glob(f"*{ext}"))
                for img_file in image_files[:10]:  # Limit to first 10 images per type
                    results_context["image_files"].append({
                        "filename": img_file.name,
                        "path": str(img_file),
                        "type": ext[1:]  # Remove dot
                    })
            
            # Build summary text
            if results_context["csv_files"]:
                csv_summary = f"发现 {len(results_context['csv_files'])} 个CSV结果文件:\n"
                for csv_info in results_context["csv_files"]:
                    csv_summary += f"- {csv_info['filename']}: {csv_info['shape']}, 列: {', '.join(csv_info['columns'][:5])}\n"
                results_context["summary"] += csv_summary
            
            if results_context["image_files"]:
                img_summary = f"\n发现 {len(results_context['image_files'])} 个图片文件:\n"
                for img_info in results_context["image_files"][:5]:
                    img_summary += f"- {img_info['filename']}\n"
                results_context["summary"] += img_summary
            
            logger.info(f"✅ [ResultReader] 读取执行结果完成: {len(results_context['csv_files'])} CSV, {len(results_context['image_files'])} 图片")
            
        except Exception as e:
            logger.error(f"❌ [ResultReader] 读取执行结果失败: {e}", exc_info=True)
        
        return results_context
    
    async def _generate_analysis_summary(
        self,
        steps_results: List[Dict[str, Any]],
        omics_type: str = "Metabolomics",
        workflow_name: str = "Analysis Pipeline",
        summary_context: Optional[Dict[str, Any]] = None,
        output_dir: Optional[str] = None  # 🔥 TASK 3: Add output_dir parameter
    ) -> Optional[str]:
        """
        基于工作流执行结果生成分析摘要（AI Expert Diagnosis）
        
        Args:
            steps_results: 步骤执行结果列表（来自 ExecutionLayer）
            omics_type: 组学类型（"Metabolomics", "scRNA", 等）
            workflow_name: 工作流名称
            summary_context: 可选的上下文信息（包含失败步骤等）
            output_dir: 输出目录路径（用于读取生成的文件）
        
        Returns:
            Markdown 格式的分析摘要，如果失败返回 None
        """
        import json
        
        try:
            logger.info(f"📝 [AnalysisSummary] 开始生成分析摘要 - 组学类型: {omics_type}")
            
            # 🔥 TASK 3: Read actual execution results from files
            execution_results = self._read_execution_results(output_dir)
            
            # 🔥 CRITICAL FIX: Extract failure information from context
            has_failures = False
            failed_steps_info = []
            if summary_context:
                has_failures = summary_context.get("has_failures", False)
                failed_steps_info = summary_context.get("failed_steps", [])
            
            # 提取关键结果
            results_summary = {
                "workflow_name": workflow_name,
                "steps_completed": len(steps_results),
                "steps": [],
                "has_failures": has_failures,
                "failed_steps": failed_steps_info
            }
            
            # 🔥 CRITICAL FIX: Process both successful and failed steps
            successful_steps = []
            failed_steps = []
            
            for step_result in steps_results:
                step_data = step_result.get("data", {})
                step_name = step_result.get("step_name", "Unknown Step")
                step_status = step_result.get("status", "unknown")
                
                if step_status == "success":
                    successful_steps.append(step_result)
                else:
                    failed_steps.append({
                        "name": step_name,
                        "status": step_status,
                        "error": step_result.get("error") or step_result.get("message", "未知错误")
                    })
                    logger.debug(f"⚠️ [AnalysisSummary] 记录失败的步骤: {step_name} (status: {step_status})")
            
            # 🔥 修复：只提取核心字段，不包含完整step_data，避免JSON过长
            # 解析成功的步骤结果
            for step_result in successful_steps:
                step_data = step_result.get("data", {})
                step_name = step_result.get("step_name", "Unknown Step")
                step_status = step_result.get("status", "unknown")
                
                # 🔥 修复：只创建最小化的step_info，不包含任何大型数据
                step_info = {
                    "name": step_name,
                    "status": step_status
                }
                
                # 🔥 修复：只从summary中提取关键指标，不包含data中的其他字段（如file_path、preview等）
                # 🔥 注意：summary可能是字符串（RNA分析）或字典（代谢组学分析）
                summary_raw = step_data.get("summary", {})
                if isinstance(summary_raw, str):
                    # 如果summary是字符串（RNA分析），将其转换为空字典，从step_data中直接提取指标
                    summary = {}
                else:
                    # 如果summary是字典（代谢组学分析），直接使用
                    summary = summary_raw if isinstance(summary_raw, dict) else {}
                
                # 根据不同的工具类型提取关键指标
                if "inspect_data" in step_name.lower() or "inspection" in step_name.lower():
                    # 🔥 修复：从step_data中提取关键指标（因为summary可能是字符串）
                    # 优先从step_data中提取，如果summary是字典，也可以从summary中提取
                    step_info["n_samples"] = step_data.get("n_samples", step_data.get("n_cells", summary.get("n_samples", summary.get("n_cells", "N/A"))))
                    step_info["n_features"] = step_data.get("n_features", step_data.get("n_genes", summary.get("n_features", summary.get("n_genes", "N/A"))))
                    step_info["missing_rate"] = step_data.get("missing_rate", summary.get("missing_rate", "N/A"))
                    # 🔥 RNA分析特定指标
                    if "n_cells" in step_data or "n_cells" in summary:
                        step_info["n_cells"] = step_data.get("n_cells", summary.get("n_cells", "N/A"))
                    if "n_genes" in step_data or "n_genes" in summary:
                        step_info["n_genes"] = step_data.get("n_genes", summary.get("n_genes", "N/A"))
                    if "mitochondrial_percentage" in step_data or "mitochondrial_percentage" in summary:
                        step_info["mitochondrial_percentage"] = step_data.get("mitochondrial_percentage", summary.get("mitochondrial_percentage", "N/A"))
                    # 🔥 不包含file_path、preview等大型数据
                
                elif "differential" in step_name.lower():
                    # 🔥 修复：只从summary提取关键指标，不包含完整的results列表
                    if summary:
                        # Use summary dict if available (from Phase 2 enhancement)
                        step_info["significant_count"] = summary.get("significant_count", summary.get("sig_count", "N/A"))
                        step_info["total_count"] = summary.get("total_metabolites", "N/A")
                        step_info["method"] = summary.get("method", "N/A")
                        step_info["case_group"] = summary.get("case_group", "N/A")
                        step_info["control_group"] = summary.get("control_group", "N/A")
                        # 🔥 修复：只保留top标记物名称，不包含完整数据
                        top_markers = summary.get("top_markers", [])
                        if top_markers:
                            step_info["top_markers"] = [
                                {
                                    "name": m.get("name", m.get("metabolite", "Unknown")),
                                    "log2fc": round(m.get("log2fc", 0), 3),  # 只保留3位小数
                                    "fdr": round(m.get("fdr", m.get("fdr_corrected_pvalue", 1.0)), 4)  # 只保留4位小数
                                }
                                for m in top_markers[:3]  # 🔥 只保留top 3
                            ]
                        else:
                            # Fallback: 从top_up/top_down提取名称（只保留名称，不包含其他数据）
                            top_up = summary.get("top_up", [])[:3]
                            top_down = summary.get("top_down", [])[:3]
                            step_info["top_up_names"] = [str(m) for m in top_up] if top_up else []
                            step_info["top_down_names"] = [str(m) for m in top_down] if top_down else []
                    else:
                        # Fallback: 只提取关键计数，不包含完整results列表
                        step_info["significant_count"] = "N/A"
                        step_info["total_count"] = "N/A"
                        step_info["method"] = "N/A"
                        step_info["case_group"] = "N/A"
                        step_info["control_group"] = "N/A"
                    
                    # 🔥 修复：不提取完整的results_list，只从summary中获取top标记物
                    # 如果summary中没有top_markers，则跳过（避免处理大型results列表）
                
                elif "plsda" in step_name.lower() or "pls-da" in step_name.lower():
                    # 🔥 修复：只从summary提取关键指标，不包含完整的vip_scores列表
                    if summary and isinstance(summary, dict):
                        # Use summary dict if available (from Phase 2 enhancement)
                        top_vip_markers = summary.get("top_vip_markers", [])
                        if top_vip_markers:
                            # 🔥 只保留top 3，只包含名称和VIP分数（保留2位小数）
                            step_info["top_vip_markers"] = [
                                {
                                    "name": v.get("name", "Unknown"),
                                    "vip": round(v.get("vip", v.get("vip_score", 0)), 2)
                                }
                                for v in top_vip_markers[:3]  # 🔥 只保留top 3
                            ]
                        else:
                            step_info["top_vip_markers"] = []
                        step_info["n_components"] = summary.get("n_components", "N/A")
                        step_info["comp1_variance"] = f"{summary.get('comp1_variance', 0):.1f}%"
                        step_info["comp2_variance"] = f"{summary.get('comp2_variance', 0):.1f}%"
                    else:
                        # 🔥 修复：不处理完整的vip_scores列表，只设置默认值
                        step_info["top_vip_markers"] = []
                        step_info["n_components"] = "N/A"
                        step_info["comp1_variance"] = "N/A"
                        step_info["comp2_variance"] = "N/A"
                
                elif "pathway" in step_name.lower() or "enrichment" in step_name.lower():
                    # 🔥 修复：只从summary提取关键指标，不包含完整的enriched_pathways列表
                    if summary and isinstance(summary, dict):
                        # Use summary dict if available (from Phase 2 enhancement)
                        step_info["enriched_pathway_count"] = summary.get("n_significant", 0)
                        top_pathways_list = summary.get("top_pathways", [])
                        if top_pathways_list:
                            # 🔥 只保留top 3通路名称，不包含完整数据
                            step_info["top_pathways"] = [
                                {"name": p if isinstance(p, str) else p.get("name", "Unknown")}
                                for p in top_pathways_list[:3]  # 🔥 只保留top 3
                            ]
                        else:
                            step_info["top_pathways"] = []
                    else:
                        # 🔥 修复：不处理完整的enriched_pathways列表，只设置默认值
                        step_info["enriched_pathway_count"] = 0
                        step_info["top_pathways"] = []
                
                elif "pca" in step_name.lower() and "visualize" not in step_name.lower():
                    # 🔥 修复：从data中提取关键指标（explained_variance在data中，不在summary中）
                    explained_variance = step_data.get("explained_variance", {})
                    if explained_variance and isinstance(explained_variance, dict):
                        pc1_var = explained_variance.get("PC1", 0)
                        pc2_var = explained_variance.get("PC2", 0)
                        # 计算总方差（前10个PC的累计）
                        total_var = sum(explained_variance.get(f"PC{i+1}", 0) for i in range(min(10, len(explained_variance))))
                        step_info["pc1_variance"] = f"{pc1_var * 100:.1f}%"
                        step_info["pc2_variance"] = f"{pc2_var * 100:.1f}%"
                        step_info["total_variance"] = f"{total_var * 100:.1f}%"
                        # RNA分析的PCA通常不评估分离质量，但可以基于方差判断
                        if total_var > 0.3:
                            step_info["separation_quality"] = "clear"
                        elif total_var > 0.15:
                            step_info["separation_quality"] = "moderate"
                        else:
                            step_info["separation_quality"] = "weak"
                    else:
                        # 🔥 修复：不处理完整的explained_variance，只设置默认值
                        step_info["pc1_variance"] = "N/A"
                        step_info["pc2_variance"] = "N/A"
                        step_info["total_variance"] = "N/A"
                        step_info["separation_quality"] = "unknown"
                
                elif "preprocess" in step_name.lower():
                    # 🔥 修复：只从summary或shape提取关键信息，不包含完整数据
                    if summary and isinstance(summary, dict):
                        shape = summary.get("shape", {})
                    else:
                        shape = step_data.get("shape", {})
                    step_info["preprocessed_rows"] = shape.get("rows", "N/A")
                    step_info["preprocessed_cols"] = shape.get("columns", "N/A")
                    # 🔥 不包含完整的shape数据或其他冗余信息
                
                # 🔥 RNA分析特定步骤：质量控制
                elif "qc" in step_name.lower() or "quality" in step_name.lower() or "filter" in step_name.lower():
                    # 🔥 修复：从data中提取关键指标（因为summary可能是字符串）
                    step_info["cells_before_qc"] = step_data.get("n_obs_before", step_data.get("cells_before_qc", "N/A"))
                    step_info["cells_after_qc"] = step_data.get("n_obs_after", step_data.get("cells_after_qc", "N/A"))
                    step_info["genes_before_qc"] = step_data.get("n_vars_before", step_data.get("genes_before_qc", "N/A"))
                    step_info["genes_after_qc"] = step_data.get("n_vars_after", step_data.get("genes_after_qc", "N/A"))
                    step_info["mitochondrial_percentage"] = step_data.get("mitochondrial_percentage", "N/A")
                
                # 🔥 RNA分析特定步骤：标记基因识别
                elif "marker" in step_name.lower() or "find_markers" in step_name.lower():
                    # 🔥 修复：从data中提取关键指标
                    step_info["n_markers"] = step_data.get("n_genes_per_cluster", step_data.get("n_markers", "N/A"))
                    step_info["n_clusters"] = step_data.get("n_clusters", "N/A")
                    # 🔥 修复：从markers_table中提取top标记基因（只保留前3个cluster的前3个基因）
                    markers_table = step_data.get("markers_table", [])
                    if markers_table and isinstance(markers_table, list) and len(markers_table) > 0:
                        # markers_table是一个列表，每个元素是一个字典，包含多个cluster的标记基因
                        top_markers = []
                        for i, marker_row in enumerate(markers_table[:3]):  # 只处理前3行
                            if isinstance(marker_row, dict):
                                # 提取每个cluster的top基因（从列名中提取cluster编号）
                                for key, value in marker_row.items():
                                    if "_names" in key and value:
                                        cluster_num = key.replace("_names", "")
                                        gene_name = value if isinstance(value, str) else str(value)
                                        if gene_name and gene_name != "None":
                                            top_markers.append({
                                                "gene": gene_name,
                                                "cluster": cluster_num,
                                                "log2fc": "N/A"  # markers_table中没有log2fc
                                            })
                                            if len(top_markers) >= 5:  # 只保留top 5
                                                break
                            if len(top_markers) >= 5:
                                break
                        if top_markers:
                            step_info["top_markers"] = top_markers[:5]
                
                # 🔥 RNA分析特定步骤：聚类
                elif "cluster" in step_name.lower() and "marker" not in step_name.lower():
                    # 🔥 修复：从data中提取关键指标
                    step_info["n_clusters"] = step_data.get("n_clusters", "N/A")
                    step_info["resolution"] = step_data.get("resolution", "N/A")
                    step_info["algorithm"] = step_data.get("algorithm", "N/A")
                
                # 🔥 RNA分析特定步骤：UMAP
                elif "umap" in step_name.lower():
                    # 🔥 修复：从data中提取关键指标（UMAP通常没有summary字段）
                    step_info["n_neighbors"] = step_data.get("n_neighbors", "N/A")
                    step_info["min_dist"] = step_data.get("min_dist", "N/A")
                    # UMAP本身不产生聚类，但可能从之前的聚类步骤获取
                    step_info["n_clusters"] = step_data.get("n_clusters", "N/A")
                
                # 🔥 RNA分析特定步骤：细胞类型注释
                elif "annotation" in step_name.lower() or "cell_type" in step_name.lower():
                    if summary:
                        step_info["cell_types"] = summary.get("cell_types", summary.get("annotated_types", []))[:5]  # 只保留前5个
                        step_info["annotation_method"] = summary.get("method", "N/A")
                
                # 🔥 修复：确保step_info不包含任何大型数据（如file_path、preview、完整results列表等）
                # 移除可能存在的冗余字段
                step_info.pop("file_path", None)
                step_info.pop("output_path", None)
                step_info.pop("preview", None)
                step_info.pop("results", None)  # 完整的results列表
                step_info.pop("vip_scores", None)  # 完整的VIP分数列表
                step_info.pop("enriched_pathways", None)  # 完整的通路列表
                step_info.pop("explained_variance", None)  # 完整的方差数据
                
                results_summary["steps"].append(step_info)
            
            # 格式化结果摘要
            # 🔥 修复：限制summary_json的长度，避免prompt过长
            # 只保留关键信息，移除冗余数据
            compact_summary = {
                "total_steps": len(steps_results),
                "successful_steps": len(successful_steps),
                "failed_steps": len(failed_steps),
                "steps": []
            }
            
            # 只提取每个步骤的关键信息，限制数据量
            for step_info in results_summary.get("steps", []):
                compact_step = {
                    "name": step_info.get("name", "Unknown"),
                    "status": step_info.get("status", "unknown")
                }
                # 🔥 修复：保留RNA分析的关键指标（用于key_findings提取）
                if "cells_after_qc" in step_info:
                    compact_step["cells_after_qc"] = step_info.get("cells_after_qc", "N/A")
                if "genes_after_qc" in step_info:
                    compact_step["genes_after_qc"] = step_info.get("genes_after_qc", "N/A")
                if "n_clusters" in step_info:
                    compact_step["n_clusters"] = step_info.get("n_clusters", "N/A")
                if "top_markers" in step_info:
                    compact_step["top_markers"] = step_info.get("top_markers", [])
                if "n_cells" in step_info:
                    compact_step["n_cells"] = step_info.get("n_cells", "N/A")
                if "n_genes" in step_info:
                    compact_step["n_genes"] = step_info.get("n_genes", "N/A")
                
                # 🔥 修复：只保留最核心的指标，移除所有冗余数据
                step_name_lower = step_info.get("name", "").lower()
                if "pca" in step_name_lower:
                    compact_step["pc1_variance"] = step_info.get("pc1_variance", "N/A")
                    compact_step["pc2_variance"] = step_info.get("pc2_variance", "N/A")
                    compact_step["separation_quality"] = step_info.get("separation_quality", "N/A")
                    # 🔥 不包含pc1_var、pc2_var、total_variance等冗余字段
                elif "differential" in step_name_lower:
                    compact_step["significant_count"] = step_info.get("significant_count", "N/A")
                    compact_step["total_count"] = step_info.get("total_count", "N/A")
                    # 🔥 只保留top 3标记物名称，不包含log2fc、fdr等详细数据
                    top_markers = step_info.get("top_markers", [])[:3]
                    if top_markers:
                        compact_step["top_marker_names"] = [m.get("name", "Unknown") for m in top_markers]
                    # 🔥 不包含top_up、top_down、method、case_group等冗余字段
                elif "plsda" in step_name_lower or "pls-da" in step_name_lower:
                    # 🔥 只保留top 3 VIP代谢物名称，不包含VIP分数
                    top_vip = step_info.get("top_vip_markers", [])[:3]
                    if top_vip:
                        compact_step["top_vip_names"] = [v.get("name", "Unknown") for v in top_vip]
                    # 🔥 不包含n_components、comp1_variance、comp2_variance等冗余字段
                elif "pathway" in step_name_lower or "enrichment" in step_name_lower:
                    compact_step["enriched_pathway_count"] = step_info.get("enriched_pathway_count", 0)
                    # 🔥 只保留top 3通路名称，不包含p_value、enrichment_score等详细数据
                    top_pathways = step_info.get("top_pathways", [])[:3]
                    if top_pathways:
                        compact_step["top_pathway_names"] = [
                            p.get("name", "Unknown") if isinstance(p, dict) else str(p) 
                            for p in top_pathways
                        ]
                elif "inspect" in step_name_lower or "inspection" in step_name_lower:
                    compact_step["n_samples"] = step_info.get("n_samples", "N/A")
                    compact_step["n_features"] = step_info.get("n_features", "N/A")
                    # 🔥 不包含missing_rate、file_path等冗余字段
                elif "preprocess" in step_name_lower:
                    compact_step["preprocessed_rows"] = step_info.get("preprocessed_rows", "N/A")
                    compact_step["preprocessed_cols"] = step_info.get("preprocessed_cols", "N/A")
                
                compact_summary["steps"].append(compact_step)
            
            summary_json = json.dumps(compact_summary, ensure_ascii=False, indent=2)
            logger.info(f"📊 [AnalysisSummary] summary_json长度: {len(summary_json)}字符")
            
            # 构建提示词
            if omics_type.lower() in ["metabolomics", "metabolomic", "metabonomics"]:
                expert_role = "代谢组学分析专家"
                domain_context = """
- 代谢物数据预处理（缺失值处理、Log2转换、标准化）
- 主成分分析（PCA）用于降维和可视化
- 差异代谢物分析（t-test/Wilcoxon）用于发现组间差异
- 火山图可视化展示差异分析结果
"""
            elif omics_type.lower() in ["scrna", "scrna-seq", "single_cell", "single-cell"]:
                expert_role = "单细胞转录组分析专家"
                domain_context = """
- 质量控制（QC）过滤低质量细胞
- 数据标准化和特征选择
- 降维分析（PCA、UMAP）
- 细胞聚类和标记基因识别
"""
            elif omics_type.lower() in ["spatial", "visium", "spatial transcriptomics"]:
                expert_role = "空间转录组学分析专家"
                domain_context = """
- 10x Visium 数据：Spot 级基因表达与空间坐标（obsm['spatial']）
- 质量控制与标准化、PCA 降维、Leiden 聚类
- 空间邻域图、空间自相关（Moran's I）识别空间可变基因（SVGs）
- 空间图按聚类/基因着色
**重要**：使用空间组学术语（Spots、Clusters、Spatial Domains、Gene Expression、Moran's I、SVGs）。不要使用代谢组学术语（代谢物、LC-MS、差异代谢物等）。
"""
            else:
                expert_role = "生物信息学分析专家"
                domain_context = "通用组学数据分析流程"
            
            # 🔥 CRITICAL FIX: Build failure information for prompt
            failure_info = ""
            if has_failures and failed_steps:
                failure_info = f"\n\n**⚠️ Failed Steps ({len(failed_steps)}/{len(steps_results)}):**\n"
                for failed_step in failed_steps:
                    failure_info += f"- **{failed_step.get('name', 'Unknown')}**: {failed_step.get('error', 'Unknown error')}\n"
                failure_info += "\n**IMPORTANT**: Some steps failed, but you should still summarize the successful steps. Explain what was accomplished and note the failures."
            
            # 🔥 CRITICAL FIX: Build rule 2 text separately to avoid f-string backslash issue
            rule2_text = ""
            if has_failures:
                rule2_text = "2. **Partial Success Handling**: ⚠️ **IMPORTANT**: Some steps failed during execution. You MUST still summarize the successful steps (e.g., PCA, PLS-DA, Volcano plots) and explain what insights can be drawn from them. For failed steps, briefly note what went wrong and why it might have failed."
            else:
                rule2_text = "2. **Complete Success**: All steps completed successfully. Provide a comprehensive analysis of all results."
            
            # 🔥 CRITICAL FIX: Build rule 3 text separately to avoid f-string backslash issue
            rule3_text = ""
            if has_failures:
                rule3_text = "- ⚠️ Some steps failed during execution. You MUST still summarize the successful steps (e.g., PCA, PLS-DA, Volcano plots) and explain what insights can be drawn from them.\n- For failed steps, briefly note what went wrong (e.g., 'Pathway enrichment failed due to missing gseapy library') but focus on interpreting the successful results."
            
            # 🔥 TASK 3: Extract Key Findings with SPECIFIC metrics (Feed the Brain)
            # 🔥 修复：根据omics_type初始化不同的key_findings结构
            is_rna_analysis = omics_type.lower() in ["scrna", "scrna-seq", "single_cell", "single-cell", "rna", "rna-seq"]
            is_spatial_analysis = omics_type.lower() in ["spatial", "visium", "spatial transcriptomics"]
            
            if is_spatial_analysis:
                key_findings = {
                    "n_spots": "N/A",
                    "n_genes": "N/A",
                    "pca_variance": {"PC1": "N/A", "PC2": "N/A"},
                    "n_clusters": "N/A",
                    "spatial_domains": "N/A",
                    "top_svgs": [],  # 空间可变基因
                    "moran_i_summary": "N/A"
                }
            elif is_rna_analysis:
                key_findings = {
                    "n_cells": "N/A",
                    "n_genes": "N/A",
                    "mitochondrial_percentage": "N/A",
                    "pca_variance": {"PC1": "N/A", "PC2": "N/A"},
                    "n_clusters": "N/A",
                    "top_marker_genes": [],  # RNA分析：标记基因
                    "cell_types": []  # RNA分析：细胞类型
                }
            else:
                key_findings = {
                    "pca_separation": "N/A",
                    "pca_variance": {"PC1": "N/A", "PC2": "N/A"},
                    "differential_count": "N/A",
                    "differential_up_down": {"up": 0, "down": 0},
                    "top_pathways": [],
                    "top_vip_metabolites": [],  # Names only for biological interpretation
                    "top_differential_metabolites": []  # Names only for biological interpretation
                }
            
            for step_info in results_summary.get("steps", []):
                step_name = step_info.get("name", "").lower()
                
                # 🔥 RNA分析特定指标提取
                if is_rna_analysis:
                    # 提取细胞和基因数量（从QC步骤，因为QC步骤有最准确的数据）
                    if "qc" in step_name or "quality" in step_name or "filter" in step_name:
                        if "cells_after_qc" in step_info:
                            key_findings["n_cells"] = step_info.get("cells_after_qc", "N/A")
                        if "genes_after_qc" in step_info:
                            key_findings["n_genes"] = step_info.get("genes_after_qc", "N/A")
                        if "mitochondrial_percentage" in step_info:
                            key_findings["mitochondrial_percentage"] = step_info.get("mitochondrial_percentage", "N/A")
                    # 如果QC步骤没有，从inspect步骤提取
                    elif ("inspect" in step_name or "inspection" in step_name) and key_findings.get("n_cells") == "N/A":
                        if "n_cells" in step_info:
                            key_findings["n_cells"] = step_info.get("n_cells", "N/A")
                        if "n_genes" in step_info:
                            key_findings["n_genes"] = step_info.get("n_genes", "N/A")
                    
                    # 提取标记基因
                    if "marker" in step_name or "find_markers" in step_name:
                        top_markers = step_info.get("top_markers", [])
                        if top_markers:
                            key_findings["top_marker_genes"] = [
                                m.get("gene", m.get("name", "Unknown")) for m in top_markers[:5]
                            ]
                    
                    # 提取聚类信息
                    if "cluster" in step_name and "marker" not in step_name:
                        if "n_clusters" in step_info:
                            key_findings["n_clusters"] = step_info.get("n_clusters", "N/A")
                    
                    # 提取细胞类型
                    if "annotation" in step_name or "cell_type" in step_name:
                        cell_types = step_info.get("cell_types", [])
                        if cell_types:
                            key_findings["cell_types"] = cell_types[:5]  # 只保留前5个
                
                # 🔥 TASK 3: Extract PCA metrics (Explained Variance) - 适用于所有分析类型
                if "pca" in step_name and "visualize" not in step_name:
                    separation = step_info.get("separation_quality", "unknown")
                    pc1_var = step_info.get("pc1_variance", step_info.get("pc1_var", "N/A"))
                    pc2_var = step_info.get("pc2_variance", step_info.get("pc2_var", "N/A"))
                    # Extract numeric values
                    if isinstance(pc1_var, str) and "%" in pc1_var:
                        pc1_val = pc1_var.replace("%", "").strip()
                    else:
                        pc1_val = str(pc1_var)
                    if isinstance(pc2_var, str) and "%" in pc2_var:
                        pc2_val = pc2_var.replace("%", "").strip()
                    else:
                        pc2_val = str(pc2_var)
                    
                    key_findings["pca_variance"] = {"PC1": pc1_val, "PC2": pc2_val}
                    # 🔥 修复：pca_separation只在代谢组学分析中存在（非 RNA、非 Spatial）
                    if not is_rna_analysis and not is_spatial_analysis:
                        if separation == "clear":
                            key_findings["pca_separation"] = f"清晰分离 (PC1: {pc1_var}, PC2: {pc2_var})"
                        else:
                            key_findings["pca_separation"] = f"中等分离 (PC1: {pc1_var}, PC2: {pc2_var})"
                
                # 🔥 空间组学：从 load_data/qc_norm 提取 n_spots、n_genes；从 cluster 提取 n_clusters；从 spatial_autocorr 提取 Moran's I、SVGs
                if is_spatial_analysis:
                    if "load" in step_name or "qc" in step_name or "norm" in step_name:
                        if "n_obs" in step_info or "n_spots" in step_info:
                            key_findings["n_spots"] = step_info.get("n_spots", step_info.get("n_obs", "N/A"))
                        if "n_vars" in step_info or "n_genes" in step_info:
                            key_findings["n_genes"] = step_info.get("n_genes", step_info.get("n_vars", "N/A"))
                    if "cluster" in step_name and "marker" not in step_name:
                        if "n_clusters" in step_info:
                            key_findings["n_clusters"] = step_info.get("n_clusters", "N/A")
                    if "spatial_autocorr" in step_name or "moran" in step_name:
                        key_findings["moran_i_summary"] = step_info.get("summary", step_info.get("n_svgs", "N/A"))
                        svgs = step_info.get("top_svgs", step_info.get("var_names", []))
                        if isinstance(svgs, list) and svgs:
                            key_findings["top_svgs"] = list(svgs)[:5]
                
                # 🔥 修复：只提取关键计数，不包含完整的top_up/top_down列表（仅代谢组学，非 Spatial）
                if not is_rna_analysis and not is_spatial_analysis and "differential" in step_name:
                    sig_count = step_info.get("significant_count", "N/A")
                    total_count = step_info.get("total_count", "N/A")
                    # 🔥 不提取top_up和top_down列表（可能很长），只使用计数
                    top_up_count = len(step_info.get("top_up", [])) if isinstance(step_info.get("top_up"), list) else 0
                    top_down_count = len(step_info.get("top_down", [])) if isinstance(step_info.get("top_down"), list) else 0
                    
                    key_findings["differential_count"] = f"发现 {sig_count} 个显著差异代谢物（共 {total_count} 个）"
                    key_findings["differential_up_down"] = {
                        "up": top_up_count,
                        "down": top_down_count
                    }
                    
                    # 🔥 修复：只提取top标记物名称，不包含完整数据
                    top_markers = step_info.get("top_markers", [])
                    if top_markers:
                        key_findings["top_differential_metabolites"] = [
                            m.get('name', 'Unknown') for m in top_markers[:3]  # 🔥 只保留top 3
                        ]
                    elif step_info.get("top_up_names"):
                        # Fallback to top_up_names（如果存在）
                        key_findings["top_differential_metabolites"] = step_info.get("top_up_names", [])[:3]
                
                # 🔥 修复：只提取top VIP代谢物名称，不包含完整数据（仅代谢组学，非 Spatial）
                if not is_rna_analysis and not is_spatial_analysis and ("plsda" in step_name or "pls-da" in step_name):
                    top_vip = step_info.get("top_vip_markers", [])
                    if top_vip:
                        # Extract metabolite NAMES only (for biological interpretation)
                        key_findings["top_vip_metabolites"] = [
                            v.get('name', 'Unknown') for v in top_vip[:3]  # 🔥 只保留top 3
                        ]
                
                # 🔥 修复：只提取top通路名称，不包含完整数据（仅代谢组学，非 Spatial）
                if not is_rna_analysis and not is_spatial_analysis and ("pathway" in step_name or "enrichment" in step_name):
                    top_pathways = step_info.get("top_pathways", [])
                    if top_pathways:
                        # Extract pathway NAMES only (for biological interpretation)
                        key_findings["top_pathways"] = [
                            p.get('name', 'Unknown') if isinstance(p, dict) else str(p)
                            for p in top_pathways[:3]  # 🔥 只保留top 3
                        ]
            
            key_findings_json = json.dumps(key_findings, ensure_ascii=False, indent=2)
            logger.info(f"📊 [AnalysisSummary] key_findings_json长度: {len(key_findings_json)}字符")
            
            # 🔥 修复：限制execution_results_text的长度，避免prompt过长
            execution_results_text = ""
            if execution_results and (execution_results.get("csv_files") or execution_results.get("image_files")):
                execution_results_text = "\n**Actual Analysis Results (From Generated Files):**\n"
                
                if execution_results.get("csv_files"):
                    execution_results_text += "\n**CSV Results Files:**\n"
                    # 🔥 限制：只处理前3个CSV文件，避免数据过多
                    for csv_info in execution_results["csv_files"][:3]:
                        execution_results_text += f"\n- **{csv_info['filename']}**: {csv_info['shape']}\n"
                        # 只显示前5个列名
                        columns_preview = ', '.join(csv_info['columns'][:5])
                        if len(csv_info['columns']) > 5:
                            columns_preview += f" ... (共{len(csv_info['columns'])}列)"
                        execution_results_text += f"  Columns: {columns_preview}\n"
                        # 🔥 限制：不包含完整的行数据，只显示统计信息
                        if csv_info.get("statistics"):
                            # 只保留关键统计信息，限制长度
                            stats = csv_info['statistics']
                            compact_stats = {}
                            for key in ['mean', 'median', 'std', 'min', 'max']:
                                if key in stats:
                                    compact_stats[key] = stats[key]
                            if compact_stats:
                                execution_results_text += f"  Statistics: {json.dumps(compact_stats, ensure_ascii=False)}\n"
                
                if execution_results.get("image_files"):
                    execution_results_text += f"\n**Generated Images ({len(execution_results['image_files'])} files):**\n"
                    # 只显示前3个图片
                    for img_info in execution_results["image_files"][:3]:
                        execution_results_text += f"- {img_info['filename']} ({img_info['type']})\n"
            else:
                execution_results_text = "\n**Note**: No generated files found in output directory. Analysis results are based on step summaries only.\n"
            
            logger.info(f"📊 [AnalysisSummary] execution_results_text长度: {len(execution_results_text)}字符")
            
            if is_spatial_analysis:
                critical_instruction_text = "Based on the provided metrics above, interpret the **spatial transcriptomics** (10x Visium) results. Use terminology: Spots, Clusters, Spatial Domains, Gene Expression, Moran's I, Spatially Variable Genes (SVGs). Do NOT mention metabolites, LC-MS, or metabolomics. Generate a structured Markdown report with spatial biology interpretation."
            else:
                critical_instruction_text = "Based on the provided metrics above, interpret the biological significance. Use your internal knowledge base (PubMed/Literature) to explain **WHY** these specific metabolites/pathways might be altered in this context. Generate a structured Markdown report with deep biological interpretation."
            
            # 空间组学使用专用输出结构，禁止出现代谢物/LC-MS 等术语
            if is_spatial_analysis:
                output_structure_section = """
### 1. 统计概览 (Statistical Overview) — 空间转录组
- 定量总结：Spot 数量、基因数量、Leiden 聚类数、PCA 方差解释（PC1/PC2）
- 数据质量与组织结构（H&E）简要评估
- 空间域（Spatial Domains）与整体数据特征

### 2. 空间聚类与基因表达 (Spatial Clusters & Gene Expression)
- **空间聚类**：Leiden 聚类结果与空间分布
- **空间可变基因（SVGs）**：基于 Moran's I 识别的基因及其空间模式
- **基因表达**：讨论关键基因在组织中的空间表达模式（勿使用代谢物、LC-MS 等代谢组学术语）

### 3. 空间生物学解读 (Spatial Biology Interpretation)
- 组织区域与空间域的关系
- 空间自相关（Moran's I）结果的生物学意义
- 与已知空间转录组学文献的关联

### 4. 结论与建议 (Conclusions & Recommendations)
- 主要发现总结（使用 Spots、Clusters、SVGs、Spatial Domains 等术语）
- 后续空间分析或实验验证建议

**禁止**：不要提及代谢物（metabolites）、LC-MS、代谢组学或差异代谢物。仅使用空间转录组学术语。
"""
            else:
                output_structure_section = """
### 1. 统计概览 (Statistical Overview)
- Quantitative summary: PCA separation quality, PC1/PC2 variance explained, differential analysis counts (up/down regulated)
- Data quality assessment based on PCA results
- Overall data characteristics and key statistics

### 2. 关键生物标志物 (Key Biomarkers)
- **VIP代谢物**: Discuss the top VIP metabolites from PLS-DA analysis (names: """ + (', '.join(key_findings.get('top_vip_metabolites', [])[:5]) if key_findings.get('top_vip_metabolites') else 'see data') + """)
- **差异代谢物**: Discuss the top differentially expressed metabolites (names: """ + (', '.join(key_findings.get('top_differential_metabolites', [])[:5]) if key_findings.get('top_differential_metabolites') else 'see data') + """)
- **生物学功能**: Use your internal knowledge base (PubMed/Literature) to explain the potential functions and biological significance of these metabolites
- **标志物潜力**: Discuss the potential of these metabolites as biomarkers

### 3. 通路机制解读 (Pathway Mechanism Interpretation)
- **富集通路**: Deep dive into the enriched pathways (names: """ + (', '.join(key_findings.get('top_pathways', [])[:5]) if key_findings.get('top_pathways') else 'see data') + """)
- **通路功能**: Explain the biological functions of these pathways and their significance in the current research context
- **机制讨论**: Relate findings to potential biological mechanisms, disease processes, or physiological states
- **功能意义**: Discuss what the differentially expressed metabolites mean in terms of biological function

### 4. 结论与建议 (Conclusions & Recommendations)
- **主要发现总结**: Summarize key findings and their biological significance
- **验证实验建议**: Suggest validation experiments (e.g., targeted metabolomics, qPCR validation)
- **后续研究**: Propose follow-up studies based on the findings
"""
            
            prompt = f"""You are a Senior Bioinformatics Scientist writing a Results & Discussion section for a top-tier journal (Nature Medicine). Your role is to interpret biological data and provide deep scientific insights, connecting findings to biological mechanisms and literature knowledge.

**User Goal:**
{workflow_name}

**Execution Results (Successful Steps):**
{summary_json}

**Key Findings Extracted (Specific Metrics):**
{key_findings_json}
{execution_results_text}
{failure_info}

**CRITICAL INSTRUCTION:**
{critical_instruction_text}

**Domain Context:**
{domain_context}

**CRITICAL RULES:**

1. **Reasoning Process (DeepSeek-R1)**: 
   - Use the `<think>` tag to show your reasoning process before generating the final report
   - Inside `<think>`, analyze the data metrics, connect metabolites to pathways, and reason about biological mechanisms
   - After reasoning, output the final report outside the `<think>` tags

2. **Scientific Persona**: You are a Senior Bioinformatics Scientist writing a publication-quality results section for Nature Medicine. Write as if you are describing results in a Methods/Results section of a high-impact research paper.

3. **NO Technical Debugging**: 
   - DO NOT mention step names, tool names, file paths, or technical errors
   - DO NOT say "Step X failed" or "Tool Y encountered an error"
   - DO NOT mention Python errors, missing libraries, or code issues
   - If a step failed, simply state the biological limitation (e.g., "Pathway enrichment analysis could not be performed due to insufficient significant features" or "Functional annotation was not available for this dataset")

4. **Deep Biological Interpretation**:
   - Connect metabolites/pathways to biological functions using your internal knowledge base (PubMed/Literature)
   - Explain the MECHANISM, not just the numbers
   - Discuss how the identified metabolites/pathways relate to biological processes, disease mechanisms, or physiological states
   - Interpret findings in the context of known metabolic pathways and their roles

5. **Professional Language**:
   - Use scientific terminology appropriate for Nature Medicine
   - Write in Simplified Chinese (简体中文)
   - Be precise, detailed, and academically rigorous
   - Minimum 800 words, aim for comprehensive coverage

6. **Output Structure (MUST FOLLOW):**
{output_structure_section}

**Output Format:**
- Use Simplified Chinese (简体中文)
- Use Markdown format with proper headings (###)
- Be professional, academic, and detailed
- Minimum 800 words, aim for comprehensive coverage
- Include specific numbers, percentages, and statistical values from the results
- Reference biological mechanisms and pathways explicitly

**Tone**: Professional, Academic, Detailed, Nature Medicine style. Focus on deep biological interpretation and scientific insights, connecting findings to mechanisms.

**CRITICAL**: You MUST provide a detailed Biological Interpretation and Mechanism Analysis. Do NOT just list steps or metrics. Explain the biological meaning, connect findings to known pathways, and discuss mechanisms.

**IMPORTANT**: Use `<think>` tags to show your reasoning process. Analyze the data deeply, then output the final report.

现在生成全面的分析报告（遵循上述结构，详细且专业，包含深度生物学机制解读）："""
            
            messages = [
                {
                    "role": "system",
                    "content": f"""You are a Senior Bioinformatics Scientist writing a publication-quality results section for a {omics_type} research paper. You are NOT a software engineer, IT support, or debugger.

**Your Scientific Persona:**
- You interpret biological data and provide scientific insights
- You write as if describing results in a Methods/Results section of a research paper
- You focus on biological meaning, statistical significance, and scientific interpretation

**What You MUST DO:**
- Describe sample characteristics, group comparisons, and statistical findings
- Interpret clustering patterns, separation between groups, and biological significance
- Explain what the data reveals about the biological system under study
- Discuss functional implications and potential biological mechanisms
- Use scientific terminology appropriate for {omics_type} research

**What You MUST NOT DO:**
- Do NOT mention step names, tool names, file paths, or technical implementation details
- Do NOT say "Step X failed" or "Tool Y encountered an error"
- Do NOT mention Python errors, missing libraries, code issues, or debugging information
- Do NOT act like IT support or a software engineer
- If a step failed, state it as a biological limitation (e.g., "Pathway enrichment could not be performed due to insufficient significant features")

**Output Style:**
- Write in Simplified Chinese (简体中文)
- Use Markdown format with proper headings
- Be precise, detailed, and academically rigorous
- Focus on biological interpretation and scientific insights"""
                },
                {"role": "user", "content": prompt}
            ]
            
            # 🔥 TASK 2: Force LLM call - ALWAYS call LLM, never return simple list
            logger.info(f"📞 [AnalysisSummary] 调用 LLM 生成深度生物学解释...")
            logger.info(f"📊 [AnalysisSummary] 提取的关键指标: {key_findings_json}")
            logger.info(f"📊 [AnalysisSummary] 成功步骤数: {len(successful_steps)}/{len(steps_results)}")
            logger.info(f"📊 [AnalysisSummary] 失败步骤数: {len(failed_steps)}")
            
            # 🔥 TASK 2: Debug logging - Log metrics being sent to LLM
            if not key_findings_json or key_findings_json == "{}":
                logger.warning(f"⚠️ [AnalysisSummary] 警告：关键指标为空，LLM可能无法生成有意义的报告")
            else:
                logger.info(f"✅ [AnalysisSummary] 关键指标已提取，包含数据，准备发送给LLM")
            
            try:
                # 🔥 TASK 3: 统一LLM客户端获取逻辑（与规划阶段一致）
                llm_client_to_use = self.llm_client
                if not llm_client_to_use:
                    logger.warning("⚠️ [AnalysisSummary] self.llm_client 不可用，使用 LLMClientFactory.create_default()")
                    from gibh_agent.core.llm_client import LLMClientFactory
                    llm_client_to_use = LLMClientFactory.create_default()
                    logger.info(f"✅ [AnalysisSummary] 已创建默认LLM客户端: {llm_client_to_use.base_url}")
                
                # 🔥 修复：检查prompt总长度，避免超过API限制
                system_message_length = len(messages[0]["content"]) if messages else 0
                user_message_length = len(messages[1]["content"]) if len(messages) > 1 else 0
                total_prompt_length = system_message_length + user_message_length
                
                logger.info(f"📊 [AnalysisSummary] Prompt长度检查:")
                logger.info(f"  - System message: {system_message_length}字符")
                logger.info(f"  - User message: {user_message_length}字符")
                logger.info(f"  - 总长度: {total_prompt_length}字符")
                
                # 🔥 修复：如果prompt过长（>100k字符），截断或简化
                MAX_PROMPT_LENGTH = 100000  # 100k字符限制
                if total_prompt_length > MAX_PROMPT_LENGTH:
                    logger.warning(f"⚠️ [AnalysisSummary] Prompt过长（{total_prompt_length}字符），超过限制（{MAX_PROMPT_LENGTH}字符），进行截断")
                    # 截断execution_results_text
                    if len(execution_results_text) > 5000:
                        execution_results_text = execution_results_text[:5000] + "\n... (内容已截断)"
                        logger.warning(f"⚠️ [AnalysisSummary] execution_results_text已截断到5000字符")
                    # 简化summary_json（如果仍然过长）
                    if len(summary_json) > 20000:
                        # 只保留最关键的步骤信息
                        import json
                        compact_summary = {
                            "total_steps": len(steps_results),
                            "successful_steps": len(successful_steps),
                            "failed_steps": len(failed_steps),
                            "key_metrics": {
                                "pca_variance": key_findings.get("pca_variance", {}),
                                "differential_count": key_findings.get("differential_count", "N/A"),
                                "top_pathways": key_findings.get("top_pathways", [])[:3],
                                "top_vip_metabolites": key_findings.get("top_vip_metabolites", [])[:3]
                            }
                        }
                        summary_json = json.dumps(compact_summary, ensure_ascii=False, indent=2)
                        logger.warning(f"⚠️ [AnalysisSummary] summary_json已简化，新长度: {len(summary_json)}字符")
                    
                    # 重新构建prompt（使用简化后的数据）
                    prompt = f"""You are a Senior Bioinformatics Scientist writing a Results & Discussion section for a top-tier journal (Nature Medicine). Your role is to interpret biological data and provide deep scientific insights, connecting findings to biological mechanisms and literature knowledge.

**User Goal:**
{workflow_name}

**Execution Results Summary:**
{summary_json}

**Key Findings Extracted (Specific Metrics):**
{key_findings_json}
{execution_results_text}
{failure_info}

**CRITICAL INSTRUCTION:**
{critical_instruction_text}

**Domain Context:**
{domain_context}

**CRITICAL RULES:**

1. **Reasoning Process (DeepSeek-R1)**: 
   - Use the `<think>` tag to show your reasoning process before generating the final report
   - Inside `<think>`, analyze the data metrics and reason about biological mechanisms
   - After reasoning, output the final report outside the `<think>` tags

2. **Scientific Persona**: You are a Senior Bioinformatics Scientist writing a publication-quality results section for Nature Medicine. Write as if you are describing results in a Methods/Results section of a high-impact research paper.

3. **NO Technical Debugging**: 
   - DO NOT mention step names, tool names, file paths, or technical errors
   - DO NOT say "Step X failed" or "Tool Y encountered an error"
   - DO NOT mention Python errors, missing libraries, or code issues
   - If a step failed, simply state the biological limitation

4. **Deep Biological Interpretation**:
   - Connect findings to biological functions using your internal knowledge base (PubMed/Literature)
   - Explain the MECHANISM, not just the numbers
   - Discuss how the identified results relate to biological processes, disease mechanisms, or physiological states

5. **Professional Language**:
   - Use scientific terminology appropriate for Nature Medicine
   - Write in Simplified Chinese (简体中文)
   - Be precise, detailed, and academically rigorous
   - Minimum 800 words, aim for comprehensive coverage

6. **Output Structure (MUST FOLLOW):**
{output_structure_section}

**Output Format:**
- Use Simplified Chinese (简体中文)
- Use Markdown format with proper headings (###)
- Be professional, academic, and detailed

**Tone**: Professional, Academic, Detailed, Nature Medicine style. Focus on deep biological interpretation and scientific insights, connecting findings to mechanisms.

**CRITICAL**: You MUST provide a detailed Biological Interpretation and Mechanism Analysis. Do NOT just list steps or metrics. Explain the biological meaning and discuss mechanisms.

**IMPORTANT**: Use `<think>` tags to show your reasoning process. Analyze the data deeply, then output the final report.

现在生成全面的分析报告（遵循上述结构，详细且专业，包含深度生物学机制解读）："""
                    
                    messages[1]["content"] = prompt
                    
                    # 重新计算长度
                    system_message_length = len(messages[0]["content"])
                    user_message_length = len(messages[1]["content"])
                    total_prompt_length = system_message_length + user_message_length
                    logger.info(f"📊 [AnalysisSummary] 截断后Prompt长度: {total_prompt_length}字符")
                
                def _clean_report_content(s: str):
                    """Strip <<<SUGGESTIONS>>> block and store in context; return cleaned text."""
                    if not s:
                        return s
                    cleaned, sug = strip_suggestions_from_text(s)
                    if sug:
                        self.context["report_suggestions"] = sug
                    return cleaned

                logger.info(f"📞 [AnalysisSummary] 开始LLM调用，max_tokens=2500...")
                completion = await llm_client_to_use.achat(messages, temperature=0.3, max_tokens=2500)  # 🔥 TASK 2: Increase tokens for comprehensive report
                logger.info(f"✅ [AnalysisSummary] LLM调用完成，开始解析响应...")
                original_content = completion.choices[0].message.content or ""
                logger.info(f"🔍 [AnalysisSummary] 原始内容长度: {len(original_content)}")
                
                think_content, response = llm_client_to_use.extract_think_and_content(completion)
                logger.info(f"🔍 [AnalysisSummary] 内容提取结果: think_length={len(think_content) if think_content else 0}, response_length={len(response) if response else 0}")
                logger.debug(f"🔍 [AnalysisSummary] original_content 预览: {original_content[:300]}...")
                logger.debug(f"🔍 [AnalysisSummary] response 预览: {response[:300] if response else 'N/A'}...")
                
                # 🔥 修复：如果提取后的response太短，但original_content很长，说明主要内容可能在标签内
                # 在这种情况下，应该使用original_content（前端会解析标签）
                if response and len(response.strip()) > 100:  # Ensure meaningful response
                    logger.info(f"✅ [AnalysisSummary] 深度生物学解释生成成功，长度: {len(response)}")
                    logger.debug(f"📝 [DEBUG] Summary preview: {response[:200]}...")
                    # Return original content with tags so frontend can parse and display reasoning
                    has_think_tags = any(tag in original_content for tag in ['<think>', '<think>', '<reasoning>', '<thought>', '<thinking>'])
                    out = original_content if has_think_tags else response
                    return _clean_report_content(out)
                elif original_content and len(original_content.strip()) > 100:
                    # 🔥 修复：如果response太短但original_content很长，使用original_content
                    logger.warning(f"⚠️ [AnalysisSummary] 提取后的内容过短，但原始内容较长，使用原始内容（长度: {len(original_content)}）")
                    return _clean_report_content(original_content)
                else:
                    logger.warning(f"⚠️ [AnalysisSummary] LLM 返回内容过短（response: {len(response) if response else 0}字符, original: {len(original_content)}字符），尝试重新生成...")
                    # Retry with simpler prompt if first attempt failed
                    retry_prompt = f"""Based on these analysis metrics: {key_findings_json}

Write a comprehensive biological interpretation report in Simplified Chinese. Include:
1. 结果摘要 (quantitative findings)
2. 生物学机制解读 (connect metabolites/pathways to biological functions)
3. 潜在标志物 (discuss VIP molecules)
4. 下一步建议 (validation experiments)

Minimum 500 words. Be scientific and detailed."""
                    
                    retry_messages = [
                        {"role": "system", "content": "You are a Senior Bioinformatics Scientist. Write detailed biological interpretations."},
                        {"role": "user", "content": retry_prompt}
                    ]
                    
                    retry_completion = await llm_client_to_use.achat(retry_messages, temperature=0.3, max_tokens=2000)
                    retry_think, retry_response = llm_client_to_use.extract_think_and_content(retry_completion)
                    
                    # 🔥 FEATURE: Return original content with tags for frontend parsing
                    retry_original_content = retry_completion.choices[0].message.content or ""
                    logger.info(f"🔍 [AnalysisSummary] 重试内容提取结果: think_length={len(retry_think) if retry_think else 0}, response_length={len(retry_response) if retry_response else 0}, original_length={len(retry_original_content)}")
                    
                    if retry_response and len(retry_response.strip()) > 100:
                        logger.info(f"✅ [AnalysisSummary] 重试成功，生成深度解释，长度: {len(retry_response)}")
                        # Return original content with tags so frontend can parse and display reasoning
                        has_think_tags = any(tag in retry_original_content for tag in ['<think>', '<think>', '<reasoning>', '<thought>', '<thinking>', '<think>'])
                        out = retry_original_content if has_think_tags else retry_response
                        return _clean_report_content(out)
                    elif retry_original_content and len(retry_original_content.strip()) > 100:
                        # 🔥 修复：如果重试后的response太短但original_content很长，使用original_content
                        logger.warning(f"⚠️ [AnalysisSummary] 重试后提取的内容过短，但原始内容较长，使用原始内容（长度: {len(retry_original_content)}）")
                        return _clean_report_content(retry_original_content)
                    else:
                        logger.error(f"❌ [AnalysisSummary] 重试后仍无法生成有效内容（response: {len(retry_response) if retry_response else 0}字符, original: {len(retry_original_content)}字符）")
                        # 🔥 TASK 3: Return user-friendly error message instead of raw traceback
                        return f"""## ⚠️ 分析报告生成失败


**已完成的步骤**: {len(successful_steps)}/{len(steps_results)}

**关键指标**:
{key_findings_json if key_findings_json != "{}" else "暂无可用指标"}

**建议**: 请查看上方的详细图表和统计结果以获取分析信息。"""
            except Exception as llm_error:
                # 🔥 TASK 2: 增强错误日志输出
                import traceback
                error_type = type(llm_error).__name__
                error_msg = (
                    f"❌ [AnalysisSummary] LLM 调用失败\n"
                    f"━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"
                    f"错误类型: {error_type}\n"
                    f"错误信息: {str(llm_error)}\n"
                    f"━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"
                    f"调用上下文:\n"
                    f"  - LLM Client base_url: {llm_client_to_use.base_url if hasattr(llm_client_to_use, 'base_url') else 'N/A'}\n"
                    f"  - LLM Client model: {llm_client_to_use.model if hasattr(llm_client_to_use, 'model') else 'N/A'}\n"
                    f"  - API Key: {'已设置' if hasattr(llm_client_to_use, 'api_key') and llm_client_to_use.api_key else '未设置'}\n"
                    f"━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"
                    f"调用参数:\n"
                    f"  - messages数量: {len(messages)}\n"
                    f"  - system message长度: {len(messages[0]['content']) if messages else 0} 字符\n"
                    f"  - user message长度: {len(messages[1]['content']) if len(messages) > 1 else 0} 字符\n"
                    f"  - temperature: 0.3\n"
                    f"  - max_tokens: 2500\n"
                    f"━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"
                    f"执行上下文:\n"
                    f"  - 成功步骤数: {len(successful_steps)}/{len(steps_results)}\n"
                    f"  - 失败步骤数: {len(failed_steps)}\n"
                    f"  - 关键指标: {key_findings_json[:200] if key_findings_json else 'N/A'}...\n"
                    f"━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"
                    f"可能原因:\n"
                    f"  - API密钥无效或过期\n"
                    f"  - 网络连接问题\n"
                    f"  - API服务暂时不可用\n"
                    f"  - 请求超时\n"
                    f"  - 请求内容过长\n"
                    f"━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"
                    f"完整堆栈:\n{traceback.format_exc()}\n"
                    f"━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
                )
                logger.error(error_msg)
                
                # 🔥 TASK: 将详细错误信息存储到context，供orchestrator通过SSE发送到前端
                self.context["last_llm_error"] = {
                    "error_type": error_type,
                    "error_message": str(llm_error),
                    "error_details": error_msg,
                    "context": {
                        "llm_client_base_url": llm_client_to_use.base_url if hasattr(llm_client_to_use, 'base_url') else 'N/A',
                        "llm_client_model": llm_client_to_use.model if hasattr(llm_client_to_use, 'model') else 'N/A',
                        "api_key_set": bool(hasattr(llm_client_to_use, 'api_key') and llm_client_to_use.api_key),
                        "messages_count": len(messages),
                        "temperature": 0.3,
                        "max_tokens": 2500,
                        "successful_steps": len(successful_steps),
                        "failed_steps": len(failed_steps),
                        "key_findings": key_findings_json[:200] if key_findings_json else 'N/A'
                    },
                    "possible_causes": [
                        "API密钥无效或过期",
                        "网络连接问题",
                        "API服务暂时不可用",
                        "请求超时",
                        "请求内容过长"
                    ]
                }
                
                # 🔥 修复：去除降级处理，LLM失败时返回用户友好的错误信息
                logger.error(f"❌ [AnalysisSummary] LLM调用失败，返回用户友好的错误信息")
                
                # 返回用户友好的错误信息（面向用户，非技术性）
                user_friendly_error = f"""## ⚠️ AI专家分析报告生成失败

很抱歉，AI专家分析报告生成时遇到了技术问题，请稍后再试。

**分析执行情况**:
- 已完成的步骤: {len(successful_steps)}/{len(steps_results)}
- 失败的步骤: {len(failed_steps)}

**建议**:
- 请查看上方的详细图表和统计结果获取分析信息
- 如果问题持续存在，请联系技术支持

**错误信息**: {error_type} - {str(llm_error)[:100]}{'...' if len(str(llm_error)) > 100 else ''}
"""
                
                logger.info(f"✅ [AnalysisSummary] 已返回用户友好的错误信息，长度: {len(user_friendly_error)}")
                return user_friendly_error
                
        except Exception as e:
            logger.error(f"❌ [AnalysisSummary] 生成分析摘要失败: {e}", exc_info=True)
            return None
    
    async def _evaluate_analysis_quality(
        self,
        steps_results: List[Dict[str, Any]],
        diagnosis: str,
        workflow_name: str = "Unknown"
    ) -> Dict[str, Any]:
        """
        🔥 PHASE 2: Evaluate analysis quality using LLM
        
        Args:
            steps_results: List of step execution results
            diagnosis: Generated diagnosis report
            workflow_name: Name of the workflow
        
        Returns:
            Dictionary with score (0-100) and critique
        """
        try:
            # Count successful/failed/warning steps
            successful_steps = [s for s in steps_results if s.get("status") == "success"]
            failed_steps = [s for s in steps_results if s.get("status") == "error"]
            warning_steps = [s for s in steps_results if s.get("status") == "warning"]
            
            # Extract key metrics
            metrics = {
                "total_steps": len(steps_results),
                "successful_steps": len(successful_steps),
                "failed_steps": len(failed_steps),
                "warning_steps": len(warning_steps),
                "completion_rate": len(successful_steps) / len(steps_results) * 100 if steps_results else 0
            }
            
            # Check for key analysis outputs
            has_pca = any("pca" in str(s.get("step_name", "")).lower() for s in successful_steps)
            has_diff = any("differential" in str(s.get("step_name", "")).lower() for s in successful_steps)
            has_visualization = any("visualize" in str(s.get("step_name", "")).lower() or "plot" in str(s.get("data", {})).lower() for s in successful_steps)
            has_pathway = any("pathway" in str(s.get("step_name", "")).lower() or "enrichment" in str(s.get("step_name", "")).lower() for s in successful_steps)
            
            # Build evaluation prompt
            evaluation_prompt = f"""You are a Senior Bioinformatics Quality Assurance Expert. Evaluate the quality of this {workflow_name} analysis.

**Execution Metrics:**
- Total Steps: {metrics['total_steps']}
- Successful: {metrics['successful_steps']}
- Failed: {metrics['failed_steps']}
- Warnings: {metrics['warning_steps']}
- Completion Rate: {metrics['completion_rate']:.1f}%

**Analysis Components:**
- PCA Analysis: {'✅ Present' if has_pca else '❌ Missing'}
- Differential Analysis: {'✅ Present' if has_diff else '❌ Missing'}
- Visualization: {'✅ Present' if has_visualization else '❌ Missing'}
- Pathway Enrichment: {'✅ Present' if has_pathway else '❌ Missing'}

**Generated Diagnosis Report:**
{diagnosis[:1000]}...

**Evaluation Criteria:**
1. **Completeness (0-30 points)**: Did all critical steps complete? Are key analyses present?
2. **Data Quality (0-25 points)**: Were data quality issues handled? Missing values? Outliers?
3. **Statistical Rigor (0-25 points)**: Were appropriate statistical methods used? Are results significant?
4. **Biological Interpretation (0-20 points)**: Is the diagnosis report scientifically sound? Does it provide biological insights?

**Output Format (JSON only, no markdown):**
{{
    "score": <integer 0-100>,
    "critique": "<brief critique in Simplified Chinese, 2-3 sentences>",
    "strengths": ["<strength 1>", "<strength 2>"],
    "weaknesses": ["<weakness 1>", "<weakness 2>"],
    "recommendations": ["<recommendation 1>", "<recommendation 2>"]
}}

Evaluate and return ONLY the JSON object:"""
            
            messages = [
                {
                    "role": "system",
                    "content": "You are a Senior Bioinformatics Quality Assurance Expert. Evaluate analysis quality and provide constructive feedback. Output ONLY valid JSON, no markdown, no explanations."
                },
                {"role": "user", "content": evaluation_prompt}
            ]
            
            logger.info(f"📞 [QualityEvaluation] 调用 LLM 评估分析质量...")
            completion = await self.llm_client.achat(messages, temperature=0.2, max_tokens=500)
            think_content, response = self.llm_client.extract_think_and_content(completion)
            
            if response:
                # Try to parse JSON from response
                import json
                import re
                
                # Extract JSON from response (handle markdown code blocks)
                json_match = re.search(r'\{[^{}]*"score"[^{}]*\}', response, re.DOTALL)
                if json_match:
                    try:
                        evaluation = json.loads(json_match.group())
                        logger.info(f"✅ [QualityEvaluation] 质量评估完成，得分: {evaluation.get('score', 'N/A')}")
                        return evaluation
                    except json.JSONDecodeError:
                        logger.warning(f"⚠️ [QualityEvaluation] JSON 解析失败，使用默认评估")
                
                # Fallback: generate basic evaluation
                base_score = int(metrics['completion_rate'])
                if has_pca and has_diff:
                    base_score += 10
                if has_visualization:
                    base_score += 5
                if has_pathway:
                    base_score += 5
                
                return {
                    "score": min(100, base_score),
                    "critique": f"分析完成率 {metrics['completion_rate']:.1f}%，{'包含关键分析步骤' if has_pca and has_diff else '缺少部分关键分析'}。",
                    "strengths": ["执行了主要分析步骤"] if has_pca or has_diff else [],
                    "weaknesses": ["部分步骤未完成"] if metrics['failed_steps'] > 0 else [],
                    "recommendations": ["建议检查失败步骤"] if metrics['failed_steps'] > 0 else []
                }
            else:
                logger.warning(f"⚠️ [QualityEvaluation] LLM 响应为空，使用默认评估")
                return {
                    "score": int(metrics['completion_rate']),
                    "critique": "无法生成详细评估",
                    "strengths": [],
                    "weaknesses": [],
                    "recommendations": []
                }
                
        except Exception as e:
            logger.error(f"❌ [QualityEvaluation] 质量评估失败: {e}", exc_info=True)
            # Return default evaluation
            return {
                "score": 50,
                "critique": f"评估过程出错: {str(e)}",
                "strengths": [],
                "weaknesses": [],
                "recommendations": []
            }
    
    def _extract_parameter_recommendations(
        self,
        diagnosis_report: str,
        omics_type: str,
        stats: Dict[str, Any]
    ) -> Optional[Dict[str, Any]]:
        """
        从诊断报告中提取参数推荐
        
        🔥 TASK 5: 解析 Markdown 表格中的参数推荐
        
        Args:
            diagnosis_report: 诊断报告 Markdown 文本
            omics_type: 组学类型
            stats: 统计数据
        
        Returns:
            参数推荐字典，格式：
            {
                "summary": "推荐摘要",
                "params": {
                    "param_name": {
                        "value": "推荐值",
                        "reason": "推荐理由"
                    }
                }
            }
        """
        if not diagnosis_report:
            return None
        
        try:
            import re
            
            # 查找参数推荐表格（Markdown 表格格式）
            # 表格格式：| 参数名 | 默认值 | **推荐值** | 推荐理由 |
            table_pattern = r'###\s*💡\s*参数推荐.*?\n(.*?)(?=\n###|\n##|$)'
            table_match = re.search(table_pattern, diagnosis_report, re.DOTALL | re.IGNORECASE)
            
            if not table_match:
                logger.debug("⚠️ [ParameterRecommendation] 未找到参数推荐表格")
                return None
            
            table_content = table_match.group(1)
            
            # 解析表格行（跳过表头）
            lines = table_content.strip().split('\n')
            params = {}
            
            for line in lines:
                line = line.strip()
                if not line or not line.startswith('|'):
                    continue
                
                # 跳过表头分隔行（如 | :--- | :--- | :--- | :--- |）
                if re.match(r'^\|[\s:---]+\|', line):
                    continue
                
                # 解析表格行：| 参数名 | 默认值 | **推荐值** | 推荐理由 |
                cells = [cell.strip() for cell in line.split('|')[1:-1]]  # 去掉首尾空元素
                
                if len(cells) >= 4:
                    param_name = cells[0].strip()
                    default_value = cells[1].strip()
                    recommended_value = cells[2].strip()
                    reason = cells[3].strip()
                    
                    # 清理推荐值（移除 Markdown 加粗标记）
                    recommended_value = re.sub(r'\*\*|\*', '', recommended_value).strip()
                    
                    # 尝试转换推荐值为合适的类型
                    try:
                        # 尝试转换为数字
                        if '.' in recommended_value:
                            recommended_value = float(recommended_value)
                        else:
                            recommended_value = int(recommended_value)
                    except ValueError:
                        # 保持字符串
                        pass
                    
                    params[param_name] = {
                        "default": default_value,
                        "value": recommended_value,
                        "reason": reason
                    }
            
            if not params:
                logger.debug("⚠️ [ParameterRecommendation] 表格解析成功但未找到参数")
                return None
            
            # 生成推荐摘要
            summary = f"基于数据特征，AI 已为您推荐 {len(params)} 个参数的优化值。"
            
            recommendation = {
                "summary": summary,
                "params": params
            }
            
            logger.info(f"✅ [ParameterRecommendation] 成功提取 {len(params)} 个参数推荐")
            return recommendation
            
        except Exception as e:
            logger.warning(f"⚠️ [ParameterRecommendation] 提取参数推荐失败: {e}", exc_info=True)
            return None

