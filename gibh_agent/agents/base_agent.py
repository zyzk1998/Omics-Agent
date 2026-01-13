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
        
        Args:
            file_metadata: FileInspector 返回的文件元数据
            omics_type: 组学类型（"scRNA", "Metabolomics", "BulkRNA", "default"）
            dataframe: 可选的数据预览（DataFrame 或 AnnData）
            system_instruction: 领域特定的系统指令（由各个 Agent 提供）
        
        Returns:
            Markdown 格式的诊断报告，如果失败返回 None
        """
        try:
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
            logger.info(f"✅ [DataDiagnostician] 统计计算完成: {len(stats)} 个指标")
            
            # Step 2: 构建 LLM Prompt
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
            
            # 使用 PromptManager 获取诊断模板
            try:
                # 🔥 确保只传递字符串给模板，不传递字典
                prompt = self.prompt_manager.get_prompt(
                    "data_diagnosis",
                    {
                        "inspection_data": stats_json,  # 字符串
                        "head_preview": head_preview[:500] if head_preview else ""  # 字符串，截断到 500 字符
                    },
                    fallback=DATA_DIAGNOSIS_PROMPT.format(inspection_data=stats_json)
                )
                logger.debug(f"📝 [DEBUG] Prompt length: {len(prompt)}")
            except Exception as prompt_err:
                logger.warning(f"⚠️ 获取诊断模板失败，使用默认模板: {prompt_err}")
                try:
                    # 🔥 安全地格式化 prompt，避免 format 错误
                    # 确保 stats_json 是字符串
                    if not isinstance(stats_json, str):
                        stats_json = json.dumps(stats_json, ensure_ascii=False)
                    prompt = DATA_DIAGNOSIS_PROMPT.format(inspection_data=stats_json)
                except Exception as format_err:
                    logger.error(f"❌ [DataDiagnostician] Prompt 格式化失败: {format_err}")
                    # 使用简单的 prompt
                    # 🔥 确保 stats_json 是字符串
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
            # 🔥 架构重构：使用策略模式，从 Agent 传入 system_instruction
            if system_instruction:
                # 使用 Agent 提供的领域特定指令
                system_prompt = system_instruction
                logger.debug(f"✅ [DataDiagnostician] Using domain-specific system instruction (length: {len(system_instruction)})")
            else:
                # 回退到通用指令（向后兼容）
                logger.warning(f"⚠️ [DataDiagnostician] No system_instruction provided, using generic prompt")
                system_prompt = "You are a Senior Bioinformatician. Generate data diagnosis and parameter recommendations in Simplified Chinese."
            
            # 🔥 架构重构：将 system_instruction 前置到用户 prompt（确保上下文隔离）
            if system_instruction:
                # 在用户 prompt 前添加系统指令，确保 LLM 理解领域约束
                prompt = f"""{system_instruction}

{prompt}"""
            
            messages = [
                {"role": "system", "content": system_prompt},
                {"role": "user", "content": prompt}
            ]
            
            # 🔥 Step 3: 调用 LLM 生成 Markdown 报告
            # 🔥 CRITICAL DEBUGGING: 包装在详细的 try-except 中
            try:
                logger.info(f"📞 [DataDiagnostician] 调用 LLM 生成报告...")
                logger.debug(f"📝 [DEBUG] LLM Client type: {type(self.llm_client)}")
                logger.debug(f"📝 [DEBUG] LLM Client methods: {dir(self.llm_client)}")
                
                completion = await self.llm_client.achat(messages, temperature=0.3, max_tokens=1500)
                
                logger.debug(f"📝 [DEBUG] LLM completion type: {type(completion)}")
                logger.debug(f"📝 [DEBUG] LLM completion: {completion}")
                
                think_content, response = self.llm_client.extract_think_and_content(completion)
                
                # 🔥 DEBUG: 打印诊断报告信息
                if response:
                    logger.info(f"✅ [DataDiagnostician] 诊断报告生成成功，长度: {len(response)}")
                    logger.debug(f"📝 [DEBUG] Diagnosis report preview: {response[:200]}...")
                else:
                    logger.warning(f"⚠️ [DataDiagnostician] 诊断报告为空")
                    logger.warning(f"⚠️ [DEBUG] Think content: {think_content[:200] if think_content else 'None'}")
                
                # Step 4: 保存到上下文（供 UI 和后续步骤使用）
                self.context["diagnosis_report"] = response
                self.context["diagnosis_stats"] = stats
                
                return response
                
            except AttributeError as attr_err:
                # LLM 客户端方法不存在
                import traceback
                error_msg = (
                    f"LLM 客户端方法调用失败: {str(attr_err)}\n"
                    f"LLM Client type: {type(self.llm_client)}\n"
                    f"Available methods: {[m for m in dir(self.llm_client) if not m.startswith('_')]}\n"
                    f"Stack trace:\n{traceback.format_exc()}"
                )
                logger.error(f"❌ [DataDiagnostician] {error_msg}")
                return f"⚠️ **诊断报告生成失败**\n\nLLM 客户端错误: {str(attr_err)}\n\n请检查服务器日志获取详细信息。"
                
            except Exception as llm_err:
                # LLM 调用失败
                import traceback
                error_msg = (
                    f"LLM 调用失败: {str(llm_err)}\n"
                    f"Error type: {type(llm_err).__name__}\n"
                    f"Stack trace:\n{traceback.format_exc()}"
                )
                logger.error(f"❌ [DataDiagnostician] {error_msg}")
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

