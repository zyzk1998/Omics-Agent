"""
基础智能体抽象类
所有领域智能体都继承此类
"""
from abc import ABC, abstractmethod
from typing import Dict, Any, Optional, List, AsyncIterator
import logging
from openai import AuthenticationError, APIError
from ..core.llm_client import LLMClient
from ..core.prompt_manager import PromptManager

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
                # DeepSeek 的 think 过程会以 <think>...</think> 标签形式返回
                has_yielded = False
                try:
                    async for chunk in self.llm_client.astream(messages):
                        if chunk.choices and chunk.choices[0].delta.content:
                            content = chunk.choices[0].delta.content
                            if content:
                                # 直接传递内容，前端会检测和处理 think 标签
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
        从上传文件列表中提取文件路径
        
        核心原则：智能体只处理文件路径（字符串），不处理二进制数据
        
        Args:
            uploaded_files: 文件列表
        
        Returns:
            文件路径列表
        """
        paths = []
        for file_info in uploaded_files:
            if isinstance(file_info, dict):
                path = file_info.get("path") or file_info.get("name")
            else:
                path = getattr(file_info, "path", None) or getattr(file_info, "name", None)
            
            if path:
                paths.append(path)
        
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

