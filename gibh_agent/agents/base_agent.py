"""
基础智能体抽象类
所有领域智能体都继承此类
"""
from abc import ABC, abstractmethod
from typing import Dict, Any, Optional, List, AsyncIterator
from ..core.llm_client import LLMClient
from ..core.prompt_manager import PromptManager


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
        
        if stream:
            async for chunk in self.llm_client.astream(messages):
                async for content in self.llm_client.get_stream_content([chunk]):
                    yield content
        else:
            completion = await self.llm_client.achat(messages)
            yield self.llm_client.get_content(completion)
    
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
            file_path: 文件路径
        
        Returns:
            文件类型（如 "fastq", "bam", "h5ad"）
        """
        ext = file_path.split('.')[-1].lower()
        
        type_mapping = {
            "fastq": ["fastq", "fq"],
            "bam": ["bam"],
            "h5ad": ["h5ad"],
            "mtx": ["mtx"],
            "vcf": ["vcf"],
            "bed": ["bed"],
            "bw": ["bw", "bigwig"],
            "bam": ["bam", "sam"]
        }
        
        for file_type, extensions in type_mapping.items():
            if ext in extensions:
                return file_type
        
        return "unknown"

