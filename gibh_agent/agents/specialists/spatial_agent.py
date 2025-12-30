"""空间组学智能体（Spatial Omics Agent）"""
from typing import Dict, Any, List
from ..base_agent import BaseAgent
from ...core.llm_client import LLMClient
from ...core.prompt_manager import PromptManager


class SpatialAgent(BaseAgent):
    """空间组学智能体（占位符，待实现）"""
    
    def __init__(self, llm_client: LLMClient, prompt_manager: PromptManager):
        super().__init__(llm_client, prompt_manager, "spatial_expert")
    
    async def process_query(
        self,
        query: str,
        history: List[Dict[str, str]] = None,
        uploaded_files: List[Dict[str, str]] = None,
        **kwargs
    ) -> Dict[str, Any]:
        """处理查询（待实现）"""
        return {
            "type": "chat",
            "response": "Spatial Omics Agent is under development. Coming soon!"
        }

