"""
GIBH-AGENT 主入口
整合所有组件：路由智能体 + 领域智能体
"""
import os
import yaml
from typing import Dict, Any, Optional
from .core.llm_client import LLMClient, LLMClientFactory
from .core.prompt_manager import PromptManager, create_default_prompt_manager
from .core.dispatcher import TaskDispatcher, create_dispatcher_from_config
from .agents.router_agent import RouterAgent
from .agents.specialists.rna_agent import RNAAgent
from .agents.specialists.dna_agent import DNAAgent
from .agents.specialists.epigenomics_agent import EpigenomicsAgent
from .agents.specialists.metabolomics_agent import MetabolomicsAgent
from .agents.specialists.proteomics_agent import ProteomicsAgent
from .agents.specialists.spatial_agent import SpatialAgent
from .agents.specialists.imaging_agent import ImagingAgent


class GIBHAgent:
    """
    GIBH-AGENT 主类
    
    整合路由智能体和所有领域智能体
    """
    
    def __init__(self, config_path: str = "config/settings.yaml"):
        """
        初始化 GIBH-AGENT
        
        Args:
            config_path: 配置文件路径
        """
        # 加载配置
        self.config = self._load_config(config_path)
        
        # 初始化核心组件
        self.llm_clients = self._init_llm_clients()
        self.prompt_manager = self._init_prompt_manager()
        self.dispatcher = self._init_dispatcher()
        
        # 初始化路由智能体
        self.router = RouterAgent(
            llm_client=self.llm_clients["logic"],
            prompt_manager=self.prompt_manager
        )
        
        # 初始化领域智能体
        self.agents = self._init_domain_agents()
    
    def _load_config(self, config_path: str) -> Dict[str, Any]:
        """加载配置文件"""
        if os.path.exists(config_path):
            with open(config_path, 'r') as f:
                return yaml.safe_load(f)
        return {}
    
    def _init_llm_clients(self) -> Dict[str, LLMClient]:
        """初始化 LLM 客户端"""
        llm_config = self.config.get("llm", {})
        default_type = llm_config.get("default", "local")
        
        clients = {}
        
        if default_type == "local":
            local_config = llm_config.get("local", {})
            clients["logic"] = LLMClientFactory.create_from_config(local_config.get("logic", {}))
            clients["vision"] = LLMClientFactory.create_from_config(local_config.get("vision", {}))
        else:
            cloud_config = llm_config.get("cloud", {})
            # 默认使用 DeepSeek
            clients["logic"] = LLMClientFactory.create_cloud_deepseek()
            clients["vision"] = LLMClientFactory.create_cloud_deepseek()
        
        return clients
    
    def _init_prompt_manager(self) -> PromptManager:
        """初始化提示管理器"""
        template_dir = os.path.join(os.path.dirname(__file__), "..", "config", "prompts")
        if os.path.exists(template_dir):
            return PromptManager(template_dir)
        return create_default_prompt_manager()
    
    def _init_dispatcher(self) -> Optional[TaskDispatcher]:
        """初始化任务分发器"""
        dispatcher_config = self.config.get("dispatcher", {})
        if dispatcher_config:
            return TaskDispatcher(dispatcher_config)
        return None
    
    def _init_domain_agents(self) -> Dict[str, BaseAgent]:
        """初始化所有领域智能体"""
        agents = {}
        
        # RNA Agent（转录组）
        agents["rna_agent"] = RNAAgent(
            llm_client=self.llm_clients["logic"],
            prompt_manager=self.prompt_manager,
            dispatcher=self.dispatcher,
            cellranger_config=self.config.get("tools", {}).get("cellranger", {}),
            scanpy_config=self.config.get("tools", {}).get("scanpy", {})
        )
        
        # DNA Agent（基因组）
        agents["dna_agent"] = DNAAgent(
            llm_client=self.llm_clients["logic"],
            prompt_manager=self.prompt_manager
        )
        
        # 其他智能体（占位符）
        agents["epigenomics_agent"] = EpigenomicsAgent(
            llm_client=self.llm_clients["logic"],
            prompt_manager=self.prompt_manager
        )
        
        agents["metabolomics_agent"] = MetabolomicsAgent(
            llm_client=self.llm_clients["logic"],
            prompt_manager=self.prompt_manager
        )
        
        agents["proteomics_agent"] = ProteomicsAgent(
            llm_client=self.llm_clients["logic"],
            prompt_manager=self.prompt_manager
        )
        
        agents["spatial_agent"] = SpatialAgent(
            llm_client=self.llm_clients["logic"],
            prompt_manager=self.prompt_manager
        )
        
        agents["imaging_agent"] = ImagingAgent(
            llm_client=self.llm_clients["logic"],
            prompt_manager=self.prompt_manager
        )
        
        return agents
    
    async def process_query(
        self,
        query: str,
        history: list = None,
        uploaded_files: list = None
    ):
        """
        处理用户查询（主入口）
        
        Args:
            query: 用户查询文本
            history: 对话历史
            uploaded_files: 上传的文件列表
        
        Returns:
            处理结果（可能是字典或异步生成器）
        """
        # 1. 路由到对应的领域智能体
        route_result = await self.router.process_query(query, history, uploaded_files)
        
        # 2. 获取目标智能体
        routing = route_result.get("routing", "rna_agent")
        target_agent = self.agents.get(routing, self.agents["rna_agent"])
        
        # 3. 处理查询
        result = await target_agent.process_query(query, history, uploaded_files)
        
        # 4. 添加路由信息
        result["routing_info"] = route_result
        
        return result


# 便捷函数
def create_agent(config_path: str = "config/settings.yaml") -> GIBHAgent:
    """创建 GIBH-AGENT 实例"""
    return GIBHAgent(config_path)


if __name__ == "__main__":
    # 示例使用
    import asyncio
    
    async def main():
        agent = create_agent()
        
        # 测试路由
        result = await agent.process_query(
            query="帮我分析一下这个单细胞数据",
            uploaded_files=[{"name": "sample.h5ad", "path": "/data/sample.h5ad"}]
        )
        
        print("Routing result:", result.get("routing_info"))
        print("Agent response:", result)
    
    asyncio.run(main())

