"""
工作流注册表（单例模式）

管理所有有效的工作流，提供路由和查询功能。
"""
import logging
from typing import Dict, Optional, Type, Any
from .base import BaseWorkflow

logger = logging.getLogger(__name__)


class WorkflowRegistry:
    """
    工作流注册表（单例）
    
    职责：
    1. 注册和存储所有有效的工作流
    2. 根据域名（domain_name）路由到对应的工作流
    3. 提供严格的域绑定检查（只支持 Metabolomics 和 RNA）
    """
    
    _instance = None
    _workflows: Dict[str, BaseWorkflow] = {}
    
    def __new__(cls):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
            cls._instance._initialized = False
        return cls._instance
    
    def __init__(self):
        """初始化注册表"""
        if self._initialized:
            return
        
        self._workflows = {}
        self._initialized = True
        
        # 自动注册所有工作流
        self._auto_register()
    
    def _auto_register(self):
        """自动注册所有工作流"""
        from .metabolomics import MetabolomicsWorkflow
        from .rna import RNAWorkflow
        from .spatial_workflow import SpatialWorkflow
        from .radiomics_workflow import RadiomicsWorkflow
        
        # 注册代谢组学工作流
        metabolomics = MetabolomicsWorkflow()
        self.register(metabolomics)
        
        # 注册 RNA 工作流
        rna = RNAWorkflow()
        self.register(rna)
        
        # 注册空间转录组工作流
        spatial = SpatialWorkflow()
        self.register(spatial)
        
        # 注册影像组学工作流
        radiomics = RadiomicsWorkflow()
        self.register(radiomics)
        
        logger.info(f"✅ [WorkflowRegistry] 已注册 {len(self._workflows)} 个工作流: {list(self._workflows.keys())}")
    
    def register(self, workflow: BaseWorkflow):
        """
        注册工作流
        
        Args:
            workflow: 工作流实例
        """
        domain_name = workflow.get_name()
        if domain_name in self._workflows:
            logger.warning(f"⚠️ 工作流 '{domain_name}' 已存在，将被覆盖")
        
        self._workflows[domain_name] = workflow
        logger.info(f"📝 [WorkflowRegistry] 注册工作流: {domain_name}")
    
    def get_workflow(self, domain_name: str) -> Optional[BaseWorkflow]:
        """
        获取工作流实例
        
        Args:
            domain_name: 域名（如 "Metabolomics", "RNA"）
            
        Returns:
            工作流实例，如果不存在返回 None
            
        Note:
            严格域绑定：只支持已注册的域名
        """
        workflow = self._workflows.get(domain_name)
        if workflow is None:
            logger.warning(f"⚠️ [WorkflowRegistry] 未找到工作流: {domain_name}")
            logger.info(f"   可用工作流: {list(self._workflows.keys())}")
        
        return workflow
    
    def is_supported(self, domain_name: str) -> bool:
        """
        检查域名是否受支持
        
        Args:
            domain_name: 域名
            
        Returns:
            True 如果支持，False 否则
        """
        return domain_name in self._workflows
    
    def list_workflows(self) -> Dict[str, str]:
        """
        列出所有已注册的工作流
        
        Returns:
            字典，键为域名，值为描述
        """
        return {
            name: workflow.get_description()
            for name, workflow in self._workflows.items()
        }
    
    def get_unsupported_error(self, domain_name: str) -> Dict[str, Any]:
        """
        生成"不支持域名"的错误响应
        
        Args:
            domain_name: 不支持的域名
            
        Returns:
            错误响应字典
        """
        supported = list(self._workflows.keys())
        return {
            "type": "error",
            "error": f"不支持的域名: {domain_name}",
            "message": f"系统目前只支持以下域名: {', '.join(supported)}。如需添加新域名（如 Proteomics），请联系管理员。",
            "supported_domains": supported
        }


# 全局注册表实例（单例）
registry = WorkflowRegistry()

