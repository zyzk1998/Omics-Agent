"""
基础工作流抽象类

定义所有工作流必须实现的接口和通用依赖解析逻辑。

完全体扩展法则 (Future-Proofing Manifesto)
-----------------------------------------
未来扩展任何新模态（如基因组、蛋白组），必须严格遵循「完全体原则」：
1. 新增底层 Tool -> 2. 注册 DAG -> 3. 更新 Agent Prompt 认知 ->
4. 补全前端 Step 映射 -> 5. 更新 UI 引导提示词 -> 6. 纳入全局容错循环。
缺一不可。详见各 workflow 与 planner 中的步骤语义与 STEP_NAME_MAP。
"""
import logging
from abc import ABC, abstractmethod
from typing import Dict, Any, List, Optional, Set
from collections import deque

logger = logging.getLogger(__name__)


class BaseWorkflow(ABC):
    """
    工作流抽象基类
    
    所有具体工作流（MetabolomicsWorkflow, RNAWorkflow）都必须继承此类。
    提供依赖解析、模板生成等核心功能。
    """
    
    def __init__(self):
        """初始化工作流"""
        self.name = self.get_name()
        self.description = self.get_description()
        self.steps_dag = self.get_steps_dag()
    
    @abstractmethod
    def get_name(self) -> str:
        """
        获取工作流名称（用于注册表路由）
        
        Returns:
            工作流名称，如 "Metabolomics" 或 "RNA"
        """
        pass
    
    @abstractmethod
    def get_description(self) -> str:
        """
        获取工作流描述
        
        Returns:
            工作流描述文本
        """
        pass
    
    @abstractmethod
    def get_steps_dag(self) -> Dict[str, List[str]]:
        """
        获取步骤依赖图（DAG）
        
        Returns:
            字典，键为步骤ID，值为该步骤依赖的步骤ID列表
            例如: {"pca": ["preprocess"], "preprocess": ["inspect"]}
            
        Note:
            - 如果步骤没有依赖，使用空列表 []
            - 依赖关系应该是传递闭包（transitive closure）
        """
        pass
    
    @abstractmethod
    def get_step_metadata(self, step_id: str) -> Dict[str, Any]:
        """
        获取步骤的元数据（名称、描述、默认参数等）
        
        Args:
            step_id: 步骤ID
            
        Returns:
            包含步骤元数据的字典，格式：
            {
                "name": "步骤显示名称",
                "description": "步骤描述",
                "tool_id": "工具ID",
                "default_params": {...}  # 可选
            }
        """
        pass
    
    def resolve_dependencies(self, target_steps: List[str]) -> List[str]:
        """
        解析依赖关系，返回完整的步骤链
        
        使用拓扑排序算法，确保依赖步骤在目标步骤之前。
        
        Args:
            target_steps: 用户请求的步骤ID列表
            
        Returns:
            完整的步骤列表（按依赖顺序排序）
            
        Example:
            DAG: {"pca": ["preprocess"], "preprocess": ["inspect"]}
            target_steps: ["pca"]
            Returns: ["inspect", "preprocess", "pca"]
        """
        if not target_steps:
            return []
        
        # 验证所有目标步骤都存在
        all_steps = set(self.steps_dag.keys())
        invalid_steps = set(target_steps) - all_steps
        if invalid_steps:
            logger.warning(f"⚠️ 无效的步骤ID: {invalid_steps}")
            target_steps = [s for s in target_steps if s in all_steps]
        
        if not target_steps:
            return []
        
        # 使用 BFS 收集所有依赖步骤
        visited: Set[str] = set()
        queue = deque(target_steps)
        
        # 先添加所有目标步骤
        for step in target_steps:
            visited.add(step)
        
        # BFS 遍历依赖图
        while queue:
            current_step = queue.popleft()
            dependencies = self.steps_dag.get(current_step, [])
            
            for dep in dependencies:
                if dep not in visited:
                    visited.add(dep)
                    queue.append(dep)
        
        # 拓扑排序：确保依赖步骤在目标步骤之前
        result = self._topological_sort(list(visited))
        
        logger.info(f"✅ 依赖解析完成: {target_steps} -> {result}")
        return result
    
    def _topological_sort(self, steps: List[str]) -> List[str]:
        """
        拓扑排序：确保依赖步骤在目标步骤之前
        
        Args:
            steps: 需要排序的步骤列表
            
        Returns:
            排序后的步骤列表
        """
        # 计算入度
        in_degree: Dict[str, int] = {step: 0 for step in steps}
        for step in steps:
            dependencies = self.steps_dag.get(step, [])
            for dep in dependencies:
                if dep in in_degree:
                    in_degree[step] += 1
        
        # Kahn 算法
        queue = deque([step for step in steps if in_degree[step] == 0])
        result = []
        
        while queue:
            current = queue.popleft()
            result.append(current)
            
            # 更新依赖此步骤的其他步骤的入度
            for step in steps:
                dependencies = self.steps_dag.get(step, [])
                if current in dependencies:
                    in_degree[step] -= 1
                    if in_degree[step] == 0:
                        queue.append(step)
        
        # 如果还有步骤未处理，说明存在循环依赖（不应该发生）
        if len(result) != len(steps):
            logger.warning(f"⚠️ 拓扑排序不完整，可能存在循环依赖。已排序: {result}, 剩余: {set(steps) - set(result)}")
            # 将未排序的步骤追加到末尾
            result.extend([s for s in steps if s not in result])
        
        return result
    
    def generate_template(
        self,
        target_steps: Optional[List[str]] = None,
        file_metadata: Optional[Dict[str, Any]] = None
    ) -> Dict[str, Any]:
        """
        生成工作流模板
        
        Args:
            target_steps: 用户请求的步骤列表（如果为 None，返回完整工作流）
            file_metadata: 文件元数据（可选，用于填充参数）
            
        Returns:
            符合前端格式的工作流配置字典
        """
        # 如果没有指定目标步骤，返回完整工作流
        if target_steps is None:
            target_steps = list(self.steps_dag.keys())
        
        # 解析依赖
        resolved_steps = self.resolve_dependencies(target_steps)
        
        # 🔥 CRITICAL FIX: 确保 resolved_steps 不为空
        if not resolved_steps:
            logger.warning(f"⚠️ [BaseWorkflow] resolve_dependencies 返回空列表，使用完整工作流")
            resolved_steps = list(self.steps_dag.keys())
        
        # 🔥 CRITICAL FIX: 再次确保不为空（如果 DAG 为空，至少返回一个占位步骤）
        if not resolved_steps:
            logger.error(f"❌ [BaseWorkflow] steps_dag 为空，无法生成模板")
            raise ValueError("工作流 DAG 为空，无法生成模板")
        
        # 生成步骤配置
        steps = []
        for step_id in resolved_steps:
            step_meta = self.get_step_metadata(step_id)
            
            # 构建步骤配置
            step_config = {
                "id": step_id,
                "step_id": step_id,
                "tool_id": step_meta.get("tool_id", step_id),
                "name": step_meta.get("name", step_id),
                "step_name": step_meta.get("name", step_id),
                "description": step_meta.get("description", ""),
                "desc": step_meta.get("description", "")[:100],
                "selected": True,
                "params": step_meta.get("default_params", {})
            }
            
            # 🔥 ARCHITECTURAL REFACTOR: Plan-First Support
            # 如果提供了文件元数据，填充 file_path 参数；否则使用占位符
            if file_metadata:
                file_path = file_metadata.get("file_path")
                if file_path:
                    # 检查步骤是否需要 file_path 参数
                    if "file_path" in step_config["params"] or "adata_path" in step_config["params"]:
                        param_name = "adata_path" if "adata_path" in step_config["params"] else "file_path"
                        step_config["params"][param_name] = file_path
                    elif not step_config["params"]:
                        # 如果没有默认参数，添加 file_path
                        step_config["params"]["file_path"] = file_path
            else:
                # 🔥 Plan-First: 如果没有文件，使用占位符
                if "file_path" in step_config["params"] or "adata_path" in step_config["params"]:
                    param_name = "adata_path" if "adata_path" in step_config["params"] else "file_path"
                    step_config["params"][param_name] = "<PENDING_UPLOAD>"
                elif step_id in ["metabo_data_validation", "inspect_data", "preprocess_data", "pca_analysis", "differential_analysis",
                                 "metabo_model_comparison", "metabolomics_plsda", "metabolomics_pathway_enrichment"]:
                    # 对于需要文件路径的步骤，添加占位符
                    step_config["params"]["file_path"] = "<PENDING_UPLOAD>"
            
            steps.append(step_config)
        
        # 构建工作流配置
        workflow_name = self._generate_workflow_name(target_steps, file_metadata)
        
        return {
            "type": "workflow_config",
            "workflow_data": {
                "workflow_name": workflow_name,
                "name": workflow_name,
                "steps": steps
            },
            "file_paths": [file_metadata.get("file_path")] if file_metadata and file_metadata.get("file_path") else []
        }
    
    def _generate_workflow_name(
        self,
        target_steps: List[str],
        file_metadata: Optional[Dict[str, Any]]
    ) -> str:
        """
        生成工作流名称
        
        Args:
            target_steps: 目标步骤列表
            file_metadata: 文件元数据
            
        Returns:
            工作流名称
        """
        # 默认使用工作流名称
        base_name = self.name
        
        # 如果只选择了部分步骤，添加说明
        all_steps = set(self.steps_dag.keys())
        if set(target_steps) != all_steps:
            return f"{base_name} 分析流程（部分步骤）"
        
        return f"{base_name} 标准分析流程"

