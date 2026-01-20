"""
Omics-Flow Visualizer: RAGFlow-style Workflow Builder

将 Omics-Agent 的工作流转换为可视化节点图，支持：
1. 节点拖拽和连接
2. 参数在线编辑
3. 实时执行监控
4. 执行日志显示
"""

import uuid
import json
from typing import List, Dict, Callable, Any, Optional, Tuple
from pathlib import Path
from dataclasses import dataclass, field
from enum import Enum
import logging

from pydantic import BaseModel, Field

logger = logging.getLogger(__name__)


class NodeStatus(str, Enum):
    """节点状态"""
    IDLE = "idle"
    RUNNING = "running"
    SUCCESS = "success"
    ERROR = "error"
    SKIPPED = "skipped"


class NodeType(str, Enum):
    """节点类型"""
    DATA_LOADER = "DataLoader"
    PREPROCESSOR = "Preprocessor"
    ANALYZER = "Analyzer"
    VISUALIZER = "Visualizer"
    ANNOTATOR = "Annotator"
    GENERIC = "Generic"


@dataclass
class NodeConfig:
    """节点配置（Pydantic 兼容）"""
    name: str
    node_type: str = NodeType.GENERIC
    tool_id: str = ""
    description: str = ""
    params: Dict[str, Any] = field(default_factory=dict)
    inputs: List[str] = field(default_factory=list)  # 输入参数名列表
    outputs: List[str] = field(default_factory=list)  # 输出参数名列表


class WorkflowNode:
    """
    工作流节点（映射到 UI 节点）
    
    对应关系：
    - Tool (e.g., DESeq2, KeggSearch) -> Node
    - Parameters -> Sidebar (可编辑)
    - Execution Status -> Node Status (运行中/完成/报错)
    """
    
    def __init__(self, node_id: str, config: NodeConfig, tool_func: Optional[Callable] = None):
        self.id = node_id
        self.config = config
        self.tool_func = tool_func  # 实际的工具函数
        self.status = NodeStatus.IDLE
        self.outputs: Dict[str, Any] = {}
        self.error: Optional[str] = None
        self.execution_time: float = 0.0
        self.logs: List[str] = []
        
    def add_log(self, message: str):
        """添加执行日志"""
        self.logs.append(message)
        logger.info(f"[Node {self.id}] {message}")
    
    def run(self, context: Dict[str, Any] = None, **kwargs) -> Dict[str, Any]:
        """
        执行节点逻辑
        
        Args:
            context: 前序节点的输出上下文
            **kwargs: 额外的执行参数
            
        Returns:
            节点执行结果
        """
        import time
        start_time = time.time()
        
        self.status = NodeStatus.RUNNING
        self.add_log(f"开始执行: {self.config.name}")
        
        try:
            # 合并参数：config.params + context + kwargs
            execution_params = {}
            execution_params.update(self.config.params)
            if context:
                execution_params.update(context)
            execution_params.update(kwargs)
            
            # 调用工具函数
            if self.tool_func:
                result = self.tool_func(**execution_params)
                self.outputs = result if isinstance(result, dict) else {"result": result}
            else:
                # 如果没有工具函数，返回模拟结果
                self.outputs = {"status": "simulated", "message": f"Node {self.config.name} executed"}
            
            self.status = NodeStatus.SUCCESS
            self.execution_time = time.time() - start_time
            self.add_log(f"执行完成 (耗时: {self.execution_time:.2f}s)")
            
            return self.outputs
            
        except Exception as e:
            self.status = NodeStatus.ERROR
            self.error = str(e)
            self.execution_time = time.time() - start_time
            self.add_log(f"执行失败: {str(e)}")
            logger.error(f"Node {self.id} execution failed: {e}", exc_info=True)
            return {"error": str(e), "status": "error"}


class OmicsGraph:
    """
    Omics 工作流图（映射到 UI 流程图）
    
    对应关系：
    - Analysis Pipeline -> Edges (连线)
    - Steps DAG -> Graph Structure
    """
    
    def __init__(self, workflow_name: str = "Omics Workflow"):
        self.workflow_name = workflow_name
        self.nodes: Dict[str, WorkflowNode] = {}
        self.edges: List[Tuple[str, str]] = []  # (source_id, target_id)
        self.execution_order: List[str] = []  # 拓扑排序后的执行顺序
        
    def add_node(
        self,
        node_id: Optional[str] = None,
        name: str = "",
        node_type: str = NodeType.GENERIC,
        tool_id: str = "",
        description: str = "",
        params: Dict[str, Any] = None,
        tool_func: Optional[Callable] = None
    ) -> str:
        """
        添加节点到图
        
        Args:
            node_id: 节点 ID（如果为 None，自动生成）
            name: 节点名称
            node_type: 节点类型
            tool_id: 工具 ID（用于查找实际函数）
            description: 节点描述
            params: 节点参数
            tool_func: 工具函数（可选）
            
        Returns:
            节点 ID
        """
        if node_id is None:
            node_id = str(uuid.uuid4())[:8]
        
        config = NodeConfig(
            name=name,
            node_type=node_type,
            tool_id=tool_id,
            description=description,
            params=params or {}
        )
        
        node = WorkflowNode(node_id, config, tool_func)
        self.nodes[node_id] = node
        
        logger.info(f"Added node: {node_id} ({name})")
        return node_id
    
    def add_edge(self, source_id: str, target_id: str):
        """
        添加边（连接两个节点）
        
        Args:
            source_id: 源节点 ID
            target_id: 目标节点 ID
        """
        if source_id not in self.nodes:
            raise ValueError(f"Source node {source_id} not found")
        if target_id not in self.nodes:
            raise ValueError(f"Target node {target_id} not found")
        
        self.edges.append((source_id, target_id))
        logger.info(f"Added edge: {source_id} -> {target_id}")
    
    def get_execution_order(self) -> List[str]:
        """
        计算拓扑排序后的执行顺序
        
        Returns:
            节点 ID 列表（按执行顺序）
        """
        if not self.nodes:
            return []
        
        # 构建入度图
        in_degree = {node_id: 0 for node_id in self.nodes}
        for source_id, target_id in self.edges:
            in_degree[target_id] += 1
        
        # 拓扑排序（Kahn's algorithm）
        queue = [node_id for node_id, degree in in_degree.items() if degree == 0]
        execution_order = []
        
        while queue:
            node_id = queue.pop(0)
            execution_order.append(node_id)
            
            # 更新依赖节点的入度
            for source_id, target_id in self.edges:
                if source_id == node_id:
                    in_degree[target_id] -= 1
                    if in_degree[target_id] == 0:
                        queue.append(target_id)
        
        # 检查是否有循环
        if len(execution_order) != len(self.nodes):
            logger.warning("Graph contains cycles, execution order may be incomplete")
        
        self.execution_order = execution_order
        return execution_order
    
    def execute(self, initial_context: Dict[str, Any] = None) -> Dict[str, Any]:
        """
        执行工作流图
        
        Args:
            initial_context: 初始上下文（包含文件路径等）
            
        Returns:
            最终执行结果
        """
        execution_order = self.get_execution_order()
        context = initial_context or {}
        
        logger.info(f"开始执行工作流: {self.workflow_name}")
        logger.info(f"执行顺序: {execution_order}")
        
        for node_id in execution_order:
            node = self.nodes[node_id]
            
            # 执行节点
            result = node.run(context=context)
            
            # 更新上下文（将节点输出传递给后续节点）
            if node.status == NodeStatus.SUCCESS:
                # 使用节点 ID 作为键，避免冲突
                context[f"{node_id}_output"] = result
                # 也更新通用键（如果节点有定义 outputs）
                for output_key in node.config.outputs:
                    if output_key in result:
                        context[output_key] = result[output_key]
            elif node.status == NodeStatus.ERROR:
                logger.error(f"节点 {node_id} 执行失败，停止工作流")
                break
        
        return {
            "workflow_name": self.workflow_name,
            "nodes_executed": len([n for n in self.nodes.values() if n.status != NodeStatus.IDLE]),
            "final_context": context
        }
    
    def to_dict(self) -> Dict[str, Any]:
        """
        转换为字典（用于前端序列化）
        
        Returns:
            图的字典表示
        """
        return {
            "workflow_name": self.workflow_name,
            "nodes": {
                node_id: {
                    "id": node_id,
                    "name": node.config.name,
                    "type": node.config.node_type,
                    "tool_id": node.config.tool_id,
                    "description": node.config.description,
                    "params": node.config.params,
                    "status": node.status.value,
                    "position": {"x": 0, "y": 0}  # 占位符，实际由前端设置
                }
                for node_id, node in self.nodes.items()
            },
            "edges": [
                {
                    "source": source_id,
                    "target": target_id,
                    "id": f"{source_id}-{target_id}"
                }
                for source_id, target_id in self.edges
            ]
        }
    
    @classmethod
    def from_workflow_data(cls, workflow_data: Dict[str, Any], tool_registry=None) -> "OmicsGraph":
        """
        从工作流数据创建图（集成现有 Omics-Agent）
        
        Args:
            workflow_data: 工作流配置（来自 SOPPlanner）
            tool_registry: 工具注册表（用于查找工具函数）
            
        Returns:
            OmicsGraph 实例
        """
        workflow_name = workflow_data.get("workflow_name", "Omics Workflow")
        steps = workflow_data.get("steps", [])
        
        graph = cls(workflow_name=workflow_name)
        
        # 映射节点类型
        type_mapping = {
            "inspect_data": NodeType.DATA_LOADER,
            "preprocess_data": NodeType.PREPROCESSOR,
            "pca_analysis": NodeType.ANALYZER,
            "differential_analysis": NodeType.ANALYZER,
            "visualize_volcano": NodeType.VISUALIZER,
            "pathway_enrichment": NodeType.ANALYZER,
        }
        
        # 添加节点
        for step in steps:
            step_id = step.get("step_id") or step.get("id", "")
            tool_id = step.get("tool_id", step_id)
            step_name = step.get("name") or step.get("step_name", tool_id)
            params = step.get("params", {})
            
            # 查找工具函数
            tool_func = None
            if tool_registry:
                tool_func = tool_registry.get_tool(tool_id)
            
            # 确定节点类型
            node_type = type_mapping.get(tool_id, NodeType.GENERIC)
            
            # 添加节点
            node_id = graph.add_node(
                node_id=step_id,
                name=step_name,
                node_type=node_type.value,
                tool_id=tool_id,
                description=step.get("description", ""),
                params=params,
                tool_func=tool_func
            )
        
        # 添加边（基于步骤顺序）
        node_ids = list(graph.nodes.keys())
        for i in range(len(node_ids) - 1):
            graph.add_edge(node_ids[i], node_ids[i + 1])
        
        return graph


class WorkflowVisualizer:
    """
    工作流可视化器（主接口）
    
    将 Omics-Agent 的工作流转换为可视化格式
    """
    
    def __init__(self, tool_registry=None):
        self.tool_registry = tool_registry
    
    def visualize_workflow(self, workflow_data: Dict[str, Any]) -> Dict[str, Any]:
        """
        可视化工作流
        
        Args:
            workflow_data: 工作流配置
            
        Returns:
            可视化数据（用于前端渲染）
        """
        graph = OmicsGraph.from_workflow_data(workflow_data, self.tool_registry)
        return graph.to_dict()
    
    def create_graph_from_steps(self, steps: List[Dict[str, Any]], workflow_name: str = "Omics Workflow") -> OmicsGraph:
        """
        从步骤列表创建图
        
        Args:
            steps: 步骤列表
            workflow_name: 工作流名称
            
        Returns:
            OmicsGraph 实例
        """
        workflow_data = {
            "workflow_name": workflow_name,
            "steps": steps
        }
        return OmicsGraph.from_workflow_data(workflow_data, self.tool_registry)

