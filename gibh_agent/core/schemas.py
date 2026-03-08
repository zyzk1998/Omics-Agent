"""
工作流规划/诊断输出的标准 Schema（Pydantic）。

用于约束大模型输出的规划结构，包含 recommended_params 以支持动态参数智能推荐与注入。
"""
from typing import Dict, Any, List, Optional
from pydantic import BaseModel, Field


class WorkflowPlanResponse(BaseModel):
    """大模型规划/诊断输出的标准结构，供解析与校验使用。"""

    workflow_name: str = Field(default="Generated Workflow", description="工作流名称")
    steps: List[Dict[str, Any]] = Field(default_factory=list, description="步骤列表")
    recommended_params: Dict[str, Any] = Field(
        default_factory=dict,
        description="大模型推荐的底层工具执行参数键值对。可为全局或按 step_id/tool_id 分组，例如 "
        "{'n_pcs': 30, 'resolution': 0.8, 'min_genes': 200} 或 "
        "{'rna_clustering': {'resolution': 0.8}, 'rna_pca': {'n_comps': 30}}",
    )

    class Config:
        extra = "allow"
