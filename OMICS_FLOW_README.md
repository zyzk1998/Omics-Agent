# Omics-Flow Visualizer: RAGFlow-style Workflow Builder

## 🎯 概述

Omics-Flow 是一个可视化工作流构建器，将 Omics-Agent 的后端逻辑转换为类似 RAGFlow 的拖拽式界面。

### 核心映射逻辑

| Omics-Agent 概念 | UI 表现 | 功能 |
|-----------------|---------|------|
| Agent / Tool (e.g., DESeq2, KeggSearch) | 节点 (Node) | 显示名称、状态（运行中/完成/报错）、耗时 |
| Analysis Pipeline (e.g., RNA-Seq Flow) | 连线 (Edge) | 显示数据流向（Count Matrix -> Diff Expr -> Volcano Plot） |
| Parameters (e.g., p_value < 0.05) | 侧边栏 (Sidebar) | 点击节点后弹出，允许在线修改参数 |
| Execution Logs | 底部控制台 (Console) | 实时显示 Agent 的思考过程（Stream Output） |

## 🚀 快速开始

### 1. 安装依赖

```bash
pip install streamlit networkx matplotlib
```

### 2. 运行可视化应用

```bash
streamlit run app_visualizer.py
```

应用将在 `http://localhost:8501` 启动。

## 📖 使用指南

### 方式一：手动构建工作流

1. **添加节点**
   - 在侧边栏选择节点类型（Data Loader, Preprocessor, Differential Analysis 等）
   - 配置参数（文件路径、p-value 阈值等）
   - 点击 "➕ Add Node" 添加到画布

2. **自动连接**
   - 节点会自动按添加顺序连接
   - 形成数据流：前一个节点的输出 → 下一个节点的输入

3. **执行工作流**
   - 点击 "🚀 Execute Workflow" 按钮
   - 查看每个节点的执行状态和日志

### 方式二：加载标准工作流

1. 在侧边栏选择工作流类型（Metabolomics 或 RNA）
2. 点击 "📥 Load Standard Workflow"
3. 系统会自动加载该类型的标准分析流程

### 方式三：从现有工作流数据创建

```python
from gibh_agent.core.visualizer import WorkflowVisualizer, OmicsGraph
from gibh_agent.core.tool_registry import registry

# 创建可视化器
visualizer = WorkflowVisualizer(tool_registry=registry)

# 从工作流数据创建图
workflow_data = {
    "workflow_name": "Metabolomics Analysis",
    "steps": [
        {"step_id": "inspect_data", "tool_id": "inspect_data", "name": "数据检查", "params": {"file_path": "/app/uploads/data.csv"}},
        {"step_id": "preprocess_data", "tool_id": "preprocess_data", "name": "数据预处理", "params": {}}
    ]
}

graph = visualizer.visualize_workflow(workflow_data)
```

## 🏗️ 架构说明

### 核心组件

1. **`WorkflowNode`** - 工作流节点
   - 映射到 UI 节点
   - 包含工具函数、参数、状态、输出

2. **`OmicsGraph`** - 工作流图
   - 管理节点和边
   - 计算拓扑排序
   - 执行工作流

3. **`WorkflowVisualizer`** - 可视化器
   - 将 Omics-Agent 工作流转换为可视化格式
   - 集成 ToolRegistry 和 WorkflowRegistry

### 数据流

```
SOPPlanner.generate_plan()
    ↓
workflow_data (steps + params)
    ↓
OmicsGraph.from_workflow_data()
    ↓
可视化节点图 (JSON)
    ↓
Streamlit UI 渲染
```

## 🎨 UI 功能

### 侧边栏（Node Builder）

- **节点类型选择**：Data Loader, Preprocessor, Analyzer, Visualizer
- **参数配置**：根据节点类型动态显示参数输入框
- **添加节点**：一键添加到工作流

### 主画布（Workflow Canvas）

- **节点列表**：显示所有节点及其状态
- **参数查看**：点击节点展开查看/编辑参数
- **执行日志**：实时显示每个节点的执行日志
- **执行历史**：记录每次执行的详细信息

### 导出功能

- **JSON 导出**：导出工作流配置为 JSON
- **JSON 导入**：从 JSON 文件导入工作流

## 🔧 集成现有系统

### 与 Orchestrator 集成

```python
from gibh_agent.core.orchestrator import AgentOrchestrator
from gibh_agent.core.visualizer import WorkflowVisualizer

# 创建 Orchestrator
orchestrator = AgentOrchestrator(agent, upload_dir="/app/uploads")

# 生成工作流
result = await orchestrator.stream_process(query="分析数据", files=[...])

# 可视化工作流
visualizer = WorkflowVisualizer(tool_registry=registry)
visualization = visualizer.visualize_workflow(result.get("workflow_data"))
```

### 与 Planner 集成

```python
from gibh_agent.core.planner import SOPPlanner
from gibh_agent.core.visualizer import WorkflowVisualizer

# 生成计划
planner = SOPPlanner(tool_retriever, llm_client)
plan = await planner.generate_plan(query="代谢组分析", file_metadata=metadata)

# 可视化计划
visualizer = WorkflowVisualizer(tool_registry=registry)
graph = OmicsGraph.from_workflow_data(plan.get("workflow_data"), tool_registry=registry)
```

## 📝 示例工作流

### Metabolomics 标准流程

1. **数据检查** (inspect_data)
   - 输入：CSV 文件
   - 输出：数据统计信息

2. **数据预处理** (preprocess_data)
   - 输入：原始数据
   - 参数：log_transform, standardize
   - 输出：预处理后的数据

3. **PCA 分析** (pca_analysis)
   - 输入：预处理数据
   - 参数：n_components=2
   - 输出：PCA 结果

4. **差异分析** (differential_analysis)
   - 输入：预处理数据
   - 参数：p_value_threshold, fold_change_threshold
   - 输出：差异代谢物列表

5. **火山图可视化** (visualize_volcano)
   - 输入：差异分析结果
   - 输出：火山图

## 🎯 未来扩展

- [ ] 拖拽式节点编辑（使用 react-flow 或类似库）
- [ ] 实时执行监控（WebSocket）
- [ ] 节点参数验证（Pydantic schema）
- [ ] 工作流模板库
- [ ] 执行结果可视化（图表、表格）

## 📚 相关文档

- [Omics-Agent Architecture](./ARCHITECTURE_REFACTORING.md)
- [Tool Registry](./gibh_agent/core/tool_registry.py)
- [Workflow Registry](./gibh_agent/core/workflows/base.py)

