# 文档中心（`docs/`）

> 原仓库根目录下的指导与复用类 Markdown 已统一迁入本目录（**保留根目录 `README.md`**）。新增规范文档请优先放在此处并在本索引登记。

## 架构与开发规范

| 文档 | 说明 |
| --- | --- |
| [设计与模块.md](./设计与模块.md) | 九大模块分层与 Mermaid 简图；导出图见 [exports/](./exports/README.md) |
| [持久化层构建规范.md](./持久化层构建规范.md) | 会话、消息入库、`state_snapshot` |
| [思考过程和历史快照.md](./思考过程和历史快照.md) | 时光机四联索引 |
| [思考过程提取和输出.md](./思考过程提取和输出.md) | SSE `thought` / 快照 `reasoning` 实现细节 |
| [组学分析task构建规范文档.md](./组学分析task构建规范文档.md) | 新模态 Task / DAG / 工具注册 SOP |
| [STED_EC_细胞类型探针与报告展示规范.md](./STED_EC_细胞类型探针与报告展示规范.md) | STED-EC 细胞类型探针、自动标注、报告 Obs 表折叠 |
| [技能扩展规范文档.md](./技能扩展规范文档.md) | 技能种子、路由、广场扩容 |
| [PARAM_RENDER_LOGIC.md](./PARAM_RENDER_LOGIC.md) | 工作流参数表单渲染逻辑 |

## 索引与导出（脚本生成）

| 文档 | 生成命令 |
| --- | --- |
| [工具库.md](./工具库.md) | `PYTHONPATH=. python3 scripts/export_tool_library_md.py` |
| [技能库.md](./技能库.md) | `PYTHONPATH=. python3 scripts/export_skill_library_md.py` |
| [generated/bioinformatics_tools_matrix.md](./generated/bioinformatics_tools_matrix.md) | `PYTHONPATH=. python3 scripts/generate_bioinformatics_tools_matrix_md.py` |
| [OMICS_THREE_MODALITIES_BACKEND_RUN_REPORT.md](./OMICS_THREE_MODALITIES_BACKEND_RUN_REPORT.md) | `PYTHONPATH=. python3 scripts/generate_omics_three_modalities_backend_report.py` |

## 其它

| 文档 | 说明 |
| --- | --- |
| [试岗期技术沉淀路线.md](./试岗期技术沉淀路线.md) | 技术沉淀路线 |
| [转正汇报PPT素材.md](./转正汇报PPT素材.md) | 汇报素材 |
| [hpc_agent_architecture.html](./hpc_agent_architecture.html) | 系统级数据流交互演示 |
| [hpc_mcp_tools_catalog.json](./hpc_mcp_tools_catalog.json) | 超算 MCP 工具固定目录 |

根目录 [README.md](../README.md) 仍为主入口；API 手册见 [api/README.md](../api/README.md) 与 `services/nginx/html/API.md`（文档中心静态页）。
