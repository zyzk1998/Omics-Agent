# 超算 MCP 侧情况简报

> 面向非 HPC 背景读者：用最少术语说明「c060」「root」差异与当前状态。

---

## 1. 系统里 MCP 是怎么接上的

- 用户/前端/API 不直接握 MCP 二进制协议，而是经 **`mcp-gateway`（本机 Docker，常见端口 8002）** 转发。
- 网关按环境变量 **`HPC_MCP_URL`**（默认 `http://192.168.32.31:8001/mcp`）连接 **上游 MCP 服务**（Streamable HTTP）。
- 对话里执行的 shell 类指令，对应上游工具 **`execute_hpc_command`**：**命令在哪个环境跑，完全由「谁在 `HPC_MCP_URL` 那个地址上提供服务」决定**，不是由前端写死。

---

## 2. 最初「直连」看到 **c060** 是怎么回事

- **`c060`** 是 **超算登录环境里的 Linux 用户名**（或同类账号标识），不是本机 Docker 默认用户。
- 当时在 **`HPC_MCP_URL` 所指的那台服务**上，`execute_hpc_command` 是在 **超算侧（或已与超算同账号、同命令环境）** 里执行的，所以 `whoami`、站点脚本 **`pestat`** 等会返回 **与真实集群一致** 的结果。
- 从链路看仍是：**应用 → 网关 → 上游 MCP**；「直连」指的是 **上游执行环境在超算一侧**，而不是浏览器跳过网关。

---

## 3. 后来变成 **root** 是怎么回事

- 部署用的 **`gpu-server` 在局域网中的地址就是 `192.168.32.31`**。
- 配置里 **`http://192.168.32.31:8001/mcp` 在路由上等于访问本机 `8001` 端口**。
- 曾加入的 **`hpc-mcp-stub`（开发用 Docker 桩）** 将 **宿主 `8001`** 映射进容器；该容器内进程用户为 **`root`**，也没有超算上的 Slurm/站点工具环境。
- 因此：**URL 未改、网关仍显示 `connected`**，但实际连到的是 **本机桩** → `whoami` 为 **`root`**，工具描述里出现「桩容器」等字样；**不是**超算 MCP「逻辑坏了」，而是 **端口上站错了程序**。

---

## 4. 当前 MCP 侧整体情况（截至桩已移除）

| 组件 | 状态 |
|------|------|
| **`mcp-gateway`** | **保留**。负责连接 `HPC_MCP_URL`、暴露 `/status`、`/call`、`/reconnect` 等。 |
| **`hpc-mcp-stub`** | **已从仓库与 compose 移除**；不再占用宿主 **`8001`**。 |
| **默认 `HPC_MCP_URL`** | 仍为 **`http://192.168.32.31:8001/mcp`**，语义为：**由你们在「该地址」上自行运行真实超算 MCP**。 |
| **恢复 c060 类效果的条件** | 在 **`192.168.32.31:8001`**（或你改为其它 IP/端口后的 `HPC_MCP_URL`）上 **启动并维持与以前相同的超算侧 MCP 服务**；然后 **重启/重连 `mcp-gateway`**（或 `POST /api/config/mcp` / `/reconnect`）。 |

---

## 5. 一句话结论

- **c060**：上游是 **超算环境**，命令在 **超算账号** 下执行。  
- **root**：上游实际是 **本机 Docker 桩** 占了 **`8001`**，命令在 **容器 root** 下执行。  
- **现在**：桩已撤，**`8001` 留给真实超算 MCP**；本仓库 **不包含** 超算登录节点上的 MCP 实现，**能否再看到 c060 取决于你们是否在约定地址上跑通该服务**。

---

## 6. 工具清单、Prompt 与 `工具库.md`（运维指引）

| 内容 | 说明 |
|------|------|
| **网关实时列表** | `GET http://127.0.0.1:8002/tools`（宿主机映射 `mcp-gateway:8002`），返回 `connected` + `tools[]`。与远端 `tools/list` 一致时，条数一般为 **22**（HPC + 工作站）。 |
| **仓库固定快照** | `docs/hpc_mcp_tools_catalog.json`：每条含 `mcp_name`、`openai_function_name`（`hpc_mcp_*`）、`description`。上游工具改名/增删后应 **更新该 JSON** 并同步改 `gibh_agent/core/openai_tools.py` 内 **`HPC_MCP_TOOL_ROUTING_CHEAT_SHEET`**（若路由描述需对齐）。 |
| **大模型选工具** | 启用 `compute_scheduler` 时，`hpc_chat_system_suffix()` 会在 system 末尾追加 **工具速查**，引导直接使用 `hpc_mcp_*`，减少滥用 `execute_hpc_command`。 |
| **根目录 `工具库.md`** | 由脚本生成：`PYTHONPATH=. python3 scripts/export_tool_library_md.py`。内含 Registry 全量 JSON + 附录 MCP 22 条 + DeepReAct 虚拟工具说明。 |
| **Chroma 工具向量库** | `ToolRetriever.sync_tools()` 会读取 `docs/hpc_mcp_tools_catalog.json` 并 **额外写入** MCP 条目（`category: HPC MCP`），便于语义检索命中超算工具名。变更目录后需 **重新执行同步**。 |
