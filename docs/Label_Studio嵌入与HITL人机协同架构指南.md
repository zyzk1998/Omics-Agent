# Label Studio 嵌入与 HITL 人机协同架构指南

> **文档性质**：GIBH-AGENT-V2 官方架构文档 · 与代码对齐（2026-06 · 同源反代 + LS 1.23 认证 + 科学语料硬 HITL）  
> **适用范围**：Label Studio 深度嵌入、Human-in-the-loop（HITL）闭环、三轨合一数据入库、科学语料数据加工技能  
> **关联规范**：`docs/ARCHITECTURE_CORE_PRINCIPLES.md`、`.cursor/rules/gibh-architecture-constitution.mdc`、`docs/思考过程和历史快照.md`

---

## 目录

1. [架构愿景](#1-架构愿景)
2. [核心状态机：软 HITL vs 硬 HITL](#2-核心状态机软-hitl-vs-硬-hitl)
3. [程序化触发（转录组 / 语料技能）](#3-程序化触发转录组--语料技能)
4. [前端 UI：独立 LS 卡 + iframe 模态 + 手动确认唤醒](#4-前端-ui独立-ls-卡--iframe-模态--手动确认唤醒)
5. [同源反代与浏览器会话桥接（2026-06 关键）](#5-同源反代与浏览器会话桥接2026-06-关键)
6. [Label Studio 1.23+ API 认证](#6-label-studio-123-api-认证)
7. [唤醒后端与持久化](#7-唤醒后端与持久化)
8. [数据打包与三轨入库](#8-数据打包与三轨入库)
9. [关键模块索引](#9-关键模块索引)
10. [环境变量与运维](#10-环境变量与运维)
11. [开发测试指南](#11-开发测试指南)
12. [常见失败排查](#12-常见失败排查)

---

## 1. 架构愿景

Label Studio 在本系统中是 **Human-in-the-loop 金标准出口**：专家在 LS 中修正 AI 对细胞类型、区域边界、语料图文标注等判断，系统再通过 `resume_from_hitl` 拉取标注并生成《专家分析报告（最终版）》或 SFT JSON。

**设计原则**：

- 工具内 **禁止** 阻塞轮询等待标注（TaaS 铁律）。
- **不依赖** LLM「凭心情」调用 `Trigger_Expert_Annotation`；标准流程在 **Executor / 专属 Agent 层程序化触发**。
- 采用 **软 HITL**（转录组默认）：初稿专家报告照常生成 → 下发 LS 入口 → 用户标注后 **手动确认** 唤醒 → 最终版报告。
- **硬 HITL**（科学语料等）：挂起等待标注，唤醒后导出结构化语料。

**浏览器嵌入铁律（2026-06）**：

- iframe **必须**与 Omics Agent **同源**（同 `host:port`），经 `/label-studio/` 反代访问 LS，禁止 iframe 直连 `:8082`（跨端口导致 CSRF 403 / 第三方 Cookie 失效）。
- 打开 iframe 前须经 **session 桥接** 写入 LS Django Cookie，并 **改写** LS 登录重定向路径。

---

## 2. 核心状态机：软 HITL vs 硬 HITL

### 2.1 软 HITL（转录组标准全流程 · 默认）

```
工作流执行 → rna_cell_annotation 成功
    → Executor 程序化 Trigger_Expert_Annotation
    → report_data.hitl + hitl_soft=true（不 break 后续步骤）
    → 继续执行并生成《专家分析报告（初稿）》
    → Orchestrator emit result
    → emit hitl_action（含 ls_url / project_id）
    → done{success} + state_snapshot.hitl_pending
    → 会话 status → waiting_for_hitl
    → 用户 LS 标注 → 点击「我已完成标注，继续生成报告」
    → POST resume_from_hitl → export_annotations → 《最终版》→ completed
```

| 字段 | 含义 |
| --- | --- |
| `report_data.hitl_soft` | `true` 表示软 HITL；Orchestrator **不会** 在初稿前 hard break |
| `report_data.hitl` | 含 `status: hitl_required`、`ls_url`、`project_id` |
| `state_snapshot.hitl_pending` | 前端渲染 LS 独立入口卡 |

### 2.2 硬 HITL（科学语料 / 步骤级 / Deep ReAct）

当某路径直接返回 `status: hitl_required` 且 **无** `hitl_soft` 时：

- Executor / 专属 Agent **break** 或挂起后续步骤；
- Orchestrator `done{waiting_for_hitl}`；
- 用户完成标注后走 `resume_from_hitl`。

**科学语料数据加工**（`skill_corpus_data_processing`）：

```
用户上传图像/文本
    → CorpusProcessingAgent.stream_skill_flow()
    → 程序化 Trigger_Expert_Annotation(scenario_type=generic_corpus_processing)
    → 硬挂起 WAITING_FOR_HITL + hitl_action
    → 用户 LS 标注 → resume_from_hitl
    → LLM 清洗为 {instruction, input, output}
    → 写入 /results/corpus_hitl/sft_corpus_*.json
```

| 组件 | 路径 |
| --- | --- |
| 专属 Agent | `gibh_agent/agents/corpus_processing_agent.py` |
| 技能 | `gibh_agent/skills/skill_corpus_data_processing.py` |
| 白名单场景 | `generic_corpus_processing`（`hitl_tools.py`） |
| 编排路由 | `orchestrator.py` 拦截 `[Skill_Route: skill_corpus_data_processing]` |

---

## 3. 程序化触发（转录组 / 语料技能）

### 3.1 转录组（Executor）

**实现位置**：`gibh_agent/core/executor.py` · `_provision_scrna_hitl_after_cell_annotation`

- 触发条件：`tool_id` / `step_id` 含 `cell_annotation` 且步骤 `success`
- 动作：`Trigger_Expert_Annotation(scenario_type=scrna_cell_type_annotation, image_path=...)`
- 图片 URL：`resolve_ls_accessible_image_url`
- 浏览器 LS 链接：`resolve_ls_browser_project_url` → 默认 `/label-studio/projects/{id}/data`
- 结果：`report_data.hitl` + `report_data.hitl_soft=true`

### 3.2 语料技能（CorpusProcessingAgent）

- 收到文件后调用 `trigger_corpus_hitl()`，不依赖 LLM 工具调用
- 仅图像时 `_normalize_tasks_for_scenario()` 自动补齐 `text` 占位字段（避免 LS 建项 HTTP 400）
- 唤醒后 `_generate_analysis_summary()` 导出 SFT JSON

### 3.3 编排器

`orchestrator.py`：

- 若 `results.hitl_soft` → **跳过** 硬挂起分支；
- `result` 事件之后优先使用 `results.hitl` 发送 `hitl_action`。

**白名单与 XML**：`gibh_agent/tools/hitl_tools.py` · `SCENARIO_TO_LS_XML` 覆盖 6 种 `scenario_type`（含 `generic_corpus_processing`）。

---

## 4. 前端 UI：独立 LS 卡 + iframe 模态 + 手动确认唤醒

**模块**：`services/nginx/html/js/components/hitl_ui.js`、`css/hitl_phase2.css`

### 4.1 布局

```
.hitl-artifact-row  (flex; space-between)
├── __left   → .analysis-artifact-card（查看工作台）
└── __right  → .hitl-entry-card（专家复核 / Label Studio）
```

### 4.2 独立 LS 卡（`.hitl-entry-card`）

| 元素 | 行为 |
| --- | --- |
| 说明文案 | 小白友好介绍 LS |
| **进入手动标注** | `bridgeLabelStudioSession()` → `openLabelStudioModal(ls_url)` |
| **跳过** | `POST resume_from_hitl {skip:true}`，保留初稿 |

> **已废弃（勿再文档化）**：嵌入「查看工作台」卡内的按钮、`.hitl-ls-entry-dock`、黄色 `.hitl-action-card` 三按钮卡。

### 4.3 模态框与唤醒（可靠路径）

1. 点击「进入手动标注」→ `POST /api/hitl/ls-session-bridge`（写入 LS Cookie）
2. iframe 加载 `ls_url`（默认 `/label-studio/projects/{id}/data`）
3. 用户在 iframe 内完成 LS 标注并 Submit
4. **必须** 点击模态框底部 **「我已完成标注，继续生成报告」**
5. 前端 `POST /api/sessions/{id}/resume_from_hitl`，SSE 消费后续事件
6. 右栏报告切换为 🟢 **专家分析报告 (最终版)** 或语料 SFT 产出；LS 卡移除

**URL 规范化**（`hitl_ui.js`）：

- `window.__LABEL_STUDIO_PUBLIC_URL` 默认 `/label-studio`
- `normalizeLsIframeUrl()`：将跨源绝对 URL 改写为同源 `/label-studio/projects/{id}/data`

**postMessage**：仅作辅助提示（高亮确认按钮 + toast），**不会** 自动调用 `resume_from_hitl`。

---

## 5. 同源反代与浏览器会话桥接（2026-06 关键）

### 5.1 架构示意

```
浏览器 (http://<host>:8028 或 :8018)
    │
    ├─ /api/*              → api-server:8028
    ├─ /label-studio/*     → 反代 → label-studio:8080（Docker 内网）
    │                         或 host.docker.internal:8082（宿主机 pip LS）
    └─ 静态页 /css /js     → nginx 或 api-server 静态挂载
```

**双通道反代**（任选其一对浏览器生效，建议同时配置）：

| 通道 | 配置位置 | 上游 |
| --- | --- | --- |
| api-server 直连 :8028 | `gibh_agent/api/label_studio_proxy.py` | `LABEL_STUDIO_URL` |
| nginx :8018 | `services/nginx/conf.d/default.conf` `location /label-studio/` | `label-studio:8080` |

nginx 额外配置：`proxy_redirect ~^/(.*)$ /label-studio/$1;`（配合 Location 头改写）。

### 5.2 `ls_url` 生成规则

`gibh_agent/utils/ls_public_url.py`：

- `LABEL_STUDIO_SAME_ORIGIN_PROXY=1`（默认）→ 相对路径 `/label-studio/projects/{id}/data`
- 显式 `LABEL_STUDIO_PUBLIC_URL` 非空时优先使用

### 5.3 登录重定向改写（修复 iframe `{"detail":"Not Found"}`）

**根因**：LS 未登录时 302 到 `/user/login/?next=/projects/...`（**无** `/label-studio` 前缀），浏览器跟随到 FastAPI 路由返回 JSON 404。

**修复**：`label_studio_proxy.py` · `rewrite_ls_proxy_location()`：

- `/user/login/` → `/label-studio/user/login/`
- `next=/projects/8/data/` → `next=/label-studio/projects/8/data/`

### 5.4 图像 Base64 内嵌（规避 CSP / 内网 DNS）

**根因**：`import_task` 若传 `http://nginx/uploads/...`，浏览器 iframe 无法解析 Docker 主机名，且 LS CSP `img-src` 拦截外部 HTTP。

**方案**（`hitl_tools.py`）：

- `resolve_local_image_file()`：将 `/uploads/`、`/results/`、`http://nginx/...` 解析为容器内物理路径
- `get_image_base64_data_uri()`：读文件 → `data:image/png;base64,...`
- `Trigger_Expert_Annotation` 导入前经 `_embed_task_image_fields_as_base64()` 写入 payload

```python
{"data": {"image": "data:image/png;base64,iVBORw0KGgo...", "text": "..."}}
```

浏览器渲染 `src="data:..."` 为同源内联数据，**无网络请求**，根除 `ERR_NAME_NOT_RESOLVED` 与 CSP 拦截。

### 5.5 Session 桥接（免二次登录）

| API | 说明 |
| --- | --- |
| `POST /api/hitl/ls-session-bridge` | api-server 用 `LABEL_STUDIO_USERNAME/PASSWORD` session 登录 LS，将 `sessionid` / `csrftoken` 写入浏览器 Cookie（`Path=/`） |

前端在 `openLabelStudioModal()` 内 **先** 调用 bridge，**再** 设置 `iframe.src`。

**验收**：

```bash
python3 -c "import httpx; r=httpx.post('http://127.0.0.1:8028/api/hitl/ls-session-bridge'); print(r.status_code, r.json())"
# 期望 200；随后带 Cookie 访问 /label-studio/projects/<id>/data 返回 text/html
```

---

## 6. Label Studio 1.23+ API 认证

LS **1.23.0** 起默认 `legacy_api_tokens_enabled=false`，`Authorization: Token <legacy>` 将返回：

```
HTTP 401: legacy token authentication has been disabled for this organization
```

### 6.1 api-server 侧（建项 / 导入 / 导出）

`gibh_agent/utils/ls_client.py`：

| 模式 | 条件 | 行为 |
| --- | --- | --- |
| **Session bootstrap**（默认） | `LABEL_STUDIO_API_KEY` 留空 | Django 登录 → API 请求带 Cookie |
| **PAT Bearer** | `LABEL_STUDIO_API_KEY` 为 JWT（`eyJ…`） | `POST /api/token/refresh` → `Authorization: Bearer <access>` |
| **Legacy Token** | 显式非 JWT 字符串 | `Authorization: Token …`；若 401 legacy disabled 则 **自动回退 session** |

### 6.2 浏览器侧（iframe 标注）

- 经 §5 session 桥接获得 LS Cookie，与 api-server bootstrap **独立**但使用同一管理员账号。
- `LABEL_STUDIO_HOST` 须与浏览器实际访问 origin 一致（含 `/label-studio` 子路径），用于 LS 生成静态资源绝对 URL。
- `CSRF_TRUSTED_ORIGINS` 须包含浏览器 origin（**无尾斜杠**），逗号分隔。

### 6.3 可选：启用 legacy Token（不推荐）

LS 容器环境变量 `LABEL_STUDIO_ENABLE_LEGACY_API_TOKEN=true`（或组织设置开启 Legacy Tokens）。新部署优先 session + PAT。

---

## 7. 唤醒后端与持久化

| API | 说明 |
| --- | --- |
| `POST /api/sessions/{session_id}/resume_from_hitl` | 前端手动确认；SSE 流式返回最终版 |
| `POST /api/hitl/ls-session-bridge` | iframe 打开前桥接 LS 会话（§5.4） |
| `POST /api/hitl/webhook` | LS Webhook（可选）；`LABEL_STUDIO_WEBHOOK_SECRET` |

**流水线**（`gibh_agent/core/hitl_resume.py`）：

1. 校验会话（`waiting_for_hitl` 或软 HITL：`completed/running + hitl_pending`）
2. `LabelStudioClient.export_annotations(project_id)`
3. `BaseAgent._generate_analysis_summary` + `hitl_annotations_json`
4. SSE：`diagnosis{hitl_final:true}` → `result` → `done`
5. 持久化 `state_snapshot` / `execution_snapshot.hitl_final` / `hitl_annotations`

语料路径：`corpus_processing_agent` 在唤醒后额外写入 `/results/corpus_hitl/sft_corpus_*.json`。

---

## 8. 数据打包与三轨入库

见 `gibh_agent/core/data_packager.py`、`ingestion_router.py`、`workspace_action_bar.js`。

- 入库 API：`POST /api/ingestion/trigger`（`skip_hitl` 可选）
- **跳过 LS 复核** 入口在 `.hitl-entry-card` 的「跳过」按钮。

---

## 9. 关键模块索引

| 层级 | 路径 | 职责 |
| --- | --- | --- |
| LS 客户端 | `gibh_agent/utils/ls_client.py` | session/PAT 认证；create / import / export |
| LS 浏览器 URL | `gibh_agent/utils/ls_public_url.py` | 同源 `/label-studio/` 基址解析 |
| LS 反代 | `gibh_agent/api/label_studio_proxy.py` | `/label-studio/*` → `LABEL_STUDIO_URL`；Location 改写 |
| HITL 工具 | `gibh_agent/tools/hitl_tools.py` | 白名单、XML、任务补齐、Trigger、URL 规范化 |
| 程序化触发 | `gibh_agent/core/executor.py` | 转录组 `_provision_scrna_hitl_after_cell_annotation` |
| 语料 Agent | `gibh_agent/agents/corpus_processing_agent.py` | 硬 HITL + SFT 导出 |
| 编排 | `gibh_agent/core/orchestrator.py` | 软/硬 HITL、`hitl_action` SSE、语料路由 |
| 唤醒 | `gibh_agent/core/hitl_resume.py` | 拉标注 + 最终版报告 |
| HITL API | `gibh_agent/api/routers/hitl_api.py` | bridge + webhook + resume |
| 会话映射 | `gibh_agent/core/hitl_session_registry.py` | session ↔ project_id |
| 前端 HITL | `services/nginx/html/js/components/hitl_ui.js` | bridge + LS 卡 + 模态 + resume |
| nginx 反代 | `services/nginx/conf.d/default.conf` | `/label-studio/` + `proxy_redirect` |
| Compose LS | `docker-compose.labelstudio.yml` | LS 容器、CSRF/HOST 环境变量 |
| 上传白名单 | `gibh_agent/core/omics_io_registry.py` | 图片格式允许上传（语料 PNG 等） |
| 测试 | 见 §11.3 | 回归单测列表 |

---

## 10. 环境变量与运维

| 变量 | 说明 |
| --- | --- |
| `LABEL_STUDIO_URL` | api-server → LS 内部 REST（compose：`http://label-studio:8080`） |
| `LABEL_STUDIO_PUBLIC_URL` | 浏览器 iframe 基址；默认同源 `/label-studio`（留空即可） |
| `LABEL_STUDIO_SAME_ORIGIN_PROXY` | `1`（默认）→ `ls_url` 为 `/label-studio/projects/{id}/data` |
| `LABEL_STUDIO_HOST` | LS 生成链接用，如 `http://127.0.0.1:8028/label-studio`（**须与浏览器 origin 一致**） |
| `LABEL_STUDIO_CSRF_TRUSTED_ORIGINS` | Django CSRF 白名单，逗号分隔、**无尾斜杠** |
| `LABEL_STUDIO_API_KEY` | **留空**（推荐）；或完整 PAT（JWT）；勿填已失效 legacy Token |
| `LABEL_STUDIO_USERNAME` / `LABEL_STUDIO_PASSWORD` | 零配置 bootstrap 管理员（与 LS 容器一致） |
| `OMICS_AGENT_WEB_URL` | LS 拉取 `/results/` 图片用（容器内可用 `http://nginx`） |
| `LABEL_STUDIO_WEBHOOK_SECRET` | Webhook 可选鉴权 |
| `GIBH_HITL_REGISTRY_DIR` | session ↔ project 映射目录 |

### 10.1 启动（Docker 叠加 LS）

```bash
cd /home/ubuntu/GIBH-AGENT-V2
# 内网部署：在 docker-compose.override.yml（已 gitignore）配置 HOST/CSRF/OMICS_AGENT_WEB_URL
# 注意：显式 -f 时不会自动加载 override，须写全三个文件：
docker compose -f docker-compose.yml -f docker-compose.labelstudio.yml -f docker-compose.override.yml up -d label-studio api-server nginx
```

**内网 IP 示例**（仅写在 `docker-compose.override.yml`，勿提交 Git）：

```yaml
services:
  api-server:
    environment:
      - OMICS_AGENT_WEB_URL=http://<浏览器IP>:8028
      - LABEL_STUDIO_URL=http://label-studio:8080
  label-studio:
    environment:
      - LABEL_STUDIO_HOST=http://<浏览器IP>:8028/label-studio
      - CSRF_TRUSTED_ORIGINS=http://<浏览器IP>:8028,http://<浏览器IP>:8018,http://127.0.0.1:8028,http://127.0.0.1:8018
```

### 10.2 宿主机 pip LS 与 Docker 端口冲突

若宿主机已占用 `:8082`（`label-studio start --port 8082`），勿再映射 Docker `8082:8080`。可选：

- 停止宿主机 LS，仅用 Docker `gibh_label_studio`；或
- nginx 上游改为 `host.docker.internal:8082`（须 `extra_hosts: host-gateway`）

`ls_data` 目录须对容器用户可写（`chmod -R a+rwX ls_data` 或 `chown 1001:0`）。

### 10.3 代码变更后重启

```bash
# 后端 Python 变更（挂载卷，重启即可）
docker restart gibh_v2_api

# nginx 配置变更
docker restart gibh_v2_nginx

# 仅前端静态：硬刷新浏览器 Ctrl+Shift+R
```

### 10.4 存量技能入库（语料技能）

```bash
cd /home/ubuntu/GIBH-AGENT-V2
PYTHONPATH=. python3 scripts/patch_corpus_skill.py
```

---

## 11. 开发测试指南

### 11.1 环境准备（本机 127.0.0.1）

```bash
LABEL_STUDIO_URL=http://label-studio:8080
LABEL_STUDIO_HOST=http://127.0.0.1:8018/label-studio
LABEL_STUDIO_CSRF_TRUSTED_ORIGINS=http://127.0.0.1:8018,http://127.0.0.1:8028
LABEL_STUDIO_API_KEY=          # 留空，走 session bootstrap
OMICS_AGENT_WEB_URL=http://127.0.0.1:8018
```

### 11.2 转录组 LS 闭环

1. 技能广场 → **转录组学标准全流程**
2. 上传 10x 三件套，执行工作流至 `rna_cell_annotation` 完成
3. 验收 §11.4 检查点 1–6

### 11.3 科学语料数据加工（硬 HITL）

1. 技能广场 → **科学语料数据加工**（`skill_corpus_data_processing`）
2. 上传 `test_corpus_umap.png` 或任意图像/文本
3. 等待 `hitl_action`（含 `project_id`、`ls_url`）
4. 「进入手动标注」→ iframe 显示 LS 图文模板
5. 标注 Submit → 「我已完成标注，继续生成报告」
6. 检查 `/results/corpus_hitl/sft_corpus_*.json`

### 11.4 验收检查点

| # | 检查项 | 预期 |
| --- | --- | --- |
| 1 | SSE `hitl_action` | 含 `ls_url`（`/label-studio/projects/...`）、`project_id` |
| 2 | 中栏双卡 | 左工作台 + 右 LS 入口卡 |
| 3 | `POST /api/hitl/ls-session-bridge` | 200，Set-Cookie `sessionid` |
| 4 | iframe | 显示 LS 标注 UI（非 JSON 404/403） |
| 5 | LS Submit + 模态确认按钮 | `resume_from_hitl` SSE 正常 |
| 6 | 产出 | 转录组：最终版报告；语料：SFT JSON |

### 11.5 单测

```bash
cd /home/ubuntu/GIBH-AGENT-V2
PYTHONPATH=. python3 -m pytest \
  tests/test_hitl_phase1.py \
  tests/test_hitl_phase2.py \
  tests/test_hitl_phase3.py \
  tests/test_executor_scrna_hitl.py \
  tests/test_ls_public_url.py \
  tests/test_label_studio_proxy.py \
  tests/test_corpus_processing_skill.py \
  tests/test_upload_image_extensions.py \
  -q
```

---

## 12. 常见失败排查

| 现象 | 原因 | 处理 |
| --- | --- | --- |
| 「进入手动标注」灰色 | LS 建项失败 / 认证失败 | 查 api-server 日志；确认 `label-studio` 容器 healthy |
| 建项 **401 legacy token disabled** | 使用旧 Token 头 | 清空 `LABEL_STUDIO_API_KEY`，重启 `gibh_v2_api` |
| iframe **CSRF 403** | 跨端口 iframe（8028 vs 8082） | 使用 `/label-studio/` 同源反代；配置 `CSRF_TRUSTED_ORIGINS` |
| iframe **`{"detail":"Not Found"}`** | LS 302 到 `/user/login/` 未加 `/label-studio` 前缀 | 确认反代 Location 改写 + session bridge；硬刷新 |
| iframe 空白 / 静态 404 | `LABEL_STUDIO_HOST` 与浏览器 origin 不一致 | override 中 HOST 改为实际 `http://<IP>:<port>/label-studio` |
| bridge **503 sessionid 未获取** | LS 未就绪或账号密码错误 | 检查 `LABEL_STUDIO_USERNAME/PASSWORD` 与 LS 容器日志 |
| 建项 **400 缺 text 键** | 仅图像未补齐语料字段 | 确认 `hitl_tools._normalize_tasks_for_scenario` 已部署 |
| PNG 上传被拒 | 上传白名单过窄 | `omics_io_registry.is_upload_allowed_filename` |
| resume 409 | 无 `hitl_pending` / 已 resume | 新开会话重跑 |
| LS 无图 / `http://nginx` ERR | 任务 image 仍为内网 HTTP URL | 确认 `resolve_ls_import_image_payload` Base64 逻辑已部署；重建 api-server |
| LS CSP 拦截图片 | 外部 HTTP 图源 | 同上，应变为 `data:image/...;base64,...` |
| LS 无图（拉取侧） | 图片 URL 不可达 | 设 `OMICS_AGENT_WEB_URL`；导入侧已优先 Base64，一般无需 `LS_IMAGE_FETCH_BASE_URL` |
| compose override 不生效 | 显式 `-f` 未包含 override 文件 | 命令行加 `-f docker-compose.override.yml` |

---

*本文档与 `tests/test_hitl_phase*.py`、`tests/test_executor_scrna_hitl.py`、`tests/test_corpus_processing_skill.py`、`tests/test_label_studio_proxy.py` 同步维护；若实现变更，须先改代码再回修本文档。*
