# 部署后重启说明（注册 404 / 500 与前端报错）

## Table 'omics_agent.xxx' doesn't exist（表不存在）

若 500 响应中 `error_type` 为 `ProgrammingError` 且 `error_message` 包含 **Table 'omics_agent.sessions' doesn't exist**（或 assets/users/messages/workflow_templates），说明库已连上但表未建。

**处理：** 在浏览器或 curl 调用一次建表接口（幂等，只创建缺失表）：

```text
GET http://<你的后端地址>:8028/api/db/init
```

例如：`http://192.168.32.31:8028/api/db/init`  

返回 `{"detail": "表已就绪", "ok": true}` 即成功，再刷新前端即可。

---

## MySQL 1030 / error 168（建表时存储引擎报错）

若调用 `/api/db/init` 或首次建表时返回 **OperationalError**，且 `error_message` 包含 **1030** 或 **error 168**、**Unknown (generic) error from engine**，说明是 **MySQL 存储引擎层** 报错，不是应用代码问题。

**常见原因：**

- **数据目录权限不足**：MySQL 进程对数据目录（如 `/var/lib/mysql`）无写权限，InnoDB 无法创建/写入表文件。
- **磁盘空间不足**：数据目录所在分区已满。
- **Docker 卷权限**：MySQL 容器挂载的数据卷权限或归属不对。

**排查步骤：**

1. 查看 MySQL 错误日志（宿主机或容器内），常有更具体的 OS 错误号，例如：
   - `[InnoDB] Operating system error number 13` → 权限拒绝（Permission denied）
   - `[InnoDB] Cannot create file` → 无法创建文件（权限或磁盘）

2. 检查数据目录权限（示例路径，按实际环境改）：
   ```bash
   ls -la /var/lib/mysql
   # 应为 mysql 用户可写，例如：
   # drwxr-xr-x  mysql mysql  ...
   ```

3. Docker 部署时，确认 MySQL 数据卷挂载正确且容器内 MySQL 用户对该目录可写：
   ```yaml
   volumes:
     - mysql_data:/var/lib/mysql
   ```
   必要时在宿主机执行：`chown -R 999:999 <卷路径>`（999 为常见 mysql 容器 UID，以实际为准）。

4. 检查磁盘空间：`df -h` 确保 MySQL 数据所在分区有足够空间。

修复权限或磁盘后，再次访问 `GET /api/db/init` 或刷新页面触发自动建表即可。

---

## 如何在浏览器中查看 500 的真实错误详情

后端已接入**全局异常显微镜**：所有未捕获异常都会以 JSON 返回，包含 `error_type`、`error_message`、`traceback`。

1. 打开开发者工具：`F12` 或 右键 → 检查。
2. 切到 **Network（网络）** 面板。
3. 触发会 500 的请求（如注册、打开侧栏拉会话/资产）。
4. 在请求列表里点选该请求（状态码为 500 的那条）。
5. 在右侧 **Response（响应）** 或 **Preview** 中查看响应体，例如：
   - `error_type`: 异常类型（如 `RuntimeError`、`OperationalError`）。
   - `error_message`: 异常信息。
   - `traceback`: 完整堆栈，可直接定位到出错文件和行号。

根据 `error_type` 与 `traceback` 即可精确定位根因（如数据库不可用、表不存在、bcrypt 未安装等）。

---

## 1. 404：新路由未加载

在 Docker 部署环境下，**FastAPI 进程在容器启动时加载路由**。若新增了 `auth.py` 等路由但**没有重启后端容器**，请求会返回 **404 Not Found**。

**解决：** 在项目根目录执行：

```bash
docker-compose restart api-server
```

## 2. 500 / “Unexpected token 'I'... is not valid JSON” 的修复说明

若后端返回 **500 Internal Server Error** 且响应体为纯文本（如 "Internal Server Error"），前端在 `response.json()` 时会抛 **“Unexpected token 'I'... is not valid JSON”**。当前版本已做如下处理：

- **后端**
  - 数据库不可用时，依赖会抛出 `RuntimeError`，由全局异常处理统一返回 **503 + JSON** `{"detail": "数据库暂不可用，请稍后重试"}`。
  - 表不存在/连接失败等 `OperationalError` 时，同样返回 **503 + JSON**，避免 500 纯文本。
  - 启动时若 MySQL 可用会自动执行建表（`Base.metadata.create_all`），建表失败仅打日志，不阻止进程启动。
- **前端**
  - 注册/登录请求使用 `parseJsonOrReject(r)`：先 `r.text()` 再解析，失败时统一 `reject({ detail })`，弹窗显示 `err.detail`。
  - 历史会话/数据资产接口在非 2xx 时也先读 body 再解析，侧栏可显示 `detail`（如“数据库暂不可用，请稍后重试”）。

**建议自测：**

1. 保证 MySQL 已启动且库存在，必要时执行建表（应用启动时会自动建表，若存储异常可检查 MySQL 日志）。
2. 重启后端使新异常处理与建表逻辑生效：`docker-compose restart api-server`。
3. 打开前端：注册/登录应正常；若 DB 不可用，应看到“数据库暂不可用，请稍后重试”等提示，而不再出现 “Unexpected token” 的 JSON 解析错误。

## 3. 首次部署或改动了 Dockerfile/依赖时

若修改了 `requirements.txt` 或 `Dockerfile`，需要重新构建并启动：

```bash
docker-compose up -d --build api-server
```
