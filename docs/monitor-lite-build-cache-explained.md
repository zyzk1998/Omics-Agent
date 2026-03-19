# monitor-lite.sh「有缓存 / 无缓存」说明

## 结论：脚本行为与名称一致

| 选项 | 脚本实际执行 | 含义 |
|------|--------------|------|
| **4) 有缓存重建** | `docker compose up -d --build` | 使用 **Docker 层缓存**：未改动的层直接复用，不传 `--no-cache`。 |
| **5) 无缓存彻底重建** | `docker compose build --no-cache` 再 `up -d` | **Docker 忽略所有层缓存**，每层都重新执行。 |

因此「有缓存」= Docker 构建缓存，「无缓存」= Docker 不用缓存，与脚本菜单描述一致。

---

## 为何选 4 仍会长时间跑 pip？

你看到的 `pip install --no-cache-dir` 来自 **Dockerfile**（`services/api/Dockerfile` 第 42–43 行），不是 monitor-lite 控制的：

- **Docker 的「有缓存」**：只决定「这一层要不要重新执行」。
- **pip 的 `--no-cache-dir`**：表示 pip 安装时**不把下载的包写入本地缓存目录**（减小镜像体积、避免把缓存打进镜像）。

一旦 Docker 判定「安装依赖」这一层失效（例如 `requirements.txt` 或更早的 `COPY` 变了），就会重新执行 `RUN pip install ...`，此时 pip 仍会带 `--no-cache-dir`，从 PyPI/清华源重新下载所有包，所以会耗时长（例如 643s）。

终端里出现 `[api-server 7/12] RUN pip install ... 643.7s` 且**没有** `CACHED`，说明这一层**没有命中 Docker 缓存**，常见原因包括：

- 修改过 `requirements.txt`
- 构建上下文（如 `.dockerignore` 或被 COPY 的文件）有变化，导致前面某层失效，进而导致本层也重建

---

## 可选优化：让 pip 在多次构建间复用下载（BuildKit 缓存挂载）

若希望「即使 Docker 层失效、需要重新执行 pip 安装」时也更快，可以在 Dockerfile 里为 pip 使用 **BuildKit 缓存挂载**，让 pip 的下载缓存保存在主机上，下次构建复用。

在 `services/api/Dockerfile` 中，将：

```dockerfile
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir -r requirements.txt
```

改为（需 Docker BuildKit）：

```dockerfile
# syntax=docker/dockerfile:1
# 使用 BuildKit 缓存挂载，加速重复构建时的 pip 安装
RUN --mount=type=cache,target=/root/.cache/pip \
    pip install --upgrade pip && \
    pip install -r requirements.txt
```

构建时需开启 BuildKit，例如：

```bash
DOCKER_BUILDKIT=1 docker compose build
```

monitor-lite 若通过 `docker compose up -d --build` 调用，需在脚本中导出 `DOCKER_BUILDKIT=1`，或用户本机已默认开启 BuildKit，才会生效。

---

## 总结

- **monitor-lite 的 4/5**：实现正确，4 = 用 Docker 缓存，5 = 不用 Docker 缓存。
- **长时间 pip**：是「Docker 层失效 + Dockerfile 里 pip 使用 `--no-cache-dir`」导致每次重建该层都全量下载，与选项 4 的“有缓存”不矛盾。
- 若需在「层失效时」也加快 pip，可加上方 BuildKit 缓存挂载并启用 BuildKit。
