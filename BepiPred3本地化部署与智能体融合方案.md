# BepiPred-3.0 本地化部署与 GIBH-AGENT-V2 智能体融合方案

**文档性质**：面向 RTX 6000（CUDA）生产服务器的运维部署与后端架构设计说明。  
**适用范围**：在国内网络环境下完成源码获取、Conda 环境、Hugging Face 模型离线缓存、本地推理服务与智能体工具链对接。  
**当前代码基线**：`gibh_agent/tools/protein_tools.py` 中 `bepipred3_prediction` **默认**通过 `asyncio.create_subprocess_exec` 调用本地 BepiPred-3.0 CLI，将产物写入 `RESULTS_DIR/bepipred3/<run_id>/` 并返回 `html_url` / `csv_url` / `zip_url`（与历史远程 API 契约一致）。设置 **`BEPIPRED3_USE_REMOTE_API=1`** 时可回退为磐石 HTTP。  
**说明**：`gibh_agent/core/utils.py` 已提供 **`@safe_tool_execution`**；下文保留的封装脚本与 env 表可与 **`BEPIPRED3_ROOT` / `BEPIPRED3_PYTHON` / `BEPIPRED3_SUBPROCESS_TIMEOUT`** 对照阅读（实现以代码为准）。

---

## 1. 环境与源码准备（国内源加速）

### 1.1 获取 BepiPred-3.0 源码（GitHub 镜像）

官方仓库 URL 请以论文 / 项目主页 / 作者公开仓库为准（不同发布渠道可能为 GitHub 或 Zenodo）。以下命令中的 `OWNER/REPO` 需替换为实际路径。

**方式 A：ghproxy 前缀（常用）**

```bash
# 示例：通过 ghproxy 加速克隆（若镜像不可用，请更换其他公益镜像或自建代理）
export GIT_REPO="https://github.com/OWNER/BepiPred-3"   # 替换为真实地址
git clone "https://mirror.ghproxy.com/${GIT_REPO}.git" /opt/bepipred3/BepiPred-3
cd /opt/bepipred3/BepiPred-3
```

**方式 B：Git 配置代理（企业内网 HTTP/SOCKS）**

```bash
git config --global http.proxy  http://127.0.0.1:7890
git config --global https.proxy http://127.0.0.1:7890
git clone https://github.com/OWNER/BepiPred-3.git /opt/bepipred3/BepiPred-3
```

**方式 C：离线包传输**

在可访问 GitHub 的环境完成 `git clone` 后，打包 `tar czf bepipred3-src.tgz BepiPred-3`，经内网堡垒机或 U 盘同步至 RTX 6000 服务器再解压。

**校验**：进入源码目录，确认 `README`、`requirements.txt` 或 `environment.yml`、`pyproject.toml` 等与官方文档一致；记录 **Python 版本下限** 与 **CUDA 对应 PyTorch 版本**。

---

### 1.2 Conda 独立环境（建议）

使用 **Miniconda / Mambaforge**，并将 **conda / pip** 指向国内镜像，避免构建环境时超时。

```bash
# 安装 Miniconda（示例：清华镜像安装脚本，以官网最新指引为准）
# wget https://mirrors.tuna.tsinghua.edu.cn/anaconda/miniconda/Miniconda3-latest-Linux-x86_64.sh
# bash Miniconda3-latest-Linux-x86_64.sh

conda create -n bepipred3-local python=3.10 -y
conda activate bepipred3-local
```

**Conda 频道镜像（清华示例）**

```bash
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/r
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/pytorch
conda config --set show_channel_urls yes
```

**PyTorch（CUDA 版，需与驱动匹配）**

RTX 6000 通常使用 **CUDA 11.8 或 12.x** 对应的官方 wheel。示例（以 PyTorch 官网生成命令为准，勿照抄版本号）：

```bash
# 示例：CUDA 12.1
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu121
```

**pip 国内镜像（清华）**

```bash
pip config set global.index-url https://pypi.tuna.tsinghua.edu.cn/simple
```

**Transformers、accelerate、Biopython 等**

```bash
pip install "transformers>=4.36" accelerate biopython pandas numpy scipy
```

若 BepiPred-3.0 仓库提供 `requirements.txt`，应 **优先** 使用：

```bash
pip install -r requirements.txt
```

---

## 2. 模型权重离线下载（Hugging Face 国内镜像）

### 2.1 推荐锁定的 ESM-2 权重

BepiPred-3.0 依赖 **ESM-2** 作为序列编码骨干。实践中最常见、与公开实现一致的检查点为：

| 用途 | Hugging Face `model_id` | 说明 |
|------|-------------------------|------|
| 标准高精度 | `facebook/esm2_t33_650M_UR50D` | 约 650M 参数，33 层，需 **约 16GB+ GPU 显存余量**（视 batch 与实现而定） |
| 资源紧张时 | `facebook/esm2_t30_150M_UR50D` 等 | 需确认 BepiPred-3.0 代码是否写死 `t33`；若写死则需改源码或保持 t33 |

**生产建议**：在方案评审中 **固定单一 model_id**，并在配置文件中显式记录 **revision（commit hash）**，避免线上静默升级导致结果漂移。

---

### 2.2 使用 HF 镜像全量下载到本地目录

**环境变量方式（推荐与现有脚本兼容）**

```bash
export HF_ENDPOINT=https://hf-mirror.com
export HF_HOME=/opt/models/huggingface
mkdir -p "$HF_HOME"
```

**使用 `huggingface-cli`**

```bash
pip install -U "huggingface_hub[cli]"
export HF_ENDPOINT=https://hf-mirror.com
huggingface-cli download facebook/esm2_t33_650M_UR50D \
  --local-dir /opt/models/huggingface/hub/models--facebook--esm2_t33_650M_UR50D \
  --local-dir-use-symlinks False
```

**纯离线推理时的环境变量（智能体进程侧同样建议设置）**

```bash
export TRANSFORMERS_OFFLINE=1
export HF_DATASETS_OFFLINE=1
```

---

### 2.3 修改 BepiPred-3.0 源码使模型指向本地路径

目标：**禁止** 在推理时访问外网，一律从 `/opt/models/...` 加载。

**思路一（推荐）**：在加载处使用 **本地目录** 作为 `model_name_or_path`：

```python
# 伪代码：在 BepiPred-3.0 内原 AutoModel.from_pretrained(...) 处
import os
local_esm = os.environ.get("BEPIPRED3_ESM_LOCAL", "/opt/models/huggingface/hub/models--facebook--esm2_t33_650M_UR50D")
model = AutoModel.from_pretrained(
    local_esm,
    local_files_only=True,
    trust_remote_code=False,
)
```

**思路二**：保留 `facebook/esm2_t33_650M_UR50D` 字符串，但依赖 Hugging Face 缓存目录已预填充，并设置：

```python
AutoModel.from_pretrained(
    "facebook/esm2_t33_650M_UR50D",
    local_files_only=True,
)
```

**运维检查清单**

- [ ] `config.json`、`pytorch_model.bin` 或 `model.safetensors`、`tokenizer` 相关文件完整。  
- [ ] 服务用户对工作目录与缓存目录有 **读权限**。  
- [ ] 首次联调时 **关闭** `TRANSFORMERS_OFFLINE`，确认缓存命中后再开启离线。

---

## 3. 智能体融合方案（代码级设计）

### 3.1 总体架构

| 组件 | 职责 |
|------|------|
| **本地推理子进程** | 在独立 Conda 环境中执行 BepiPred-3.0 官方/封装 CLI，输入 FASTA 路径，输出 JSON 或结构化文件 |
| **GIBH-AGENT-2 工具** | `bepipred3_prediction` 异步起子进程、限时长、读 stdout、吞并 stderr、映射为 `{"status","message",...}` |
| **装饰器** | `@registry.register` 注册工具；`@safe_tool_execution`（建议新增）统一捕获未处理异常，防止事件循环/主进程因未捕获错误崩溃 |

当前 `protein_tools.py` 已使用 `@registry.register`；`utils.py` 尚无 `@safe_tool_execution`，下文给出 **可粘贴实现的骨架**。

---

### 3.2 建议新增：`gibh_agent/core/utils.py` 中的 `@safe_tool_execution`

**设计目标**

- 将工具函数未捕获异常转为 `{"status": "error", "message": "..."}`，与 `SkillAgent`、前端展示约定一致。  
- 可选：记录结构化日志，不向外泄露堆栈细节（或仅 `logger.exception` 在服务端）。

```python
# 建议追加至 gibh_agent/core/utils.py（节选骨架）
from functools import wraps
import logging
import asyncio

logger = logging.getLogger(__name__)


def safe_tool_execution(func):
    """异步工具装饰器：统一吞掉异常，返回 Agent 可序列化的 error dict。"""
    if asyncio.iscoroutinefunction(func):

        @wraps(func)
        async def async_wrapper(*args, **kwargs):
            try:
                return await func(*args, **kwargs)
            except Exception as e:
                logger.exception("工具执行失败: %s", func.__name__)
                return {"status": "error", "message": f"{func.__name__} 执行异常: {e}"}

        return async_wrapper

    @wraps(func)
    def sync_wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except Exception as e:
            logger.exception("工具执行失败: %s", func.__name__)
            return {"status": "error", "message": f"{func.__name__} 执行异常: {e}"}

    return sync_wrapper
```

---

### 3.3 本地封装 CLI 脚本（与智能体解耦）

建议在 `/opt/bepipred3/BepiPred-3/scripts/run_bepipred3_local.py`（路径自定）实现：

- 参数：`--fasta`、`--top`、`--smoothing`、`--out-json`  
- 内部调用 BepiPred-3.0 推理管线，**仅向 stdout 打印一行 JSON** 或写入 `--out-json`，exit code 0 表示成功。

智能体 **不直接 import** 重型依赖，只负责 **传路径与读结果**，降低 Agent 进程与 CUDA 上下文耦合。

---

### 3.4 `bepipred3_prediction` 重构骨架（`asyncio.create_subprocess_exec`）

以下骨架展示：**异步子进程**、**明确 Conda 环境 Python 解释器**、**FASTA 落盘路径传递**、**异步读取 stdout/stderr**、**超时 kill**，并与 `@registry.register` + `@safe_tool_execution` 组合。

**环境变量约定（运维配置）**

| 变量 | 含义 |
|------|------|
| `BEPIPRED3_LOCAL_PYTHON` | 例如 `/opt/miniconda3/envs/bepipred3-local/bin/python` |
| `BEPIPRED3_LOCAL_RUNNER` | 封装脚本绝对路径，如 `/opt/bepipred3/BepiPred-3/scripts/run_bepipred3_local.py` |
| `BEPIPRED3_SUBPROCESS_TIMEOUT` | 秒，默认如 `600` |
| `BEPIPRED3_ESM_LOCAL` | 上一节本地 ESM 目录，供子进程内读取 |

```python
# 建议重构方向：gibh_agent/tools/protein_tools.py（骨架示例，非直接替换文件）
from __future__ import annotations

import asyncio
import json
import logging
import os
import tempfile
from typing import Any, Dict, Literal

from gibh_agent.core.tool_registry import registry
from gibh_agent.core.utils import safe_tool_execution  # 需先在 utils.py 实现

logger = logging.getLogger(__name__)

TopEpitopePercentageCutoff = Literal["top_20", "top_50", "top_70", "all"]


def _resolve_fasta_content(sequence_or_path: str) -> str:
    raw = (sequence_or_path or "").strip()
    if not raw:
        raise ValueError("sequence_or_path 为空")
    if os.path.isfile(raw):
        with open(raw, "r", encoding="utf-8", errors="replace") as f:
            return f.read()
    return raw


@registry.register(
    name="bepipred3_prediction",
    description=(
        "BepiPred-3.0：基于 ESM-2 的 B 细胞线性/构象表位预测（本地子进程执行）。"
        "输入可为 FASTA 文本或服务器上已存在的 FASTA 路径。"
    ),
    category="Proteomics",
    output_type="json",
)
@safe_tool_execution
async def bepipred3_prediction(
    sequence_or_path: str,
    top_epitope_percentage_cutoff: TopEpitopePercentageCutoff = "top_20",
    use_sequential_smoothing: bool = False,
) -> Dict[str, Any]:
    py_exe = os.getenv("BEPIPRED3_LOCAL_PYTHON", "").strip()
    runner = os.getenv("BEPIPRED3_LOCAL_RUNNER", "").strip()
    if not py_exe or not os.path.isfile(py_exe):
        return {"status": "error", "message": "未配置 BEPIPRED3_LOCAL_PYTHON 或文件不存在"}
    if not runner or not os.path.isfile(runner):
        return {"status": "error", "message": "未配置 BEPIPRED3_LOCAL_RUNNER 或文件不存在"}

    try:
        fasta_text = _resolve_fasta_content(sequence_or_path)
    except Exception as e:
        return {"status": "error", "message": f"输入解析失败: {e}"}

    if not fasta_text.strip().startswith(">"):
        fasta_text = f">query\n{fasta_text.strip()}"

    timeout = float(os.getenv("BEPIPRED3_SUBPROCESS_TIMEOUT", "600"))

    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False, encoding="utf-8") as tmp:
        tmp.write(fasta_text)
        fasta_path = tmp.name

    out_json_path = fasta_path + ".out.json"
    env = os.environ.copy()
    env.setdefault("BEPIPRED3_ESM_LOCAL", "/opt/models/huggingface/hub/models--facebook--esm2_t33_650M_UR50D")
    env.setdefault("TRANSFORMERS_OFFLINE", "1")
    env.setdefault("HF_HOME", "/opt/models/huggingface")

    proc: asyncio.subprocess.Process | None = None
    try:
        proc = await asyncio.create_subprocess_exec(
            py_exe,
            runner,
            "--fasta",
            fasta_path,
            "--top",
            top_epitope_percentage_cutoff,
            "--smoothing",
            "1" if use_sequential_smoothing else "0",
            "--out-json",
            out_json_path,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE,
            env=env,
        )
        try:
            stdout_b, stderr_b = await asyncio.wait_for(proc.communicate(), timeout=timeout)
        except asyncio.TimeoutError:
            proc.kill()
            await proc.communicate()
            return {"status": "error", "message": f"BepiPred3 子进程超时（>{timeout}s）"}

        stderr_text = (stderr_b or b"").decode("utf-8", errors="replace").strip()
        if stderr_text:
            # stderr 仅记录与回传摘要，避免主进程因未读管道导致阻塞；不当作成功依据
            logger.warning("BepiPred3 stderr: %s", stderr_text[:4000])

        if proc.returncode != 0:
            return {
                "status": "error",
                "message": f"BepiPred3 进程退出码 {proc.returncode}；stderr: {stderr_text[:2000]}",
            }

        if not os.path.isfile(out_json_path):
            out_preview = (stdout_b or b"").decode("utf-8", errors="replace")[:2000]
            return {
                "status": "error",
                "message": f"未生成输出文件；stdout 摘要: {out_preview}",
            }

        with open(out_json_path, "r", encoding="utf-8", errors="replace") as f:
            payload = json.load(f)

        # 按现有智能体约定归一化：例如 status / html_url / csv_url 等
        return {"status": "success", "message": "预测完成", **payload}

    finally:
        for p in (fasta_path, out_json_path):
            try:
                if p and os.path.isfile(p):
                    os.unlink(p)
            except OSError:
                logger.debug("临时文件删除失败: %s", p)
```

**要点说明**

1. **`asyncio.create_subprocess_exec`**：避免阻塞事件循环；与 FastAPI/异步 Orchestrator 同栈时必选。  
2. **独立 Conda Python**：通过 `BEPIPRED3_LOCAL_PYTHON` 固定解释器，避免与 Agent 主环境混用。  
3. **FASTA 文件路径**：智能体将用户粘贴的序列写入 **受控临时文件**，再把路径传给子进程，避免命令行长度限制与注入。  
4. **`stderr` 管道**：必须 `PIPE` 并 `communicate()` 消费，否则子进程可能 **写满管道阻塞**；错误时 **不抛未捕获异常**，而是返回 `status: error`。  
5. **超时**：`asyncio.wait_for` + `proc.kill()`，防止 GPU 死锁拖垮 Worker。  
6. **`@safe_tool_execution`**：最后一道防线，防止 `asyncio` / `OSError` 等遗漏导致 **主进程崩溃**。

---

### 3.5 与现有 HTTP 版本的迁移策略

| 阶段 | 行为 |
|------|------|
| 灰度 | `BEPIPRED3_LOCAL_PYTHON` 非空时走子进程；否则回退 `httpx` 调 `BEPIPRED3_PREDICT_URL` |
| 全量本地 | 删除或禁用远程 URL，强制本地子进程，未配置则 **启动时 fail-fast** |

---

## 4. 运维与安全建议

- **GPU 独占**：同一容器内避免多个 Worker 并发抢占 RTX 6000；可用 `CUDA_VISIBLE_DEVICES` + 队列或单 worker。  
- **资源限额**：systemd / cgroups 限制子进程 CPU、内存，防止 OOM 拖垮 Agent。  
- **输入大小**：限制 FASTA 总长度与文件大小，防止恶意超大输入。  
- **日志**：`stderr` 全文进日志文件，对用户仅返回摘要。  
- **复现性**：固定 `model_id` + `revision` + `conda env export` 导出环境描述。

---

## 5. 验收标准

- [ ] 无外网条件下，子进程可完成单次预测并返回约定 JSON。  
- [ ] 人为制造 GPU OOM / 脚本异常时，Agent 返回 `status: error`，**主服务不崩溃**。  
- [ ] 超时场景子进程被终止，无僵尸进程。  
- [ ] 与 `Skill_Route: bepipred3_prediction` 技能联调通过。

---

*文档版本：与仓库 `gibh_agent/tools/protein_tools.py`、`gibh_agent/core/utils.py` 当前实现对齐说明；落地时请以实际 BepiPred-3.0 仓库结构与官方依赖为准进行微调。*
