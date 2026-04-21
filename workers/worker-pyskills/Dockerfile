# =============================================================================
# worker-pyskills — ubuntu:22.04 + Micromamba（BuildKit 缓存 + 清华源）
# -----------------------------------------------------------------------------
# 设计说明（给总指挥 / CI）：
# - 未采用多阶段构建：单阶段即可满足「小基底 + 依赖缓存 + 业务层分离」；
#   多阶段在此场景收益有限，反而增加层间传递复杂度。
# - 构建前请在宿主机启用 BuildKit，否则 RUN --mount=type=cache 会报错：
#     export DOCKER_BUILDKIT=1
#     docker compose build worker-pyskills
# =============================================================================
FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive \
    # Micromamba / conda 环境根前缀（conda 包缓存见 Step 4：独立挂载 + CONDA_PKGS_DIRS）
    MAMBA_ROOT_PREFIX=/opt/conda \
    QT_QPA_PLATFORM=offscreen \
    MPLBACKEND=Agg \
    PATH=/opt/conda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

WORKDIR /app

# -----------------------------------------------------------------------------
# Step 1 — 极简 APT + 国内镜像（清华；可改为阿里云 mirrors.aliyun.com/ubuntu）
# -----------------------------------------------------------------------------
# 仅安装 curl（拉 micromamba 安装脚本）、bzip2（部分包解压依赖）、ca-certificates。
RUN set -eux; \
    if [ -f /etc/apt/sources.list ]; then \
      sed -i 's@http://archive.ubuntu.com/ubuntu@http://mirrors.tuna.tsinghua.edu.cn/ubuntu@g' /etc/apt/sources.list; \
      sed -i 's@http://security.ubuntu.com/ubuntu@http://mirrors.tuna.tsinghua.edu.cn/ubuntu@g' /etc/apt/sources.list; \
    fi; \
    apt-get update; \
    apt-get install -y --no-install-recommends curl bzip2 ca-certificates; \
    rm -rf /var/lib/apt/lists/*

# -----------------------------------------------------------------------------
# Step 2 — micromamba 静态二进制（GitHub Release，避免 install.sh 行为/网络导致未写入 BIN_FOLDER）
# -----------------------------------------------------------------------------
# 标签与 https://github.com/mamba-org/micromamba-releases/releases 一致；多架构可扩展 micromamba-linux-aarch64 等。
ARG MICROMAMBA_VERSION=2.4.0-0
ENV HOME=/root

RUN set -eux; \
    mkdir -p /usr/local/bin; \
    curl -fSL "https://github.com/mamba-org/micromamba-releases/releases/download/${MICROMAMBA_VERSION}/micromamba-linux-64" \
      -o /usr/local/bin/micromamba; \
    chmod +x /usr/local/bin/micromamba; \
    /usr/local/bin/micromamba --version

# -----------------------------------------------------------------------------
# Step 3 — 国内 Conda 频道（在安装任何 conda 包之前生效）
# -----------------------------------------------------------------------------
COPY .condarc /root/.condarc

# -----------------------------------------------------------------------------
# Step 4 — 先复制依赖清单，再安装（业务代码延后 COPY，最大化缓存命中）
# -----------------------------------------------------------------------------
# 说明：总指挥模板中的 `micromamba install -n base -f` 在「空 root」场景下易与
# MAMBA_ROOT_PREFIX=/opt/conda 冲突（预建空目录会触发 libmamba「禁止覆盖 root prefix」）。
# 此处采用官方推荐的显式前缀引导：`create -p /opt/conda`，等价得到位于该前缀下的
# 单一 base 式环境，且与 `install … -n base` 一样由 environment.yml 驱动依赖集合。
# BuildKit 缓存注意：切勿把 cache 直接挂在 /opt/conda/pkgs —— 挂载会先创建父目录 /opt/conda，
# libmamba 会误判为「已有 root prefix」而报错 Overwriting root prefix is not permitted。
# 正确做法：缓存目录与 root prefix 分离，通过 CONDA_PKGS_DIRS 让 micromamba 写入独立挂载点。
COPY environment.yml /tmp/environment.yml

# third_party/drug-similarity：与 find_similar.py 同源使用 requests + RDKit +（可选）scipy/numpy；
# 依赖集合以 environment.yml 为准，保证容器内联网检索与本地脚本「同一土壤」。
RUN --mount=type=cache,target=/var/cache/mamba/pkgs \
    set -eux; \
    export CONDA_PKGS_DIRS=/var/cache/mamba/pkgs; \
    micromamba create -y -p /opt/conda -f /tmp/environment.yml && \
    micromamba clean --all --yes

# pyradiomics：清华 PyPI 上常无适用 manylinux wheel，需本地编译 C 扩展；不在全局基底永久安装 gcc，
# 仅在本层临时 apt 安装 build-essential，pip 成功后 purge，兼顾「基底包尽量精简」与可构建性。
RUN --mount=type=cache,target=/root/.cache/pip \
    set -eux; \
    if [ -f /etc/apt/sources.list ]; then \
      sed -i 's@http://archive.ubuntu.com/ubuntu@http://mirrors.tuna.tsinghua.edu.cn/ubuntu@g' /etc/apt/sources.list; \
      sed -i 's@http://security.ubuntu.com/ubuntu@http://mirrors.tuna.tsinghua.edu.cn/ubuntu@g' /etc/apt/sources.list; \
    fi; \
    apt-get update; \
    apt-get install -y --no-install-recommends build-essential; \
    /opt/conda/bin/python -m pip install --no-cache-dir --no-build-isolation \
      -i "https://pypi.tuna.tsinghua.edu.cn/simple" \
      "pyradiomics==3.0.1"; \
    apt-get purge -y build-essential; \
    apt-get autoremove -y --purge; \
    rm -rf /var/lib/apt/lists/*

# -----------------------------------------------------------------------------
# Step 5 — 业务代码（变更频繁，置于依赖层之后）
# -----------------------------------------------------------------------------
COPY . /app/

RUN /opt/conda/bin/python -c "import RNA; import pymol; from rdkit import Chem; import gseapy; from radiomics import featureextractor; import sklearn; from sa_score_compat import load_calculate_score; load_calculate_score(); print('pyskills imports OK')"

EXPOSE 8000

CMD ["/opt/conda/bin/uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8000"]
