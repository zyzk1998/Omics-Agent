# syntax=docker/dockerfile:1
# 启用 BuildKit 以支持 RUN --mount=type=cache（pip 缓存加速重建）
# HPC MCP 由独立服务 mcp-gateway（Python 3.11 + 官方 mcp SDK）承载；本镜像保持 3.9 稳定栈
FROM python:3.9-slim

WORKDIR /app

# 使用清华大学镜像源加速 apt-get（更稳定，支持更多版本）
RUN if [ -f /etc/apt/sources.list.d/debian.sources ]; then \
        sed -i 's|deb.debian.org|mirrors.tuna.tsinghua.edu.cn|g' /etc/apt/sources.list.d/debian.sources && \
        sed -i 's|security.debian.org|mirrors.tuna.tsinghua.edu.cn|g' /etc/apt/sources.list.d/debian.sources; \
    elif [ -f /etc/apt/sources.list ]; then \
        sed -i 's|deb.debian.org|mirrors.tuna.tsinghua.edu.cn|g' /etc/apt/sources.list && \
        sed -i 's|security.debian.org|mirrors.tuna.tsinghua.edu.cn|g' /etc/apt/sources.list; \
    else \
        # 检测 Debian 版本并设置对应的镜像源
        DEB_VERSION=$(cat /etc/os-release | grep VERSION_CODENAME | cut -d= -f2 || echo "bookworm") && \
        echo "deb https://mirrors.tuna.tsinghua.edu.cn/debian/ ${DEB_VERSION} main" > /etc/apt/sources.list && \
        echo "deb https://mirrors.tuna.tsinghua.edu.cn/debian/ ${DEB_VERSION}-updates main" >> /etc/apt/sources.list && \
        echo "deb https://mirrors.tuna.tsinghua.edu.cn/debian-security ${DEB_VERSION}-security main" >> /etc/apt/sources.list; \
    fi

# 安装系统级依赖（编译生信库 + 组学下游 CLI）
# bedtools 为 pybedtools (Epigenomics) 所需
# bwa/samtools/bcftools/bowtie2/fastp 为 Genomics 真实管线所需（见 scripts/install_omics_real_env.sh apt 段）
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    gcc \
    g++ \
    wget \
    ca-certificates \
    bzip2 \
    libhdf5-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    gosu \
    bedtools \
    bwa \
    samtools \
    bcftools \
    bowtie2 \
    fastp \
    default-jdk-headless \
    libxrender1 \
    libxext6 \
    libfontconfig1 \
    libxft2 \
    libcairo2 \
    openbabel \
    hmmer \
    picard-tools \
    trimmomatic \
    && rm -rf /var/lib/apt/lists/*

# 组学 Conda 栈：GATK4 / MACS2（独立 env，避免 base 与 Python 3.13 求解冲突）
RUN wget -qO /tmp/miniforge.sh https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh \
    && bash /tmp/miniforge.sh -b -p /opt/omics-conda \
    && /opt/omics-conda/bin/conda create -y -n omics-real -c conda-forge -c bioconda \
        python=3.11 gatk4 macs2 cnvkit delly snpeff \
    && rm -f /tmp/miniforge.sh
# 注意：组学 Conda 的 PATH 须在 BepiPred/chem venv 创建之后再 prepend，否则 `python3` 会变成 3.11 导致 pandas 编译失败

# 配置 pip 使用国内镜像源（清华大学镜像，速度更快）
RUN /usr/local/bin/python3 -m pip config set global.index-url https://pypi.tuna.tsinghua.edu.cn/simple && \
    /usr/local/bin/python3 -m pip config set global.trusted-host pypi.tuna.tsinghua.edu.cn && \
    /usr/local/bin/python3 -m pip config set global.timeout 60

# 先复制依赖文件（利用 Docker 层缓存）
COPY requirements.txt .

# 安装 Python 依赖（层缓存失效时用 BuildKit 缓存挂载加速 pip 下载；需 DOCKER_BUILDKIT=1）
RUN --mount=type=cache,target=/root/.cache/pip \
    /usr/local/bin/python3 -m pip install --upgrade pip && \
    /usr/local/bin/python3 -m pip install -r requirements.txt && \
    /usr/local/bin/python3 -m pip install --no-cache-dir 'biopython>=1.83' 'promb>=0.1.0' && \
    /usr/local/bin/python3 -c "from Bio import SeqIO; import promb; print('biopython+promb ok')" && \
    /usr/local/bin/python3 -m pip install --no-cache-dir 'anarci' || true

# RDKit 化学原子工具（chem_*）：独立 venv；须固定 numpy<2（与仓库 requirements 一致，避免 rdInfoTheory _ARRAY_API / NumPy 2 ABI 崩溃）
# 系统库 libxrender 等：MolToImage / Draw 依赖 libXrender.so.1（无头容器亦需安装）
RUN /usr/local/bin/python3 -m venv /opt/chem-rdkit-venv && \
    /opt/chem-rdkit-venv/bin/pip install -U pip && \
    /opt/chem-rdkit-venv/bin/pip install 'numpy>=1.24.0,<2' && \
    /opt/chem-rdkit-venv/bin/pip install --no-cache-dir 'rdkit-pypi>=2022.3.1' && \
    /opt/chem-rdkit-venv/bin/python -c \
    "import numpy as np; import rdkit; from rdkit.Chem import AllChem, DataStructs, Descriptors, Draw; print('chem-rdkit-venv ok', rdkit.__version__, 'numpy', np.__version__)"

ENV CHEM_RDKIT_PYTHON=/opt/chem-rdkit-venv/bin/python

# BepiPred-3.0 独立 venv（须用 /usr/local/bin/python3 创建，避免 Conda PATH 污染导致 pandas 编译失败）
# 默认 CPU 版 torch；GPU：--build-arg BEPIPRED_TORCH_INDEX=https://download.pytorch.org/whl/cu121
# INSTALL_BEPIPRED3=0 仅用于极快重建（生产/技能请保持 1；开发可挂载 third_party/BepiPred-3.0/.venv）
ARG INSTALL_BEPIPRED3=1
ARG BEPIPRED_TORCH_INDEX=https://download.pytorch.org/whl/cpu
COPY third_party/BepiPred-3.0/requirements.txt /opt/bepipred3/requirements.txt
RUN --mount=type=cache,target=/root/.cache/pip \
    if [ "${INSTALL_BEPIPRED3}" = "1" ]; then \
        /usr/local/bin/python3 -m venv /opt/bepipred3-venv && \
        /opt/bepipred3-venv/bin/pip install --upgrade pip 'setuptools>=65,<81' wheel && \
        /opt/bepipred3-venv/bin/pip install torch torchvision --index-url "${BEPIPRED_TORCH_INDEX}" && \
        /opt/bepipred3-venv/bin/pip install --only-binary=:all: 'numpy==1.22.4' 'pandas==1.4.3' && \
        /opt/bepipred3-venv/bin/pip install --no-build-isolation -r /opt/bepipred3/requirements.txt && \
        /opt/bepipred3-venv/bin/python -c "import torch; import esm; print('bepipred3-venv ok')"; \
    else \
        mkdir -p /opt/bepipred3-venv/bin && \
        printf '%s\n' 'INSTALL_BEPIPRED3=0: use mounted /app/third_party/BepiPred-3.0/.venv or entrypoint bootstrap.' \
            > /opt/bepipred3-venv/SKIPPED.txt; \
    fi

ENV BEPIPRED3_PYTHON=/opt/bepipred3-venv/bin/python

# 运行时：组学 CLI（cnvkit/delly/snpEff/gatk）优先；API 主进程仍为 /usr/local/bin/python3
ENV PATH="/opt/omics-conda/envs/omics-real/bin:/opt/omics-conda/bin:${PATH}"

# 最后复制项目代码（代码变化时只需要重新构建这一层）
# 使用 .dockerignore 排除大文件，只复制必要的代码
COPY gibh_agent/ /app/gibh_agent/
COPY server.py /app/
COPY scripts/docker_api_entrypoint.sh /entrypoint.sh
COPY scripts/init_omics_mock_references.py /app/scripts/init_omics_mock_references.py

# 🔒 安全：创建非 root 用户（带 home 目录）
# 注意：UID 1000 应该与主机用户匹配，以便挂载目录有正确的权限
RUN groupadd -r appuser && \
    useradd -r -g appuser -u 1000 -m -d /home/appuser appuser

# 创建必要的目录并设置权限
RUN mkdir -p /app/uploads /app/results /app/data /home/appuser/.cache /home/appuser/.config && \
    chown -R appuser:appuser /app /home/appuser /opt/bepipred3-venv /opt/chem-rdkit-venv /opt/omics-conda

# 入口脚本：权限修复 +（可选）运行时补齐 BepiPred venv
RUN chmod +x /entrypoint.sh

# 设置环境变量（解决 Numba 和 Matplotlib 权限问题）
ENV NUMBA_CACHE_DIR=/home/appuser/.cache/numba
ENV MPLCONFIGDIR=/home/appuser/.config/matplotlib
ENV HOME=/home/appuser
# 禁用 Numba 缓存（在容器环境中可能有问题，使用环境变量在运行时设置）
# 注意：NUMBA_DISABLE_CACHING 不是标准环境变量，需要在代码中设置

# 🔧 修复：entrypoint 脚本需要 root 权限来修复挂载目录权限
# 注意：entrypoint.sh 会在修复权限后使用 gosu 切换到 appuser
# 不要在这里切换用户，让 entrypoint 脚本处理
# USER appuser  # 注释掉，让 entrypoint 脚本处理用户切换

# 暴露端口
EXPOSE 8028

# 默认命令（会被 docker-compose 覆盖）
CMD ["/usr/local/bin/python3", "server.py"]

