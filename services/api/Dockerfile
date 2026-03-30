# syntax=docker/dockerfile:1
# 启用 BuildKit 以支持 RUN --mount=type=cache（pip 缓存加速重建）
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

# 安装系统级依赖（编译生信库需要）
# bedtools 为 pybedtools (Epigenomics) 所需
RUN apt-get update && apt-get install -y \
    build-essential \
    gcc \
    g++ \
    libhdf5-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    gosu \
    bedtools \
    && rm -rf /var/lib/apt/lists/*

# 配置 pip 使用国内镜像源（清华大学镜像，速度更快）
RUN pip config set global.index-url https://pypi.tuna.tsinghua.edu.cn/simple && \
    pip config set global.trusted-host pypi.tuna.tsinghua.edu.cn && \
    pip config set global.timeout 60

# 先复制依赖文件（利用 Docker 层缓存）
COPY requirements.txt .

# 安装 Python 依赖（层缓存失效时用 BuildKit 缓存挂载加速 pip 下载；需 DOCKER_BUILDKIT=1）
RUN --mount=type=cache,target=/root/.cache/pip \
    pip install --upgrade pip && \
    pip install -r requirements.txt

# BepiPred-3.0 独立 venv（必须放在 /opt，避免 compose 挂载 .:/app 覆盖 /app 下 .venv）
# 默认 CPU 版 torch；GPU 构建示例：--build-arg BEPIPRED_TORCH_INDEX=https://download.pytorch.org/whl/cu121
ARG BEPIPRED_TORCH_INDEX=https://download.pytorch.org/whl/cpu
COPY third_party/BepiPred-3.0/requirements.txt /opt/bepipred3/requirements.txt
RUN --mount=type=cache,target=/root/.cache/pip \
    python3 -m venv /opt/bepipred3-venv && \
    /opt/bepipred3-venv/bin/pip install --upgrade pip && \
    /opt/bepipred3-venv/bin/pip install torch torchvision --index-url "${BEPIPRED_TORCH_INDEX}" && \
    /opt/bepipred3-venv/bin/pip install -r /opt/bepipred3/requirements.txt

ENV BEPIPRED3_PYTHON=/opt/bepipred3-venv/bin/python

# 最后复制项目代码（代码变化时只需要重新构建这一层）
# 使用 .dockerignore 排除大文件，只复制必要的代码
COPY gibh_agent/ /app/gibh_agent/
COPY server.py /app/
COPY scripts/docker_api_entrypoint.sh /entrypoint.sh

# 🔒 安全：创建非 root 用户（带 home 目录）
# 注意：UID 1000 应该与主机用户匹配，以便挂载目录有正确的权限
RUN groupadd -r appuser && \
    useradd -r -g appuser -u 1000 -m -d /home/appuser appuser

# 创建必要的目录并设置权限
RUN mkdir -p /app/uploads /app/results /app/data /home/appuser/.cache /home/appuser/.config && \
    chown -R appuser:appuser /app /home/appuser /opt/bepipred3-venv

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
CMD ["python", "server.py"]

