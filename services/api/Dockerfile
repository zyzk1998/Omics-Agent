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

# 安装 Python 依赖（这一层会被缓存，除非 requirements.txt 改变）
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir -r requirements.txt

# 最后复制项目代码（代码变化时只需要重新构建这一层）
# 使用 .dockerignore 排除大文件，只复制必要的代码
COPY gibh_agent/ /app/gibh_agent/
COPY server.py /app/

# 🔒 安全：创建非 root 用户（带 home 目录）
# 注意：UID 1000 应该与主机用户匹配，以便挂载目录有正确的权限
RUN groupadd -r appuser && \
    useradd -r -g appuser -u 1000 -m -d /home/appuser appuser

# 创建必要的目录并设置权限
RUN mkdir -p /app/uploads /app/results /app/data /home/appuser/.cache /home/appuser/.config && \
    chown -R appuser:appuser /app /home/appuser

# 🔧 修复：创建一个启动脚本，在运行时修复挂载目录的权限
# 注意：这个脚本需要以 root 运行来修复权限，然后使用 gosu 切换到 appuser
RUN echo '#!/bin/sh\n\
set -e\n\
# 修复挂载目录的权限（如果挂载的目录权限不正确）\n\
if [ -d /app/uploads ]; then\n\
    chown -R appuser:appuser /app/uploads 2>/dev/null || true\n\
fi\n\
if [ -d /app/results ]; then\n\
    chown -R appuser:appuser /app/results 2>/dev/null || true\n\
fi\n\
if [ -d /app/data ]; then\n\
    chown -R appuser:appuser /app/data 2>/dev/null || true\n\
fi\n\
# 切换到 appuser 并执行原始命令\n\
exec gosu appuser "$@"' > /entrypoint.sh && \
    chmod +x /entrypoint.sh

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

