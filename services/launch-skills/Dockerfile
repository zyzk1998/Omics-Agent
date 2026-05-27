# syntax=docker/dockerfile:1
# 首发 10 项 BaseSkill 专用镜像：不继承 api-server 组学/Conda/BepiPred 层，单独重建约 2–5 分钟。
FROM python:3.9-slim

WORKDIR /app

RUN if [ -f /etc/apt/sources.list.d/debian.sources ]; then \
        sed -i 's|deb.debian.org|mirrors.tuna.tsinghua.edu.cn|g' /etc/apt/sources.list.d/debian.sources && \
        sed -i 's|security.debian.org|mirrors.tuna.tsinghua.edu.cn|g' /etc/apt/sources.list.d/debian.sources; \
    elif [ -f /etc/apt/sources.list ]; then \
        sed -i 's|deb.debian.org|mirrors.tuna.tsinghua.edu.cn|g' /etc/apt/sources.list && \
        sed -i 's|security.debian.org|mirrors.tuna.tsinghua.edu.cn|g' /etc/apt/sources.list; \
    fi

RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    ncbi-blast+ \
    libxrender1 \
    libxext6 \
    libfontconfig1 \
    libxft2 \
    libcairo2 \
    && rm -rf /var/lib/apt/lists/*

RUN blastn -version && blastp -version

RUN python -m pip config set global.index-url https://pypi.tuna.tsinghua.edu.cn/simple && \
    python -m pip config set global.trusted-host pypi.tuna.tsinghua.edu.cn

COPY services/launch-skills/requirements.txt /tmp/launch-skills-requirements.txt

RUN --mount=type=cache,target=/root/.cache/pip \
    python -m pip install --upgrade pip && \
    python -m pip install --no-cache-dir -r /tmp/launch-skills-requirements.txt && \
    python -c "import httpx, chembl_webresource_client; from rdkit import Chem; print('launch-skills venv ok')"

COPY gibh_agent/ /app/gibh_agent/
COPY services/launch-skills/main.py /app/main.py
COPY scripts/verify_launch_skill_deps.sh /app/scripts/verify_launch_skill_deps.sh

ENV PYTHONPATH=/app \
    LAUNCH_SKILLS_IN_WORKER=1 \
    HOME=/tmp \
    MPLCONFIGDIR=/tmp/matplotlib

EXPOSE 8010

HEALTHCHECK --interval=30s --timeout=5s --start-period=15s --retries=3 \
    CMD python -c "import urllib.request; urllib.request.urlopen('http://127.0.0.1:8010/health', timeout=3)"

CMD ["python", "-m", "uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8010"]
