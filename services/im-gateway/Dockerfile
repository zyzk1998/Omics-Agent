# IM Gateway：独立接收各 IM 平台 Webhook，与主 API 解耦（Python 3.10）。
FROM python:3.10-slim

WORKDIR /app

RUN pip config set global.index-url https://pypi.tuna.tsinghua.edu.cn/simple && \
    pip config set global.trusted-host pypi.tuna.tsinghua.edu.cn && \
    pip config set global.timeout 120

COPY services/im-gateway/requirements.txt /app/requirements.txt
RUN pip install --no-cache-dir -r /app/requirements.txt

COPY services/im-gateway/main.py /app/main.py

EXPOSE 8030
CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8030"]
