# GIBH-AGENT-V2 访问指南

## 📋 当前服务状态

### 服务端口配置

- **API 服务器端口**: `8028` (直接暴露)
- **Nginx 端口**: `80` (如果启用，需要额外配置)

### 服务器 IP 地址

```bash
# 查看服务器 IP
hostname -I
# 或
ip addr show
```

常见 IP 地址：
- 内网 IP: `192.168.32.31`
- 其他内网 IP: `172.20.12.178`

---

## 🌐 访问方式

### 方式 1: 直接访问 API 服务器（当前配置）

**适用场景**: 开发环境、内网访问

**访问地址**:
- **前端界面**: `http://<服务器IP>:8028`
- **API 端点**: `http://<服务器IP>:8028/api/chat`
- **API 文档**: `http://<服务器IP>:8028/api/docs` (Swagger UI)
- **文件上传**: `http://<服务器IP>:8028/api/upload`

**示例**:
```bash
# 如果服务器 IP 是 192.168.32.31
前端: http://192.168.32.31:8028
API:  http://192.168.32.31:8028/api/chat
```

### 方式 2: 通过 Nginx 反向代理（推荐生产环境）

**适用场景**: 生产环境、需要 HTTPS、统一端口

**配置步骤**:

1. **启用 Nginx 服务**（在 `docker-compose.yml` 中添加）:
```yaml
nginx:
  image: nginx:alpine
  container_name: gibh_v2_nginx
  restart: always
  ports:
    - "80:80"
    - "443:443"  # 如果需要 HTTPS
  volumes:
    - ./services/nginx/conf.d:/etc/nginx/conf.d:ro
    - ./services/nginx/html:/usr/share/nginx/html:ro
  depends_on:
    - api-server
  networks:
    - gibh-network
```

2. **访问地址**:
- **前端界面**: `http://<服务器IP>` 或 `http://<域名>`
- **API 端点**: `http://<服务器IP>/api/chat`
- **API 文档**: `http://<服务器IP>/api/docs`

---

## 🔧 启动服务

### 方式 1: 使用 Docker Compose（推荐）

```bash
cd /home/ubuntu/GIBH-AGENT-V2

# 启动所有服务
docker compose up -d

# 查看服务状态
docker compose ps

# 查看日志
docker compose logs -f api-server
```

### 方式 2: 直接运行 Python 服务器

```bash
cd /home/ubuntu/GIBH-AGENT-V2

# 设置环境变量
export SILICONFLOW_API_KEY="your_api_key"
export SILICONFLOW_MODEL="deepseek-ai/DeepSeek-R1"

# 运行服务器
python3 server.py
# 或使用 gunicorn
gunicorn server:app -w 2 -k uvicorn.workers.UvicornWorker --bind 0.0.0.0:8028
```

---

## 🔒 安全配置

### 1. 防火墙设置

**Ubuntu/Debian (ufw)**:
```bash
# 允许 8028 端口
sudo ufw allow 8028/tcp

# 如果使用 Nginx，允许 80 和 443 端口
sudo ufw allow 80/tcp
sudo ufw allow 443/tcp

# 查看防火墙状态
sudo ufw status
```

**CentOS/RHEL (firewalld)**:
```bash
# 允许 8028 端口
sudo firewall-cmd --permanent --add-port=8028/tcp
sudo firewall-cmd --reload

# 如果使用 Nginx
sudo firewall-cmd --permanent --add-service=http
sudo firewall-cmd --permanent --add-service=https
sudo firewall-cmd --reload
```

### 2. CORS 配置

当前配置允许所有来源 (`ALLOWED_ORIGINS=*`)，生产环境应限制：

```bash
# 在 docker-compose.yml 或环境变量中设置
ALLOWED_ORIGINS=https://yourdomain.com,https://www.yourdomain.com
```

### 3. HTTPS 配置（生产环境推荐）

如果需要 HTTPS，需要：
1. 配置 SSL 证书
2. 修改 Nginx 配置支持 HTTPS
3. 使用 Let's Encrypt 或其他证书服务

---

## 📊 验证服务状态

### 检查服务是否运行

```bash
# 检查端口监听
netstat -tlnp | grep 8028
# 或
ss -tlnp | grep 8028

# 检查进程
ps aux | grep -E "(gunicorn|uvicorn|server.py)"

# 检查 Docker 容器
docker ps | grep gibh_v2
```

### 测试访问

```bash
# 测试前端
curl http://localhost:8028/

# 测试 API
curl http://localhost:8028/api/docs

# 从其他机器测试（替换为实际 IP）
curl http://<服务器IP>:8028/
```

---

## 🌍 外网访问配置

### 1. 云服务器配置

如果服务器在云平台（阿里云、腾讯云、AWS 等）：

1. **安全组配置**:
   - 开放端口 `8028` (TCP)
   - 如果使用 Nginx，开放端口 `80` 和 `443`

2. **获取公网 IP**:
   ```bash
   curl ifconfig.me
   # 或
   curl ip.sb
   ```

3. **访问地址**:
   - `http://<公网IP>:8028`

### 2. 内网穿透（如果需要）

如果服务器在内网，可以使用：
- **frp**: https://github.com/fatedier/frp
- **ngrok**: https://ngrok.com/
- **花生壳**: https://hsk.oray.com/

---

## 📝 常见问题

### Q1: 无法从其他机器访问

**检查项**:
1. ✅ 服务是否绑定到 `0.0.0.0`（不是 `127.0.0.1`）
2. ✅ 防火墙是否开放端口
3. ✅ 服务器网络是否可访问
4. ✅ 安全组（云服务器）是否配置正确

### Q2: 端口被占用

```bash
# 查看端口占用
sudo lsof -i :8028
# 或
sudo netstat -tlnp | grep 8028

# 停止占用端口的进程
sudo kill -9 <PID>
```

### Q3: 服务启动失败

```bash
# 查看日志
docker compose logs api-server
# 或
tail -f gibh_agent.log

# 检查环境变量
echo $SILICONFLOW_API_KEY
```

---

## 🎯 快速访问命令

```bash
# 获取服务器 IP
SERVER_IP=$(hostname -I | awk '{print $1}')
echo "访问地址: http://$SERVER_IP:8028"

# 测试服务
curl -I http://localhost:8028/
```

---

## 📞 技术支持

如有问题，请检查：
1. 服务器日志: `docker compose logs -f api-server`
2. 系统日志: `journalctl -u docker` (如果使用 systemd)
3. 项目文档: `README.md`, `API.md`

