# GIBH-AGENT-V2 访问信息摘要

## 📍 当前服务状态

### 服务器信息
- **服务器 IP**: 192.168.32.31 (主 IP)
- **服务端口**: 8028
- **服务绑定**: 0.0.0.0:8028 (允许外部访问)

### 访问地址

#### 方式 1: 直接访问 API 服务器（当前配置）

**内网访问**:
- 前端界面: `http://192.168.32.31:8028`
- API 文档: `http://192.168.32.31:8028/api/docs`
- API 端点: `http://192.168.32.31:8028/api/chat`

**本地访问**:
- 前端界面: `http://localhost:8028`
- API 文档: `http://localhost:8028/api/docs`

#### 方式 2: 通过 Nginx 反向代理（需要配置）

如果启用了 Nginx 服务（端口 80）:
- 前端界面: `http://192.168.32.31`
- API 端点: `http://192.168.32.31/api/chat`

---

## 🚀 启动服务

### 检查服务状态

```bash
# 检查端口是否监听
netstat -tlnp | grep 8028
# 或
ss -tlnp | grep 8028

# 检查进程
ps aux | grep -E "(gunicorn|uvicorn|server.py)"

# 测试访问
curl http://localhost:8028/
```

### 启动服务

**方式 1: Docker Compose（推荐）**
```bash
cd /home/ubuntu/GIBH-AGENT-V2
docker compose up -d
docker compose logs -f api-server
```

**方式 2: 直接运行**
```bash
cd /home/ubuntu/GIBH-AGENT-V2
export SILICONFLOW_API_KEY="your_api_key"
python3 server.py
```

---

## 🔒 防火墙配置

### Ubuntu/Debian (ufw)
```bash
sudo ufw allow 8028/tcp
sudo ufw status
```

### CentOS/RHEL (firewalld)
```bash
sudo firewall-cmd --permanent --add-port=8028/tcp
sudo firewall-cmd --reload
```

---

## 🌍 外网访问（云服务器）

如果服务器在云平台，需要：

1. **配置安全组**:
   - 开放端口 `8028` (TCP)
   - 源地址: `0.0.0.0/0` (或限制特定 IP)

2. **获取公网 IP**:
   ```bash
   curl ifconfig.me
   ```

3. **访问地址**:
   - `http://<公网IP>:8028`

---

## ⚠️ 注意事项

1. **安全**: 当前 CORS 配置允许所有来源 (`*`)，生产环境应限制
2. **HTTPS**: 生产环境建议配置 HTTPS（需要 SSL 证书）
3. **认证**: 当前无认证机制，建议添加 API 密钥或用户认证
4. **防火墙**: 确保防火墙开放相应端口

---

## 📞 故障排查

1. **无法访问**: 检查防火墙、安全组、服务状态
2. **端口被占用**: `sudo lsof -i :8028`
3. **服务启动失败**: 查看日志 `docker compose logs api-server`

详细文档请参考: `ACCESS_GUIDE.md`
