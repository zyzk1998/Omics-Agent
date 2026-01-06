# GIBH-AGENT-V2 可视化监控运维脚本使用说明

## 📋 概述

`monitor.sh` 是一个综合性的可视化监控运维脚本，整合了所有 Docker、服务、日志监控等功能，方便调试和 review。

## 🚀 快速开始

```bash
# 直接运行（交互式菜单）
./monitor.sh

# 或使用 bash
bash monitor.sh
```

## 🎯 功能模块

### 1. Docker 服务管理
- **查看服务状态**：显示所有容器的运行状态和资源使用情况
- **启动服务**：创建必要目录并启动所有 Docker 服务
- **停止服务**：停止所有 Docker 服务
- **重启服务**：重启所有 Docker 服务
- **重新构建镜像**：强制重新构建 Docker 镜像（可选删除旧镜像）

### 2. 健康检查
- 检查 Docker 容器状态
- 检查 API 服务器响应
- 检查端口监听状态
- 检查直接运行的服务器进程

### 3. 日志监控
- **实时日志**：实时监控 API 服务器或 Worker 的日志输出
- **最近日志**：查看指定服务的最近 N 行日志
- **所有服务日志摘要**：查看所有服务的日志摘要
- **错误日志分析**：自动筛选和显示错误日志
- **本地日志文件**：查看本地日志文件（gibh_agent.log, server.log, debug.log）

### 4. 数据监控
- **数据状态**：
  - 上传文件统计（数量、大小、最近文件）
  - 结果文件统计（数量、大小、最近文件）
  - 磁盘使用情况
  - 目录大小排序
- **数据清理**：
  - 清理上传文件（保留元数据）
  - 清理结果文件
  - 清理所有数据

### 5. 错误诊断
- **诊断 502 错误**：
  - 容器状态检查
  - 日志分析
  - 依赖检查
  - 网络连接测试
- **自动修复 502 错误**：
  - 自动停止服务
  - 检查并修复 requirements.txt
  - 重新构建镜像
  - 启动服务并验证

### 6. 智能分析
- **LLM 解读 JSON 数据**：
  - 支持从文件读取或 stdin 输入
  - 自动调用 DeepSeek API 解读 JSON
  - 提供结构化说明和健康度评估

### 7. 综合监控面板
- 实时刷新监控面板（每 5 秒）
- 显示服务状态、API 健康、数据统计、最近错误
- 按 Ctrl+C 退出

## 📖 使用示例

### 示例 1：启动服务并查看状态
```bash
./monitor.sh
# 选择 2（启动服务）
# 选择 1（查看服务状态）
```

### 示例 2：监控实时日志
```bash
./monitor.sh
# 选择 7（实时日志 - API 服务器）
# 按 Ctrl+C 退出
```

### 示例 3：诊断 502 错误
```bash
./monitor.sh
# 选择 15（诊断 502 错误）
# 如果发现问题，选择 16（自动修复 502 错误）
```

### 示例 4：查看数据状态
```bash
./monitor.sh
# 选择 13（数据状态）
```

### 示例 5：LLM 解读 JSON 数据
```bash
# 方式 1：从文件读取
./monitor.sh
# 选择 17，然后输入文件路径

# 方式 2：从 stdin 输入
echo '{"status": "success", "data": {...}}' | ./monitor.sh
# 选择 17，直接按 Enter，然后粘贴 JSON
```

### 示例 6：实时监控面板
```bash
./monitor.sh
# 选择 18（实时监控面板）
# 按 Ctrl+C 退出
```

## 🔧 配置说明

脚本中的配置项（可在脚本开头修改）：

```bash
PROJECT_DIR="/home/ubuntu/GIBH-AGENT-V2"  # 项目目录
UPLOAD_DIR="${PROJECT_DIR}/data/uploads"   # 上传目录
RESULTS_DIR="${PROJECT_DIR}/results"      # 结果目录
LOG_DIR="${PROJECT_DIR}"                   # 日志目录
API_PORT=8028                              # Docker API 端口
DIRECT_PORT=8018                           # 直接运行端口
```

## 📝 注意事项

1. **权限要求**：
   - 需要 Docker 和 Docker Compose 权限
   - 某些操作可能需要 sudo（如删除镜像）

2. **API 密钥**：
   - LLM 解读功能需要 `SILICONFLOW_API_KEY` 环境变量
   - 或从 `docker-compose.yml` 中自动读取

3. **网络要求**：
   - 健康检查需要能够访问 `http://localhost:${API_PORT}`
   - LLM 解读需要能够访问 `https://api.siliconflow.cn`

4. **日志文件**：
   - 脚本会自动查找以下日志文件：
     - `gibh_agent.log`
     - `server.log`
     - `.cursor/debug.log`

## 🐛 故障排查

### 问题 1：无法获取 Docker 状态
- 检查 Docker 是否运行：`docker ps`
- 检查 Docker Compose 是否安装：`docker compose version`

### 问题 2：API 无响应
- 检查服务是否启动：选择菜单项 1
- 查看日志：选择菜单项 9 或 11
- 尝试重启：选择菜单项 4

### 问题 3：LLM 解读失败
- 检查 API 密钥是否正确
- 检查网络连接
- 查看错误信息

## 🔄 更新日志

- **v1.0** (2024-12-XX)
  - 初始版本
  - 整合所有现有脚本功能
  - 添加可视化输出
  - 添加 LLM JSON 解读功能
  - 添加实时监控面板

## 📞 支持

如有问题或建议，请查看项目文档或联系维护人员。

