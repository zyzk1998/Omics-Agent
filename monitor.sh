#!/bin/bash
# ============================================
# GIBH-AGENT-V2 可视化监控运维脚本
# 整合 Docker、服务、日志监控等功能
# ============================================

# 不使用 set -e，以便更好地处理错误和权限问题

# 颜色定义
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
MAGENTA='\033[0;35m'
CYAN='\033[0;36m'
WHITE='\033[1;37m'
NC='\033[0m' # No Color
BOLD='\033[1m'

# 项目配置
PROJECT_DIR="/home/ubuntu/GIBH-AGENT-V2"
UPLOAD_DIR="${PROJECT_DIR}/data/uploads"
RESULTS_DIR="${PROJECT_DIR}/results"
LOG_DIR="${PROJECT_DIR}"
API_PORT=8028
DIRECT_PORT=8018

# Docker 命令前缀（用于处理权限）
DOCKER_CMD_PREFIX=""

# Ollama 配置
OLLAMA_MODEL="qwen3-coder:30b"
OLLAMA_URL="http://localhost:11434"

# 检查命令是否存在
check_command() {
    if ! command -v $1 &> /dev/null; then
        echo -e "${RED}❌ $1 未安装${NC}"
        return 1
    fi
    return 0
}

# 检查 Docker 权限并设置命令前缀
check_docker_permission() {
    # 先尝试普通权限
    if docker info > /dev/null 2>&1; then
        DOCKER_CMD_PREFIX=""
        return 0
    fi
    
    # 检查是否在 docker 组中
    if groups | grep -q docker; then
        # 在 docker 组中但可能 Docker 未运行
        print_status "warning" "Docker 可能未运行，或需要重启会话以应用 docker 组权限"
        return 1
    fi
    
    # 需要 sudo 权限
    echo -e "${YELLOW}需要 sudo 权限访问 Docker${NC}"
    echo -e "${YELLOW}请输入 sudo 密码（如果需要），或按 Ctrl+C 取消：${NC}"
    
    # 测试 sudo 权限（使用 -S 从 stdin 读取密码，但这里先测试）
    if sudo -v 2>/dev/null; then
        DOCKER_CMD_PREFIX="sudo "
        print_status "ok" "已获取 sudo 权限"
        return 0
    else
        # 尝试使用 sudo -S 从 stdin 读取
        echo -e "${YELLOW}尝试使用 sudo...${NC}"
        if echo "" | sudo -S -v 2>/dev/null; then
            DOCKER_CMD_PREFIX="sudo "
            print_status "ok" "已配置 sudo 权限"
            return 0
        else
            print_status "error" "无法获取 sudo 权限"
            echo -e "${YELLOW}提示：可以手动运行 'sudo ./monitor.sh' 或添加用户到 docker 组${NC}"
            return 1
        fi
    fi
}

# 执行 Docker Compose 命令（自动处理权限）
docker_compose_cmd() {
    local cmd="$1"
    shift
    
    # 如果已设置前缀，直接使用
    if [ -n "${DOCKER_CMD_PREFIX}" ]; then
        if ${DOCKER_CMD_PREFIX}docker compose $cmd "$@" 2>/dev/null; then
            return 0
        elif ${DOCKER_CMD_PREFIX}docker-compose $cmd "$@" 2>/dev/null; then
            return 0
        else
            return 1
        fi
    fi
    
    # 尝试普通权限
    if docker compose $cmd "$@" 2>/dev/null; then
        return 0
    elif docker-compose $cmd "$@" 2>/dev/null; then
        return 0
    else
        return 1
    fi
}

# 执行 Docker 命令（自动处理权限）
docker_cmd() {
    local cmd="$1"
    shift
    
    if [ -z "${DOCKER_CMD_PREFIX}" ]; then
        docker $cmd "$@" 2>&1
    else
        ${DOCKER_CMD_PREFIX}docker $cmd "$@" 2>&1
    fi
}

# 打印分隔线
print_separator() {
    echo -e "${CYAN}${BOLD}============================================${NC}"
}

# 打印标题
print_title() {
    echo -e "\n${MAGENTA}${BOLD}$1${NC}"
    print_separator
}

# 打印状态
print_status() {
    if [ "$1" = "ok" ]; then
        echo -e "${GREEN}✅ $2${NC}"
    elif [ "$1" = "error" ]; then
        echo -e "${RED}❌ $2${NC}"
    elif [ "$1" = "warning" ]; then
        echo -e "${YELLOW}⚠️  $2${NC}"
    elif [ "$1" = "info" ]; then
        echo -e "${BLUE}ℹ️  $2${NC}"
    fi
}

# 打开浏览器
open_browser() {
    local url="http://localhost:${API_PORT}"
    local max_attempts=10
    local attempt=0
    
    echo ""
    echo -e "${CYAN}🔍 检查服务是否就绪...${NC}"
    
    # 等待服务启动
    while [ $attempt -lt $max_attempts ]; do
        if curl -s -o /dev/null -w "%{http_code}" "${url}" 2>/dev/null | grep -q "200\|301\|302"; then
            print_status "ok" "服务已就绪！"
            echo ""
            echo -e "${GREEN}${BOLD}🌐 正在打开浏览器...${NC}"
            echo -e "${WHITE}访问地址: ${url}${NC}"
            
            # 尝试使用不同的方式打开浏览器
            if command -v xdg-open &> /dev/null; then
                xdg-open "${url}" 2>/dev/null &
            elif command -v gnome-open &> /dev/null; then
                gnome-open "${url}" 2>/dev/null &
            elif command -v kde-open &> /dev/null; then
                kde-open "${url}" 2>/dev/null &
            elif [ -n "$DISPLAY" ] && command -v firefox &> /dev/null; then
                firefox "${url}" 2>/dev/null &
            elif [ -n "$DISPLAY" ] && command -v google-chrome &> /dev/null; then
                google-chrome "${url}" 2>/dev/null &
            elif [ -n "$DISPLAY" ] && command -v chromium-browser &> /dev/null; then
                chromium-browser "${url}" 2>/dev/null &
            else
                echo -e "${YELLOW}⚠️  无法自动打开浏览器，请手动访问: ${url}${NC}"
                return 1
            fi
            
            sleep 1
            print_status "ok" "浏览器已打开！"
            return 0
        fi
        
        attempt=$((attempt + 1))
        echo -e "${YELLOW}⏳ 等待服务启动... (${attempt}/${max_attempts})${NC}"
        sleep 2
    done
    
    print_status "warning" "服务可能尚未完全启动，请稍后手动访问: ${url}"
    return 1
}

# ============================================
# 1. Docker 服务管理
# ============================================

docker_status() {
    print_title "🐳 Docker 服务状态"
    
    if ! check_command docker; then
        return 1
    fi
    
    # 检查 Docker 权限（不强制要求，只是尝试）
    check_docker_permission 2>/dev/null || true
    
    echo -e "${WHITE}容器状态：${NC}"
    docker_compose_cmd ps || {
        print_status "error" "无法获取容器状态，请检查 Docker Compose"
        return 1
    }
    
    echo ""
    echo -e "${WHITE}容器资源使用：${NC}"
    $(docker_cmd stats --no-stream --format "table {{.Container}}\t{{.CPUPerc}}\t{{.MemUsage}}\t{{.NetIO}}") | head -10 || echo "无法获取资源使用情况"
}

docker_start() {
    print_title "🚀 启动 Docker 服务"
    
    # 检查 Docker 权限
    if ! check_docker_permission; then
        return 1
    fi
    
    # 创建必要的目录
    echo "📁 创建必要的目录..."
    mkdir -p ${UPLOAD_DIR} ${RESULTS_DIR} ${PROJECT_DIR}/data/redis
    
    # 启动服务
    echo "🚀 启动所有服务..."
    docker_compose_cmd up -d || {
        print_status "error" "启动服务失败"
        return 1
    }
    
    echo ""
    echo "⏳ 等待服务启动（5秒）..."
    sleep 5
    
    docker_status
    
    # 自动打开浏览器
    open_browser
}

docker_stop() {
    print_title "🛑 停止 Docker 服务"
    
    # 检查 Docker 权限
    if ! check_docker_permission; then
        return 1
    fi
    
    docker_compose_cmd down || {
        print_status "error" "停止服务失败"
        return 1
    }
    print_status "ok" "服务已停止"
}

docker_restart() {
    print_title "🔄 重启 Docker 服务"
    docker_stop
    sleep 2
    docker_start
}

docker_rebuild() {
    print_title "🔨 重新构建 Docker 镜像"
    
    read -p "是否删除旧镜像以强制重新构建？(y/N): " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        echo "🗑️  删除旧镜像..."
        docker rmi gibh-v2-api:latest 2>/dev/null || echo "镜像不存在，跳过删除"
    fi
    
    echo "🔨 重新构建镜像..."
    docker compose build --no-cache 2>/dev/null || docker-compose build --no-cache 2>/dev/null
    
    print_status "ok" "构建完成"
    
    read -p "是否立即启动服务？(Y/n): " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Nn]$ ]]; then
        docker_start
        # docker_start 内部已经会调用 open_browser
    fi
}

# ============================================
# 2. 服务健康检查
# ============================================

health_check() {
    print_title "🏥 服务健康检查"
    
    # 检查 Docker 权限（不强制要求）
    check_docker_permission 2>/dev/null || true
    
    # 检查 Docker 服务
    echo -e "${WHITE}1. Docker 容器状态：${NC}"
    if $(docker_cmd ps) 2>/dev/null | grep -q gibh_v2_api; then
        print_status "ok" "API 服务器容器运行中"
    else
        print_status "warning" "API 服务器容器未运行（可能使用直接运行模式）"
    fi
    
    if $(docker_cmd ps) 2>/dev/null | grep -q gibh_v2_redis; then
        print_status "ok" "Redis 容器运行中"
    else
        print_status "warning" "Redis 容器未运行（可能使用直接运行模式）"
    fi
    
    if $(docker_cmd ps) 2>/dev/null | grep -q gibh_v2_worker; then
        print_status "ok" "Worker 容器运行中"
    else
        print_status "warning" "Worker 容器未运行（可选服务）"
    fi
    
    echo ""
    echo -e "${WHITE}2. API 服务器响应：${NC}"
    api_found=false
    
    # 检查 Docker 端口 (8028)
    if curl -s -o /dev/null -w "%{http_code}" http://localhost:${API_PORT}/ 2>/dev/null | grep -q "200\|301\|302"; then
        print_status "ok" "Docker API 服务器响应正常 (http://localhost:${API_PORT})"
        api_found=true
    elif curl -s -o /dev/null -w "%{http_code}" http://localhost:${API_PORT}/api/docs 2>/dev/null | grep -q "200\|301\|302"; then
        print_status "ok" "Docker API 服务器响应正常 (http://localhost:${API_PORT}/api/docs)"
        api_found=true
    fi
    
    # 检查直接运行端口 (8018)
    if curl -s -o /dev/null -w "%{http_code}" http://localhost:${DIRECT_PORT}/ 2>/dev/null | grep -q "200\|301\|302"; then
        print_status "ok" "直接运行服务器响应正常 (http://localhost:${DIRECT_PORT})"
        api_found=true
    elif curl -s -o /dev/null -w "%{http_code}" http://localhost:${DIRECT_PORT}/api/docs 2>/dev/null | grep -q "200\|301\|302"; then
        print_status "ok" "直接运行服务器响应正常 (http://localhost:${DIRECT_PORT}/api/docs)"
        api_found=true
    fi
    
    if [ "$api_found" = false ]; then
        print_status "error" "API 服务器无响应"
        echo -e "${YELLOW}已检查端口：${NC}"
        echo "  - ${API_PORT} (Docker)"
        echo "  - ${DIRECT_PORT} (直接运行)"
    fi
    
    echo ""
    echo -e "${WHITE}3. 端口监听：${NC}"
    if netstat -tlnp 2>/dev/null | grep -q ":${API_PORT} "; then
        print_status "ok" "端口 ${API_PORT} 正在监听"
    elif ss -tlnp 2>/dev/null | grep -q ":${API_PORT} "; then
        print_status "ok" "端口 ${API_PORT} 正在监听"
    else
        print_status "warning" "端口 ${API_PORT} 未监听（可能在使用 Docker 网络）"
    fi
    
    echo ""
    echo -e "${WHITE}4. 直接运行服务器状态：${NC}"
    if pgrep -f "python.*server.py" > /dev/null; then
        print_status "info" "检测到直接运行的服务器进程"
        ps aux | grep -E "python.*server.py" | grep -v grep | head -3
    else
        print_status "info" "未检测到直接运行的服务器"
    fi
}

# ============================================
# 3. 日志监控
# ============================================

logs_realtime() {
    print_title "📡 实时日志监控"
    
    service=${1:-api-server}
    echo -e "${WHITE}监控服务: ${service}${NC}"
    
    # 检查 Docker 权限
    if ! check_docker_permission; then
        return 1
    fi
    
    # 检查容器是否存在
    container_name="gibh_v2_${service//-/_}"
    if ! $(docker_cmd ps -a --format "{{.Names}}") | grep -q "^${container_name}$"; then
        print_status "warning" "容器 ${container_name} 不存在"
        echo -e "${YELLOW}提示：容器可能未启动，请先启动 Docker 服务（菜单项 2）${NC}"
        return 1
    fi
    
    echo -e "${YELLOW}按 Ctrl+C 退出${NC}\n"
    
    # 获取实时日志
    docker_compose_cmd logs -f ${service} || {
        print_status "error" "无法获取日志，请检查 Docker Compose"
        return 1
    }
}

logs_recent() {
    print_title "📋 最近日志"
    
    service=${1:-api-server}
    lines=${2:-50}
    
    echo -e "${WHITE}服务: ${service} | 行数: ${lines}${NC}\n"
    
    # 检查 Docker 权限
    if ! check_docker_permission; then
        return 1
    fi
    
    # 检查容器是否存在
    container_name="gibh_v2_${service//-/_}"
    if ! $(docker_cmd ps -a --format "{{.Names}}") | grep -q "^${container_name}$"; then
        print_status "warning" "容器 ${container_name} 不存在"
        echo -e "${YELLOW}提示：容器可能未启动，请先启动 Docker 服务（菜单项 2）${NC}"
        return 1
    fi
    
    # 获取日志
    docker_compose_cmd logs --tail ${lines} ${service} | head -100 || {
        print_status "error" "无法获取日志"
        return 1
    }
}

logs_all() {
    print_title "📚 所有服务日志摘要"
    
    # 检查 Docker 权限
    if ! check_docker_permission; then
        return 1
    fi
    
    services=("api-server" "worker" "redis")
    
    for service in "${services[@]}"; do
        echo -e "\n${CYAN}${BOLD}--- ${service} (最近 20 行) ---${NC}"
        
        # 检查容器是否存在
        container_name="gibh_v2_${service//-/_}"
        if ! $(docker_cmd ps -a --format "{{.Names}}") | grep -q "^${container_name}$"; then
            echo -e "${YELLOW}  容器不存在或未运行${NC}"
            continue
        fi
        
        # 获取日志
        docker_compose_cmd logs --tail 20 ${service} | head -25 || echo -e "${YELLOW}  无法获取日志${NC}"
    done
}

logs_errors() {
    print_title "🔍 错误日志分析"
    
    service=${1:-api-server}
    lines=${2:-100}
    
    echo -e "${WHITE}在 ${service} 的最近 ${lines} 行日志中查找错误：${NC}\n"
    
    # 检查 Docker 权限
    if ! check_docker_permission; then
        return 1
    fi
    
    # 获取日志
    temp_log=$(mktemp)
    if docker_compose_cmd logs --tail ${lines} ${service} > "${temp_log}" 2>&1; then
        all_logs=$(cat "${temp_log}")
    else
        print_status "error" "无法获取日志"
        rm -f "${temp_log}"
        return 1
    fi
    
    # 过滤错误日志
    error_logs=$(echo "${all_logs}" | grep -i -E "error|exception|failed|traceback|❌" || echo "")
    
    if [ -n "${error_logs}" ]; then
        echo -e "${RED}${BOLD}错误日志（已过滤）：${NC}"
        echo "${error_logs}"
    else
        print_status "info" "未发现错误日志"
    fi
    
    echo ""
    echo -e "${YELLOW}是否显示完整的原始日志（最近 ${lines} 行）？(y/N): ${NC}"
    read -n 1 -r show_full
    echo
    if [[ $show_full =~ ^[Yy]$ ]]; then
        echo -e "\n${WHITE}${BOLD}完整原始日志：${NC}\n"
        echo "${all_logs}"
    fi
    
    rm -f "${temp_log}"
}

logs_file() {
    print_title "📄 本地日志文件"
    
    log_files=(
        "${LOG_DIR}/gibh_agent.log"
        "${LOG_DIR}/server.log"
        "${LOG_DIR}/.cursor/debug.log"
    )
    
    for log_file in "${log_files[@]}"; do
        if [ -f "${log_file}" ]; then
            size=$(du -h "${log_file}" | cut -f1)
            lines=$(wc -l < "${log_file}" 2>/dev/null || echo "0")
            echo -e "${WHITE}${log_file}${NC}"
            echo -e "  大小: ${size} | 行数: ${lines}"
            
            if [ "${lines}" -gt 0 ]; then
                echo -e "  最后 5 行："
                tail -5 "${log_file}" | sed 's/^/  /'
                
                # 如果是 debug.log，使用 Ollama 解读
                if [[ "${log_file}" == *"debug.log" ]]; then
                    echo ""
                    echo -e "${YELLOW}是否使用 Ollama 解读所有 JSON 记录？(y/N): ${NC}"
                    read -n 1 -r use_ollama
                    echo
                    if [[ $use_ollama =~ ^[Yy]$ ]]; then
                        interpret_debug_log_with_ollama "${log_file}"
                    fi
                fi
            fi
            echo ""
        fi
    done
}

# 使用 Ollama 解读 debug.log
interpret_debug_log_with_ollama() {
    local log_file="$1"
    
    print_title "🤖 使用 Ollama 解读 debug.log"
    
    # 检查 Ollama 是否可用
    if ! command -v ollama &> /dev/null; then
        print_status "error" "Ollama 未安装或不在 PATH 中"
        return 1
    fi
    
    # 检查模型是否存在
    if ! ollama list 2>/dev/null | grep -q "${OLLAMA_MODEL}"; then
        print_status "error" "模型 ${OLLAMA_MODEL} 不存在"
        echo -e "${YELLOW}可用模型：${NC}"
        ollama list 2>/dev/null || echo "无法列出模型"
        return 1
    fi
    
    # 读取所有 JSON 记录
    echo -e "${WHITE}正在读取 JSON 记录...${NC}"
    json_content=$(cat "${log_file}")
    
    if [ -z "${json_content}" ]; then
        print_status "warning" "日志文件为空"
        return 1
    fi
    
    # 统计记录数
    record_count=$(echo "${json_content}" | grep -c "^{" || echo "0")
    echo -e "${WHITE}找到 ${record_count} 条 JSON 记录${NC}\n"
    
    # 构建提示词（包含项目背景信息）
    prompt="## 项目背景

你正在分析 **GIBH-AGENT-V2** 项目的调试日志。这是一个基于多模态大模型的生物信息学分析智能体平台，主要功能包括：

1. **多智能体架构**：
   - RouterAgent：路由智能体，识别用户查询的组学类型（转录组、基因组、代谢组等）
   - Domain Agents：领域智能体（如 RNAAgent、MetabolomicsAgent），处理特定类型的分析任务
   - 工具类：生成分析脚本（如 Scanpy、Cell Ranger）

2. **技术栈**：
   - FastAPI 服务器（server.py）
   - 使用 DeepSeek API（硅基流动）作为 LLM
   - 支持文件上传和工作流执行

3. **工作流程**：
   - 用户通过 Web 界面输入自然语言查询（如\"我要做代谢组分析\"）
   - 系统通过 RouterAgent 识别意图，路由到对应的 Domain Agent
   - Domain Agent 生成工作流配置（workflow_config）或直接回复（chat）
   - 执行分析任务并返回结果

## 调试日志说明

以下 JSON 日志记录的是 FastAPI 服务器（server.py）处理用户查询时的关键执行点：

- **location**: 代码位置（文件:行号），如 \"server.py:1161\" 表示在 server.py 第 1161 行
- **message**: 调试消息，描述当前执行的操作
  - \"Before process_query\": 调用智能体处理查询之前
  - \"After process_query\": 处理查询之后，显示返回结果类型
  - \"chat_endpoint entry\": 进入聊天接口
- **data**: 关键数据
  - \"query\": 用户输入的查询内容
  - \"uploaded_files_count\": 上传的文件数量
  - \"result_type\": 返回结果类型（\"dict\" 表示字典）
  - \"result_keys\": 返回字典的键列表（如 [\"type\", \"workflow_data\", \"file_paths\", \"routing_info\"]）
  - \"result_type_value\": 结果类型值（\"workflow_config\" 表示工作流配置，\"chat\" 表示聊天回复）
- **timestamp**: Unix 时间戳（毫秒）
- **sessionId/runId/hypothesisId**: 会话和运行标识，用于追踪调试流程

## 你的任务

作为调试助手，请直接分析以下 JSON 日志，重点关注：

1. **执行流程问题**：
   - 是否有异常的执行路径？
   - 是否有重复的请求？
   - 流程是否按预期执行？

2. **数据异常**：
   - 返回结果类型是否符合预期？
   - 是否有数据丢失或不一致？
   - 上传文件处理是否正确？

3. **性能问题**：
   - 请求处理时间是否过长？
   - 是否有阻塞或延迟？

4. **潜在 Bug**：
   - 是否有错误或异常？
   - 是否有逻辑问题？

请直接开始分析，不需要解释 JSON 结构，直接指出问题和建议。

---

JSON 日志内容：
\`\`\`json
${json_content}
\`\`\`"
    
    # 使用临时文件存储提示词（避免命令行长度限制）
    temp_prompt=$(mktemp)
    echo "${prompt}" > "${temp_prompt}"
    
    echo -e "${WHITE}正在调用 Ollama (${OLLAMA_MODEL}) 解读...${NC}"
    echo -e "${YELLOW}这可能需要一些时间，请耐心等待...${NC}\n"
    
    # 调用 Ollama（使用临时文件）
    # 使用 ollama run 的 stdin 输入方式
    response=$(cat "${temp_prompt}" | ollama run ${OLLAMA_MODEL} 2>&1)
    exit_code=$?
    
    # 清理临时文件
    rm -f "${temp_prompt}"
    
    if [ ${exit_code} -eq 0 ] && [ -n "${response}" ]; then
        echo -e "\n${GREEN}${BOLD}Ollama 解读结果：${NC}\n"
        echo "${response}"
    else
        print_status "error" "Ollama 调用失败 (退出码: ${exit_code})"
        if [ -n "${response}" ]; then
            echo -e "${YELLOW}错误信息：${NC}"
            echo "${response}"
        fi
        echo -e "\n${YELLOW}提示：${NC}"
        echo "  1. 检查 Ollama 服务是否运行: ollama serve"
        echo "  2. 检查模型是否存在: ollama list"
        echo "  3. 尝试手动运行: ollama run ${OLLAMA_MODEL}"
    fi
}

# ============================================
# 4. 数据监控
# ============================================

data_status() {
    print_title "📊 数据状态监控"
    
    # 上传文件统计
    echo -e "${WHITE}1. 上传文件统计：${NC}"
    if [ -d "${UPLOAD_DIR}" ]; then
        upload_count=$(find "${UPLOAD_DIR}" -type f ! -name "*.meta.json" | wc -l)
        upload_size=$(du -sh "${UPLOAD_DIR}" 2>/dev/null | cut -f1 || echo "0")
        echo -e "  文件数: ${upload_count}"
        echo -e "  总大小: ${upload_size}"
        
        if [ "${upload_count}" -gt 0 ]; then
            echo -e "  最近上传的文件："
            find "${UPLOAD_DIR}" -type f ! -name "*.meta.json" -printf "  %Tb %Td %TH:%TM  %s  %f\n" | sort -r | head -5
        fi
    else
        print_status "warning" "上传目录不存在"
    fi
    
    echo ""
    echo -e "${WHITE}2. 结果文件统计：${NC}"
    if [ -d "${RESULTS_DIR}" ]; then
        result_count=$(find "${RESULTS_DIR}" -type f | wc -l)
        result_size=$(du -sh "${RESULTS_DIR}" 2>/dev/null | cut -f1 || echo "0")
        echo -e "  文件数: ${result_count}"
        echo -e "  总大小: ${result_size}"
        
        if [ "${result_count}" -gt 0 ]; then
            echo -e "  最近生成的结果："
            find "${RESULTS_DIR}" -type f -printf "  %Tb %Td %TH:%TM  %s  %f\n" | sort -r | head -5
        fi
    else
        print_status "warning" "结果目录不存在"
    fi
    
    echo ""
    echo -e "${WHITE}3. 磁盘使用情况：${NC}"
    df -h "${PROJECT_DIR}" | tail -1 | awk '{print "  总空间: " $2 " | 已用: " $3 " (" $5 ") | 可用: " $4}'
    
    echo ""
    echo -e "${WHITE}4. 目录大小排序（前 10）：${NC}"
    du -sh "${PROJECT_DIR}"/* 2>/dev/null | sort -hr | head -10 | sed 's/^/  /'
}

data_cleanup() {
    print_title "🧹 数据清理"
    
    echo -e "${YELLOW}警告：此操作将删除数据文件${NC}\n"
    
    echo "1. 清理上传文件（保留元数据）"
    echo "2. 清理结果文件"
    echo "3. 清理所有数据"
    echo "4. 取消"
    
    read -p "请选择 (1-4): " choice
    
    case $choice in
        1)
            read -p "确认删除上传文件？(y/N): " -n 1 -r
            echo
            if [[ $REPLY =~ ^[Yy]$ ]]; then
                find "${UPLOAD_DIR}" -type f ! -name "*.meta.json" -delete
                print_status "ok" "上传文件已清理"
            fi
            ;;
        2)
            read -p "确认删除结果文件？(y/N): " -n 1 -r
            echo
            if [[ $REPLY =~ ^[Yy]$ ]]; then
                rm -rf "${RESULTS_DIR}"/*
                print_status "ok" "结果文件已清理"
            fi
            ;;
        3)
            read -p "确认删除所有数据？(y/N): " -n 1 -r
            echo
            if [[ $REPLY =~ ^[Yy]$ ]]; then
                find "${UPLOAD_DIR}" -type f ! -name "*.meta.json" -delete
                rm -rf "${RESULTS_DIR}"/*
                print_status "ok" "所有数据已清理"
            fi
            ;;
        4)
            print_status "info" "已取消"
            ;;
        *)
            print_status "error" "无效选择"
            ;;
    esac
}

# ============================================
# 5. 错误诊断和修复
# ============================================

diagnose_502() {
    print_title "🔍 502 错误诊断"
    
    echo -e "${WHITE}1. 容器状态检查：${NC}"
    docker compose ps
    
    echo ""
    echo -e "${WHITE}2. API 服务器日志（最近 50 行）：${NC}"
    docker compose logs --tail 50 api-server 2>/dev/null | tail -20
    
    echo ""
    echo -e "${WHITE}3. Worker 日志（最近 30 行）：${NC}"
    docker compose logs --tail 30 worker 2>/dev/null | tail -20 || echo "Worker 未运行"
    
    echo ""
    echo -e "${WHITE}4. 检查依赖：${NC}"
    if grep -q "paramiko" "${PROJECT_DIR}/requirements.txt"; then
        print_status "ok" "paramiko 已在 requirements.txt 中"
    else
        print_status "error" "paramiko 缺失"
    fi
    
    echo ""
    echo -e "${WHITE}5. 网络连接测试：${NC}"
    if docker compose exec -T api-server curl -s http://localhost:${API_PORT}/api/docs > /dev/null 2>&1; then
        print_status "ok" "容器内部 API 响应正常"
    else
        print_status "error" "容器内部 API 无响应"
    fi
    
    echo ""
    print_status "info" "诊断完成，请查看上述信息"
}

fix_502() {
    print_title "🔧 自动修复 502 错误"
    
    echo "1️⃣ 停止所有容器..."
    docker compose down
    
    echo ""
    echo "2️⃣ 检查 requirements.txt..."
    if ! grep -q "paramiko" "${PROJECT_DIR}/requirements.txt"; then
        print_status "warning" "paramiko 缺失，正在添加..."
        echo "paramiko>=3.0.0" >> "${PROJECT_DIR}/requirements.txt"
        print_status "ok" "已添加 paramiko"
    fi
    
    echo ""
    echo "3️⃣ 重新构建镜像..."
    docker compose build --no-cache
    
    echo ""
    echo "4️⃣ 启动服务..."
    docker compose up -d
    
    echo ""
    echo "5️⃣ 等待服务启动（10秒）..."
    sleep 10
    
    echo ""
    echo "6️⃣ 检查服务状态..."
    docker compose ps
    
    echo ""
    echo "7️⃣ 测试 API..."
    if curl -s -f http://localhost:${API_PORT}/api/docs > /dev/null 2>&1; then
        print_status "ok" "修复成功！API 服务器响应正常"
    else
        print_status "error" "修复后仍无法访问，请查看日志："
        echo "  docker compose logs api-server"
    fi
}

# ============================================
# 6. JSON 数据查看（原始数据）
# ============================================

interpret_json() {
    print_title "📄 JSON 数据查看（原始数据）"
    
    # 获取 JSON 内容
    if [ -z "$1" ]; then
        echo -e "${YELLOW}请输入 JSON 文件路径，或直接按 Enter 从 stdin 输入（Ctrl+D 结束）：${NC}"
        read json_file
        if [ -z "$json_file" ]; then
            echo -e "${WHITE}请输入 JSON 内容（Ctrl+D 结束）：${NC}"
            json_content=$(cat)
        else
            if [ ! -f "$json_file" ]; then
                print_status "error" "文件不存在: $json_file"
                return 1
            fi
            json_content=$(cat "$json_file")
        fi
    else
        if [ ! -f "$1" ]; then
            print_status "error" "文件不存在: $1"
            return 1
        fi
        json_content=$(cat "$1")
    fi
    
    # 验证并格式化 JSON
    echo -e "${WHITE}正在格式化 JSON 数据...${NC}\n"
    
    formatted_json=$(echo "${json_content}" | python3 -m json.tool 2>&1)
    
    if [ $? -eq 0 ]; then
        echo -e "${GREEN}${BOLD}格式化后的 JSON 数据：${NC}\n"
        echo "${formatted_json}"
    else
        print_status "warning" "JSON 格式验证失败，显示原始数据："
        echo -e "\n${YELLOW}原始数据：${NC}\n"
        echo "${json_content}"
    fi
    
    # 显示统计信息
    echo -e "\n${CYAN}${BOLD}数据统计：${NC}"
    char_count=$(echo "${json_content}" | wc -c)
    line_count=$(echo "${json_content}" | wc -l)
    echo -e "  字符数: ${char_count}"
    echo -e "  行数: ${line_count}"
    
    # 尝试提取基本信息（如果可能）
    if command -v python3 &> /dev/null; then
        echo -e "\n${CYAN}${BOLD}数据结构信息：${NC}"
        # 使用临时文件避免转义问题
        temp_json=$(mktemp)
        echo "${json_content}" > "${temp_json}"
        
        python3 << EOF
import json
import sys

try:
    with open("${temp_json}", "r", encoding="utf-8") as f:
        data = json.load(f)
    
    def analyze_json(obj, indent=0):
        prefix = "  " * indent
        if isinstance(obj, dict):
            print(f"{prefix}类型: 对象 (dict)")
            print(f"{prefix}键数量: {len(obj)}")
            if len(obj) > 0:
                print(f"{prefix}键列表:")
                for key in list(obj.keys())[:10]:  # 只显示前10个键
                    print(f"{prefix}  - {key}")
                if len(obj) > 10:
                    print(f"{prefix}  ... (还有 {len(obj) - 10} 个键)")
        elif isinstance(obj, list):
            print(f"{prefix}类型: 数组 (list)")
            print(f"{prefix}元素数量: {len(obj)}")
            if len(obj) > 0:
                print(f"{prefix}第一个元素类型: {type(obj[0]).__name__}")
        elif isinstance(obj, str):
            print(f"{prefix}类型: 字符串 (str)")
            print(f"{prefix}长度: {len(obj)}")
        elif isinstance(obj, (int, float)):
            print(f"{prefix}类型: 数字 ({type(obj).__name__})")
            print(f"{prefix}值: {obj}")
        elif isinstance(obj, bool):
            print(f"{prefix}类型: 布尔值 (bool)")
            print(f"{prefix}值: {obj}")
        elif obj is None:
            print(f"{prefix}类型: 空值 (None)")
        else:
            print(f"{prefix}类型: {type(obj).__name__}")
    
    analyze_json(data)
except Exception as e:
    print(f"  无法解析 JSON: {e}")
EOF
        rm -f "${temp_json}"
    fi
}

# ============================================
# 7. 综合监控面板
# ============================================

monitor_dashboard() {
    print_title "📊 综合监控面板"
    
    while true; do
        clear
        echo -e "${MAGENTA}${BOLD}╔════════════════════════════════════════════════╗${NC}"
        echo -e "${MAGENTA}${BOLD}║   GIBH-AGENT-V2 实时监控面板                  ║${NC}"
        echo -e "${MAGENTA}${BOLD}╚════════════════════════════════════════════════╝${NC}"
        echo ""
        
        # 服务状态
        echo -e "${CYAN}${BOLD}🐳 Docker 服务状态：${NC}"
        if check_docker_permission 2>/dev/null; then
            docker_compose_cmd ps 2>/dev/null | tail -n +2 | while read line; do
                if echo "$line" | grep -q "Up"; then
                    echo -e "  ${GREEN}✅ $line${NC}"
                else
                    echo -e "  ${RED}❌ $line${NC}"
                fi
            done || echo -e "  ${YELLOW}⚠️  无法获取 Docker 状态（可能使用直接运行模式）${NC}"
        else
            echo -e "  ${YELLOW}⚠️  Docker 未运行或无权限（可能使用直接运行模式）${NC}"
        fi
        
        echo ""
        
        # API 健康（检查多个端点）
        echo -e "${CYAN}${BOLD}🏥 API 健康状态：${NC}"
        api_healthy=false
        
        # 检查 Docker 端口 (8028)
        if curl -s -o /dev/null -w "%{http_code}" http://localhost:${API_PORT}/ 2>/dev/null | grep -q "200\|301\|302"; then
            echo -e "  ${GREEN}✅ Docker API 服务器正常 (http://localhost:${API_PORT})${NC}"
            api_healthy=true
        elif curl -s -o /dev/null -w "%{http_code}" http://localhost:${API_PORT}/api/docs 2>/dev/null | grep -q "200\|301\|302"; then
            echo -e "  ${GREEN}✅ Docker API 服务器正常 (http://localhost:${API_PORT}/api/docs)${NC}"
            api_healthy=true
        fi
        
        # 检查直接运行端口 (8018)
        if curl -s -o /dev/null -w "%{http_code}" http://localhost:${DIRECT_PORT}/ 2>/dev/null | grep -q "200\|301\|302"; then
            echo -e "  ${GREEN}✅ 直接运行服务器正常 (http://localhost:${DIRECT_PORT})${NC}"
            api_healthy=true
        elif curl -s -o /dev/null -w "%{http_code}" http://localhost:${DIRECT_PORT}/api/docs 2>/dev/null | grep -q "200\|301\|302"; then
            echo -e "  ${GREEN}✅ 直接运行服务器正常 (http://localhost:${DIRECT_PORT}/api/docs)${NC}"
            api_healthy=true
        fi
        
        if [ "$api_healthy" = false ]; then
            echo -e "  ${RED}❌ API 服务器无响应${NC}"
            echo -e "  ${YELLOW}提示：检查端口 ${API_PORT} (Docker) 或 ${DIRECT_PORT} (直接运行)${NC}"
        fi
        
        echo ""
        
        # 数据统计
        echo -e "${CYAN}${BOLD}📊 数据统计：${NC}"
        if [ -d "${UPLOAD_DIR}" ]; then
            upload_count=$(find "${UPLOAD_DIR}" -type f ! -name "*.meta.json" 2>/dev/null | wc -l)
            echo -e "  上传文件: ${upload_count} 个"
        fi
        if [ -d "${RESULTS_DIR}" ]; then
            result_count=$(find "${RESULTS_DIR}" -type f 2>/dev/null | wc -l)
            echo -e "  结果文件: ${result_count} 个"
        fi
        
        echo ""
        
        # 最近错误
        echo -e "${CYAN}${BOLD}🔍 最近错误（API 服务器，最后 5 条）：${NC}"
        docker compose logs --tail 100 api-server 2>/dev/null | grep -i -E "error|exception|failed|❌" | tail -5 || echo "  无错误"
        
        echo ""
        echo -e "${YELLOW}按 Ctrl+C 退出监控${NC}"
        echo -e "${YELLOW}每 5 秒自动刷新...${NC}"
        
        sleep 5
    done
}

# ============================================
# 7. API 请求/响应调试
# ============================================

api_debug_realtime() {
    print_title "🔍 实时监控 API 请求/响应"
    
    echo -e "${WHITE}正在监控 API 请求和响应...${NC}"
    echo -e "${YELLOW}按 Ctrl+C 退出${NC}\n"
    
    # 监控多个日志源
    (
        # 监控服务器日志（包含请求/响应）
        tail -f "${LOG_DIR}/gibh_agent.log" 2>/dev/null | grep -E "(收到|返回|请求|响应|POST|GET|JSON|result|response)" --line-buffered &
        PID1=$!
        
        # 监控 debug.log（包含详细的请求/响应数据）
        tail -f "${LOG_DIR}/.cursor/debug.log" 2>/dev/null | while read line; do
            if echo "$line" | python3 -m json.tool > /dev/null 2>&1; then
                echo -e "${CYAN}[DEBUG]${NC} $line" | python3 -m json.tool 2>/dev/null || echo "$line"
            fi
        done &
        PID2=$!
        
        # 监控 Docker 日志（如果有）
        if check_docker_permission 2>/dev/null; then
            docker_compose_cmd logs -f api-server 2>/dev/null | grep -E "(收到|返回|请求|响应|POST|GET|JSON|result|response)" --line-buffered &
            PID3=$!
        fi
        
        # 等待用户中断
        trap "kill $PID1 $PID2 $PID3 2>/dev/null; exit" INT TERM
        wait
    )
}

api_debug_recent() {
    print_title "📋 最近的 API 请求/响应日志"
    
    echo -e "${WHITE}分析最近的 API 交互...${NC}\n"
    
    # 从 debug.log 提取请求/响应信息
    debug_log="${LOG_DIR}/.cursor/debug.log"
    
    if [ -f "${debug_log}" ]; then
        echo -e "${CYAN}${BOLD}1. 完整的请求-响应流程：${NC}\n"
        
        # 提取最近的完整请求-响应对
        python3 << EOF
import json
import sys
from datetime import datetime

try:
    with open("${debug_log}", "r", encoding="utf-8") as f:
        lines = f.readlines()
    
    # 找到最近的请求-响应对
    entries = []
    for line in lines[-50:]:  # 只检查最后50行
        line = line.strip()
        if not line:
            continue
        try:
            entry = json.loads(line)
            if entry.get("message") in ["chat_endpoint entry", "Before process_query", "After process_query"]:
                entries.append(entry)
        except:
            continue
    
    if not entries:
        print("  未找到请求/响应数据")
        sys.exit(0)
    
    # 按时间戳排序
    entries.sort(key=lambda x: x.get("timestamp", 0))
    
    # 显示最近的3个完整流程
    for i, entry in enumerate(entries[-6:], 1):
        msg = entry.get("message", "")
        data = entry.get("data", {})
        ts = entry.get("timestamp", 0)
        
        if ts:
            dt = datetime.fromtimestamp(ts / 1000)
            time_str = dt.strftime("%H:%M:%S.%f")[:-3]
        else:
            time_str = "N/A"
        
        print(f"\n{'='*60}")
        print(f"📌 流程 #{i} - {time_str}")
        print(f"📍 位置: {entry.get('location', 'N/A')}")
        print(f"💬 消息: {msg}")
        print(f"{'='*60}")
        
        if msg == "chat_endpoint entry":
            print("📥 请求数据（前端发送）：")
            print(f"  - 消息: {data.get('req_message', 'N/A')[:100]}")
            print(f"  - Agent 状态: {data.get('agent_is_none', 'N/A')}")
        
        elif msg == "Before process_query":
            print("🔄 处理前数据：")
            print(f"  - 查询内容: {data.get('query', 'N/A')[:100]}")
            print(f"  - 上传文件数: {data.get('uploaded_files_count', 'N/A')}")
            print(f"  - 测试数据集: {data.get('test_dataset_id', 'N/A')}")
        
        elif msg == "After process_query":
            print("📤 处理后数据（准备返回给前端）：")
            print(f"  - 结果类型: {data.get('result_type', 'N/A')}")
            if data.get('result_type') == 'dict':
                print(f"  - 返回键: {data.get('result_keys', [])}")
                print(f"  - 结果值类型: {data.get('result_type_value', 'N/A')}")
                if data.get('result_type_value') == 'workflow_config':
                    print("  ✅ 这是工作流配置，会返回完整的 workflow_data")
                elif data.get('result_type_value') == 'chat':
                    print("  ✅ 这是聊天回复，会返回流式响应")
            print(f"  - 完整数据: {json.dumps(data, ensure_ascii=False, indent=2)}")
        
        print()
    
    print(f"\n💡 提示：")
    print(f"  - 'chat_endpoint entry': 前端发送的原始请求")
    print(f"  - 'Before process_query': 处理前的数据（已解析）")
    print(f"  - 'After process_query': 处理后的数据（准备返回给前端）")
    print(f"  - 如果 result_type_value 是 'workflow_config'，会返回 JSON")
    print(f"  - 如果 result_type_value 是 'chat'，会返回流式文本")

except Exception as e:
    print(f"  解析错误: {e}")
    import traceback
    traceback.print_exc()
EOF
    fi
    
    # 从服务器日志提取请求/响应
    server_log="${LOG_DIR}/gibh_agent.log"
    if [ -f "${server_log}" ]; then
        echo -e "\n${CYAN}${BOLD}2. 服务器日志中的请求/响应摘要：${NC}\n"
        
        echo -e "${WHITE}最近的请求：${NC}"
        grep -E "(收到聊天请求|处理文件|返回|result|response|JSONResponse|StreamingResponse)" "${server_log}" | tail -20 | sed 's/^/  /'
    fi
    
    # 从 Docker 日志提取（如果有）
    if check_docker_permission 2>/dev/null; then
        echo -e "\n${CYAN}${BOLD}3. Docker 容器日志中的请求/响应：${NC}\n"
        docker_compose_cmd logs --tail 50 api-server 2>/dev/null | grep -E "(收到|返回|请求|响应|POST|GET|JSON|JSONResponse|StreamingResponse)" | tail -20 | sed 's/^/  /'
    fi
    
    echo ""
    echo -e "${YELLOW}是否使用 Ollama 分析这些请求/响应数据？(y/N): ${NC}"
    read -n 1 -r use_ollama
    echo
    if [[ $use_ollama =~ ^[Yy]$ ]]; then
        if [ -f "${debug_log}" ]; then
            # 提取最近的完整请求-响应对
            recent_entries=$(tail -30 "${debug_log}" | grep -E "(chat_endpoint entry|Before process_query|After process_query)")
            if [ -n "${recent_entries}" ]; then
                temp_data=$(mktemp)
                echo "${recent_entries}" > "${temp_data}"
                interpret_debug_log_with_ollama "${temp_data}"
                rm -f "${temp_data}"
            else
                print_status "warning" "没有找到足够的请求/响应数据"
            fi
        fi
    fi
}

# ============================================
# 8. 主菜单
# ============================================

show_menu() {
    clear
    echo -e "${MAGENTA}${BOLD}"
    echo "╔════════════════════════════════════════════════╗"
    echo "║     GIBH-AGENT-V2 可视化监控运维脚本          ║"
    echo "╚════════════════════════════════════════════════╝"
    echo -e "${NC}"
    echo ""
    echo -e "${CYAN}${BOLD}🐳 Docker 服务管理：${NC}"
    echo "  1) 查看服务状态"
    echo "  2) 启动服务"
    echo "  3) 停止服务"
    echo "  4) 重启服务"
    echo "  5) 重新构建镜像"
    echo ""
    echo -e "${CYAN}${BOLD}🏥 健康检查：${NC}"
    echo "  6) 服务健康检查"
    echo ""
    echo -e "${CYAN}${BOLD}📋 日志监控：${NC}"
    echo "  7) 实时日志（API 服务器）"
    echo "  8) 实时日志（Worker）"
    echo "  9) 最近日志（API 服务器）"
    echo "  10) 所有服务日志摘要"
    echo "  11) 错误日志分析"
    echo "  12) 本地日志文件"
    echo ""
    echo -e "${CYAN}${BOLD}📊 数据监控：${NC}"
    echo "  13) 数据状态"
    echo "  14) 数据清理"
    echo ""
    echo -e "${CYAN}${BOLD}🔧 错误诊断：${NC}"
    echo "  15) 诊断 502 错误"
    echo "  16) 自动修复 502 错误"
    echo ""
    echo -e "${CYAN}${BOLD}📄 数据查看：${NC}"
    echo "  17) 查看 JSON 数据（原始格式）"
    echo ""
    echo -e "${CYAN}${BOLD}🔍 API 调试：${NC}"
    echo "  19) 实时监控 API 请求/响应"
    echo "  20) 查看最近的 API 请求/响应日志"
    echo ""
    echo -e "${CYAN}${BOLD}📊 综合监控：${NC}"
    echo "  18) 实时监控面板"
    echo ""
    echo -e "${CYAN}${BOLD}其他：${NC}"
    echo "  0) 退出"
    echo ""
    echo -e "${YELLOW}请选择操作 (0-20): ${NC}"
}

main() {
    cd "${PROJECT_DIR}" || {
        print_status "error" "无法进入项目目录: ${PROJECT_DIR}"
        exit 1
    }
    
    while true; do
        show_menu
        read -r choice
        
        case $choice in
            1) docker_status; read -p "按 Enter 继续..."; ;;
            2) docker_start; read -p "按 Enter 继续..."; ;;
            3) docker_stop; read -p "按 Enter 继续..."; ;;
            4) docker_restart; read -p "按 Enter 继续..."; ;;
            5) docker_rebuild; read -p "按 Enter 继续..."; ;;
            6) health_check; read -p "按 Enter 继续..."; ;;
            7) logs_realtime "api-server"; ;;
            8) logs_realtime "worker"; ;;
            9) logs_recent "api-server" 50; read -p "按 Enter 继续..."; ;;
            10) logs_all; read -p "按 Enter 继续..."; ;;
            11) logs_errors "api-server" 100; read -p "按 Enter 继续..."; ;;
            12) logs_file; read -p "按 Enter 继续..."; ;;
            13) data_status; read -p "按 Enter 继续..."; ;;
            14) data_cleanup; read -p "按 Enter 继续..."; ;;
            15) diagnose_502; read -p "按 Enter 继续..."; ;;
            16) fix_502; read -p "按 Enter 继续..."; ;;
            17)
                read -p "请输入 JSON 文件路径（或直接按 Enter 从 stdin 输入）: " json_file
                if [ -z "$json_file" ]; then
                    interpret_json
                else
                    interpret_json "$json_file"
                fi
                read -p "按 Enter 继续..."; ;;
            18) monitor_dashboard; ;;
            19) api_debug_realtime; ;;
            20) api_debug_recent; read -p "按 Enter 继续..."; ;;
            0)
                print_status "info" "再见！"
                exit 0
                ;;
            *)
                print_status "error" "无效选择，请重试"
                sleep 1
                ;;
        esac
    done
}

# 如果直接运行脚本，显示主菜单
if [ "${BASH_SOURCE[0]}" = "${0}" ]; then
    main "$@"
fi

