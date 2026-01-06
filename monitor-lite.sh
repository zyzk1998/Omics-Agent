#!/bin/bash
# ============================================
# GIBH-AGENT-V2 简化监控脚本（Lite Version）
# 专注于智能体逻辑调试
# ============================================

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
API_PORT=8028
DIRECT_PORT=8018

# Docker 命令前缀（用于处理权限）
DOCKER_CMD_PREFIX=""

# 检查 Docker 权限并设置命令前缀
check_docker_permission() {
    # 先尝试普通权限
    if docker info > /dev/null 2>&1; then
        DOCKER_CMD_PREFIX=""
        return 0
    fi
    
    # 检查是否在 docker 组中
    if groups | grep -q docker; then
        return 1
    fi
    
    # 需要 sudo 权限
    echo -e "${YELLOW}需要 sudo 权限访问 Docker${NC}"
    if sudo -v 2>/dev/null; then
        DOCKER_CMD_PREFIX="sudo "
        return 0
    else
        return 1
    fi
}

# 执行 Docker Compose 命令
docker_compose_cmd() {
    local cmd="$1"
    shift
    
    if [ -n "${DOCKER_CMD_PREFIX}" ]; then
        ${DOCKER_CMD_PREFIX}docker compose $cmd "$@" 2>/dev/null || \
        ${DOCKER_CMD_PREFIX}docker-compose $cmd "$@" 2>/dev/null
    else
        docker compose $cmd "$@" 2>/dev/null || \
        docker-compose $cmd "$@" 2>/dev/null
    fi
}

# 执行 Docker 命令
docker_cmd() {
    local cmd="$1"
    shift
    
    if [ -n "${DOCKER_CMD_PREFIX}" ]; then
        ${DOCKER_CMD_PREFIX}docker $cmd "$@" 2>&1
    else
        docker $cmd "$@" 2>&1
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
# 1. 服务管理（启动/重启/停止）
# ============================================

manage_services() {
    print_title "🚀 服务管理"
    
    echo "1) 启动服务"
    echo "2) 停止服务"
    echo "3) 重启服务"
    echo "4) 重新构建并启动"
    echo "5) 查看服务状态"
    echo "6) 返回主菜单"
    echo ""
    echo -e "${YELLOW}请选择 (1-6): ${NC}"
    read -r choice
    
    case $choice in
        1)
            if ! check_docker_permission; then
                print_status "error" "无法访问 Docker"
                return 1
            fi
            echo "🚀 启动服务..."
            mkdir -p ${PROJECT_DIR}/data/uploads ${PROJECT_DIR}/results ${PROJECT_DIR}/data/redis
            docker_compose_cmd up -d
            sleep 3
            docker_compose_cmd ps
            open_browser
            ;;
        2)
            if ! check_docker_permission; then
                print_status "error" "无法访问 Docker"
                return 1
            fi
            echo "🛑 停止服务..."
            docker_compose_cmd down
            print_status "ok" "服务已停止"
            ;;
        3)
            if ! check_docker_permission; then
                print_status "error" "无法访问 Docker"
                return 1
            fi
            echo "🔄 重启服务..."
            docker_compose_cmd restart
            sleep 3
            docker_compose_cmd ps
            open_browser
            ;;
        4)
            if ! check_docker_permission; then
                print_status "error" "无法访问 Docker"
                return 1
            fi
            echo "🔨 重新构建并启动..."
            docker_compose_cmd build --no-cache
            docker_compose_cmd up -d
            sleep 5
            docker_compose_cmd ps
            open_browser
            ;;
        5)
            if ! check_docker_permission; then
                print_status "error" "无法访问 Docker"
                return 1
            fi
            docker_compose_cmd ps
            ;;
        6)
            return 0
            ;;
        *)
            print_status "error" "无效选择"
            ;;
    esac
    
    read -p "按 Enter 继续..."
}

# ============================================
# 2. Agent Logic Trace（核心功能）
# ============================================

agent_trace() {
    clear
    echo -e "${MAGENTA}${BOLD}"
    echo "╔════════════════════════════════════════════════╗"
    echo "║      🕵️  Agent Logic Trace Mode (God Mode)     ║"
    echo "║      全栈日志监控 - 高亮模式，不隐藏任何日志   ║"
    echo "╚════════════════════════════════════════════════╝"
    echo -e "${NC}"
    echo ""
    echo -e "${CYAN}${BOLD}高亮规则：${NC}"
    echo -e "  ${RED}${BOLD}ERROR / Exception / Traceback${NC} - 错误信息（红色粗体）"
    echo -e "  ${GREEN}User Query / Process Query${NC} - 用户查询"
    echo -e "  ${CYAN}Router${NC} - 路由决策"
    echo -e "  ${YELLOW}Thought / <think>${NC} - LLM 思考过程"
    echo -e "  ${MAGENTA}Action / Tool Call${NC} - 工具调用"
    echo -e "  ${BLUE}Observation / Tool Output${NC} - 工具输出"
    echo -e "  ${WHITE}其他所有日志${NC} - 显示为灰色/白色（不隐藏）"
    echo ""
    echo -e "${YELLOW}正在监听全栈日志（api-server + worker）...${NC}"
    echo -e "${YELLOW}按 Ctrl+C 退出${NC}"
    echo ""
    print_separator
    echo ""
    
    # 检查 Docker 权限
    if ! check_docker_permission; then
        print_status "error" "无法访问 Docker，尝试使用本地日志..."
        # 降级到本地日志，使用 Python 脚本进行高亮处理
        tail -f ${PROJECT_DIR}/gibh_agent.log 2>/dev/null | python3 -u -c "
import sys
import re

# 只过滤真正无用的噪音（严格限制）
hard_noise_patterns = [
    r'^GET /health',
    r'^GET /static',
    r'^200 OK$',
    r'^$'  # 空行
]

# 定义关键词和颜色（按优先级排序，错误最高优先级）
keywords = [
    # 错误相关（最高优先级，红色粗体）
    ('Traceback', '\033[1;31m'),  # Bold Red
    ('SyntaxError', '\033[1;31m'),
    ('ImportError', '\033[1;31m'),
    ('IndentationError', '\033[1;31m'),
    ('ModuleNotFoundError', '\033[1;31m'),
    ('FileNotFoundError', '\033[1;31m'),
    ('413', '\033[1;31m'),  # Nginx 413 Request Entity Too Large
    ('Request Entity Too Large', '\033[1;31m'),
    ('ERROR', '\033[1;31m'),
    ('Exception', '\033[1;31m'),
    ('❌', '\033[1;31m'),
    ('Failed', '\033[1;31m'),
    ('failed', '\033[1;31m'),
    ('Error', '\033[1;31m'),
    # AI 相关（绿色）
    ('收到聊天请求', '\033[0;32m'),
    ('处理查询', '\033[0;32m'),
    ('chat_endpoint entry', '\033[0;32m'),
    ('Before process_query', '\033[0;32m'),
    ('After process_query', '\033[0;32m'),
    ('workflow_config', '\033[0;32m'),
    ('Final Answer', '\033[0;32m'),
    ('✅', '\033[0;32m'),
    # 路由相关（青色）
    ('路由', '\033[0;36m'),
    ('Router', '\033[0;36m'),
    ('routing', '\033[0;36m'),
    ('🎯', '\033[0;36m'),
    # 思考过程（黄色）
    ('Thought', '\033[1;33m'),
    ('<think>', '\033[1;33m'),
    ('reasoning', '\033[1;33m'),
    # 工具调用（洋红色）
    ('执行步骤', '\033[0;35m'),
    ('execute_workflow', '\033[0;35m'),
    ('Tool Call', '\033[0;35m'),
    ('Action', '\033[0;35m'),
    ('🔧', '\033[0;35m'),
    # 工具输出（蓝色）
    ('Tool Output', '\033[0;34m'),
    ('Observation', '\033[0;34m'),
    ('📊', '\033[0;34m'),
    ('💬', '\033[0;32m')
]
NC = '\033[0m'
DEFAULT_COLOR = '\033[0;37m'  # 灰色（默认显示所有日志）

try:
    for line in sys.stdin:
        if not line:
            continue
        line = line.rstrip()
        
        # 硬过滤：只隐藏真正无用的噪音
        is_hard_noise = any(re.search(pattern, line, re.IGNORECASE) for pattern in hard_noise_patterns)
        if is_hard_noise:
            continue
        
        # 检查是否包含关键词（按优先级）
        matched = False
        for keyword, color in keywords:
            if keyword in line:
                # 高亮关键词
                highlighted = re.sub(
                    f'({re.escape(keyword)})',
                    f'{color}\\1{NC}',
                    line,
                    flags=re.IGNORECASE
                )
                print(highlighted, flush=True)
                matched = True
                break
        
        # 如果没有匹配关键词，显示原始日志（灰色）- 不隐藏！
        if not matched:
            print(f'{DEFAULT_COLOR}{line}{NC}', flush=True)
except KeyboardInterrupt:
    sys.exit(0)
except Exception as e:
    # 如果 Python 脚本出错，降级到原始日志
    print(f'\033[1;31m过滤脚本错误: {e}\033[0m', file=sys.stderr)
    sys.exit(1)
" || tail -f ${PROJECT_DIR}/gibh_agent.log 2>/dev/null
        return
    fi
    
    # 实时监控所有容器日志（api-server + worker），使用 Python 脚本进行高亮处理
    # 注意：如果 nginx 服务不存在，docker compose 会忽略它
    docker_compose_cmd logs -f api-server worker 2>/dev/null 2>&1 | python3 -u -c "
import sys
import re

# 只过滤真正无用的噪音（严格限制）
hard_noise_patterns = [
    r'^GET /health',
    r'^GET /static',
    r'^200 OK$',
    r'^$'  # 空行
]

# 定义关键词和颜色（按优先级排序，错误最高优先级）
keywords = [
    # 错误相关（最高优先级，红色粗体）
    ('Traceback', '\033[1;31m'),  # Bold Red
    ('SyntaxError', '\033[1;31m'),
    ('ImportError', '\033[1;31m'),
    ('IndentationError', '\033[1;31m'),
    ('ModuleNotFoundError', '\033[1;31m'),
    ('FileNotFoundError', '\033[1;31m'),
    ('413', '\033[1;31m'),  # Nginx 413 Request Entity Too Large
    ('Request Entity Too Large', '\033[1;31m'),
    ('ERROR', '\033[1;31m'),
    ('Exception', '\033[1;31m'),
    ('❌', '\033[1;31m'),
    ('Failed', '\033[1;31m'),
    ('failed', '\033[1;31m'),
    ('Error', '\033[1;31m'),
    # AI 相关（绿色）
    ('收到聊天请求', '\033[0;32m'),
    ('处理查询', '\033[0;32m'),
    ('chat_endpoint entry', '\033[0;32m'),
    ('Before process_query', '\033[0;32m'),
    ('After process_query', '\033[0;32m'),
    ('workflow_config', '\033[0;32m'),
    ('Final Answer', '\033[0;32m'),
    ('✅', '\033[0;32m'),
    # 路由相关（青色）
    ('路由', '\033[0;36m'),
    ('Router', '\033[0;36m'),
    ('routing', '\033[0;36m'),
    ('🎯', '\033[0;36m'),
    # 思考过程（黄色）
    ('Thought', '\033[1;33m'),
    ('<think>', '\033[1;33m'),
    ('reasoning', '\033[1;33m'),
    # 工具调用（洋红色）
    ('执行步骤', '\033[0;35m'),
    ('execute_workflow', '\033[0;35m'),
    ('Tool Call', '\033[0;35m'),
    ('Action', '\033[0;35m'),
    ('🔧', '\033[0;35m'),
    # 工具输出（蓝色）
    ('Tool Output', '\033[0;34m'),
    ('Observation', '\033[0;34m'),
    ('📊', '\033[0;34m'),
    ('💬', '\033[0;32m')
]
NC = '\033[0m'
DEFAULT_COLOR = '\033[0;37m'  # 灰色（默认显示所有日志）

try:
    for line in sys.stdin:
        if not line:
            continue
        line = line.rstrip()
        
        # 硬过滤：只隐藏真正无用的噪音
        is_hard_noise = any(re.search(pattern, line, re.IGNORECASE) for pattern in hard_noise_patterns)
        if is_hard_noise:
            continue
        
        # 检查是否包含关键词（按优先级）
        matched = False
        for keyword, color in keywords:
            if keyword in line:
                # 高亮关键词
                highlighted = re.sub(
                    f'({re.escape(keyword)})',
                    f'{color}\\1{NC}',
                    line,
                    flags=re.IGNORECASE
                )
                print(highlighted, flush=True)
                matched = True
                break
        
        # 如果没有匹配关键词，显示原始日志（灰色）- 不隐藏！
        if not matched:
            print(f'{DEFAULT_COLOR}{line}{NC}', flush=True)
except KeyboardInterrupt:
    sys.exit(0)
except Exception as e:
    # 如果 Python 脚本出错，降级到原始日志
    print(f'\033[1;31m过滤脚本错误: {e}\033[0m', file=sys.stderr)
    sys.exit(1)
" || {
    # 如果 Python 脚本失败，使用原始日志（不做任何过滤）
    print_status "warning" "Python 过滤脚本失败，显示原始日志..."
    docker_compose_cmd logs -f api-server worker 2>/dev/null
}
}

# ============================================
# 3. 系统日志（原始）
# ============================================

system_logs() {
    print_title "📋 系统日志（原始）"
    
    echo "1) 实时日志（API 服务器）"
    echo "2) 实时日志（Worker）"
    echo "3) 最近日志（API 服务器，50行）"
    echo "4) 错误日志（API 服务器）"
    echo "5) 本地日志文件（gibh_agent.log）"
    echo "6) Debug 日志（.cursor/debug.log）"
    echo "7) 返回主菜单"
    echo ""
    echo -e "${YELLOW}请选择 (1-7): ${NC}"
    read -r choice
    
    if ! check_docker_permission 2>/dev/null; then
        print_status "warning" "Docker 未运行或无权限，只能查看本地日志"
    fi
    
    case $choice in
        1)
            if check_docker_permission 2>/dev/null; then
                echo -e "${YELLOW}按 Ctrl+C 退出${NC}\n"
                docker_compose_cmd logs -f api-server 2>/dev/null
            else
                echo -e "${YELLOW}按 Ctrl+C 退出${NC}\n"
                tail -f ${PROJECT_DIR}/gibh_agent.log 2>/dev/null
            fi
            ;;
        2)
            if check_docker_permission 2>/dev/null; then
                echo -e "${YELLOW}按 Ctrl+C 退出${NC}\n"
                docker_compose_cmd logs -f worker 2>/dev/null
            else
                print_status "warning" "Worker 日志需要 Docker"
            fi
            ;;
        3)
            if check_docker_permission 2>/dev/null; then
                echo -e "${CYAN}${BOLD}原始日志（最近 100 行，无过滤）${NC}\n"
                docker_compose_cmd logs --tail 100 api-server worker 2>/dev/null
            else
                echo -e "${CYAN}${BOLD}本地日志（最近 50 行，无过滤）${NC}\n"
                tail -50 ${PROJECT_DIR}/gibh_agent.log 2>/dev/null
            fi
            read -p "按 Enter 继续..."
            ;;
        4)
            if check_docker_permission 2>/dev/null; then
                docker_compose_cmd logs --tail 100 api-server 2>/dev/null | \
                    grep -i -E "error|exception|failed|traceback|❌" | tail -20
            else
                grep -i -E "error|exception|failed|traceback|❌" ${PROJECT_DIR}/gibh_agent.log 2>/dev/null | tail -20
            fi
            read -p "按 Enter 继续..."
            ;;
        5)
            if [ -f "${PROJECT_DIR}/gibh_agent.log" ]; then
                echo -e "${WHITE}文件: ${PROJECT_DIR}/gibh_agent.log${NC}"
                echo -e "${WHITE}大小: $(du -h ${PROJECT_DIR}/gibh_agent.log | cut -f1)${NC}"
                echo -e "${WHITE}最后 20 行：${NC}\n"
                tail -20 ${PROJECT_DIR}/gibh_agent.log
            else
                print_status "warning" "日志文件不存在"
            fi
            read -p "按 Enter 继续..."
            ;;
        6)
            if [ -f "${PROJECT_DIR}/.cursor/debug.log" ]; then
                echo -e "${WHITE}文件: ${PROJECT_DIR}/.cursor/debug.log${NC}"
                echo -e "${WHITE}大小: $(du -h ${PROJECT_DIR}/.cursor/debug.log | cut -f1)${NC}"
                echo -e "${WHITE}最后 10 条 JSON 记录：${NC}\n"
                tail -10 ${PROJECT_DIR}/.cursor/debug.log | while read line; do
                    if command -v python3 &> /dev/null; then
                        echo "$line" | python3 -m json.tool 2>/dev/null || echo "$line"
                    else
                        echo "$line"
                    fi
                    echo ""
                done
            else
                print_status "warning" "Debug 日志文件不存在"
            fi
            read -p "按 Enter 继续..."
            ;;
        7)
            return 0
            ;;
        *)
            print_status "error" "无效选择"
            ;;
    esac
}

# ============================================
# 4. 高级工具（子菜单）
# ============================================

advanced_tools() {
    print_title "🛠️ 高级工具"
    
    echo "1) 数据清理"
    echo "2) 修复 502 错误"
    echo "3) 健康检查"
    echo "4) 查看数据状态"
    echo "5) 查看 JSON 数据"
    echo "6) 返回主菜单"
    echo ""
    echo -e "${YELLOW}请选择 (1-6): ${NC}"
    read -r choice
    
    case $choice in
        1)
            echo -e "${YELLOW}警告：此操作将删除数据文件${NC}\n"
            echo "1. 清理上传文件"
            echo "2. 清理结果文件"
            echo "3. 清理所有数据"
            echo "4. 取消"
            read -p "请选择 (1-4): " sub_choice
            
            case $sub_choice in
                1)
                    read -p "确认删除上传文件？(y/N): " -n 1 -r
                    echo
                    if [[ $REPLY =~ ^[Yy]$ ]]; then
                        find ${PROJECT_DIR}/data/uploads -type f ! -name "*.meta.json" -delete 2>/dev/null
                        print_status "ok" "上传文件已清理"
                    fi
                    ;;
                2)
                    read -p "确认删除结果文件？(y/N): " -n 1 -r
                    echo
                    if [[ $REPLY =~ ^[Yy]$ ]]; then
                        rm -rf ${PROJECT_DIR}/results/* 2>/dev/null
                        print_status "ok" "结果文件已清理"
                    fi
                    ;;
                3)
                    read -p "确认删除所有数据？(y/N): " -n 1 -r
                    echo
                    if [[ $REPLY =~ ^[Yy]$ ]]; then
                        find ${PROJECT_DIR}/data/uploads -type f ! -name "*.meta.json" -delete 2>/dev/null
                        rm -rf ${PROJECT_DIR}/results/* 2>/dev/null
                        print_status "ok" "所有数据已清理"
                    fi
                    ;;
            esac
            read -p "按 Enter 继续..."
            ;;
        2)
            print_title "🔧 自动修复 502 错误"
            if ! check_docker_permission; then
                print_status "error" "无法访问 Docker"
                read -p "按 Enter 继续..."
                return 1
            fi
            
            echo "1️⃣ 停止所有容器..."
            docker_compose_cmd down
            
            echo ""
            echo "2️⃣ 检查 requirements.txt..."
            if ! grep -q "paramiko" ${PROJECT_DIR}/requirements.txt; then
                echo "paramiko>=3.0.0" >> ${PROJECT_DIR}/requirements.txt
                print_status "ok" "已添加 paramiko"
            fi
            
            echo ""
            echo "3️⃣ 重新构建镜像..."
            docker_compose_cmd build --no-cache
            
            echo ""
            echo "4️⃣ 启动服务..."
            docker_compose_cmd up -d
            
            echo ""
            echo "5️⃣ 等待服务启动（10秒）..."
            sleep 10
            
            echo ""
            echo "6️⃣ 检查服务状态..."
            docker_compose_cmd ps
            
            read -p "按 Enter 继续..."
            ;;
        3)
            print_title "🏥 健康检查"
            if ! check_docker_permission 2>/dev/null; then
                print_status "warning" "Docker 未运行或无权限"
            else
                echo -e "${WHITE}容器状态：${NC}"
                docker_compose_cmd ps
            fi
            
            echo ""
            echo -e "${WHITE}API 服务器响应：${NC}"
            api_found=false
            if curl -s -o /dev/null -w "%{http_code}" http://localhost:${API_PORT}/ 2>/dev/null | grep -q "200\|301\|302"; then
                print_status "ok" "Docker API 服务器正常 (http://localhost:${API_PORT})"
                api_found=true
            fi
            if curl -s -o /dev/null -w "%{http_code}" http://localhost:${DIRECT_PORT}/ 2>/dev/null | grep -q "200\|301\|302"; then
                print_status "ok" "直接运行服务器正常 (http://localhost:${DIRECT_PORT})"
                api_found=true
            fi
            if [ "$api_found" = false ]; then
                print_status "error" "API 服务器无响应"
            fi
            
            read -p "按 Enter 继续..."
            ;;
        4)
            print_title "📊 数据状态"
            echo -e "${WHITE}上传文件：${NC}"
            if [ -d "${PROJECT_DIR}/data/uploads" ]; then
                upload_count=$(find ${PROJECT_DIR}/data/uploads -type f ! -name "*.meta.json" 2>/dev/null | wc -l)
                upload_size=$(du -sh ${PROJECT_DIR}/data/uploads 2>/dev/null | cut -f1 || echo "0")
                echo -e "  文件数: ${upload_count} | 大小: ${upload_size}"
            fi
            
            echo ""
            echo -e "${WHITE}结果文件：${NC}"
            if [ -d "${PROJECT_DIR}/results" ]; then
                result_count=$(find ${PROJECT_DIR}/results -type f 2>/dev/null | wc -l)
                result_size=$(du -sh ${PROJECT_DIR}/results 2>/dev/null | cut -f1 || echo "0")
                echo -e "  文件数: ${result_count} | 大小: ${result_size}"
            fi
            
            read -p "按 Enter 继续..."
            ;;
        5)
            print_title "📄 查看 JSON 数据"
            read -p "请输入 JSON 文件路径: " json_file
            if [ -f "$json_file" ]; then
                if command -v python3 &> /dev/null; then
                    echo ""
                    cat "$json_file" | python3 -m json.tool 2>/dev/null || cat "$json_file"
                else
                    cat "$json_file"
                fi
            else
                print_status "error" "文件不存在: $json_file"
            fi
            read -p "按 Enter 继续..."
            ;;
        6)
            return 0
            ;;
        *)
            print_status "error" "无效选择"
            ;;
    esac
}

# ============================================
# 主菜单
# ============================================

show_menu() {
    clear
    echo -e "${MAGENTA}${BOLD}"
    echo "╔════════════════════════════════════════════════╗"
    echo "║     GIBH-AGENT-V2 简化监控脚本 (Lite)         ║"
    echo "║     专注于智能体逻辑调试                       ║"
    echo "╚════════════════════════════════════════════════╝"
    echo -e "${NC}"
    echo ""
    echo -e "${CYAN}${BOLD}核心功能：${NC}"
    echo "  1) 🚀 服务管理（启动/重启/停止）"
    echo "  2) 🕵️  Agent Logic Trace（实时监听智能体逻辑流程）"
    echo "  3) 📋 系统日志（原始日志，用于深度调试）"
    echo "  4) 🛠️  高级工具（数据清理、502修复、健康检查等）"
    echo ""
    echo -e "${CYAN}${BOLD}其他：${NC}"
    echo "  0) ❌ 退出"
    echo ""
    echo -e "${YELLOW}请选择操作 (0-4): ${NC}"
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
            1) manage_services; ;;
            2) agent_trace; ;;
            3) system_logs; ;;
            4) advanced_tools; ;;
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

