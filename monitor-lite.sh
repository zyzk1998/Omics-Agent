#!/bin/bash
# ============================================
# GIBH-AGENT-V2 极简监控脚本（Minimal Version）
# ============================================

# 颜色定义
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
MAGENTA='\033[0;35m'
CYAN='\033[0;36m'
WHITE='\033[1;37m'
NC='\033[0m'
BOLD='\033[1m'

# 项目配置
PROJECT_DIR="/home/ubuntu/GIBH-AGENT-V2"
API_PORT=8028

# Docker 命令前缀
DOCKER_CMD_PREFIX=""

check_docker_permission() {
    if docker info > /dev/null 2>&1; then
        DOCKER_CMD_PREFIX=""
        return 0
    fi
    if groups | grep -q docker; then
        return 1
    fi
    if sudo -v 2>/dev/null; then
        DOCKER_CMD_PREFIX="sudo "
        return 0
    fi
    return 1
}

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

print_status() {
    if [ "$1" = "ok" ]; then
        echo -e "${GREEN}✅ $2${NC}"
    elif [ "$1" = "error" ]; then
        echo -e "${RED}❌ $2${NC}"
    elif [ "$1" = "warning" ]; then
        echo -e "${YELLOW}⚠️  $2${NC}"
    fi
}

# ============================================
# 1. 服务管理
# ============================================

manage_services() {
    clear
    echo -e "${CYAN}${BOLD}🚀 服务管理${NC}\n"
    echo "1) 启动服务"
    echo "2) 停止服务"
    echo "3) 重启服务"
    echo "0) 返回"
    echo ""
    read -p "请选择: " choice
    
    if ! check_docker_permission; then
        print_status "error" "无法访问 Docker"
        read -p "按 Enter 继续..."
        return
    fi
    
    case $choice in
        1)
            echo "🚀 启动服务..."
            mkdir -p ${PROJECT_DIR}/data/uploads ${PROJECT_DIR}/results ${PROJECT_DIR}/data/redis
            docker_compose_cmd up -d
            sleep 3
            docker_compose_cmd ps
            ;;
        2)
            echo "🛑 停止服务..."
            docker_compose_cmd down
            print_status "ok" "服务已停止"
            ;;
        3)
            echo "🔄 重启服务..."
            docker_compose_cmd restart
            sleep 3
            docker_compose_cmd ps
            ;;
    esac
    read -p "按 Enter 继续..."
}

# ============================================
# 2. Agent Logic Trace (God Mode)
# ============================================

agent_trace() {
    clear
    echo -e "${MAGENTA}${BOLD}"
    echo "╔════════════════════════════════════════════════╗"
    echo "║      🕵️  Agent Logic Trace Mode (God Mode)     ║"
    echo "╚════════════════════════════════════════════════╝"
    echo -e "${NC}\n"
    echo -e "${CYAN}${BOLD}高亮规则：${NC}"
    echo -e "  ${RED}${BOLD}ERROR / Exception${NC} - 错误（红色）"
    echo -e "  ${GREEN}User Query / Process Query${NC} - 用户查询（绿色）"
    echo -e "  ${CYAN}Router${NC} - 路由决策（青色）"
    echo -e "  ${YELLOW}Thought${NC} - LLM 思考（黄色）"
    echo -e "  ${MAGENTA}Tool Call${NC} - 工具调用（洋红色）"
    echo -e "  ${BLUE}Tool Output${NC} - 工具输出（蓝色）"
    echo -e "  ${CYAN}🔥 [LLM_RAW_DUMP]${NC} - LLM 原始 JSON（青色，美化打印）"
    echo ""
    echo -e "${YELLOW}按 Ctrl+C 退出${NC}\n"
    
    if ! check_docker_permission; then
        print_status "warning" "无法访问 Docker，使用本地日志..."
        tail -f ${PROJECT_DIR}/gibh_agent.log 2>/dev/null | python3 -u -c "
import sys, re, json

hard_noise = [r'^GET /health', r'^GET /static', r'^200 OK$', r'^$']
keywords = [
    ('[LLM_RAW_DUMP]', '\033[1;95m'),
    ('Traceback', '\033[1;31m'), ('ERROR', '\033[1;31m'), ('Exception', '\033[1;31m'),
    ('收到聊天请求', '\033[0;32m'), ('处理查询', '\033[0;32m'), ('✅', '\033[0;32m'),
    ('路由', '\033[0;36m'), ('Router', '\033[0;36m'), ('🎯', '\033[0;36m'),
    ('Thought', '\033[1;33m'), ('reasoning', '\033[1;33m'),
    ('Tool Call', '\033[0;35m'), ('Action', '\033[0;35m'),
    ('Tool Output', '\033[0;34m'), ('Observation', '\033[0;34m')
]
NC = '\033[0m'

for line in sys.stdin:
    if not line.strip():
        continue
    line = line.rstrip()
    
    # 🔥 Task 4: 特殊处理 LLM_RAW_DUMP - 美化 JSON (Cyan 颜色)
    if '[LLM_RAW_DUMP]' in line:
        parts = line.split('[LLM_RAW_DUMP]', 1)
        if len(parts) > 1:
            json_str = parts[1].strip()
            try:
                # 尝试解析并美化 JSON
                parsed = json.loads(json_str)
                pretty = json.dumps(parsed, indent=2, ensure_ascii=False)
                print(f'\033[0;36m🔥 [LLM_RAW_DUMP]\033[0m\n\033[0;36m{pretty}\033[0m', flush=True)
            except:
                # 如果不是 JSON，直接显示
                print(f'\033[0;36m🔥 [LLM_RAW_DUMP]\033[0m {json_str}', flush=True)
            continue
    
    if any(re.search(p, line, re.I) for p in hard_noise):
        continue
    
    matched = False
    for kw, color in keywords:
        if kw in line:
            print(f'{color}{line}{NC}', flush=True)
            matched = True
            break
    if not matched:
        print(f'\033[0;37m{line}{NC}', flush=True)
" || tail -f ${PROJECT_DIR}/gibh_agent.log 2>/dev/null
        return
    fi
    
    docker_compose_cmd logs -f api-server worker 2>/dev/null 2>&1 | python3 -u -c "
import sys, re, json

hard_noise = [r'^GET /health', r'^GET /static', r'^200 OK$', r'^$']
keywords = [
    ('[LLM_RAW_DUMP]', '\033[1;95m'),
    ('Traceback', '\033[1;31m'), ('ERROR', '\033[1;31m'), ('Exception', '\033[1;31m'),
    ('收到聊天请求', '\033[0;32m'), ('处理查询', '\033[0;32m'), ('✅', '\033[0;32m'),
    ('路由', '\033[0;36m'), ('Router', '\033[0;36m'), ('🎯', '\033[0;36m'),
    ('Thought', '\033[1;33m'), ('reasoning', '\033[1;33m'),
    ('Tool Call', '\033[0;35m'), ('Action', '\033[0;35m'),
    ('Tool Output', '\033[0;34m'), ('Observation', '\033[0;34m')
]
NC = '\033[0m'

for line in sys.stdin:
    if not line.strip():
        continue
    line = line.rstrip()
    
    # 🔥 Task 4: 特殊处理 LLM_RAW_DUMP - 美化 JSON (Cyan 颜色)
    if '[LLM_RAW_DUMP]' in line:
        parts = line.split('[LLM_RAW_DUMP]', 1)
        if len(parts) > 1:
            json_str = parts[1].strip()
            try:
                # 尝试解析并美化 JSON
                parsed = json.loads(json_str)
                pretty = json.dumps(parsed, indent=2, ensure_ascii=False)
                print(f'\033[0;36m🔥 [LLM_RAW_DUMP]\033[0m\n\033[0;36m{pretty}\033[0m', flush=True)
            except:
                # 如果不是 JSON，直接显示
                print(f'\033[0;36m🔥 [LLM_RAW_DUMP]\033[0m {json_str}', flush=True)
            continue
    
    if any(re.search(p, line, re.I) for p in hard_noise):
        continue
    
    matched = False
    for kw, color in keywords:
        if kw in line:
            print(f'{color}{line}{NC}', flush=True)
            matched = True
            break
    if not matched:
        print(f'\033[0;37m{line}{NC}', flush=True)
" || docker_compose_cmd logs -f api-server worker 2>/dev/null
}

# ============================================
# 3. 系统日志（原始）
# ============================================

system_logs() {
    clear
    echo -e "${CYAN}${BOLD}📋 系统日志（原始）${NC}\n"
    echo "1) 实时日志（API + Worker）"
    echo "2) 最近 100 行"
    echo "3) 错误日志"
    echo "0) 返回"
    echo ""
    read -p "请选择: " choice
    
    if ! check_docker_permission 2>/dev/null; then
        print_status "warning" "Docker 未运行，只能查看本地日志"
    fi
    
    case $choice in
        1)
            if check_docker_permission 2>/dev/null; then
                echo -e "${YELLOW}按 Ctrl+C 退出${NC}\n"
                docker_compose_cmd logs -f api-server worker 2>/dev/null
            else
                echo -e "${YELLOW}按 Ctrl+C 退出${NC}\n"
                tail -f ${PROJECT_DIR}/gibh_agent.log 2>/dev/null
            fi
            ;;
        2)
            if check_docker_permission 2>/dev/null; then
                docker_compose_cmd logs --tail 100 api-server worker 2>/dev/null
            else
                tail -100 ${PROJECT_DIR}/gibh_agent.log 2>/dev/null
            fi
            read -p "按 Enter 继续..."
            ;;
        3)
            if check_docker_permission 2>/dev/null; then
                docker_compose_cmd logs --tail 200 api-server worker 2>/dev/null | \
                    grep -i -E "error|exception|failed|traceback|❌" | tail -30
            else
                grep -i -E "error|exception|failed|traceback|❌" ${PROJECT_DIR}/gibh_agent.log 2>/dev/null | tail -30
            fi
            read -p "按 Enter 继续..."
            ;;
    esac
}

# ============================================
# 4. 清理数据
# ============================================

clean_data() {
    clear
    echo -e "${CYAN}${BOLD}🧹 清理数据${NC}\n"
    echo "1) 清理上传文件"
    echo "2) 清理结果文件"
    echo "3) 清理所有数据"
    echo "0) 返回"
    echo ""
    read -p "请选择: " choice
    
    case $choice in
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
}

# ============================================
# 5. 实时监听 LLM JSON
# ============================================

llm_json_monitor() {
    clear
    echo -e "${CYAN}${BOLD}🧠 实时监听 LLM JSON${NC}\n"
    echo -e "${YELLOW}按 Ctrl+C 退出${NC}\n"
    
    # 🔥 Task 2: 检查 jq 是否安装
    if command -v jq >/dev/null 2>&1; then
        USE_JQ=true
    else
        USE_JQ=false
        echo -e "${YELLOW}⚠️  jq 未安装，将使用 Python 格式化 JSON${NC}\n"
    fi
    
    if ! check_docker_permission; then
        print_status "warning" "无法访问 Docker，使用本地日志..."
        if [ "$USE_JQ" = true ]; then
            tail -f ${PROJECT_DIR}/gibh_agent.log 2>/dev/null | grep --line-buffered "LLM_RAW_DUMP" | sed 's/.*\[LLM_RAW_DUMP\] //' | jq . 2>/dev/null || \
            tail -f ${PROJECT_DIR}/gibh_agent.log 2>/dev/null | grep --line-buffered "LLM_RAW_DUMP" | sed 's/.*\[LLM_RAW_DUMP\] //' | sed 's/^/\033[0;36m🔥 [LLM_RAW_DUMP]\033[0m /'
        else
            tail -f ${PROJECT_DIR}/gibh_agent.log 2>/dev/null | grep --line-buffered "LLM_RAW_DUMP" | python3 -u -c "
import sys, json

for line in sys.stdin:
    line = line.rstrip()
    if '[LLM_RAW_DUMP]' in line:
        # 提取 JSON 部分（可能在行中的任何位置）
        idx = line.find('[LLM_RAW_DUMP]')
        if idx >= 0:
            json_str = line[idx + len('[LLM_RAW_DUMP]'):].strip()
            if json_str:
                try:
                    parsed = json.loads(json_str)
                    pretty = json.dumps(parsed, indent=2, ensure_ascii=False)
                    print(f'\033[0;36m🔥 [LLM_RAW_DUMP]\033[0m\n\033[0;36m{pretty}\033[0m', flush=True)
                except json.JSONDecodeError:
                    # 如果不是完整 JSON，尝试提取并显示
                    print(f'\033[0;36m🔥 [LLM_RAW_DUMP]\033[0m {json_str[:500]}...', flush=True)
                except:
                    print(f'\033[0;36m🔥 [LLM_RAW_DUMP]\033[0m {json_str}', flush=True)
" || echo "无法读取日志文件"
        fi
    else
        # 🔥 Task 2: 使用更简单、更健壮的命令
        if [ "$USE_JQ" = true ]; then
            docker_compose_cmd logs -f api-server worker 2>/dev/null | grep --line-buffered "LLM_RAW_DUMP" | sed 's/.*\[LLM_RAW_DUMP\] //' | jq . 2>/dev/null || \
            docker_compose_cmd logs -f api-server worker 2>/dev/null | grep --line-buffered "LLM_RAW_DUMP" | sed 's/.*\[LLM_RAW_DUMP\] //' | sed 's/^/\033[0;36m🔥 [LLM_RAW_DUMP]\033[0m /'
        else
            docker_compose_cmd logs -f api-server worker 2>/dev/null | grep --line-buffered "LLM_RAW_DUMP" | python3 -u -c "
import sys, json

for line in sys.stdin:
    line = line.rstrip()
    if '[LLM_RAW_DUMP]' in line:
        # 提取 JSON 部分（可能在行中的任何位置，Docker logs 可能包含时间戳和容器名）
        idx = line.find('[LLM_RAW_DUMP]')
        if idx >= 0:
            json_str = line[idx + len('[LLM_RAW_DUMP]'):].strip()
            # 移除可能的前导空格
            json_str = json_str.lstrip()
            if json_str:
                try:
                    parsed = json.loads(json_str)
                    pretty = json.dumps(parsed, indent=2, ensure_ascii=False)
                    print(f'\033[0;36m🔥 [LLM_RAW_DUMP]\033[0m\n\033[0;36m{pretty}\033[0m', flush=True)
                except json.JSONDecodeError as e:
                    # 如果不是完整 JSON，显示原始内容（可能被截断）
                    print(f'\033[0;36m🔥 [LLM_RAW_DUMP]\033[0m {json_str[:500]}...', flush=True)
                    print(f'\033[0;33m⚠️  JSON 解析失败: {e}\033[0m', flush=True)
                except Exception as e:
                    print(f'\033[0;36m🔥 [LLM_RAW_DUMP]\033[0m {json_str}', flush=True)
                    print(f'\033[0;33m⚠️  处理失败: {e}\033[0m', flush=True)
" || echo "无法读取 Docker 日志"
        fi
    fi
}

# ============================================
# 主菜单
# ============================================

show_menu() {
    clear
    echo -e "${MAGENTA}${BOLD}"
    echo "╔════════════════════════════════════════════════╗"
    echo "║     GIBH-AGENT-V2 极简监控脚本                ║"
    echo "╚════════════════════════════════════════════════╝"
    echo -e "${NC}\n"
    echo "  1) 🚀 服务管理"
    echo "  2) 🕵️  Agent Logic Trace"
    echo "  3) 📋 系统日志"
    echo "  4) 🧹 清理数据"
    echo "  5) 🧠 实时监听 LLM JSON"
    echo ""
    echo "  0) ❌ 退出"
    echo ""
    echo -e "${YELLOW}请选择 (0-5): ${NC}"
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
            4) clean_data; ;;
            5) llm_json_monitor; ;;
            0)
                print_status "ok" "再见！"
                exit 0
                ;;
            *)
                print_status "error" "无效选择"
                sleep 1
                ;;
        esac
    done
}

if [ "${BASH_SOURCE[0]}" = "${0}" ]; then
    main "$@"
fi
