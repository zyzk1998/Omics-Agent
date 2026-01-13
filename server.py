"""
GIBH-AGENT-V2 测试服务器
提供简单的 Web 接口用于测试功能
"""
import os
import sys
import json
import logging
import asyncio
import re
import secrets
from pathlib import Path
from typing import List, Optional, Set
from datetime import datetime
from collections import deque

from fastapi import FastAPI, UploadFile, File, HTTPException, Request
from typing import Optional
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import StreamingResponse, JSONResponse, HTMLResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel

# 添加项目路径
sys.path.insert(0, str(Path(__file__).parent))

from gibh_agent import create_agent
from gibh_agent.core.file_inspector import FileInspector

# 配置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler('gibh_agent.log', encoding='utf-8')
    ]
)
logger = logging.getLogger(__name__)

# 创建 FastAPI 应用
app = FastAPI(title="GIBH-AGENT-V2 Test Server")

# 🔥 Step 2: Tool-RAG 架构 - Vector Database Integration
# 初始化工具检索器（在启动时同步工具）
tool_retriever = None
workflow_planner = None
try:
    from gibh_agent.core.tool_retriever import ToolRetriever
    # 🔥 Step 4: 模块化工具系统 - 自动发现和加载所有工具
    from gibh_agent.tools import load_all_tools
    
    # 初始化工具检索器
    chroma_dir = os.getenv("CHROMA_PERSIST_DIR", "./data/chroma_tools")
    embedding_model = os.getenv("OLLAMA_EMBEDDING_MODEL", "nomic-embed-text")
    ollama_url = os.getenv("OLLAMA_BASE_URL", "http://localhost:11434")
    
    logger.info(f"🔧 初始化工具检索器...")
    logger.info(f"   ChromaDB 目录: {chroma_dir}")
    logger.info(f"   Embedding 模型: {embedding_model}")
    logger.info(f"   Ollama URL: {ollama_url}")
    
    tool_retriever = ToolRetriever(
        persist_directory=chroma_dir,
        embedding_model=embedding_model,
        ollama_base_url=ollama_url
    )
    
    logger.info("✅ 工具检索器初始化成功")
except ImportError as e:
    logger.warning(f"⚠️ 工具检索器依赖未安装: {e}")
    logger.warning("   跳过工具检索器初始化（需要: pip install langchain-chroma langchain-ollama）")
except Exception as e:
    logger.error(f"❌ 工具检索器初始化失败: {e}", exc_info=True)
    logger.warning("   继续启动，但工具检索功能将不可用")

# 🔥 Step 3: 初始化工作流规划器（需要 agent 初始化后才能获取 LLM client）
# 这将在 agent 初始化后设置

# 配置 CORS（安全配置）
# 生产环境应该限制为特定域名
ALLOWED_ORIGINS = os.getenv("ALLOWED_ORIGINS", "*").split(",")
if ALLOWED_ORIGINS == ["*"]:
    logger.warning("⚠️  CORS 配置允许所有来源，生产环境应限制为特定域名")

app.add_middleware(
    CORSMiddleware,
    allow_origins=ALLOWED_ORIGINS,
    allow_credentials=True,
    allow_methods=["GET", "POST", "OPTIONS"],
    allow_headers=["Content-Type", "Authorization"],
)

# 创建上传目录（使用绝对路径，兼容容器环境）
# 🔧 修复：优先使用环境变量，否则使用容器内的默认路径
# 注意：在容器环境中，默认路径应该是 /app/uploads，而不是相对路径
UPLOAD_DIR = Path(os.getenv("UPLOAD_DIR", "/app/uploads"))
RESULTS_DIR = Path(os.getenv("RESULTS_DIR", "/app/results"))

# 确保目录存在且可写
try:
    UPLOAD_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    
    # 检查目录权限
    if not os.access(UPLOAD_DIR, os.W_OK):
        logger.warning(f"⚠️ 上传目录不可写: {UPLOAD_DIR}")
    if not os.access(RESULTS_DIR, os.W_OK):
        logger.warning(f"⚠️ 结果目录不可写: {RESULTS_DIR}")
    
    logger.info(f"📁 上传目录: {UPLOAD_DIR.absolute()} (可写: {os.access(UPLOAD_DIR, os.W_OK)})")
    logger.info(f"📁 结果目录: {RESULTS_DIR.absolute()} (可写: {os.access(RESULTS_DIR, os.W_OK)})")
except Exception as e:
    logger.error(f"❌ 创建目录失败: {e}", exc_info=True)
    raise

# 安全配置
MAX_FILE_SIZE = int(os.getenv("MAX_FILE_SIZE", 100 * 1024 * 1024))  # 默认 100MB
ALLOWED_EXTENSIONS = {'.h5ad', '.mtx', '.tsv', '.csv', '.txt', '.gz', '.tar', '.zip'}
ALLOWED_MIME_TYPES = {
    'text/csv', 'text/tab-separated-values', 'text/plain',
    'application/gzip', 'application/x-gzip',
    'application/zip', 'application/x-tar'
}

def sanitize_filename(filename: str) -> str:
    """
    清理文件名，防止路径遍历攻击
    
    Args:
        filename: 原始文件名
    
    Returns:
        清理后的安全文件名
    """
    if not filename:
        # 如果文件名为空，生成随机名称
        return f"file_{secrets.token_hex(8)}"
    
    # 移除路径分隔符和危险字符
    filename = os.path.basename(filename)  # 移除路径部分
    filename = re.sub(r'[<>:"|?*\x00-\x1f]', '', filename)  # 移除危险字符
    filename = filename.strip('. ')  # 移除开头和结尾的点/空格
    
    # 如果清理后为空，生成随机名称
    if not filename:
        filename = f"file_{secrets.token_hex(8)}"
    
    # 限制文件名长度
    if len(filename) > 255:
        name, ext = os.path.splitext(filename)
        filename = name[:255-len(ext)] + ext
    
    return filename

def validate_file_path(file_path: Path, base_dir: Path) -> Path:
    """
    验证文件路径是否在允许的目录内（防止路径遍历）
    
    Args:
        file_path: 要验证的路径
        base_dir: 基础目录
    
    Returns:
        规范化的安全路径
    
    Raises:
        HTTPException: 如果路径不安全
    """
    try:
        # 解析并规范化路径
        resolved_path = file_path.resolve()
        resolved_base = base_dir.resolve()
        
        # 检查路径是否在基础目录内
        if not str(resolved_path).startswith(str(resolved_base)):
            raise HTTPException(
                status_code=403,
                detail="文件路径不安全：不允许访问基础目录外的文件"
            )
        
        return resolved_path
    except (ValueError, OSError) as e:
        raise HTTPException(
            status_code=400,
            detail=f"无效的文件路径: {str(e)}"
        )

# 初始化文件检测器
file_inspector = FileInspector(str(UPLOAD_DIR))

# 添加静态文件服务（用于访问结果图片）
from fastapi.staticfiles import StaticFiles
app.mount("/results", StaticFiles(directory="results"), name="results")
app.mount("/uploads", StaticFiles(directory="uploads"), name="uploads")

# 初始化智能体
agent = None
try:
    # 尝试从当前目录加载配置
    config_path = Path(__file__).parent / "gibh_agent" / "config" / "settings.yaml"
    logger.info(f"🔍 查找配置文件: {config_path}")
    logger.info(f"📂 配置文件存在: {config_path.exists()}")
    
    if not config_path.exists():
        # 如果不存在，尝试其他路径
        alt_path = Path(__file__).parent / "config" / "settings.yaml"
        logger.info(f"🔍 尝试备用路径: {alt_path}")
        if alt_path.exists():
            config_path = alt_path
        else:
            config_path = "gibh_agent/config/settings.yaml"
            logger.info(f"🔍 使用默认路径: {config_path}")
    
    logger.info(f"📄 使用配置文件: {config_path}")
    
    # 设置 scanpy 工具的默认输出目录（使用相对路径）
    import os
    scanpy_output_dir = os.path.join(os.getcwd(), "results")
    logger.info(f"📁 Scanpy 输出目录: {scanpy_output_dir}")
    
    # 创建智能体
    agent = create_agent(str(config_path))
    
    # 更新 scanpy 工具的输出目录
    if agent and hasattr(agent, 'agents') and 'rna_agent' in agent.agents:
        rna_agent = agent.agents['rna_agent']
        if hasattr(rna_agent, 'scanpy_tool'):
            rna_agent.scanpy_tool.output_dir = scanpy_output_dir
            os.makedirs(scanpy_output_dir, exist_ok=True)
            logger.info(f"✅ 已设置 Scanpy 输出目录: {scanpy_output_dir}")
    
    logger.info("✅ GIBH-AGENT 初始化成功")
    
    # 🔥 Step 3: 初始化工作流规划器（需要 agent 和 tool_retriever）
    if agent and tool_retriever:
        try:
            from gibh_agent.core.planner import WorkflowPlanner
            # 获取 LLM client（从 agent 的某个智能体中获取）
            llm_client = None
            if hasattr(agent, 'agents') and agent.agents:
                # 尝试从第一个智能体获取 LLM client
                first_agent = list(agent.agents.values())[0]
                if hasattr(first_agent, 'llm_client'):
                    llm_client = first_agent.llm_client
            
            if llm_client:
                workflow_planner = WorkflowPlanner(
                    tool_retriever=tool_retriever,
                    llm_client=llm_client
                )
                logger.info("✅ 工作流规划器初始化成功")
            else:
                logger.warning("⚠️ 无法获取 LLM client，跳过工作流规划器初始化")
        except Exception as e:
            logger.error(f"❌ 工作流规划器初始化失败: {e}", exc_info=True)
            logger.warning("   继续启动，但动态规划功能将不可用")
    
except Exception as e:
    import traceback
    error_msg = f"❌ GIBH-AGENT 初始化失败: {e}"
    logger.error(error_msg, exc_info=True)
    logger.error(f"详细错误:\n{traceback.format_exc()}")
    agent = None

# 🔥 Step 2: 启动时同步工具到 ChromaDB
@app.on_event("startup")
async def sync_tools_on_startup():
    """
    启动时同步工具到 Vector Database
    
    确保 ChromaDB 中的工具定义与代码中的 @register 装饰器保持一致。
    """
    # 🔥 Step 4: 首先加载所有工具模块（自动发现）
    try:
        logger.info("🔍 启动时自动发现和加载工具模块...")
        load_result = load_all_tools()
        logger.info(f"✅ 工具模块加载完成: {load_result['loaded']} 个成功, {load_result['failed']} 个失败")
    except Exception as e:
        logger.error(f"❌ 工具模块加载失败: {e}", exc_info=True)
        logger.warning("   继续启动，但工具可能未完全加载")
    
    # 然后同步到 ChromaDB
    if tool_retriever is None:
        logger.warning("⚠️ 工具检索器未初始化，跳过工具同步")
        return
    
    try:
        logger.info("🔄 启动时同步工具到 ChromaDB...")
        synced_count = tool_retriever.sync_tools(clear_existing=True)
        logger.info(f"✅ 工具同步完成: {synced_count} 个工具已同步到 ChromaDB")
    except Exception as e:
        logger.error(f"❌ 工具同步失败: {e}", exc_info=True)
        logger.warning("   继续启动，但工具检索功能可能不可用")


# 请求模型
class ChatRequest(BaseModel):
    message: str = ""
    history: List[dict] = []
    uploaded_files: List[dict] = []
    workflow_data: Optional[dict] = None
    test_dataset_id: Optional[str] = None


# 日志缓冲区（保留用于未来扩展）
log_buffer = deque(maxlen=1000)
log_listeners: Set[asyncio.Queue] = set()


def log_handler(record):
    """日志处理器，将日志发送到所有监听者"""
    log_entry = {
        "timestamp": datetime.now().isoformat(),
        "level": record.levelname,
        "message": record.getMessage(),
        "module": record.name
    }
    log_buffer.append(log_entry)
    
    # 通知所有监听者
    for listener in list(log_listeners):
        try:
            listener.put_nowait(log_entry)
        except:
            # 如果队列已满或已关闭，移除监听者
            log_listeners.discard(listener)


# 添加自定义日志处理器
class StreamLogHandler(logging.Handler):
    def emit(self, record):
        try:
            # 确保记录被格式化
            self.format(record)
            log_handler(record)
        except Exception as e:
            # 避免日志处理器本身出错，但记录错误
            print(f"日志处理器错误: {e}")


stream_handler = StreamLogHandler()
stream_handler.setLevel(logging.DEBUG)  # 降低级别以捕获更多日志
stream_handler.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))

# 添加到根日志记录器，捕获所有模块的日志
root_logger = logging.getLogger()
root_logger.setLevel(logging.DEBUG)  # 降低级别
# 移除现有的处理器，避免重复
for handler in root_logger.handlers[:]:
    root_logger.removeHandler(handler)
root_logger.addHandler(stream_handler)

# 也添加到当前logger
if stream_handler not in logger.handlers:
    logger.addHandler(stream_handler)

# 测试日志
logger.info("📋 日志系统初始化完成")
logger.info("🔍 测试日志输出 - 这应该出现在前端")


@app.get("/", response_class=HTMLResponse)
async def index():
    """返回前端页面"""
    # 🔥 优先读取外部 HTML 文件（如果存在）
    html_file_path = Path(__file__).parent / "services" / "nginx" / "html" / "index.html"
    if html_file_path.exists():
        try:
            with open(html_file_path, "r", encoding="utf-8") as f:
                html_content = f.read()
            logger.info(f"✅ 已加载外部前端文件: {html_file_path}")
            return HTMLResponse(content=html_content)
        except Exception as e:
            logger.warning(f"⚠️ 读取外部 HTML 文件失败，使用内嵌版本: {e}")
    
    # 如果外部文件不存在，使用内嵌的 HTML
    html_content = """
<!DOCTYPE html>
<html lang="zh-CN">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>GIBH-AGENT-V2 测试界面</title>
    <style>
        * {
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }
        body {
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            background: #f5f5f5;
            padding: 20px;
        }
        .container {
            max-width: 1200px;
            margin: 0 auto;
            display: flex;
            flex-direction: column;
            gap: 20px;
            height: calc(100vh - 40px);
        }
        .panel {
            background: white;
            border-radius: 8px;
            padding: 20px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
            display: flex;
            flex-direction: column;
        }
        .panel h2 {
            margin-bottom: 15px;
            color: #333;
            border-bottom: 2px solid #4CAF50;
            padding-bottom: 10px;
        }
        .chat-panel {
            flex: 1;
        }
        .chat-area {
            flex: 1;
            border: 1px solid #ddd;
            border-radius: 4px;
            padding: 15px;
            overflow-y: auto;
            overflow-x: hidden;
            margin-bottom: 15px;
            background: #fafafa;
            min-height: 300px;
            word-wrap: break-word;
            overflow-wrap: break-word;
        }
        .message {
            margin-bottom: 10px;
            padding: 8px;
            border-radius: 4px;
            word-wrap: break-word;
            overflow-wrap: break-word;
            max-width: 100%;
            box-sizing: border-box;
            overflow-x: hidden;
        }
        .message.user {
            background: #e3f2fd;
            text-align: right;
        }
        .message.assistant {
            background: #f1f8e9;
        }
        .message.error {
            background: #ffebee;
            color: #c62828;
        }
        .input-area {
            display: flex;
            gap: 10px;
        }
        input[type="text"], input[type="file"] {
            flex: 1;
            padding: 10px;
            border: 1px solid #ddd;
            border-radius: 4px;
            font-size: 14px;
        }
        button {
            padding: 10px 20px;
            background: #4CAF50;
            color: white;
            border: none;
            border-radius: 4px;
            cursor: pointer;
            font-size: 14px;
        }
        button:hover {
            background: #45a049;
        }
        button:disabled {
            background: #ccc;
            cursor: not-allowed;
        }
        .file-info {
            margin-top: 10px;
            padding: 10px;
            background: #fff3cd;
            border-radius: 4px;
            font-size: 12px;
        }
        .analysis-result {
            background: #f1f8e9 !important;
            word-wrap: break-word;
            overflow-wrap: break-word;
            max-width: 100%;
            box-sizing: border-box;
            overflow-x: hidden;
        }
        .analysis-summary {
            padding: 15px;
            max-width: 100%;
            box-sizing: border-box;
            overflow-x: hidden;
        }
        .analysis-summary h3 {
            margin-top: 0;
            color: #4CAF50;
        }
        .analysis-summary h4 {
            margin-top: 15px;
            margin-bottom: 10px;
            color: #333;
            border-bottom: 1px solid #ddd;
            padding-bottom: 5px;
        }
        .analysis-summary ul {
            margin: 10px 0;
            padding-left: 20px;
            max-width: 100%;
            box-sizing: border-box;
        }
        .analysis-summary li {
            margin: 5px 0;
            word-wrap: break-word;
            overflow-wrap: break-word;
            max-width: 100%;
        }
        .qc-metrics, .steps-details, .visualization, .step-plots, .markers-table, .diagnosis {
            max-width: 100%;
            box-sizing: border-box;
            overflow-x: hidden;
            word-wrap: break-word;
            overflow-wrap: break-word;
        }
        .diagnosis div {
            max-width: 100%;
            word-wrap: break-word;
            overflow-wrap: break-word;
            white-space: pre-wrap;
        }
        .visualization, .step-plots {
            max-width: 100%;
            box-sizing: border-box;
            overflow-x: hidden;
        }
        .visualization img, .step-plots img {
            max-width: 100%;
            height: auto;
            border: 1px solid #ddd;
            border-radius: 4px;
            margin: 10px 0;
            display: block;
        }
        .step-plots > div {
            max-width: 100%;
            box-sizing: border-box;
            overflow-x: hidden;
            word-wrap: break-word;
        }
        .markers-table {
            overflow-x: auto;
            max-width: 100%;
            box-sizing: border-box;
            margin: 10px 0;
        }
        .markers-table table {
            width: 100%;
            max-width: 100%;
            border-collapse: collapse;
            margin: 0;
            table-layout: auto;
        }
        .markers-table th, .markers-table td {
            word-wrap: break-word;
            overflow-wrap: break-word;
            max-width: 200px;
        }
        .markers-table th, .markers-table td {
            border: 1px solid #ddd;
            padding: 8px;
            text-align: left;
        }
        .markers-table th {
            background: #f5f5f5;
            font-weight: bold;
        }
        .think-card {
            background: #f1f8e9 !important;
        }
        .think-process {
            margin-bottom: 10px;
        }
        .think-header {
            background: #e8f5e9;
            padding: 10px 15px;
            border-radius: 4px;
            cursor: pointer;
            display: flex;
            align-items: center;
            gap: 10px;
            user-select: none;
            transition: background 0.2s;
        }
        .think-header:hover {
            background: #c8e6c9;
        }
        .think-icon {
            font-size: 18px;
        }
        .think-title {
            flex: 1;
            font-weight: bold;
            color: #2e7d32;
        }
        .think-toggle {
            color: #666;
            font-size: 12px;
        }
        .think-content {
            margin-top: 10px;
            padding: 15px;
            background: #fff;
            border: 1px solid #ddd;
            border-radius: 4px;
            white-space: pre-wrap;
            font-family: 'Courier New', monospace;
            font-size: 13px;
            line-height: 1.6;
            color: #333;
            max-height: 500px;
            overflow-y: auto;
        }
        .final-answer {
            margin-top: 10px;
            padding: 10px;
        }
        .test-data-selection {
            background: #f1f8e9 !important;
            max-width: 100%;
            box-sizing: border-box;
            overflow-x: hidden;
            word-wrap: break-word;
            overflow-wrap: break-word;
        }
        .test-data-selection h3 {
            margin-top: 0;
            color: #4CAF50;
            word-wrap: break-word;
        }
        .test-data-selection div[onclick] {
            border: 1px solid #ddd;
            border-radius: 4px;
            padding: 15px;
            cursor: pointer;
            transition: background 0.2s;
            max-width: 100%;
            box-sizing: border-box;
            word-wrap: break-word;
            overflow-wrap: break-word;
        }
        .test-data-selection div[onclick]:hover {
            background: #f5f5f5;
            border-color: #4CAF50;
        }
        .dataset-card {
            max-width: 100%;
            box-sizing: border-box;
            word-wrap: break-word;
            overflow-wrap: break-word;
        }
    </style>
</head>
<body>
    <div class="container">
        <div class="panel chat-panel">
            <h2>💬 对话界面</h2>
            <div id="chatArea" class="chat-area"></div>
            <div class="input-area">
                <input type="text" id="messageInput" placeholder="输入消息或上传文件进行分析..." />
                <input type="file" id="fileInput" accept=".h5ad,.mtx,.tsv,.csv" multiple />
                <button id="sendBtn" onclick="sendMessage()">发送</button>
            </div>
            <div id="fileInfo" class="file-info" style="display:none;"></div>
        </div>
    </div>

    <script>
        // 文件上下文管理（记住已上传的文件）
        let uploadedFilesContext = [];
        
        // 文件选择（支持多文件）
        let selectedFiles = [];
        document.getElementById('fileInput').addEventListener('change', function(e) {
            const files = Array.from(e.target.files);
            if (files.length > 0) {
                selectedFiles = files;
                const fileList = files.map(f => `${f.name} (${(f.size / 1024 / 1024).toFixed(2)} MB)`).join('<br>');
                document.getElementById('fileInfo').style.display = 'block';
                document.getElementById('fileInfo').innerHTML = `📁 已选择 ${files.length} 个文件:<br>${fileList}`;
            }
        });

        // 发送消息
        async function sendMessage() {
            const input = document.getElementById('messageInput');
            const message = input.value.trim();
            const btn = document.getElementById('sendBtn');
            
            if (!message && selectedFiles.length === 0) {
                alert('请输入消息或选择文件');
                return;
            }

            btn.disabled = true;
            const fileNames = selectedFiles.length > 0 ? selectedFiles.map(f => f.name).join(', ') : '';
            addMessage('user', message || (fileNames ? `上传文件: ${fileNames}` : ''));

            try {
                let uploadedFiles = [];
                
                // 如果有新选择的文件，先上传所有文件
                if (selectedFiles.length > 0) {
                    for (const file of selectedFiles) {
                        const formData = new FormData();
                        formData.append('file', file);
                        
                        const uploadRes = await fetch('/api/upload', {
                            method: 'POST',
                            body: formData
                        });
                        
                        if (!uploadRes.ok) {
                            throw new Error(`文件上传失败: ${file.name}`);
                        }
                        
                        const uploadData = await uploadRes.json();
                        uploadedFiles.push(uploadData);
                        // 添加到上下文
                        uploadedFilesContext.push(uploadData);
                        addMessage('assistant', `✅ 文件上传成功: ${uploadData.file_name}`);
                    }
                } else if (uploadedFilesContext.length > 0) {
                    // 如果没有新文件，使用上下文中的文件
                    uploadedFiles = uploadedFilesContext;
                    addMessage('assistant', `📁 使用已上传的文件: ${uploadedFiles.map(f => f.file_name).join(', ')}`);
                }

                // 发送聊天请求
                const response = await fetch('/api/chat', {
                    method: 'POST',
                    headers: {
                        'Content-Type': 'application/json',
                    },
                    body: JSON.stringify({
                        message: message || (uploadedFiles.length > 0 ? '分析这个文件' : ''),
                        history: [],
                        uploaded_files: uploadedFiles
                    })
                });

                if (!response.ok) {
                    throw new Error(`请求失败: ${response.status}`);
                }

                const contentType = response.headers.get('content-type');
                
                if (contentType && contentType.includes('application/json')) {
                    const data = await response.json();
                    
                    if (data.type === 'workflow_config') {
                        // 执行工作流
                        addMessage('assistant', '🚀 开始执行分析流程...');
                        await executeWorkflow(data.workflow_data, data.file_paths);
                    } else {
                        addMessage('assistant', JSON.stringify(data, null, 2));
                    }
                } else {
                    // 流式响应（支持 think 过程提取）
                    const reader = response.body.getReader();
                    const decoder = new TextDecoder();
                    let fullText = '';
                    let thinkBuffer = '';
                    let isThinking = false;
                    let hasThinkBlock = false;
                    let finalAnswer = '';
                    let thinkStartIndex = -1;
                    let datasetsJson = null;
                    
                    while (true) {
                        const { done, value } = await reader.read();
                        if (done) break;
                        
                        const chunk = decoder.decode(value);
                        fullText += chunk;
                        
                        // 检查是否包含数据集 JSON（测试数据选择响应）
                        // 使用非贪婪匹配，但需要匹配多行（因为 JSON 可能跨行）
                        const datasetsMatch = fullText.match(/<!-- DATASETS_JSON: (\[[\s\S]*?\]) -->/);
                        if (datasetsMatch && !datasetsJson) {
                            try {
                                // JSON 中的换行符已被替换为空格，直接解析即可
                                datasetsJson = JSON.parse(datasetsMatch[1]);
                            } catch (e) {
                                console.error('解析数据集JSON失败:', e, datasetsMatch[1].substring(0, 100));
                            }
                        }
                        
                        // 检测 think 开始标签（支持多种格式）
                        const thinkStartPatterns = [
                            /<think>/i,
                            /<think>/i,
                            /<reasoning>/i,
                            /<thought>/i,
                            /<thinking>/i
                        ];
                        
                        for (const pattern of thinkStartPatterns) {
                            const match = fullText.match(pattern);
                            if (match && !hasThinkBlock) {
                            isThinking = true;
                            hasThinkBlock = true;
                                thinkStartIndex = match.index + match[0].length;
                            // 创建 think 卡片
                            if (!document.querySelector('.think-card:last-child .think-process')) {
                                createThinkCard();
                                }
                                break;
                            }
                        }
                        
                        // 检测 think 结束标签
                        const thinkEndPatterns = [
                            /<\/think>/i,
                            /<\/redacted_reasoning>/i,
                            /<\/reasoning>/i,
                            /<\/thought>/i,
                            /<\/thinking>/i
                        ];
                        
                        for (const pattern of thinkEndPatterns) {
                            const match = fullText.match(pattern);
                            if (match && isThinking) {
                                // 提取 think 内容
                                thinkBuffer = fullText.substring(thinkStartIndex, match.index);
                            updateThinkContent(thinkBuffer);
                            isThinking = false;
                            
                                // 提取 think 标签之后的内容作为最终答案
                                const afterThinkIndex = match.index + match[0].length;
                                finalAnswer = fullText.substring(afterThinkIndex);
                                if (finalAnswer.trim()) {
                                    updateLastMessage('assistant', finalAnswer.trim());
                            }
                                break;
                        }
                        }
                        
                        // 更新显示
                        if (isThinking) {
                            // 在 think 块中，更新 think 内容
                            if (thinkStartIndex >= 0) {
                                thinkBuffer = fullText.substring(thinkStartIndex);
                            updateThinkContent(thinkBuffer);
                            }
                        } else if (hasThinkBlock && !isThinking) {
                            // think 块已结束，更新最终答案
                            if (finalAnswer) {
                                updateLastMessage('assistant', finalAnswer);
                            }
                        } else {
                            // 没有 think 块，直接更新消息
                            // 在流式响应过程中，先显示文本内容（去除 JSON 注释）
                            const cleanText = fullText.replace(/<!-- DATASETS_JSON: \[[\s\S]*?\] -->/g, '').trim();
                            updateLastMessage('assistant', cleanText);
                        }
                    }
                    
                    // 流式响应结束后，如果检测到数据集信息，替换为选择界面
                    if (datasetsJson && datasetsJson.length > 0) {
                        // 移除 JSON 注释，只保留用户友好的文本
                        const cleanText = fullText.replace(/<!-- DATASETS_JSON: \[[\s\S]*?\] -->/g, '').trim();
                        // 移除之前的普通消息
                        const lastMessage = chatArea.querySelector('.message.assistant:last-child');
                        if (lastMessage && !lastMessage.classList.contains('test-data-selection')) {
                            lastMessage.remove();
                        }
                        // 显示选择界面
                        if (!document.querySelector('.test-data-selection')) {
                            displayTestDataSelection(cleanText, datasetsJson);
                        }
                    }
                }
            } catch (error) {
                addMessage('error', `❌ 错误: ${error.message}`);
                console.error(error);
            } finally {
                btn.disabled = false;
                input.value = '';
                // 不清空 selectedFiles，保留文件选择
                // 但清空文件输入框，允许用户重新选择
                document.getElementById('fileInput').value = '';
                // 如果有上下文文件，显示提示
                if (uploadedFilesContext.length > 0) {
                    document.getElementById('fileInfo').style.display = 'block';
                    document.getElementById('fileInfo').innerHTML = `📁 已上传 ${uploadedFilesContext.length} 个文件，可直接输入需求继续分析`;
                } else {
                    document.getElementById('fileInfo').style.display = 'none';
                }
            }
        }

        // 执行工作流
        async function executeWorkflow(workflowData, filePaths) {
            try {
                const response = await fetch('/api/execute', {
                    method: 'POST',
                    headers: {
                        'Content-Type': 'application/json',
                    },
                    body: JSON.stringify({
                        workflow_data: workflowData,
                        file_paths: filePaths
                    })
                });

                const data = await response.json();
                
                if (data.status === 'success') {
                    // 美化显示分析结果
                    displayAnalysisResult(data);
                } else {
                    addMessage('error', `❌ 分析失败: ${data.error || '未知错误'}`);
                }
            } catch (error) {
                addMessage('error', `❌ 执行错误: ${error.message}`);
            }
        }
        
        // 美化显示分析结果
        function displayAnalysisResult(data) {
            const resultDiv = document.createElement('div');
            resultDiv.className = 'message assistant analysis-result';
            
            let html = '<div class="analysis-summary">';
            html += '<h3>✅ 分析完成</h3>';
            
            // QC 指标
            if (data.qc_metrics) {
                html += '<div class="qc-metrics">';
                html += '<h4>📊 质量控制指标</h4>';
                html += '<ul>';
                html += `<li>原始细胞数: <strong>${data.qc_metrics.raw_cells || 'N/A'}</strong></li>`;
                html += `<li>原始基因数: <strong>${data.qc_metrics.raw_genes || 'N/A'}</strong></li>`;
                if (data.qc_metrics.filtered_cells) {
                    html += `<li>过滤后细胞数: <strong>${data.qc_metrics.filtered_cells}</strong></li>`;
                }
                if (data.qc_metrics.filtered_genes) {
                    html += `<li>过滤后基因数: <strong>${data.qc_metrics.filtered_genes}</strong></li>`;
                }
                html += '</ul>';
                html += '</div>';
            }
            
            // 步骤详情
            if (data.steps_details && data.steps_details.length > 0) {
                html += '<div class="steps-details">';
                html += '<h4>📋 执行步骤</h4>';
                html += '<ul>';
                data.steps_details.forEach(step => {
                    const stepName = step.name || step.tool_id || '未知步骤';
                    const stepSummary = step.summary || '完成';
                    html += `<li><strong>${stepName}</strong>: ${stepSummary}</li>`;
                });
                html += '</ul>';
                html += '</div>';
            }
            
            // 可视化图片（只显示步骤的图片，避免与 final_plot 重复）
            if (data.steps_details) {
                const plotSteps = data.steps_details.filter(s => s.plot);
                if (plotSteps.length > 0) {
                    html += '<div class="step-plots">';
                    html += '<h4>📈 可视化结果</h4>';
                    plotSteps.forEach(step => {
                        let plotUrl = step.plot;
                        if (!plotUrl.startsWith('http') && !plotUrl.startsWith('/')) {
                            // 如果路径包含 results，直接使用
                            if (plotUrl.includes('results/')) {
                                plotUrl = '/' + plotUrl;
                            } else {
                                plotUrl = '/results/' + plotUrl;
                            }
                        }
                        html += `<div style="margin: 10px 0;">`;
                        html += `<strong>${step.name || step.tool_id}</strong><br>`;
                        html += `<img src="${plotUrl}" alt="${step.name}" style="max-width: 100%; border-radius: 4px; margin-top: 10px;" onerror="this.style.display='none';">`;
                        html += `</div>`;
                    });
                    html += '</div>';
                } else if (data.final_plot) {
                    // 如果没有步骤图片，使用 final_plot（向后兼容）
                    html += '<div class="visualization">';
                    html += '<h4>📈 可视化结果</h4>';
                    let plotUrl = data.final_plot;
                    if (!plotUrl.startsWith('http') && !plotUrl.startsWith('/')) {
                        if (plotUrl.includes('results/')) {
                            plotUrl = '/' + plotUrl;
                        } else {
                            plotUrl = '/results/' + plotUrl;
                        }
                    }
                    html += `<img src="${plotUrl}" alt="Visualization" style="max-width: 100%; border-radius: 4px; margin-top: 10px;" onerror="this.style.display='none'; this.nextElementSibling.style.display='block';"><p style="display:none; color: #999;">图片加载失败: ${plotUrl}</p>`;
                    html += '</div>';
                }
            } else if (data.final_plot) {
                // 如果没有 steps_details，使用 final_plot
                html += '<div class="visualization">';
                html += '<h4>📈 可视化结果</h4>';
                let plotUrl = data.final_plot;
                if (!plotUrl.startsWith('http') && !plotUrl.startsWith('/')) {
                    if (plotUrl.includes('results/')) {
                        plotUrl = '/' + plotUrl;
                    } else {
                        plotUrl = '/results/' + plotUrl;
                    }
                }
                html += `<img src="${plotUrl}" alt="Visualization" style="max-width: 100%; border-radius: 4px; margin-top: 10px;" onerror="this.style.display='none'; this.nextElementSibling.style.display='block';"><p style="display:none; color: #999;">图片加载失败: ${plotUrl}</p>`;
                html += '</div>';
            }
            
            // Marker 基因表格（如果有）
            const markersStep = data.steps_details?.find(s => s.name === 'local_markers' || s.tool_id === 'local_markers');
            if (markersStep && markersStep.details) {
                html += '<div class="markers-table">';
                html += '<h4>🧬 Marker 基因</h4>';
                // 直接显示 HTML 表格
                html += markersStep.details;
                html += '</div>';
            }
            
            // 诊断信息
            if (data.diagnosis) {
                html += '<div class="diagnosis">';
                html += '<h4>💡 分析诊断</h4>';
                html += `<div style="white-space: pre-wrap;">${data.diagnosis}</div>`;
                html += '</div>';
            }
            
            html += '</div>';
            resultDiv.innerHTML = html;
            
            const chatArea = document.getElementById('chatArea');
            chatArea.appendChild(resultDiv);
            chatArea.scrollTop = chatArea.scrollHeight;
        }

        // 添加消息
        function addMessage(role, content) {
            const chatArea = document.getElementById('chatArea');
            const msgDiv = document.createElement('div');
            msgDiv.className = `message ${role}`;
            msgDiv.textContent = content;
            chatArea.appendChild(msgDiv);
            chatArea.scrollTop = chatArea.scrollHeight;
        }

        // 更新最后一条消息
        function updateLastMessage(role, content) {
            const chatArea = document.getElementById('chatArea');
            const messages = chatArea.querySelectorAll('.message');
            if (messages.length > 0 && messages[messages.length - 1].classList.contains(role)) {
                const lastMsg = messages[messages.length - 1];
                // 如果已经有 think 卡片，更新最终答案部分
                const finalAnswerDiv = lastMsg.querySelector('.final-answer');
                if (finalAnswerDiv) {
                    finalAnswerDiv.textContent = content;
                } else {
                    lastMsg.textContent = content;
                }
            } else {
                addMessage(role, content);
            }
            chatArea.scrollTop = chatArea.scrollHeight;
        }
        
        // 创建 think 卡片
        function createThinkCard() {
            const chatArea = document.getElementById('chatArea');
            const thinkCard = document.createElement('div');
            thinkCard.className = 'message assistant think-card';
            thinkCard.innerHTML = `
                <div class="think-process">
                    <div class="think-header" onclick="toggleThink(this)">
                        <span class="think-icon">🤔</span>
                        <span class="think-title">思考过程</span>
                        <span class="think-toggle">▼</span>
                    </div>
                    <div class="think-content" style="display: none;"></div>
                </div>
                <div class="final-answer"></div>
            `;
            chatArea.appendChild(thinkCard);
            chatArea.scrollTop = chatArea.scrollHeight;
        }
        
        // 更新 think 内容
        function updateThinkContent(content) {
            const chatArea = document.getElementById('chatArea');
            const thinkCards = chatArea.querySelectorAll('.think-card');
            if (thinkCards.length > 0) {
                const lastCard = thinkCards[thinkCards.length - 1];
                const thinkContentDiv = lastCard.querySelector('.think-content');
                if (thinkContentDiv) {
                    thinkContentDiv.textContent = content;
                }
            }
        }
        
        // 切换 think 卡片展开/折叠
        function toggleThink(header) {
            const thinkCard = header.closest('.think-process');
            const content = thinkCard.querySelector('.think-content');
            const toggle = header.querySelector('.think-toggle');
            
            if (content.style.display === 'none') {
                content.style.display = 'block';
                toggle.textContent = '▲';
            } else {
                content.style.display = 'none';
                toggle.textContent = '▼';
            }
        }
        
        // 显示测试数据选择界面
        function displayTestDataSelection(messageText, datasets) {
            const chatArea = document.getElementById('chatArea');
            
            // 检查是否已经显示过选择界面
            const existing = document.querySelector('.test-data-selection');
            if (existing) {
                // 如果已存在，更新它而不是创建新的
                return;
            }
            
            // 移除之前的普通消息（如果有）
            const lastMessage = chatArea.querySelector('.message.assistant:last-child');
            if (lastMessage && !lastMessage.classList.contains('test-data-selection')) {
                lastMessage.remove();
            }
            
            const selectionDiv = document.createElement('div');
            selectionDiv.className = 'message assistant test-data-selection';
            
            let html = '<div style="padding: 15px; max-width: 100%; box-sizing: border-box;">';
            html += '<h3 style="margin-top: 0; color: #4CAF50; margin-bottom: 10px; word-wrap: break-word;">📊 选择测试数据集</h3>';
            
            // 显示消息文本（去除数据集列表部分）
            const lines = messageText.split('\\n');
            const messageLine = lines.find(line => line.includes('检测到') || line.includes('请选择'));
            if (messageLine) {
                html += '<p style="margin-bottom: 15px; color: #333; word-wrap: break-word; max-width: 100%;">' + messageLine + '</p>';
            }
            
            // 显示数据集选择卡片
            html += '<div style="display: flex; flex-direction: column; gap: 10px; margin-bottom: 15px; max-width: 100%; box-sizing: border-box;">';
            datasets.forEach(dataset => {
                html += `<div class="dataset-card" 
                             style="border: 2px solid #ddd; border-radius: 8px; padding: 15px; cursor: pointer; transition: all 0.2s; background: white; max-width: 100%; box-sizing: border-box; word-wrap: break-word; overflow-wrap: break-word;" 
                             onmouseover="this.style.borderColor='#4CAF50'; this.style.boxShadow='0 2px 8px rgba(76,175,80,0.2)'" 
                             onmouseout="this.style.borderColor='#ddd'; this.style.boxShadow='none'"
                             onclick="selectTestDataset('${dataset.id}', '${dataset.name}')">`;
                html += `<div style="display: flex; align-items: center; gap: 10px; margin-bottom: 8px; max-width: 100%; flex-wrap: wrap;">`;
                html += `<span style="font-size: 24px; flex-shrink: 0;">📦</span>`;
                html += `<strong style="color: #4CAF50; font-size: 18px; word-wrap: break-word; flex: 1; min-width: 0;">${dataset.name}</strong>`;
                html += `</div>`;
                html += `<p style="margin: 0; color: #666; font-size: 14px; word-wrap: break-word; max-width: 100%;">${dataset.description}</p>`;
                html += `<div style="margin-top: 8px; font-size: 12px; color: #999; word-wrap: break-word; max-width: 100%;">ID: <code style="word-break: break-all;">${dataset.id}</code></div>`;
                html += '</div>';
            });
            html += '</div>';
            
            html += '<p style="margin-top: 10px; color: #666; font-size: 14px; font-style: italic; word-wrap: break-word; max-width: 100%;">💡 点击上面的数据集卡片选择，或上传您自己的数据文件。</p>';
            html += '</div>';
            
            selectionDiv.innerHTML = html;
            chatArea.appendChild(selectionDiv);
            chatArea.scrollTop = chatArea.scrollHeight;
        }
        
        // 选择测试数据集
        async function selectTestDataset(datasetId, datasetName) {
            addMessage('user', `使用测试数据集: ${datasetName} (${datasetId})`);
            addMessage('assistant', `正在使用测试数据集 ${datasetName} 执行分析...`);
            
            try {
                const response = await fetch('/api/chat', {
                    method: 'POST',
                    headers: {
                        'Content-Type': 'application/json',
                    },
                    body: JSON.stringify({
                        message: `使用测试数据集 ${datasetId} 执行完整的单细胞转录组分析流程`,
                        history: [],
                        uploaded_files: [],
                        test_dataset_id: datasetId
                    })
                });
                
                // 处理响应
                const contentType = response.headers.get('content-type');
                if (contentType && contentType.includes('application/json')) {
                    const data = await response.json();
                    if (data.type === 'workflow_config') {
                        addMessage('assistant', '🚀 开始执行分析流程...');
                        await executeWorkflow(data.workflow_data, data.file_paths || []);
                    } else {
                        addMessage('assistant', JSON.stringify(data, null, 2));
                    }
                } else {
                    // 流式响应
                    const reader = response.body.getReader();
                    const decoder = new TextDecoder();
                    let fullText = '';
                    
                    while (true) {
                        const { done, value } = await reader.read();
                        if (done) break;
                        const chunk = decoder.decode(value);
                        fullText += chunk;
                        updateLastMessage('assistant', fullText);
                    }
                }
            } catch (error) {
                addMessage('error', `❌ 错误: ${error.message}`);
                console.error(error);
            }
        }
        
        // 全局函数，供 HTML 调用
        window.toggleThink = toggleThink;
        window.selectTestDataset = selectTestDataset;

        // 日志功能已移除

        // 回车发送
        document.getElementById('messageInput').addEventListener('keypress', function(e) {
            if (e.key === 'Enter') {
                sendMessage();
            }
        });

        // 初始化完成
    </script>
</body>
</html>
    """
    return HTMLResponse(content=html_content)


@app.post("/api/upload")
async def upload_file(files: List[UploadFile] = File(...)):
    """文件上传接口（支持多文件上传）"""
    try:
        if not files or len(files) == 0:
            raise HTTPException(status_code=400, detail="No files provided")
        
        # 限制文件数量
        if len(files) > 20:
            raise HTTPException(status_code=400, detail="一次最多上传20个文件")
        
        logger.info(f"📤 收到文件上传: {len(files)} 个文件")
        
        # 检测是否是10x Genomics文件（matrix.mtx, barcodes.tsv, features.tsv）
        is_10x_data = False
        tenx_files = []
        other_files = []
        
        for file in files:
            # 🔒 安全：清理文件名
            if not file.filename:
                raise HTTPException(status_code=400, detail="文件名不能为空")
            
            original_filename = file.filename
            safe_filename = sanitize_filename(original_filename)
            
            # 🔒 安全：验证文件扩展名
            file_ext = Path(safe_filename).suffix.lower()
            if file_ext and file_ext not in ALLOWED_EXTENSIONS:
                raise HTTPException(
                    status_code=400,
                    detail=f"不允许的文件类型: {file_ext}。允许的类型: {', '.join(ALLOWED_EXTENSIONS)}"
                )
            
            # 更新文件名为安全版本
            file.filename = safe_filename
            
            filename_lower = safe_filename.lower()
            if any(pattern in filename_lower for pattern in ['matrix.mtx', 'barcodes.tsv', 'features.tsv']):
                is_10x_data = True
                tenx_files.append(file)
            else:
                other_files.append(file)
        
        uploaded_results = []
        
        # 如果是10x数据，创建子目录并保存
        if is_10x_data and len(tenx_files) >= 2:  # 至少需要2个文件（通常是matrix + barcodes/features）
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            tenx_dir = UPLOAD_DIR / f"10x_data_{timestamp}"
            tenx_dir.mkdir(exist_ok=True)
            
            logger.info(f"📁 检测到10x数据，创建目录: {tenx_dir}")
            
            # 保存10x文件到子目录
            for file in tenx_files:
                # 🔒 安全：验证文件路径
                file_path = tenx_dir / file.filename
                try:
                    file_path = validate_file_path(file_path, UPLOAD_DIR)
                except HTTPException as e:
                    logger.error(f"❌ 文件路径验证失败: {file.filename} -> {e.detail}")
                    raise
                
                # 🔒 安全：检查文件大小
                content = await file.read()
                if len(content) > MAX_FILE_SIZE:
                    raise HTTPException(
                        status_code=413,
                        detail=f"文件 {file.filename} 超过最大大小限制 ({MAX_FILE_SIZE / 1024 / 1024:.0f}MB)"
                    )
                
                # 🔧 修复：确保父目录存在
                file_path.parent.mkdir(parents=True, exist_ok=True)
                
                try:
                    with open(file_path, "wb") as f:
                        f.write(content)
                except PermissionError as e:
                    logger.error(f"❌ 文件写入权限错误: {file_path} -> {e}")
                    raise HTTPException(status_code=500, detail=f"文件保存失败：权限不足 ({file.filename})")
                except OSError as e:
                    logger.error(f"❌ 文件写入系统错误: {file_path} -> {e}")
                    raise HTTPException(status_code=500, detail=f"文件保存失败：{str(e)} ({file.filename})")
                
                logger.info(f"✅ 10x文件保存成功: {file_path}")
                
                # 生成元数据
                try:
                    metadata = file_inspector.generate_metadata(str(file_path.relative_to(UPLOAD_DIR)))
                except Exception as e:
                    logger.warning(f"⚠️ 生成文件元数据失败: {e}")
                    metadata = None
                
                uploaded_results.append({
                    "file_id": str(tenx_dir.relative_to(UPLOAD_DIR)),
                    "file_name": file.filename,
                    "file_path": str(file_path),
                    "file_size": len(content),
                    "metadata": metadata,
                    "is_10x": True,
                    "group_dir": str(tenx_dir.relative_to(UPLOAD_DIR))
                })
            
            # 返回10x目录路径（而不是单个文件路径）
            file_paths = [str(tenx_dir.relative_to(UPLOAD_DIR))]
            return {
                "status": "success",
                "is_10x_data": True,
                "group_dir": str(tenx_dir.relative_to(UPLOAD_DIR)),
                "files": uploaded_results,
                "file_paths": file_paths,  # 🔥 添加 file_paths 数组
                "message": f"10x数据已保存到: {tenx_dir.relative_to(UPLOAD_DIR)}"
            }
        
        # 处理其他文件（非10x或单独的10x文件）
        # 🔧 修复：如果只有1个10x文件，也当作普通文件处理
        files_to_process = other_files
        if is_10x_data and len(tenx_files) == 1:
            # 只有1个10x文件，当作普通文件处理
            logger.info(f"⚠️ 只有1个10x文件，当作普通文件处理: {tenx_files[0].filename}")
            files_to_process = other_files + tenx_files
        elif not is_10x_data:
            # 不是10x数据，处理所有文件
            files_to_process = other_files + tenx_files
        
        for file in files_to_process:
            # 🔒 安全：验证文件路径
            file_path = UPLOAD_DIR / file.filename
            try:
                file_path = validate_file_path(file_path, UPLOAD_DIR)
            except HTTPException as e:
                logger.error(f"❌ 文件路径验证失败: {file.filename} -> {e.detail}")
                raise
            
            # 🔒 安全：检查文件大小
            content = await file.read()
            if len(content) > MAX_FILE_SIZE:
                raise HTTPException(
                    status_code=413,
                    detail=f"文件 {file.filename} 超过最大大小限制 ({MAX_FILE_SIZE / 1024 / 1024:.0f}MB)"
                )
            
            # 🔧 修复：确保父目录存在
            file_path.parent.mkdir(parents=True, exist_ok=True)
            
            try:
                with open(file_path, "wb") as f:
                    f.write(content)
            except PermissionError as e:
                logger.error(f"❌ 文件写入权限错误: {file_path} -> {e}")
                raise HTTPException(status_code=500, detail=f"文件保存失败：权限不足 ({file.filename})")
            except OSError as e:
                logger.error(f"❌ 文件写入系统错误: {file_path} -> {e}")
                raise HTTPException(status_code=500, detail=f"文件保存失败：{str(e)} ({file.filename})")
            
            logger.info(f"✅ 文件保存成功: {file_path}")
            
            # 生成元数据
            try:
                metadata = file_inspector.generate_metadata(file.filename)
                if metadata:
                    logger.info(f"📊 文件元数据已生成: {metadata.get('file_type', 'unknown')}")
            except Exception as e:
                logger.warning(f"⚠️ 生成文件元数据失败: {e}")
                metadata = None
            
            uploaded_results.append({
                "file_id": file.filename,
                "file_name": file.filename,
                "file_path": str(file_path),
                "file_size": len(content),
                "metadata": metadata,
                "is_10x": False
            })
        
        # 🔥 统一返回格式：始终返回 file_paths 数组和 file_info 数组（用于前端发送聊天请求）
        # 注意：使用相对路径，因为前端需要相对于 UPLOAD_DIR 的路径
        file_paths = []
        file_info = []
        for result in uploaded_results:
            # 转换为相对路径（相对于 UPLOAD_DIR）
            file_path_abs = result["file_path"]
            if isinstance(file_path_abs, str) and str(UPLOAD_DIR) in file_path_abs:
                # 提取相对路径
                rel_path = str(Path(file_path_abs).relative_to(UPLOAD_DIR))
            else:
                rel_path = result["file_id"]  # 回退到 file_id
            file_paths.append(rel_path)
            
            # 构建 file_info 条目
            file_info.append({
                "name": result["file_name"],
                "size": result["file_size"],
                "path": rel_path  # 使用相对路径
            })
        
        # 🔥 统一返回格式：始终返回一致的 JSON 结构
        response = {
            "status": "success",
            "file_paths": file_paths,  # 文件路径数组（相对路径）
            "file_info": file_info,    # 文件信息数组
            "count": len(uploaded_results)
        }
        
        # 如果只有一个文件，添加单个文件的详细信息（向后兼容）
        if len(uploaded_results) == 1:
            result = uploaded_results[0]
            response.update({
                "file_id": result["file_id"],
                "file_name": result["file_name"],
                "file_path": result["file_path"],  # 绝对路径（向后兼容）
                "file_size": result["file_size"],
                "metadata": result["metadata"]
            })
        else:
            # 多个文件时，添加 files 数组（向后兼容）
            response["files"] = uploaded_results
        
        return response
        
    except HTTPException:
        # 重新抛出 HTTP 异常（保持状态码和详细信息）
        raise
    except Exception as e:
        # 🔧 改进：记录详细错误信息，但返回用户友好的错误消息
        import traceback
        error_detail = f"{type(e).__name__}: {str(e)}"
        logger.error(f"❌ 文件上传失败: {error_detail}", exc_info=True)
        logger.error(f"详细堆栈:\n{traceback.format_exc()}")
        
        # 根据错误类型返回更具体的错误信息
        if "Permission" in error_detail or "permission" in error_detail.lower():
            raise HTTPException(status_code=500, detail="文件上传失败：权限不足，请检查服务器配置")
        elif "No such file" in error_detail or "directory" in error_detail.lower():
            raise HTTPException(status_code=500, detail="文件上传失败：目录不存在，请检查服务器配置")
        elif "disk" in error_detail.lower() or "space" in error_detail.lower():
            raise HTTPException(status_code=500, detail="文件上传失败：磁盘空间不足")
        else:
            # 开发环境返回更详细的错误信息，生产环境返回通用错误
            import os
            if os.getenv("DEBUG", "false").lower() == "true":
                raise HTTPException(status_code=500, detail=f"文件上传失败: {error_detail}")
            else:
                raise HTTPException(status_code=500, detail="文件上传失败，请稍后重试")


@app.post("/api/chat")
async def chat_endpoint(req: ChatRequest):
    """聊天接口"""
    # #region debug log - entry point
    import json
    import traceback
    # 🔧 修复：使用容器内的日志路径（统一使用 /app/debug.log）
    debug_log_path = Path("/app/debug.log")
    try:
        # 确保目录存在
        debug_log_path.parent.mkdir(parents=True, exist_ok=True)
        with open(debug_log_path, 'a') as f:
            f.write(json.dumps({"location":"server.py:1112","message":"chat_endpoint entry","data":{"agent_is_none":agent is None,"req_message":req.message[:100] if req.message else None},"timestamp":int(__import__('time').time()*1000),"sessionId":"debug-session","runId":"run1","hypothesisId":"ENTRY"})+"\n")
    except Exception as log_err:
        pass  # 即使日志写入失败也不影响主流程
    # #endregion
    
    if not agent:
        error_msg = "智能体未初始化，请检查配置和日志。可能的原因：1) 配置文件路径错误 2) API Key未设置 3) 依赖包缺失"
        logger.error(error_msg)
        logger.error("请检查终端日志中的详细错误信息")
        return JSONResponse(
            status_code=500,
            content={
                "type": "error",
                "error": error_msg,
                "message": "智能体初始化失败，请查看服务器日志获取详细信息"
            }
        )
    
    try:
        logger.info(f"💬 收到聊天请求: {req.message}")
        logger.info(f"📁 上传文件数: {len(req.uploaded_files)}")
        logger.info(f"🔄 工作流数据: {req.workflow_data is not None}")
        
        # 🔧 修复：如果包含工作流数据，直接执行工作流（而不是通过 agent.process_query）
        if req.workflow_data:
            logger.info("🚀 检测到工作流执行请求，直接调用 execute_workflow")
            try:
                # 🔧 修复：优先使用 workflow_data 中的 file_paths（前端已经设置好）
                file_paths = req.workflow_data.get("file_paths", [])
                logger.info(f"📁 从 workflow_data 获取的文件路径: {file_paths}")
                
                # 如果 workflow_data 中没有 file_paths，再从 uploaded_files 中提取
                if not file_paths:
                    logger.info("⚠️ workflow_data 中没有 file_paths，从 uploaded_files 中提取")
                    for file_info in req.uploaded_files:
                        file_name = file_info.get("file_name", "")
                        file_path_str = file_info.get("file_path", "")
                        
                        if file_path_str:
                            file_path = Path(file_path_str)
                        else:
                            file_path = UPLOAD_DIR / file_name if file_name else None
                        
                        if file_path and file_path.exists():
                            file_paths.append(str(file_path))
                
                logger.info(f"📂 最终文件路径列表: {file_paths}")
                
                # 验证文件路径是否存在
                valid_file_paths = []
                for fp in file_paths:
                    fp_path = Path(fp)
                    if fp_path.exists():
                        valid_file_paths.append(str(fp_path))
                    else:
                        logger.warning(f"⚠️ 文件不存在，跳过: {fp}")
                
                if not valid_file_paths:
                    raise ValueError("没有找到有效的输入文件。请确保文件已正确上传。")
                
                # 直接调用 execute_workflow 函数（不通过 HTTP）
                execute_request = {
                    "workflow_data": req.workflow_data,
                    "file_paths": valid_file_paths
                }
                # 调用 execute_workflow 函数（定义在下方）
                result = await execute_workflow(execute_request)
                return result
            except Exception as e:
                logger.error(f"❌ 工作流执行失败: {e}", exc_info=True)
                return JSONResponse(
                    status_code=500,
                    content={
                        "type": "error",
                        "error": str(e),
                        "message": f"工作流执行失败: {str(e)}"
                    }
                )
        
        # 🔥 Step 3: 尝试使用动态规划器（如果可用且查询看起来是工作流规划请求）
        if workflow_planner and not req.workflow_data:
            # 简单的启发式检测：如果查询包含分析相关的关键词，尝试使用规划器
            query_lower = req.message.lower()
            workflow_keywords = [
                "analyze", "analysis", "pca", "differential", "preprocess",
                "分析", "处理", "降维", "差异", "预处理"
            ]
            
            # 如果有上传文件或包含关键词，尝试使用规划器
            has_files = len(req.uploaded_files) > 0
            has_keywords = any(keyword in query_lower for keyword in workflow_keywords)
            
            if has_files or has_keywords:
                try:
                    logger.info("🧠 尝试使用动态规划器生成工作流...")
                    
                    # 提取文件路径（先转换 uploaded_files）
                    file_paths = []
                    for file_info in req.uploaded_files:
                        file_path = file_info.get("path") or file_info.get("file_name")
                        if file_path:
                            # 如果是相对路径，转换为绝对路径
                            if not Path(file_path).is_absolute():
                                file_path = str(UPLOAD_DIR / Path(file_path).name)
                            file_paths.append(file_path)
                    
                    # 检测类别（简单启发式）
                    category_filter = None
                    if any(keyword in query_lower for keyword in ["metabolite", "代谢", "metabolomics"]):
                        category_filter = "Metabolomics"
                    elif any(keyword in query_lower for keyword in ["rna", "gene", "transcript", "转录"]):
                        category_filter = "scRNA-seq"
                    
                    # 调用规划器
                    plan_result = await workflow_planner.plan(
                        user_query=req.message,
                        context_files=file_paths,
                        category_filter=category_filter
                    )
                    
                    # 如果规划成功，返回结果
                    if plan_result.get("type") == "workflow_config":
                        logger.info("✅ 动态规划器成功生成工作流")
                        return JSONResponse(content=plan_result)
                    else:
                        logger.info(f"⚠️ 动态规划器返回: {plan_result.get('type')}，继续使用传统流程")
                        # 继续使用传统流程
                except Exception as planner_err:
                    logger.warning(f"⚠️ 动态规划器失败，回退到传统流程: {planner_err}")
                    # 继续使用传统流程
        
        # 🔥 转换文件路径：支持多种前端格式
        uploaded_files = []
        logger.info(f"📥 收到 uploaded_files: {len(req.uploaded_files)} 个文件")
        
        for file_info in req.uploaded_files:
            # 支持多种字段名：file_name/name, file_path/path
            file_name = file_info.get("file_name") or file_info.get("name", "")
            file_path_str = file_info.get("file_path") or file_info.get("path", "")
            
            # 🔒 安全：清理文件名
            if file_name:
                file_name = sanitize_filename(file_name)
            
            # 🔥 构建文件路径：优先使用 file_path，如果是相对路径则拼接 UPLOAD_DIR
            if file_path_str:
                file_path = Path(file_path_str)
                # 如果是相对路径，拼接 UPLOAD_DIR
                if not file_path.is_absolute():
                    file_path = UPLOAD_DIR / file_path
            elif file_name:
                # 如果没有路径，使用文件名在 UPLOAD_DIR 中查找
                file_path = UPLOAD_DIR / file_name
            else:
                logger.warning(f"⚠️ 无法确定文件路径，跳过: {file_info}")
                continue
            
            # 🔒 安全：验证路径在允许的目录内
            try:
                file_path = validate_file_path(file_path, UPLOAD_DIR)
            except HTTPException:
                logger.warning(f"⚠️ 不安全的文件路径，跳过: {file_path}")
                continue
            
            # 检查文件是否存在
            if not file_path.exists():
                logger.warning(f"⚠️ 文件不存在: {file_path}，尝试查找...")
                # 尝试在 UPLOAD_DIR 中查找同名文件
                if file_name:
                    alt_path = UPLOAD_DIR / file_name
                    if alt_path.exists():
                        file_path = alt_path
                        logger.info(f"✅ 找到文件: {file_path}")
                    else:
                        logger.warning(f"⚠️ 文件不存在，跳过: {file_path}")
                        continue
                else:
                    continue
            
            uploaded_files.append({
                "name": file_name or os.path.basename(str(file_path)),
                "path": str(file_path)
            })
        
        logger.info(f"📂 处理文件: {len(uploaded_files)} 个有效文件")
        logger.info(f"📂 文件路径列表: {[f['path'] for f in uploaded_files]}")
        
        # 处理查询
        # #region debug log
        try:
            debug_log_path = Path("/app/debug.log")
            debug_log_path.parent.mkdir(parents=True, exist_ok=True)
            with open(debug_log_path, 'a') as f:
                f.write(json.dumps({"location":"server.py:1161","message":"Before process_query","data":{"query":req.message,"uploaded_files_count":len(uploaded_files),"test_dataset_id":req.test_dataset_id},"timestamp":int(__import__('time').time()*1000),"sessionId":"debug-session","runId":"run1","hypothesisId":"A"})+"\n")
        except Exception as log_err:
            pass  # 日志写入失败不影响主流程
        # #endregion
        try:
            result = await agent.process_query(
                query=req.message,
                history=req.history,
                uploaded_files=uploaded_files,
                test_dataset_id=req.test_dataset_id
            )
        except Exception as process_err:
            # #region debug log
            try:
                debug_log_path = Path("/app/debug.log")
                debug_log_path.parent.mkdir(parents=True, exist_ok=True)
                with open(debug_log_path, 'a') as f:
                    f.write(json.dumps({"location":"server.py:1156","message":"process_query exception","data":{"error_type":type(process_err).__name__,"error_message":str(process_err),"traceback":traceback.format_exc()},"timestamp":int(__import__('time').time()*1000),"sessionId":"debug-session","runId":"run1","hypothesisId":"PROCESS_QUERY"})+"\n")
            except:
                pass
            # #endregion
            raise  # 重新抛出异常，让外层异常处理捕获
        
        # #region debug log
        try:
            debug_log_path = Path("/app/debug.log")
            debug_log_path.parent.mkdir(parents=True, exist_ok=True)
            with open(debug_log_path, 'a') as f:
                f.write(json.dumps({"location":"server.py:1168","message":"After process_query","data":{"result_type":type(result).__name__,"result_keys":list(result.keys()) if isinstance(result,dict) else None,"result_type_value":result.get('type') if isinstance(result,dict) else None},"timestamp":int(__import__('time').time()*1000),"sessionId":"debug-session","runId":"run1","hypothesisId":"A"})+"\n")
        except:
            pass
        # #endregion
        
        logger.info(f"✅ 处理完成，返回类型: {result.get('type', 'unknown')}")
        
        # 如果是工作流配置，返回 JSON
        if result.get("type") == "workflow_config":
            # 优先使用 result 中的 file_paths（可能来自测试数据集）
            # 如果没有，则使用 uploaded_files
            result_file_paths = result.get("file_paths", [])
            if not result_file_paths:
                result_file_paths = [f["path"] for f in uploaded_files]
            
            response_content = {
                "type": "workflow_config",
                "workflow_data": result.get("workflow_data"),
                "file_paths": result_file_paths
            }
            
            # 🔧 修复：如果包含诊断报告，也返回给前端
            if "diagnosis_report" in result:
                response_content["diagnosis_report"] = result["diagnosis_report"]
            
            # 🔧 修复：如果包含推荐信息，也返回给前端（代谢组学）
            if "recommendation" in result:
                response_content["recommendation"] = result["recommendation"]
            
            logger.info(f"📤 返回工作流配置: 包含推荐={('recommendation' in response_content)}, 包含诊断={('diagnosis_report' in response_content)}")
            
            return JSONResponse(content=response_content)
        
        # 如果是测试数据选择请求，格式化为用户友好的文本
        if result.get("type") == "test_data_selection":
            # #region debug log
            import json
            try:
                debug_log_path = Path("/app/debug.log")
                debug_log_path.parent.mkdir(parents=True, exist_ok=True)
                with open(debug_log_path, 'a') as f:
                    f.write(json.dumps({"location":"server.py:1178","message":"Entering test_data_selection handler","data":{"has_message":"message" in result,"has_options":"options" in result,"has_datasets_display":"datasets_display" in result,"has_datasets_json":"datasets_json" in result},"timestamp":int(__import__('time').time()*1000),"sessionId":"debug-session","runId":"run1","hypothesisId":"B"})+"\n")
            except:
                pass  # 日志写入失败不影响主流程
            # #endregion
            async def generate():
                try:
                    # #region debug log
                    try:
                        debug_log_path = Path("/app/debug.log")
                        debug_log_path.parent.mkdir(parents=True, exist_ok=True)
                        with open(debug_log_path, 'a') as f:
                            f.write(json.dumps({"location":"server.py:1181","message":"Inside generate()","data":{},"timestamp":int(__import__('time').time()*1000),"sessionId":"debug-session","runId":"run1","hypothesisId":"C"})+"\n")
                    except:
                        pass
                    # #endregion
                    # 构建用户友好的消息
                    message = result.get("message", "检测到您没有上传相关数据。请选择：")
                    options = result.get("options", [])
                    datasets_display = result.get("datasets_display", "")
                    # #region debug log
                    try:
                        debug_log_path = Path("/app/debug.log")
                        debug_log_path.parent.mkdir(parents=True, exist_ok=True)
                        with open(debug_log_path, 'a') as f:
                            f.write(json.dumps({"location":"server.py:1187","message":"Before datasets_json processing","data":{"message_type":type(message).__name__,"options_type":type(options).__name__,"datasets_display_type":type(datasets_display).__name__,"datasets_display_len":len(str(datasets_display)) if datasets_display else 0},"timestamp":int(__import__('time').time()*1000),"sessionId":"debug-session","runId":"run1","hypothesisId":"D"})+"\n")
                    except:
                        pass
                    # #endregion
                    
                    response_text = f"{message}\n\n"
                    for option in options:
                        response_text += f"  {option}\n"
                    
                    if datasets_display:
                        response_text += f"\n{datasets_display}\n"
                    
                    response_text += "\n💡 提示：回复数据集ID（如：pbmc_1k_v3）或上传您自己的数据文件。\n"
                    
                    # 同时保存数据集信息到响应中（用于前端处理）
                    # 这里我们通过特殊标记来传递 JSON 数据
                    # 将 JSON 中的换行符替换为空格，避免破坏 HTML 注释
                    datasets_json_raw = result.get('datasets_json', '[]')
                    # #region debug log
                    try:
                        debug_log_path = Path("/app/debug.log")
                        debug_log_path.parent.mkdir(parents=True, exist_ok=True)
                        with open(debug_log_path, 'a') as f:
                            f.write(json.dumps({"location":"server.py:1200","message":"Before datasets_json replace","data":{"datasets_json_type":type(datasets_json_raw).__name__,"datasets_json_is_none":datasets_json_raw is None,"datasets_json_len":len(str(datasets_json_raw)) if datasets_json_raw else 0},"timestamp":int(__import__('time').time()*1000),"sessionId":"debug-session","runId":"run1","hypothesisId":"B"})+"\n")
                    except:
                        pass
                    # #endregion
                    datasets_json = str(datasets_json_raw).replace('\n', ' ').replace('\r', '') if datasets_json_raw else '[]'
                    # #region debug log
                    try:
                        debug_log_path = Path("/app/debug.log")
                        debug_log_path.parent.mkdir(parents=True, exist_ok=True)
                        with open(debug_log_path, 'a') as f:
                            f.write(json.dumps({"location":"server.py:1203","message":"After datasets_json replace","data":{"datasets_json_len":len(datasets_json),"response_text_len":len(response_text)},"timestamp":int(__import__('time').time()*1000),"sessionId":"debug-session","runId":"run1","hypothesisId":"B"})+"\n")
                    except:
                        pass
                    # #endregion
                    response_text += f"\n<!-- DATASETS_JSON: {datasets_json} -->\n"
                    
                    # #region debug log
                    try:
                        debug_log_path = Path("/app/debug.log")
                        debug_log_path.parent.mkdir(parents=True, exist_ok=True)
                        with open(debug_log_path, 'a') as f:
                            f.write(json.dumps({"location":"server.py:1207","message":"Before yield","data":{"final_response_text_len":len(response_text)},"timestamp":int(__import__('time').time()*1000),"sessionId":"debug-session","runId":"run1","hypothesisId":"C"})+"\n")
                    except:
                        pass
                    # #endregion
                    yield response_text
                    # #region debug log
                    try:
                        debug_log_path = Path("/app/debug.log")
                        debug_log_path.parent.mkdir(parents=True, exist_ok=True)
                        with open(debug_log_path, 'a') as f:
                            f.write(json.dumps({"location":"server.py:1209","message":"After yield","data":{},"timestamp":int(__import__('time').time()*1000),"sessionId":"debug-session","runId":"run1","hypothesisId":"C"})+"\n")
                    except:
                        pass
                    # #endregion
                except Exception as e:
                    # #region debug log
                    import traceback
                    try:
                        debug_log_path = Path("/app/debug.log")
                        debug_log_path.parent.mkdir(parents=True, exist_ok=True)
                        with open(debug_log_path, 'a') as f:
                            f.write(json.dumps({"location":"server.py:1212","message":"Exception in generate()","data":{"error_type":type(e).__name__,"error_message":str(e),"traceback":traceback.format_exc()},"timestamp":int(__import__('time').time()*1000),"sessionId":"debug-session","runId":"run1","hypothesisId":"C"})+"\n")
                    except:
                        pass
                    # #endregion
                    logger.error(f"❌ 格式化测试数据选择响应错误: {e}", exc_info=True)
                    yield f"\n\n❌ 错误: {str(e)}"
            
            # #region debug log
            try:
                debug_log_path = Path("/app/debug.log")
                debug_log_path.parent.mkdir(parents=True, exist_ok=True)
                with open(debug_log_path, 'a') as f:
                    f.write(json.dumps({"location":"server.py:1218","message":"Before StreamingResponse","data":{},"timestamp":int(__import__('time').time()*1000),"sessionId":"debug-session","runId":"run1","hypothesisId":"E"})+"\n")
            except:
                pass
            # #endregion
            return StreamingResponse(generate(), media_type="text/plain")
        
        # 如果是聊天响应，返回流式
        if result.get("type") == "chat":
            async def generate():
                try:
                    response_iter = result.get("response")
                    if response_iter:
                        # 确保 response_iter 是异步迭代器
                        async for chunk in response_iter:
                            if chunk:
                                yield chunk
                    else:
                        logger.warning("⚠️ 聊天响应中没有 response 迭代器")
                        yield "❌ 错误: 无法获取响应"
                except Exception as e:
                    logger.error(f"❌ 流式响应错误: {e}", exc_info=True)
                    import traceback
                    logger.error(f"详细错误: {traceback.format_exc()}")
                    yield f"\n\n❌ 错误: {str(e)}"
            
            return StreamingResponse(
                generate(), 
                media_type="text/plain",
                headers={
                    "Cache-Control": "no-cache",
                    "Connection": "keep-alive",
                    "X-Accel-Buffering": "no"
                }
            )
        
        # 其他情况返回 JSON
        return JSONResponse(content=result)
        
    except Exception as e:
        # #region debug log
        import traceback
        try:
            debug_log_path = Path("/app/debug.log")
            debug_log_path.parent.mkdir(parents=True, exist_ok=True)
            with open(debug_log_path, 'a') as f:
                f.write(json.dumps({"location":"server.py:1210","message":"Exception in chat_endpoint","data":{"error_type":type(e).__name__,"error_message":str(e),"traceback":traceback.format_exc()},"timestamp":int(__import__('time').time()*1000),"sessionId":"debug-session","runId":"run1","hypothesisId":"ALL"})+"\n")
        except:
            pass  # 日志写入失败不影响主流程
        # #endregion
        error_detail = f"{type(e).__name__}: {str(e)}"
        logger.error(f"❌ 处理失败: {error_detail}", exc_info=True)
        logger.error(f"详细错误: {traceback.format_exc()}")
        raise HTTPException(status_code=500, detail=error_detail)


@app.post("/api/execute")
async def execute_workflow(request: dict):
    """执行工作流接口"""
    if not agent:
        raise HTTPException(status_code=500, detail="智能体未初始化")
    
    try:
        workflow_data = request.get("workflow_data")
        file_paths = request.get("file_paths", [])
        
        logger.info(f"🚀 开始执行工作流: {len(file_paths)} 个文件")
        
        # 🔧 修复：优先检查 workflow_name 中是否包含代谢组关键词
        workflow_name = workflow_data.get("workflow_name", "")
        routing = None
        target_agent = None
        route_query = None
        
        # 方法1: 如果有 workflow_name，优先检查是否包含代谢组关键词
        if workflow_name:
            workflow_name_lower = workflow_name.lower()
            # 如果 workflow_name 包含代谢组关键词，直接路由到 metabolomics_agent
            if any(kw in workflow_name_lower for kw in ["metabolomics", "代谢组", "代谢"]):
                logger.info(f"✅ 根据 workflow_name 直接路由到 metabolomics_agent: {workflow_name}")
                routing = "metabolomics_agent"
                target_agent = agent.agents.get(routing)
                if not target_agent:
                    logger.warning(f"⚠️ metabolomics_agent 不存在，使用默认 rna_agent")
                    target_agent = agent.agents.get("rna_agent")
                    routing = "rna_agent"
            else:
                route_query = workflow_name
        # 方法2: 根据文件类型构建查询
        elif file_paths:
            file_path = file_paths[0]
            file_ext = os.path.splitext(file_path)[1].lower()
            if file_ext == ".csv":
                route_query = "metabolomics analysis"
            elif file_ext in [".h5ad", ".h5"]:
                route_query = "single cell transcriptomics analysis"
            elif "fastq" in file_path.lower():
                route_query = "single cell RNA-seq analysis"
            else:
                route_query = "bioinformatics analysis"
        else:
            route_query = "bioinformatics analysis"
        
        # 准备上传文件列表（用于 RouterAgent）
        uploaded_files_for_router = []
        for file_path in file_paths:
            uploaded_files_for_router.append({
                "name": os.path.basename(file_path),
                "path": file_path
            })
        
        # 🔧 修复：如果还没有路由，使用 RouterAgent 进行路由决策
        if not routing or not target_agent:
            try:
                route_result = await agent.router.process_query(
                    query=route_query,
                    history=[],
                    uploaded_files=uploaded_files_for_router
                )
                
                routing = route_result.get("routing", "rna_agent")
                target_agent = agent.agents.get(routing)
                
                # 如果路由的智能体不存在，使用默认的 RNA Agent
                if not target_agent:
                    logger.warning(f"⚠️ 路由的智能体不存在: {routing}，使用默认 rna_agent")
                    target_agent = agent.agents.get("rna_agent")
                    routing = "rna_agent"
                
                if not target_agent:
                    raise HTTPException(status_code=500, detail="RNA Agent 未找到")
                
                logger.info(f"✅ RouterAgent 路由结果: {routing} (confidence: {route_result.get('confidence', 0):.2f}, modality: {route_result.get('modality', 'unknown')})")
                
            except Exception as e:
                logger.error(f"❌ RouterAgent 路由失败: {e}，使用默认 rna_agent", exc_info=True)
                # 降级到默认 Agent
                target_agent = agent.agents.get("rna_agent")
                routing = "rna_agent"
                if not target_agent:
                    raise HTTPException(status_code=500, detail="RNA Agent 未找到")
        
        # 🔥 Step 4: 使用通用执行器（动态执行，不依赖硬编码逻辑）
        try:
            from gibh_agent.core.executor import WorkflowExecutor
            
            logger.info("🔧 使用通用执行器执行工作流...")
            
            # 设置输出目录
            output_dir = str(RESULTS_DIR / f"run_{datetime.now().strftime('%Y%m%d_%H%M%S')}")
            
            # 创建执行器并执行
            executor = WorkflowExecutor(output_dir=output_dir)
            report_data = executor.execute_workflow(
                workflow_data=workflow_data,
                file_paths=file_paths,
                output_dir=output_dir
            )
            
            logger.info("✅ 通用执行器执行完成")
            
            # 构建返回结果（符合前端格式）
            return JSONResponse(content={
                "type": "analysis_report",
                "status": "success",
                "report_data": report_data,
                "reply": "✅ 工作流执行完成（使用动态执行引擎）",
                "thought": "[THOUGHT] 使用 ToolRegistry 动态执行，工具无关"
            })
        
        except ImportError:
            logger.warning("⚠️ 通用执行器未找到，回退到传统执行方式")
            # 回退到传统执行方式
            output_dir = str(RESULTS_DIR / f"run_{datetime.now().strftime('%Y%m%d_%H%M%S')}")
            os.makedirs(output_dir, exist_ok=True)
            
            report = await target_agent.execute_workflow(
                workflow_config=workflow_data,
                file_paths=file_paths,
                output_dir=output_dir
            )
            
            logger.info(f"✅ 工作流执行完成: {report.get('status')}")
            
            # 处理图片路径（传统方式）
            if report.get("final_plot"):
                plot_path = report["final_plot"]
                if not plot_path.startswith("/results/"):
                    if plot_path.startswith("results/"):
                        plot_path = "/" + plot_path
                    elif "/" in plot_path:
                        plot_path = f"/results/{plot_path}"
                    else:
                        run_name = os.path.basename(output_dir)
                        plot_path = f"/results/{run_name}/{plot_path}"
                report["final_plot"] = plot_path
            
            # 处理步骤中的图片路径
            run_name = os.path.basename(output_dir)
            if report.get("steps_details"):
                for step in report["steps_details"]:
                    if step.get("plot"):
                        plot_path = step["plot"]
                        if not plot_path.startswith("/results/"):
                            if plot_path.startswith("results/"):
                                plot_path = "/" + plot_path
                            elif "/" in plot_path:
                                plot_path = f"/results/{plot_path}"
                            else:
                                plot_path = f"/results/{run_name}/{plot_path}"
                        step["plot"] = plot_path
            
            # 返回传统格式的结果
            return JSONResponse(content={
                "type": "analysis_report",
                "status": report.get("status", "success"),
                "report_data": report
            })
        
        # 处理通用执行器返回的图片路径
        logger.info(f"✅ 工作流执行完成: {report_data.get('status')}")
        
        # 处理图片路径，转换为可访问的 URL（在返回之前）
        # 图片保存在 results/run_xxx/ 目录，需要转换为 /results/run_xxx/filename
        if report.get("final_plot"):
            plot_path = report["final_plot"]
            # 确保路径以 /results/ 开头
            if not plot_path.startswith("/results/"):
                if plot_path.startswith("results/"):
                    plot_path = "/" + plot_path
                elif "/" in plot_path:
                    # 如果包含 run_xxx/filename 格式，添加 results 前缀
                    plot_path = f"/results/{plot_path}"
                else:
                    # 如果只是文件名，需要找到对应的 run 目录
                    # 从 output_dir 中提取 run_xxx
                    run_name = os.path.basename(output_dir)
                    plot_path = f"/results/{run_name}/{plot_path}"
            report_data["final_plot"] = plot_path
        
        # 处理步骤中的图片路径（steps_details）
        run_name = os.path.basename(output_dir)
        if report_data.get("steps_details"):
            for step in report_data["steps_details"]:
                if step.get("plot"):
                    plot_path = step["plot"]
                    # 确保路径以 /results/ 开头
                    if not plot_path.startswith("/results/"):
                        if plot_path.startswith("results/"):
                            plot_path = "/" + plot_path
                        elif "/" in plot_path:
                            plot_path = f"/results/{plot_path}"
                        else:
                            plot_path = f"/results/{run_name}/{plot_path}"
                    step["plot"] = plot_path
                
                # 处理 step_result 中的图片路径
                if step.get("step_result") and step["step_result"].get("data", {}).get("images"):
                    images = step["step_result"]["data"]["images"]
                    fixed_images = []
                    for img_path in images:
                        if not img_path.startswith("/results/"):
                            if img_path.startswith("results/"):
                                img_path = "/" + img_path
                            elif "/" in img_path:
                                img_path = f"/results/{img_path}"
                            else:
                                img_path = f"/results/{run_name}/{img_path}"
                        fixed_images.append(img_path)
                    step["step_result"]["data"]["images"] = fixed_images
        
        # 确保 steps_results 存在（前端可直接使用）
        if "steps_results" not in report and "steps_details" in report:
            steps_results = []
            for step_detail in report.get("steps_details", []):
                if "step_result" in step_detail:
                    step_result = step_detail["step_result"].copy()
                    # 确保图片路径正确
                    if step_result.get("data", {}).get("images"):
                        images = step_result["data"]["images"]
                        fixed_images = []
                        for img_path in images:
                            if not img_path.startswith("/results/"):
                                if img_path.startswith("results/"):
                                    img_path = "/" + img_path
                                elif "/" in img_path:
                                    img_path = f"/results/{img_path}"
                                else:
                                    img_path = f"/results/{run_name}/{img_path}"
                                fixed_images.append(img_path)
                            else:
                                fixed_images.append(img_path)
                        step_result["data"]["images"] = fixed_images
                    steps_results.append(step_result)
                else:
                    # 兼容旧格式
                    step_result = {
                        "step_name": step_detail.get("name", "Unknown"),
                        "status": step_detail.get("status", "success"),
                        "logs": step_detail.get("summary", ""),
                        "data": {}
                    }
                    # 如果有 plot，添加到 data.images
                    if step_detail.get("plot"):
                        plot_path = step_detail["plot"]
                        if not plot_path.startswith("/results/"):
                            if plot_path.startswith("results/"):
                                plot_path = "/" + plot_path
                            elif "/" in plot_path:
                                plot_path = f"/results/{plot_path}"
                            else:
                                plot_path = f"/results/{run_name}/{plot_path}"
                        step_result["data"]["images"] = [plot_path]
                    steps_results.append(step_result)
            report["steps_results"] = steps_results
        
        # 处理 steps_results 中的图片路径（如果存在）
        if report.get("steps_results"):
            for step_result in report["steps_results"]:
                if step_result.get("data", {}).get("images"):
                    images = step_result["data"]["images"]
                    fixed_images = []
                    for img_path in images:
                        if not img_path.startswith("/results/"):
                            if img_path.startswith("results/"):
                                img_path = "/" + img_path
                            elif "/" in img_path:
                                img_path = f"/results/{img_path}"
                            else:
                                img_path = f"/results/{run_name}/{img_path}"
                            fixed_images.append(img_path)
                        else:
                            fixed_images.append(img_path)
                    step_result["data"]["images"] = fixed_images
        
        # 🔧 修复：返回正确的工作流执行结果格式
        return JSONResponse(content={
            "type": "analysis_report",
            "status": "success",
            "report_data": report
        })
        
    except Exception as e:
        import traceback
        error_detail = f"{type(e).__name__}: {str(e)}"
        error_traceback = traceback.format_exc()
        logger.error(f"❌ 工作流执行失败: {error_detail}", exc_info=True)
        logger.error(f"详细错误: {error_traceback}")
        # 返回更详细的错误信息
        return JSONResponse(
            status_code=500,
            content={
                "status": "error",
                "error": error_detail,
                "error_detail": error_traceback,
                "message": f"工作流执行失败: {error_detail}"
            }
        )


@app.get("/api/logs/stream")
async def stream_logs():
    """实时日志流接口（Server-Sent Events）"""
    logger.info("📡 新的日志流连接")
    
    async def event_generator():
        q = asyncio.Queue(maxsize=100)
        log_listeners.add(q)
        
        try:
            # 先发送历史日志
            history_logs = list(log_buffer)[-100:]  # 最近100条
            logger.info(f"📤 发送历史日志: {len(history_logs)} 条")
            for entry in history_logs:
                yield f"data: {json.dumps(entry, ensure_ascii=False)}\\n\\n"
            
            # 实时发送新日志
            while True:
                try:
                    entry = await asyncio.wait_for(q.get(), timeout=1.0)
                    yield f"data: {json.dumps(entry, ensure_ascii=False)}\\n\\n"
                except asyncio.TimeoutError:
                    # 发送心跳保持连接
                    yield f"data: {json.dumps({'type': 'heartbeat', 'timestamp': datetime.now().isoformat()})}\\n\\n"
        except asyncio.CancelledError:
            logger.info("📡 日志流连接已取消")
        except Exception as e:
            logger.error(f"❌ 日志流错误: {e}", exc_info=True)
        finally:
            log_listeners.discard(q)
            logger.info("📡 日志流连接已关闭")
    
    return StreamingResponse(
        event_generator(),
        media_type="text/event-stream",
        headers={
            "Cache-Control": "no-cache",
            "Connection": "keep-alive",
            "X-Accel-Buffering": "no"
        }
    )


@app.get("/api/logs")
async def get_logs(limit: int = 100):
    """获取历史日志"""
    return JSONResponse(content={
        "logs": list(log_buffer)[-limit:],
        "total": len(log_buffer)
    })


# 🔥 Step 2: Tool-RAG API - 工具检索端点
@app.get("/api/tools/search")
async def search_tools(
    query: str,
    top_k: int = 5,
    category: Optional[str] = None
):
    """
    语义搜索工具
    
    Args:
        query: 查询文本（自然语言）
        top_k: 返回前 k 个最相关的工具（默认 5）
        category: 可选的类别过滤器
    
    Returns:
        相关工具的 JSON Schema 列表
    """
    if tool_retriever is None:
        raise HTTPException(
            status_code=503,
            detail="工具检索器未初始化。请检查 Ollama 服务和依赖是否已安装。"
        )
    
    try:
        tools = tool_retriever.retrieve(query=query, top_k=top_k, category_filter=category)
        return {
            "status": "success",
            "query": query,
            "count": len(tools),
            "tools": tools
        }
    except Exception as e:
        logger.error(f"❌ 工具搜索失败: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=f"工具搜索失败: {str(e)}")


@app.get("/api/tools/list")
async def list_tools():
    """
    列出所有已注册的工具
    
    Returns:
        工具名称列表
    """
    if tool_retriever is None:
        raise HTTPException(
            status_code=503,
            detail="工具检索器未初始化"
        )
    
    try:
        tools = tool_retriever.list_all_tools()
        return {
            "status": "success",
            "count": len(tools),
            "tools": tools
        }
    except Exception as e:
        logger.error(f"❌ 列出工具失败: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=f"列出工具失败: {str(e)}")


@app.get("/api/tools/{tool_name}")
async def get_tool_schema(tool_name: str):
    """
    获取特定工具的完整 Schema
    
    Args:
        tool_name: 工具名称
    
    Returns:
        工具的完整 JSON Schema
    """
    if tool_retriever is None:
        raise HTTPException(
            status_code=503,
            detail="工具检索器未初始化"
        )
    
    try:
        tool_schema = tool_retriever.get_tool_by_name(tool_name)
        if tool_schema is None:
            raise HTTPException(status_code=404, detail=f"工具 '{tool_name}' 不存在")
        
        return {
            "status": "success",
            "tool": tool_schema
        }
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"❌ 获取工具 Schema 失败: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=f"获取工具 Schema 失败: {str(e)}")


@app.get("/api/workflow/status/{run_id}")
async def get_workflow_status(run_id: str):
    """
    查询工作流状态（兼容旧架构）
    如果使用 Celery，查询 Celery 任务状态
    如果使用同步执行，返回 not_found（因为同步执行没有任务ID）
    """
    try:
        # 尝试从 Celery 查询任务状态
        from celery.result import AsyncResult
        from gibh_agent.core.celery_app import celery_app
        
        task_result = AsyncResult(run_id, app=celery_app)
        
        response = {
            "status": "running",
            "completed": False,
            "steps_status": [],
            "error": None
        }
        
        if task_result.state == 'PENDING':
            response["status"] = "running"
        elif task_result.state == 'SUCCESS':
            response["status"] = "success"
            response["completed"] = True
            result_data = task_result.result
            if result_data:
                response["report_data"] = result_data
                if "steps_details" in result_data:
                    response["steps_status"] = result_data["steps_details"]
                elif "steps" in result_data:
                    response["steps_status"] = result_data["steps"]
        elif task_result.state == 'FAILURE':
            response["status"] = "failed"
            response["completed"] = True
            response["error"] = str(task_result.result)
        elif task_result.state == 'PROGRESS':
            info = task_result.info
            if isinstance(info, dict):
                response["steps_status"] = info.get("steps", [])
        
        return JSONResponse(content=response)
        
    except ImportError:
        # Celery 未安装或未配置，返回 not_found
        return JSONResponse(
            status_code=404,
            content={
                "status": "not_found",
                "message": "工作流状态查询需要 Celery 支持，当前使用同步执行模式"
            }
        )
    except Exception as e:
        logger.error(f"❌ 查询工作流状态失败: {e}", exc_info=True)
        return JSONResponse(
            status_code=500,
            content={
                "status": "error",
                "error": str(e)
            }
        )


if __name__ == "__main__":
    import uvicorn
    import json
    
    port = int(os.getenv("PORT", 8018))
    logger.info(f"🚀 启动服务器，端口: {port}")
    logger.info(f"📁 上传目录: {UPLOAD_DIR.absolute()}")
    logger.info(f"📁 结果目录: {RESULTS_DIR.absolute()}")
    
    uvicorn.run(
        "server:app",
        host="0.0.0.0",
        port=port,
        log_level="info",
        reload=True
    )

