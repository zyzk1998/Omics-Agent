"""
工具检索器 - Vector Database Integration

将 ToolRegistry 中的工具同步到 ChromaDB，并提供语义搜索功能。
使用 LangChain ChromaDB 和 Ollama Embeddings。
"""
import os
import json
import logging
from typing import List, Dict, Any, Optional
from pathlib import Path

try:
    from langchain_chroma import Chroma
    from langchain_ollama import OllamaEmbeddings
    from langchain_core.documents import Document
except ImportError as e:
    logging.warning(f"⚠️ LangChain 依赖未安装: {e}")
    Chroma = None
    OllamaEmbeddings = None
    Document = None

from .tool_registry import registry

logger = logging.getLogger(__name__)


class ToolRetriever:
    """
    工具检索器
    
    职责：
    1. 将 ToolRegistry 中的工具同步到 ChromaDB
    2. 提供基于语义搜索的工具检索功能
    """
    
    def __init__(
        self,
        persist_directory: Optional[str] = None,
        embedding_model: str = "nomic-embed-text",
        ollama_base_url: Optional[str] = None
    ):
        """
        初始化工具检索器
        
        Args:
            persist_directory: ChromaDB 持久化目录
            embedding_model: Ollama embedding 模型名称（默认 "nomic-embed-text"）
            ollama_base_url: Ollama 服务 URL（默认从环境变量读取或使用 http://localhost:11434）
        """
        if Chroma is None or OllamaEmbeddings is None:
            raise ImportError(
                "❌ 缺少必要的依赖。请安装: pip install langchain-chroma langchain-ollama"
            )
        
        _pd = persist_directory or os.getenv("CHROMA_TOOLS_PERSIST_DIR", "./data/chroma_tools")
        self.persist_directory = Path(_pd)
        self.persist_directory.mkdir(parents=True, exist_ok=True)
        
        # 初始化 Embeddings
        ollama_url = ollama_base_url or os.getenv("OLLAMA_BASE_URL", "http://localhost:11434")
        logger.info(f"🔗 初始化 Ollama Embeddings: {ollama_url}, 模型: {embedding_model}")
        
        try:
            self.embeddings = OllamaEmbeddings(
                model=embedding_model,
                base_url=ollama_url
            )
        except Exception as e:
            logger.error(f"❌ 初始化 Ollama Embeddings 失败: {e}")
            logger.warning("⚠️ 尝试使用备用配置...")
            # 备用：尝试不同的配置
            try:
                self.embeddings = OllamaEmbeddings(model=embedding_model)
            except Exception as e2:
                logger.error(f"❌ 备用配置也失败: {e2}")
                raise
        
        # 初始化 ChromaDB Vector Store
        collection_name = "tools"
        try:
            self.vector_store = Chroma(
                collection_name=collection_name,
                embedding_function=self.embeddings,
                persist_directory=str(self.persist_directory)
            )
            logger.info(f"✅ ChromaDB 初始化成功: {self.persist_directory}")
        except Exception as e:
            logger.error(f"❌ ChromaDB 初始化失败: {e}", exc_info=True)
            raise
    
    def sync_tools(self, clear_existing: bool = True) -> int:
        """
        同步 ToolRegistry 中的工具到 ChromaDB
        
        🔥 关键功能：确保 VDB 与代码中的工具定义保持一致
        
        Args:
            clear_existing: 是否清除现有集合（默认 True，确保一致性）
        
        Returns:
            同步的工具数量
        """
        logger.info("🔄 开始同步工具到 ChromaDB...")
        
        # 获取所有工具
        tools_json = registry.get_all_tools_json()
        
        if not tools_json:
            logger.warning("⚠️ ToolRegistry 中没有工具，跳过同步")
            return 0
        
        logger.info(f"📋 发现 {len(tools_json)} 个工具需要同步")
        
        # 清除现有集合（如果需要）
        if clear_existing:
            try:
                # 删除现有集合并重新创建
                # 注意：Chroma 的 delete_collection 方法可能不存在，我们通过删除并重建来实现
                logger.info("🗑️  清除现有工具集合...")
                # 创建一个新的集合（会自动覆盖）
                self.vector_store = Chroma(
                    collection_name="tools",
                    embedding_function=self.embeddings,
                    persist_directory=str(self.persist_directory)
                )
            except Exception as e:
                logger.warning(f"⚠️ 清除集合时出错（可能集合不存在）: {e}")
        
        # 准备文档
        documents = []
        metadatas = []
        
        for tool in tools_json:
            # page_content: 用于搜索的文本（description + name + category）
            page_content = f"""
工具名称: {tool['name']}
类别: {tool['category']}
描述: {tool['description']}
输出类型: {tool['output_type']}
""".strip()
            
            # metadata: 完整的工具 schema（用于后续传递给 LLM）
            metadata = {
                "name": tool['name'],
                "category": tool['category'],
                "output_type": tool['output_type'],
                "description": tool['description'],
                "args_schema": json.dumps(tool['args_schema'], ensure_ascii=False)  # JSON 字符串
            }
            
            documents.append(Document(page_content=page_content, metadata=metadata))
            metadatas.append(metadata)

        # 附录：超算/工作站 MCP 固定清单（与 docs/hpc_mcp_tools_catalog.json 一致），便于向量库命中 hpc_mcp_*
        repo_root = Path(__file__).resolve().parents[2]
        mcp_catalog_path = repo_root / "docs" / "hpc_mcp_tools_catalog.json"
        if mcp_catalog_path.is_file():
            try:
                raw_cat = json.loads(mcp_catalog_path.read_text(encoding="utf-8"))
            except Exception as e:
                logger.warning("⚠️ 读取 hpc_mcp_tools_catalog.json 失败，跳过 MCP 向量条目: %s", e)
                raw_cat = []
            if isinstance(raw_cat, list):
                for row in raw_cat:
                    if not isinstance(row, dict):
                        continue
                    oname = (row.get("openai_function_name") or "").strip()
                    mcp_n = (row.get("mcp_name") or "").strip()
                    desc = (row.get("description") or "").strip()
                    if not oname:
                        continue
                    page_content = (
                        f"工具名称: {oname}\n"
                        f"远端 MCP 名: {mcp_n}\n"
                        f"类别: HPC MCP\n"
                        f"描述: {desc}\n"
                        f"输出类型: json"
                    ).strip()
                    metadata = {
                        "name": oname,
                        "category": "HPC MCP",
                        "output_type": "json",
                        "description": desc,
                        "args_schema": json.dumps(
                            {"type": "object", "properties": {}, "x_mcp_tool": mcp_n},
                            ensure_ascii=False,
                        ),
                    }
                    documents.append(Document(page_content=page_content, metadata=metadata))
                logger.info("📎 已附加 %s 条 MCP 目录工具供向量检索", len(raw_cat))
        
        # 添加到 ChromaDB
        try:
            # 使用 add_documents 方法
            ids = [f"tool_{i}" for i in range(len(documents))]
            self.vector_store.add_documents(documents=documents, ids=ids)
            
            # 持久化（新版本 ChromaDB 自动持久化，如果 persist 方法存在则调用）
            if hasattr(self.vector_store, 'persist'):
                self.vector_store.persist()
            else:
                # 新版本 ChromaDB 使用 persist_directory 时自动持久化
                logger.debug("ChromaDB 自动持久化（无需手动调用 persist）")
            
            logger.info(f"✅ 成功同步 {len(documents)} 个工具到 ChromaDB")
            return len(documents)
        
        except Exception as e:
            logger.error(f"❌ 同步工具失败: {e}", exc_info=True)
            raise
    
    def retrieve(
        self,
        query: str,
        top_k: int = 5,
        category_filter: Optional[str] = None
    ) -> List[Dict[str, Any]]:
        """
        检索相关工具
        
        Args:
            query: 查询文本（自然语言）
            top_k: 返回前 k 个最相关的工具（默认 5）
            category_filter: 可选的类别过滤器（如 "Metabolomics"）
        
        Returns:
            工具 schema 列表，每个元素包含完整的工具信息（name, description, args_schema 等）
        """
        try:
            # 构建搜索查询
            search_kwargs = {}
            
            # 如果指定了类别过滤，添加到 metadata 过滤中
            if category_filter:
                # ChromaDB 支持 metadata 过滤
                search_kwargs["filter"] = {"category": category_filter}
            
            # 执行相似度搜索（k 参数直接传递，不放在 search_kwargs 中）
            results = self.vector_store.similarity_search_with_score(
                query,
                k=top_k,
                **search_kwargs
            )
            
            # 转换为工具 schema 格式
            tools = []
            for doc, score in results:
                metadata = doc.metadata
                
                # 解析 args_schema（从 JSON 字符串）
                try:
                    args_schema = json.loads(metadata.get("args_schema", "{}"))
                except Exception:
                    args_schema = {}
                
                tool_schema = {
                    "name": metadata.get("name"),
                    "description": metadata.get("description"),
                    "category": metadata.get("category"),
                    "output_type": metadata.get("output_type"),
                    "args_schema": args_schema,
                    "similarity_score": float(score)  # 相似度分数（越低越相似）
                }
                
                tools.append(tool_schema)
            
            logger.info(f"🔍 检索到 {len(tools)} 个相关工具 (查询: '{query}')")
            return tools
        
        except Exception as e:
            logger.error(f"❌ 工具检索失败: {e}", exc_info=True)
            return []
    
    def get_tool_by_name(self, name: str) -> Optional[Dict[str, Any]]:
        """
        根据工具名称获取工具 schema
        
        Args:
            name: 工具名称
        
        Returns:
            工具 schema，如果不存在返回 None
        """
        try:
            # 使用精确匹配搜索
            results = self.vector_store.similarity_search_with_score(
                f"工具名称: {name}",
                k=1
            )
            
            if results:
                doc, score = results[0]
                metadata = doc.metadata
                
                if metadata.get("name") == name:
                    try:
                        args_schema = json.loads(metadata.get("args_schema", "{}"))
                    except Exception:
                        args_schema = {}
                    
                    return {
                        "name": metadata.get("name"),
                        "description": metadata.get("description"),
                        "category": metadata.get("category"),
                        "output_type": metadata.get("output_type"),
                        "args_schema": args_schema
                    }
            
            return None
        
        except Exception as e:
            logger.error(f"❌ 获取工具失败: {e}", exc_info=True)
            return None
    
    def list_all_tools(self) -> List[str]:
        """
        列出所有已注册的工具名称
        
        Returns:
            工具名称列表
        """
        return registry.list_tools()

