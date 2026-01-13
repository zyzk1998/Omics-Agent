#!/usr/bin/env python3
"""
测试工具检索器

测试脚本：验证 Tool Retriever 的语义搜索功能
"""
import sys
import os
from pathlib import Path

# 添加项目路径
sys.path.insert(0, str(Path(__file__).parent))

def test_tool_retrieval():
    """测试工具检索功能"""
    print("🔍 测试工具检索器")
    print("=" * 60)
    
    try:
        # 1. 导入工具定义（触发注册）
        print("\n1️⃣ 导入工具定义...")
        from gibh_agent.tools.definitions import *
        from gibh_agent.core.tool_registry import registry
        
        tools = registry.list_tools()
        print(f"   ✅ 已注册 {len(tools)} 个工具: {tools}")
        
        # 2. 初始化工具检索器
        print("\n2️⃣ 初始化工具检索器...")
        from gibh_agent.core.tool_retriever import ToolRetriever
        
        chroma_dir = "./data/chroma_tools_test"
        embedding_model = os.getenv("OLLAMA_EMBEDDING_MODEL", "nomic-embed-text")
        ollama_url = os.getenv("OLLAMA_BASE_URL", "http://localhost:11434")
        
        print(f"   ChromaDB 目录: {chroma_dir}")
        print(f"   Embedding 模型: {embedding_model}")
        print(f"   Ollama URL: {ollama_url}")
        
        retriever = ToolRetriever(
            persist_directory=chroma_dir,
            embedding_model=embedding_model,
            ollama_base_url=ollama_url
        )
        print("   ✅ 工具检索器初始化成功")
        
        # 3. 同步工具到 ChromaDB
        print("\n3️⃣ 同步工具到 ChromaDB...")
        synced_count = retriever.sync_tools(clear_existing=True)
        print(f"   ✅ 成功同步 {synced_count} 个工具")
        
        # 4. 测试语义搜索
        print("\n4️⃣ 测试语义搜索...")
        test_queries = [
            "I want to reduce dimensions",
            "perform differential analysis",
            "preprocess metabolite data",
            "inspect a file"
        ]
        
        for query in test_queries:
            print(f"\n   🔍 查询: '{query}'")
            results = retriever.retrieve(query=query, top_k=3)
            
            if results:
                print(f"   ✅ 找到 {len(results)} 个相关工具:")
                for i, tool in enumerate(results, 1):
                    print(f"      {i}. {tool['name']} (相似度: {tool['similarity_score']:.4f})")
                    print(f"         描述: {tool['description'][:60]}...")
            else:
                print("   ⚠️ 未找到相关工具")
        
        # 5. 测试按名称获取工具
        print("\n5️⃣ 测试按名称获取工具...")
        tool_name = "metabolomics_pca"
        tool_schema = retriever.get_tool_by_name(tool_name)
        
        if tool_schema:
            print(f"   ✅ 找到工具: {tool_schema['name']}")
            print(f"      类别: {tool_schema['category']}")
            print(f"      参数: {list(tool_schema['args_schema'].get('properties', {}).keys())}")
        else:
            print(f"   ❌ 未找到工具: {tool_name}")
        
        print("\n" + "=" * 60)
        print("✅ 测试完成！")
        
    except ImportError as e:
        print(f"\n❌ 导入失败: {e}")
        print("\n请安装依赖:")
        print("  pip install langchain-chroma langchain-ollama langchain-core")
        print("\n确保 Ollama 服务正在运行:")
        print("  ollama serve")
        print("  ollama pull nomic-embed-text")
        sys.exit(1)
    except Exception as e:
        print(f"\n❌ 测试失败: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    test_tool_retrieval()

