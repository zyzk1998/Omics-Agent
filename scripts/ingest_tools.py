#!/usr/bin/env python3
"""
Ingest all registered tools (including 7 omics atomic tools) into ChromaDB.
Run from project root: python scripts/ingest_tools.py
or: PYTHONPATH=. python scripts/ingest_tools.py
"""
import os
import sys
import logging
from pathlib import Path

# Ensure project root is on path
project_root = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(project_root))
os.chdir(project_root)

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


def main():
    # 1. Load all tools (triggers @registry.register for every module)
    from gibh_agent.tools import load_all_tools
    from gibh_agent.core.tool_registry import registry
    from gibh_agent.core.tool_retriever import ToolRetriever

    logger.info("Loading all tools (including atomic omics modules)...")
    load_all_tools()
    n_registered = len(registry._tools)
    logger.info("Total tools registered: %s", n_registered)

    # 2. Initialize retriever and sync to ChromaDB
    chroma_dir = os.getenv("CHROMA_PERSIST_DIR", "./data/chroma_tools")
    embedding_model = os.getenv("OLLAMA_EMBEDDING_MODEL", "nomic-embed-text")
    ollama_url = os.getenv("OLLAMA_BASE_URL", "http://localhost:11434")

    try:
        retriever = ToolRetriever(
            persist_directory=chroma_dir,
            embedding_model=embedding_model,
            ollama_base_url=ollama_url,
        )
    except Exception as e:
        logger.error("ToolRetriever init failed (need langchain-chroma, langchain-ollama, Ollama): %s", e)
        return 1

    n_synced = retriever.sync_tools(clear_existing=True)
    logger.info("Tools ingested into ChromaDB: %s", n_synced)
    print("\n--- Summary ---")
    print(f"Tools registered: {n_registered}")
    print(f"Tools ingested (ChromaDB): {n_synced}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
