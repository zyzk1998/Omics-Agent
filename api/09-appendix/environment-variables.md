# 环境变量配置

| 变量名 | 说明 | 默认值 |
|--------|------|--------|
| `UPLOAD_DIR` | 上传文件目录 | `/app/uploads` |
| `RESULTS_DIR` | 结果输出目录 | `/app/results` |
| `MAX_FILE_SIZE` | 最大文件大小（字节） | `104857600` (100MB) |
| `ALLOWED_ORIGINS` | CORS 允许的来源 | `*` |
| `SILICONFLOW_API_KEY` | SiliconFlow API Key | - |
| `SILICONFLOW_MODEL` | SiliconFlow 模型名称 | - |
| `OLLAMA_BASE_URL` | Ollama 服务地址 | `http://localhost:11434` |
| `OLLAMA_EMBEDDING_MODEL` | Ollama Embedding 模型 | `nomic-embed-text` |
| `CHROMA_PERSIST_DIR` | ChromaDB 持久化目录 | `./data/chroma_tools` |

---

**返回**: [附录目录](../README.md) | [API 手册首页](../../README.md)
