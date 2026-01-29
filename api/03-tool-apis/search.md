# 语义搜索工具

## `GET /api/tools/search`

**说明**: 语义搜索工具（基于 ChromaDB + Embeddings）

**请求参数**:
- `query` (string, 必需): 查询文本（自然语言）
- `top_k` (int, 可选): 返回前 k 个最相关的工具，默认 5
- `category` (string, 可选): 类别过滤器（如 "Metabolomics", "scRNA-seq"）

**响应示例**:

```json
{
  "status": "success",
  "query": "数据预处理",
  "count": 3,
  "tools": [
    {
      "name": "preprocess_data",
      "description": "数据预处理工具",
      "category": "Metabolomics",
      "parameters": { ... }
    }
  ]
}
```

**错误响应** (503 Service Unavailable):

```json
{
  "detail": "工具检索器未初始化。请检查 Ollama 服务和依赖是否已安装。"
}
```

---

**返回**: [工具 API 目录](../README.md) | [API 手册首页](../../README.md)
