# 列出所有工具

## `GET /api/tools/list`

**说明**: 列出所有已注册的工具

**请求参数**: 无

**响应示例**:

```json
{
  "status": "success",
  "count": 15,
  "tools": [
    "inspect_data",
    "preprocess_data",
    "pca_analysis",
    ...
  ]
}
```

---

**返回**: [工具 API 目录](../README.md) | [API 手册首页](../../README.md)
