# 获取工具 Schema

## `GET /api/tools/{tool_name}`

**说明**: 获取特定工具的完整 Schema

**路径参数**:
- `tool_name` (string, 必需): 工具名称

**响应示例**:

```json
{
  "status": "success",
  "tool": {
    "name": "preprocess_data",
    "description": "数据预处理工具",
    "category": "Metabolomics",
    "parameters": {
      "file_path": {
        "type": "string",
        "description": "文件路径",
        "required": true
      },
      "missing_threshold": {
        "type": "string",
        "description": "缺失值阈值",
        "default": "0.5"
      }
    }
  }
}
```

**错误响应** (404 Not Found):

```json
{
  "detail": "工具 'unknown_tool' 不存在"
}
```

---

**返回**: [工具 API 目录](../README.md) | [API 手册首页](../../README.md)
