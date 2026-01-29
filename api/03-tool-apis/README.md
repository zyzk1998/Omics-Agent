# 工具 API

工具 API 提供工具检索和管理功能，支持语义搜索、列出工具、获取工具 Schema 等。

## 接口列表

1. [语义搜索工具](search.md) - `GET /api/tools/search`
2. [列出所有工具](list.md) - `GET /api/tools/list`
3. [获取工具 Schema](get-schema.md) - `GET /api/tools/{tool_name}`

---

## 使用场景

- **语义搜索**: 根据自然语言查询，找到相关的分析工具
- **工具发现**: 列出所有可用的工具
- **工具详情**: 获取工具的完整 Schema（参数定义、描述等）

---

**返回**: [API 手册首页](../README.md)
