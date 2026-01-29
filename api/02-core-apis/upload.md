# 文件上传接口

## `POST /api/upload`

**说明**: 上传一个或多个文件，支持 10x Genomics 数据（自动识别并分组）

**请求格式**: `multipart/form-data`

**请求参数**:
- `files` (File[], 必需): 文件列表（支持多文件上传，最多 20 个）
- `user_id` (string, 可选): 用户ID，默认 "guest"
- `session_id` (string, 可选): 会话ID，未提供时自动生成（格式: `YYYYMMDD_HHMMSS`）

**支持的文件类型**:
- `.h5ad` - AnnData 格式（单细胞数据）
- `.mtx` - Matrix Market 格式
- `.tsv`, `.csv` - 表格数据
- `.txt` - 文本文件
- `.gz`, `.tar`, `.zip` - 压缩文件

**文件大小限制**: 默认 100MB（可通过环境变量 `MAX_FILE_SIZE` 配置）

**10x Genomics 数据自动识别**:
- 如果上传的文件包含 `matrix.mtx`、`barcodes.tsv`、`features.tsv`（或 `genes.tsv`），系统会自动识别为 10x Genomics 数据
- 10x 数据会被保存到独立的子目录中（格式: `10x_data_YYYYMMDD_HHMMSS`）
- 返回的 `file_paths` 将指向该子目录，而不是单个文件

**成功响应** (200 OK):

```json
{
  "status": "success",
  "file_paths": [
    "guest/20250128_120000/example.csv",
    "guest/20250128_120000/10x_data_20250128_120000"
  ],
  "file_info": [
    {
      "name": "example.csv",
      "size": 1024000,
      "path": "guest/20250128_120000/example.csv"
    }
  ],
  "count": 2,
  "user_id": "guest",
  "session_id": "20250128_120000",
  "is_10x_data": true,
  "group_dir": "guest/20250128_120000/10x_data_20250128_120000",
  "files": [...]
}
```

**响应字段说明**:
- `status`: 操作状态（"success" 表示成功）
- `file_paths`: 文件路径数组（相对路径，相对于 `/app/uploads`）
- `file_info`: 文件信息数组，包含 `name`、`size`、`path`
- `count`: 上传的文件数量
- `user_id`: 用户ID
- `session_id`: 会话ID
- `is_10x_data`: 是否为 10x Genomics 数据（仅当检测到 10x 数据时存在）
- `group_dir`: 10x 数据组目录路径（仅当检测到 10x 数据时存在）

**错误响应**:
- `400 Bad Request`: 请求参数错误
- `403 Forbidden`: 文件路径不安全
- `413 Payload Too Large`: 文件大小超限

**前端集成示例**:

```javascript
const formData = new FormData();
for (let file of fileInput.files) {
    formData.append('files', file);
}
formData.append('user_id', 'guest');
formData.append('session_id', '20250128_120000');

const response = await fetch('/api/upload', {
    method: 'POST',
    body: formData
});

const result = await response.json();
if (result.status === 'success') {
    console.log('上传成功:', result.file_paths);
    uploadedFiles = result.file_paths;
}
```

---

**返回**: [核心 API 目录](../README.md) | [API 手册首页](../../README.md)
