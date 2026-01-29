# 10x Genomics 数据特殊处理

## 自动识别

系统自动识别 10x Genomics 数据文件：
- `matrix.mtx` - 表达矩阵
- `barcodes.tsv` - 细胞条形码
- `features.tsv` 或 `genes.tsv` - 基因/特征信息

## 自动分组

当检测到 10x 数据文件时：
- 自动创建子目录: `10x_data_YYYYMMDD_HHMMSS`
- 将所有 10x 文件保存到该子目录
- 返回的 `file_paths` 指向该子目录，而不是单个文件

## 响应格式

```json
{
  "status": "success",
  "is_10x_data": true,
  "group_dir": "guest/20250128_120000/10x_data_20250128_120000",
  "file_paths": ["guest/20250128_120000/10x_data_20250128_120000"],
  ...
}
```

---

**返回**: [附录目录](../README.md) | [API 手册首页](../../README.md)
