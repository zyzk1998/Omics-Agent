# 文件路径说明

## 上传文件路径

- **格式**: `{user_id}/{session_id}/{filename}`
- **相对路径**: 相对于 `UPLOAD_DIR`（默认 `/app/uploads`）
- **示例**: `guest/20250128_120000/example.csv`

## 结果文件路径

- **格式**: `run_{timestamp}/{filename}`
- **相对路径**: 相对于 `RESULTS_DIR`（默认 `/app/results`）
- **示例**: `run_20250128_120000/pca_plot.png`

## 访问结果文件

通过静态文件服务访问：
- **URL 格式**: `/results/{path}`
- **示例**: `/results/run_20250128_120000/pca_plot.png`

---

**返回**: [附录目录](../README.md) | [API 手册首页](../../README.md)
