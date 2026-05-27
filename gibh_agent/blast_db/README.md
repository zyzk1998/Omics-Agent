# 本地 BLAST 核酸库挂载目录

将 NCBI `blastdb` 下载的 **nt**（或 **core_nt**）库文件置于此目录，例如：

- `nt.nhr`, `nt.nin`, `nt.nsq`（及可选 `.nog` 等）

或通过环境变量指定其它路径：

```bash
export LOCAL_BLASTDB_PATH=/data/blastdb
```

`launch-skills` / `api-server` 检测到本地库后，核酸比对将**自动走本地引擎**（无 `-remote`），响应通常为秒级至数分钟级，且耗时稳定。

未挂载时回退 **NCBI 远程**，可能排队 1–10 分钟，属正常现象。
