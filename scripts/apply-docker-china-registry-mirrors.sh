#!/usr/bin/env bash
# 一键合并 Docker Hub 国内 registry-mirrors（缓解 FROM docker.io/... 极慢）。
# 用法：sudo bash scripts/apply-docker-china-registry-mirrors.sh
# 然后：sudo systemctl restart docker   （脚本末尾会尝试执行）
set -euo pipefail

if [[ "${EUID:-$(id -u)}" -ne 0 ]]; then
  echo "请使用 root 执行：sudo bash $0" >&2
  exit 1
fi

TARGET="/etc/docker/daemon.json"
BACKUP="/etc/docker/daemon.json.bak.gibh-$(date +%Y%m%d%H%M%S)"
# 多镜像：前者失败时可由 Docker 尝试列表中其它项（视引擎版本行为而定）
MIRRORS='["https://docker.m.daocloud.io","https://docker.1ms.run"]'

python3 - "$TARGET" "$BACKUP" "$MIRRORS" <<'PY'
import json
import sys
from pathlib import Path

target = Path(sys.argv[1])
backup = Path(sys.argv[2])
extra = json.loads(sys.argv[3])

if target.exists():
    backup.write_text(target.read_text(encoding="utf-8"), encoding="utf-8")
    data = json.loads(target.read_text(encoding="utf-8"))
else:
    data = {}

cur = list(data.get("registry-mirrors") or [])
seen = set(cur)
for u in extra:
    if u not in seen:
        cur.append(u)
        seen.add(u)
data["registry-mirrors"] = cur
# 多 blob 并行拉取（默认常为 3）；对大镜像可减轻「首层飞快、后续层像卡住」体感
if "max-concurrent-downloads" not in data:
    data["max-concurrent-downloads"] = 10
if "max-concurrent-uploads" not in data:
    data["max-concurrent-uploads"] = 5
target.parent.mkdir(parents=True, exist_ok=True)
target.write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
print("已写入:", target)
if backup.exists():
    print("已备份:", backup)
print("registry-mirrors:", cur)
print("max-concurrent-downloads:", data.get("max-concurrent-downloads"))
PY

if systemctl is-active --quiet docker 2>/dev/null; then
  systemctl restart docker
  echo "已执行: systemctl restart docker"
else
  echo "未检测到 systemd active 的 docker 服务，请手动重启 Docker 后再 build。" >&2
fi

echo "下一步: cd 到项目根目录后执行 sudo docker compose build worker-pyskills"
