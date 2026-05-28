#!/usr/bin/env bash
# 将仓库内 prompt_specs 同步到运行时 SKILLS_ASSETS_ROOT/prompt/（供运维/容器挂载后读取）
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SRC="${ROOT}/gibh_agent/skills/prompt_specs"
DEST="${SKILLS_ASSETS_ROOT:-${UPLOAD_DIR:-${ROOT}/data/uploads}/skills_assets}/prompt"
mkdir -p "${DEST}"
for d in "${SRC}"/*/; do
  [[ -d "${d}" ]] || continue
  id="$(basename "${d}")"
  if [[ -f "${d}/SKILL.md" ]]; then
    mkdir -p "${DEST}/${id}"
    cp -f "${d}/SKILL.md" "${DEST}/${id}/SKILL.md"
    echo "synced ${id} -> ${DEST}/${id}/SKILL.md"
  fi
done
# 文档镜像
DOCS="${ROOT}/docs/skills/prompt"
mkdir -p "${DOCS}"
for d in "${SRC}"/*/; do
  id="$(basename "${d}")"
  [[ -f "${d}/SKILL.md" ]] && cp -f "${d}/SKILL.md" "${DOCS}/${id}.md"
done
echo "Done. Runtime prompt dir: ${DEST}"
