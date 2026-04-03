#!/usr/bin/env bash
# 兼容入口：转发至 scripts/清理后台结果数据.sh（请勿删除，旧文档/CI 可能仍引用本路径）。
exec "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/清理后台结果数据.sh" "$@"
