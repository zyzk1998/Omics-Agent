#!/bin/bash
# 环境变量和挂载点验证脚本

echo "📋 环境变量和挂载点检查"
echo ""

# 检查 API 容器
if docker ps | grep -q "gibh_v2_api"; then
    echo "=== API 容器 (gibh_v2_api) ==="
    echo ""
    echo "环境变量:"
    docker exec gibh_v2_api env | grep -E "UPLOAD_DIR|RESULTS_DIR" || echo "  ⚠️ 未找到环境变量"
    echo ""
    echo "挂载点检查:"
    docker exec gibh_v2_api ls -ld /app/uploads 2>/dev/null && echo "  ✅ /app/uploads 存在" || echo "  ❌ /app/uploads 不存在"
    docker exec gibh_v2_api test -w /app/uploads 2>/dev/null && echo "  ✅ /app/uploads 可写" || echo "  ❌ /app/uploads 不可写"
    docker exec gibh_v2_api ls -ld /app/results 2>/dev/null && echo "  ✅ /app/results 存在" || echo "  ❌ /app/results 不存在"
    docker exec gibh_v2_api test -w /app/results 2>/dev/null && echo "  ✅ /app/results 可写" || echo "  ❌ /app/results 不可写"
    echo ""
else
    echo "⚠️ API 容器 (gibh_v2_api) 未运行"
    echo ""
fi

# 检查 Worker 容器
if docker ps | grep -q "gibh_v2_worker"; then
    echo "=== Worker 容器 (gibh_v2_worker) ==="
    echo ""
    echo "环境变量:"
    docker exec gibh_v2_worker env | grep -E "UPLOAD_DIR|RESULTS_DIR" || echo "  ⚠️ 未找到环境变量"
    echo ""
    echo "挂载点检查:"
    docker exec gibh_v2_worker ls -ld /app/uploads 2>/dev/null && echo "  ✅ /app/uploads 存在" || echo "  ❌ /app/uploads 不存在"
    docker exec gibh_v2_worker test -w /app/uploads 2>/dev/null && echo "  ✅ /app/uploads 可写" || echo "  ❌ /app/uploads 不可写"
    docker exec gibh_v2_worker ls -ld /app/results 2>/dev/null && echo "  ✅ /app/results 存在" || echo "  ❌ /app/results 不存在"
    docker exec gibh_v2_worker test -w /app/results 2>/dev/null && echo "  ✅ /app/results 可写" || echo "  ❌ /app/results 不可写"
    echo ""
else
    echo "⚠️ Worker 容器 (gibh_v2_worker) 未运行"
    echo ""
fi

# 检查主机挂载目录
echo "=== 主机挂载目录 ==="
echo ""
if [ -d "./data/uploads" ]; then
    echo "  ✅ ./data/uploads 存在"
    ls -ld ./data/uploads | awk '{print "     权限: " $1 " 所有者: " $3 ":" $4}'
else
    echo "  ❌ ./data/uploads 不存在"
fi

if [ -d "./results" ]; then
    echo "  ✅ ./results 存在"
    ls -ld ./results | awk '{print "     权限: " $1 " 所有者: " $3 ":" $4}'
else
    echo "  ❌ ./results 不存在"
fi

echo ""
echo "=== Docker Compose 配置 ==="
echo ""
echo "环境变量配置:"
grep -A 1 "UPLOAD_DIR" docker-compose.yml | head -2
echo ""
echo "挂载点配置:"
grep -A 1 "volumes:" docker-compose.yml | grep "uploads" | head -1

echo ""
echo "✅ 检查完成"

