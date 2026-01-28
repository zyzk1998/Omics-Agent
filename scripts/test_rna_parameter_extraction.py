#!/usr/bin/env python3
"""
测试RNA分析参数提取逻辑
模拟从steps_details中提取steps_results的过程
"""

import json
import sys
import os

# 添加项目根目录到路径
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# 模拟从日志中提取的steps_details数据
mock_steps_details = [
    {
        "step_id": "rna_qc_filter",
        "tool_id": "rna_qc_filter",
        "name": "质量控制过滤",
        "status": "success",
        "summary": "步骤 质量控制过滤 执行完成",
        "step_result": {
            "step_name": "质量控制过滤",
            "status": "success",
            "logs": "步骤 质量控制过滤 执行完成",
            "data": {
                "status": "success",
                "n_obs_before": 2866,
                "n_obs_after": 2528,
                "n_vars_before": 33538,
                "n_vars_after": 20728,
                "plot_path": "results/run_20260128_163805/qc_violin_1769589485.png",
                "output_h5ad": "results/run_20260128_163805/filtered.h5ad",
                "summary": "过滤后剩余 2528 个细胞，20728 个基因"
            }
        }
    },
    {
        "step_id": "rna_pca",
        "tool_id": "rna_pca",
        "name": "主成分分析 (PCA)",
        "status": "success",
        "summary": "步骤 主成分分析 (PCA) 执行完成",
        "step_result": {
            "step_name": "主成分分析 (PCA)",
            "status": "success",
            "logs": "步骤 主成分分析 (PCA) 执行完成",
            "data": {
                "status": "success",
                "n_comps": 50,
                "explained_variance": {
                    "PC1": 0.0739179328083992,
                    "PC2": 0.03615681827068329
                },
                "plot_path": "results/run_20260128_163805/pca_variance_1769589488.png",
                "output_h5ad": "results/run_20260128_163805/pca.h5ad",
                "summary": "PCA 降维完成"
            }
        }
    },
    {
        "step_id": "rna_clustering",
        "tool_id": "rna_clustering",
        "name": "Leiden 聚类",
        "status": "success",
        "summary": "步骤 Leiden 聚类 执行完成",
        "step_result": {
            "step_name": "Leiden 聚类",
            "status": "success",
            "logs": "步骤 Leiden 聚类 执行完成",
            "data": {
                "status": "success",
                "algorithm": "leiden",
                "resolution": 0.5,
                "n_clusters": 13,
                "cluster_key": "leiden",
                "output_h5ad": "results/run_20260128_163805/leiden_clustered.h5ad",
                "summary": "Leiden 聚类 (Res=0.5): 13 个簇"
            }
        }
    },
    {
        "step_id": "rna_find_markers",
        "tool_id": "rna_find_markers",
        "name": "Marker 基因检测",
        "status": "success",
        "summary": "步骤 Marker 基因检测 执行完成",
        "step_result": {
            "step_name": "Marker 基因检测",
            "status": "success",
            "logs": "步骤 Marker 基因检测 执行完成",
            "data": {
                "status": "success",
                "method": "t-test",
                "n_clusters": 13,
                "n_genes_per_cluster": 5,
                "markers_table": [
                    {"0_names": "S100A8", "0_pvals": 0.0},
                    {"1_names": "SRGN", "1_pvals": 4.84041e-318},
                    {"2_names": "NKG7", "2_pvals": 1.0933693412268645e-168}
                ],
                "output_csv": "results/run_20260128_163805/markers.csv",
                "summary": "Marker 基因鉴定完成"
            }
        }
    }
]

def test_parameter_extraction():
    """测试参数提取逻辑"""
    print("=" * 80)
    print("🧪 测试RNA分析参数提取逻辑")
    print("=" * 80)
    
    # 模拟orchestrator中的提取逻辑
    steps_results = []
    for step_detail in mock_steps_details:
        if "step_result" in step_detail:
            steps_results.append(step_detail["step_result"])
        elif "status" in step_detail:
            # 如果没有step_result，构建一个基本的step_result
            steps_results.append({
                "step_name": step_detail.get("name", step_detail.get("step_id", "Unknown")),
                "status": step_detail.get("status", "unknown"),
                "data": step_detail.get("data", {})
            })
    
    print(f"\n✅ 提取到 {len(steps_results)} 个步骤结果")
    print(f"\n步骤详情:")
    for i, step_result in enumerate(steps_results, 1):
        step_name = step_result.get("step_name", "Unknown")
        status = step_result.get("status", "unknown")
        data = step_result.get("data", {})
        print(f"\n{i}. {step_name} (状态: {status})")
        print(f"   - data keys: {list(data.keys())}")
        
        # 检查是否有summary
        if "summary" in data:
            print(f"   - summary: {data['summary']}")
        
        # 检查RNA特定字段
        if "n_obs_after" in data:
            print(f"   - n_obs_after: {data['n_obs_after']}")
        if "n_vars_after" in data:
            print(f"   - n_vars_after: {data['n_vars_after']}")
        if "n_clusters" in data:
            print(f"   - n_clusters: {data['n_clusters']}")
        if "explained_variance" in data:
            print(f"   - explained_variance: PC1={data['explained_variance'].get('PC1', 'N/A')}, PC2={data['explained_variance'].get('PC2', 'N/A')}")
        if "markers_table" in data:
            print(f"   - markers_table: {len(data['markers_table'])} 行")
    
    # 检查提取的数据结构是否符合_generate_analysis_summary的期望
    print("\n" + "=" * 80)
    print("📊 数据结构检查")
    print("=" * 80)
    
    required_fields = ["step_name", "status", "data"]
    for i, step_result in enumerate(steps_results, 1):
        missing_fields = [field for field in required_fields if field not in step_result]
        if missing_fields:
            print(f"❌ 步骤 {i} 缺少字段: {missing_fields}")
        else:
            print(f"✅ 步骤 {i} 数据结构完整")
    
    # 检查data中是否有summary字段（用于指标提取）
    print("\n" + "=" * 80)
    print("📋 Summary字段检查")
    print("=" * 80)
    
    for i, step_result in enumerate(steps_results, 1):
        data = step_result.get("data", {})
        if "summary" in data:
            print(f"✅ 步骤 {i} ({step_result.get('step_name')}) 有summary字段")
        else:
            print(f"⚠️ 步骤 {i} ({step_result.get('step_name')}) 没有summary字段，但可能有其他字段")
            print(f"   - data keys: {list(data.keys())}")
    
    return steps_results

if __name__ == "__main__":
    steps_results = test_parameter_extraction()
    print("\n" + "=" * 80)
    print("✅ 测试完成")
    print("=" * 80)
    print(f"\n提取的steps_results数量: {len(steps_results)}")
    print(f"JSON长度: {len(json.dumps(steps_results, ensure_ascii=False))} 字符")
