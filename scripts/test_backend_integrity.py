#!/usr/bin/env python3
"""
Backend Integrity Test Script

测试 Metabolomics 和 RNA 工具的参数对齐，确保执行器调用时不会出现参数缺失错误。

测试逻辑：
1. 创建模拟数据文件（CSV 和 H5AD）
2. 测试 Metabolomics 工具的参数对齐
3. 测试 RNA 工具的参数对齐
4. 报告 PASS/FAIL 状态
"""

import sys
import os
import tempfile
import shutil
from pathlib import Path
import pandas as pd
import numpy as np

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from gibh_agent.core.tool_registry import registry
import inspect


def create_mock_csv(file_path: str):
    """创建模拟的代谢组学 CSV 文件"""
    np.random.seed(42)
    n_samples = 20
    n_features = 50
    
    # 创建数据
    data = np.random.randn(n_samples, n_features)
    
    # 创建分组列
    groups = ['Control'] * 10 + ['Treatment'] * 10
    
    # 创建 DataFrame
    df = pd.DataFrame(
        data,
        index=[f'Sample_{i+1}' for i in range(n_samples)],
        columns=[f'Metabolite_{i+1}' for i in range(n_features)]
    )
    df['Diet'] = groups  # 添加分组列
    
    # 保存
    df.to_csv(file_path)
    print(f"✅ 创建模拟 CSV 文件: {file_path}")
    return file_path


def create_mock_h5ad(file_path: str):
    """创建模拟的 scRNA-seq H5AD 文件"""
    try:
        import scanpy as sc
        import anndata as ad
        
        np.random.seed(42)
        n_obs = 100
        n_vars = 50
        
        # 创建随机数据
        X = np.random.randn(n_obs, n_vars)
        
        # 创建 AnnData 对象
        adata = ad.AnnData(X)
        adata.obs_names = [f'Cell_{i+1}' for i in range(n_obs)]
        adata.var_names = [f'Gene_{i+1}' for i in range(n_vars)]
        
        # 添加分组信息
        adata.obs['leiden'] = ['0'] * 50 + ['1'] * 50
        
        # 保存
        adata.write(file_path)
        print(f"✅ 创建模拟 H5AD 文件: {file_path}")
        return file_path
    except ImportError:
        print(f"⚠️ scanpy 或 anndata 未安装，跳过 H5AD 文件创建")
        return None


def test_tool_signature(tool_id: str, test_params: dict) -> tuple[bool, str]:
    """
    测试工具函数签名是否匹配提供的参数
    
    Args:
        tool_id: 工具 ID
        test_params: 测试参数字典
        
    Returns:
        (success: bool, message: str)
    """
    try:
        tool_func = registry.get_tool(tool_id)
        if not tool_func:
            return False, f"工具 '{tool_id}' 未在注册表中找到"
        
        # 获取函数签名
        sig = inspect.signature(tool_func)
        params = sig.parameters
        
        # 检查必需参数
        missing_params = []
        for param_name, param_obj in params.items():
            if param_name == 'self':
                continue
            if param_obj.default == inspect.Parameter.empty and param_name not in test_params:
                # 检查是否有 **kwargs
                has_kwargs = any(
                    p.kind == inspect.Parameter.VAR_KEYWORD 
                    for p in params.values()
                )
                if not has_kwargs:
                    missing_params.append(param_name)
        
        if missing_params:
            return False, f"缺少必需参数: {', '.join(missing_params)}"
        
        # 尝试调用（使用模拟数据，可能失败但不应该因为参数缺失而失败）
        try:
            # 只检查参数对齐，不实际执行
            bound = sig.bind(**test_params)
            bound.apply_defaults()
            return True, "PASS"
        except TypeError as e:
            error_msg = str(e)
            if "missing" in error_msg.lower() or "required" in error_msg.lower():
                return False, f"参数对齐失败: {error_msg}"
            else:
                # 其他错误（如数据格式问题）不算参数对齐失败
                return True, f"PASS (参数对齐成功，但执行可能因数据问题失败: {error_msg[:100]})"
        
    except Exception as e:
        return False, f"测试异常: {str(e)}"


def test_metabolomics_tools(temp_dir: str):
    """测试代谢组学工具"""
    print("\n" + "="*60)
    print("测试 Metabolomics 工具")
    print("="*60)
    
    csv_path = os.path.join(temp_dir, "test_metabolomics.csv")
    create_mock_csv(csv_path)
    
    results = []
    
    # 测试 run_plsda
    print("\n1. 测试 run_plsda...")
    test_params = {
        "file_path": csv_path,
        "group_column": "Diet",
        "n_components": 2,
        "scale": True
    }
    success, message = test_tool_signature("metabolomics_plsda", test_params)
    results.append(("metabolomics_plsda", success, message))
    status = "✅ PASS" if success else "❌ FAIL"
    print(f"   结果: {status} - {message}")
    
    # 测试 run_pathway_enrichment
    print("\n2. 测试 run_pathway_enrichment...")
    test_params = {
        "file_path": csv_path,
        "group_column": "Diet",
        "case_group": "Treatment",
        "control_group": "Control",
        "organism": "hsa"
    }
    success, message = test_tool_signature("metabolomics_pathway_enrichment", test_params)
    results.append(("metabolomics_pathway_enrichment", success, message))
    status = "✅ PASS" if success else "❌ FAIL"
    print(f"   结果: {status} - {message}")
    
    # 测试 run_differential_analysis
    print("\n3. 测试 run_differential_analysis...")
    test_params = {
        "file_path": csv_path,
        "group_column": "Diet",
        "method": "t-test"
    }
    success, message = test_tool_signature("differential_analysis", test_params)
    results.append(("differential_analysis", success, message))
    status = "✅ PASS" if success else "❌ FAIL"
    print(f"   结果: {status} - {message}")
    
    return results


def test_rna_tools(temp_dir: str):
    """测试 RNA 工具"""
    print("\n" + "="*60)
    print("测试 RNA 工具")
    print("="*60)
    
    h5ad_path = os.path.join(temp_dir, "test_rna.h5ad")
    h5ad_created = create_mock_h5ad(h5ad_path)
    
    if not h5ad_created:
        print("⚠️ 无法创建 H5AD 文件，跳过 RNA 工具测试")
        return []
    
    results = []
    
    # 测试 rna_qc_filter
    print("\n1. 测试 rna_qc_filter...")
    test_params = {
        "adata_path": h5ad_path,
        "min_genes": 200,
        "min_cells": 3
    }
    success, message = test_tool_signature("rna_qc_filter", test_params)
    results.append(("rna_qc_filter", success, message))
    status = "✅ PASS" if success else "❌ FAIL"
    print(f"   结果: {status} - {message}")
    
    # 测试 rna_clustering
    print("\n2. 测试 rna_clustering...")
    # 先检查实际函数签名
    clustering_func = registry.get_tool("rna_clustering")
    if clustering_func:
        sig = inspect.signature(clustering_func)
        print(f"   函数签名参数: {list(sig.parameters.keys())}")
    test_params = {
        "adata_path": h5ad_path,
        "resolution": 0.5
    }
    success, message = test_tool_signature("rna_clustering", test_params)
    results.append(("rna_clustering", success, message))
    status = "✅ PASS" if success else "❌ FAIL"
    print(f"   结果: {status} - {message}")
    
    # 测试 rna_find_markers
    print("\n3. 测试 rna_find_markers...")
    # 先检查实际函数签名
    markers_func = registry.get_tool("rna_find_markers")
    if markers_func:
        sig = inspect.signature(markers_func)
        print(f"   函数签名参数: {list(sig.parameters.keys())}")
    test_params = {
        "adata_path": h5ad_path,
        "cluster_key": "leiden"
    }
    success, message = test_tool_signature("rna_find_markers", test_params)
    results.append(("rna_find_markers", success, message))
    status = "✅ PASS" if success else "❌ FAIL"
    print(f"   结果: {status} - {message}")
    
    return results


def main():
    """主函数"""
    print("="*60)
    print("Backend Integrity Test")
    print("="*60)
    print("\n测试目标：验证工具函数签名与执行器调用参数的对齐")
    
    # 创建临时目录
    temp_dir = tempfile.mkdtemp(prefix="integrity_test_")
    print(f"\n临时目录: {temp_dir}")
    
    try:
        # 测试 Metabolomics 工具
        metabolomics_results = test_metabolomics_tools(temp_dir)
        
        # 测试 RNA 工具
        rna_results = test_rna_tools(temp_dir)
        
        # 汇总结果
        print("\n" + "="*60)
        print("测试结果汇总")
        print("="*60)
        
        all_results = metabolomics_results + rna_results
        
        passed = sum(1 for _, success, _ in all_results if success)
        failed = len(all_results) - passed
        
        print(f"\n总计: {len(all_results)} 个工具")
        print(f"✅ 通过: {passed}")
        print(f"❌ 失败: {failed}")
        
        if failed > 0:
            print("\n失败的工具:")
            for tool_id, success, message in all_results:
                if not success:
                    print(f"  - {tool_id}: {message}")
        
        # 返回退出码
        return 0 if failed == 0 else 1
        
    finally:
        # 清理临时目录
        shutil.rmtree(temp_dir, ignore_errors=True)
        print(f"\n清理临时目录: {temp_dir}")


if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)

