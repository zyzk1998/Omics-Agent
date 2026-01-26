#!/usr/bin/env python3
"""
验证 DeepSeek-R1 流程和 Expert Report 生成
"""
import os
import sys
import asyncio
import json
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from gibh_agent.core.llm_client import LLMClient
from gibh_agent.agents.base_agent import BaseAgent
from gibh_agent.agents.specialists.metabolomics_agent import MetabolomicsAgent


async def main():
    print("=" * 80)
    print("DeepSeek-R1 流程验证脚本")
    print("=" * 80)
    print()
    
    # Phase 1: Initialize LLM Client with R1
    print("📋 初始化 LLM 客户端...")
    api_key = os.getenv("DEEPSEEK_API_KEY", os.getenv("SILICONFLOW_API_KEY", os.getenv("LLM_API_KEY", "EMPTY")))
    
    if api_key == "EMPTY":
        print("⚠️  警告: API Key 未设置，使用 'EMPTY' (可能失败)")
        print("   设置方法: export DEEPSEEK_API_KEY='your_key'")
        print()
    
    # Use SiliconFlow DeepSeek-R1 endpoint
    base_url = "https://api.siliconflow.cn/v1"
    model = "deepseek-ai/DeepSeek-R1"
    
    llm_client = LLMClient(
        base_url=base_url,
        api_key=api_key,
        model=model,
        temperature=0.3,
        max_tokens=2500
    )
    
    print(f"   ✅ 使用 LLM: {base_url} (model: {model})")
    print()
    
    # Phase 2: Mock tool results with summary dictionaries
    print("📋 模拟工具输出（包含 summary 字典）...")
    steps_results = [
        {
            "step_name": "pca_analysis",
            "status": "success",
            "data": {
                "explained_variance": {"PC1": 0.452, "PC2": 0.187},
                "plot_path": "/app/results/pca_plot.png",
                "summary": {
                    "pc1_var": 0.452,
                    "pc2_var": 0.187,
                    "separation": "observed",
                    "total_variance_explained": 0.639
                }
            }
        },
        {
            "step_name": "differential_analysis",
            "status": "success",
            "data": {
                "results": [
                    {"metabolite": "Glucose", "log2fc": 2.3, "fdr": 0.001},
                    {"metabolite": "Lactate", "log2fc": -1.8, "fdr": 0.003},
                    {"metabolite": "Pyruvate", "log2fc": 1.5, "fdr": 0.005}
                ],
                "summary": {
                    "total_metabolites": 150,
                    "significant_count": 23,
                    "sig_count": 23,
                    "method": "t-test",
                    "case_group": "Treatment",
                    "control_group": "Control",
                    "top_up": ["Glucose", "Pyruvate", "Citrate"],
                    "top_down": ["Lactate", "Alanine", "Glutamate"]
                }
            }
        },
        {
            "step_name": "pathway_enrichment",
            "status": "success",
            "data": {
                "enriched_pathways": [
                    {"Term": "Glycolysis / Gluconeogenesis", "Adjusted P-value": 0.001},
                    {"Term": "TCA Cycle", "Adjusted P-value": 0.003},
                    {"Term": "Amino Acid Metabolism", "Adjusted P-value": 0.005}
                ],
                "summary": {
                    "n_significant": 3,
                    "n_total": 50,
                    "top_pathways": ["Glycolysis / Gluconeogenesis", "TCA Cycle", "Amino Acid Metabolism"],
                    "p_value_threshold": 0.05
                }
            }
        }
    ]
    
    print(f"   ✅ 模拟了 {len(steps_results)} 个步骤的结果")
    print(f"      - PCA: PC1={steps_results[0]['data']['summary']['pc1_var']:.1%}, PC2={steps_results[0]['data']['summary']['pc2_var']:.1%}")
    print(f"      - Differential: {steps_results[1]['data']['summary']['sig_count']} 个显著差异代谢物")
    print(f"      - Pathways: {len(steps_results[2]['data']['summary']['top_pathways'])} 个富集通路")
    print()
    
    # Phase 3: Initialize Agent and call _generate_analysis_summary
    print("📋 初始化 Agent...")
    agent = MetabolomicsAgent(llm_client=llm_client)
    print(f"   ✅ 使用 {agent.__class__.__name__}")
    print()
    
    # Phase 4: Call _generate_analysis_summary
    print("📞 调用 _generate_analysis_summary...")
    print()
    
    try:
        summary = await agent._generate_analysis_summary(
            steps_results=steps_results,
            omics_type="Metabolomics",
            workflow_name="Metabolomics Analysis",
            summary_context={
                "has_failures": False,
                "has_warnings": False,
                "failed_steps": [],
                "warning_steps": [],
                "successful_steps": steps_results
            }
        )
        
        print("=" * 80)
        print("生成的报告内容：")
        print("=" * 80)
        print(summary)
        print("=" * 80)
        print()
        
        # Phase 5: Assertions
        print("🔍 检查报告质量...")
        
        # Check 1: Contains reasoning tags (proof of R1)
        has_reasoning = False
        reasoning_patterns = ['<think>', '<think>', '<reasoning>']
        for pattern in reasoning_patterns:
            if pattern in summary:
                has_reasoning = True
                print(f"   ✅ 找到 reasoning 标签: {pattern}")
                break
        
        if not has_reasoning:
            print("   ⚠️  未找到 reasoning 标签（可能 R1 未启用或模型未使用 reasoning）")
        
        # Check 2: Contains data injection keywords
        data_keywords = ['Glycolysis', 'Glucose', 'Lactate', 'Pyruvate', 'TCA', '代谢物', '通路', 'PCA', 'PC1', 'PC2']
        found_keywords = [kw for kw in data_keywords if kw in summary]
        print(f"   ✅ 找到 {len(found_keywords)} 个数据关键词: {found_keywords[:10]}")
        
        # Check 3: Report length
        report_length = len(summary)
        print(f"   📏 报告长度: {report_length} 字符")
        
        if report_length < 500:
            print("   ⚠️  报告长度不足 500 字符")
        else:
            print("   ✅ 报告长度充足")
        
        # Check 4: Contains scientific interpretation keywords
        scientific_keywords = ['分析', '生物学', '机制', '显著', '差异', 'pathway', 'PCA', 'VIP', '代谢']
        found_scientific = [kw for kw in scientific_keywords if kw in summary]
        print(f"   ✅ 找到 {len(found_scientific)} 个科学关键词: {found_scientific[:10]}")
        
        # Check 5: Not a fallback/error message
        error_indicators = ['LLM 生成失败', '错误信息', 'Error', '失败', 'fallback']
        has_error = any(indicator in summary for indicator in error_indicators)
        if has_error:
            print("   ⚠️  检测到错误或 fallback 信息")
        else:
            print("   ✅ 未检测到错误或 fallback 信息")
        
        # Check 6: Contains narrative structure
        narrative_indicators = ['通过', '发现', '表明', '显示', '说明', '分析', '解读']
        has_narrative = any(indicator in summary for indicator in narrative_indicators)
        if has_narrative:
            print("   ✅ 检测到叙述性文本结构（科学解释）")
        else:
            print("   ⚠️  未检测到叙述性文本结构")
        
        print()
        print("=" * 80)
        
        # Final verdict
        all_checks = [
            report_length >= 500,
            len(found_keywords) >= 3,
            len(found_scientific) >= 5,
            not has_error,
            has_narrative
        ]
        
        if all(all_checks):
            print("✅ 报告质量验证通过！")
            print("   - 长度充足")
            print("   - 包含数据关键词")
            print("   - 包含科学关键词")
            print("   - 无错误信息")
            print("   - 具有叙述性结构")
            if has_reasoning:
                print("   - ✅ 包含 reasoning 标签（R1 工作正常）")
            else:
                print("   - ⚠️  未检测到 reasoning 标签（可能需要检查模型配置）")
        else:
            print("❌ 报告质量验证未完全通过")
            print("   请检查上述警告项")
        
        print("=" * 80)
        
    except Exception as e:
        print(f"❌ 错误: {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    return 0


if __name__ == "__main__":
    exit_code = asyncio.run(main())
    sys.exit(exit_code)

