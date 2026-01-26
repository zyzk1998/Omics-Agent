#!/usr/bin/env python3
"""
验证报告质量脚本
测试_generate_analysis_summary是否生成高质量的科学解释报告
"""

import sys
import asyncio
import os
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from gibh_agent.agents.base_agent import BaseAgent
from gibh_agent.core.llm_client import LLMClient
from gibh_agent.core.prompt_manager import PromptManager


async def test_report_generation():
    """测试报告生成"""
    print("=" * 80)
    print("报告质量验证脚本")
    print("=" * 80)
    
    # Mock steps_results data (simulating successful Metabolomics run)
    mock_steps_results = [
        {
            "step_name": "metabolomics_preprocess_data",
            "status": "success",
            "data": {
                "shape": {"rows": 100, "columns": 500},
                "summary": {
                    "n_samples": 100,
                    "n_features": 500,
                    "missing_rate": 0.05
                }
            }
        },
        {
            "step_name": "metabolomics_pca_analysis",
            "status": "success",
            "data": {
                "explained_variance": {
                    "PC1": 0.452,
                    "PC2": 0.187
                },
                "summary": {
                    "pc1_variance": "45.2%",
                    "pc2_variance": "18.7%"
                }
            }
        },
        {
            "step_name": "metabolomics_differential_analysis",
            "status": "success",
            "data": {
                "results": [
                    {
                        "metabolite": "Glutamate",
                        "log2fc": 2.34,
                        "fdr": 0.001,
                        "pvalue": 0.0001
                    },
                    {
                        "metabolite": "Alanine",
                        "log2fc": 1.89,
                        "fdr": 0.003,
                        "pvalue": 0.0005
                    },
                    {
                        "metabolite": "Leucine",
                        "log2fc": 1.67,
                        "fdr": 0.005,
                        "pvalue": 0.001
                    }
                ],
                "summary": {
                    "significant_count": 23,
                    "total_metabolites": 500,
                    "method": "t-test",
                    "case_group": "Treatment",
                    "control_group": "Control"
                }
            }
        },
        {
            "step_name": "metabolomics_plsda",
            "status": "success",
            "data": {
                "vip_scores": [
                    {
                        "metabolite": "Glutamate",
                        "vip_score": 2.45
                    },
                    {
                        "metabolite": "Alanine",
                        "vip_score": 2.12
                    },
                    {
                        "metabolite": "Leucine",
                        "vip_score": 1.98
                    }
                ]
            }
        },
        {
            "step_name": "metabolomics_pathway_enrichment",
            "status": "success",
            "data": {
                "enriched_pathways": [
                    {
                        "pathway": "Amino acid metabolism",
                        "p_value": 0.0001,
                        "enrichment_score": 0.85
                    },
                    {
                        "pathway": "Fatty acid synthesis",
                        "p_value": 0.0005,
                        "enrichment_score": 0.72
                    },
                    {
                        "pathway": "TCA cycle",
                        "p_value": 0.001,
                        "enrichment_score": 0.68
                    }
                ]
            }
        }
    ]
    
    # Initialize LLM client and agent
    print("\n📋 初始化LLM客户端和Agent...")
    try:
        # Get LLM config from environment or use defaults
        import os
        base_url = os.getenv("LLM_BASE_URL", "https://api.deepseek.com/v1")
        api_key = os.getenv("DEEPSEEK_API_KEY", os.getenv("LLM_API_KEY", "EMPTY"))
        model = os.getenv("LLM_MODEL", "deepseek-chat")
        
        print(f"  使用 LLM: {base_url} (model: {model})")
        
        llm_client = LLMClient(
            base_url=base_url,
            api_key=api_key,
            model=model
        )
        prompt_manager = PromptManager()
        
        # Create a mock agent (we'll use MetabolomicsAgent)
        try:
            from gibh_agent.agents.specialists.metabolomics_agent import MetabolomicsAgent
            agent = MetabolomicsAgent(llm_client, prompt_manager)
            print("✅ 使用 MetabolomicsAgent")
        except ImportError as e:
            print(f"⚠️ 无法导入 MetabolomicsAgent: {e}")
            # Create a minimal mock agent that has _generate_analysis_summary
            class MockAgent(BaseAgent):
                async def process_query(self, query, **kwargs):
                    return {"type": "chat", "response": "Mock"}
            agent = MockAgent(llm_client, prompt_manager, "metabolomics_expert")
            print("✅ 使用 MockAgent")
    except Exception as e:
        print(f"❌ 初始化失败: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    # Call _generate_analysis_summary
    print("\n📞 调用 _generate_analysis_summary...")
    try:
            summary_context = {
                "has_failures": False,
                "has_warnings": False,
                "failed_steps": [],
                "warning_steps": [],
                "successful_steps": mock_steps_results,
                "workflow_status": "success"
            }
            
            # Create results dict in the format expected by _generate_analysis_summary
            # The method expects steps_results to be a list of step result dicts
            results = {
                "steps_results": mock_steps_results,  # This should be a list of dicts
                "status": "success",
                "steps_details": [
                    {
                        "step_id": f"step_{i}",
                        "step_name": step["step_name"],
                        "status": step["status"],
                        "step_result": step  # Wrap in step_result
                    }
                    for i, step in enumerate(mock_steps_results)
                ]
            }
            
            # Extract steps_results from results dict (as _generate_analysis_summary expects a list)
            steps_results_list = results.get("steps_results", mock_steps_results)
            
            report = await agent._generate_analysis_summary(
                steps_results_list,  # Pass list of step results
                "Metabolomics",
                summary_context=summary_context
            )
    except Exception as e:
        print(f"❌ 报告生成异常: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    if not report:
        print("❌ 报告生成失败：返回 None")
        return False
    
    print("\n" + "=" * 80)
    print("生成的报告内容：")
    print("=" * 80)
    print(report)
    print("=" * 80)
    
    # Assert Content: Check for keywords
    print("\n🔍 检查报告质量...")
    
    # Keywords that indicate scientific interpretation
    scientific_keywords = [
        "分析", "分析结果", "结果摘要",
        "生物学", "机制", "机制解读", "生物学机制",
        "显著", "差异", "差异分析",
        "通路", "代谢通路", "pathway",
        "标志物", "潜在标志物", "biomarker",
        "建议", "下一步", "验证",
        "PCA", "主成分", "PC1", "PC2",
        "VIP", "差异表达", "富集"
    ]
    
    found_keywords = []
    report_lower = report.lower()
    for keyword in scientific_keywords:
        if keyword.lower() in report_lower:
            found_keywords.append(keyword)
    
    print(f"  ✅ 找到 {len(found_keywords)} 个科学关键词: {found_keywords[:10]}")
    
    # Assert Length: Should be substantial
    report_length = len(report)
    print(f"  📏 报告长度: {report_length} 字符")
    
    if report_length < 200:
        print(f"  ❌ 报告过短（< 200字符），可能只是步骤列表")
        return False
    
    # Check for fallback indicators (bad signs)
    fallback_indicators = [
        "✅ 成功步骤",
        "成功步骤",
        "失败步骤",
        "跳过步骤",
        "步骤列表",
        "LLM 生成失败",
        "Error",
        "错误信息"
    ]
    
    has_fallback = any(indicator in report for indicator in fallback_indicators)
    if has_fallback:
        print(f"  ⚠️  检测到fallback或错误信息")
        # Check if it's an error message (which is OK) or a lazy fallback (which is bad)
        if "LLM 生成失败" in report or "错误信息" in report:
            print(f"  ⚠️  这是LLM失败的错误信息（可以接受，但需要修复）")
        else:
            print(f"  ❌ 这是lazy fallback（不应该出现）")
            return False
    
    # Check for narrative structure (good signs)
    narrative_indicators = [
        "通过",
        "发现",
        "表明",
        "显示",
        "识别",
        "富集",
        "相关",
        "可能",
        "建议"
    ]
    
    has_narrative = any(indicator in report for indicator in narrative_indicators)
    if has_narrative:
        print(f"  ✅ 检测到叙述性文本结构（科学解释）")
    else:
        print(f"  ⚠️  缺少叙述性文本结构")
    
    # Final verdict
    print("\n" + "=" * 80)
    if report_length >= 200 and len(found_keywords) >= 5 and has_narrative:
        print("✅ 报告质量验证通过！")
        print("   - 长度充足")
        print("   - 包含科学关键词")
        print("   - 具有叙述性结构")
        return True
    else:
        print("❌ 报告质量验证失败")
        print(f"   - 长度: {report_length} {'✅' if report_length >= 200 else '❌'}")
        print(f"   - 关键词: {len(found_keywords)} {'✅' if len(found_keywords) >= 5 else '❌'}")
        print(f"   - 叙述性: {'✅' if has_narrative else '❌'}")
        return False


async def main():
    """主函数"""
    success = await test_report_generation()
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    asyncio.run(main())

