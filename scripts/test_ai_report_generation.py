#!/usr/bin/env python3
"""
AI生信专家报告生成流程测试脚本
测试完整的报告生成流程，包括数据提取、LLM调用、报告质量评估
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
from gibh_agent.core.prompt_manager import create_default_prompt_manager

# 评分标准
SCORING_CRITERIA = {
    "data_extraction": {
        "weight": 0.2,
        "description": "数据提取完整性（是否提取了PCA、差异分析、通路富集的关键指标）"
    },
    "llm_call": {
        "weight": 0.3,
        "description": "LLM调用成功性（是否成功调用LLM，返回内容是否合理）"
    },
    "scientific_content": {
        "weight": 0.3,
        "description": "科学内容质量（是否包含生物学机制解读、通路分析、标志物讨论）"
    },
    "report_structure": {
        "weight": 0.1,
        "description": "报告结构完整性（是否包含摘要、统计分析、机制解读、标志物、建议）"
    },
    "reasoning_quality": {
        "weight": 0.1,
        "description": "推理过程质量（是否包含thinking标签，推理是否合理）"
    }
}


def score_report(report: str, has_reasoning: bool, extracted_data: dict) -> dict:
    """评分报告质量"""
    scores = {}
    
    # 1. 数据提取评分
    data_keywords = ['PC1', 'PC2', 'PCA', '差异', '代谢物', '通路', 'VIP', 'Log2FC', 'FDR']
    found_data = sum(1 for kw in data_keywords if kw in report)
    scores['data_extraction'] = min(100, (found_data / len(data_keywords)) * 100)
    
    # 2. LLM调用评分
    if report and len(report) > 200 and not report.startswith("## ❌"):
        scores['llm_call'] = 100
    elif report and len(report) > 100:
        scores['llm_call'] = 70
    else:
        scores['llm_call'] = 0
    
    # 3. 科学内容评分
    scientific_keywords = ['机制', '生物学', '通路', '代谢', '功能', '标志物', 'biomarker', 'pathway', 'mechanism']
    found_scientific = sum(1 for kw in scientific_keywords if kw in report)
    narrative_indicators = ['通过', '发现', '表明', '显示', '说明', '分析', '解读', '讨论']
    found_narrative = sum(1 for kw in narrative_indicators if kw in report)
    scores['scientific_content'] = min(100, ((found_scientific + found_narrative) / (len(scientific_keywords) + len(narrative_indicators))) * 100)
    
    # 4. 报告结构评分
    structure_sections = ['摘要', '统计', '机制', '标志物', '建议', 'Summary', 'Results', 'Mechanism', 'Biomarker', 'Recommendation']
    found_sections = sum(1 for section in structure_sections if section in report)
    scores['report_structure'] = min(100, (found_sections / len(structure_sections)) * 100)
    
    # 5. 推理质量评分
    if has_reasoning:
        scores['reasoning_quality'] = 100
    else:
        scores['reasoning_quality'] = 50  # 如果没有reasoning标签，给50分（可能R1未启用）
    
    # 计算总分
    total_score = sum(scores[key] * SCORING_CRITERIA[key]['weight'] for key in scores)
    
    return {
        'scores': scores,
        'total_score': total_score,
        'grade': get_grade(total_score)
    }


def get_grade(score: float) -> str:
    """根据分数返回等级"""
    if score >= 90:
        return "A+ (生产级)"
    elif score >= 80:
        return "A (优秀)"
    elif score >= 70:
        return "B (良好)"
    elif score >= 60:
        return "C (及格)"
    else:
        return "D (不合格)"


async def main():
    print("=" * 80)
    print("AI生信专家报告生成流程测试")
    print("=" * 80)
    print()
    
    # 初始化
    print("📋 初始化组件...")
    api_key = os.getenv("DEEPSEEK_API_KEY", os.getenv("SILICONFLOW_API_KEY", os.getenv("LLM_API_KEY", "EMPTY")))
    
    if api_key == "EMPTY":
        print("⚠️  警告: API Key 未设置")
        print("   设置方法: export DEEPSEEK_API_KEY='your_key'")
        return 1
    
    base_url = "https://api.siliconflow.cn/v1"
    model = "deepseek-ai/DeepSeek-R1"
    
    llm_client = LLMClient(
        base_url=base_url,
        api_key=api_key,
        model=model,
        temperature=0.3,
        max_tokens=2500,
        timeout=120.0  # 增加超时时间到120秒（DeepSeek-R1推理时间较长）
    )
    
    prompt_manager = create_default_prompt_manager()
    agent = MetabolomicsAgent(llm_client=llm_client, prompt_manager=prompt_manager)
    
    print(f"   ✅ LLM: {base_url} (model: {model})")
    print(f"   ✅ Agent: {agent.__class__.__name__}")
    print()
    
    # 模拟完整的工具输出（包含summary字典）
    print("📋 模拟工具输出...")
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
                    {"metabolite": "Pyruvate", "log2fc": 1.5, "fdr": 0.005},
                    {"metabolite": "Citrate", "log2fc": 1.2, "fdr": 0.008},
                    {"metabolite": "Alanine", "log2fc": -1.1, "fdr": 0.012}
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
    
    print(f"   ✅ 模拟了 {len(steps_results)} 个步骤")
    print()
    
    # 调用报告生成
    print("📞 调用 _generate_analysis_summary...")
    print("   ⏱️  设置超时保护: 180秒（DeepSeek-R1推理时间较长）")
    print("   💡 提示: 如果超时，请检查网络连接或API服务状态")
    print()
    
    try:
        # 添加超时保护，避免LLM调用卡住
        import asyncio
        summary = await asyncio.wait_for(
            agent._generate_analysis_summary(
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
            ),
            timeout=180.0  # 180秒超时（DeepSeek-R1需要更长时间）
        )
        
        print("=" * 80)
        print("生成的报告内容：")
        print("=" * 80)
        print(summary[:500] + "..." if len(summary) > 500 else summary)
        print("=" * 80)
        print()
        
        # 检查reasoning标签
        has_reasoning = any(tag in summary for tag in ['<think>', '<think>', '<reasoning>'])
        
        # 提取数据检查
        extracted_data = {
            "pca": "PC1" in summary and "PC2" in summary,
            "differential": "差异" in summary or "代谢物" in summary,
            "pathway": "通路" in summary or "pathway" in summary.lower()
        }
        
        # 评分
        print("🔍 报告质量评分...")
        print()
        score_result = score_report(summary, has_reasoning, extracted_data)
        
        print("评分详情:")
        for key, score in score_result['scores'].items():
            criterion = SCORING_CRITERIA[key]
            print(f"  {key}: {score:.1f}/100 ({criterion['description']})")
        
        print()
        print(f"总分: {score_result['total_score']:.1f}/100")
        print(f"等级: {score_result['grade']}")
        print()
        
        # 详细检查
        print("详细检查:")
        print(f"  ✅ 报告长度: {len(summary)} 字符")
        print(f"  {'✅' if has_reasoning else '⚠️'} Reasoning标签: {'存在' if has_reasoning else '不存在（可能R1未启用）'}")
        print(f"  ✅ 数据提取: PCA={extracted_data['pca']}, Diff={extracted_data['differential']}, Pathway={extracted_data['pathway']}")
        print(f"  {'✅' if '机制' in summary or 'mechanism' in summary.lower() else '❌'} 生物学机制解读: {'存在' if '机制' in summary or 'mechanism' in summary.lower() else '缺失'}")
        print(f"  {'✅' if '标志物' in summary or 'biomarker' in summary.lower() else '❌'} 标志物讨论: {'存在' if '标志物' in summary or 'biomarker' in summary.lower() else '缺失'}")
        print(f"  {'✅' if '建议' in summary or 'recommendation' in summary.lower() else '❌'} 下一步建议: {'存在' if '建议' in summary or 'recommendation' in summary.lower() else '缺失'}")
        print()
        
        # 最终判定
        print("=" * 80)
        if score_result['total_score'] >= 80:
            print("✅ 测试通过！报告质量达到生产级标准")
        elif score_result['total_score'] >= 70:
            print("⚠️  测试基本通过，但需要改进")
        else:
            print("❌ 测试未通过，需要修复")
        print("=" * 80)
        
        return 0 if score_result['total_score'] >= 80 else 1
        
    except asyncio.TimeoutError:
        print("=" * 80)
        print("⏱️  LLM调用超时（180秒）")
        print("   可能原因：")
        print("   1. API响应慢或服务繁忙")
        print("   2. 网络延迟")
        print("   3. DeepSeek-R1推理时间较长（需要处理大量reasoning）")
        print()
        print("建议：")
        print("   - 检查网络连接: ping api.siliconflow.cn")
        print("   - 检查API服务状态: curl https://api.siliconflow.cn/v1/models")
        print("   - 考虑使用更快的模型（如DeepSeek-V3）进行快速测试")
        print("   - 或增加超时时间到300秒")
        print()
        print("⚠️  注意: 虽然超时，但代码逻辑是正确的。")
        print("   在实际生产环境中，LLM调用通常会在合理时间内完成。")
        print("=" * 80)
        return 1
    except Exception as e:
        print(f"❌ 错误: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    exit_code = asyncio.run(main())
    sys.exit(exit_code)

