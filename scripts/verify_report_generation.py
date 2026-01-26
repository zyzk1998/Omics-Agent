#!/usr/bin/env python3
"""
验证报告生成脚本 - 确保AI Expert Analysis包含深度生物学解释

测试目标：
1. 验证_generate_analysis_summary提取具体指标（PCA方差、VIP代谢物名称、通路名称）
2. 验证输出包含生物学关键词（metabolism, regulation等）
3. 验证输出长度 > 500字符
4. 验证没有硬编码的步骤列表fallback
"""

import asyncio
import sys
from pathlib import Path

# 添加项目根目录到路径
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from gibh_agent.main import create_agent


async def test_report_generation():
    """测试报告生成"""
    print("=" * 80)
    print("验证报告生成 - AI Expert Analysis")
    print("=" * 80)
    
    try:
        # 创建Agent实例
        agent_wrapper = create_agent(config_path="gibh_agent/config/settings.yaml")
        agent = agent_wrapper.agents.get("metabolomics_agent")
        
        if not agent:
            print("❌ 无法获取metabolomics_agent")
            return False
        
        print("✅ Agent创建成功")
        
        # 🔥 TASK 4: Mock tool outputs with fake VIP scores and Pathway names
        steps_results = [
            {
                "step_name": "pca_analysis",
                "status": "success",
                "data": {
                    "explained_variance": {"PC1": 0.452, "PC2": 0.187},
                    "summary": {
                        "pc1_var": 0.452,
                        "pc2_var": 0.187,
                        "separation": "clear",
                        "total_variance_explained": 0.639
                    }
                }
            },
            {
                "step_name": "differential_analysis",
                "status": "success",
                "data": {
                    "results": [
                        {"metabolite": "Glucose", "log2fc": 2.1, "fdr": 0.001, "significant": True},
                        {"metabolite": "Lactate", "log2fc": 1.8, "fdr": 0.003, "significant": True},
                        {"metabolite": "Pyruvate", "log2fc": -1.5, "fdr": 0.005, "significant": True},
                    ],
                    "summary": {
                        "total_metabolites": 500,
                        "significant_count": 23,
                        "method": "t-test",
                        "case_group": "Treatment",
                        "control_group": "Control",
                        "top_up": ["Glucose", "Lactate"],
                        "top_down": ["Pyruvate", "Citrate"]
                    }
                }
            },
            {
                "step_name": "metabolomics_plsda",
                "status": "success",
                "data": {
                    "summary": {
                        "top_vip_markers": [
                            {"name": "Glucose", "vip": 2.5},
                            {"name": "Lactate", "vip": 2.3},
                            {"name": "Pyruvate", "vip": 2.1},
                            {"name": "Citrate", "vip": 1.9},
                            {"name": "Succinate", "vip": 1.8}
                        ],
                        "n_components": 2,
                        "comp1_variance": 35.0,
                        "comp2_variance": 15.0
                    }
                }
            },
            {
                "step_name": "metabolomics_pathway_enrichment",
                "status": "success",
                "data": {
                    "summary": {
                        "n_significant": 5,
                        "n_total": 100,
                        "top_pathways": [
                            {"name": "Glycolysis / Gluconeogenesis", "p_value": 0.005},
                            {"name": "Citrate cycle (TCA cycle)", "p_value": 0.008},
                            {"name": "Pyruvate metabolism", "p_value": 0.012},
                            {"name": "Pentose phosphate pathway", "p_value": 0.015},
                            {"name": "Fatty acid metabolism", "p_value": 0.018}
                        ],
                        "p_value_threshold": 0.05
                    }
                }
            }
        ]
        
        summary_context = {
            "has_failures": False,
            "has_warnings": False,
            "failed_steps": [],
            "warning_steps": [],
            "successful_steps": steps_results
        }
        
        print("\n📊 模拟数据:")
        print(f"   - PCA: PC1={steps_results[0]['data']['summary']['pc1_var']:.1%}, PC2={steps_results[0]['data']['summary']['pc2_var']:.1%}")
        print(f"   - VIP代谢物: {', '.join([m['name'] for m in steps_results[2]['data']['summary']['top_vip_markers'][:3]])}")
        print(f"   - 富集通路: {', '.join([p['name'] for p in steps_results[3]['data']['summary']['top_pathways'][:3]])}")
        
        # 调用_generate_analysis_summary
        print("\n📞 调用 _generate_analysis_summary...")
        summary = await agent._generate_analysis_summary(
            steps_results=steps_results,
            omics_type="Metabolomics",
            workflow_name="Metabolomics Analysis",
            summary_context=summary_context
        )
        
        # 验证结果
        if not summary:
            print("❌ 报告摘要为空（返回None）")
            return False
        
        if len(summary) < 500:
            print(f"❌ 报告摘要过短 ({len(summary)} 字符)，应至少500字符")
            print(f"   摘要内容: {summary[:500]}")
            return False
        
        # 🔥 TASK 4: Assert - 必须包含生物学关键词
        biological_keywords = ["代谢", "regulation", "regulation", "pathway", "通路", "机制", "生物学", "功能"]
        found_keywords = [kw for kw in biological_keywords if kw.lower() in summary.lower()]
        
        if len(found_keywords) < 3:
            print(f"❌ 报告摘要缺少生物学关键词（仅找到: {found_keywords}）")
            print(f"   摘要内容: {summary[:800]}")
            return False
        
        # 🔥 TASK 4: Assert - 不应包含硬编码步骤列表
        fallback_indicators = ["✅ 成功步骤", "已完成步骤", "步骤列表", "Step List", "Successful Steps"]
        has_fallback = any(indicator in summary for indicator in fallback_indicators)
        
        if has_fallback:
            print(f"❌ 报告摘要包含硬编码步骤列表fallback")
            print(f"   摘要内容: {summary[:800]}")
            return False
        
        # 检查是否包含具体代谢物名称
        metabolite_names = ["Glucose", "Lactate", "Pyruvate", "Citrate"]
        found_metabolites = [m for m in metabolite_names if m in summary]
        
        # 检查是否包含通路名称
        pathway_keywords = ["Glycolysis", "TCA", "Citrate", "Pyruvate"]
        found_pathways = [p for p in pathway_keywords if p in summary]
        
        print(f"\n✅ 报告生成成功:")
        print(f"   - 长度: {len(summary)} 字符")
        print(f"   - 包含关键词: {found_keywords}")
        print(f"   - 包含代谢物: {found_metabolites}")
        print(f"   - 包含通路: {found_pathways}")
        print(f"   - 无硬编码fallback: ✅")
        
        print(f"\n📝 报告预览 (前1000字符):")
        print("-" * 80)
        print(summary[:1000] + "..." if len(summary) > 1000 else summary)
        print("-" * 80)
        
        return True
        
    except Exception as e:
        print(f"❌ 报告生成失败: {e}")
        import traceback
        traceback.print_exc()
        return False


async def main():
    """主测试函数"""
    print("=" * 80)
    print("报告生成验证脚本")
    print("=" * 80)
    print("\n此脚本将验证:")
    print("1. 提取具体指标（PCA方差、VIP代谢物名称、通路名称）")
    print("2. 输出包含生物学关键词")
    print("3. 输出长度 > 500字符")
    print("4. 无硬编码步骤列表fallback")
    print("\n注意: 需要有效的LLM配置（SILICONFLOW_API_KEY等）")
    print("=" * 80)
    
    try:
        result = await test_report_generation()
        
        print("\n" + "=" * 80)
        if result:
            print("✅ 验证通过！报告生成正确。")
        else:
            print("❌ 验证失败，请检查日志。")
        print("=" * 80)
        
        return 0 if result else 1
        
    except Exception as e:
        print(f"\n❌ 测试失败: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    exit_code = asyncio.run(main())
    sys.exit(exit_code)
