#!/usr/bin/env python3
"""
LLM真实连接诊断脚本 - 测试实际的LLM连接和报告生成

目标：
1. 测试真实的LLM连接（不mock）
2. 诊断LLM失败原因（Auth? Context Window? Timeout?）
3. 模拟报告生成流程
"""

import os
import sys
import asyncio
import json
from pathlib import Path
from dotenv import load_dotenv

# 添加项目根目录到路径
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

# 🔥 TASK 4: Load .env explicitly
env_path = project_root / ".env"
if env_path.exists():
    load_dotenv(env_path)
    print(f"✅ 已加载 .env 文件: {env_path}")
else:
    print(f"⚠️ .env 文件不存在: {env_path}")
    print("   将从环境变量读取配置")

from gibh_agent.core.llm_client import LLMClient
from gibh_agent.agents.base_agent import BaseAgent
from gibh_agent.core.prompt_manager import create_default_prompt_manager
from gibh_agent.main import create_agent


async def test_llm_connection():
    """测试1: 测试基本的LLM连接"""
    print("=" * 80)
    print("测试 1: LLM 基本连接测试")
    print("=" * 80)
    
    # 🔥 TASK 4: Initialize Client with deepseek-ai/DeepSeek-R1
    base_url = os.getenv("SILICONFLOW_BASE_URL", "https://api.siliconflow.cn/v1")
    api_key = os.getenv("SILICONFLOW_API_KEY", "")
    model = os.getenv("SILICONFLOW_MODEL", "deepseek-ai/DeepSeek-R1")
    
    print(f"\n📊 LLM配置:")
    print(f"   - Base URL: {base_url}")
    print(f"   - Model: {model}")
    print(f"   - API Key: {'✅ 已设置' if api_key else '❌ 未设置'}")
    if api_key:
        print(f"   - API Key长度: {len(api_key)} 字符")
        print(f"   - API Key前缀: {api_key[:10]}...")
    
    if not api_key:
        print("\n❌ API Key未设置！请设置环境变量 SILICONFLOW_API_KEY")
        return False
    
    try:
        llm_client = LLMClient(
            base_url=base_url,
            api_key=api_key,
            model=model,
            timeout=180.0  # 增加超时时间，DeepSeek-R1可能需要更长时间
        )
        print(f"\n✅ LLM客户端创建成功")
        print(f"   - Timeout: {llm_client.timeout}秒")
    except Exception as e:
        print(f"\n❌ LLM客户端创建失败: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    # 🔥 TASK 4: Test Call - Send simple prompt
    print(f"\n📞 发送测试消息: 'Hello, are you R1?'")
    try:
        messages = [
            {"role": "user", "content": "Hello, are you R1? Please respond briefly."}
        ]
        
        print("   等待响应...")
        completion = await llm_client.achat(messages, temperature=0.3, max_tokens=100)
        
        # Extract response
        think_content, response = llm_client.extract_think_and_content(completion)
        
        print(f"\n✅ LLM响应成功:")
        print(f"   - 响应长度: {len(response)} 字符")
        print(f"   - 思考内容长度: {len(think_content)} 字符")
        print(f"\n📝 响应内容:")
        print("-" * 80)
        print(response[:500] + "..." if len(response) > 500 else response)
        print("-" * 80)
        
        if think_content:
            print(f"\n💭 思考过程 (前200字符):")
            print("-" * 80)
            print(think_content[:200] + "..." if len(think_content) > 200 else think_content)
            print("-" * 80)
        
        return True
        
    except Exception as e:
        print(f"\n❌ LLM调用失败: {e}")
        print(f"\n🔍 错误类型: {type(e).__name__}")
        
        # 🔥 TASK 4: Print debug info
        error_str = str(e).lower()
        if "timeout" in error_str or "timed out" in error_str:
            print("\n💡 诊断: 请求超时")
            print("   - 可能原因: 网络连接慢或LLM服务响应慢")
            print("   - 建议: 增加timeout时间或检查网络连接")
        elif "auth" in error_str or "unauthorized" in error_str or "401" in error_str:
            print("\n💡 诊断: 认证失败")
            print("   - 可能原因: API Key无效或过期")
            print("   - 建议: 检查SILICONFLOW_API_KEY是否正确")
        elif "rate limit" in error_str or "429" in error_str:
            print("\n💡 诊断: 请求频率限制")
            print("   - 可能原因: API调用频率过高")
            print("   - 建议: 等待一段时间后重试")
        elif "context" in error_str or "token" in error_str:
            print("\n💡 诊断: 上下文窗口问题")
            print("   - 可能原因: 输入内容过长或模型不支持")
            print("   - 建议: 减少输入内容长度")
        else:
            print("\n💡 诊断: 未知错误")
            print("   - 请查看完整错误信息")
        
        import traceback
        print("\n📋 完整错误堆栈:")
        traceback.print_exc()
        
        return False


async def test_report_generation():
    """测试2: 模拟报告生成（包含失败步骤）"""
    print("\n" + "=" * 80)
    print("测试 2: 模拟报告生成（包含失败步骤）")
    print("=" * 80)
    
    try:
        # 创建Agent实例
        agent_wrapper = create_agent(config_path="gibh_agent/config/settings.yaml")
        agent = agent_wrapper.agents.get("metabolomics_agent")
        
        if not agent:
            print("❌ 无法获取metabolomics_agent")
            return False
        
        print("✅ Agent创建成功")
        
        # 🔥 TASK 4: Construct dummy steps_results with a failed step
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
                        {"metabolite": "M1", "log2fc": 2.1, "fdr": 0.001, "significant": True},
                    ],
                    "summary": {
                        "total_metabolites": 500,
                        "significant_count": 23,
                        "method": "t-test",
                        "case_group": "Treatment",
                        "control_group": "Control",
                        "top_up": ["M1", "M2"],
                        "top_down": ["M3", "M4"]
                    }
                }
            },
            {
                "step_name": "metabolomics_pathway_enrichment",
                "status": "error",  # 🔥 TASK 4: Simulate failed step
                "error": "Group column 'Sample' not found",
                "data": {}
            }
        ]
        
        summary_context = {
            "has_failures": True,  # 🔥 TASK 4: Has failures
            "has_warnings": False,
            "failed_steps": [
                {
                    "name": "metabolomics_pathway_enrichment",
                    "status": "error",
                    "error": "Group column 'Sample' not found"
                }
            ],
            "warning_steps": [],
            "successful_steps": [s for s in steps_results if s.get("status") == "success"]
        }
        
        print(f"\n📊 模拟数据:")
        print(f"   - 成功步骤: {len([s for s in steps_results if s.get('status') == 'success'])}")
        print(f"   - 失败步骤: {len([s for s in steps_results if s.get('status') == 'error'])}")
        print(f"   - 失败步骤名称: {summary_context['failed_steps'][0]['name']}")
        
        # 🔥 TASK 4: Call _generate_analysis_summary logic directly
        print(f"\n📞 调用 _generate_analysis_summary...")
        print("   等待LLM生成报告...")
        
        summary = await agent._generate_analysis_summary(
            steps_results=steps_results,
            omics_type="Metabolomics",
            workflow_name="Metabolomics Analysis",
            summary_context=summary_context
        )
        
        # 验证结果
        if not summary:
            print("\n❌ 报告摘要为空（返回None）")
            return False
        
        print(f"\n✅ 报告生成成功:")
        print(f"   - 长度: {len(summary)} 字符")
        
        # 检查是否包含失败信息
        if "失败" in summary or "failed" in summary.lower():
            print(f"   - ✅ 包含失败步骤信息")
        
        # 检查是否包含成功步骤信息
        if "PCA" in summary or "差异" in summary or "分析" in summary:
            print(f"   - ✅ 包含成功步骤分析")
        
        print(f"\n📝 报告预览 (前800字符):")
        print("-" * 80)
        print(summary[:800] + "..." if len(summary) > 800 else summary)
        print("-" * 80)
        
        return True
        
    except Exception as e:
        print(f"\n❌ 报告生成失败: {e}")
        import traceback
        traceback.print_exc()
        return False


async def main():
    """主测试函数"""
    print("=" * 80)
    print("LLM真实连接诊断脚本")
    print("=" * 80)
    print("\n此脚本将:")
    print("1. 测试真实的LLM连接（不mock）")
    print("2. 诊断LLM失败原因")
    print("3. 模拟报告生成流程（包含失败步骤）")
    print("\n注意: 需要有效的LLM配置（SILICONFLOW_API_KEY等）")
    print("=" * 80)
    
    results = []
    
    # 测试1: LLM连接
    try:
        result1 = await test_llm_connection()
        results.append(("LLM连接测试", result1))
    except Exception as e:
        print(f"❌ 测试1失败: {e}")
        import traceback
        traceback.print_exc()
        results.append(("LLM连接测试", False))
    
    # 测试2: 报告生成（仅在连接成功时测试）
    if results[0][1]:  # 如果连接测试通过
        try:
            result2 = await test_report_generation()
            results.append(("报告生成测试", result2))
        except Exception as e:
            print(f"❌ 测试2失败: {e}")
            import traceback
            traceback.print_exc()
            results.append(("报告生成测试", False))
    else:
        print("\n⚠️ 跳过报告生成测试（LLM连接失败）")
        results.append(("报告生成测试", False))
    
    # 汇总结果
    print("\n" + "=" * 80)
    print("测试结果汇总")
    print("=" * 80)
    for test_name, passed in results:
        status = "✅ 通过" if passed else "❌ 失败"
        print(f"{test_name}: {status}")
    
    all_passed = all(result for _, result in results)
    print("\n" + "=" * 80)
    if all_passed:
        print("✅ 所有测试通过！LLM连接正常。")
    else:
        print("❌ 部分测试失败，请查看上述诊断信息。")
    print("=" * 80)
    
    return 0 if all_passed else 1


if __name__ == "__main__":
    exit_code = asyncio.run(main())
    sys.exit(exit_code)
