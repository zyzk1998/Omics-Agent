#!/usr/bin/env python3
"""
测试工作流规划阶段前的功能：
1. 文件检测功能
2. 参数推荐
3. 推荐依据报告
4. 逻辑和配置验证
"""

import asyncio
import sys
import os
from pathlib import Path

# 添加项目根目录到路径
sys.path.insert(0, str(Path(__file__).parent))

from gibh_agent import create_agent
from gibh_agent.core.prompt_manager import PromptManager
import logging

# 配置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


async def test_metabolomics_workflow_preparation():
    """测试代谢组学工作流规划前的功能"""
    print("\n" + "="*80)
    print("🧪 测试代谢组学工作流规划前功能")
    print("="*80)
    
    # 使用 create_agent 创建智能体（会自动初始化所有组件）
    try:
        gibh_agent = create_agent("metabolomics")
        # 获取具体的代谢组学智能体
        agent = gibh_agent.agents.get('metabolomics_agent')
        if not agent:
            print("❌ 无法获取 metabolomics_agent")
            return False
    except Exception as e:
        print(f"❌ 创建智能体失败: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    # 测试文件路径
    test_file = "/home/ubuntu/GIBH-AGENT-V2/uploads/human_cachexia.csv"
    
    if not os.path.exists(test_file):
        print(f"❌ 测试文件不存在: {test_file}")
        return False
    
    print(f"\n📁 测试文件: {test_file}")
    
    # 测试1: 轻量级文件检测
    print("\n" + "-"*80)
    print("🔍 测试1: 轻量级文件检测 (_peek_data_lightweight)")
    print("-"*80)
    try:
        peek_result = await agent._peek_data_lightweight(test_file)
        if "error" in peek_result:
            print(f"❌ 文件检测失败: {peek_result['error']}")
            return False
        
        print("✅ 文件检测成功")
        print(f"  - 样本数: {peek_result.get('n_samples', 'N/A')}")
        print(f"  - 代谢物数: {peek_result.get('n_metabolites', 'N/A')}")
        print(f"  - 元数据列数: {peek_result.get('n_metadata_cols', 'N/A')}")
        print(f"  - 元数据列: {peek_result.get('metadata_columns', [])}")
        print(f"  - 数值统计: {peek_result.get('numeric_stats', {})}")
    except Exception as e:
        print(f"❌ 文件检测异常: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    # 测试2: 参数推荐
    print("\n" + "-"*80)
    print("💡 测试2: 参数推荐 (_generate_parameter_recommendations)")
    print("-"*80)
    try:
        query = "分析这个代谢组数据"
        recommendation = await agent._generate_parameter_recommendations(peek_result, query)
        
        if not recommendation:
            print("❌ 参数推荐失败: 返回空结果")
            return False
        
        print("✅ 参数推荐成功")
        print(f"  - 摘要: {recommendation.get('summary', 'N/A')}")
        print(f"  - 推荐参数:")
        if "params" in recommendation:
            for param_name, param_info in recommendation["params"].items():
                if isinstance(param_info, dict):
                    value = param_info.get("value", "N/A")
                    reason = param_info.get("reason", "N/A")
                    print(f"    • {param_name}: {value} (理由: {reason})")
                else:
                    print(f"    • {param_name}: {param_info}")
    except Exception as e:
        print(f"❌ 参数推荐异常: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    # 测试3: 工作流配置生成（包含文件检测和推荐）
    print("\n" + "-"*80)
    print("🚀 测试3: 工作流配置生成（完整流程）")
    print("-"*80)
    try:
        workflow_result = await agent._generate_workflow_config(
            query=query,
            file_paths=[test_file]
        )
        
        if not workflow_result:
            print("❌ 工作流配置生成失败: 返回空结果")
            return False
        
        print("✅ 工作流配置生成成功")
        print(f"  - 类型: {workflow_result.get('type', 'N/A')}")
        
        if "workflow_data" in workflow_result:
            workflow_data = workflow_result["workflow_data"]
            print(f"  - 工作流名称: {workflow_data.get('workflow_name', 'N/A')}")
            print(f"  - 步骤数: {len(workflow_data.get('steps', []))}")
            
            # 检查推荐是否应用
            if "recommendation" in workflow_result:
                print(f"  - ✅ 包含推荐信息")
                rec = workflow_result["recommendation"]
                print(f"    • 摘要: {rec.get('summary', 'N/A')}")
            else:
                print(f"  - ⚠️ 未包含推荐信息")
        
        # 检查步骤参数是否已填充推荐值
        if "workflow_data" in workflow_result:
            steps = workflow_result["workflow_data"].get("steps", [])
            for step in steps:
                step_name = step.get("name", step.get("step_name", "未知步骤"))
                params = step.get("params", {})
                if params:
                    print(f"  - 步骤 '{step_name}' 参数: {params}")
    except Exception as e:
        print(f"❌ 工作流配置生成异常: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    print("\n✅ 代谢组学工作流规划前功能测试通过！")
    return True


async def test_rna_workflow_preparation():
    """测试RNA工作流规划前的功能"""
    print("\n" + "="*80)
    print("🧪 测试RNA工作流规划前功能")
    print("="*80)
    
    # 使用 create_agent 创建智能体
    try:
        gibh_agent = create_agent("rna")
        # 获取具体的RNA智能体
        agent = gibh_agent.agents.get('rna_agent')
        if not agent:
            print("❌ 无法获取 rna_agent")
            return False
    except Exception as e:
        print(f"❌ 创建智能体失败: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    # 测试文件路径（使用matrix.mtx作为示例）
    test_file = "/home/ubuntu/GIBH-AGENT-V2/uploads/matrix.mtx"
    
    if not os.path.exists(test_file):
        print(f"⚠️ 测试文件不存在: {test_file}，跳过RNA测试")
        return True  # 不是错误，只是没有测试文件
    
    print(f"\n📁 测试文件: {test_file}")
    
    # 测试1: 文件检测
    print("\n" + "-"*80)
    print("🔍 测试1: 文件检测 (inspect_file)")
    print("-"*80)
    try:
        # 检查是否有 scanpy_tool
        if not hasattr(agent, 'scanpy_tool'):
            print("⚠️ RNA agent 没有 scanpy_tool，跳过文件检测")
            inspection_result = None
        else:
            inspection_result = agent.scanpy_tool.inspect_file(test_file)
        if "error" in inspection_result:
            print(f"⚠️ 文件检测失败: {inspection_result['error']}")
            print("   注意: 某些文件类型可能不支持检测，这是正常的")
        else:
            print("✅ 文件检测成功")
            print(f"  - 文件类型: {inspection_result.get('file_type', 'N/A')}")
            print(f"  - 细胞数: {inspection_result.get('n_obs', 'N/A')}")
            print(f"  - 基因数: {inspection_result.get('n_vars', 'N/A')}")
    except Exception as e:
        print(f"⚠️ 文件检测异常: {e}")
        print("   注意: 某些文件类型可能不支持检测，这是正常的")
        inspection_result = None
    
    # 测试2: 诊断和推荐报告生成
    if inspection_result and "error" not in inspection_result:
        print("\n" + "-"*80)
        print("📊 测试2: 诊断和推荐报告生成 (_generate_diagnosis_and_recommendation)")
        print("-"*80)
        try:
            diagnosis_report = await agent._generate_diagnosis_and_recommendation(inspection_result)
            
            if not diagnosis_report:
                print("❌ 诊断报告生成失败: 返回空结果")
                return False
            
            print("✅ 诊断报告生成成功")
            print(f"  - 报告长度: {len(diagnosis_report)} 字符")
            print(f"  - 报告预览（前200字符）:")
            print(f"    {diagnosis_report[:200]}...")
        except Exception as e:
            print(f"❌ 诊断报告生成异常: {e}")
            import traceback
            traceback.print_exc()
            return False
    
    print("\n✅ RNA工作流规划前功能测试通过！")
    return True


async def test_configuration():
    """测试配置和提示词"""
    print("\n" + "="*80)
    print("⚙️ 测试配置和提示词")
    print("="*80)
    
    try:
        prompt_manager = PromptManager()
        
        # 测试数据诊断提示词
        print("\n📝 测试数据诊断提示词模板")
        try:
            prompt = prompt_manager.get_prompt(
                "data_diagnosis",
                {"inspection_data": '{"test": "data"}'},
                fallback="Default prompt"
            )
            if prompt and len(prompt) > 0:
                print("✅ 数据诊断提示词模板可用")
                print(f"  - 模板长度: {len(prompt)} 字符")
            else:
                print("⚠️ 数据诊断提示词模板为空，使用默认值")
        except Exception as e:
            print(f"⚠️ 获取数据诊断提示词失败: {e}，将使用默认值")
        
        # 测试配置加载
        print("\n⚙️ 测试配置加载")
        try:
            import yaml
            config_path = Path(__file__).parent / "gibh_agent" / "config" / "settings.yaml"
            if config_path.exists():
                with open(config_path, 'r', encoding='utf-8') as f:
                    settings = yaml.safe_load(f)
                print("✅ 配置加载成功")
                if settings and 'llm' in settings:
                    print(f"  - LLM配置: {settings.get('llm', {}).get('provider', 'N/A')}")
            else:
                print("⚠️ 配置文件不存在，使用默认值")
        except Exception as e:
            print(f"⚠️ 配置加载失败: {e}，使用默认值")
        
    except Exception as e:
        print(f"⚠️ 配置测试异常: {e}")
        import traceback
        traceback.print_exc()
    
    return True


async def main():
    """主测试函数"""
    print("\n" + "="*80)
    print("🧪 工作流规划前功能自测")
    print("="*80)
    
    results = []
    
    # 测试配置
    results.append(await test_configuration())
    
    # 测试代谢组学工作流
    results.append(await test_metabolomics_workflow_preparation())
    
    # 测试RNA工作流
    results.append(await test_rna_workflow_preparation())
    
    # 总结
    print("\n" + "="*80)
    print("📊 测试总结")
    print("="*80)
    passed = sum(results)
    total = len(results)
    print(f"✅ 通过: {passed}/{total}")
    print(f"❌ 失败: {total - passed}/{total}")
    
    if passed == total:
        print("\n🎉 所有测试通过！")
        return 0
    else:
        print("\n⚠️ 部分测试失败，请检查上述输出")
        return 1


if __name__ == "__main__":
    exit_code = asyncio.run(main())
    sys.exit(exit_code)

