#!/usr/bin/env python3
"""
测试程序：验证执行和预览模式的严格分离

验证点：
1. Path A (有文件) -> is_template=False -> template_mode=False -> 步骤不为空
2. Path B (无文件) -> is_template=True -> template_mode=True -> 步骤不为空
3. 文件路径正确传递和填充
"""

import asyncio
import sys
import os
import json
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent))

from gibh_agent.core.planner import SOPPlanner
from gibh_agent.core.tool_retriever import ToolRetriever

# Mock LLM client for testing
class MockLLMClient:
    async def chat(self, messages, **kwargs):
        # Return a simple response for intent classification
        if "Intent Classifier" in messages[0]["content"]:
            return {
                "choices": [{
                    "message": {
                        "content": json.dumps({
                            "domain_name": "Metabolomics",
                            "mode": "EXECUTION",
                            "target_steps": []
                        })
                    }
                }]
            }
        return {"choices": [{"message": {"content": "OK"}}]}

async def test_execution_mode_with_file():
    """测试 Path A: 有文件 -> 执行模式"""
    print("\n" + "="*60)
    print("测试 1: Path A - 执行模式（有文件）")
    print("="*60)
    
    # Setup
    tool_retriever = ToolRetriever()
    llm_client = MockLLMClient()
    planner = SOPPlanner(tool_retriever, llm_client)
    
    # Mock file metadata
    file_metadata = {
        "file_path": "/app/uploads/cow_diet.csv",
        "status": "success",
        "n_samples": 39,
        "n_features": 100,
        "file_type": "CSV"
    }
    
    # Test: Generate plan with file (is_template=False)
    result = await planner.generate_plan(
        user_query="分析这个文件",
        file_metadata=file_metadata,
        domain_name="Metabolomics",
        target_steps=None,
        is_template=False  # 🔥 CRITICAL: Execution mode
    )
    
    # Verify
    print(f"\n✅ 结果类型: {result.get('type', 'N/A')}")
    print(f"✅ template_mode: {result.get('template_mode', 'N/A')}")
    
    workflow_data = result.get("workflow_data") or result
    steps = workflow_data.get("steps", [])
    print(f"✅ 步骤数量: {len(steps)}")
    
    # Check file_path in params
    if steps:
        first_step = steps[0]
        params = first_step.get("params", {})
        file_path = params.get("file_path") or params.get("adata_path")
        print(f"✅ 第一个步骤的文件路径: {file_path}")
    
    # Assertions
    assert result.get("template_mode") == False, f"❌ template_mode 应该是 False，但得到 {result.get('template_mode')}"
    assert len(steps) > 0, f"❌ 步骤列表不应该为空，但得到 {len(steps)} 个步骤"
    
    if steps:
        first_step = steps[0]
        params = first_step.get("params", {})
        file_path = params.get("file_path") or params.get("adata_path")
        assert file_path not in ["<待上传数据>", "<PENDING_UPLOAD>", ""], \
            f"❌ 文件路径不应该是占位符，但得到 {file_path}"
        assert file_path == file_metadata["file_path"], \
            f"❌ 文件路径应该匹配，期望 {file_metadata['file_path']}，但得到 {file_path}"
    
    print("\n✅ 测试 1 通过：执行模式正确")

async def test_preview_mode_without_file():
    """测试 Path B: 无文件 -> 预览模式"""
    print("\n" + "="*60)
    print("测试 2: Path B - 预览模式（无文件）")
    print("="*60)
    
    # Setup
    tool_retriever = ToolRetriever()
    llm_client = MockLLMClient()
    planner = SOPPlanner(tool_retriever, llm_client)
    
    # Test: Generate plan without file (is_template=True)
    result = await planner.generate_plan(
        user_query="显示代谢组分析流程",
        file_metadata=None,  # 🔥 CRITICAL: No file
        domain_name="Metabolomics",
        target_steps=None,
        is_template=True  # 🔥 CRITICAL: Template mode
    )
    
    # Verify
    print(f"\n✅ 结果类型: {result.get('type', 'N/A')}")
    print(f"✅ template_mode: {result.get('template_mode', 'N/A')}")
    
    workflow_data = result.get("workflow_data") or result
    steps = workflow_data.get("steps", [])
    print(f"✅ 步骤数量: {len(steps)}")
    
    # Check file_path in params (should be placeholder)
    if steps:
        first_step = steps[0]
        params = first_step.get("params", {})
        file_path = params.get("file_path") or params.get("adata_path")
        print(f"✅ 第一个步骤的文件路径: {file_path}")
    
    # Assertions
    assert result.get("template_mode") == True, f"❌ template_mode 应该是 True，但得到 {result.get('template_mode')}"
    assert len(steps) > 0, f"❌ 步骤列表不应该为空，但得到 {len(steps)} 个步骤"
    
    print("\n✅ 测试 2 通过：预览模式正确")

async def test_file_path_validation():
    """测试文件路径验证"""
    print("\n" + "="*60)
    print("测试 3: 文件路径验证")
    print("="*60)
    
    # Setup
    tool_retriever = ToolRetriever()
    llm_client = MockLLMClient()
    planner = SOPPlanner(tool_retriever, llm_client)
    
    # Mock file metadata
    file_metadata = {
        "file_path": "/app/uploads/test_data.csv",
        "status": "success",
        "n_samples": 10,
        "n_features": 50,
        "file_type": "CSV"
    }
    
    # Test: Generate plan with file
    result = await planner.generate_plan(
        user_query="分析数据",
        file_metadata=file_metadata,
        domain_name="Metabolomics",
        target_steps=None,
        is_template=False
    )
    
    # Verify all steps have correct file_path
    workflow_data = result.get("workflow_data") or result
    steps = workflow_data.get("steps", [])
    
    placeholder_count = 0
    correct_path_count = 0
    
    for step in steps:
        params = step.get("params", {})
        for param_name in ["file_path", "adata_path"]:
            if param_name in params:
                param_value = params[param_name]
                if param_value in ["<待上传数据>", "<PENDING_UPLOAD>", ""]:
                    placeholder_count += 1
                    print(f"⚠️ 步骤 {step.get('id')} 的参数 {param_name} 仍是占位符: {param_value}")
                elif param_value == file_metadata["file_path"]:
                    correct_path_count += 1
                    print(f"✅ 步骤 {step.get('id')} 的参数 {param_name} 正确: {param_value}")
    
    print(f"\n✅ 正确路径数量: {correct_path_count}")
    print(f"⚠️ 占位符数量: {placeholder_count}")
    
    # Assertions
    assert placeholder_count == 0, f"❌ 不应该有占位符，但发现 {placeholder_count} 个"
    assert correct_path_count > 0, f"❌ 应该有正确的文件路径，但发现 {correct_path_count} 个"
    
    print("\n✅ 测试 3 通过：文件路径验证正确")

async def main():
    """运行所有测试"""
    print("\n" + "="*60)
    print("开始测试：执行和预览模式的严格分离")
    print("="*60)
    
    try:
        await test_execution_mode_with_file()
        await test_preview_mode_without_file()
        await test_file_path_validation()
        
        print("\n" + "="*60)
        print("✅ 所有测试通过！")
        print("="*60)
        return 0
    except AssertionError as e:
        print(f"\n❌ 测试失败: {e}")
        return 1
    except Exception as e:
        print(f"\n❌ 测试出错: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    exit_code = asyncio.run(main())
    sys.exit(exit_code)

