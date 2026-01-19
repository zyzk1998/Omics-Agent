#!/usr/bin/env python3
"""
完整流程测试：验证规划、执行阶段的核心功能

测试场景：
1. 未上传文件 - 完整工作流规划
2. 未上传文件 - 部分工作流规划（PCA）
3. 已上传文件 - 完整工作流规划
4. 已上传文件 - 部分工作流规划（PCA）
5. 验证 SSE 事件格式是否符合前端期望
"""
import sys
import os
import asyncio
import json
from pathlib import Path
from typing import Dict, Any, List

# 添加项目根目录到路径
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from gibh_agent.core.planner import SOPPlanner
from gibh_agent.core.workflows import WorkflowRegistry
from gibh_agent.core.llm_client import LLMClient


class MockLLMClient:
    """Mock LLM 客户端，用于测试"""
    def __init__(self):
        self.call_count = 0
    
    async def achat(self, messages, **kwargs):
        """模拟 LLM 响应"""
        from openai.types.chat import ChatCompletion
        from openai.types.chat.chat_completion import ChatCompletionMessage
        from openai.types.chat.chat_completion import Choice as ChatCompletionChoice
        
        self.call_count += 1
        
        # 根据消息内容返回不同的响应
        user_message = messages[-1]["content"] if messages else ""
        system_message = messages[0]["content"] if messages else ""
        
        # 意图分类响应（SOPPlanner._classify_intent）
        if "Intent Classifier" in system_message and "domain_name" in system_message:
            if "PCA" in user_message or "pca" in user_message.lower():
                response_text = '{"domain_name": "Metabolomics", "target_steps": ["pca_analysis"]}'
            elif "full" in user_message.lower() or "完整" in user_message:
                response_text = '{"domain_name": "Metabolomics", "target_steps": []}'
            else:
                response_text = '{"domain_name": "Metabolomics", "target_steps": []}'
        
        # 意图分析响应（SOPPlanner._analyze_user_intent）
        elif "Intent Analyzer" in system_message and "Target Steps" in system_message:
            if "PCA" in user_message or "pca" in user_message.lower():
                response_text = '["pca_analysis"]'
            elif "full" in user_message.lower() or "完整" in user_message:
                response_text = '[]'
            else:
                response_text = '[]'
        
        # 查询重写响应
        elif "Query Rewriter" in system_message:
            response_text = user_message  # 返回原查询
        
        # 意图检查响应（Clarifier._check_intent）
        elif "Intent Classifier" in system_message and "PLAN" in system_message:
            if "I want PCA" in user_message or "Plan" in user_message or "plan" in user_message.lower():
                response_text = "planning"
            else:
                response_text = "execution"
        
        # 澄清检查响应
        elif "Clarification Assistant" in system_message:
            response_text = "OK"
        
        # 反思响应
        elif "Workflow Reflector" in system_message:
            response_text = '{"is_valid": true, "corrected_plan": null}'
        
        else:
            response_text = "OK"
        
        # 创建 ChatCompletion 对象
        try:
            message = ChatCompletionMessage(role="assistant", content=response_text)
            choice = ChatCompletionChoice(index=0, message=message, finish_reason="stop")
            completion = ChatCompletion(
                id=f"mock-{self.call_count}",
                model="mock-model",
                choices=[choice],
                created=1234567890,
                object="chat.completion"
            )
        except Exception as e:
            # 如果导入失败，创建一个简单的对象
            class SimpleCompletion:
                def __init__(self, content):
                    class Choice:
                        def __init__(self, content):
                            class Message:
                                def __init__(self, content):
                                    self.content = content
                                    self.role = "assistant"
                            self.message = Message(content)
                            self.index = 0
                            self.finish_reason = "stop"
                    self.choices = [Choice(content)]
                    self.id = f"mock-{self.call_count}"
                    self.model = "mock-model"
            
            completion = SimpleCompletion(response_text)
        
        return completion


class TestResult:
    """测试结果"""
    def __init__(self, scenario_name: str):
        self.scenario_name = scenario_name
        self.passed = True
        self.errors = []
        self.warnings = []
        self.events = []
    
    def fail(self, message: str):
        self.passed = False
        self.errors.append(message)
        print(f"  ❌ FAIL: {message}")
    
    def warn(self, message: str):
        self.warnings.append(message)
        print(f"  ⚠️  WARN: {message}")
    
    def success(self, message: str):
        print(f"  ✅ PASS: {message}")
    
    def add_event(self, event_type: str, data: Dict[str, Any]):
        self.events.append({"type": event_type, "data": data})


async def test_scenario_1_no_file_full_workflow():
    """场景1：未上传文件 - 完整工作流规划"""
    result = TestResult("场景1: 未上传文件 - 完整工作流规划")
    
    print(f"\n{'='*60}")
    print(f"测试: {result.scenario_name}")
    print(f"{'='*60}")
    
    try:
        # 创建 Mock LLM
        mock_llm = MockLLMClient()
        
        # 创建 planner
        from gibh_agent.core.tool_retriever import ToolRetriever
        tool_retriever = ToolRetriever()
        planner = SOPPlanner(tool_retriever=tool_retriever, llm_client=mock_llm)
        
        # 调用 generate_plan
        query = "完整分析"
        file_metadata = None
        
        print(f"\n📤 输入:")
        print(f"  Query: {query}")
        print(f"  Files: None")
        
        workflow_config = await planner.generate_plan(
            user_query=query,
            file_metadata=file_metadata
        )
        
        print(f"\n📥 返回结构:")
        print(f"  Type: {type(workflow_config)}")
        print(f"  Keys: {list(workflow_config.keys()) if isinstance(workflow_config, dict) else 'N/A'}")
        
        # 验证结构
        if workflow_config.get("type") == "error":
            result.fail(f"返回错误: {workflow_config.get('error')}")
            return result
        
        # 检查 workflow_data
        workflow_data = workflow_config.get("workflow_data", {})
        if not workflow_data:
            result.fail("缺少 'workflow_data' 字段")
        else:
            result.success("workflow_data 存在")
            
            # 检查 steps
            steps = workflow_data.get("steps", [])
            if not steps:
                result.fail("workflow_data.steps 为空")
            else:
                result.success(f"steps 包含 {len(steps)} 个步骤")
                
                # 验证步骤结构
                for i, step in enumerate(steps):
                    step_id = step.get('step_id') or step.get('id') or step.get('tool_id')
                    if not step_id:
                        result.fail(f"步骤 {i} 缺少 step_id/id/tool_id")
                    
                    # 检查占位符
                    params = step.get('params', {})
                    file_path = params.get('file_path') or params.get('adata_path')
                    if file_path and ('<待上传数据>' not in str(file_path) and '<PENDING_UPLOAD>' not in str(file_path)):
                        result.warn(f"步骤 {i} 的 file_path 不是占位符: {file_path}")
        
        # 检查 template_mode
        template_mode = workflow_config.get("template_mode")
        if template_mode is not True:
            result.fail(f"template_mode 应为 True，实际为: {template_mode}")
        else:
            result.success("template_mode = True")
        
        # 检查 diagnosis
        diagnosis = workflow_config.get("diagnosis")
        if not diagnosis:
            result.fail("缺少 'diagnosis' 字段")
        else:
            result.success("diagnosis 存在")
            if isinstance(diagnosis, dict):
                message = diagnosis.get("message", "")
                if "方案已生成" in message or "Template Ready" in message:
                    result.success("diagnosis 包含模板信息")
                else:
                    result.warn(f"diagnosis 消息可能不正确: {message[:100]}")
        
        result.add_event("workflow", workflow_config)
        
    except Exception as e:
        result.fail(f"测试执行失败: {e}")
        import traceback
        print(f"  详细错误:\n{traceback.format_exc()}")
    
    return result


async def test_scenario_2_no_file_partial_workflow():
    """场景2：未上传文件 - 部分工作流规划（PCA）"""
    result = TestResult("场景2: 未上传文件 - 部分工作流规划（PCA）")
    
    print(f"\n{'='*60}")
    print(f"测试: {result.scenario_name}")
    print(f"{'='*60}")
    
    try:
        # 创建 Mock LLM
        mock_llm = MockLLMClient()
        
        # 创建 planner
        from gibh_agent.core.tool_retriever import ToolRetriever
        tool_retriever = ToolRetriever()
        planner = SOPPlanner(tool_retriever=tool_retriever, llm_client=mock_llm)
        
        # 调用 generate_plan
        query = "I want PCA"
        file_metadata = None
        
        print(f"\n📤 输入:")
        print(f"  Query: {query}")
        print(f"  Files: None")
        
        workflow_config = await planner.generate_plan(
            user_query=query,
            file_metadata=file_metadata
        )
        
        # 验证结构
        if workflow_config.get("type") == "error":
            result.fail(f"返回错误: {workflow_config.get('error')}")
            return result
        
        # 检查 steps 数量（应该是3个：inspect, preprocess, pca）
        workflow_data = workflow_config.get("workflow_data", {})
        steps = workflow_data.get("steps", [])
        
        if len(steps) != 3:
            result.fail(f"部分工作流应包含3个步骤，实际为: {len(steps)}")
        else:
            result.success(f"部分工作流包含 {len(steps)} 个步骤（正确）")
        
        # 检查是否包含 PCA
        has_pca = any(
            step.get('step_id') == 'pca_analysis' or 
            step.get('id') == 'pca_analysis' or 
            step.get('tool_id') == 'pca_analysis'
            for step in steps
        )
        if not has_pca:
            result.fail("工作流中未找到 pca_analysis 步骤")
        else:
            result.success("工作流包含 pca_analysis 步骤")
        
        # 检查 template_mode
        if workflow_config.get("template_mode") is not True:
            result.fail("template_mode 应为 True")
        else:
            result.success("template_mode = True")
        
        result.add_event("workflow", workflow_config)
        
    except Exception as e:
        result.fail(f"测试执行失败: {e}")
        import traceback
        print(f"  详细错误:\n{traceback.format_exc()}")
    
    return result


async def test_scenario_3_with_file_full_workflow(upload_dir: str):
    """场景3：已上传文件 - 完整工作流规划"""
    result = TestResult("场景3: 已上传文件 - 完整工作流规划")
    
    print(f"\n{'='*60}")
    print(f"测试: {result.scenario_name}")
    print(f"{'='*60}")
    
    try:
        # 创建模拟文件
        mock_file_path = Path(upload_dir) / "test_metabolomics.csv"
        mock_file_path.parent.mkdir(parents=True, exist_ok=True)
        mock_file_path.write_text("sample,group,metabolite1,metabolite2\n1,A,1.0,2.0\n2,B,1.5,2.5\n")
        
        # 创建 Mock LLM
        mock_llm = MockLLMClient()
        
        # 创建 planner
        from gibh_agent.core.tool_retriever import ToolRetriever
        tool_retriever = ToolRetriever()
        planner = SOPPlanner(tool_retriever=tool_retriever, llm_client=mock_llm)
        
        # 创建文件元数据
        file_metadata = {
            "file_path": str(mock_file_path),
            "file_type": "csv",
            "semantic_map": {
                "group_cols": ["group"]
            }
        }
        
        # 调用 generate_plan
        query = "完整分析"
        
        print(f"\n📤 输入:")
        print(f"  Query: {query}")
        print(f"  Files: {file_metadata['file_path']}")
        
        workflow_config = await planner.generate_plan(
            user_query=query,
            file_metadata=file_metadata
        )
        
        # 验证结构
        if workflow_config.get("type") == "error":
            result.fail(f"返回错误: {workflow_config.get('error')}")
            return result
        
        # 检查 template_mode 应为 False 或不存在
        template_mode = workflow_config.get("template_mode")
        if template_mode is True:
            result.fail(f"有文件时 template_mode 应为 False 或不存在，实际为: {template_mode}")
        else:
            result.success(f"template_mode = {template_mode} (正确)")
        
        # 检查 file_path 应为真实路径
        workflow_data = workflow_config.get("workflow_data", {})
        steps = workflow_data.get("steps", [])
        has_real_path = False
        for step in steps:
            params = step.get('params', {})
            file_path = params.get('file_path') or params.get('adata_path')
            if file_path and '<待上传数据>' not in str(file_path) and '<PENDING_UPLOAD>' not in str(file_path):
                has_real_path = True
                break
        
        if not has_real_path:
            result.warn("未找到真实 file_path（可能仍为占位符）")
        else:
            result.success("找到真实 file_path")
        
        # 清理
        if mock_file_path.exists():
            mock_file_path.unlink()
        
        result.add_event("workflow", workflow_config)
        
    except Exception as e:
        result.fail(f"测试执行失败: {e}")
        import traceback
        print(f"  详细错误:\n{traceback.format_exc()}")
    
    return result


def verify_frontend_contract(events: List[Dict[str, Any]]) -> TestResult:
    """验证前端契约：检查 SSE 事件格式"""
    result = TestResult("前端契约验证")
    
    print(f"\n{'='*60}")
    print(f"测试: {result.scenario_name}")
    print(f"{'='*60}")
    
    if not events:
        result.fail("没有事件可验证")
        return result
    
    for event in events:
        event_type = event.get("type")
        data = event.get("data", {})
        
        if event_type == "workflow":
            # 检查 workflow 事件结构
            workflow_config = data.get("workflow_config") or data.get("workflow_data") or data
            
            # 检查嵌套结构
            if isinstance(workflow_config, dict) and "workflow_data" in workflow_config:
                nested = workflow_config.get("workflow_data", {})
                if isinstance(nested, dict) and "workflow_data" in nested:
                    result.fail("检测到 workflow_data.workflow_data 嵌套结构（错误）")
                else:
                    result.success("workflow_data 结构正确（无嵌套）")
            
            # 检查步骤结构
            steps = workflow_config.get("steps") or workflow_config.get("workflow_data", {}).get("steps", [])
            for i, step in enumerate(steps):
                if not step.get('step_id') and not step.get('id') and not step.get('tool_id'):
                    result.fail(f"步骤 {i} 缺少 step_id/id/tool_id")
    
    return result


async def main():
    """主函数"""
    print("="*60)
    print("完整流程测试：验证规划、执行阶段的核心功能")
    print("="*60)
    
    # 获取上传目录
    upload_dir = os.getenv("UPLOAD_DIR", str(project_root / "uploads"))
    Path(upload_dir).mkdir(parents=True, exist_ok=True)
    print(f"📁 上传目录: {upload_dir}")
    
    # 运行测试场景
    results = []
    
    # 场景1：未上传文件 - 完整工作流
    result1 = await test_scenario_1_no_file_full_workflow()
    results.append(result1)
    
    # 场景2：未上传文件 - 部分工作流
    result2 = await test_scenario_2_no_file_partial_workflow()
    results.append(result2)
    
    # 场景3：已上传文件 - 完整工作流
    result3 = await test_scenario_3_with_file_full_workflow(upload_dir)
    results.append(result3)
    
    # 前端契约验证
    all_events = []
    for r in results:
        all_events.extend(r.events)
    result4 = verify_frontend_contract(all_events)
    results.append(result4)
    
    # 汇总结果
    print(f"\n{'='*60}")
    print("测试结果汇总")
    print(f"{'='*60}")
    
    all_passed = True
    for result in results:
        status = "✅ PASS" if result.passed else "❌ FAIL"
        print(f"\n{status} {result.scenario_name}")
        if result.errors:
            print(f"  错误 ({len(result.errors)}):")
            for error in result.errors:
                print(f"    - {error}")
        if result.warnings:
            print(f"  警告 ({len(result.warnings)}):")
            for warning in result.warnings:
                print(f"    - {warning}")
        
        if not result.passed:
            all_passed = False
    
    print(f"\n{'='*60}")
    if all_passed:
        print("✅ 所有测试通过！")
        return 0
    else:
        print("❌ 部分测试失败，请检查上述错误")
        return 1


if __name__ == "__main__":
    exit_code = asyncio.run(main())
    sys.exit(exit_code)

