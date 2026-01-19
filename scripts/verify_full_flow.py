#!/usr/bin/env python3
"""
集成测试脚本：验证后端逻辑是否符合前端预期

模拟前端API调用，验证JSON响应结构。
"""
import sys
import os
import asyncio
import json
import re
from pathlib import Path
from typing import Dict, Any, List, Optional

# 添加项目根目录到路径
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from gibh_agent.main import GIBHAgent
from gibh_agent.core.orchestrator import AgentOrchestrator


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


def parse_sse_stream(stream_text: str) -> List[Dict[str, Any]]:
    """
    解析 SSE 流文本，提取事件
    
    Args:
        stream_text: SSE 格式的文本流
        
    Returns:
        事件列表，每个事件包含 type 和 data
    """
    events = []
    lines = stream_text.split('\n')
    
    current_event_type = None
    current_data_lines = []
    
    for line in lines:
        line = line.strip()
        if not line:
            # 空行表示事件结束
            if current_event_type and current_data_lines:
                try:
                    data_json = '\n'.join(current_data_lines)
                    data = json.loads(data_json)
                    events.append({
                        "type": current_event_type,
                        "data": data
                    })
                except json.JSONDecodeError as e:
                    print(f"  ⚠️  JSON 解析错误: {e}, 数据: {data_json[:100]}")
            current_event_type = None
            current_data_lines = []
            continue
        
        if line.startswith('event: '):
            current_event_type = line[7:].strip()
        elif line.startswith('data: '):
            current_data_lines.append(line[6:])
        elif current_data_lines:
            # 多行数据
            current_data_lines.append(line)
    
    # 处理最后一个事件
    if current_event_type and current_data_lines:
        try:
            data_json = '\n'.join(current_data_lines)
            data = json.loads(data_json)
            events.append({
                "type": current_event_type,
                "data": data
            })
        except json.JSONDecodeError as e:
            print(f"  ⚠️  JSON 解析错误: {e}, 数据: {data_json[:100]}")
    
    return events


async def test_scenario_1_plan_first(agent: GIBHAgent, upload_dir: str) -> TestResult:
    """
    场景1：Plan-First（无文件）
    
    输入: query="I want PCA", files=[]
    验证: workflow 事件，template_mode=True，steps 不为空
    """
    result = TestResult("Scenario 1: Plan-First (No File)")
    
    print(f"\n{'='*60}")
    print(f"测试场景: {result.scenario_name}")
    print(f"{'='*60}")
    
    try:
        # 创建编排器
        orchestrator = AgentOrchestrator(agent, upload_dir=upload_dir)
        
        # 调用 stream_process
        query = "I want PCA"
        files = []
        
        print(f"\n📤 输入:")
        print(f"  Query: {query}")
        print(f"  Files: {files}")
        
        # 收集所有 SSE 事件
        stream_chunks = []
        async for chunk in orchestrator.stream_process(
            query=query,
            files=files,
            session_id="test-session-1",
            user_id="test-user"
        ):
            stream_chunks.append(chunk)
        
        # 解析 SSE 流
        stream_text = ''.join(stream_chunks)
        events = parse_sse_stream(stream_text)
        
        print(f"\n📥 收到 {len(events)} 个事件:")
        for i, event in enumerate(events):
            print(f"  [{i+1}] {event['type']}: {json.dumps(event['data'], ensure_ascii=False)[:100]}...")
            result.add_event(event['type'], event['data'])
        
        # 验证 workflow 事件
        workflow_events = [e for e in events if e['type'] == 'workflow']
        if not workflow_events:
            result.fail("未找到 'workflow' 事件")
        else:
            result.success(f"找到 {len(workflow_events)} 个 'workflow' 事件")
            workflow_data = workflow_events[0]['data']
            
            # 验证 workflow_config 结构
            workflow_config = workflow_data.get('workflow_config') or workflow_data.get('workflow_data') or workflow_data
            
            # 检查 steps
            steps = workflow_config.get('steps') or workflow_config.get('workflow_data', {}).get('steps', [])
            if not steps or len(steps) == 0:
                result.fail("workflow_config.steps 为空")
            else:
                result.success(f"workflow_config.steps 包含 {len(steps)} 个步骤")
                
                # 验证每个步骤都有 step_id
                for i, step in enumerate(steps):
                    step_id = step.get('step_id') or step.get('id') or step.get('tool_id')
                    if not step_id:
                        result.fail(f"步骤 {i} 缺少 step_id/id/tool_id: {step}")
            
            # 验证 template_mode
            template_mode = workflow_config.get('template_mode') or workflow_data.get('template_mode')
            if template_mode is not True:
                result.fail(f"template_mode 应为 True，实际为: {template_mode}")
            else:
                result.success("template_mode = True")
            
            # 验证 file_path 占位符
            has_placeholder = False
            for step in steps:
                params = step.get('params', {})
                file_path = params.get('file_path') or params.get('adata_path')
                if file_path and ('<待上传数据>' in str(file_path) or '<PENDING_UPLOAD>' in str(file_path)):
                    has_placeholder = True
                    break
            
            if not has_placeholder:
                result.warn("未找到 file_path 占位符（<待上传数据> 或 <PENDING_UPLOAD>）")
            else:
                result.success("找到 file_path 占位符")
        
        # 验证 diagnosis 事件
        diagnosis_events = [e for e in events if e['type'] == 'diagnosis']
        if diagnosis_events:
            diagnosis_data = diagnosis_events[0]['data']
            diagnosis_text = diagnosis_data.get('message') or str(diagnosis_data)
            if '方案已生成' in diagnosis_text or 'Template Ready' in diagnosis_text or 'template_ready' in str(diagnosis_data):
                result.success("诊断报告包含 'Template Ready' 信息")
            else:
                result.warn(f"诊断报告可能不包含模板信息: {diagnosis_text[:100]}")
        
        # 验证 result 事件
        result_events = [e for e in events if e['type'] == 'result']
        if result_events:
            result_data = result_events[0]['data']
            
            # 检查结构嵌套
            if 'workflow_data' in result_data and 'workflow_data' in result_data.get('workflow_data', {}):
                result.fail("检测到嵌套的 workflow_data.workflow_data（结构错误）")
            
            # 检查 diagnosis_report 和 workflow_config
            has_diagnosis = 'diagnosis_report' in result_data or 'report_data' in result_data
            has_workflow = 'workflow_config' in result_data or 'report_data' in result_data
            
            if has_diagnosis and has_workflow:
                result.success("result 事件包含 diagnosis_report 和 workflow_config")
            else:
                result.warn(f"result 事件结构: diagnosis={has_diagnosis}, workflow={has_workflow}")
        
        # 验证 done 事件
        done_events = [e for e in events if e['type'] == 'done']
        if not done_events:
            result.warn("未找到 'done' 事件")
        else:
            result.success("找到 'done' 事件")
        
    except Exception as e:
        result.fail(f"测试执行失败: {e}")
        import traceback
        print(f"  详细错误:\n{traceback.format_exc()}")
    
    return result


async def test_scenario_2_execution(agent: GIBHAgent, upload_dir: str) -> TestResult:
    """
    场景2：执行模式（有文件）
    
    输入: query="Analyze this", files=["/mock/path/cow_diet.csv"]
    验证: workflow 事件，template_mode=False，file_path 为真实路径
    """
    result = TestResult("Scenario 2: Execution (With File)")
    
    print(f"\n{'='*60}")
    print(f"测试场景: {result.scenario_name}")
    print(f"{'='*60}")
    
    try:
        # 创建编排器
        orchestrator = AgentOrchestrator(agent, upload_dir=upload_dir)
        
        # 创建模拟文件
        mock_file_path = Path(upload_dir) / "test_cow_diet.csv"
        mock_file_path.parent.mkdir(parents=True, exist_ok=True)
        mock_file_path.write_text("sample,group,metabolite1,metabolite2\n1,A,1.0,2.0\n2,B,1.5,2.5\n")
        
        # 调用 stream_process
        query = "Analyze this"
        files = [{"name": "test_cow_diet.csv", "path": str(mock_file_path)}]
        
        print(f"\n📤 输入:")
        print(f"  Query: {query}")
        print(f"  Files: {files}")
        
        # 收集所有 SSE 事件
        stream_chunks = []
        async for chunk in orchestrator.stream_process(
            query=query,
            files=files,
            session_id="test-session-2",
            user_id="test-user"
        ):
            stream_chunks.append(chunk)
        
        # 解析 SSE 流
        stream_text = ''.join(stream_chunks)
        events = parse_sse_stream(stream_text)
        
        print(f"\n📥 收到 {len(events)} 个事件:")
        for i, event in enumerate(events):
            print(f"  [{i+1}] {event['type']}: {json.dumps(event['data'], ensure_ascii=False)[:100]}...")
            result.add_event(event['type'], event['data'])
        
        # 验证 workflow 事件
        workflow_events = [e for e in events if e['type'] == 'workflow']
        if not workflow_events:
            result.fail("未找到 'workflow' 事件")
        else:
            result.success(f"找到 {len(workflow_events)} 个 'workflow' 事件")
            workflow_data = workflow_events[0]['data']
            
            # 验证 workflow_config 结构
            workflow_config = workflow_data.get('workflow_config') or workflow_data.get('workflow_data') or workflow_data
            
            # 检查 steps
            steps = workflow_config.get('steps') or workflow_config.get('workflow_data', {}).get('steps', [])
            if not steps or len(steps) == 0:
                result.fail("workflow_config.steps 为空")
            else:
                result.success(f"workflow_config.steps 包含 {len(steps)} 个步骤")
            
            # 验证 template_mode 应为 False 或不存在
            template_mode = workflow_config.get('template_mode') or workflow_data.get('template_mode')
            if template_mode is True:
                result.fail(f"template_mode 应为 False 或不存在，实际为: {template_mode}")
            else:
                result.success(f"template_mode = {template_mode} (正确)")
            
            # 验证 file_path 为真实路径（不是占位符）
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
        
        # 验证 result 事件
        result_events = [e for e in events if e['type'] == 'result']
        if not result_events:
            result.warn("未找到 'result' 事件")
        else:
            result.success("找到 'result' 事件")
        
        # 清理模拟文件
        if mock_file_path.exists():
            mock_file_path.unlink()
        
    except Exception as e:
        result.fail(f"测试执行失败: {e}")
        import traceback
        print(f"  详细错误:\n{traceback.format_exc()}")
    
    return result


async def test_scenario_3_frontend_contract(agent: GIBHAgent, upload_dir: str) -> TestResult:
    """
    场景3：前端契约检查
    
    检查 JSON 结构是否符合前端期望
    """
    result = TestResult("Scenario 3: Frontend Contract Check")
    
    print(f"\n{'='*60}")
    print(f"测试场景: {result.scenario_name}")
    print(f"{'='*60}")
    
    try:
        # 使用场景1的结果进行验证
        scenario1_result = await test_scenario_1_plan_first(agent, upload_dir)
        
        if not scenario1_result.events:
            result.fail("场景1未产生事件，无法进行契约检查")
            return result
        
        # 检查所有 workflow 事件的结构
        workflow_events = [e for e in scenario1_result.events if e['type'] == 'workflow']
        
        for event in workflow_events:
            data = event['data']
            
            # 检查嵌套结构
            if 'workflow_data' in data and isinstance(data.get('workflow_data'), dict):
                if 'workflow_data' in data['workflow_data']:
                    result.fail("检测到 workflow_data.workflow_data 嵌套结构（错误）")
                else:
                    result.success("workflow_data 结构正确（无嵌套）")
            
            # 检查步骤结构
            workflow_config = data.get('workflow_config') or data.get('workflow_data') or data
            steps = workflow_config.get('steps') or workflow_config.get('workflow_data', {}).get('steps', [])
            
            for i, step in enumerate(steps):
                # 检查必需字段
                if not step.get('step_id') and not step.get('id') and not step.get('tool_id'):
                    result.fail(f"步骤 {i} 缺少 step_id/id/tool_id: {step}")
                
                # 检查 params 结构
                if 'params' not in step:
                    result.warn(f"步骤 {i} 缺少 params 字段")
        
        # 检查 result 事件结构
        result_events = [e for e in scenario1_result.events if e['type'] == 'result']
        for event in result_events:
            data = event['data']
            
            # 检查是否有 diagnosis_report 和 workflow_config
            has_diagnosis = 'diagnosis_report' in data or 'report_data' in data
            has_workflow = 'workflow_config' in data or 'report_data' in data
            
            if not has_diagnosis:
                result.warn("result 事件缺少 diagnosis_report")
            if not has_workflow:
                result.warn("result 事件缺少 workflow_config")
            
            if has_diagnosis and has_workflow:
                result.success("result 事件结构完整")
        
    except Exception as e:
        result.fail(f"测试执行失败: {e}")
        import traceback
        print(f"  详细错误:\n{traceback.format_exc()}")
    
    return result


async def main():
    """主函数"""
    print("="*60)
    print("集成测试：验证后端逻辑是否符合前端预期")
    print("="*60)
    
    # 初始化 Agent
    config_path = "config/settings.yaml"
    if not os.path.exists(config_path):
        print(f"❌ 配置文件不存在: {config_path}")
        print("   使用默认配置...")
        config_path = None
    
    try:
        agent = GIBHAgent(config_path) if config_path else GIBHAgent()
        print("✅ GIBHAgent 初始化成功")
    except Exception as e:
        print(f"❌ GIBHAgent 初始化失败: {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    # 获取上传目录（使用项目内的目录避免权限问题）
    upload_dir = os.getenv("UPLOAD_DIR", str(project_root / "uploads"))
    print(f"📁 上传目录: {upload_dir}")
    # 确保目录存在
    Path(upload_dir).mkdir(parents=True, exist_ok=True)
    
    # 运行测试场景
    results = []
    
    # 场景1：Plan-First
    result1 = await test_scenario_1_plan_first(agent, upload_dir)
    results.append(result1)
    
    # 场景2：Execution
    result2 = await test_scenario_2_execution(agent, upload_dir)
    results.append(result2)
    
    # 场景3：Frontend Contract
    result3 = await test_scenario_3_frontend_contract(agent, upload_dir)
    results.append(result3)
    
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

