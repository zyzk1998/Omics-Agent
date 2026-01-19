#!/usr/bin/env python3
"""
测试 Orchestrator SSE 事件格式

验证 orchestrator 发送的 SSE 事件格式是否符合前端期望。
"""
import sys
import os
import asyncio
import json
import re
from pathlib import Path

# 添加项目根目录到路径
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from gibh_agent.main import GIBHAgent
from gibh_agent.core.orchestrator import AgentOrchestrator


class MockLLMClient:
    """Mock LLM 客户端"""
    def __init__(self):
        self.call_count = 0
    
    async def achat(self, messages, **kwargs):
        from openai.types.chat import ChatCompletion
        from openai.types.chat.chat_completion import ChatCompletionMessage
        
        self.call_count += 1
        user_message = messages[-1]["content"] if messages else ""
        system_message = messages[0]["content"] if messages else ""
        
        # 根据消息类型返回不同响应
        if "Intent Classifier" in system_message and "domain_name" in system_message:
            response_text = '{"domain_name": "Metabolomics", "target_steps": []}'
        elif "Intent Analyzer" in system_message:
            response_text = '[]'
        elif "Query Rewriter" in system_message:
            response_text = user_message
        elif "Intent Classifier" in system_message and "PLAN" in system_message:
            response_text = "planning"
        elif "Clarification Assistant" in system_message:
            response_text = "OK"
        elif "Workflow Reflector" in system_message:
            response_text = '{"is_valid": true}'
        else:
            response_text = "OK"
        
        # 创建简单的 ChatCompletion 对象
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
        
        return SimpleCompletion(response_text)


def parse_sse_stream(stream_text: str) -> list:
    """解析 SSE 流"""
    events = []
    lines = stream_text.split('\n')
    
    current_event_type = None
    current_data_lines = []
    
    for line in lines:
        line = line.strip()
        if not line:
            if current_event_type and current_data_lines:
                try:
                    data_json = '\n'.join(current_data_lines)
                    data = json.loads(data_json)
                    events.append({"type": current_event_type, "data": data})
                except json.JSONDecodeError:
                    pass
            current_event_type = None
            current_data_lines = []
            continue
        
        if line.startswith('event: '):
            current_event_type = line[7:].strip()
        elif line.startswith('data: '):
            current_data_lines.append(line[6:])
        elif current_data_lines:
            current_data_lines.append(line)
    
    if current_event_type and current_data_lines:
        try:
            data_json = '\n'.join(current_data_lines)
            data = json.loads(data_json)
            events.append({"type": current_event_type, "data": data})
        except json.JSONDecodeError:
            pass
    
    return events


async def test_orchestrator_sse_no_file():
    """测试：无文件时的 SSE 事件格式"""
    print("\n" + "="*60)
    print("测试：Orchestrator SSE 事件格式（无文件）")
    print("="*60)
    
    try:
        # 创建 Mock Agent
        class MockAgent:
            def __init__(self):
                self.agents = {}
            
            async def process_query(self, query, history=None, uploaded_files=None, **kwargs):
                # 返回模拟结果
                return {
                    "type": "workflow_config",
                    "report_data": {
                        "diagnosis": {
                            "status": "template_ready",
                            "message": "### 📋 分析方案已生成\n\n根据您的需求，我为您规划了以下流程...",
                            "template_mode": True
                        },
                        "workflow": {
                            "workflow_data": {
                                "steps": [
                                    {"id": "inspect_data", "name": "数据检查", "params": {"file_path": "<待上传数据>"}},
                                    {"id": "preprocess_data", "name": "数据预处理", "params": {"file_path": "<待上传数据>"}},
                                    {"id": "pca_analysis", "name": "PCA分析", "params": {"file_path": "<待上传数据>"}}
                                ],
                                "template_mode": True
                            },
                            "template_mode": True
                        }
                    }
                }
        
        mock_agent = MockAgent()
        upload_dir = str(project_root / "uploads")
        Path(upload_dir).mkdir(parents=True, exist_ok=True)
        
        orchestrator = AgentOrchestrator(mock_agent, upload_dir=upload_dir)
        
        # Mock LLM client
        mock_llm = MockLLMClient()
        orchestrator.query_rewriter = type('obj', (object,), {'rewrite': lambda self, q, h: asyncio.coroutine(lambda: q)()})()
        orchestrator.clarifier = type('obj', (object,), {'check_and_clarify': lambda self, q, f, d: asyncio.coroutine(lambda: None)()})()
        
        # 调用 stream_process
        query = "I want PCA"
        files = []
        
        print(f"\n📤 输入:")
        print(f"  Query: {query}")
        print(f"  Files: {files}")
        
        # 收集 SSE 事件
        stream_chunks = []
        async for chunk in orchestrator.stream_process(
            query=query,
            files=files,
            session_id="test-sse-1",
            user_id="test-user"
        ):
            stream_chunks.append(chunk)
        
        # 解析 SSE 流
        stream_text = ''.join(stream_chunks)
        events = parse_sse_stream(stream_text)
        
        print(f"\n📥 收到 {len(events)} 个 SSE 事件:")
        for i, event in enumerate(events):
            event_type = event['type']
            data = event['data']
            print(f"  [{i+1}] {event_type}: {json.dumps(data, ensure_ascii=False)[:100]}...")
        
        # 验证事件格式
        errors = []
        warnings = []
        
        # 检查是否有 workflow 事件
        workflow_events = [e for e in events if e['type'] == 'workflow']
        if not workflow_events:
            errors.append("未找到 'workflow' 事件")
        else:
            print(f"\n✅ 找到 {len(workflow_events)} 个 'workflow' 事件")
            
            workflow_data = workflow_events[0]['data']
            
            # 检查结构
            workflow_config = workflow_data.get('workflow_config') or workflow_data.get('workflow_data') or workflow_data
            
            # 检查 steps
            steps = workflow_config.get('steps') or workflow_config.get('workflow_data', {}).get('steps', [])
            if not steps:
                errors.append("workflow_config.steps 为空")
            else:
                print(f"  ✅ steps 包含 {len(steps)} 个步骤")
            
            # 检查 template_mode
            template_mode = workflow_config.get('template_mode') or workflow_data.get('template_mode')
            if template_mode is not True:
                errors.append(f"template_mode 应为 True，实际为: {template_mode}")
            else:
                print(f"  ✅ template_mode = True")
        
        # 检查是否有 diagnosis 事件
        diagnosis_events = [e for e in events if e['type'] == 'diagnosis']
        if diagnosis_events:
            print(f"✅ 找到 {len(diagnosis_events)} 个 'diagnosis' 事件")
        
        # 检查是否有 result 事件
        result_events = [e for e in events if e['type'] == 'result']
        if result_events:
            print(f"✅ 找到 {len(result_events)} 个 'result' 事件")
            result_data = result_events[0]['data']
            
            # 检查 result 事件结构
            has_diagnosis = 'diagnosis_report' in result_data or 'report_data' in result_data
            has_workflow = 'workflow_config' in result_data or 'report_data' in result_data
            
            if has_diagnosis and has_workflow:
                print(f"  ✅ result 事件包含 diagnosis_report 和 workflow_config")
            else:
                warnings.append(f"result 事件结构: diagnosis={has_diagnosis}, workflow={has_workflow}")
        
        # 检查是否有 done 事件
        done_events = [e for e in events if e['type'] == 'done']
        if not done_events:
            warnings.append("未找到 'done' 事件")
        else:
            print(f"✅ 找到 'done' 事件")
        
        # 输出结果
        print(f"\n{'='*60}")
        if errors:
            print("❌ 发现错误:")
            for error in errors:
                print(f"  - {error}")
            return False
        else:
            print("✅ SSE 事件格式验证通过")
            if warnings:
                print("\n⚠️  警告:")
                for warning in warnings:
                    print(f"  - {warning}")
            return True
        
    except Exception as e:
        print(f"❌ 测试失败: {e}")
        import traceback
        traceback.print_exc()
        return False


async def main():
    """主函数"""
    print("="*60)
    print("测试：Orchestrator SSE 事件格式验证")
    print("="*60)
    
    success = await test_orchestrator_sse_no_file()
    
    print(f"\n{'='*60}")
    if success:
        print("✅ 所有测试通过！")
        return 0
    else:
        print("❌ 测试失败")
        return 1


if __name__ == "__main__":
    exit_code = asyncio.run(main())
    sys.exit(exit_code)


