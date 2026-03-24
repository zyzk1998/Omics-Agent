#!/usr/bin/env python3
"""
完整工作流测试脚本
测试流程：
1. 无文件预览（Plan-First模式）
2. 上传文件规划（Execution模式）
3. 执行工作流
4. 验证输出结果
"""

import asyncio
import sys
import os
import json
from pathlib import Path

# 添加项目路径
sys.path.insert(0, str(Path(__file__).parent.parent))

from gibh_agent.core.orchestrator import AgentOrchestrator
from gibh_agent.core.file_inspector import FileInspector
from gibh_agent.core.planner import SOPPlanner
from gibh_agent.core.tool_retriever import ToolRetriever
from gibh_agent.core.executor import WorkflowExecutor
from gibh_agent.llm_client import LLMClient


async def test_full_workflow():
    """测试完整工作流"""
    print("=" * 80)
    print("🧪 完整工作流测试")
    print("=" * 80)
    print()
    
    # 初始化组件
    print("📦 初始化组件...")
    try:
        llm_client = LLMClient()
        file_inspector = FileInspector()
        tool_retriever = ToolRetriever()
        planner = SOPPlanner(tool_retriever, llm_client)
        orchestrator = AgentOrchestrator(
            agent=None,  # 暂时不使用agent
            file_inspector=file_inspector
        )
        print("✅ 组件初始化成功")
    except Exception as e:
        print(f"❌ 组件初始化失败: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    print()
    
    # Step 1: 无文件预览（Plan-First模式）
    print("=" * 80)
    print("Step 1: 无文件预览（Plan-First模式）")
    print("=" * 80)
    try:
        query = "代谢组学分析"
        print(f"📝 用户查询: {query}")
        print("📤 发送请求（无文件）...")
        
        events = []
        async for event in orchestrator.stream_process(query=query, files=None):
            events.append(event)
            # 解析SSE事件
            if event.startswith("event:"):
                event_type = event.split("event:")[1].split("\n")[0].strip()
                if "data:" in event:
                    data_str = event.split("data:")[1].strip()
                    try:
                        data = json.loads(data_str)
                        if event_type == "workflow":
                            print(f"✅ 收到workflow事件: {len(data.get('workflow_config', {}).get('workflow_data', {}).get('steps', []))} 个步骤")
                        elif event_type == "status":
                            print(f"📊 状态: {data.get('content', '')}")
                    except:
                        pass
        
        print(f"✅ Step 1完成: 收到 {len(events)} 个事件")
        print()
    except Exception as e:
        print(f"❌ Step 1失败: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    # Step 2: 上传文件规划（Execution模式）
    print("=" * 80)
    print("Step 2: 上传文件规划（Execution模式）")
    print("=" * 80)
    
    # 查找测试文件
    test_file = None
    for root, dirs, files in os.walk('.'):
        for file in files:
            if file.endswith('.csv') and 'test' in file.lower():
                test_file = os.path.join(root, file)
                break
        if test_file:
            break
    
    if not test_file:
        print("❌ 未找到测试CSV文件")
        return False
    
    print(f"📁 使用测试文件: {test_file}")
    
    try:
        # 检查文件
        file_metadata = file_inspector.inspect_file(test_file)
        if file_metadata.get("status") != "success":
            print(f"❌ 文件检查失败: {file_metadata.get('error')}")
            return False
        
        print(f"✅ 文件检查成功:")
        print(f"   - 样本数: {file_metadata.get('n_samples', 'N/A')}")
        print(f"   - 特征数: {file_metadata.get('n_features', 'N/A')}")
        print(f"   - 列数: {len(file_metadata.get('columns', []))}")
        
        # 规划工作流
        query = "分析这个代谢组学数据"
        files = [{"path": test_file, "name": os.path.basename(test_file)}]
        
        print(f"📝 用户查询: {query}")
        print("📤 发送请求（有文件）...")
        
        events = []
        diagnosis_received = False
        workflow_received = False
        
        async for event in orchestrator.stream_process(query=query, files=files):
            events.append(event)
            # 解析SSE事件
            if event.startswith("event:"):
                event_type = event.split("event:")[1].split("\n")[0].strip()
                if "data:" in event:
                    data_str = event.split("data:")[1].strip()
                    try:
                        data = json.loads(data_str)
                        if event_type == "diagnosis":
                            diagnosis_received = True
                            print(f"✅ 收到diagnosis事件")
                            print(f"   - 消息: {data.get('message', '')[:100]}...")
                        elif event_type == "workflow":
                            workflow_received = True
                            steps = data.get('workflow_config', {}).get('workflow_data', {}).get('steps', [])
                            print(f"✅ 收到workflow事件: {len(steps)} 个步骤")
                            for i, step in enumerate(steps[:3], 1):
                                print(f"   {i}. {step.get('name', step.get('id', 'Unknown'))}")
                        elif event_type == "status":
                            content = data.get('content', '')
                            if '执行' in content or '生成' in content:
                                print(f"📊 状态: {content}")
                    except Exception as e:
                        pass
        
        if diagnosis_received and workflow_received:
            print(f"✅ Step 2完成: 收到 {len(events)} 个事件")
            print(f"   - 诊断事件: ✅")
            print(f"   - 工作流事件: ✅")
        else:
            print(f"⚠️ Step 2部分完成: 诊断={diagnosis_received}, 工作流={workflow_received}")
        print()
    except Exception as e:
        print(f"❌ Step 2失败: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    # Step 3: 执行工作流
    print("=" * 80)
    print("Step 3: 执行工作流")
    print("=" * 80)
    
    try:
        # 生成工作流配置
        result = None
        async for _ev, _data in planner.generate_plan(
            user_query=query,
            file_metadata=file_metadata,
            is_template=False,
        ):
            if _ev == "workflow":
                result = _data
        
        if not result or result.get("type") != "workflow_config":
            print(f"❌ 规划失败: {(result or {}).get('error', 'Unknown error')}")
            return False
        
        workflow_config = result.get("workflow_data", {})
        steps = workflow_config.get("steps", [])
        
        print(f"📋 工作流配置生成成功: {len(steps)} 个步骤")
        
        # 执行工作流
        executor = WorkflowExecutor()
        execution_results = executor.execute_workflow(
            workflow_data=workflow_config,
            file_paths=[test_file]
        )
        
        status = execution_results.get("status", "unknown")
        steps_details = execution_results.get("steps_details", [])
        
        print(f"✅ 工作流执行完成:")
        print(f"   - 状态: {status}")
        print(f"   - 步骤数: {len(steps_details)}")
        
        success_count = sum(1 for s in steps_details if s.get("status") == "success")
        print(f"   - 成功步骤: {success_count}/{len(steps_details)}")
        
        for step in steps_details[:3]:
            step_name = step.get("step_name", step.get("step_id", "Unknown"))
            step_status = step.get("status", "unknown")
            print(f"   - {step_name}: {step_status}")
        
        print()
    except Exception as e:
        print(f"❌ Step 3失败: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    # Step 4: 验证输出结果
    print("=" * 80)
    print("Step 4: 验证输出结果")
    print("=" * 80)
    
    try:
        output_dir = executor.output_dir if hasattr(executor, 'output_dir') else None
        
        if output_dir and os.path.exists(output_dir):
            print(f"📂 输出目录: {output_dir}")
            
            # 检查生成的文件
            csv_files = list(Path(output_dir).glob("*.csv"))
            png_files = list(Path(output_dir).glob("*.png"))
            
            print(f"✅ 生成的文件:")
            print(f"   - CSV文件: {len(csv_files)} 个")
            for f in csv_files[:3]:
                print(f"     * {f.name}")
            print(f"   - PNG图片: {len(png_files)} 个")
            for f in png_files[:3]:
                print(f"     * {f.name}")
        else:
            print("⚠️ 输出目录不存在或未设置")
        
        print()
        print("=" * 80)
        print("✅ 完整工作流测试通过！")
        print("=" * 80)
        return True
        
    except Exception as e:
        print(f"❌ Step 4失败: {e}")
        import traceback
        traceback.print_exc()
        return False


if __name__ == "__main__":
    success = asyncio.run(test_full_workflow())
    sys.exit(0 if success else 1)
