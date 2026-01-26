#!/usr/bin/env python3
"""
前端模拟器 - 验证工作流规划和执行逻辑
模拟前端的两阶段交互：规划 -> 执行
"""
import os
import sys
import asyncio
import json
from pathlib import Path
from typing import Dict, Any, List, Optional

# 添加项目路径
sys.path.insert(0, str(Path(__file__).parent.parent))

from gibh_agent import create_agent
from gibh_agent.core.orchestrator import AgentOrchestrator
from gibh_agent.core.file_inspector import FileInspector

# 配置
TEST_DATA_DIR = Path(__file__).parent.parent / "test_data"
TEST_CSV_PATH = TEST_DATA_DIR / "cow_diet.csv"
UPLOAD_DIR = Path(os.getenv("UPLOAD_DIR", "/app/uploads"))


class FrontendSimulator:
    """前端模拟器"""
    
    def __init__(self):
        self.agent = create_agent()
        # 使用测试数据目录作为上传目录（避免权限问题）
        test_upload_dir = str(TEST_DATA_DIR)
        self.orchestrator = AgentOrchestrator(self.agent, upload_dir=test_upload_dir)
        self.file_inspector = FileInspector(test_upload_dir)
        self.test_file_path = str(TEST_CSV_PATH)
        
    async def simulate_plan_request(self) -> Dict[str, Any]:
        """
        模拟规划请求
        
        Returns:
            工作流配置字典
        """
        print("=" * 80)
        print("📋 Step 1: 模拟规划请求")
        print("=" * 80)
        
        # 准备文件
        if not Path(self.test_file_path).exists():
            print(f"❌ 测试文件不存在: {self.test_file_path}")
            return None
        
        # 检查文件
        file_metadata = self.file_inspector.inspect_file(self.test_file_path)
        print(f"✅ 文件检查完成: {file_metadata.get('file_type', 'unknown')}")
        
        # 准备文件列表
        files = [{
            "name": Path(self.test_file_path).name,
            "path": self.test_file_path
        }]
        
        # 发送规划请求
        query = "Analyze this"
        print(f"\n💬 发送规划请求: '{query}'")
        print(f"   文件: {files[0]['name']}")
        
        # 收集事件
        events = []
        workflow_config = None
        result_data = None
        has_step_result = False
        current_event_type = None
        
        async for event_str in self.orchestrator.stream_process(
            query=query,
            files=files,
            history=[],
            session_id="test-session",
            user_id="test_user"
        ):
            # 解析 SSE 事件（格式: "event: type\ndata: {...}\n\n"）
            # event_str 可能包含多行，需要按行解析
            lines = event_str.strip().split('\n')
            for line in lines:
                line = line.strip()
                if not line:
                    continue
                if line.startswith("event:"):
                    current_event_type = line.split(":", 1)[1].strip()
                elif line.startswith("data:"):
                    data_str = line.split(":", 1)[1].strip()
                    try:
                        data = json.loads(data_str)
                        if current_event_type:
                            events.append({"type": current_event_type, "data": data})
                            
                            if current_event_type == "workflow":
                                workflow_config = data.get("workflow_config") or data.get("workflow_data")
                                print(f"   ✅ 收到 workflow 事件")
                                print(f"      template_mode: {data.get('template_mode', 'N/A')}")
                                if workflow_config:
                                    steps = workflow_config.get("steps", [])
                                    print(f"      步骤数: {len(steps)}")
                            
                            elif current_event_type == "result":
                                result_data = data
                                print(f"   ✅ 收到 result 事件")
                                print(f"      键: {list(data.keys())}")
                                # 也从 result 事件中提取 workflow_config
                                if not workflow_config:
                                    workflow_config = data.get("workflow_config") or data.get("workflow_data")
                            
                            elif current_event_type == "status":
                                content = data.get("content", "")
                                state = data.get("state", "")
                                print(f"   📊 [{state}] {content}")
                            
                            elif current_event_type == "step_result":
                                has_step_result = True
                                print(f"   ⚠️ 收到 step_result 事件（不应该在规划阶段出现）")
                            
                            elif current_event_type == "done":
                                print(f"   🏁 流式传输完成")
                    except json.JSONDecodeError as e:
                        print(f"   ⚠️ JSON 解析失败: {e} - {data_str[:100]}")
        
        # 验证结果
        print("\n" + "=" * 80)
        print("🔍 验证规划结果")
        print("=" * 80)
        
        assertions = []
        
        # 断言 1: 收到了 workflow 事件
        if workflow_config:
            assertions.append(("收到 workflow 事件", True))
            print("✅ 断言 1: 收到 workflow 事件")
        else:
            assertions.append(("收到 workflow 事件", False))
            print("❌ 断言 1: 未收到 workflow 事件")
        
        # 断言 2: template_mode == False
        if result_data:
            template_mode = result_data.get("template_mode", True)
            if template_mode == False:
                assertions.append(("template_mode == False", True))
                print("✅ 断言 2: template_mode == False")
            else:
                assertions.append(("template_mode == False", False))
                print(f"❌ 断言 2: template_mode == {template_mode} (期望 False)")
        else:
            assertions.append(("template_mode == False", False))
            print("❌ 断言 2: 未收到 result 事件，无法验证 template_mode")
        
        # 断言 3: steps 列表不为空
        if workflow_config:
            steps = workflow_config.get("steps", [])
            if steps and len(steps) > 0:
                assertions.append(("steps 列表不为空", True))
                print(f"✅ 断言 3: steps 列表不为空 ({len(steps)} 个步骤)")
            else:
                assertions.append(("steps 列表不为空", False))
                print(f"❌ 断言 3: steps 列表为空")
        else:
            assertions.append(("steps 列表不为空", False))
            print("❌ 断言 3: 无 workflow_config，无法验证 steps")
        
        # 断言 4: metabolomics_plsda 有 group_column 参数
        if workflow_config:
            steps = workflow_config.get("steps", [])
            plsda_step = next((s for s in steps if s.get("id") == "metabolomics_plsda"), None)
            if plsda_step:
                params = plsda_step.get("params", {})
                group_column = params.get("group_column")
                if group_column:
                    assertions.append(("metabolomics_plsda 有 group_column", True))
                    print(f"✅ 断言 4: metabolomics_plsda 有 group_column = '{group_column}'")
                else:
                    assertions.append(("metabolomics_plsda 有 group_column", False))
                    print(f"❌ 断言 4: metabolomics_plsda 没有 group_column")
                    print(f"   参数: {list(params.keys())}")
            else:
                assertions.append(("metabolomics_plsda 有 group_column", False))
                print("❌ 断言 4: 未找到 metabolomics_plsda 步骤")
        else:
            assertions.append(("metabolomics_plsda 有 group_column", False))
            print("❌ 断言 4: 无 workflow_config，无法验证")
        
        # 断言 5: 没有 step_result 事件（证明停止了）
        if not has_step_result:
            assertions.append(("没有 step_result 事件", True))
            print("✅ 断言 5: 没有 step_result 事件（正确停止在规划阶段）")
        else:
            assertions.append(("没有 step_result 事件", False))
            print("❌ 断言 5: 有 step_result 事件（不应该在规划阶段出现）")
        
        # 总结
        print("\n" + "=" * 80)
        passed = all(assertion[1] for assertion in assertions)
        if passed:
            print("✅ 所有断言通过！")
        else:
            print("❌ 部分断言失败：")
            for name, result in assertions:
                if not result:
                    print(f"   - {name}")
        
        return workflow_config if passed else None
    
    async def simulate_execute_request(self, workflow_config: Dict[str, Any]) -> bool:
        """
        模拟执行请求
        
        Args:
            workflow_config: 工作流配置
            
        Returns:
            是否成功
        """
        print("\n" + "=" * 80)
        print("🚀 Step 2: 模拟执行请求")
        print("=" * 80)
        
        if not workflow_config:
            print("❌ 工作流配置为空，无法执行")
            return False
        
        steps = workflow_config.get("steps", [])
        print(f"   工作流步骤数: {len(steps)}")
        
        # 发送执行请求
        print(f"\n💬 发送执行请求（带 workflow_data）")
        
        # 收集事件
        events = []
        steps_details = []
        has_step_result = False
        
        current_event_type = None
        async for event_str in self.orchestrator.stream_process(
            query="",
            files=[],
            history=[],
            session_id="test-session",
            user_id="test_user",
            workflow_data=workflow_config
        ):
            # 解析 SSE 事件（格式: "event: type\ndata: {...}\n\n"）
            lines = event_str.strip().split('\n')
            for line in lines:
                line = line.strip()
                if not line:
                    continue
                if line.startswith("event:"):
                    current_event_type = line.split(":", 1)[1].strip()
                elif line.startswith("data:"):
                    data_str = line.split(":", 1)[1].strip()
                    try:
                        data = json.loads(data_str)
                        if current_event_type:
                            events.append({"type": current_event_type, "data": data})
                            
                            if current_event_type == "status":
                                content = data.get("content", "")
                                state = data.get("state", "")
                                print(f"   📊 [{state}] {content}")
                            
                            elif current_event_type == "step_result":
                                has_step_result = True
                                step_id = data.get("step_id", "unknown")
                                step_name = data.get("step_name", step_id)
                                status = data.get("status", "unknown")
                                print(f"   🔧 [{status}] {step_name}")
                                
                                steps_details.append({
                                    "step_id": step_id,
                                    "step_name": step_name,
                                    "status": status,
                                    "result": data.get("result"),
                                    "error": data.get("error")
                                })
                            
                            elif current_event_type == "result":
                                result_data = data
                                report_data = result_data.get("report_data", {})
                                if report_data.get("steps_details"):
                                    steps_details = report_data.get("steps_details", [])
                                print(f"   ✅ 收到最终结果")
                                print(f"      步骤详情数: {len(steps_details)}")
                                print(f"      结果键: {list(data.keys())}")
                            
                            elif current_event_type == "done":
                                print(f"   🏁 流式传输完成")
                            
                            elif current_event_type == "error":
                                error_msg = data.get("error", "未知错误")
                                print(f"   ❌ 错误: {error_msg}")
                            
                    except json.JSONDecodeError as e:
                        print(f"   ⚠️ JSON 解析失败: {e} - {data_str[:100]}")
        
        # 验证结果
        print("\n" + "=" * 80)
        print("🔍 验证执行结果")
        print("=" * 80)
        
        assertions = []
        
        # 断言 1: 有 step_result 事件
        if has_step_result or len(steps_details) > 0:
            assertions.append(("有 step_result 事件", True))
            print(f"✅ 断言 1: 有 step_result 事件 ({len(steps_details)} 个步骤)")
        else:
            assertions.append(("有 step_result 事件", False))
            print("❌ 断言 1: 没有 step_result 事件")
        
        # 断言 2: metabolomics_plsda 返回 success
        plsda_result = next((s for s in steps_details if s.get("step_id") == "metabolomics_plsda"), None)
        if plsda_result:
            status = plsda_result.get("status")
            if status == "success":
                assertions.append(("metabolomics_plsda 返回 success", True))
                print(f"✅ 断言 2: metabolomics_plsda 返回 success")
            else:
                assertions.append(("metabolomics_plsda 返回 success", False))
                error = plsda_result.get("error", "未知错误")
                print(f"❌ 断言 2: metabolomics_plsda 返回 {status}")
                print(f"   错误: {error}")
        else:
            assertions.append(("metabolomics_plsda 返回 success", False))
            print("❌ 断言 2: 未找到 metabolomics_plsda 执行结果")
        
        # 断言 3: visualize_volcano 返回 success
        volcano_result = next((s for s in steps_details if s.get("step_id") == "visualize_volcano"), None)
        if volcano_result:
            status = volcano_result.get("status")
            if status == "success":
                assertions.append(("visualize_volcano 返回 success", True))
                print(f"✅ 断言 3: visualize_volcano 返回 success")
            else:
                assertions.append(("visualize_volcano 返回 success", False))
                error = volcano_result.get("error", "未知错误")
                print(f"❌ 断言 3: visualize_volcano 返回 {status}")
                print(f"   错误: {error}")
        else:
            assertions.append(("visualize_volcano 返回 success", False))
            print("❌ 断言 3: 未找到 visualize_volcano 执行结果")
        
        # 总结
        print("\n" + "=" * 80)
        passed = all(assertion[1] for assertion in assertions)
        if passed:
            print("✅ 所有断言通过！")
        else:
            print("❌ 部分断言失败：")
            for name, result in assertions:
                if not result:
                    print(f"   - {name}")
        
        return passed


async def main():
    """主函数"""
    print("=" * 80)
    print("🚀 前端模拟器 - 验证工作流规划和执行逻辑")
    print("=" * 80)
    
    simulator = FrontendSimulator()
    
    # Step 1: 模拟规划请求
    workflow_config = await simulator.simulate_plan_request()
    
    if not workflow_config:
        print("\n❌ 规划阶段失败，无法继续执行")
        return False
    
    # Step 2: 模拟执行请求
    success = await simulator.simulate_execute_request(workflow_config)
    
    # 最终总结
    print("\n" + "=" * 80)
    if success:
        print("✅ 所有测试通过！")
    else:
        print("❌ 测试失败！")
    print("=" * 80)
    
    return success


if __name__ == "__main__":
    success = asyncio.run(main())
    sys.exit(0 if success else 1)

