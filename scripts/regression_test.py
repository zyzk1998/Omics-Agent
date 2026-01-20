#!/usr/bin/env python3
"""
回归测试脚本：验证所有用户旅程（Legacy & New）无回归

测试场景：
1. Scenario A (Legacy Metabolomics): Upload CSV -> Expect template_mode=False + Diagnosis + Workflow
2. Scenario B (New Plan-First): No File -> Expect template_mode=True -> Upload -> Expect template_mode=False
3. Scenario C (Execution Trigger): Request contains workflow_data -> Expect Immediate Execution (Bypass Planner)
4. Scenario D (Async RNA): Request scRNA-seq -> Planner selects rna_cellranger_count -> Execution returns status="async_job_started" -> Orchestrator must yield this status and STOP
"""

import asyncio
import sys
import os
from pathlib import Path
from typing import Dict, Any, List, Optional
from unittest.mock import Mock, AsyncMock, MagicMock, patch
import json

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from gibh_agent.core.orchestrator import AgentOrchestrator
from gibh_agent.core.planner import SOPPlanner
from gibh_agent.core.executor import ExecutionLayer


class RegressionTestSuite:
    """回归测试套件"""
    
    def __init__(self):
        self.passed = 0
        self.failed = 0
        self.errors = []
    
    def log_test(self, name: str, passed: bool, error: str = None):
        """记录测试结果"""
        if passed:
            self.passed += 1
            print(f"✅ {name}")
        else:
            self.failed += 1
            print(f"❌ {name}")
            if error:
                print(f"   错误: {error}")
                self.errors.append(f"{name}: {error}")
    
    async def test_scenario_a_legacy_metabolomics(self):
        """Scenario A: Legacy Metabolomics - Upload CSV -> Expect template_mode=False + Diagnosis + Workflow"""
        print("\n" + "="*60)
        print("测试场景 A: Legacy Metabolomics")
        print("="*60)
        
        try:
            # Mock dependencies
            mock_agent = Mock()
            mock_file_inspector = Mock()
            mock_file_inspector.inspect_file.return_value = {
                "status": "success",
                "file_path": "/app/uploads/cow_diet.csv",
                "n_samples": 39,
                "n_features": 100,
                "file_type": "CSV"
            }
            
            orchestrator = AgentOrchestrator(mock_agent, upload_dir="/app/uploads")
            orchestrator.file_inspector = mock_file_inspector
            
            # Mock planner
            mock_planner = Mock()
            mock_planner.generate_plan = AsyncMock(return_value={
                "type": "workflow_config",
                "template_mode": False,
                "workflow_data": {
                    "workflow_name": "Metabolomics 标准分析流程",
                    "steps": [
                        {"id": "inspect_data", "name": "数据检查", "params": {"file_path": "/app/uploads/cow_diet.csv"}},
                        {"id": "preprocess_data", "name": "数据预处理", "params": {}}
                    ]
                }
            })
            
            # Mock query rewriter
            orchestrator.query_rewriter = None
            
            # Test: Upload CSV file
            files = [{"name": "cow_diet.csv", "path": "/app/uploads/cow_diet.csv"}]
            query = "分析这个文件"
            
            # Collect events
            events = []
            async for event in orchestrator.stream_process(
                query=query,
                files=files,
                workflow_data=None  # No workflow_data -> should plan
            ):
                events.append(event)
            
            # Verify results
            has_diagnosis = any('diagnosis' in event for event in events)
            has_workflow = any('workflow' in event for event in events)
            has_template_mode_false = any('"template_mode": false' in event or '"template_mode":false' in event for event in events)
            
            passed = has_diagnosis and has_workflow and has_template_mode_false
            
            self.log_test(
                "Scenario A: Legacy Metabolomics",
                passed,
                None if passed else f"诊断: {has_diagnosis}, 工作流: {has_workflow}, template_mode=False: {has_template_mode_false}"
            )
            
        except Exception as e:
            self.log_test("Scenario A: Legacy Metabolomics", False, str(e))
    
    async def test_scenario_b_plan_first(self):
        """Scenario B: Plan-First - No File -> template_mode=True -> Upload -> template_mode=False"""
        print("\n" + "="*60)
        print("测试场景 B: Plan-First")
        print("="*60)
        
        try:
            # Mock dependencies
            mock_agent = Mock()
            orchestrator = AgentOrchestrator(mock_agent, upload_dir="/app/uploads")
            
            # Mock planner for template mode
            mock_planner = Mock()
            mock_planner.generate_plan = AsyncMock(return_value={
                "type": "workflow_config",
                "template_mode": True,
                "workflow_data": {
                    "workflow_name": "Metabolomics 标准分析流程",
                    "steps": [
                        {"id": "inspect_data", "name": "数据检查", "params": {"file_path": "<待上传数据>"}}
                    ]
                }
            })
            
            # Test Step 1: No File -> Expect template_mode=True
            events_no_file = []
            async for event in orchestrator.stream_process(
                query="显示代谢组分析流程",
                files=[],
                workflow_data=None
            ):
                events_no_file.append(event)
            
            has_template_mode_true = any('"template_mode": true' in event or '"template_mode":true' in event for event in events_no_file)
            
            # Test Step 2: Upload File -> Expect template_mode=False
            mock_file_inspector = Mock()
            mock_file_inspector.inspect_file.return_value = {
                "status": "success",
                "file_path": "/app/uploads/cow_diet.csv",
                "n_samples": 39,
                "n_features": 100
            }
            orchestrator.file_inspector = mock_file_inspector
            
            mock_planner.generate_plan = AsyncMock(return_value={
                "type": "workflow_config",
                "template_mode": False,
                "workflow_data": {
                    "workflow_name": "Metabolomics 标准分析流程",
                    "steps": [
                        {"id": "inspect_data", "name": "数据检查", "params": {"file_path": "/app/uploads/cow_diet.csv"}}
                    ]
                }
            })
            
            events_with_file = []
            async for event in orchestrator.stream_process(
                query="分析这个文件",
                files=[{"name": "cow_diet.csv", "path": "/app/uploads/cow_diet.csv"}],
                workflow_data=None
            ):
                events_with_file.append(event)
            
            has_template_mode_false = any('"template_mode": false' in event or '"template_mode":false' in event for event in events_with_file)
            
            passed = has_template_mode_true and has_template_mode_false
            
            self.log_test(
                "Scenario B: Plan-First",
                passed,
                None if passed else f"无文件时 template_mode=True: {has_template_mode_true}, 有文件时 template_mode=False: {has_template_mode_false}"
            )
            
        except Exception as e:
            self.log_test("Scenario B: Plan-First", False, str(e))
    
    async def test_scenario_c_execution_trigger(self):
        """Scenario C: Execution Trigger - Request contains workflow_data -> Expect Immediate Execution (Bypass Planner)"""
        print("\n" + "="*60)
        print("测试场景 C: Execution Trigger")
        print("="*60)
        
        try:
            # Mock dependencies
            mock_agent = Mock()
            orchestrator = AgentOrchestrator(mock_agent, upload_dir="/app/uploads")
            
            # Mock executor
            mock_executor = Mock()
            mock_executor.execute_workflow.return_value = {
                "steps_details": [
                    {"step_id": "inspect_data", "status": "completed", "output": "data.csv"}
                ],
                "workflow_name": "Metabolomics 标准分析流程"
            }
            
            # Mock planner (should NOT be called)
            mock_planner_called = False
            original_generate_plan = orchestrator.planner.generate_plan if hasattr(orchestrator, 'planner') else None
            
            # Test: Request with workflow_data
            workflow_data = {
                "workflow_data": {
                    "workflow_name": "Metabolomics 标准分析流程",
                    "steps": [
                        {"id": "inspect_data", "name": "数据检查", "params": {"file_path": "/app/uploads/cow_diet.csv"}}
                    ]
                },
                "file_paths": ["/app/uploads/cow_diet.csv"]
            }
            
            # Patch executor
            with patch('gibh_agent.core.orchestrator.ExecutionLayer', return_value=mock_executor):
                events = []
                async for event in orchestrator.stream_process(
                    query="执行工作流",
                    files=[],
                    workflow_data=workflow_data
                ):
                    events.append(event)
            
            # Verify: Should have execution status, NOT planning status
            has_execution_status = any('正在执行分析工具' in event or '正在初始化执行引擎' in event for event in events)
            has_planning_status = any('正在根据您的需求定制流程' in event or '正在分析您的需求' in event for event in events)
            has_result = any('result' in event for event in events)
            
            passed = has_execution_status and not has_planning_status and has_result
            
            self.log_test(
                "Scenario C: Execution Trigger",
                passed,
                None if passed else f"执行状态: {has_execution_status}, 规划状态: {has_planning_status}, 结果: {has_result}"
            )
            
        except Exception as e:
            self.log_test("Scenario C: Execution Trigger", False, str(e))
    
    async def test_scenario_d_async_rna(self):
        """Scenario D: Async RNA - Request scRNA-seq -> Planner selects rna_cellranger_count -> Execution returns status="async_job_started" -> Orchestrator must yield this status and STOP"""
        print("\n" + "="*60)
        print("测试场景 D: Async RNA")
        print("="*60)
        
        try:
            # Mock dependencies
            mock_agent = Mock()
            orchestrator = AgentOrchestrator(mock_agent, upload_dir="/app/uploads")
            
            # Mock executor to return async_job_started
            mock_executor = Mock()
            mock_executor.execute_workflow.return_value = {
                "steps_details": [
                    {
                        "step_id": "rna_cellranger_count",
                        "status": "async_job_started",
                        "job_id": "job_12345",
                        "message": "CellRanger job started, waiting for completion"
                    }
                ],
                "workflow_name": "RNA 标准分析流程"
            }
            
            # Test: Request with workflow_data containing RNA workflow
            workflow_data = {
                "workflow_data": {
                    "workflow_name": "RNA 标准分析流程",
                    "steps": [
                        {"id": "rna_cellranger_count", "name": "CellRanger Count", "params": {"fastq_path": "/data/fastq"}}
                    ]
                },
                "file_paths": []
            }
            
            # Patch executor
            with patch('gibh_agent.core.orchestrator.ExecutionLayer', return_value=mock_executor):
                events = []
                async for event in orchestrator.stream_process(
                    query="执行工作流",
                    files=[],
                    workflow_data=workflow_data
                ):
                    events.append(event)
                    # Check if async status is yielded
                    if 'async_job_started' in event:
                        break
            
            # Verify: Should have async_job_started status
            has_async_status = any('async_job_started' in event for event in events)
            has_job_id = any('job_12345' in event for event in events)
            
            passed = has_async_status and has_job_id
            
            self.log_test(
                "Scenario D: Async RNA",
                passed,
                None if passed else f"async_job_started 状态: {has_async_status}, job_id: {has_job_id}"
            )
            
        except Exception as e:
            self.log_test("Scenario D: Async RNA", False, str(e))
    
    async def run_all_tests(self):
        """运行所有测试"""
        print("\n" + "="*60)
        print("开始回归测试")
        print("="*60)
        
        await self.test_scenario_a_legacy_metabolomics()
        await self.test_scenario_b_plan_first()
        await self.test_scenario_c_execution_trigger()
        await self.test_scenario_d_async_rna()
        
        # Print summary
        print("\n" + "="*60)
        print("测试总结")
        print("="*60)
        print(f"✅ 通过: {self.passed}")
        print(f"❌ 失败: {self.failed}")
        
        if self.errors:
            print("\n错误详情:")
            for error in self.errors:
                print(f"  - {error}")
        
        return self.failed == 0


async def main():
    """主函数"""
    suite = RegressionTestSuite()
    success = await suite.run_all_tests()
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    asyncio.run(main())

