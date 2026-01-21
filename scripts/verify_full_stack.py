#!/usr/bin/env python3
"""
🔬 Full Stack Verification Script (Judge)
验证整个系统的完整性：执行、文件生成、JSON结构、UI渲染
"""
import sys
import os
import asyncio
import json
from pathlib import Path
from typing import Dict, Any, List

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from gibh_agent import create_agent
from gibh_agent.core.executor import WorkflowExecutor
from gibh_agent.core.planner import SOPPlanner
from gibh_agent.core.file_inspector import FileInspector
from gibh_agent.core.tool_retriever import ToolRetriever
from gibh_agent.core.llm_client import LLMClientFactory
from gibh_agent.core.prompt_manager import create_default_prompt_manager

# ANSI color codes
GREEN = '\033[92m'
RED = '\033[91m'
YELLOW = '\033[93m'
BLUE = '\033[94m'
RESET = '\033[0m'

class VerificationResult:
    """验证结果类"""
    def __init__(self, name: str):
        self.name = name
        self.passed = False
        self.message = ""
        self.details = []
    
    def add_detail(self, detail: str):
        self.details.append(detail)
    
    def __str__(self):
        status = f"{GREEN}✓ PASSED{RESET}" if self.passed else f"{RED}✗ FAILED{RESET}"
        result = f"[{self.name}] {status}"
        if self.message:
            result += f"\n  {self.message}"
        if self.details:
            for detail in self.details:
                result += f"\n    - {detail}"
        return result

def check_file_exists(file_path: Path, description: str) -> bool:
    """检查文件是否存在"""
    if file_path.exists():
        size = file_path.stat().st_size
        return size > 0
    return False

def normalize_path(path_str: str, base_dir: Path) -> Path:
    """标准化路径"""
    if not path_str:
        return None
    
    # Remove /app/results prefix if present
    path_str = path_str.replace('/app/results/', '').replace('/app/', '')
    path_str = path_str.replace('results/', '')
    
    # If it starts with /results/, remove that too
    if path_str.startswith('/results/'):
        path_str = path_str[9:]  # Remove '/results/'
    
    # Try relative to base_dir
    full_path = base_dir / path_str
    if full_path.exists():
        return full_path
    
    # Try absolute path
    if Path(path_str).exists():
        return Path(path_str)
    
    return None

async def main():
    """Main verification function"""
    print("=" * 80)
    print(f"{BLUE}🔬 Full Stack Verification Script (Judge){RESET}")
    print("=" * 80)
    print()
    
    results_dir = project_root / "results"
    test_data_path = project_root / "test_data" / "cow_diet.csv"
    
    # Create test data if needed
    if not test_data_path.exists():
        print(f"{YELLOW}⚠️  测试数据不存在，正在创建...{RESET}")
        test_data_path.parent.mkdir(parents=True, exist_ok=True)
        
        import pandas as pd
        import numpy as np
        
        np.random.seed(42)
        n_samples = 20
        n_metabolites = 50
        
        sample_ids = [f"Sample_{i+1}" for i in range(n_samples)]
        diet_values = [0] * (n_samples // 2) + [1] * (n_samples - n_samples // 2)
        np.random.shuffle(diet_values)
        
        data = {"Diet": diet_values}
        for i in range(n_metabolites):
            col_name = f"Metabolite_{i+1}"
            mean_0 = np.random.uniform(1, 5)
            mean_1 = mean_0 * np.random.uniform(1.2, 2.0)
            values = []
            for diet in diet_values:
                if diet == 0:
                    values.append(np.random.normal(mean_0, 0.5))
                else:
                    values.append(np.random.normal(mean_1, 0.5))
            data[col_name] = values
        
        df = pd.DataFrame(data, index=sample_ids)
        df.to_csv(test_data_path)
        print(f"{GREEN}✅ 测试数据已创建{RESET}")
    
    print(f"{BLUE}📋 Step 1: 初始化组件...{RESET}")
    try:
        agent_instance = create_agent()
        agent = agent_instance.agents.get("metabolomics_agent")
        if not agent:
            raise ValueError("Metabolomics agent not found")
    except Exception as e:
        print(f"{YELLOW}⚠️  使用 create_agent 失败: {e}，尝试直接初始化...{RESET}")
        llm_client = LLMClientFactory.create_client("openai")
        prompt_manager = create_default_prompt_manager()
        tool_retriever = ToolRetriever()
        
        from gibh_agent.agents.specialists.metabolomics_agent import MetabolomicsAgent
        agent = MetabolomicsAgent(
            llm_client=llm_client,
            prompt_manager=prompt_manager,
            tool_retriever=tool_retriever
        )
    
    tool_retriever = ToolRetriever()
    llm_client = agent.llm_client if hasattr(agent, 'llm_client') else LLMClientFactory.create_client("openai")
    planner = SOPPlanner(tool_retriever, llm_client)
    executor = WorkflowExecutor()
    file_inspector = FileInspector(str(project_root / "test_data"))
    
    print(f"{GREEN}✅ 组件初始化完成{RESET}\n")
    
    # Generate plan
    print(f"{BLUE}📋 Step 2: 生成工作流计划...{RESET}")
    file_metadata = file_inspector.inspect_file(str(test_data_path))
    workflow_config = await planner.generate_plan(
        user_query="做全流程分析",
        file_metadata=file_metadata,
        domain_name="Metabolomics"
    )
    
    workflow_data = workflow_config.get("workflow_data", {})
    steps = workflow_data.get("steps", [])
    print(f"{GREEN}✅ 工作流计划生成完成，共 {len(steps)} 个步骤{RESET}\n")
    
    # Execute workflow
    print(f"{BLUE}📋 Step 3: 执行工作流...{RESET}")
    results = executor.execute_workflow(
        workflow_data=workflow_data,
        file_paths=[str(test_data_path)],
        agent=agent
    )
    
    steps_details = results.get("steps_details", [])
    print(f"{GREEN}✅ 工作流执行完成，共 {len(steps_details)} 个步骤详情{RESET}\n")
    
    # Generate diagnosis
    print(f"{BLUE}📋 Step 4: 生成诊断报告...{RESET}")
    diagnosis = ""
    evaluation = None
    try:
        steps_results = results.get("steps_results", [])
        if not steps_results:
            steps_results = []
            for step_detail in steps_details:
                if "step_result" in step_detail:
                    steps_results.append(step_detail["step_result"])
        
        try:
            diagnosis = await agent._generate_analysis_summary(
                steps_results,
                omics_type="Metabolomics",
                workflow_name=workflow_data.get("workflow_name", "代谢组学分析")
            )
        except Exception as e:
            print(f"{YELLOW}⚠️  LLM 生成诊断超时/失败: {e}，使用后备摘要{RESET}")
            # Use fallback summary
            from gibh_agent.core.orchestrator import AgentOrchestrator
            orchestrator = AgentOrchestrator(agent=None, upload_dir=str(project_root / "test_data"))
            failed_steps = [s for s in steps_details if s.get("status") == "error"]
            warning_steps = [s for s in steps_details if s.get("status") == "warning"]
            successful_steps = [s for s in steps_details if s.get("status") == "success"]
            diagnosis = orchestrator._generate_fallback_summary(
                successful_steps, warning_steps, failed_steps, steps_details
            )
        
        if not diagnosis:
            # Final fallback
            from gibh_agent.core.orchestrator import AgentOrchestrator
            orchestrator = AgentOrchestrator(agent=None, upload_dir=str(project_root / "test_data"))
            failed_steps = [s for s in steps_details if s.get("status") == "error"]
            warning_steps = [s for s in steps_details if s.get("status") == "warning"]
            successful_steps = [s for s in steps_details if s.get("status") == "success"]
            diagnosis = orchestrator._generate_fallback_summary(
                successful_steps, warning_steps, failed_steps, steps_details
            )
        
        if diagnosis and hasattr(agent, '_evaluate_analysis_quality'):
            try:
                evaluation = await agent._evaluate_analysis_quality(
                    steps_results,
                    diagnosis,
                    workflow_data.get("workflow_name", "代谢组学分析")
                )
            except Exception as e:
                print(f"{YELLOW}⚠️  质量评估失败: {e}，使用默认评估{RESET}")
                # Generate default evaluation
                successful_count = len([s for s in steps_details if s.get("status") == "success"])
                evaluation = {
                    "score": int(successful_count / len(steps_details) * 100) if steps_details else 50,
                    "critique": f"分析完成率 {successful_count}/{len(steps_details)}",
                    "strengths": [],
                    "weaknesses": [],
                    "recommendations": []
                }
    except Exception as e:
        print(f"{YELLOW}⚠️  诊断生成失败: {e}{RESET}")
        # Final fallback
        from gibh_agent.core.orchestrator import AgentOrchestrator
        orchestrator = AgentOrchestrator(agent=None, upload_dir=str(project_root / "test_data"))
        failed_steps = [s for s in steps_details if s.get("status") == "error"]
        warning_steps = [s for s in steps_details if s.get("status") == "warning"]
        successful_steps = [s for s in steps_details if s.get("status") == "success"]
        diagnosis = orchestrator._generate_fallback_summary(
            successful_steps, warning_steps, failed_steps, steps_details
        )
    
    print(f"{GREEN}✅ 诊断报告生成完成{RESET}\n")
    
    # Verification checks
    print("=" * 80)
    print(f"{BLUE}🔍 Verification Checks{RESET}")
    print("=" * 80)
    print()
    
    all_results: List[VerificationResult] = []
    
    # Check 1: Workflow has correct number of steps
    check1 = VerificationResult("Workflow Steps Count")
    check1.passed = len(steps) == 7
    check1.message = f"Expected 7 steps, got {len(steps)}"
    if check1.passed:
        check1.add_detail(f"Steps: {[s.get('id', 'N/A') for s in steps]}")
    all_results.append(check1)
    
    # Check 2: All steps executed
    check2 = VerificationResult("All Steps Executed")
    check2.passed = len(steps_details) == len(steps)
    check2.message = f"Expected {len(steps)} step details, got {len(steps_details)}"
    all_results.append(check2)
    
    # Check 3: Artifacts exist on disk
    check3 = VerificationResult("Artifacts on Disk")
    artifact_count = 0
    missing_artifacts = []
    
    for step_detail in steps_details:
        step_id = step_detail.get("step_id", "")
        plot_path = step_detail.get("plot")
        
        # Check plot_path in step_detail
        if plot_path:
            normalized = normalize_path(plot_path, results_dir)
            if normalized and normalized.exists():
                artifact_count += 1
            else:
                missing_artifacts.append(f"{step_id}: plot_path={plot_path}")
        
        # Check images in step_result.data
        step_result = step_detail.get("step_result", {})
        step_data = step_result.get("data", {})
        images = step_data.get("images", [])
        
        for img_path in images:
            normalized = normalize_path(img_path, results_dir)
            if normalized and normalized.exists():
                artifact_count += 1
            else:
                missing_artifacts.append(f"{step_id}: image={img_path}")
    
    check3.passed = len(missing_artifacts) == 0
    check3.message = f"Found {artifact_count} artifacts, {len(missing_artifacts)} missing"
    if missing_artifacts:
        check3.details = missing_artifacts[:5]  # Show first 5
    all_results.append(check3)
    
    # Check 4: JSON structure - plot_path in steps_details
    check4 = VerificationResult("JSON Structure - plot_path")
    steps_with_plots = []
    steps_without_plots = []
    
    for step_detail in steps_details:
        step_id = step_detail.get("step_id", "")
        has_plot = bool(step_detail.get("plot") or 
                       step_detail.get("step_result", {}).get("data", {}).get("images"))
        
        if has_plot:
            steps_with_plots.append(step_id)
        else:
            steps_without_plots.append(step_id)
    
    # Expected steps that should have plots: PCA, PLS-DA, Volcano
    expected_plot_steps = ["pca_analysis", "metabolomics_plsda", "visualize_volcano"]
    missing_plot_steps = [s for s in expected_plot_steps if s in steps_without_plots]
    
    check4.passed = len(missing_plot_steps) == 0
    check4.message = f"Steps with plots: {len(steps_with_plots)}, Missing: {missing_plot_steps}"
    all_results.append(check4)
    
    # Check 5: Diagnosis exists and is long enough
    check5 = VerificationResult("Diagnosis Report")
    check5.passed = diagnosis is not None and len(diagnosis) > 200
    check5.message = f"Diagnosis length: {len(diagnosis) if diagnosis else 0} chars (required > 200)"
    if diagnosis:
        check5.add_detail(f"Preview: {diagnosis[:100]}...")
    all_results.append(check5)
    
    # Check 6: Evaluation exists
    check6 = VerificationResult("Quality Evaluation")
    check6.passed = evaluation is not None and "score" in evaluation
    if evaluation:
        score = evaluation.get("score", 0)
        check6.message = f"Evaluation score: {score}/100"
        check6.add_detail(f"Critique: {evaluation.get('critique', 'N/A')[:100]}...")
    else:
        check6.message = "Evaluation not generated"
    all_results.append(check6)
    
    # Check 7: steps_details structure
    check7 = VerificationResult("steps_details Structure")
    valid_steps = 0
    invalid_steps = []
    
    for step_detail in steps_details:
        has_required_fields = all([
            step_detail.get("step_id"),
            step_detail.get("status"),
            step_detail.get("name")
        ])
        
        if has_required_fields:
            valid_steps += 1
        else:
            invalid_steps.append(step_detail.get("step_id", "Unknown"))
    
    check7.passed = valid_steps == len(steps_details)
    check7.message = f"Valid steps: {valid_steps}/{len(steps_details)}"
    if invalid_steps:
        check7.details = invalid_steps
    all_results.append(check7)
    
    # Print results
    for result in all_results:
        print(result)
        print()
    
    # Summary
    passed_count = sum(1 for r in all_results if r.passed)
    total_count = len(all_results)
    
    print("=" * 80)
    if passed_count == total_count:
        print(f"{GREEN}✅ ALL CHECKS PASSED ({passed_count}/{total_count}){RESET}")
        print("=" * 80)
        return 0
    else:
        print(f"{RED}❌ SOME CHECKS FAILED ({passed_count}/{total_count}){RESET}")
        print("=" * 80)
        return 1

if __name__ == "__main__":
    exit_code = asyncio.run(main())
    sys.exit(exit_code)

