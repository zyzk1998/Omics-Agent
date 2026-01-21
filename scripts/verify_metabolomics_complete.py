#!/usr/bin/env python3
"""
验证脚本：验证代谢组学完整流程
确保所有步骤正确执行，即使某些步骤返回 warning
"""
import sys
import os
import asyncio
from pathlib import Path

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
RESET = '\033[0m'

def print_check(check_num: int, description: str, passed: bool, details: str = ""):
    """Print check result with color"""
    status = f"{GREEN}✓ PASSED{RESET}" if passed else f"{RED}✗ FAILED{RESET}"
    print(f"[Check {check_num}] {description}: {status}")
    if details:
        print(f"  {details}")
    return passed

async def main():
    """Main verification function"""
    print("=" * 80)
    print("🔬 代谢组学完整流程验证脚本")
    print("=" * 80)
    print()
    
    # Setup
    print("📋 Step 1: 初始化组件...")
    
    # Initialize components
    try:
        # Try to create agent using factory
        agent_instance = create_agent()
        # Get metabolomics agent from the main agent
        agent = agent_instance.agents.get("metabolomics_agent")
        if not agent:
            raise ValueError("Metabolomics agent not found")
    except Exception as e:
        print(f"⚠️  使用 create_agent 失败: {e}，尝试直接初始化...")
        # Fallback: initialize directly
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
    
    # Test data path
    test_data_path = project_root / "test_data" / "cow_diet.csv"
    
    # Create test data if it doesn't exist
    if not test_data_path.exists():
        print(f"⚠️  测试数据不存在，正在创建: {test_data_path}")
        test_data_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Create dummy cow_diet.csv with Diet column (0/1)
        import pandas as pd
        import numpy as np
        
        np.random.seed(42)
        n_samples = 20
        n_metabolites = 50
        
        # Create sample IDs
        sample_ids = [f"Sample_{i+1}" for i in range(n_samples)]
        
        # Create Diet column (0 or 1)
        diet_values = [0] * (n_samples // 2) + [1] * (n_samples - n_samples // 2)
        np.random.shuffle(diet_values)
        
        # Create metabolite data
        data = {
            "Diet": diet_values
        }
        
        # Add metabolite columns
        for i in range(n_metabolites):
            col_name = f"Metabolite_{i+1}"
            # Different means for different diets
            mean_0 = np.random.uniform(1, 5)
            mean_1 = mean_0 * np.random.uniform(1.2, 2.0)  # Some metabolites differ between diets
            values = []
            for diet in diet_values:
                if diet == 0:
                    values.append(np.random.normal(mean_0, 0.5))
                else:
                    values.append(np.random.normal(mean_1, 0.5))
            data[col_name] = values
        
        df = pd.DataFrame(data, index=sample_ids)
        df.to_csv(test_data_path)
        print(f"✅ 测试数据已创建: {test_data_path}")
    
    print(f"✅ 使用测试数据: {test_data_path}")
    print()
    
    # Inspect file
    print("📋 Step 2: 检查文件...")
    file_metadata = file_inspector.inspect_file(str(test_data_path))
    print(f"✅ 文件检查完成")
    print()
    
    # Generate plan
    print("📋 Step 3: 生成工作流计划...")
    workflow_config = await planner.generate_plan(
        user_query="做全流程分析",
        file_metadata=file_metadata,
        domain_name="Metabolomics"
    )
    
    # Extract steps from workflow_config
    if "workflow_config" in workflow_config:
        workflow_data = workflow_config["workflow_config"]
    elif "workflow_data" in workflow_config:
        workflow_data = workflow_config["workflow_data"]
    else:
        workflow_data = workflow_config
    
    steps = workflow_data.get("steps", [])
    
    print(f"✅ 工作流计划生成完成，共 {len(steps)} 个步骤")
    if len(steps) > 0:
        print(f"   步骤列表: {[s.get('id', 'N/A') for s in steps]}")
    print()
    
    # Execute workflow
    print("📋 Step 4: 执行工作流...")
    results = executor.execute_workflow(
        workflow_data=workflow_data,
        file_paths=[str(test_data_path)],
        agent=agent
    )
    
    steps_details = results.get("steps_details", [])
    print(f"✅ 工作流执行完成，共 {len(steps_details)} 个步骤详情")
    print()
    
    # Generate summary
    print("📋 Step 5: 生成 AI 诊断报告...")
    diagnosis = ""
    try:
        steps_results = results.get("steps_results", [])
        if not steps_results:
            # Extract from steps_details
            steps_results = []
            for step_detail in steps_details:
                if "step_result" in step_detail:
                    steps_results.append(step_detail["step_result"])
        
        # Build summary context
        failed_steps = [s for s in steps_details if s.get("status") == "error"]
        warning_steps = [s for s in steps_details if s.get("status") == "warning"]
        successful_steps = [s for s in steps_details if s.get("status") == "success"]
        
        summary_context = {
            "has_failures": len(failed_steps) > 0,
            "has_warnings": len(warning_steps) > 0,
            "failed_steps": failed_steps,
            "warning_steps": warning_steps,
            "successful_steps": successful_steps,
            "workflow_status": results.get("status", "unknown")
        }
        
        try:
            diagnosis = await agent._generate_analysis_summary(
                steps_results,
                omics_type="Metabolomics",
                workflow_name=workflow_data.get("workflow_name", "代谢组学分析"),
                summary_context=summary_context
            )
        except Exception as e:
            print(f"⚠️  LLM 生成诊断失败: {e}，使用后备摘要")
            # Use fallback summary
            from gibh_agent.core.orchestrator import AgentOrchestrator
            orchestrator = AgentOrchestrator(agent=None, upload_dir=str(project_root / "test_data"))
            diagnosis = orchestrator._generate_fallback_summary(
                successful_steps, warning_steps, failed_steps, steps_details
            )
        
        if not diagnosis:
            # Final fallback: check results
            diagnosis = results.get("report_data", {}).get("diagnosis", "")
            if not diagnosis:
                # Generate basic fallback
                from gibh_agent.core.orchestrator import AgentOrchestrator
                orchestrator = AgentOrchestrator(agent=None, upload_dir=str(project_root / "test_data"))
                diagnosis = orchestrator._generate_fallback_summary(
                    successful_steps, warning_steps, failed_steps, steps_details
                )
        
        print(f"✅ 诊断报告生成完成，长度: {len(diagnosis) if diagnosis else 0} 字符")
    except Exception as e:
        print(f"⚠️  诊断报告生成失败: {e}")
        # Final fallback
        from gibh_agent.core.orchestrator import AgentOrchestrator
        orchestrator = AgentOrchestrator(agent=None, upload_dir=str(project_root / "test_data"))
        failed_steps = [s for s in steps_details if s.get("status") == "error"]
        warning_steps = [s for s in steps_details if s.get("status") == "warning"]
        successful_steps = [s for s in steps_details if s.get("status") == "success"]
        diagnosis = orchestrator._generate_fallback_summary(
            successful_steps, warning_steps, failed_steps, steps_details
        )
    print()
    
    # Verification checks
    print("=" * 80)
    print("🔍 验证检查")
    print("=" * 80)
    print()
    
    all_passed = True
    
    # Check 1: Workflow has 7 steps
    check1 = len(steps) == 7
    all_passed = print_check(1, f"工作流有 7 个步骤", check1, f"实际: {len(steps)} 个步骤") and all_passed
    
    # Check 2: PLS-DA status is success
    plsda_step = next((s for s in steps_details if s.get("step_id") == "metabolomics_plsda"), None)
    check2 = plsda_step is not None and plsda_step.get("status") == "success"
    plsda_status = plsda_step.get("status") if plsda_step else "not found"
    all_passed = print_check(2, f"metabolomics_plsda 状态为 success", check2, f"实际状态: {plsda_status}") and all_passed
    
    # Check 3: Pathway enrichment status is success OR warning (NOT error)
    enrichment_step = next((s for s in steps_details if s.get("step_id") == "metabolomics_pathway_enrichment"), None)
    check3 = enrichment_step is not None and enrichment_step.get("status") in ["success", "warning"]
    enrichment_status = enrichment_step.get("status") if enrichment_step else "not found"
    all_passed = print_check(3, f"metabolomics_pathway_enrichment 状态为 success 或 warning", check3, f"实际状态: {enrichment_status}") and all_passed
    
    # Check 4: Result contains diagnosis string (must be > 50 chars)
    check4 = diagnosis is not None and len(diagnosis) > 50
    diagnosis_len = len(diagnosis) if diagnosis else 0
    all_passed = print_check(4, f"结果包含 diagnosis 字符串（长度 > 50）", check4, f"实际长度: {diagnosis_len} 字符") and all_passed
    
    # Check 5: Result contains steps_details
    check5 = len(steps_details) > 0
    all_passed = print_check(5, f"结果包含 steps_details", check5, f"实际步骤数: {len(steps_details)}") and all_passed
    
    print()
    print("=" * 80)
    if all_passed:
        print(f"{GREEN}✅ ALL CHECKS PASSED{RESET}")
        print("=" * 80)
        return 0
    else:
        print(f"{RED}❌ SOME CHECKS FAILED{RESET}")
        print("=" * 80)
        return 1

if __name__ == "__main__":
    exit_code = asyncio.run(main())
    sys.exit(exit_code)

