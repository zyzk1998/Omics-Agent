#!/usr/bin/env python3
"""
RNA分析全流程自测程序

测试目标：
1. 验证RNA分析全流程（从10x数据到最终报告）
2. 搞清楚数据诊断报告生成成功的逻辑
3. 为代谢组分析流程数据诊断过程失败做好充足准备

测试数据：test_data文件夹内的10x数据文件
"""
import os
import sys
import asyncio
import logging
from pathlib import Path
from typing import Dict, Any, List

# 添加项目根目录到路径
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from gibh_agent.core.file_inspector import FileInspector
from gibh_agent.core.prompt_manager import PromptManager
from gibh_agent.core.tool_retriever import ToolRetriever
from gibh_agent.agents.specialists.rna_agent import RNAAgent
from gibh_agent.core.orchestrator import Orchestrator
from gibh_agent.core.llm_client import LLMClientFactory
from gibh_agent.core.llm_cloud_providers import get_default_chat_model

# 配置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('test_rna_workflow.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)


async def test_data_diagnosis():
    """测试1: 数据诊断报告生成"""
    logger.info("=" * 80)
    logger.info("测试1: 数据诊断报告生成")
    logger.info("=" * 80)
    
    # 查找test_data中的10x数据文件
    test_data_dir = project_root / "test_data"
    h5ad_file = test_data_dir / "pbmc_1k_v3_filtered.h5ad"
    
    if not h5ad_file.exists():
        logger.error(f"❌ 测试数据文件不存在: {h5ad_file}")
        return False
    
    logger.info(f"📁 使用测试数据: {h5ad_file}")
    
    # 初始化组件
    upload_dir = str(test_data_dir)
    prompt_manager = PromptManager()
    tool_retriever = ToolRetriever()
    
    # 初始化FileInspector
    file_inspector = FileInspector(upload_dir)
    
    # 初始化RNAAgent
    llm_client = LLMClientFactory.create_for_model(get_default_chat_model())
    rna_agent = RNAAgent(
        prompt_manager=prompt_manager,
        tool_retriever=tool_retriever,
        upload_dir=upload_dir,
        llm_client=llm_client
    )
    
    # 执行数据诊断
    logger.info("🔍 开始数据诊断...")
    try:
        file_metadata = file_inspector.inspect_file(str(h5ad_file))
        logger.info(f"✅ 文件检查完成: {file_metadata.get('status')}")
        
        # 调用数据诊断
        diagnosis = await rna_agent._perform_data_diagnosis(
            file_metadata=file_metadata,
            omics_type="scRNA",
        )
        
        if diagnosis:
            logger.info(f"✅ 数据诊断报告生成成功，长度: {len(diagnosis)}")
            logger.info(f"📝 诊断报告预览:\n{diagnosis[:500]}...")
            return True
        else:
            logger.error("❌ 数据诊断报告生成失败（返回None）")
            return False
            
    except Exception as e:
        logger.error(f"❌ 数据诊断失败: {e}", exc_info=True)
        return False


async def test_full_workflow():
    """测试2: 完整工作流执行"""
    logger.info("=" * 80)
    logger.info("测试2: 完整工作流执行")
    logger.info("=" * 80)
    
    # 查找test_data中的10x数据文件
    test_data_dir = project_root / "test_data"
    h5ad_file = test_data_dir / "pbmc_1k_v3_filtered.h5ad"
    
    if not h5ad_file.exists():
        logger.error(f"❌ 测试数据文件不存在: {h5ad_file}")
        return False
    
    logger.info(f"📁 使用测试数据: {h5ad_file}")
    
    # 初始化组件
    upload_dir = str(test_data_dir)
    prompt_manager = PromptManager()
    tool_retriever = ToolRetriever()
    llm_client = LLMClientFactory.create_for_model(get_default_chat_model())
    
    # 初始化RNAAgent
    rna_agent = RNAAgent(
        prompt_manager=prompt_manager,
        tool_retriever=tool_retriever,
        upload_dir=upload_dir,
        llm_client=llm_client
    )
    
    # 初始化Orchestrator
    orchestrator = Orchestrator(rna_agent)
    
    # 构建工作流配置（简化版，只测试关键步骤）
    workflow_config = {
        "workflow_name": "RNA分析测试流程",
        "workflow_data": {
            "steps": [
                {
                    "step_id": "rna_qc_filter",
                    "tool_id": "rna_qc_filter",
                    "name": "质量控制过滤",
                    "params": {
                        "adata_path": str(h5ad_file),
                        "min_genes": 200,
                        "max_mt": 20
                    }
                },
                {
                    "step_id": "rna_normalize",
                    "tool_id": "rna_normalize",
                    "name": "数据标准化",
                    "params": {
                        "adata_path": "<rna_qc_filter_output>",
                        "target_sum": 10000
                    }
                },
                {
                    "step_id": "rna_pca",
                    "tool_id": "rna_pca",
                    "name": "主成分分析",
                    "params": {
                        "adata_path": "<rna_normalize_output>",
                        "n_comps": 50
                    }
                }
            ]
        },
        "file_paths": [str(h5ad_file)]
    }
    
    logger.info("🚀 开始执行工作流...")
    try:
        # 执行工作流（使用异步生成器）
        results = None
        async for event in orchestrator.execute_workflow_stream(workflow_config):
            event_type = event.get("type", "")
            event_data = event.get("data", {})
            
            if event_type == "status":
                logger.info(f"📊 状态更新: {event_data.get('content', '')}")
            elif event_type == "step_result":
                logger.info(f"✅ 步骤结果: {event_data.get('step_name', '')}")
            elif event_type == "diagnosis":
                logger.info(f"💡 AI专家报告生成: {len(event_data.get('report_data', {}).get('summary', ''))} 字符")
            elif event_type == "result":
                results = event_data
                logger.info("✅ 工作流执行完成")
        
        if results:
            logger.info(f"✅ 工作流执行成功，共 {len(results.get('steps_details', []))} 个步骤")
            return True
        else:
            logger.error("❌ 工作流执行失败（未返回结果）")
            return False
            
    except Exception as e:
        logger.error(f"❌ 工作流执行失败: {e}", exc_info=True)
        return False


async def test_ai_report_generation():
    """测试3: AI专家分析报告生成"""
    logger.info("=" * 80)
    logger.info("测试3: AI专家分析报告生成")
    logger.info("=" * 80)
    
    # 初始化组件
    upload_dir = str(project_root / "test_data")
    prompt_manager = PromptManager()
    tool_retriever = ToolRetriever()
    llm_client = LLMClientFactory.create_for_model(get_default_chat_model())
    
    # 初始化RNAAgent
    rna_agent = RNAAgent(
        prompt_manager=prompt_manager,
        tool_retriever=tool_retriever,
        upload_dir=upload_dir,
        llm_client=llm_client
    )
    
    # 模拟执行结果
    steps_results = [
        {
            "step_name": "质量控制过滤",
            "status": "success",
            "data": {
                "n_cells": 1000,
                "n_genes": 2000,
                "summary": "质量控制完成"
            }
        },
        {
            "step_name": "数据标准化",
            "status": "success",
            "data": {
                "summary": "标准化完成"
            }
        },
        {
            "step_name": "主成分分析",
            "status": "success",
            "data": {
                "n_components": 50,
                "variance_explained": [0.3, 0.2, 0.1],
                "summary": "PCA分析完成"
            }
        }
    ]
    
    logger.info("📝 开始生成AI专家分析报告...")
    try:
        summary = await rna_agent._generate_analysis_summary(
            steps_results=steps_results,
            omics_type="scRNA",
            workflow_name="RNA分析测试",
            summary_context={
                "has_failures": False,
                "has_warnings": False,
                "failed_steps": [],
                "warning_steps": [],
                "successful_steps": steps_results
            }
        )
        
        if summary:
            logger.info(f"✅ AI专家分析报告生成成功，长度: {len(summary)}")
            logger.info(f"📝 报告预览:\n{summary[:500]}...")
            
            # 检查报告质量
            has_biology = "生物" in summary or "biology" in summary.lower() or "机制" in summary
            has_analysis = "分析" in summary or "analysis" in summary.lower()
            is_programmer_log = "错误" in summary or "error" in summary.lower() or "失败" in summary or "failed" in summary.lower()
            
            logger.info(f"  - 包含生物学内容: {has_biology}")
            logger.info(f"  - 包含分析内容: {has_analysis}")
            logger.info(f"  - 是否为程序员日志: {is_programmer_log}")
            
            if has_biology and has_analysis and not is_programmer_log:
                logger.info("✅ 报告质量检查通过：包含生物学分析，非程序员日志")
                return True
            else:
                logger.warning("⚠️ 报告质量可能不足：缺少生物学分析或包含程序员日志")
                return False
        else:
            logger.error("❌ AI专家分析报告生成失败（返回None）")
            return False
            
    except Exception as e:
        logger.error(f"❌ AI专家分析报告生成失败: {e}", exc_info=True)
        return False


async def main():
    """主测试函数"""
    logger.info("=" * 80)
    logger.info("RNA分析全流程自测程序")
    logger.info("=" * 80)
    
    results = {
        "数据诊断": False,
        "完整工作流": False,
        "AI专家报告": False
    }
    
    # 测试1: 数据诊断
    try:
        results["数据诊断"] = await test_data_diagnosis()
    except Exception as e:
        logger.error(f"❌ 测试1失败: {e}", exc_info=True)
    
    # 测试2: 完整工作流
    try:
        results["完整工作流"] = await test_full_workflow()
    except Exception as e:
        logger.error(f"❌ 测试2失败: {e}", exc_info=True)
    
    # 测试3: AI专家报告
    try:
        results["AI专家报告"] = await test_ai_report_generation()
    except Exception as e:
        logger.error(f"❌ 测试3失败: {e}", exc_info=True)
    
    # 输出测试结果
    logger.info("=" * 80)
    logger.info("测试结果汇总")
    logger.info("=" * 80)
    for test_name, result in results.items():
        status = "✅ 通过" if result else "❌ 失败"
        logger.info(f"{test_name}: {status}")
    
    total_passed = sum(1 for r in results.values() if r)
    logger.info(f"\n总计: {total_passed}/{len(results)} 测试通过")
    
    return all(results.values())


if __name__ == "__main__":
    success = asyncio.run(main())
    sys.exit(0 if success else 1)
