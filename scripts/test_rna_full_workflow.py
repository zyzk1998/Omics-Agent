#!/usr/bin/env python3
"""
RNA分析全流程自测程序

测试目标：
1. 验证10x格式数据检测和自动跳过cellranger步骤
2. 验证数据诊断报告生成逻辑
3. 验证AI专家分析报告生成（生信专家版本，非程序员日志版本）
4. 验证所有工具执行流程

测试数据：test_data文件夹中的10x数据文件
"""
import os
import sys
import asyncio
import logging
from pathlib import Path
from typing import Dict, Any, List
from datetime import datetime

# 添加项目根目录到路径
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from gibh_agent.core.file_inspector import FileInspector
from gibh_agent.core.llm_client import LLMClientFactory
from gibh_agent.core.prompt_manager import PromptManager
from gibh_agent.core.tool_retriever import ToolRetriever
from gibh_agent.agents.specialists.rna_agent import RNAAgent
from gibh_agent.core.orchestrator import Orchestrator
from gibh_agent.core.executor import WorkflowExecutor

# 配置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(f'test_rna_workflow_{datetime.now().strftime("%Y%m%d_%H%M%S")}.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)


def find_10x_data_directories(test_data_dir: Path) -> List[Path]:
    """
    在test_data目录中查找10x数据目录
    
    Returns:
        10x数据目录路径列表
    """
    tenx_dirs = []
    
    # 递归搜索包含matrix.mtx的目录
    for root, dirs, files in os.walk(test_data_dir):
        if 'matrix.mtx' in files or 'matrix.mtx.gz' in files:
            # 检查是否包含barcodes或features文件
            has_barcodes = any('barcodes' in f for f in files)
            has_features = any(f in files for f in ['features.tsv', 'features.tsv.gz', 'genes.tsv', 'genes.tsv.gz'])
            
            if has_barcodes or has_features:
                tenx_dirs.append(Path(root))
                logger.info(f"✅ 找到10x数据目录: {root}")
    
    return tenx_dirs


async def test_data_diagnosis(file_path: str, upload_dir: str) -> Dict[str, Any]:
    """
    测试数据诊断报告生成
    
    Args:
        file_path: 10x数据文件路径
        upload_dir: 上传目录
    
    Returns:
        诊断报告结果
    """
    logger.info("=" * 80)
    logger.info("📊 测试1: 数据诊断报告生成")
    logger.info("=" * 80)
    
    try:
        # 初始化FileInspector
        file_inspector = FileInspector(upload_dir)
        
        # 检查文件
        logger.info(f"🔍 检查文件: {file_path}")
        inspection_result = file_inspector.inspect_file(file_path)
        
        if inspection_result.get("status") != "success":
            logger.error(f"❌ 文件检查失败: {inspection_result.get('error')}")
            return {"status": "error", "error": inspection_result.get("error")}
        
        logger.info(f"✅ 文件检查成功: {inspection_result.get('file_type')}")
        logger.info(f"   - 细胞数: {inspection_result.get('n_obs', 'N/A')}")
        logger.info(f"   - 基因数: {inspection_result.get('n_vars', 'N/A')}")
        
        # 初始化RNAAgent进行数据诊断
        prompt_manager = PromptManager()
        tool_retriever = ToolRetriever()
        
        rna_agent = RNAAgent(
            prompt_manager=prompt_manager,
            tool_retriever=tool_retriever,
            upload_dir=upload_dir
        )
        
        # 执行数据诊断
        logger.info("🔬 执行LLM数据诊断...")
        diagnosis_result = await rna_agent._perform_data_diagnosis(
            file_metadata=inspection_result,
            omics_type="scRNA"
        )
        
        if diagnosis_result:
            logger.info("✅ 数据诊断报告生成成功")
            logger.info(f"   报告长度: {len(diagnosis_result)} 字符")
            logger.info(f"   报告预览: {diagnosis_result[:200]}...")
            return {
                "status": "success",
                "diagnosis": diagnosis_result,
                "inspection": inspection_result
            }
        else:
            logger.error("❌ 数据诊断报告生成失败（返回None）")
            return {"status": "error", "error": "诊断报告生成失败"}
            
    except Exception as e:
        logger.error(f"❌ 数据诊断测试失败: {e}", exc_info=True)
        return {"status": "error", "error": str(e)}


async def test_workflow_execution(file_path: str, output_dir: str, upload_dir: str) -> Dict[str, Any]:
    """
    测试工作流执行（包括10x格式检测和自动跳过cellranger步骤）
    
    Args:
        file_path: 10x数据文件路径
        output_dir: 输出目录
        upload_dir: 上传目录
    
    Returns:
        执行结果
    """
    logger.info("=" * 80)
    logger.info("🚀 测试2: 工作流执行（10x格式检测和自动跳过cellranger）")
    logger.info("=" * 80)
    
    try:
        # 创建标准RNA工作流配置
        workflow_data = {
            "workflow_name": "RNA分析全流程测试",
            "steps": [
                {
                    "step_id": "rna_qc_filter",
                    "tool_id": "rna_qc_filter",
                    "name": "质量控制过滤",
                    "params": {
                        "min_genes": 200,
                        "max_mt": 20
                    }
                },
                {
                    "step_id": "rna_cellranger_count",
                    "tool_id": "rna_cellranger_count",
                    "name": "Cell Ranger 计数（异步）",
                    "params": {
                        "localcores": 8,
                        "localmem": 32,
                        "create_bam": False
                    }
                },
                {
                    "step_id": "rna_convert_cellranger_to_h5ad",
                    "tool_id": "rna_convert_cellranger_to_h5ad",
                    "name": "转换为 H5AD 格式",
                    "params": {}
                },
                {
                    "step_id": "rna_normalize",
                    "tool_id": "rna_normalize",
                    "name": "数据标准化",
                    "params": {
                        "target_sum": 10000
                    }
                }
            ]
        }
        
        # 初始化执行器
        executor = WorkflowExecutor(output_dir=output_dir, upload_dir=upload_dir)
        
        # 检测是否为10x格式
        is_10x = executor._is_10x_format(file_path)
        logger.info(f"🔍 10x格式检测结果: {is_10x}")
        
        if is_10x:
            logger.info("✅ 检测到10x格式，cellranger和convert步骤应该被自动跳过")
        else:
            logger.warning("⚠️ 未检测到10x格式，cellranger步骤将正常执行")
        
        # 执行工作流
        logger.info("🚀 开始执行工作流...")
        execution_result = executor.execute_workflow(
            workflow_data=workflow_data,
            file_paths=[file_path],
            output_dir=output_dir
        )
        
        # 检查执行结果
        if execution_result.get("status") == "success":
            logger.info("✅ 工作流执行成功")
            steps_details = execution_result.get("steps_details", [])
            
            # 检查cellranger步骤是否被跳过
            for step in steps_details:
                step_id = step.get("step_id")
                step_status = step.get("status")
                if step_id == "rna_cellranger_count":
                    if step_status == "skipped":
                        logger.info("✅ Cell Ranger步骤已正确跳过")
                    else:
                        logger.warning(f"⚠️ Cell Ranger步骤状态: {step_status}（预期: skipped）")
                elif step_id == "rna_convert_cellranger_to_h5ad":
                    if step_status == "success":
                        logger.info("✅ Convert步骤已成功执行（10x格式自动转换）")
                    else:
                        logger.warning(f"⚠️ Convert步骤状态: {step_status}")
            
            return {
                "status": "success",
                "result": execution_result,
                "is_10x_detected": is_10x
            }
        else:
            logger.error(f"❌ 工作流执行失败: {execution_result.get('error')}")
            return {"status": "error", "error": execution_result.get("error")}
            
    except Exception as e:
        logger.error(f"❌ 工作流执行测试失败: {e}", exc_info=True)
        return {"status": "error", "error": str(e)}


async def test_ai_report_generation(execution_results: Dict[str, Any], omics_type: str = "scRNA") -> Dict[str, Any]:
    """
    测试AI专家分析报告生成（验证是否为生信专家版本，而非程序员日志版本）
    
    Args:
        execution_results: 工作流执行结果
        omics_type: 组学类型
    
    Returns:
        报告生成结果
    """
    logger.info("=" * 80)
    logger.info("📝 测试3: AI专家分析报告生成（验证生信专家版本）")
    logger.info("=" * 80)
    
    try:
        # 初始化RNAAgent
        prompt_manager = PromptManager()
        tool_retriever = ToolRetriever()
        upload_dir = os.getenv("UPLOAD_DIR", "/app/uploads")
        
        rna_agent = RNAAgent(
            prompt_manager=prompt_manager,
            tool_retriever=tool_retriever,
            upload_dir=upload_dir
        )
        
        # 生成分析报告
        logger.info("🔬 生成AI专家分析报告...")
        report = await rna_agent._generate_analysis_summary(
            steps_results=execution_results.get("steps_details", []),
            omics_type=omics_type
        )
        
        if report:
            logger.info("✅ AI专家分析报告生成成功")
            logger.info(f"   报告长度: {len(report)} 字符")
            
            # 验证报告内容（检查是否为生信专家版本）
            is_programmer_log = any(keyword in report.lower() for keyword in [
                "step", "tool", "file_path", "error", "failed", "python", "code",
                "execution", "workflow", "parameter", "function"
            ])
            
            is_bioinformatics = any(keyword in report.lower() for keyword in [
                "细胞", "基因", "表达", "差异", "通路", "功能", "生物学",
                "样本", "聚类", "注释", "标志物", "机制", "代谢"
            ])
            
            logger.info(f"   程序员日志关键词检测: {is_programmer_log}")
            logger.info(f"   生信专家关键词检测: {is_bioinformatics}")
            
            if is_bioinformatics and not is_programmer_log:
                logger.info("✅ 报告内容验证通过：生信专家版本（非程序员日志版本）")
            elif is_programmer_log:
                logger.warning("⚠️ 报告内容验证失败：包含程序员日志关键词")
            else:
                logger.warning("⚠️ 报告内容验证：未检测到明显的生信专家关键词")
            
            logger.info(f"   报告预览: {report[:500]}...")
            
            return {
                "status": "success",
                "report": report,
                "is_bioinformatics": is_bioinformatics,
                "is_programmer_log": is_programmer_log
            }
        else:
            logger.error("❌ AI专家分析报告生成失败（返回None）")
            return {"status": "error", "error": "报告生成失败"}
            
    except Exception as e:
        logger.error(f"❌ AI报告生成测试失败: {e}", exc_info=True)
        return {"status": "error", "error": str(e)}


async def main():
    """
    主测试函数
    """
    logger.info("=" * 80)
    logger.info("🧪 RNA分析全流程自测程序")
    logger.info("=" * 80)
    
    # 设置路径
    project_root = Path(__file__).parent.parent
    test_data_dir = project_root / "test_data"
    upload_dir = os.getenv("UPLOAD_DIR", str(project_root / "uploads"))
    output_dir = project_root / "test_results" / f"rna_test_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    logger.info(f"📂 测试数据目录: {test_data_dir}")
    logger.info(f"📂 上传目录: {upload_dir}")
    logger.info(f"📂 输出目录: {output_dir}")
    
    # 查找10x数据目录
    tenx_dirs = find_10x_data_directories(test_data_dir)
    
    if not tenx_dirs:
        logger.error("❌ 未找到10x数据目录，请确保test_data文件夹中包含10x格式数据")
        return
    
    # 使用第一个找到的10x数据目录
    test_file_path = str(tenx_dirs[0])
    logger.info(f"📁 使用测试数据: {test_file_path}")
    
    # 测试1: 数据诊断报告生成
    diagnosis_result = await test_data_diagnosis(test_file_path, upload_dir)
    
    # 测试2: 工作流执行
    execution_result = await test_workflow_execution(test_file_path, str(output_dir), upload_dir)
    
    # 测试3: AI专家分析报告生成
    if execution_result.get("status") == "success":
        report_result = await test_ai_report_generation(execution_result)
    else:
        logger.warning("⚠️ 跳过AI报告生成测试（工作流执行失败）")
        report_result = {"status": "skipped"}
    
    # 汇总测试结果
    logger.info("=" * 80)
    logger.info("📊 测试结果汇总")
    logger.info("=" * 80)
    logger.info(f"1. 数据诊断报告生成: {diagnosis_result.get('status', 'unknown')}")
    logger.info(f"2. 工作流执行: {execution_result.get('status', 'unknown')}")
    logger.info(f"3. AI专家分析报告生成: {report_result.get('status', 'unknown')}")
    
    if report_result.get("status") == "success":
        logger.info(f"   - 生信专家版本: {report_result.get('is_bioinformatics', False)}")
        logger.info(f"   - 程序员日志版本: {report_result.get('is_programmer_log', False)}")
    
    # 保存测试结果到文件
    result_file = output_dir / "test_summary.txt"
    with open(result_file, 'w', encoding='utf-8') as f:
        f.write("RNA分析全流程测试结果汇总\n")
        f.write("=" * 80 + "\n\n")
        f.write(f"1. 数据诊断报告生成: {diagnosis_result.get('status', 'unknown')}\n")
        f.write(f"2. 工作流执行: {execution_result.get('status', 'unknown')}\n")
        f.write(f"3. AI专家分析报告生成: {report_result.get('status', 'unknown')}\n\n")
        
        if report_result.get("status") == "success":
            f.write(f"AI报告验证:\n")
            f.write(f"  - 生信专家版本: {report_result.get('is_bioinformatics', False)}\n")
            f.write(f"  - 程序员日志版本: {report_result.get('is_programmer_log', False)}\n\n")
            f.write(f"报告内容:\n{report_result.get('report', 'N/A')}\n")
    
    logger.info(f"✅ 测试结果已保存到: {result_file}")


if __name__ == "__main__":
    asyncio.run(main())
