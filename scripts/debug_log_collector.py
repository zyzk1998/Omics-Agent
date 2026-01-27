#!/usr/bin/env python3
"""
专业调试日志收集器

捕获以下关键过程的详细日志：
1. 有文件模式下规划阶段的LLM数据诊断过程
2. 代谢组分析过程所有工具的参数、输入输出
3. 执行结果回滚LLM生成分析报告的全部过程

使用方法：
    python3 scripts/debug_log_collector.py
"""

import os
import sys
import json
import logging
import asyncio
from pathlib import Path
from datetime import datetime
from typing import Dict, Any, List, Optional
from contextlib import contextmanager

# 添加项目路径
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

# 导入必要的模块
from dotenv import load_dotenv
load_dotenv(project_root / ".env")

from gibh_agent.core.file_inspector import FileInspector
from gibh_agent.core.orchestrator import AgentOrchestrator
from gibh_agent.core.executor import WorkflowExecutor
from gibh_agent.core.llm_client import LLMClient
from gibh_agent.core.prompt_manager import PromptManager, create_default_prompt_manager
from gibh_agent.agents.specialists.metabolomics_agent import MetabolomicsAgent


# ============================================
# 自定义日志处理器
# ============================================

class CategoryFilter(logging.Filter):
    """按类别过滤日志"""
    
    def __init__(self, category: str, allowed_modules: List[str]):
        super().__init__()
        self.category = category
        self.allowed_modules = allowed_modules
    
    def filter(self, record):
        """检查日志记录是否属于此类别"""
        module_name = record.name
        
        # 检查模块是否匹配
        for allowed in self.allowed_modules:
            if module_name.startswith(allowed):
                return True
        
        return False


class DetailedFileHandler(logging.Handler):
    """详细的文件日志处理器，记录函数调用、参数、返回值等"""
    
    def __init__(self, log_file: str, category: str):
        super().__init__()
        self.log_file = log_file
        self.category = category
        self.log_dir = project_root / "debug_logs"
        self.log_dir.mkdir(exist_ok=True)
        
        # 创建文件处理器
        file_path = self.log_dir / log_file
        self.file_handler = logging.FileHandler(file_path, encoding='utf-8')
        self.file_handler.setLevel(logging.DEBUG)
        
        # 详细格式
        formatter = logging.Formatter(
            '%(asctime)s | %(levelname)-8s | %(name)-40s | %(funcName)s:%(lineno)d | %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S.%f'
        )
        self.file_handler.setFormatter(formatter)
    
    def emit(self, record):
        """发送日志记录到文件"""
        try:
            self.file_handler.emit(record)
        except Exception:
            self.handleError(record)
    
    def close(self):
        """关闭文件处理器"""
        self.file_handler.close()
        super().close()


class LogCollector:
    """日志收集器，管理多个日志处理器"""
    
    def __init__(self):
        self.log_dir = project_root / "debug_logs"
        self.log_dir.mkdir(exist_ok=True)
        
        # 创建时间戳
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        self.timestamp = timestamp
        
        # 三个日志文件
        self.log_files = {
            "diagnosis": f"01_diagnosis_planning_{timestamp}.log",
            "execution": f"02_tool_execution_{timestamp}.log",
            "report": f"03_llm_report_generation_{timestamp}.log"
        }
        
        # 日志处理器字典
        self.handlers: Dict[str, DetailedFileHandler] = {}
        
        # 根日志记录器
        self.root_logger = logging.getLogger()
        self.original_level = self.root_logger.level
        
    def setup(self):
        """设置日志收集"""
        print("=" * 80)
        print("🔍 设置专业调试日志收集器")
        print("=" * 80)
        
        # 设置根日志级别为DEBUG
        self.root_logger.setLevel(logging.DEBUG)
        
        # 创建三个专门的日志处理器
        for category, log_file in self.log_files.items():
            handler = DetailedFileHandler(log_file, category)
            handler.setLevel(logging.DEBUG)
            self.handlers[category] = handler
            self.root_logger.addHandler(handler)
            print(f"✅ 已创建日志处理器: {category} -> {log_file}")
        
        # 添加控制台输出（INFO级别）
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(logging.INFO)
        console_formatter = logging.Formatter('%(asctime)s | %(levelname)-8s | %(message)s')
        console_handler.setFormatter(console_formatter)
        self.root_logger.addHandler(console_handler)
        
        print(f"\n📁 日志文件将保存到: {self.log_dir}")
        print("=" * 80)
    
    def add_module_filter(self, module_names: List[str], category: str):
        """为特定模块添加过滤器"""
        if category not in self.handlers:
            return
        
        handler = self.handlers[category]
        filter_obj = CategoryFilter(category, module_names)
        handler.addFilter(filter_obj)
        print(f"✅ 已为类别 '{category}' 添加模块过滤器: {module_names}")
    
    def cleanup(self):
        """清理日志处理器"""
        for handler in self.handlers.values():
            handler.close()
            self.root_logger.removeHandler(handler)
        
        # 恢复原始日志级别
        self.root_logger.setLevel(self.original_level)
        
        print("\n" + "=" * 80)
        print("✅ 日志收集完成")
        print("=" * 80)
        print(f"\n📁 日志文件位置:")
        for category, log_file in self.log_files.items():
            file_path = self.log_dir / log_file
            if file_path.exists():
                size = file_path.stat().st_size
                print(f"   {category:12s}: {file_path} ({size:,} bytes)")
        print("=" * 80)


# ============================================
# 函数包装器（用于捕获参数和返回值）
# ============================================

def log_function_call(logger_name: str, category: str):
    """装饰器：记录函数调用、参数和返回值"""
    def decorator(func):
        async def async_wrapper(*args, **kwargs):
            logger = logging.getLogger(logger_name)
            
            # 记录函数调用
            logger.info(f"🔵 [CALL] {func.__name__}()")
            logger.debug(f"   📥 Args: {args}")
            logger.debug(f"   📥 Kwargs: {kwargs}")
            
            try:
                # 执行函数
                result = await func(*args, **kwargs)
                
                # 记录返回值
                if isinstance(result, (str, dict, list)):
                    result_str = json.dumps(result, ensure_ascii=False, indent=2) if isinstance(result, dict) else str(result)
                    if len(result_str) > 1000:
                        result_str = result_str[:1000] + "... (truncated)"
                    logger.debug(f"   📤 Return: {result_str}")
                else:
                    logger.debug(f"   📤 Return: {type(result).__name__}")
                
                logger.info(f"✅ [SUCCESS] {func.__name__}()")
                return result
            except Exception as e:
                logger.error(f"❌ [ERROR] {func.__name__}(): {e}", exc_info=True)
                raise
        
        def sync_wrapper(*args, **kwargs):
            logger = logging.getLogger(logger_name)
            
            # 记录函数调用
            logger.info(f"🔵 [CALL] {func.__name__}()")
            logger.debug(f"   📥 Args: {args}")
            logger.debug(f"   📥 Kwargs: {kwargs}")
            
            try:
                # 执行函数
                result = func(*args, **kwargs)
                
                # 记录返回值
                if isinstance(result, (str, dict, list)):
                    result_str = json.dumps(result, ensure_ascii=False, indent=2) if isinstance(result, dict) else str(result)
                    if len(result_str) > 1000:
                        result_str = result_str[:1000] + "... (truncated)"
                    logger.debug(f"   📤 Return: {result_str}")
                else:
                    logger.debug(f"   📤 Return: {type(result).__name__}")
                
                logger.info(f"✅ [SUCCESS] {func.__name__}()")
                return result
            except Exception as e:
                logger.error(f"❌ [ERROR] {func.__name__}(): {e}", exc_info=True)
                raise
        
        if asyncio.iscoroutinefunction(func):
            return async_wrapper
        else:
            return sync_wrapper
    
    return decorator


# ============================================
# 测试工作流
# ============================================

async def test_metabolomics_workflow(log_collector: LogCollector):
    """测试代谢组学工作流，收集所有日志"""
    
    print("\n" + "=" * 80)
    print("🧪 开始测试代谢组学工作流")
    print("=" * 80)
    
    # 1. 初始化组件
    print("\n📋 Step 1: 初始化组件...")
    
    # LLM客户端
    base_url = os.getenv("LLM_BASE_URL", "http://localhost:8000/v1")
    api_key = os.getenv("LLM_API_KEY", "sk-test")
    llm_client = LLMClient(base_url=base_url, api_key=api_key)
    
    # Prompt Manager
    prompt_manager = create_default_prompt_manager()
    
    # Agent
    agent = MetabolomicsAgent(
        llm_client=llm_client,
        prompt_manager=prompt_manager
    )
    
    # File Inspector
    upload_dir = project_root / "test_data"
    upload_dir.mkdir(exist_ok=True)
    file_inspector = FileInspector(upload_dir=str(upload_dir))
    
    # Orchestrator
    orchestrator = AgentOrchestrator(agent=agent, upload_dir=str(upload_dir))
    
    # Executor
    executor = WorkflowExecutor(upload_dir=str(upload_dir))
    
    # 2. 查找测试文件
    print("\n📋 Step 2: 查找测试文件...")
    test_files = list(upload_dir.glob("*.csv"))
    if not test_files:
        print("❌ 未找到测试文件，请将CSV文件放入 test_data/ 目录")
        print("   提示: 可以使用示例数据或创建测试文件")
        return
    
    test_file = test_files[0]
    print(f"✅ 找到测试文件: {test_file}")
    
    # 记录文件信息
    diagnosis_logger = logging.getLogger("gibh_agent.core.file_inspector")
    diagnosis_logger.info(f"📁 测试文件路径: {test_file}")
    diagnosis_logger.info(f"📁 文件大小: {test_file.stat().st_size:,} bytes")
    
    # 3. 文件检查和诊断（规划阶段）
    print("\n" + "=" * 80)
    print("📊 Phase 1: 文件检查和LLM数据诊断（规划阶段）")
    print("=" * 80)
    
    diagnosis_logger = logging.getLogger("gibh_agent.core.file_inspector")
    diagnosis_logger.info("=" * 80)
    diagnosis_logger.info("PHASE 1: 文件检查和LLM数据诊断（规划阶段）")
    diagnosis_logger.info("=" * 80)
    
    file_metadata = None
    diagnosis = None
    
    try:
        # 文件检查
        diagnosis_logger.info(f"🔍 开始文件检查: {test_file}")
        file_metadata = file_inspector.inspect_file(str(test_file))
        diagnosis_logger.info(f"✅ 文件检查完成: {file_metadata.get('status', 'unknown')}")
        
        # 详细记录文件元数据
        if file_metadata:
            diagnosis_logger.info("📊 文件元数据摘要:")
            diagnosis_logger.info(f"   - 状态: {file_metadata.get('status')}")
            diagnosis_logger.info(f"   - 文件类型: {file_metadata.get('file_type')}")
            diagnosis_logger.info(f"   - 行数: {file_metadata.get('shape', {}).get('rows', 'N/A')}")
            diagnosis_logger.info(f"   - 列数: {file_metadata.get('shape', {}).get('cols', 'N/A')}")
            diagnosis_logger.info(f"   - 缺失率: {file_metadata.get('missing_rate', 'N/A')}%")
            
            # 完整元数据（DEBUG级别）
            diagnosis_logger.debug(f"📋 完整文件元数据:\n{json.dumps(file_metadata, ensure_ascii=False, indent=2)}")
        
        # LLM诊断
        if file_metadata and file_metadata.get("status") == "success":
            diagnosis_logger = logging.getLogger("gibh_agent.agents.base_agent")
            diagnosis_logger.info("=" * 80)
            diagnosis_logger.info("🔍 开始LLM数据诊断...")
            diagnosis_logger.info("=" * 80)
            
            # 加载数据预览
            import pandas as pd
            df = pd.read_csv(test_file, nrows=100)
            diagnosis_logger.info(f"📊 数据预览加载完成: {len(df)} 行 × {len(df.columns)} 列")
            
            # 记录LLM调用参数
            diagnosis_logger.info("📞 LLM调用参数:")
            diagnosis_logger.info(f"   - omics_type: Metabolomics")
            diagnosis_logger.info(f"   - file_metadata: {json.dumps(file_metadata, ensure_ascii=False)[:500]}...")
            
            diagnosis = await agent._perform_data_diagnosis(
                file_metadata=file_metadata,
                omics_type="Metabolomics",
                dataframe=df,
                system_instruction=None
            )
            
            diagnosis_logger.info(f"✅ LLM诊断完成，长度: {len(diagnosis) if diagnosis else 0}")
            if diagnosis:
                diagnosis_logger.info("📝 诊断报告预览:")
                diagnosis_logger.info(f"{diagnosis[:1000]}...")
                # 完整报告（DEBUG级别）
                diagnosis_logger.debug(f"📝 完整诊断报告:\n{diagnosis}")
    
    except Exception as e:
        diagnosis_logger.error(f"❌ Phase 1 失败: {e}", exc_info=True)
        import traceback
        diagnosis_logger.error(f"❌ 完整错误堆栈:\n{traceback.format_exc()}")
    
    # 4. 工作流执行（工具调用）
    print("\n" + "=" * 80)
    print("🔧 Phase 2: 工作流执行（工具参数和输入输出）")
    print("=" * 80)
    
    execution_logger = logging.getLogger("gibh_agent.core.executor")
    execution_logger.info("=" * 80)
    execution_logger.info("PHASE 2: 工作流执行（工具参数和输入输出）")
    execution_logger.info("=" * 80)
    
    results = None
    
    try:
        # 创建简单的工作流配置
        workflow_config = {
            "workflow_name": "代谢组学分析测试",
            "steps": [
                {
                    "id": "preprocess",
                    "name": "数据预处理",
                    "tool_id": "metabolomics_preprocess_data",
                    "parameters": {
                        "file_path": str(test_file),
                        "log_transform": True,
                        "standardize": True
                    }
                }
            ]
        }
        
        execution_logger.info("📋 工作流配置:")
        execution_logger.info(f"{json.dumps(workflow_config, ensure_ascii=False, indent=2)}")
        
        # 执行工作流（同步方法）
        execution_logger.info("🚀 开始执行工作流...")
        results = executor.execute_workflow(
            workflow_data=workflow_config,
            file_paths=[str(test_file)],
            output_dir=None,
            agent=agent
        )
        
        execution_logger.info(f"✅ 工作流执行完成")
        
        # 详细记录执行结果
        if results:
            execution_logger.info("📊 执行结果摘要:")
            execution_logger.info(f"   - 状态: {results.get('status', 'unknown')}")
            execution_logger.info(f"   - 步骤数: {len(results.get('steps_details', []))}")
            
            # 记录每个步骤的详细信息
            steps_details = results.get('steps_details', [])
            for i, step_detail in enumerate(steps_details, 1):
                execution_logger.info(f"\n📋 步骤 {i}:")
                execution_logger.info(f"   - ID: {step_detail.get('step_id', 'N/A')}")
                execution_logger.info(f"   - 名称: {step_detail.get('step_name', 'N/A')}")
                execution_logger.info(f"   - 状态: {step_detail.get('status', 'N/A')}")
                
                # 记录工具调用参数
                step_result = step_detail.get('step_result', {})
                if step_result:
                    execution_logger.info(f"   - 工具输出状态: {step_result.get('status', 'N/A')}")
                    execution_logger.debug(f"   - 完整步骤结果:\n{json.dumps(step_result, ensure_ascii=False, indent=6)}")
            
            # 完整结果（DEBUG级别）
            execution_logger.debug(f"📋 完整执行结果:\n{json.dumps(results, ensure_ascii=False, indent=2)}")
    
    except Exception as e:
        execution_logger.error(f"❌ Phase 2 失败: {e}", exc_info=True)
        import traceback
        execution_logger.error(f"❌ 完整错误堆栈:\n{traceback.format_exc()}")
    
    # 5. LLM生成分析报告
    print("\n" + "=" * 80)
    print("📝 Phase 3: LLM生成分析报告")
    print("=" * 80)
    
    report_logger = logging.getLogger("gibh_agent.agents.base_agent")
    report_logger.info("=" * 80)
    report_logger.info("PHASE 3: LLM生成分析报告（执行结果回滚）")
    report_logger.info("=" * 80)
    
    summary = None
    
    try:
        # 从执行结果中提取步骤结果
        steps_results = []
        if results and results.get('steps_details'):
            for step_detail in results.get('steps_details', []):
                step_result = step_detail.get('step_result', {})
                if step_result:
                    steps_results.append({
                        "step_name": step_detail.get('step_name', 'Unknown'),
                        "status": step_result.get('status', 'unknown'),
                        "data": step_result.get('data', {})
                    })
        
        # 如果没有执行结果，使用模拟数据
        if not steps_results:
            report_logger.warning("⚠️ 未找到执行结果，使用模拟数据")
            steps_results = [
                {
                    "step_name": "数据预处理",
                    "status": "success",
                    "data": {
                        "summary": {
                            "n_samples": 10,
                            "n_features": 50
                        }
                    }
                }
            ]
        
        report_logger.info(f"📊 步骤结果数量: {len(steps_results)}")
        report_logger.debug(f"📊 步骤结果详情:\n{json.dumps(steps_results, ensure_ascii=False, indent=2)}")
        
        # 记录LLM调用参数
        report_logger.info("📞 LLM调用参数:")
        report_logger.info(f"   - omics_type: Metabolomics")
        report_logger.info(f"   - workflow_name: 代谢组学分析测试")
        report_logger.info(f"   - steps_results数量: {len(steps_results)}")
        
        # 生成报告
        report_logger.info("🚀 开始调用LLM生成分析报告...")
        summary = await agent._generate_analysis_summary(
            steps_results=steps_results,
            omics_type="Metabolomics",
            workflow_name="代谢组学分析测试",
            summary_context=None,
            output_dir=None
        )
        
        report_logger.info(f"✅ LLM报告生成完成，长度: {len(summary) if summary else 0}")
        if summary:
            report_logger.info("📝 报告预览:")
            report_logger.info(f"{summary[:1000]}...")
            # 完整报告（DEBUG级别）
            report_logger.debug(f"📝 完整报告:\n{summary}")
        else:
            report_logger.warning("⚠️ LLM报告为空")
    
    except Exception as e:
        report_logger.error(f"❌ Phase 3 失败: {e}", exc_info=True)
        import traceback
        report_logger.error(f"❌ 完整错误堆栈:\n{traceback.format_exc()}")
    
    print("\n" + "=" * 80)
    print("✅ 测试完成")
    print("=" * 80)


async def main():
    """主函数"""
    # 创建日志收集器
    log_collector = LogCollector()
    log_collector.setup()
    
    # 添加模块过滤器
    log_collector.add_module_filter([
        "gibh_agent.agents.base_agent",
        "gibh_agent.core.file_inspector",
        "gibh_agent.core.orchestrator",
        "gibh_agent.core.llm_client"
    ], "diagnosis")
    
    log_collector.add_module_filter([
        "gibh_agent.core.executor",
        "gibh_agent.tools",
        "gibh_agent.core.tool_registry"
    ], "execution")
    
    log_collector.add_module_filter([
        "gibh_agent.agents.base_agent",
        "gibh_agent.core.llm_client"
    ], "report")
    
    try:
        # 运行测试
        await test_metabolomics_workflow(log_collector)
    
    finally:
        # 清理
        log_collector.cleanup()


if __name__ == "__main__":
    asyncio.run(main())
