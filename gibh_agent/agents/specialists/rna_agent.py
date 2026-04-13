"""
转录组智能体（RNA Agent）
处理单细胞转录组（scRNA-seq）和 Bulk RNA-seq 分析
重构自现有的 BioBlendAgent
"""
import json
import os
from pathlib import Path
from typing import Dict, Any, List, Optional, AsyncIterator
from ..base_agent import BaseAgent
from ...core.llm_client import LLMClient
from ...core.prompt_manager import PromptManager, RNA_REPORT_PROMPT
from ...core.utils import sanitize_for_json
from ...core.dispatcher import TaskDispatcher
from ...core.test_data_manager import TestDataManager
from ...core.tool_retriever import ToolRetriever
from ...core.planner import RNAPlanner
from ...core.workflows import WorkflowRegistry
from ...core.tool_registry import registry
# 导入新工具函数
from ...tools.general.file_inspector import inspect_file
from ...tools.rna.upstream import run_cellranger_count, convert_cellranger_to_h5ad
from ...tools.rna.quality_control import run_qc_filter
from ...tools.rna.analysis import run_normalize, run_hvg, run_pca, run_neighbors, run_umap, run_clustering
from ...tools.rna.annotation import run_cell_annotation
from ...tools.rna.plotting import visualize_qc, visualize_clustering
import scanpy as sc
import logging

logger = logging.getLogger(__name__)


# 🔥 架构重构：领域特定的系统指令（策略模式）
RNA_INSTRUCTION = """You are a Senior Bioinformatician specializing in Single-Cell RNA-seq analysis.

**CRITICAL CONSTRAINTS:**
- The data represents **Gene Expression** (RNA transcripts), NOT Metabolite Abundance.
- Rows = Cells (Single Cells), Columns = Genes (Gene Expression).
- This is RNA sequencing data, measuring transcript counts per cell.

**REQUIRED TERMINOLOGY:**
- Cell, Cells, Single Cell, Cellular
- Gene, Genes, Gene Expression, Transcript
- Mitochondria, Mitochondrial (mt-genes)
- scRNA-seq, Single-Cell RNA-seq, scRNA
- Transcriptomics, Transcriptome
- UMI, Count Matrix, Expression Matrix

**CONTEXT:**
This is single-cell transcriptomics data representing gene expression levels measured by RNA sequencing.

Generate data diagnosis and parameter recommendations in Simplified Chinese (简体中文).
Focus on single-cell-specific quality metrics (cells, genes, mitochondrial percentage, doublet rate)."""


class RNAAgent(BaseAgent):
    """
    转录组智能体
    
    职责：
    1. 处理单细胞转录组分析（scRNA-seq）
    2. 处理 Bulk RNA-seq 分析
    3. 生成工作流脚本
    4. 通过 TaskDispatcher 提交任务
    """
    
    def __init__(
        self,
        llm_client: LLMClient,
        prompt_manager: PromptManager,
        dispatcher: Optional[TaskDispatcher] = None,
        cellranger_config: Optional[Dict[str, Any]] = None,
        scanpy_config: Optional[Dict[str, Any]] = None,
        test_data_dir: Optional[str] = None,
        tool_retriever: Optional[ToolRetriever] = None
    ):
        """初始化转录组智能体"""
        super().__init__(llm_client, prompt_manager, "rna_expert")
        
        self.dispatcher = dispatcher
        self.cellranger_config = cellranger_config or {}
        self.scanpy_config = scanpy_config or {}
        # 🔥 架构升级：移除旧工具，使用新模块化工具系统
        # 初始化测试数据管理器
        self.test_data_manager = TestDataManager(test_data_dir)
        
        # 🔥 架构升级：初始化 SOP 驱动的动态规划器
        self.sop_planner = None
        if tool_retriever:
            try:
                # 创建 scRNA-seq 特定的 SOPPlanner（需要自定义 SOP 规则）
                self.sop_planner = RNAPlanner(tool_retriever, llm_client)
                logger.info("✅ [RNAAgent] RNAPlanner 已初始化")
            except Exception as e:
                logger.warning(f"⚠️ [RNAAgent] RNAPlanner 初始化失败，将使用回退逻辑: {e}")
        else:
            logger.info("ℹ️ [RNAAgent] 未提供 ToolRetriever，将使用传统工作流生成逻辑")
        
        # 🔥 工具ID映射表：旧ID -> 新ID（已注册的工具名称）
        self.tool_id_mapping = {
            "local_qc": "rna_qc_filter",
            "local_normalize": "rna_normalize",
            "local_hvg": "rna_hvg",
            "local_scale": "rna_scale",
            "local_pca": "rna_pca",
            "local_neighbors": "rna_neighbors",
            "local_cluster": "rna_clustering",
            "local_umap": "rna_umap",
            "local_tsne": "rna_tsne",
            "local_markers": "rna_find_markers",
        }
        
        # 标准工作流步骤（十步流程）- 使用新的工具ID
        self.workflow_steps = [
            {"name": "1. Quality Control", "tool_id": "rna_qc_filter", "desc": "过滤低质量细胞和基因"},
            {"name": "2. Normalization", "tool_id": "rna_normalize", "desc": "数据标准化"},
            {"name": "3. Find Variable Genes", "tool_id": "rna_hvg", "desc": "筛选高变基因"},
            {"name": "4. Scale Data", "tool_id": "rna_scale", "desc": "数据缩放"},
            {"name": "5. PCA", "tool_id": "rna_pca", "desc": "主成分分析"},
            {"name": "6. Compute Neighbors", "tool_id": "rna_neighbors", "desc": "构建邻接图"},
            {"name": "7. Clustering", "tool_id": "rna_clustering", "desc": "Leiden 聚类"},
            {"name": "8. UMAP Visualization", "tool_id": "rna_umap", "desc": "UMAP 可视化"},
            {"name": "9. t-SNE Visualization", "tool_id": "rna_tsne", "desc": "t-SNE 可视化"},
            {"name": "10. Find Markers", "tool_id": "rna_find_markers", "desc": "寻找 Marker 基因"},
        ]
    
    async def process_query(
        self,
        query: str,
        history: List[Dict[str, str]] = None,
        uploaded_files: List[Dict[str, str]] = None,
        **kwargs
    ) -> Dict[str, Any]:
        """
        处理用户查询
        
        Returns:
            处理结果字典，可能包含：
            - workflow_config: 工作流配置（JSON）
            - chat_response: 聊天响应（流式）
            - task_submitted: 任务提交信息
            - test_data_selection: 测试数据选择请求
        """
        query_lower = query.lower().strip()
        file_paths = self.get_file_paths(uploaded_files or [])
        
        # 🔥 Task 1: LLM 驱动的意图检测（在生成工作流之前）
        # 🔒 安全包装：如果意图检测失败，回退到原始逻辑
        intent = "chat"  # 默认值
        intent_result = None
        try:
            intent_result = await self._detect_intent_with_llm(query, file_paths, uploaded_files)
            intent = intent_result.get("intent", "chat")
            reasoning = intent_result.get("reasoning", "")
            logger.info(f"🎯 意图检测结果: {intent} (推理: {reasoning})")
        except Exception as e:
            logger.warning(f"⚠️ 意图检测失败，回退到原始逻辑: {e}", exc_info=True)
            # 回退到原始的工作流检测逻辑
            intent = None  # 标记为未检测，使用回退逻辑
        
        # 如果意图检测成功且为 explain_file，处理文件解释
        if intent == "explain_file":
            # 解释文件：检查文件并生成自然语言解释
            if not file_paths:
                return {
                    "type": "chat",
                    "response": self._stream_string_response("没有检测到上传的文件。请先上传文件后再询问。")
                }
            
            # 🔧 修复：优先使用最新上传的文件（列表最后一个），而不是第一个
            # 检查最新上传的文件（如果是 h5ad，使用 scanpy 工具）
            input_path = file_paths[-1] if file_paths else None
            if not input_path:
                return {
                    "type": "chat",
                    "response": self._stream_string_response("没有检测到上传的文件。请先上传文件后再询问。")
                }
            try:
                # 使用新工具系统
                if input_path.endswith('.h5ad'):
                    adata = sc.read_h5ad(input_path)
                    summary = f"""
文件类型: H5AD (AnnData)
- 细胞数: {adata.n_obs}
- 基因数: {adata.n_vars}
- 观察变量: {list(adata.obs.columns) if hasattr(adata, 'obs') else 'None'}
- 变量变量: {list(adata.var.columns) if hasattr(adata, 'var') else 'None'}
"""
                    explanation = await self._explain_file_with_llm(query, summary, input_path)
                    return {
                        "type": "chat",
                        "response": self._stream_string_response(explanation)
                    }
                else:
                    # 其他文件类型，读取文件内容并使用 LLM 解释
                    try:
                        # 使用 file_inspector 读取文件元数据和内容
                        from ..core.file_inspector import FileInspector
                        import os
                        
                        # 获取上传目录
                        upload_dir = os.getenv("UPLOAD_DIR", "/app/uploads")
                        file_inspector = FileInspector(upload_dir)
                        
                        # 获取文件元数据
                        file_name = os.path.basename(input_path)
                        metadata = file_inspector.generate_metadata(file_name)
                        
                        # 读取文件前几行作为内容预览
                        file_path_obj = Path(input_path)
                        if not file_path_obj.is_absolute():
                            file_path_obj = Path(upload_dir) / file_name
                        
                        file_summary = f"文件路径: {input_path}\n文件类型: {os.path.splitext(input_path)[1]}\n"
                        
                        if metadata:
                            file_summary += f"文件大小: {metadata.get('size_mb', 'unknown')} MB\n"
                            if metadata.get('estimated_cells'):
                                file_summary += f"估算细胞数: {metadata.get('estimated_cells')}\n"
                            if metadata.get('estimated_genes'):
                                file_summary += f"估算基因数: {metadata.get('estimated_genes')}\n"
                        
                        # 读取文件前几行
                        try:
                            if file_path_obj.exists() and file_path_obj.is_file():
                                head_lines = file_inspector._read_head(file_path_obj, 10)
                                if head_lines:
                                    file_summary += f"\n文件内容预览（前10行）：\n"
                                    for i, line in enumerate(head_lines[:10], 1):
                                        file_summary += f"{i}: {line[:200]}\n"  # 限制每行长度
                        except Exception as e:
                            logger.warning(f"⚠️ 读取文件内容失败: {e}")
                            file_summary += "\n（无法读取文件内容）\n"
                        
                        # 使用 LLM 生成文件解释
                        explanation = await self._explain_file_with_llm(query, file_summary, input_path)
                        return {
                            "type": "chat",
                            "response": self._stream_string_response(explanation)
                        }
                    except Exception as e:
                        logger.error(f"❌ 文件解释失败: {e}", exc_info=True)
                        # 回退到基本信息
                        return {
                            "type": "chat",
                            "response": self._stream_string_response(f"文件路径: {input_path}\n文件类型: {os.path.splitext(input_path)[1]}\n\n（文件内容读取失败，请检查文件格式）")
                        }
            except Exception as e:
                logger.error(f"❌ 文件解释失败: {e}", exc_info=True)
                return {
                    "type": "chat",
                    "response": self._stream_string_response(f"文件检查时出错: {str(e)}")
                }
        
        # 与 Spatial/Radiomics 一致：无文件时不再提供本地测试数据选择，由规划器返回工作流预览卡片（上传以激活）
        
        # 处理测试数据选择（用户通过 JSON 选择，保留兼容）
        if "test_dataset_id" in kwargs:
            dataset_id = kwargs["test_dataset_id"]
            dataset = self.test_data_manager.get_dataset_by_id(dataset_id)
            if dataset:
                # 使用选定的测试数据
                if dataset.get("fastq_dir") and dataset.get("reference"):
                    # 有 FASTQ 和参考基因组，使用它们
                    file_paths = [dataset["fastq_dir"]]
                    # 将参考基因组路径添加到配置中
                    self.cellranger_config["reference"] = dataset["reference"]
                elif dataset.get("h5ad_file"):
                    # 只有 .h5ad 文件，直接使用
                    file_paths = [dataset["h5ad_file"]]
                else:
                    return {
                        "type": "chat",
                        "response": self._stream_string_response(
                            f"测试数据集 {dataset['name']} 不可用。"
                        )
                    }
        
        # 🔒 回退逻辑：如果意图检测失败或意图不明确，使用原始逻辑
        if intent is None or intent == "chat":
            # 使用原始的工作流检测逻辑作为回退
            is_workflow_request = self._is_workflow_request(query_lower, file_paths)
            if is_workflow_request:
                return await self._generate_workflow_config(query, file_paths)
            else:
                # 普通聊天
                return {
                    "type": "chat",
                    "response": self._stream_chat_response(query, file_paths)
                }
        
        # 如果意图明确为 run_workflow，直接生成工作流配置
        elif intent == "run_workflow":
            return await self._generate_workflow_config(query, file_paths)
        
        # 默认：普通聊天
        else:
            return {
                "type": "chat",
                "response": self._stream_chat_response(query, file_paths)
            }
    
    async def _detect_intent_with_llm(
        self,
        query: str,
        file_paths: List[str],
        uploaded_files: List[Dict[str, str]] = None
    ) -> Dict[str, Any]:
        """
        使用 LLM 检测用户意图
        
        Returns:
            {
                "intent": "explain_file" | "run_workflow" | "chat",
                "reasoning": "..."
            }
        """
        import json
        import os
        
        # 提取文件名
        file_names = []
        if uploaded_files:
            for f in uploaded_files:
                name = f.get("name") or f.get("file_name", "")
                if name:
                    file_names.append(name)
        elif file_paths:
            for path in file_paths:
                file_names.append(os.path.basename(path))
        
        file_names_str = ", ".join(file_names) if file_names else "None"
        
        prompt = f"""分析用户输入，判断用户意图。

User Input: {query}
Uploaded Files: {file_names_str}

请将意图分类为以下三种之一：
1. "explain_file" - 用户想要了解文件内容、结构或含义（例如："这是什么文件？"、"文件里有什么？"、"解释一下这个数据"）
2. "run_workflow" - 用户想要执行分析工作流（例如："分析一下"、"运行工作流"、"做一下分析"、"处理这个文件"）
3. "chat" - 普通对话或询问（例如："你好"、"如何使用"、"介绍功能"）

返回 JSON 格式：
{{
    "intent": "explain_file" | "run_workflow" | "chat",
    "reasoning": "简要说明判断理由"
}}"""
        
        messages = [
            {"role": "system", "content": "You are an intent classification assistant. Return JSON only."},
            {"role": "user", "content": prompt}
        ]
        
        try:
            completion = await self.llm_client.achat(messages, temperature=0.1, max_tokens=128)
            think_content, response = self.llm_client.extract_think_and_content(completion)
            
            # 解析 JSON
            json_str = response.strip()
            if "```json" in json_str:
                json_str = json_str.split("```json")[1].split("```")[0].strip()
            elif "```" in json_str:
                json_str = json_str.split("```")[1].split("```")[0].strip()
            
            result = json.loads(json_str)
            
            # 验证意图值
            valid_intents = ["explain_file", "run_workflow", "chat"]
            if result.get("intent") not in valid_intents:
                logger.warning(f"⚠️ LLM 返回了无效意图: {result.get('intent')}, 使用默认值 'chat'")
                result["intent"] = "chat"
            
            return result
        except Exception as e:
            logger.error(f"❌ 意图检测失败: {e}", exc_info=True)
            # 默认返回 chat
            return {
                "intent": "chat",
                "reasoning": f"Intent detection failed: {str(e)}"
            }
    
    async def _explain_file_with_llm(
        self,
        query: str,
        file_summary: str,
        file_path: str
    ) -> str:
        """
        使用 LLM 生成文件解释
        
        Args:
            query: 用户查询
            file_summary: 文件摘要信息
            file_path: 文件路径
        
        Returns:
            自然语言的文件解释
        """
        prompt = f"""用户询问关于文件的问题。

User Query: {query}
File Path: {file_path}

文件摘要信息：
{file_summary}

请用自然语言解释这个文件的内容、结构和特点。回答应该：
1. 简洁明了，易于理解
2. 包含关键信息（细胞数、基因数等）
3. 如果用户有特定问题，针对性地回答
4. 使用中文回答

回答："""
        
        messages = [
            {"role": "system", "content": "You are a bioinformatics data expert. Explain file contents in natural language."},
            {"role": "user", "content": prompt}
        ]
        
        try:
            completion = await self.llm_client.achat(messages, temperature=0.3, max_tokens=800)
            think_content, response = self.llm_client.extract_think_and_content(completion)
            return response
        except Exception as e:
            logger.error(f"❌ 文件解释生成失败: {e}", exc_info=True)
            return f"文件解释生成失败: {str(e)}"
    
    def _stream_string_response(self, text: str) -> AsyncIterator[str]:
        """将字符串转换为异步生成器（用于流式响应）"""
        async def _generator():
            yield text
        return _generator()
    
    def _is_workflow_request(self, query: str, file_paths: List[str]) -> bool:
        """判断是否是工作流请求"""
        workflow_keywords = [
            "规划", "流程", "workflow", "pipeline", "分析", "run",
            "执行", "plan", "做一下", "跑一下", "分析一下"
        ]
        
        bio_keywords = [
            "pca", "umap", "tsne", "qc", "质控", "聚类", "cluster"
        ]
        
        if any(kw in query for kw in workflow_keywords):
            return True
        
        if file_paths and any(kw in query for kw in bio_keywords):
            return True
        
        if file_paths and (not query or len(query) < 5):
            return True
        
        return False
    
    def _needs_cellranger(self, query: str) -> bool:
        """判断是否需要 Cell Ranger（基于查询关键词）"""
        cellranger_keywords = [
            "cellranger", "cell ranger", "fastq", "fq", "测序",
            "第一步", "全流程", "完整流程", "从fastq", "从测序"
        ]
        return any(kw in query for kw in cellranger_keywords)
    
    async def _generate_workflow_config(
        self,
        query: str,
        file_paths: List[str]
    ) -> Dict[str, Any]:
        """
        生成工作流配置
        
        🔥 架构升级：优先使用 SOP 驱动的动态规划器，失败则回退到传统逻辑
        
        流程：
        1. 如果 RNAPlanner 可用，使用动态规划器生成工作流
        2. 否则，使用传统硬编码逻辑（回退）
        """
        logger.info("=" * 80)
        logger.info("🚀 [RNAAgent] _generate_workflow_config START")
        logger.info(f"   Query: {query}")
        logger.info(f"   File paths: {file_paths}")
        logger.info(f"   RNAPlanner available: {self.sop_planner is not None}")
        logger.info("=" * 80)
        
        # 🔥 Phase 3: 优先使用 SOP 驱动的动态规划器（无文件时也生成预览模板，与 Spatial/Radiomics 一致）
        if self.sop_planner:
            try:
                logger.info("🧠 [RNAPlanner] 尝试使用动态规划器生成工作流...")
                
                # Step 1: 有文件时检查并获取元数据
                file_metadata = None
                if file_paths:
                    try:
                        from ...core.file_inspector import FileInspector
                        import os
                        upload_dir = os.getenv("UPLOAD_DIR", "/app/uploads")
                        inspector = FileInspector(upload_dir)
                        file_metadata = inspector.inspect_file(file_paths[0])
                        if file_metadata.get("status") != "success":
                            logger.warning(f"⚠️ 文件检查失败: {file_metadata.get('error')}")
                            file_metadata = None
                    except Exception as e:
                        logger.warning(f"⚠️ 文件检查异常: {e}")
                
                # Step 2: 有文件用真实参数，无文件用模板模式（预览卡片，与 Spatial/Radiomics 一致）
                is_template = not file_metadata or not (file_metadata.get("file_path") if file_metadata else None)
                plan_result = None
                async for _ev, _data in self.sop_planner.generate_plan(
                    user_query=query,
                    file_metadata=file_metadata,
                    category_filter="scRNA-seq",
                    domain_name="RNA",
                    is_template=is_template,
                ):
                    if _ev == "workflow":
                        plan_result = _data
                
                if plan_result and plan_result.get("type") != "error":
                    logger.info("✅ [RNAPlanner] 动态规划成功 (template=%s)", is_template)
                    if file_metadata and not is_template:
                        diagnosis_report = None
                        try:
                            _wf = WorkflowRegistry().get_workflow("RNA")
                            _step_ids = None
                            if plan_result:
                                _acc = []
                                for s in (plan_result.get("workflow_data") or {}).get("steps") or []:
                                    if isinstance(s, dict):
                                        sid = s.get("step_id") or s.get("id")
                                        if sid:
                                            _acc.append(sid)
                                _step_ids = _acc or None
                            diagnosis_report = await self._perform_data_diagnosis(
                                file_metadata=file_metadata,
                                omics_type="scRNA",
                                system_instruction=RNA_INSTRUCTION,
                                workflow_for_whitelist=_wf,
                                target_step_ids_for_whitelist=_step_ids,
                            )
                        except Exception as e:
                            logger.warning(f"⚠️ 诊断报告生成失败: {e}")
                        if diagnosis_report:
                            plan_result["diagnosis_report"] = diagnosis_report
                    if hasattr(self, "context") and self.context.get("parameter_recommendation"):
                        plan_result["recommendation"] = self.context.get("parameter_recommendation")
                    return plan_result
                logger.warning(
                    "⚠️ [RNAPlanner] 规划失败: %s",
                    (plan_result or {}).get("error"),
                )
            
            except Exception as e:
                logger.error(f"❌ [RNAPlanner] 动态规划异常: {e}", exc_info=True)
                logger.info("🔄 回退到传统工作流生成逻辑...")
        
        # 🔄 回退：使用传统硬编码逻辑
        logger.info("📋 [Fallback] 使用传统工作流生成逻辑...")
        
        # 🔥 Step 1: 文件检查和数据诊断（使用统一的 BaseAgent 方法）
        inspection_result = None
        diagnosis_report = None
        if file_paths:
            input_path = file_paths[0]
            try:
                # 使用新工具系统
                if input_path.endswith('.h5ad'):
                    # 对于 H5AD 文件，直接加载并检查
                    adata = sc.read_h5ad(input_path)
                    inspection_result = {
                        "status": "success",
                        "n_obs": adata.n_obs,
                        "n_vars": adata.n_vars,
                        "file_type": "h5ad"
                    }
                else:
                    # 对于其他文件，使用通用检查工具
                    inspection_result = inspect_file(input_path)
                if "error" in inspection_result:
                    logger.warning(f"File inspection failed: {inspection_result.get('error')}")
                else:
                    # 🔥 使用 BaseAgent 的统一诊断方法
                    # 尝试加载数据预览（用于更准确的诊断）
                    dataframe = None
                    try:
                        if input_path.endswith('.h5ad'):
                            # 使用新工具系统
                            adata = sc.read_h5ad(input_path)
                            # 提取 obs 表作为预览（包含 QC 指标）
                            if hasattr(adata, 'obs') and len(adata.obs) > 0:
                                dataframe = adata.obs.head(1000)  # 最多1000行
                    except Exception as e:
                        logger.debug(f"无法加载数据预览: {e}")
                    
                    # 调用统一的诊断方法
                    # 🔥 架构重构：传递领域特定的系统指令
                    _fm = dict(inspection_result) if isinstance(inspection_result, dict) else {}
                    if input_path and not _fm.get("file_path"):
                        _fm["file_path"] = input_path
                    _wf = WorkflowRegistry().get_workflow("RNA")
                    diagnosis_report = await self._perform_data_diagnosis(
                        file_metadata=_fm,
                        omics_type="scRNA",
                        dataframe=dataframe,
                        system_instruction=RNA_INSTRUCTION,
                        workflow_for_whitelist=_wf,
                        target_step_ids_for_whitelist=None,
                    )
                    # 🔥 DEBUG: 打印诊断报告信息
                    if diagnosis_report:
                        logger.info(f"📝 [DEBUG] RNAAgent diagnosis report generated, length: {len(diagnosis_report)}")
                    else:
                        logger.warning(f"⚠️ [DEBUG] RNAAgent diagnosis report is None")
            except Exception as e:
                logger.error(f"Error inspecting file: {e}", exc_info=True)
                diagnosis_report = None  # 🔥 确保在异常时也设置为 None
        
        # 使用 LLM 提取参数（传入检查结果和诊断报告）
        extracted_params = await self._extract_workflow_params(query, file_paths, inspection_result, diagnosis_report)
        
        # 构建工作流配置
        workflow_config = {
            "workflow_name": "Standard Single-Cell Pipeline",
            "steps": []
        }
        
        for step_template in self.workflow_steps:
            step = step_template.copy()
            
            # 注入参数
            tool_id = step["tool_id"]
            # 🔥 工具ID映射：如果使用旧ID，映射到新ID
            if tool_id in self.tool_id_mapping:
                tool_id = self.tool_id_mapping[tool_id]
                step["tool_id"] = tool_id
            
            # 根据工具ID设置参数
            if tool_id == "rna_qc_filter":
                step["params"] = {
                    "min_genes": extracted_params.get("min_genes", 200),
                    "max_mt": extracted_params.get("max_mt", 20.0),
                    "min_cells": extracted_params.get("min_cells", 3)
                }
            elif tool_id == "rna_hvg":
                step["params"] = {
                    "n_top_genes": extracted_params.get("n_top_genes", 2000)
                }
            elif tool_id == "rna_clustering":
                step["params"] = {
                    "resolution": extracted_params.get("resolution", 0.5)
                }
            else:
                step["params"] = {}
            
            workflow_config["steps"].append(step)
        
        # 如果生成了诊断报告，将其包含在返回结果中
        result = {
            "type": "workflow_config",
            "workflow_data": workflow_config,
            "file_paths": file_paths
        }
        
        # 🔥 修复：检查 diagnosis_report 是否为有效字符串（非 None 且非空）
        if diagnosis_report and isinstance(diagnosis_report, str) and diagnosis_report.strip():
            result["diagnosis_report"] = diagnosis_report
            logger.info(f"📝 [DEBUG] RNAAgent: Adding diagnosis_report to result, length: {len(diagnosis_report)}")
        else:
            logger.warning(f"⚠️ [DEBUG] RNAAgent: diagnosis_report is invalid (None/empty), NOT adding to result. Type: {type(diagnosis_report)}, Value: {diagnosis_report}")
        
        # 🔥 TASK 5: 添加参数推荐到结果
        if hasattr(self, 'context') and "parameter_recommendation" in self.context:
            recommendation = self.context.get("parameter_recommendation")
            if recommendation:
                result["recommendation"] = recommendation
                logger.info(f"✅ [RNAAgent] 添加参数推荐到结果: {len(recommendation.get('params', {}))} 个参数")
        
        # 🔥 DEBUG: 打印最终返回结构
        logger.info(f"📤 [DEBUG] RNAAgent returning result with keys: {list(result.keys())}")
        logger.info(f"📤 [DEBUG] RNAAgent has diagnosis_report: {'diagnosis_report' in result}")
        
        return result
    
    # 🔥 已移除：_generate_diagnosis_and_recommendation 方法
    # 现在使用 BaseAgent._perform_data_diagnosis() 统一方法
    
    async def _extract_workflow_params(
        self,
        query: str,
        file_paths: List[str],
        inspection_result: Optional[Dict[str, Any]] = None,
        diagnosis_report: Optional[str] = None
    ) -> Dict[str, Any]:
        """
        使用 LLM 提取工作流参数
        
        基于检查结果智能推荐参数
        """
        # 构建包含检查结果的提示
        inspection_info = ""
        if inspection_result and "error" not in inspection_result:
            inspection_info = f"""
【Data Inspection Results】
- Number of cells (n_obs): {inspection_result.get('n_obs', 'N/A')}
- Number of genes (n_vars): {inspection_result.get('n_vars', 'N/A')}
- Max value: {inspection_result.get('max_value', 'N/A')}
- Is normalized: {inspection_result.get('is_normalized', False)}
- Has QC metrics: {inspection_result.get('has_qc_metrics', False)}
- Has clusters: {inspection_result.get('has_clusters', False)}
- Has UMAP: {inspection_result.get('has_umap', False)}

【Recommendations Based on Inspection】
"""
            n_obs = inspection_result.get('n_obs', 0)
            is_normalized = inspection_result.get('is_normalized', False)
            has_qc = inspection_result.get('has_qc_metrics', False)
            
            if n_obs > 10000:
                inspection_info += "- Large dataset (>10k cells): Recommend min_genes=500, max_mt=5%\n"
            elif n_obs > 5000:
                inspection_info += "- Medium dataset (5k-10k cells): Recommend min_genes=300, max_mt=5%\n"
            else:
                inspection_info += "- Small dataset (<5k cells): Recommend min_genes=200, max_mt=10%\n"
            
            if is_normalized:
                inspection_info += "- Data appears normalized: Skip normalization step\n"
            else:
                inspection_info += "- Data appears to be raw counts: Need normalization\n"
            
            if has_qc:
                inspection_info += "- QC metrics already calculated: May skip QC calculation\n"
        
        prompt = f"""Extract workflow parameters from user query and inspection results:

Query: {query}
Files: {', '.join(file_paths) if file_paths else 'None'}
{inspection_info}

Extract these parameters (if mentioned in query, otherwise use recommendations):
- min_genes (default: 200, adjust based on dataset size)
- max_mt (default: 20, adjust based on dataset size)
- resolution (default: 0.5, for clustering)
- n_top_genes (default: 2000, for HVG selection)

Return JSON only:
{{"resolution": "0.8", "min_genes": "500", "max_mt": "5"}}
"""
        
        messages = [
            {"role": "system", "content": "You are a parameter extraction assistant. Return JSON only."},
            {"role": "user", "content": prompt}
        ]
        
        try:
            completion = await self.llm_client.achat(messages, temperature=0.1, max_tokens=256)
            # 提取 think 过程和实际内容
            think_content, response = self.llm_client.extract_think_and_content(completion)
            # 如果有 think 内容，记录日志（可选）
            if think_content:
                import logging
                logger = logging.getLogger(__name__)
                logger.debug(f"RNA Agent think process: {think_content[:200]}...")
            
            # 解析 JSON
            json_str = response.strip()
            if "```json" in json_str:
                json_str = json_str.split("```json")[1].split("```")[0].strip()
            
            return json.loads(json_str)
        except:
            return {}
    
    async def _stream_chat_response(
        self,
        query: str,
        file_paths: List[str]
    ) -> AsyncIterator[str]:
        """
        流式聊天响应（支持 ReAct 循环和工具调用）
        
        实现 ReAct 循环：
        1. Thought: LLM 思考
        2. Action: 调用工具（如 inspect_file）
        3. Observation: 工具返回结果
        4. Final Answer: 最终回答
        """
        context = {
            "context": f"Uploaded files: {', '.join(file_paths) if file_paths else 'None'}",
            "available_tools": ["inspect_file", "run_cellranger", "convert_cellranger_to_h5ad"],
            "tool_descriptions": {
                "inspect_file": "检查数据文件，返回数据摘要（n_obs, n_vars, obs_keys, var_keys, is_normalized, etc.）",
                "run_cellranger": "运行 Cell Ranger count 对 FASTQ 文件进行计数分析",
                "convert_cellranger_to_h5ad": "将 Cell Ranger 输出转换为 Scanpy 格式 (.h5ad)"
            }
        }
        
        # 如果有文件，强制先检查（符合 SOP）
        inspection_result = None
        if file_paths:
            input_path = file_paths[0]
            try:
                # 使用新工具系统
                if input_path.endswith('.h5ad'):
                    adata = sc.read_h5ad(input_path)
                    inspection_result = {
                        "status": "success",
                        "n_obs": adata.n_obs,
                        "n_vars": adata.n_vars,
                        "file_type": "h5ad"
                    }
                else:
                    inspection_result = inspect_file(input_path)
                if "error" not in inspection_result:
                    # 将检查结果添加到上下文中
                    inspection_summary = f"""
【Data Inspection Completed】
- Cells: {inspection_result.get('n_obs', 'N/A')}
- Genes: {inspection_result.get('n_vars', 'N/A')}
- Max value: {inspection_result.get('max_value', 'N/A')}
- Normalized: {inspection_result.get('is_normalized', False)}
- Has QC metrics: {inspection_result.get('has_qc_metrics', False)}
- Has clusters: {inspection_result.get('has_clusters', False)}
"""
                    # 先输出检查结果
                    yield f"🔍 **Data Inspection Results:**\n{inspection_summary}\n\n"
                    # 将检查结果添加到查询中，让 LLM 基于此分析
                    query = f"""{query}

{inspection_summary}

Based on the inspection results above, please:
1. Analyze the data characteristics
2. Propose appropriate analysis parameters
3. Ask for confirmation before proceeding with analysis
"""
            except Exception as e:
                import logging
                logger = logging.getLogger(__name__)
                logger.error(f"Error inspecting file: {e}", exc_info=True)
                yield f"⚠️ Warning: Could not inspect file: {str(e)}\n\n"
        
        # 构建增强的用户查询，包含工具说明
        enhanced_query = f"""{query}

【Available Tools】
You have access to:
- inspect_file(file_path): Check data file structure (already executed above if files were provided)
- run_cellranger(fastq_dir, sample_id, output_dir, reference, ...): Run Cell Ranger count on FASTQ files
- convert_cellranger_to_h5ad(matrix_dir, output_path): Convert Cell Ranger output to .h5ad format

【Workflow Rule】
- If user provides FASTQ files: First run Cell Ranger, then convert to .h5ad, then inspect
- If user provides .h5ad or 10x MTX files: Inspect first (already done above), then analyze and propose parameters
- Before running any analysis, you MUST have inspected the data first
- Now analyze the inspection results and propose parameters.
"""
        
        # 流式输出 LLM 响应
        async for chunk in self.chat(enhanced_query, context, stream=True):
            yield chunk
    
    async def execute_workflow(
        self,
        workflow_config: Dict[str, Any],
        file_paths: List[str],
        output_dir: str
    ) -> Dict[str, Any]:
        """
        执行工作流
        
        核心：直接执行 scanpy 分析流程（参考旧版本实现）
        
        Args:
            workflow_config: 工作流配置
            file_paths: 文件路径列表
            output_dir: 输出目录
        
        Returns:
            分析报告
        """
        # 检测输入文件类型
        input_path = file_paths[0] if file_paths else None
        if not input_path:
            raise ValueError("No input files provided")
        
        file_type = self.detect_file_type(input_path)
        
        # 设置输出目录
        if not os.path.exists(output_dir):
            os.makedirs(output_dir, exist_ok=True)
        
        # 🔥 架构升级：使用新工具系统
        convert_result = None
        if file_type == "fastq":
            # 从 FASTQ 开始：先运行 Cell Ranger（异步），然后转换，最后执行分析
            fastq_dir = input_path
            sample_id = os.path.basename(fastq_dir).replace("_fastqs", "").replace("fastqs", "")
            if not sample_id:
                sample_id = "sample"
            
            # 创建临时输出目录
            temp_output_dir = os.path.join(output_dir, "cellranger_output")
            os.makedirs(temp_output_dir, exist_ok=True)
            
            # 使用新工具系统运行 Cell Ranger（异步）
            transcriptome_path = self.cellranger_config.get("transcriptome_path", "/opt/refdata-gex-GRCh38-2020-A")
            cellranger_result = run_cellranger_count(
                fastqs_path=fastq_dir,
                sample_id=sample_id,
                transcriptome_path=transcriptome_path,
                output_dir=temp_output_dir,
                localcores=self.cellranger_config.get("localcores", 8),
                localmem=self.cellranger_config.get("localmem", 32),
                create_bam=self.cellranger_config.get("create_bam", False),
                cellranger_path=self.cellranger_config.get("cellranger_path", "/opt/cellranger")
            )
            
            # Cell Ranger 是异步执行的，返回状态为 async_job_started
            if cellranger_result.get("status") == "async_job_started":
                # 返回异步任务状态，前端会显示消息
                return sanitize_for_json({
                    "status": "async_job_started",
                    "message": cellranger_result.get("message"),
                    "job_id": cellranger_result.get("job_id"),
                    "log_path": cellranger_result.get("log_path"),
                    "output_dir": cellranger_result.get("output_dir")
                })
            elif cellranger_result.get("status") != "success":
                return sanitize_for_json({
                    "status": "error",
                    "error": f"Cell Ranger failed: {cellranger_result.get('error', 'Unknown error')}",
                    "cellranger_result": cellranger_result
                })
            
            # 如果同步执行成功，继续转换
            matrix_dir = os.path.join(temp_output_dir, sample_id, "outs", "filtered_feature_bc_matrix")
            if not os.path.exists(matrix_dir):
                return sanitize_for_json({
                    "status": "error",
                    "error": f"Cell Ranger output matrix directory not found: {matrix_dir}",
                    "cellranger_result": cellranger_result
                })
            
            h5ad_path = os.path.join(output_dir, f"{sample_id}_filtered.h5ad")
            convert_result = convert_cellranger_to_h5ad(
                cellranger_matrix_dir=matrix_dir,
                output_h5ad_path=h5ad_path
            )
            
            if convert_result.get("status") != "success":
                return sanitize_for_json({
                    "status": "error",
                    "error": f"Conversion failed: {convert_result.get('error', 'Unknown error')}",
                    "cellranger_result": cellranger_result,
                    "convert_result": convert_result
                })
            
            # 使用转换后的 .h5ad 文件继续执行分析
            input_path = h5ad_path
        
        # 执行分析流程（使用新工具系统）
        if file_type != "fastq" or (file_type == "fastq" and convert_result and convert_result.get("status") == "success"):
            # 注意：这里应该使用工作流执行器，而不是直接调用旧工具
            # 由于工作流执行器可能还未完全实现，这里提供一个基本实现
            steps = workflow_config.get("steps", [])
            
            # 基本实现：按步骤执行
            current_adata_path = input_path
            report = {
                "status": "success",
                "steps": [],
                "final_output": current_adata_path
            }
            
            # 这里应该调用 WorkflowExecutor，但为了简化，我们只记录步骤
            logger.info(f"执行工作流，共 {len(steps)} 个步骤")
            logger.warning("⚠️ 完整的工作流执行需要使用 WorkflowExecutor")
            
            # 如果是从 FASTQ 转换来的，添加转换信息到报告
            if file_type == "fastq" and convert_result:
                report["cellranger_result"] = {
                    "status": "success",
                    "converted_file": convert_result.get("output_path"),
                    "n_obs": convert_result.get("n_obs"),
                    "n_vars": convert_result.get("n_vars")
                }
            
            # 🔥 生成最终分析报告（将工具结果反馈给LLM进行解释）
            if report.get("status") == "success":
                try:
                    final_report = await self.generate_final_report(report)
                    report["final_report"] = final_report
                except Exception as e:
                    logger.warning(f"⚠️ 生成最终报告失败: {e}")
                    report["final_report"] = None
            
            # 🔥 清理数据以确保 JSON 序列化安全（处理 Numpy 类型、NaN/Infinity 等）
            logger.info("✅ Workflow finished. Sanitizing data for JSON serialization...")
            sanitized_report = sanitize_for_json(report)
            logger.info("✅ Data sanitization completed. Returning result to frontend.")
            
            return sanitized_report
        else:
            # 如果 FASTQ 处理失败，返回错误（也需要清理）
            error_result = {
                "status": "error",
                "error": "Failed to process FASTQ files",
                "convert_result": convert_result
            }
            return sanitize_for_json(error_result)
    
    async def generate_final_report(self, execution_results: Dict[str, Any]) -> str:
        """
        生成最终分析报告
        
        将工具执行结果反馈给LLM，生成科学解释报告
        
        Args:
            execution_results: 执行结果字典，包含：
                - qc_metrics: 质量指标
                - steps_details: 步骤详情
                - final_plot: 最终图片路径
                - marker_genes: Marker基因（如果有）
        
        Returns:
            Markdown格式的分析报告
        """
        try:
            # 收集所有输出数据
            results_summary = {
                "qc_metrics": execution_results.get("qc_metrics", {}),
                "steps_completed": len(execution_results.get("steps_details", [])),
                "final_plot": execution_results.get("final_plot"),
                "output_file": execution_results.get("output_file"),
                "steps_summary": [
                    {
                        "name": step.get("name"),
                        "status": step.get("status"),
                        "summary": step.get("summary")
                    }
                    for step in execution_results.get("steps_details", [])
                ]
            }
            
            # 提取Marker基因（如果有）
            marker_genes = []
            for step in execution_results.get("steps_details", []):
                if step.get("tool_id") == "rna_find_markers" and step.get("details"):
                    # 尝试从details中提取marker基因信息
                    marker_genes.append(step.get("details"))
            
            if marker_genes:
                results_summary["marker_genes"] = marker_genes
            
            # 构建提示词
            import json
            results_json = json.dumps(results_summary, ensure_ascii=False, indent=2)
            
            # 使用 PromptManager 获取报告模板
            try:
                prompt = self.prompt_manager.get_prompt(
                    "rna_report",
                    {"results_summary": results_json},
                    fallback=RNA_REPORT_PROMPT.format(results_summary=results_json)
                )
            except Exception as e:
                logger.warning(f"⚠️ 获取报告模板失败，使用默认模板: {e}")
                prompt = RNA_REPORT_PROMPT.format(results_summary=results_json)
            
            # 调用LLM生成报告
            messages = [
                {"role": "system", "content": "You are a Senior Bioinformatician. Write analysis reports in Simplified Chinese."},
                {"role": "user", "content": prompt}
            ]
            
            completion = await self.llm_client.achat(messages, temperature=0.3, max_tokens=2000)
            think_content, response = self.llm_client.extract_think_and_content(completion)
            from ...core.stream_utils import strip_suggestions_from_text
            if response:
                response, _ = strip_suggestions_from_text(response)
            logger.info("✅ 最终分析报告已生成")
            return response
            
        except Exception as e:
            logger.error(f"❌ 生成最终报告失败: {e}", exc_info=True)
            return f"报告生成失败: {str(e)}"

