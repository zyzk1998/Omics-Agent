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
from ...core.prompt_manager import PromptManager, RNA_REPORT_PROMPT, DATA_DIAGNOSIS_PROMPT
from ...core.utils import sanitize_for_json
from ...core.dispatcher import TaskDispatcher
from ...core.test_data_manager import TestDataManager
from ...tools.cellranger_tool import CellRangerTool
from ...tools.scanpy_tool import ScanpyTool
import logging

logger = logging.getLogger(__name__)


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
        test_data_dir: Optional[str] = None
    ):
        """初始化转录组智能体"""
        super().__init__(llm_client, prompt_manager, "rna_expert")
        
        self.dispatcher = dispatcher
        self.cellranger_config = cellranger_config or {}
        self.scanpy_config = scanpy_config or {}
        self.cellranger_tool = CellRangerTool(self.cellranger_config)
        # 将 cellranger_tool 传递给 scanpy_tool，使其可以使用 Cell Ranger 功能
        self.scanpy_tool = ScanpyTool(self.scanpy_config, cellranger_tool=self.cellranger_tool)
        # 初始化测试数据管理器
        self.test_data_manager = TestDataManager(test_data_dir)
        
        # 标准工作流步骤（十步流程）
        self.workflow_steps = [
            {"name": "1. Quality Control", "tool_id": "local_qc", "desc": "过滤低质量细胞和基因"},
            {"name": "2. Normalization", "tool_id": "local_normalize", "desc": "数据标准化"},
            {"name": "3. Find Variable Genes", "tool_id": "local_hvg", "desc": "筛选高变基因"},
            {"name": "4. Scale Data", "tool_id": "local_scale", "desc": "数据缩放"},
            {"name": "5. PCA", "tool_id": "local_pca", "desc": "主成分分析"},
            {"name": "6. Compute Neighbors", "tool_id": "local_neighbors", "desc": "构建邻接图"},
            {"name": "7. Clustering", "tool_id": "local_cluster", "desc": "Leiden 聚类"},
            {"name": "8. UMAP Visualization", "tool_id": "local_umap", "desc": "UMAP 可视化"},
            {"name": "9. t-SNE Visualization", "tool_id": "local_tsne", "desc": "t-SNE 可视化"},
            {"name": "10. Find Markers", "tool_id": "local_markers", "desc": "寻找 Marker 基因"},
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
                # 使用 scanpy 工具检查文件
                if input_path.endswith('.h5ad'):
                    adata = self.scanpy_tool.load_data(input_path)
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
        
        # 智能数据检测：如果需要 Cell Ranger 但没有上传文件，提供测试数据选择
        needs_cellranger = self._needs_cellranger(query_lower)
        if needs_cellranger and not file_paths:
            # 检查是否有测试数据可用
            test_datasets = self.test_data_manager.scan_test_datasets()
            if test_datasets:
                # 返回测试数据选择请求
                return {
                    "type": "test_data_selection",
                    "message": "检测到您没有上传相关数据。请选择：",
                    "options": [
                        "1. 使用本地测试数据集",
                        "2. 上传您自己的数据"
                    ],
                    "datasets": test_datasets,
                    "datasets_json": self.test_data_manager.format_datasets_for_selection(test_datasets),
                    "datasets_display": self.test_data_manager.format_datasets_for_display(test_datasets)
                }
            else:
                # 没有测试数据，提示用户上传
                return {
                    "type": "chat",
                    "response": self._stream_string_response(
                        "检测到您没有上传相关数据，且没有可用的测试数据集。\n"
                        "请上传 FASTQ 文件或 .h5ad 文件以开始分析。"
                    )
                }
        
        # 处理测试数据选择（用户通过 JSON 选择）
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
        
        强制流程：
        1. 先检查文件（inspect_file）
        2. 基于检查结果提取参数
        3. 生成工作流配置
        """
        # 强制检查：如果有文件，先检查
        inspection_result = None
        diagnosis_report = None
        if file_paths:
            input_path = file_paths[0]
            try:
                inspection_result = self.scanpy_tool.inspect_file(input_path)
                if "error" in inspection_result:
                    # 检查失败，但仍然继续（可能是文件路径问题）
                    import logging
                    logger = logging.getLogger(__name__)
                    logger.warning(f"File inspection failed: {inspection_result.get('error')}")
                else:
                    # 🔥 生成数据诊断和参数推荐
                    diagnosis_report = await self._generate_diagnosis_and_recommendation(inspection_result)
            except Exception as e:
                import logging
                logger = logging.getLogger(__name__)
                logger.error(f"Error inspecting file: {e}", exc_info=True)
        
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
            if tool_id == "local_qc":
                step["params"] = {
                    "min_genes": extracted_params.get("min_genes", "200"),
                    "max_mt": extracted_params.get("max_mt", "20")
                }
            elif tool_id == "local_hvg":
                step["params"] = {
                    "n_top_genes": extracted_params.get("n_top_genes", "2000")
                }
            elif tool_id == "local_cluster":
                step["params"] = {
                    "resolution": extracted_params.get("resolution", "0.5")
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
        
        if diagnosis_report:
            result["diagnosis_report"] = diagnosis_report
        
        return result
    
    async def _generate_diagnosis_and_recommendation(
        self,
        inspection_result: Dict[str, Any]
    ) -> str:
        """
        生成数据诊断和参数推荐报告
        
        Args:
            inspection_result: 文件检查结果
        
        Returns:
            Markdown格式的诊断和推荐报告
        """
        try:
            import json
            # 格式化检查结果为JSON字符串
            inspection_json = json.dumps(inspection_result, ensure_ascii=False, indent=2)
            
            # 使用 PromptManager 获取诊断模板
            try:
                prompt = self.prompt_manager.get_prompt(
                    "data_diagnosis",
                    {"inspection_data": inspection_json},
                    fallback=DATA_DIAGNOSIS_PROMPT.format(inspection_data=inspection_json)
                )
            except Exception as e:
                import logging
                logger = logging.getLogger(__name__)
                logger.warning(f"⚠️ 获取诊断模板失败，使用默认模板: {e}")
                prompt = DATA_DIAGNOSIS_PROMPT.format(inspection_data=inspection_json)
            
            # 调用LLM生成诊断报告
            messages = [
                {"role": "system", "content": "You are a Senior Bioinformatician. Generate data diagnosis and parameter recommendations in Simplified Chinese."},
                {"role": "user", "content": prompt}
            ]
            
            completion = await self.llm_client.achat(messages, temperature=0.3, max_tokens=1500)
            think_content, response = self.llm_client.extract_think_and_content(completion)
            
            import logging
            logger = logging.getLogger(__name__)
            logger.info("✅ 数据诊断和参数推荐已生成")
            return response
            
        except Exception as e:
            import logging
            logger = logging.getLogger(__name__)
            logger.error(f"❌ 生成诊断报告失败: {e}", exc_info=True)
            return f"诊断报告生成失败: {str(e)}"
    
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
                inspection_result = self.scanpy_tool.inspect_file(input_path)
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
        
        # 更新 scanpy 工具的输出目录
        self.scanpy_config["output_dir"] = output_dir
        # 重新初始化 scanpy 工具以使用新的输出目录（保留 cellranger_tool）
        self.scanpy_tool = ScanpyTool(self.scanpy_config, cellranger_tool=self.cellranger_tool)
        
        # 直接执行 Scanpy 流程
        convert_result = None
        if file_type == "fastq":
            # 从 FASTQ 开始：先运行 Cell Ranger，然后转换，最后执行 Scanpy 分析
            # 提取参数
            fastq_dir = input_path
            sample_id = os.path.basename(fastq_dir).replace("_fastqs", "").replace("fastqs", "")
            if not sample_id:
                sample_id = "sample"
            
            # 创建临时输出目录
            temp_output_dir = os.path.join(output_dir, "cellranger_output")
            os.makedirs(temp_output_dir, exist_ok=True)
            
            # 运行 Cell Ranger
            cellranger_result = self.scanpy_tool.run_cellranger(
                fastq_dir=fastq_dir,
                sample_id=sample_id,
                output_dir=temp_output_dir,
                localcores=self.cellranger_config.get("localcores", 8),
                localmem=self.cellranger_config.get("localmem", 32),
                create_bam=self.cellranger_config.get("create_bam", False)
            )
            
            if cellranger_result.get("status") != "success":
                return sanitize_for_json({
                    "status": "error",
                    "error": f"Cell Ranger failed: {cellranger_result.get('error', 'Unknown error')}",
                    "cellranger_result": cellranger_result
                })
            
            # 转换 Cell Ranger 输出为 .h5ad
            matrix_dir = cellranger_result.get("matrix_dir")
            if not matrix_dir:
                return sanitize_for_json({
                    "status": "error",
                    "error": "Cell Ranger output matrix directory not found",
                    "cellranger_result": cellranger_result
                })
            
            h5ad_path = os.path.join(output_dir, f"{sample_id}_filtered.h5ad")
            convert_result = self.scanpy_tool.convert_cellranger_to_h5ad(
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
            
            # 使用转换后的 .h5ad 文件继续执行 Scanpy 分析
            input_path = h5ad_path
        
        # 执行 Scanpy 分析流程
        if file_type != "fastq" or (file_type == "fastq" and convert_result and convert_result.get("status") == "success"):
            # 直接运行 Scanpy 分析
            steps = workflow_config.get("steps", [])
            
            # 执行分析流程
            report = self.scanpy_tool.run_pipeline(
                data_input=input_path,
                steps_config=steps
            )
            
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
                if step.get("name") == "local_markers" and step.get("details"):
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
            
            logger.info("✅ 最终分析报告已生成")
            return response
            
        except Exception as e:
            logger.error(f"❌ 生成最终报告失败: {e}", exc_info=True)
            return f"报告生成失败: {str(e)}"

