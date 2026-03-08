"""代谢组学智能体（Metabolomics Agent）"""
from typing import Dict, Any, List, AsyncIterator, Optional
from ..base_agent import BaseAgent
from ...core.llm_client import LLMClient
from ...core.prompt_manager import PromptManager
from ...core.utils import sanitize_for_json
from ...core.tool_retriever import ToolRetriever
from ...core.planner import SOPPlanner
from ...core.tool_registry import registry
# 导入新工具函数
from ...tools.general.file_inspector import inspect_file
from ...tools.metabolomics.preprocessing import preprocess_metabolite_data
from ...tools.metabolomics.statistics import run_pca, run_differential_analysis
from ...tools.metabolomics.plotting import plot_volcano
import logging

logger = logging.getLogger(__name__)


# 🔥 架构重构：领域特定的系统指令（策略模式）
METABO_INSTRUCTION = """You are an expert Chemist/Metabolomics Analyst specializing in Metabolomics data analysis.

**CRITICAL CONSTRAINTS:**
- The data represents **Metabolite Abundance** (Chemical Compounds), NOT Gene Expression.
- Rows = Samples (Biological Samples), Columns = Metabolites (Chemical Compounds).
- This is Mass Spectrometry or LC-MS/GC-MS data, measuring chemical concentrations.

**STRICTLY FORBIDDEN TERMS:**
- Cell, Cells, Cellular
- Gene, Genes, Gene Expression, Transcript
- Mitochondria, Mitochondrial
- scRNA, Single-Cell RNA-seq, scRNA-seq
- Transcriptomics, Transcriptome
- RNA-seq, RNA sequencing

**REQUIRED TERMINOLOGY:**
- Metabolite, Metabolites, Metabolite Abundance
- Sample, Samples, Biological Sample
- Metabolomics, Metabolomic Analysis
- Mass Spectrometry, LC-MS, GC-MS
- Chemical Compound, Compound

**CONTEXT ISOLATION:**
This is NOT single-cell data. This is NOT transcriptomics data.
This is Metabolomics data representing metabolite abundance levels measured by mass spectrometry.

Generate data diagnosis and parameter recommendations in Simplified Chinese (简体中文).
Focus on metabolite-specific quality metrics (missing values, abundance range, normalization needs)."""


class MetabolomicsAgent(BaseAgent):
    """代谢组学智能体"""
    
    # 定义严格的步骤顺序（依赖链，与 MetabolomicsWorkflow DAG 一致）
    STEPS_ORDER = [
        "metabo_data_validation",  # 步骤0: 数据校验与稀疏性检查
        "inspect_data",          # 步骤1: 数据检查
        "preprocess_data",       # 步骤2: 数据预处理
        "pca_analysis",          # 步骤3: PCA 分析
        "metabo_model_comparison",  # 步骤3b: 多维模型对比
        "differential_analysis", # 步骤4: 差异分析
        "visualize_pca",        # 步骤5: PCA 可视化
        "visualize_volcano"     # 步骤6: 火山图可视化
    ]
    
    # 步骤映射（step_id -> 在 STEPS_ORDER 中的索引）
    STEP_INDEX_MAP = {step: idx for idx, step in enumerate(STEPS_ORDER)}
    
    def __init__(
        self,
        llm_client: LLMClient,
        prompt_manager: PromptManager,
        metabolomics_config: Dict[str, Any] = None,
        tool_retriever: Optional[ToolRetriever] = None
    ):
        super().__init__(llm_client, prompt_manager, "metabolomics_expert")
        self.metabolomics_config = metabolomics_config or {}
        # 🔥 架构升级：移除旧工具，使用新模块化工具系统
        
        # 🔥 架构升级：初始化 SOP 驱动的动态规划器
        self.sop_planner = None
        if tool_retriever:
            try:
                self.sop_planner = SOPPlanner(tool_retriever, llm_client)
                logger.info("✅ [MetabolomicsAgent] SOPPlanner 已初始化")
            except Exception as e:
                logger.warning(f"⚠️ [MetabolomicsAgent] SOPPlanner 初始化失败，将使用回退逻辑: {e}")
        else:
            logger.info("ℹ️ [MetabolomicsAgent] 未提供 ToolRetriever，将使用传统工作流生成逻辑")
        
        # 标准工作流步骤（代谢组学分析流程）- 保留作为回退
        self.workflow_steps = [
            {"name": "0. 数据校验", "step_id": "metabo_data_validation", "tool_id": "metabo_data_validation", "desc": "丰度矩阵与样本分组校验：shape、缺失率、零值比例"},
            {"name": "1. 数据检查", "step_id": "inspect_data", "tool_id": "inspect_data", "desc": "检查数据文件的基本信息（样本数、代谢物数、缺失值、分组信息等）"},
            {"name": "2. 数据预处理", "step_id": "preprocess_data", "tool_id": "preprocess_data", "desc": "数据预处理：处理缺失值、标准化、缩放"},
            {"name": "3. 主成分分析", "step_id": "pca_analysis", "tool_id": "pca_analysis", "desc": "执行主成分分析 (PCA)，降维并提取主要变异"},
            {"name": "3b. 多维模型对比", "step_id": "metabo_model_comparison", "tool_id": "metabo_model_comparison", "desc": "PCA + PLS-DA + VIP 1x3 对比图"},
            {"name": "4. 差异代谢物分析", "step_id": "differential_analysis", "tool_id": "differential_analysis", "desc": "执行差异代谢物分析（两组比较），识别显著差异的代谢物"},
            {"name": "5. PCA 可视化", "step_id": "visualize_pca", "tool_id": "visualize_pca", "desc": "生成 PCA 可视化图，展示样本在主成分空间的分布"},
            {"name": "6. 火山图可视化", "step_id": "visualize_volcano", "tool_id": "visualize_volcano", "desc": "生成火山图 (Volcano Plot)，展示差异代谢物的统计显著性"},
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
            - chat: 聊天响应（流式）
            - workflow_config: 工作流配置
        """
        # 🔥 架构重构：会话级文件注册表
        query_lower = query.lower().strip()
        file_paths = self.get_file_paths(uploaded_files or [])
        
        # 🔥 修复：当有新文件上传时，清除旧上下文，确保使用新文件
        if uploaded_files and len(uploaded_files) > 0:
            # 检查是否有新文件（通过比较文件名）
            current_active_file = self.context.get("active_file")
            new_file_names = []
            for file_info in uploaded_files:
                if isinstance(file_info, dict):
                    filename = file_info.get("name") or file_info.get("path") or file_info.get("file_id", "unknown")
                else:
                    filename = getattr(file_info, "name", None) or getattr(file_info, "path", None) or "unknown"
                new_file_names.append(filename)
            
            # 如果当前活动文件不在新文件列表中，说明有新文件上传
            has_new_file = current_active_file not in new_file_names if current_active_file else True
            
            if has_new_file:
                logger.info(f"🔄 [FileRegistry] 检测到新文件上传，清除旧上下文")
                logger.info(f"   旧活动文件: {current_active_file}")
                logger.info(f"   新文件列表: {new_file_names}")
                # 清除旧上下文
                self._refresh_context_for_new_files(uploaded_files)
            
            # 注册所有上传的文件
            for file_info in uploaded_files:
                if isinstance(file_info, dict):
                    filename = file_info.get("name") or file_info.get("path") or file_info.get("file_id", "unknown")
                    file_path = file_info.get("path") or file_info.get("file_id", filename)
                else:
                    filename = getattr(file_info, "name", None) or getattr(file_info, "path", None) or "unknown"
                    file_path = getattr(file_info, "path", None) or filename
                
                # 找到对应的绝对路径
                if file_paths:
                    # 尝试匹配路径
                    absolute_path = None
                    for abs_path in file_paths:
                        if file_path in abs_path or abs_path.endswith(file_path.split('/')[-1]):
                            absolute_path = abs_path
                            break
                    if not absolute_path:
                        absolute_path = file_paths[0]  # 使用第一个路径作为回退
                else:
                    absolute_path = file_path
                
                # 注册文件
                self.register_file(filename, absolute_path, file_metadata=None)
            
            # 🔥 关键修复：优先使用最新上传的文件（列表最后一个）作为活动文件
            if new_file_names:
                latest_file = new_file_names[-1]
                self.set_active_file(latest_file)
                logger.info(f"✅ [FileRegistry] 设置最新文件为活动文件: {latest_file}")
                # 确保 file_paths 使用最新文件
                if file_paths:
                    file_paths = [file_paths[-1]]  # 只使用最后一个文件
                    logger.info(f"📂 [FileRegistry] 使用最新文件路径: {file_paths[0]}")
        
        # Scenario B: 没有新文件 - 使用当前活动文件
        if not file_paths:
            active_file_info = self.get_active_file_info()
            if active_file_info:
                file_paths = [active_file_info["path"]]
                logger.info(f"📂 [FileRegistry] Using active file: {active_file_info['filename']}")
            else:
                logger.warning("⚠️ [FileRegistry] No files available (no uploads and no active file)")
        
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
            # 检查最新上传的文件
            input_path = file_paths[-1] if file_paths else None
            if not input_path:
                return {
                    "type": "chat",
                    "response": self._stream_string_response("没有检测到上传的文件。请先上传文件后再询问。")
                }
            try:
                # 使用新工具系统
                inspection_result = inspect_file(input_path)
                if "error" in inspection_result:
                    return {
                        "type": "chat",
                        "response": self._stream_string_response(f"文件检查失败: {inspection_result.get('error')}")
                    }
                
                # 使用 LLM 生成文件解释
                explanation = await self._explain_file_with_llm(query, inspection_result, input_path)
                return {
                    "type": "chat",
                    "response": self._stream_string_response(explanation)
                }
            except Exception as e:
                logger.error(f"❌ 文件解释失败: {e}", exc_info=True)
                return {
                    "type": "chat",
                    "response": self._stream_string_response(f"文件检查时出错: {str(e)}")
                }
        
        # 🔒 回退逻辑：如果意图检测失败或意图不明确，使用原始逻辑
        if intent is None or intent == "chat":
            # 使用原始的工作流检测逻辑作为回退
            is_workflow_request = self._is_workflow_request(query_lower, file_paths)
            if is_workflow_request:
                # 工作流请求：生成工作流配置
                return await self._generate_workflow_config(query, file_paths)
            else:
                # 普通聊天：流式响应
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
        inspection_result: Dict[str, Any],
        file_path: str
    ) -> str:
        """
        使用 LLM 生成文件解释
        
        Args:
            query: 用户查询
            inspection_result: 文件检查结果
            file_path: 文件路径
        
        Returns:
            自然语言的文件解释
        """
        import json
        
        # 格式化检查结果
        inspection_summary = json.dumps(inspection_result, ensure_ascii=False, indent=2)
        
        prompt = f"""用户询问关于文件的问题。

User Query: {query}
File Path: {file_path}

文件检查结果：
{inspection_summary}

请用自然语言解释这个文件的内容、结构和特点。回答应该：
1. 简洁明了，易于理解
2. 包含关键信息（样本数、变量数、缺失值等）
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
        """判断是否是工作流请求（保留用于向后兼容）"""
        workflow_keywords = [
            "规划", "流程", "workflow", "pipeline", "分析", "run",
            "执行", "plan", "做一下", "跑一下", "分析一下", "全流程",
            "代谢组", "代谢组分析", "metabolomics"
        ]
        
        # 🔧 修复：如果查询包含工作流关键词，返回 True
        if query and any(kw in query for kw in workflow_keywords):
            return True
        
        # 🔧 修复：如果只有文件没有文本（或文本很短），且不是自我介绍等常见查询，返回 True
        if file_paths and (not query or len(query.strip()) < 5):
            # 排除一些常见的非工作流查询
            non_workflow_queries = ["你好", "hello", "hi", "介绍", "自我介绍", "你是谁", "who are you"]
            if not query or query.strip().lower() not in [q.lower() for q in non_workflow_queries]:
                return True
        
        return False
    
    async def _generate_workflow_config(
        self,
        query: str,
        file_paths: List[str]
    ) -> Dict[str, Any]:
        """
        生成工作流配置
        
        🔥 架构升级：优先使用 SOP 驱动的动态规划器，失败则回退到传统逻辑
        
        流程：
        1. 如果 SOPPlanner 可用，使用动态规划器生成工作流
        2. 否则，使用传统硬编码逻辑（回退）
        """
        logger.info("=" * 80)
        logger.info("🚀 [CHECKPOINT] _generate_workflow_config START")
        logger.info(f"   Query: {query}")
        logger.info(f"   File paths: {file_paths}")
        logger.info(f"   SOPPlanner available: {self.sop_planner is not None}")
        
        # 🔥 修复：确保使用最新文件路径，清除缓存的旧文件元数据
        if file_paths:
            # 使用最新文件（列表最后一个）
            current_file = file_paths[-1]
            logger.info(f"   Using latest file: {current_file}")
            
            # 清除旧的元数据缓存，强制重新检查文件
            self.context.pop("file_metadata", None)
            self.context.pop("diagnosis_report", None)
            self.context.pop("diagnosis_stats", None)
            logger.info("   ✅ Cleared cached file metadata, will re-inspect file")
        else:
            # 如果没有文件路径，尝试使用活动文件
            active_file_info = self.get_active_file_info()
            if active_file_info:
                current_file = active_file_info["path"]
                file_paths = [current_file]
                logger.info(f"   Using active file: {current_file}")
            else:
                logger.warning("   ⚠️ No file paths provided and no active file")
        
        logger.info("=" * 80)
        
        # 🔥 Phase 3: 优先使用 SOP 驱动的动态规划器
        if self.sop_planner:
            try:
                logger.info("🧠 [SOPPlanner] 尝试使用动态规划器生成工作流...")
                
                # Step 1: 检查文件获取元数据
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
                
                # Step 2: 使用 SOPPlanner 生成计划
                if file_metadata:
                    plan_result = await self.sop_planner.generate_plan(
                        user_query=query,
                        file_metadata=file_metadata,
                        category_filter="Metabolomics"
                    )
                    
                    # 检查是否成功
                    if plan_result.get("type") != "error":
                        logger.info("✅ [SOPPlanner] 动态规划成功")
                        
                        # 生成诊断报告（可选）
                        diagnosis_report = None
                        try:
                            diagnosis_report = await self._perform_data_diagnosis(
                                file_metadata=file_metadata,
                                system_instruction=METABO_INSTRUCTION
                            )
                        except Exception as e:
                            logger.warning(f"⚠️ 诊断报告生成失败: {e}")
                        
                        # 添加诊断报告到结果
                        if diagnosis_report:
                            plan_result["diagnosis_report"] = diagnosis_report
                        
                        return plan_result
                    else:
                        logger.warning(f"⚠️ [SOPPlanner] 规划失败: {plan_result.get('error')}")
                else:
                    logger.warning("⚠️ [SOPPlanner] 文件元数据不可用，回退到传统逻辑")
            
            except Exception as e:
                logger.error(f"❌ [SOPPlanner] 动态规划异常: {e}", exc_info=True)
                logger.info("🔄 回退到传统工作流生成逻辑...")
        
        # 🔄 回退：使用传统硬编码逻辑
        logger.info("📋 [Fallback] 使用传统工作流生成逻辑...")
        
        # 🔥 Step 1: 立即检查文件（修复 N/A 问题）
        file_metadata = None
        stats = {"n_samples": "N/A", "n_features": "N/A"}
        diagnosis_report = None
        recommendation = None
        
        if not file_paths:
            logger.warning("⚠️ 没有提供文件路径")
            return {
                "type": "workflow_config",
                "workflow_data": {
                    "workflow_name": "Metabolomics Analysis Pipeline",
                    "steps": []
                },
                "file_paths": [],
                "diagnosis_report": "⚠️ 未提供数据文件，无法生成诊断报告。"
            }
        
        current_file = file_paths[0]
        logger.info(f"🔍 [CHECKPOINT] Inspecting file IMMEDIATELY: {current_file}")
        
        try:
            # 🔥 使用 FileInspector 立即检查文件
            from ...core.file_inspector import FileInspector
            import os
            upload_dir = os.getenv("UPLOAD_DIR", "/app/uploads")
            inspector = FileInspector(upload_dir)
            
            file_metadata = inspector.inspect_file(current_file)
            
            if file_metadata.get("status") != "success" or not file_metadata.get("success", True):
                error_msg = file_metadata.get("error", "未知错误")
                logger.warning(f"⚠️ File inspection failed: {error_msg}")
                
                # 🔥 修复：使用详细的错误信息，而不是硬编码消息
                # 这样用户可以看到系统在哪里查找文件
                diagnosis_report = f"⚠️ **文件读取失败**\n\n{error_msg}"
                
                # Fallback: 使用默认值
                stats = {"n_samples": "N/A", "n_features": "N/A"}
            else:
                logger.info(f"✅ [CHECKPOINT] File inspection successful")
                
                # 🔥 CRITICAL FIX: 映射 shape.rows/cols 到 n_samples/n_features
                shape = file_metadata.get("shape", {})
                stats = {
                    "n_samples": shape.get("rows", file_metadata.get("n_samples", "N/A")),
                    "n_features": shape.get("cols", file_metadata.get("n_features", "N/A"))
                }
                
                # 如果 shape 中没有，尝试从 file_metadata 直接获取
                if stats["n_samples"] == "N/A" or stats["n_samples"] == 0:
                    stats["n_samples"] = file_metadata.get("n_samples", "N/A")
                if stats["n_features"] == "N/A" or stats["n_features"] == 0:
                    stats["n_features"] = file_metadata.get("n_features", "N/A")
                
                logger.info(f"📊 [CHECKPOINT] Stats mapped: n_samples={stats['n_samples']}, n_features={stats['n_features']}")
                
                # 🔥 Step 2: 立即生成诊断报告（修复 UI 缺失报告问题）
                try:
                    logger.info(f"🔍 [CHECKPOINT] Generating diagnosis report IMMEDIATELY...")
                    
                    # 尝试加载数据预览（用于更准确的诊断）
                    dataframe = None
                    try:
                        import pandas as pd
                        head_data = file_metadata.get("head", {})
                        if head_data and isinstance(head_data, dict) and "json" in head_data:
                            dataframe = pd.DataFrame(head_data["json"])
                    except Exception as e:
                        logger.debug(f"无法构建数据预览: {e}")
                    
                    # 调用统一的诊断方法（在 planning 阶段）
                    # 🔥 架构重构：传递领域特定的系统指令
                    diagnosis_report = await self._perform_data_diagnosis(
                        file_metadata=file_metadata,
                        omics_type="Metabolomics",
                        dataframe=dataframe,
                        system_instruction=METABO_INSTRUCTION
                    )
                    
                    if diagnosis_report:
                        logger.info(f"✅ [CHECKPOINT] Diagnosis report generated, length: {len(diagnosis_report)}")
                    else:
                        logger.warning(f"⚠️ [CHECKPOINT] Diagnosis report is None")
                        diagnosis_report = "⚠️ 诊断报告生成失败，但可以继续进行分析。"
                    
                except Exception as diag_err:
                    logger.error(f"❌ [CHECKPOINT] Diagnosis generation failed: {diag_err}", exc_info=True)
                    diagnosis_report = "⚠️ 诊断报告生成失败，但可以继续进行分析。"
                
                # 从诊断结果中提取推荐参数（如果可用）
                if diagnosis_report and hasattr(self, 'context') and self.context.get("diagnosis_stats"):
                    stats_context = self.context.get("diagnosis_stats", {})
                    recommendations = stats_context.get("recommendations", {})
                    if recommendations:
                        # 转换为 MetabolomicsAgent 期望的格式
                        recommendation = {
                            "params": {
                                "normalization": {
                                    "value": recommendations.get("normalization", {}).get("recommended", "log2")
                                },
                                "missing_threshold": {
                                    "value": "0.5"  # 默认值
                                },
                                "scale": {
                                    "value": True
                                },
                                "n_components": {
                                    "value": "10"
                                }
                            }
                        }
                    else:
                        recommendation = None
                else:
                    recommendation = None
                    
        except Exception as e:
            logger.error(f"❌ [CHECKPOINT] Error inspecting file: {e}", exc_info=True)
            stats = {"n_samples": "N/A", "n_features": "N/A"}
            diagnosis_report = "⚠️ 文件检查失败，请检查文件路径和格式。"
        
        # 使用 file_metadata 作为 inspection_result
        inspection_result = file_metadata
        
        # 使用 LLM 提取目标结束步骤（例如："做到PCA" -> "pca_analysis"）
        target_end_step = None
        try:
            logger.info(f"🔍 [CHECKPOINT] Extracting target end step from query...")
            target_end_step = await self._extract_target_end_step(query, inspection_result)
            logger.info(f"✅ [CHECKPOINT] Target end step extracted: {target_end_step}")
        except Exception as e:
            logger.error(f"❌ [CHECKPOINT] Error extracting target end step: {e}", exc_info=True)
            target_end_step = None  # 使用默认值（所有步骤）
        
        # 🔥 Task 1: 使用推荐值或 LLM 提取参数（优先使用推荐值）
        extracted_params = {}
        if recommendation and "params" in recommendation:
            # 优先使用推荐值
            rec_params = recommendation["params"]
            extracted_params = {
                "normalization": rec_params.get("normalization", {}).get("value", "log2"),
                "missing_threshold": rec_params.get("missing_threshold", {}).get("value", "0.5"),
                "scale": str(rec_params.get("scale", {}).get("value", True)).lower(),
                "n_components": rec_params.get("n_components", {}).get("value", "10")
            }
            logger.info(f"✅ [CHECKPOINT] Using recommended parameters: {extracted_params}")
        else:
            # 如果没有推荐，使用 LLM 提取
            try:
                logger.info(f"🔍 [CHECKPOINT] Extracting workflow parameters with LLM...")
                extracted_params = await self._extract_workflow_params(query, file_paths, inspection_result, None)
                logger.info(f"✅ [CHECKPOINT] Workflow parameters extracted: {list(extracted_params.keys())}")
            except Exception as e:
                logger.error(f"❌ [CHECKPOINT] Error extracting workflow params: {e}", exc_info=True)
                extracted_params = {}  # 使用默认值
        
        # 🔥 修复 2: 启发式检测分组列（如果未指定）
        if not extracted_params.get("group_column"):
            detected_group_col = self._detect_group_column_heuristic(file_metadata)
            if detected_group_col:
                extracted_params["group_column"] = detected_group_col
                logger.info(f"✅ [Heuristic] 自动检测到分组列: {detected_group_col}")
            else:
                # 回退到默认值（但记录警告）
                extracted_params["group_column"] = "Group"  # 默认值
                logger.warning(f"⚠️ 未检测到分组列，使用默认值: Group")
        
        # 定义所有可用步骤（包含友好的中文名称）
        all_steps = [
            {
                "step_id": "inspect_data",
                "tool_id": "inspect_data",
                "name": "数据检查",  # 🔧 修复：添加 name 字段
                "step_name": "数据检查",  # 🔧 修复：添加 step_name 字段（兼容前端）
                "desc": f"检查数据文件的基本信息（样本数: {stats.get('n_samples', 'N/A')}, 代谢物数: {stats.get('n_features', 'N/A')}）",
                "params": {"file_path": current_file if file_paths else ""}
            },
            {
                "step_id": "preprocess_data",
                "tool_id": "preprocess_data",
                "name": "数据预处理",  # 🔧 修复：添加 name 字段
                "step_name": "数据预处理",  # 🔧 修复：添加 step_name 字段（兼容前端）
                "desc": "数据预处理：处理缺失值、标准化、缩放",
                "params": {
                    "file_path": current_file if file_paths else "",
                    "missing_threshold": extracted_params.get("missing_threshold", "0.5"),
                    "normalization": extracted_params.get("normalization", "log2"),
                    "scale": extracted_params.get("scale", "true")
                }
            },
            {
                "step_id": "pca_analysis",
                "tool_id": "pca_analysis",
                "name": "主成分分析",  # 🔧 修复：添加 name 字段
                "step_name": "主成分分析",  # 🔧 修复：添加 step_name 字段（兼容前端）
                "desc": "执行主成分分析 (PCA)，降维并提取主要变异",
                "params": {
                    "n_components": extracted_params.get("n_components", "10")
                }
            },
            {
                "step_id": "differential_analysis",
                "tool_id": "differential_analysis",
                "name": "差异代谢物分析",  # 🔧 修复：添加 name 字段
                "step_name": "差异代谢物分析",  # 🔧 修复：添加 step_name 字段（兼容前端）
                "desc": "执行差异代谢物分析（两组比较），识别显著差异的代谢物",
                "params": {
                    "group_column": extracted_params.get("group_column", "Group"),  # 🔥 修复：使用启发式检测的值
                    "method": extracted_params.get("method", "t-test"),
                    "p_value_threshold": extracted_params.get("p_value_threshold", "0.05"),
                    "fold_change_threshold": extracted_params.get("fold_change_threshold", "1.5"),
                    "group1": extracted_params.get("group1"),
                    "group2": extracted_params.get("group2")
                }
            },
            {
                "step_id": "visualize_pca",
                "tool_id": "visualize_pca",
                "name": "PCA 可视化",  # 🔧 修复：添加 name 字段
                "step_name": "PCA 可视化",  # 🔧 修复：添加 step_name 字段（兼容前端）
                "desc": "生成 PCA 可视化图，展示样本在主成分空间的分布",
                "params": {
                    "group_column": extracted_params.get("group_column", "Group"),  # 🔥 修复：使用启发式检测的值
                    "pc1": "1",
                    "pc2": "2"
                }
            },
            {
                "step_id": "visualize_volcano",
                "tool_id": "visualize_volcano",
                "name": "火山图可视化",  # 🔧 修复：添加 name 字段
                "step_name": "火山图可视化",  # 🔧 修复：添加 step_name 字段（兼容前端）
                "desc": "生成火山图 (Volcano Plot)，展示差异代谢物的统计显著性",
                "params": {
                    "p_value_threshold": extracted_params.get("p_value_threshold", "0.05"),
                    "fold_change_threshold": extracted_params.get("fold_change_threshold", "1.5")
                }
            }
        ]
        
        # 根据目标结束步骤自动包含所有前置步骤
        if target_end_step and target_end_step in self.STEP_INDEX_MAP:
            # 找到目标步骤的索引
            target_index = self.STEP_INDEX_MAP[target_end_step]
            # 包含从开始到目标步骤的所有步骤（包括目标步骤）
            required_step_ids = self.STEPS_ORDER[:target_index + 1]
            logger.info(f"🎯 目标步骤: {target_end_step}, 将执行: {required_step_ids}")
            
            # 构建步骤映射并筛选
            step_map = {step["step_id"]: step for step in all_steps}
            selected_steps = [step_map[s] for s in required_step_ids if s in step_map]
        else:
            # 如果没有指定或无效，使用所有步骤
            selected_steps = all_steps
        
        # 构建工作流配置
        workflow_config = {
            "workflow_name": "Metabolomics Analysis Pipeline",
            "steps": selected_steps
        }
        
        # 🔥 Task 1: 构建返回结果，包含推荐信息和诊断报告
        result = {
            "type": "workflow_config",
            "workflow_data": workflow_config,
            "file_paths": file_paths
        }
        
        # 添加诊断报告（如果生成成功）
        # 🔥 修复：检查 diagnosis_report 是否为有效字符串（非 None 且非空）
        if diagnosis_report and isinstance(diagnosis_report, str) and diagnosis_report.strip():
            result["diagnosis_report"] = diagnosis_report
            logger.info(f"📝 [DEBUG] Adding diagnosis_report to result, length: {len(diagnosis_report)}")
        else:
            logger.warning(f"⚠️ [DEBUG] diagnosis_report is invalid (None/empty), NOT adding to result. Type: {type(diagnosis_report)}, Value: {diagnosis_report}")
        
        # 添加推荐信息（如果生成成功）
        if recommendation:
            result["recommendation"] = recommendation
            # 自动填充推荐值到步骤参数
            self._apply_recommendations_to_steps(workflow_config["steps"], recommendation)
        
        logger.info("=" * 80)
        logger.info("✅ [CHECKPOINT] _generate_workflow_config SUCCESS")
        logger.info(f"   Workflow name: {workflow_config.get('workflow_name')}")
        logger.info(f"   Steps count: {len(workflow_config.get('steps', []))}")
        logger.info(f"   Has recommendation: {recommendation is not None}")
        logger.info(f"   Has diagnosis_report: {diagnosis_report is not None}")
        
        # 🔥 DEBUG: 打印最终返回结构
        logger.info(f"📤 [DEBUG] MetabolomicsAgent returning result with keys: {list(result.keys())}")
        logger.info(f"📤 [DEBUG] MetabolomicsAgent has diagnosis_report: {'diagnosis_report' in result}")
        if 'diagnosis_report' in result:
            logger.info(f"📤 [DEBUG] MetabolomicsAgent diagnosis_report length: {len(result['diagnosis_report'])}")
        logger.info("=" * 80)
        
        return result
    
    def _detect_group_column_heuristic(self, file_metadata: Dict[str, Any]) -> Optional[str]:
        """
        启发式检测分组列
        
        🔥 修复 2: 工具健壮性 - 自动检测分组列，避免硬编码 "Group"
        
        Args:
            file_metadata: FileInspector 返回的文件元数据
        
        Returns:
            检测到的分组列名，如果未找到返回 None
        """
        # 优先级关键词列表
        priority_keywords = ['Diet', 'diet', 'Group', 'group', 'Condition', 'condition', 
                            'Treatment', 'treatment', 'Class', 'class', 'Category', 'category',
                            'Type', 'type', 'Label', 'label', 'Status', 'status']
        
        # 从 file_metadata 获取列信息
        columns = file_metadata.get("columns", [])
        if not columns:
            logger.warning("⚠️ 无法获取列信息，无法检测分组列")
            return None
        
        # 方法1: 检查列名是否包含优先级关键词
        for col in columns:
            if any(keyword in col for keyword in priority_keywords):
                logger.info(f"✅ [Heuristic] 检测到分组列（关键词匹配）: {col}")
                return col
        
        # 方法2: 检查 potential_groups（FileInspector 可能已经检测到）
        potential_groups = file_metadata.get("potential_groups", {})
        if isinstance(potential_groups, dict) and len(potential_groups) > 0:
            # 返回第一个潜在分组列
            first_group_col = list(potential_groups.keys())[0]
            logger.info(f"✅ [Heuristic] 检测到分组列（potential_groups）: {first_group_col}")
            return first_group_col
        
        # 方法3: 查找第一个非数值的分类列（唯一值数量 < 样本数的50%）
        try:
            import pandas as pd
            head_data = file_metadata.get("head", {})
            if head_data and isinstance(head_data, dict) and "json" in head_data:
                df_preview = pd.DataFrame(head_data["json"])
                n_samples = len(df_preview)
                
                for col in columns:
                    if col in df_preview.columns:
                        # 检查是否为数值类型
                        if not pd.api.types.is_numeric_dtype(df_preview[col]):
                            # 检查唯一值数量
                            unique_count = df_preview[col].nunique()
                            if 2 <= unique_count <= max(2, n_samples * 0.5):
                                logger.info(f"✅ [Heuristic] 检测到分组列（分类列）: {col} (唯一值: {unique_count})")
                                return col
        except Exception as e:
            logger.warning(f"⚠️ 启发式检测失败: {e}")
        
        logger.warning("⚠️ 未检测到分组列")
        return None
    
    # 🔥 已移除：_generate_diagnosis_and_recommendation 方法
    # 现在使用 BaseAgent._perform_data_diagnosis() 统一方法
    
    async def _extract_target_end_step(
        self,
        query: str,
        inspection_result: Dict[str, Any] = None
    ) -> Optional[str]:
        """
        使用 LLM 从用户查询中提取目标结束步骤
        
        支持：
        - "做到PCA" / "up to PCA" -> "pca_analysis" (会自动包含 inspect_data, preprocess_data)
        - "只做预处理" -> "preprocess_data" (会自动包含 inspect_data)
        - "做到差异分析" -> "differential_analysis" (会自动包含所有前置步骤)
        - 默认返回 None（执行所有步骤）
        
        Returns:
            目标结束步骤的 step_id，或 None（执行所有步骤）
        """
        # 定义步骤关键词映射（用于匹配"做到XXX"）
        step_keywords = {
            "inspect_data": ["检查", "inspect", "检查数据", "数据检查", "步骤1", "step 1", "第一步", "到检查"],
            "preprocess_data": ["预处理", "preprocess", "数据预处理", "步骤2", "step 2", "第二步", "到预处理"],
            "pca_analysis": ["pca", "主成分", "主成分分析", "步骤3", "step 3", "第三步", "做到pca", "到pca", "up to pca"],
            "differential_analysis": ["差异", "differential", "差异分析", "步骤4", "step 4", "第四步", "到差异分析"],
            "visualize_pca": ["pca图", "pca可视化", "pca plot", "步骤5", "step 5", "第五步", "到pca图"],
            "visualize_volcano": ["火山图", "volcano", "volcano plot", "步骤6", "step 6", "第六步", "到火山图"]
        }
        
        # 先进行关键词匹配（优先匹配"做到"、"up to"等表达）
        query_lower = query.lower()
        
        # 检查"做到"、"up to"等表达
        for step_id, keywords in step_keywords.items():
            for kw in keywords:
                if kw in query_lower:
                    # 特别检查"做到"、"up to"等表达
                    if any(phrase in query_lower for phrase in ["做到", "up to", "until", "到"]):
                        logger.info(f"🎯 检测到目标结束步骤: {step_id} (关键词: {kw})")
                        return step_id
        
        # 如果没有匹配到"做到"表达，使用 LLM 提取
        prompt = f"""
Extract the TARGET END STEP from the user query. This is the LAST step the user wants to run.

User Query: {query}

Available steps (in order):
1. inspect_data - Check data file
2. preprocess_data - Preprocess data
3. pca_analysis - PCA analysis
4. differential_analysis - Differential analysis
5. visualize_pca - PCA visualization
6. visualize_volcano - Volcano plot visualization

Important: If the user says "do analysis up to PCA" or "做到PCA", they want to run steps 1, 2, and 3 (all prerequisites + PCA).
If they say "only preprocessing", they want steps 1 and 2.

Examples:
- "做到PCA" / "up to PCA" -> "pca_analysis"
- "只做预处理" -> "preprocess_data"
- "做到差异分析" -> "differential_analysis"
- "做完整分析" / "全部" -> null (run all steps)

Return JSON only (single step_id string or null):
"""
        
        messages = [
            {"role": "system", "content": "You are a step extraction assistant. Return JSON only (single string or null)."},
            {"role": "user", "content": prompt}
        ]
        
        try:
            logger.info(f"🔍 [CHECKPOINT] Calling LLM to extract target end step...")
            completion = await self.llm_client.achat(messages, temperature=0.1, max_tokens=64)
            think_content, response = self.llm_client.extract_think_and_content(completion)
            logger.info(f"✅ [CHECKPOINT] LLM response received: {response[:100]}...")
            
            import json
            json_str = response.strip()
            if "```json" in json_str:
                json_str = json_str.split("```json")[1].split("```")[0].strip()
            elif "```" in json_str:
                json_str = json_str.split("```")[1].split("```")[0].strip()
            
            result = json.loads(json_str)
            if result is None or result == "":
                logger.info(f"✅ [CHECKPOINT] LLM returned null, will run all steps")
                return None
            # 验证结果是否在有效步骤列表中
            if isinstance(result, str) and result in self.STEPS_ORDER:
                logger.info(f"🎯 LLM 提取的目标结束步骤: {result}")
                return result
            logger.warning(f"⚠️ LLM returned invalid step: {result}, will run all steps")
            return None
        except json.JSONDecodeError as e:
            logger.error(f"❌ [CHECKPOINT] JSON decode error: {e}, response: {response[:200]}")
            return None  # 默认执行所有步骤
        except Exception as e:
            logger.error(f"❌ [CHECKPOINT] Error extracting target end step: {e}", exc_info=True)
            return None  # 默认执行所有步骤
    
    async def _extract_workflow_params(
        self,
        query: str,
        file_paths: List[str],
        inspection_result: Dict[str, Any] = None,
        diagnosis_report: Optional[str] = None
    ) -> Dict[str, Any]:
        """
        使用 LLM 提取工作流参数
        
        基于检查结果和诊断报告智能推荐参数
        """
        # 构建包含检查结果和诊断报告的提示
        inspection_info = ""
        if diagnosis_report:
            inspection_info = f"""
【Data Diagnosis & Recommendations】
{diagnosis_report}

"""
        elif inspection_result and "error" not in inspection_result:
            inspection_info = f"""
【Data Diagnosis & Recommendations】
{diagnosis_report}

"""
        elif inspection_result and "error" not in inspection_result:
            inspection_info = f"""
【Data Inspection Results】
- Number of samples: {inspection_result.get('n_samples', 'N/A')}
- Number of metabolites: {inspection_result.get('n_metabolites', 'N/A')}
- Missing values: {inspection_result.get('missing_values', {}).get('percentage', 'N/A')}%
- Group column: {inspection_result.get('group_info', {}).get('column', 'N/A')}
- Groups: {inspection_result.get('group_info', {}).get('groups', {})}
"""
        
        prompt = f"""
Extract workflow parameters for metabolomics analysis from the user query.

User Query: {query}
File Paths: {file_paths}

{inspection_info}

Based on the inspection results and user query, extract the following parameters:
- missing_threshold: Threshold for removing metabolites with high missing values (default: 0.5)
- normalization: Normalization method - "log2", "zscore", or "none" (default: "log2")
- scale: Whether to apply StandardScaler (default: true)
- n_components: Number of PCA components (default: 10)
- group_column: Column name for group comparison (default: "Muscle loss" or first metadata column)
- method: Statistical method for differential analysis - "t-test" or "mann-whitney" (default: "t-test")
- p_value_threshold: P-value threshold for significance (default: 0.05)
- fold_change_threshold: Fold change threshold (default: 1.5)

Return JSON only:
{{"missing_threshold": "0.5", "normalization": "log2", "scale": "true", "n_components": "10", "group_column": "Muscle loss", "method": "t-test", "p_value_threshold": "0.05", "fold_change_threshold": "1.5"}}
"""
        
        messages = [
            {"role": "system", "content": "You are a parameter extraction assistant. Return JSON only."},
            {"role": "user", "content": prompt}
        ]
        
        try:
            logger.info(f"🔍 [CHECKPOINT] Calling LLM to extract workflow parameters...")
            completion = await self.llm_client.achat(messages, temperature=0.1, max_tokens=256)
            think_content, response = self.llm_client.extract_think_and_content(completion)
            logger.info(f"✅ [CHECKPOINT] LLM response received: {response[:200]}...")
            
            # 解析 JSON
            import json
            json_str = response.strip()
            if "```json" in json_str:
                json_str = json_str.split("```json")[1].split("```")[0].strip()
            
            params = json.loads(json_str)
            logger.info(f"✅ [CHECKPOINT] Parameters extracted: {list(params.keys())}")
            return params
        except json.JSONDecodeError as e:
            logger.error(f"❌ [CHECKPOINT] JSON decode error: {e}, response: {response[:200]}")
            return {}  # 返回空字典，使用默认值
        except Exception as e:
            logger.error(f"❌ [CHECKPOINT] Error extracting parameters: {e}", exc_info=True)
            return {}  # 返回空字典，使用默认值
    
    async def _peek_data_lightweight(self, file_path: str) -> Dict[str, Any]:
        """
        轻量级数据预览（使用 FileInspector）
        
        🔧 升级：委托给 FileInspector，获得准确的统计信息
        
        Args:
            file_path: 文件路径
        
        Returns:
            包含基本信息的字典（样本数、列数、数值范围等）
        """
        try:
            # 🔧 使用 FileInspector（Universal Eyes）
            from ...core.file_inspector import FileInspector
            upload_dir = os.getenv("UPLOAD_DIR", "/app/uploads")
            inspector = FileInspector(upload_dir)
            
            # 使用通用检查器
            result = inspector.inspect_file(file_path)
            
            if result.get("status") == "success" and result.get("file_type") == "tabular":
                # 转换为兼容格式
                summary = result.get("data", {}).get("summary", {})
                data_range = result.get("data_range", {})
                
                # 构建兼容格式
                peek_result = {
                    "n_samples": summary.get("n_samples", "N/A"),
                    "n_metabolites": summary.get("n_features", 0),
                    "n_metadata_cols": result.get("n_metadata_cols", 0),
                    "metadata_columns": result.get("metadata_columns", []),
                    "data_range": data_range,  # 🔧 添加数据范围（用于 Log2 判断）
                    "missing_rate": summary.get("missing_rate", 0),
                    "numeric_stats": {
                        "min": data_range.get("min", 0),
                        "max": data_range.get("max", 0),
                        "mean": data_range.get("mean", 0),
                        "median": data_range.get("median", 0),
                        "has_large_values": data_range.get("max", 0) > 1000 if isinstance(data_range.get("max"), (int, float)) else False,
                        "has_negative": data_range.get("min", 0) < 0 if isinstance(data_range.get("min"), (int, float)) else False
                    },
                    "is_sampled": summary.get("is_sampled", False),
                    "file_path": file_path
                }
                
                return peek_result
            else:
                # 检查失败
                logger.warning(f"⚠️ File inspection failed: {result.get('error', 'Unknown error')}")
                return {
                    "error": result.get("error", "File inspection failed"),
                    "n_samples": "N/A",
                    "n_metabolites": 0
                }
                
        except Exception as e:
            logger.error(f"❌ Error in _peek_data_lightweight: {e}", exc_info=True)
            return {
                "error": str(e),
                "n_samples": "N/A",
                "n_metabolites": 0
            }
    
    async def _generate_parameter_recommendations(
        self,
        peek_result: Dict[str, Any],
        query: str
    ) -> Dict[str, Any]:
        """
        基于轻量级预览生成参数推荐
        
        Args:
            peek_result: 轻量级预览结果
            query: 用户查询
        
        Returns:
            推荐字典，包含 summary 和 params
        """
        import json
        
        try:
            # 构建预览摘要
            preview_summary = f"""
数据预览结果：
- 样本数: {peek_result.get('n_samples', 'N/A')}
- 代谢物数: {peek_result.get('n_metabolites', 'N/A')}
- 元数据列数: {peek_result.get('n_metadata_cols', 'N/A')}
- 元数据列: {', '.join(peek_result.get('metadata_columns', []))}
- 数值范围: {peek_result.get('numeric_stats', {})}
"""
            
            # 🔧 升级：只传递统计信息，不传递原始数据行
            summary = peek_result.get("data", {}).get("summary", {})
            data_range = peek_result.get("data_range", {})
            
            # 🔥 Step 1: 使用新的元数据字段（head, columns, separator）
            columns = peek_result.get('columns', [])
            head_data = peek_result.get('head', {})
            separator = peek_result.get('separator', ',')
            file_path = peek_result.get('file_path', 'N/A')
            
            # 构建列信息摘要
            columns_summary = ""
            if columns:
                metadata_cols = peek_result.get('metadata_columns', [])
                feature_cols = [col for col in columns if col not in metadata_cols]
                columns_summary = f"""
- 总列数: {len(columns)}
- 元数据列 ({len(metadata_cols)}): {', '.join(metadata_cols[:5])}{'...' if len(metadata_cols) > 5 else ''}
- 特征列 ({len(feature_cols)}): {', '.join(feature_cols[:10])}{'...' if len(feature_cols) > 10 else ''}
"""
            
            # 构建数据预览摘要（使用 head）
            head_summary = ""
            if head_data:
                head_markdown = head_data.get('markdown', '')
                if head_markdown:
                    # 只显示前3行，避免 prompt 过长
                    head_lines = head_markdown.split('\n')[:4]  # 表头 + 前3行数据
                    head_summary = f"""
- 数据预览（前3行）:
{chr(10).join(head_lines)}
"""
            
            stats_summary = f"""
数据统计信息（基于完整文件或大文件采样）：
- 文件路径: {file_path}
- 分隔符: {separator}
- 样本数: {summary.get('n_samples', 'N/A')}
- 特征数: {summary.get('n_features', 'N/A')}
- 缺失率: {summary.get('missing_rate', 0):.2f}%
- 数据范围:
  * 最小值: {data_range.get('min', 'N/A')}
  * 最大值: {data_range.get('max', 'N/A')}
  * 平均值: {data_range.get('mean', 'N/A'):.2f if isinstance(data_range.get('mean'), (int, float)) else 'N/A'}
  * 中位数: {data_range.get('median', 'N/A'):.2f if isinstance(data_range.get('median'), (int, float)) else 'N/A'}
- 是否采样: {summary.get('is_sampled', False)}
{columns_summary}{head_summary}
"""
            
            prompt = f"""基于数据统计信息，生成参数推荐。

用户查询: {query}

{stats_summary}

请分析数据特征并推荐合适的参数。返回 JSON 格式：
{{
    "summary": "数据特征摘要（1-2句话）",
    "params": {{
        "normalization": {{"value": "log2" | "zscore" | "none", "reason": "推荐理由"}},
        "missing_threshold": {{"value": "0.5", "reason": "推荐理由"}},
        "scale": {{"value": true | false, "reason": "推荐理由"}},
        "n_components": {{"value": "10", "reason": "推荐理由"}}
    }}
}}

重要判断规则：
- **Log2 变换判断**：如果最大值 > 1000 且最小值 >= 0，推荐 "log2"（数据跨度大，需要对数变换）
- **Z-score 标准化**：如果数据已标准化（均值接近0，标准差接近1）或包含负值，推荐 "zscore"
- **缺失值阈值**：根据缺失率推荐，如果缺失率 > 50%，推荐更高的阈值（如 0.7）
- **PCA 主成分数**：根据样本数推荐，通常为 min(10, 样本数/2)
- **缩放（Scale）**：如果数据范围差异大，推荐 true

注意：只基于统计信息（最大值、最小值、缺失率等）进行推荐，不查看原始数据行。
"""
            
            messages = [
                {"role": "system", "content": "You are a bioinformatics expert. Return JSON only."},
                {"role": "user", "content": prompt}
            ]
            
            completion = await self.llm_client.achat(messages, temperature=0.2, max_tokens=800)
            think_content, response = self.llm_client.extract_think_and_content(completion)
            
            # 解析 JSON
            json_str = response.strip()
            if "```json" in json_str:
                json_str = json_str.split("```json")[1].split("```")[0].strip()
            elif "```" in json_str:
                json_str = json_str.split("```")[1].split("```")[0].strip()
            
            recommendation = json.loads(json_str)
            
            # 验证和补充推荐
            if "summary" not in recommendation:
                recommendation["summary"] = f"检测到数据包含 {peek_result.get('n_samples', 'N/A')} 个样本，{peek_result.get('n_metabolites', 'N/A')} 个代谢物。"
            
            if "params" not in recommendation:
                recommendation["params"] = {}
            
            # 确保关键参数存在
            numeric_stats = peek_result.get("numeric_stats", {})
            if "normalization" not in recommendation["params"]:
                if numeric_stats.get("has_large_values", False) and not numeric_stats.get("has_negative", False):
                    recommendation["params"]["normalization"] = {"value": "log2", "reason": "数值跨度大，建议 Log 变换以符合正态分布"}
                else:
                    recommendation["params"]["normalization"] = {"value": "zscore", "reason": "标准 Z-score 标准化"}
            
            if "missing_threshold" not in recommendation["params"]:
                recommendation["params"]["missing_threshold"] = {"value": "0.5", "reason": "标准质控阈值"}
            
            if "scale" not in recommendation["params"]:
                recommendation["params"]["scale"] = {"value": True, "reason": "标准化有助于后续分析"}
            
            if "n_components" not in recommendation["params"]:
                n_samples = peek_result.get("n_samples", 100)
                n_comp = min(10, max(2, n_samples // 2))
                recommendation["params"]["n_components"] = {"value": str(n_comp), "reason": f"根据样本数 ({n_samples}) 推荐"}
            
            return recommendation
            
        except Exception as e:
            logger.error(f"❌ 生成推荐失败: {e}", exc_info=True)
            # 返回默认推荐
            return {
                "summary": f"检测到数据包含 {peek_result.get('n_samples', 'N/A')} 个样本。",
                "params": {
                    "normalization": {"value": "log2", "reason": "默认推荐"},
                    "missing_threshold": {"value": "0.5", "reason": "标准质控阈值"},
                    "scale": {"value": True, "reason": "标准化有助于后续分析"},
                    "n_components": {"value": "10", "reason": "默认值"}
                }
            }
    
    def _apply_recommendations_to_steps(
        self,
        steps: List[Dict[str, Any]],
        recommendation: Dict[str, Any]
    ):
        """
        将推荐值自动填充到步骤参数中
        
        Args:
            steps: 工作流步骤列表
            recommendation: 推荐字典
        """
        if not recommendation or "params" not in recommendation:
            return
        
        rec_params = recommendation["params"]
        
        for step in steps:
            if step.get("step_id") == "preprocess_data":
                # 填充预处理参数
                if "normalization" in rec_params:
                    step["params"]["normalization"] = rec_params["normalization"]["value"]
                if "missing_threshold" in rec_params:
                    step["params"]["missing_threshold"] = rec_params["missing_threshold"]["value"]
                if "scale" in rec_params:
                    step["params"]["scale"] = str(rec_params["scale"]["value"]).lower()
            
            elif step.get("step_id") == "pca_analysis":
                # 填充 PCA 参数
                if "n_components" in rec_params:
                    step["params"]["n_components"] = rec_params["n_components"]["value"]
    
    async def _generate_final_diagnosis(
        self,
        steps_details: List[Dict[str, Any]],
        workflow_config: Dict[str, Any]
    ) -> Optional[str]:
        """
        基于工作流执行结果生成最终诊断报告
        
        Args:
            steps_details: 步骤执行详情列表
            workflow_config: 工作流配置
        
        Returns:
            Markdown 格式的诊断报告，如果失败返回 None
        """
        import json
        
        try:
            # 提取关键结果
            inspection_result = None
            differential_result = None
            pca_result = None
            
            for step_detail in steps_details:
                tool_id = step_detail.get("tool_id")
                step_result = step_detail.get("step_result", {})
                
                if tool_id == "inspect_data":
                    inspection_result = step_result.get("data", {})
                elif tool_id == "differential_analysis":
                    differential_result = step_result.get("data", {})
                elif tool_id == "pca_analysis":
                    pca_result = step_result.get("data", {})
            
            # 🔧 修复：构建结果摘要（修复字段名不匹配问题）
            # 差异分析：工具返回 n_significant 和 n_total，不是 significant_count 和 total_count
            # PCA：工具返回 explained_variance 在顶层，不在 data.summary 中
            results_summary = {
                "workflow_name": workflow_config.get("workflow_name", "Metabolomics Analysis"),
                "steps_completed": len(steps_details),
                "inspection": inspection_result.get("summary", {}) if inspection_result else None,
                "differential_analysis": {
                    "significant_metabolites": differential_result.get("summary", {}).get("n_significant", "N/A") if differential_result else "N/A",
                    "total_metabolites": differential_result.get("summary", {}).get("n_total", "N/A") if differential_result else "N/A"
                } if differential_result else None,
                "pca": {
                    # 🔧 修复：PCA 结果在 step_result.data 中，但 explained_variance 在顶层的 result 中
                    # 需要从步骤详情中获取完整的 result
                    "variance_explained": self._extract_pca_variance_explained(steps_details) if pca_result else "N/A"
                } if pca_result else None
            }
            
            # 格式化结果摘要
            summary_json = json.dumps(results_summary, ensure_ascii=False, indent=2)
            
            # 🔥 修复：严格的数据驱动诊断 prompt（防止幻觉）
            prompt = f"""You are a strict Data Analyst. Generate a concise diagnosis report based ONLY on the execution results below.

执行结果摘要：
{summary_json}

**CRITICAL RULES:**

1. **Fact-Check First**: Look at `differential_analysis.significant_metabolites`.
   - If it is 0 or "N/A", state clearly: "本次分析未发现显著差异代谢物。"
   - DO NOT invent hypotheses or excuses (like "technical noise", "metabolic homeostasis", "biological similarity") unless there is explicit evidence in the QC metrics.
   - DO NOT write long essays about why there might be no differences.

2. **Interpret PCA**: Look at `pca.variance_explained`.
   - If PC1 is very high (>50%), mention it might indicate a strong batch effect or dominant biological factor.
   - If PC1 is low (<20%), mention the data might be highly heterogeneous.

3. **Concise Conclusion**: Keep it short (3-5 sentences max). Do not write a thesis.
   - Focus on what the data shows, not what it might mean theoretically.

4. **Actionable Advice**: If 0 differences found, suggest:
   - "尝试放宽 P 值阈值（如 0.1）"
   - "检查分组标签是否正确"
   - "考虑增加样本量"
   - DO NOT suggest complex biological interpretations without evidence.

**Output Format:**
- Use Simplified Chinese (简体中文)
- Use Markdown format
- Be direct and factual
- Maximum 200 words

**Example of Good Output (when n_significant = 0):**
"本次分析未发现显著差异代谢物（FDR < 0.05, |Log2FC| > 1）。建议：1) 尝试放宽 P 值阈值至 0.1；2) 检查分组标签是否正确；3) 考虑增加样本量以提高统计功效。"

**Example of Bad Output (DO NOT DO THIS):**
"虽然未发现显著差异，但这可能反映了代谢稳态的维持机制，表明两组样本在代谢水平上保持了高度的生物学相似性..." (This is speculation without evidence!)

现在生成诊断报告："""
            
            messages = [
                {"role": "system", "content": "You are a strict Data Analyst. You must base your diagnosis ONLY on the provided data. Do not invent hypotheses or write speculative essays. Be concise and factual. Use Simplified Chinese."},
                {"role": "user", "content": prompt}
            ]
            
            # 🔥 修复：降低 max_tokens 以匹配简洁性要求（最多 200 字）
            completion = await self.llm_client.achat(messages, temperature=0.2, max_tokens=500)
            think_content, response = self.llm_client.extract_think_and_content(completion)
            from ...core.stream_utils import strip_suggestions_from_text
            if response:
                response, _ = strip_suggestions_from_text(response)
            logger.info(f"📝 Generating diagnosis... Result length: {len(response)}")
            logger.info("📝 Diagnosis generated successfully.")
            return response
            
        except Exception as e:
            logger.error(f"❌ 生成最终诊断失败: {e}", exc_info=True)
            return None
    
    def _extract_pca_variance_explained(self, steps_details: List[Dict[str, Any]]) -> str:
        """
        从步骤详情中提取 PCA 解释方差
        
        Args:
            steps_details: 步骤详情列表
        
        Returns:
            解释方差字符串，格式如 "PC1: 45.23%, PC2: 12.56%"
        """
        for step_detail in steps_details:
            if step_detail.get("tool_id") == "pca_analysis":
                step_result = step_detail.get("step_result", {})
                # 优先从 _full_result 中获取
                full_result = step_result.get("_full_result", {})
                if full_result and "explained_variance" in full_result:
                    pc1_var = full_result["explained_variance"].get("PC1", 0) * 100
                    pc2_var = full_result["explained_variance"].get("PC2", 0) * 100
                    return f"PC1: {pc1_var:.2f}%, PC2: {pc2_var:.2f}%"
                # 如果没有，尝试从 data.tables.variance_table 中提取
                elif step_result.get("data", {}).get("tables", {}).get("variance_table"):
                    variance_table = step_result["data"]["tables"]["variance_table"]
                    if variance_table and len(variance_table) > 0:
                        pc1_var = variance_table[0].get("解释方差", variance_table[0].get("Explained Variance", "N/A"))
                        pc2_var = variance_table[1].get("解释方差", variance_table[1].get("Explained Variance", "N/A")) if len(variance_table) > 1 else "N/A"
                        return f"PC1: {pc1_var}, PC2: {pc2_var}"
        return "N/A"
    
    async def _stream_chat_response(
        self,
        query: str,
        file_paths: List[str]
    ) -> AsyncIterator[str]:
        """
        流式聊天响应（支持 ReAct 循环和工具调用）
        """
        context = {
            "context": f"Uploaded files: {', '.join(file_paths) if file_paths else 'None'}",
            "available_tools": ["file_inspect", "metabolomics_preprocess", "metabolomics_pca", "metabolomics_differential_analysis", "metabolomics_volcano"],
            "tool_descriptions": {
                "inspect_data": "检查代谢组学数据文件，返回数据摘要（样本数、代谢物数、缺失值、分组信息等）",
                "preprocess_data": "预处理数据：处理缺失值、标准化、缩放",
                "pca_analysis": "执行主成分分析 (PCA)",
                "differential_analysis": "执行差异代谢物分析（两组比较）",
                "visualize_pca": "生成 PCA 可视化图",
                "visualize_volcano": "生成火山图（Volcano Plot）"
            }
        }
        
        # 如果有文件，强制先检查（符合 SOP）
        inspection_result = None
        if file_paths:
            input_path = file_paths[0]
            try:
                # 使用新工具系统
                inspection_result = inspect_file(input_path)
                if "error" not in inspection_result:
                    # 将检查结果添加到上下文中
                    inspection_summary = f"""
【Data Inspection Completed】
- Samples: {inspection_result.get('n_samples', 'N/A')}
- Metabolites: {inspection_result.get('n_metabolites', 'N/A')}
- Missing values: {inspection_result.get('missing_values', {}).get('percentage', 'N/A')}%
- Group column: {inspection_result.get('group_info', {}).get('column', 'N/A')}
- Groups: {inspection_result.get('group_info', {}).get('groups', {})}
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
                logger.error(f"Error inspecting file: {e}", exc_info=True)
                yield f"⚠️ Warning: Could not inspect file: {str(e)}\n\n"
        
        # 构建增强的用户查询，包含工具说明
        enhanced_query = f"""{query}

【Available Tools】
You have access to:
- inspect_data(file_path): Check data file structure (already executed above if files were provided)
- preprocess_data(file_path, missing_threshold, normalization, scale): Preprocess metabolomics data
- pca_analysis(n_components, file_path): Perform PCA analysis
- differential_analysis(group_column, file_path, method, p_value_threshold, fold_change_threshold): Find differential metabolites
- visualize_pca(group_column, pca_file, pc1, pc2): Generate PCA plot
- visualize_volcano(diff_file, p_value_threshold, fold_change_threshold): Generate volcano plot

【Workflow Rule】
- If user provides CSV files: Inspect first (already done above), then analyze and propose parameters
- Before running any analysis, you MUST have inspected the data first
- Now analyze the inspection results and propose parameters.
"""
        
        # 先输出一个提示，让用户知道系统正在工作
        yield "💭 正在分析您的需求，请稍候...\n\n"
        
        # 流式输出 LLM 响应
        try:
            has_content = False
            async for chunk in self.chat(enhanced_query, context, stream=True):
                if chunk:
                    has_content = True
                    yield chunk
            # 如果没有任何内容输出，说明可能出错了
            if not has_content:
                yield "\n\n⚠️ 抱歉，响应生成出现问题。请重试或检查日志。"
        except Exception as e:
            logger.error(f"❌ 流式响应错误: {e}", exc_info=True)
            yield f"\n\n❌ 错误: {str(e)}\n\n请检查服务器日志获取详细信息。"
    
    async def execute_workflow(
        self,
        workflow_config: Dict[str, Any],
        file_paths: List[str],
        output_dir: str
    ) -> Dict[str, Any]:
        """
        执行代谢组学工作流
        
        Args:
            workflow_config: 工作流配置
            file_paths: 文件路径列表
            output_dir: 输出目录
        
        Returns:
            分析报告
        """
        import os
        
        logger.info("=" * 80)
        logger.info("🚀 [CHECKPOINT] execute_workflow START")
        logger.info(f"📁 Input file paths: {file_paths}")
        logger.info(f"📂 Output directory: {output_dir}")
        logger.info(f"⚙️  Workflow config: {workflow_config.get('workflow_name', 'Unknown')}")
        logger.info("=" * 80)
        
        input_path = file_paths[0] if file_paths else None
        if not input_path:
            error_msg = "No input files provided"
            logger.error(f"❌ [CHECKPOINT] execute_workflow FAILED: {error_msg}")
            raise ValueError(error_msg)
        
        # 检查输入文件是否存在
        logger.info(f"🔍 [CHECKPOINT] Checking input file: {input_path}")
        logger.info(f"   File exists? {os.path.exists(input_path)}")
        if not os.path.exists(input_path):
            error_msg = f"Input file not found: {input_path}"
            logger.error(f"❌ [CHECKPOINT] execute_workflow FAILED: {error_msg}")
            raise FileNotFoundError(error_msg)
        logger.info(f"   File size: {os.path.getsize(input_path)} bytes")
        
        # 设置输出目录
        logger.info(f"📂 [CHECKPOINT] Setting up output directory: {output_dir}")
        if not os.path.exists(output_dir):
            os.makedirs(output_dir, exist_ok=True)
            logger.info(f"   Created output directory: {output_dir}")
        else:
            logger.info(f"   Output directory already exists: {output_dir}")
        
        # 🔥 架构升级：不再需要设置旧工具的输出目录
        logger.info(f"✅ [CHECKPOINT] Output directory set: {output_dir}")
        
        # 执行工作流步骤
        steps = workflow_config.get("steps", [])
        steps_details = []
        final_plot = None
        
        logger.info(f"📋 [CHECKPOINT] Workflow has {len(steps)} steps to execute")
        for idx, step in enumerate(steps, 1):
            logger.info(f"   Step {idx}: {step.get('step_id')} ({step.get('tool_id')})")
        
        try:
            for step in steps:
                step_id = step.get("step_id")
                tool_id = step.get("tool_id")
                params = step.get("params", {})
                
                logger.info(f"🔧 [CHECKPOINT] Executing step {len(steps_details) + 1}/{len(steps)}: {step_id} ({tool_id})")
                logger.info(f"   Step params: {params}")
                
                if tool_id == "inspect_data":
                    # 🔥 Step 3: 确保使用绝对路径
                    file_path_to_inspect = params.get("file_path", input_path)
                    
                    # 如果路径不是绝对路径，转换为绝对路径
                    if not os.path.isabs(file_path_to_inspect):
                        from pathlib import Path
                        upload_dir = os.getenv("UPLOAD_DIR", "/app/uploads")
                        file_path_to_inspect = str(Path(upload_dir) / file_path_to_inspect)
                    
                    # 确保路径存在
                    file_path_to_inspect = os.path.abspath(file_path_to_inspect)
                    
                    logger.info(f"🔍 [CHECKPOINT] inspect_data: Using absolute path: {file_path_to_inspect}")
                    logger.info(f"   File exists? {os.path.exists(file_path_to_inspect)}")
                    if os.path.exists(file_path_to_inspect):
                        logger.info(f"   File size: {os.path.getsize(file_path_to_inspect)} bytes")
                    else:
                        logger.error(f"❌ File not found: {file_path_to_inspect}")
                        logger.error(f"   Original path: {params.get('file_path', input_path)}")
                        logger.error(f"   Upload dir: {os.getenv('UPLOAD_DIR', '/app/uploads')}")
                    
                    # 使用新工具系统
                    result = inspect_file(file_path_to_inspect)
                    logger.info(f"✅ [CHECKPOINT] inspect_data completed: {result.get('status', 'unknown')}")
                    step_result = {
                        "step_name": step.get("desc", step_id),
                        "status": result.get("status", "success"),
                        "logs": f"检查完成: {result.get('n_samples', 'N/A')} 个样本, {result.get('n_metabolites', 'N/A')} 个代谢物",
                        "data": result.get("data", {})  # 包含 preview 和 summary
                    }
                    steps_details.append({
                        "step_id": step_id,
                        "tool_id": tool_id,
                        "name": step.get("desc", step_id),
                        "summary": step_result["logs"],
                        "status": step_result["status"],
                        "step_result": step_result  # 完整的步骤结果
                    })
                
                elif tool_id == "preprocess_data":
                    # 🔥 修复：使用智能路径解析，确保文件能被找到
                    file_path_to_preprocess = params.get("file_path", input_path)
                    
                    # 如果路径不是绝对路径或文件不存在，尝试智能路径解析
                    if not os.path.isabs(file_path_to_preprocess) or not os.path.exists(file_path_to_preprocess):
                        from ...core.file_inspector import FileInspector
                        upload_dir = os.getenv("UPLOAD_DIR", "/app/uploads")
                        inspector = FileInspector(upload_dir)
                        resolved_path, _ = inspector._resolve_actual_path(file_path_to_preprocess)
                        if resolved_path:
                            file_path_to_preprocess = resolved_path
                            logger.info(f"✅ [CHECKPOINT] preprocess_data: Resolved path to: {file_path_to_preprocess}")
                    
                    logger.info(f"🔍 [CHECKPOINT] preprocess_data: Trying to read file at: {file_path_to_preprocess}")
                    logger.info(f"   File exists? {os.path.exists(file_path_to_preprocess)}")
                    logger.info(f"   Parameters: missing_threshold={params.get('missing_threshold', '0.5')}, normalization={params.get('normalization', 'log2')}, scale={params.get('scale', 'true')}")
                    # 使用新工具系统
                    missing_imputation = "min"  # 默认使用最小值填充
                    if float(params.get("missing_threshold", "0.5")) > 0.5:
                        missing_imputation = "zero"  # 如果缺失值超过50%，使用零填充
                    
                    result = preprocess_metabolite_data(
                        file_path=file_path_to_preprocess,
                        missing_imputation=missing_imputation,
                        log_transform=(params.get("normalization", "log2") == "log2"),
                        standardize=params.get("scale", "true").lower() == "true"
                    )
                    
                    # 适配返回格式（保存预处理后的文件）
                    if result.get("status") == "success" and output_dir:
                        import pandas as pd
                        preprocessed_df = pd.DataFrame(result.get("preprocessed_data", {}))
                        output_path = os.path.join(output_dir, "preprocessed_data.csv")
                        preprocessed_df.to_csv(output_path)
                        result["output_path"] = output_path
                    logger.info(f"✅ [CHECKPOINT] preprocess_data completed: {result.get('status', 'unknown')}")
                    step_result = {
                        "step_name": step.get("desc", step_id),
                        "status": result.get("status", "success"),
                        "logs": result.get("message", "预处理完成"),
                        "data": result.get("data", {})  # 包含 preview 和 summary
                    }
                    steps_details.append({
                        "step_id": step_id,
                        "tool_id": tool_id,
                        "name": step.get("desc", step_id),
                        "summary": step_result["logs"],
                        "status": step_result["status"],
                        "step_result": step_result
                    })
                
                elif tool_id == "pca_analysis":
                    # 如果上一步是预处理，使用预处理后的文件
                    preprocessed_file = None
                    for prev_step in steps_details:
                        if prev_step.get("tool_id") == "preprocess_data":
                            # 从预处理结果中获取文件路径
                            preprocessed_file = os.path.join(output_dir, "preprocessed_data.csv")
                            break
                    
                    # 使用新工具系统
                    result = run_pca(
                        file_path=preprocessed_file or params.get("file_path", input_path),
                        n_components=int(params.get("n_components", "10")),
                        scale=True,
                        output_dir=output_dir
                    )
                    # 保存 PCA 结果供后续使用
                    if result.get("status") == "success":
                        self._last_pca_result = result
                    step_result = {
                        "step_name": step.get("desc", step_id),
                        "status": result.get("status", "success"),
                        "logs": result.get("message", "PCA 分析完成"),
                        "data": result.get("data", {}),  # 包含 preview 和 tables
                        "_full_result": result  # 🔧 修复：保存完整结果以便后续提取
                    }
                    steps_details.append({
                        "step_id": step_id,
                        "tool_id": tool_id,
                        "name": step.get("desc", step_id),
                        "summary": step_result["logs"],
                        "status": step_result["status"],
                        "step_result": step_result
                    })
                
                elif tool_id == "differential_analysis":
                    # 使用预处理后的文件
                    preprocessed_file = os.path.join(output_dir, "preprocessed_data.csv")
                    if not os.path.exists(preprocessed_file):
                        preprocessed_file = input_path
                    
                    # 使用新工具系统
                    import pandas as pd
                    df = pd.read_csv(preprocessed_file, index_col=0)
                    group_col = params.get("group_column", "Muscle loss")
                    
                    # 确定分组
                    if group_col in df.columns:
                        groups = df[group_col].unique()
                        case_group = params.get("group1") or (groups[0] if len(groups) > 0 else "Group1")
                        control_group = params.get("group2") or (groups[1] if len(groups) > 1 else groups[0])
                    else:
                        case_group = params.get("group1", "Group1")
                        control_group = params.get("group2", "Group2")
                    
                    result = run_differential_analysis(
                        file_path=preprocessed_file,
                        group_column=group_col,
                        case_group=case_group,
                        control_group=control_group,
                        fdr_method="fdr_bh"
                    )
                    
                    # 保存差异分析结果到文件
                    if result.get("status") == "success" and output_dir:
                        results_df = pd.DataFrame(result.get("results", []))
                        diff_output_path = os.path.join(output_dir, "differential_analysis.csv")
                        results_df.to_csv(diff_output_path, index=False)
                        result["output_path"] = diff_output_path
                    
                    step_result = {
                        "step_name": step.get("desc", step_id),
                        "status": result.get("status", "success"),
                        "logs": result.get("message", "差异分析完成"),
                        "data": result.get("data", {})  # 包含 tables 和 summary
                    }
                    steps_details.append({
                        "step_id": step_id,
                        "tool_id": tool_id,
                        "name": step.get("desc", step_id),
                        "summary": step_result["logs"],
                        "status": step_result["status"],
                        "details": result.get("summary", ""),
                        "step_result": step_result
                    })
                
                elif tool_id == "visualize_pca":
                    # PCA 工具已经返回 plot_path，这里只需要从之前的步骤中提取
                    plot_path = None
                    for prev_step in steps_details:
                        if prev_step.get("tool_id") == "pca_analysis":
                            prev_result = prev_step.get("step_result", {}).get("_full_result", {})
                            plot_path = prev_result.get("plot_path")
                            break
                    
                    if plot_path:
                        # 转换为相对路径（相对于 output_dir）
                        if os.path.isabs(plot_path):
                            relative_plot_path = os.path.relpath(plot_path, output_dir)
                        else:
                            relative_plot_path = plot_path
                        relative_plot_path = relative_plot_path.replace("\\", "/")
                        final_plot = relative_plot_path
                    else:
                        relative_plot_path = None
                        logger.warning("⚠️ PCA 可视化：未找到 PCA 图路径")
                    
                    step_result = {
                        "step_name": step.get("desc", step_id),
                        "status": "success" if plot_path else "warning",
                        "logs": "PCA 可视化完成" if plot_path else "PCA 图未生成",
                        "data": {"images": [relative_plot_path]} if relative_plot_path else {}
                    }
                    
                    steps_details.append({
                        "step_id": step_id,
                        "tool_id": tool_id,
                        "name": step.get("desc", step_id),
                        "summary": "PCA 可视化完成",
                        "status": step_result["status"],
                        "plot": relative_plot_path,
                        "step_result": step_result
                    })
                
                elif tool_id == "visualize_volcano":
                    # 使用差异分析结果文件
                    diff_file = params.get("diff_file") or os.path.join(output_dir, "differential_analysis.csv")
                    if not os.path.exists(diff_file):
                        # 如果差异分析文件不存在，尝试从步骤详情中获取
                        for prev_step in steps_details:
                            if prev_step.get("tool_id") == "differential_analysis":
                                diff_file = os.path.join(output_dir, "differential_analysis.csv")
                                break
                    
                    # 使用新工具系统
                    # 读取差异分析结果
                    import pandas as pd
                    diff_df = pd.read_csv(diff_file)
                    diff_results = {"results": diff_df.to_dict(orient="records")}
                    
                    # 生成火山图
                    volcano_output_path = os.path.join(output_dir, "volcano_plot.png")
                    result = plot_volcano(
                        diff_results=diff_results,
                        output_path=volcano_output_path,
                        fdr_threshold=float(params.get("p_value_threshold", "0.05")),
                        log2fc_threshold=float(params.get("fold_change_threshold", "1.5"))
                    )
                    plot_path = result.get("plot_path") or result.get("plot_file")
                    relative_plot_path = None
                    if plot_path:
                        # 转换为相对路径（相对于 output_dir）
                        if os.path.isabs(plot_path):
                            relative_plot_path = os.path.relpath(plot_path, output_dir)
                        else:
                            relative_plot_path = plot_path
                        # 确保路径使用正斜杠（Web 兼容）
                        relative_plot_path = relative_plot_path.replace("\\", "/")
                        if not final_plot:
                            final_plot = relative_plot_path
                    
                    step_result = {
                        "step_name": step.get("desc", step_id),
                        "status": result.get("status", "success"),
                        "logs": "火山图可视化完成",
                        "data": result.get("data", {})  # 包含 images 数组
                    }
                    # 如果 result.data 中没有 images，添加
                    if "images" not in step_result["data"] and relative_plot_path:
                        step_result["data"]["images"] = [relative_plot_path]
                    
                    steps_details.append({
                        "step_id": step_id,
                        "tool_id": tool_id,
                        "name": step.get("desc", step_id),
                        "summary": "火山图可视化完成",
                        "status": result.get("status", "success"),
                        "plot": relative_plot_path,
                        "step_result": step_result
                    })
            
            # 构建 steps_results 列表（前端可直接使用）
            steps_results = []
            for step_detail in steps_details:
                if "step_result" in step_detail:
                    steps_results.append(step_detail["step_result"])
                else:
                    # 兼容旧格式
                    steps_results.append({
                        "step_name": step_detail.get("name", "Unknown"),
                        "status": step_detail.get("status", "success"),
                        "logs": step_detail.get("summary", ""),
                        "data": {}
                    })
            
            logger.info("=" * 80)
            logger.info("✅ [CHECKPOINT] execute_workflow SUCCESS")
            logger.info(f"📊 Completed {len(steps_details)} steps")
            logger.info(f"📂 Output directory: {output_dir}")
            if final_plot:
                logger.info(f"🖼️  Final plot: {final_plot}")
            logger.info("=" * 80)
            
            # 🔥 Task 1: 生成 AI 诊断（在所有步骤完成后）
            diagnosis = None
            try:
                logger.info("📝 [CHECKPOINT] Generating AI diagnosis from workflow results...")
                diagnosis = await self._generate_final_diagnosis(steps_details, workflow_config)
                if diagnosis:
                    logger.info(f"📝 Diagnosis generated successfully. Result length: {len(diagnosis)}")
                    logger.info(f"📝 Diagnosis preview: {diagnosis[:200]}...")  # 显示前200字符
                else:
                    logger.warning("⚠️ [CHECKPOINT] Diagnosis generation returned None")
            except Exception as diag_err:
                logger.error(f"❌ [CHECKPOINT] Diagnosis generation failed: {diag_err}", exc_info=True)
                diagnosis = "⚠️ 诊断生成失败，但分析已完成。"
            
            # 🔥 构建返回结果
            workflow_result = {
                "status": "success",
                "workflow_name": workflow_config.get("workflow_name", "Metabolomics Analysis"),
                "steps_details": steps_details,  # 保留旧格式以兼容
                "steps_results": steps_results,  # 新的格式，前端可直接使用
                "final_plot": final_plot,
                "output_dir": output_dir,
                "diagnosis": diagnosis or "✅ **分析成功完成！**"  # 🔥 Task 1: 添加诊断字段
            }
            
            # 🔥 清理数据以确保 JSON 序列化安全（处理 Numpy 类型、NaN/Infinity 等）
            logger.info("✅ Workflow finished. Sanitizing data for JSON serialization...")
            sanitized_result = sanitize_for_json(workflow_result)
            logger.info("✅ Data sanitization completed. Returning result to frontend.")
            
            return sanitized_result
            
        except Exception as e:
            import traceback
            error_traceback = traceback.format_exc()
            logger.error("=" * 80)
            logger.error("❌ [CHECKPOINT] execute_workflow FAILED")
            logger.error(f"❌ Error type: {type(e).__name__}")
            logger.error(f"❌ Error message: {str(e)}")
            logger.error(f"❌ Completed steps: {len(steps_details)}/{len(steps)}")
            logger.error(f"📂 Output directory: {output_dir}")
            logger.error("❌ Full traceback:")
            logger.error(error_traceback)
            logger.error("=" * 80)
            
            # 返回错误信息给前端（也需要清理）
            error_result = {
                "status": "error",
                "error": str(e),
                "error_type": type(e).__name__,
                "error_traceback": error_traceback,
                "steps_details": steps_details,
                "output_dir": output_dir,
                "message": f"❌ 工作流执行出错: {str(e)}\n\n(请查看后台日志获取详细堆栈)"
            }
            
            # 清理错误结果（虽然错误结果通常不包含 Numpy 数据，但为了安全起见）
            return sanitize_for_json(error_result)

