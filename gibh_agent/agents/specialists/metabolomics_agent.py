"""代谢组学智能体（Metabolomics Agent）"""
from typing import Dict, Any, List, AsyncIterator, Optional
from ..base_agent import BaseAgent
from ...core.llm_client import LLMClient
from ...core.prompt_manager import PromptManager, DATA_DIAGNOSIS_PROMPT
from ...core.utils import sanitize_for_json
from ...tools.metabolomics_tool import MetabolomicsTool
import logging

logger = logging.getLogger(__name__)


class MetabolomicsAgent(BaseAgent):
    """代谢组学智能体"""
    
    # 定义严格的步骤顺序（依赖链）
    STEPS_ORDER = [
        "inspect_data",      # 步骤1: 数据检查
        "preprocess_data",   # 步骤2: 数据预处理
        "pca_analysis",      # 步骤3: PCA 分析
        "differential_analysis",  # 步骤4: 差异分析
        "visualize_pca",     # 步骤5: PCA 可视化
        "visualize_volcano"  # 步骤6: 火山图可视化
    ]
    
    # 步骤映射（step_id -> 在 STEPS_ORDER 中的索引）
    STEP_INDEX_MAP = {step: idx for idx, step in enumerate(STEPS_ORDER)}
    
    def __init__(
        self,
        llm_client: LLMClient,
        prompt_manager: PromptManager,
        metabolomics_config: Dict[str, Any] = None
    ):
        super().__init__(llm_client, prompt_manager, "metabolomics_expert")
        self.metabolomics_config = metabolomics_config or {}
        self.metabolomics_tool = MetabolomicsTool(self.metabolomics_config)
    
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
        """
        query_lower = query.lower().strip()
        file_paths = self.get_file_paths(uploaded_files or [])
        
        # 判断是否是工作流请求
        is_workflow_request = self._is_workflow_request(query_lower, file_paths)
        
        if is_workflow_request:
            # 工作流请求：先检查数据，然后生成工作流配置
            return await self._generate_workflow_config(query, file_paths)
        else:
            # 普通聊天：流式响应
            return {
                "type": "chat",
                "response": self._stream_chat_response(query, file_paths)
            }
    
    def _is_workflow_request(self, query: str, file_paths: List[str]) -> bool:
        """判断是否是工作流请求"""
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
        
        流程：
        1. 先检查数据（inspect_data）
        2. 使用 LLM 提取目标步骤（支持用户指定"只运行步骤1"等）
        3. 基于检查结果提取参数
        4. 生成工作流配置（只包含目标步骤）
        """
        logger.info("=" * 80)
        logger.info("🚀 [CHECKPOINT] _generate_workflow_config START")
        logger.info(f"   Query: {query}")
        logger.info(f"   File paths: {file_paths}")
        logger.info("=" * 80)
        
        # 强制检查：如果有文件，先检查
        inspection_result = None
        diagnosis_report = None
        if file_paths:
            input_path = file_paths[0]
            logger.info(f"🔍 [CHECKPOINT] Inspecting file: {input_path}")
            try:
                inspection_result = self.metabolomics_tool.inspect_data(input_path)
                if "error" in inspection_result:
                    logger.warning(f"⚠️ File inspection failed: {inspection_result.get('error')}")
                else:
                    logger.info(f"✅ [CHECKPOINT] File inspection successful")
                    # 🔥 生成数据诊断和参数推荐
                    try:
                        logger.info(f"🔍 [CHECKPOINT] Generating diagnosis report...")
                        diagnosis_report = await self._generate_diagnosis_and_recommendation(inspection_result)
                        logger.info(f"✅ [CHECKPOINT] Diagnosis report generated")
                    except Exception as diag_err:
                        logger.error(f"❌ [CHECKPOINT] Diagnosis report generation failed: {diag_err}", exc_info=True)
                        diagnosis_report = None  # 继续执行，不阻塞
            except Exception as e:
                logger.error(f"❌ [CHECKPOINT] Error inspecting file: {e}", exc_info=True)
        
        # 使用 LLM 提取目标结束步骤（例如："做到PCA" -> "pca_analysis"）
        target_end_step = None
        try:
            logger.info(f"🔍 [CHECKPOINT] Extracting target end step from query...")
            target_end_step = await self._extract_target_end_step(query, inspection_result)
            logger.info(f"✅ [CHECKPOINT] Target end step extracted: {target_end_step}")
        except Exception as e:
            logger.error(f"❌ [CHECKPOINT] Error extracting target end step: {e}", exc_info=True)
            target_end_step = None  # 使用默认值（所有步骤）
        
        # 使用 LLM 提取参数（传入检查结果和诊断报告）
        extracted_params = {}
        try:
            logger.info(f"🔍 [CHECKPOINT] Extracting workflow parameters...")
            extracted_params = await self._extract_workflow_params(query, file_paths, inspection_result, diagnosis_report)
            logger.info(f"✅ [CHECKPOINT] Workflow parameters extracted: {list(extracted_params.keys())}")
        except Exception as e:
            logger.error(f"❌ [CHECKPOINT] Error extracting workflow params: {e}", exc_info=True)
            extracted_params = {}  # 使用默认值
        
        # 定义所有可用步骤（包含友好的中文名称）
        all_steps = [
            {
                "step_id": "inspect_data",
                "tool_id": "inspect_data",
                "name": "数据检查",  # 🔧 修复：添加 name 字段
                "step_name": "数据检查",  # 🔧 修复：添加 step_name 字段（兼容前端）
                "desc": "检查数据文件的基本信息（样本数、代谢物数、缺失值、分组信息等）",
                "params": {"file_path": file_paths[0] if file_paths else ""}
            },
            {
                "step_id": "preprocess_data",
                "tool_id": "preprocess_data",
                "name": "数据预处理",  # 🔧 修复：添加 name 字段
                "step_name": "数据预处理",  # 🔧 修复：添加 step_name 字段（兼容前端）
                "desc": "数据预处理：处理缺失值、标准化、缩放",
                "params": {
                    "file_path": file_paths[0] if file_paths else "",
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
                    "group_column": extracted_params.get("group_column", "Muscle loss"),
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
                    "group_column": extracted_params.get("group_column", "Muscle loss"),
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
        
        # 如果生成了诊断报告，将其包含在返回结果中
        result = {
            "type": "workflow_config",
            "workflow_data": workflow_config,
            "file_paths": file_paths
        }
        
        if diagnosis_report:
            result["diagnosis_report"] = diagnosis_report
        
        logger.info("=" * 80)
        logger.info("✅ [CHECKPOINT] _generate_workflow_config SUCCESS")
        logger.info(f"   Workflow name: {workflow_config.get('workflow_name')}")
        logger.info(f"   Steps count: {len(workflow_config.get('steps', []))}")
        logger.info("=" * 80)
        
        return result
    
    async def _generate_diagnosis_and_recommendation(
        self,
        inspection_result: Dict[str, Any]
    ) -> Optional[str]:
        """
        生成数据诊断和参数推荐报告
        
        Args:
            inspection_result: 文件检查结果
        
        Returns:
            Markdown格式的诊断和推荐报告，如果失败返回 None
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
                logger.warning(f"⚠️ 获取诊断模板失败，使用默认模板: {e}")
                prompt = DATA_DIAGNOSIS_PROMPT.format(inspection_data=inspection_json)
            
            # 调用LLM生成诊断报告
            messages = [
                {"role": "system", "content": "You are a Senior Bioinformatician specializing in Metabolomics. Generate data diagnosis and parameter recommendations in Simplified Chinese."},
                {"role": "user", "content": prompt}
            ]
            
            completion = await self.llm_client.achat(messages, temperature=0.3, max_tokens=1500)
            think_content, response = self.llm_client.extract_think_and_content(completion)
            
            logger.info("✅ 数据诊断和参数推荐已生成")
            return response
            
        except Exception as e:
            logger.error(f"❌ 生成诊断报告失败: {e}", exc_info=True)
            return None  # 返回 None，不阻塞工作流生成
    
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
            "available_tools": list(self.metabolomics_tool.tool_map.keys()),
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
                inspection_result = self.metabolomics_tool.inspect_data(input_path)
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
        
        # 更新代谢组工具的输出目录
        from pathlib import Path
        self.metabolomics_config["output_dir"] = output_dir
        self.metabolomics_tool.output_dir = Path(output_dir)
        self.metabolomics_tool.output_dir.mkdir(parents=True, exist_ok=True)
        logger.info(f"✅ [CHECKPOINT] Metabolomics tool output_dir set to: {self.metabolomics_tool.output_dir}")
        
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
                    file_path_to_inspect = params.get("file_path", input_path)
                    logger.info(f"🔍 [CHECKPOINT] inspect_data: Trying to read file at: {file_path_to_inspect}")
                    logger.info(f"   File exists? {os.path.exists(file_path_to_inspect)}")
                    if os.path.exists(file_path_to_inspect):
                        logger.info(f"   File size: {os.path.getsize(file_path_to_inspect)} bytes")
                    result = self.metabolomics_tool.inspect_data(file_path_to_inspect)
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
                    file_path_to_preprocess = params.get("file_path", input_path)
                    logger.info(f"🔍 [CHECKPOINT] preprocess_data: Trying to read file at: {file_path_to_preprocess}")
                    logger.info(f"   File exists? {os.path.exists(file_path_to_preprocess)}")
                    logger.info(f"   Parameters: missing_threshold={params.get('missing_threshold', '0.5')}, normalization={params.get('normalization', 'log2')}, scale={params.get('scale', 'true')}")
                    result = self.metabolomics_tool.preprocess_data(
                        file_path=file_path_to_preprocess,
                        missing_threshold=float(params.get("missing_threshold", "0.5")),
                        normalization=params.get("normalization", "log2"),
                        scale=params.get("scale", "true").lower() == "true"
                    )
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
                    
                    result = self.metabolomics_tool.pca_analysis(
                        n_components=int(params.get("n_components", "10")),
                        file_path=preprocessed_file or params.get("file_path", input_path)
                    )
                    step_result = {
                        "step_name": step.get("desc", step_id),
                        "status": result.get("status", "success"),
                        "logs": result.get("message", "PCA 分析完成"),
                        "data": result.get("data", {})  # 包含 preview 和 tables
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
                    
                    result = self.metabolomics_tool.differential_analysis(
                        group_column=params.get("group_column", "Muscle loss"),
                        file_path=preprocessed_file,
                        method=params.get("method", "t-test"),
                        p_value_threshold=float(params.get("p_value_threshold", "0.05")),
                        fold_change_threshold=float(params.get("fold_change_threshold", "1.5")),
                        group1=params.get("group1"),
                        group2=params.get("group2")
                    )
                    
                    # 检查是否需要用户选择分组
                    if result.get("status") == "need_selection":
                        return {
                            "status": "need_selection",
                            "message": result.get("message", "需要选择比较的分组"),
                            "groups": result.get("groups", []),
                            "steps_details": steps_details
                        }
                    
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
                    # 使用 PCA 分析结果文件
                    pca_file = params.get("pca_file") or os.path.join(output_dir, "pca_results.csv")
                    if not os.path.exists(pca_file):
                        # 如果 PCA 文件不存在，尝试从步骤详情中获取
                        for prev_step in steps_details:
                            if prev_step.get("tool_id") == "pca_analysis":
                                pca_file = os.path.join(output_dir, "pca_results.csv")
                                break
                    
                    result = self.metabolomics_tool.visualize_pca(
                        group_column=params.get("group_column", "Muscle loss"),
                        pca_file=pca_file,
                        pc1=int(params.get("pc1", "1")),
                        pc2=int(params.get("pc2", "2"))
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
                        final_plot = relative_plot_path
                    
                    step_result = {
                        "step_name": step.get("desc", step_id),
                        "status": result.get("status", "success"),
                        "logs": "PCA 可视化完成",
                        "data": result.get("data", {})  # 包含 images 数组
                    }
                    # 如果 result.data 中没有 images，添加
                    if "images" not in step_result["data"] and relative_plot_path:
                        step_result["data"]["images"] = [relative_plot_path]
                    
                    steps_details.append({
                        "step_id": step_id,
                        "tool_id": tool_id,
                        "name": step.get("desc", step_id),
                        "summary": "PCA 可视化完成",
                        "status": result.get("status", "success"),
                        "plot": relative_plot_path,
                        "step_result": step_result
                    })
                
                elif tool_id == "visualize_volcano":
                    # 使用差异分析结果文件
                    diff_file = params.get("diff_file") or os.path.join(output_dir, "differential_results.csv")
                    if not os.path.exists(diff_file):
                        # 如果差异分析文件不存在，尝试从步骤详情中获取
                        for prev_step in steps_details:
                            if prev_step.get("tool_id") == "differential_analysis":
                                diff_file = os.path.join(output_dir, "differential_results.csv")
                                break
                    
                    result = self.metabolomics_tool.visualize_volcano(
                        diff_file=diff_file,
                        p_value_threshold=float(params.get("p_value_threshold", "0.05")),
                        fold_change_threshold=float(params.get("fold_change_threshold", "1.5"))
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
            
            # 🔥 构建返回结果
            workflow_result = {
                "status": "success",
                "workflow_name": workflow_config.get("workflow_name", "Metabolomics Analysis"),
                "steps_details": steps_details,  # 保留旧格式以兼容
                "steps_results": steps_results,  # 新的格式，前端可直接使用
                "final_plot": final_plot,
                "output_dir": output_dir
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

