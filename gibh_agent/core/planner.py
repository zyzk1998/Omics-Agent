"""
动态工作流规划器 - The Brain

使用 Tool-RAG 架构动态生成可执行的工作流计划。
结合工具检索和 LLM 推理，生成符合前端格式的工作流配置。

包含两个规划器：
1. WorkflowPlanner: 通用工作流规划器
2. SOPPlanner: 基于 SOP（标准操作程序）的领域特定规划器（已升级为使用 WorkflowRegistry）
"""
import json
import logging
from typing import Dict, Any, List, Optional
from pathlib import Path

from .tool_retriever import ToolRetriever
from .tool_registry import registry
from .llm_client import LLMClient
from .workflows import WorkflowRegistry

logger = logging.getLogger(__name__)


def _is_flat_params(d: dict) -> bool:
    """判断 recommended_params 是否为扁平键值对（参数名->值），而非 step_id->params 的嵌套。"""
    if not d or not isinstance(d, dict):
        return False
    return not any(isinstance(v, dict) for v in d.values())


# Step semantics for Dynamic Planning: LLM uses these to map "skip X" / "no Y" to step IDs.
SOP_PLANNER_STEP_SEMANTICS = {
    "Spatial": [
        ("load_data", "Load Visium data from Space Ranger output directory"),
        ("spatial_data_validation", "Image & coordinate validation (秀肌肉 first step)"),
        ("qc_norm", "QC and normalization of spots/genes"),
        ("dimensionality_reduction", "PCA dimensionality reduction"),
        ("clustering", "Leiden clustering on the graph"),
        ("spatial_clustering_comparison", "Multi-resolution spatial domain 1x3 comparison (秀肌肉)"),
        ("spatial_neighbors", "Spatial neighborhood graph (for Moran's I)"),
        ("spatial_autocorr", "Spatially variable genes (Moran's I)"),
        ("functional_enrichment", "Pathway enrichment of top SVGs (pathway enrichment)"),
        ("plot_clusters", "Spatial scatter plot colored by clusters"),
        ("plot_genes", "Spatial scatter plot colored by gene expression"),
    ],
    "Radiomics": [
        ("radiomics_data_validation", "Data & ROI validation (秀肌肉 first step)"),
        ("load_image", "Load NIfTI/DICOM medical image"),
        ("preprocess", "Resampling and intensity normalization (preprocessing)"),
        ("preview_slice", "Export mid-slice PNG for preview"),
        ("extract_features", "PyRadiomics texture/shape feature extraction"),
        ("radiomics_model_comparison", "Multi-algorithm LR/SVM/RF ROC comparison (秀肌肉)"),
        ("calc_score", "Rad-Score and risk probability (scoring)"),
        ("viz_score", "Rad-Score visualization (bar/gauge chart)"),
    ],
}


class WorkflowPlanner:
    """
    工作流规划器
    
    职责：
    1. 从用户查询中检索相关工具
    2. 使用 LLM 生成工作流计划
    3. 验证和转换输出格式
    """
    
    def __init__(
        self,
        tool_retriever: ToolRetriever,
        llm_client: LLMClient
    ):
        """
        初始化工作流规划器
        
        Args:
            tool_retriever: 工具检索器实例
            llm_client: LLM 客户端实例
        """
        self.tool_retriever = tool_retriever
        self.llm_client = llm_client
    
    async def plan(
        self,
        user_query: str,
        context_files: List[str] = None,
        category_filter: Optional[str] = None
    ) -> Dict[str, Any]:
        """
        生成工作流计划
        
        Args:
            user_query: 用户查询文本
            context_files: 可用的文件路径列表
            category_filter: 可选的类别过滤器（如 "Metabolomics"）
        
        Returns:
            工作流配置字典，符合前端格式
        """
        try:
            logger.info(f"🧠 开始规划工作流: '{user_query}'")
            
            # Step 1: 检索相关工具
            logger.info("🔍 Step 1: 检索相关工具...")
            retrieved_tools = self.tool_retriever.retrieve(
                query=user_query,
                top_k=10,
                category_filter=category_filter
            )
            
            if not retrieved_tools:
                logger.warning("⚠️ 未检索到相关工具")
                return {
                    "type": "error",
                    "error": "未找到相关工具，请检查查询或工具注册",
                    "message": "无法生成工作流计划"
                }
            
            logger.info(f"✅ 检索到 {len(retrieved_tools)} 个相关工具")
            for tool in retrieved_tools[:3]:  # 只打印前3个
                logger.info(f"   - {tool['name']} (相似度: {tool['similarity_score']:.4f})")
            
            # Step 2: 构建 LLM Prompt
            logger.info("📝 Step 2: 构建 LLM Prompt...")
            system_prompt = self._build_system_prompt()
            user_prompt = self._build_user_prompt(
                user_query=user_query,
                retrieved_tools=retrieved_tools,
                context_files=context_files or []
            )
            
            # Step 3: 调用 LLM 生成计划
            logger.info("🤖 Step 3: 调用 LLM 生成计划...")
            messages = [
                {"role": "system", "content": system_prompt},
                {"role": "user", "content": user_prompt}
            ]
            
            # 尝试使用 JSON mode（如果支持）
            response = await self.llm_client.achat(
                messages=messages,
                temperature=0.1,  # 低温度确保一致性
                max_tokens=2048
            )
            
            # Step 4: 解析 LLM 响应
            logger.info("🔧 Step 4: 解析 LLM 响应...")
            # 🔥 FIX: 提取 ChatCompletion 对象的内容
            if hasattr(response, 'choices') and response.choices:
                response_text = response.choices[0].message.content or ""
            else:
                response_text = str(response)
            workflow_plan = self._parse_llm_response(response_text)
            
            # Step 5: 验证工具存在性
            logger.info("✅ Step 5: 验证工具...")
            validated_plan = self._validate_and_adapt(workflow_plan, context_files or [])
            
            logger.info(f"✅ 工作流规划完成: {len(validated_plan.get('steps', []))} 个步骤")
            return validated_plan
        
        except Exception as e:
            logger.error(f"❌ 工作流规划失败: {e}", exc_info=True)
            return {
                "type": "error",
                "error": str(e),
                "message": f"工作流规划失败: {str(e)}"
            }
    
    def _build_system_prompt(self) -> str:
        """
        构建系统提示词
        
        Returns:
            系统提示词文本
        """
        return """You are a Senior Bioinformatics Workflow Architect.

Your task is to generate executable workflow plans based on user queries and available tools.

**CRITICAL RULES:**
1. Output MUST be a valid JSON object (no markdown code blocks, no extra text).
2. Structure: {"workflow_name": "...", "steps": [...], "recommended_params": {...} (optional)}.
3. **Data Flow**: The output of Step N must match the input of Step N+1 (e.g., file paths).
4. If a parameter is a file path, use placeholders like `<step1_output>` or match the uploaded filename.
5. Use tool names EXACTLY as provided in the tool schemas.
6. Only include parameters that exist in the tool's args_schema.
7. For file paths, prefer using the actual uploaded filename if available.
8. **recommended_params**: You MUST recommend execution parameters based on data scale and output them as a JSON object in "recommended_params". Example: {"n_pcs": 30, "resolution": 0.8, "min_genes": 200} or per-step: {"rna_clustering": {"resolution": 0.8}, "rna_pca": {"n_comps": 30}}. Only include parameter names that exist in the tool's signature; values will be type-coerced by the executor.

**Step Structure:**
- "id": Unique step identifier (e.g., "step1", "step2")
- "tool_name": The exact tool name from the retrieved tools
- "params": Dictionary of parameters matching the tool's args_schema
- "dependency": Optional, the step ID this step depends on (e.g., "step1")

**Example Output:**
{
  "workflow_name": "PCA Analysis Pipeline",
  "steps": [
    {
      "id": "step1",
      "tool_name": "metabolomics_pca",
      "params": {
        "file_path": "cow_diet.csv",
        "n_components": 2,
        "scale": true
      },
      "dependency": null
    }
  ],
  "recommended_params": {"n_components": 10, "resolution": 0.6}
}

Generate ONLY the JSON object, no additional text."""
    
    def _build_user_prompt(
        self,
        user_query: str,
        retrieved_tools: List[Dict[str, Any]],
        context_files: List[str]
    ) -> str:
        """
        构建用户提示词
        
        Args:
            user_query: 用户查询
            retrieved_tools: 检索到的工具列表
            context_files: 可用文件列表
        
        Returns:
            用户提示词文本
        """
        # 格式化工具信息
        tools_text = []
        for i, tool in enumerate(retrieved_tools, 1):
            tool_info = f"""
Tool {i}: {tool['name']}
  Description: {tool['description']}
  Category: {tool['category']}
  Output Type: {tool['output_type']}
  Parameters Schema:
{json.dumps(tool['args_schema'], indent=2, ensure_ascii=False)}
"""
            tools_text.append(tool_info)
        
        # 格式化文件信息
        files_text = ""
        if context_files:
            files_text = "\n**Available Files:**\n"
            for file_path in context_files:
                filename = Path(file_path).name
                files_text += f"- {filename} (path: {file_path})\n"
        else:
            files_text = "\n**Available Files:** None (user may upload files later)\n"
        
        prompt = f"""**User Query:**
{user_query}

**Retrieved Tools:**
{''.join(tools_text)}

{files_text}

**Task:**
Generate a workflow plan that fulfills the user's request. Use the retrieved tools and available files.

**Instructions:**
1. Select the most relevant tools from the retrieved list.
2. Arrange them in a logical order (data flow: output of step N → input of step N+1).
3. Fill in parameters using:
   - Actual file names from "Available Files" if they match the query
   - Default values from the tool's args_schema
   - Placeholders like `<step1_output>` if a step depends on another step's output
4. Keep the workflow concise - only include necessary steps.

**Output Format:**
Return ONLY a valid JSON object with this structure:
{{
  "workflow_name": "Descriptive workflow name",
  "steps": [
    {{
      "id": "step1",
      "tool_name": "exact_tool_name",
      "params": {{"param1": "value1", "param2": "value2"}},
      "dependency": null
    }},
    {{
      "id": "step2",
      "tool_name": "another_tool_name",
      "params": {{"file_path": "<step1_output>", "other_param": "value"}},
      "dependency": "step1"
    }}
  ]
}}

Remember: Output ONLY the JSON object, no markdown, no code blocks, no explanations."""
        
        return prompt
    
    def _parse_llm_response(self, response: str) -> Dict[str, Any]:
        """
        解析 LLM 响应
        
        Args:
            response: LLM 返回的文本
        
        Returns:
            解析后的工作流计划字典
        """
        # 尝试提取 JSON（可能被 markdown 代码块包裹）
        response = response.strip()
        
        # 移除 markdown 代码块标记
        if response.startswith("```json"):
            response = response[7:]
        elif response.startswith("```"):
            response = response[3:]
        
        if response.endswith("```"):
            response = response[:-3]
        
        response = response.strip()
        
        try:
            plan = json.loads(response)
            return plan
        except json.JSONDecodeError as e:
            logger.error(f"❌ JSON 解析失败: {e}")
            logger.error(f"响应内容: {response[:500]}")
            
            # 尝试提取 JSON 对象（使用正则表达式）
            import re
            json_match = re.search(r'\{[^{}]*(?:\{[^{}]*\}[^{}]*)*\}', response, re.DOTALL)
            if json_match:
                try:
                    plan = json.loads(json_match.group())
                    logger.info("✅ 使用正则表达式成功提取 JSON")
                    return plan
                except:
                    pass
            
            raise ValueError(f"无法解析 LLM 响应为 JSON: {e}")
    
    def _validate_and_adapt(
        self,
        workflow_plan: Dict[str, Any],
        context_files: List[str]
    ) -> Dict[str, Any]:
        """
        验证工作流计划并适配为前端格式
        
        Args:
            workflow_plan: LLM 生成的原始计划
            context_files: 可用文件列表
        
        Returns:
            验证和适配后的工作流配置（符合前端格式）
        """
        # 验证基本结构
        if "workflow_name" not in workflow_plan:
            workflow_plan["workflow_name"] = "Generated Workflow"
        
        if "steps" not in workflow_plan or not isinstance(workflow_plan["steps"], list):
            raise ValueError("工作流计划必须包含 'steps' 数组")
        
        # 适配步骤格式
        adapted_steps = []
        for i, step in enumerate(workflow_plan["steps"], 1):
            # 验证工具存在性
            tool_name = step.get("tool_name") or step.get("tool_id")
            if not tool_name:
                logger.warning(f"⚠️ 步骤 {i} 缺少 tool_name，跳过")
                continue
            
            # 检查工具是否在注册表中
            tool_metadata = registry.get_metadata(tool_name)
            if not tool_metadata:
                logger.warning(f"⚠️ 工具 '{tool_name}' 不在注册表中，跳过")
                continue
            
            # 获取工具描述
            tool_desc = tool_metadata.description
            
            # 处理文件路径参数
            params = step.get("params", {})
            adapted_params = {}
            
            for param_name, param_value in params.items():
                # 如果参数值是文件路径占位符，尝试匹配实际文件
                if isinstance(param_value, str) and param_value.startswith("<") and param_value.endswith(">"):
                    # 占位符，保持原样（执行时会处理）
                    adapted_params[param_name] = param_value
                elif param_name == "file_path" and context_files:
                    # 如果是 file_path 参数，尝试匹配上传的文件
                    if isinstance(param_value, str):
                        # 尝试匹配文件名
                        matched_file = None
                        for file_path in context_files:
                            if param_value.lower() in Path(file_path).name.lower():
                                matched_file = file_path
                                break
                        
                        if matched_file:
                            adapted_params[param_name] = matched_file
                        else:
                            # 使用第一个文件作为默认值
                            adapted_params[param_name] = context_files[0]
                    else:
                        adapted_params[param_name] = param_value
                else:
                    adapted_params[param_name] = param_value
            
            # 构建适配后的步骤（符合前端格式）
            # 🔥 Frontend Contract: 必须包含 step_id, tool_id, name, step_name, desc
            step_display_name = self._get_step_display_name(tool_name, tool_desc)
            step_desc = tool_desc[:100] if tool_desc else ""
            
            adapted_step = {
                "step_id": tool_name,  # 🔥 Frontend requires step_id (must match tool_id)
                "id": tool_name,  # 保留 id 作为兼容字段
                "tool_id": tool_name,  # Frontend requires tool_id
                "name": step_display_name,  # Frontend requires name
                "step_name": step_display_name,  # Frontend requires step_name (compatibility)
                "description": tool_desc if tool_desc else "",  # 完整描述
                "desc": step_desc,  # Frontend requires desc (truncated to 100 chars)
                "selected": True,  # Frontend may use this
                "params": adapted_params
            }
            # 动态参数推荐：若 LLM 输出了 recommended_params（按 step_id/tool_id 或全局），则挂到步骤上供 Executor 注入
            rec = workflow_plan.get("recommended_params") or {}
            if isinstance(rec, dict):
                step_rec = rec.get(tool_name) or rec.get(step.get("id")) or (rec if _is_flat_params(rec) else {})
                if step_rec and isinstance(step_rec, dict):
                    adapted_step["recommended_params"] = step_rec
            adapted_steps.append(adapted_step)
        
        workflow_data_inner = {
            "workflow_name": workflow_plan["workflow_name"],
            "steps": adapted_steps
        }
        if workflow_plan.get("recommended_params"):
            workflow_data_inner["recommended_params"] = workflow_plan["recommended_params"]
        # 构建最终的工作流配置（符合前端格式）
        workflow_config = {
            "type": "workflow_config",
            "workflow_data": workflow_data_inner,
            "file_paths": context_files
        }
        # 前端“参数推荐”卡片：从 recommended_params 生成 recommendation（前端展示与后端执行一致）
        rec = workflow_plan.get("recommended_params")
        if isinstance(rec, dict) and rec and _is_flat_params(rec):
            workflow_config["recommendation"] = {
                "summary": "基于数据特征生成的参数推荐",
                "params": {
                    k: {"value": v, "reason": "基于数据特征推荐"}
                    for k, v in rec.items()
                },
            }
        return workflow_config
    
    def _get_step_display_name(self, tool_name: str, tool_desc: str) -> str:
        """
        获取步骤的显示名称
        
        Args:
            tool_name: 工具名称
            tool_desc: 工具描述
        
        Returns:
            显示名称
        """
        # 简单的名称映射（含各模态新工具）
        name_mapping = {
            "metabolomics_pca": "主成分分析 (PCA)",
            "metabolomics_differential_analysis": "差异分析",
            "metabolomics_preprocess": "数据预处理",
            "file_inspect": "文件检查",
            "rna_data_validation": "正在执行底层数据结构与内存预检...",
            "rna_clustering_comparison": "正在生成多分辨率聚类鲁棒性对比图...",
            "metabo_data_validation": "正在执行底层数据与稀疏性校验...",
            "metabo_model_comparison": "正在执行多维机器学习模型 (PCA/PLS-DA) 效能对比...",
            "spatial_data_validation": "正在执行空间影像与坐标校验...",
            "spatial_clustering_comparison": "正在生成多分辨率空间域物理映射对比图...",
            "radiomics_data_validation": "正在执行数据与 ROI 校验...",
            "radiomics_model_comparison": "正在执行多算法 (LR/SVM/RF) ROC 诊断效能对比...",
        }
        
        if tool_name in name_mapping:
            return name_mapping[tool_name]
        
        # 从描述中提取（简单启发式）
        if "PCA" in tool_desc or "principal component" in tool_desc.lower():
            return "主成分分析"
        elif "differential" in tool_desc.lower():
            return "差异分析"
        elif "preprocess" in tool_desc.lower() or "预处理" in tool_desc:
            return "数据预处理"
        elif "inspect" in tool_desc.lower() or "检查" in tool_desc:
            return "文件检查"
        
        # 默认：使用工具名称（美化）
        return tool_name.replace("_", " ").title()



class SOPPlanner:
    """
    SOP 驱动的动态规划器（已升级）
    
    基于标准操作程序（SOP）规则，使用 LLM 生成符合专业流程的工作流计划。
    
    🔥 ARCHITECTURAL UPGRADE:
    - 使用 WorkflowRegistry 进行严格的域绑定（只支持 Metabolomics 和 RNA）
    - 使用 DAG 依赖解析（代码逻辑，非 LLM 幻觉）
    - 支持计划优先（plan-first）：可以在没有文件的情况下生成工作流
    """
    
    def __init__(
        self,
        tool_retriever: ToolRetriever,
        llm_client: LLMClient
    ):
        """
        初始化 SOP 规划器
        
        Args:
            tool_retriever: 工具检索器实例
            llm_client: LLM 客户端实例
        """
        self.tool_retriever = tool_retriever
        self.llm_client = llm_client
        self.workflow_registry = WorkflowRegistry()
    
    async def generate_plan(
        self,
        user_query: str,
        file_metadata: Optional[Dict[str, Any]] = None,
        category_filter: Optional[str] = None,
        domain_name: Optional[str] = None,
        target_steps: Optional[List[str]] = None,
        is_template: bool = False  # 🔥 ARCHITECTURAL RESET: Explicit template flag
    ) -> Dict[str, Any]:
        """
        生成基于 SOP 规则的工作流计划（架构重置版 - 严格分离执行和预览）
        
        🔥 ARCHITECTURAL RESET: Strict Separation of Execution and Preview
        - If is_template=False: MUST use _fill_parameters, MUST return template_mode=False
        - If is_template=True: MUST use _fill_placeholders, MUST return template_mode=True
        
        Args:
            user_query: 用户查询文本
            file_metadata: FileInspector 返回的文件元数据（可选）
            category_filter: 工具类别过滤器（可选）
            domain_name: 可选的域名（如果提供，跳过意图分类）
            target_steps: 可选的目标步骤列表（如果提供，跳过意图分析）
            is_template: 是否为模板模式（True=预览，False=执行）
        
        Returns:
            符合前端格式的工作流配置字典
        """
        try:
            logger.info(f"🧠 [SOPPlanner] 开始生成计划: '{user_query}'")
            
            # Step 1: Intent Classification (LLM) - 识别域名和模式（如果未提供）
            execution_mode = None  # 🔥 NEW: Track execution mode from intent classification
            
            # 🔥 CRITICAL FIX: 如果有 file_metadata，默认应该是 EXECUTION 模式（除非用户明确要求预览）
            has_file_metadata = file_metadata is not None
            
            if not domain_name:
                logger.info("🔍 [SOPPlanner] Step 1: 意图分类（识别域名和模式）...")
                logger.info(f"🔍 [SOPPlanner] file_metadata 存在: {has_file_metadata}")
                intent_result = await self._classify_intent(user_query, file_metadata)
                domain_name = intent_result.get("domain_name")
                execution_mode = intent_result.get("mode", "PLANNING")  # 🔥 NEW: Extract mode
                
                # 🔥 CRITICAL FIX: 如果 file_metadata 存在但模式是 PLANNING，检查是否是明确的预览请求
                if has_file_metadata and execution_mode == "PLANNING":
                    query_lower = user_query.lower()
                    preview_keywords = ["preview", "预览", "show", "显示", "查看", "plan", "规划", "what", "什么"]
                    if not any(kw in query_lower for kw in preview_keywords):
                        # 没有明确的预览关键词，但有文件 -> 强制 EXECUTION
                        logger.warning(f"⚠️ [SOPPlanner] 有文件但模式是 PLANNING，且无预览关键词，强制设置为 EXECUTION")
                        execution_mode = "EXECUTION"
                
                logger.info(f"✅ [SOPPlanner] 意图分类结果: domain={domain_name}, mode={execution_mode}")
            else:
                logger.info(f"✅ [SOPPlanner] 使用提供的域名: {domain_name}")
                # 🔥 CRITICAL: If domain_name provided but no mode, infer from file_metadata
                if has_file_metadata:
                    # Default to EXECUTION if file exists and domain provided
                    execution_mode = "EXECUTION"
                    logger.info(f"✅ [SOPPlanner] 有 file_metadata，默认设置为 EXECUTION 模式")
                else:
                    execution_mode = "PLANNING"
                    logger.info(f"✅ [SOPPlanner] 无 file_metadata，设置为 PLANNING 模式")
            
            # Step 2: 严格域绑定检查
            if not self.workflow_registry.is_supported(domain_name):
                logger.warning(f"⚠️ [SOPPlanner] 不支持的域名: {domain_name}")
                return self.workflow_registry.get_unsupported_error(domain_name)
            
            # Step 3: 获取工作流实例
            workflow = self.workflow_registry.get_workflow(domain_name)
            if not workflow:
                return {
                    "type": "error",
                    "error": f"无法获取工作流: {domain_name}",
                    "message": "工作流注册表错误"
                }
            
            # 🔥 ARCHITECTURAL FIX: 优先运行意图分析（Plan-First）
            # Step 4: Analyze User Intent (LLM) - 从可用工具集中选择目标步骤（如果未提供）
            steps_to_skip: List[str] = []
            if target_steps is None:
                logger.info("🔍 [SOPPlanner] Step 2: 分析用户意图（选择目标步骤）...")
                intent_result = await self._analyze_user_intent(user_query, workflow)
                if isinstance(intent_result, dict):
                    steps_to_skip = intent_result.get("skip_steps") or []
                    target_steps = intent_result.get("target_steps") or []
                else:
                    target_steps = intent_result if isinstance(intent_result, list) else []
            else:
                logger.info(f"✅ [SOPPlanner] 使用提供的目标步骤: {target_steps}")
            
            # 🔥 CRITICAL FIX: 确保 target_steps 不为空（Plan-First 必须返回完整标准流程）
            # 如果查询模糊（如"Analyze this", "Full analysis", "完整分析"），使用完整 SOP
            # 否则，使用用户明确请求的步骤（即使只有一个）
            # 🔥 URGENT: 如果没有文件且查询是规划类（"Plan", "预览", "show me"），默认使用完整 SOP
            if not target_steps:
                query_lower = user_query.lower()
                vague_keywords = ["analyze this", "full analysis", "完整分析", "全部", "all", "complete", "help me analyze", "帮我分析"]
                planning_keywords = ["plan", "预览", "show me", "显示", "生成", "规划", "workflow", "流程"]
                
                if any(kw in query_lower for kw in vague_keywords):
                    logger.info("ℹ️ [SOPPlanner] 查询明确要求完整分析，使用完整 SOP")
                    target_steps = list(workflow.steps_dag.keys())
                elif not file_metadata and any(kw in query_lower for kw in planning_keywords):
                    # 🔥 CRITICAL FIX: Plan-First 模式（无文件），默认返回完整标准流程
                    logger.info("ℹ️ [SOPPlanner] Plan-First 模式（无文件），使用完整标准流程")
                    target_steps = list(workflow.steps_dag.keys())
                else:
                    # 如果意图分析失败，尝试回退关键词匹配
                    logger.info("ℹ️ [SOPPlanner] 意图分析未返回步骤，尝试关键词匹配...")
                    target_steps = self._fallback_intent_analysis(user_query, list(workflow.steps_dag.keys()))
                    if not target_steps:
                        # 🔥 CRITICAL FIX: 如果所有回退都失败，使用完整 SOP（确保不为空）
                        logger.info("ℹ️ [SOPPlanner] 关键词匹配也失败，使用完整 SOP（确保 Plan-First 返回完整流程）")
                        target_steps = list(workflow.steps_dag.keys())
            
            # 🔥 CRITICAL FIX: 再次确保 target_steps 不为空
            if not target_steps:
                logger.warning("⚠️ [SOPPlanner] target_steps 仍然为空，强制使用完整 SOP")
                target_steps = list(workflow.steps_dag.keys())
            
            logger.info(f"✅ [SOPPlanner] 目标步骤: {target_steps} (共 {len(target_steps)} 个)")
            
            # Step 5: Resolve Dependencies (Code) - 使用硬编码 DAG 解析依赖
            logger.info("🔍 [SOPPlanner] Step 3: 解析依赖关系...")
            resolved_steps = self._resolve_dependencies(target_steps, workflow)
            
            # 🔥 CRITICAL FIX: 确保 resolved_steps 不为空
            if not resolved_steps:
                logger.error(f"❌ [SOPPlanner] resolve_dependencies 返回空列表！target_steps: {target_steps}")
                logger.warning("⚠️ [SOPPlanner] 强制使用完整 SOP")
                resolved_steps = list(workflow.steps_dag.keys())
            
            logger.info(f"✅ [SOPPlanner] 依赖解析完成: {target_steps} -> {resolved_steps} (共 {len(resolved_steps)} 个)")
            
            # Step 6: Generate Template - 生成工作流模板（支持占位符）
            logger.info("🔍 [SOPPlanner] Step 4: 生成工作流模板...")
            workflow_config = workflow.generate_template(
                target_steps=resolved_steps,
                file_metadata=file_metadata
            )
            
            # 🔥 Dynamic Planning: Uncheck steps the user asked to skip (selected: false)
            if steps_to_skip:
                wd = workflow_config.get("workflow_data") or {}
                for step in (wd.get("steps") or []):
                    if step.get("id") in steps_to_skip:
                        step["selected"] = False
                        logger.info(f"✅ [SOPPlanner] 已取消勾选步骤: {step.get('id')}")
            
            # 🔥 URGENT FIX: 验证生成的模板包含步骤
            steps_count = len(workflow_config.get('workflow_data', {}).get('steps', []))
            if steps_count == 0:
                logger.error(f"❌ [SOPPlanner] generate_template 返回空步骤！resolved_steps: {resolved_steps}")
                logger.error(f"❌ [SOPPlanner] workflow_config: {workflow_config}")
                # 尝试修复：如果 resolved_steps 不为空但模板为空，可能是 generate_template 实现有问题
                if resolved_steps:
                    logger.warning(f"⚠️ [SOPPlanner] resolved_steps 不为空但模板为空，尝试手动构建步骤...")
                    # 这里不应该手动构建，应该修复 generate_template
                    # 但为了不破坏流程，我们返回错误
                    return {
                        "type": "error",
                        "error": "工作流模板生成失败",
                        "message": f"无法生成工作流模板：步骤列表为空。resolved_steps: {resolved_steps}",
                        "workflow_data": workflow_config.get("workflow_data", {})
                    }
            
            logger.info(f"✅ [SOPPlanner] 模板生成成功: {steps_count} 个步骤")
            
            # 🔥 ARCHITECTURAL RESET: Step 7 - Strict Separation Based on is_template Flag
            logger.info("🔍 [SOPPlanner] Step 5: 处理元数据...")
            logger.info(f"🔍 [SOPPlanner] is_template={is_template}, file_metadata存在={file_metadata is not None}")
            
            # 🔥 CRITICAL: Remove ambiguity - Use is_template flag explicitly
            if is_template:
                # TEMPLATE MODE: Use placeholders, set template_mode = True
                workflow_config = self._fill_placeholders(workflow_config, user_query)
                logger.info("✅ [SOPPlanner] TEMPLATE 模式：已使用占位符，template_mode = True")
            else:
                # EXECUTION MODE: MUST use _fill_parameters, MUST return template_mode = False
                if not file_metadata:
                    logger.error("❌ [SOPPlanner] EXECUTION 模式但 file_metadata 不存在！这是逻辑错误。")
                    # Fallback: Use placeholders but log error
                    workflow_config = self._fill_placeholders(workflow_config, user_query)
                    logger.warning("⚠️ [SOPPlanner] 回退到占位符模式（但这是错误的）")
                else:
                    workflow_config = self._fill_parameters(workflow_config, file_metadata, workflow, template_mode=False)
                    logger.info("✅ [SOPPlanner] EXECUTION 模式：已填充真实参数，template_mode = False")
            
                    # 🔥 CRITICAL: Validate that file_path in params is NOT <PENDING_UPLOAD>
                    steps = workflow_config.get("workflow_data", {}).get("steps", [])
                    for step in steps:
                        params = step.get("params", {})
                        for param_name in ["file_path", "adata_path"]:
                            if param_name in params:
                                param_value = params[param_name]
                                if param_value in ["<待上传数据>", "<PENDING_UPLOAD>", ""]:
                                    logger.error(f"❌ [SOPPlanner] EXECUTION 模式但步骤 {step.get('id')} 的参数 {param_name} 仍是占位符: {param_value}")
                                    # Try to fix: use file_path from metadata
                                    if file_metadata.get("file_path"):
                                        params[param_name] = file_metadata.get("file_path")
                                        logger.warning(f"⚠️ [SOPPlanner] 已修复：将 {param_name} 设置为 {file_metadata.get('file_path')}")
            
            # 🔥 ARCHITECTURAL RESET: Final Validation - Enforce is_template flag
            final_steps_count = len(workflow_config.get('workflow_data', {}).get('steps', []))
            template_mode = workflow_config.get("template_mode", False)
            
            # 🔥 CRITICAL: Force validation based on is_template flag
            if is_template:
                # TEMPLATE MODE: MUST be True
                if not template_mode:
                    logger.warning(f"⚠️ [SOPPlanner] is_template=True 但 template_mode=False，强制设置为 True")
                    template_mode = True
                    workflow_config["template_mode"] = True
                    if "workflow_data" in workflow_config:
                        workflow_config["workflow_data"]["template_mode"] = True
            else:
                # EXECUTION MODE: MUST be False
                if template_mode:
                    logger.error(f"❌ [SOPPlanner] is_template=False 但 template_mode=True，这是逻辑错误！强制设置为 False")
                    template_mode = False
                    workflow_config["template_mode"] = False
                    if "workflow_data" in workflow_config:
                        workflow_config["workflow_data"]["template_mode"] = False
                
                # 🔥 CRITICAL: Validate file paths are NOT placeholders
                if file_metadata:
                    steps = workflow_config.get("workflow_data", {}).get("steps", [])
                    for step in steps:
                        params = step.get("params", {})
                        for param_name in ["file_path", "adata_path"]:
                            if param_name in params:
                                param_value = params[param_name]
                                if param_value in ["<待上传数据>", "<PENDING_UPLOAD>", ""]:
                                    logger.error(f"❌ [SOPPlanner] EXECUTION 模式但步骤 {step.get('id')} 的参数 {param_name} 仍是占位符")
                                    # Try to fix
                                    if file_metadata.get("file_path"):
                                        params[param_name] = file_metadata.get("file_path")
                                        logger.warning(f"⚠️ [SOPPlanner] 已修复：将 {param_name} 设置为 {file_metadata.get('file_path')}")
            
            logger.info(f"✅ [SOPPlanner] 工作流规划完成: {final_steps_count} 个步骤, template_mode = {template_mode}")
            logger.info(f"✅ [SOPPlanner] file_metadata 存在: {file_metadata is not None}")
            
            if final_steps_count == 0:
                logger.error(f"❌ [SOPPlanner] 最终工作流配置步骤为空！")
                return {
                    "type": "error",
                    "error": "工作流步骤为空",
                    "message": "工作流规划失败：最终步骤列表为空。",
                    "workflow_data": workflow_config.get("workflow_data", {})
                }
            
            # 🔥 CRITICAL FIX: 构建返回结果，清理冗余字段
            result = {
                "type": "workflow_config",
                "workflow_data": workflow_config.get("workflow_data"),
                "template_mode": template_mode
            }
            
            # 🔥 CRITICAL: 只在模板模式时包含诊断消息
            # 如果有文件，不包含诊断（让 Orchestrator 从 file_metadata 生成真实诊断）
            if template_mode and "diagnosis" in workflow_config:
                result["diagnosis"] = workflow_config["diagnosis"]
            
            return result
        
        except Exception as e:
            import traceback
            tb = traceback.format_exc()
            logger.error("❌ [SOPPlanner] 工作流规划失败 | 根因: %s | traceback:\n%s", e, tb, exc_info=False)
            return {
                "type": "error",
                "error": str(e),
                "message": f"工作流规划失败: {str(e)}",
            }
    
    def _get_domain_knowledge(self, workflow: "BaseWorkflow") -> str:
        """
        返回当前工作流各步骤的领域描述，供 LLM 理解“跳过/仅运行”等意图。
        
        Spatial / Radiomics 步骤包含详细语义，便于 Dynamic Planning（如 skip preprocessing、no pathway enrichment）。
        """
        steps_dag = workflow.steps_dag if hasattr(workflow, "steps_dag") else getattr(workflow, "get_steps_dag", lambda: {})()
        steps_dag = steps_dag if isinstance(steps_dag, dict) else {}
        domain_name = getattr(workflow, "get_name", lambda: "Unknown")()
        lines = [f"**Domain:** {domain_name}", ""]
        for step_id in steps_dag.keys():
            meta = workflow.get_step_metadata(step_id) if hasattr(workflow, "get_step_metadata") else {}
            name = meta.get("name", step_id)
            desc = (meta.get("description") or "")[:160]
            lines.append(f"- **{step_id}**: {name}. {desc}")
        # Append explicit step semantics for Spatial/Radiomics so LLM knows what to skip (e.g. "pathway enrichment" -> functional_enrichment, "scoring" -> calc_score/viz_score)
        semantics = SOP_PLANNER_STEP_SEMANTICS.get(domain_name)
        if semantics:
            lines.append("")
            lines.append("**Step semantics (for skip/keep intent):**")
            for sid, hint in semantics:
                if sid in steps_dag:
                    lines.append(f"- {sid}: {hint}")
        return "\n".join(lines)

    async def _analyze_user_intent(
        self,
        user_query: str,
        workflow: "BaseWorkflow"
    ):
        """
        分析用户意图，从可用工具集中选择目标步骤，并识别需要“跳过”的步骤。
        
        支持 Dynamic Planning：用户说 "skip preprocessing" / "不要通路富集" 时，
        返回 skip_steps，规划器会将对应步骤设为 selected: false。
        
        Returns:
            dict: {"target_steps": ["step1", ...], "skip_steps": ["step_x"]}
            或 list: 兼容旧版，表示 target_steps，skip_steps 视为 []。
        """
        # 获取可用步骤列表
        available_steps = list(workflow.steps_dag.keys())
        domain_knowledge = self._get_domain_knowledge(workflow)
        
        # 构建提示词（含领域步骤说明 + 跳过/仅运行 few-shot）
        system_prompt = """You are an Intent Analyzer for Bioinformatics Workflows.

Your task is to output TWO lists:
1. **target_steps**: Steps to include (same as before). Use [] for "full workflow".
2. **skip_steps**: Steps the user explicitly asked to SKIP (e.g. "skip preprocessing", "不要聚类", "no pathway enrichment").

**CRITICAL - User Intent:**
- "Analyze this but skip X" / "全部分析但跳过X" -> target_steps = [] (full), skip_steps = ["X"]
- "Just do A and B, no C" -> target_steps = ["A", "B"] (or minimal set including deps), skip_steps = ["C"] if C is in workflow
- Specific only (e.g. "PCA only") -> target_steps = ["pca_analysis"], skip_steps = []
- Vague ("Analyze this", "完整分析") -> target_steps = [], skip_steps = []

**Output Format (strict JSON object only):**
{"target_steps": ["step1", "step2", ...], "skip_steps": ["step_to_skip", ...]}

- target_steps: only use step IDs from the Available Steps list. Empty [] = full workflow.
- skip_steps: only use step IDs from the Available Steps list. Empty [] = no skip.

**Few-Shot (Spatial):**
- User: "Run spatial analysis but I don't need pathway enrichment." -> {"target_steps": [], "skip_steps": ["functional_enrichment"]}
- User: "Analyze this Visium data, skip clustering." -> {"target_steps": [], "skip_steps": ["clustering"]}

**Few-Shot (Radiomics):**
- User: "Just extract features from this CT, no scoring needed." -> {"target_steps": ["load_image", "preprocess", "preview_slice", "extract_features"], "skip_steps": ["calc_score", "viz_score"]}
- User: "Analyze this Radiomics data but skip the preprocessing step." -> {"target_steps": [], "skip_steps": ["preprocess"]}

**Few-Shot (Metabolomics/RNA):**
- "Do PCA only." -> {"target_steps": ["pca_analysis"], "skip_steps": []}
- "Full analysis." -> {"target_steps": [], "skip_steps": []}"""

        user_prompt = f"""**User Query:**
{user_query}

**Step Definitions (use these IDs in target_steps and skip_steps):**
{domain_knowledge}

**Available Step IDs:**
{json.dumps(available_steps, ensure_ascii=False, indent=2)}

**Task:**
From the user query, output a JSON object with "target_steps" and "skip_steps". Use only step IDs from the list above.

**Output (JSON only):**
{{"target_steps": [...], "skip_steps": [...]}}"""

        messages = [
            {"role": "system", "content": system_prompt},
            {"role": "user", "content": user_prompt}
        ]
        
        response = await self.llm_client.achat(
            messages=messages,
            temperature=0.1,
            max_tokens=512
        )
        
        # 🔥 FIX: 提取 ChatCompletion 对象的内容
        if hasattr(response, 'choices') and response.choices:
            response_text = response.choices[0].message.content or ""
        else:
            response_text = str(response)
        
        # 解析响应：支持 {"target_steps": [...], "skip_steps": [...]} 或旧版纯数组
        try:
            raw = json.loads(response_text.strip())
            if isinstance(raw, dict):
                target_steps = raw.get("target_steps")
                skip_steps = raw.get("skip_steps")
                if not isinstance(target_steps, list):
                    target_steps = []
                if not isinstance(skip_steps, list):
                    skip_steps = []
            else:
                target_steps = raw if isinstance(raw, list) else []
                skip_steps = []
            
            valid_steps = [s for s in target_steps if s in available_steps]
            valid_skip = [s for s in skip_steps if s in available_steps]
            invalid = set(target_steps) - set(available_steps) | (set(skip_steps) - set(available_steps))
            if invalid:
                logger.warning(f"⚠️ LLM 返回了无效的步骤ID: {invalid}，已过滤")
            
            logger.info(f"✅ [SOPPlanner] 意图分析完成: target_steps={valid_steps}, skip_steps={valid_skip}")
            return {"target_steps": valid_steps, "skip_steps": valid_skip}
        except json.JSONDecodeError as e:
            logger.error(f"❌ 意图分析 JSON 解析失败: {e}")
            logger.error(f"响应内容: {response_text[:200] if 'response_text' in locals() else str(response)[:200]}")
            fallback_list = self._fallback_intent_analysis(user_query, available_steps)
            return {"target_steps": fallback_list, "skip_steps": []}
    
    def _fallback_intent_analysis(self, user_query: str, available_steps: List[str]) -> List[str]:
        """
        回退的意图分析（基于关键词），含 Spatial / Radiomics 步骤关键词。
        
        Returns:
            目标步骤列表（供 plan 使用；意图分析异常时由 _analyze_user_intent 包装为 dict）
        """
        query_lower = user_query.lower()
        
        step_keywords = {
            "pca_analysis": ["pca", "主成分", "主成分分析", "principal component"],
            "differential_analysis": ["differential", "差异", "diff", "差异分析"],
            "metabolomics_plsda": ["plsda", "pls-da", "pls da", "监督分析"],
            "visualize_volcano": ["volcano", "火山图"],
            "metabolomics_pathway_enrichment": ["pathway", "通路", "enrichment", "富集", "通路富集"],
            "preprocess_data": ["preprocess", "预处理"],
            "inspect_data": ["inspect", "检查", "数据检查"],
            "metabo_data_validation": ["validation", "校验", "数据校验"],
            "metabo_model_comparison": ["model comparison", "模型对比", "pca plsda"],
            "rna_data_validation": ["validation", "预检", "数据校验"],
            "rna_clustering_comparison": ["clustering comparison", "多分辨率", "聚类对比"],
            # Spatial
            "load_data": ["load", "加载"],
            "spatial_data_validation": ["spatial validation", "空间校验", "坐标校验"],
            "spatial_clustering_comparison": ["spatial comparison", "空间域对比", "多分辨率空间"],
            "qc_norm": ["qc", "quality", "标准化", "norm"],
            "dimensionality_reduction": ["pca", "dimensionality", "降维"],
            "clustering": ["clustering", "cluster", "聚类", "leiden"],
            "spatial_neighbors": ["neighbors", "邻域"],
            "spatial_autocorr": ["autocorr", "moran", "svg", "spatially variable"],
            "functional_enrichment": ["pathway", "enrichment", "富集", "通路"],
            "plot_clusters": ["plot cluster", "聚类图"],
            "plot_genes": ["plot gene", "基因图"],
            # Radiomics
            "load_image": ["load", "加载"],
            "preprocess": ["preprocess", "预处理", "resample", "重采样"],
            "preview_slice": ["preview", "预览"],
            "extract_features": ["extract", "feature", "特征提取", "texture"],
            "radiomics_data_validation": ["radiomics validation", "roi 校验", "数据校验"],
            "radiomics_model_comparison": ["roc", "model comparison", "lr svm rf", "诊断效能"],
            "calc_score": ["score", "rad-score", "评分"],
            "viz_score": ["viz", "visualize score", "评分图"],
        }
        
        matched_steps = []
        for step_id, keywords in step_keywords.items():
            if step_id in available_steps:
                if any(kw in query_lower for kw in keywords):
                    matched_steps.append(step_id)
        
        return matched_steps if matched_steps else []
    
    def _resolve_dependencies(
        self,
        target_steps: List[str],
        workflow: "BaseWorkflow"
    ) -> List[str]:
        """
        解析依赖关系（使用硬编码 DAG）
        
        🔥 ARCHITECTURAL REFACTOR: Hardcoded DAG Logic
        
        使用工作流的 DAG 递归解析依赖，确保所有前置步骤都被包含。
        
        Args:
            target_steps: 用户请求的目标步骤列表
            workflow: 工作流实例（包含 steps_dag）
        
        Returns:
            完整的步骤列表（按依赖顺序排序）
        """
        # 使用 BaseWorkflow 的 resolve_dependencies 方法
        resolved_steps = workflow.resolve_dependencies(target_steps)
        return resolved_steps
    
    def _fill_placeholders(
        self,
        workflow_config: Dict[str, Any],
        user_query: str
    ) -> Dict[str, Any]:
        """
        填充占位符（当没有文件时）
        
        🔥 ARCHITECTURAL FIX: Plan-First Interactive Workflow
        
        如果没有文件元数据，使用占位符填充参数，并生成结构化的诊断内容。
        
        Args:
            workflow_config: 工作流配置
            user_query: 用户查询（用于生成诊断信息）
        
        Returns:
            填充占位符后的工作流配置
        """
        steps = workflow_config.get("workflow_data", {}).get("steps", [])
        
        # 🔥 URGENT FIX: 确保 steps 不为空
        if not steps or len(steps) == 0:
            logger.error(f"❌ [SOPPlanner] _fill_placeholders: steps 为空！workflow_config: {workflow_config}")
            # 尝试从 workflow_data 的其他位置获取
            workflow_data = workflow_config.get("workflow_data", {})
            if isinstance(workflow_data, dict):
                # 检查是否有其他字段包含步骤信息
                for key in ["steps", "workflow_steps", "pipeline_steps"]:
                    if key in workflow_data and workflow_data[key]:
                        steps = workflow_data[key]
                        logger.warning(f"⚠️ [SOPPlanner] 从 {key} 获取步骤: {len(steps)} 个")
                        break
            
            # 如果仍然为空，返回错误
            if not steps or len(steps) == 0:
                logger.error(f"❌ [SOPPlanner] _fill_placeholders: 无法获取步骤，返回错误配置")
                return {
                    "type": "error",
                    "error": "工作流步骤为空",
                    "message": "无法生成工作流：步骤列表为空。请检查工作流配置。",
                    "workflow_data": workflow_config.get("workflow_data", {})
                }
        
        logger.info(f"✅ [SOPPlanner] _fill_placeholders: 处理 {len(steps)} 个步骤")
        
        # 🔥 使用中文占位符
        placeholder_text = "<待上传数据>"
        
        for step in steps:
            params = step.get("params", {})
            
            # 填充 file_path 占位符（使用中文）
            if "file_path" in params or "adata_path" in params:
                param_name = "adata_path" if "adata_path" in params else "file_path"
                params[param_name] = placeholder_text
            
            # 填充 group_column 占位符（如果步骤需要）
            step_id = step.get("id") or step.get("tool_id")
            if step_id in ["metabolomics_plsda", "differential_analysis", "metabolomics_pathway_enrichment"]:
                if "group_column" in params:
                    params["group_column"] = "自动检测"
        
        # 🔥 ARCHITECTURAL FIX: 生成结构化的诊断内容（而不是错误）
        # 构建步骤列表文本
        step_names = []
        for step in steps:
            step_name = step.get("name") or step.get("step_name") or step.get("id", "未知步骤")
            step_names.append(step_name)
        
        step_list_text = "\n".join([f"{i+1}. {name}" for i, name in enumerate(step_names)])
        
        guide_message = f"""### 📋 分析方案已生成

根据您的需求 **'{user_query}'**，我为您规划了以下流程：

{step_list_text}

**当前状态**：等待数据。
**下一步**：请点击下方按钮上传文件。"""
        
        workflow_config["diagnosis"] = {
            "status": "template_ready",
            "message": guide_message,
            "steps_count": len(steps),
            "template_mode": True
        }
        
        # 确保 workflow_data 中包含模板标记
        if "workflow_data" in workflow_config:
            workflow_config["workflow_data"]["template_mode"] = True
        
        # 🔥 CRITICAL: 设置顶层 template_mode 标记
        workflow_config["template_mode"] = True
        
        return workflow_config
    
    async def _classify_intent(
        self,
        user_query: str,
        file_metadata: Optional[Dict[str, Any]]
    ) -> Dict[str, Any]:
        """
        使用 LLM 进行意图分类
        
        识别：
        1. 域名（domain_name）："Metabolomics" 或 "RNA"
        2. 目标步骤（target_steps）：用户请求的具体步骤列表
        
        Args:
            user_query: 用户查询
            file_metadata: 文件元数据（可选）
            
        Returns:
            {
                "domain_name": "Metabolomics" | "RNA",
                "target_steps": ["step1", "step2", ...]  # 如果为空，表示完整工作流
            }
        """
        # 🔥 CONTEXT-AWARE INTENT CLASSIFICATION: Determine execution mode based on query + file context
        has_file = file_metadata is not None
        
        # 构建意图分类提示词
        system_prompt = """You are an Intent Classifier for Bioinformatics Workflows.

Your task is to classify user queries into:
1. Domain Name: "Metabolomics", "RNA", "Spatial", or "Radiomics" (strictly one of these four)
2. Mode: "EXECUTION" or "PLANNING" (determines if user wants to run or preview)
3. Target Steps: List of specific steps the user wants (e.g., ["pca_analysis"], ["load_image", "preview_slice", "extract_features"])

**Available Domains:**
- Metabolomics: For metabolite data analysis (CSV files with metabolite measurements)
- RNA: For single-cell RNA-seq analysis (H5AD files, FASTQ files)
- Spatial: For spatial transcriptomics / 10x Visium (Visium directories, spatial coordinates, spots, Moran's I, SVGs)
- Radiomics: For medical imaging (NIfTI, DICOM), CT/MRI, texture analysis, biomarker discovery

**Available Steps (Metabolomics):**
- inspect_data, preprocess_data, pca_analysis, metabolomics_plsda, differential_analysis, visualize_volcano, metabolomics_pathway_enrichment

**Available Steps (RNA):**
- rna_qc_filter, rna_normalize, rna_pca, rna_clustering, rna_find_markers, etc.

**Available Steps (Spatial):**
- load_data, qc_norm, dimensionality_reduction, clustering, spatial_neighbors, spatial_autocorr, plot_clusters, plot_genes

**Available Steps (Radiomics):**
- load_image, preprocess, preview_slice, extract_features, calc_score, viz_score

**Mode Classification Rules (CRITICAL):**
1. IF File is **False** (no file uploaded): ALWAYS return "PLANNING".
2. IF File is **True** (file uploaded):
   - Query implies ACTION ("analyze", "run", "do", "start", "执行", "分析", "运行", "开始"): -> Return "EXECUTION"
   - Query implies INQUIRY ("show me the plan", "what steps?", "preview", "预览", "显示", "查看"): -> Return "PLANNING"
   - Query is VAGUE ("metabolomics", "RNA", "Spatial", "Radiomics", "代谢组", "转录组", "空间", "影像组学"): -> Default to "EXECUTION"

**Output Format:**
Return ONLY a JSON object:
{
  "domain_name": "Metabolomics" | "RNA" | "Spatial" | "Radiomics",
  "mode": "EXECUTION" | "PLANNING",
  "target_steps": ["step1", "step2", ...]  // Empty array [] means full workflow
}

**Few-Shot Example (Spatial):**
- User: "Analyze spatial autocorrelation." -> {"domain_name": "Spatial", "mode": "EXECUTION" or "PLANNING" (depends on file), "target_steps": ["spatial_autocorr", "plot_genes"]}
- User: "Show me the spatial workflow." (no file) -> {"domain_name": "Spatial", "mode": "PLANNING", "target_steps": []}

**Few-Shot Example (Radiomics):**
- User: "Extract features from this CT." -> {"domain_name": "Radiomics", "mode": "EXECUTION", "target_steps": ["load_image", "preprocess", "preview_slice", "extract_features", "calc_score", "viz_score"]}
- User: "Run radiomics on this image." (with file) -> {"domain_name": "Radiomics", "mode": "EXECUTION", "target_steps": []}

**Rules:**
- **User Intent Priority**: If user asks for specific analysis (e.g., "PCA", "spatial autocorrelation", "extract features"), target_steps MUST include only those steps or the full Radiomics pipeline [load_image, preview_slice, extract_features]. Do NOT return full pipeline unless vague.
- If user asks for "full analysis" or "完整分析" or "analyze this file" (vague), use empty array []
- Domain name MUST be exactly "Metabolomics", "RNA", "Spatial", or "Radiomics" (case-sensitive)
- Mode MUST be exactly "EXECUTION" or "PLANNING" (case-sensitive)"""

        # 🔥 TASK 3 FIX: 增强用户查询和文件格式的关联
        query_lower = user_query.lower()
        rna_keywords = ["rna", "scrna", "single cell", "单细胞", "转录组", "cellranger", "cell ranger", "fastq", "测序", "sequencing"]
        metabolomics_keywords = ["metabolomics", "metabolite", "代谢组", "代谢物", "代谢"]
        spatial_keywords = ["visium", "spatial", "slice", "spot", "moran", "spatial transcriptomics", "空间转录组", "空间组学"]
        radiomics_keywords = ["ct", "mri", "radiomics", "texture", "nifti", "dicom", "影像组学", "放射组学", "纹理特征"]

        # 检测用户查询中的领域关键词
        has_rna_keyword = any(kw in query_lower for kw in rna_keywords)
        has_metabolomics_keyword = any(kw in query_lower for kw in metabolomics_keywords)
        has_spatial_keyword = any(kw in query_lower for kw in spatial_keywords)
        has_radiomics_keyword = any(kw in query_lower for kw in radiomics_keywords)

        user_prompt = f"""**User Query:**
{user_query}

**File Context:**
File Uploaded: {has_file} ({'True' if has_file else 'False'})

**File Metadata (if available):**
{json.dumps(file_metadata, ensure_ascii=False, indent=2) if file_metadata else "No file metadata available"}

**CRITICAL ROUTING RULES (User Intent + File Format):**
1. **Visium/Spatial files** (file_type="visium" or domain="Spatial") MUST route to "Spatial" domain.
2. **Medical imaging / Radiomics files** (file_type="medical_image" or domain="Radiomics" or extension .nii/.nii.gz/.dcm) MUST route to "Radiomics" domain.
3. **FASTQ files** (file_type="fastq" or extension=".fastq"/".fq") MUST route to "RNA" domain, regardless of query.
4. **H5AD/10x files** (file_type="h5ad" or "10x_mtx"): If file_type is "visium" or domain is "Spatial" → "Spatial"; else → "RNA".
5. **CSV/Tabular files** (file_type="tabular" or extension=".csv") MUST route to "Metabolomics" domain, unless user explicitly mentions RNA.
6. **User Query Keywords:**
   - If query contains Spatial keywords (visium, spatial, slice, spot, moran): Prefer "Spatial" domain
   - If query contains Radiomics keywords (ct, mri, radiomics, texture, nifti, dicom): Prefer "Radiomics" domain
   - If query contains RNA keywords ({', '.join(rna_keywords[:5])}): Prefer "RNA" domain
   - If query contains Metabolomics keywords ({', '.join(metabolomics_keywords[:3])}): Prefer "Metabolomics" domain
7. **Priority Order:** File format (visium→Spatial, medical_image→Radiomics) > User query keywords > LLM inference

**Task:**
Classify the intent and return JSON only. Remember:
- If File=False: mode MUST be "PLANNING"
- If File=True + Action words: mode = "EXECUTION"
- If File=True + Inquiry words: mode = "PLANNING"
- If File=True + Vague: mode = "EXECUTION" (default)
- **CRITICAL**: File format takes precedence over LLM inference for domain classification."""

        messages = [
            {"role": "system", "content": system_prompt},
            {"role": "user", "content": user_prompt}
        ]
        
        response = await self.llm_client.achat(
            messages=messages,
            temperature=0.1,
            max_tokens=512
        )
        
        # 🔥 FIX: 提取 ChatCompletion 对象的内容
        if hasattr(response, 'choices') and response.choices:
            response_text = response.choices[0].message.content or ""
        else:
            response_text = str(response)
        
        # 解析响应
        try:
            intent_result = json.loads(response_text.strip())
            domain_name = intent_result.get("domain_name", "Metabolomics")
            mode = intent_result.get("mode", "PLANNING")  # 🔥 NEW: Extract mode
            target_steps = intent_result.get("target_steps", [])
            
            # 🔥 TASK 3 FIX: 增强文件格式和用户查询的关联路由（优先级：文件格式 > 用户查询关键词 > LLM推理）
            logger.info(f"🔍 [SOPPlanner] 意图分类后检查: domain_name={domain_name}, file_metadata exists={file_metadata is not None}")
            
            # 检测用户查询中的领域关键词
            query_lower = user_query.lower()
            rna_keywords = ["rna", "scrna", "single cell", "单细胞", "转录组", "cellranger", "cell ranger", "fastq", "测序", "sequencing"]
            metabolomics_keywords = ["metabolomics", "metabolite", "代谢组", "代谢物", "代谢"]
            has_rna_keyword = any(kw in query_lower for kw in rna_keywords)
            has_metabolomics_keyword = any(kw in query_lower for kw in metabolomics_keywords)
            
            if file_metadata:
                file_path = file_metadata.get("file_path", "")
                file_type = file_metadata.get("file_type", "")
                file_domain = file_metadata.get("domain", "")
                logger.info(f"🔍 [SOPPlanner] 文件元数据: file_path={file_path}, file_type={file_type}, domain={file_domain}")
                
                # 🔥 Spatial: Visium / domain=Spatial → Spatial
                if file_type == "visium" or file_domain == "Spatial":
                    if domain_name != "Spatial":
                        logger.warning(f"⚠️ [SOPPlanner] 检测到 Visium/Spatial 文件，强制覆盖域名: {domain_name} → Spatial")
                        domain_name = "Spatial"
                    else:
                        logger.info("✅ [SOPPlanner] Visium/Spatial 文件已正确分类为 Spatial")
                # 🔥 Radiomics: medical_image / domain=Radiomics → Radiomics
                elif file_type == "medical_image" or file_domain == "Radiomics":
                    if domain_name != "Radiomics":
                        logger.warning(f"⚠️ [SOPPlanner] 检测到医学影像/Radiomics 文件，强制覆盖域名: {domain_name} → Radiomics")
                        domain_name = "Radiomics"
                    else:
                        logger.info("✅ [SOPPlanner] 医学影像已正确分类为 Radiomics")
                # 🔥 TASK 3 FIX: 文件格式优先路由（最高优先级）
                # 1. FASTQ文件 → 强制路由到RNA
                elif file_type == "fastq" or (file_path and any(ext in file_path.lower() for ext in [".fastq", ".fq", "fastq"])):
                    if domain_name != "RNA":
                        logger.warning(f"⚠️ [SOPPlanner] 检测到FASTQ文件，强制覆盖域名: {domain_name} → RNA")
                        domain_name = "RNA"
                    else:
                        logger.info(f"✅ [SOPPlanner] FASTQ文件已正确分类为 RNA")
                
                # 2. H5AD/10x文件 → 强制路由到RNA
                elif file_type in ["h5ad", "10x_mtx", "anndata"] or (file_path and ".h5ad" in file_path.lower()):
                    if domain_name != "RNA":
                        logger.warning(f"⚠️ [SOPPlanner] 检测到{file_type}文件，强制覆盖域名: {domain_name} → RNA")
                        domain_name = "RNA"
                    else:
                        logger.info(f"✅ [SOPPlanner] {file_type}文件已正确分类为 RNA")
                
                # 3. CSV/Tabular文件 → 优先路由到Metabolomics（除非用户明确提到RNA）
                elif file_type == "tabular" or (file_path and file_path.lower().endswith(".csv")):
                    if domain_name == "RNA" and not has_rna_keyword:
                        # 如果LLM分类为RNA但用户查询中没有RNA关键词，强制覆盖为Metabolomics
                        logger.warning(f"⚠️ [SOPPlanner] 检测到CSV/Tabular文件且用户未明确提到RNA，强制覆盖域名: {domain_name} → Metabolomics")
                        domain_name = "Metabolomics"
                    elif domain_name == "Metabolomics":
                        logger.info(f"✅ [SOPPlanner] CSV/Tabular文件已正确分类为 Metabolomics")
                
                # 4. 检查文件扩展名（作为补充）
                if file_path:
                    file_ext = file_path.lower().split('.')[-1] if '.' in file_path else ""
                    logger.info(f"🔍 [SOPPlanner] 文件扩展名: {file_ext}")
                    
                    if file_ext == "csv" and domain_name == "RNA" and not has_rna_keyword:
                        logger.warning(f"⚠️ [SOPPlanner] CSV扩展名检测，强制覆盖域名: {domain_name} → Metabolomics")
                        domain_name = "Metabolomics"
                    elif file_ext in ["fastq", "fq"] and domain_name == "Metabolomics":
                        logger.warning(f"⚠️ [SOPPlanner] FASTQ扩展名检测，强制覆盖域名: {domain_name} → RNA")
                        domain_name = "RNA"
                    elif file_ext == "h5ad" and domain_name == "Metabolomics":
                        logger.warning(f"⚠️ [SOPPlanner] H5AD扩展名检测，强制覆盖域名: {domain_name} → RNA")
                        domain_name = "RNA"
            else:
                logger.warning(f"⚠️ [SOPPlanner] 没有文件元数据，使用用户查询关键词进行路由")
                # 如果没有文件，使用用户查询关键词进行路由
                if has_spatial_keyword and not has_rna_keyword and not has_metabolomics_keyword:
                    if domain_name != "Spatial":
                        logger.info(f"ℹ️ [SOPPlanner] 用户查询包含空间组学关键词，调整域名: {domain_name} → Spatial")
                        domain_name = "Spatial"
                elif has_rna_keyword and not has_metabolomics_keyword and not has_spatial_keyword:
                    if domain_name != "RNA":
                        logger.info(f"ℹ️ [SOPPlanner] 用户查询包含RNA关键词，调整域名: {domain_name} → RNA")
                        domain_name = "RNA"
                elif has_metabolomics_keyword and not has_rna_keyword and not has_spatial_keyword and not has_radiomics_keyword:
                    if domain_name != "Metabolomics":
                        logger.info(f"ℹ️ [SOPPlanner] 用户查询包含代谢组关键词，调整域名: {domain_name} → Metabolomics")
                        domain_name = "Metabolomics"
                elif has_radiomics_keyword and not has_spatial_keyword:
                    if domain_name != "Radiomics":
                        logger.info(f"ℹ️ [SOPPlanner] 用户查询包含影像组学关键词，调整域名: {domain_name} → Radiomics")
                        domain_name = "Radiomics"

            # 验证域名（支持 Spatial, Radiomics）
            if domain_name not in ["Metabolomics", "RNA", "Spatial", "Radiomics"]:
                logger.warning(f"⚠️ LLM 返回了无效的域名: {domain_name}，使用默认值 Metabolomics")
                domain_name = "Metabolomics"
            
            # 🔥 CRITICAL: Validate and enforce mode rules
            if not has_file and mode == "EXECUTION":
                logger.warning(f"⚠️ LLM 返回了不一致的模式: 无文件但 mode=EXECUTION，强制设置为 PLANNING")
                mode = "PLANNING"
            
            # 验证模式
            if mode not in ["EXECUTION", "PLANNING"]:
                logger.warning(f"⚠️ LLM 返回了无效的模式: {mode}，使用默认值 PLANNING")
                mode = "PLANNING"
            
            logger.info(f"✅ [SOPPlanner] 意图分类结果: domain={domain_name}, mode={mode}, target_steps={len(target_steps)}")
            
            return {
                "domain_name": domain_name,
                "mode": mode,  # 🔥 NEW: Include mode in result
                "target_steps": target_steps if isinstance(target_steps, list) else []
            }
        except json.JSONDecodeError as e:
            logger.error(f"❌ 意图分类 JSON 解析失败: {e}")
            logger.error(f"响应内容: {response_text[:200] if 'response_text' in locals() else str(response)[:200]}")
            # 回退：尝试从查询和文件元数据中推断
            return self._fallback_intent_classification(user_query, has_file, file_metadata)
    
    def _fallback_intent_classification(self, user_query: str, has_file: bool = False, file_metadata: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
        """
        回退的意图分类（基于文件内容和关键词）
        
        🔥 STRATEGIC FIX: Content-Aware Routing
        - 优先使用 FileInspector 元数据（列名、数据类型）进行智能路由
        - 回退到文件扩展名
        - 最后使用关键词匹配
        
        Args:
            user_query: 用户查询
            has_file: 是否有文件（用于决定模式）
            file_metadata: 文件元数据（可选）
            
        Returns:
            意图分类结果（包含 mode）
        """
        query_lower = user_query.lower()
        domain_name = "Metabolomics"  # 默认值
        spatial_keywords_fallback = ["visium", "spatial", "slice", "spot", "moran", "空间转录组", "空间组学"]

        # 🔥 STRATEGIC FIX: Content-Aware Routing - 基于文件内容
        if file_metadata:
            file_path = file_metadata.get("file_path", "")
            file_type = file_metadata.get("file_type", "")
            file_domain = file_metadata.get("domain", "")
            columns = file_metadata.get("columns", [])
            # Visium/Spatial 优先
            if file_type == "visium" or file_domain == "Spatial":
                domain_name = "Spatial"
                logger.info("✅ [SOPPlanner] Fallback: 检测到 Visium/Spatial 文件，使用 Spatial 域名")
            elif file_type == "medical_image" or file_domain == "Radiomics":
                domain_name = "Radiomics"
                logger.info("✅ [SOPPlanner] Fallback: 检测到医学影像/Radiomics 文件，使用 Radiomics 域名")
            elif file_type == "anndata" or (file_path and file_path.lower().endswith(('.h5ad', '.h5', '.loom'))):
                logger.info("✅ [SOPPlanner] 检测到 RNA 文件类型（anndata/h5ad），使用 RNA 域名")
                domain_name = "RNA"
            elif file_type == "tabular" or (file_path and file_path.lower().endswith('.csv')):
                feature_columns = file_metadata.get("feature_columns", [])
                metadata_columns = file_metadata.get("metadata_columns", [])
                # 策略2: 检查列名模式（内容感知）
                # RNA 特征：包含 "gene", "barcode", "cell", "UMAP", "tSNE", "cluster" 等
                # Metabolomics 特征：包含 "Metabolite", "Compound", "m/z", "RT" 等，或数值列很多
                rna_column_keywords = ["gene", "barcode", "cell", "umap", "tsne", "cluster", "leiden", "pca"]
                metabolomics_column_keywords = ["metabolite", "compound", "m/z", "rt", "retention", "mass"]
                
                column_names_lower = [col.lower() for col in columns] if columns else []
                
                # 检查是否有 RNA 特征列名
                has_rna_keywords = any(kw in ' '.join(column_names_lower) for kw in rna_column_keywords)
                # 检查是否有 Metabolomics 特征列名
                has_metabolomics_keywords = any(kw in ' '.join(column_names_lower) for kw in metabolomics_column_keywords)
                
                # 策略3: 检查数据结构
                # RNA: 通常有大量特征列（基因），少量元数据列
                # Metabolomics: 通常有中等数量的特征列（代谢物），可能有分组列
                is_likely_rna = False
                is_likely_metabolomics = False
                
                if feature_columns and metadata_columns:
                    # 如果特征列数量 > 1000，很可能是 RNA
                    if len(feature_columns) > 1000:
                        is_likely_rna = True
                        logger.info(f"✅ [SOPPlanner] 检测到大量特征列 ({len(feature_columns)})，推断为 RNA")
                    # 如果特征列数量在 10-500 之间，且没有 RNA 关键词，很可能是 Metabolomics
                    elif 10 <= len(feature_columns) <= 500 and not has_rna_keywords:
                        is_likely_metabolomics = True
                        logger.info(f"✅ [SOPPlanner] 检测到中等数量特征列 ({len(feature_columns)})，推断为 Metabolomics")
                
                # 决策逻辑
                if has_rna_keywords or is_likely_rna:
                    domain_name = "RNA"
                    logger.info("✅ [SOPPlanner] 基于列名/数据结构，使用 RNA 域名")
                elif has_metabolomics_keywords or is_likely_metabolomics:
                    domain_name = "Metabolomics"
                    logger.info("✅ [SOPPlanner] 基于列名/数据结构，使用 Metabolomics 域名")
                else:
                    # 回退到文件扩展名
                    if file_path.lower().endswith('.csv'):
                        logger.info("✅ [SOPPlanner] 检测到 CSV 文件，默认使用 Metabolomics 域名（回退）")
                        domain_name = "Metabolomics"
                    else:
                        # 最后使用关键词匹配
                        rna_keywords = ["rna", "scrna", "single cell", "单细胞", "转录组", "cellranger"]
                        if any(kw in query_lower for kw in spatial_keywords_fallback):
                            domain_name = "Spatial"
                        else:
                            domain_name = "RNA" if any(kw in query_lower for kw in rna_keywords) else "Metabolomics"
            else:
                # 未知文件类型，使用关键词匹配
                rna_keywords = ["rna", "scrna", "single cell", "单细胞", "转录组", "cellranger", "h5ad"]
                radiomics_keywords_fallback = ["ct", "mri", "radiomics", "texture", "nifti", "dicom", "影像组学"]
                if any(kw in query_lower for kw in spatial_keywords_fallback):
                    domain_name = "Spatial"
                elif any(kw in query_lower for kw in radiomics_keywords_fallback):
                    domain_name = "Radiomics"
                else:
                    domain_name = "RNA" if any(kw in query_lower for kw in rna_keywords) else "Metabolomics"
        else:
            # No file metadata - check query keywords
            rna_keywords = ["rna", "scrna", "single cell", "单细胞", "转录组", "cellranger", "h5ad"]
            radiomics_keywords_fallback = ["ct", "mri", "radiomics", "texture", "nifti", "dicom", "影像组学"]
            if any(kw in query_lower for kw in spatial_keywords_fallback):
                domain_name = "Spatial"
            elif any(kw in query_lower for kw in radiomics_keywords_fallback):
                domain_name = "Radiomics"
            else:
                domain_name = "RNA" if any(kw in query_lower for kw in rna_keywords) else "Metabolomics"
        
        # 🔥 CRITICAL: Determine mode based on query and file presence
        # Action keywords
        action_keywords = ["analyze", "run", "do", "start", "执行", "分析", "运行", "开始"]
        # Inquiry keywords
        inquiry_keywords = ["show", "preview", "plan", "what", "预览", "显示", "查看", "规划"]
        
        if not has_file:
            mode = "PLANNING"  # No file -> always planning
        elif any(kw in query_lower for kw in action_keywords):
            mode = "EXECUTION"  # Action words + file -> execution
        elif any(kw in query_lower for kw in inquiry_keywords):
            mode = "PLANNING"  # Inquiry words + file -> planning
        else:
            mode = "EXECUTION"  # Vague query + file -> default to execution
        
        return {
            "domain_name": domain_name,
            "mode": mode,
            "target_steps": []
        }
    
    def _fill_parameters(
        self,
        workflow_config: Dict[str, Any],
        file_metadata: Dict[str, Any],
        workflow: "BaseWorkflow",
        template_mode: bool = False  # 🔥 NEW: Allow PLANNING mode with file paths filled
    ) -> Dict[str, Any]:
        """
        填充工作流参数（基于文件元数据）
        
        Args:
            workflow_config: 工作流配置
            file_metadata: 文件元数据
            workflow: 工作流实例
            
        Returns:
            填充参数后的工作流配置
        """
        # 🔥 CRITICAL FIX: Ensure we work on a copy to avoid modifying the original
        workflow_config = workflow_config.copy() if isinstance(workflow_config, dict) else workflow_config
        
        # 🔥 CRITICAL FIX: Ensure workflow_data exists and is a dict
        if "workflow_data" not in workflow_config:
            workflow_config["workflow_data"] = {}
        elif not isinstance(workflow_config["workflow_data"], dict):
            workflow_config["workflow_data"] = {}
        
        steps = workflow_config.get("workflow_data", {}).get("steps", [])
        
        # 🔥 CRITICAL FIX: Ensure steps is a list (not None)
        if not isinstance(steps, list):
            logger.warning(f"⚠️ [SOPPlanner] steps 不是列表类型: {type(steps)}，初始化为空列表")
            steps = []
            workflow_config["workflow_data"]["steps"] = steps
        
        file_path = file_metadata.get("file_path")
        real_data_path = file_metadata.get("real_data_path") or file_path  # Visium root for spatial workflow
        # 🔥 Spatial: 若 real_data_path 为文件路径（如 tissue_positions_list.csv），推导为 Visium 根目录
        if real_data_path:
            rp = Path(real_data_path)
            if rp.is_file():
                parent = rp.parent
                if rp.name == "tissue_positions_list.csv" or "tissue_positions" in rp.name.lower():
                    real_data_path = str(parent.parent if parent.name == "spatial" else parent)
                elif "spatial" in rp.parts:
                    for i, part in enumerate(rp.parts):
                        if part == "spatial" and i > 0:
                            real_data_path = str(Path(*rp.parts[:i]))
                            break
                    else:
                        real_data_path = str(parent)
                else:
                    real_data_path = str(parent)
                logger.info("✅ [SOPPlanner] 从文件路径推导 data_dir -> %s", real_data_path)

        # 🔥 CRITICAL FIX: 检测分组列（用于代谢组学）
        # 首先尝试从 semantic_map 获取
        semantic_map = file_metadata.get("semantic_map", {})
        group_cols = semantic_map.get("group_cols", [])
        
        # 如果没有，使用启发式方法检测
        if not group_cols:
            detected_group_col = self._detect_group_column_heuristic(file_metadata)
            if detected_group_col:
                group_cols = [detected_group_col]
                logger.info(f"✅ [SOPPlanner] 启发式检测到分组列: {detected_group_col}")
                # 更新 semantic_map
                if "semantic_map" not in file_metadata:
                    file_metadata["semantic_map"] = {}
                file_metadata["semantic_map"]["group_cols"] = group_cols
        
        for step in steps:
            params = step.get("params", {})
            step_id = step.get("id")
            
            # 🔥 CRITICAL FIX: 填充 data_dir（Spatial Visium load_data 步骤）
            if "data_dir" in params and real_data_path:
                old_val = params.get("data_dir", "")
                params["data_dir"] = real_data_path
                if old_val and old_val != real_data_path:
                    logger.info(f"✅ [SOPPlanner] 填充 data_dir: {old_val} -> {real_data_path} ({step_id})")
            # 🔥 Spatial Consumer: pass Sensor's detected matrix file name (dynamic H5 loading)
            if "counts_file" in params and file_metadata.get("matrix_file"):
                params["counts_file"] = file_metadata["matrix_file"]
            # 🔥 Radiomics: fill image_path and mask_path from file_metadata
            if "image_path" in params and file_path:
                params["image_path"] = file_path
            if "mask_path" in params and file_metadata.get("mask_path"):
                params["mask_path"] = file_metadata["mask_path"]

            # 🔥 CRITICAL FIX: 填充 file_path 或 adata_path（覆盖占位符）
            if "file_path" in params or "adata_path" in params:
                param_name = "adata_path" if "adata_path" in params else "file_path"
                if file_path:
                    # 🔥 CRITICAL: 覆盖占位符（如 "<待上传数据>"）
                    old_value = params.get(param_name, "")
                    params[param_name] = file_path
                    if old_value and old_value != file_path:
                        logger.info(f"✅ [SOPPlanner] 覆盖占位符: {old_value} -> {file_path}")
                elif params.get(param_name) in ["<待上传数据>", "<PENDING_UPLOAD>", ""]:
                    # 如果没有文件路径但参数是占位符，记录警告
                    logger.warning(f"⚠️ [SOPPlanner] 步骤 {step_id} 的参数 {param_name} 仍然是占位符，但 file_metadata 存在")
            
            # 🔥 CRITICAL FIX: 强制填充 group_column（对于需要分组列的步骤）
            if step_id in ["metabolomics_plsda", "differential_analysis", "metabolomics_pathway_enrichment"]:
                if group_cols:
                    # 强制设置 group_column 参数
                    if "group_column" not in params:
                        params["group_column"] = group_cols[0]
                        logger.info(f"✅ [SOPPlanner] 强制填充 group_column: {group_cols[0]} -> {step_id}")
                    elif params.get("group_column") != group_cols[0]:
                        # 如果已存在但值不同，更新它
                        params["group_column"] = group_cols[0]
                        logger.info(f"✅ [SOPPlanner] 更新 group_column: {params.get('group_column')} -> {group_cols[0]} ({step_id})")
                else:
                    # 如果没有检测到分组列，标记为需要用户输入
                    logger.warning(f"⚠️ [SOPPlanner] 步骤 {step_id} 需要分组列，但未检测到")
                    step["status"] = "waiting_for_upload"
                    step["description"] += " ⚠️ 需要分组信息"
            elif "group_column" in params and group_cols:
                # 对于其他步骤，如果有 group_column 参数且有分组列，填充它
                params["group_column"] = group_cols[0]
        
        # 🔥 CONTEXT-AWARE FIX: Set template_mode based on parameter (allows PLANNING mode with file)
        workflow_config["template_mode"] = template_mode
        
        # 🔥 TASK 2: 清除诊断字段（仅在 EXECUTION 模式）
        # 如果 diagnosis 存在，且是 EXECUTION 模式，则完全移除它
        # Reason: Orchestrator 已经发送了真实的 diagnosis 事件（从 FileInspector）
        # 如果 Planner 返回 diagnosis: null，可能会覆盖 UI
        if not template_mode and "diagnosis" in workflow_config:
            # 🔥 TASK 2: Remove diagnosis key entirely in execution mode
            workflow_config.pop("diagnosis", None)
            logger.info("✅ [SOPPlanner] EXECUTION 模式：已移除 diagnosis 字段，避免覆盖 Orchestrator 的真实诊断")
        
        # 确保 workflow_data 中包含模式标记
        if "workflow_data" in workflow_config:
            workflow_config["workflow_data"]["template_mode"] = template_mode
        
        return workflow_config
    
    def _detect_group_column_heuristic(self, file_metadata: Dict[str, Any]) -> Optional[str]:
        """
        启发式检测分组列
        
        🔥 CRITICAL FIX: 即使列是数值型（int/float），如果唯一值 <= 5，也当作分类变量
        
        Args:
            file_metadata: FileInspector 返回的文件元数据
        
        Returns:
            检测到的分组列名，如果未找到返回 None
        """
        # 优先级关键词列表
        priority_keywords = ['Diet', 'diet', 'Group', 'group', 'Condition', 'condition', 
                            'Treatment', 'treatment', 'Class', 'class', 'Category', 'category',
                            'Type', 'type', 'Label', 'label', 'Status', 'status']
        
        # 方法1: 检查 metadata_columns（FileInspector 可能已经检测到）
        metadata_cols = file_metadata.get("metadata_columns", [])
        if metadata_cols:
            # 优先检查关键词匹配
            for col in metadata_cols:
                if any(keyword in col for keyword in priority_keywords):
                    logger.info(f"✅ [Heuristic] 检测到分组列（关键词匹配）: {col}")
                    return col
            # 如果没有关键词匹配，返回第一个元数据列
            logger.info(f"✅ [Heuristic] 检测到分组列（metadata_columns）: {metadata_cols[0]}")
            return metadata_cols[0]
        
        # 方法2: 检查 potential_groups
        potential_groups = file_metadata.get("potential_groups", {})
        if isinstance(potential_groups, dict) and len(potential_groups) > 0:
            first_group_col = list(potential_groups.keys())[0]
            logger.info(f"✅ [Heuristic] 检测到分组列（potential_groups）: {first_group_col}")
            return first_group_col
        
        # 方法3: 🔥 CRITICAL FIX - 检查数值列，如果唯一值 <= 5，当作分类变量
        columns = file_metadata.get("columns", [])
        head_data = file_metadata.get("head", {})
        
        if columns and head_data:
            try:
                import pandas as pd
                head_json = head_data.get("json", [])
                if isinstance(head_json, list) and len(head_json) > 0:
                    df_preview = pd.DataFrame(head_json)
                    
                    # 首先检查关键词匹配的列
                    for col in columns:
                        if col in df_preview.columns:
                            if any(keyword in col for keyword in priority_keywords):
                                # 检查唯一值数量
                                unique_count = df_preview[col].nunique()
                                if 2 <= unique_count <= 5:  # 🔥 关键：即使数值型，唯一值 <= 5 也当作分类
                                    logger.info(f"✅ [Heuristic] 检测到分组列（数值型关键词匹配）: {col} (唯一值: {unique_count})")
                                    return col
                    
                    # 然后检查所有列（包括数值列）
                    for col in columns:
                        if col in df_preview.columns:
                            unique_count = df_preview[col].nunique()
                            # 🔥 CRITICAL: 即使列是数值型，如果唯一值 <= 5，当作分类变量
                            if 2 <= unique_count <= 5:
                                # 检查列名是否包含分组关键词
                                if any(keyword in col for keyword in priority_keywords):
                                    logger.info(f"✅ [Heuristic] 检测到分组列（数值型，唯一值 <= 5）: {col} (唯一值: {unique_count})")
                                    return col
                                # 或者如果唯一值正好是 2（典型的二元分组）
                                elif unique_count == 2:
                                    logger.info(f"✅ [Heuristic] 检测到分组列（二元数值型）: {col} (唯一值: {unique_count})")
                                    return col
            except Exception as e:
                logger.warning(f"⚠️ [Heuristic] 数值列检测失败: {e}")
        
        # 方法4: 检查所有列名（关键词匹配）
        if columns:
            for col in columns:
                if any(keyword in col for keyword in priority_keywords):
                    logger.info(f"✅ [Heuristic] 检测到分组列（列名关键词匹配）: {col}")
                    return col
        
        logger.info("⚠️ [Heuristic] 未检测到分组列")
        return None
    
    def _generate_metabolomics_plan(self, file_metadata: Optional[Dict[str, Any]]) -> Dict[str, Any]:
        """
        生成确定性的代谢组学 SOP 流程
        
        严格按照线性流程图逻辑：
        1. inspect_data (Always)
        2. preprocess_data (Always)
        3. pca_analysis (Always)
        4. Decision: 如果有分组列
           - metabolomics_plsda
           - differential_analysis
           - visualize_volcano
           - metabolomics_pathway_enrichment
        5. 如果没有分组列：停止
        
        Args:
            file_metadata: FileInspector 返回的文件元数据（None表示模板模式）
        
        Returns:
            符合前端格式的工作流配置字典
        """
        logger.info("📋 [SOPPlanner] 生成确定性代谢组学 SOP 流程")
        
        # 🔥 CRITICAL FIX: 移除文件缺失检查逻辑
        # Planner应该只接受file_metadata作为输入（None或Dict）并相应地输出计划
        # Orchestrator已经处理了分支逻辑，这里不需要再次检查
        
        # 获取文件路径（如果有）
        file_path = file_metadata.get("file_path", "") if file_metadata else ""
        
        # 检测分组列（如果有文件元数据）
        group_column = None
        has_groups = False
        if file_metadata:
            group_column = self._detect_group_column_heuristic(file_metadata)
            has_groups = group_column is not None
        
        logger.info(f"🔍 [SOPPlanner] 分组检测结果: {group_column if has_groups else '无分组列'}")
        
        # 构建步骤列表
        steps = []
        
        # Step 1: Data Inspection (Always)
        steps.append({
            "id": "inspect_data",
            "name": "数据检查",
            "description": "SOP规则：必须首先进行数据质量评估，检查缺失值、数据范围等",
            "selected": True,
            "params": {
                "file_path": file_path if file_path else "<PENDING_UPLOAD>"
            }
        })
        
        # Step 2: Preprocessing (Always)
        steps.append({
            "id": "preprocess_data",
            "name": "数据预处理",
            "description": "SOP规则：必须进行Log2转换和标准化，缺失值处理",
            "selected": True,
            "params": {
                "file_path": file_path if file_path else "<PENDING_UPLOAD>",  # 将自动更新为预处理后的文件
                "log_transform": True,
                "standardize": True,
                "missing_imputation": "min"
            }
        })
        
        # Step 3: Unsupervised Analysis - PCA (Always)
        steps.append({
            "id": "pca_analysis",
            "name": "主成分分析 (PCA)",
            "description": "SOP规则：必须进行PCA分析以探索数据结构和降维",
            "selected": True,
            "params": {
                "file_path": "<preprocess_data_output>",  # 将自动更新
                "n_components": 2,
                "scale": True
            }
        })
        
        # Decision Node: Check for Group Column
        if has_groups:
            logger.info(f"✅ [SOPPlanner] 检测到分组列 '{group_column}'，添加监督分析步骤")
            
            # Step 4: Supervised Analysis - PLS-DA
            steps.append({
                "id": "metabolomics_plsda",
                "name": "PLS-DA 分析",
                "description": f"SOP规则：检测到分组列 '{group_column}'，必须进行监督分析（PLS-DA）以识别组间差异",
                "selected": True,
                "params": {
                    "file_path": "<preprocess_data_output>",
                    "group_column": group_column,
                    "n_components": 2
                }
            })
            
            # Step 5: Differential Analysis
            steps.append({
                "id": "differential_analysis",
                "name": "差异代谢物分析",
                "description": f"SOP规则：必须进行差异分析以识别显著差异的代谢物（分组列: {group_column}）",
                "selected": True,
                "params": {
                    "file_path": "<preprocess_data_output>",
                    "group_column": group_column,
                    "method": "t-test",
                    "p_value_threshold": 0.05,
                    "fold_change_threshold": 1.5
                }
            })
            
            # Step 6: Visualization - Volcano Plot
            # 注意：visualize_volcano 需要 diff_results 字典，包含 results 列表
            steps.append({
                "id": "visualize_volcano",
                "name": "火山图可视化",
                "description": "SOP规则：必须可视化差异分析结果，展示显著差异代谢物",
                "selected": True,
                "params": {
                    "diff_results": "<differential_analysis_output>",  # 将自动从上一个步骤获取
                    "output_path": None,  # 将自动生成
                    "fdr_threshold": 0.05,
                    "log2fc_threshold": 1.0
                }
            })
            
            # Step 7: Functional Analysis - Pathway Enrichment
            # 注意：metabolomics_pathway_enrichment 需要 file_path, group_column, case_group, control_group
            # 需要从数据中自动检测分组值
            potential_groups = file_metadata.get("potential_groups", {}) if file_metadata else {}
            case_group = None
            control_group = None
            
            # 方法1: 从 potential_groups 字典中提取
            if isinstance(potential_groups, dict) and group_column in potential_groups:
                group_values = potential_groups[group_column]
                if isinstance(group_values, list) and len(group_values) >= 2:
                    # 使用前两个分组值
                    case_group = group_values[0]
                    control_group = group_values[1]
                    logger.info(f"✅ [SOPPlanner] 从 potential_groups 检测到分组值: {case_group} vs {control_group}")
            
            # 方法2: 如果 potential_groups 是列表格式（旧格式兼容）
            if (not case_group or not control_group) and isinstance(potential_groups, dict):
                # 尝试从其他键中查找
                for key, values in potential_groups.items():
                    if isinstance(values, list) and len(values) >= 2:
                        case_group = values[0]
                        control_group = values[1]
                        logger.info(f"✅ [SOPPlanner] 从其他分组列检测到分组值: {case_group} vs {control_group}")
                        break
            
            # 方法3: 如果仍然没有检测到，尝试从 head 数据中推断
            if not case_group or not control_group:
                head_data = file_metadata.get("head", {})
                if isinstance(head_data, dict):
                    head_json = head_data.get("json", [])
                    if isinstance(head_json, list) and len(head_json) > 0:
                        # 尝试从第一行数据中获取分组列的值
                        first_row = head_json[0] if isinstance(head_json[0], dict) else {}
                        if group_column in first_row:
                            # 需要读取更多数据来获取唯一值，这里先使用占位符
                            # 实际执行时 differential_analysis 会自动检测
                            logger.info(f"⚠️ [SOPPlanner] 无法从元数据中确定分组值，将在执行时自动检测")
            
            # 如果没有检测到分组值，设置为 None（让工具自动检测）
            # 注意：differential_analysis 会自动检测，但 pathway_enrichment 需要明确的值
            # 所以我们从 differential_analysis 的结果中提取，或者让用户指定
            if not case_group or not control_group:
                # 使用占位符，ExecutionLayer 会尝试从前一步骤的结果中提取
                # 或者工具会在执行时自动检测
                case_group = None  # 让工具自动检测
                control_group = None  # 让工具自动检测
                logger.warning(f"⚠️ [SOPPlanner] 未检测到分组值，将在执行时自动检测或从 differential_analysis 结果中提取")
            
            steps.append({
                "id": "metabolomics_pathway_enrichment",
                "name": "通路富集分析",
                "description": f"SOP规则：必须进行通路富集分析以理解差异代谢物的生物学意义（分组列: {group_column}）",
                "selected": True,
                "params": {
                    "file_path": "<preprocess_data_output>",
                    "group_column": group_column,
                    "case_group": case_group if case_group else "<differential_analysis_case_group>",  # 占位符，ExecutionLayer 会处理
                    "control_group": control_group if control_group else "<differential_analysis_control_group>",  # 占位符
                    "organism": "hsa",  # 默认人类，可根据需要调整
                    "p_value_threshold": 0.05
                }
            })
        else:
            logger.info("⚠️ [SOPPlanner] 未检测到分组列，仅执行无监督分析（PCA）")
            # 添加警告描述到最后一个步骤
            steps[-1]["description"] += " ⚠️ 注意：未检测到分组列，无法进行监督分析和差异分析。"
        
        # 构建工作流配置（符合前端格式）
        workflow_name = "代谢组学标准分析流程" + ("（含分组分析）" if has_groups else "（无监督分析）")
        
        # 适配步骤格式（符合前端格式）
        adapted_steps = []
        for i, step in enumerate(steps, 1):
            tool_id = step["id"]
            
            # 验证工具是否在注册表中
            tool_metadata_obj = registry.get_metadata(tool_id)
            if not tool_metadata_obj:
                logger.warning(f"⚠️ [SOPPlanner] 工具 '{tool_id}' 不在注册表中，跳过")
                continue
            
            # 获取显示名称
            step_display_name = self._get_step_display_name(tool_id, tool_metadata_obj.description)
            
            # 构建适配后的步骤
            adapted_step = {
                "step_id": tool_id,
                "id": tool_id,
                "tool_id": tool_id,
                "name": step.get("name", step_display_name),
                "step_name": step.get("name", step_display_name),
                "description": step.get("description", tool_metadata_obj.description[:100]),
                "desc": step.get("description", tool_metadata_obj.description[:100])[:100],
                "selected": step.get("selected", True),
                "params": step.get("params", {})
            }
            
            adapted_steps.append(adapted_step)
        
        # 构建最终的工作流配置
        workflow_config = {
            "type": "workflow_config",
            "workflow_data": {
                "workflow_name": workflow_name,
                "name": workflow_name,
                "steps": adapted_steps
            },
            "file_paths": [file_path] if file_path else []
        }
        
        logger.info(f"✅ [SOPPlanner] 确定性流程生成完成: {len(adapted_steps)} 个步骤")
        logger.info(f"   步骤列表: {[s['id'] for s in adapted_steps]}")
        
        return workflow_config
    
    def _build_sop_system_prompt(self) -> str:
        """
        构建 SOP 驱动的系统提示词
        
        Returns:
            包含 SOP 规则的系统提示词
        """
        return """You are an expert Bioinformatics Pipeline Architect specializing in Metabolomics data analysis.

Your task is to generate executable workflow plans that STRICTLY follow the Metabolomics Standard Operating Procedure (SOP).

**CRITICAL SOP RULES (MUST FOLLOW):**

0. **Data Validation (MANDATORY FIRST STEP - 秀肌肉):**
   - ALWAYS include metabo_data_validation as the very first step when planning the full workflow. It checks matrix/metadata shape, missing and zero value ratio.
   - Then proceed to inspect_data and preprocess_data.

1. **Data Quality Assessment (AFTER VALIDATION):**
   - IF missing_values > 50% → MUST use a dropping/imputation tool.
   - IF missing_values > 0% AND missing_values <= 50% → MUST use an imputation tool.
   - ALWAYS perform data inspection (inspect_data) after metabo_data_validation.

2. **Data Preprocessing (ALWAYS REQUIRED):**
   - ALWAYS perform Normalization (Log2 transformation + Scaling).
   - Use preprocess_data tool for this step.
   - Parameters should be auto-filled based on file metadata.

3. **Multi-Algorithm Model Comparison (RECOMMENDED - 秀肌肉):**
   - Include metabo_model_comparison after preprocess_data when planning the full pipeline. It runs PCA + PLS-DA + VIP and produces a 1x3 comparison plot. This step is optional for user-skip but MUST be in the default DAG.

4. **Analysis Type Selection (CONDITIONAL):**
   - IF group_columns exist (metadata_columns detected OR numeric columns with ≤5 unique values like 0/1) → MUST perform:
     a. Unsupervised Analysis: PCA (pca_analysis) - ALWAYS perform PCA first
     b. Supervised Analysis: PLS-DA (metabolomics_plsda) - MUST add if groups detected
     c. Differential Analysis (differential_analysis) - MUST add if groups detected
     d. Pathway Enrichment (metabolomics_pathway_enrichment) - MUST be added if differential analysis is planned
   - IF NO group_columns → Perform Unsupervised Analysis only:
     a. PCA (pca_analysis)

5. **Visualization Rules (CRITICAL - MUST FOLLOW):**
   - 🔥 ANTI-REDUNDANCY RULE: The tool `pca_analysis` ALREADY generates a plot. If you select `pca_analysis`, you MUST NOT select `visualize_pca`. They are mutually exclusive. Adding both will cause errors.
   - 🔥 DATA FLOW RULE: If you perform `differential_analysis`, you MUST follow it with `visualize_volcano` to plot the results.
   - 🔥 GROUP DETECTION: If a column (like 'Diet') has few unique values (e.g., 0 and 1, or 2-5 unique values), treat it as a Grouping Column. In this case, you MUST add `metabolomics_plsda` and `metabolomics_pathway_enrichment`.

**OUTPUT CONTRACT (CRITICAL - MUST MATCH EXACTLY):**

You MUST output a JSON object that matches the existing Frontend UI structure EXACTLY:

```json
{
  "name": "Generated Pipeline Name",
  "steps": [
    {
      "id": "tool_name_from_registry",
      "name": "Human Readable Step Name",
      "description": "Why this step is needed (based on SOP rules)",
      "selected": true,
      "params": {
        "file_path": "<auto-filled from metadata>",
        "param1": "value1",
        "param2": "value2"
      }
    }
  ]
}
```

**Step Structure Requirements:**
- "id": MUST be the exact tool name from the retrieved tools (e.g., "pca_analysis", "differential_analysis", "metabolomics_plsda")
- "name": Human-readable name in Chinese (e.g., "主成分分析", "PLS-DA 分析")
- "description": Brief explanation of why this step is needed, referencing SOP rules
- "selected": Always true (all steps are selected by default)
- "params": Dictionary matching the tool's args_schema, with file_path auto-filled from metadata

**Parameter Auto-Filling Rules:**
- file_path: Use the file_path from file_metadata
- group_column: 🔥 CRITICAL - MUST be one of the strings in file_metadata['semantic_map']['group_cols']. If semantic_map['group_cols'] is empty, DO NOT add supervised steps (PLS-DA, Differential Analysis).
- n_components: Default to 2 for PCA/PLS-DA
- scale: Default to true for normalization
- Other parameters: Use sensible defaults from the tool's args_schema

**🔥 SEMANTIC MAPPING CONSTRAINT (CRITICAL):**
- The group_column parameter MUST be selected from semantic_map['group_cols'] list.
- If semantic_map['group_cols'] is empty, SKIP all supervised analysis steps (metabolomics_plsda, differential_analysis, visualize_volcano, metabolomics_pathway_enrichment).
- Do NOT hallucinate or guess group column names. Only use what is explicitly provided in semantic_map['group_cols'].

**Example Output:**
{
  "name": "代谢组学标准分析流程",
  "steps": [
    {
      "id": "inspect_data",
      "name": "数据检查",
      "description": "SOP规则：必须首先进行数据质量评估，检查缺失值、数据范围等",
      "selected": true,
      "params": {
        "file_path": "/app/uploads/data.csv"
      }
    },
    {
      "id": "preprocess_data",
      "name": "数据预处理",
      "description": "SOP规则：必须进行Log2转换和标准化，缺失值率5%需要插补",
      "selected": true,
      "params": {
        "file_path": "/app/uploads/data.csv",
        "log_transform": true,
        "standardize": true,
        "missing_imputation": "min"
      }
    }
  ]
}

Generate ONLY the JSON object, no markdown code blocks, no additional text."""
    
    def _build_sop_user_prompt(
        self,
        user_query: str,
        file_metadata: Dict[str, Any],
        retrieved_tools: List[Dict[str, Any]]
    ) -> str:
        """
        构建 SOP 驱动的用户提示词
        
        Args:
            user_query: 用户查询
            file_metadata: 文件元数据
            retrieved_tools: 检索到的工具列表
        
        Returns:
            用户提示词文本
        """
        # 格式化文件元数据
        metadata_text = self._format_file_metadata(file_metadata)
        
        # 格式化工具信息
        tools_text = []
        for i, tool in enumerate(retrieved_tools, 1):
            tool_info = f"""
Tool {i}: {tool['name']}
  Description: {tool['description']}
  Category: {tool['category']}
  Parameters Schema:
{json.dumps(tool['args_schema'], indent=2, ensure_ascii=False)}
"""
            tools_text.append(tool_info)
        
        # 🔥 ARCHITECTURAL UPGRADE: Phase 2 - Extract semantic_map for fast constraint
        semantic_map = file_metadata.get("semantic_map", {})
        group_cols = semantic_map.get("group_cols", [])
        
        # Fail-Fast: If no group columns, skip supervised steps
        has_groups = len(group_cols) > 0
        
        prompt = f"""**User Query:**
{user_query}

**File Metadata:**
{metadata_text}

**Available Tools:**
{''.join(tools_text)}

**Task (SIMPLIFIED - NO BIOINFORMATICS REASONING):**
Map the user's intent to the available group_cols. If user says 'analyze cachexia', and group_cols contains 'Muscle loss', pick 'Muscle loss'.

**🔥 SEMANTIC MAPPING CONSTRAINT (CRITICAL):**
- semantic_map['group_cols'] = {group_cols}
- The group_column parameter in your JSON MUST be one of these strings: {group_cols if group_cols else '[]'}
- If group_cols is empty ([]), SKIP all supervised steps (metabolomics_plsda, differential_analysis, visualize_volcano, metabolomics_pathway_enrichment)
- Do NOT hallucinate or guess group column names. Only use what is explicitly provided.

**Workflow Decision:**
- Has groups: {has_groups}
- If has_groups=True: Add PCA + PLS-DA + Differential + Volcano + Pathway
- If has_groups=False: Add PCA only (unsupervised)

**CRITICAL: Visualization Rules (MUST FOLLOW STRICTLY)**
- 🔥 ANTI-REDUNDANCY RULE: The tool `pca_analysis` ALREADY generates a plot. If you select `pca_analysis`, you MUST NOT select `visualize_pca`. They are mutually exclusive.
- 🔥 DATA FLOW RULE: If you perform `differential_analysis`, you MUST follow it with `visualize_volcano` to plot the results.

**Output:**
Return ONLY a valid JSON object matching the structure:
{{
  "name": "Pipeline Name",
  "steps": [
    {{
      "id": "tool_name",
      "name": "Step Name",
      "description": "Why needed",
      "selected": true,
      "params": {{"file_path": "...", ...}}
    }}
  ]
}}

Remember: Output ONLY the JSON object, no markdown, no code blocks, no explanations."""
        
        return prompt
    
    def _format_file_metadata(self, file_metadata: Dict[str, Any]) -> str:
        """
        格式化文件元数据为可读文本
        
        Args:
            file_metadata: 文件元数据字典
        
        Returns:
            格式化的元数据文本
        """
        if not file_metadata or file_metadata.get("status") != "success":
            return "File metadata not available or invalid."
        
        shape = file_metadata.get("shape", {})
        missing_rate = file_metadata.get("missing_rate", 0)
        metadata_cols = file_metadata.get("metadata_columns", [])
        feature_cols = file_metadata.get("feature_columns", [])
        file_path = file_metadata.get("file_path", "N/A")
        
        # 🔥 TASK 2: Ensure ALL column names are included for LLM context
        all_columns = file_metadata.get("columns", [])
        columns_text = ', '.join(all_columns) if all_columns else 'None'
        
        text = f"""File Path: {file_path}
Shape: {shape.get('rows', 'N/A')} rows × {shape.get('cols', 'N/A')} columns
Missing Rate: {missing_rate}%
**ALL COLUMNS (CRITICAL - Use these exact names):** {columns_text}
Metadata Columns: {', '.join(metadata_cols) if metadata_cols else 'None'}
Feature Columns (first 10): {', '.join(feature_cols[:10]) if feature_cols else 'None'}
Total Features: {file_metadata.get('total_feature_columns', 'N/A')}
"""
        
        # 🔥 ARCHITECTURAL UPGRADE: Phase 2 - Use semantic_map for clear constraints
        semantic_map = file_metadata.get("semantic_map", {})
        group_cols = semantic_map.get("group_cols", [])
        id_col = semantic_map.get("id_col", "N/A")
        feature_count = semantic_map.get("feature_count", "N/A")
        
        text += f"""
**Semantic Map (CRITICAL - USE THIS FOR group_column):**
- ID Column: {id_col}
- Group Columns: {group_cols if group_cols else '[] (NO GROUPS - SKIP SUPERVISED STEPS)'}
- Feature Count: {feature_count}

🔥 CONSTRAINT: The group_column parameter MUST be one of: {group_cols if group_cols else '[]'}
🔥 FAIL-FAST: If group_cols is empty, DO NOT add supervised steps (PLS-DA, Differential, Volcano, Pathway).
"""
        
        # 保留旧格式以兼容
        potential_groups = file_metadata.get("potential_groups", {})
        if isinstance(potential_groups, dict) and len(potential_groups) > 0:
            text += "\n**Grouping Columns (Legacy Format - Use semantic_map instead):**\n"
            for col_name, col_info in potential_groups.items():
                if isinstance(col_info, dict):
                    n_unique = col_info.get("n_unique", "?")
                    values = col_info.get("values", [])
                    values_str = ", ".join([str(v) for v in values[:10]])
                    text += f"  - {col_name}: {n_unique} unique values ({values_str})\n"
        
        # 添加所有列名（用于检测数值型分组列）
        all_columns = file_metadata.get("columns", [])
        if all_columns:
            # 检查是否有数值型列可能包含分组信息
            head_data = file_metadata.get("head", {})
            if head_data and isinstance(head_data, dict):
                head_json = head_data.get("json", [])
                if isinstance(head_json, list) and len(head_json) > 0:
                    try:
                        import pandas as pd
                        df_preview = pd.DataFrame(head_json)
                        text += "\n**Column Analysis (for Group Detection):**\n"
                        priority_keywords = ['Diet', 'diet', 'Group', 'group', 'Condition', 'condition']
                        for col in all_columns:
                            if col in df_preview.columns:
                                unique_count = df_preview[col].nunique()
                                # 如果列名包含分组关键词且唯一值 <= 5，标记为潜在分组列
                                if any(kw in col for kw in priority_keywords) and 2 <= unique_count <= 5:
                                    unique_values = sorted(df_preview[col].unique().tolist())
                                    text += f"  - {col}: {unique_count} unique values {unique_values} ⚠️ POTENTIAL GROUPING COLUMN\n"
                    except Exception as e:
                        logger.debug(f"无法分析列的唯一值: {e}")
        
        # 添加数据范围信息
        data_range = file_metadata.get("data_range", {})
        if data_range:
            text += f"\nData Range: min={data_range.get('min', 'N/A')}, max={data_range.get('max', 'N/A')}\n"
        
        return text
    
    def _parse_llm_response(self, response: str) -> Dict[str, Any]:
        """
        解析 LLM 响应（复用 WorkflowPlanner 的逻辑）
        
        Args:
            response: LLM 返回的文本
        
        Returns:
            解析后的工作流计划字典
        """
        # 尝试提取 JSON（可能被 markdown 代码块包裹）
        response = response.strip()
        
        # 移除 markdown 代码块标记
        if response.startswith("```json"):
            response = response[7:]
        elif response.startswith("```"):
            response = response[3:]
        
        if response.endswith("```"):
            response = response[:-3]
        
        response = response.strip()
        
        try:
            plan = json.loads(response)
            return plan
        except json.JSONDecodeError as e:
            logger.error(f"❌ [SOPPlanner] JSON 解析失败: {e}")
            logger.error(f"响应内容: {response[:500]}")
            
            # 尝试提取 JSON 对象（使用正则表达式）
            import re
            json_match = re.search(r'\{[^{}]*(?:\{[^{}]*\}[^{}]*)*\}', response, re.DOTALL)
            if json_match:
                try:
                    plan = json.loads(json_match.group())
                    logger.info("✅ [SOPPlanner] 使用正则表达式成功提取 JSON")
                    return plan
                except:
                    pass
            
            raise ValueError(f"无法解析 LLM 响应为 JSON: {e}")
    
    def _validate_and_adapt_sop_plan(
        self,
        workflow_plan: Dict[str, Any],
        file_metadata: Dict[str, Any],
        retrieved_tools: List[Dict[str, Any]]
    ) -> Dict[str, Any]:
        """
        验证 SOP 计划并适配为前端格式
        
        Args:
            workflow_plan: LLM 生成的原始计划
            file_metadata: 文件元数据
            retrieved_tools: 检索到的工具列表
        
        Returns:
            验证和适配后的工作流配置（符合前端格式）
        """
        # 验证基本结构
        if "name" not in workflow_plan:
            workflow_plan["name"] = "代谢组学分析流程"
        
        if "steps" not in workflow_plan or not isinstance(workflow_plan["steps"], list):
            raise ValueError("工作流计划必须包含 'steps' 数组")
        
        # 获取文件路径
        file_path = file_metadata.get("file_path", "") if file_metadata else ""
        
        # 🔥 ARCHITECTURAL UPGRADE: Phase 2 - Fail-Fast Logic
        # 检查 semantic_map，如果 group_cols 为空，移除所有监督步骤
        semantic_map = file_metadata.get("semantic_map", {}) if file_metadata else {}
        group_cols = semantic_map.get("group_cols", [])
        has_groups = len(group_cols) > 0
        
        if not has_groups:
            logger.warning("⚠️ [SOPPlanner] Fail-Fast: group_cols 为空，将移除所有监督步骤")
        
        # 创建工具名称到工具信息的映射
        tool_map = {tool['name']: tool for tool in retrieved_tools}
        
        # 🔥 CRITICAL FIX: 硬删除所有 visualize_pca 步骤（pca_analysis 已包含可视化）
        # 不再需要检查，直接移除所有 visualize_pca
        original_steps = workflow_plan.get("steps", [])
        workflow_plan["steps"] = [
            step for step in original_steps
            if step.get("id") != "visualize_pca" and 
               step.get("tool_id") != "visualize_pca" and 
               step.get("tool_name") != "visualize_pca"
        ]
        removed_count = len(original_steps) - len(workflow_plan["steps"])
        if removed_count > 0:
            logger.info(f"✅ [SOPPlanner] 硬删除 {removed_count} 个 visualize_pca 步骤（pca_analysis 已包含可视化）")
        
        # 🔥 Fail-Fast: 如果 group_cols 为空，移除所有监督步骤
        if not has_groups:
            supervised_tools = ["metabolomics_plsda", "differential_analysis", "visualize_volcano", "metabolomics_pathway_enrichment"]
            before_failfast = len(workflow_plan["steps"])
            workflow_plan["steps"] = [
                step for step in workflow_plan["steps"]
                if step.get("id") not in supervised_tools and 
                   step.get("tool_id") not in supervised_tools and 
                   step.get("tool_name") not in supervised_tools
            ]
            removed_supervised = before_failfast - len(workflow_plan["steps"])
            if removed_supervised > 0:
                logger.info(f"✅ [SOPPlanner] Fail-Fast: 移除 {removed_supervised} 个监督步骤（无分组列）")
        
        # 适配步骤格式
        adapted_steps = []
        for i, step in enumerate(workflow_plan["steps"], 1):
            # 验证工具存在性
            tool_id = step.get("id") or step.get("tool_id") or step.get("tool_name")
            if not tool_id:
                logger.warning(f"⚠️ [SOPPlanner] 步骤 {i} 缺少 tool_id，跳过")
                continue
            
            # 🔥 CRITICAL FIX: 硬删除 visualize_pca（pca_analysis 已包含可视化）
            if tool_id == "visualize_pca":
                logger.warning(f"⚠️ [SOPPlanner] 硬删除 visualize_pca 步骤（pca_analysis 已包含可视化）")
                continue
            
            # 检查工具是否在注册表中
            tool_metadata_obj = registry.get_metadata(tool_id)
            if not tool_metadata_obj:
                logger.warning(f"⚠️ [SOPPlanner] 工具 '{tool_id}' 不在注册表中，跳过")
                continue
            
            # 获取工具信息
            tool_info = tool_map.get(tool_id, {})
            tool_desc = tool_metadata_obj.description
            
            # 处理参数
            params = step.get("params", {})
            adapted_params = {}
            
            # 自动填充 file_path
            if "file_path" not in params and file_path:
                adapted_params["file_path"] = file_path
            elif "file_path" in params:
                adapted_params["file_path"] = params["file_path"]
            
            # 🔥 ARCHITECTURAL UPGRADE: Phase 2 - Use semantic_map for group_column
            # 自动填充 group_column（优先使用 semantic_map['group_cols']）
            if "group_column" not in params and file_metadata:
                semantic_map = file_metadata.get("semantic_map", {})
                group_cols = semantic_map.get("group_cols", [])
                if group_cols:
                    # 使用第一个分组列
                    adapted_params["group_column"] = group_cols[0]
                    logger.info(f"✅ [SOPPlanner] 从 semantic_map 自动填充 group_column: {group_cols[0]}")
                else:
                    # 回退到旧逻辑（兼容性）
                    metadata_cols = file_metadata.get("metadata_columns", [])
                    if metadata_cols:
                        adapted_params["group_column"] = metadata_cols[0]
            
            # 🔥 CRITICAL: 验证 group_column 是否在 semantic_map['group_cols'] 中
            if "group_column" in adapted_params and file_metadata:
                semantic_map = file_metadata.get("semantic_map", {})
                group_cols = semantic_map.get("group_cols", [])
                planned_group_col = adapted_params.get("group_column")
                if group_cols and planned_group_col not in group_cols:
                    logger.warning(f"⚠️ [SOPPlanner] group_column '{planned_group_col}' 不在 semantic_map['group_cols'] 中，自动替换为: {group_cols[0]}")
                    adapted_params["group_column"] = group_cols[0]
            
            # 复制其他参数
            for param_name, param_value in params.items():
                if param_name not in ["file_path", "group_column"]:  # 避免重复
                    adapted_params[param_name] = param_value
            
            # 构建适配后的步骤（符合前端格式）
            # 🔥 Frontend Contract: 必须包含 step_id, tool_id, name, step_name, desc
            step_display_name = step.get("name", self._get_step_display_name(tool_id, tool_desc))
            step_description = step.get("description", tool_desc[:100] if tool_desc else "")
            
            adapted_step = {
                "step_id": tool_id,  # 🔥 Frontend requires step_id (not just id)
                "id": tool_id,  # 保留 id 作为兼容字段
                "tool_id": tool_id,  # Frontend requires tool_id
                "name": step_display_name,  # Frontend requires name
                "step_name": step_display_name,  # Frontend requires step_name (compatibility)
                "description": step_description,  # 完整描述
                "desc": step_description[:100] if len(step_description) > 100 else step_description,  # Frontend requires desc (truncated)
                "selected": step.get("selected", True),  # Frontend may use this
                "params": adapted_params
            }
            
            adapted_steps.append(adapted_step)
        
        # 构建最终的工作流配置（符合前端格式）
        workflow_config = {
            "type": "workflow_config",
            "workflow_data": {
                "workflow_name": workflow_plan.get("name", "代谢组学分析流程"),
                "name": workflow_plan.get("name", "代谢组学分析流程"),  # 兼容字段
                "steps": adapted_steps
            },
            "file_paths": [file_path] if file_path else []
        }
        
        return workflow_config
    
    def _get_step_display_name(self, tool_name: str, tool_desc: str) -> str:
        """
        获取步骤的显示名称
        
        Args:
            tool_name: 工具名称
            tool_desc: 工具描述
        
        Returns:
            显示名称（中文）
        """
        # 名称映射（匹配新的工具 ID，含数据校验与多维对比）
        name_mapping = {
            "metabo_data_validation": "正在执行底层数据与稀疏性校验...",
            "metabo_model_comparison": "正在执行多维机器学习模型 (PCA/PLS-DA) 效能对比...",
            "inspect_data": "数据检查",
            "preprocess_data": "数据预处理",
            "pca_analysis": "主成分分析 (PCA)",
            "visualize_pca": "PCA 可视化",
            "differential_analysis": "差异代谢物分析",
            "visualize_volcano": "火山图",
            "metabolomics_plsda": "PLS-DA 分析",
            "metabolomics_pathway_enrichment": "通路富集分析",
            "metabolomics_heatmap": "热图"
        }
        
        if tool_name in name_mapping:
            return name_mapping[tool_name]
        
        # 默认：使用工具名称（美化）
        return tool_name.replace("_", " ").title()


class RNAPlanner(SOPPlanner):
    """
    scRNA-seq 特定的 SOP 规划器
    
    继承自 SOPPlanner，但使用 scRNA-seq 特定的 SOP 规则。
    """
    
    def _build_sop_system_prompt(self) -> str:
        """
        构建 scRNA-seq 特定的 SOP 系统提示词
        
        Returns:
            包含 scRNA-seq SOP 规则的系统提示词
        """
        return """You are an expert Bioinformatics Pipeline Architect specializing in Single-Cell RNA-seq (scRNA-seq) data analysis.

Your task is to generate executable workflow plans that STRICTLY follow the scRNA-seq Standard Operating Procedure (SOP).

**CRITICAL SOP RULES (MUST FOLLOW):**

0. **Data Validation (MANDATORY FIRST STEP - 秀肌肉):**
   - ALWAYS include rna_data_validation as the very first step when planning the full workflow. It performs fast shape/obs/var checks and ensures downstream steps receive valid input.
   - If input is FASTQ → after CellRanger/convert, the first analytical step is still rna_data_validation (or QC); if input is H5AD/matrix → start with rna_data_validation then rna_qc_filter.

1. **Input Type Detection (MANDATORY AFTER DATA VALIDATION):**
   - IF input is FASTQ files (.fastq, .fq) → MUST start with CellRanger (rna_cellranger_count) - This runs ASYNCHRONOUSLY
   - IF input is Matrix/H5AD (.h5ad, .mtx, 10x directory) → Start directly with QC (rna_qc_filter)

2. **Quality Control (ALWAYS REQUIRED AFTER INPUT PROCESSING):**
   - ALWAYS perform QC filtering (rna_qc_filter)
   - Filter cells based on: min_genes (default 200), max_mt (default 20%)
   - ALWAYS perform Doublet Detection (rna_doublet_detection) after QC
   - Generate QC visualization (rna_visualize_qc)

3. **Preprocessing Pipeline (MANDATORY SEQUENCE):**
   - ALWAYS perform Normalization (rna_normalize) - LogNormalize with target_sum=1e4
   - ALWAYS find Highly Variable Genes (rna_hvg) - default n_top_genes=2000
   - ALWAYS scale data (rna_scale) - for PCA preparation

4. **Dimensionality Reduction (REQUIRED):**
   - ALWAYS perform PCA (rna_pca) - default n_comps=50
   - ALWAYS compute Neighbors (rna_neighbors) - default n_neighbors=15
   - ALWAYS perform UMAP (rna_umap) - for visualization
   - OPTIONAL: t-SNE (rna_tsne) - alternative visualization

5. **Clustering (REQUIRED):**
   - ALWAYS perform Leiden Clustering (rna_clustering) - default resolution=0.5
   - Include rna_clustering_comparison after rna_clustering when planning the full pipeline: it runs Leiden at 0.3/0.5/0.8 and produces a 1x3 UMAP comparison (multi-resolution robustness). This step is optional for user-skip but MUST be in the default DAG.
   - Generate clustering visualization (rna_visualize_clustering)

6. **Marker Detection & Annotation (REQUIRED - FINAL):**
   - ALWAYS find Marker Genes (rna_find_markers) - for each cluster
   - ALWAYS perform Cell Type Annotation (rna_cell_annotation) - using markers
   - Generate marker visualization (rna_visualize_markers)
   - Workflow ends after rna_cell_annotation (no export step).

**OUTPUT CONTRACT (CRITICAL - MUST MATCH EXACTLY):**

You MUST output a JSON object that matches the existing Frontend UI structure EXACTLY:

```json
{
  "name": "Generated Pipeline Name",
  "steps": [
    {
      "id": "tool_name_from_registry",
      "name": "Human Readable Step Name",
      "description": "Why this step is needed (based on SOP rules)",
      "selected": true,
      "params": {
        "adata_path": "<auto-filled from metadata>",
        "param1": "value1",
        "param2": "value2"
      }
    }
  ]
}
```

**Step Structure Requirements:**
- "id": MUST be the exact tool name from the retrieved tools (e.g., "rna_qc_filter", "rna_normalize")
- "name": Human-readable name in Chinese (e.g., "质量控制", "数据标准化")
- "description": Brief explanation of why this step is needed, referencing SOP rules
- "selected": Always true (all steps are selected by default)
- "params": Dictionary matching the tool's args_schema

**Parameter Auto-Filling Rules:**
- adata_path: Use the file_path from file_metadata (or output from previous step)
- For CellRanger: fastqs_path, transcriptome_path, output_dir should be auto-filled
- min_genes: Default to 200
- max_mt: Default to 20.0
- n_top_genes: Default to 2000
- resolution: Default to 0.5 for clustering
- Other parameters: Use sensible defaults from the tool's args_schema

**Example Output:**
{
  "name": "单细胞转录组标准分析流程",
  "steps": [
    {
      "id": "rna_qc_filter",
      "name": "质量控制过滤",
      "description": "SOP规则：必须首先进行质量控制，过滤低质量细胞和线粒体基因",
      "selected": true,
      "params": {
        "adata_path": "/app/uploads/data.h5ad",
        "min_genes": 200,
        "max_mt": 20.0
      }
    },
    {
      "id": "rna_normalize",
      "name": "数据标准化",
      "description": "SOP规则：必须进行LogNormalize标准化，为后续分析做准备",
      "selected": true,
      "params": {
        "adata_path": "<step1_output>",
        "target_sum": 10000
      }
    }
  ]
}

Generate ONLY the JSON object, no markdown code blocks, no additional text."""
    
    def _get_step_display_name(self, tool_name: str, tool_desc: str) -> str:
        """
        获取步骤的显示名称（scRNA-seq 特定）
        
        Args:
            tool_name: 工具名称
            tool_desc: 工具描述
        
        Returns:
            显示名称（中文）
        """
        # scRNA-seq 特定的名称映射（含数据校验与多维对比）
        name_mapping = {
            "rna_data_validation": "正在执行底层数据结构与内存预检...",
            "rna_cellranger_count": "Cell Ranger 计数（异步）",
            "rna_convert_cellranger_to_h5ad": "转换为 H5AD 格式",
            "rna_qc_filter": "质量控制过滤",
            "rna_doublet_detection": "双联体检测",
            "rna_normalize": "数据标准化",
            "rna_hvg": "高变基因筛选",
            "rna_scale": "数据缩放",
            "rna_pca": "主成分分析 (PCA)",
            "rna_neighbors": "构建邻接图",
            "rna_umap": "UMAP 降维",
            "rna_tsne": "t-SNE 降维",
            "rna_clustering": "Leiden 聚类",
            "rna_clustering_comparison": "正在生成多分辨率聚类鲁棒性对比图...",
            "rna_find_markers": "Marker 基因检测",
            "rna_cell_annotation": "细胞类型注释",
            "rna_visualize_qc": "QC 可视化",
            "rna_visualize_clustering": "聚类可视化",
            "rna_visualize_markers": "Marker 可视化",
        }
        
        if tool_name in name_mapping:
            return name_mapping[tool_name]
        
        # 默认：使用工具名称（美化）
        return tool_name.replace("_", " ").title()
