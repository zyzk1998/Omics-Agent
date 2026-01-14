"""
动态工作流规划器 - The Brain

使用 Tool-RAG 架构动态生成可执行的工作流计划。
结合工具检索和 LLM 推理，生成符合前端格式的工作流配置。

包含两个规划器：
1. WorkflowPlanner: 通用工作流规划器
2. SOPPlanner: 基于 SOP（标准操作程序）的领域特定规划器
"""
import json
import logging
from typing import Dict, Any, List, Optional
from pathlib import Path

from .tool_retriever import ToolRetriever
from .tool_registry import registry
from .llm_client import LLMClient

logger = logging.getLogger(__name__)


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
            workflow_plan = self._parse_llm_response(response)
            
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
2. Structure: {"workflow_name": "...", "steps": [{"id": "step1", "tool_name": "...", "params": {...}, "dependency": "..."}]}
3. **Data Flow**: The output of Step N must match the input of Step N+1 (e.g., file paths).
4. If a parameter is a file path, use placeholders like `<step1_output>` or match the uploaded filename.
5. Use tool names EXACTLY as provided in the tool schemas.
6. Only include parameters that exist in the tool's args_schema.
7. For file paths, prefer using the actual uploaded filename if available.

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
  ]
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
            
            adapted_steps.append(adapted_step)
        
        # 构建最终的工作流配置（符合前端格式）
        workflow_config = {
            "type": "workflow_config",
            "workflow_data": {
                "workflow_name": workflow_plan["workflow_name"],
                "steps": adapted_steps
            },
            "file_paths": context_files
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
        # 简单的名称映射
        name_mapping = {
            "metabolomics_pca": "主成分分析 (PCA)",
            "metabolomics_differential_analysis": "差异分析",
            "metabolomics_preprocess": "数据预处理",
            "file_inspect": "文件检查"
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
    SOP 驱动的动态规划器
    
    基于标准操作程序（SOP）规则，使用 LLM 生成符合专业流程的工作流计划。
    专门为代谢组学分析设计，严格遵循 SOP 规则，确保输出符合前端 UI 格式。
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
    
    async def generate_plan(
        self,
        user_query: str,
        file_metadata: Dict[str, Any],
        category_filter: str = "Metabolomics"
    ) -> Dict[str, Any]:
        """
        生成基于 SOP 规则的工作流计划
        
        Args:
            user_query: 用户查询文本
            file_metadata: FileInspector 返回的文件元数据
            category_filter: 工具类别过滤器（默认 "Metabolomics"）
        
        Returns:
            符合前端格式的工作流配置字典
        """
        try:
            logger.info(f"🧠 [SOPPlanner] 开始生成计划: '{user_query}'")
            
            # 🔥 CRITICAL: 使用 LLM 生成工作流计划（保持智能性）
            # 对于所有类型（包括代谢组学），都使用 LLM 规划
            # Step 1: 检索相关工具
            logger.info("🔍 [SOPPlanner] Step 1: 检索相关工具...")
            retrieved_tools = self.tool_retriever.retrieve(
                query=user_query,
                top_k=15,  # 获取更多工具以支持 SOP 规则
                category_filter=category_filter
            )
            
            if not retrieved_tools:
                logger.warning("⚠️ [SOPPlanner] 未检索到相关工具")
                return {
                    "type": "error",
                    "error": "未找到相关工具，请检查查询或工具注册",
                    "message": "无法生成工作流计划"
                }
            
            logger.info(f"✅ [SOPPlanner] 检索到 {len(retrieved_tools)} 个相关工具")
            
            # Step 2: 构建 SOP 驱动的系统提示词
            logger.info("📝 [SOPPlanner] Step 2: 构建 SOP 提示词...")
            system_prompt = self._build_sop_system_prompt()
            user_prompt = self._build_sop_user_prompt(
                user_query=user_query,
                file_metadata=file_metadata,
                retrieved_tools=retrieved_tools
            )
            
            # Step 3: 调用 LLM 生成计划
            logger.info("🤖 [SOPPlanner] Step 3: 调用 LLM 生成计划...")
            messages = [
                {"role": "system", "content": system_prompt},
                {"role": "user", "content": user_prompt}
            ]
            
            response = await self.llm_client.achat(
                messages=messages,
                temperature=0.1,  # 低温度确保遵循 SOP 规则
                max_tokens=2048
            )
            
            # Step 4: 解析 LLM 响应
            logger.info("🔧 [SOPPlanner] Step 4: 解析 LLM 响应...")
            workflow_plan = self._parse_llm_response(response)
            
            # Step 5: 验证和适配为前端格式
            logger.info("✅ [SOPPlanner] Step 5: 验证和适配...")
            validated_plan = self._validate_and_adapt_sop_plan(
                workflow_plan,
                file_metadata,
                retrieved_tools
            )
            
            logger.info(f"✅ [SOPPlanner] 工作流规划完成: {len(validated_plan.get('steps', []))} 个步骤")
            return validated_plan
        
        except Exception as e:
            logger.error(f"❌ [SOPPlanner] 工作流规划失败: {e}", exc_info=True)
            return {
                "type": "error",
                "error": str(e),
                "message": f"工作流规划失败: {str(e)}"
            }
    
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
    
    def _generate_metabolomics_plan(self, file_metadata: Dict[str, Any]) -> Dict[str, Any]:
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
            file_metadata: FileInspector 返回的文件元数据
        
        Returns:
            符合前端格式的工作流配置字典
        """
        logger.info("📋 [SOPPlanner] 生成确定性代谢组学 SOP 流程")
        
        # 获取文件路径
        file_path = file_metadata.get("file_path", "")
        if not file_path:
            raise ValueError("文件元数据中缺少 file_path")
        
        # 检测分组列
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
                "file_path": file_path
            }
        })
        
        # Step 2: Preprocessing (Always)
        steps.append({
            "id": "preprocess_data",
            "name": "数据预处理",
            "description": "SOP规则：必须进行Log2转换和标准化，缺失值处理",
            "selected": True,
            "params": {
                "file_path": file_path,  # 将自动更新为预处理后的文件
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
            potential_groups = file_metadata.get("potential_groups", {})
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

1. **Data Quality Assessment (MANDATORY FIRST STEP):**
   - IF missing_values > 50% → MUST use a dropping/imputation tool.
   - IF missing_values > 0% AND missing_values <= 50% → MUST use an imputation tool.
   - ALWAYS perform data inspection first (inspect_data).

2. **Data Preprocessing (ALWAYS REQUIRED):**
   - ALWAYS perform Normalization (Log2 transformation + Scaling).
   - Use preprocess_data tool for this step.
   - Parameters should be auto-filled based on file metadata.

3. **Analysis Type Selection (CONDITIONAL):**
   - IF group_columns exist (metadata_columns detected OR numeric columns with ≤5 unique values like 0/1) → MUST perform:
     a. Unsupervised Analysis: PCA (pca_analysis) - ALWAYS perform PCA first
     b. Supervised Analysis: PLS-DA (metabolomics_plsda) - MUST add if groups detected
     c. Differential Analysis (differential_analysis) - MUST add if groups detected
     d. Pathway Enrichment (metabolomics_pathway_enrichment) - MUST be added if differential analysis is planned
   - IF NO group_columns → Perform Unsupervised Analysis only:
     a. PCA (pca_analysis)

4. **Visualization Rules (CRITICAL - MUST FOLLOW):**
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
        
        text = f"""File Path: {file_path}
Shape: {shape.get('rows', 'N/A')} rows × {shape.get('cols', 'N/A')} columns
Missing Rate: {missing_rate}%
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
        # 名称映射（匹配新的工具 ID）
        name_mapping = {
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

1. **Input Type Detection (MANDATORY FIRST STEP):**
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
   - Generate clustering visualization (rna_visualize_clustering)

6. **Marker Detection & Annotation (REQUIRED):**
   - ALWAYS find Marker Genes (rna_find_markers) - for each cluster
   - ALWAYS perform Cell Type Annotation (rna_cell_annotation) - using markers
   - Generate marker visualization (rna_visualize_markers)

7. **Export (FINAL STEP):**
   - ALWAYS export results (rna_export_results) - H5AD, CSVs, figures ZIP

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
        # scRNA-seq 特定的名称映射
        name_mapping = {
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
            "rna_find_markers": "Marker 基因检测",
            "rna_cell_annotation": "细胞类型注释",
            "rna_visualize_qc": "QC 可视化",
            "rna_visualize_clustering": "聚类可视化",
            "rna_visualize_markers": "Marker 可视化",
            "rna_export_results": "结果导出"
        }
        
        if tool_name in name_mapping:
            return name_mapping[tool_name]
        
        # 默认：使用工具名称（美化）
        return tool_name.replace("_", " ").title()
