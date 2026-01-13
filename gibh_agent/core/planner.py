"""
动态工作流规划器 - The Brain

使用 Tool-RAG 架构动态生成可执行的工作流计划。
结合工具检索和 LLM 推理，生成符合前端格式的工作流配置。
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
            adapted_step = {
                "step_id": step.get("id", f"step{i}"),
                "tool_id": tool_name,  # 前端使用 tool_id
                "name": self._get_step_display_name(tool_name, tool_desc),
                "step_name": self._get_step_display_name(tool_name, tool_desc),  # 兼容字段
                "desc": tool_desc[:100] if tool_desc else "",  # 描述截断
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

