"""
Agentic 组件 - 智能查询处理

包含：
1. QueryRewriter: 查询重写器（将模糊查询转换为技术术语）
2. Clarifier: 澄清器（主动询问缺失信息）
3. Reflector: 反思器（自我检查和纠正工作流计划）
"""
import json
import logging
from typing import Dict, Any, List, Optional
from .llm_client import LLMClient

logger = logging.getLogger(__name__)


class QueryRewriter:
    """
    查询重写器
    
    将用户的模糊查询转换为精确的生物信息学术语。
    例如："run that cell thing" -> "Execute scRNA-seq pipeline starting from CellRanger"
    """
    
    def __init__(self, llm_client: LLMClient):
        """
        初始化查询重写器
        
        Args:
            llm_client: LLM 客户端实例
        """
        self.llm_client = llm_client
    
    async def rewrite(self, raw_query: str, history: List[Dict[str, str]] = None) -> str:
        """
        重写查询
        
        Args:
            raw_query: 用户的原始查询
            history: 对话历史（可选，用于上下文）
            
        Returns:
            重写后的查询
        """
        try:
            logger.info(f"🔄 [QueryRewriter] 重写查询: '{raw_query}'")
            
            system_prompt = """You are a Bioinformatics Query Translator.

Your task is to translate vague user queries into precise bioinformatics technical terms.

**Examples:**
- "run that cell thing" -> "Execute scRNA-seq pipeline starting from CellRanger"
- "do pca" -> "Perform Principal Component Analysis (PCA)"
- "analyze this file" -> "Analyze the uploaded file using standard metabolomics pipeline"
- "那个细胞分析" -> "执行单细胞转录组分析流程"

**Rules:**
1. Translate vague terms into specific bioinformatics operations
2. Preserve the user's intent
3. Use standard bioinformatics terminology
4. If the query is already clear, return it as-is
5. Output ONLY the refined query, no explanations

**Output Format:**
Return ONLY the refined query text, no markdown, no code blocks."""

            user_prompt = f"""**User Query:**
{raw_query}

**Context:**
{self._format_history(history) if history else "No conversation history"}

**Task:**
Translate this query into precise bioinformatics technical terms."""

            messages = [
                {"role": "system", "content": system_prompt},
                {"role": "user", "content": user_prompt}
            ]
            
            response = await self.llm_client.achat(
                messages=messages,
                temperature=0.3,
                max_tokens=256
            )
            
            # 🔥 FIX: 提取 ChatCompletion 对象的内容
            if hasattr(response, 'choices') and response.choices:
                refined_query = response.choices[0].message.content or ""
            else:
                refined_query = str(response)
            refined_query = refined_query.strip()
            logger.info(f"✅ [QueryRewriter] 重写完成: '{raw_query}' -> '{refined_query}'")
            
            return refined_query
        
        except Exception as e:
            logger.error(f"❌ [QueryRewriter] 重写失败: {e}", exc_info=True)
            # 失败时返回原始查询
            return raw_query
    
    def _format_history(self, history: List[Dict[str, str]]) -> str:
        """格式化对话历史"""
        if not history:
            return ""
        
        formatted = []
        for item in history[-5:]:  # 只使用最近5条
            role = item.get("role", "user")
            content = item.get("content", "")
            formatted.append(f"{role.capitalize()}: {content}")
        
        return "\n".join(formatted)


class Clarifier:
    """
    澄清器
    
    主动检查缺失的关键信息，并在需要时询问用户。
    例如：代谢组学分析需要分组列，如果缺失则询问用户。
    """
    
    def __init__(self, llm_client: LLMClient):
        """
        初始化澄清器
        
        Args:
            llm_client: LLM 客户端实例
        """
        self.llm_client = llm_client
    
    async def check_and_clarify(
        self,
        refined_query: str,
        file_metadata: Optional[Dict[str, Any]] = None,
        domain: Optional[str] = None
    ) -> Optional[str]:
        """
        检查是否需要澄清
        
        🔥 ARCHITECTURAL MERGE: 修复 Plan-First 逻辑
        - 区分规划和执行意图
        - 如果是规划意图，不要求文件
        
        Args:
            refined_query: 重写后的查询
            file_metadata: 文件元数据（可选）
            domain: 域名（Metabolomics 或 RNA）
            
        Returns:
            如果需要澄清，返回澄清问题；否则返回 None
        """
        try:
            logger.info(f"🔍 [Clarifier] 检查是否需要澄清: '{refined_query}'")
            
            # 🔥 ARCHITECTURAL MERGE: 首先检查用户意图（规划 vs 执行）
            intent = await self._check_intent(refined_query, file_metadata)
            
            if intent == "planning":
                # 用户只是想规划工作流，不需要文件
                logger.info("✅ [Clarifier] 用户意图：规划工作流，无需澄清")
                return None
            
            # 构建上下文信息
            context_info = self._build_context(file_metadata, domain)
            
            system_prompt = """You are a Bioinformatics Clarification Assistant.

Your task is to check if critical information is missing for executing the SOP (Standard Operating Procedure).

**Critical Checks:**
1. **Metabolomics Analysis:**
   - Requires a group column for differential analysis
   - If group column is missing, ask: "请问您想使用哪一列作为分组列？"
   
2. **scRNA-seq Analysis:**
   - Requires input file type (FASTQ vs H5AD)
   - If ambiguous, ask: "请问您的输入数据是 FASTQ 文件还是 H5AD 文件？"

3. **General:**
   - Check if file is uploaded (if file operations are needed)
   - Check if parameters are specified (if tool requires parameters)

**Output Format:**
- If clarification is needed: Return ONLY the clarification question in Chinese
- If no clarification needed: Return "OK"

**Rules:**
- Be polite and concise
- Ask ONE question at a time
- Use Simplified Chinese
- Only ask if information is truly critical"""

            user_prompt = f"""**Refined Query:**
{refined_query}

**Context:**
{context_info}

**Task:**
Check if critical information is missing. If yes, generate a polite clarification question."""

            messages = [
                {"role": "system", "content": system_prompt},
                {"role": "user", "content": user_prompt}
            ]
            
            response = await self.llm_client.achat(
                messages=messages,
                temperature=0.2,
                max_tokens=128
            )
            
            # 🔥 FIX: 提取 ChatCompletion 对象的内容
            if hasattr(response, 'choices') and response.choices:
                clarification = response.choices[0].message.content or ""
            else:
                clarification = str(response)
            clarification = clarification.strip()
            
            # 检查是否需要澄清
            if clarification.upper() == "OK" or clarification.lower() == "ok":
                logger.info("✅ [Clarifier] 无需澄清")
                return None
            else:
                logger.info(f"❓ [Clarifier] 需要澄清: {clarification}")
                return clarification
        
        except Exception as e:
            logger.error(f"❌ [Clarifier] 澄清检查失败: {e}", exc_info=True)
            # 失败时假设不需要澄清
            return None
    
    async def _check_intent(self, query: str, file_metadata: Optional[Dict[str, Any]]) -> str:
        """
        检查用户意图：规划工作流 vs 执行工作流
        
        Returns:
            "planning" 或 "execution"
        """
        try:
            # 如果已经有文件，通常是执行意图
            if file_metadata:
                return "execution"
            
            # 🔥 URGENT FIX: 使用 LLM 判断意图 - 更严格识别 planning
            system_prompt = """You are an Intent Classifier.

Your task is to determine if the user wants to PLAN a workflow or EXECUTE a workflow.

**Planning Intent Keywords (STRICT):**
- "plan", "规划", "设计", "生成工作流", "create workflow", "生成方案"
- "what steps", "哪些步骤", "how to", "怎么做"
- "without file", "先规划", "先看看"
- "I want PCA", "I want to do PCA", "我想做PCA" (when no file)
- "show me", "显示", "预览"
- Any query that asks for a workflow template WITHOUT immediate execution

**Execution Intent Keywords:**
- "run", "执行", "分析", "analyze", "process", "开始分析"
- "do", "做", "开始", "立即"
- "proceed", "继续", "go ahead" (when plan already exists)
- Implies immediate action WITH a file

**CRITICAL RULE:**
- If user says "I want PCA" or "Do PCA" WITHOUT a file, classify as "planning"
- If user says "Plan PCA" or "规划PCA", classify as "planning"
- If user says "Proceed" or "继续" AFTER a plan exists, classify as "execution"

**Output Format:**
Return ONLY one word: "planning" or "execution"."""

            user_prompt = f"""**User Query:**
{query}

**File Available:** {"Yes" if file_metadata else "No"}

**Task:**
Classify the user's intent."""

            messages = [
                {"role": "system", "content": system_prompt},
                {"role": "user", "content": user_prompt}
            ]
            
            response = await self.llm_client.achat(
                messages=messages,
                temperature=0.1,
                max_tokens=16
            )
            
            # 🔥 FIX: 提取 ChatCompletion 对象的内容
            if hasattr(response, 'choices') and response.choices:
                intent = response.choices[0].message.content or ""
            else:
                intent = str(response)
            intent = intent.strip().lower()
            if "plan" in intent:
                return "planning"
            else:
                return "execution"
        
        except Exception as e:
            logger.warning(f"⚠️ [Clarifier] 意图检查失败: {e}，默认返回 execution")
            return "execution"
    
    def _build_context(self, file_metadata: Optional[Dict[str, Any]], domain: Optional[str]) -> str:
        """构建上下文信息"""
        context_parts = []
        
        if domain:
            context_parts.append(f"Domain: {domain}")
        
        if file_metadata:
            context_parts.append("File Metadata Available:")
            
            # 检查分组列
            semantic_map = file_metadata.get("semantic_map", {})
            group_cols = semantic_map.get("group_cols", [])
            if group_cols:
                context_parts.append(f"  - Group columns: {', '.join(group_cols)}")
            else:
                context_parts.append("  - Group columns: None")
            
            # 检查文件类型
            file_type = file_metadata.get("file_type", "unknown")
            context_parts.append(f"  - File type: {file_type}")
        else:
            context_parts.append("File Metadata: Not available")
        
        return "\n".join(context_parts)


class Reflector:
    """
    反思器
    
    自我检查和纠正工作流计划，确保符合 SOP 规则。
    """
    
    def __init__(self, llm_client: LLMClient):
        """
        初始化反思器
        
        Args:
            llm_client: LLM 客户端实例
        """
        self.llm_client = llm_client
    
    async def reflect_and_correct(
        self,
        workflow_plan: Dict[str, Any],
        domain: str,
        file_metadata: Optional[Dict[str, Any]] = None,
        dag_issues: Optional[List[str]] = None
    ) -> Dict[str, Any]:
        """
        反思并纠正工作流计划
        
        Args:
            workflow_plan: 生成的工作流计划
            domain: 域名（Metabolomics 或 RNA）
            file_metadata: 文件元数据（可选）
            
        Returns:
            纠正后的工作流计划（如果发现问题）或原始计划
        """
        try:
            logger.info(f"🔍 [Reflector] 反思工作流计划: {domain}")
            
            # 获取 SOP 规则
            sop_rules = self._get_sop_rules(domain)
            
            # 处理 DAG 问题
            dag_issues_text = ""
            if dag_issues:
                dag_issues_text = "\n".join([f"- {issue}" for issue in dag_issues])
            else:
                dag_issues_text = "None (DAG validation passed)"
            
            system_prompt = """You are a Bioinformatics Workflow Validator.

Your task is to review a generated workflow plan against SOP (Standard Operating Procedure) rules and identify any flaws.

**SOP Rules:**
{sop_rules}

            **Common Flaws to Check:**
1. Missing required steps (e.g., Normalization before PCA)
2. Incorrect step order (e.g., Clustering before Neighbors)
3. Missing prerequisites (e.g., QC before Normalization)
4. Incompatible parameters

**DAG Issues (if provided):**
{dag_issues_text}

**Output Format:**
Return JSON only:
{{
  "is_valid": true/false,
  "issues": ["issue1", "issue2", ...],
  "corrected_steps": [{{"id": "...", "name": "...", ...}}, ...]  // Only if is_valid=false
}}

If the plan is valid, return:
{{
  "is_valid": true,
  "issues": []
}}"""

            user_prompt = f"""**Workflow Plan:**
{json.dumps(workflow_plan, ensure_ascii=False, indent=2)}

**Domain:** {domain}

**File Metadata:**
{json.dumps(file_metadata, ensure_ascii=False, indent=2) if file_metadata else "Not available"}

**Task:**
Review this workflow plan against SOP rules. If flawed, provide corrected steps."""

            messages = [
                {"role": "system", "content": system_prompt.format(sop_rules=sop_rules)},
                {"role": "user", "content": user_prompt}
            ]
            
            response = await self.llm_client.achat(
                messages=messages,
                temperature=0.1,
                max_tokens=1024
            )
            
            # 🔥 FIX: 提取 ChatCompletion 对象的内容
            if hasattr(response, 'choices') and response.choices:
                response_text = response.choices[0].message.content or ""
            else:
                response_text = str(response)
            
            # 解析响应
            try:
                reflection_result = json.loads(response_text.strip())
                
                if reflection_result.get("is_valid", True):
                    logger.info("✅ [Reflector] 工作流计划有效")
                    return workflow_plan
                else:
                    issues = reflection_result.get("issues", [])
                    logger.warning(f"⚠️ [Reflector] 发现 {len(issues)} 个问题: {issues}")
                    
                    # 如果有纠正后的步骤，使用它们
                    corrected_steps = reflection_result.get("corrected_steps")
                    if corrected_steps:
                        # 处理不同的工作流计划格式
                        if "workflow_data" in workflow_plan:
                            workflow_plan["workflow_data"]["steps"] = corrected_steps
                        elif "steps" in workflow_plan:
                            workflow_plan["steps"] = corrected_steps
                        else:
                            # 如果格式不匹配，尝试创建 workflow_data 结构
                            workflow_plan["workflow_data"] = {
                                "workflow_name": workflow_plan.get("name", "Corrected Workflow"),
                                "steps": corrected_steps
                            }
                        logger.info("✅ [Reflector] 已应用纠正")
                    
                    return workflow_plan
            
            except json.JSONDecodeError as e:
                logger.error(f"❌ [Reflector] JSON 解析失败: {e}")
                logger.error(f"响应内容: {response[:200]}")
                # 解析失败时返回原始计划
                return workflow_plan
        
        except Exception as e:
            logger.error(f"❌ [Reflector] 反思失败: {e}", exc_info=True)
            # 失败时返回原始计划
            return workflow_plan
    
    def _get_sop_rules(self, domain: str) -> str:
        """获取 SOP 规则"""
        if domain == "Metabolomics":
            return """**Metabolomics SOP Rules:**
1. MUST start with inspect_data (data quality assessment)
2. MUST perform preprocess_data (Log2 transformation + Scaling)
3. MUST perform pca_analysis (unsupervised analysis)
4. IF group columns exist:
   - MUST perform metabolomics_plsda (supervised analysis)
   - MUST perform differential_analysis
   - MUST perform visualize_volcano
   - MUST perform metabolomics_pathway_enrichment
5. Step order: inspect -> preprocess -> pca -> (if groups) plsda -> diff -> volcano -> pathway"""
        
        elif domain == "RNA":
            return """**scRNA-seq SOP Rules:**
1. IF input is FASTQ: MUST start with rna_cellranger_count
2. IF input is FASTQ: MUST convert with rna_convert_cellranger_to_h5ad
3. MUST perform rna_qc_filter (quality control)
4. MUST perform rna_doublet_detection (after QC)
5. MUST perform rna_normalize (LogNormalize)
6. MUST perform rna_hvg (highly variable genes)
7. MUST perform rna_scale (for PCA)
8. MUST perform rna_pca (before neighbors)
9. MUST perform rna_neighbors (before UMAP/clustering)
10. MUST perform rna_umap (visualization)
11. MUST perform rna_clustering (Leiden)
12. MUST perform rna_find_markers (after clustering)
13. MUST perform rna_cell_annotation (after markers)
14. Step order: QC -> Doublet -> Normalize -> HVG -> Scale -> PCA -> Neighbors -> UMAP/Clustering -> Markers -> Annotation"""
        
        else:
            return f"**{domain} SOP Rules:**\n(No specific rules defined)"

