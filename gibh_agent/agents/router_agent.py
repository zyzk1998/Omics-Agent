"""
路由智能体
分析用户自然语言输入，识别意图和组学类型，路由到对应的领域智能体
"""
import json
from typing import Dict, Any, List, Optional
from .base_agent import BaseAgent
from ..core.llm_client import LLMClient
from ..core.prompt_manager import PromptManager


class RouterAgent(BaseAgent):
    """
    路由智能体
    
    职责：
    1. 分析用户查询，识别组学类型
    2. 识别用户意图（分析、可视化、解释等）
    3. 路由到对应的领域智能体
    """
    
    # 组学类型关键词映射
    MODALITY_KEYWORDS = {
        "transcriptomics": [
            "rna", "transcript", "expression", "gene", "scRNA", "scrna",
            "single cell", "bulk rna", "rna-seq", "转录组", "单细胞", "基因表达"
        ],
        "genomics": [
            "genome", "wgs", "wes", "variant", "snp", "indel",
            "gatk", "bwa", "基因组", "变异", "突变"
        ],
        "epigenomics": [
            "chip", "atac", "methylation", "epigenetic", "histone",
            "ChIP-seq", "ATAC-seq", "表观遗传"
        ],
        "metabolomics": [
            "metabolite", "metabolism", "lc-ms", "gc-ms", "xcms",
            "代谢组", "代谢物", "代谢组学", "代谢组分析", "代谢分析",
            "metabolomics", "metabolomic", "metabonomics"
        ],
        "proteomics": [
            "protein", "proteome", "mass spec", "maxquant",
            "蛋白质组", "质谱"
        ],
        "spatial_omics": [
            "spatial", "visium", "st", "spatial transcriptomics",
            "空间转录组"
        ],
        "imaging": [
            "image", "microscopy", "histology", "病理", "影像"
        ]
    }
    
    # 路由映射
    ROUTING_MAP = {
        "transcriptomics": "rna_agent",
        "genomics": "dna_agent",
        "epigenomics": "epigenomics_agent",
        "metabolomics": "metabolomics_agent",
        "proteomics": "proteomics_agent",
        "spatial_omics": "spatial_agent",
        "imaging": "imaging_agent"
    }
    
    def __init__(self, llm_client: LLMClient, prompt_manager: PromptManager):
        """初始化路由智能体"""
        super().__init__(llm_client, prompt_manager, "router")
    
    async def process_query(
        self,
        query: str,
        history: List[Dict[str, str]] = None,
        uploaded_files: List[Dict[str, str]] = None,
        **kwargs
    ) -> Dict[str, Any]:
        """
        处理查询并返回路由决策
        
        Returns:
            {
                "modality": "transcriptomics",
                "intent": "single_cell_analysis",
                "confidence": 0.95,
                "routing": "rna_agent",
                "reasoning": "..."
            }
        """
        import logging
        logger = logging.getLogger(__name__)
        
        logger.info(f"🔀 RouterAgent: 开始路由决策 - 查询: {query[:50]}...")
        
        # 方法1：基于关键词的快速路由
        quick_route = self._quick_route(query, uploaded_files)
        # 降低置信度阈值，让更多情况能快速路由（避免 LLM 调用延迟）
        if quick_route and quick_route.get("confidence", 0) > 0.7:  # 从 0.8 降低到 0.7
            logger.info(f"✅ RouterAgent: 快速路由成功 - {quick_route.get('routing')} (confidence: {quick_route.get('confidence', 0):.2f})")
            return quick_route
        
        logger.info(f"⚠️ RouterAgent: 快速路由失败或置信度低，使用 LLM 路由...")
        
        # 方法2：使用 LLM 进行深度分析（带超时保护）
        try:
            import asyncio
            # 设置 10 秒超时，避免 LLM 调用卡住
            llm_route = await asyncio.wait_for(
                self._llm_route(query, uploaded_files),
                timeout=10.0
            )
            if llm_route:
                logger.info(f"✅ RouterAgent: LLM 路由成功 - {llm_route.get('routing')} (confidence: {llm_route.get('confidence', 0):.2f})")
                return llm_route
        except asyncio.TimeoutError:
            logger.warning(f"⏱️ RouterAgent: LLM 路由超时（10秒），使用默认路由")
        except Exception as e:
            logger.error(f"❌ RouterAgent: LLM 路由失败: {e}", exc_info=True)
        
        # 默认路由到转录组（向后兼容）
        logger.warning(f"⚠️ RouterAgent: 使用默认路由 (rna_agent)")
        return {
            "modality": "transcriptomics",
            "intent": "analysis",
            "confidence": 0.5,
            "routing": "rna_agent",
            "reasoning": "Default routing to transcriptomics"
        }
    
    def _quick_route(
        self,
        query: str,
        uploaded_files: List[Dict[str, str]] = None
    ) -> Optional[Dict[str, Any]]:
        """基于关键词的快速路由"""
        import logging
        import os
        logger = logging.getLogger(__name__)
        
        query_lower = query.lower()
        file_paths = self.get_file_paths(uploaded_files or [])
        
        logger.debug(f"🔍 快速路由: 查询='{query_lower[:50]}...', 文件数={len(file_paths)}")
        
        # 🔥 文件优先启发式：在调用 LLM 之前，根据文件扩展名强制路由
        if file_paths:
            # 提取所有文件扩展名
            file_extensions = set()
            for path in file_paths:
                # 处理相对路径和绝对路径
                if os.path.isabs(path):
                    ext = os.path.splitext(path)[1].lower()
                else:
                    ext = os.path.splitext(path)[1].lower()
                if ext:
                    file_extensions.add(ext)
            
            logger.info(f"📁 检测到的文件扩展名: {file_extensions}")
            
            # RNA/转录组文件扩展名
            rna_extensions = {'.h5ad', '.mtx', '.fastq', '.fq', '.bam', '.sam'}
            # 代谢组文件扩展名
            metabolomics_extensions = {'.csv', '.txt', '.xlsx', '.xls', '.tsv'}
            
            # 🔧 修复：优先使用最新上传的文件（列表最后一个）
            # 检查是否有 RNA 文件
            if file_extensions & rna_extensions:
                # 如果同时有 RNA 和代谢组文件，优先使用最新上传的文件（最后一个）
                if file_extensions & metabolomics_extensions:
                    # 检查最后一个文件的扩展名
                    last_file_path = file_paths[-1] if file_paths else ""
                    last_ext = os.path.splitext(last_file_path)[1].lower()
                    if last_ext in metabolomics_extensions:
                        logger.info(f"✅ 文件优先路由: 检测到最新文件是代谢组文件 {last_ext} → metabolomics_agent")
                        return {
                            "modality": "metabolomics",
                            "intent": self._detect_intent(query) if query else "analysis",
                            "confidence": 0.95,
                            "routing": "metabolomics_agent",
                            "reasoning": f"Latest file extension-based routing: {last_ext}"
                        }
                
                logger.info(f"✅ 文件优先路由: 检测到 RNA 文件扩展名 {file_extensions & rna_extensions} → rna_agent")
                return {
                    "modality": "transcriptomics",
                    "intent": self._detect_intent(query) if query else "analysis",
                    "confidence": 0.95,
                    "routing": "rna_agent",
                    "reasoning": f"File extension-based routing: {file_extensions & rna_extensions}"
                }
            
            # 检查是否有代谢组文件
            if file_extensions & metabolomics_extensions:
                logger.info(f"✅ 文件优先路由: 检测到代谢组文件扩展名 {file_extensions & metabolomics_extensions} → metabolomics_agent")
                return {
                    "modality": "metabolomics",
                    "intent": self._detect_intent(query) if query else "analysis",
                    "confidence": 0.95,
                    "routing": "metabolomics_agent",
                    "reasoning": f"File extension-based routing: {file_extensions & metabolomics_extensions}"
                }
        
        # 检查文件扩展名（使用 detect_file_type）
        file_types = set()
        for path in file_paths:
            file_type = self.detect_file_type(path)
            if file_type != "unknown":
                file_types.add(file_type)
        
        logger.debug(f"📁 检测到的文件类型: {file_types}")
        
        # 🔧 修复：优先检查查询中的代谢组关键词（即使没有文件）
        metabolomics_keywords = [
            "代谢组", "代谢物", "代谢组学", "代谢组分析", "代谢分析",
            "metabolite", "metabolism", "metabolomics", "metabolomic", "metabonomics",
            "lc-ms", "gc-ms", "xcms"
        ]
        if any(kw in query_lower for kw in metabolomics_keywords):
            logger.info(f"✅ 快速路由: 查询包含代谢组关键词 → metabolomics_agent")
            return {
                "modality": "metabolomics",
                "intent": self._detect_intent(query),
                "confidence": 0.95,  # 高置信度
                "routing": "metabolomics_agent",
                "reasoning": "Query contains metabolomics keywords"
            }
        
        # 特殊规则：CSV 文件 + 代谢组关键词 = 高置信度快速路由
        if "csv" in file_types:
            # 检查是否有代谢组相关关键词
            if any(kw in query_lower for kw in metabolomics_keywords):
                logger.info(f"✅ 快速路由: CSV 文件 + 代谢组关键词 → metabolomics_agent")
                return {
                    "modality": "metabolomics",
                    "intent": self._detect_intent(query),
                    "confidence": 0.98,  # 极高置信度
                    "routing": "metabolomics_agent",
                    "reasoning": "CSV file with metabolomics keywords"
                }
            # 即使没有明确关键词，CSV 文件也优先考虑代谢组（除非有明确的其他关键词）
            elif not any(kw in query_lower for kw in ["rna", "转录", "单细胞", "scrna", "gene", "表达", "transcript"]):
                logger.info(f"✅ 快速路由: CSV 文件（无其他关键词）→ metabolomics_agent")
                return {
                    "modality": "metabolomics",
                    "intent": self._detect_intent(query),
                    "confidence": 0.85,  # 较高置信度
                    "routing": "metabolomics_agent",
                    "reasoning": "CSV file (likely metabolomics data)"
                }
        
        # 匹配组学类型
        scores = {}
        for modality, keywords in self.MODALITY_KEYWORDS.items():
            score = 0
            # 查询文本匹配
            matched_keywords = []
            for keyword in keywords:
                if keyword in query_lower:
                    score += 1
                    matched_keywords.append(keyword)
            
            # 🔧 修复：提高代谢组关键词的权重
            if modality == "metabolomics" and matched_keywords:
                # 代谢组关键词匹配给予更高权重（每个匹配的关键词额外 +2 分）
                score += len(matched_keywords) * 2
            
            # 文件类型匹配（给予更高权重）
            if modality == "transcriptomics" and ("fastq" in file_types or "h5ad" in file_types):
                score += 3  # 提高权重
            elif modality == "genomics" and ("bam" in file_types or "vcf" in file_types):
                score += 3
            elif modality == "metabolomics" and "csv" in file_types:
                score += 5  # 🔧 修复：CSV 文件强烈暗示代谢组数据，提高权重
            
            if score > 0:
                scores[modality] = score
                logger.debug(f"  {modality}: score={score} (matched: {matched_keywords[:3]})")
        
        if scores:
            best_modality = max(scores.items(), key=lambda x: x[1])
            modality, score = best_modality
            
            # 如果分数足够高，直接返回（不需要 LLM）
            # 降低置信度阈值，让更多情况能快速路由
            confidence = min(0.95, 0.5 + score * 0.15)  # 调整公式，让置信度更容易达到 0.8
            
            result = {
                "modality": modality,
                "intent": self._detect_intent(query),
                "confidence": confidence,
                "routing": self.ROUTING_MAP.get(modality, "rna_agent"),
                "reasoning": f"Matched keywords and file types (score: {score})"
            }
            
            logger.debug(f"✅ 快速路由结果: {result}")
            return result
        
        logger.debug("❌ 快速路由: 无匹配")
        return None
    
    async def _llm_route(
        self,
        query: str,
        uploaded_files: List[Dict[str, str]] = None
    ) -> Optional[Dict[str, Any]]:
        """使用 LLM 进行路由决策"""
        file_paths = self.get_file_paths(uploaded_files or [])
        
        routing_prompt = self.prompt_manager.get_prompt(
            "router",
            {
                "user_query": query,
                "uploaded_files": ", ".join(file_paths) if file_paths else "None"
            },
            fallback=f"""Analyze this query and determine routing:

Query: {query}
Files: {', '.join(file_paths) if file_paths else 'None'}

Return JSON:
{{
    "modality": "transcriptomics|genomics|epigenomics|metabolomics|proteomics|spatial_omics|imaging",
    "intent": "analysis|visualization|interpretation|workflow",
    "confidence": 0.0-1.0,
    "routing": "rna_agent|dna_agent|...",
    "reasoning": "brief explanation"
}}"""
        )
        
        messages = [
            {"role": "system", "content": routing_prompt},
            {"role": "user", "content": query}
        ]
        
        try:
            completion = await self.llm_client.achat(messages, temperature=0.1, max_tokens=256)
            # 提取 think 过程和实际内容
            think_content, response = self.llm_client.extract_think_and_content(completion)
            # 如果有 think 内容，记录日志（可选）
            if think_content:
                import logging
                logger = logging.getLogger(__name__)
                logger.debug(f"Router think process: {think_content[:200]}...")
            
            # 解析 JSON
            json_str = response.strip()
            if "```json" in json_str:
                json_str = json_str.split("```json")[1].split("```")[0].strip()
            elif "```" in json_str:
                json_str = json_str.split("```")[1].split("```")[0].strip()
            
            route_result = json.loads(json_str)
            
            # 确保路由字段存在
            if "routing" not in route_result:
                modality = route_result.get("modality", "transcriptomics")
                route_result["routing"] = self.ROUTING_MAP.get(modality, "rna_agent")
            
            return route_result
        except json.JSONDecodeError as e:
            # JSON 解析失败，使用默认路由
            import logging
            logger = logging.getLogger(__name__)
            logger.warning(f"LLM routing JSON parse failed: {e}, response: {response[:200]}")
            return None
        except Exception as e:
            # 其他错误，使用默认路由
            import logging
            logger = logging.getLogger(__name__)
            logger.error(f"LLM routing failed: {e}", exc_info=True)
            return None
    
    def _detect_intent(self, query: str) -> str:
        """检测用户意图"""
        query_lower = query.lower()
        
        if any(kw in query_lower for kw in ["可视化", "visualize", "plot", "图", "画"]):
            return "visualization"
        elif any(kw in query_lower for kw in ["解释", "interpret", "说明", "含义"]):
            return "interpretation"
        elif any(kw in query_lower for kw in ["流程", "workflow", "pipeline", "分析", "run"]):
            return "workflow"
        else:
            return "analysis"

