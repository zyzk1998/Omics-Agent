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
            "代谢组", "代谢物"
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
        # 方法1：基于关键词的快速路由
        quick_route = self._quick_route(query, uploaded_files)
        if quick_route and quick_route.get("confidence", 0) > 0.8:
            return quick_route
        
        # 方法2：使用 LLM 进行深度分析
        llm_route = await self._llm_route(query, uploaded_files)
        
        # 合并结果
        if llm_route:
            return llm_route
        
        # 默认路由到转录组（向后兼容）
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
        query_lower = query.lower()
        file_paths = self.get_file_paths(uploaded_files or [])
        
        # 检查文件扩展名
        file_types = set()
        for path in file_paths:
            file_type = self.detect_file_type(path)
            if file_type != "unknown":
                file_types.add(file_type)
        
        # 匹配组学类型
        scores = {}
        for modality, keywords in self.MODALITY_KEYWORDS.items():
            score = 0
            # 查询文本匹配
            for keyword in keywords:
                if keyword in query_lower:
                    score += 1
            # 文件类型匹配
            if modality == "transcriptomics" and ("fastq" in file_types or "h5ad" in file_types):
                score += 2
            elif modality == "genomics" and ("bam" in file_types or "vcf" in file_types):
                score += 2
            
            if score > 0:
                scores[modality] = score
        
        if scores:
            best_modality = max(scores.items(), key=lambda x: x[1])
            modality, score = best_modality
            
            return {
                "modality": modality,
                "intent": self._detect_intent(query),
                "confidence": min(0.9, 0.5 + score * 0.1),
                "routing": self.ROUTING_MAP.get(modality, "rna_agent"),
                "reasoning": f"Matched keywords and file types (score: {score})"
            }
        
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
            response = self.llm_client.get_content(completion)
            
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
        except Exception as e:
            print(f"LLM routing failed: {e}")
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

