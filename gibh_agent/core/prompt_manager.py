"""
提示管理器
使用 Jinja2 模板引擎动态生成提示词
支持专家角色（Expert Personas）模板化
"""
from typing import Dict, Any, Optional
from pathlib import Path
import yaml
from jinja2 import Environment, FileSystemLoader, Template


class PromptManager:
    """
    提示管理器，统一管理所有提示词模板
    
    使用方式：
        manager = PromptManager(template_dir="config/prompts")
        
        # 加载模板
        prompt = manager.get_prompt("rna_expert", {
            "file_path": "/data/sample.h5ad",
            "user_intent": "进行单细胞分析"
        })
    """
    
    def __init__(self, template_dir: str = "config/prompts"):
        """
        初始化提示管理器
        
        Args:
            template_dir: 模板文件目录
        """
        self.template_dir = Path(template_dir)
        self.env = Environment(
            loader=FileSystemLoader(str(self.template_dir)),
            trim_blocks=True,
            lstrip_blocks=True
        )
        self._templates: Dict[str, Template] = {}
        self._load_templates()
    
    def _load_templates(self):
        """加载所有 YAML 模板文件"""
        if not self.template_dir.exists():
            return
        
        for yaml_file in self.template_dir.glob("*.yaml"):
            try:
                with open(yaml_file, 'r', encoding='utf-8') as f:
                    data = yaml.safe_load(f)
                    if data and 'template' in data:
                        template_name = yaml_file.stem
                        self._templates[template_name] = self.env.from_string(data['template'])
            except Exception as e:
                print(f"Warning: Failed to load template {yaml_file}: {e}")
    
    def get_prompt(
        self,
        template_name: str,
        context: Dict[str, Any],
        fallback: Optional[str] = None
    ) -> str:
        """
        获取渲染后的提示词
        
        Args:
            template_name: 模板名称（不含 .yaml 扩展名）
            context: 上下文变量字典
            fallback: 如果模板不存在，使用的备用模板字符串
        
        Returns:
            渲染后的提示词字符串
        """
        if template_name in self._templates:
            return self._templates[template_name].render(**context)
        
        if fallback:
            template = Template(fallback)
            return template.render(**context)
        
        raise ValueError(f"Template '{template_name}' not found and no fallback provided")
    
    def register_template(self, name: str, template_str: str):
        """动态注册模板"""
        self._templates[name] = self.env.from_string(template_str)
    
    def get_system_prompt(
        self,
        expert_role: str,
        context: Dict[str, Any] = None
    ) -> str:
        """
        获取专家角色的系统提示词
        
        Args:
            expert_role: 专家角色名称（如 "rna_expert", "dna_expert"）
            context: 上下文变量
        
        Returns:
            系统提示词
        """
        context = context or {}
        return self.get_prompt(
            f"{expert_role}_system",
            context,
            fallback=f"You are a {expert_role} expert. Please help the user."
        )
    
    def get_user_prompt(
        self,
        template_name: str,
        context: Dict[str, Any]
    ) -> str:
        """
        获取用户提示词
        
        Args:
            template_name: 模板名称
            context: 上下文变量
        
        Returns:
            用户提示词
        """
        return self.get_prompt(template_name, context)


# 预定义的专家角色模板
EXPERT_ROLES = {
    "rna_expert": """You are a Senior Transcriptomics Bioinformatics Expert.

【Your Expertise】
- Single-cell RNA-seq (scRNA-seq) analysis
- Bulk RNA-seq differential expression analysis
- Quality control, normalization, dimensionality reduction
- Cell type annotation, trajectory analysis
- Tools: Cell Ranger, Scanpy, Seurat, DESeq2

【Your Approach】
- Always start with quality control metrics
- Explain each step clearly
- Provide code examples when needed
- Consider batch effects and normalization strategies

【Current Context】
{{ context }}
""",
    
    "dna_expert": """You are a Senior Genomics Bioinformatics Expert.

【Your Expertise】
- Whole Genome Sequencing (WGS)
- Whole Exome Sequencing (WES)
- Variant calling, annotation
- Tools: GATK, BWA, Samtools, VEP

【Your Approach】
- Follow GATK best practices
- Ensure proper quality filtering
- Provide variant annotation and interpretation

【Current Context】
{{ context }}
""",
    
    "router": """You are a Bioinformatics Task Router.

【Your Task】
Analyze user's natural language input and determine:
1. Which omics modality is involved (Transcriptomics, Genomics, Epigenomics, etc.)
2. What is the user's intent (analysis, visualization, interpretation, etc.)
3. Route to the appropriate specialist agent

【Available Modalities】
- Transcriptomics (RNA-seq, scRNA-seq)
- Genomics (WGS, WES)
- Epigenomics (ChIP-seq, ATAC-seq)
- Metabolomics (LC-MS, GC-MS)
- Proteomics (Mass Spec)
- Spatial Omics
- Imaging

【Output Format】
Return JSON:
{
    "modality": "transcriptomics",
    "intent": "single_cell_analysis",
    "confidence": 0.95,
    "routing": "rna_agent"
}

【User Query】
{{ user_query }}

【Uploaded Files】
{{ uploaded_files }}
"""
}


def create_default_prompt_manager() -> PromptManager:
    """创建默认的提示管理器（使用内置模板）"""
    manager = PromptManager()
    
    # 注册内置模板
    for role, template_str in EXPERT_ROLES.items():
        manager.register_template(f"{role}_system", template_str)
    
    return manager

