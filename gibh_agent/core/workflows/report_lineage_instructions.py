"""
专家报告 / 工作流叙事：数据血缘（原始文件名）强制披露片段。

供 BaseAgent._generate_analysis_summary、各专科 Agent 提示词复用。
"""

EXPERT_REPORT_DATA_LINEAGE_RULE_ZH = (
    "在报告的第一段，你必须用一两句话明确指出正在分析的原始数据文件名"
    "（从下方 context_data.primary_source_files 与步骤结果中的路径提取 basename），"
    "以增强报告的针对性和可信度；禁止省略为「该数据」「上传文件」等模糊指代。"
)
