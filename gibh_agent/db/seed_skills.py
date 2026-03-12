"""
7 大核心组学技能种子数据：供 server 启动注入与 API 按需补种使用。
保证 status=approved、author_id=system、main_category=多模态组学。
"""
from sqlalchemy.orm import Session

from gibh_agent.db.models import Skill as SkillModel

CORE_OMICS_SKILLS = [
    {
        "name": "转录组学标准全流程",
        "main_category": "多模态组学",
        "sub_category": "转录组学",
        "description": "一键生成 scRNA-seq 质控、降维聚类、细胞注释、差异分析与通路富集的标准工作流。",
        "prompt_template": """【角色设定】
你是一位拥有 10 年以上经验的资深单细胞生物信息学科学家，精通肿瘤微环境与免疫学，熟练掌握 Scanpy/Seurat 底层逻辑与参数调优。

【任务目标】
接收单细胞表达矩阵，构建标准化的分析工作流，揭示组织中的细胞异质性，并输出高质量图表与专业解读。

【输入规范与前置检查】
在规划工作流前，请确认：1. 物种信息（人/鼠）；2. 数据格式（h5ad/10x）；3. 是否需要调整默认质控阈值（如线粒体比例）。

【标准分析流程】
1. [数据校验]：检查数据格式与内存预分配。
2. [质控过滤]：执行线粒体/核糖体基因过滤，输出 QC 小提琴图。
3. [降维聚类]：执行 PCA 与多分辨率 Leiden 聚类，输出 1x3 UMAP 对比图。
4. [细胞注释]：基于 Marker 基因或预训练模型进行细胞类型鉴定。
5.[差异与富集]：寻找亚群特异性 Marker，输出气泡图(DotPlot)与通路富集结果。

【约束与安全护栏】
1. 绝对禁止幻觉：若 Marker 基因不符合已知细胞类型，请标记为 Unknown，严禁捏造生物学结论。
2. 必须严格调用系统提供的底层分析工具，不可自行编造 Python 脚本输出。

【输出要求】
请输出包含 QC 统计表、多分辨率 UMAP 对比图、核心 Marker 气泡图的完整数据诊断与专家解读报告。""",
    },
    {
        "name": "空间域识别与聚类",
        "main_category": "多模态组学",
        "sub_category": "空间组学",
        "description": "聚焦空间域识别：坐标/影像校验、空间高变基因、降维聚类、多分辨率空间域物理映射对比。",
        "prompt_template": """【角色设定】
你是一位顶尖的空间转录组学(Spatial Transcriptomics)生信专家，精通 Squidpy 空间数据分析框架，擅长结合基因表达与物理空间坐标进行联合建模。

【任务目标】
接收包含空间坐标的表达矩阵，执行空间域识别与聚类，输出多分辨率空间域物理映射图及生物学解释。

【输入规范与前置检查】
请提示用户上传包含空间坐标的 h5ad 文件，并确认是否包含 H&E 组织切片图像数据 (uns['spatial'])。

【标准分析流程】
1. [空间数据校验]：检查坐标矩阵与 H&E 染色底图的完整性与对齐情况。
2. [空间高变基因]：计算 Moran's I 指数，提取空间特异性基因(SVG)。
3.[空间域聚类]：构建空间邻域图并执行多分辨率聚类。
4. [双屏联动可视化]：输出 UMAP 与真实组织切片的双屏联动对比图。

【约束与安全护栏】
1. 若 H&E 底图缺失，必须安全降级为纯坐标散点图，严禁抛出 KeyError 导致流程中断。
2. 空间域的生物学注释必须结合组织解剖学常识（如肿瘤核心区、浸润边缘区）。

【输出要求】
请输出包含空间高变基因列表、双屏联动聚类图、以及 Top SVG 物理映射热图的深度分析报告。""",
    },
    {
        "name": "差异标志物发现",
        "main_category": "多模态组学",
        "sub_category": "代谢组学",
        "description": "面向代谢物 biomarker：缺失值与归一化、PCA/PLS-DA、多维模型效能对比、VIP 与显著性筛选。",
        "prompt_template": """【角色设定】
你是一位资深的非靶向代谢组学(Untargeted Metabolomics)分析专家，精通高维代谢特征矩阵的缺失值插补、数据归一化及多元统计学分析。

【任务目标】
接收代谢物丰度矩阵与样本分组信息，寻找具有极高临床诊断价值的关键生物标志物(Biomarkers)。

【输入规范与前置检查】
请引导用户上传代谢物丰度矩阵与样本分组 Metadata，并明确询问预期的对比组别（如 Disease vs Control）。

【标准分析流程】
1. [数据预检]：检查特征矩阵的稀疏度与零值比例。
2. [数据预处理]：执行 KNN 插补与 Log2 转换/Auto-scaling。
3. [多维模型对比]：同时执行 PCA 与 PLS-DA，输出 1x3 模型效能对比拼图。
4. [标志物筛选]：提取 VIP > 1 且 P-value < 0.05 的差异代谢物。
5. [高级可视化]：绘制 VIP Score 棒棒糖图与特征聚类热图(Clustermap)。

【约束与安全护栏】
1. 执行有监督学习(PLS-DA)前，必须严格校验样本分组标签的有效性。若分组不足两类，必须主动报错拦截，严禁静默降级。
2. 严禁捏造代谢物的 KEGG 通路信息，必须基于真实数据库检索。

【输出要求】
请输出包含模型对比图、VIP 棒棒糖图、火山图的完整报告，并重点解读 Top 5 代谢物的潜在病理意义。""",
    },
    {
        "name": "多算法诊断模型构建",
        "main_category": "多模态组学",
        "sub_category": "医学影像组学",
        "description": "影像组学建模：ROI/数据校验、特征提取、LR/SVM/RF 多算法 ROC 对比、特征相关性聚类热图。",
        "prompt_template": """【角色设定】
你是一位资深的医学影像组学(Radiomics)与 AI 医疗算法专家，精通 PyRadiomics 特征提取与高维数据降维，擅长构建具有极高临床可解释性的诊断模型。

【任务目标】
接收医学影像特征矩阵，筛选核心特征，构建多算法分类模型，并评估其临床诊断效能。

【输入规范与前置检查】
请提示用户上传影像特征矩阵(CSV)，并明确指出用于分类预测的临床标签列(Label Column)。

【标准分析流程】
1.[影像数据校验]：检查特征矩阵是否存在 NaN/Inf 异常值。
2. [特征降维]：执行方差过滤与 LASSO 回归，筛选核心影像特征。
3. [多算法建模]：基于同一随机种子，同时训练 LR、SVM 和 Random Forest 模型。
4. [效能评估]：在同一坐标系内绘制多算法 ROC 曲线对比图，并标注 AUC 值。
5.[可解释性分析]：输出特征权重棒棒糖图与相关性热图。

【约束与安全护栏】
1. 必须保证训练集与测试集的严格隔离，防止数据泄露(Data Leakage)。
2. 必须提供多算法的横向对比，打破 AI 黑盒，提供临床可解释性。

【输出要求】
请输出包含 LASSO 系数路径图、多算法 ROC 对比图、特征权重图的专业诊断模型评估报告。""",
    },
    {
        "name": "基因组变异检测",
        "main_category": "多模态组学",
        "sub_category": "基因组学",
        "description": "从原始测序数据到变异 calling、注释与致病性评估，输出质控与变异解读报告。",
        "prompt_template": """【角色设定】
你是一位资深的基因组学(Genomics)与临床遗传学专家，精通 WES/WGS 数据分析流程及 GATK 最佳实践。

【任务目标】
接收原始测序数据或 BAM/VCF 文件，执行变异检测与功能注释，评估突变的临床致病性。

【输入规范与前置检查】
请确认用户输入的数据类型（FASTQ/BAM/VCF），以及参考基因组版本（如 hg19/hg38）。

【标准分析流程】
1. [测序质控]：执行 FastQC 质量评估。
2. [序列比对]：使用 BWA-MEM 将 Reads 比对至参考基因组，并进行 PCR 去重。
3. [变异检测]：执行 GATK HaplotypeCaller 进行 SNP/Indel calling。
4. [变异注释]：使用 ANNOVAR 或 SnpEff 进行功能注释与致病性(ClinVar)评估。
5. [报告生成]：输出高频突变基因全景图(OncoPlot)与临床解读报告。

【约束与安全护栏】
1. 变异致病性解读必须严格依据 ACMG 指南，严禁主观臆断突变的临床意义。

【输出要求】
请输出包含突变频谱图(Mutational Signatures)、OncoPlot 及关键致病变异列表的临床级基因组报告。""",
    },
    {
        "name": "表观遗传峰与 motif 分析",
        "main_category": "多模态组学",
        "sub_category": "表观遗传组学",
        "description": "峰 calling、差异分析、motif 富集与转录组/基因组整合解读。",
        "prompt_template": """【角色设定】
你是一位资深的表观遗传学(Epigenomics)专家，精通 ChIP-seq 与 ATAC-seq 分析，擅长解析染色质开放性与转录因子调控网络。

【任务目标】
接收表观遗传测序比对文件，识别富集区域(Peaks)，寻找差异调控位点及核心转录因子结合基序。

【输入规范与前置检查】
请确认实验类型（ChIP-seq 靶点或 ATAC-seq），以及是否包含 Input/Control 对照样本。

【标准分析流程】
1. [数据质控]：评估测序深度、文库复杂度及 TSS 富集得分。
2.[Peak Calling]：使用 MACS2 识别染色质开放区或转录因子结合位点。
3. [差异 Peak 分析]：寻找组间显著变化的调控区域。
4.[Motif 富集]：执行 Homer 分析，寻找核心转录因子结合基序(Motifs)。
5.[多组学整合]：将 Peak 关联至最近的靶基因，并尝试结合转录组数据进行联合解释。

【约束与安全护栏】
1. 必须严格校验 Peak 的假阳性率(FDR)，确保调控位点的可靠性。

【输出要求】
请输出包含 Peak 基因组分布饼图、差异 Peak 火山图、Top Motif 序列 Logo 图的表观调控分析报告。""",
    },
    {
        "name": "蛋白质互作网络分析",
        "main_category": "多模态组学",
        "sub_category": "蛋白质组学",
        "description": "质谱定量、差异蛋白、富集分析与蛋白质互作网络构建与通路解读。",
        "prompt_template": """【角色设定】
你是一位资深的蛋白质组学(Proteomics)与系统生物学专家，精通质谱数据处理与复杂生物网络拓扑分析。

【任务目标】
接收蛋白质定量矩阵，筛选差异表达蛋白，构建并解析蛋白质-蛋白质相互作用(PPI)网络，识别核心枢纽(Hub)蛋白。

【输入规范与前置检查】
请确认定量数据的类型（如 TMT, Label-free），以及是否已进行过 Log 转换和中位数归一化。

【标准分析流程】
1. [数据预处理]：质谱定量数据的缺失值插补与标准化。
2. [差异蛋白筛选]：执行统计学检验(如 limma/t-test)，输出火山图。
3.[PPI 网络构建]：基于 STRING 数据库映射差异蛋白，构建蛋白质互作网络。
4. [Hub 蛋白识别]：使用网络拓扑算法(如 Degree, Betweenness)提取核心枢纽蛋白。
5. [功能富集]：对核心网络模块执行 KEGG/GO 深度富集分析。

【约束与安全护栏】
1. PPI 网络的边权重(Edge Confidence)必须设置合理的阈值（如 >0.4），过滤低质量互作证据。

【输出要求】
请输出包含差异蛋白火山图、高分辨率 PPI 网络拓扑图、Hub 蛋白列表及核心通路富集气泡图的系统生物学报告。""",
    },
]

# 统一的占位符 Prompt（生物医药 UI 占位技能使用）
PLACEHOLDER_PROMPT = "【系统提示】该技能正在建设中，敬请期待..."

# 生物医药大类技能矩阵 (UI 占位)
BIOMEDICINE_SKILLS = [
    {"name": "BepiPred3", "sub_category": "预测与建模", "description": "基于蛋白语言模型的 B 细胞表位预测工具，可高效识别蛋白序列中的潜在线性与构象表位。"},
    {"name": "Evo2", "sub_category": "预测与建模", "description": "生物学基础模型，能够整合长基因组序列的信息，同时保持对单核苷酸变化的敏感性。"},
    {"name": "ESM3", "sub_category": "预测与建模", "description": "模拟蛋白质进化，多模态生成新型蛋白。"},
    {"name": "AlphaFold2", "sub_category": "预测与建模", "description": "根据蛋白质的氨基酸序列预测其三维结构。"},
    {"name": "AlphaFold2-Multimer", "sub_category": "预测与建模", "description": "根据蛋白质的氨基酸序列预测其三维结构，预测精度高，适合需要高准确度的结构研究。"},
    {"name": "DiffSBDD", "sub_category": "预测与建模", "description": "基于扩散模型的结构驱动药物设计工具，可自动生成与蛋白结合口袋高度匹配的小分子配体。"},
    {"name": "OligoFormer", "sub_category": "预测与建模", "description": "可以自动设计并推荐高效靶向特定 mRNA 的 siRNA 分子。"},
    {"name": "MSA search", "sub_category": "数据分析", "description": "通过将查询序列与蛋白质序列数据库进行比对，生成多序列比对（MSA）。"},
    {"name": "ESMFold", "sub_category": "预测与建模", "description": "根据蛋白质的氨基酸序列预测其三维结构，预测速度快，适合快速分析蛋白质。"},
    {"name": "DiffDock", "sub_category": "预测与建模", "description": "预测分子与蛋白质相互作用的 3D 结构。"},
    {"name": "ProtGPT2", "sub_category": "预测与建模", "description": "能够理解蛋白质语言，可用于从头设计和构建蛋白质。"},
    {"name": "RNAFold", "sub_category": "预测与建模", "description": "基于最小自由能原理，预测 RNA 分子的二级结构。"},
    {"name": "RFdiffusion", "sub_category": "预测与建模", "description": "用于蛋白质结合剂设计的蛋白质骨架生成模型。"},
    {"name": "BioGPT", "sub_category": "文本处理", "description": "可用于生物医学命名实体识别、关系提取、文本摘要、对话生成等任务。"},
    {"name": "BioGraph", "sub_category": "数据可视化", "description": "专为基因表达数据和生物信息分析设计的可视化工具，可用于生成基因表达热图、PCA降维图等。"},
    {"name": "基因蛋白信息查询器", "sub_category": "信息检索", "description": "根据基因名查询并返回其在不同物种中的蛋白质注释、序列、亚细胞定位等详细信息。"},
    {"name": "蛋白质资料提取工具", "sub_category": "信息检索", "description": "基于蛋白质ID快速获取其功能注释、亚细胞定位、组织表达及疾病关联信息。"},
    {"name": "蛋白同源结构评估器", "sub_category": "数据分析", "description": "集成 BLAST 序列比对与结构相似性分析，用于精准判定蛋白质同源关系与三维结构差异。"},
    {"name": "蛋白质结构渲染工具", "sub_category": "数据可视化", "description": "支持从PDB格式导入并可视化蛋白质三维结构。"},
    {"name": "RNA二级结构可视化工具", "sub_category": "数据可视化", "description": "基于RNA碱基序列生成其对应的二级结构图示。"},
    {"name": "残基互作分析器", "sub_category": "数据分析", "description": "基于PyMOL分析PDB文件中的残基互作关系，支持复合物结构解析。"},
    {"name": "抗体人源化", "sub_category": "预测与建模", "description": "使用对天然抗体库（Sapiens）或 CDR 移植的深度学习来实现抗体人源化。"},
    {"name": "分析细菌生长曲线", "sub_category": "数据分析", "description": "基于OD600数据拟合细菌生长曲线，提供参数分析、倍增时间和延迟期等结果。"},
    {"name": "UniProt数据库查询", "sub_category": "信息检索", "description": "包含蛋白质序列、功能信息和研究论文索引的蛋白质数据库查询。"},
    {"name": "PubMed数据库查询", "sub_category": "信息检索", "description": "生物医学文献信息检索系统。"},
    {"name": "GWAS catalog 数据库查询", "sub_category": "信息检索", "description": "遗传学研究的重要资源，收录了许多基因组关联数据。"},
    {"name": "dbSNP 数据库查询", "sub_category": "信息检索", "description": "单核苷酸多态性（SNP）数据库，快速检索SNP信息。"},
    {"name": "CRISPR-Cas9 基因编辑工具", "sub_category": "预测与建模", "description": "模拟 CRISPR-Cas9 基因组编辑流程，包含向导RNA验证、目标位点识别等。"},
    {"name": "合成可行性分析工具", "sub_category": "数据分析", "description": "对输入的候选分子 SDF 文件逐一分子计算合成可行性分数（SA Score）。"},
    {"name": "ADMET性质预测工具", "sub_category": "预测与建模", "description": "预测一组化合物的 ADMET 药代动力学属性，包括溶解性、吸收、代谢、毒性等。"},
    {"name": "LigandMPNN", "sub_category": "预测与建模", "description": "深度学习驱动的蛋白质序列设计模型，能够显式考虑小分子、核酸等非蛋白质环境的作用。"},
]

# 自动补全通用字段
for skill in BIOMEDICINE_SKILLS:
    skill["main_category"] = "生物医药"
    skill["prompt_template"] = PLACEHOLDER_PROMPT

# 补全遗漏的生物医药技能 (Phase 2)
ADDITIONAL_BIOMED_SKILLS = [
    {"name": "ESM-Variants", "sub_category": "预测与建模", "description": "交互式可视化蛋白质序列中的氨基酸变化。"},
    {"name": "人源性评估", "sub_category": "数据分析", "description": "使用天然抗体库（OASis）中的肽搜索和种系序列同一性来评估抗体的人性。"},
    {"name": "抗体序列生成", "sub_category": "预测与建模", "description": "对抗体序列进行突变，同时实时监测序列特性。"},
    {"name": "AlphaFold数据库查询", "sub_category": "信息检索", "description": "查询 AlphaFold 数据库并可选下载结构文件。"},
    {"name": "GEO 数据库查询", "sub_category": "信息检索", "description": "快速检索基因表达数据的重要数据库。"},
    {"name": "ClinVar 数据库查询", "sub_category": "信息检索", "description": "快速检索遗传变异的临床相关信息。"},
    {"name": "UCSC 基因组浏览器查询", "sub_category": "信息检索", "description": "基因组学研究中的重要工具，快速检索基因组数据。"},
    {"name": "Reactome 数据库查询", "sub_category": "信息检索", "description": "经过手动筛选和同行评审的生物分子通路知识数据库。"},
    {"name": "InterPro 数据库查询", "sub_category": "信息检索", "description": "整合了多个蛋白质家族、结构域和功能位点的数据库系统。"},
    {"name": "获取mRNA序列工具", "sub_category": "信息检索", "description": "根据目标疾病靶点基因名称检索 mRNA 序列，并返回 FASTA 文件。"},
    {"name": "圆二色谱分析工具", "sub_category": "数据分析", "description": "用于分析圆二色性(CD)光谱数据以确定二级结构和热稳定性。"},
    {"name": "蛋白质序列保守性分析工具", "sub_category": "数据分析", "description": "进行蛋白质多序列比对、系统发育树构建与保守性分析。"},
    {"name": "ITC结合热力学分析工具", "sub_category": "数据分析", "description": "分析等温滴定量热 (ITC) 数据，返回 Kd、ΔH、ΔS 等热力学参数。"},
    {"name": "蛋白酶动力学分析工具", "sub_category": "数据分析", "description": "基于荧光读数拟合 Michaelis-Menten 模型分析蛋白酶动力学数据。"},
    {"name": "酶动力学分析工具", "sub_category": "数据分析", "description": "执行酶动力学实验分析，模拟并拟合酶在不同底物浓度下的反应速率。"},
    {"name": "RNA二级结构分析工具", "sub_category": "数据分析", "description": "分析 RNA 的二级结构特征，包括碱基配对、茎区、环区数量与大小。"},
    {"name": "基因集富集分析工具", "sub_category": "数据分析", "description": "对基因列表执行基因集富集分析。"},
    {"name": "免疫细胞分离与纯化", "sub_category": "临床应用", "description": "模拟免疫细胞的分离与纯化流程。"},
    {"name": "估计细胞周期各阶段持续时间", "sub_category": "数据分析", "description": "基于双核苷脉冲标记的流式细胞术数据，估算细胞周期各阶段时长。"},
    {"name": "查询GSEA支持的数据库工具", "sub_category": "信息检索", "description": "返回 gene set enrichment analysis 支持的数据库名称列表。"},
]
for skill in ADDITIONAL_BIOMED_SKILLS:
    skill["main_category"] = "生物医药"
    skill["prompt_template"] = PLACEHOLDER_PROMPT

# 化学大类技能 (UI 占位，可后续扩充；先有基础列表，再追加海量资产)
CHEMISTRY_SKILLS = []

ADDITIONAL_CHEMISTRY_SKILLS = [
    {"name": "Open Babel", "sub_category": "数据处理", "description": "支持分子文件格式转换、分子对接和虚拟筛选。"},
    {"name": "分子分析工具", "sub_category": "数据分析", "description": "基于RDKit，实现分子性质计算、子结构匹配与相似性分析等功能。"},
    {"name": "分子格式转换工具", "sub_category": "数据处理", "description": "基于RDKit，实现SMILES、SDF、Mol、InChI等化学结构格式的相互转换与标准化。"},
    {"name": "分子可视化工具", "sub_category": "数据可视化", "description": "基于RDKit.js，支持从多种格式读取分子结构并生成2D图像。"},
    {"name": "3D分子结构渲染工具", "sub_category": "数据可视化", "description": "用于渲染通用小分子的三维结构，支持MOL、SDF格式。"},
    {"name": "LAMMPS", "sub_category": "预测与建模", "description": "开源的分子动力学模拟软件，常被用于模拟液体、固体或气态的粒子集合。"},
    {"name": "分子量计算", "sub_category": "数据分析", "description": "根据化学式计算分子质量和元素组成，可用于辅助化学、制药与生物分析研究。"},
    {"name": "药物相似性评估工具", "sub_category": "预测与建模", "description": "根据 Lipinski Rule of Five 检查分子是否具有药物相似性。"},
    {"name": "分子胃肠道吸收能力评估工具", "sub_category": "预测与建模", "description": "估算分子在胃肠道的吸收能力。"},
    {"name": "分子官能团识别工具", "sub_category": "数据分析", "description": "识别分子中的常见官能团 (functional groups)。"},
    {"name": "化学元素查询", "sub_category": "信息检索", "description": "查询指定化学元素的详细信息。"},
    {"name": "分子相似性评估工具", "sub_category": "数据分析", "description": "计算两个分子的 Tanimoto 相似度。"},
    {"name": "分子假阳性片段检测工具", "sub_category": "数据分析", "description": "使用 PAINS filter 检查分子中是否含有可能导致假阳性的片段。"},
    {"name": "分子图像生成工具", "sub_category": "数据可视化", "description": "根据 SMILES 生成分子结构图像 (PNG)，上传到 MinIO 并返回文件 URL。"},
    {"name": "分子Kekulization 转换工具", "sub_category": "数据处理", "description": "对分子进行 Kekulization 转换，将芳香键转为交替单/双键。"},
    {"name": "分子芳香性感知操作工具", "sub_category": "数据处理", "description": "对分子进行芳香性感知操作。"},
    {"name": "分子Pattern Fingerprint生成工具", "sub_category": "数据处理", "description": "生成分子的 Pattern Fingerprint。"},
    {"name": "分子Morgan Fingerprint生成工具", "sub_category": "数据处理", "description": "生成分子的 Morgan Fingerprint（局部化学环境编码的比特向量）。"},
    {"name": "Brenk filter分子毒性检查工具", "sub_category": "预测与建模", "description": "使用 Brenk filter 检查分子是否含有潜在毒性或不良片段。"},
    {"name": "分子BBB评估工具", "sub_category": "预测与建模", "description": "估算分子是否可能穿过血脑屏障 (BBB)。"},
    {"name": "分子量计算工具", "sub_category": "数据分析", "description": "根据 SMILES 表达式计算分子量。"},
    {"name": "Tanimoto 距离矩阵计算工具", "sub_category": "数据分析", "description": "基于分子指纹 (Morgan Fingerprint) 计算 Tanimoto 距离矩阵。"},
    {"name": "CP2K", "sub_category": "预测与建模", "description": "量子化学与固体物理模拟软件包，主要用于原子级别的分子动力学模拟。"},
    {"name": "DRTtools", "sub_category": "数据分析", "description": "基于弛豫时间分布方法的电化学阻抗谱分析工具。"},
]
for skill in ADDITIONAL_CHEMISTRY_SKILLS:
    skill["main_category"] = "化学"
    skill["prompt_template"] = PLACEHOLDER_PROMPT
CHEMISTRY_SKILLS.extend(ADDITIONAL_CHEMISTRY_SKILLS)


def run_seed_core_skills(db: Session) -> int:
    """向当前 Session 插入 7 大核心组学技能（不 commit）。返回插入条数。"""
    for core in CORE_OMICS_SKILLS:
        db.add(SkillModel(
            name=core["name"],
            description=core["description"],
            main_category=core.get("main_category", "多模态组学"),
            sub_category=core["sub_category"],
            prompt_template=core["prompt_template"],
            author_id="system",
            status="approved",
        ))
    return len(CORE_OMICS_SKILLS)


def run_seed_all_system_skills(db: Session) -> int:
    """合并四部分后统一插入：CORE 最前，其次 BIOMEDICINE、ADDITIONAL_BIOMED、CHEMISTRY；所有技能 status=approved、author_id=system。"""
    all_skills = CORE_OMICS_SKILLS + BIOMEDICINE_SKILLS + ADDITIONAL_BIOMED_SKILLS + CHEMISTRY_SKILLS
    for skill in all_skills:
        db.add(SkillModel(
            name=skill["name"],
            description=skill["description"],
            main_category=skill.get("main_category", "多模态组学"),
            sub_category=skill["sub_category"],
            prompt_template=skill["prompt_template"],
            author_id="system",
            status="approved",
        ))
    return len(all_skills)
