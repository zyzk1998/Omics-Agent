# -*- coding: utf-8 -*-
"""
六个 RDKit 化学原子技能的 prompt_template（与 seed_skills / patch_chem_skills_batch 同源）。

铁律：首行 [Skill_Route: <tool_id>] 与 gibh_agent.tools.chem_rdkit_tools 注册名完全一致。

原始技能包（2026.04.20 skills 内六个子 ZIP）仅有 SKILL.md + 脚本，无独立 CSV 测试文件；
注入模板与仓库 **`docs/chem_skills_test_inputs/`** 中的 ASCII 示例，均以各子包 **`SKILL.md`** 规定的 **`--smiles` / `--file`（每行一条 SMILES）** 为输入依据，**仅供测试与验收**。
批量或超出会话上下文长度时，改用该目录下的原件或通过会话上传附件。
"""

from __future__ import annotations

# 技能广场卡片「名称」Skill.name -> 注入模板（快车道暗号 + 业务引导 + SMILES 占位）
CHEM_PROMPTS_BY_SKILL_NAME: dict[str, str] = {
    "分子相似性评估工具": """[Skill_Route: chem_tanimoto_similarity]
您好。本工具基于 **Morgan 指纹（半径 2，2048-bit）** 计算 **双分子 Tanimoto 相似度**，适用于先导化合物结构簇对比与排序筛选。

**技能包数据说明**：压缩包内 **无独立「原始数据文件」**，`tanimoto/SKILL.md` 仅以命令行示例给出 SMILES；下表为可从文档复现的一对示例（**可整段替换**）。双行文件形态见 **`docs/chem_skills_test_inputs/pair_smiles_two_lines.txt`**。

| 分子 | 示例 SMILES | 来源说明 |
|------|-------------|----------|
| A（阿司匹林） | CC(=O)Oc1ccccc1C(=O)O | 与分子量文档示例一致 |
| B（咖啡因） | Cn1cnc2c1c(=O)[nH]c(=O)n2C | 经典小分子示例 |

**亦可改用文档中的小矩阵示例**：`CCO` · `c1ccccc1`（乙醇 · 苯），写成一行 `CCO|c1ccccc1` 即可。

**提交方式（任选其一，不要求必须上传文件）**
1. **直接在本提示下方改正文表格中的 SMILES**，再发送；或在与助手对话中用一行写 `SMILES_A|SMILES_B`。
2. 在输入框另写两行 SMILES，由助手填入 **smiles_text** / **smiles2_text**。
3. 上传包含 **两行** SMILES 的 `.csv`/`.txt`（每行一个分子），路径写入 **file_path**。

（助手侧：优先会话附件 **file_path**；否则解析竖线分隔或正文双字段；路径须真实来自附件列表。）
""",
    "分子假阳性片段检测工具": """[Skill_Route: chem_pains_filter]
您好。本工具基于 RDKit **PAINS FilterCatalog（A/B/C）** 检测常见 **泛测定干扰片段**，用于苗头化合物质量把关。

**技能包数据说明**：`pains-filter.zip` 内 **无单独测试 CSV**，文档以占位符示意；下表给出一条可直接跑的示例结构（**可替换**）。单行文件示例见 **`docs/chem_skills_test_inputs/single_smiles_first_line.txt`**。

| 字段 | 示例 SMILES |
|------|-------------|
| 待测分子 | CC(=O)Oc1ccccc1C(=O)O |

**提交方式（任选其一，不要求必须上传文件）**
1. **修改本提示中的示例 SMILES** 后发送；或在对话正文中写出新的单行 SMILES。
2. 上传含 **首行 SMILES** 的 `.csv`/`.txt`，路径写入 **file_path**。

（助手侧：**smiles_text** 与 **file_path** 二选一；路径须来自当前会话附件。）
""",
    "分子Morgan Fingerprint生成工具": """[Skill_Route: chem_morgan_fingerprint]
您好。本工具生成 **Morgan（ECFP 类）环形指纹**，输出 **位向量摘要（哈希预览、开位数）** 并附带 **2D 结构图**，便于对接相似性检索与机器学习特征。

**技能包数据说明**：`morgan-fingerprint/SKILL.md` 使用命令行示例 SMILES；下表摘录文档中的典型值（**可替换**）。首行 SMILES 文件示例同上目录 **`single_smiles_first_line.txt`**。

| 示例标签 | 示例 SMILES |
|----------|-------------|
| 乙醇 | CCO |
| 阿司匹林 | CC(=O)OC1=CC=CC=C1C(=O)O |
| 苯 | c1ccccc1 |

**提交方式（任选其一，不要求必须上传文件）**
1. **指定上表中某一 SMILES**，或在本提示内改写后再发送。
2. 在对话中给出任意新的单行 SMILES。
3. 上传单行 SMILES 文本文件，路径写入 **file_path**。

（助手侧：**smiles_text** 或 **file_path**；禁止臆造路径。）
""",
    "Brenk filter分子毒性检查工具": """[Skill_Route: chem_brenk_filter]
您好。本工具执行 **Brenk 不良片段/毒性预警筛查**（优先 RDKit **BRENK FilterCatalog**，环境缺失时回退至内置 SMARTS 子集），用于早期淘汰高风险母核或官能团。

**技能包数据说明**：`brenk-filter/SKILL.md` 中的命令行示例如下表（**可替换**）；包内 **无独立数据表文件**。首行 SMILES 文件示例见 **`docs/chem_skills_test_inputs/single_smiles_first_line.txt`**（可改为下表 SMILES）。

| 字段 | 示例 SMILES |
|------|-------------|
| 待测分子 | CC(=O)NS(=O)(=O)c1ccc(N)cc1 |

**提交方式（任选其一，不要求必须上传文件）**
1. **修改本提示中的 SMILES** 后发送；或在对话正文给出新的 SMILES。
2. 上传含首行 SMILES 的 `.csv`/`.txt`，路径写入 **file_path**。

（助手侧：**smiles_text** 或 **file_path**。）
""",
    "分子BBB评估工具": """[Skill_Route: chem_bbb_assessment]
您好。本工具基于 **理化描述符组合**给出 **血脑屏障透过倾向的启发式评分**（非临床 PBPK），用于虚拟筛选阶段的 **CNS 暴露风险粗筛**。

**技能包数据说明**：`bbb/SKILL.md` 命令行示例如下表（**可替换**）；包内 **无独立 SMILES 附件**。首行 SMILES 文件示例见 **`docs/chem_skills_test_inputs/single_smiles_first_line.txt`**。

| 字段 | 示例 SMILES |
|------|-------------|
| 待测分子（咖啡因骨架示例） | CN1C=NC2=C1C(=O)N(C(=O)N2C)C |

**提交方式（任选其一，不要求必须上传文件）**
1. **修改本提示中的 SMILES** 后发送；或在对话中给出新的 SMILES。
2. 上传含首行 SMILES 的 `.csv`/`.txt`，路径写入 **file_path**。

（助手侧：**smiles_text** 或 **file_path**。）
""",
    "分子量计算工具": """[Skill_Route: chem_molecular_weight]
您好。本工具使用 RDKit **Descriptors** 计算 **平均分子量与精确分子量**，用于理化性质核对与报表。

**技能包数据说明**：`molecular-weight/SKILL.md` 文档示例（下表）；包内 **无单独 molecules.txt**，表格供用户在会话中 **直接改正文**。多行批量见 **`docs/chem_skills_test_inputs/batch_smiles_one_per_line.txt`**。

| 化合物 | 示例 SMILES |
|--------|-------------|
| 阿司匹林 | CC(=O)OC1=CC=CC=C1C(=O)O |
| 乙醇 | CCO |
| 苯 | c1ccccc1 |

**提交方式（任选其一，不要求必须上传文件）**
1. **在上表中选一列或改写 SMILES** 后发送；或在对话中另写单行 SMILES。
2. 上传含首行 SMILES 的 `.csv`/`.txt`，路径写入 **file_path**。

（助手侧：**smiles_text** 或 **file_path**。）
""",
}
