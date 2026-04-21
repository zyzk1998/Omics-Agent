# 化学 RDKit 六技能 · 可测试原始数据（与 `SKILL.md` 输入形态一致）

> 与 `技能扩展规范文档.md` 同处 **`docs/`** 目录。本目录下为 **纯文本、可测试** 的输入样例，**依据**为指挥包内各子技能 **`SKILL.md`** 对 `--smiles`、**`--file`（每行一条 SMILES）** 的说明，**非**生产保密数据。

## 原则

1. **生成方向**：仅服务联调、冒烟、回归；SMILES 均从各包 `SKILL.md` 命令行示例或等价的公开结构抄录/派生。  
2. **形态约定**：与 `gibh_agent/assets/chem_rdkit/run_chem_tool.py` 及 `SKILL.md` 一致 — 单分子取 **首行非空、非 `#` 注释**；Tanimoto 取 **前两行** 各一条。  
3. **何时用文件而不用长 Prompt**：当行数多、会占满输入框或超过模型上下文时，**不要**把整表贴进对话，应使用本目录 **`.txt` 原件** 或经 **`/uploads` 上传的附件**（会话里由助手写 `file_path`）。

## 文件说明

| 文件 | 对应 `SKILL.md` 输入 | 适用工具（Tool ID） |
|------|----------------------|----------------------|
| `single_smiles_first_line.txt` | 单行 / 文件首行一条 | `chem_pains_filter`, `chem_morgan_fingerprint`, `chem_brenk_filter`, `chem_bbb_assessment`, `chem_molecular_weight`（单分子） |
| `pair_smiles_two_lines.txt` | 文件内两行各一条 | `chem_tanimoto_similarity` |
| `batch_smiles_one_per_line.txt` | `molecules.txt` 风格，多行 | 单分子类工具**仅消费首行**；本文件用于**压测/长列表**时走附件，或供外部分批脚本；Tanimoto 仍只读前两行 |

## 与 DB 模板的关系

`gibh_agent/db/chem_skill_prompt_templates.py` 中内嵌表格与上述文件**同一数据口径**；更新任一方时建议同步另一方并执行 `scripts/patch_chem_skills_batch.py`（或受控 SQL）更新 `skills.prompt_template`。
