# 化学 RDKit 闭环验收用例（与指挥包对齐说明）

本目录用于 **`run_chem_tool.py`** 与仓库内 **黄金 JSON** 的一致性验收。

## 样例文件

| 文件 | 用途 |
|------|------|
| `test_smiles_single.txt` | 单分子工具（BBB / PAINS / Brenk / Morgan / MolWt）：首行 SMILES |
| `test_smiles_pair.txt` | Tanimoto：两行，每行一个 SMILES |

当前 SMILES 为 **阿司匹林**（`CC(=O)Oc1ccccc1C(=O)O`）与 **咖啡因**（第二行），用于双分子相似度；若指挥侧压缩包 **`2026.04.20 skills`** 内含 `test_smiles.csv` / README 样例，可将该文件 **覆盖本目录同名输入** 后重新生成黄金基准（见下方命令）。

## 黄金基准（`expected/*.chem_summary.json`）

由 **`gibh_agent/assets/chem_rdkit/run_chem_tool.py`** 在同环境 RDKit 下生成。**数值与指纹摘要随 RDKit 版本可能略有浮点差异**；验收测试以「JSON 深度相等」为准，升级 RDKit 后若失败请按下面命令刷新 `expected/` 并提交 MR 说明版本。

### 重新生成 expected（需在已安装 `rdkit-pypi` 的解释器中执行）

```bash
cd /path/to/GIBH-AGENT-V2
FIX=tests/fixtures/chem_2026_04_20
RUNNER=gibh_agent/assets/chem_rdkit/run_chem_tool.py
for t in bbb pains brenk morgan_fp mol_weight; do
  python "$RUNNER" --tool "$t" --file "$FIX/test_smiles_single.txt" --output-dir "$FIX/_tmp" && \
  cp "$FIX/_tmp/chem_summary.json" "$FIX/expected/${t}.chem_summary.json"
done
python "$RUNNER" --tool tanimoto --file "$FIX/test_smiles_pair.txt" --output-dir "$FIX/_tmp" && \
  cp "$FIX/_tmp/chem_summary.json" "$FIX/expected/tanimoto.chem_summary.json"
rm -rf "$FIX/_tmp"
```

## 与「2026.04.20 skills」压缩包的关系

若已将 **`2026.04.20 skills.zip`** 解压到服务器（例如 `data/uploads/skills_assets/script/inbound/` 下），请将其中的 **化学类测试 CSV/TXT** 复制为本目录 `test_smiles_*.txt`（或并列保留 `test_smiles.csv` 并改写 `run_chem_tool.py` 的 `--file` 指向），再执行 **`pytest tests/test_chem_rdkit_golden.py`**。  
**当前仓库未包含该 zip**；下列黄金文件即本次在 **RDKit 2022.09.x** 上对本目录输入跑出的 **权威输出**，用于 CI 回归。
