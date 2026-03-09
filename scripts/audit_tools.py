#!/usr/bin/env python3
"""
Full Tool Inventory Audit.
Read-only: does NOT call sync_tools(clear_existing=True). Just loads tools and reports by category.
"""
import os
import sys
from pathlib import Path
from collections import defaultdict

# Project root
project_root = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(project_root))
os.chdir(project_root)

# Optional: avoid numba cache errors when importing full agent stack
os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba_cache")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl")

# 1) Get registry reference first (singleton)
from gibh_agent.core.tool_registry import registry

# 2) Import gibh_agent.tools so __init__.py runs and loads ALL 7 modality modules
import gibh_agent.tools  # noqa: F401

# 3) Access registry._tools (the single source of truth)
tools_dict = registry._tools
total = len(tools_dict)

# 4) Group by category
by_category = defaultdict(list)
for name, meta in tools_dict.items():
    cat = meta.category
    by_category[cat].append(name)

# 5) Expected 7 modalities (must have at least 1 tool each)
EXPECTED_MODALITIES = [
    "scRNA-seq",
    "Metabolomics",
    "Genomics",
    "Transcriptomics",
    "Epigenomics",
    "Proteomics",
    "Spatial",
    "Radiomics",
]

# 6) Print table: Category | Count | Example Tool Names
print("=" * 80)
print("TOOL INVENTORY AUDIT (read-only, no ChromaDB overwrite)")
print("=" * 80)
print(f"Total tools in registry: {total}")
print()

# Sort categories for stable output (expected first, then rest)
all_cats = list(by_category.keys())
ordered_cats = [c for c in EXPECTED_MODALITIES if c in all_cats]
ordered_cats += [c for c in sorted(all_cats) if c not in EXPECTED_MODALITIES]

max_cat_len = max(len(c) for c in all_cats) if all_cats else 10
max_cat_len = max(max_cat_len, 14)
fmt = f"{{:<{max_cat_len}}} | {{:>5}} | {{}}"

print(fmt.format("Category", "Count", "Example Tool Names"))
print("-" * 80)

missing_modalities = []
for cat in ordered_cats:
    names = by_category[cat]
    count = len(names)
    examples = ", ".join(sorted(names)[:5])
    if len(names) > 5:
        examples += ", ..."
    print(fmt.format(cat, count, examples))
    if cat in EXPECTED_MODALITIES and count == 0:
        missing_modalities.append(cat)

# 7) Flag missing modalities
print()
if missing_modalities:
    print("CRITICAL ERROR: The following modalities have 0 tools registered:")
    for m in missing_modalities:
        print(f"  - {m}")
    sys.exit(1)

missing_from_list = [m for m in EXPECTED_MODALITIES if m not in all_cats]
if missing_from_list:
    print("CRITICAL ERROR: The following modalities have 0 tools (category not present):")
    for m in missing_from_list:
        print(f"  - {m}")
    sys.exit(1)

print("All 7 modalities have at least one tool registered.")
print("=" * 80)
