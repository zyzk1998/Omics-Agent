"""
RDKit SA_Score 在不同发行版中的路径不一致：
- pip / 部分构建：rdkit.Chem.SA_Score
- conda-forge 常见：rdkit.Contrib.SA_Score.sascorer 或 RDContribDir/SA_Score/sascorer.py
"""
from __future__ import annotations

import os
import sys
from typing import Any, Callable


def load_calculate_score() -> Callable[[Any], float]:
    """返回与 sascorer.calculateScore(mol) 等价的可调用对象。"""
    try:
        from rdkit.Chem.SA_Score import calculateScore

        return calculateScore
    except ImportError:
        pass
    try:
        from rdkit.Contrib.SA_Score import sascorer

        return sascorer.calculateScore
    except ImportError:
        pass
    try:
        from rdkit.Contrib.SA_Score.sascorer import calculateScore

        return calculateScore
    except ImportError:
        pass
    try:
        from rdkit.Chem import RDConfig

        d = os.path.join(RDConfig.RDContribDir, "SA_Score")
        p = os.path.join(d, "sascorer.py")
        if os.path.isfile(p):
            if d not in sys.path:
                sys.path.insert(0, d)
            import sascorer as _sc

            return _sc.calculateScore
    except Exception:
        pass
    raise ImportError(
        "无法加载 RDKit SA_Score（已尝试 rdkit.Chem.SA_Score、"
        "rdkit.Contrib.SA_Score、RDContribDir/SA_Score/sascorer.py）"
    ) from None
