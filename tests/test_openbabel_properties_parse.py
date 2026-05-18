# -*- coding: utf-8 -*-
"""Open Babel --append 性质解析：MW/rotors 映射与 ?? 过滤。"""
from __future__ import annotations

import importlib.util
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[1]
_CONVERTER = _ROOT / "gibh_agent" / "assets" / "openbabel_converter" / "openbabel_converter.py"


def _load_converter():
    spec = importlib.util.spec_from_file_location("obabel_conv", _CONVERTER)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_parse_append_line_tab_separated_six_columns() -> None:
    mod = _load_converter()
    line = "CC(=O)Oc1ccccc1C(=O)O\t180.157 1.3101 63.6 4 1 3"
    vals = mod._parse_obabel_append_line(line, 6)
    assert vals == ["180.157", "1.3101", "63.6", "4", "1", "3"]


def test_parse_append_line_legacy_broken_molwt_nrot() -> None:
    """旧版 molwt/nrot 描述符失败时，性质挤在第二列空格分隔。"""
    mod = _load_converter()
    line = "CC(=O)Oc1ccccc1C(=O)O\t?? 1.3101 63.6 4 1 ??"
    vals = mod._parse_obabel_append_line(line, 6)
    assert vals == ["??", "1.3101", "63.6", "4", "1", "??"]


def test_calculate_properties_aspirin_integration() -> None:
    import shutil

    if not shutil.which("obabel"):
        return
    mod = _load_converter()
    props = mod.calculate_properties("CC(=O)Oc1ccccc1C(=O)O", "smi")
    assert "error" not in props
    assert props.get("molwt") not in (None, "", "??")
    assert props.get("nrot") not in (None, "", "??")
    mw = float(props["molwt"])
    assert 179.0 < mw < 181.5
    assert int(props["nrot"]) == 3
