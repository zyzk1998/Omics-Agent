# -*- coding: utf-8 -*-
"""BioPython Entrez 解析辅助（DictionaryElement / ListElement / StringElement）。"""
from __future__ import annotations

from typing import Any, Iterable, List, Sequence


def _is_entrez_list(obj: Any) -> bool:
    if isinstance(obj, (list, tuple)):
        return True
    cls = type(obj).__name__
    return cls in ("ListElement",)


def entrez_as_list(obj: Any) -> List[Any]:
    """将 IdList / DocumentSummary 等转为 Python list。"""
    if obj is None:
        return []
    if isinstance(obj, list):
        return obj
    if isinstance(obj, tuple):
        return list(obj)
    if _is_entrez_list(obj):
        return list(obj)
    return [obj]


def entrez_get(obj: Any, key: str, default: Any = None) -> Any:
    """安全读取 Entrez 记录字段（兼容 DictionaryElement）。"""
    if obj is None:
        return default
    if isinstance(obj, dict):
        return obj.get(key, default)
    if hasattr(obj, "keys"):
        try:
            return obj.get(key, default)  # type: ignore[attr-defined]
        except AttributeError:
            try:
                return obj[key]
            except (KeyError, TypeError):
                return default
    return default


def entrez_text(val: Any, default: str = "") -> str:
    if val is None:
        return default
    if isinstance(val, str):
        return val.strip()
    return str(val).strip()


def format_entrez_author_list(author_list: Any) -> str:
    """AuthorList 可能是字符串列表或含 Name 的字典列表。"""
    names: List[str] = []
    for item in entrez_as_list(author_list):
        if isinstance(item, str):
            s = item.strip()
            if s:
                names.append(s)
        elif isinstance(item, dict) or hasattr(item, "keys"):
            nm = entrez_text(entrez_get(item, "Name") or entrez_get(item, "name"))
            if nm:
                names.append(nm)
        elif item is not None:
            s = str(item).strip()
            if s:
                names.append(s)
    return ", ".join(names)[:200]


def iter_pubmed_esummary_docs(summaries: Any) -> Iterable[Any]:
    """
    解析 esummary 返回结构。

    当前 PubMed esummary 常直接返回 ListElement[DocumentSummary, ...]；
    旧版/其它库可能为 DocumentSummarySet → DocumentSummary。
    """
    if summaries is None:
        return

    if _is_entrez_list(summaries):
        for doc in summaries:
            yield doc
        return

    if isinstance(summaries, list):
        for doc in summaries:
            yield doc
        return

    if hasattr(summaries, "keys"):
        dss = entrez_get(summaries, "DocumentSummarySet")
        if dss is not None:
            if hasattr(dss, "keys"):
                raw = entrez_get(dss, "DocumentSummary")
            else:
                raw = dss
            for doc in entrez_as_list(raw):
                yield doc
            return
        raw = entrez_get(summaries, "DocumentSummary")
        for doc in entrez_as_list(raw):
            yield doc


def parse_pubmed_esummary_doc(doc: Any) -> dict[str, str]:
    """单条 DocumentSummary → 表格行。"""
    return {
        "PMID": entrez_text(entrez_get(doc, "Id")),
        "Title": entrez_text(entrez_get(doc, "Title")),
        "Journal": entrez_text(
            entrez_get(doc, "FullJournalName") or entrez_get(doc, "Source")
        ),
        "PubDate": entrez_text(entrez_get(doc, "PubDate")),
        "Authors": format_entrez_author_list(entrez_get(doc, "AuthorList")),
    }
