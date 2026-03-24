"""
轻量级公开 REST API 查询工具（MyGene、UniProt、Reactome、AlphaFold）。
返回 Markdown 摘要，供工作流步骤与 LLM 直接展示。
"""
from __future__ import annotations

import logging
import re
from typing import Any, Dict, List, Optional
from urllib.parse import quote

import requests

from ..core.tool_registry import registry

logger = logging.getLogger(__name__)

_DEFAULT_TIMEOUT = 25
_USER_AGENT = "GIBH-Agent/2.0 (research-api-query; +https://github.com)"


def _md_escape(text: str) -> str:
    if not text:
        return ""
    return str(text).replace("|", "\\|").replace("\n", " ")


def _strip_html(text: str) -> str:
    if not text:
        return ""
    s = re.sub(r"<[^>]+>", "", str(text))
    return re.sub(r"\s+", " ", s).strip()


def _http_get_json(url: str, params: Optional[Dict[str, Any]] = None) -> Any:
    headers = {"User-Agent": _USER_AGENT, "Accept": "application/json"}
    r = requests.get(url, params=params, headers=headers, timeout=_DEFAULT_TIMEOUT)
    r.raise_for_status()
    return r.json()


def _http_get(url: str, params: Optional[Dict[str, Any]] = None) -> requests.Response:
    headers = {"User-Agent": _USER_AGENT}
    return requests.get(url, params=params, headers=headers, timeout=_DEFAULT_TIMEOUT)


@registry.register(
    name="query_gene_info",
    description=(
        "Query MyGene.info for a human gene symbol: full name, genomic location (chromosome/position), "
        "aliases, and NCBI gene summary. Use for gene cards and quick annotation."
    ),
    category="Public API",
    output_type="json",
)
def query_gene_info(gene: str) -> Dict[str, Any]:
    """
    Args:
        gene: Gene symbol or query string (e.g. TP53).
    """
    gene = (gene or "").strip()
    if not gene:
        return {"status": "error", "markdown": "**错误**：请提供基因符号或查询词。", "error": "empty gene"}

    q = f"symbol:{gene}" if ":" not in gene and " " not in gene else gene
    url = "https://mygene.info/v3/query"
    try:
        data = _http_get_json(
            url,
            params={
                "q": q,
                "species": "human",
                "fields": "symbol,name,alias,genomic_pos,summary,entrezgene,ensemblgene,type_of_gene",
                "size": 5,
            },
        )
    except requests.RequestException as e:
        logger.warning("query_gene_info network error: %s", e)
        return {
            "status": "error",
            "markdown": f"**网络错误**（MyGene）：{_md_escape(str(e))}",
            "error": str(e),
        }
    except Exception as e:
        logger.exception("query_gene_info failed")
        return {"status": "error", "markdown": f"**查询失败**：{_md_escape(str(e))}", "error": str(e)}

    hits: List[Dict[str, Any]] = data.get("hits") or []
    if not hits:
        return {
            "status": "success",
            "markdown": f"未在 MyGene 中找到与 **{_md_escape(gene)}** 匹配的条目。",
            "hits": [],
        }

    lines: List[str] = [f"## 基因信息（MyGene）：{_md_escape(gene)}\n"]
    for i, h in enumerate(hits[:3], 1):
        sym = h.get("symbol") or ""
        name = h.get("name") or ""
        aliases = h.get("alias")
        if isinstance(aliases, list):
            alias_s = ", ".join(str(a) for a in aliases[:12])
        elif aliases:
            alias_s = str(aliases)
        else:
            alias_s = ""
        gp = h.get("genomic_pos") or {}
        if isinstance(gp, list) and gp:
            gp = gp[0]
        loc = ""
        if isinstance(gp, dict):
            chrom = gp.get("chr")
            start = gp.get("start")
            end = gp.get("end")
            strand = gp.get("strand")
            if chrom is not None:
                loc = f"chr{chrom}:{start}-{end}"
                if strand is not None:
                    loc += f" (strand {strand})"
        summary = (h.get("summary") or "").strip()
        if len(summary) > 1200:
            summary = summary[:1197] + "..."

        lines.append(f"### 命中 {i}: {_md_escape(sym)} — {_md_escape(name)}\n")
        if loc:
            lines.append(f"- **染色体定位**：{loc}\n")
        if alias_s:
            lines.append(f"- **别名**：{_md_escape(alias_s)}\n")
        if h.get("entrezgene"):
            lines.append(f"- **Entrez Gene ID**：{h.get('entrezgene')}\n")
        if summary:
            lines.append(f"- **摘要**：{summary}\n")
        lines.append("\n")

    md = "".join(lines).strip()
    return {"status": "success", "markdown": md, "hits": hits[:3]}


@registry.register(
    name="query_uniprot_protein",
    description=(
        "Search UniProt KB for a protein (gene name, accession, or keyword): molecular function hints, "
        "protein name, organism, and subcellular location from reviewed entries when available."
    ),
    category="Public API",
    output_type="json",
)
def query_uniprot_protein(protein: str) -> Dict[str, Any]:
    protein = (protein or "").strip()
    if not protein:
        return {"status": "error", "markdown": "**错误**：请提供蛋白名称、基因名或 UniProt ID。", "error": "empty protein"}

    url = "https://rest.uniprot.org/uniprotkb/search"
    pu = protein.strip().upper()
    # accession:基因符号 会触发 UniProt 400，故按形态分流
    if re.match(r"^[A-Z][0-9][A-Z0-9]{3,9}$", pu):
        query = f"(accession:{pu}) AND (organism_id:9606) AND (reviewed:true)"
    else:
        query = f"(gene_exact:{pu}) AND (organism_id:9606) AND (reviewed:true)"
    try:
        data = _http_get_json(
            url,
            params={"query": query, "format": "json", "size": 5},
        )
    except requests.RequestException as e:
        logger.warning("query_uniprot_protein network error: %s", e)
        return {
            "status": "error",
            "markdown": f"**网络错误**（UniProt）：{_md_escape(str(e))}",
            "error": str(e),
        }
    except Exception as e:
        logger.exception("query_uniprot_protein failed")
        return {"status": "error", "markdown": f"**查询失败**：{_md_escape(str(e))}", "error": str(e)}

    results = (data.get("results") or []) if isinstance(data, dict) else []
    if not results:
        return {
            "status": "success",
            "markdown": f"未在 UniProt 中找到与 **{_md_escape(protein)}** 匹配的条目。",
            "results": [],
        }

    lines: List[str] = [f"## UniProt 检索：{_md_escape(protein)}\n\n"]
    for entry in results[:2]:
        acc = entry.get("primaryAccession", "")
        uni = entry.get("uniProtkbId", "")
        desc = (entry.get("proteinDescription") or {}).get("recommendedUniProtName") or {}
        full_name = ""
        if isinstance(desc, dict):
            full = desc.get("fullName") or {}
            full_name = full.get("value") or ""
        org = ((entry.get("organism") or {}).get("scientificName")) or ""
        lines.append(f"### {acc} ({uni})\n")
        lines.append(f"- **名称**：{_md_escape(full_name)}\n")
        if org:
            lines.append(f"- **物种**：{_md_escape(org)}\n")

        func_texts: List[str] = []
        for c in entry.get("comments") or []:
            if c.get("commentType") == "FUNCTION":
                for t in c.get("texts") or []:
                    val = (t or {}).get("value")
                    if val:
                        func_texts.append(val)
        if func_texts:
            ft = " ".join(func_texts)[:1500]
            lines.append(f"- **功能注释**：{ft}\n")

        locs: List[str] = []
        for c in entry.get("comments") or []:
            if c.get("commentType") == "SUBCELLULAR LOCATION":
                for loc in c.get("subcellularLocations") or []:
                    l = loc.get("location") or {}
                    val = l.get("value")
                    if val:
                        locs.append(val)
        if locs:
            lines.append(f"- **亚细胞定位**：{_md_escape('; '.join(locs))}\n")
        lines.append("\n")

    md = "".join(lines).strip()
    return {"status": "success", "markdown": md, "accessions": [r.get("primaryAccession") for r in results[:2]]}


@registry.register(
    name="query_reactome_pathway",
    description=(
        "Search Reactome (human) for pathways by keyword or name; returns short descriptions and stable IDs "
        "for apoptosis, signaling, metabolism, etc."
    ),
    category="Public API",
    output_type="json",
)
def query_reactome_pathway(pathway_query: str) -> Dict[str, Any]:
    pathway_query = (pathway_query or "").strip()
    if not pathway_query:
        return {"status": "error", "markdown": "**错误**：请提供通路名称或关键词。", "error": "empty pathway_query"}

    base = "https://reactome.org/ContentService/search/query"
    try:
        data = _http_get_json(
            base,
            params={"query": pathway_query, "species": "9606", "types": "Pathway"},
        )
    except requests.RequestException as e:
        logger.warning("query_reactome_pathway network error: %s", e)
        return {
            "status": "error",
            "markdown": f"**网络错误**（Reactome）：{_md_escape(str(e))}",
            "error": str(e),
        }
    except Exception as e:
        logger.exception("query_reactome_pathway failed")
        return {"status": "error", "markdown": f"**查询失败**：{_md_escape(str(e))}", "error": str(e)}

    entries: List[Dict[str, Any]] = []
    if isinstance(data, dict):
        for block in data.get("results") or []:
            if isinstance(block, dict):
                for e in block.get("entries") or []:
                    if isinstance(e, dict):
                        entries.append(e)

    if not entries:
        return {
            "status": "success",
            "markdown": f"未在 Reactome（人类）中找到与 **{_md_escape(pathway_query)}** 相关的通路条目。",
            "entries": [],
        }

    lines: List[str] = [f"## Reactome 通路检索：{_md_escape(pathway_query)}\n\n"]
    for item in entries[:5]:
        name = _strip_html(item.get("name") or item.get("displayName") or "")
        st = item.get("stableIdentifier") or item.get("stId") or item.get("id") or ""
        summ_raw = item.get("summation") or item.get("summary") or ""
        if isinstance(summ_raw, dict):
            summ_raw = str(summ_raw)
        summ = _strip_html(str(summ_raw).strip())
        if len(summ) > 800:
            summ = summ[:797] + "..."
        lines.append(f"### {_md_escape(name)}\n")
        if st:
            lines.append(f"- **Stable ID**：`{st}`\n")
        if summ:
            lines.append(f"- **描述**：{summ}\n")
        lines.append("\n")

    md = "".join(lines).strip()
    return {"status": "success", "markdown": md, "count": len(entries)}


@registry.register(
    name="query_alphafold_db",
    description=(
        "Query AlphaFold DB (EBI) by UniProt accession: whether a structure prediction exists and links to "
        "PDB/mmCIF and coverage metadata."
    ),
    category="Public API",
    output_type="json",
)
def query_alphafold_db(uniprot_id: str) -> Dict[str, Any]:
    uniprot_id = (uniprot_id or "").strip().upper()
    if not uniprot_id:
        return {"status": "error", "markdown": "**错误**：请提供 UniProt Accession（如 P04637）。", "error": "empty uniprot_id"}

    url = f"https://alphafold.ebi.ac.uk/api/prediction/{quote(uniprot_id, safe='')}"
    try:
        r = _http_get(url)
        if r.status_code == 404:
            return {
                "status": "success",
                "markdown": f"AlphaFold 数据库中**暂无** UniProt **{uniprot_id}** 的预测条目（404）。",
                "predictions": [],
            }
        r.raise_for_status()
        data = r.json()
    except requests.RequestException as e:
        logger.warning("query_alphafold_db network error: %s", e)
        return {
            "status": "error",
            "markdown": f"**网络错误**（AlphaFold）：{_md_escape(str(e))}",
            "error": str(e),
        }
    except Exception as e:
        logger.exception("query_alphafold_db failed")
        return {"status": "error", "markdown": f"**查询失败**：{_md_escape(str(e))}", "error": str(e)}

    if not isinstance(data, list):
        data = [data] if isinstance(data, dict) else []

    if not data:
        return {
            "status": "success",
            "markdown": f"未返回 **{uniprot_id}** 的 AlphaFold 预测数据。",
            "predictions": [],
        }

    lines: List[str] = [f"## AlphaFold 结构预测：{uniprot_id}\n\n"]
    for pred in data[:3]:
        if not isinstance(pred, dict):
            continue
        entry = pred.get("entryId") or pred.get("uniprotAccession") or ""
        gene = pred.get("gene") or ""
        org = pred.get("organismScientificName") or ""
        pdb = pred.get("pdbUrl") or ""
        cif = pred.get("cifUrl") or ""
        bcif = pred.get("bcifUrl") or ""
        cov = pred.get("coverage") or pred.get("uniprotEnd") or ""
        lines.append(f"### 条目 {entry}\n")
        if gene:
            lines.append(f"- **基因**：{_md_escape(str(gene))}\n")
        if org:
            lines.append(f"- **物种**：{_md_escape(org)}\n")
        if cov:
            lines.append(f"- **覆盖/长度信息**：{_md_escape(str(cov))}\n")
        if pdb:
            lines.append(f"- **PDB 下载/查看**：[链接]({pdb})\n")
        if cif:
            lines.append(f"- **mmCIF**：[链接]({cif})\n")
        if bcif:
            lines.append(f"- **BCIF**：[链接]({bcif})\n")
        lines.append("\n")

    md = "".join(lines).strip()
    return {"status": "success", "markdown": md, "predictions": data[:3]}
