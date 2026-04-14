"""
传统生信微服务适配器：通过 HTTP 调用 worker-pyskills（ViennaRNA RNAfold、RDKit SA_Score 等）。

与《底层算子融合智能体系统标准方法论》对齐：
- 入参优先为已上传文件的绝对路径；若无上传，支持 fasta_content / smiles_text / pdb_content / table_content 内联正文，
  适配器写入 `UPLOAD_DIR/.pyskills_scratch/`（与 worker-pyskills 共享挂载）再调算子；
- 成功时写入 RESULTS_DIR 并返回 csv_url / json_url，供 SkillAgent Markdown 卡片展示。

容器内默认基址：http://worker-pyskills:8000；本地可设环境变量 PYSKILLS_BASE_URL 覆盖。
"""
from __future__ import annotations

import csv
import io
import json
import logging
import os
import re
import uuid
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import httpx

from ..core.executor import WorkflowExecutor
from ..core.tool_registry import registry
from ..core.utils import safe_tool_execution

logger = logging.getLogger(__name__)

PYSKILLS_BASE_URL = (os.getenv("PYSKILLS_BASE_URL") or "http://worker-pyskills:8000").rstrip("/")
_HTTP_TIMEOUT = float(os.getenv("PYSKILLS_HTTP_TIMEOUT", "120"))


def _pyskills_url(path: str) -> str:
    p = path if path.startswith("/") else f"/{path}"
    return f"{PYSKILLS_BASE_URL}{p}"


def _results_root() -> Path:
    raw = (os.getenv("RESULTS_DIR") or "").strip()
    if raw:
        return Path(raw).expanduser().resolve()
    p = Path("/app/results")
    return p.resolve() if p.parent.is_dir() else Path.cwd() / "results"


def _public_results_url(rel_under_results: str) -> str:
    """与 protein_tools._public_path_url 一致：默认 /results/..."""
    rel = rel_under_results.lstrip("/")
    base = (os.getenv("PUBLIC_RESULTS_BASE_URL") or "").rstrip("/")
    path = f"/results/{rel}"
    return f"{base}{path}" if base else path


def _resolve_tool_input_path(file_path: str) -> str:
    ex = WorkflowExecutor()
    try:
        return ex._resolve_file_path(file_path)
    except FileNotFoundError as e:
        raise ValueError(
            f"找不到输入文件：{file_path}。"
            "请确认已通过对话上传并使用列表中的绝对路径；或未上传时在参数中提供内联正文字段（如 fasta_content），"
            "勿编造 `/app/uploads/demo_*.fa` 等示例路径。"
        ) from e


def _uploads_root() -> Path:
    raw = (os.getenv("UPLOAD_DIR") or "/app/uploads").strip()
    return Path(raw).expanduser().resolve()


def _write_scratch_file(content: str, suffix: str) -> Path:
    """
    写入与 worker-pyskills 共享的挂载目录（须为 UPLOAD_DIR 下子路径），
    使容器内 `http://worker-pyskills` 能按同一绝对路径读盘。
    """
    d = _uploads_root() / ".pyskills_scratch"
    d.mkdir(parents=True, exist_ok=True)
    p = d / f"{uuid.uuid4().hex}{suffix}"
    p.write_text(content, encoding="utf-8")
    return p.resolve()


def _try_resolve_path(file_path: str) -> Optional[str]:
    if not (file_path or "").strip():
        return None
    try:
        return _resolve_tool_input_path(file_path.strip())
    except ValueError:
        return None


def _normalize_fasta_inline(raw: str) -> str:
    t = raw.strip().replace("\r\n", "\n")
    if not t:
        raise ValueError("fasta_content 为空")
    if not t.startswith(">"):
        seq = "".join(
            ln.strip()
            for ln in t.splitlines()
            if ln.strip() and not ln.strip().startswith("#")
        )
        seq = re.sub(r"[^AUCGTNaucgtn]", "", seq)
        if len(seq) < 1:
            raise ValueError("fasta_content 中未识别到 RNA 碱基（A/U/G/C/T/N）")
        return f">sequence\n{seq.upper()}\n"
    return t if t.endswith("\n") else t + "\n"


def _resolve_rnafold_input(file_path: str, fasta_content: str) -> str:
    p = _try_resolve_path(file_path)
    if p:
        return p
    if (fasta_content or "").strip():
        body = _normalize_fasta_inline(fasta_content)
        return str(_write_scratch_file(body, ".fa"))
    raise ValueError(
        "请提供已上传文件的 file_path（列表中的真实绝对路径），或在参数 fasta_content 中提供 FASTA 正文"
        "（可含 `>id` 头行；或仅一行 AUCG 序列）。禁止编造不存在的磁盘路径。"
    )


def _normalize_smiles_inline(raw: str) -> str:
    for line in raw.splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        return s.split()[0].strip()
    raise ValueError("smiles_text 中未找到 SMILES")


def _resolve_sascore_input(file_path: str, smiles_text: str) -> str:
    p = _try_resolve_path(file_path)
    if p:
        return p
    if (smiles_text or "").strip():
        line = _normalize_smiles_inline(smiles_text) + "\n"
        return str(_write_scratch_file(line, ".smi"))
    raise ValueError(
        "请提供 file_path（已上传 SMILES 文本文件）或 smiles_text（单行 SMILES）。"
        "禁止编造路径。"
    )


def _resolve_pymol_input(file_path: str, pdb_content: str) -> str:
    p = _try_resolve_path(file_path)
    if p:
        return p
    if (pdb_content or "").strip():
        body = pdb_content.strip().replace("\r\n", "\n")
        if "ATOM" not in body and "HETATM" not in body:
            raise ValueError("pdb_content 需包含 ATOM/HETATM 记录")
        return str(_write_scratch_file(body if body.endswith("\n") else body + "\n", ".pdb"))
    raise ValueError("请提供 file_path（已上传 PDB/mmCIF）或 pdb_content（正文）。")


def _resolve_gseapy_input(file_path: str, table_content: str) -> str:
    p = _try_resolve_path(file_path)
    if p:
        return p
    if (table_content or "").strip():
        body = table_content.strip().replace("\r\n", "\n")
        return str(_write_scratch_file(body if body.endswith("\n") else body + "\n", ".csv"))
    raise ValueError("请提供 file_path（已上传 csv/gmt 等）或 table_content（CSV 正文）。")


def _sanitize_for_client_storage(obj: Union[Dict[str, Any], List[Any], Any]) -> Any:
    """落盘 result.json 前剔除易泄露服务器路径的字段与字符串，禁止将 /app/... 暴露给前端。"""
    if isinstance(obj, dict):
        skip = {"path", "file_path", "absolute_path", "input_path", "output_path", "full_path"}
        return {k: _sanitize_for_client_storage(v) for k, v in obj.items() if k not in skip}
    if isinstance(obj, list):
        return [_sanitize_for_client_storage(x) for x in obj]
    if isinstance(obj, str):
        if obj.startswith("/app/") or obj.startswith("/root/") or "/site-packages/" in obj:
            return "[redacted]"
    return obj


def _write_pyskills_artifacts(
    tool_tag: str,
    payload: Dict[str, Any],
) -> tuple[Optional[str], Optional[str]]:
    """
    将算子 JSON 落盘，返回 (csv_url, json_url)。
    csv_url 供技能卡片主链接；json_url 供完整结果下载。
    worker 已生成的 /uploads/results/... 制品 URL 直接复用，不在此重复拷贝二进制。
    """
    run_id = uuid.uuid4().hex[:12]
    out_dir = _results_root() / "pyskills" / run_id
    out_dir.mkdir(parents=True, exist_ok=True)

    safe_payload = _sanitize_for_client_storage(payload)
    json_path = out_dir / "result.json"
    json_path.write_text(
        json.dumps(safe_payload, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    rel_json = f"pyskills/{run_id}/result.json"
    json_url = _public_results_url(rel_json)

    csv_url: Optional[str] = None
    data = safe_payload.get("data") if isinstance(safe_payload.get("data"), dict) else {}
    wc = data.get("csv_url")
    if isinstance(wc, str) and wc.startswith("/uploads/"):
        csv_url = wc
    if tool_tag == "rnafold" and csv_url is None and "structure" in data and "mfe" in data:
        csv_path = out_dir / "rnafold.csv"
        with csv_path.open("w", encoding="utf-8", newline="") as f:
            w = csv.writer(f)
            w.writerow(
                [
                    "structure",
                    "mfe_kcal_per_mol",
                    "temperature_celsius",
                    "use_temperature_control",
                    "partition_free_energy",
                    "pair_prob_structure",
                ]
            )
            w.writerow(
                [
                    data.get("structure", ""),
                    data.get("mfe"),
                    data.get("temperature_celsius", ""),
                    data.get("use_temperature_control", ""),
                    data.get("partition_free_energy", ""),
                    (data.get("pair_prob_structure") or "")[:8000],
                ]
            )
        csv_url = _public_results_url(f"pyskills/{run_id}/rnafold.csv")
    elif tool_tag == "sascore" and csv_url is None and "sa_score" in data:
        csv_path = out_dir / "sascore.csv"
        with csv_path.open("w", encoding="utf-8", newline="") as f:
            w = csv.writer(f)
            w.writerow(["sa_score"])
            w.writerow([data.get("sa_score")])
        csv_url = _public_results_url(f"pyskills/{run_id}/sascore.csv")
    elif tool_tag == "lipinski" and isinstance(data, dict) and "MW" in data:
        csv_path = out_dir / "lipinski.csv"
        with csv_path.open("w", encoding="utf-8", newline="") as f:
            w = csv.writer(f)
            w.writerow(["MW", "logP", "HBD", "HBA", "lipinski_violations", "is_druglike"])
            w.writerow(
                [
                    data.get("MW"),
                    data.get("logP"),
                    data.get("HBD"),
                    data.get("HBA"),
                    data.get("lipinski_violations", ""),
                    data.get("is_druglike"),
                ]
            )
        csv_url = _public_results_url(f"pyskills/{run_id}/lipinski.csv")
    # drug_similarity：真实 Top-N 表由 worker + drug_similarity.py 写入 similarity_results.csv，经 data.csv_url 复用；不再生成 query/html 摘要 CSV。

    if (
        not csv_url
        and not (tool_tag == "pymol" and isinstance(data.get("image_url"), str))
        and tool_tag != "drug_similarity"
    ):
        # 无表格字段时导出 key-value 摘要；PyMOL 仅图像时可不生成无意义 CSV
        buf = io.StringIO()
        w = csv.writer(buf)
        w.writerow(["key", "value"])
        for k, v in (data or {}).items():
            w.writerow([k, v])
        fallback = out_dir / "summary.csv"
        fallback.write_text(buf.getvalue(), encoding="utf-8")
        csv_url = _public_results_url(f"pyskills/{run_id}/summary.csv")

    return csv_url, json_url


def _normalize_worker_success(
    tool_tag: str,
    body: Dict[str, Any],
) -> Dict[str, Any]:
    """合并 worker 响应与方法论约定的 html/csv/json 链接；禁止向下游透传服务器物理路径。"""
    if body.get("status") != "success":
        return body
    safe = _sanitize_for_client_storage(body)
    if not isinstance(safe, dict):
        safe = dict(body)
    csv_u, json_u = _write_pyskills_artifacts(tool_tag, safe)
    if tool_tag == "rnafold":
        msg = "RNA 二级结构预测完成，结构以点括号表示，MFE 单位为 kcal/mol。"
        if isinstance(safe.get("data"), dict) and safe["data"].get("pair_prob_structure"):
            msg += " 已包含分区函数碱基配对概率串（可用于 dot plot）；详见 JSON / CSV。"
    elif tool_tag == "sascore":
        msg = "合成可行性（SA Score）计算完成，数值越低通常越易合成。"
    elif tool_tag == "lipinski":
        msg = "Lipinski 五规则评估完成（MW、LogP、氢键供体/受体与违反条数）。"
    elif tool_tag == "drug_similarity":
        msg = (
            (safe.get("message") if isinstance(safe.get("message"), str) else None)
            or "结构相似性检索完成，详见 HTML 报告链接。"
        )
    elif tool_tag == "pymol":
        msg = (safe.get("message") if isinstance(safe.get("message"), str) else None) or "蛋白质三维结构卡通渲染已完成。"
    elif tool_tag == "gseapy":
        msg = (safe.get("message") if isinstance(safe.get("message"), str) else None) or "基因集富集分析已完成，详见结果表与条形图。"
    else:
        msg = "PySkills 执行完成。"
    out = dict(safe)
    out["message"] = msg
    if csv_u:
        out["csv_url"] = csv_u
    if json_u:
        out["json_url"] = json_u
    d = out.get("data")
    if isinstance(d, dict):
        for ukey in ("image_url", "csv_url", "html_url", "pdf_url"):
            if isinstance(d.get(ukey), str) and d[ukey].startswith("/uploads/") and not out.get(ukey):
                out[ukey] = d[ukey]
    # drug_similarity：仅 similarity_report.html + similarity_results.csv（由 DrugSimilarityTool
    # 与 find_similar 同一检索面写入）；此处只做 URL 对齐，不存在「精简版 HTML」第二路径。
    if tool_tag == "drug_similarity" and isinstance(out.get("data"), dict):
        dd = out["data"]
        for ukey in ("html_url", "csv_url"):
            v = dd.get(ukey)
            if isinstance(v, str) and v.strip().startswith("/uploads/"):
                out[ukey] = v.strip()
    return out


@registry.register(
    name="rnafold_analysis",
    description=(
        "调用 ViennaRNA 对 RNA FASTA 进行二级结构预测（MFE + 点括号）。"
        "输入二选一：file_path（已上传文件的真实绝对路径）或 fasta_content（FASTA 正文，可仅一行碱基序列）。"
        "可选：temperature_celsius、use_temperature_control、include_pairing_probabilities。"
    ),
    category="Bioinformatics",
    output_type="json",
)
@safe_tool_execution
async def rnafold_analysis(
    file_path: str = "",
    fasta_content: str = "",
    temperature_celsius: float = 37.0,
    use_temperature_control: bool = True,
    include_pairing_probabilities: bool = False,
) -> dict:
    try:
        resolved = _resolve_rnafold_input(file_path, fasta_content)
    except ValueError as e:
        return {"status": "error", "message": str(e)}

    try:
        async with httpx.AsyncClient(timeout=_HTTP_TIMEOUT) as client:
            resp = await client.post(
                _pyskills_url("/api/rnafold"),
                json={
                    "file_path": resolved,
                    "temperature_celsius": temperature_celsius,
                    "use_temperature_control": use_temperature_control,
                    "include_pairing_probabilities": include_pairing_probabilities,
                },
            )
            resp.raise_for_status()
    except httpx.HTTPStatusError as e:
        r = e.response
        detail = _http_error_detail(r)
        logger.warning("rnafold_analysis HTTP %s: %s", r.status_code, detail[:500])
        return {"status": "error", "message": _friendly_http_error(detail), "http_status": r.status_code}
    except httpx.RequestError as e:
        logger.warning("rnafold_analysis 连接失败: %s", e)
        return {
            "status": "error",
            "message": (
                "无法连接 PySkills 微服务（请确认已启动 worker-pyskills 容器，且 "
                f"PYSKILLS_BASE_URL={PYSKILLS_BASE_URL} 可达）。详情：{e}"
            ),
        }

    try:
        body = resp.json()
    except json.JSONDecodeError:
        return {"status": "error", "message": "PySkills 返回非 JSON 响应"}

    if isinstance(body, dict) and body.get("status") == "success":
        return _normalize_worker_success("rnafold", body)
    if isinstance(body, dict) and body.get("status") != "success":
        return {"status": "error", "message": body.get("message") or str(body)}
    return body if isinstance(body, dict) else {"status": "error", "message": str(body)}


@registry.register(
    name="sascore_analysis",
    description=(
        "使用 RDKit SA Score 评估合成可行性。"
        "输入二选一：file_path（已上传文件）或 smiles_text（单行 SMILES 字符串）。"
    ),
    category="Bioinformatics",
    output_type="json",
)
@safe_tool_execution
async def sascore_analysis(file_path: str = "", smiles_text: str = "") -> dict:
    try:
        resolved = _resolve_sascore_input(file_path, smiles_text)
    except ValueError as e:
        return {"status": "error", "message": str(e)}

    try:
        async with httpx.AsyncClient(timeout=_HTTP_TIMEOUT) as client:
            resp = await client.post(
                _pyskills_url("/api/sascore"),
                json={"file_path": resolved},
            )
            resp.raise_for_status()
    except httpx.HTTPStatusError as e:
        r = e.response
        detail = _http_error_detail(r)
        logger.warning("sascore_analysis HTTP %s: %s", r.status_code, detail[:500])
        return {"status": "error", "message": _friendly_http_error(detail), "http_status": r.status_code}
    except httpx.RequestError as e:
        logger.warning("sascore_analysis 连接失败: %s", e)
        return {
            "status": "error",
            "message": (
                "无法连接 PySkills 微服务（请确认已启动 worker-pyskills 容器，且 "
                f"PYSKILLS_BASE_URL={PYSKILLS_BASE_URL} 可达）。详情：{e}"
            ),
        }

    try:
        body = resp.json()
    except json.JSONDecodeError:
        return {"status": "error", "message": "PySkills 返回非 JSON 响应"}

    if isinstance(body, dict) and body.get("status") == "success":
        return _normalize_worker_success("sascore", body)
    if isinstance(body, dict) and body.get("status") != "success":
        return {"status": "error", "message": body.get("message") or str(body)}
    return body if isinstance(body, dict) else {"status": "error", "message": str(body)}


_DRUG_SIM_HTTP_TIMEOUT = float(os.getenv("PYSKILLS_DRUG_SIM_TIMEOUT", "920"))


@registry.register(
    name="lipinski_druglikeness",
    description=(
        "基于 Lipinski Rule of Five 评估小分子口服成药潜势。"
        "输入二选一：file_path（已上传 SMILES 文本文件）或 smiles_text（单行 SMILES）。"
        "输出 MW、logP、HBD、HBA、is_druglike（违反规则条数≤1 为 true）。"
    ),
    category="Bioinformatics",
    output_type="json",
)
@safe_tool_execution
async def lipinski_druglikeness(file_path: str = "", smiles_text: str = "") -> dict:
    try:
        resolved = _resolve_sascore_input(file_path, smiles_text)
    except ValueError as e:
        return {"status": "error", "message": str(e)}

    try:
        async with httpx.AsyncClient(timeout=_HTTP_TIMEOUT) as client:
            resp = await client.post(_pyskills_url("/api/lipinski"), json={"file_path": resolved})
            resp.raise_for_status()
    except httpx.HTTPStatusError as e:
        r = e.response
        detail = _http_error_detail(r)
        logger.warning("lipinski_druglikeness HTTP %s: %s", r.status_code, detail[:500])
        return {"status": "error", "message": _friendly_http_error(detail), "http_status": r.status_code}
    except httpx.RequestError as e:
        logger.warning("lipinski_druglikeness 连接失败: %s", e)
        return {
            "status": "error",
            "message": (
                "无法连接 PySkills 微服务（请确认已启动 worker-pyskills 容器，且 "
                f"PYSKILLS_BASE_URL={PYSKILLS_BASE_URL} 可达）。详情：{e}"
            ),
        }

    try:
        body = resp.json()
    except json.JSONDecodeError:
        return {"status": "error", "message": "PySkills 返回非 JSON 响应"}

    if isinstance(body, dict) and body.get("status") == "success":
        return _normalize_worker_success("lipinski", body)
    if isinstance(body, dict) and body.get("status") != "success":
        return {"status": "error", "message": body.get("message") or str(body)}
    return body if isinstance(body, dict) else {"status": "error", "message": str(body)}


@registry.register(
    name="drug_similarity_search",
    description=(
        "在 PubChem/ChEMBL 等公共库中检索与查询分子结构相似的化合物，生成交互式 HTML 报告（需联网，耗时较长）。"
        "输入二选一：file_path（已上传 SMILES 文本）或 smiles_text（单行 SMILES）。"
        "可选：top_n（默认 10）、fingerprint（默认 morgan_ecfp）、metric（默认 tanimoto）、"
        "databases_csv（逗号分隔库名，如 pubchem,chembl）；"
        "similarity_threshold（默认 0.5，与原始包 find_similar.py 中 get_similar_from_all 一致）、"
        "max_results_per_db（默认 100，每库相似检索候选上限）。"
    ),
    category="Bioinformatics",
    output_type="json",
)
@safe_tool_execution
async def drug_similarity_search(
    file_path: str = "",
    smiles_text: str = "",
    top_n: int = 10,
    fingerprint: str = "morgan_ecfp",
    metric: str = "tanimoto",
    databases_csv: str = "pubchem,chembl",
    similarity_threshold: float = 0.5,
    max_results_per_db: int = 100,
) -> dict:
    try:
        resolved = _resolve_sascore_input(file_path, smiles_text)
    except ValueError as e:
        return {"status": "error", "message": str(e)}

    dbs = [x.strip().lower() for x in (databases_csv or "").split(",") if x.strip()]
    if not dbs:
        dbs = ["pubchem", "chembl"]
    payload: Dict[str, Any] = {
        "file_path": resolved,
        "top_n": int(top_n),
        "fingerprint": (fingerprint or "morgan_ecfp").strip().lower(),
        "metric": (metric or "tanimoto").strip().lower(),
        "databases": dbs,
        # 与 worker CLI `--similarity-threshold` / `--max-results-per-db` 同源，禁止在编排层静默改写
        "similarity_threshold": float(similarity_threshold),
        "max_results_per_db": int(max_results_per_db),
    }

    try:
        async with httpx.AsyncClient(timeout=_DRUG_SIM_HTTP_TIMEOUT) as client:
            resp = await client.post(_pyskills_url("/api/drug_similarity"), json=payload)
            resp.raise_for_status()
    except httpx.HTTPStatusError as e:
        r = e.response
        detail = _http_error_detail(r)
        logger.warning("drug_similarity_search HTTP %s: %s", r.status_code, detail[:500])
        return {"status": "error", "message": _friendly_http_error(detail), "http_status": r.status_code}
    except httpx.RequestError as e:
        logger.warning("drug_similarity_search 连接失败: %s", e)
        return {
            "status": "error",
            "message": (
                "无法连接 PySkills 微服务（请确认已启动 worker-pyskills 容器，且 "
                f"PYSKILLS_BASE_URL={PYSKILLS_BASE_URL} 可达）。详情：{e}"
            ),
        }

    try:
        body = resp.json()
    except json.JSONDecodeError:
        return {"status": "error", "message": "PySkills 返回非 JSON 响应"}

    if isinstance(body, dict) and body.get("status") == "success":
        return _normalize_worker_success("drug_similarity", body)
    if isinstance(body, dict) and body.get("status") != "success":
        return {"status": "error", "message": body.get("message") or str(body)}
    return body if isinstance(body, dict) else {"status": "error", "message": str(body)}


def _http_error_detail(resp: httpx.Response) -> str:
    try:
        body = resp.json()
        if isinstance(body, dict) and "detail" in body:
            d = body["detail"]
            return d if isinstance(d, str) else str(d)
    except Exception:
        pass
    return (resp.text or "")[:4000] or resp.reason_phrase


def _friendly_http_error(detail: str) -> str:
    """将 FastAPI traceback 等转为更短的用户可读说明。"""
    d = (detail or "").strip()
    if len(d) > 800:
        head = d[:800]
        if "Traceback" in d or "File \"" in d:
            return (
                "PySkills 服务端处理失败（详见后台日志）。"
                "请确认输入文件格式正确（RNAfold 需 FASTA；SA Score 需 SMILES 文本）。"
            )
        return head + "…"
    return d or "PySkills 请求失败"


@registry.register(
    name="pymol_analysis",
    description=(
        "调用 PySkills：PyMOL 无头渲染蛋白质结构（白底、cartoon、二级结构着色），输出 PNG 预览。"
        "输入二选一：file_path（已上传 PDB/mmCIF）或 pdb_content（PDB 正文）。"
    ),
    category="Bioinformatics",
    output_type="json",
)
@safe_tool_execution
async def pymol_analysis(file_path: str = "", pdb_content: str = "") -> dict:
    try:
        resolved = _resolve_pymol_input(file_path, pdb_content)
    except ValueError as e:
        return {"status": "error", "message": str(e)}
    try:
        async with httpx.AsyncClient(timeout=_HTTP_TIMEOUT) as client:
            resp = await client.post(_pyskills_url("/api/pymol"), json={"file_path": resolved})
            resp.raise_for_status()
    except httpx.HTTPStatusError as e:
        r = e.response
        detail = _http_error_detail(r)
        return {"status": "error", "message": _friendly_http_error(detail), "http_status": r.status_code}
    except httpx.RequestError as e:
        return {
            "status": "error",
            "message": f"无法连接 PySkills 微服务：{e}",
        }
    try:
        body = resp.json()
    except json.JSONDecodeError:
        return {"status": "error", "message": "PySkills 返回非 JSON 响应"}
    if isinstance(body, dict) and body.get("status") == "success":
        return _normalize_worker_success("pymol", body)
    if isinstance(body, dict) and body.get("status") != "success":
        return {"status": "error", "message": body.get("message") or str(body)}
    return body if isinstance(body, dict) else {"status": "error", "message": str(body)}


@registry.register(
    name="gseapy_analysis",
    description=(
        "调用 PySkills：基于 gseapy.enrichr 对人类基因列表做 KEGG_2021_Human 富集，输出结果 CSV 与条形图 PNG。"
        "输入二选一：file_path（基因列表 TXT/CSV）或 table_content（含 gene 列或首列为基因名的 CSV）。"
    ),
    category="Bioinformatics",
    output_type="json",
)
@safe_tool_execution
async def gseapy_analysis(file_path: str = "", table_content: str = "") -> dict:
    try:
        resolved = _resolve_gseapy_input(file_path, table_content)
    except ValueError as e:
        return {"status": "error", "message": str(e)}
    try:
        async with httpx.AsyncClient(timeout=_HTTP_TIMEOUT) as client:
            resp = await client.post(_pyskills_url("/api/gseapy"), json={"file_path": resolved})
            resp.raise_for_status()
    except httpx.HTTPStatusError as e:
        r = e.response
        detail = _http_error_detail(r)
        return {"status": "error", "message": _friendly_http_error(detail), "http_status": r.status_code}
    except httpx.RequestError as e:
        return {
            "status": "error",
            "message": f"无法连接 PySkills 微服务：{e}",
        }
    try:
        body = resp.json()
    except json.JSONDecodeError:
        return {"status": "error", "message": "PySkills 返回非 JSON 响应"}
    if isinstance(body, dict) and body.get("status") == "success":
        return _normalize_worker_success("gseapy", body)
    if isinstance(body, dict) and body.get("status") != "success":
        return {"status": "error", "message": body.get("message") or str(body)}
    return body if isinstance(body, dict) else {"status": "error", "message": str(body)}
