"""
统一资源定位符（Unified Asset Locator）

将「新上传 / 数据资产 / 历史会话路径」归一为容器内可验证的绝对路径。
"""
from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)


class AssetLocatorError(Exception):
    """无法解析或定位持久化/上传资产。"""

    def __init__(
        self,
        message: str,
        *,
        raw_path: str = "",
        attempted: Optional[List[str]] = None,
        asset_id: Optional[int] = None,
    ) -> None:
        super().__init__(message)
        self.raw_path = raw_path
        self.attempted = attempted or []
        self.asset_id = asset_id


def _default_upload_dir() -> Path:
    return Path(os.getenv("UPLOAD_DIR", "/app/uploads")).expanduser()


def _lookup_asset_path_from_db(
    asset_id: int,
    owner_id: Optional[str],
    db: Any,
) -> Optional[str]:
    if db is None or not asset_id:
        return None
    try:
        from gibh_agent.db.models import Asset as AssetModel

        q = db.query(AssetModel).filter(AssetModel.id == int(asset_id))
        if owner_id:
            q = q.filter(AssetModel.owner_id == owner_id)
        row = q.first()
        if row and row.file_path:
            return str(row.file_path).strip()
        # owner_id 与入库不一致时（游客/用户名漂移），放宽 owner 再查一次
        if owner_id:
            row = (
                db.query(AssetModel)
                .filter(AssetModel.id == int(asset_id))
                .first()
            )
            if row and row.file_path:
                logger.warning(
                    "[AssetLocator] asset_id=%s owner 不匹配(%s vs %s)，仍采用 DB 路径",
                    asset_id,
                    owner_id,
                    row.owner_id,
                )
                return str(row.file_path).strip()
    except Exception as exc:
        logger.warning("[AssetLocator] DB 查询 asset_id=%s 失败: %s", asset_id, exc)
    return None


def _lookup_asset_by_name_from_db(
    file_name: str,
    owner_id: Optional[str],
    db: Any,
) -> Optional[str]:
    """无 asset_id 时按 file_name + owner 回退查库。"""
    if db is None or not file_name or not str(file_name).strip():
        return None
    try:
        from gibh_agent.db.models import Asset as AssetModel

        name = str(file_name).strip()
        q = db.query(AssetModel).filter(AssetModel.file_name == name)
        if owner_id:
            q = q.filter(AssetModel.owner_id == owner_id)
        row = q.order_by(AssetModel.created_at.desc()).first()
        if row and row.file_path:
            return str(row.file_path).strip()
    except Exception as exc:
        logger.warning("[AssetLocator] DB 按文件名查询失败 name=%s: %s", file_name, exc)
    return None


def is_data_asset_entry(file_info: Dict[str, Any]) -> bool:
    """是否为「数据资产库」选中条目（含 asset_id 或 source 标记）。"""
    if not isinstance(file_info, dict):
        return False
    raw_id = file_info.get("asset_id")
    try:
        has_id = raw_id is not None and str(raw_id).strip() != ""
    except Exception:
        has_id = False
    src = str(file_info.get("source") or "").lower()
    return has_id or src in ("data_asset", "asset", "data")


def bridge_data_asset_to_physical_path(
    file_info: Dict[str, Any],
    *,
    upload_dir: Optional[str] = None,
    owner_id: Optional[str] = None,
    db: Any = None,
    context_paths: Optional[List[str]] = None,
    uploaded_entries: Optional[List[Dict[str, Any]]] = None,
) -> Tuple[str, List[str]]:
    """
    数据资产实体映射：优先 asset_id 查库得到权威物理路径，再前置 fuzzy heal。

    Returns:
        (resolved_path, warnings)
    """
    warnings: List[str] = []
    name = (
        file_info.get("file_name")
        or file_info.get("name")
        or ""
    )
    zombie_path = str(file_info.get("file_path") or file_info.get("path") or "").strip()
    asset_id = file_info.get("asset_id")
    try:
        asset_id = int(asset_id) if asset_id is not None and str(asset_id).strip() else None
    except (TypeError, ValueError):
        asset_id = None

    authoritative = ""
    if asset_id and db is not None:
        authoritative = _lookup_asset_path_from_db(asset_id, owner_id, db) or ""
        if authoritative and zombie_path and authoritative.rstrip("/") != zombie_path.rstrip("/"):
            logger.info(
                "[AssetBridge] asset_id=%s 替换僵尸路径: %s -> %s",
                asset_id,
                zombie_path,
                authoritative,
            )
            warnings.append(
                f"数据资产 #{asset_id} 已映射为服务器当前物理路径（忽略前端旧路径）。"
            )
    if not authoritative and is_data_asset_entry(file_info) and db is not None:
        authoritative = _lookup_asset_by_name_from_db(name or Path(zombie_path).name, owner_id, db) or ""
        if authoritative:
            warnings.append(f"数据资产按文件名 {name!r} 查库映射成功。")

    raw_for_heal = authoritative or zombie_path
    if not raw_for_heal:
        return "", warnings

    ud = Path(upload_dir) if upload_dir else _default_upload_dir()
    try:
        resolved = resolve_asset_path(
            raw_for_heal,
            upload_dir=str(ud),
            asset_id=asset_id,
            owner_id=owner_id,
            db=db,
            file_name_hint=name or None,
            context_paths=context_paths,
            uploaded_entries=uploaded_entries,
            must_exist=False,
        )
        p = Path(resolved)
        if p.exists() or (not authoritative and resolved):
            if resolved != zombie_path and zombie_path:
                logger.info("[AssetBridge] 路径已更新: %s -> %s", zombie_path, resolved)
            return resolved, warnings
    except AssetLocatorError as exc:
        warnings.append(str(exc))

    from gibh_agent.core.omics_io_registry import fuzzy_heal_path_by_basename

    healed, _ = fuzzy_heal_path_by_basename(
        raw_for_heal,
        upload_dir=ud,
        owner_id=owner_id,
        context_paths=context_paths,
        uploaded_entries=uploaded_entries,
        allow_dir=True,
    )
    if healed:
        if healed != zombie_path:
            warnings.append(f"Basename 模糊自愈: {Path(zombie_path).name} -> {healed}")
        return healed, warnings

    return raw_for_heal, warnings


def bridge_uploaded_files_at_entry(
    files: List[Any],
    *,
    upload_dir: Optional[str] = None,
    owner_id: Optional[str] = None,
    db: Any = None,
    strict: bool = False,
) -> Tuple[List[Dict[str, Any]], List[str]]:
    """
    系统第一入口：数据资产 DB 映射 + 全量 preemptive fuzzy heal。

    在 FileInspector / Planner / Executor 之前调用；strict 仅影响是否抛错，默认软放行。
    """
    all_warnings: List[str] = []
    raw_entries: List[Dict[str, Any]] = []
    for fi in files or []:
        if isinstance(fi, str):
            fi = {"path": fi, "name": os.path.basename(fi)}
        if isinstance(fi, dict):
            raw_entries.append(dict(fi))

    context_paths: List[str] = []
    for entry in raw_entries:
        for k in ("path", "file_path", "group_dir"):
            v = entry.get(k)
            if v:
                context_paths.append(str(v))

    out: List[Dict[str, Any]] = []
    for fi in raw_entries:
        name = fi.get("file_name") or fi.get("name") or ""
        if is_data_asset_entry(fi) or fi.get("asset_id"):
            resolved, warns = bridge_data_asset_to_physical_path(
                fi,
                upload_dir=upload_dir,
                owner_id=owner_id,
                db=db,
                context_paths=context_paths,
                uploaded_entries=raw_entries,
            )
            all_warnings.extend(warns)
        else:
            resolved, warns = bridge_data_asset_to_physical_path(
                {**fi, "source": fi.get("source")},
                upload_dir=upload_dir,
                owner_id=owner_id,
                db=db,
                context_paths=context_paths,
                uploaded_entries=raw_entries,
            )
            all_warnings.extend(warns)

        if not resolved:
            if strict:
                raise AssetLocatorError(
                    f"无法解析上传条目: {fi!r}",
                    raw_path=str(fi.get("path") or ""),
                    asset_id=fi.get("asset_id"),
                )
            all_warnings.append(f"未能定位文件条目: {name or fi}")
            if fi.get("path") or fi.get("file_path"):
                out.append(
                    {
                        "name": name or os.path.basename(str(fi.get("path") or "")),
                        "path": str(fi.get("path") or fi.get("file_path")),
                        "source": fi.get("source"),
                        "asset_id": fi.get("asset_id"),
                        "locator_warning": "路径未能验证存在，已软放行",
                    }
                )
            continue

        entry: Dict[str, Any] = {
            "name": name or Path(resolved).name,
            "path": resolved,
        }
        if fi.get("source"):
            entry["source"] = fi["source"]
        if fi.get("asset_id") is not None:
            entry["asset_id"] = fi["asset_id"]
        if fi.get("is_10x") is not None:
            entry["is_10x"] = bool(fi.get("is_10x"))
        if fi.get("group_dir"):
            try:
                entry["group_dir"] = resolve_asset_path(
                    str(fi["group_dir"]),
                    upload_dir=upload_dir,
                    owner_id=owner_id,
                    db=db,
                    context_paths=context_paths,
                    uploaded_entries=raw_entries,
                    must_exist=False,
                )
            except AssetLocatorError:
                entry["group_dir"] = fi["group_dir"]
        if warns:
            entry["locator_warnings"] = warns
        out.append(entry)

    return out, all_warnings


def _resolve_with_container_rules(
    raw: str,
    upload_dir: Path,
    owner_id: Optional[str] = None,
) -> Tuple[Optional[str], List[str]]:
    from gibh_agent.utils.path_resolver import (
        normalize_duplicate_tail_filename,
        resolve_upload_path_for_container,
    )

    attempted: List[str] = []
    s = normalize_duplicate_tail_filename(str(raw or "").strip())
    if not s:
        return None, attempted

    try:
        resolved = resolve_upload_path_for_container(s, str(upload_dir))
        attempted.append(resolved)
        p = Path(resolved)
        if p.exists() and (p.is_file() or p.is_dir()):
            return str(p.resolve()), attempted
    except Exception:
        pass

    # 绝对路径直查
    p0 = Path(s)
    if p0.is_absolute():
        attempted.append(str(p0))
        try:
            rp = p0.resolve()
            if rp.exists():
                return str(rp), attempted
        except OSError:
            pass

    # 相对 uploads 子路径
    rel = s.lstrip("/")
    if rel.startswith("uploads/"):
        rel = rel[len("uploads/") :]
    cand = (upload_dir / rel).resolve()
    attempted.append(str(cand))
    if cand.exists():
        return str(cand), attempted

    # basename 仅在 owner 目录下递归搜索（严禁全局 UPLOAD_DIR 盲搜，防跨用户串读）
    base = p0.name
    if base and len(base) >= 3 and owner_id and str(owner_id).strip():
        owner_root = (upload_dir / str(owner_id).strip()).resolve()
        attempted.append(str(owner_root / "**" / base))
        try:
            if owner_root.exists():
                for hit in owner_root.rglob(base):
                    if hit.is_file() or hit.is_dir():
                        try:
                            hit.resolve().relative_to(owner_root)
                        except ValueError:
                            continue
                        hit_s = str(hit.resolve())
                        attempted.append(hit_s)
                        return hit_s, attempted
        except OSError:
            pass

    return None, attempted


def resolve_asset_path(
    raw_path: str = "",
    *,
    upload_dir: Optional[str] = None,
    asset_id: Optional[int] = None,
    owner_id: Optional[str] = None,
    db: Any = None,
    file_name_hint: Optional[str] = None,
    context_paths: Optional[List[str]] = None,
    uploaded_entries: Optional[List[Dict[str, Any]]] = None,
    must_exist: bool = True,
) -> str:
    """
    解析单条资产路径为容器内绝对路径。

    优先级：asset_id → DB.file_path → raw_path → owner 目录下 basename 搜索。
    """
    ud = Path(upload_dir) if upload_dir else _default_upload_dir()
    attempted: List[str] = []
    candidates: List[str] = []

    if asset_id:
        db_path = _lookup_asset_path_from_db(int(asset_id), owner_id, db)
        if db_path:
            candidates.append(db_path)

    if raw_path and str(raw_path).strip():
        candidates.append(str(raw_path).strip())

    if file_name_hint and file_name_hint.strip():
        hint = file_name_hint.strip()
        if owner_id and str(owner_id).strip():
            candidates.append(str(ud / str(owner_id).strip() / hint))
        elif not owner_id:
            candidates.append(str(ud / hint))

    seen: set[str] = set()
    for cand in candidates:
        key = cand.rstrip("/")
        if not key or key in seen:
            continue
        seen.add(key)
        resolved, tries = _resolve_with_container_rules(cand, ud, owner_id=owner_id)
        attempted.extend(t for t in tries if t not in attempted)
        if resolved:
            logger.info("[AssetLocator] 解析成功: %s -> %s", cand, resolved)
            return resolved

    # Basename 模糊自愈（上下文 + owner 隔离）— 对每个候选路径的 basename 尝试
    from gibh_agent.core.omics_io_registry import fuzzy_heal_path_by_basename

    heal_sources = list(candidates)
    if raw_path and str(raw_path).strip():
        heal_sources.append(str(raw_path).strip())
    seen_heal: set[str] = set()
    for src in heal_sources:
        key = src.rstrip("/")
        if not key or key in seen_heal:
            continue
        seen_heal.add(key)
        healed, heal_tries = fuzzy_heal_path_by_basename(
            src,
            upload_dir=ud,
            owner_id=owner_id,
            context_paths=context_paths,
            uploaded_entries=uploaded_entries,
            allow_dir=True,
        )
        attempted.extend(t for t in heal_tries if t not in attempted)
        if healed:
            logger.info("[AssetLocator] 模糊自愈: %s -> %s", src, healed)
            return healed

    if not must_exist and candidates:
        return candidates[0]

    msg = (
        f"无法在容器内定位文件：raw={raw_path!r}"
        + (f", asset_id={asset_id}" if asset_id else "")
        + (f"；已尝试 {attempted[:8]}" if attempted else "")
    )
    raise AssetLocatorError(msg, raw_path=raw_path or "", attempted=attempted, asset_id=asset_id)


def normalize_uploaded_files_list(
    files: List[Any],
    *,
    upload_dir: Optional[str] = None,
    owner_id: Optional[str] = None,
    db: Any = None,
    strict: bool = False,
) -> List[Dict[str, Any]]:
    """第一入口包装：资产 DB 桥接 + preemptive heal，再返回归一化列表。"""
    resolved, _warnings = bridge_uploaded_files_at_entry(
        files,
        upload_dir=upload_dir,
        owner_id=owner_id,
        db=db,
        strict=strict,
    )
    return resolved


def normalize_uploaded_file_entry(
    file_info: Dict[str, Any],
    *,
    upload_dir: Optional[str] = None,
    owner_id: Optional[str] = None,
    db: Any = None,
    strict: bool = False,
) -> Dict[str, Any]:
    """单条目归一化：走资产桥接 + preemptive heal。"""
    items, warnings = bridge_uploaded_files_at_entry(
        [file_info],
        upload_dir=upload_dir,
        owner_id=owner_id,
        db=db,
        strict=strict,
    )
    if not items:
        raise AssetLocatorError("无法归一化 uploaded_files 条目", raw_path=str(file_info))
    out = items[0]
    if warnings:
        out.setdefault("locator_warnings", []).extend(warnings)
    return out
