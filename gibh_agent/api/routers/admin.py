"""
管理员 API（统一挂载在 app.include_router(admin.router)，与 skills 路由解耦）

- GET  /api/admin/skills — 待审核技能（UGC + dynamic_skill_plugins pending/accepted）
- POST /api/admin/plugins/{id}/review — 状态机 accept | publish | reject
"""
from __future__ import annotations

import logging
from datetime import datetime
from typing import Any, Dict, List, Optional

from fastapi import APIRouter, Depends, HTTPException, status
from pydantic import BaseModel, ConfigDict, Field
from sqlalchemy.orm import Session

from gibh_agent.core.deps import get_current_admin_user
from gibh_agent.core.user_notifications import create_user_notification
from gibh_agent.db.connection import get_db_session
from gibh_agent.db.models import Skill as SkillModel, User
from gibh_agent.plugin_system.registry import DynamicSkillPlugin

logger = logging.getLogger(__name__)

router = APIRouter(prefix="/api", tags=["Admin"])

PLUGIN_REVIEW_ACTIVE_STATUSES = ("pending", "accepted")
PLUGIN_PLAZA_VISIBLE_STATUSES = ("approved", "published")


class SkillStatusUpdate(BaseModel):
    status: str = Field(..., pattern="^(approved|rejected)$")


class SkillReviewAction(BaseModel):
    action: str = Field(..., pattern="^(approve|reject)$")


class PluginReviewAction(BaseModel):
    """用户安装技能审核：pending → accepted → published / rejected。"""

    action: str = Field(..., pattern="^(accept|publish|reject)$")


class PendingUserOut(BaseModel):
    model_config = ConfigDict(from_attributes=True)

    id: int
    username: str
    email: Optional[str] = None
    created_at: Optional[datetime] = None


def _plugin_display_name(row: DynamicSkillPlugin) -> str:
    return (row.display_name or row.name or "未命名技能").strip()


def _notify_plugin_review(
    db: Session,
    row: DynamicSkillPlugin,
    *,
    ntype: str,
    title: str,
    content: str,
) -> None:
    try:
        create_user_notification(
            db,
            user_id=row.author_id,
            ntype=ntype,
            title=title,
            content=content,
            commit=False,
        )
    except Exception as e:
        logger.warning("写入用户通知失败 plugin_id=%s: %s", row.id, e)


def _apply_skill_review_status(
    db: Session, skill_id: int, new_status: str, admin_username: str
) -> SkillModel:
    skill = db.query(SkillModel).filter(SkillModel.id == skill_id).first()
    if not skill:
        raise HTTPException(status_code=404, detail="技能不存在")
    try:
        skill.status = new_status
        db.commit()
        db.refresh(skill)
        logger.info(
            "管理员审批技能: id=%s status=%s by=%s", skill_id, skill.status, admin_username
        )
        return skill
    except Exception as e:
        db.rollback()
        logger.exception("更新技能状态失败: %s", e)
        raise HTTPException(status_code=500, detail="更新失败")


@router.get("/admin/skills", response_model=List[dict])
def list_pending_skills(
    current_user: User = Depends(get_current_admin_user),
    db: Session = Depends(get_db_session),
):
    """待审核技能：UGC skills 表 pending + 用户上传 dynamic_skill_plugins（pending/accepted）。"""
    ugc_rows = (
        db.query(SkillModel)
        .filter(SkillModel.status == "pending")
        .order_by(SkillModel.created_at.desc())
        .all()
    )
    plugin_rows = (
        db.query(DynamicSkillPlugin)
        .filter(DynamicSkillPlugin.status.in_(PLUGIN_REVIEW_ACTIVE_STATUSES))
        .order_by(DynamicSkillPlugin.created_at.desc())
        .all()
    )
    items: List[dict] = []
    for r in ugc_rows:
        items.append(
            {
                "id": r.id,
                "source_type": "ugc",
                "name": r.name,
                "description": r.description,
                "main_category": r.main_category,
                "sub_category": r.sub_category,
                "prompt_template": r.prompt_template,
                "author_id": r.author_id,
                "status": r.status,
                "created_at": r.created_at.isoformat() if r.created_at else None,
            }
        )
    for r in plugin_rows:
        review_meta = {}
        if isinstance(r.parameters_schema, dict):
            review_meta = r.parameters_schema.get("_review_meta") or {}
        driver = review_meta.get("driver_type") or r.skill_type or "prompt"
        items.append(
            {
                "id": r.id,
                "source_type": "plugin",
                "name": r.display_name or r.name,
                "description": r.description,
                "main_category": "用户上传",
                "sub_category": driver,
                "prompt_template": None,
                "author_id": r.author_id,
                "status": r.status,
                "created_at": r.created_at.isoformat() if r.created_at else None,
                "extract_dir": r.extract_dir,
                "skill_type": r.skill_type,
                "submission_uid": review_meta.get("submission_uid"),
            }
        )
    items.sort(key=lambda x: x.get("created_at") or "", reverse=True)
    return items


@router.put("/admin/skills/{skill_id}/status")
def update_skill_status(
    skill_id: int,
    body: SkillStatusUpdate,
    current_user: User = Depends(get_current_admin_user),
    db: Session = Depends(get_db_session),
):
    st = body.status.strip().lower()
    skill = _apply_skill_review_status(db, skill_id, st, current_user.username)
    return {"skill_id": skill.id, "status": skill.status, "message": "已更新"}


@router.post("/admin/skills/{skill_id}/review")
def review_skill_by_action(
    skill_id: int,
    body: SkillReviewAction,
    current_user: User = Depends(get_current_admin_user),
    db: Session = Depends(get_db_session),
):
    new_status = "approved" if body.action == "approve" else "rejected"
    skill = _apply_skill_review_status(db, skill_id, new_status, current_user.username)
    return {"skill_id": skill.id, "status": skill.status, "message": "已更新"}


def _transition_plugin_status(
    row: DynamicSkillPlugin, action: str
) -> str:
    cur = (row.status or "pending").strip().lower()
    act = (action or "").strip().lower()
    if act == "accept":
        if cur != "pending":
            raise HTTPException(status_code=400, detail=f"当前状态 {cur!r} 不可受理，仅 pending 可受理")
        return "accepted"
    if act == "publish":
        if cur not in ("accepted", "approved"):
            raise HTTPException(
                status_code=400,
                detail=f"当前状态 {cur!r} 不可上架，须先受理（accepted）",
            )
        return "published"
    if act == "reject":
        if cur in ("rejected", "published"):
            raise HTTPException(status_code=400, detail=f"当前状态 {cur!r} 不可驳回")
        return "rejected"
    raise HTTPException(status_code=400, detail="action 必须是 accept、publish 或 reject")


@router.post("/admin/plugins/{plugin_id}/review")
def review_plugin_by_action(
    plugin_id: int,
    body: PluginReviewAction,
    current_user: User = Depends(get_current_admin_user),
    db: Session = Depends(get_db_session),
):
    """用户「安装技能」审核状态机 + 用户通知。"""
    row = db.query(DynamicSkillPlugin).filter(DynamicSkillPlugin.id == plugin_id).first()
    if not row:
        raise HTTPException(status_code=404, detail="动态插件不存在")

    new_status = _transition_plugin_status(row, body.action)
    display = _plugin_display_name(row)

    try:
        row.status = new_status
        if body.action == "accept":
            _notify_plugin_review(
                db,
                row,
                ntype="skill_review",
                title="技能已受理",
                content=(
                    f"🎉 您的技能「{display}」已被平台受理！"
                    "感谢您的贡献，我们的开发团队正在进行底层的安全审查与融合。"
                    "请关注后续上架进度。"
                ),
            )
        elif body.action == "publish":
            _notify_plugin_review(
                db,
                row,
                ntype="skill_published",
                title="技能已上架",
                content=(
                    f"🚀 好消息！您提交的技能「{display}」现已在 Omics Agent 技能广场正式上线！"
                    "感谢您为生物信息生态做出的卓越贡献！"
                ),
            )
        elif body.action == "reject":
            _notify_plugin_review(
                db,
                row,
                ntype="skill_review",
                title="技能未通过审核",
                content=(
                    f"很遗憾，您提交的技能「{display}」暂未通过平台审核。"
                    "可能原因包括：规范不符、安全风险或重复条目。"
                    "您可对照技能扩展规范修订后重新提交。"
                ),
            )
        db.commit()
        db.refresh(row)
        logger.info(
            "管理员审批动态插件: id=%s action=%s status=%s by=%s",
            plugin_id,
            body.action,
            row.status,
            current_user.username,
        )
    except HTTPException:
        raise
    except Exception as e:
        db.rollback()
        logger.exception("动态插件审核失败: %s", e)
        raise HTTPException(status_code=500, detail="更新失败") from e

    return {
        "plugin_id": row.id,
        "status": row.status,
        "message": "已更新",
        "name": display,
        "action": body.action,
    }


@router.get("/admin/users/pending", response_model=List[PendingUserOut])
def list_pending_users(
    db: Session = Depends(get_db_session),
    _admin: User = Depends(get_current_admin_user),
):
    rows = (
        db.query(User)
        .filter(User.approval_status == "pending")
        .order_by(User.created_at.asc())
        .all()
    )
    return rows


@router.post("/admin/users/{user_id}/approve")
def approve_user(
    user_id: int,
    db: Session = Depends(get_db_session),
    _admin: User = Depends(get_current_admin_user),
) -> Dict[str, Any]:
    user = db.query(User).filter(User.id == user_id).first()
    if not user:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="用户不存在")
    user.approval_status = "approved"
    db.commit()
    logger.info("管理员审核通过用户: id=%s username=%s by=%s", user_id, user.username, _admin.username)
    return {"ok": True, "user_id": user_id, "approval_status": "approved"}


@router.post("/admin/users/{user_id}/reject")
def reject_user(
    user_id: int,
    db: Session = Depends(get_db_session),
    _admin: User = Depends(get_current_admin_user),
) -> Dict[str, Any]:
    user = db.query(User).filter(User.id == user_id).first()
    if not user:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="用户不存在")
    user.approval_status = "rejected"
    db.commit()
    logger.info("管理员拒绝用户: id=%s username=%s by=%s", user_id, user.username, _admin.username)
    return {"ok": True, "user_id": user_id, "approval_status": "rejected"}
