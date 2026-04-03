"""
管理员 API（统一挂载在 app.include_router(admin.router)，与 skills 路由解耦）

- GET  /api/admin/skills — 待审核技能（原在 skills 路由；skills 模块加载失败时此处仍可用）
- PUT  /api/admin/skills/{skill_id}/status
- POST /api/admin/skills/{skill_id}/review
- GET  /api/admin/users/pending
- POST /api/admin/users/{user_id}/approve | /reject
"""
from __future__ import annotations

import logging
from datetime import datetime
from typing import Any, List

from fastapi import APIRouter, Depends, HTTPException, status
from pydantic import BaseModel, ConfigDict, Field
from sqlalchemy.orm import Session

from gibh_agent.core.deps import get_current_admin_user
from gibh_agent.db.connection import get_db_session
from gibh_agent.db.models import Skill as SkillModel, User

logger = logging.getLogger(__name__)

router = APIRouter(prefix="/api", tags=["Admin"])


class SkillStatusUpdate(BaseModel):
    status: str = Field(..., pattern="^(approved|rejected)$")


class SkillReviewAction(BaseModel):
    action: str = Field(..., pattern="^(approve|reject)$")


class PendingUserOut(BaseModel):
    model_config = ConfigDict(from_attributes=True)

    id: int
    username: str
    email: str | None
    created_at: datetime | None


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
    rows = (
        db.query(SkillModel)
        .filter(SkillModel.status == "pending")
        .order_by(SkillModel.created_at.desc())
        .all()
    )
    return [
        {
            "id": r.id,
            "name": r.name,
            "description": r.description,
            "main_category": r.main_category,
            "sub_category": r.sub_category,
            "prompt_template": r.prompt_template,
            "author_id": r.author_id,
            "status": r.status,
            "created_at": r.created_at.isoformat() if r.created_at else None,
        }
        for r in rows
    ]


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
) -> dict[str, Any]:
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
) -> dict[str, Any]:
    user = db.query(User).filter(User.id == user_id).first()
    if not user:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="用户不存在")
    user.approval_status = "rejected"
    db.commit()
    logger.info("管理员拒绝用户: id=%s username=%s by=%s", user_id, user.username, _admin.username)
    return {"ok": True, "user_id": user_id, "approval_status": "rejected"}
