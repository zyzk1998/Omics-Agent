"""
UGC 技能上传与管理员审核 API

路由规范：APIRouter(prefix="/api")，路径不再带 /api，避免嵌套重复。
- GET  /api/skills: 公开分页查询（正交过滤 main_cat, sub_cat），仅返回 approved
- POST /api/skills: 用户上传技能，status 强制 pending
- GET  /api/admin/skills: 管理员拉取待审核列表（严格 role 校验）
- PUT  /api/admin/skills/{skill_id}/status: 管理员审批
"""
import logging
from typing import List, Optional

from fastapi import APIRouter, Depends, HTTPException, Query, status
from pydantic import BaseModel, Field
from sqlalchemy.orm import Session

from gibh_agent.core.deps import get_current_user, get_current_admin_user
from gibh_agent.db.connection import get_db_session
from gibh_agent.db.models import User, Skill as SkillModel

logger = logging.getLogger(__name__)

# 统一 prefix="/api"，路径写相对路径，最终 URL 与前端 fetch 一致
router = APIRouter(prefix="/api", tags=["Skills"])

# Pydantic 强校验：限制长度与类型，防御注入与超大 payload
MAX_NAME = 256
MAX_MAIN_CATEGORY = 128
MAX_SUB_CATEGORY = 128
MAX_DESCRIPTION = 10000
MAX_PROMPT = 50000


class SkillCreate(BaseModel):
    """用户上传技能表单（强校验，防 SQL 注入 / XSS 通过长度与类型约束）。"""
    name: str = Field(..., min_length=1, max_length=MAX_NAME)
    description: str | None = Field(None, max_length=MAX_DESCRIPTION)
    main_category: str | None = Field(None, max_length=MAX_MAIN_CATEGORY)
    sub_category: str | None = Field(None, max_length=MAX_SUB_CATEGORY)
    prompt_template: str | None = Field(None, max_length=MAX_PROMPT)


class SkillStatusUpdate(BaseModel):
    """管理员审批请求体。"""
    status: str = Field(..., pattern="^(approved|rejected)$")


@router.get("/skills")
def list_skills_public(
    main_cat: Optional[str] = Query(None, description="大类筛选"),
    sub_cat: Optional[str] = Query(None, description="小类筛选"),
    page: int = Query(1, ge=1, description="页码"),
    size: int = Query(12, ge=1, le=50, description="每页条数"),
    db: Session = Depends(get_db_session),
):
    """公开分页查询：仅 status=approved，正交过滤 main_category + sub_category。"""
    q = db.query(SkillModel).filter(SkillModel.status == "approved")
    if main_cat and (main_cat := main_cat.strip()):
        q = q.filter(SkillModel.main_category == main_cat)
    if sub_cat and (sub_cat := sub_cat.strip()):
        q = q.filter(SkillModel.sub_category == sub_cat)
    total = q.count()
    offset = (page - 1) * size
    rows = q.order_by(SkillModel.created_at.desc()).offset(offset).limit(size).all()
    items = [
        {
            "id": r.id,
            "name": r.name,
            "description": r.description,
            "main_category": r.main_category,
            "sub_category": r.sub_category,
            "prompt_template": r.prompt_template,
            "author_id": r.author_id,
            "created_at": r.created_at.isoformat() if r.created_at else None,
        }
        for r in rows
    ]
    return {"items": items, "total": total}


@router.post("/skills")
def create_skill(
    body: SkillCreate,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db_session),
):
    """用户上传技能；status 强制为 pending，author_id 为当前用户。"""
    name = (body.name or "").strip()
    if not name:
        raise HTTPException(status_code=400, detail="技能名称不能为空")
    try:
        skill = SkillModel(
            name=name,
            description=(body.description or "").strip() or None,
            main_category=(body.main_category or "").strip() or None,
            sub_category=(body.sub_category or "").strip() or None,
            prompt_template=(body.prompt_template or "").strip() or None,
            author_id=current_user.username,
            status="pending",
        )
        db.add(skill)
        db.commit()
        db.refresh(skill)
        logger.info("UGC 技能提交: id=%s name=%s author=%s", skill.id, skill.name, skill.author_id)
        return {"skill_id": skill.id, "message": "提交成功，等待管理员审核", "status": "pending"}
    except Exception as e:
        db.rollback()
        logger.exception("创建技能失败: %s", e)
        raise HTTPException(status_code=500, detail="提交失败")


@router.get("/admin/skills", response_model=List[dict])
def list_pending_skills(
    current_user: User = Depends(get_current_admin_user),
    db: Session = Depends(get_db_session),
):
    """管理员拉取待审核技能；必须 role == admin。"""
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
    """管理员审批；必须 role == admin。"""
    skill = db.query(SkillModel).filter(SkillModel.id == skill_id).first()
    if not skill:
        raise HTTPException(status_code=404, detail="技能不存在")
    try:
        skill.status = body.status.strip().lower()
        db.commit()
        db.refresh(skill)
        logger.info("管理员审批技能: id=%s status=%s by=%s", skill_id, skill.status, current_user.username)
        return {"skill_id": skill.id, "status": skill.status, "message": "已更新"}
    except Exception as e:
        db.rollback()
        logger.exception("更新技能状态失败: %s", e)
        raise HTTPException(status_code=500, detail="更新失败")
