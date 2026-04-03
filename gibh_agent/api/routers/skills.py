"""
UGC 技能上传与管理员审核 API

路由规范：APIRouter(prefix="/api")，路径不再带 /api，避免嵌套重复。
- GET    /api/skills: 公开分页查询（正交过滤 main_cat, sub_cat），仅返回 approved；saved_only=True 时返回当前用户收藏
- POST   /api/skills: 用户上传；普通用户 pending，管理员 approved
- POST   /api/skills/{skill_id}/bookmark: 收藏技能（防重复）
- DELETE /api/skills/{skill_id}/bookmark: 取消收藏
- 管理员审核技能接口已迁至 gibh_agent.api.routers.admin（与 /api/admin/users/* 同路由挂载）
"""
import logging
from typing import Optional

from fastapi import APIRouter, Depends, HTTPException, Query, Request, status
from pydantic import BaseModel, Field
from sqlalchemy import case, or_
from sqlalchemy.orm import Session

from gibh_agent.core.deps import (
    get_current_user,
    get_current_owner_id,
    get_optional_owner_id,
)
from gibh_agent.db.connection import get_db_session
from gibh_agent.db.models import User, Skill as SkillModel, UserSavedSkill

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


@router.get("/skills")
def list_skills_public(
    request: Request,
    main_cat: Optional[str] = Query(None, description="大类筛选"),
    sub_cat: Optional[str] = Query(None, description="小类筛选"),
    saved_only: bool = Query(False, description="仅返回当前用户收藏的技能"),
    page: int = Query(1, ge=1, description="页码"),
    size: int = Query(12, ge=1, le=50, description="每页条数"),
    db: Session = Depends(get_db_session),
):
    """公开分页查询：仅 status=approved；可无 Token/X-Guest（橱窗）。saved_only=True 时需 owner_id。若当前无技能则自动补种后重查。"""
    owner_id: Optional[str] = get_optional_owner_id(request)
    if saved_only and not owner_id:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="查看「我的」收藏需要登录或提供身份",
            headers={"WWW-Authenticate": "Bearer"},
        )

    def _query():
        q = db.query(SkillModel).filter(SkillModel.status == "approved")
        if saved_only and owner_id:
            q = q.join(
                UserSavedSkill,
                (UserSavedSkill.skill_id == SkillModel.id) & (UserSavedSkill.owner_id == owner_id),
            )
        if main_cat and (mc := main_cat.strip()):
            q = q.filter(SkillModel.main_category == mc)
        if sub_cat and (sc := sub_cat.strip()):
            q = q.filter(SkillModel.sub_category == sc)
        return q

    q = _query()
    total = q.count()
    if total == 0 and not saved_only:
        try:
            from gibh_agent.db.seed_skills import run_upsert_system_skills
            run_upsert_system_skills(db)
            logger.info("[Skills] 已按需幂等补种核心组学 + 生物医药 + 化学技能")
        except Exception as e:
            db.rollback()
            logger.warning("[Skills] 按需补种失败: %s", e)
        q = _query()
        total = q.count()
    offset = (page - 1) * size
    # 快车道暗号置顶：全局顺序优先于分页，故在 ORDER BY 中体现（等价于全量按暗号+时间排序后 slice）
    _pt = SkillModel.prompt_template
    _fast_lane_first = case(
        (or_(_pt.contains("[Skill_Route:"), _pt.contains("[Omics_Route:")), 0),
        else_=1,
    )
    rows = (
        q.order_by(_fast_lane_first.asc(), SkillModel.created_at.desc())
        .offset(offset)
        .limit(size)
        .all()
    )
    saved_ids: set = set()
    if owner_id:
        saved_rows = (
            db.query(UserSavedSkill.skill_id)
            .filter(UserSavedSkill.owner_id == owner_id)
            .all()
        )
        saved_ids = {r[0] for r in saved_rows}
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
            "saved": r.id in saved_ids,
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
    """用户上传技能；管理员直接 approved，普通用户 pending。"""
    name = (body.name or "").strip()
    if not name:
        raise HTTPException(status_code=400, detail="技能名称不能为空")
    is_admin = (current_user.role or "").strip().lower() == "admin"
    initial_status = "approved" if is_admin else "pending"
    try:
        skill = SkillModel(
            name=name,
            description=(body.description or "").strip() or None,
            main_category=(body.main_category or "").strip() or None,
            sub_category=(body.sub_category or "").strip() or None,
            prompt_template=(body.prompt_template or "").strip() or None,
            author_id=current_user.username,
            status=initial_status,
        )
        db.add(skill)
        db.commit()
        db.refresh(skill)
        logger.info(
            "UGC 技能提交: id=%s name=%s author=%s status=%s",
            skill.id,
            skill.name,
            skill.author_id,
            initial_status,
        )
        msg = "已发布（管理员上传直接通过）" if is_admin else "提交成功，等待管理员审核"
        return {"skill_id": skill.id, "message": msg, "status": initial_status}
    except Exception as e:
        db.rollback()
        logger.exception("创建技能失败: %s", e)
        raise HTTPException(status_code=500, detail="提交失败")


@router.post("/skills/{skill_id}/bookmark")
def bookmark_skill(
    skill_id: int,
    owner_id: str = Depends(get_current_owner_id),
    db: Session = Depends(get_db_session),
):
    """收藏技能：防重复插入。"""
    skill = db.query(SkillModel).filter(SkillModel.id == skill_id, SkillModel.status == "approved").first()
    if not skill:
        raise HTTPException(status_code=404, detail="技能不存在或未通过审核")
    existing = (
        db.query(UserSavedSkill)
        .filter(UserSavedSkill.owner_id == owner_id, UserSavedSkill.skill_id == skill_id)
        .first()
    )
    if existing:
        return {"status": "success", "skill_id": skill_id, "message": "已在收藏中", "saved": True}
    try:
        db.add(UserSavedSkill(owner_id=owner_id, skill_id=skill_id))
        db.commit()
        logger.info("用户收藏技能: owner=%s skill_id=%s", owner_id, skill_id)
        return {"status": "success", "skill_id": skill_id, "message": "已添加到我的工具", "saved": True}
    except Exception as e:
        db.rollback()
        logger.exception("收藏技能失败: %s", e)
        raise HTTPException(status_code=500, detail="收藏失败")


@router.delete("/skills/{skill_id}/bookmark")
def unbookmark_skill(
    skill_id: int,
    owner_id: str = Depends(get_current_owner_id),
    db: Session = Depends(get_db_session),
):
    """取消收藏。"""
    deleted = (
        db.query(UserSavedSkill)
        .filter(UserSavedSkill.owner_id == owner_id, UserSavedSkill.skill_id == skill_id)
        .delete()
    )
    db.commit()
    if deleted:
        logger.info("用户取消收藏技能: owner=%s skill_id=%s", owner_id, skill_id)
    return {"status": "success", "skill_id": skill_id, "message": "已取消收藏", "saved": False}
