"""
UGC 技能上传与管理员审核 API

路由规范：APIRouter(prefix="/api")，路径不再带 /api，避免嵌套重复。
- GET    /api/skills: 公开分页查询（正交过滤 main_cat, sub_cat），仅返回 approved；saved_only=True 时返回当前用户收藏
- POST   /api/skills: 用户上传；普通用户 pending，管理员 approved
- POST   /api/skills/{skill_id}/bookmark: 收藏技能（防重复）
- DELETE /api/skills/{skill_id}/bookmark: 取消收藏
- GET    /api/admin/skills: 管理员拉取待审核列表（严格 role 校验）
- PUT    /api/admin/skills/{skill_id}/status: 管理员审批（body: {status}）
- POST   /api/admin/skills/{skill_id}/review: 管理员审批（body: {action: approve|reject}）
"""
import logging
from typing import List, Optional

from fastapi import APIRouter, Depends, HTTPException, Query, Request, status
from pydantic import BaseModel, Field
from sqlalchemy.orm import Session

from gibh_agent.core.deps import get_current_user, get_current_admin_user, get_current_owner_id
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


class SkillStatusUpdate(BaseModel):
    """管理员审批请求体。"""
    status: str = Field(..., pattern="^(approved|rejected)$")


class SkillReviewAction(BaseModel):
    """管理员审批（POST 别名）：action → status。"""
    action: str = Field(..., pattern="^(approve|reject)$")


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
    """公开分页查询：仅 status=approved。saved_only=True 时需鉴权并只返回当前用户收藏。若当前无技能则自动补种 7 大核心组学后重查。"""
    owner_id: Optional[str] = None
    try:
        owner_id = get_current_owner_id(request)
    except HTTPException:
        pass
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
    rows = q.order_by(SkillModel.created_at.desc()).offset(offset).limit(size).all()
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


def _apply_skill_review_status(db: Session, skill_id: int, new_status: str, admin_username: str) -> SkillModel:
    skill = db.query(SkillModel).filter(SkillModel.id == skill_id).first()
    if not skill:
        raise HTTPException(status_code=404, detail="技能不存在")
    try:
        skill.status = new_status
        db.commit()
        db.refresh(skill)
        logger.info("管理员审批技能: id=%s status=%s by=%s", skill_id, skill.status, admin_username)
        return skill
    except Exception as e:
        db.rollback()
        logger.exception("更新技能状态失败: %s", e)
        raise HTTPException(status_code=500, detail="更新失败")


@router.put("/admin/skills/{skill_id}/status")
def update_skill_status(
    skill_id: int,
    body: SkillStatusUpdate,
    current_user: User = Depends(get_current_admin_user),
    db: Session = Depends(get_db_session),
):
    """管理员审批；必须 role == admin。"""
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
    """管理员审批（POST）：action 为 approve / reject。"""
    new_status = "approved" if body.action == "approve" else "rejected"
    skill = _apply_skill_review_status(db, skill_id, new_status, current_user.username)
    return {"skill_id": skill.id, "status": skill.status, "message": "已更新"}
