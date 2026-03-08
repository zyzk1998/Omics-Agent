"""
Phase 2 - 鉴权与游客资产继承 API

- POST /register: 注册
- POST /login: 登录，返回 JWT
- POST /merge_guest_data: 已登录用户将 guest_uuid 下资产合并到当前用户（事务内 3 表 UPDATE）
"""
import logging
from datetime import timedelta

from fastapi import APIRouter, Depends, HTTPException, status
from fastapi.security import OAuth2PasswordRequestForm
from pydantic import BaseModel
from sqlalchemy.orm import Session

from gibh_agent.core.security import (
    ACCESS_TOKEN_EXPIRE_MINUTES,
    create_access_token,
    get_password_hash,
    verify_password,
)
from gibh_agent.core.deps import get_current_user
from gibh_agent.db.connection import get_db_session
from gibh_agent.db.models import User, Session as SessionModel, Asset, WorkflowTemplate

logger = logging.getLogger(__name__)

router = APIRouter(tags=["Authentication"])


# ----- 请求/响应模型 -----

class RegisterBody(BaseModel):
    username: str
    password: str


class MergeGuestBody(BaseModel):
    guest_uuid: str


# ----- 路由 -----

@router.post("/register")
def register(body: RegisterBody, db: Session = Depends(get_db_session)):
    """注册：检查用户名是否存在，哈希密码后写入 User 表。"""
    if not body.username or not body.password:
        raise HTTPException(status_code=400, detail="用户名和密码不能为空")
    existing = db.query(User).filter(User.username == body.username).first()
    if existing:
        raise HTTPException(status_code=400, detail="用户名已存在")
    user = User(
        username=body.username,
        hashed_password=get_password_hash(body.password),
        role="user",
    )
    db.add(user)
    db.commit()
    db.refresh(user)
    logger.info("用户注册: %s", body.username)
    return {"username": user.username, "message": "注册成功"}


@router.post("/login")
def login(
    form: OAuth2PasswordRequestForm = Depends(),
    db: Session = Depends(get_db_session),
):
    """登录：校验密码，成功后返回 access_token、token_type、username。"""
    user = db.query(User).filter(User.username == form.username).first()
    if not user or not verify_password(form.password, user.hashed_password):
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="用户名或密码错误",
            headers={"WWW-Authenticate": "Bearer"},
        )
    expires = timedelta(minutes=ACCESS_TOKEN_EXPIRE_MINUTES)
    token = create_access_token(data={"sub": user.username}, expires_delta=expires)
    return {
        "access_token": token,
        "token_type": "bearer",
        "username": user.username,
    }


@router.post("/merge_guest_data")
def merge_guest_data(
    body: MergeGuestBody,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db_session),
):
    """
    魔法继承：将 guest_uuid 下的会话/资产/工作流模板归属改为当前登录用户。
    在事务内执行 3 条 UPDATE（Session, Asset, WorkflowTemplate）；Message 通过 session 归属自然继承。
    """
    guest_uuid = (body.guest_uuid or "").strip()
    if not guest_uuid:
        raise HTTPException(status_code=400, detail="guest_uuid 不能为空")

    try:
        # 事务内 3 表 UPDATE
        count_sessions = (
            db.query(SessionModel)
            .filter(SessionModel.owner_id == guest_uuid)
            .update(
                {SessionModel.owner_id: current_user.username},
                synchronize_session=False,
            )
        )
        count_assets = (
            db.query(Asset)
            .filter(Asset.owner_id == guest_uuid)
            .update(
                {Asset.owner_id: current_user.username},
                synchronize_session=False,
            )
        )
        count_templates = (
            db.query(WorkflowTemplate)
            .filter(WorkflowTemplate.owner_id == guest_uuid)
            .update(
                {WorkflowTemplate.owner_id: current_user.username},
                synchronize_session=False,
            )
        )
        db.commit()
        logger.info(
            "merge_guest_data: guest_uuid=%s -> user=%s, sessions=%s, assets=%s, templates=%s",
            guest_uuid,
            current_user.username,
            count_sessions,
            count_assets,
            count_templates,
        )
        return {
            "message": "游客数据已合并到当前账户",
            "sessions_updated": count_sessions,
            "assets_updated": count_assets,
            "workflow_templates_updated": count_templates,
        }
    except Exception as e:
        db.rollback()
        logger.exception("merge_guest_data 失败: %s", e)
        raise HTTPException(status_code=500, detail="合并失败，请稍后重试")
