"""
Phase 2 - 鉴权与游客资产继承 API

- POST /register: 注册
- POST /login: 登录，返回 JWT
- POST /merge_guest_data: 已登录用户将 guest_uuid 下资产合并到当前用户（事务内 3 表 UPDATE）
- INITIAL_ADMINS 冷启动：环境变量逗号分隔 username，注册/登录时若在其中则强制 role=admin
"""
from __future__ import annotations

import logging
import os
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

# bcrypt 算法限制：输入不得超过 72 字节
BCRYPT_MAX_PASSWORD_BYTES = 72

# 管理员冷启动：环境变量 INITIAL_ADMINS 逗号分隔的 username 列表
INITIAL_ADMINS = frozenset(
    s.strip() for s in (os.environ.get("INITIAL_ADMINS") or "").split(",") if s.strip()
)


def _sync_admin_role(db: Session, user: User) -> None:
    """若 username 在 INITIAL_ADMINS 中，强制 admin + 审核通过（可登录）。"""
    if not user or user.username not in INITIAL_ADMINS:
        return
    changed = False
    if user.role != "admin":
        user.role = "admin"
        changed = True
    if getattr(user, "approval_status", None) != "approved":
        user.approval_status = "approved"
        changed = True
    if not changed:
        return
    try:
        db.commit()
        db.refresh(user)
        logger.info("INITIAL_ADMINS 冷启动: 已将用户 %s 设为 admin 且审核通过", user.username)
    except Exception as e:
        db.rollback()
        logger.warning("sync admin role 失败: %s", e)


def _approval_status_of(user: User) -> str | None:
    v = getattr(user, "approval_status", None)
    if v is None:
        return None
    s = str(v).strip().lower()
    return s if s else None


@router.post("/register")
def register(body: RegisterBody, db: Session = Depends(get_db_session)):
    """注册：检查用户名是否存在，哈希密码后写入 User 表。"""
    if not body.username or not body.password:
        raise HTTPException(status_code=400, detail="用户名和密码不能为空")
    plain = body.password if isinstance(body.password, str) else str(body.password)
    pwd_bytes = plain.encode("utf-8")
    if len(pwd_bytes) > BCRYPT_MAX_PASSWORD_BYTES:
        logger.warning(
            "register 收到超长 password: username=%r, password 类型=%s, 字节长=%d",
            body.username,
            type(body.password).__name__,
            len(pwd_bytes),
        )
        raise HTTPException(
            status_code=400,
            detail="密码长度不能超过 72 字节（当前 %d 字节），请检查输入或请求体是否被错误绑定" % len(pwd_bytes),
        )
    existing = db.query(User).filter(User.username == body.username).first()
    if existing:
        raise HTTPException(status_code=400, detail="用户名已存在")
    approval = "approved" if body.username in INITIAL_ADMINS else "pending"
    user = User(
        username=body.username,
        hashed_password=get_password_hash(plain),
        role="user",
        approval_status=approval,
    )
    db.add(user)
    db.commit()
    db.refresh(user)
    _sync_admin_role(db, user)
    db.refresh(user)
    logger.info("用户注册: %s approval=%s", body.username, user.approval_status)
    if (user.approval_status or "") == "pending":
        return {
            "username": user.username,
            "message": "注册已提交，请等待管理员审核",
            "approval_status": user.approval_status,
        }
    return {
        "username": user.username,
        "message": "注册成功（管理员账户已自动通过审核）",
        "role": user.role,
        "approval_status": user.approval_status,
    }


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
    _sync_admin_role(db, user)
    db.refresh(user)
    st = _approval_status_of(user)
    if st == "pending":
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="ACCOUNT_PENDING",
        )
    if st == "rejected":
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="ACCOUNT_REJECTED",
        )
    expires = timedelta(minutes=ACCESS_TOKEN_EXPIRE_MINUTES)
    token = create_access_token(data={"sub": user.username}, expires_delta=expires)
    return {
        "access_token": token,
        "token_type": "bearer",
        "username": user.username,
        "role": user.role,
    }


@router.get("/me")
def me(current_user: User = Depends(get_current_user)):
    """当前用户信息（含 role），供前端显隐管理员入口。"""
    return {"username": current_user.username, "role": current_user.role or "user"}


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
