"""
Phase 2 - FastAPI 鉴权依赖注入

get_current_user: 从 Bearer JWT 解码 username，查 User 表，无效或不存在则 401。
get_current_owner_id: Phase 4 身份解析，支持 JWT 或 X-Guest-UUID，返回 owner_id 字符串。
"""
from typing import Generator

from fastapi import Depends, HTTPException, Request, status
from fastapi.security import OAuth2PasswordBearer
from sqlalchemy.orm import Session

from gibh_agent.core.security import decode_access_token
from gibh_agent.db.connection import get_db_session
from gibh_agent.db.models import User

# Bearer 方案，tokenUrl 仅用于 OpenAPI 文档「Authorize」
oauth2_scheme = OAuth2PasswordBearer(tokenUrl="/api/auth/login", auto_error=True)


def get_current_owner_id(request: Request) -> str:
    """
    Phase 4 全局身份解析：优先 JWT 的 username，否则 X-Guest-UUID。
    两者都无或无效则 401。不查库，仅解析身份。
    Bearer null/格式错误/JWTError 等一律按无效 token 处理，抛 401 而非 500。
    """
    auth = request.headers.get("Authorization")
    if auth and auth.startswith("Bearer "):
        token = auth[7:].strip()
        if token and token.lower() != "null":
            try:
                payload = decode_access_token(token)
                if payload and payload.get("sub"):
                    return str(payload["sub"]).strip()
            except Exception:
                pass
    guest = request.headers.get("X-Guest-UUID")
    if guest and guest.strip():
        return guest.strip()
    raise HTTPException(
        status_code=status.HTTP_401_UNAUTHORIZED,
        detail="缺少身份标识：请提供 Authorization Bearer Token 或 X-Guest-UUID",
        headers={"WWW-Authenticate": "Bearer"},
    )


def get_current_user(
    token: str = Depends(oauth2_scheme),
    db: Session = Depends(get_db_session),
) -> User:
    """
    依赖：解码 JWT，按 username 查 User；无效或不存在则 401。
    仅当 MySQL 可用且 get_db_session 注入成功时调用。
    """
    payload = decode_access_token(token)
    if payload is None:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="无效或过期的 Token",
            headers={"WWW-Authenticate": "Bearer"},
        )
    username = payload.get("sub")
    if not username:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Token 缺少用户标识",
            headers={"WWW-Authenticate": "Bearer"},
        )
    user = db.query(User).filter(User.username == username).first()
    if user is None:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="用户不存在",
            headers={"WWW-Authenticate": "Bearer"},
        )
    return user
