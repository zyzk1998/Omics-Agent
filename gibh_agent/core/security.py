"""
Phase 2 - 鉴权安全工具：密码哈希与 JWT 签发

- passlib (bcrypt) 密码哈希
- python-jose JWT 创建与解码
- SECRET_KEY / ALGORITHM 从环境变量读取，默认 ALGORITHM=HS256
"""
import os
from datetime import datetime, timedelta
from typing import Any, Optional

# 从环境变量读取，无则使用默认（生产环境务必设置 SECRET_KEY）
SECRET_KEY = os.environ.get("SECRET_KEY", "change-me-in-production-gibh-agent-v2")
ALGORITHM = os.environ.get("JWT_ALGORITHM", "HS256")
ACCESS_TOKEN_EXPIRE_MINUTES = int(os.environ.get("ACCESS_TOKEN_EXPIRE_MINUTES", "60"))

# passlib 上下文：bcrypt
_pwd_context: Optional[Any] = None


def _get_pwd_context():
    global _pwd_context
    if _pwd_context is None:
        from passlib.context import CryptContext
        _pwd_context = CryptContext(schemes=["bcrypt"], deprecated="auto")
    return _pwd_context


def get_password_hash(password: str) -> str:
    """对明文密码做 bcrypt 哈希。Docker 未装 bcrypt C 扩展时会抛异常，此处统一包装便于排查。"""
    try:
        return _get_pwd_context().hash(password)
    except Exception as e:
        raise RuntimeError(
            "密码哈希不可用，请检查 bcrypt 安装（如 Docker 中需安装 build 依赖）: %s" % e
        ) from e


def verify_password(plain_password: str, hashed_password: str) -> bool:
    """校验明文密码与哈希是否一致。"""
    try:
        return _get_pwd_context().verify(plain_password, hashed_password)
    except Exception as e:
        raise RuntimeError(
            "密码校验不可用，请检查 bcrypt 安装: %s" % e
        ) from e


def create_access_token(
    data: dict,
    expires_delta: Optional[timedelta] = None,
) -> str:
    """
    签发 JWT。
    data 通常包含 sub (username) 等声明；过期时间由 expires_delta 或默认配置决定。
    """
    from jose import jwt

    to_encode = data.copy()
    if expires_delta is not None:
        expire = datetime.utcnow() + expires_delta
    else:
        expire = datetime.utcnow() + timedelta(minutes=ACCESS_TOKEN_EXPIRE_MINUTES)
    to_encode.update({"exp": expire})
    return jwt.encode(to_encode, SECRET_KEY, algorithm=ALGORITHM)


def decode_access_token(token: str) -> Optional[dict]:
    """解码 JWT，失败返回 None。"""
    from jose import jwt
    from jose.exceptions import JWTError

    try:
        payload = jwt.decode(token, SECRET_KEY, algorithms=[ALGORITHM])
        return payload
    except JWTError:
        return None
