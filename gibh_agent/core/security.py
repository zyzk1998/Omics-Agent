"""
Phase 2 - 鉴权安全工具：密码哈希与 JWT 签发

- pwdlib (bcrypt) 密码哈希，替代 passlib，兼容最新 bcrypt
- python-jose JWT 创建与解码
- SECRET_KEY / ALGORITHM 从环境变量读取，默认 ALGORITHM=HS256
"""
import os
from datetime import datetime, timedelta
from typing import Optional

from pwdlib import PasswordHash
from pwdlib.hashers.bcrypt import BcryptHasher

def _load_secret_key() -> str:
    """从环境或宿主机持久化来源读取 SECRET_KEY，不依赖容器生命周期，重启后不变。"""
    key = (os.environ.get("SECRET_KEY") or "").strip()
    if key:
        return key
    path = (os.environ.get("SECRET_KEY_FILE") or "").strip()
    if path and os.path.isfile(path):
        try:
            with open(path, "r", encoding="utf-8") as f:
                key = (f.read() or "").strip()
            if key:
                return key
        except Exception:
            pass
    return "09d25e094faa6ca2556c818166b7a9563b93f7099f6f0f4caa6cf63b88e8d3e7"


SECRET_KEY = _load_secret_key()
ALGORITHM = os.environ.get("JWT_ALGORITHM", "HS256")
# 默认 7 天（10080 分钟）；可通过环境变量 ACCESS_TOKEN_EXPIRE_MINUTES 覆盖
ACCESS_TOKEN_EXPIRE_MINUTES = int(os.environ.get("ACCESS_TOKEN_EXPIRE_MINUTES", "10080"))

# 仅使用 bcrypt，不依赖 argon2（避免 HasherNotAvailable）
pwd_context = PasswordHash((BcryptHasher(),))


def get_password_hash(password: str) -> str:
    """对明文密码做 bcrypt 哈希。"""
    return pwd_context.hash(password)


def verify_password(plain_password: str, hashed_password: str) -> bool:
    """校验明文密码与哈希是否一致。"""
    if not hashed_password:
        return False
    try:
        return pwd_context.verify(plain_password, hashed_password)
    except Exception:
        return False


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
