"""
Phase 4 - 侧栏数据读取：会话列表、消息历史、资产列表

- GET /api/sessions: 当前用户（owner_id）的历史会话，按 created_at 倒序
- GET /api/sessions/{session_id}/messages: 指定会话的消息列表（校验归属）
- GET /api/assets: 当前用户的数据资产列表
"""
from typing import List

from fastapi import APIRouter, Depends, HTTPException, status
from sqlalchemy.orm import Session

from gibh_agent.core.deps import get_current_owner_id
from gibh_agent.db.connection import get_db_session
from gibh_agent.db.models import Session as SessionModel, Message as MessageModel, Asset as AssetModel

router = APIRouter(tags=["User Data"])


@router.get("/api/sessions")
def list_sessions(
    owner_id: str = Depends(get_current_owner_id),
    db: Session = Depends(get_db_session),
) -> List[dict]:
    """按 owner_id 返回历史会话列表，按创建时间倒序。"""
    rows = (
        db.query(SessionModel)
        .filter(SessionModel.owner_id == owner_id)
        .order_by(SessionModel.created_at.desc())
        .all()
    )
    return [
        {
            "id": r.id,
            "owner_id": r.owner_id,
            "title": r.title,
            "created_at": r.created_at.isoformat() if r.created_at else None,
        }
        for r in rows
    ]


@router.get("/api/sessions/{session_id}/messages")
def list_messages(
    session_id: str,
    owner_id: str = Depends(get_current_owner_id),
    db: Session = Depends(get_db_session),
) -> List[dict]:
    """返回指定会话的消息列表；仅当会话属于当前 owner_id 时允许访问。"""
    session = db.query(SessionModel).filter(SessionModel.id == session_id).first()
    if not session:
        raise HTTPException(status_code=404, detail="会话不存在")
    if session.owner_id != owner_id:
        raise HTTPException(status_code=status.HTTP_403_FORBIDDEN, detail="无权访问该会话")
    rows = (
        db.query(MessageModel)
        .filter(MessageModel.session_id == session_id)
        .order_by(MessageModel.created_at.asc())
        .all()
    )
    return [
        {
            "id": r.id,
            "session_id": r.session_id,
            "role": r.role,
            "content": r.content,
            "created_at": r.created_at.isoformat() if r.created_at else None,
        }
        for r in rows
    ]


@router.get("/api/assets")
def list_assets(
    owner_id: str = Depends(get_current_owner_id),
    db: Session = Depends(get_db_session),
) -> List[dict]:
    """按 owner_id 返回数据资产列表。"""
    rows = (
        db.query(AssetModel)
        .filter(AssetModel.owner_id == owner_id)
        .order_by(AssetModel.created_at.desc())
        .all()
    )
    return [
        {
            "id": r.id,
            "owner_id": r.owner_id,
            "file_name": r.file_name,
            "file_path": r.file_path,
            "modality": r.modality,
            "created_at": r.created_at.isoformat() if r.created_at else None,
        }
        for r in rows
    ]
