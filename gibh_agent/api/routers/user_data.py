"""
Phase 4 - 侧栏数据读取：会话列表、消息历史、资产列表、工作流收藏

- GET /api/sessions: 当前用户（owner_id）的历史会话，按 created_at 倒序
- GET /api/sessions/{session_id}/messages: 指定会话的消息列表（校验归属）
- GET /api/assets: 当前用户的数据资产列表
- DELETE /api/sessions/{session_id}: 删除会话及该会话下所有消息（校验 owner_id）
- DELETE /api/assets/{asset_id}: 删除资产记录（校验 owner_id，可选删除物理文件）
- GET /api/workflow_templates: 当前用户的工作流收藏列表
- POST /api/workflow_templates: 新建工作流收藏（name + config_json）
- DELETE /api/workflow_templates/{template_id}: 删除工作流收藏（校验 owner_id）
"""
import os
import logging
from typing import List, Optional

from fastapi import APIRouter, Depends, HTTPException, status
from pydantic import BaseModel
from sqlalchemy.orm import Session

from gibh_agent.core.deps import get_current_owner_id
from gibh_agent.db.connection import get_db_session
from gibh_agent.db.models import (
    Session as SessionModel,
    Message as MessageModel,
    Asset as AssetModel,
    WorkflowTemplate as WorkflowTemplateModel,
)


class WorkflowTemplateCreate(BaseModel):
    """POST /api/workflow_templates 请求体"""
    name: str
    config_json: dict = {}


class SessionRenameBody(BaseModel):
    """PUT /api/sessions/{session_id} 请求体"""
    title: str


class AssetRenameBody(BaseModel):
    """PUT /api/assets/{asset_id} 请求体"""
    file_name: str


class WorkflowTemplateRenameBody(BaseModel):
    """PUT /api/workflow_templates/{template_id} 请求体"""
    name: str

logger = logging.getLogger(__name__)

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
    logger.info("🔍 [DB] 查询到 Owner: %s 共有 %s 条历史会话", owner_id, len(rows))
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


@router.delete("/api/messages/{message_id}")
def delete_message(
    message_id: int,
    owner_id: str = Depends(get_current_owner_id),
    db: Session = Depends(get_db_session),
) -> dict:
    """删除单条消息；仅当该消息所属会话属于当前 owner_id 时允许。"""
    msg = db.query(MessageModel).filter(MessageModel.id == message_id).first()
    if not msg:
        raise HTTPException(status_code=404, detail="消息不存在")
    session = db.query(SessionModel).filter(SessionModel.id == msg.session_id).first()
    if not session or session.owner_id != owner_id:
        raise HTTPException(status_code=status.HTTP_403_FORBIDDEN, detail="无权删除该消息")
    try:
        db.delete(msg)
        db.commit()
        return {"status": "success"}
    except Exception as e:
        db.rollback()
        logger.exception("删除消息失败: %s", e)
        raise HTTPException(status_code=500, detail="删除失败")


def _infer_modality_from_file_name(file_name: str) -> Optional[str]:
    """根据文件名推断组学类型，与前端 FIXED_MODALITY_LABELS 及 orchestrator 一致。"""
    if not file_name or not isinstance(file_name, str):
        return None
    name = file_name.strip().lower()
    # 10x 转录组三件套
    if name in ("barcodes.tsv", "features.tsv", "matrix.mtx", "genes.tsv"):
        return "rna"
    if name.endswith(".h5ad"):
        return "rna"
    if name.endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz")):
        return "rna"
    # 代谢组
    if name.endswith(".csv"):
        return "metabolomics"
    # 影像组
    if name.endswith((".nii", ".nii.gz", ".dcm")):
        return "radiomics"
    # 空间（常见后缀，保守归类）
    if name.endswith(".h5") or "spatial" in name or "visium" in name:
        return "spatial"
    return None


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


@router.post("/api/assets/reclassify")
def reclassify_assets(
    owner_id: str = Depends(get_current_owner_id),
    db: Session = Depends(get_db_session),
) -> dict:
    """对当前用户未分类的资产按文件名推断 modality，并刷新列表。"""
    rows = (
        db.query(AssetModel)
        .filter(AssetModel.owner_id == owner_id, AssetModel.modality.is_(None))
        .all()
    )
    updated = 0
    for r in rows:
        modality = _infer_modality_from_file_name(r.file_name)
        if modality:
            r.modality = modality
            updated += 1
    if updated:
        try:
            db.commit()
            logger.info("🔁 [Assets] 重新分类: owner=%s, 更新 %s 条", owner_id, updated)
        except Exception as e:
            db.rollback()
            logger.exception("重新分类提交失败: %s", e)
            raise HTTPException(status_code=500, detail="重新分类保存失败")
    return {"status": "success", "updated": updated}


@router.delete("/api/sessions/{session_id}")
def delete_session(
    session_id: str,
    owner_id: str = Depends(get_current_owner_id),
    db: Session = Depends(get_db_session),
) -> dict:
    """删除会话；仅当会话属于当前 owner_id 时允许。同时删除该会话下所有 Message。"""
    session = db.query(SessionModel).filter(SessionModel.id == session_id).first()
    if not session:
        raise HTTPException(status_code=404, detail="会话不存在")
    if session.owner_id != owner_id:
        raise HTTPException(status_code=status.HTTP_403_FORBIDDEN, detail="无权删除该会话")
    try:
        db.query(MessageModel).filter(MessageModel.session_id == session_id).delete()
        db.delete(session)
        db.commit()
        return {"status": "success"}
    except Exception as e:
        db.rollback()
        logger.exception("删除会话失败: %s", e)
        raise HTTPException(status_code=500, detail="删除失败")


@router.put("/api/sessions/{session_id}")
def rename_session(
    session_id: str,
    body: SessionRenameBody,
    owner_id: str = Depends(get_current_owner_id),
    db: Session = Depends(get_db_session),
) -> dict:
    """重命名会话；仅当会话属于当前 owner_id 时允许。"""
    session = (
        db.query(SessionModel)
        .filter(SessionModel.id == session_id, SessionModel.owner_id == owner_id)
        .first()
    )
    if not session:
        raise HTTPException(status_code=404, detail="会话不存在或无权修改")
    try:
        session.title = (body.title or "").strip() or session.title
        db.commit()
        db.refresh(session)
        return {"status": "success", "title": session.title}
    except Exception as e:
        db.rollback()
        logger.exception("重命名会话失败: %s", e)
        raise HTTPException(status_code=500, detail="更新失败")


@router.delete("/api/assets/{asset_id}")
def delete_asset(
    asset_id: int,
    owner_id: str = Depends(get_current_owner_id),
    db: Session = Depends(get_db_session),
) -> dict:
    """删除资产记录；仅当资产属于当前 owner_id 时允许。可选：尝试删除服务器上的物理文件。"""
    asset = db.query(AssetModel).filter(AssetModel.id == asset_id).first()
    if not asset:
        raise HTTPException(status_code=404, detail="资产不存在")
    if asset.owner_id != owner_id:
        raise HTTPException(status_code=status.HTTP_403_FORBIDDEN, detail="无权删除该资产")
    try:
        file_path = (asset.file_path or "").strip()
        if file_path and os.path.isfile(file_path):
            try:
                os.remove(file_path)
            except OSError as e:
                logger.warning("删除物理文件失败（已忽略）: %s", e)
        db.delete(asset)
        db.commit()
        return {"status": "success"}
    except Exception as e:
        db.rollback()
        logger.exception("删除资产失败: %s", e)
        raise HTTPException(status_code=500, detail="删除失败")


@router.put("/api/assets/{asset_id}")
def rename_asset(
    asset_id: int,
    body: AssetRenameBody,
    owner_id: str = Depends(get_current_owner_id),
    db: Session = Depends(get_db_session),
) -> dict:
    """重命名资产（仅更新 file_name 显示名）；仅当资产属于当前 owner_id 时允许。"""
    asset = (
        db.query(AssetModel)
        .filter(AssetModel.id == asset_id, AssetModel.owner_id == owner_id)
        .first()
    )
    if not asset:
        raise HTTPException(status_code=404, detail="资产不存在或无权修改")
    try:
        asset.file_name = (body.file_name or "").strip() or asset.file_name
        db.commit()
        db.refresh(asset)
        return {"status": "success", "file_name": asset.file_name}
    except Exception as e:
        db.rollback()
        logger.exception("重命名资产失败: %s", e)
        raise HTTPException(status_code=500, detail="更新失败")


@router.get("/api/workflow_templates")
def list_workflow_templates(
    owner_id: str = Depends(get_current_owner_id),
    db: Session = Depends(get_db_session),
) -> List[dict]:
    """按 owner_id 返回工作流收藏列表，按创建时间倒序。"""
    rows = (
        db.query(WorkflowTemplateModel)
        .filter(WorkflowTemplateModel.owner_id == owner_id)
        .order_by(WorkflowTemplateModel.created_at.desc())
        .all()
    )
    return [
        {
            "id": r.id,
            "owner_id": r.owner_id,
            "name": r.name,
            "config_json": r.config_json,
            "created_at": r.created_at.isoformat() if r.created_at else None,
        }
        for r in rows
    ]


@router.post("/api/workflow_templates")
def create_workflow_template(
    body: WorkflowTemplateCreate,
    owner_id: str = Depends(get_current_owner_id),
    db: Session = Depends(get_db_session),
) -> dict:
    """新建工作流收藏；返回 template_id。"""
    try:
        template = WorkflowTemplateModel(
            owner_id=owner_id,
            name=(body.name or "未命名工作流").strip() or "未命名工作流",
            config_json=body.config_json if body.config_json is not None else {},
        )
        db.add(template)
        db.commit()
        db.refresh(template)
        return {"template_id": template.id, "status": "success"}
    except Exception as e:
        db.rollback()
        logger.exception("创建工作流收藏失败: %s", e)
        raise HTTPException(status_code=500, detail="创建失败")


@router.delete("/api/workflow_templates/{template_id}")
def delete_workflow_template(
    template_id: int,
    owner_id: str = Depends(get_current_owner_id),
    db: Session = Depends(get_db_session),
) -> dict:
    """删除工作流收藏；仅当该记录属于当前 owner_id 时允许。"""
    template = db.query(WorkflowTemplateModel).filter(WorkflowTemplateModel.id == template_id).first()
    if not template:
        raise HTTPException(status_code=404, detail="工作流收藏不存在")
    if template.owner_id != owner_id:
        raise HTTPException(status_code=status.HTTP_403_FORBIDDEN, detail="无权删除该工作流收藏")
    try:
        db.delete(template)
        db.commit()
        return {"status": "success"}
    except Exception as e:
        db.rollback()
        logger.exception("删除工作流收藏失败: %s", e)
        raise HTTPException(status_code=500, detail="删除失败")


@router.put("/api/workflow_templates/{template_id}")
def rename_workflow_template(
    template_id: int,
    body: WorkflowTemplateRenameBody,
    owner_id: str = Depends(get_current_owner_id),
    db: Session = Depends(get_db_session),
) -> dict:
    """重命名工作流收藏；仅当该记录属于当前 owner_id 时允许。"""
    template = (
        db.query(WorkflowTemplateModel)
        .filter(
            WorkflowTemplateModel.id == template_id,
            WorkflowTemplateModel.owner_id == owner_id,
        )
        .first()
    )
    if not template:
        raise HTTPException(status_code=404, detail="工作流收藏不存在或无权修改")
    try:
        template.name = (body.name or "").strip() or template.name
        db.commit()
        db.refresh(template)
        return {"status": "success", "name": template.name}
    except Exception as e:
        db.rollback()
        logger.exception("重命名工作流收藏失败: %s", e)
        raise HTTPException(status_code=500, detail="更新失败")
