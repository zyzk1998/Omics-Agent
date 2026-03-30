"""
GIBH-AGENT-V2 测试服务器
提供简单的 Web 接口用于测试功能
"""
import os
import sys

import json
import logging
import traceback
import time
import asyncio
import re
import secrets
from pathlib import Path
from typing import List, Optional, Set, Dict, Any
from datetime import datetime
from collections import deque

from fastapi import FastAPI, UploadFile, File, HTTPException, Request, Depends, Query
from fastapi import status
from typing import Optional
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import StreamingResponse, JSONResponse, HTMLResponse, FileResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel

# 添加项目路径
sys.path.insert(0, str(Path(__file__).parent))

# 加载 .env（INITIAL_ADMINS 等），若未安装 python-dotenv 则跳过
try:
    from dotenv import load_dotenv
    load_dotenv()
except ImportError:
    pass

from gibh_agent import create_agent
from gibh_agent.core.file_inspector import FileInspector
from gibh_agent.core.orchestrator import AgentOrchestrator
from gibh_agent.core.file_handlers.structure_normalizer import normalize_session_directory
from gibh_agent.core.file_handlers.universal_normalizer import normalize_session_directory as universal_unpack
from gibh_agent.core.file_handlers.modality_sniffer import detect_dominant_modality, paths_for_response

# 配置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler('gibh_agent.log', encoding='utf-8')
    ]
)
logger = logging.getLogger(__name__)

# 创建 FastAPI 应用
app = FastAPI(title="GIBH-AGENT-V2 Test Server")


@app.exception_handler(Exception)
async def global_exception_handler(request: Request, exc: Exception):
    """全链路异常透传：将 Python Traceback 打包进 500 响应体，供前端控制台高亮打印。"""
    tb_str = traceback.format_exc()
    logger.error("🔥 [全局异常捕获] 500 Internal Server Error:\n%s", tb_str)
    return JSONResponse(
        status_code=500,
        content={
            "detail": "Internal Server Error",
            "error_message": str(exc),
            "traceback": tb_str,
        },
    )


# 🔥 Step 2: Tool-RAG 架构 - Vector Database Integration
# 初始化工具检索器（在启动时同步工具）
tool_retriever = None
workflow_planner = None
try:
    from gibh_agent.core.tool_retriever import ToolRetriever
    # 🔥 Step 4: 模块化工具系统 - 自动发现和加载所有工具
    from gibh_agent.tools import load_all_tools
    
    # 初始化工具检索器
    chroma_dir = os.getenv("CHROMA_PERSIST_DIR", "./data/chroma_tools")
    embedding_model = os.getenv("OLLAMA_EMBEDDING_MODEL", "nomic-embed-text")
    ollama_url = os.getenv("OLLAMA_BASE_URL", "http://localhost:11434")
    
    logger.info(f"🔧 初始化工具检索器...")
    logger.info(f"   ChromaDB 目录: {chroma_dir}")
    logger.info(f"   Embedding 模型: {embedding_model}")
    logger.info(f"   Ollama URL: {ollama_url}")
    
    tool_retriever = ToolRetriever(
        persist_directory=chroma_dir,
        embedding_model=embedding_model,
        ollama_base_url=ollama_url
    )
    
    logger.info("✅ 工具检索器初始化成功")
except ImportError as e:
    logger.warning(f"⚠️ 工具检索器依赖未安装: {e}")
    logger.warning("   跳过工具检索器初始化（需要: pip install langchain-chroma langchain-ollama）")
except Exception as e:
    logger.error(f"❌ 工具检索器初始化失败: {e}", exc_info=True)
    logger.warning("   继续启动，但工具检索功能将不可用")

# 🔥 Step 3: 初始化工作流规划器（需要 agent 初始化后才能获取 LLM client）
# 这将在 agent 初始化后设置

# 配置 CORS（安全配置）
# 生产环境应该限制为特定域名
ALLOWED_ORIGINS = os.getenv("ALLOWED_ORIGINS", "*").split(",")
if ALLOWED_ORIGINS == ["*"]:
    logger.warning("⚠️  CORS 配置允许所有来源，生产环境应限制为特定域名")

app.add_middleware(
    CORSMiddleware,
    allow_origins=ALLOWED_ORIGINS,
    allow_credentials=True,
    allow_methods=["GET", "POST", "PUT", "DELETE", "OPTIONS"],
    allow_headers=["Content-Type", "Authorization", "X-Guest-UUID"],
)

# Phase 2: 鉴权与游客资产继承（仅注册路由，不修改现有业务）
try:
    from gibh_agent.api.routers.auth import router as auth_router
    app.include_router(auth_router, prefix="/api/auth")
    logger.info("✅ Auth 路由已注册: /api/auth/register, /api/auth/login, /api/auth/merge_guest_data")
except Exception as e:
    logger.warning("⚠️ Auth 路由注册失败（鉴权功能不可用）: %s", e)

# Phase 4: 侧栏数据读取（会话/消息/资产）
try:
    from gibh_agent.api.routers.user_data import router as user_data_router
    app.include_router(user_data_router)
    logger.info("✅ User Data 路由已注册: GET /api/sessions, /api/sessions/{id}/messages, /api/assets")
except Exception as e:
    logger.warning("⚠️ User Data 路由注册失败: %s", e)

# UGC 技能上传与管理员审核
try:
    from gibh_agent.api.routers.skills import router as skills_router
    app.include_router(skills_router)
    logger.info("✅ Skills 路由已注册: POST /api/skills, GET /api/admin/skills, PUT status, POST review")
except Exception as e:
    logger.warning("⚠️ Skills 路由注册失败: %s", e)

# Phase 4: 身份解析与 DB 依赖（upload/chat 持久化）
from gibh_agent.core.deps import get_current_owner_id, get_current_admin_user
from gibh_agent.db.connection import get_db_session, engine, is_available, Base
from sqlalchemy.orm import Session
from sqlalchemy import text, case, or_
from sqlalchemy.exc import DataError, OperationalError

# 注册 ORM 模型以便 create_all 建表
try:
    from gibh_agent.db import models  # noqa: F401
    from gibh_agent.db.models import Skill as SkillModel, UserSavedSkill
except Exception:
    SkillModel = None
    UserSavedSkill = None


def _run_create_all():
    """执行建表（幂等）。确保 models 已挂到 Base 后再调用 create_all。"""
    import gibh_agent.db.models as _  # noqa: F401 强制注册所有表
    Base.metadata.create_all(bind=engine)


# 7 大核心组学技能：CORE_OMICS_SKILLS 字典列表定义在 gibh_agent/db/seed_skills.py，修改 prompt_template 请编辑该文件。此处仅调用 run_seed_core_skills 注入。
def _insert_core_skills_only(db: Session) -> None:
    """向当前 Session 仅插入 7 大核心组学技能（不 commit，由调用方提交）。"""
    from gibh_agent.db.seed_skills import run_seed_core_skills
    run_seed_core_skills(db)


def _bootstrap_skills_mock():
    """幂等注入：由 seed_skills.run_upsert_system_skills 内部先暴力清洗 author_id=system，再逐条 Upsert（含 IntegrityError 防并发）。
    【红线】启动时禁止对 users/sessions/messages/assets/workflow_templates 做任何 delete/drop/truncate。"""
    if not is_available() or engine is None or SkillModel is None:
        logger.warning("[DB] 技能注入跳过: 数据库或 SkillModel 不可用")
        return
    from gibh_agent.db.connection import SessionLocal
    from gibh_agent.db.seed_skills import run_upsert_system_skills
    db = SessionLocal()
    try:
        run_upsert_system_skills(db)
        logger.info("✅ [DB] 7 大核心组学 + 生物医药 + 化学技能已幂等注入并就绪！")
    except Exception as e:
        db.rollback()
        logger.error("❌ 技能注入失败: %s", e, exc_info=True)
    finally:
        db.close()


@app.on_event("startup")
def _ensure_db_tables():
    """等待 MySQL 就绪（最多 30 秒）后建表，避免 Docker 下 API 先于 MySQL 启动导致假死。"""
    if not is_available() or engine is None or Base is None:
        return
    max_wait = 30
    for attempt in range(max_wait):
        try:
            with engine.connect() as conn:
                conn.execute(text("SELECT 1"))
            break
        except Exception as e:
            if attempt == 0:
                logger.info("[DB] 等待 MySQL 就绪...")
            if attempt >= max_wait - 1:
                logger.warning(
                    "[DB] MySQL 未在 %d 秒内就绪，跳过建表。可稍后调用 GET /api/db/init 建表: %s",
                    max_wait,
                    e,
                )
                return
            time.sleep(1)
    try:
        _run_create_all()
        logger.info("✅ [DB] 表结构已就绪（如需建表已自动创建）")
        _bootstrap_skills_mock()
    except Exception as e:
        logger.error(
            "❌ [DB] 建表失败，后续查询可能 500。请调用 GET /api/db/init 重试或检查 MySQL 权限/版本。原因: %s\n%s",
            e,
            traceback.format_exc(),
        )


@app.get("/api/db/init")
def init_db_tables():
    """
    手动触发建表（startup 未执行或失败时使用）。幂等，仅创建缺失表；
    会同时触发系统技能强制刷新（仅清理 author_id=system 后注入 7 大核心组学技能）。
    部署后若出现 Table 'xxx' doesn't exist 或技能广场为空，访问: GET http://<host>:8028/api/db/init
    """
    if not is_available() or engine is None or Base is None:
        return JSONResponse(status_code=503, content={"detail": "数据库不可用，无法建表"})
    try:
        _run_create_all()
        _bootstrap_skills_mock()
        return JSONResponse(status_code=200, content={"detail": "表已就绪，技能表已检查/注入", "ok": True})
    except Exception as e:
        logger.error("建表失败: %s\n%s", e, traceback.format_exc())
        detail = _detail_for_db_error(e, str(e))
        return JSONResponse(
            status_code=500,
            content={
                "detail": detail,
                "error_type": type(e).__name__,
                "error_message": str(e),
                "traceback": traceback.format_exc(),
            },
        )


@app.get("/api/skills")
def list_skills_public(
    request: Request,
    main_cat: Optional[str] = Query(None, description="大类筛选"),
    sub_cat: Optional[str] = Query(None, description="小类筛选"),
    saved_only: bool = Query(False, description="仅返回当前用户收藏的技能"),
    page: int = Query(1, ge=1, description="页码"),
    size: int = Query(12, ge=1, le=50, description="每页条数"),
    db: Session = Depends(get_db_session),
):
    """技能广场公开分页查询：仅 status=approved；saved_only=True 时需鉴权并只返回收藏。无数据时自动补种 7 大核心组学。"""
    if SkillModel is None:
        return JSONResponse(status_code=503, content={"detail": "技能服务不可用", "items": [], "total": 0})
    owner_id = None
    try:
        owner_id = get_current_owner_id(request)
    except HTTPException:
        pass
    if saved_only and not owner_id:
        raise HTTPException(status_code=status.HTTP_401_UNAUTHORIZED, detail="查看「我的」收藏需要登录或提供身份", headers={"WWW-Authenticate": "Bearer"})
    mc = (main_cat or "").strip()
    sc = (sub_cat or "").strip()
    q = db.query(SkillModel).filter(SkillModel.status == "approved")
    if saved_only and owner_id and UserSavedSkill is not None:
        q = q.join(UserSavedSkill, (UserSavedSkill.skill_id == SkillModel.id) & (UserSavedSkill.owner_id == owner_id))
    if mc:
        q = q.filter(SkillModel.main_category == mc)
    if sc:
        q = q.filter(SkillModel.sub_category == sc)
    total = q.count()
    if total == 0 and not saved_only:
        try:
            from gibh_agent.db.seed_skills import run_upsert_system_skills
            run_upsert_system_skills(db)
            logger.info("[Skills] 已按需幂等补种核心组学 + 生物医药 + 化学技能")
        except Exception as e:
            db.rollback()
            logger.warning("[Skills] 按需补种失败: %s", e)
        q = db.query(SkillModel).filter(SkillModel.status == "approved")
        if saved_only and owner_id and UserSavedSkill is not None:
            q = q.join(UserSavedSkill, (UserSavedSkill.skill_id == SkillModel.id) & (UserSavedSkill.owner_id == owner_id))
        if mc:
            q = q.filter(SkillModel.main_category == mc)
        if sc:
            q = q.filter(SkillModel.sub_category == sc)
        total = q.count()
    offset = (page - 1) * size
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
    saved_ids = set()
    if owner_id and UserSavedSkill is not None:
        saved_rows = db.query(UserSavedSkill.skill_id).filter(UserSavedSkill.owner_id == owner_id).all()
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


@app.post("/api/skills/{skill_id}/bookmark")
def bookmark_skill(
    skill_id: int,
    owner_id: str = Depends(get_current_owner_id),
    db: Session = Depends(get_db_session),
):
    """收藏技能：防重复。与 skills 路由重复注册以保证 router 未加载时仍可用。"""
    if SkillModel is None or UserSavedSkill is None:
        raise HTTPException(status_code=503, detail="技能服务不可用")
    skill = db.query(SkillModel).filter(SkillModel.id == skill_id, SkillModel.status == "approved").first()
    if not skill:
        raise HTTPException(status_code=404, detail="技能不存在")
    existing = (
        db.query(UserSavedSkill)
        .filter(UserSavedSkill.owner_id == owner_id, UserSavedSkill.skill_id == skill_id)
        .first()
    )
    if not existing:
        db.add(UserSavedSkill(owner_id=owner_id, skill_id=skill_id))
        db.commit()
    return {"status": "success", "message": "已添加到我的工具", "saved": True}


@app.delete("/api/skills/{skill_id}/bookmark")
def remove_skill_bookmark(
    skill_id: int,
    owner_id: str = Depends(get_current_owner_id),
    db: Session = Depends(get_db_session),
):
    """取消收藏。与 skills 路由重复注册以保证 router 未加载时仍可用。"""
    if UserSavedSkill is None:
        raise HTTPException(status_code=503, detail="技能服务不可用")
    bookmark = (
        db.query(UserSavedSkill)
        .filter(UserSavedSkill.owner_id == owner_id, UserSavedSkill.skill_id == skill_id)
        .first()
    )
    if bookmark:
        db.delete(bookmark)
        db.commit()
    return {"status": "success", "message": "已取消收藏", "saved": False}


@app.post("/api/admin/bootstrap-skills")
def admin_bootstrap_skills(
    db: Session = Depends(get_db_session),
    _admin=Depends(get_current_admin_user),
):
    """
    管理员专用：仅清理 author_id=system 的旧系统技能后，重新注入 7 大核心组学技能，保留用户上传技能。
    技能广场需要恢复默认技能时，以管理员身份调用此接口即可。
    """
    if not is_available() or SkillModel is None:
        raise HTTPException(status_code=503, detail="数据库不可用")
    try:
        from gibh_agent.db.seed_skills import run_upsert_system_skills
        deleted = db.query(SkillModel).filter(SkillModel.author_id == "system").delete(synchronize_session=False)
        db.commit()
        run_upsert_system_skills(db)
        logger.info("✅ [DB] 管理员已触发：清理 %d 条旧系统技能并幂等注入核心组学 + 生物医药 + 化学", deleted)
        return JSONResponse(status_code=200, content={"detail": "已清理旧系统技能并重新注入 7 大核心组学技能", "ok": True})
    except Exception as e:
        db.rollback()
        logger.exception("管理员注入技能失败: %s", e)
        raise HTTPException(status_code=500, detail="注入失败: " + str(e))


def _detail_for_db_error(exc: Exception, err_msg: str) -> str:
    """对常见 MySQL 错误返回可读说明，便于运维排查。"""
    try:
        from sqlalchemy.exc import OperationalError
        if isinstance(exc, OperationalError) and ("1030" in err_msg or "error 168" in err_msg or "error from engine" in err_msg.lower()):
            return (
                "MySQL 存储引擎错误(1030/168)：多为数据目录权限或磁盘问题。"
                "请检查 MySQL 数据目录(如 /var/lib/mysql)权限、磁盘空间，并查看 MySQL 错误日志。"
            )
    except ImportError:
        pass
    return "服务器内部致命错误"


# 全局异常显微镜：捕获所有未处理异常；若为「表不存在」则尝试自动建表并返回提示
@app.exception_handler(Exception)
async def global_exception_handler(request: Request, exc: Exception):
    traceback.print_exc()
    err_msg = str(exc)
    # 表不存在时尝试自动建表，避免 500 循环
    if "doesn't exist" in err_msg and ("Table " in err_msg or "table " in err_msg.lower()):
        try:
            from sqlalchemy.exc import ProgrammingError
            if isinstance(exc, ProgrammingError) and is_available() and engine is not None:
                _run_create_all()
                logger.info("[DB] 已根据表不存在错误自动建表")
                return JSONResponse(
                    status_code=503,
                    content={"detail": "表已自动创建，请刷新页面后重试", "ok": True},
                )
        except Exception as e:
            logger.warning("[DB] 自动建表失败: %s", e)
            detail = _detail_for_db_error(e, str(e)) if "1030" in str(e) or "168" in str(e) else "自动建表失败，请检查 MySQL 权限与磁盘"
            return JSONResponse(
                status_code=503,
                content={"detail": detail, "error_type": type(e).__name__, "error_message": str(e)},
            )
    detail = _detail_for_db_error(exc, err_msg)
    return JSONResponse(
        status_code=500,
        content={
            "detail": detail,
            "error_type": type(exc).__name__,
            "error_message": err_msg,
            "traceback": traceback.format_exc(),
        },
    )

# 创建上传目录（使用绝对路径，兼容容器环境）
# 🔧 修复：优先使用环境变量，否则使用容器内的默认路径
# 注意：在容器环境中，默认路径应该是 /app/uploads，而不是相对路径
UPLOAD_DIR = Path(os.getenv("UPLOAD_DIR", "/app/uploads"))
RESULTS_DIR = Path(os.getenv("RESULTS_DIR", "/app/results"))


def _ensure_absolute_upload_path(path: str, base: Optional[Path] = None) -> str:
    """将相对路径转为基于 base 的绝对路径；已是绝对路径则直接 resolve 返回，避免重复拼接。"""
    if not path or not str(path).strip():
        return path
    base = base or UPLOAD_DIR
    p = Path(path)
    if p.is_absolute():
        return str(p.resolve())
    return str((base / p).resolve())

# 确保目录存在且可写
try:
    UPLOAD_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    
    # 检查目录权限
    if not os.access(UPLOAD_DIR, os.W_OK):
        logger.warning(f"⚠️ 上传目录不可写: {UPLOAD_DIR}")
    if not os.access(RESULTS_DIR, os.W_OK):
        logger.warning(f"⚠️ 结果目录不可写: {RESULTS_DIR}")
    
    logger.info(f"📁 上传目录: {UPLOAD_DIR.absolute()} (可写: {os.access(UPLOAD_DIR, os.W_OK)})")
    logger.info(f"📁 结果目录: {RESULTS_DIR.absolute()} (可写: {os.access(RESULTS_DIR, os.W_OK)})")
except Exception as e:
    logger.error(f"❌ 创建目录失败: {e}", exc_info=True)
    raise

# ---------------------------------------------------------------------------
# 存储守护：APScheduler 定时清理 + 管理员手动触发 API（依赖 UPLOAD_DIR / RESULTS_DIR）
# ---------------------------------------------------------------------------
_storage_scheduler = None


def _scheduled_storage_cleanup_job() -> None:
    """后台任务：清理过期重型文件，保留 .csv/.md/.png/.json 等报告。"""
    try:
        from gibh_agent.utils.storage_manager import StorageManager

        sm = StorageManager(results_dir=str(RESULTS_DIR), upload_dir=str(UPLOAD_DIR))
        r = sm.auto_cleanup(
            max_days=int(os.getenv("STORAGE_CLEANUP_MAX_DAYS", "7")),
            max_size_gb=float(os.getenv("STORAGE_CLEANUP_MAX_SIZE_GB", "100")),
            deep=False,
            dry_run=False,
        )
        logger.info(
            "[StorageDaemon] 定时清理完成 deleted=%s freed_gb=%.4f errors=%s",
            r.get("deleted_count"),
            r.get("freed_gb", 0),
            len(r.get("errors") or []),
        )
    except Exception:
        logger.exception("[StorageDaemon] 定时清理失败")


@app.on_event("startup")
def _mount_storage_cleanup_scheduler() -> None:
    global _storage_scheduler
    if os.getenv("STORAGE_CLEANUP_DISABLED", "").lower() in ("1", "true", "yes"):
        logger.info("📦 存储自动清理已禁用（STORAGE_CLEANUP_DISABLED）")
        return
    try:
        from apscheduler.schedulers.background import BackgroundScheduler
    except ImportError:
        logger.warning("⚠️ 未安装 APScheduler，跳过存储定时清理。请: pip install APScheduler")
        return
    _storage_scheduler = BackgroundScheduler()
    hours = int(os.getenv("STORAGE_CLEANUP_INTERVAL_HOURS", "12"))
    _storage_scheduler.add_job(
        _scheduled_storage_cleanup_job,
        "interval",
        hours=hours,
        id="storage_auto_cleanup",
        replace_existing=True,
        max_instances=1,
        coalesce=True,
    )
    _storage_scheduler.start()
    logger.info("📦 存储守护已挂载：每 %s 小时执行 StorageManager.auto_cleanup", hours)


@app.on_event("shutdown")
def _shutdown_storage_cleanup_scheduler() -> None:
    global _storage_scheduler
    if _storage_scheduler is not None:
        try:
            _storage_scheduler.shutdown(wait=False)
        except Exception:
            pass
        _storage_scheduler = None


@app.delete("/api/admin/storage/cleanup")
def admin_storage_cleanup(
    max_days: int = Query(7, ge=1, le=3650),
    max_size_gb: float = Query(100.0, ge=1.0, le=100000.0),
    deep: bool = Query(False),
    dry_run: bool = Query(False),
    _admin=Depends(get_current_admin_user),
):
    """
    管理员：强制触发存储清理。返回删除数量与释放空间（MB/GB）。
    deep=true 时有效保留期缩短；dry_run=true 仅统计不删除。
    """
    try:
        from gibh_agent.utils.storage_manager import StorageManager

        sm = StorageManager(results_dir=str(RESULTS_DIR), upload_dir=str(UPLOAD_DIR))
        r = sm.auto_cleanup(
            max_days=max_days,
            max_size_gb=max_size_gb,
            deep=deep,
            dry_run=dry_run,
        )
        return JSONResponse(status_code=200, content={"ok": True, **r})
    except Exception as e:
        logger.exception("管理员存储清理失败: %s", e)
        raise HTTPException(status_code=500, detail=str(e))


# 安全配置
MAX_FILE_SIZE = int(os.getenv("MAX_FILE_SIZE", 200 * 1024 * 1024))  # 默认 200MB
# .h5 = 10x Visium / scRNA Feature Barcode Matrix (industry standard); .h5ad = AnnData
ALLOWED_EXTENSIONS = {'.h5', '.h5ad', '.mtx', '.tsv', '.csv', '.txt', '.gz', '.tar', '.zip'}
ALLOWED_MIME_TYPES = {
    'text/csv', 'text/tab-separated-values', 'text/plain',
    'application/gzip', 'application/x-gzip',
    'application/zip', 'application/x-tar'
}

def sanitize_filename(filename: str) -> str:
    """
    清理文件名，防止路径遍历攻击
    
    Args:
        filename: 原始文件名
    
    Returns:
        清理后的安全文件名
    """
    if not filename:
        # 如果文件名为空，生成随机名称
        return f"file_{secrets.token_hex(8)}"
    
    # 移除路径分隔符和危险字符
    filename = os.path.basename(filename)  # 移除路径部分
    filename = re.sub(r'[<>:"|?*\x00-\x1f]', '', filename)  # 移除危险字符
    filename = filename.strip('. ')  # 移除开头和结尾的点/空格
    
    # 如果清理后为空，生成随机名称
    if not filename:
        filename = f"file_{secrets.token_hex(8)}"
    
    # 限制文件名长度
    if len(filename) > 255:
        name, ext = os.path.splitext(filename)
        filename = name[:255-len(ext)] + ext
    
    return filename

def validate_file_path(file_path: Path, base_dir: Path) -> Path:
    """
    验证文件路径是否在允许的目录内（防止路径遍历）
    
    Args:
        file_path: 要验证的路径
        base_dir: 基础目录
    
    Returns:
        规范化的安全路径
    
    Raises:
        HTTPException: 如果路径不安全
    """
    try:
        # 解析并规范化路径
        resolved_path = file_path.resolve()
        resolved_base = base_dir.resolve()
        
        # 检查路径是否在基础目录内
        if not str(resolved_path).startswith(str(resolved_base)):
            raise HTTPException(
                status_code=403,
                detail="文件路径不安全：不允许访问基础目录外的文件"
            )
        
        return resolved_path
    except (ValueError, OSError) as e:
        raise HTTPException(
            status_code=400,
            detail=f"无效的文件路径: {str(e)}"
        )

# 初始化文件检测器
file_inspector = FileInspector(str(UPLOAD_DIR))

# 添加静态文件服务（用于访问结果图片、前端 Logo 等）；使用绝对路径与 RESULTS_DIR/UPLOAD_DIR 一致，避免打印/预览时图片 404
from fastapi.staticfiles import StaticFiles
_results_dir = str(RESULTS_DIR.resolve()) if RESULTS_DIR.exists() else "results"
_uploads_dir = str(UPLOAD_DIR.resolve()) if UPLOAD_DIR.exists() else "uploads"
app.mount("/results", StaticFiles(directory=_results_dir), name="results")
app.mount("/uploads", StaticFiles(directory=_uploads_dir), name="uploads")
# 前端静态资源（Logo 等）：与 index.html 同级的 static 目录
_static_dir = Path(__file__).parent / "services" / "nginx" / "html" / "static"
if _static_dir.is_dir():
    app.mount("/static", StaticFiles(directory=str(_static_dir)), name="static")

# 帮助文档与 API 文档（与 index.html 同级，供用户菜单点击访问）
_html_dir = Path(__file__).parent / "services" / "nginx" / "html"
# 直连 API 端口（如 :8028）时首页由本进程返回 HTML，但浏览器仍请求 /css、/js；需挂载否则 404
_css_dir = _html_dir / "css"
_js_dir = _html_dir / "js"
if _css_dir.is_dir():
    app.mount("/css", StaticFiles(directory=str(_css_dir)), name="css")
if _js_dir.is_dir():
    app.mount("/js", StaticFiles(directory=str(_js_dir)), name="js")


@app.get("/whotowork.html")
async def serve_whotowork():
    """帮助文档：如何工作 / 使用说明"""
    path = _html_dir / "whotowork.html"
    if not path.exists():
        raise HTTPException(status_code=404, detail="Not Found")
    return FileResponse(path, media_type="text/html; charset=utf-8")


@app.get("/API.md")
async def serve_api_md():
    """API 文档（原始 Markdown）"""
    path = _html_dir / "API.md"
    if not path.exists():
        path = Path(__file__).parent / "API.md"
    if not path.exists():
        raise HTTPException(status_code=404, detail="Not Found")
    return FileResponse(path, media_type="text/markdown; charset=utf-8")


@app.get("/api-doc.html")
async def serve_api_doc_html():
    """API 文档渲染页：优先返回静态文件，不存在则返回内联 HTML 查看器（fetch /API.md + marked 渲染）"""
    path = _html_dir / "api-doc.html"
    if path.exists():
        return FileResponse(path, media_type="text/html; charset=utf-8")
    # 回退：内联最小可用查看器，确保路由不 404
    inline_html = """<!DOCTYPE html><html lang="zh-CN"><head><meta charset="UTF-8"><title>API 文档</title>
<script src="https://cdn.jsdelivr.net/npm/marked/marked.min.js"></script></head><body style="font-family:sans-serif;max-width:960px;margin:0 auto;padding:24px;">
<div id="c">正在加载 /API.md…</div>
<script>
fetch('/API.md').then(function(r){if(!r.ok)throw new Error(r.status);return r.text();})
.then(function(md){var c=document.getElementById('c');if(typeof marked!=='undefined'){marked.setOptions({gfm:true});c.innerHTML=marked.parse(md);}else c.textContent=md;})
.catch(function(e){document.getElementById('c').textContent='加载失败: '+e.message;});
</script></body></html>"""
    return HTMLResponse(inline_html)

# 初始化智能体
agent = None
try:
    # 尝试从当前目录加载配置
    config_path = Path(__file__).parent / "gibh_agent" / "config" / "settings.yaml"
    logger.info(f"🔍 查找配置文件: {config_path}")
    logger.info(f"📂 配置文件存在: {config_path.exists()}")
    
    if not config_path.exists():
        # 如果不存在，尝试其他路径
        alt_path = Path(__file__).parent / "config" / "settings.yaml"
        logger.info(f"🔍 尝试备用路径: {alt_path}")
        if alt_path.exists():
            config_path = alt_path
        else:
            config_path = "gibh_agent/config/settings.yaml"
            logger.info(f"🔍 使用默认路径: {config_path}")
    
    logger.info(f"📄 使用配置文件: {config_path}")
    
    # 设置 scanpy 工具的默认输出目录（使用相对路径）
    import os
    scanpy_output_dir = os.path.join(os.getcwd(), "results")
    logger.info(f"📁 Scanpy 输出目录: {scanpy_output_dir}")
    
    # 创建智能体
    agent = create_agent(str(config_path))
    
    # 更新 scanpy 工具的输出目录
    if agent and hasattr(agent, 'agents') and 'rna_agent' in agent.agents:
        rna_agent = agent.agents['rna_agent']
        if hasattr(rna_agent, 'scanpy_tool'):
            rna_agent.scanpy_tool.output_dir = scanpy_output_dir
            os.makedirs(scanpy_output_dir, exist_ok=True)
            logger.info(f"✅ 已设置 Scanpy 输出目录: {scanpy_output_dir}")
    
    logger.info("✅ GIBH-AGENT 初始化成功")
    
    # 🔥 Step 3: 初始化工作流规划器（需要 agent 和 tool_retriever）
    if agent and tool_retriever:
        try:
            from gibh_agent.core.planner import WorkflowPlanner
            # 获取 LLM client（从 agent 的某个智能体中获取）
            llm_client = None
            if hasattr(agent, 'agents') and agent.agents:
                # 尝试从第一个智能体获取 LLM client
                first_agent = list(agent.agents.values())[0]
                if hasattr(first_agent, 'llm_client'):
                    llm_client = first_agent.llm_client
            
            if llm_client:
                workflow_planner = WorkflowPlanner(
                    tool_retriever=tool_retriever,
                    llm_client=llm_client
                )
                logger.info("✅ 工作流规划器初始化成功")
            else:
                logger.warning("⚠️ 无法获取 LLM client，跳过工作流规划器初始化")
        except Exception as e:
            logger.error(f"❌ 工作流规划器初始化失败: {e}", exc_info=True)
            logger.warning("   继续启动，但动态规划功能将不可用")
    
except Exception as e:
    import traceback
    error_msg = f"❌ GIBH-AGENT 初始化失败: {e}"
    logger.error(error_msg, exc_info=True)
    logger.error(f"详细错误:\n{traceback.format_exc()}")
    agent = None

# 🔥 Step 2: 启动时同步工具到 ChromaDB
@app.on_event("startup")
async def sync_tools_on_startup():
    """
    启动时同步工具到 Vector Database
    
    确保 ChromaDB 中的工具定义与代码中的 @register 装饰器保持一致。
    """
    # 🔥 Step 4: 首先加载所有工具模块（自动发现）
    try:
        logger.info("🔍 启动时自动发现和加载工具模块...")
        load_result = load_all_tools()
        logger.info(f"✅ 工具模块加载完成: {load_result['loaded']} 个成功, {load_result['failed']} 个失败")
    except Exception as e:
        logger.error(f"❌ 工具模块加载失败: {e}", exc_info=True)
        logger.warning("   继续启动，但工具可能未完全加载")

    # 🔌 MCP 插件（与 tools/ 物理隔离，注册到同一 ToolRegistry）
    try:
        from gibh_agent.mcp import load_mcp_plugins
        load_mcp_plugins()
        logger.info("✅ MCP 插件模块已加载")
    except Exception as e:
        logger.warning("⚠️ MCP 插件加载失败: %s", e, exc_info=True)
    
    # 然后同步到 ChromaDB
    if tool_retriever is None:
        logger.warning("⚠️ 工具检索器未初始化，跳过工具同步")
        return
    
    try:
        logger.info("🔄 启动时同步工具到 ChromaDB...")
        synced_count = tool_retriever.sync_tools(clear_existing=True)
        logger.info(f"✅ 工具同步完成: {synced_count} 个工具已同步到 ChromaDB")
    except Exception as e:
        logger.error(f"❌ 工具同步失败: {e}", exc_info=True)
        logger.warning("   继续启动，但工具检索功能可能不可用")


# 请求模型
class ChatRequest(BaseModel):
    message: str = ""
    history: List[dict] = []
    uploaded_files: List[dict] = []
    workflow_data: Optional[dict] = None
    test_dataset_id: Optional[str] = None
    stream: Optional[bool] = False  # 🔥 SSE 流式传输开关
    session_id: Optional[str] = None  # 🔥 BUG FIX: 添加 session_id 字段
    user_id: Optional[str] = "guest"  # 🔥 BUG FIX: 添加 user_id 字段，默认为 guest
    model_name: Optional[str] = "qwen3.5-plus"  # 🔥 前端模型切换：默认阿里百炼 Qwen3.5-Plus，qwen 走 DashScope，其余走 SiliconFlow
    target_domain: Optional[str] = None  # 🔥 快车道硬路由：技能广场点击时传入 'rna'|'metabolomics'|'radiomics'|'spatial'，跳过意图识别
    enabled_mcps: Optional[List[str]] = None  # 🔌 前端 MCP 开关：如 ["web_search"]，与系统设置联动


# 日志缓冲区（保留用于未来扩展）
log_buffer = deque(maxlen=1000)
log_listeners: Set[asyncio.Queue] = set()


def log_handler(record):
    """日志处理器，将日志发送到所有监听者"""
    log_entry = {
        "timestamp": datetime.now().isoformat(),
        "level": record.levelname,
        "message": record.getMessage(),
        "module": record.name
    }
    log_buffer.append(log_entry)
    
    # 通知所有监听者
    for listener in list(log_listeners):
        try:
            listener.put_nowait(log_entry)
        except:
            # 如果队列已满或已关闭，移除监听者
            log_listeners.discard(listener)


# 添加自定义日志处理器
class StreamLogHandler(logging.Handler):
    def emit(self, record):
        try:
            # 确保记录被格式化
            self.format(record)
            log_handler(record)
        except Exception as e:
            # 避免日志处理器本身出错，但记录错误
            print(f"日志处理器错误: {e}")


stream_handler = StreamLogHandler()
stream_handler.setLevel(logging.DEBUG)  # 降低级别以捕获更多日志
stream_handler.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))

# 添加到根日志记录器，捕获所有模块的日志
root_logger = logging.getLogger()
root_logger.setLevel(logging.DEBUG)  # 降低级别
# 移除现有的处理器，避免重复
for handler in root_logger.handlers[:]:
    root_logger.removeHandler(handler)
root_logger.addHandler(stream_handler)

# 也添加到当前logger
if stream_handler not in logger.handlers:
    logger.addHandler(stream_handler)

# 测试日志
logger.info("📋 日志系统初始化完成")
logger.info("🔍 测试日志输出 - 这应该出现在前端")


@app.get("/health")
async def health_check():
    """健康检查端点"""
    return {"status": "ok", "service": "GIBH-AGENT-V2"}


@app.get("/api/health")
async def api_health_check():
    """API 健康检查端点"""
    return {
        "status": "ok",
        "service": "GIBH-AGENT-V2",
        "agent_initialized": agent is not None,
        "tool_retriever_initialized": tool_retriever is not None
    }


@app.get("/", response_class=HTMLResponse)
async def index():
    """返回前端页面"""
    # 🔥 优先读取外部 HTML 文件（如果存在）
    html_file_path = Path(__file__).parent / "services" / "nginx" / "html" / "index.html"
    if html_file_path.exists():
        try:
            with open(html_file_path, "r", encoding="utf-8") as f:
                html_content = f.read()
            logger.info(f"✅ 已加载外部前端文件: {html_file_path}")
            return HTMLResponse(content=html_content)
        except Exception as e:
            logger.warning(f"⚠️ 读取外部 HTML 文件失败，使用内嵌版本: {e}")
    
    # 如果外部文件不存在，使用内嵌的 HTML
    html_content = """
<!DOCTYPE html>
<html lang="zh-CN">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>GIBH-AGENT-V2 测试界面</title>
    <style>
        * {
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }
        body {
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            background: #f5f5f5;
            padding: 20px;
        }
        .container {
            max-width: 1200px;
            margin: 0 auto;
            display: flex;
            flex-direction: column;
            gap: 20px;
            height: calc(100vh - 40px);
        }
        .panel {
            background: white;
            border-radius: 8px;
            padding: 20px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
            display: flex;
            flex-direction: column;
        }
        .panel h2 {
            margin-bottom: 15px;
            color: #333;
            border-bottom: 2px solid #4CAF50;
            padding-bottom: 10px;
        }
        .chat-panel {
            flex: 1;
        }
        .chat-area {
            flex: 1;
            border: 1px solid #ddd;
            border-radius: 4px;
            padding: 15px;
            overflow-y: auto;
            overflow-x: hidden;
            margin-bottom: 15px;
            background: #fafafa;
            min-height: 300px;
            word-wrap: break-word;
            overflow-wrap: break-word;
        }
        .message {
            margin-bottom: 10px;
            padding: 8px;
            border-radius: 4px;
            word-wrap: break-word;
            overflow-wrap: break-word;
            max-width: 100%;
            box-sizing: border-box;
            overflow-x: hidden;
        }
        .message.user {
            background: #e3f2fd;
            text-align: right;
        }
        .message.assistant {
            background: #f1f8e9;
        }
        .message.error {
            background: #ffebee;
            color: #c62828;
        }
        .input-area {
            display: flex;
            gap: 10px;
        }
        input[type="text"], input[type="file"] {
            flex: 1;
            padding: 10px;
            border: 1px solid #ddd;
            border-radius: 4px;
            font-size: 14px;
        }
        button {
            padding: 10px 20px;
            background: #4CAF50;
            color: white;
            border: none;
            border-radius: 4px;
            cursor: pointer;
            font-size: 14px;
        }
        button:hover {
            background: #45a049;
        }
        button:disabled {
            background: #ccc;
            cursor: not-allowed;
        }
        .file-info {
            margin-top: 10px;
            padding: 10px;
            background: #fff3cd;
            border-radius: 4px;
            font-size: 12px;
        }
        .analysis-result {
            background: #f1f8e9 !important;
            word-wrap: break-word;
            overflow-wrap: break-word;
            max-width: 100%;
            box-sizing: border-box;
            overflow-x: hidden;
        }
        .analysis-summary {
            padding: 15px;
            max-width: 100%;
            box-sizing: border-box;
            overflow-x: hidden;
        }
        .analysis-summary h3 {
            margin-top: 0;
            color: #4CAF50;
        }
        .analysis-summary h4 {
            margin-top: 15px;
            margin-bottom: 10px;
            color: #333;
            border-bottom: 1px solid #ddd;
            padding-bottom: 5px;
        }
        .analysis-summary ul {
            margin: 10px 0;
            padding-left: 20px;
            max-width: 100%;
            box-sizing: border-box;
        }
        .analysis-summary li {
            margin: 5px 0;
            word-wrap: break-word;
            overflow-wrap: break-word;
            max-width: 100%;
        }
        .qc-metrics, .steps-details, .visualization, .step-plots, .markers-table, .diagnosis {
            max-width: 100%;
            box-sizing: border-box;
            overflow-x: hidden;
            word-wrap: break-word;
            overflow-wrap: break-word;
        }
        .diagnosis div {
            max-width: 100%;
            word-wrap: break-word;
            overflow-wrap: break-word;
            white-space: pre-wrap;
        }
        .visualization, .step-plots {
            max-width: 100%;
            box-sizing: border-box;
            overflow-x: hidden;
        }
        .visualization img, .step-plots img {
            max-width: 100%;
            height: auto;
            border: 1px solid #ddd;
            border-radius: 4px;
            margin: 10px 0;
            display: block;
        }
        .step-plots > div {
            max-width: 100%;
            box-sizing: border-box;
            overflow-x: hidden;
            word-wrap: break-word;
        }
        .markers-table {
            overflow-x: auto;
            max-width: 100%;
            box-sizing: border-box;
            margin: 10px 0;
        }
        .markers-table table {
            width: 100%;
            max-width: 100%;
            border-collapse: collapse;
            margin: 0;
            table-layout: auto;
        }
        .markers-table th, .markers-table td {
            word-wrap: break-word;
            overflow-wrap: break-word;
            max-width: 200px;
        }
        .markers-table th, .markers-table td {
            border: 1px solid #ddd;
            padding: 8px;
            text-align: left;
        }
        .markers-table th {
            background: #f5f5f5;
            font-weight: bold;
        }
        .think-card {
            background: #f1f8e9 !important;
        }
        .think-process {
            margin-bottom: 10px;
        }
        .think-header {
            background: #e8f5e9;
            padding: 10px 15px;
            border-radius: 4px;
            cursor: pointer;
            display: flex;
            align-items: center;
            gap: 10px;
            user-select: none;
            transition: background 0.2s;
        }
        .think-header:hover {
            background: #c8e6c9;
        }
        .think-icon {
            font-size: 18px;
        }
        .think-title {
            flex: 1;
            font-weight: bold;
            color: #2e7d32;
        }
        .think-toggle {
            color: #666;
            font-size: 12px;
        }
        .think-content {
            margin-top: 10px;
            padding: 15px;
            background: #fff;
            border: 1px solid #ddd;
            border-radius: 4px;
            white-space: pre-wrap;
            font-family: 'Courier New', monospace;
            font-size: 13px;
            line-height: 1.6;
            color: #333;
            max-height: 500px;
            overflow-y: auto;
        }
        .final-answer {
            margin-top: 10px;
            padding: 10px;
        }
        .test-data-selection {
            background: #f1f8e9 !important;
            max-width: 100%;
            box-sizing: border-box;
            overflow-x: hidden;
            word-wrap: break-word;
            overflow-wrap: break-word;
        }
        .test-data-selection h3 {
            margin-top: 0;
            color: #4CAF50;
            word-wrap: break-word;
        }
        .test-data-selection div[onclick] {
            border: 1px solid #ddd;
            border-radius: 4px;
            padding: 15px;
            cursor: pointer;
            transition: background 0.2s;
            max-width: 100%;
            box-sizing: border-box;
            word-wrap: break-word;
            overflow-wrap: break-word;
        }
        .test-data-selection div[onclick]:hover {
            background: #f5f5f5;
            border-color: #4CAF50;
        }
        .dataset-card {
            max-width: 100%;
            box-sizing: border-box;
            word-wrap: break-word;
            overflow-wrap: break-word;
        }
    </style>
</head>
<body>
    <div class="container">
        <div class="panel chat-panel">
            <h2>💬 对话界面</h2>
            <div id="chatArea" class="chat-area"></div>
            <div class="input-area">
                <input type="text" id="messageInput" placeholder="输入消息或上传文件进行分析..." />
                <input type="file" id="fileInput" accept=".h5ad,.mtx,.tsv,.csv" multiple />
                <button id="sendBtn" onclick="sendMessage()">发送</button>
            </div>
            <div id="fileInfo" class="file-info" style="display:none;"></div>
        </div>
    </div>

    <script>
        // 文件上下文管理（记住已上传的文件）
        let uploadedFilesContext = [];
        
        // 文件选择（支持多文件）
        let selectedFiles = [];
        document.getElementById('fileInput').addEventListener('change', function(e) {
            const files = Array.from(e.target.files);
            if (files.length > 0) {
                selectedFiles = files;
                const fileList = files.map(f => `${f.name} (${(f.size / 1024 / 1024).toFixed(2)} MB)`).join('<br>');
                document.getElementById('fileInfo').style.display = 'block';
                document.getElementById('fileInfo').innerHTML = `📁 已选择 ${files.length} 个文件:<br>${fileList}`;
            }
        });

        // 发送消息
        async function sendMessage() {
            const input = document.getElementById('messageInput');
            const message = input.value.trim();
            const btn = document.getElementById('sendBtn');
            
            if (!message && selectedFiles.length === 0) {
                alert('请输入消息或选择文件');
                return;
            }

            btn.disabled = true;
            const fileNames = selectedFiles.length > 0 ? selectedFiles.map(f => f.name).join(', ') : '';
            addMessage('user', message || (fileNames ? `上传文件: ${fileNames}` : ''));

            try {
                let uploadedFiles = [];
                
                // 如果有新选择的文件，先上传所有文件
                if (selectedFiles.length > 0) {
                    for (const file of selectedFiles) {
                        const formData = new FormData();
                        formData.append('file', file);
                        
                        const uploadRes = await fetch('/api/upload', {
                            method: 'POST',
                            body: formData
                        });
                        
                        if (!uploadRes.ok) {
                            throw new Error(`文件上传失败: ${file.name}`);
                        }
                        
                        const uploadData = await uploadRes.json();
                        uploadedFiles.push(uploadData);
                        // 添加到上下文
                        uploadedFilesContext.push(uploadData);
                        addMessage('assistant', `✅ 文件上传成功: ${uploadData.file_name}`);
                    }
                } else if (uploadedFilesContext.length > 0) {
                    // 如果没有新文件，使用上下文中的文件
                    uploadedFiles = uploadedFilesContext;
                    addMessage('assistant', `📁 使用已上传的文件: ${uploadedFiles.map(f => f.file_name).join(', ')}`);
                }

                // 发送聊天请求
                const response = await fetch('/api/chat', {
                    method: 'POST',
                    headers: {
                        'Content-Type': 'application/json',
                    },
                    body: JSON.stringify({
                        message: message || (uploadedFiles.length > 0 ? '分析这个文件' : ''),
                        history: [],
                        uploaded_files: uploadedFiles
                    })
                });

                if (!response.ok) {
                    throw new Error(`请求失败: ${response.status}`);
                }

                const contentType = response.headers.get('content-type');
                
                if (contentType && contentType.includes('application/json')) {
                    const data = await response.json();
                    
                    if (data.type === 'workflow_config') {
                        // 执行工作流
                        addMessage('assistant', '🚀 开始执行分析流程...');
                        await executeWorkflow(data.workflow_data, data.file_paths);
                    } else {
                        addMessage('assistant', JSON.stringify(data, null, 2));
                    }
                } else {
                    // 流式响应（支持 think 过程提取）
                    const reader = response.body.getReader();
                    const decoder = new TextDecoder();
                    let fullText = '';
                    let thinkBuffer = '';
                    let isThinking = false;
                    let hasThinkBlock = false;
                    let finalAnswer = '';
                    let thinkStartIndex = -1;
                    let datasetsJson = null;
                    
                    while (true) {
                        const { done, value } = await reader.read();
                        if (done) break;
                        
                        const chunk = decoder.decode(value);
                        fullText += chunk;
                        
                        // 检查是否包含数据集 JSON（测试数据选择响应）
                        // 使用非贪婪匹配，但需要匹配多行（因为 JSON 可能跨行）
                        const datasetsMatch = fullText.match(/<!-- DATASETS_JSON: (\[[\s\S]*?\]) -->/);
                        if (datasetsMatch && !datasetsJson) {
                            try {
                                // JSON 中的换行符已被替换为空格，直接解析即可
                                datasetsJson = JSON.parse(datasetsMatch[1]);
                            } catch (e) {
                                console.error('解析数据集JSON失败:', e, datasetsMatch[1].substring(0, 100));
                            }
                        }
                        
                        // 检测 think 开始标签（支持多种格式）
                        const thinkStartPatterns = [
                            /<think>/i,
                            /<think>/i,
                            /<reasoning>/i,
                            /<thought>/i,
                            /<thinking>/i
                        ];
                        
                        for (const pattern of thinkStartPatterns) {
                            const match = fullText.match(pattern);
                            if (match && !hasThinkBlock) {
                            isThinking = true;
                            hasThinkBlock = true;
                                thinkStartIndex = match.index + match[0].length;
                            // 创建 think 卡片
                            if (!document.querySelector('.think-card:last-child .think-process')) {
                                createThinkCard();
                                }
                                break;
                            }
                        }
                        
                        // 检测 think 结束标签
                        const thinkEndPatterns = [
                            /<\/think>/i,
                            /<\/redacted_reasoning>/i,
                            /<\/reasoning>/i,
                            /<\/thought>/i,
                            /<\/thinking>/i
                        ];
                        
                        for (const pattern of thinkEndPatterns) {
                            const match = fullText.match(pattern);
                            if (match && isThinking) {
                                // 提取 think 内容
                                thinkBuffer = fullText.substring(thinkStartIndex, match.index);
                            updateThinkContent(thinkBuffer);
                            isThinking = false;
                            
                                // 提取 think 标签之后的内容作为最终答案
                                const afterThinkIndex = match.index + match[0].length;
                                finalAnswer = fullText.substring(afterThinkIndex);
                                if (finalAnswer.trim()) {
                                    updateLastMessage('assistant', finalAnswer.trim());
                            }
                                break;
                        }
                        }
                        
                        // 更新显示
                        if (isThinking) {
                            // 在 think 块中，更新 think 内容
                            if (thinkStartIndex >= 0) {
                                thinkBuffer = fullText.substring(thinkStartIndex);
                            updateThinkContent(thinkBuffer);
                            }
                        } else if (hasThinkBlock && !isThinking) {
                            // think 块已结束，更新最终答案
                            if (finalAnswer) {
                                updateLastMessage('assistant', finalAnswer);
                            }
                        } else {
                            // 没有 think 块，直接更新消息
                            // 在流式响应过程中，先显示文本内容（去除 JSON 注释）
                            const cleanText = fullText.replace(/<!-- DATASETS_JSON: \[[\s\S]*?\] -->/g, '').trim();
                            updateLastMessage('assistant', cleanText);
                        }
                    }
                    
                    // 流式响应结束后，如果检测到数据集信息，替换为选择界面
                    if (datasetsJson && datasetsJson.length > 0) {
                        // 移除 JSON 注释，只保留用户友好的文本
                        const cleanText = fullText.replace(/<!-- DATASETS_JSON: \[[\s\S]*?\] -->/g, '').trim();
                        // 移除之前的普通消息
                        const lastMessage = chatArea.querySelector('.message.assistant:last-child');
                        if (lastMessage && !lastMessage.classList.contains('test-data-selection')) {
                            lastMessage.remove();
                        }
                        // 显示选择界面
                        if (!document.querySelector('.test-data-selection')) {
                            displayTestDataSelection(cleanText, datasetsJson);
                        }
                    }
                }
            } catch (error) {
                addMessage('error', `❌ 错误: ${error.message}`);
                console.error(error);
            } finally {
                btn.disabled = false;
                input.value = '';
                // 不清空 selectedFiles，保留文件选择
                // 但清空文件输入框，允许用户重新选择
                document.getElementById('fileInput').value = '';
                // 如果有上下文文件，显示提示
                if (uploadedFilesContext.length > 0) {
                    document.getElementById('fileInfo').style.display = 'block';
                    document.getElementById('fileInfo').innerHTML = `📁 已上传 ${uploadedFilesContext.length} 个文件，可直接输入需求继续分析`;
                } else {
                    document.getElementById('fileInfo').style.display = 'none';
                }
            }
        }

        // 执行工作流
        async function executeWorkflow(workflowData, filePaths) {
            try {
                const response = await fetch('/api/execute', {
                    method: 'POST',
                    headers: {
                        'Content-Type': 'application/json',
                    },
                    body: JSON.stringify({
                        workflow_data: workflowData,
                        file_paths: filePaths
                    })
                });

                const data = await response.json();
                
                if (data.status === 'success') {
                    // 美化显示分析结果
                    displayAnalysisResult(data);
                } else {
                    addMessage('error', `❌ 分析失败: ${data.error || '未知错误'}`);
                }
            } catch (error) {
                addMessage('error', `❌ 执行错误: ${error.message}`);
            }
        }
        
        // 美化显示分析结果
        function displayAnalysisResult(data) {
            const resultDiv = document.createElement('div');
            resultDiv.className = 'message assistant analysis-result';
            
            let html = '<div class="analysis-summary">';
            html += '<h3>✅ 分析完成</h3>';
            
            // QC 指标
            if (data.qc_metrics) {
                html += '<div class="qc-metrics">';
                html += '<h4>📊 质量控制指标</h4>';
                html += '<ul>';
                html += `<li>原始细胞数: <strong>${data.qc_metrics.raw_cells || 'N/A'}</strong></li>`;
                html += `<li>原始基因数: <strong>${data.qc_metrics.raw_genes || 'N/A'}</strong></li>`;
                if (data.qc_metrics.filtered_cells) {
                    html += `<li>过滤后细胞数: <strong>${data.qc_metrics.filtered_cells}</strong></li>`;
                }
                if (data.qc_metrics.filtered_genes) {
                    html += `<li>过滤后基因数: <strong>${data.qc_metrics.filtered_genes}</strong></li>`;
                }
                html += '</ul>';
                html += '</div>';
            }
            
            // 步骤详情
            if (data.steps_details && data.steps_details.length > 0) {
                html += '<div class="steps-details">';
                html += '<h4>📋 执行步骤</h4>';
                html += '<ul>';
                data.steps_details.forEach(step => {
                    const stepName = step.name || step.tool_id || '未知步骤';
                    const stepSummary = step.summary || '完成';
                    html += `<li><strong>${stepName}</strong>: ${stepSummary}</li>`;
                });
                html += '</ul>';
                html += '</div>';
            }
            
            // 可视化图片（只显示步骤的图片，避免与 final_plot 重复）
            if (data.steps_details) {
                const plotSteps = data.steps_details.filter(s => s.plot);
                if (plotSteps.length > 0) {
                    html += '<div class="step-plots">';
                    html += '<h4>📈 可视化结果</h4>';
                    plotSteps.forEach(step => {
                        let plotUrl = step.plot;
                        if (!plotUrl.startsWith('http') && !plotUrl.startsWith('/')) {
                            // 如果路径包含 results，直接使用
                            if (plotUrl.includes('results/')) {
                                plotUrl = '/' + plotUrl;
                            } else {
                                plotUrl = '/results/' + plotUrl;
                            }
                        }
                        html += `<div style="margin: 10px 0;">`;
                        html += `<strong>${step.name || step.tool_id}</strong><br>`;
                        html += `<img src="${plotUrl}" alt="${step.name}" style="max-width: 100%; border-radius: 4px; margin-top: 10px;" onerror="this.style.display='none';">`;
                        html += `</div>`;
                    });
                    html += '</div>';
                } else if (data.final_plot) {
                    // 如果没有步骤图片，使用 final_plot（向后兼容）
                    html += '<div class="visualization">';
                    html += '<h4>📈 可视化结果</h4>';
                    let plotUrl = data.final_plot;
                    if (!plotUrl.startsWith('http') && !plotUrl.startsWith('/')) {
                        if (plotUrl.includes('results/')) {
                            plotUrl = '/' + plotUrl;
                        } else {
                            plotUrl = '/results/' + plotUrl;
                        }
                    }
                    html += `<img src="${plotUrl}" alt="Visualization" style="max-width: 100%; border-radius: 4px; margin-top: 10px;" onerror="this.style.display='none'; this.nextElementSibling.style.display='block';"><p style="display:none; color: #999;">图片加载失败: ${plotUrl}</p>`;
                    html += '</div>';
                }
            } else if (data.final_plot) {
                // 如果没有 steps_details，使用 final_plot
                html += '<div class="visualization">';
                html += '<h4>📈 可视化结果</h4>';
                let plotUrl = data.final_plot;
                if (!plotUrl.startsWith('http') && !plotUrl.startsWith('/')) {
                    if (plotUrl.includes('results/')) {
                        plotUrl = '/' + plotUrl;
                    } else {
                        plotUrl = '/results/' + plotUrl;
                    }
                }
                html += `<img src="${plotUrl}" alt="Visualization" style="max-width: 100%; border-radius: 4px; margin-top: 10px;" onerror="this.style.display='none'; this.nextElementSibling.style.display='block';"><p style="display:none; color: #999;">图片加载失败: ${plotUrl}</p>`;
                html += '</div>';
            }
            
            // Marker 基因表格（如果有）
            const markersStep = data.steps_details?.find(s => s.name === 'local_markers' || s.tool_id === 'local_markers');
            if (markersStep && markersStep.details) {
                html += '<div class="markers-table">';
                html += '<h4>🧬 Marker 基因</h4>';
                // 直接显示 HTML 表格
                html += markersStep.details;
                html += '</div>';
            }
            
            // 诊断信息
            if (data.diagnosis) {
                html += '<div class="diagnosis">';
                html += '<h4>💡 分析诊断</h4>';
                html += `<div style="white-space: pre-wrap;">${data.diagnosis}</div>`;
                html += '</div>';
            }
            
            html += '</div>';
            resultDiv.innerHTML = html;
            
            const chatArea = document.getElementById('chatArea');
            chatArea.appendChild(resultDiv);
            chatArea.scrollTop = chatArea.scrollHeight;
        }

        // 添加消息
        function addMessage(role, content) {
            const chatArea = document.getElementById('chatArea');
            const msgDiv = document.createElement('div');
            msgDiv.className = `message ${role}`;
            msgDiv.textContent = content;
            chatArea.appendChild(msgDiv);
            chatArea.scrollTop = chatArea.scrollHeight;
        }

        // 更新最后一条消息
        function updateLastMessage(role, content) {
            const chatArea = document.getElementById('chatArea');
            const messages = chatArea.querySelectorAll('.message');
            if (messages.length > 0 && messages[messages.length - 1].classList.contains(role)) {
                const lastMsg = messages[messages.length - 1];
                // 如果已经有 think 卡片，更新最终答案部分
                const finalAnswerDiv = lastMsg.querySelector('.final-answer');
                if (finalAnswerDiv) {
                    finalAnswerDiv.textContent = content;
                } else {
                    lastMsg.textContent = content;
                }
            } else {
                addMessage(role, content);
            }
            chatArea.scrollTop = chatArea.scrollHeight;
        }
        
        // 创建 think 卡片
        function createThinkCard() {
            const chatArea = document.getElementById('chatArea');
            const thinkCard = document.createElement('div');
            thinkCard.className = 'message assistant think-card';
            thinkCard.innerHTML = `
                <div class="think-process">
                    <div class="think-header" onclick="toggleThink(this)">
                        <span class="think-icon">🤔</span>
                        <span class="think-title">思考过程</span>
                        <span class="think-toggle">▼</span>
                    </div>
                    <div class="think-content" style="display: none;"></div>
                </div>
                <div class="final-answer"></div>
            `;
            chatArea.appendChild(thinkCard);
            chatArea.scrollTop = chatArea.scrollHeight;
        }
        
        // 更新 think 内容
        function updateThinkContent(content) {
            const chatArea = document.getElementById('chatArea');
            const thinkCards = chatArea.querySelectorAll('.think-card');
            if (thinkCards.length > 0) {
                const lastCard = thinkCards[thinkCards.length - 1];
                const thinkContentDiv = lastCard.querySelector('.think-content');
                if (thinkContentDiv) {
                    thinkContentDiv.textContent = content;
                }
            }
        }
        
        // 切换 think 卡片展开/折叠
        function toggleThink(header) {
            const thinkCard = header.closest('.think-process');
            const content = thinkCard.querySelector('.think-content');
            const toggle = header.querySelector('.think-toggle');
            
            if (content.style.display === 'none') {
                content.style.display = 'block';
                toggle.textContent = '▲';
            } else {
                content.style.display = 'none';
                toggle.textContent = '▼';
            }
        }
        
        // 显示测试数据选择界面
        function displayTestDataSelection(messageText, datasets) {
            const chatArea = document.getElementById('chatArea');
            
            // 检查是否已经显示过选择界面
            const existing = document.querySelector('.test-data-selection');
            if (existing) {
                // 如果已存在，更新它而不是创建新的
                return;
            }
            
            // 移除之前的普通消息（如果有）
            const lastMessage = chatArea.querySelector('.message.assistant:last-child');
            if (lastMessage && !lastMessage.classList.contains('test-data-selection')) {
                lastMessage.remove();
            }
            
            const selectionDiv = document.createElement('div');
            selectionDiv.className = 'message assistant test-data-selection';
            
            let html = '<div style="padding: 15px; max-width: 100%; box-sizing: border-box;">';
            html += '<h3 style="margin-top: 0; color: #4CAF50; margin-bottom: 10px; word-wrap: break-word;">📊 选择测试数据集</h3>';
            
            // 显示消息文本（去除数据集列表部分）
            const lines = messageText.split('\\n');
            const messageLine = lines.find(line => line.includes('检测到') || line.includes('请选择'));
            if (messageLine) {
                html += '<p style="margin-bottom: 15px; color: #333; word-wrap: break-word; max-width: 100%;">' + messageLine + '</p>';
            }
            
            // 显示数据集选择卡片
            html += '<div style="display: flex; flex-direction: column; gap: 10px; margin-bottom: 15px; max-width: 100%; box-sizing: border-box;">';
            datasets.forEach(dataset => {
                html += `<div class="dataset-card" 
                             style="border: 2px solid #ddd; border-radius: 8px; padding: 15px; cursor: pointer; transition: all 0.2s; background: white; max-width: 100%; box-sizing: border-box; word-wrap: break-word; overflow-wrap: break-word;" 
                             onmouseover="this.style.borderColor='#4CAF50'; this.style.boxShadow='0 2px 8px rgba(76,175,80,0.2)'" 
                             onmouseout="this.style.borderColor='#ddd'; this.style.boxShadow='none'"
                             onclick="selectTestDataset('${dataset.id}', '${dataset.name}')">`;
                html += `<div style="display: flex; align-items: center; gap: 10px; margin-bottom: 8px; max-width: 100%; flex-wrap: wrap;">`;
                html += `<span style="font-size: 24px; flex-shrink: 0;">📦</span>`;
                html += `<strong style="color: #4CAF50; font-size: 18px; word-wrap: break-word; flex: 1; min-width: 0;">${dataset.name}</strong>`;
                html += `</div>`;
                html += `<p style="margin: 0; color: #666; font-size: 14px; word-wrap: break-word; max-width: 100%;">${dataset.description}</p>`;
                html += `<div style="margin-top: 8px; font-size: 12px; color: #999; word-wrap: break-word; max-width: 100%;">ID: <code style="word-break: break-all;">${dataset.id}</code></div>`;
                html += '</div>';
            });
            html += '</div>';
            
            html += '<p style="margin-top: 10px; color: #666; font-size: 14px; font-style: italic; word-wrap: break-word; max-width: 100%;">💡 点击上面的数据集卡片选择，或上传您自己的数据文件。</p>';
            html += '</div>';
            
            selectionDiv.innerHTML = html;
            chatArea.appendChild(selectionDiv);
            chatArea.scrollTop = chatArea.scrollHeight;
        }
        
        // 选择测试数据集
        async function selectTestDataset(datasetId, datasetName) {
            addMessage('user', `使用测试数据集: ${datasetName} (${datasetId})`);
            addMessage('assistant', `正在使用测试数据集 ${datasetName} 执行分析...`);
            
            try {
                const response = await fetch('/api/chat', {
                    method: 'POST',
                    headers: {
                        'Content-Type': 'application/json',
                    },
                    body: JSON.stringify({
                        message: `使用测试数据集 ${datasetId} 执行完整的单细胞转录组分析流程`,
                        history: [],
                        uploaded_files: [],
                        test_dataset_id: datasetId
                    })
                });
                
                // 处理响应
                const contentType = response.headers.get('content-type');
                if (contentType && contentType.includes('application/json')) {
                    const data = await response.json();
                    if (data.type === 'workflow_config') {
                        addMessage('assistant', '🚀 开始执行分析流程...');
                        await executeWorkflow(data.workflow_data, data.file_paths || []);
                    } else {
                        addMessage('assistant', JSON.stringify(data, null, 2));
                    }
                } else {
                    // 流式响应
                    const reader = response.body.getReader();
                    const decoder = new TextDecoder();
                    let fullText = '';
                    
                    while (true) {
                        const { done, value } = await reader.read();
                        if (done) break;
                        const chunk = decoder.decode(value);
                        fullText += chunk;
                        updateLastMessage('assistant', fullText);
                    }
                }
            } catch (error) {
                addMessage('error', `❌ 错误: ${error.message}`);
                console.error(error);
            }
        }
        
        // 全局函数，供 HTML 调用
        window.toggleThink = toggleThink;
        window.selectTestDataset = selectTestDataset;

        // 日志功能已移除

        // 回车发送
        document.getElementById('messageInput').addEventListener('keypress', function(e) {
            if (e.key === 'Enter') {
                sendMessage();
            }
        });

        // 初始化完成
    </script>
</body>
</html>
    """
    return HTMLResponse(content=html_content)


@app.post("/api/upload")
async def upload_file(
    files: List[UploadFile] = File(...),  # 🔥 必须 List：单数 file 只会接收第一个，多文件会被静默丢弃
    owner_id: str = Depends(get_current_owner_id),
    db: Session = Depends(get_db_session),
):
    """
    文件上传接口（支持多文件上传）
    Phase 4: 身份由 Authorization / X-Guest-UUID 解析为 owner_id，路径按 owner_id 隔离并写入 Asset 表。
    """
    from gibh_agent.db.models import Asset as AssetModel

    try:
        if not files or len(files) == 0:
            raise HTTPException(status_code=400, detail="No files provided")
        
        # 限制文件数量
        if len(files) > 20:
            raise HTTPException(status_code=400, detail="一次最多上传20个文件")
        
        # Phase 4: 按 owner_id 隔离目录，用时间戳子目录避免重名
        batch_id = datetime.now().strftime("%Y%m%d_%H%M%S")
        user_dir = UPLOAD_DIR / owner_id / batch_id
        user_dir.mkdir(parents=True, exist_ok=True)
        user_id = owner_id  # 兼容下方日志与返回
        session_id = batch_id
        
        logger.info(f"📤 收到文件上传: {len(files)} 个文件 (owner_id: {owner_id})")
        
        # 检测是否是10x Genomics文件（matrix.mtx, barcodes.tsv, features.tsv）
        is_10x_data = False
        tenx_files = []
        other_files = []
        
        for file in files:
            # 🔒 安全：清理文件名
            if not file.filename:
                raise HTTPException(status_code=400, detail="文件名不能为空")
            
            original_filename = file.filename
            safe_filename = sanitize_filename(original_filename)
            
            # 🔒 安全：验证文件扩展名
            file_ext = Path(safe_filename).suffix.lower()
            if file_ext and file_ext not in ALLOWED_EXTENSIONS:
                raise HTTPException(
                    status_code=400,
                    detail=f"不允许的文件类型: {file_ext}。允许的类型: {', '.join(ALLOWED_EXTENSIONS)}"
                )
            
            # 更新文件名为安全版本
            file.filename = safe_filename
            
            filename_lower = safe_filename.lower()
            # 🔥 修复：支持 genes.tsv（旧版10x格式）和 features.tsv（新版10x格式）
            if any(pattern in filename_lower for pattern in ['matrix.mtx', 'barcodes.tsv', 'features.tsv', 'genes.tsv']):
                is_10x_data = True
                tenx_files.append(file)
            else:
                other_files.append(file)
        
        uploaded_results = []
        
        # Phase 4: user_dir 已在上面按 owner_id 创建
        logger.info(f"📁 用户目录: {user_dir}")
        
        # 如果是10x数据，创建子目录并保存
        if is_10x_data and len(tenx_files) >= 2:  # 至少需要2个文件（通常是matrix + barcodes/features）
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            tenx_dir = user_dir / f"10x_data_{timestamp}"
            tenx_dir.mkdir(exist_ok=True)
            
            logger.info(f"📁 检测到10x数据，创建目录: {tenx_dir}")
            
            # 保存10x文件到子目录
            for file in tenx_files:
                # 🔒 安全：验证文件路径
                file_path = tenx_dir / file.filename
                try:
                    file_path = validate_file_path(file_path, UPLOAD_DIR)
                except HTTPException as e:
                    logger.error(f"❌ 文件路径验证失败: {file.filename} -> {e.detail}")
                    raise
                
                # 🔒 安全：检查文件大小
                content = await file.read()
                if len(content) > MAX_FILE_SIZE:
                    raise HTTPException(
                        status_code=413,
                        detail=f"文件 {file.filename} 超过最大大小限制 ({MAX_FILE_SIZE / 1024 / 1024:.0f}MB)"
                    )
                
                # 🔧 修复：确保父目录存在
                file_path.parent.mkdir(parents=True, exist_ok=True)
                
                try:
                    with open(file_path, "wb") as f:
                        f.write(content)
                except PermissionError as e:
                    logger.error(f"❌ 文件写入权限错误: {file_path} -> {e}")
                    raise HTTPException(status_code=500, detail=f"文件保存失败：权限不足 ({file.filename})")
                except OSError as e:
                    logger.error(f"❌ 文件写入系统错误: {file_path} -> {e}")
                    raise HTTPException(status_code=500, detail=f"文件保存失败：{str(e)} ({file.filename})")
                
                logger.info(f"✅ 10x文件保存成功: {file_path}")
                
                # 生成元数据
                try:
                    metadata = file_inspector.generate_metadata(str(file_path.relative_to(UPLOAD_DIR)))
                except Exception as e:
                    logger.warning(f"⚠️ 生成文件元数据失败: {e}")
                    metadata = None
                
                uploaded_results.append({
                    "file_id": str(tenx_dir.relative_to(UPLOAD_DIR)),
                    "file_name": file.filename,
                    "file_path": str(file_path),
                    "file_size": len(content),
                    "metadata": metadata,
                    "is_10x": True,
                    "group_dir": str(tenx_dir.relative_to(UPLOAD_DIR))
                })
                
                # 🔒 Async signing (anti-regression: do not block 200 on failure)
                try:
                    from gibh_agent.core.tasks import sign_uploaded_file_task
                    sign_uploaded_file_task.delay(str(file_path))
                except Exception as e:
                    logger.warning("⚠️ 签名任务入队失败（文件已保存，不影响上传）: %s", e)
            
            # SpatialStructureNormalizer: reorganize loose spatial archives into spatial/ before inspection
            try:
                normalize_session_directory(Path(tenx_dir))
            except Exception as e:
                logger.warning("⚠️ 目录结构规范化失败（不影响上传）: %s", e)
            
            # Phase 4: 10x 文件入库
            if db is not None:
                try:
                    for r in uploaded_results:
                        db.add(AssetModel(owner_id=owner_id, file_name=r.get("file_name", ""), file_path=r.get("file_path", ""), modality=None))
                    db.commit()
                except Exception as e:
                    logger.warning("⚠️ [Upload] 10x Asset 入库失败: %s", e)
            # 返回10x目录路径（强制绝对路径，供下游与前端一致使用）
            abs_tenx = str(tenx_dir.resolve())
            for r in uploaded_results:
                r["file_path"] = r.get("file_path") and _ensure_absolute_upload_path(r["file_path"]) or abs_tenx
            return {
                "status": "success",
                "is_10x_data": True,
                "group_dir": str(tenx_dir.relative_to(UPLOAD_DIR)),
                "files": uploaded_results,
                "file_paths": [abs_tenx],
                "message": f"10x数据已保存到: {tenx_dir.relative_to(UPLOAD_DIR)}"
            }
        
        # 处理其他文件（非10x或单独的10x文件）
        # 🔧 修复：如果只有1个10x文件，也当作普通文件处理
        files_to_process = other_files
        if is_10x_data and len(tenx_files) == 1:
            # 只有1个10x文件，当作普通文件处理
            logger.info(f"⚠️ 只有1个10x文件，当作普通文件处理: {tenx_files[0].filename}")
            files_to_process = other_files + tenx_files
        elif not is_10x_data:
            # 不是10x数据，处理所有文件
            files_to_process = other_files + tenx_files
        
        for file in files_to_process:
            # 🔥 多用户支持：文件保存到用户目录
            file_path = user_dir / file.filename
            try:
                file_path = validate_file_path(file_path, UPLOAD_DIR)
            except HTTPException as e:
                logger.error(f"❌ 文件路径验证失败: {file.filename} -> {e.detail}")
                raise
            
            # 🔒 安全：检查文件大小
            content = await file.read()
            if len(content) > MAX_FILE_SIZE:
                raise HTTPException(
                    status_code=413,
                    detail=f"文件 {file.filename} 超过最大大小限制 ({MAX_FILE_SIZE / 1024 / 1024:.0f}MB)"
                )
            
            # 🔧 修复：确保父目录存在
            file_path.parent.mkdir(parents=True, exist_ok=True)
            
            try:
                with open(file_path, "wb") as f:
                    f.write(content)
            except PermissionError as e:
                logger.error(f"❌ 文件写入权限错误: {file_path} -> {e}")
                raise HTTPException(status_code=500, detail=f"文件保存失败：权限不足 ({file.filename})")
            except OSError as e:
                logger.error(f"❌ 文件写入系统错误: {file_path} -> {e}")
                raise HTTPException(status_code=500, detail=f"文件保存失败：{str(e)} ({file.filename})")
            
            logger.info(f"✅ 文件保存成功: {file_path}")
            
            # 生成元数据
            try:
                metadata = file_inspector.generate_metadata(file.filename)
                if metadata:
                    logger.info(f"📊 文件元数据已生成: {metadata.get('file_type', 'unknown')}")
            except Exception as e:
                logger.warning(f"⚠️ 生成文件元数据失败: {e}")
                metadata = None
            
            uploaded_results.append({
                "file_id": file.filename,
                "file_name": file.filename,
                "file_path": str(file_path),
                "file_size": len(content),
                "metadata": metadata,
                "is_10x": False
            })
            
            # 🔒 Async signing (anti-regression: do not block 200 on failure)
            try:
                from gibh_agent.core.tasks import sign_uploaded_file_task
                sign_uploaded_file_task.delay(str(file_path))
            except Exception as e:
                logger.warning("⚠️ 签名任务入队失败（文件已保存，不影响上传）: %s", e)
        
        # Universal Ingestion: unpack archives (zip/tar.gz) and get effective root
        effective_root = None
        replaced_archive_paths: List[Path] = []
        try:
            unpacked, effective_root, replaced_archive_paths = universal_unpack(
                Path(user_dir), remove_archive=False
            )
            if unpacked and effective_root is not None:
                # Visium 分体上传：.h5 在 user_dir，spatial 在 effective_root；复制 .h5 进 effective_root 以便识别为 Spatial
                effective_path = Path(effective_root)
                has_h5 = any(
                    f.is_file() and (f.name.endswith("_feature_bc_matrix.h5") or f.name == "filtered_feature_bc_matrix.h5" or f.name == "raw_feature_bc_matrix.h5")
                    for f in effective_path.iterdir()
                )
                if not has_h5:
                    for f in Path(user_dir).iterdir():
                        if f.is_file() and (f.name.endswith("_feature_bc_matrix.h5") or f.name in ("filtered_feature_bc_matrix.h5", "raw_feature_bc_matrix.h5")):
                            dest = effective_path / f.name
                            if not dest.exists():
                                try:
                                    import shutil
                                    shutil.copy2(str(f), str(dest))
                                    logger.info("✅ [UniversalIngestion] 已将 matrix .h5 复制到解压目录: %s -> %s", f.name, effective_path)
                                except Exception as e:
                                    logger.warning("⚠️ 复制 .h5 到解压目录失败: %s", e)
                            break
                # 先规范化布局（创建 spatial/、从父目录补 .h5），再检测模态
                try:
                    normalize_session_directory(Path(effective_root))
                except Exception as e:
                    logger.warning("⚠️ 解压目录规范化失败: %s", e)
                modality, payload = detect_dominant_modality(effective_root)
                if modality != "unknown" and payload is not None:
                    if modality == "Spatial" and payload is not None:
                        try:
                            normalize_session_directory(Path(payload))
                        except Exception as e:
                            logger.warning("⚠️ Spatial 目录规范化失败: %s", e)
                    inner_paths = paths_for_response(modality, payload, UPLOAD_DIR)
                    if inner_paths:
                        # 🔥 合并而非覆盖：保留同批上传的独立文件（如 Visium 的 .h5），再追加嗅探/解压得到的路径
                        unpack_entries = []
                        for p in inner_paths:
                            try:
                                size = p.stat().st_size if p.is_file() else 0
                            except OSError:
                                size = 0
                            unpack_entries.append({
                                "file_id": p.name,
                                "file_name": p.name,
                                "file_path": str(p.resolve()),
                                "file_size": size,
                                "metadata": None,
                                "is_10x": False,
                            })
                        merged = list(uploaded_results)
                        seen_norm = {
                            _ensure_absolute_upload_path(r.get("file_path") or "").rstrip("/")
                            for r in merged
                            if r.get("file_path")
                        }
                        for ent in unpack_entries:
                            rawp = ent.get("file_path") or ""
                            ap = _ensure_absolute_upload_path(rawp).rstrip("/") if rawp else ""
                            if ap and ap not in seen_norm:
                                seen_norm.add(ap)
                                merged.append(ent)
                        effective_path = Path(effective_root)
                        abs_eff = str(effective_path.resolve()).rstrip("/")
                        if abs_eff and abs_eff not in seen_norm:
                            merged.append({
                                "file_id": effective_path.name,
                                "file_name": effective_path.name,
                                "file_path": str(effective_path.resolve()),
                                "file_size": 0,
                                "metadata": None,
                                "is_10x": False,
                            })
                            seen_norm.add(abs_eff)
                        uploaded_results = merged
                        logger.info(
                            "✅ [UniversalIngestion] 已合并解压/嗅探路径与原始上传项，共 %s 条",
                            len(uploaded_results),
                        )
                elif modality == "unknown" and effective_root is not None:
                    # 🔥 多文件修复：解压目录追加到列表，不覆盖（保留已上传的 .h5 等，避免多文件变单文件）
                    effective_path = Path(effective_root)
                    abs_eff = str(effective_path.resolve())
                    if not any((r.get("file_path") or "").rstrip("/") == abs_eff.rstrip("/") for r in uploaded_results):
                        uploaded_results.append({
                            "file_id": effective_path.name,
                            "file_name": effective_path.name,
                            "file_path": abs_eff,
                            "file_size": 0,
                            "metadata": None,
                            "is_10x": False,
                        })
                    try:
                        rel = str(effective_path.relative_to(UPLOAD_DIR))
                    except ValueError:
                        rel = effective_path.name
                    logger.info("✅ [UniversalIngestion] 已追加解压目录到返回列表（供下游使用）: %s", rel)
        except Exception as e:
            logger.warning("⚠️ Universal unpack / Modality Sniffer 失败（不影响上传）: %s", e)

        # 若已解压但列表里尚未包含解压目录（例如 modality 已识别而 paths_for_response 为空），补一条以免后续剔除归档后无有效路径
        if replaced_archive_paths and effective_root is not None:
            eff_p = Path(effective_root).resolve()
            abs_eff = str(eff_p)

            def _upload_results_contain_path(rows: List[Dict[str, Any]], target_abs: str) -> bool:
                for r in rows:
                    fp = r.get("file_path") or ""
                    if not fp:
                        continue
                    try:
                        if str(Path(_ensure_absolute_upload_path(fp)).resolve()) == target_abs:
                            return True
                    except Exception:
                        continue
                return False

            if not _upload_results_contain_path(uploaded_results, abs_eff):
                uploaded_results.append({
                    "file_id": eff_p.name,
                    "file_name": eff_p.name,
                    "file_path": abs_eff,
                    "file_size": 0,
                    "metadata": None,
                    "is_10x": False,
                })
                logger.info("✅ [UniversalIngestion] 已补充解压目录到上传结果（paths_for_response 未返回子路径时）: %s", abs_eff)

        # 归档替换：成功解压后从返回列表剔除原始压缩包路径，避免前端/上下文出现「幽灵资产」
        if replaced_archive_paths:
            rep_norm = {str(p.resolve()) for p in replaced_archive_paths}

            def _upload_row_resolved(row: Dict[str, Any]) -> str:
                fp = row.get("file_path") or ""
                if not fp:
                    return ""
                try:
                    return str(Path(_ensure_absolute_upload_path(fp)).resolve())
                except Exception:
                    return str(Path(_ensure_absolute_upload_path(fp)))

            before_n = len(uploaded_results)
            uploaded_results = [r for r in uploaded_results if _upload_row_resolved(r) not in rep_norm]
            if len(uploaded_results) != before_n:
                logger.info(
                    "✅ [UniversalIngestion] 已从上传结果剔除已解压归档 %s 条（幽灵资产）",
                    before_n - len(uploaded_results),
                )
        
        # SpatialStructureNormalizer: reorganize loose spatial.tar.gz + .h5 into Visium layout (spatial/ + matrix)
        try:
            normalize_session_directory(Path(effective_root if effective_root is not None else user_dir))
        except Exception as e:
            logger.warning("⚠️ 目录结构规范化失败（不影响上传）: %s", e)
        
        # 🔥 统一返回格式：强制使用绝对路径，避免底层工具收到相对路径导致 [Errno 2] No such file or directory
        file_paths = []
        file_info = []
        for result in uploaded_results:
            raw = result.get("file_path") or result.get("file_id") or ""
            abs_path = _ensure_absolute_upload_path(raw) if raw else ""
            file_paths.append(abs_path)
            file_info.append({
                "name": result.get("file_name", ""),
                "size": result.get("file_size", 0),
                "path": abs_path
            })
        
        # Phase 4: 非 10x 文件入库（Asset 表与返回前端均使用绝对路径）
        if db is not None:
            try:
                for result in uploaded_results:
                    fp = result.get("file_path") or ""
                    abs_fp = _ensure_absolute_upload_path(fp) if fp else ""
                    db.add(AssetModel(owner_id=owner_id, file_name=result.get("file_name", ""), file_path=abs_fp, modality=None))
                db.commit()
            except Exception as e:
                logger.warning("⚠️ [Upload] Asset 入库失败: %s", e)
        # 🔥 统一返回格式：始终返回一致的 JSON 结构
        response = {
            "status": "success",
            "file_paths": file_paths,  # 文件路径数组（相对路径，包含 owner_id/batch_id）
            "file_info": file_info,    # 文件信息数组
            "count": len(uploaded_results),
            "user_id": user_id,        # 兼容前端
            "session_id": session_id   # 兼容前端
        }
        
        # 如果只有一个文件，添加单个文件的详细信息（向后兼容，统一绝对路径）
        if len(uploaded_results) == 1:
            result = uploaded_results[0]
            response.update({
                "file_id": result["file_id"],
                "file_name": result["file_name"],
                "file_path": file_paths[0] if file_paths else _ensure_absolute_upload_path(result.get("file_path", "")),
                "file_size": result["file_size"],
                "metadata": result.get("metadata")
            })
        else:
            # 多个文件时，添加 files 数组（向后兼容）
            response["files"] = uploaded_results
        
        return response
        
    except HTTPException:
        # 重新抛出 HTTP 异常（保持状态码和详细信息）
        raise
    except Exception as e:
        # 🔧 改进：记录详细错误信息，但返回用户友好的错误消息
        import traceback
        error_detail = f"{type(e).__name__}: {str(e)}"
        logger.error(f"❌ 文件上传失败: {error_detail}", exc_info=True)
        logger.error(f"详细堆栈:\n{traceback.format_exc()}")
        
        # 根据错误类型返回更具体的错误信息
        if "Permission" in error_detail or "permission" in error_detail.lower():
            raise HTTPException(status_code=500, detail="文件上传失败：权限不足，请检查服务器配置")
        elif "No such file" in error_detail or "directory" in error_detail.lower():
            raise HTTPException(status_code=500, detail="文件上传失败：目录不存在，请检查服务器配置")
        elif "disk" in error_detail.lower() or "space" in error_detail.lower():
            raise HTTPException(status_code=500, detail="文件上传失败：磁盘空间不足")
        else:
            # 开发环境返回更详细的错误信息，生产环境返回通用错误
            import os
            if os.getenv("DEBUG", "false").lower() == "true":
                raise HTTPException(status_code=500, detail=f"文件上传失败: {error_detail}")
            else:
                raise HTTPException(status_code=500, detail="文件上传失败，请稍后重试")


@app.post("/api/chat")
async def chat_endpoint(
    req: ChatRequest,
    owner_id: str = Depends(get_current_owner_id),
    db: Session = Depends(get_db_session),
):
    """聊天接口。Phase 4: 身份解析、Session/Message 持久化，AI 消息在流结束后写入，不阻塞 SSE。"""
    _chat_start = time.time()
    logger.info("[Profiler] /api/chat 请求进入 - 耗时: 0.00s")

    from gibh_agent.db.models import Session as SessionModel, Message as MessageModel

    # 🔥 CRITICAL DEBUG: 打印接收到的文件列表
    print(f"🔍 API RECEIVED FILES: {req.uploaded_files}")
    logger.info(f"🔍 [ChatEndpoint] 接收到的文件列表: {req.uploaded_files}")
    
    # #region debug log - entry point
    import json
    import traceback
    import uuid
    # 🔧 修复：使用容器内的日志路径（统一使用 /app/debug.log）
    debug_log_path = Path("/app/debug.log")
    try:
        # 确保目录存在
        debug_log_path.parent.mkdir(parents=True, exist_ok=True)
        with open(debug_log_path, 'a') as f:
            f.write(json.dumps({"location":"server.py:1112","message":"chat_endpoint entry","data":{"agent_is_none":agent is None,"req_message":req.message[:100] if req.message else None,"uploaded_files":req.uploaded_files},"timestamp":int(__import__('time').time()*1000),"sessionId":"debug-session","runId":"run1","hypothesisId":"ENTRY"})+"\n")
    except Exception as log_err:
        pass  # 即使日志写入失败也不影响主流程
    # #endregion
    
    # Phase 4: 会话与身份（新建会话必须立即 commit + refresh，避免“数据黑洞”）
    if not req.session_id:
        req.session_id = str(uuid.uuid4())
        try:
            title = (req.message or "新会话")[:20]
            new_session = SessionModel(id=req.session_id, owner_id=owner_id, title=title)
            db.add(new_session)
            db.commit()
            db.refresh(new_session)
            logger.info(f"✅ [DB] 成功创建新会话: {new_session.id}, Owner: {owner_id}, Title: {new_session.title}")
        except Exception as e:
            logger.warning("⚠️ [Chat] Session 入库失败: %s", e)
            db.rollback()
        logger.info(f"🔑 [ChatEndpoint] 新建 session_id: {req.session_id}")
    else:
        # 校验已有 session 归属
        existing = db.query(SessionModel).filter(SessionModel.id == req.session_id).first()
        if existing and existing.owner_id != owner_id:
            raise HTTPException(status_code=403, detail="无权使用该会话")
    if not req.user_id:
        req.user_id = owner_id
        logger.debug(f"🔑 [ChatEndpoint] user_id = owner_id: {req.user_id}")
    
    if not agent:
        error_msg = "智能体未初始化，请检查配置和日志。可能的原因：1) 配置文件路径错误 2) API Key未设置 3) 依赖包缺失"
        logger.error(error_msg)
        logger.error("请检查终端日志中的详细错误信息")
        return JSONResponse(
            status_code=500,
            content={
                "type": "error",
                "error": error_msg,
                "message": "智能体初始化失败，请查看服务器日志获取详细信息"
            }
        )
    
    # Phase 4: 用户消息入库（大模型生成前）
    try:
        db.add(MessageModel(session_id=req.session_id, role="user", content={"text": req.message or ""}))
        db.commit()
    except Exception as e:
        logger.warning("⚠️ [Chat] 用户消息入库失败: %s", e)
        db.rollback()
    
    # 🔥 SSE 流式传输模式
    if req.stream:
        logger.info(f"🔥 [SSE] 启用流式传输模式")
        try:
            # 🔥 CRITICAL FIX: 转换文件格式（健壮处理）
            uploaded_files = []
            logger.info(f"🔍 [ChatEndpoint] 原始 uploaded_files: {req.uploaded_files}")
            logger.info(f"🔍 [ChatEndpoint] uploaded_files 类型: {type(req.uploaded_files)}, 长度: {len(req.uploaded_files) if req.uploaded_files else 0}")
            
            for i, file_info in enumerate(req.uploaded_files):
                logger.info(f"🔍 [ChatEndpoint] 处理文件 [{i}]: {file_info}, 类型: {type(file_info)}")
                
                # 🔥 CRITICAL: 支持多种格式（dict, Pydantic model, string）
                file_name = ""
                file_path_str = ""
                
                if isinstance(file_info, dict):
                    file_name = file_info.get("file_name") or file_info.get("name", "")
                    file_path_str = file_info.get("file_path") or file_info.get("path", "")
                elif hasattr(file_info, "path"):
                    # Pydantic model
                    file_path_str = file_info.path
                    file_name = getattr(file_info, "name", "") or getattr(file_info, "file_name", "")
                elif isinstance(file_info, str):
                    # 直接是路径字符串
                    file_path_str = file_info
                    file_name = os.path.basename(file_path_str)
                else:
                    logger.warning(f"⚠️ [ChatEndpoint] 未知的文件格式: {type(file_info)}")
                    continue
                
                logger.info(f"🔍 [ChatEndpoint] 提取的文件名: {file_name}, 路径: {file_path_str}")
                
                # 🔥 CRITICAL: 构建文件路径（支持绝对路径和相对路径）
                if file_path_str:
                    file_path = Path(file_path_str)
                    if not file_path.is_absolute():
                        # 相对路径：相对于 UPLOAD_DIR
                        file_path = UPLOAD_DIR / file_path
                    # 如果路径不存在，尝试使用文件名查找
                    if not file_path.exists() and file_name:
                        # 尝试在 UPLOAD_DIR 中查找文件名
                        potential_path = UPLOAD_DIR / file_name
                        if potential_path.exists():
                            file_path = potential_path
                            logger.info(f"✅ [ChatEndpoint] 在 UPLOAD_DIR 中找到文件: {file_path}")
                        else:
                            logger.warning(f"⚠️ [ChatEndpoint] 文件不存在: {file_path}, 尝试: {potential_path}")
                elif file_name:
                    file_path = UPLOAD_DIR / file_name
                else:
                    logger.warning(f"⚠️ [ChatEndpoint] 无法确定文件路径，跳过")
                    continue
                
                # 始终加入文件项（路径存在则用绝对路径，不存在也保留以便进入执行分支，体检时再报错）
                file_dict = {
                    "name": file_name or os.path.basename(str(file_path)),
                    "path": str(file_path.resolve()) if file_path.exists() else str(file_path)
                }
                uploaded_files.append(file_dict)
                if file_path.exists():
                    logger.info(f"✅ [ChatEndpoint] 添加文件: {file_dict}")
                else:
                    logger.warning(f"⚠️ [ChatEndpoint] 路径暂不存在，仍加入列表以便执行分支: {file_path}")
            
            logger.info(f"✅ [ChatEndpoint] 转换后的 uploaded_files: {uploaded_files}")
            logger.info(f"✅ [ChatEndpoint] uploaded_files 数量: {len(uploaded_files)}")
            
            # 🔥 任务2：从 prompt 中提取「系统注入」的历史资产路径，并入 uploaded_files；强制绝对路径
            injected = re.findall(r'\[系统注入：用户选择了历史资产文件：(.*?)\]', (req.message or ""))
            for p in injected:
                p = (p or "").strip()
                if not p:
                    continue
                abs_p = _ensure_absolute_upload_path(p)
                if any((f.get("path") or "").rstrip("/") == abs_p.rstrip("/") for f in uploaded_files):
                    continue
                path_obj = Path(abs_p)
                if path_obj.exists():
                    uploaded_files.append({"name": path_obj.name, "path": abs_p})
                    logger.info("✅ [ChatEndpoint] 注入历史资产路径(绝对): %s", abs_p)
                else:
                    uploaded_files.append({"name": os.path.basename(p), "path": abs_p})
                    logger.warning("⚠️ [ChatEndpoint] 注入路径不存在，仍加入列表: %s", abs_p)
            
            # 创建编排器（传递 upload_dir）
            orchestrator = AgentOrchestrator(agent, upload_dir=str(UPLOAD_DIR))
            
            # 返回 SSE 流式响应（Phase 4: 首包下发 session_id，流结束后写 agent 消息，不阻塞 SSE）
            async def generate_sse():
                original_llm_clients = {}
                override_llm = None
                state_snapshot_for_db = None  # 🔥 有状态应用切片：从流中解析 event: state_snapshot 作为入库 content
                # Phase 4: 第一个事件下发 session_id，供前端持久化
                yield f"event: session\ndata: {json.dumps({'session_id': req.session_id}, ensure_ascii=False)}\n\n"
                try:
                    # 🔥 动态模型路由：前端 model_name 传入 stream_process，每次请求在调用处传入，避免单例竞态
                    model_name = (getattr(req, "model_name", None) or "").strip() or "qwen3.5-plus"
                    if model_name and hasattr(agent, "agents") and agent.agents:
                        try:
                            from gibh_agent.core.llm_client import LLMClientFactory
                            override_llm = LLMClientFactory.create_for_model(model_name)
                            for name, a in agent.agents.items():
                                if hasattr(a, "llm_client") and a.llm_client is not None:
                                    original_llm_clients[name] = a.llm_client
                                    a.llm_client = override_llm
                            logger.info(f"🔀 [ChatEndpoint] 已切换模型: {model_name}")
                        except Exception as swap_err:
                            logger.warning(f"⚠️ [ChatEndpoint] 模型切换失败，使用默认: {swap_err}")
                    async for event in orchestrator.stream_process(
                        query=req.message,
                        files=uploaded_files,
                        history=req.history or [],
                        session_id=req.session_id or "default",
                        test_dataset_id=req.test_dataset_id,
                        workflow_data=req.workflow_data,
                        user_id=req.user_id or "guest",
                        owner_id=owner_id,
                        db=db,
                        model_name=model_name,
                        target_domain=req.target_domain,
                        enabled_mcps=req.enabled_mcps or [],
                    ):
                        if isinstance(event, str) and "event: state_snapshot" in event:
                            for line in event.split("\n"):
                                if line.startswith("data:"):
                                    raw_line = line[5:].strip()
                                    try:
                                        state_snapshot_for_db = json.loads(raw_line)
                                    except Exception as parse_err:
                                        logger.error(
                                            "❌ [Chat] state_snapshot JSON 解析失败，将原始数据备份入库: %s",
                                            parse_err,
                                            exc_info=True,
                                        )
                                        raw_preview = raw_line[:1000] if len(raw_line) > 1000 else raw_line
                                        logger.error("❌ [Chat] 原始 data 前 1000 字符: %s", raw_preview)
                                        raw_backup = raw_line[:50000] if len(raw_line) > 50000 else raw_line
                                        state_snapshot_for_db = {
                                            "text": "⚠️ [系统警告] 工作流数据解析失败，原始数据备份:\n" + raw_backup,
                                            "workflow": None,
                                            "steps": [],
                                            "report": None,
                                        }
                                    break
                        yield event
                except Exception as e:
                    logger.error(f"❌ SSE 流式传输错误: {e}", exc_info=True)
                    error_event = f"event: error\ndata: {json.dumps({'error': str(e)}, ensure_ascii=False)}\n\n"
                    yield error_event
                finally:
                    # Phase 4: 流结束后写 agent 消息（有状态切片入库，废弃 content.events）
                    # 🔥 强制快照同步：前端传来的 workflow_data 含真实 enabled/selected，必须原样写入快照
                    if req.workflow_data and state_snapshot_for_db is not None:
                        state_snapshot_for_db["workflow"] = req.workflow_data
                    try:
                        content = {"state_snapshot": state_snapshot_for_db if state_snapshot_for_db is not None else {"text": "", "workflow": None, "steps": [], "report": None}}
                        msg = MessageModel(session_id=req.session_id, role="agent", content=content)
                        db.add(msg)
                        db.commit()
                        # 下发 message_saved 事件，前端可绑定到当前 AI 气泡的 dataset.messageId（commit 后 msg.id 已填充）
                        if getattr(msg, "id", None) is not None:
                            yield f"event: message_saved\ndata: {json.dumps({'message_id': msg.id}, ensure_ascii=False)}\n\n"
                    except (DataError, OperationalError) as db_err:
                        logger.error(
                            "❌ [Chat] Agent 消息入库失败（数据库错误，可能超出 max_allowed_packet 或字段限制）: %s",
                            db_err,
                            exc_info=True,
                        )
                        db.rollback()
                    except Exception as persist_err:
                        logger.error("❌ [Chat] Agent 消息入库失败: %s", persist_err, exc_info=True)
                        db.rollback()
                    # 恢复各 agent 原有 LLM 客户端
                    if original_llm_clients and hasattr(agent, "agents"):
                        for name, a in agent.agents.items():
                            if name in original_llm_clients:
                                a.llm_client = original_llm_clients[name]
                        logger.debug("🔀 [ChatEndpoint] 已恢复默认模型")
            
            return StreamingResponse(
                generate_sse(),
                media_type="text/event-stream",
                headers={
                    "Cache-Control": "no-cache",
                    "Connection": "keep-alive",
                    "X-Accel-Buffering": "no",  # 禁用 Nginx 缓冲
                    "Content-Type": "text/event-stream; charset=utf-8"
                }
            )
        except Exception as e:
            logger.error(f"❌ SSE 流式传输初始化失败: {e}", exc_info=True)
            return JSONResponse(
                status_code=500,
                content={
                    "type": "error",
                    "error": str(e),
                    "message": f"流式传输失败: {str(e)}"
                }
            )
    
    # 传统同步模式（保持向后兼容）
    try:
        logger.info(f"💬 收到聊天请求: {req.message}")
        logger.info(f"📁 上传文件数: {len(req.uploaded_files)}")
        logger.info(f"🔄 工作流数据: {req.workflow_data is not None}")
        
        # 🔧 修复：如果包含工作流数据，直接执行工作流（而不是通过 agent.process_query）
        if req.workflow_data:
            logger.info("🚀 检测到工作流执行请求，直接调用 execute_workflow")
            try:
                # 🔧 修复：优先使用 workflow_data 中的 file_paths（前端已经设置好）
                file_paths = req.workflow_data.get("file_paths", [])
                # 🔥 拆开逗号/分号拼接的路径（历史或旧前端可能传 "path1; path2" 单字符串），保证每条独立
                _flat = []
                for p in file_paths or []:
                    s = (p if isinstance(p, str) else str(p)).strip()
                    if not s:
                        continue
                    if ";" in s or "," in s:
                        for part in s.replace(";", ",").split(","):
                            t = part.strip()
                            if t:
                                _flat.append(t)
                    else:
                        _flat.append(s)
                file_paths = _flat
                logger.info(f"📁 从 workflow_data 获取的文件路径: {file_paths}")
                
                # 如果 workflow_data 中没有 file_paths，再从 uploaded_files 中提取
                if not file_paths:
                    logger.info("⚠️ workflow_data 中没有 file_paths，从 uploaded_files 中提取")
                    for file_info in req.uploaded_files:
                        file_name = file_info.get("file_name", "")
                        file_path_str = file_info.get("file_path", "")
                        
                        if file_path_str:
                            file_path = Path(file_path_str)
                        else:
                            file_path = UPLOAD_DIR / file_name if file_name else None
                        
                        if file_path and file_path.exists():
                            file_paths.append(str(file_path))
                
                logger.info(f"📂 最终文件路径列表: {file_paths}")
                
                # 验证文件路径是否存在
                valid_file_paths = []
                for fp in file_paths:
                    fp_path = Path(fp)
                    if fp_path.exists():
                        valid_file_paths.append(str(fp_path))
                    else:
                        logger.warning(f"⚠️ 文件不存在，跳过: {fp}")
                
                if not valid_file_paths:
                    raise ValueError("没有找到有效的输入文件。请确保文件已正确上传。")
                
                # 直接调用 execute_workflow 函数（不通过 HTTP）
                execute_request = {
                    "workflow_data": req.workflow_data,
                    "file_paths": valid_file_paths
                }
                # 调用 execute_workflow 函数（定义在下方）
                result = await execute_workflow(execute_request)
                return result
            except Exception as e:
                logger.error(f"❌ 工作流执行失败: {e}", exc_info=True)
                return JSONResponse(
                    status_code=500,
                    content={
                        "type": "error",
                        "error": str(e),
                        "message": f"工作流执行失败: {str(e)}"
                    }
                )
        
        # 🔥 Step 3: 尝试使用动态规划器（如果可用且查询看起来是工作流规划请求）
        if workflow_planner and not req.workflow_data:
            # 简单的启发式检测：如果查询包含分析相关的关键词，尝试使用规划器
            query_lower = req.message.lower()
            workflow_keywords = [
                "analyze", "analysis", "pca", "differential", "preprocess",
                "分析", "处理", "降维", "差异", "预处理"
            ]
            
            # 如果有上传文件或包含关键词，尝试使用规划器
            has_files = len(req.uploaded_files) > 0
            has_keywords = any(keyword in query_lower for keyword in workflow_keywords)
            
            if has_files or has_keywords:
                try:
                    logger.info("🧠 尝试使用动态规划器生成工作流...")
                    
                    # 提取文件路径（先转换 uploaded_files）
                    file_paths = []
                    for file_info in req.uploaded_files:
                        file_path = file_info.get("path") or file_info.get("file_name")
                        if file_path:
                            # 如果是相对路径，转换为绝对路径（保留完整相对结构，如 guest/session/10x_data_xxx）
                            if not Path(file_path).is_absolute():
                                file_path = str(UPLOAD_DIR / file_path)
                            file_paths.append(file_path)
                    
                    # 检测类别（简单启发式）
                    category_filter = None
                    if any(keyword in query_lower for keyword in ["metabolite", "代谢", "metabolomics"]):
                        category_filter = "Metabolomics"
                    elif any(keyword in query_lower for keyword in ["rna", "gene", "transcript", "转录"]):
                        category_filter = "scRNA-seq"
                    
                    # 调用规划器
                    plan_result = await workflow_planner.plan(
                        user_query=req.message,
                        context_files=file_paths,
                        category_filter=category_filter
                    )
                    
                    # 如果规划成功，返回结果
                    if plan_result.get("type") == "workflow_config":
                        logger.info("✅ 动态规划器成功生成工作流")
                        return JSONResponse(content=plan_result)
                    else:
                        logger.info(f"⚠️ 动态规划器返回: {plan_result.get('type')}，继续使用传统流程")
                        # 继续使用传统流程
                except Exception as planner_err:
                    logger.warning(f"⚠️ 动态规划器失败，回退到传统流程: {planner_err}")
                    # 继续使用传统流程
        
        # 🔥 转换文件路径：支持多种前端格式
        uploaded_files = []
        logger.info(f"📥 收到 uploaded_files: {len(req.uploaded_files)} 个文件")
        
        for file_info in req.uploaded_files:
            # 支持多种字段名：file_name/name, file_path/path
            file_name = file_info.get("file_name") or file_info.get("name", "")
            file_path_str = file_info.get("file_path") or file_info.get("path", "")
            
            # 🔒 安全：清理文件名
            if file_name:
                file_name = sanitize_filename(file_name)
            
            # 🔥 构建文件路径：优先使用 file_path，如果是相对路径则拼接 UPLOAD_DIR
            if file_path_str:
                file_path = Path(file_path_str)
                # 如果是相对路径，拼接 UPLOAD_DIR
                if not file_path.is_absolute():
                    file_path = UPLOAD_DIR / file_path
            elif file_name:
                # 如果没有路径，使用文件名在 UPLOAD_DIR 中查找
                file_path = UPLOAD_DIR / file_name
            else:
                logger.warning(f"⚠️ 无法确定文件路径，跳过: {file_info}")
                continue
            
            # 🔒 安全：验证路径在允许的目录内
            try:
                file_path = validate_file_path(file_path, UPLOAD_DIR)
            except HTTPException:
                logger.warning(f"⚠️ 不安全的文件路径，跳过: {file_path}")
                continue
            
            # 检查文件是否存在
            if not file_path.exists():
                logger.warning(f"⚠️ 文件不存在: {file_path}，尝试查找...")
                # 尝试在 UPLOAD_DIR 中查找同名文件
                if file_name:
                    alt_path = UPLOAD_DIR / file_name
                    if alt_path.exists():
                        file_path = alt_path
                        logger.info(f"✅ 找到文件: {file_path}")
                    else:
                        logger.warning(f"⚠️ 文件不存在，跳过: {file_path}")
                        continue
                else:
                    continue
            
            uploaded_files.append({
                "name": file_name or os.path.basename(str(file_path)),
                "path": str(file_path)
            })
        
        logger.info(f"📂 处理文件: {len(uploaded_files)} 个有效文件")
        logger.info(f"📂 文件路径列表: {[f['path'] for f in uploaded_files]}")
        
        # 处理查询
        # #region debug log
        try:
            debug_log_path = Path("/app/debug.log")
            debug_log_path.parent.mkdir(parents=True, exist_ok=True)
            with open(debug_log_path, 'a') as f:
                f.write(json.dumps({"location":"server.py:1161","message":"Before process_query","data":{"query":req.message,"uploaded_files_count":len(uploaded_files),"test_dataset_id":req.test_dataset_id},"timestamp":int(__import__('time').time()*1000),"sessionId":"debug-session","runId":"run1","hypothesisId":"A"})+"\n")
        except Exception as log_err:
            pass  # 日志写入失败不影响主流程
        # #endregion
        try:
            result = await agent.process_query(
                query=req.message,
                history=req.history or [],
                uploaded_files=uploaded_files,
                test_dataset_id=req.test_dataset_id,
                workflow_data=req.workflow_data,
                user_id=req.user_id,
                session_id=req.session_id
            )
        except Exception as process_err:
            # #region debug log
            try:
                debug_log_path = Path("/app/debug.log")
                debug_log_path.parent.mkdir(parents=True, exist_ok=True)
                with open(debug_log_path, 'a') as f:
                    f.write(json.dumps({"location":"server.py:1156","message":"process_query exception","data":{"error_type":type(process_err).__name__,"error_message":str(process_err),"traceback":traceback.format_exc()},"timestamp":int(__import__('time').time()*1000),"sessionId":"debug-session","runId":"run1","hypothesisId":"PROCESS_QUERY"})+"\n")
            except:
                pass
            # #endregion
            raise  # 重新抛出异常，让外层异常处理捕获
        
        # #region debug log
        try:
            debug_log_path = Path("/app/debug.log")
            debug_log_path.parent.mkdir(parents=True, exist_ok=True)
            with open(debug_log_path, 'a') as f:
                f.write(json.dumps({"location":"server.py:1168","message":"After process_query","data":{"result_type":type(result).__name__,"result_keys":list(result.keys()) if isinstance(result,dict) else None,"result_type_value":result.get('type') if isinstance(result,dict) else None},"timestamp":int(__import__('time').time()*1000),"sessionId":"debug-session","runId":"run1","hypothesisId":"A"})+"\n")
        except:
            pass
        # #endregion
        
        logger.info(f"✅ 处理完成，返回类型: {result.get('type', 'unknown')}")
        
        # 如果是工作流配置，返回 JSON
        if result.get("type") == "workflow_config":
            # 优先使用 result 中的 file_paths（可能来自测试数据集）
            # 如果没有，则使用 uploaded_files
            result_file_paths = result.get("file_paths", [])
            if not result_file_paths:
                result_file_paths = [f["path"] for f in uploaded_files]
            
            response_content = {
                "type": "workflow_config",
                "workflow_data": result.get("workflow_data"),
                "file_paths": result_file_paths
            }
            
            # 🔧 修复：如果包含诊断报告，也返回给前端
            if "diagnosis_report" in result:
                response_content["diagnosis_report"] = result["diagnosis_report"]
            
            # 🔧 修复：如果包含推荐信息，也返回给前端（代谢组学）
            if "recommendation" in result:
                response_content["recommendation"] = result["recommendation"]
            
            logger.info(f"📤 返回工作流配置: 包含推荐={('recommendation' in response_content)}, 包含诊断={('diagnosis_report' in response_content)}")
            
            return JSONResponse(content=response_content)
        
        # 如果是测试数据选择请求，格式化为用户友好的文本
        if result.get("type") == "test_data_selection":
            # #region debug log
            import json
            try:
                debug_log_path = Path("/app/debug.log")
                debug_log_path.parent.mkdir(parents=True, exist_ok=True)
                with open(debug_log_path, 'a') as f:
                    f.write(json.dumps({"location":"server.py:1178","message":"Entering test_data_selection handler","data":{"has_message":"message" in result,"has_options":"options" in result,"has_datasets_display":"datasets_display" in result,"has_datasets_json":"datasets_json" in result},"timestamp":int(__import__('time').time()*1000),"sessionId":"debug-session","runId":"run1","hypothesisId":"B"})+"\n")
            except:
                pass  # 日志写入失败不影响主流程
            # #endregion
            async def generate():
                try:
                    # #region debug log
                    try:
                        debug_log_path = Path("/app/debug.log")
                        debug_log_path.parent.mkdir(parents=True, exist_ok=True)
                        with open(debug_log_path, 'a') as f:
                            f.write(json.dumps({"location":"server.py:1181","message":"Inside generate()","data":{},"timestamp":int(__import__('time').time()*1000),"sessionId":"debug-session","runId":"run1","hypothesisId":"C"})+"\n")
                    except:
                        pass
                    # #endregion
                    # 构建用户友好的消息
                    message = result.get("message", "检测到您没有上传相关数据。请选择：")
                    options = result.get("options", [])
                    datasets_display = result.get("datasets_display", "")
                    # #region debug log
                    try:
                        debug_log_path = Path("/app/debug.log")
                        debug_log_path.parent.mkdir(parents=True, exist_ok=True)
                        with open(debug_log_path, 'a') as f:
                            f.write(json.dumps({"location":"server.py:1187","message":"Before datasets_json processing","data":{"message_type":type(message).__name__,"options_type":type(options).__name__,"datasets_display_type":type(datasets_display).__name__,"datasets_display_len":len(str(datasets_display)) if datasets_display else 0},"timestamp":int(__import__('time').time()*1000),"sessionId":"debug-session","runId":"run1","hypothesisId":"D"})+"\n")
                    except:
                        pass
                    # #endregion
                    
                    response_text = f"{message}\n\n"
                    for option in options:
                        response_text += f"  {option}\n"
                    
                    if datasets_display:
                        response_text += f"\n{datasets_display}\n"
                    
                    response_text += "\n💡 提示：回复数据集ID（如：pbmc_1k_v3）或上传您自己的数据文件。\n"
                    
                    # 同时保存数据集信息到响应中（用于前端处理）
                    # 这里我们通过特殊标记来传递 JSON 数据
                    # 将 JSON 中的换行符替换为空格，避免破坏 HTML 注释
                    datasets_json_raw = result.get('datasets_json', '[]')
                    # #region debug log
                    try:
                        debug_log_path = Path("/app/debug.log")
                        debug_log_path.parent.mkdir(parents=True, exist_ok=True)
                        with open(debug_log_path, 'a') as f:
                            f.write(json.dumps({"location":"server.py:1200","message":"Before datasets_json replace","data":{"datasets_json_type":type(datasets_json_raw).__name__,"datasets_json_is_none":datasets_json_raw is None,"datasets_json_len":len(str(datasets_json_raw)) if datasets_json_raw else 0},"timestamp":int(__import__('time').time()*1000),"sessionId":"debug-session","runId":"run1","hypothesisId":"B"})+"\n")
                    except:
                        pass
                    # #endregion
                    datasets_json = str(datasets_json_raw).replace('\n', ' ').replace('\r', '') if datasets_json_raw else '[]'
                    # #region debug log
                    try:
                        debug_log_path = Path("/app/debug.log")
                        debug_log_path.parent.mkdir(parents=True, exist_ok=True)
                        with open(debug_log_path, 'a') as f:
                            f.write(json.dumps({"location":"server.py:1203","message":"After datasets_json replace","data":{"datasets_json_len":len(datasets_json),"response_text_len":len(response_text)},"timestamp":int(__import__('time').time()*1000),"sessionId":"debug-session","runId":"run1","hypothesisId":"B"})+"\n")
                    except:
                        pass
                    # #endregion
                    response_text += f"\n<!-- DATASETS_JSON: {datasets_json} -->\n"
                    
                    # #region debug log
                    try:
                        debug_log_path = Path("/app/debug.log")
                        debug_log_path.parent.mkdir(parents=True, exist_ok=True)
                        with open(debug_log_path, 'a') as f:
                            f.write(json.dumps({"location":"server.py:1207","message":"Before yield","data":{"final_response_text_len":len(response_text)},"timestamp":int(__import__('time').time()*1000),"sessionId":"debug-session","runId":"run1","hypothesisId":"C"})+"\n")
                    except:
                        pass
                    # #endregion
                    yield response_text
                    # #region debug log
                    try:
                        debug_log_path = Path("/app/debug.log")
                        debug_log_path.parent.mkdir(parents=True, exist_ok=True)
                        with open(debug_log_path, 'a') as f:
                            f.write(json.dumps({"location":"server.py:1209","message":"After yield","data":{},"timestamp":int(__import__('time').time()*1000),"sessionId":"debug-session","runId":"run1","hypothesisId":"C"})+"\n")
                    except:
                        pass
                    # #endregion
                except Exception as e:
                    # #region debug log
                    import traceback
                    try:
                        debug_log_path = Path("/app/debug.log")
                        debug_log_path.parent.mkdir(parents=True, exist_ok=True)
                        with open(debug_log_path, 'a') as f:
                            f.write(json.dumps({"location":"server.py:1212","message":"Exception in generate()","data":{"error_type":type(e).__name__,"error_message":str(e),"traceback":traceback.format_exc()},"timestamp":int(__import__('time').time()*1000),"sessionId":"debug-session","runId":"run1","hypothesisId":"C"})+"\n")
                    except:
                        pass
                    # #endregion
                    logger.error(f"❌ 格式化测试数据选择响应错误: {e}", exc_info=True)
                    yield f"\n\n❌ 错误: {str(e)}"
            
            # #region debug log
            try:
                debug_log_path = Path("/app/debug.log")
                debug_log_path.parent.mkdir(parents=True, exist_ok=True)
                with open(debug_log_path, 'a') as f:
                    f.write(json.dumps({"location":"server.py:1218","message":"Before StreamingResponse","data":{},"timestamp":int(__import__('time').time()*1000),"sessionId":"debug-session","runId":"run1","hypothesisId":"E"})+"\n")
            except:
                pass
            # #endregion
            return StreamingResponse(generate(), media_type="text/plain")
        
        # 如果是聊天响应，返回流式
        if result.get("type") == "chat":
            async def generate():
                try:
                    response_iter = result.get("response")
                    if response_iter:
                        # 确保 response_iter 是异步迭代器
                        async for chunk in response_iter:
                            if chunk:
                                yield chunk
                    else:
                        logger.warning("⚠️ 聊天响应中没有 response 迭代器")
                        yield "❌ 错误: 无法获取响应"
                except Exception as e:
                    logger.error(f"❌ 流式响应错误: {e}", exc_info=True)
                    import traceback
                    logger.error(f"详细错误: {traceback.format_exc()}")
                    yield f"\n\n❌ 错误: {str(e)}"
            
            return StreamingResponse(
                generate(), 
                media_type="text/plain",
                headers={
                    "Cache-Control": "no-cache",
                    "Connection": "keep-alive",
                    "X-Accel-Buffering": "no"
                }
            )
        
        # 其他情况返回 JSON
        return JSONResponse(content=result)
        
    except Exception as e:
        # #region debug log
        import traceback
        try:
            debug_log_path = Path("/app/debug.log")
            debug_log_path.parent.mkdir(parents=True, exist_ok=True)
            with open(debug_log_path, 'a') as f:
                f.write(json.dumps({"location":"server.py:1210","message":"Exception in chat_endpoint","data":{"error_type":type(e).__name__,"error_message":str(e),"traceback":traceback.format_exc()},"timestamp":int(__import__('time').time()*1000),"sessionId":"debug-session","runId":"run1","hypothesisId":"ALL"})+"\n")
        except:
            pass  # 日志写入失败不影响主流程
        # #endregion
        error_detail = f"{type(e).__name__}: {str(e)}"
        logger.error(f"❌ 处理失败: {error_detail}", exc_info=True)
        logger.error(f"详细错误: {traceback.format_exc()}")
        raise HTTPException(status_code=500, detail=error_detail)


@app.post("/api/execute")
async def execute_workflow(request: dict):
    """执行工作流接口"""
    if not agent:
        raise HTTPException(status_code=500, detail="智能体未初始化")
    
    try:
        workflow_data = request.get("workflow_data")
        file_paths = request.get("file_paths", [])
        enabled_mcps = request.get("enabled_mcps") or []
        if isinstance(enabled_mcps, str):
            try:
                enabled_mcps = json.loads(enabled_mcps)
            except Exception:
                enabled_mcps = []
        if not isinstance(enabled_mcps, list):
            enabled_mcps = []
        enabled_mcps = [str(x).strip() for x in enabled_mcps if x]
        
        logger.info(f"🚀 开始执行工作流: {len(file_paths)} 个文件")
        
        # 🔧 修复：优先检查 workflow_name 中是否包含代谢组关键词
        workflow_name = workflow_data.get("workflow_name", "")
        routing = None
        target_agent = None
        route_query = None
        
        # 方法1: 如果有 workflow_name，优先检查是否包含代谢组关键词
        if workflow_name:
            workflow_name_lower = workflow_name.lower()
            # 如果 workflow_name 包含代谢组关键词，直接路由到 metabolomics_agent
            if any(kw in workflow_name_lower for kw in ["metabolomics", "代谢组", "代谢"]):
                logger.info(f"✅ 根据 workflow_name 直接路由到 metabolomics_agent: {workflow_name}")
                routing = "metabolomics_agent"
                target_agent = agent.agents.get(routing)
                if not target_agent:
                    logger.warning(f"⚠️ metabolomics_agent 不存在，使用默认 rna_agent")
                    target_agent = agent.agents.get("rna_agent")
                    routing = "rna_agent"
            else:
                route_query = workflow_name
        # 方法2: 根据文件类型构建查询
        elif file_paths:
            file_path = file_paths[0]
            file_ext = os.path.splitext(file_path)[1].lower()
            if file_ext == ".csv":
                route_query = "metabolomics analysis"
            elif file_ext in [".h5ad", ".h5"]:
                route_query = "single cell transcriptomics analysis"
            elif "fastq" in file_path.lower():
                route_query = "single cell RNA-seq analysis"
            else:
                route_query = "bioinformatics analysis"
        else:
            route_query = "bioinformatics analysis"
        
        # 准备上传文件列表（用于 RouterAgent）
        uploaded_files_for_router = []
        for file_path in file_paths:
            uploaded_files_for_router.append({
                "name": os.path.basename(file_path),
                "path": file_path
            })
        
        # 🔧 修复：如果还没有路由，使用 RouterAgent 进行路由决策
        if not routing or not target_agent:
            try:
                route_result = await agent.router.process_query(
                    query=route_query,
                    history=[],
                    uploaded_files=uploaded_files_for_router
                )
                
                routing = route_result.get("routing", "rna_agent")
                modality = route_result.get("modality", "")
                # HITL：模态未知时返回「请选择数据类型」，由前端展示按钮
                if routing == "ask_modality" or modality == "unknown":
                    logger.info("📋 RouterAgent 返回 ask_modality，要求用户选择数据类型")
                    return JSONResponse(content={
                        "type": "ask_modality",
                        "message": "检测到数据类型不明确，请选择分析类型：",
                        "action": "show_modality_selector",
                        "modalities": route_result.get("available_modalities", [
                            {"id": "RNA", "label": "单细胞/转录组 (RNA)"},
                            {"id": "Metabolomics", "label": "代谢组学 (Metabolomics)"},
                            {"id": "Spatial", "label": "空间转录组 (Spatial)"},
                            {"id": "Radiomics", "label": "影像组学 (Radiomics)"},
                        ]),
                        "file_paths": file_paths,
                    })
                target_agent = agent.agents.get(routing)
                
                # 如果路由的智能体不存在，使用默认的 RNA Agent
                if not target_agent:
                    logger.warning(f"⚠️ 路由的智能体不存在: {routing}，使用默认 rna_agent")
                    target_agent = agent.agents.get("rna_agent")
                    routing = "rna_agent"
                
                if not target_agent:
                    raise HTTPException(status_code=500, detail="RNA Agent 未找到")
                
                logger.info(f"✅ RouterAgent 路由结果: {routing} (confidence: {route_result.get('confidence', 0):.2f}, modality: {route_result.get('modality', 'unknown')})")
                
            except Exception as e:
                logger.error(f"❌ RouterAgent 路由失败: {e}，使用默认 rna_agent", exc_info=True)
                # 降级到默认 Agent
                target_agent = agent.agents.get("rna_agent")
                routing = "rna_agent"
                if not target_agent:
                    raise HTTPException(status_code=500, detail="RNA Agent 未找到")
        
        # 🔥 Step 4: 使用通用执行器（动态执行，不依赖硬编码逻辑）
        try:
            from gibh_agent.core.executor import WorkflowExecutor
            
            logger.info("🔧 使用通用执行器执行工作流...")
            
            # 设置输出目录
            output_dir = str(RESULTS_DIR / f"run_{datetime.now().strftime('%Y%m%d_%H%M%S')}")
            
            # 创建执行器并执行（传递 agent 实例以生成诊断）
            executor = WorkflowExecutor(output_dir=output_dir)
            report_data = executor.execute_workflow(
                workflow_data=workflow_data,
                file_paths=file_paths,
                output_dir=output_dir,
                agent=target_agent,  # 🔥 传递 agent 实例以生成 AI Expert Diagnosis
                enabled_mcps=enabled_mcps,
            )
            
            logger.info("✅ 通用执行器执行完成")
            
            # 🔥 生成 AI Expert Diagnosis（如果提供了 Agent 实例）
            if target_agent and hasattr(target_agent, '_generate_analysis_summary'):
                try:
                    logger.info("📝 [Server] 生成 AI Expert Diagnosis...")
                    
                    # 检测组学类型
                    omics_type = "Metabolomics"  # 默认
                    steps = workflow_data.get("steps", [])
                    if any("rna" in step.get("id", "").lower() or "rna" in step.get("tool_id", "").lower() for step in steps):
                        omics_type = "scRNA"
                    elif any("metabolomics" in step.get("id", "").lower() or "metabolomics" in step.get("tool_id", "").lower() for step in steps):
                        omics_type = "Metabolomics"
                    
                    # 调用异步方法生成诊断
                    steps_results = report_data.get("steps_results", [])
                    workflow_name = report_data.get("workflow_name", "Analysis Pipeline")
                    diagnosis = await target_agent._generate_analysis_summary(
                        steps_results, 
                        omics_type, 
                        workflow_name
                    )
                    
                    if diagnosis:
                        logger.info(f"✅ [Server] AI Expert Diagnosis 生成成功，长度: {len(diagnosis)}")
                        report_data["diagnosis"] = diagnosis
                    else:
                        logger.warning("⚠️ [Server] AI Expert Diagnosis 生成失败或返回空")
                except Exception as diag_err:
                    logger.error(f"❌ [Server] 生成 AI Expert Diagnosis 失败: {diag_err}", exc_info=True)
                    # 不中断工作流，继续执行
            
            # 构建返回结果（符合前端格式）
            return JSONResponse(content={
                "type": "analysis_report",
                "status": "success",
                "report_data": report_data,
                "reply": "✅ 工作流执行完成（使用动态执行引擎）",
                "thought": "[THOUGHT] 使用 ToolRegistry 动态执行，工具无关"
            })
        
        except ImportError:
            logger.warning("⚠️ 通用执行器未找到，回退到传统执行方式")
            # 回退到传统执行方式
            output_dir = str(RESULTS_DIR / f"run_{datetime.now().strftime('%Y%m%d_%H%M%S')}")
            os.makedirs(output_dir, exist_ok=True)
            
            report = await target_agent.execute_workflow(
                workflow_config=workflow_data,
                file_paths=file_paths,
                output_dir=output_dir
            )
            
            logger.info(f"✅ 工作流执行完成: {report.get('status')}")
            
            # 处理图片路径（传统方式）
            if report.get("final_plot"):
                plot_path = report["final_plot"]
                if not plot_path.startswith("/results/"):
                    if plot_path.startswith("results/"):
                        plot_path = "/" + plot_path
                    elif "/" in plot_path:
                        plot_path = f"/results/{plot_path}"
                    else:
                        run_name = os.path.basename(output_dir)
                        plot_path = f"/results/{run_name}/{plot_path}"
                report["final_plot"] = plot_path
            
            # 处理步骤中的图片路径
            run_name = os.path.basename(output_dir)
            if report.get("steps_details"):
                for step in report["steps_details"]:
                    if step.get("plot"):
                        plot_path = step["plot"]
                        if not plot_path.startswith("/results/"):
                            if plot_path.startswith("results/"):
                                plot_path = "/" + plot_path
                            elif "/" in plot_path:
                                plot_path = f"/results/{plot_path}"
                            else:
                                plot_path = f"/results/{run_name}/{plot_path}"
                        step["plot"] = plot_path
            
            # 返回传统格式的结果
            return JSONResponse(content={
                "type": "analysis_report",
                "status": report.get("status", "success"),
                "report_data": report
            })
        
        # 处理通用执行器返回的图片路径
        logger.info(f"✅ 工作流执行完成: {report_data.get('status')}")
        
        # 处理图片路径，转换为可访问的 URL（在返回之前）
        # 图片保存在 results/run_xxx/ 目录，需要转换为 /results/run_xxx/filename
        if report.get("final_plot"):
            plot_path = report["final_plot"]
            # 确保路径以 /results/ 开头
            if not plot_path.startswith("/results/"):
                if plot_path.startswith("results/"):
                    plot_path = "/" + plot_path
                elif "/" in plot_path:
                    # 如果包含 run_xxx/filename 格式，添加 results 前缀
                    plot_path = f"/results/{plot_path}"
                else:
                    # 如果只是文件名，需要找到对应的 run 目录
                    # 从 output_dir 中提取 run_xxx
                    run_name = os.path.basename(output_dir)
                    plot_path = f"/results/{run_name}/{plot_path}"
            report_data["final_plot"] = plot_path
        
        # 处理步骤中的图片路径（steps_details）
        run_name = os.path.basename(output_dir)
        if report_data.get("steps_details"):
            for step in report_data["steps_details"]:
                if step.get("plot"):
                    plot_path = step["plot"]
                    # 确保路径以 /results/ 开头
                    if not plot_path.startswith("/results/"):
                        if plot_path.startswith("results/"):
                            plot_path = "/" + plot_path
                        elif "/" in plot_path:
                            plot_path = f"/results/{plot_path}"
                        else:
                            plot_path = f"/results/{run_name}/{plot_path}"
                    step["plot"] = plot_path
                
                # 处理 step_result 中的图片路径（含 STED-EC report_data.images 展平后的绝对路径）
                if step.get("step_result") and step["step_result"].get("data", {}).get("images"):
                    images = step["step_result"]["data"]["images"]
                    fixed_images = []
                    _results_prefix = str(RESULTS_DIR).rstrip("/")
                    for img_path in images:
                        if not isinstance(img_path, str):
                            continue
                        if img_path.startswith("/results/"):
                            fixed_images.append(img_path)
                        elif img_path.startswith(_results_prefix + "/") or img_path.startswith(_results_prefix):
                            fixed_images.append("/results/" + img_path[len(_results_prefix):].lstrip("/"))
                        elif img_path.startswith("results/"):
                            fixed_images.append("/" + img_path)
                        elif "/" in img_path:
                            fixed_images.append(f"/results/{img_path}")
                        else:
                            fixed_images.append(f"/results/{run_name}/{img_path}")
                    step["step_result"]["data"]["images"] = fixed_images
        
        # 确保 steps_results 存在（前端可直接使用）
        if "steps_results" not in report and "steps_details" in report:
            steps_results = []
            for step_detail in report.get("steps_details", []):
                if "step_result" in step_detail:
                    step_result = step_detail["step_result"].copy()
                    # 确保图片路径正确（含 STED-EC 绝对路径规范化）
                    if step_result.get("data", {}).get("images"):
                        images = step_result["data"]["images"]
                        fixed_images = []
                        _results_prefix = str(RESULTS_DIR).rstrip("/")
                        for img_path in images:
                            if not isinstance(img_path, str):
                                continue
                            if img_path.startswith("/results/"):
                                fixed_images.append(img_path)
                            elif img_path.startswith(_results_prefix + "/") or img_path.startswith(_results_prefix):
                                fixed_images.append("/results/" + img_path[len(_results_prefix):].lstrip("/"))
                            elif img_path.startswith("results/"):
                                fixed_images.append("/" + img_path)
                            elif "/" in img_path:
                                fixed_images.append(f"/results/{img_path}")
                            else:
                                fixed_images.append(f"/results/{run_name}/{img_path}")
                        step_result["data"]["images"] = fixed_images
                    steps_results.append(step_result)
                else:
                    # 兼容旧格式
                    step_result = {
                        "step_name": step_detail.get("name", "Unknown"),
                        "status": step_detail.get("status", "success"),
                        "logs": step_detail.get("summary", ""),
                        "data": {}
                    }
                    # 如果有 plot，添加到 data.images
                    if step_detail.get("plot"):
                        plot_path = step_detail["plot"]
                        if not plot_path.startswith("/results/"):
                            if plot_path.startswith("results/"):
                                plot_path = "/" + plot_path
                            elif "/" in plot_path:
                                plot_path = f"/results/{plot_path}"
                            else:
                                plot_path = f"/results/{run_name}/{plot_path}"
                        step_result["data"]["images"] = [plot_path]
                    steps_results.append(step_result)
            report["steps_results"] = steps_results
        
        # 处理 steps_results 中的图片路径（如果存在，含 STED-EC 绝对路径规范化）
        if report.get("steps_results"):
            _results_prefix = str(RESULTS_DIR).rstrip("/")
            for step_result in report["steps_results"]:
                if step_result.get("data", {}).get("images"):
                    images = step_result["data"]["images"]
                    fixed_images = []
                    for img_path in images:
                        if not isinstance(img_path, str):
                            continue
                        if img_path.startswith("/results/"):
                            fixed_images.append(img_path)
                        elif img_path.startswith(_results_prefix + "/") or img_path.startswith(_results_prefix):
                            fixed_images.append("/results/" + img_path[len(_results_prefix):].lstrip("/"))
                        elif img_path.startswith("results/"):
                            fixed_images.append("/" + img_path)
                        elif "/" in img_path:
                            fixed_images.append(f"/results/{img_path}")
                        else:
                            fixed_images.append(f"/results/{run_name}/{img_path}")
                    step_result["data"]["images"] = fixed_images
        
        # 🔧 修复：返回正确的工作流执行结果格式
        return JSONResponse(content={
            "type": "analysis_report",
            "status": "success",
            "report_data": report
        })
        
    except Exception as e:
        import traceback
        error_detail = f"{type(e).__name__}: {str(e)}"
        error_traceback = traceback.format_exc()
        logger.error(f"❌ 工作流执行失败: {error_detail}", exc_info=True)
        logger.error(f"详细错误: {error_traceback}")
        # 返回更详细的错误信息
        return JSONResponse(
            status_code=500,
            content={
                "status": "error",
                "error": error_detail,
                "error_detail": error_traceback,
                "message": f"工作流执行失败: {error_detail}"
            }
        )


@app.get("/api/logs/stream")
async def stream_logs():
    """实时日志流接口（Server-Sent Events）"""
    logger.info("📡 新的日志流连接")
    
    async def event_generator():
        q = asyncio.Queue(maxsize=100)
        log_listeners.add(q)
        
        try:
            # 先发送历史日志
            history_logs = list(log_buffer)[-100:]  # 最近100条
            logger.info(f"📤 发送历史日志: {len(history_logs)} 条")
            for entry in history_logs:
                yield f"data: {json.dumps(entry, ensure_ascii=False)}\\n\\n"
            
            # 实时发送新日志
            while True:
                try:
                    entry = await asyncio.wait_for(q.get(), timeout=1.0)
                    yield f"data: {json.dumps(entry, ensure_ascii=False)}\\n\\n"
                except asyncio.TimeoutError:
                    # 发送心跳保持连接
                    yield f"data: {json.dumps({'type': 'heartbeat', 'timestamp': datetime.now().isoformat()})}\\n\\n"
        except asyncio.CancelledError:
            logger.info("📡 日志流连接已取消")
        except Exception as e:
            logger.error(f"❌ 日志流错误: {e}", exc_info=True)
        finally:
            log_listeners.discard(q)
            logger.info("📡 日志流连接已关闭")
    
    return StreamingResponse(
        event_generator(),
        media_type="text/event-stream",
        headers={
            "Cache-Control": "no-cache",
            "Connection": "keep-alive",
            "X-Accel-Buffering": "no"
        }
    )


@app.get("/api/logs")
async def get_logs(limit: int = 100):
    """获取历史日志"""
    return JSONResponse(content={
        "logs": list(log_buffer)[-limit:],
        "total": len(log_buffer)
    })


# 🔥 Step 2: Tool-RAG API - 工具检索端点
@app.get("/api/tools/search")
async def search_tools(
    query: str,
    top_k: int = 5,
    category: Optional[str] = None
):
    """
    语义搜索工具
    
    Args:
        query: 查询文本（自然语言）
        top_k: 返回前 k 个最相关的工具（默认 5）
        category: 可选的类别过滤器
    
    Returns:
        相关工具的 JSON Schema 列表
    """
    if tool_retriever is None:
        raise HTTPException(
            status_code=503,
            detail="工具检索器未初始化。请检查 Ollama 服务和依赖是否已安装。"
        )
    
    try:
        tools = tool_retriever.retrieve(query=query, top_k=top_k, category_filter=category)
        return {
            "status": "success",
            "query": query,
            "count": len(tools),
            "tools": tools
        }
    except Exception as e:
        logger.error(f"❌ 工具搜索失败: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=f"工具搜索失败: {str(e)}")


@app.get("/api/tools/list")
async def list_tools():
    """
    列出所有已注册的工具
    
    Returns:
        工具名称列表
    """
    if tool_retriever is None:
        raise HTTPException(
            status_code=503,
            detail="工具检索器未初始化"
        )
    
    try:
        tools = tool_retriever.list_all_tools()
        return {
            "status": "success",
            "count": len(tools),
            "tools": tools
        }
    except Exception as e:
        logger.error(f"❌ 列出工具失败: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=f"列出工具失败: {str(e)}")


@app.get("/api/tools/{tool_name}")
async def get_tool_schema(tool_name: str):
    """
    获取特定工具的完整 Schema
    
    Args:
        tool_name: 工具名称
    
    Returns:
        工具的完整 JSON Schema
    """
    if tool_retriever is None:
        raise HTTPException(
            status_code=503,
            detail="工具检索器未初始化"
        )
    
    try:
        tool_schema = tool_retriever.get_tool_by_name(tool_name)
        if tool_schema is None:
            raise HTTPException(status_code=404, detail=f"工具 '{tool_name}' 不存在")
        
        return {
            "status": "success",
            "tool": tool_schema
        }
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"❌ 获取工具 Schema 失败: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=f"获取工具 Schema 失败: {str(e)}")


# ==================== 🔥 ARCHITECTURAL UPGRADE: 新的工作流 API 端点 ====================

class WorkflowPlanRequest(BaseModel):
    """工作流规划请求"""
    query: str
    file_metadata: Optional[Dict[str, Any]] = None
    user_id: Optional[str] = None


class WorkflowSaveRequest(BaseModel):
    """保存工作流请求"""
    name: str
    workflow_json: Dict[str, Any]
    user_id: Optional[str] = None


@app.post("/api/workflows/plan")
async def plan_workflow(req: WorkflowPlanRequest):
    """
    规划工作流（plan-first：可以在没有文件的情况下生成工作流）
    
    🔥 ARCHITECTURAL UPGRADE:
    - 支持 plan-first：可以在没有文件的情况下生成工作流
    - 使用 WorkflowRegistry 进行严格的域绑定
    - 使用 DAG 依赖解析（代码逻辑，非 LLM 幻觉）
    """
    try:
        user_id = req.user_id or "guest"
        logger.info(f"📋 [WorkflowPlan] 用户查询: '{req.query}' (User: {user_id})")
        
        # 初始化 SOPPlanner（如果还没有）
        if not workflow_planner:
            from gibh_agent.core.llm_client import LLMClient
            from gibh_agent.core.planner import SOPPlanner
            # 🔥 TASK 2: 使用统一的LLM客户端创建方法
            from gibh_agent.core.llm_client import LLMClientFactory
            llm_client = LLMClientFactory.create_default() if agent else None
            if not llm_client:
                raise HTTPException(status_code=500, detail="LLM 客户端未初始化")
            
            planner = SOPPlanner(tool_retriever, llm_client)
        else:
            planner = workflow_planner
        
        # 生成工作流计划（SOPPlanner.generate_plan 为异步生成器，取最终 workflow 事件）
        workflow_config = None
        async for _ev, _data in planner.generate_plan(
            user_query=req.query,
            file_metadata=req.file_metadata,
        ):
            if _ev == "workflow":
                workflow_config = _data
        if workflow_config is None:
            raise HTTPException(status_code=500, detail="工作流规划未返回结果")

        return {
            "status": "success",
            "workflow": workflow_config,
            "user_id": user_id
        }
    
    except Exception as e:
        logger.error(f"❌ 工作流规划失败: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=f"工作流规划失败: {str(e)}")


@app.post("/api/workflows/save")
async def save_workflow(req: WorkflowSaveRequest):
    """
    保存工作流（书签）
    
    🔥 ARCHITECTURAL UPGRADE: 多用户支持
    """
    try:
        user_id = req.user_id or "guest"
        logger.info(f"💾 [WorkflowSave] 保存工作流: '{req.name}' (User: {user_id})")
        
        from gibh_agent.db import get_db
        db = get_db()
        workflow_id = db.save_workflow(
            user_id=user_id,
            name=req.name,
            workflow_json=req.workflow_json
        )
        
        return {
            "status": "success",
            "workflow_id": workflow_id,
            "message": f"工作流 '{req.name}' 已保存"
        }
    
    except Exception as e:
        logger.error(f"❌ 保存工作流失败: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=f"保存工作流失败: {str(e)}")


@app.get("/api/workflows/list")
async def list_workflows(user_id: Optional[str] = None):
    """
    列出用户的所有工作流（书签）
    
    🔥 ARCHITECTURAL UPGRADE: 多用户支持
    """
    try:
        user_id = user_id or "guest"
        logger.info(f"📋 [WorkflowList] 列出工作流 (User: {user_id})")
        
        from gibh_agent.db import get_db
        db = get_db()
        workflows = db.list_workflows(user_id=user_id)
        
        return {
            "status": "success",
            "workflows": workflows,
            "count": len(workflows)
        }
    
    except Exception as e:
        logger.error(f"❌ 列出工作流失败: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=f"列出工作流失败: {str(e)}")


@app.delete("/api/workflows/{workflow_id}")
async def delete_workflow(workflow_id: int, user_id: Optional[str] = None):
    """
    删除工作流
    
    🔥 ARCHITECTURAL UPGRADE: 多用户支持
    """
    try:
        user_id = user_id or "guest"
        logger.info(f"🗑️ [WorkflowDelete] 删除工作流: {workflow_id} (User: {user_id})")
        
        from gibh_agent.db import get_db
        db = get_db()
        deleted = db.delete_workflow(workflow_id=workflow_id, user_id=user_id)
        
        if deleted:
            return {
                "status": "success",
                "message": f"工作流 {workflow_id} 已删除"
            }
        else:
            raise HTTPException(status_code=404, detail="工作流不存在或无权删除")
    
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"❌ 删除工作流失败: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=f"删除工作流失败: {str(e)}")


@app.get("/api/jobs/history")
async def get_job_history(
    user_id: Optional[str] = None,
    status: Optional[str] = None,
    limit: int = 50
):
    """
    获取任务执行历史
    
    🔥 ARCHITECTURAL UPGRADE: 多用户支持
    """
    try:
        user_id = user_id or "guest"
        logger.info(f"📜 [JobHistory] 获取任务历史 (User: {user_id}, Status: {status or 'all'})")
        
        from gibh_agent.db import get_db
        db = get_db()
        jobs = db.list_jobs(user_id=user_id, status=status, limit=limit)
        
        return {
            "status": "success",
            "jobs": jobs,
            "count": len(jobs)
        }
    
    except Exception as e:
        logger.error(f"❌ 获取任务历史失败: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=f"获取任务历史失败: {str(e)}")


@app.get("/api/workflow/status/{run_id}")
async def get_workflow_status(run_id: str):
    """
    查询工作流状态（兼容旧架构）
    如果使用 Celery，查询 Celery 任务状态
    如果使用同步执行，返回 not_found（因为同步执行没有任务ID）
    """
    try:
        # 尝试从 Celery 查询任务状态
        from celery.result import AsyncResult
        from gibh_agent.core.celery_app import celery_app
        
        task_result = AsyncResult(run_id, app=celery_app)
        
        response = {
            "status": "running",
            "completed": False,
            "steps_status": [],
            "error": None
        }
        
        if task_result.state == 'PENDING':
            response["status"] = "running"
        elif task_result.state == 'SUCCESS':
            response["status"] = "success"
            response["completed"] = True
            result_data = task_result.result
            if result_data:
                response["report_data"] = result_data
                if "steps_details" in result_data:
                    response["steps_status"] = result_data["steps_details"]
                elif "steps" in result_data:
                    response["steps_status"] = result_data["steps"]
        elif task_result.state == 'FAILURE':
            response["status"] = "failed"
            response["completed"] = True
            response["error"] = str(task_result.result)
        elif task_result.state == 'PROGRESS':
            info = task_result.info
            if isinstance(info, dict):
                response["steps_status"] = info.get("steps", [])
        
        return JSONResponse(content=response)
        
    except ImportError:
        # Celery 未安装或未配置，返回 not_found
        return JSONResponse(
            status_code=404,
            content={
                "status": "not_found",
                "message": "工作流状态查询需要 Celery 支持，当前使用同步执行模式"
            }
        )
    except Exception as e:
        logger.error(f"❌ 查询工作流状态失败: {e}", exc_info=True)
        return JSONResponse(
            status_code=500,
            content={
                "status": "error",
                "error": str(e)
            }
        )


if __name__ == "__main__":
    import uvicorn
    import json
    
    port = int(os.getenv("PORT", 8018))
    logger.info(f"🚀 启动服务器，端口: {port}")
    logger.info(f"📁 上传目录: {UPLOAD_DIR.absolute()}")
    logger.info(f"📁 结果目录: {RESULTS_DIR.absolute()}")
    
    uvicorn.run(
        "server:app",
        host="0.0.0.0",
        port=port,
        log_level="info",
        reload=True
    )

