"""
工作流数据库模块

使用 SQLite 存储：
1. saved_workflows: 保存的工作流（书签）
2. job_history: 任务执行历史
"""
import sqlite3
import json
import logging
from pathlib import Path
from typing import Dict, Any, List, Optional
from datetime import datetime
from contextlib import contextmanager

logger = logging.getLogger(__name__)


class WorkflowDB:
    """
    工作流数据库管理器
    
    提供工作流和任务历史的 CRUD 操作。
    """
    
    def __init__(self, db_path: str = "workflows.db"):
        """
        初始化数据库
        
        Args:
            db_path: 数据库文件路径
        """
        self.db_path = Path(db_path)
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        self._init_tables()
    
    def _init_tables(self):
        """初始化数据库表"""
        with self._get_connection() as conn:
            cursor = conn.cursor()
            
            # 表1: saved_workflows（保存的工作流）
            cursor.execute("""
                CREATE TABLE IF NOT EXISTS saved_workflows (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    user_id TEXT NOT NULL DEFAULT 'guest',
                    name TEXT NOT NULL,
                    workflow_json TEXT NOT NULL,
                    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
                )
            """)
            
            # 表2: job_history（任务执行历史）
            cursor.execute("""
                CREATE TABLE IF NOT EXISTS job_history (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    user_id TEXT NOT NULL DEFAULT 'guest',
                    workflow_snapshot TEXT NOT NULL,
                    status TEXT NOT NULL DEFAULT 'pending',
                    result_summary TEXT,
                    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                    completed_at TIMESTAMP
                )
            """)
            
            # 创建索引
            cursor.execute("CREATE INDEX IF NOT EXISTS idx_saved_workflows_user_id ON saved_workflows(user_id)")
            cursor.execute("CREATE INDEX IF NOT EXISTS idx_job_history_user_id ON job_history(user_id)")
            cursor.execute("CREATE INDEX IF NOT EXISTS idx_job_history_status ON job_history(status)")
            
            conn.commit()
            logger.info("✅ [WorkflowDB] 数据库表初始化完成")
    
    @contextmanager
    def _get_connection(self):
        """获取数据库连接（上下文管理器）"""
        conn = sqlite3.connect(str(self.db_path))
        conn.row_factory = sqlite3.Row  # 返回字典格式的行
        try:
            yield conn
        finally:
            conn.close()
    
    # ==================== saved_workflows 操作 ====================
    
    def save_workflow(
        self,
        user_id: str,
        name: str,
        workflow_json: Dict[str, Any]
    ) -> int:
        """
        保存工作流（书签）
        
        Args:
            user_id: 用户ID（默认 "guest"）
            name: 工作流名称
            workflow_json: 工作流配置（字典）
            
        Returns:
            保存的工作流ID
        """
        with self._get_connection() as conn:
            cursor = conn.cursor()
            cursor.execute("""
                INSERT INTO saved_workflows (user_id, name, workflow_json, updated_at)
                VALUES (?, ?, ?, ?)
            """, (
                user_id,
                name,
                json.dumps(workflow_json, ensure_ascii=False),
                datetime.now().isoformat()
            ))
            conn.commit()
            workflow_id = cursor.lastrowid
            logger.info(f"✅ [WorkflowDB] 保存工作流: {name} (ID: {workflow_id}, User: {user_id})")
            return workflow_id
    
    def get_workflow(self, workflow_id: int, user_id: Optional[str] = None) -> Optional[Dict[str, Any]]:
        """
        获取工作流
        
        Args:
            workflow_id: 工作流ID
            user_id: 用户ID（可选，用于权限检查）
            
        Returns:
            工作流字典，如果不存在返回 None
        """
        with self._get_connection() as conn:
            cursor = conn.cursor()
            if user_id:
                cursor.execute("""
                    SELECT * FROM saved_workflows
                    WHERE id = ? AND user_id = ?
                """, (workflow_id, user_id))
            else:
                cursor.execute("""
                    SELECT * FROM saved_workflows
                    WHERE id = ?
                """, (workflow_id,))
            
            row = cursor.fetchone()
            if row:
                return {
                    "id": row["id"],
                    "user_id": row["user_id"],
                    "name": row["name"],
                    "workflow_json": json.loads(row["workflow_json"]),
                    "created_at": row["created_at"],
                    "updated_at": row["updated_at"]
                }
            return None
    
    def list_workflows(self, user_id: str = "guest") -> List[Dict[str, Any]]:
        """
        列出用户的所有工作流
        
        Args:
            user_id: 用户ID
            
        Returns:
            工作流列表
        """
        with self._get_connection() as conn:
            cursor = conn.cursor()
            cursor.execute("""
                SELECT * FROM saved_workflows
                WHERE user_id = ?
                ORDER BY updated_at DESC
            """, (user_id,))
            
            workflows = []
            for row in cursor.fetchall():
                workflows.append({
                    "id": row["id"],
                    "user_id": row["user_id"],
                    "name": row["name"],
                    "workflow_json": json.loads(row["workflow_json"]),
                    "created_at": row["created_at"],
                    "updated_at": row["updated_at"]
                })
            
            logger.info(f"✅ [WorkflowDB] 列出工作流: {len(workflows)} 个 (User: {user_id})")
            return workflows
    
    def delete_workflow(self, workflow_id: int, user_id: str = "guest") -> bool:
        """
        删除工作流
        
        Args:
            workflow_id: 工作流ID
            user_id: 用户ID（用于权限检查）
            
        Returns:
            True 如果删除成功，False 否则
        """
        with self._get_connection() as conn:
            cursor = conn.cursor()
            cursor.execute("""
                DELETE FROM saved_workflows
                WHERE id = ? AND user_id = ?
            """, (workflow_id, user_id))
            conn.commit()
            
            deleted = cursor.rowcount > 0
            if deleted:
                logger.info(f"✅ [WorkflowDB] 删除工作流: {workflow_id} (User: {user_id})")
            else:
                logger.warning(f"⚠️ [WorkflowDB] 工作流不存在或无权删除: {workflow_id} (User: {user_id})")
            
            return deleted
    
    # ==================== job_history 操作 ====================
    
    def create_job(
        self,
        user_id: str,
        workflow_snapshot: Dict[str, Any],
        status: str = "pending"
    ) -> int:
        """
        创建任务记录
        
        Args:
            user_id: 用户ID
            workflow_snapshot: 工作流快照（字典）
            status: 任务状态（pending, running, completed, failed）
            
        Returns:
            任务ID
        """
        with self._get_connection() as conn:
            cursor = conn.cursor()
            cursor.execute("""
                INSERT INTO job_history (user_id, workflow_snapshot, status)
                VALUES (?, ?, ?)
            """, (
                user_id,
                json.dumps(workflow_snapshot, ensure_ascii=False),
                status
            ))
            conn.commit()
            job_id = cursor.lastrowid
            logger.info(f"✅ [WorkflowDB] 创建任务: {job_id} (User: {user_id}, Status: {status})")
            return job_id
    
    def update_job(
        self,
        job_id: int,
        status: Optional[str] = None,
        result_summary: Optional[Dict[str, Any]] = None
    ) -> bool:
        """
        更新任务状态
        
        Args:
            job_id: 任务ID
            status: 新状态（可选）
            result_summary: 结果摘要（可选）
            
        Returns:
            True 如果更新成功，False 否则
        """
        with self._get_connection() as conn:
            cursor = conn.cursor()
            
            updates = []
            params = []
            
            if status:
                updates.append("status = ?")
                params.append(status)
            
            if result_summary:
                updates.append("result_summary = ?")
                params.append(json.dumps(result_summary, ensure_ascii=False))
            
            if status in ["completed", "failed"]:
                updates.append("completed_at = ?")
                params.append(datetime.now().isoformat())
            
            if not updates:
                return False
            
            params.append(job_id)
            cursor.execute(f"""
                UPDATE job_history
                SET {', '.join(updates)}
                WHERE id = ?
            """, params)
            conn.commit()
            
            updated = cursor.rowcount > 0
            if updated:
                logger.info(f"✅ [WorkflowDB] 更新任务: {job_id} (Status: {status})")
            return updated
    
    def get_job(self, job_id: int) -> Optional[Dict[str, Any]]:
        """
        获取任务
        
        Args:
            job_id: 任务ID
            
        Returns:
            任务字典，如果不存在返回 None
        """
        with self._get_connection() as conn:
            cursor = conn.cursor()
            cursor.execute("""
                SELECT * FROM job_history
                WHERE id = ?
            """, (job_id,))
            
            row = cursor.fetchone()
            if row:
                return {
                    "id": row["id"],
                    "user_id": row["user_id"],
                    "workflow_snapshot": json.loads(row["workflow_snapshot"]),
                    "status": row["status"],
                    "result_summary": json.loads(row["result_summary"]) if row["result_summary"] else None,
                    "created_at": row["created_at"],
                    "completed_at": row["completed_at"]
                }
            return None
    
    def list_jobs(
        self,
        user_id: str = "guest",
        status: Optional[str] = None,
        limit: int = 50
    ) -> List[Dict[str, Any]]:
        """
        列出任务历史
        
        Args:
            user_id: 用户ID
            status: 状态过滤（可选）
            limit: 返回数量限制
            
        Returns:
            任务列表
        """
        with self._get_connection() as conn:
            cursor = conn.cursor()
            
            if status:
                cursor.execute("""
                    SELECT * FROM job_history
                    WHERE user_id = ? AND status = ?
                    ORDER BY created_at DESC
                    LIMIT ?
                """, (user_id, status, limit))
            else:
                cursor.execute("""
                    SELECT * FROM job_history
                    WHERE user_id = ?
                    ORDER BY created_at DESC
                    LIMIT ?
                """, (user_id, limit))
            
            jobs = []
            for row in cursor.fetchall():
                jobs.append({
                    "id": row["id"],
                    "user_id": row["user_id"],
                    "workflow_snapshot": json.loads(row["workflow_snapshot"]),
                    "status": row["status"],
                    "result_summary": json.loads(row["result_summary"]) if row["result_summary"] else None,
                    "created_at": row["created_at"],
                    "completed_at": row["completed_at"]
                })
            
            logger.info(f"✅ [WorkflowDB] 列出任务: {len(jobs)} 个 (User: {user_id}, Status: {status or 'all'})")
            return jobs


# 全局数据库实例（单例）
_db_instance: Optional[WorkflowDB] = None


def get_db(db_path: str = "workflows.db") -> WorkflowDB:
    """
    获取数据库实例（单例）
    
    Args:
        db_path: 数据库文件路径
        
    Returns:
        WorkflowDB 实例
    """
    global _db_instance
    if _db_instance is None:
        _db_instance = WorkflowDB(db_path)
    return _db_instance

