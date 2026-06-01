"""
Phase 1 - 用户与资产中台：核心数据表 ORM 模型

5 张表：User, Session, Message, Asset, WorkflowTemplate。
- 所有归属字段 owner_id 均为 String，不设 FK，支持游客 guest_xxx 及后续一句 SQL 继承给正式用户。
- Message.content / WorkflowTemplate.config_json 使用 JSON 类型，无损兼容前端复杂嵌套结构。
"""
from datetime import datetime
from sqlalchemy import Column, Integer, String, DateTime, Text
from sqlalchemy.types import JSON

# 使用 Part 1 建立的容错连接中的 Base（同包内引用）
from .connection import Base


class User(Base):
    """用户表：正式注册用户。"""

    __tablename__ = "users"

    id = Column(Integer, primary_key=True, autoincrement=True)
    username = Column(String(255), unique=True, nullable=False, index=True)
    hashed_password = Column(String(255), nullable=False)
    role = Column(String(64), nullable=False, default="user")
    # 注册审核：pending | approved | rejected；老库迁移前可为 NULL，登录侧视 NULL 为已通过
    approval_status = Column(String(32), nullable=True, default="pending")
    email = Column(String(512), nullable=True)
    created_at = Column(DateTime, nullable=False, default=datetime.utcnow)


class Session(Base):
    """历史会话表：归属用 owner_id（username 或 guest_uuid），不设外键便于游客继承。"""

    __tablename__ = "sessions"

    id = Column(String(64), primary_key=True)  # 前端生成的会话 UUID
    owner_id = Column(String(255), nullable=False, index=True)  # username 或 guest_xxx
    title = Column(String(512), nullable=True)
    # 主页面 Composer 草稿：{ input_draft_text, input_draft_attachments }，与 execution_snapshot 同键名
    composer_draft = Column(JSON, nullable=True)
    # 会话执行状态：idle | running | completed | failed（后台任务与 SSE 重连）
    status = Column(String(32), nullable=False, default="idle")
    created_at = Column(DateTime, nullable=False, default=datetime.utcnow)


class Message(Base):
    """消息记录表：content 为 JSON，兼容前端 steps_details 等复杂嵌套结构。

    时光机全量执行快照写入 content.execution_snapshot（steps_details 含 step_result.markdown/data），
    与 content.state_snapshot 双写；无需单独 ALTER 表（JSON 列容量由库配置决定）。
    """

    __tablename__ = "messages"

    id = Column(Integer, primary_key=True, autoincrement=True)
    session_id = Column(String(64), nullable=False, index=True)
    role = Column(String(32), nullable=False)  # user / agent
    content = Column(JSON, nullable=True)  # 复杂嵌套 JSON，无损存取
    created_at = Column(DateTime, nullable=False, default=datetime.utcnow)


class Asset(Base):
    """数据资产表：上传文件元数据，owner_id 为 String 支持游客与继承。"""

    __tablename__ = "assets"

    id = Column(Integer, primary_key=True, autoincrement=True)
    owner_id = Column(String(255), nullable=False, index=True)
    file_name = Column(String(512), nullable=False)
    file_path = Column(String(1024), nullable=False)  # 服务器物理路径
    modality = Column(String(64), nullable=True)  # rna, spatial, metabolomics, radiomics 等
    created_at = Column(DateTime, nullable=False, default=datetime.utcnow)


class WorkflowTemplate(Base):
    """工作流收藏表：收藏的流程配置，config_json 为完整参数配置。"""

    __tablename__ = "workflow_templates"

    id = Column(Integer, primary_key=True, autoincrement=True)
    owner_id = Column(String(255), nullable=False, index=True)
    name = Column(String(512), nullable=False)
    config_json = Column(JSON, nullable=True)  # 工作流参数配置，复杂嵌套无损
    created_at = Column(DateTime, nullable=False, default=datetime.utcnow)


class Skill(Base):
    """UGC 技能表：用户上传技能，待管理员审核后展示。系统技能 name 唯一约束防并发重复。"""

    __tablename__ = "skills"

    id = Column(Integer, primary_key=True, autoincrement=True)
    name = Column(String(255), unique=True, index=True, nullable=False)
    description = Column(Text, nullable=True)
    main_category = Column(String(128), nullable=True)
    sub_category = Column(String(128), nullable=True)
    prompt_template = Column(Text, nullable=True)
    author_id = Column(String(255), nullable=False, index=True)  # username
    status = Column(String(32), nullable=False, default="pending")  # pending | approved | rejected
    created_at = Column(DateTime, nullable=False, default=datetime.utcnow)


class UserSavedSkill(Base):
    """用户收藏技能关联表：轻量级，owner_id 为 String 支持游客与继承。"""

    __tablename__ = "user_saved_skills"

    id = Column(Integer, primary_key=True, autoincrement=True)
    owner_id = Column(String(255), nullable=False, index=True)
    skill_id = Column(Integer, nullable=False, index=True)  # 逻辑关联 skills.id，不设 FK 保持与项目风格一致
    created_at = Column(DateTime, nullable=False, default=datetime.utcnow)


class UserFeedback(Base):
    """用户反馈（帮助台）：支持游客与注册用户 owner_id。"""

    __tablename__ = "user_feedbacks"

    id = Column(Integer, primary_key=True, autoincrement=True)
    owner_id = Column(String(255), nullable=False, index=True)
    feedback_type = Column(String(64), nullable=False, default="other")
    content = Column(Text, nullable=False)
    error_context = Column(Text, nullable=True)
    client_timestamp = Column(String(64), nullable=True)
    created_at = Column(DateTime, nullable=False, default=datetime.utcnow)


class UserNotification(Base):
    """用户站内通知：技能审核进度、上架结果等。"""

    __tablename__ = "user_notifications"

    id = Column(Integer, primary_key=True, autoincrement=True)
    user_id = Column(String(255), nullable=False, index=True)  # username / owner_id
    type = Column(String(64), nullable=False, default="system_notice")  # skill_review | skill_published | system_notice
    title = Column(String(512), nullable=False)
    content = Column(Text, nullable=False)
    is_read = Column(Integer, nullable=False, default=0)  # 0/1 兼容 SQLite/MySQL
    created_at = Column(DateTime, nullable=False, default=datetime.utcnow)
