# -*- coding: utf-8 -*-
"""双边站内通知：用户提交 → 管理员；管理员处理 → 用户。"""
from __future__ import annotations

import tempfile
from pathlib import Path

import pytest

from gibh_agent.db.connection import Base, SessionLocal, engine
from gibh_agent.db.models import Skill as SkillModel, User, UserFeedback, UserNotification
from gibh_agent.plugin_system.registry import DynamicSkillPlugin, persist_plugin


@pytest.fixture(autouse=True)
def _ensure_tables():
    Base.metadata.create_all(bind=engine)
    yield


def _add_admin(db, username="admin_notify"):
    u = User(username=username, hashed_password="x", role="admin", approval_status="approved")
    db.add(u)
    db.commit()
    return u


def test_registration_notifies_admins():
    db = SessionLocal()
    try:
        _add_admin(db)
        from gibh_agent.api.routers.auth import RegisterBody, register

        register(RegisterBody(username="newbie@test.com", password="secret123"), db=db)
        notes = (
            db.query(UserNotification)
            .filter(UserNotification.user_id == "admin_notify")
            .all()
        )
        assert any(n.type == "admin_registration_pending" for n in notes)
    finally:
        db.query(UserNotification).delete()
        db.query(User).filter(User.username.in_(["admin_notify", "newbie@test.com"])).delete()
        db.commit()
        db.close()


def test_ugc_skill_review_notifies_author():
    db = SessionLocal()
    try:
        admin = _add_admin(db, "admin_ugc")
        skill = SkillModel(
            name="ugc_notify_skill",
            author_id="author@test.com",
            status="pending",
        )
        db.add(skill)
        db.commit()
        db.refresh(skill)
        from gibh_agent.api.routers.admin import SkillReviewAction, review_skill_by_action

        review_skill_by_action(
            skill.id,
            SkillReviewAction(action="approve"),
            current_user=admin,
            db=db,
        )
        notes = (
            db.query(UserNotification)
            .filter(UserNotification.user_id == "author@test.com")
            .all()
        )
        assert len(notes) == 1
        assert notes[0].type == "ugc_skill_approved"
    finally:
        db.query(UserNotification).delete()
        db.query(SkillModel).filter(SkillModel.name == "ugc_notify_skill").delete()
        db.query(User).filter(User.username == "admin_ugc").delete()
        db.commit()
        db.close()


def test_feedback_submit_notifies_admins():
    db = SessionLocal()
    try:
        _add_admin(db, "admin_fb")
        from gibh_agent.api.routers.feedback import FeedbackCreate, submit_feedback

        submit_feedback(
            FeedbackCreate(type="issue", content="按钮点不动"),
            owner_id="guest-uuid-1",
            db=db,
        )
        notes = (
            db.query(UserNotification)
            .filter(UserNotification.user_id == "admin_fb")
            .all()
        )
        assert any(n.type == "admin_feedback_new" for n in notes)
    finally:
        db.query(UserNotification).delete()
        db.query(UserFeedback).delete()
        db.query(User).filter(User.username == "admin_fb").delete()
        db.commit()
        db.close()


def test_feedback_resolve_notifies_user():
    db = SessionLocal()
    try:
        _add_admin(db, "admin_fb2")
        row = UserFeedback(
            owner_id="user_fb@test.com",
            feedback_type="issue",
            content="test",
            status="open",
        )
        db.add(row)
        db.commit()
        db.refresh(row)
        from gibh_agent.api.routers.feedback import FeedbackResolveAction, resolve_feedback

        resolve_feedback(
            row.id,
            FeedbackResolveAction(action="resolve"),
            db=db,
            _admin=db.query(User).filter(User.username == "admin_fb2").first(),
        )
        notes = (
            db.query(UserNotification)
            .filter(UserNotification.user_id == "user_fb@test.com")
            .all()
        )
        assert len(notes) == 1
        assert notes[0].type == "feedback_resolved"
    finally:
        db.query(UserNotification).delete()
        db.query(UserFeedback).delete()
        db.query(User).filter(User.username == "admin_fb2").delete()
        db.commit()
        db.close()


def test_plugin_upload_notifies_admins():
    db = SessionLocal()
    try:
        _add_admin(db, "admin_plg")
        user = User(username="uploader@test.com", hashed_password="x", role="user", approval_status="approved")
        db.add(user)
        db.commit()
        staging = Path(tempfile.mkdtemp()) / "staging"
        staging.mkdir(parents=True)
        row = persist_plugin(
            db,
            name="plugin_notify_test",
            display_name="Plugin Notify",
            description="d",
            parameters_schema={"type": "object", "properties": {}},
            skill_type="prompt",
            script_path=None,
            extract_dir=str(staging),
            author_id=user.username,
            status="pending",
        )
        from gibh_agent.core.user_notifications import notify_all_admins

        notify_all_admins(
            db,
            ntype="admin_plugin_submitted",
            title="新安装技能待审核",
            content=f"plugin_id={row.id}",
            commit=True,
        )
        notes = (
            db.query(UserNotification)
            .filter(UserNotification.user_id == "admin_plg")
            .all()
        )
        assert any(n.type == "admin_plugin_submitted" for n in notes)
    finally:
        db.query(UserNotification).delete()
        db.query(DynamicSkillPlugin).filter(DynamicSkillPlugin.name == "plugin_notify_test").delete()
        db.query(User).filter(User.username.in_(["admin_plg", "uploader@test.com"])).delete()
        db.commit()
        db.close()
