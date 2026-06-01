# -*- coding: utf-8 -*-
"""插件审核状态机与用户通知。"""
from __future__ import annotations

import tempfile
from pathlib import Path

import pytest
from fastapi import HTTPException

from gibh_agent.db.connection import Base, SessionLocal, engine
from gibh_agent.db.models import User, UserNotification
from gibh_agent.plugin_system.registry import DynamicSkillPlugin, persist_plugin


@pytest.fixture(autouse=True)
def _ensure_tables():
    Base.metadata.create_all(bind=engine)
    yield


def _seed_plugin(db, *, status="pending", author="alice@test.com"):
    staging = Path(tempfile.mkdtemp()) / "staging"
    staging.mkdir(parents=True)
    return persist_plugin(
        db,
        name="notify_test_skill",
        display_name="Notify Test",
        description="desc",
        parameters_schema={"type": "object", "properties": {}},
        skill_type="prompt",
        script_path=None,
        extract_dir=str(staging),
        author_id=author,
        status=status,
    )


def test_plugin_review_accept_notifies_user():
    db = SessionLocal()
    try:
        admin = User(username="admin_n", hashed_password="x", role="admin", approval_status="approved")
        db.add(admin)
        db.commit()
        row = _seed_plugin(db, status="pending", author="creator@test.com")
        from gibh_agent.api.routers.admin import PluginReviewAction, review_plugin_by_action

        out = review_plugin_by_action(
            row.id,
            PluginReviewAction(action="accept"),
            current_user=admin,
            db=db,
        )
        assert out["status"] == "accepted"
        notes = (
            db.query(UserNotification)
            .filter(UserNotification.user_id == "creator@test.com")
            .all()
        )
        assert len(notes) == 1
        assert "受理" in notes[0].content
    finally:
        db.query(UserNotification).delete()
        db.query(DynamicSkillPlugin).filter(DynamicSkillPlugin.name == "notify_test_skill").delete()
        db.query(User).filter(User.username == "admin_n").delete()
        db.commit()
        db.close()


def test_plugin_review_publish_requires_accepted():
    db = SessionLocal()
    try:
        admin = User(username="admin_n2", hashed_password="x", role="admin", approval_status="approved")
        db.add(admin)
        db.commit()
        row = _seed_plugin(db, status="pending")
        from gibh_agent.api.routers.admin import PluginReviewAction, review_plugin_by_action

        with pytest.raises(HTTPException) as exc:
            review_plugin_by_action(
                row.id,
                PluginReviewAction(action="publish"),
                current_user=admin,
                db=db,
            )
        assert exc.value.status_code == 400
    finally:
        db.query(DynamicSkillPlugin).filter(DynamicSkillPlugin.name == "notify_test_skill").delete()
        db.query(User).filter(User.username == "admin_n2").delete()
        db.commit()
        db.close()
