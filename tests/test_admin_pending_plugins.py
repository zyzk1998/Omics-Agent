# -*- coding: utf-8 -*-
"""管理员待审核列表应包含 dynamic_skill_plugins。"""
from __future__ import annotations

import tempfile
from pathlib import Path

import pytest

from gibh_agent.db.connection import Base, SessionLocal, engine
from gibh_agent.db.models import User
from gibh_agent.plugin_system.registry import DynamicSkillPlugin, persist_plugin


@pytest.fixture(autouse=True)
def _ensure_tables():
    Base.metadata.create_all(bind=engine)
    yield


def test_list_pending_includes_plugin_submissions(monkeypatch, tmp_path):
    monkeypatch.setenv("UPLOAD_DIR", str(tmp_path))
    db = SessionLocal()
    try:
        admin = User(username="admin_test", hashed_password="x", role="admin", approval_status="approved")
        db.add(admin)
        db.commit()
        staging = tmp_path / "skills_review_staging" / "alice" / "uid-1"
        staging.mkdir(parents=True)
        persist_plugin(
            db,
            name="test_pending_plugin",
            display_name="Test Pending",
            description="awaiting review",
            parameters_schema={"type": "object", "properties": {}, "_review_meta": {"driver_type": "prompt"}},
            skill_type="prompt",
            script_path=None,
            extract_dir=str(staging),
            author_id="alice",
            status="pending",
        )
        from gibh_agent.api.routers.admin import list_pending_skills

        items = list_pending_skills(current_user=admin, db=db)
        plugin_items = [x for x in items if x.get("source_type") == "plugin"]
        assert len(plugin_items) >= 1
        assert plugin_items[0]["name"] == "Test Pending"
    finally:
        db.query(DynamicSkillPlugin).filter(DynamicSkillPlugin.name == "test_pending_plugin").delete()
        db.query(User).filter(User.username == "admin_test").delete()
        db.commit()
        db.close()
