# -*- coding: utf-8 -*-
"""messages 轻量查询单元测试。"""
from __future__ import annotations

from types import SimpleNamespace

from gibh_agent.db.message_queries import get_latest_agent_message_id, load_latest_agent_snapshot
from gibh_agent.db.models import Message as MessageModel, Session as SessionModel


class _IdQuery:
    def __init__(self, row_id):
        self._row_id = row_id

    def filter(self, *args, **kwargs):
        return self

    def order_by(self, *args, **kwargs):
        return self

    def limit(self, n):
        return self

    def scalar(self):
        return self._row_id


class _FullRowQuery:
    def __init__(self, msg):
        self._msg = msg

    def filter(self, *args, **kwargs):
        return self

    def first(self):
        return self._msg


class _SessionQuery:
    def __init__(self, session):
        self._session = session

    def filter(self, *args, **kwargs):
        return self

    def first(self):
        return self._session


def test_get_latest_agent_message_id_returns_scalar():
    class _Db:
        def query(self, model):
            if model is MessageModel.id:
                return _IdQuery(42)
            raise AssertionError(f"unexpected query model: {model!r}")

    assert get_latest_agent_message_id(_Db(), "sess-1") == 42


def test_load_latest_agent_snapshot_checks_owner():
    session = SimpleNamespace(id="sess-1", owner_id="owner-a")
    msg_obj = SimpleNamespace(id=7, content={"state_snapshot": {"text": "ok"}})

    class _Db:
        def query(self, model):
            if model is SessionModel:
                return _SessionQuery(session)
            if model is MessageModel.id:
                return _IdQuery(7)
            if model is MessageModel:
                return _FullRowQuery(msg_obj)
            raise AssertionError(f"unexpected query model: {model!r}")

    msg, snap = load_latest_agent_snapshot(_Db(), "sess-1", "owner-a")
    assert msg is msg_obj
    assert snap.get("text") == "ok"

    msg2, snap2 = load_latest_agent_snapshot(_Db(), "sess-1", "owner-b")
    assert msg2 is None
    assert snap2 == {}
