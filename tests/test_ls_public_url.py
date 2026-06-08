# -*- coding: utf-8 -*-
"""Label Studio 浏览器基址 / 同源反代解析单元测试。"""
from __future__ import annotations

from gibh_agent.utils.ls_public_url import (
    resolve_ls_browser_project_url,
    resolve_ls_public_base_url,
)


def test_same_origin_proxy_default_relative_base(monkeypatch):
    monkeypatch.delenv("LABEL_STUDIO_PUBLIC_URL", raising=False)
    monkeypatch.setenv("LABEL_STUDIO_SAME_ORIGIN_PROXY", "1")
    monkeypatch.setenv("OMICS_AGENT_WEB_URL", "http://host.docker.internal:8018")
    assert resolve_ls_public_base_url() == "/label-studio"
    assert resolve_ls_browser_project_url(7) == "/label-studio/projects/7/data"


def test_explicit_public_url_wins(monkeypatch):
    monkeypatch.setenv("LABEL_STUDIO_PUBLIC_URL", "http://example.com:9000")
    monkeypatch.setenv("LABEL_STUDIO_SAME_ORIGIN_PROXY", "1")
    assert resolve_ls_public_base_url() == "http://example.com:9000"
    assert resolve_ls_browser_project_url(3) == "http://example.com:9000/projects/3/data"


def test_same_origin_proxy_always_relative_for_browser(monkeypatch):
    monkeypatch.delenv("LABEL_STUDIO_PUBLIC_URL", raising=False)
    monkeypatch.setenv("LABEL_STUDIO_SAME_ORIGIN_PROXY", "1")
    monkeypatch.setenv("OMICS_AGENT_WEB_URL", "http://127.0.0.1:8028")
    assert resolve_ls_public_base_url() == "/label-studio"
    assert resolve_ls_browser_project_url(12) == "/label-studio/projects/12/data"


def test_legacy_direct_port_when_proxy_disabled(monkeypatch):
    monkeypatch.delenv("LABEL_STUDIO_PUBLIC_URL", raising=False)
    monkeypatch.setenv("LABEL_STUDIO_SAME_ORIGIN_PROXY", "0")
    monkeypatch.setenv("LABEL_STUDIO_URL", "http://127.0.0.1:8082")
    assert resolve_ls_public_base_url() == "http://127.0.0.1:8082"
