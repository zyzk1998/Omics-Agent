# -*- coding: utf-8 -*-
"""Label Studio 同源反代重定向改写单元测试。"""
from __future__ import annotations

from unittest.mock import MagicMock

from gibh_agent.api.label_studio_proxy import rewrite_ls_proxy_location


def _req(host: str = "192.168.32.31:8028"):
    request = MagicMock()
    request.headers = {"host": host}
    return request


def test_rewrite_login_redirect_with_next():
    loc = "/user/login/?next=/projects/8/data/"
    out = rewrite_ls_proxy_location(loc, _req())
    assert out == "/label-studio/user/login/?next=%2Flabel-studio%2Fprojects%2F8%2Fdata%2F"


def test_rewrite_already_prefixed_location():
    loc = "/label-studio/projects/8/data/"
    out = rewrite_ls_proxy_location(loc, _req())
    assert out == "/label-studio/projects/8/data/"
