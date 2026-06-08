# -*- coding: utf-8 -*-
"""HTTP 上传白名单：图片与语料格式。"""
from gibh_agent.core.omics_io_registry import (
    is_upload_allowed_filename,
    upload_rejection_hint,
)


def test_png_jpg_and_corpus_formats_allowed():
    for name in (
        "test_corpus_umap.png",
        "slice.JPG",
        "photo.jpeg",
        "plot.webp",
        "pathology.TIF",
        "anim.gif",
        "icon.svg",
        "notes.md",
        "tasks.json",
        "matrix.mtx",
        "data.h5ad",
    ):
        assert is_upload_allowed_filename(name), name


def test_unknown_extension_rejected():
    assert not is_upload_allowed_filename("malware.exe")
    assert "png" in upload_rejection_hint("bad.exe") or "不允许" in upload_rejection_hint("bad.exe")
