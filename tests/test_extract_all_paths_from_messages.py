# -*- coding: utf-8 -*-
from gibh_agent.core.storage.session_file_index import extract_all_paths_from_messages


def test_extract_paths_from_snapshot_content():
    messages = [
        {
            "role": "assistant",
            "content": {
                "image_path": "/assets/images/demos/corpus/test_corpus_umap.png",
                "text": "done",
            },
        }
    ]
    paths = extract_all_paths_from_messages(messages)
    assert "/assets/images/demos/corpus/test_corpus_umap.png" in paths


def test_extract_paths_from_tool_calls_and_attachments():
    messages = [
        {
            "role": "user",
            "content": {
                "attachments": [{"path": "C:/data/matrix.mtx"}],
            },
            "tool_calls": [
                {
                    "function": {
                        "arguments": '{"file_path": "C:/data/features.tsv"}',
                    }
                }
            ],
            "metadata": {"mask_path": "C:/data/mask.nii.gz"},
        }
    ]
    paths = extract_all_paths_from_messages(messages)
    assert "C:/data/matrix.mtx" in paths
    assert "C:/data/features.tsv" in paths
    assert "C:/data/mask.nii.gz" in paths
