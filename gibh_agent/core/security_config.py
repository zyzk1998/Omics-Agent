"""
Auto-initializing signing keys for file signatures (BLAKE3 + Ed25519).

On first use: if the keys file does not exist, generate a new Ed25519 key pair,
persist to a JSON file on the mounted volume, and return the keys.
On subsequent runs: load keys from the file. No manual .env or copy-paste required.

宿主机开发：默认 `/app/data/security/` 往往不可写；会按候选路径依次回退到用户缓存目录，
避免刷屏告警并确保密钥可持久化。
"""

from __future__ import annotations

import base64
import json
import logging
import os
from pathlib import Path
from typing import List, Optional, Tuple

logger = logging.getLogger(__name__)

# 容器内挂载卷（可被 SECURITY_KEYS_PATH 覆盖）
_ENV_DEFAULT = "/app/data/security/signing_keys.json"

_keys_cache: Optional[Tuple[str, str]] = None
_resolved_keys_file: Optional[Path] = None


def _key_path_candidates() -> List[Path]:
    """Ordered: env override -> container default -> user cache (always writable on POSIX)."""
    out: List[Path] = []
    env = (os.getenv("SECURITY_KEYS_PATH") or "").strip()
    if env:
        out.append(Path(env))
    out.append(Path(_ENV_DEFAULT))
    out.append(Path.home() / ".cache" / "gibh-agent" / "security" / "signing_keys.json")
    seen: set = set()
    uniq: List[Path] = []
    for p in out:
        try:
            key = str(p.resolve())
        except OSError:
            key = str(p)
        if key not in seen:
            seen.add(key)
            uniq.append(p)
    return uniq


def _find_existing_keys_file() -> Optional[Path]:
    for p in _key_path_candidates():
        try:
            if p.is_file():
                return p.resolve()
        except OSError:
            continue
    return None


def _allocate_keys_file_for_write() -> Path:
    """First candidate whose parent directory can be created and is writable."""
    last_err: Optional[Exception] = None
    for p in _key_path_candidates():
        try:
            p.parent.mkdir(parents=True, exist_ok=True)
            return p.resolve()
        except OSError as e:
            last_err = e
            continue
    raise RuntimeError(
        "Cannot create writable directory for signing keys; set SECURITY_KEYS_PATH to a "
        f"writable JSON path. Last error: {last_err}"
    )


def _ensure_keys_file() -> Path:
    global _resolved_keys_file
    if _resolved_keys_file is not None:
        return _resolved_keys_file
    existing = _find_existing_keys_file()
    if existing is not None:
        _resolved_keys_file = existing
        return _resolved_keys_file
    _resolved_keys_file = _allocate_keys_file_for_write()
    return _resolved_keys_file


def _generate_and_persist_keys() -> Tuple[str, str]:
    try:
        from nacl.signing import SigningKey
    except ImportError:
        raise RuntimeError("pynacl is required for auto key generation. Install: pip install pynacl")
    sk = SigningKey.generate()
    private_b64 = base64.standard_b64encode(bytes(sk)).decode("ascii")
    public_b64 = base64.standard_b64encode(bytes(sk.verify_key)).decode("ascii")
    path = _ensure_keys_file()
    payload = {
        "private_key_b64": private_b64,
        "public_key_b64": public_b64,
    }
    with open(path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)
    logger.info("Auto-generated signing keys and persisted to %s", path)
    return private_b64, public_b64


def _load_keys() -> Optional[Tuple[str, str]]:
    global _resolved_keys_file
    path = _find_existing_keys_file()
    if path is None:
        return None
    _resolved_keys_file = path
    try:
        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)
        priv = data.get("private_key_b64") or data.get("private_key")
        pub = data.get("public_key_b64") or data.get("public_key")
        if priv and pub:
            return (priv.strip(), pub.strip())
    except (json.JSONDecodeError, OSError) as e:
        logger.warning("Failed to load signing keys from %s: %s", path, e)
    return None


def get_signing_keys() -> Tuple[str, str]:
    """
    Return (private_key_b64, public_key_b64). Creates and persists keys on first run.
    """
    global _keys_cache
    if _keys_cache is not None:
        return _keys_cache
    loaded = _load_keys()
    if loaded is not None:
        _keys_cache = loaded
        return _keys_cache
    _keys_cache = _generate_and_persist_keys()
    return _keys_cache


def get_signing_private_key() -> Optional[str]:
    """
    Return the private key (Base64) for signing, or None if unavailable.
    Uses auto-generated persisted keys; no env var required.
    """
    try:
        private_b64, _ = get_signing_keys()
        return private_b64
    except Exception as e:
        logger.warning("Could not load signing private key: %s", e)
        return None


def get_signing_public_key() -> Optional[str]:
    """Return the public key (Base64) for verification, or None if unavailable."""
    try:
        _, public_b64 = get_signing_keys()
        return public_b64
    except Exception as e:
        logger.warning("Could not load signing public key: %s", e)
        return None
