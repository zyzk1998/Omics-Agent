"""
Auto-initializing signing keys for file signatures (BLAKE3 + Ed25519).

On first use: if the keys file does not exist, generate a new Ed25519 key pair,
persist to a JSON file on the mounted volume, and return the keys.
On subsequent runs: load keys from the file. No manual .env or copy-paste required.
"""

from __future__ import annotations

import base64
import json
import logging
import os
from pathlib import Path
from typing import Optional, Tuple

logger = logging.getLogger(__name__)

# Persisted under mounted volume so keys survive container restarts (e.g. /app/data)
DEFAULT_KEYS_PATH = os.getenv("SECURITY_KEYS_PATH", "/app/data/security/signing_keys.json")

_keys_cache: Optional[Tuple[str, str]] = None


def _ensure_keys_file() -> Path:
    path = Path(DEFAULT_KEYS_PATH)
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


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
    path = Path(DEFAULT_KEYS_PATH)
    if not path.is_file():
        return None
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
