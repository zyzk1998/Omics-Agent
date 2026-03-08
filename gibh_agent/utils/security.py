"""
Digital signature utilities for uploaded files.

Uses BLAKE3 (chunked hashing) and Ed25519 (signing) for integrity and non-repudiation.
Sidecar format: {file_path}.sig (JSON).
"""

from __future__ import annotations

import base64
import json
import logging
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)

# Chunk size for hashing large bio-files (1 MiB)
CHUNK_SIZE = 1024 * 1024

try:
    import blake3
except ImportError:
    blake3 = None  # type: ignore[assignment]

try:
    from nacl.signing import SigningKey, VerifyKey
    from nacl.encoding import RawEncoder
except ImportError:
    SigningKey = None  # type: ignore[assignment, misc]
    VerifyKey = None  # type: ignore[assignment, misc]
    RawEncoder = None  # type: ignore[assignment, misc]


def compute_blake3_hash(file_path: Path, chunk_size: int = CHUNK_SIZE) -> str:
    """
    Compute BLAKE3 hash of a file using chunked reading (1 MiB default).

    Suitable for large bio-files without loading entire file into memory.

    Args:
        file_path: Path to the file.
        chunk_size: Read chunk size in bytes (default 1 MiB).

    Returns:
        Lowercase hex digest string. Empty string "" on blake3 unavailable or I/O error.
    """
    if blake3 is None:
        logger.warning("blake3 not installed; cannot compute hash")
        return ""
    if not file_path.is_file():
        logger.warning("Not a file or missing: %s", file_path)
        return ""
    try:
        hasher = blake3.blake3()
        with open(file_path, "rb") as f:
            while True:
                chunk = f.read(chunk_size)
                if not chunk:
                    break
                hasher.update(chunk)
        return hasher.digest().hex()
    except OSError as e:
        logger.exception("Failed to hash file %s: %s", file_path, e)
        return ""


def sign_file(
    file_path: Path,
    private_key_b64: str,
    sig_path: Optional[Path] = None,
) -> bool:
    """
    Compute BLAKE3 hash of the file, sign it with Ed25519, write sidecar .sig.

    Sidecar file: {original_filename}.sig in the same directory.
    Content (JSON): {"blake3": "...", "signature_b64": "...", "timestamp": "ISO8601"}.

    Args:
        file_path: Path to the file to sign.
        private_key_b64: Base64-encoded Ed25519 private key (32 bytes).
        sig_path: Override path for sidecar; default is {file_path}.sig.

    Returns:
        True if sidecar was written successfully; False otherwise. IO errors handled gracefully.
    """
    if blake3 is None or SigningKey is None or RawEncoder is None:
        logger.warning("blake3 or pynacl not installed")
        return False
    if not file_path.is_file():
        logger.warning("Not a file or missing: %s", file_path)
        return False
    sig_path = sig_path or Path(str(file_path) + ".sig")
    try:
        digest_hex = compute_blake3_hash(file_path)
        if not digest_hex:
            logger.warning("Failed to compute BLAKE3 hash for %s", file_path)
            return False
        key_b64 = (private_key_b64 or "").strip()
        if not key_b64:
            logger.warning("Private key is empty")
            return False
        key_bytes = base64.standard_b64decode(key_b64)
        sk = SigningKey(key_bytes)
        message = digest_hex.encode("utf-8")
        sig_bytes = sk.sign(message, encoder=RawEncoder).signature
        payload = {
            "blake3": digest_hex,
            "signature_b64": base64.standard_b64encode(sig_bytes).decode("ascii"),
            "timestamp": datetime.now(timezone.utc).isoformat(),
        }
        sig_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
        logger.info("Signed file: %s -> %s", file_path, sig_path)
        return True
    except (ValueError, TypeError) as e:
        logger.warning("Invalid key or encoding for signing %s: %s", file_path, e)
        return False
    except OSError as e:
        logger.exception("I/O error while signing %s: %s", file_path, e)
        return False
    except Exception as e:
        logger.exception("Unexpected error signing %s: %s", file_path, e)
        return False


def verify_file_signature(
    file_path: Path,
    public_key_b64: str,
    sig_path: Optional[Path] = None,
) -> bool:
    """
    Verify Ed25519 signature in sidecar against current file BLAKE3 hash.

    Args:
        file_path: Path to the file.
        public_key_b64: Base64-encoded Ed25519 public key (32 bytes).
        sig_path: Override path for sidecar; default is {file_path}.sig.

    Returns:
        True if sidecar exists, hash matches, and signature verifies; False otherwise.
    """
    if blake3 is None or VerifyKey is None:
        logger.debug("blake3 or pynacl not installed")
        return False
    if not file_path.is_file():
        return False
    sig_path = sig_path or Path(str(file_path) + ".sig")
    if not sig_path.exists():
        logger.debug("No signature file: %s", sig_path)
        return False
    try:
        payload = json.loads(sig_path.read_text(encoding="utf-8"))
        digest_hex = payload.get("blake3")
        sig_b64 = payload.get("signature_b64")
        if not digest_hex or not sig_b64:
            return False
        current_hex = compute_blake3_hash(file_path)
        if current_hex != digest_hex:
            logger.warning("Hash mismatch for %s", file_path)
            return False
        key_b64 = public_key_b64.strip()
        if not key_b64:
            return False
        vk = VerifyKey(base64.standard_b64decode(key_b64))
        vk.verify(digest_hex.encode("utf-8"), base64.standard_b64decode(sig_b64))
        return True
    except (ValueError, TypeError, KeyError) as e:
        logger.debug("Invalid signature or key for %s: %s", file_path, e)
        return False
    except OSError as e:
        logger.debug("I/O error verifying %s: %s", file_path, e)
        return False
    except Exception as e:
        logger.debug("Verification failed for %s: %s", file_path, e)
        return False
