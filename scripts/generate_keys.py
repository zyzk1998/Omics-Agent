#!/usr/bin/env python3
"""
Generate Ed25519 key pair for file signing (BLAKE3 + Ed25519 sidecar .sig).

Prints Base64-encoded private and public keys to stdout so they can be
added to .env as FILE_SIGNING_PRIVATE_KEY and FILE_SIGNING_PUBLIC_KEY.

Usage:
    python scripts/generate_keys.py

Then add to .env (or Docker/CI secrets):
    FILE_SIGNING_PRIVATE_KEY=<private_b64>
    FILE_SIGNING_PUBLIC_KEY=<public_b64>
"""

from __future__ import annotations

import base64
import sys


def main() -> int:
    try:
        from nacl.signing import SigningKey
    except ImportError:
        print("ERROR: pynacl is required. Install with: pip install pynacl", file=sys.stderr)
        return 1
    sk = SigningKey.generate()
    private_b64 = base64.standard_b64encode(bytes(sk)).decode("ascii")
    public_b64 = base64.standard_b64encode(bytes(sk.verify_key)).decode("ascii")
    print("# Add these to .env (or your environment / Docker secrets)")
    print("FILE_SIGNING_PRIVATE_KEY=" + private_b64)
    print("FILE_SIGNING_PUBLIC_KEY=" + public_b64)
    return 0


if __name__ == "__main__":
    sys.exit(main())
