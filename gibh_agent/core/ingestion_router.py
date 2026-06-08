# -*- coding: utf-8 -*-
"""三轨合一入库物理路由（Local Volume / HPC SFTP / API multipart）。"""
from __future__ import annotations

import logging
import shutil
import socket
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any, Dict, Optional

logger = logging.getLogger(__name__)


class IngestionRouterError(Exception):
    """入库路由失败；message 面向用户，details 供日志。"""

    def __init__(
        self,
        message: str,
        *,
        mount_type: Optional[str] = None,
        details: Optional[str] = None,
    ) -> None:
        super().__init__(message)
        self.message = message
        self.mount_type = mount_type
        self.details = details


def deliver_archive(
    archive_path: str,
    settings: Dict[str, Any],
    *,
    session_id: str = "",
    owner_id: str = "",
) -> Dict[str, Any]:
    """根据用户 mount 配置，将归档包推送到目标位置。"""
    return IngestionRouter().deliver(
        archive_path,
        settings,
        session_id=session_id,
        owner_id=owner_id,
    )


class IngestionStrategy(ABC):
    mount_type: str

    @abstractmethod
    def deliver(
        self,
        archive_path: Path,
        settings: Dict[str, Any],
        *,
        session_id: str,
        owner_id: str,
    ) -> Dict[str, Any]:
        ...


class LocalVolumeStrategy(IngestionStrategy):
    mount_type = "local_volume"

    def deliver(
        self,
        archive_path: Path,
        settings: Dict[str, Any],
        *,
        session_id: str,
        owner_id: str,
    ) -> Dict[str, Any]:
        cfg = settings.get("local_volume") or {}
        mount_path = str(cfg.get("mount_path") or "").strip()
        if not mount_path:
            raise IngestionRouterError(
                "local_volume.mount_path 未配置",
                mount_type=self.mount_type,
            )

        base = Path(mount_path).expanduser()
        if not base.exists():
            raise IngestionRouterError(
                f"本地挂载路径不存在: {base}",
                mount_type=self.mount_type,
            )
        if not base.is_dir():
            raise IngestionRouterError(
                f"本地挂载路径不是目录: {base}",
                mount_type=self.mount_type,
            )

        sub = session_id.strip() or (f"owner-{owner_id[:8]}" if owner_id else "anonymous")
        dest_dir = base / "gibh_ingest" / sub
        try:
            dest_dir.mkdir(parents=True, exist_ok=True)
        except PermissionError as exc:
            raise IngestionRouterError(
                f"无法在挂载目录创建子目录 {dest_dir} — 权限不足: {exc}",
                mount_type=self.mount_type,
                details=str(exc),
            ) from exc
        except OSError as exc:
            raise IngestionRouterError(
                f"创建目标目录失败 {dest_dir}: {exc}",
                mount_type=self.mount_type,
                details=str(exc),
            ) from exc

        dest_file = dest_dir / archive_path.name
        try:
            shutil.copy2(str(archive_path), str(dest_file))
        except PermissionError as exc:
            raise IngestionRouterError(
                f"复制归档到 {dest_file} 失败 — 权限不足: {exc}",
                mount_type=self.mount_type,
                details=str(exc),
            ) from exc
        except OSError as exc:
            raise IngestionRouterError(
                f"复制归档失败 {archive_path} → {dest_file}: {exc}",
                mount_type=self.mount_type,
                details=str(exc),
            ) from exc

        return {
            "strategy": self.mount_type,
            "destination": str(dest_file.resolve()),
            "destination_dir": str(dest_dir.resolve()),
        }


class HpcSlurmStrategy(IngestionStrategy):
    mount_type = "hpc_slurm"

    def deliver(
        self,
        archive_path: Path,
        settings: Dict[str, Any],
        *,
        session_id: str,
        owner_id: str,
    ) -> Dict[str, Any]:
        try:
            import paramiko
        except ImportError as exc:
            raise IngestionRouterError(
                "HPC/SFTP 入库需要 paramiko，请安装 requirements.txt 中的依赖",
                mount_type=self.mount_type,
                details=str(exc),
            ) from exc

        cfg = settings.get("hpc_slurm") or {}
        host = str(cfg.get("host") or "").strip()
        port = int(cfg.get("port") or 22)
        username = str(cfg.get("username") or "").strip()
        base_path = str(cfg.get("base_path") or "").strip().rstrip("/")

        if not host:
            raise IngestionRouterError("hpc_slurm.host 未配置", mount_type=self.mount_type)
        if not username:
            raise IngestionRouterError(
                "hpc_slurm.username 未配置（SFTP 需要 SSH 用户名）",
                mount_type=self.mount_type,
            )

        sub = session_id.strip() or (f"owner-{owner_id[:8]}" if owner_id else "anonymous")
        remote_dir = f"{base_path}/gibh_ingest/{sub}" if base_path else f"gibh_ingest/{sub}"
        remote_file = f"{remote_dir}/{archive_path.name}"

        transport = None
        sftp = None
        stat = None

        try:
            transport = paramiko.Transport((host, port))
            transport.connect(
                username=username,
                allow_agent=True,
                look_for_keys=True,
            )
        except paramiko.AuthenticationException as exc:
            raise IngestionRouterError(
                f"SSH 认证失败 ({username}@{host}:{port}) — 请确认 SSH 密钥或凭据",
                mount_type=self.mount_type,
                details=str(exc),
            ) from exc
        except paramiko.SSHException as exc:
            raise IngestionRouterError(
                f"SSH 连接失败 {host}:{port} — {exc}",
                mount_type=self.mount_type,
                details=str(exc),
            ) from exc
        except socket.timeout as exc:
            raise IngestionRouterError(
                f"SSH 连接超时 {host}:{port}",
                mount_type=self.mount_type,
                details=str(exc),
            ) from exc
        except OSError as exc:
            raise IngestionRouterError(
                f"无法连接 HPC 主机 {host}:{port} — {exc}",
                mount_type=self.mount_type,
                details=str(exc),
            ) from exc

        try:
            sftp = paramiko.SFTPClient.from_transport(transport)
            if sftp is None:
                raise IngestionRouterError(
                    "无法建立 SFTP 会话",
                    mount_type=self.mount_type,
                )
            self._mkdir_p(sftp, remote_dir)
            sftp.put(str(archive_path), remote_file)
            stat = sftp.stat(remote_file)
        except IOError as exc:
            raise IngestionRouterError(
                f"SFTP 上传失败 → {remote_file} — {exc}",
                mount_type=self.mount_type,
                details=str(exc),
            ) from exc
        except FileNotFoundError as exc:
            raise IngestionRouterError(
                f"SFTP 远程路径不可写或不存在: {remote_dir} — {exc}",
                mount_type=self.mount_type,
                details=str(exc),
            ) from exc
        except IngestionRouterError:
            raise
        except Exception as exc:
            raise IngestionRouterError(
                f"SFTP 传输异常 ({host}): {exc}",
                mount_type=self.mount_type,
                details=str(exc),
            ) from exc
        finally:
            if sftp is not None:
                try:
                    sftp.close()
                except Exception:
                    pass
            if transport is not None:
                try:
                    transport.close()
                except Exception:
                    pass

        return {
            "strategy": self.mount_type,
            "destination": remote_file,
            "destination_host": f"{username}@{host}:{port}",
            "remote_size_bytes": getattr(stat, "st_size", None) if stat else None,
        }

    @staticmethod
    def _mkdir_p(sftp, remote_dir: str) -> None:
        parts = [p for p in remote_dir.replace("\\", "/").split("/") if p]
        cur = ""
        for part in parts:
            cur = f"{cur}/{part}" if cur else part
            try:
                sftp.stat(cur)
            except IOError:
                sftp.mkdir(cur)


class ApiUrlStrategy(IngestionStrategy):
    mount_type = "api_url"

    def deliver(
        self,
        archive_path: Path,
        settings: Dict[str, Any],
        *,
        session_id: str,
        owner_id: str,
    ) -> Dict[str, Any]:
        import requests

        cfg = settings.get("api_url") or {}
        endpoint = str(cfg.get("endpoint") or "").strip()
        token = str(cfg.get("token") or "").strip()
        if not endpoint:
            raise IngestionRouterError("api_url.endpoint 未配置", mount_type=self.mount_type)

        headers: Dict[str, str] = {}
        if token:
            headers["Authorization"] = f"Bearer {token}"

        data = {
            "session_id": session_id or "",
            "owner_id": owner_id or "",
        }

        try:
            with open(archive_path, "rb") as fh:
                files = {"archive": (archive_path.name, fh, "application/gzip")}
                resp = requests.post(
                    endpoint,
                    headers=headers,
                    data=data,
                    files=files,
                    timeout=120,
                )
        except requests.RequestException as exc:
            raise IngestionRouterError(
                f"API 上传请求失败: {exc}",
                mount_type=self.mount_type,
                details=str(exc),
            ) from exc

        if resp.status_code >= 400:
            body_preview = (resp.text or "")[:500]
            raise IngestionRouterError(
                f"API 入库端点返回 HTTP {resp.status_code}: {body_preview}",
                mount_type=self.mount_type,
                details=body_preview,
            )

        try:
            resp_json = resp.json() if resp.content else {}
        except ValueError:
            resp_json = {"raw": (resp.text or "")[:2000]}

        return {
            "strategy": self.mount_type,
            "destination": endpoint,
            "http_status": resp.status_code,
            "response": resp_json if isinstance(resp_json, dict) else {},
        }


class IngestionRouter:
    _STRATEGIES: Dict[str, IngestionStrategy] = {
        "local_volume": LocalVolumeStrategy(),
        "hpc_slurm": HpcSlurmStrategy(),
        "api_url": ApiUrlStrategy(),
    }

    def deliver(
        self,
        archive_path: str,
        settings: Dict[str, Any],
        *,
        session_id: str = "",
        owner_id: str = "",
    ) -> Dict[str, Any]:
        path = Path(str(archive_path)).expanduser()
        if not path.is_file():
            raise IngestionRouterError(f"归档文件不存在: {path}")

        mount_type = str(settings.get("mount_type") or "local_volume")
        strategy = self._STRATEGIES.get(mount_type)
        if not strategy:
            raise IngestionRouterError(f"未知 mount_type: {mount_type}", mount_type=mount_type)

        logger.info(
            "[IngestionRouter] deliver mount_type=%s archive=%s session=%s",
            mount_type,
            path.name,
            session_id or "na",
        )
        result = strategy.deliver(
            path,
            settings,
            session_id=session_id,
            owner_id=owner_id,
        )
        result["mount_type"] = mount_type
        result["archive_name"] = path.name
        return result
