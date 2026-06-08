# -*- coding: utf-8 -*-
"""
Label Studio REST API 客户端（Phase 1 · HITL 基础设施）。

通过 HTTP 与本地或远程 Label Studio 实例通信，不引入 label-studio-sdk 重型依赖。
支持零配置：无 LABEL_STUDIO_API_KEY 时通过预设管理员账号 session 登录（LS 1.23+ 兼容）。
"""
from __future__ import annotations

import logging
import os
import re
import time
from typing import Any, Dict, List, Optional, Union
from urllib.parse import urljoin

import requests
from requests.exceptions import ConnectionError as RequestsConnectionError
from requests.exceptions import RequestException, Timeout

logger = logging.getLogger(__name__)

DEFAULT_LS_BASE_URL = os.getenv("LABEL_STUDIO_URL", "http://127.0.0.1:8082").rstrip("/")
DEFAULT_LS_API_KEY = os.getenv("LABEL_STUDIO_API_KEY", "")
DEFAULT_LS_USERNAME = os.getenv("LABEL_STUDIO_USERNAME", "system_admin@omics.local")
DEFAULT_LS_PASSWORD = os.getenv("LABEL_STUDIO_PASSWORD", "AutoBootstrap_O1m2i3c4s!")
DEFAULT_TIMEOUT_SEC = float(os.getenv("LABEL_STUDIO_TIMEOUT_SEC", "30"))
BOOTSTRAP_MAX_RETRIES = int(os.getenv("LABEL_STUDIO_BOOTSTRAP_RETRIES", "30"))
BOOTSTRAP_RETRY_INTERVAL_SEC = float(os.getenv("LABEL_STUDIO_BOOTSTRAP_INTERVAL_SEC", "2"))

# 按 LS 基址缓存 bootstrap session（LS 1.23+ 默认禁用 legacy Token，改走 Cookie session）
_BOOTSTRAP_SESSION_CACHE: Dict[str, requests.Session] = {}

# 通用图像标注配置（可被调用方覆盖）
DEFAULT_IMAGE_LABEL_CONFIG = """<View>
  <Image name="image" value="$image"/>
  <RectangleLabels name="label" toName="image">
    <Label value="Region" background="#FF0000"/>
  </RectangleLabels>
</View>"""

_AUTH_SESSION = "session"
_AUTH_TOKEN = "token"
_AUTH_BEARER = "bearer"


class LabelStudioClientError(Exception):
    """Label Studio API 调用失败。"""

    def __init__(self, message: str, *, status_code: Optional[int] = None, response_body: Any = None):
        super().__init__(message)
        self.status_code = status_code
        self.response_body = response_body


def _extract_csrf_from_html(html: str) -> str:
    if not html:
        return ""
    m = re.search(
        r'name=["\']csrfmiddlewaretoken["\']\s+value=["\']([^"\']+)["\']',
        html,
        flags=re.IGNORECASE,
    )
    return m.group(1) if m else ""


def _parse_token_payload(data: Any) -> str:
    if not isinstance(data, dict):
        return ""
    for key in ("token", "detail", "key", "access"):
        val = data.get(key)
        if val and len(str(val).strip()) >= 8:
            return str(val).strip()
    return ""


def _is_jwt_pat(token: str) -> bool:
    raw = str(token or "").strip()
    parts = raw.split(".")
    return len(parts) == 3 and raw.startswith("eyJ")


def _legacy_token_disabled(body: Any) -> bool:
    text = body if isinstance(body, str) else str(body)
    return "legacy token authentication has been disabled" in text.lower()


def _wait_ls_health(base_url: str, timeout_sec: float) -> bool:
    """等待 LS /health 就绪。"""
    for path in ("/health", "/api/health"):
        try:
            resp = requests.get(
                f"{base_url.rstrip('/')}{path}",
                timeout=min(timeout_sec, 5.0),
            )
            if resp.status_code < 500:
                return True
        except RequestException:
            continue
    return False


def _session_login(base_url: str, username: str, password: str, timeout_sec: float) -> requests.Session:
    """Django session 登录（POST /user/login/ + CSRF）。"""
    base = base_url.rstrip("/")
    login_url = f"{base}/user/login/"
    session = requests.Session()
    page = session.get(login_url, timeout=timeout_sec)
    if page.status_code >= 400:
        raise LabelStudioClientError(
            f"Label Studio 登录页不可达 HTTP {page.status_code}: {login_url}",
            status_code=page.status_code,
        )

    csrf = session.cookies.get("csrftoken") or _extract_csrf_from_html(page.text)
    if not csrf:
        raise LabelStudioClientError("Label Studio 登录页未返回 CSRF token")

    headers = {
        "Referer": login_url,
        "X-CSRFToken": csrf,
    }
    form = {
        "email": username.strip().lower(),
        "password": password,
        "csrfmiddlewaretoken": csrf,
        "persist_session": "on",
    }
    resp = session.post(
        login_url,
        data=form,
        headers=headers,
        timeout=timeout_sec,
        allow_redirects=True,
    )
    if "sessionid" not in session.cookies:
        raise LabelStudioClientError(
            f"Label Studio 登录失败（账号 {username!r}），HTTP {resp.status_code}"
        )
    session.headers.setdefault("Content-Type", "application/json")
    return session


def _fetch_legacy_token_with_session(session: requests.Session, base_url: str, timeout_sec: float) -> str:
    """已登录 session 下抓取 legacy API token（仅 legacy_api_tokens_enabled=true 时可用）。"""
    base = base_url.rstrip("/")
    for path in ("/api/current-user/token", "/api/current-user/token/"):
        try:
            resp = session.get(f"{base}{path}", timeout=timeout_sec)
            if resp.status_code == 200:
                tok = _parse_token_payload(resp.json())
                if tok:
                    return tok
        except (RequestException, ValueError):
            continue
    return ""


def _create_pat_with_session(session: requests.Session, base_url: str, timeout_sec: float) -> str:
    """创建 Personal Access Token（refresh JWT）；409 表示已存在有效 PAT。"""
    base = base_url.rstrip("/")
    csrf = session.cookies.get("csrftoken", "")
    headers = {"Referer": base, "X-CSRFToken": csrf} if csrf else {}
    try:
        resp = session.post(f"{base}/api/token/", headers=headers, timeout=timeout_sec)
        if resp.status_code == 201:
            return _parse_token_payload(resp.json())
    except (RequestException, ValueError):
        pass
    return ""


def _refresh_bearer_access(pat_refresh: str, base_url: str, timeout_sec: float) -> str:
    base = base_url.rstrip("/")
    resp = requests.post(
        f"{base}/api/token/refresh",
        json={"refresh": pat_refresh},
        timeout=timeout_sec,
    )
    if resp.status_code >= 400:
        raise LabelStudioClientError(
            f"Label Studio PAT refresh 失败 HTTP {resp.status_code}: {resp.text[:300]}",
            status_code=resp.status_code,
            response_body=resp.text,
        )
    access = _parse_token_payload(resp.json())
    if not access:
        raise LabelStudioClientError("Label Studio PAT refresh 未返回 access token")
    return access


def _auto_bootstrap_session(
    base_url: str,
    *,
    username: Optional[str] = None,
    password: Optional[str] = None,
    timeout_sec: Optional[float] = None,
    force_refresh: bool = False,
) -> requests.Session:
    """
    零配置认证：等待 LS 健康 → Django session 登录。
    LS 1.23+ 默认禁用 legacy Token，API 调用直接使用 session Cookie。
    """
    base = base_url.rstrip("/")
    if not force_refresh and base in _BOOTSTRAP_SESSION_CACHE:
        cached = _BOOTSTRAP_SESSION_CACHE[base]
        if cached.cookies.get("sessionid"):
            return cached

    user = (username or DEFAULT_LS_USERNAME).strip()
    pwd = password or DEFAULT_LS_PASSWORD
    timeout = timeout_sec if timeout_sec is not None else DEFAULT_TIMEOUT_SEC
    last_err = "未知错误"

    for attempt in range(1, BOOTSTRAP_MAX_RETRIES + 1):
        if attempt > 1:
            time.sleep(BOOTSTRAP_RETRY_INTERVAL_SEC)
        try:
            if not _wait_ls_health(base, timeout):
                last_err = "Label Studio 健康检查未通过（/health）"
                logger.debug("LS bootstrap [%s/%s]: %s", attempt, BOOTSTRAP_MAX_RETRIES, last_err)
                continue

            session = _session_login(base, user, pwd, timeout)
            _BOOTSTRAP_SESSION_CACHE[base] = session
            logger.info(
                "Label Studio 零配置 session bootstrap 成功 base=%s user=%s",
                base,
                user,
            )
            return session
        except LabelStudioClientError as exc:
            last_err = str(exc)
            logger.debug("LS bootstrap [%s/%s]: %s", attempt, BOOTSTRAP_MAX_RETRIES, last_err)
        except RequestException as exc:
            last_err = f"网络异常: {exc}"
            logger.debug("LS bootstrap [%s/%s]: %s", attempt, BOOTSTRAP_MAX_RETRIES, last_err)

    raise LabelStudioClientError(
        f"Label Studio 服务未就绪，自动登录失败（已重试 {BOOTSTRAP_MAX_RETRIES} 次）: {last_err}"
    )


def _auto_bootstrap_token(
    base_url: str,
    *,
    username: Optional[str] = None,
    password: Optional[str] = None,
    timeout_sec: Optional[float] = None,
    force_refresh: bool = False,
) -> str:
    """
    兼容旧调用方：返回可用于 Bearer refresh 的 PAT，或 legacy Token。
    优先 PAT；legacy 仅在组织启用时回退。
    """
    base = base_url.rstrip("/")
    timeout = timeout_sec if timeout_sec is not None else DEFAULT_TIMEOUT_SEC
    session = _auto_bootstrap_session(
        base,
        username=username,
        password=password,
        timeout_sec=timeout,
        force_refresh=force_refresh,
    )
    pat = _create_pat_with_session(session, base, timeout)
    if pat:
        return pat
    legacy = _fetch_legacy_token_with_session(session, base, timeout)
    if legacy:
        return legacy
    raise LabelStudioClientError(
        "Label Studio 登录成功但无法获取 API 凭据（legacy 已禁用且 PAT 已存在）。"
        "请在 LS Account & Settings 撤销旧 PAT 后重试，或设置 LABEL_STUDIO_API_KEY 为完整 PAT。"
    )


class LabelStudioClient:
    """
    Label Studio HTTP API 轻量封装。

    环境变量：
    - LABEL_STUDIO_URL：服务基址，默认 http://127.0.0.1:8082
    - LABEL_STUDIO_API_KEY：可选 PAT（JWT）或 legacy Token；未设置时自动 session 登录
    - LABEL_STUDIO_USERNAME / LABEL_STUDIO_PASSWORD：零配置登录凭据
    - LABEL_STUDIO_TIMEOUT_SEC：请求超时秒数，默认 30
    """

    def __init__(
        self,
        base_url: Optional[str] = None,
        api_key: Optional[str] = None,
        timeout_sec: Optional[float] = None,
        *,
        auto_bootstrap: bool = True,
    ):
        self.base_url = (base_url or DEFAULT_LS_BASE_URL).rstrip("/")
        self.timeout_sec = timeout_sec if timeout_sec is not None else DEFAULT_TIMEOUT_SEC
        self._auto_bootstrap_enabled = auto_bootstrap
        self._session = requests.Session()
        self._session.headers.setdefault("Content-Type", "application/json")

        env_key = (api_key if api_key is not None else DEFAULT_LS_API_KEY).strip()
        self.api_key = env_key
        self._explicit_api_key = bool(env_key)
        self._auth_mode: Optional[str] = None
        self._pat_refresh_token = ""
        self._bearer_access_token = ""

        if self.api_key:
            self._init_explicit_api_key(self.api_key)

    def _init_explicit_api_key(self, token: str) -> None:
        if _is_jwt_pat(token):
            self._auth_mode = _AUTH_BEARER
            self._pat_refresh_token = token
            return
        self._auth_mode = _AUTH_TOKEN
        self._apply_token_header(token)

    @staticmethod
    def _set_token_header(session: requests.Session, token: str) -> None:
        session.headers["Authorization"] = f"Token {token}"

    @staticmethod
    def _set_bearer_header(session: requests.Session, access: str) -> None:
        session.headers["Authorization"] = f"Bearer {access}"

    def _apply_token_header(self, token: str) -> None:
        self._set_token_header(self._session, token)

    def _ensure_bearer_access_token(self, *, force_refresh: bool = False) -> None:
        if self._bearer_access_token and not force_refresh:
            self._set_bearer_header(self._session, self._bearer_access_token)
            return
        self._bearer_access_token = _refresh_bearer_access(
            self._pat_refresh_token,
            self.base_url,
            self.timeout_sec,
        )
        self._set_bearer_header(self._session, self._bearer_access_token)

    def _ensure_authenticated(self, *, force_refresh: bool = False) -> None:
        if self._explicit_api_key and self._auth_mode:
            if self._auth_mode == _AUTH_BEARER:
                self._ensure_bearer_access_token(force_refresh=force_refresh)
            elif self._auth_mode == _AUTH_TOKEN:
                self._apply_token_header(self.api_key)
            return

        if not self._auto_bootstrap_enabled:
            raise LabelStudioClientError("Label Studio 服务未就绪：缺少 API 凭据且未启用自动 bootstrap")

        session = _auto_bootstrap_session(
            self.base_url,
            timeout_sec=self.timeout_sec,
            force_refresh=force_refresh,
        )
        self._session = session
        self._auth_mode = _AUTH_SESSION

    # 兼容旧私有方法名
    def _ensure_api_key(self) -> None:
        self._ensure_authenticated()

    def _fallback_to_session_bootstrap(self) -> None:
        logger.info("Label Studio 显式 Token 不可用，回退 session bootstrap…")
        self._explicit_api_key = False
        self.api_key = ""
        self._auth_mode = None
        self._pat_refresh_token = ""
        self._bearer_access_token = ""
        self._session = requests.Session()
        self._session.headers.setdefault("Content-Type", "application/json")
        _BOOTSTRAP_SESSION_CACHE.pop(self.base_url, None)
        self._ensure_authenticated(force_refresh=True)

    def _request(
        self,
        method: str,
        path: str,
        *,
        json_body: Any = None,
        params: Optional[Dict[str, Any]] = None,
        _auth_retry: bool = False,
        _session_retry: bool = False,
    ) -> Any:
        self._ensure_authenticated()
        url = urljoin(f"{self.base_url}/", path.lstrip("/"))
        try:
            resp = self._session.request(
                method=method.upper(),
                url=url,
                json=json_body,
                params=params,
                timeout=self.timeout_sec,
            )
        except Timeout as exc:
            logger.warning("Label Studio 请求超时: %s %s", method, url)
            raise LabelStudioClientError(
                f"Label Studio 请求超时（>{self.timeout_sec}s）: {url}",
            ) from exc
        except RequestsConnectionError as exc:
            logger.warning("Label Studio 连接失败: %s %s — %s", method, url, exc)
            raise LabelStudioClientError(
                f"Label Studio 服务未就绪或无法连接（{self.base_url}）。",
            ) from exc
        except RequestException as exc:
            logger.warning("Label Studio 网络异常: %s %s — %s", method, url, exc)
            raise LabelStudioClientError(f"Label Studio 网络请求失败: {exc}") from exc

        if resp.status_code == 401 and not _auth_retry:
            body_preview: Any
            try:
                body_preview = resp.json()
            except ValueError:
                body_preview = resp.text

            if (
                self._auth_mode == _AUTH_BEARER
                and self._explicit_api_key
                and not _session_retry
            ):
                logger.info("Label Studio Bearer access 过期，刷新 PAT…")
                self._ensure_bearer_access_token(force_refresh=True)
                return self._request(
                    method,
                    path,
                    json_body=json_body,
                    params=params,
                    _auth_retry=True,
                )

            if (
                self._auth_mode == _AUTH_TOKEN
                and self._explicit_api_key
                and _legacy_token_disabled(body_preview)
                and self._auto_bootstrap_enabled
            ):
                self._fallback_to_session_bootstrap()
                return self._request(
                    method,
                    path,
                    json_body=json_body,
                    params=params,
                    _auth_retry=True,
                    _session_retry=True,
                )

            if (
                self._auto_bootstrap_enabled
                and not self._explicit_api_key
                and self._auth_mode == _AUTH_SESSION
            ):
                logger.info("Label Studio session 失效，尝试重新登录…")
                _BOOTSTRAP_SESSION_CACHE.pop(self.base_url, None)
                self._auth_mode = None
                self._ensure_authenticated(force_refresh=True)
                return self._request(
                    method,
                    path,
                    json_body=json_body,
                    params=params,
                    _auth_retry=True,
                )

        if resp.status_code >= 400:
            body: Any
            try:
                body = resp.json()
            except ValueError:
                body = resp.text
            detail = body if isinstance(body, str) else str(body)
            logger.warning(
                "Label Studio HTTP %s: %s %s — %s",
                resp.status_code,
                method,
                url,
                detail[:500],
            )
            raise LabelStudioClientError(
                f"Label Studio API 错误 HTTP {resp.status_code}: {detail}",
                status_code=resp.status_code,
                response_body=body,
            )

        if resp.status_code == 204 or not resp.content:
            return {}
        try:
            return resp.json()
        except ValueError:
            return {"raw": resp.text}

    def health_check(self) -> Dict[str, Any]:
        """探测 Label Studio 是否可达（/health 或根路径）。"""
        for path in ("/health", "/api/health", "/"):
            try:
                data = self._request("GET", path)
                return {"status": "ok", "path": path, "detail": data}
            except LabelStudioClientError:
                continue
        raise LabelStudioClientError(f"Label Studio 健康检查失败: {self.base_url}")

    def create_project(
        self,
        title: str,
        *,
        label_config: Optional[str] = None,
        description: str = "",
    ) -> Dict[str, Any]:
        """
        创建标注项目。

        Returns:
            含 id、title 等字段的项目字典。
        """
        title = str(title or "").strip()
        if not title:
            raise LabelStudioClientError("项目 title 不能为空")

        payload: Dict[str, Any] = {
            "title": title,
            "label_config": label_config or DEFAULT_IMAGE_LABEL_CONFIG,
        }
        if description:
            payload["description"] = description

        result = self._request("POST", "/api/projects/", json_body=payload)
        project_id = result.get("id")
        if project_id is None:
            raise LabelStudioClientError("创建项目成功但未返回 id", response_body=result)
        logger.info("Label Studio 项目已创建: id=%s title=%s", project_id, title)
        return result

    def import_task(
        self,
        project_id: int,
        tasks: Union[List[Dict[str, Any]], Dict[str, Any]],
        *,
        predictions: Optional[List[Dict[str, Any]]] = None,
    ) -> Dict[str, Any]:
        """
        向项目导入一条或多条任务（可含模型预测 pre-annotation）。

        Args:
            project_id: Label Studio 项目 ID
            tasks: 单条 dict 或任务列表；每项可为 {"data": {...}} 或扁平字段（如 {"image": "..."}）
            predictions: 可选，与 tasks 对齐的预测列表
        """
        if isinstance(tasks, dict):
            task_list: List[Dict[str, Any]] = [tasks]
        else:
            task_list = list(tasks or [])

        if not task_list:
            raise LabelStudioClientError("import_task 需要至少一条任务数据")

        normalized: List[Dict[str, Any]] = []
        for idx, item in enumerate(task_list):
            if not isinstance(item, dict):
                raise LabelStudioClientError(f"任务 #{idx} 必须是 dict")
            if "data" in item:
                normalized.append(item)
            else:
                normalized.append({"data": item})

        if predictions:
            for idx, (task, pred) in enumerate(zip(normalized, predictions)):
                if pred:
                    task["predictions"] = pred if isinstance(pred, list) else [pred]

        path = f"/api/projects/{int(project_id)}/import"
        result = self._request("POST", path, json_body=normalized)
        logger.info(
            "Label Studio 任务已导入: project_id=%s count=%s sample=%s",
            project_id,
            len(normalized),
            normalized[0] if normalized else None,
        )
        return result

    def project_url(self, project_id: int) -> str:
        """生成可在浏览器中打开的标注项目 URL。"""
        return f"{self.base_url}/projects/{int(project_id)}/data"

    def export_annotations(self, project_id: int) -> Any:
        """
        导出项目标注（JSON）。Label Studio API: GET /api/projects/{id}/export
        """
        path = f"/api/projects/{int(project_id)}/export"
        return self._request("GET", path, params={"exportType": "JSON"})

    def list_tasks(self, project_id: int, *, page_size: int = 100) -> List[Dict[str, Any]]:
        """列出项目任务（含 annotations 字段）。"""
        path = f"/api/projects/{int(project_id)}/tasks"
        result = self._request("GET", path, params={"page_size": page_size})
        if isinstance(result, list):
            return result
        if isinstance(result, dict):
            return result.get("tasks") or result.get("results") or []
        return []
