/**
 * 工作台底部 · 三轨合一入库动态操作栏（Smart Mount Discovery + 静默入库）
 */
(function () {
    'use strict';

    var _lastIngestionResult = null;
    var _ingestInFlight = false;
    var LS_INGESTION_PATH_KEY = 'omics_ingestion_path';

    function authHeadersMerge() {
        if (typeof window.mergeJsonAuthHeaders === 'function') return window.mergeJsonAuthHeaders();
        return Object.assign(
            { 'Content-Type': 'application/json' },
            typeof getAuthHeaders === 'function' ? getAuthHeaders() : {}
        );
    }

    function getStoredIngestionPath() {
        try {
            return (localStorage.getItem(LS_INGESTION_PATH_KEY) || '').trim();
        } catch (_e) {
            return '';
        }
    }

    function setStoredIngestionPath(path) {
        try {
            localStorage.setItem(LS_INGESTION_PATH_KEY, String(path || '').trim());
        } catch (_e) { /* ignore quota */ }
    }

    function isConfigured(cfg) {
        if (getStoredIngestionPath()) return true;
        if (!cfg) return false;
        var mt = cfg.mount_type || 'local_volume';
        if (mt === 'local_volume') return !!((cfg.local_volume || {}).mount_path);
        if (mt === 'hpc_slurm') return !!((cfg.hpc_slurm || {}).host);
        if (mt === 'api_url') return !!((cfg.api_url || {}).endpoint);
        return false;
    }

    function getEls() {
        return {
            bar: document.getElementById('workspace-ingestion-action-bar'),
            hint: document.getElementById('workspace-ingestion-hint'),
            btn: document.getElementById('workspace-ingestion-primary-btn'),
            spinner: document.getElementById('workspace-ingestion-spinner'),
            err: document.getElementById('workspace-ingestion-error'),
        };
    }

    function setSpinner(visible) {
        var els = getEls();
        if (!els.spinner) return;
        if (visible) {
            els.spinner.hidden = false;
            els.spinner.setAttribute('aria-hidden', 'false');
        } else {
            els.spinner.hidden = true;
            els.spinner.setAttribute('aria-hidden', 'true');
        }
    }

    function setError(msg) {
        var els = getEls();
        if (!els.err) return;
        if (msg) {
            els.err.textContent = msg;
            els.err.hidden = false;
        } else {
            els.err.textContent = '';
            els.err.hidden = true;
        }
    }

    function resetBtnVisual(btn) {
        if (!btn) return;
        btn.disabled = false;
        btn.classList.remove(
            'workspace-ingestion-action-bar__btn--loading',
            'workspace-ingestion-action-bar__btn--success',
            'workspace-ingestion-action-bar__btn--error'
        );
    }

    function setLoadingState() {
        var els = getEls();
        if (!els.btn) return;
        resetBtnVisual(els.btn);
        els.btn.disabled = true;
        els.btn.classList.add('workspace-ingestion-action-bar__btn--loading');
        els.btn.textContent = '⏳ 正在打包并推送入库...';
        setSpinner(true);
        setError('');
    }

    function setSuccessState(res) {
        var els = getEls();
        if (!els.btn) return;
        resetBtnVisual(els.btn);
        els.btn.classList.add('workspace-ingestion-action-bar__btn--success');
        els.btn.textContent = '✅ 入库成功 (查看日志)';
        els.btn.disabled = false;
        setSpinner(false);
        setError('');
        els.btn.onclick = function () {
            showIngestionResultModal(res);
        };
    }

    function setFailureState(message) {
        var els = getEls();
        if (!els.btn) return;
        resetBtnVisual(els.btn);
        els.btn.classList.add('workspace-ingestion-action-bar__btn--error');
        els.btn.textContent = '❌ 入库失败 (重试)';
        els.btn.disabled = false;
        setSpinner(false);
        setError(message || '入库失败，请检查挂载配置与网络后重试');
    }

    function escapeHtml(s) {
        return String(s || '')
            .replace(/&/g, '&amp;')
            .replace(/</g, '&lt;')
            .replace(/>/g, '&gt;')
            .replace(/"/g, '&quot;');
    }

    function showIngestionResultModal(res) {
        var existing = document.getElementById('ingestion-result-modal');
        if (existing) existing.remove();

        var delivery = (res && res.delivery) || {};
        var manifest = (res && res.manifest) || {};
        var dest = delivery.destination || res.archive_path || '—';
        var files = Array.isArray(manifest.files) ? manifest.files : [];
        var filesPreview = files.slice(0, 12).map(function (f) {
            return '<li><code>' + escapeHtml(f) + '</code></li>';
        }).join('');
        if (files.length > 12) {
            filesPreview += '<li>… 共 ' + files.length + ' 项</li>';
        }

        var overlay = document.createElement('div');
        overlay.id = 'ingestion-result-modal';
        overlay.className = 'ingestion-result-modal-overlay';
        overlay.innerHTML =
            '<div class="ingestion-result-modal" role="dialog" aria-labelledby="ingestion-result-title">' +
            '<button type="button" class="ingestion-result-modal__close" aria-label="关闭">&times;</button>' +
            '<h3 id="ingestion-result-title">入库成功</h3>' +
            '<p class="ingestion-result-modal__dest"><strong>目标路径：</strong><code>' + escapeHtml(dest) + '</code></p>' +
            '<p><strong>挂载类型：</strong>' + escapeHtml(res.mount_type || delivery.strategy || '—') + '</p>' +
            '<p><strong>manifest.json 摘要</strong></p>' +
            '<ul class="ingestion-result-modal__files">' + (filesPreview || '<li>（无文件列表）</li>') + '</ul>' +
            '<button type="button" class="hitl-action-btn hitl-action-btn--primary ingestion-result-modal__ok">确定</button>' +
            '</div>';

        function closeModal() {
            overlay.remove();
        }
        overlay.querySelector('.ingestion-result-modal__close').addEventListener('click', closeModal);
        overlay.querySelector('.ingestion-result-modal__ok').addEventListener('click', closeModal);
        overlay.addEventListener('click', function (e) {
            if (e.target === overlay) closeModal();
        });
        document.body.appendChild(overlay);
    }

    function showMountConfirmModal(defaultPath, hint) {
        return new Promise(function (resolve, reject) {
            var existing = document.getElementById('ingestion-mount-confirm-modal');
            if (existing) existing.remove();

            var overlay = document.createElement('div');
            overlay.id = 'ingestion-mount-confirm-modal';
            overlay.className = 'ingestion-result-modal-overlay';
            overlay.innerHTML =
                '<div class="ingestion-result-modal" role="dialog" aria-labelledby="ingestion-mount-title">' +
                '<button type="button" class="ingestion-result-modal__close" aria-label="关闭">&times;</button>' +
                '<h3 id="ingestion-mount-title">确认入库挂载路径</h3>' +
                '<p class="text-muted small">系统已自动探测到 Docker 容器内可写目录。此为<strong>容器内路径</strong>，' +
                '对应 docker-compose 挂载卷；请勿填写宿主机盘符（如 D:\\project）。</p>' +
                '<p><strong>落盘路径（不可修改）：</strong></p>' +
                '<p><span class="mount-path-card" role="button" tabindex="0" data-path="' + escapeHtml(defaultPath).replace(/"/g, '&quot;') + '" title="在文件管理器中打开">' +
                '<span class="mount-path-card__icon" aria-hidden="true">📂</span>' +
                '<span class="mount-path-card__text">' + escapeHtml(defaultPath) + '</span></span></p>' +
                (hint ? ('<p class="small text-muted">' + escapeHtml(hint) + '</p>') : '') +
                '<div class="d-flex gap-2 justify-content-end mt-3">' +
                '<button type="button" class="btn btn-sm btn-outline-secondary ingestion-mount-cancel">取消</button>' +
                '<button type="button" class="btn btn-sm btn-primary ingestion-mount-confirm">确认并开始入库</button>' +
                '</div></div>';

            function closeModal(result) {
                overlay.remove();
                if (result) resolve(defaultPath);
                else reject(new Error('用户取消入库'));
            }
            overlay.querySelector('.ingestion-result-modal__close').addEventListener('click', function () {
                closeModal(false);
            });
            overlay.querySelector('.ingestion-mount-cancel').addEventListener('click', function () {
                closeModal(false);
            });
            overlay.querySelector('.ingestion-mount-confirm').addEventListener('click', function () {
                closeModal(true);
            });
            overlay.addEventListener('click', function (e) {
                if (e.target === overlay) closeModal(false);
            });
            document.body.appendChild(overlay);
        });
    }

    function discoverDefaultMountPath() {
        return fetch('/api/ingestion/discover-mount', { headers: authHeadersMerge() })
            .then(function (r) { return r.json(); })
            .then(function (body) {
                if (body.status !== 'success' || !body.default_path) {
                    throw new Error(body.message || '未能探测到可用挂载目录');
                }
                return body;
            });
    }

    function resolveMountPathForIngestion() {
        var stored = getStoredIngestionPath();
        if (stored) {
            return Promise.resolve({ path: stored, firstTime: false });
        }
        return discoverDefaultMountPath().then(function (body) {
            return showMountConfirmModal(body.default_path, body.hint).then(function (path) {
                setStoredIngestionPath(path);
                return { path: path, firstTime: true };
            });
        });
    }

    function triggerIngestion(skipHitl, mountPathOverride) {
        if (_ingestInFlight) return Promise.resolve();
        _ingestInFlight = true;
        setLoadingState();

        var mountPromise;
        if (mountPathOverride) {
            mountPromise = Promise.resolve({ path: mountPathOverride, firstTime: false });
        } else {
            mountPromise = resolveMountPathForIngestion();
        }

        return mountPromise
            .then(function (ctx) {
                return fetch('/api/ingestion/trigger', {
                    method: 'POST',
                    headers: authHeadersMerge(),
                    body: JSON.stringify({
                        session_id: typeof currentSessionId !== 'undefined' ? currentSessionId : null,
                        skip_hitl: !!skipHitl,
                        mount_path: ctx.path,
                        persist_mount_path: !!ctx.firstTime,
                    }),
                }).then(function (r) {
                    return r.json().then(function (body) {
                        if (!r.ok && body && !body.message) {
                            body.message = body.detail || ('HTTP ' + r.status);
                        }
                        return body;
                    });
                });
            })
            .then(function (res) {
                _lastIngestionResult = res;
                if (res.needs_settings && typeof window.openDatabaseSettingsPanel === 'function') {
                    window.openDatabaseSettingsPanel();
                    resetBtnVisual(getEls().btn);
                    if (getEls().btn) getEls().btn.textContent = '💾 一键入库 / 关联业务库';
                    setSpinner(false);
                    return res;
                }
                if (res.status === 'success') {
                    setSuccessState(res);
                    if (typeof showToast === 'function') {
                        showToast(res.message || '入库成功，数据已落盘至挂载目录', 'success');
                    }
                } else {
                    setFailureState(res.message || '入库失败');
                    if (typeof showToast === 'function') {
                        showToast(res.message || '入库失败', 'danger');
                    }
                    var els = getEls();
                    if (els.btn) {
                        els.btn.onclick = function () {
                            triggerIngestion(!!skipHitl);
                        };
                    }
                }
                return res;
            })
            .catch(function (e) {
                if (e && e.message === '用户取消入库') {
                    resetBtnVisual(getEls().btn);
                    if (getEls().btn) getEls().btn.textContent = '💾 一键入库';
                    setSpinner(false);
                    return;
                }
                setFailureState('入库请求失败: ' + (e.message || e));
                if (typeof showToast === 'function') showToast('入库请求失败: ' + e.message, 'danger');
                var els = getEls();
                if (els.btn) {
                    els.btn.onclick = function () {
                        triggerIngestion(!!skipHitl);
                    };
                }
            })
            .finally(function () {
                _ingestInFlight = false;
            });
    }

    function renderBar(cfg) {
        var els = getEls();
        if (!els.bar || !els.hint || !els.btn) return;

        var wp = document.getElementById('workspace-pane');
        var configured = isConfigured(cfg);
        var auto = cfg && cfg.is_auto_ingestion_enabled;
        var storedPath = getStoredIngestionPath();

        els.bar.classList.add('is-visible');
        if (wp && wp.classList.contains('workspace-active')) {
            els.bar.style.display = 'block';
        }

        if (_ingestInFlight) return;

        resetBtnVisual(els.btn);
        setSpinner(false);
        setError('');

        if (storedPath) {
            els.hint.textContent = '已记忆入库路径，点击后将静默落盘至 ' + storedPath;
        } else if (!configured) {
            els.hint.textContent = '首次入库将自动探测容器挂载路径，确认后终身免打扰。';
        } else if (auto) {
            els.hint.textContent = '已开启自动入库，可随时检查入库状态。';
        } else {
            els.hint.textContent = '业务库已关联，可手动触发三轨合一入库。';
        }

        if (auto) {
            els.btn.textContent = '✅ 检查数据是否成功入库';
        } else {
            els.btn.textContent = storedPath ? '💾 一键入库' : '💾 一键入库 / 关联业务库';
        }

        els.btn.onclick = function () {
            triggerIngestion(false);
        };
    }

    window.refreshWorkspaceIngestionBar = function () {
        var cfg = window.__userDatabaseMountConfig;
        if (cfg) {
            renderBar(cfg);
            return Promise.resolve(cfg);
        }
        if (typeof window.loadDatabaseSettings === 'function') {
            return window.loadDatabaseSettings().then(renderBar);
        }
        return fetch('/api/settings/database', { headers: authHeadersMerge() })
            .then(function (r) { return r.json(); })
            .then(function (res) {
                window.__userDatabaseMountConfig = res.config || null;
                renderBar(res.config);
            })
            .catch(function () { renderBar(null); });
    };

    window.triggerWorkspaceIngestion = triggerIngestion;
    window.getOmicsIngestionPath = getStoredIngestionPath;
    window.clearOmicsIngestionPath = function () {
        try { localStorage.removeItem(LS_INGESTION_PATH_KEY); } catch (_e) { /* ignore */ }
    };

    document.addEventListener('DOMContentLoaded', function () {
        window.refreshWorkspaceIngestionBar();
        var origOpen = window.openWorkspace;
        if (origOpen && !origOpen.__ingestWrapped) {
            window.openWorkspace = function () {
                var r = origOpen.apply(this, arguments);
                window.refreshWorkspaceIngestionBar();
                return r;
            };
            window.openWorkspace.__ingestWrapped = true;
        }
    });
})();
