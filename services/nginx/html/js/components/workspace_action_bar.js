/**
 * 工作台底部 · 三轨合一入库动态操作栏（Phase 3：进度 / 成功 / 失败态）
 */
(function () {
    'use strict';

    var _lastIngestionResult = null;
    var _ingestInFlight = false;

    function authHeadersMerge() {
        if (typeof window.mergeJsonAuthHeaders === 'function') return window.mergeJsonAuthHeaders();
        return Object.assign(
            { 'Content-Type': 'application/json' },
            typeof getAuthHeaders === 'function' ? getAuthHeaders() : {}
        );
    }

    function isConfigured(cfg) {
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

    function triggerIngestion(skipHitl) {
        if (_ingestInFlight) return Promise.resolve();
        _ingestInFlight = true;
        setLoadingState();

        return fetch('/api/ingestion/trigger', {
            method: 'POST',
            headers: authHeadersMerge(),
            body: JSON.stringify({
                session_id: typeof currentSessionId !== 'undefined' ? currentSessionId : null,
                skip_hitl: !!skipHitl,
            }),
        })
            .then(function (r) {
                return r.json().then(function (body) {
                    if (!r.ok && body && !body.message) {
                        body.message = body.detail || ('HTTP ' + r.status);
                    }
                    return body;
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
                        showToast(res.message || '入库成功', 'success');
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

        els.bar.classList.add('is-visible');
        if (wp && wp.classList.contains('workspace-active')) {
            els.bar.style.display = 'block';
        }

        if (_ingestInFlight) return;

        resetBtnVisual(els.btn);
        setSpinner(false);
        setError('');

        if (!configured) {
            els.hint.textContent = '尚未关联业务数据库，完成挂载后可一键入库。';
            els.btn.textContent = '💾 一键入库 / 关联业务库';
            els.btn.onclick = function () {
                if (typeof window.openDatabaseSettingsPanel === 'function') {
                    window.openDatabaseSettingsPanel();
                }
            };
            return;
        }

        if (auto) {
            els.hint.textContent = '已开启自动入库，可随时检查入库状态。';
            els.btn.textContent = '✅ 检查数据是否成功入库';
            els.btn.onclick = function () {
                triggerIngestion(false);
            };
            return;
        }

        els.hint.textContent = '业务库已关联，可手动触发三轨合一入库。';
        els.btn.textContent = '💾 一键入库';
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
