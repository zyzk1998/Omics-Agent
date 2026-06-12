/**
 * 工作台底部 · 三轨合一入库动态操作栏（Smart Mount Discovery + 静默入库）
 */
(function () {
    'use strict';

    var _lastIngestionResult = null;
    var _lastIngestionError = null;
    var _ingestInFlight = false;
    var LS_INGESTION_PATH_KEY = 'omics_ingestion_path';
    var LS_INGESTION_SESSION_STATE_PREFIX = 'omics_ingestion_state_';

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

    function getSessionIngestionStateKey(sessionId) {
        var sid = String(sessionId || '').trim();
        return sid ? (LS_INGESTION_SESSION_STATE_PREFIX + sid) : '';
    }

    function clearSessionIngestionStateCache(sessionId) {
        var key = getSessionIngestionStateKey(sessionId);
        if (!key) return;
        try { localStorage.removeItem(key); } catch (_e) { /* ignore */ }
    }

    function saveSessionIngestionState(sessionId, payload) {
        var key = getSessionIngestionStateKey(sessionId);
        if (!key) return;
        try {
            localStorage.setItem(key, JSON.stringify(payload || {}));
        } catch (_e) { /* ignore */ }
    }

    function loadSessionIngestionState(sessionId) {
        var key = getSessionIngestionStateKey(sessionId);
        if (!key) return null;
        try {
            var raw = localStorage.getItem(key);
            return raw ? JSON.parse(raw) : null;
        } catch (_e) {
            return null;
        }
    }

    function getCachedWorkspaceSelectionSafe() {
        if (typeof window.getCachedWorkspaceSelection === 'function') {
            return window.getCachedWorkspaceSelection();
        }
        try {
            var raw = localStorage.getItem('workspaceSelection');
            return raw ? JSON.parse(raw) : null;
        } catch (_e) {
            return null;
        }
    }

    function getLocalWorkspaceRootPath() {
        var sel = getCachedWorkspaceSelectionSafe();
        if (!sel) return '';
        return String(sel.workspace_path || '').trim();
    }

    function getExplicitLocalProjectMountPath() {
        var root = getLocalWorkspaceRootPath();
        if (root) return root;
        var sel = getCachedWorkspaceSelectionSafe();
        if (!sel) return '';
        return String(sel.result_path || sel.workspace_path || '').trim();
    }

    function isLocalProjectMounted() {
        return !!getLocalWorkspaceRootPath();
    }

    function getDatabaseMountPath() {
        var cfg = window.__userDatabaseMountConfig;
        if (!cfg) return '';
        return String((cfg.local_volume || {}).mount_path || '').trim();
    }

    function normalizeMountPathForCompare(path) {
        return String(path || '').trim().replace(/\\/g, '/').replace(/\/+$/, '').toLowerCase();
    }

    function renderIngestionPathLink(path, label) {
        var p = String(path || '').trim();
        if (!p) return '—';
        var text = label != null ? String(label) : p;
        return (
            '<a href="#" class="ingestion-path-link" data-ingestion-path="' +
            escapeHtml(p).replace(/"/g, '&quot;') +
            '" title="在文件管理器中打开">' + escapeHtml(text) + '</a>'
        );
    }

    function resolveLocalResultNavigatePath(checkBody) {
        if (!checkBody) return '';
        if (checkBody.result_dir) return checkBody.result_dir;
        var paths = checkBody.sft_corpus_paths || [];
        if (paths.length) {
            var p = String(paths[0]);
            var idx = p.lastIndexOf('/');
            return idx > 0 ? p.slice(0, idx) : p;
        }
        if (checkBody.session_root) return checkBody.session_root + '/result';
        return '';
    }

    function checkLocalSessionResultExists(mountPath) {
        var sid = typeof currentSessionId !== 'undefined' ? currentSessionId : null;
        if (!sid) return Promise.resolve({ has_local_result: false });
        var mp = String(mountPath || getExplicitLocalProjectMountPath() || getDatabaseMountPath() || getStoredIngestionPath() || '').trim();
        var url = '/api/ingestion/sessions/' + encodeURIComponent(sid) + '/local-result-check';
        if (mp) url += '?mount_path=' + encodeURIComponent(mp);
        return fetch(url, { headers: authHeadersMerge() })
            .then(function (r) { return r.json(); })
            .catch(function () { return { status: 'error', has_local_result: false }; });
    }

    function handlePrimaryIngestionClick() {
        var mountPath = getExplicitLocalProjectMountPath() || getDatabaseMountPath() || getStoredIngestionPath();
        return checkLocalSessionResultExists(mountPath).then(function (body) {
            if (body && body.status === 'success' && body.has_local_result) {
                var navPath = resolveLocalResultNavigatePath(body);
                if (typeof showToast === 'function') {
                    showToast('数据已安全落盘至本地', 'success');
                }
                if (navPath && typeof window.setRightPanelState === 'function') {
                    window.setRightPanelState('FILE_MANAGER', 'local', navPath);
                } else if (navPath && typeof window.navigateFileManagerToPath === 'function') {
                    window.navigateFileManagerToPath(navPath, { mode: 'local' });
                }
                if (body.files_notification && typeof window.dispatchEvent === 'function') {
                    window.dispatchEvent(new CustomEvent('omics:session-files-changed', { detail: body.files_notification }));
                }
                return body;
            }
            return triggerIngestion(false);
        });
    }

    function bindIngestionPathLinks(rootEl) {
        if (!rootEl) return;
        rootEl.querySelectorAll('.ingestion-path-link').forEach(function (a) {
            a.addEventListener('click', function (ev) {
                ev.preventDefault();
                ev.stopPropagation();
                var targetPath = a.getAttribute('data-ingestion-path') || '';
                if (!targetPath) return;
                if (typeof window.navigateFileManagerToPath === 'function') {
                    var mode = typeof window.inferPathNavigationMode === 'function'
                        ? window.inferPathNavigationMode(targetPath)
                        : undefined;
                    window.navigateFileManagerToPath(targetPath, { mode: mode });
                } else if (typeof showToast === 'function') {
                    showToast('文件管理器未就绪', 'warning');
                }
            });
        });
    }

    function resetIngestionUiState() {
        _lastIngestionResult = null;
        _lastIngestionError = null;
        var els = getEls();
        if (!els.btn) return;
        resetBtnVisual(els.btn);
        setSpinner(false);
        setError('');
        els.btn.textContent = '💾 一键入库';
        els.btn.onclick = function () {
            handlePrimaryIngestionClick();
        };
    }

    function resetIngestionForSessionChange(prevSessionId) {
        if (prevSessionId) clearSessionIngestionStateCache(prevSessionId);
        resetIngestionUiState();
        if (typeof window.refreshWorkspaceIngestionBar === 'function') {
            window.refreshWorkspaceIngestionBar();
        }
    }

    window.resetWorkspaceIngestionForSessionChange = resetIngestionForSessionChange;

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

    function setFailureState(message, errorPayload) {
        var els = getEls();
        if (!els.btn) return;
        resetBtnVisual(els.btn);
        els.btn.classList.add('workspace-ingestion-action-bar__btn--error');
        els.btn.textContent = '❌ 入库失败 (查看详情)';
        els.btn.disabled = false;
        setSpinner(false);
        setError(message || '入库失败，请检查挂载配置与网络后重试');
        _lastIngestionError = errorPayload || { message: message };
        els.btn.onclick = function () {
            showIngestionErrorModal(_lastIngestionError);
        };
    }

    function buildIngestionErrorText(payload, fallbackMessage) {
        if (!payload) return String(fallbackMessage || '入库失败');
        var parts = [];
        if (payload.message) parts.push('【错误摘要】\n' + payload.message);
        if (payload.error_type) parts.push('【异常类型】\n' + payload.error_type);
        var trace = payload.traceback || payload.stack_trace || payload.detail || payload.error_details;
        if (typeof trace === 'object') trace = JSON.stringify(trace, null, 2);
        if (trace) parts.push('【完整 Stack Trace】\n' + trace);
        if (payload.mount_type) parts.push('【挂载类型】\n' + payload.mount_type);
        if (payload.archive_path) parts.push('【归档路径】\n' + payload.archive_path);
        return parts.join('\n\n') || String(fallbackMessage || '入库失败');
    }

    function resetIngestionTargetAndRetry() {
        clearOmicsIngestionPath();
        var sid = typeof currentSessionId !== 'undefined' ? currentSessionId : null;
        if (sid) clearSessionIngestionStateCache(sid);
        resetIngestionUiState();
        return buildMountPathCandidates().then(function (candidates) {
            var defaultPath = candidates.length ? candidates[0].path : '';
            return showMountConfirmModal(defaultPath, '请重新选择入库落盘目录', candidates);
        }).then(function (path) {
            setStoredIngestionPath(path);
            return triggerIngestion(false, path);
        }).catch(function (e) {
            if (e && e.message === '用户取消入库') return;
            if (typeof showToast === 'function') {
                showToast('未能重新指定入库路径', 'warning');
            }
        });
    }

    window.resetIngestionTargetAndRetry = resetIngestionTargetAndRetry;

    function buildMountPathCandidates() {
        var list = [];
        var seen = {};
        function add(path, label) {
            var p = String(path || '').trim();
            if (!p) return;
            var key = normalizeMountPathForCompare(p);
            if (seen[key]) return;
            seen[key] = true;
            list.push({ path: p, label: label || p });
        }
        add(getExplicitLocalProjectMountPath(), '本地项目挂载');
        add(getDatabaseMountPath(), '数据库绑定路径');
        add(getStoredIngestionPath(), '已记忆的入库路径');
        return discoverDefaultMountPath().then(function (body) {
            if (body.default_path) add(body.default_path, '容器默认路径');
            (body.paths || []).forEach(function (item) {
                if (item && item.path) add(item.path, item.label || '探测路径');
            });
            return list;
        }).catch(function () {
            return list;
        });
    }

    function bindResetIngestionTargetButton(rootEl) {
        if (!rootEl) return;
        var btn = rootEl.querySelector('.btn-reset-ingestion-target');
        if (!btn) return;
        btn.addEventListener('click', function (ev) {
            ev.preventDefault();
            var modal = rootEl.closest('.ingestion-result-modal-overlay');
            if (modal) modal.remove();
            resetIngestionTargetAndRetry();
        });
    }

    function showIngestionErrorModal(errPayload) {
        var existing = document.getElementById('ingestion-error-modal');
        if (existing) existing.remove();

        var payload = errPayload || _lastIngestionError || {};
        var fullText = buildIngestionErrorText(payload, payload.message || '服务器内部致命错误');
        var traceHtml = escapeHtml(fullText);
        var archivePath = payload.archive_path || '';
        var archiveLink = archivePath
            ? ('<p class="ingestion-result-modal__dest"><strong>归档路径：</strong>' + renderIngestionPathLink(archivePath) + '</p>')
            : '';

        var overlay = document.createElement('div');
        overlay.id = 'ingestion-error-modal';
        overlay.className = 'ingestion-result-modal-overlay ingestion-error-modal-overlay';
        overlay.innerHTML =
            '<div class="ingestion-result-modal ingestion-error-modal" role="alertdialog" aria-labelledby="ingestion-error-title">' +
            '<button type="button" class="ingestion-result-modal__close" aria-label="关闭">&times;</button>' +
            '<h3 id="ingestion-error-title" class="ingestion-error-modal__title">入库失败 · 详细日志</h3>' +
            '<p class="ingestion-error-modal__hint text-muted small">请将下方完整报错信息复制后反馈给运维或开发同事。</p>' +
            archiveLink +
            '<pre class="ingestion-error-modal__trace"><code>' + traceHtml + '</code></pre>' +
            '<div class="ingestion-error-modal__actions">' +
            '<button type="button" class="hitl-action-btn hitl-action-btn--primary ingestion-error-modal__copy">📋 一键复制报错信息</button>' +
            '<button type="button" class="hitl-action-btn ingestion-error-modal__retry">重试入库</button>' +
            '<button type="button" class="hitl-action-btn btn-reset-ingestion-target">🔄 重新指定数据库/目录</button>' +
            '</div></div>';

        function closeModal() {
            overlay.remove();
        }

        overlay.querySelector('.ingestion-result-modal__close').addEventListener('click', closeModal);
        overlay.querySelector('.ingestion-error-modal__retry').addEventListener('click', function () {
            closeModal();
            triggerIngestion(false);
        });
        overlay.querySelector('.ingestion-error-modal__copy').addEventListener('click', function () {
            var btn = overlay.querySelector('.ingestion-error-modal__copy');
            function done(ok) {
                if (typeof showToast === 'function') {
                    showToast(ok ? '报错信息已复制到剪贴板' : '复制失败，请手动选中日志复制', ok ? 'success' : 'warning');
                }
            }
            if (navigator.clipboard && navigator.clipboard.writeText) {
                navigator.clipboard.writeText(fullText).then(function () { done(true); }).catch(function () { done(false); });
            } else {
                try {
                    var ta = document.createElement('textarea');
                    ta.value = fullText;
                    document.body.appendChild(ta);
                    ta.select();
                    document.execCommand('copy');
                    ta.remove();
                    done(true);
                } catch (_e) {
                    done(false);
                }
            }
        });
        overlay.addEventListener('click', function (e) {
            if (e.target === overlay) closeModal();
        });
        bindIngestionPathLinks(overlay);
        bindResetIngestionTargetButton(overlay);
        document.body.appendChild(overlay);
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
            '<p class="ingestion-result-modal__dest"><strong>目标路径：</strong>' + renderIngestionPathLink(dest) + '</p>' +
            '<p><strong>挂载类型：</strong>' + escapeHtml(res.mount_type || delivery.strategy || '—') + '</p>' +
            '<p><strong>manifest.json 摘要</strong></p>' +
            '<ul class="ingestion-result-modal__files">' + (filesPreview || '<li>（无文件列表）</li>') + '</ul>' +
            '<div class="ingestion-result-modal__footer-actions">' +
            '<button type="button" class="hitl-action-btn hitl-action-btn--primary ingestion-result-modal__ok">确定</button>' +
            '<button type="button" class="hitl-action-btn btn-reset-ingestion-target">🔄 重新指定数据库/目录</button>' +
            '</div></div>';

        function closeModal() {
            overlay.remove();
        }
        overlay.querySelector('.ingestion-result-modal__close').addEventListener('click', closeModal);
        overlay.querySelector('.ingestion-result-modal__ok').addEventListener('click', closeModal);
        overlay.addEventListener('click', function (e) {
            if (e.target === overlay) closeModal();
        });
        bindIngestionPathLinks(overlay);
        bindResetIngestionTargetButton(overlay);
        document.body.appendChild(overlay);
    }

    function showMountConfirmModal(defaultPath, hint, candidates) {
        return new Promise(function (resolve, reject) {
            var existing = document.getElementById('ingestion-mount-confirm-modal');
            if (existing) existing.remove();

            var list = Array.isArray(candidates) && candidates.length
                ? candidates
                : [{ path: defaultPath, label: '推荐路径' }];
            var optionsHtml = list.map(function (item, idx) {
                var p = String(item.path || '').trim();
                var lbl = String(item.label || p || ('候选 ' + (idx + 1)));
                var checked = idx === 0 ? ' checked' : '';
                return (
                    '<label class="ingestion-mount-option">' +
                    '<input type="radio" name="ingestion-mount-choice" value="' + escapeHtml(p).replace(/"/g, '&quot;') + '"' + checked + '>' +
                    '<span class="ingestion-mount-option__label">' + escapeHtml(lbl) + '</span>' +
                    '<span class="ingestion-mount-option__path">' + renderIngestionPathLink(p, p) + '</span>' +
                    '</label>'
                );
            }).join('');

            var overlay = document.createElement('div');
            overlay.id = 'ingestion-mount-confirm-modal';
            overlay.className = 'ingestion-result-modal-overlay';
            overlay.innerHTML =
                '<div class="ingestion-result-modal" role="dialog" aria-labelledby="ingestion-mount-title">' +
                '<button type="button" class="ingestion-result-modal__close" aria-label="关闭">&times;</button>' +
                '<h3 id="ingestion-mount-title">确认入库挂载路径</h3>' +
                '<p class="text-muted small">优先级：<strong>本地项目挂载</strong> → 数据库绑定 → 已记忆路径 → 容器默认路径。</p>' +
                (hint ? ('<p class="small text-muted">' + escapeHtml(hint) + '</p>') : '') +
                '<div class="ingestion-mount-options">' + optionsHtml + '</div>' +
                '<div class="d-flex gap-2 justify-content-end mt-3">' +
                '<button type="button" class="btn btn-sm btn-outline-secondary ingestion-mount-cancel">取消</button>' +
                '<button type="button" class="btn btn-sm btn-primary ingestion-mount-confirm">确认并开始入库</button>' +
                '</div></div>';

            function closeModal(result, chosenPath) {
                overlay.remove();
                if (result) resolve(chosenPath || defaultPath);
                else reject(new Error('用户取消入库'));
            }
            overlay.querySelector('.ingestion-result-modal__close').addEventListener('click', function () {
                closeModal(false);
            });
            overlay.querySelector('.ingestion-mount-cancel').addEventListener('click', function () {
                closeModal(false);
            });
            overlay.querySelector('.ingestion-mount-confirm').addEventListener('click', function () {
                var picked = overlay.querySelector('input[name="ingestion-mount-choice"]:checked');
                var chosen = picked ? picked.value : defaultPath;
                closeModal(true, chosen);
            });
            overlay.addEventListener('click', function (e) {
                if (e.target === overlay) closeModal(false);
            });
            bindIngestionPathLinks(overlay);
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
        var localProject = getExplicitLocalProjectMountPath();
        var dbPath = getDatabaseMountPath();
        var stored = getStoredIngestionPath();

        if (localProject) {
            return Promise.resolve({ path: localProject, firstTime: false, source: 'local_project' });
        }
        if (dbPath) {
            return Promise.resolve({ path: dbPath, firstTime: false, source: 'database' });
        }
        if (stored) {
            return Promise.resolve({ path: stored, firstTime: false, source: 'stored' });
        }
        return buildMountPathCandidates().then(function (candidates) {
            var defaultPath = candidates.length ? candidates[0].path : '';
            return showMountConfirmModal(defaultPath, null, candidates).then(function (path) {
                setStoredIngestionPath(path);
                return { path: path, firstTime: true, source: 'confirmed' };
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
                var wsRoot = getLocalWorkspaceRootPath();
                var payload = {
                    session_id: typeof currentSessionId !== 'undefined' ? currentSessionId : null,
                    skip_hitl: !!skipHitl,
                    mount_path: ctx.path,
                    persist_mount_path: !!ctx.firstTime,
                };
                if (wsRoot) {
                    payload.workspace_path = wsRoot;
                    payload.client_track = 'local_sidecar';
                }
                return fetch('/api/ingestion/trigger', {
                    method: 'POST',
                    headers: authHeadersMerge(),
                    body: JSON.stringify(payload),
                }).then(function (r) {
                    return r.json().then(function (body) {
                        if (!r.ok) {
                            if (!body.message) {
                                body.message = body.detail || ('HTTP ' + r.status + ' 服务器内部致命错误');
                            }
                            if (!body.traceback && !body.stack_trace) {
                                body.stack_trace = 'HTTP ' + r.status + '\n' + (typeof body.detail === 'string' ? body.detail : JSON.stringify(body, null, 2));
                            }
                        }
                        return body;
                    });
                });
            })
            .then(function (res) {
                _lastIngestionResult = res;
                var sid = typeof currentSessionId !== 'undefined' ? currentSessionId : null;
                if (res.needs_settings && typeof window.openDatabaseSettingsPanel === 'function') {
                    window.openDatabaseSettingsPanel();
                    resetBtnVisual(getEls().btn);
                    if (getEls().btn) getEls().btn.textContent = '💾 一键入库 / 关联业务库';
                    setSpinner(false);
                    return res;
                }
                if (res.status === 'client_deploy') {
                    var deployFn = typeof window.deployPackageViaLocalSidecar === 'function'
                        ? window.deployPackageViaLocalSidecar
                        : null;
                    var pkg = res.deploy_package;
                    var hasPkg = pkg && ((pkg.files && pkg.files.length) || (pkg.download_items && pkg.download_items.length));
                    if (!deployFn || !hasPkg) {
                        setFailureState('缺少本机 Sidecar 落盘能力', res);
                        showIngestionErrorModal(res);
                        return res;
                    }
                    return deployFn(res.deploy_package).then(function (sideBody) {
                        var notify = res.files_notification || { session_id: sid, changed_paths: [] };
                        var paths = sideBody.changed_paths || sideBody.copied_paths || [];
                        notify.changed_paths = paths;
                        notify.changed_files = paths.map(function (p) { return { path: p, status: 'added' }; });
                        if (sideBody.mount_tree) notify.mount_tree = sideBody.mount_tree;
                        delete notify.client_deploy_package;
                        res.delivery = {
                            strategy: 'local_sidecar',
                            deploy_mode: 'client_relay',
                            destination_dir: sideBody.result_dir || sideBody.destination_dir,
                            changed_paths: paths,
                        };
                        res.status = 'success';
                        res.message = '入库成功，数据已落盘至本地挂载目录';
                        setSuccessState(res);
                        if (sid) {
                            saveSessionIngestionState(sid, {
                                status: 'success',
                                at: Date.now(),
                                destination: paths[0] || '',
                            });
                        }
                        if (typeof showToast === 'function') {
                            showToast(res.message, 'success');
                        }
                        window.dispatchEvent(new CustomEvent('omics:session-files-changed', { detail: notify }));
                        if (typeof window.refreshWorkspaceSessionFiles === 'function') {
                            window.refreshWorkspaceSessionFiles(true);
                        }
                        return res;
                    }).catch(function (sideErr) {
                        var errMsg = typeof window.formatThrownError === 'function'
                            ? window.formatThrownError(sideErr)
                            : (sideErr && sideErr.message) || String(sideErr);
                        var errRes = Object.assign({}, res, {
                            status: 'error',
                            message: '本机 Sidecar 落盘失败: ' + errMsg,
                            stack_trace: (sideErr && sideErr.stack) ? String(sideErr.stack) : errMsg,
                        });
                        setFailureState(errRes.message, errRes);
                        showIngestionErrorModal(errRes);
                        if (typeof showToast === 'function') showToast(errRes.message, 'danger');
                        return errRes;
                    });
                }
                if (res.status === 'success') {
                    setSuccessState(res);
                    if (sid) {
                        saveSessionIngestionState(sid, {
                            status: 'success',
                            at: Date.now(),
                            destination: ((res.delivery || {}).destination || res.archive_path || ''),
                        });
                    }
                    if (typeof showToast === 'function') {
                        showToast(res.message || '入库成功，数据已落盘至挂载目录', 'success');
                    }
                    if (res.files_notification) {
                        window.dispatchEvent(new CustomEvent('omics:session-files-changed', {
                            detail: res.files_notification,
                        }));
                    }
                    if (typeof window.refreshWorkspaceSessionFiles === 'function') {
                        window.refreshWorkspaceSessionFiles(true);
                    }
                } else {
                    setFailureState(res.message || '入库失败', res);
                    showIngestionErrorModal(res);
                    if (typeof showToast === 'function') {
                        showToast(res.message || '入库失败', 'danger');
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
                setFailureState('入库请求失败: ' + (e.message || e), {
                    message: '入库请求失败: ' + (e.message || e),
                    stack_trace: (e && e.stack) ? String(e.stack) : String(e),
                });
                showIngestionErrorModal(_lastIngestionError);
                if (typeof showToast === 'function') showToast('入库请求失败: ' + e.message, 'danger');
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

        var effectivePath = getExplicitLocalProjectMountPath() || getDatabaseMountPath() || storedPath;
        if (effectivePath) {
            els.hint.innerHTML = '已记忆入库路径，点击后将静默落盘至 ' + renderIngestionPathLink(effectivePath, effectivePath);
            bindIngestionPathLinks(els.bar);
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
            els.btn.textContent = storedPath || effectivePath ? '💾 一键入库' : '💾 一键入库 / 关联业务库';
        }

        els.btn.onclick = function () {
            handlePrimaryIngestionClick();
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
    window.handlePrimaryIngestionClick = handlePrimaryIngestionClick;
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
