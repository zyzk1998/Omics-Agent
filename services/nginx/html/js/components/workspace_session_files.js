/**
 * 工作台顶部 · 会话输出文件总览（持久化壳层 + 原生右键菜单）
 */
(function () {
    'use strict';

    var _lastFetchSid = null;
    var _expanded = false;
    var _lastPanelData = null;
    var _contextMenuEl = null;

    var CTX_ACTION_CLASS = {
        open: 'ctx-open-file',
        explorer: 'ctx-open-dir',
        copy: 'ctx-copy-path',
        delete: 'ctx-delete-file',
    };

    function authHeadersMerge() {
        if (typeof window.mergeJsonAuthHeaders === 'function') return window.mergeJsonAuthHeaders();
        return Object.assign(
            { 'Content-Type': 'application/json' },
            typeof getAuthHeaders === 'function' ? getAuthHeaders() : {}
        );
    }

    function escapeHtml(s) {
        return String(s || '')
            .replace(/&/g, '&amp;')
            .replace(/</g, '&lt;')
            .replace(/>/g, '&gt;')
            .replace(/"/g, '&quot;');
    }

    function formatSize(bytes) {
        var n = Number(bytes) || 0;
        if (n < 1024) return n + ' B';
        if (n < 1024 * 1024) return (n / 1024).toFixed(2) + ' kb';
        return (n / (1024 * 1024)).toFixed(2) + ' MB';
    }

    function resolveSessionId(explicitSid) {
        if (explicitSid != null && String(explicitSid).trim()) return String(explicitSid).trim();
        return typeof currentSessionId !== 'undefined' && currentSessionId ? String(currentSessionId).trim() : '';
    }

    function sidecarUrl(pathname) {
        if (typeof localSidecarUrl === 'function') return localSidecarUrl(pathname);
        var port = typeof getLocalSidecarPort === 'function' ? getLocalSidecarPort() : '8019';
        var p = String(pathname || '');
        if (p.charAt(0) !== '/') p = '/' + p;
        return 'http://127.0.0.1:' + port + p;
    }

    function isLocalNativeFilePath(path) {
        var p = String(path || '').trim();
        if (!p) return false;
        if (typeof inferPathNavigationMode === 'function') {
            return inferPathNavigationMode(p) === 'local';
        }
        var lower = p.toLowerCase().replace(/\\/g, '/');
        if (lower.indexOf('/app/results/') === 0 || lower.indexOf('/app/uploads/') === 0) return false;
        if (lower.indexOf('/results/') === 0 || lower.indexOf('/uploads/') === 0) return false;
        return /^[a-z]:[\\/]/i.test(p) || lower.indexOf('/home/') === 0 || lower.indexOf('/users/') === 0;
    }

    function filterOutputFiles(data) {
        var list = (data && (data.outputs || data.results)) || [];
        return list.filter(function (f) {
            var zone = (f && (f.zone || f.kind)) || '';
            return zone !== 'upload';
        });
    }

    function ensureAccordionHost() {
        var scroll = document.querySelector('#workspace-pane .workspace-scroll');
        if (!scroll) return null;
        var host = document.getElementById('workspace-session-files-accordion');
        if (host && host.parentNode !== scroll) {
            host.parentNode.removeChild(host);
            host = null;
        }
        if (host) return host;
        host = document.createElement('div');
        host.id = 'workspace-session-files-accordion';
        host.className = 'workspace-session-files-accordion';
        host.setAttribute('aria-label', '会话中的文件');
        var printHeader = document.getElementById('workspace-print-header');
        if (printHeader && printHeader.parentNode === scroll) {
            scroll.insertBefore(host, printHeader);
        } else {
            scroll.insertBefore(host, scroll.firstChild);
        }
        return host;
    }

    function renderFileList(items) {
        if (!items || !items.length) {
            return '<p class="workspace-session-files-empty">暂无输出文件</p>';
        }
        return items.map(function (f) {
            var path = f.path || '';
            var name = f.name || path.split(/[/\\]/).pop() || '未命名';
            var tier = f.storage_tier === 'permanent' ? '永久' : '缓存';
            var nativeOk = isLocalNativeFilePath(path);
            return (
                '<a href="#" class="workspace-session-file-link"' +
                ' data-path="' + escapeHtml(path).replace(/"/g, '&quot;') + '"' +
                ' data-native="' + (nativeOk ? '1' : '0') + '"' +
                ' title="' + escapeHtml(path) + '">' +
                '<span class="workspace-session-file-link__icon" aria-hidden="true">📄</span>' +
                '<span class="workspace-session-file-link__meta">' +
                '<span class="workspace-session-file-link__name">' + escapeHtml(name) + '</span>' +
                '<span class="workspace-session-file-link__size">' + escapeHtml(formatSize(f.size_bytes)) + ' · ' + tier + '</span>' +
                '</span></a>'
            );
        }).join('');
    }

    function buildShellData(sessionTitle, outputs, errMsg) {
        return {
            session_title: sessionTitle || '当前会话',
            results: outputs || [],
            outputs: outputs || [],
            _error: errMsg || '',
        };
    }

    function updateBadge(count) {
        var badge = document.querySelector('#workspace-session-files-accordion .workspace-session-files-card__badge');
        if (badge) badge.textContent = String(count);
    }

    function hideContextMenu() {
        if (_contextMenuEl) {
            _contextMenuEl.remove();
            _contextMenuEl = null;
        }
    }

    function installGlobalContextMenuDismiss() {
        if (installGlobalContextMenuDismiss._bound) return;
        installGlobalContextMenuDismiss._bound = true;
        document.addEventListener('click', function (e) {
            if (!_contextMenuEl) return;
            if (e.target && _contextMenuEl.contains(e.target)) return;
            hideContextMenu();
        });
        document.addEventListener('contextmenu', function (e) {
            if (!_contextMenuEl) return;
            if (e.target && _contextMenuEl.contains(e.target)) return;
            hideContextMenu();
        });
        document.addEventListener('keydown', function (e) {
            if (e.key === 'Escape') hideContextMenu();
        });
        window.addEventListener('scroll', hideContextMenu, true);
        window.addEventListener('resize', hideContextMenu);
    }

    function installGlobalContextMenuActions() {
        if (installGlobalContextMenuActions._bound) return;
        installGlobalContextMenuActions._bound = true;

        document.addEventListener('mousedown', function (e) {
            if (!_contextMenuEl) return;
            var btn = e.target.closest('.workspace-session-file-context-menu__item[data-action]');
            if (!btn || !_contextMenuEl.contains(btn) || btn.disabled) return;
            e.preventDefault();
            e.stopPropagation();
            var action = btn.getAttribute('data-action') || '';
            var filePath = _contextMenuEl.dataset.targetPath || '';
            var nativeOk = _contextMenuEl.dataset.nativeOk === '1';
            hideContextMenu();
            handleContextAction(action, filePath, nativeOk);
        }, true);

        document.addEventListener('click', function (e) {
            var delegated = e.target.closest('.ctx-copy-path, .ctx-open-file, .ctx-open-dir, .ctx-delete-file');
            if (!delegated || !delegated.closest('#workspace-session-file-context-menu')) return;
        });
    }

    function showContextMenu(clientX, clientY, filePath, nativeOk) {
        hideContextMenu();
        var menu = document.createElement('div');
        menu.id = 'workspace-session-file-context-menu';
        menu.className = 'workspace-session-file-context-menu';
        menu.setAttribute('role', 'menu');
        menu.dataset.targetPath = String(filePath || '');
        menu.dataset.nativeOk = nativeOk ? '1' : '0';

        function addItem(label, action, disabled, title) {
            var btn = document.createElement('button');
            btn.type = 'button';
            btn.className = 'workspace-session-file-context-menu__item';
            if (CTX_ACTION_CLASS[action]) btn.classList.add(CTX_ACTION_CLASS[action]);
            btn.setAttribute('role', 'menuitem');
            btn.setAttribute('data-action', action);
            btn.dataset.targetPath = String(filePath || '');
            btn.textContent = label;
            if (disabled) {
                btn.disabled = true;
                btn.classList.add('is-disabled');
                if (title) btn.title = title;
            }
            menu.appendChild(btn);
        }

        addItem('打开文件', 'open', !nativeOk, '仅本地挂载文件支持原生打开');
        addItem('打开文件目录', 'explorer', !nativeOk, '仅本地挂载文件支持在资源管理器中定位');
        addItem('复制路径', 'copy', false, '');
        addItem('删除', 'delete', !nativeOk, '仅本地文件可删除');

        document.body.appendChild(menu);
        _contextMenuEl = menu;

        var mw = menu.offsetWidth || 200;
        var mh = menu.offsetHeight || 160;
        menu.style.left = Math.max(8, Math.min(clientX, window.innerWidth - mw - 8)) + 'px';
        menu.style.top = Math.max(8, Math.min(clientY, window.innerHeight - mh - 8)) + 'px';
    }

    function copyTextToClipboard(text, okMessage) {
        var msg = okMessage || '路径已复制';
        if (navigator.clipboard && navigator.clipboard.writeText) {
            return navigator.clipboard.writeText(String(text || '')).then(function () {
                if (typeof showToast === 'function') showToast(msg, 'success');
            });
        }
        return Promise.reject(new Error('clipboard unavailable'));
    }

    function sidecarFsPost(endpoint, filePath) {
        return fetch(sidecarUrl(endpoint), {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ path: filePath }),
        }).then(function (r) {
            return r.json().then(function (body) {
                if (!r.ok) {
                    var msg = (body && (body.detail || body.message)) || ('HTTP ' + r.status);
                    throw new Error(typeof msg === 'string' ? msg : JSON.stringify(msg));
                }
                return body;
            });
        });
    }

    function removeFileFromPanelData(filePath) {
        if (!_lastPanelData) return;
        var key = _lastPanelData.outputs ? 'outputs' : 'results';
        var list = _lastPanelData[key] || _lastPanelData.results || [];
        _lastPanelData[key] = list.filter(function (f) { return (f.path || '') !== filePath; });
        if (_lastPanelData.results) _lastPanelData.results = _lastPanelData[key];
    }

    function handleContextAction(action, filePath, nativeOk) {
        if (action === 'copy') {
            copyTextToClipboard(filePath, '文件路径已复制').catch(function () {
                if (nativeOk) {
                    sidecarFsPost('/api/fs/copy_file', filePath)
                        .then(function (body) {
                            var text = (body && body.clipboard_text) ? body.clipboard_text : filePath;
                            return copyTextToClipboard(text, '已复制到剪贴板');
                        })
                        .catch(function () {
                            copyTextToClipboard(filePath, '路径已复制');
                        });
                } else {
                    copyTextToClipboard(filePath, '路径已复制');
                }
            });
            return;
        }
        if (!nativeOk) {
            if (typeof showToast === 'function') showToast('仅本地挂载文件支持原生操作', 'info');
            return;
        }
        if (action === 'open') {
            sidecarFsPost('/api/fs/open_file', filePath)
                .then(function () {
                    if (typeof showToast === 'function') showToast('已请求打开文件', 'success');
                })
                .catch(function (err) {
                    if (typeof showToast === 'function') showToast('打开失败: ' + (err.message || err), 'danger');
                });
            return;
        }
        if (action === 'explorer') {
            sidecarFsPost('/api/fs/show_in_explorer', filePath)
                .then(function () {
                    if (typeof showToast === 'function') showToast('已在文件管理器中定位', 'success');
                })
                .catch(function (err) {
                    if (typeof showToast === 'function') showToast('定位失败: ' + (err.message || err), 'danger');
                });
            return;
        }
        if (action === 'delete') {
            if (!window.confirm('确定删除本地文件？\n' + filePath)) return;
            sidecarFsPost('/api/fs/delete_file', filePath)
                .then(function () {
                    var host = document.getElementById('workspace-session-files-accordion');
                    if (host) {
                        host.querySelectorAll('.workspace-session-file-link').forEach(function (link) {
                            if ((link.getAttribute('data-path') || '') === filePath) link.remove();
                        });
                    }
                    removeFileFromPanelData(filePath);
                    updateBadge(filterOutputFiles(_lastPanelData).length);
                    if (typeof showToast === 'function') showToast('文件已删除', 'success');
                })
                .catch(function (err) {
                    if (typeof showToast === 'function') showToast('删除失败: ' + (err.message || err), 'danger');
                });
        }
    }

    function bindPanelEvents(host, data) {
        var headerBtn = host.querySelector('.workspace-session-files-card__header');
        if (headerBtn) {
            headerBtn.addEventListener('click', function () {
                _expanded = !_expanded;
                renderAccordion(data);
            });
        }
        host.querySelectorAll('.workspace-session-file-link').forEach(function (link) {
            link.addEventListener('click', function (e) {
                e.preventDefault();
                var p = link.getAttribute('data-path') || '';
                if (!p) return;
                if (typeof window.navigateFileManagerToPath === 'function') {
                    var mode = typeof inferPathNavigationMode === 'function' ? inferPathNavigationMode(p) : undefined;
                    window.navigateFileManagerToPath(p, { mode: mode });
                } else if (typeof showToast === 'function') {
                    showToast('文件管理器未就绪', 'warning');
                }
            });
            link.addEventListener('contextmenu', function (e) {
                e.preventDefault();
                e.stopPropagation();
                showContextMenu(
                    e.clientX,
                    e.clientY,
                    link.getAttribute('data-path') || '',
                    link.getAttribute('data-native') === '1'
                );
            });
        });
    }

    function renderAccordion(data) {
        var host = ensureAccordionHost();
        if (!host) return;
        _lastPanelData = data || buildShellData('当前会话', []);
        var outputs = filterOutputFiles(_lastPanelData);
        var title = _lastPanelData.session_title || '当前会话';
        var errMsg = _lastPanelData._error || '';
        var expandedClass = _expanded ? ' is-expanded' : '';
        host.innerHTML =
            '<div class="workspace-session-files-card' + expandedClass + '">' +
            '<button type="button" class="workspace-session-files-card__header" aria-expanded="' + (_expanded ? 'true' : 'false') + '">' +
            '<span class="workspace-session-files-card__chevron" aria-hidden="true">▸</span>' +
            '<span class="workspace-session-files-card__title">会话中的文件</span>' +
            '<span class="workspace-session-files-card__badge">' + outputs.length + '</span>' +
            '</button>' +
            '<div class="workspace-session-files-card__body"' + (_expanded ? '' : ' hidden') + '>' +
            '<p class="workspace-session-files-card__subtitle text-muted small">' + escapeHtml(title) + '</p>' +
            (errMsg ? ('<p class="workspace-session-files-error text-danger small">' + escapeHtml(errMsg) + '</p>') : '') +
            '<section class="workspace-session-files-zone workspace-session-files-zone--outputs">' +
            '<h4 class="workspace-session-files-zone__title">输出文件</h4>' +
            renderFileList(outputs) +
            '</section>' +
            '</div></div>';

        bindPanelEvents(host, _lastPanelData);
    }

    function fetchAndRenderSessionFiles(explicitSid, force) {
        var sid = resolveSessionId(explicitSid);
        var host = ensureAccordionHost();
        if (!host) return Promise.resolve();

        if (!sid) {
            host.innerHTML = '';
            _lastFetchSid = null;
            return Promise.resolve();
        }

        if (String(sid) !== String(_lastFetchSid)) {
            _lastFetchSid = null;
        }

        var sessionTitle = '当前会话';
        if (window.cachedFullSessionsList) {
            for (var i = 0; i < window.cachedFullSessionsList.length; i++) {
                if (window.cachedFullSessionsList[i].id === sid) {
                    sessionTitle = window.cachedFullSessionsList[i].title || sessionTitle;
                    break;
                }
            }
        }

        if (!force && sid === _lastFetchSid && host.querySelector('.workspace-session-files-card')) {
            return Promise.resolve();
        }

        renderAccordion(buildShellData(sessionTitle, []));

        _lastFetchSid = sid;
        return fetch('/api/sessions/' + encodeURIComponent(sid) + '/files', { headers: authHeadersMerge() })
            .then(function (r) { return r.json(); })
            .then(function (body) {
                if (body.status === 'success') {
                    if (!body.session_title) body.session_title = sessionTitle;
                    renderAccordion(body);
                } else {
                    renderAccordion(buildShellData(sessionTitle, [], body.message || '文件清单加载失败'));
                }
            })
            .catch(function (err) {
                renderAccordion(buildShellData(sessionTitle, [], '文件清单加载失败: ' + (err.message || err)));
            });
    }

    window.refreshWorkspaceSessionFiles = function (force) {
        return fetchAndRenderSessionFiles(null, !!force);
    };
    window.fetchAndRenderSessionFiles = fetchAndRenderSessionFiles;

    window.addEventListener('omics:session-files-changed', function () {
        fetchAndRenderSessionFiles(null, true);
    });

    document.addEventListener('DOMContentLoaded', function () {
        installGlobalContextMenuDismiss();
        installGlobalContextMenuActions();
        var origOpen = window.openWorkspace;
        if (origOpen && !origOpen.__sessionFilesWrapped) {
            window.openWorkspace = function () {
                var r = origOpen.apply(this, arguments);
                fetchAndRenderSessionFiles(null, true);
                return r;
            };
            window.openWorkspace.__sessionFilesWrapped = true;
        }

        setTimeout(function () {
            var sid = resolveSessionId(null);
            if (sid) fetchAndRenderSessionFiles(sid, true);
        }, 800);
    });
})();
