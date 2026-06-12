/**
 * 浏览器 → Local Sidecar 落盘中继（远程 API 无法直连用户本机 8019）。
 * 小文件 inline Base64；大文件 staging URL + Sidecar 流式下载（避免浏览器 OOM）。
 */
(function () {
    'use strict';

    var _deployInflight = Object.create(null);
    var _sidecarCapsCache = null;
    var _sidecarCapsInflight = null;
    var _relayToastEl = null;

    function sidecarUrl(pathname) {
        if (typeof localSidecarUrl === 'function') return localSidecarUrl(pathname);
        return 'http://127.0.0.1:8019' + (pathname || '');
    }

    function authMerge() {
        if (typeof getAuthHeaders === 'function') return getAuthHeaders();
        return {};
    }

    function cloudBaseUrl() {
        if (typeof getCloudUploadBaseUrl === 'function') {
            return String(getCloudUploadBaseUrl()).replace(/\/+$/, '');
        }
        try {
            return String(window.location.origin || '').replace(/\/+$/, '');
        } catch (_e) {
            return 'http://127.0.0.1:8018';
        }
    }

    function resolveAbsoluteDownloadUrl(relativeOrAbsolute) {
        var u = String(relativeOrAbsolute || '').trim();
        if (/^https?:\/\//i.test(u)) return u;
        return new URL(u, cloudBaseUrl() + '/').href;
    }

    function formatSidecarApiError(body, httpStatus) {
        var status = httpStatus != null ? httpStatus : '';
        if (!body) return 'HTTP ' + status;
        var d = body.detail;
        if (typeof d === 'string' && d.trim()) return d.trim();
        if (Array.isArray(d)) {
            var parts = d.map(function (item) {
                if (typeof item === 'string') return item;
                if (item && typeof item === 'object') {
                    var loc = Array.isArray(item.loc) ? item.loc.filter(function (x) { return x !== 'body'; }).join('.') : '';
                    var msg = item.msg || item.message || '';
                    return (loc ? loc + ': ' : '') + msg;
                }
                return String(item);
            }).filter(Boolean);
            if (parts.length) return parts.join('; ');
        }
        if (d && typeof d === 'object') {
            try { return JSON.stringify(d); } catch (_e) { return String(d); }
        }
        if (typeof body.message === 'string' && body.message.trim()) return body.message.trim();
        try { return JSON.stringify(body); } catch (_e2) { return 'HTTP ' + status; }
    }

    function sidecarErrorFromResponse(out, context) {
        var msg = formatSidecarApiError(out.body, out.res && out.res.status);
        var status = out.res && out.res.status;
        if (status === 404) {
            return new Error('Sidecar 接口不存在 (404): ' + (context || '/api/workspace/silent_deploy') + '。' + sidecarUpgradeHint());
        }
        if (status === 422) {
            return new Error('Sidecar 请求参数被拒绝 (422): ' + msg + '。' + sidecarUpgradeHint());
        }
        if (status === 403) {
            return new Error(msg || 'Sidecar 拒绝访问 (403)，请确认通过桌面客户端打开且 Sidecar 监听 127.0.0.1:8019');
        }
        return new Error(msg || ('HTTP ' + status));
    }

    function sidecarUpgradeHint() {
        return '请完全退出并重启 Omics Agent 桌面客户端，打开 http://127.0.0.1:8019/health 确认 capabilities 含 silent_deploy_from_url 或 write_bytes。';
    }

    function formatThrownError(err) {
        if (!err) return '未知错误';
        if (typeof err === 'string') return err;
        if (err.message && typeof err.message === 'string' && err.message !== '[object Object]') return err.message;
        if (err.detail) return formatSidecarApiError(err, err.status);
        try { return JSON.stringify(err); } catch (_e) { return String(err); }
    }

    function ensureRelaySpinnerStyle() {
        if (document.getElementById('sidecar-relay-spinner-style')) return;
        var st = document.createElement('style');
        st.id = 'sidecar-relay-spinner-style';
        st.textContent = '@keyframes sidecar-relay-spin{to{transform:rotate(360deg)}}'
            + '.sidecar-relay-spinner{display:inline-block;width:14px;height:14px;margin-right:8px;'
            + 'border:2px solid rgba(255,255,255,.35);border-top-color:#fff;border-radius:50%;'
            + 'animation:sidecar-relay-spin .8s linear infinite;vertical-align:-2px}';
        document.head.appendChild(st);
    }

    function showRelayProgress(message) {
        ensureRelaySpinnerStyle();
        dismissRelayProgress(false);
        var el = document.createElement('div');
        el.id = 'sidecar-relay-progress-toast';
        el.className = 'omics-toast omics-toast--info';
        el.style.cssText = 'position:fixed;bottom:24px;left:50%;transform:translateX(-50%);background:#1d4ed8;color:#eff6ff;'
            + 'padding:12px 20px;border-radius:8px;font-size:14px;z-index:10003;box-shadow:0 4px 16px rgba(0,0,0,.35);'
            + 'max-width:min(92vw,520px);display:flex;align-items:center;line-height:1.45;';
        el.innerHTML = '<span class="sidecar-relay-spinner" aria-hidden="true"></span><span></span>';
        el.querySelector('span:last-child').textContent = message || '正在通过中继高速落盘到本地…';
        document.body.appendChild(el);
        _relayToastEl = el;
    }

    function dismissRelayProgress() {
        if (_relayToastEl && _relayToastEl.parentNode) {
            _relayToastEl.parentNode.removeChild(_relayToastEl);
        }
        _relayToastEl = null;
    }

    function finishRelayProgress(success, message) {
        dismissRelayProgress();
        if (typeof showToast === 'function') {
            showToast(message, success ? 'success' : 'danger');
        }
    }

    function hasDeployContent(pkg) {
        if (!pkg) return false;
        var files = pkg.files && pkg.files.length;
        var urls = pkg.download_items && pkg.download_items.length;
        return !!(files || urls);
    }

    function deployKey(pkg) {
        var sid = (pkg && pkg.session_id) ? String(pkg.session_id) : 'na';
        var host = (pkg && pkg.host_mount_path) ? String(pkg.host_mount_path) : '';
        var n = (pkg && pkg.files && pkg.files.length) ? pkg.files.length : 0;
        var u = (pkg && pkg.download_items && pkg.download_items.length) ? pkg.download_items.length : 0;
        var tok = pkg && pkg.staging_token ? String(pkg.staging_token) : '';
        return sid + '|' + host + '|' + n + '|' + u + '|' + tok;
    }

    function safeSessionFolderName(sessionTitle, sessionId, folderTimestamp, sessionFolder) {
        if (sessionFolder && String(sessionFolder).trim()) return String(sessionFolder).trim();
        var sid = String(sessionId || '').trim() || 'unknown_session';
        var raw = String(sessionTitle || '').trim();
        if (!raw) return sid;
        var cleaned = raw.replace(/[<>:"/\\|?*\x00-\x1f]/g, '_').replace(/\s+/g, ' ').replace(/^[\s.]+|[\s.]+$/g, '');
        if (!cleaned) return sid;
        var ts = String(folderTimestamp || '').trim();
        if (!ts) {
            var now = new Date();
            ts = String(now.getFullYear())
                + String(now.getMonth() + 1).padStart(2, '0')
                + String(now.getDate()).padStart(2, '0')
                + String(now.getHours()).padStart(2, '0')
                + String(now.getMinutes()).padStart(2, '0');
        }
        var suffix = '-' + ts;
        var maxBase = Math.max(1, 80 - suffix.length);
        if (cleaned.length > maxBase) cleaned = cleaned.slice(0, maxBase).replace(/[\s.]+$/g, '');
        return cleaned ? cleaned + suffix : sid;
    }

    function resolveDeploySessionFolder(deployPackage) {
        var pkg = deployPackage || {};
        if (pkg.session_folder && String(pkg.session_folder).trim()) return String(pkg.session_folder).trim();
        return safeSessionFolderName(pkg.session_title, pkg.session_id, pkg.folder_timestamp, null);
    }

    function joinHostPath(hostMount, relParts) {
        var host = String(hostMount || '').replace(/[/\\]+$/, '');
        var useBackslash = host.indexOf('\\') >= 0 && host.indexOf('/') < 0;
        var sep = useBackslash ? '\\' : '/';
        return [host].concat(relParts || []).join(sep);
    }

    function buildResultFilePath(hostMount, sessionTitle, sessionId, destName, deployPackage) {
        var folder = resolveDeploySessionFolder(deployPackage || {
            session_title: sessionTitle,
            session_id: sessionId,
        });
        var rel = String(destName || '').replace(/\\/g, '/').replace(/^\/+/, '');
        if (!rel) return joinHostPath(hostMount, [folder, 'result']);
        return joinHostPath(hostMount, [folder, 'result'].concat(rel.split('/')));
    }

    function withSessionFolderPayload(deployPackage) {
        var payload = Object.assign({}, deployPackage || {});
        if (!payload.session_folder) payload.session_folder = resolveDeploySessionFolder(deployPackage);
        return payload;
    }

    function mergeDeployResult(detail, body, hostMount, sessionTitle, sessionId, deployPackage) {
        var merged = Object.assign({}, detail || {});
        var paths = body.changed_paths || body.copied_paths || merged.changed_paths || [];
        if (!paths.length && hostMount && body && body.status === 'success') {
            var pkg = deployPackage || (detail && (detail.client_deploy_package || detail.deploy_package)) || {};
            var specs = (pkg.files || []).concat(pkg.download_items || []);
            paths = specs.map(function (spec) {
                return buildResultFilePath(hostMount, sessionTitle, sessionId, spec.dest_name, deployPackage);
            });
        }
        merged.changed_paths = paths;
        merged.changed_files = paths.map(function (p) {
            return { path: p, status: 'added' };
        });
        if (body.mount_tree) merged.mount_tree = body.mount_tree;
        if (body.result_dir) merged.result_dir = body.result_dir;
        if (body.session_folder) merged.session_folder = body.session_folder;
        else if (deployPackage && deployPackage.session_folder) merged.session_folder = deployPackage.session_folder;
        delete merged.client_deploy_package;
        delete merged.deploy_package;
        return merged;
    }

    function fetchSidecarCapabilities() {
        if (_sidecarCapsCache) return Promise.resolve(_sidecarCapsCache);
        if (_sidecarCapsInflight) return _sidecarCapsInflight;
        _sidecarCapsInflight = fetch(sidecarUrl('/health'), { method: 'GET' })
            .then(function (res) { return res.json().catch(function () { return {}; }); })
            .then(function (body) {
                _sidecarCapsCache = Array.isArray(body && body.capabilities) ? body.capabilities : [];
                return _sidecarCapsCache;
            })
            .catch(function () {
                _sidecarCapsCache = [];
                return _sidecarCapsCache;
            })
            .finally(function () {
                _sidecarCapsInflight = null;
            });
        return _sidecarCapsInflight;
    }

    function postJson(path, payload) {
        return fetch(sidecarUrl(path), {
            method: 'POST',
            headers: Object.assign({}, authMerge(), { 'Content-Type': 'application/json' }),
            body: JSON.stringify(payload),
        }).then(function (res) {
            return res.json().catch(function () { return {}; }).then(function (body) {
                return { res: res, body: body };
            });
        });
    }

    function buildUrlDeployPayload(deployPackage) {
        var auth = authMerge();
        var items = (deployPackage.download_items || []).map(function (it) {
            return {
                dest_name: it.dest_name,
                download_url: resolveAbsoluteDownloadUrl(it.download_url),
                authorization: auth.Authorization || auth.authorization || null,
                x_guest_uuid: auth['X-Guest-UUID'] || null,
            };
        });
        return {
            session_id: deployPackage.session_id,
            session_title: deployPackage.session_title,
            session_folder: deployPackage.session_folder || resolveDeploySessionFolder(deployPackage),
            folder_timestamp: deployPackage.folder_timestamp,
            host_mount_path: deployPackage.host_mount_path,
            items: items,
            files: deployPackage.files || [],
        };
    }

    function deployViaUrlMode(deployPackage) {
        var payload = buildUrlDeployPayload(deployPackage);
        var label = (deployPackage.download_items || []).length > 1
            ? '正在经 Sidecar 流式下载 ' + deployPackage.download_items.length + ' 个大文件到本地…'
            : '正在经 Sidecar 流式下载归档到本地（浏览器不经手大文件）…';
        showRelayProgress(label);
        return postJson('/api/workspace/silent_deploy_from_url', payload).then(function (out) {
            if (!out.res.ok || !out.body || out.body.status !== 'success') {
                throw sidecarErrorFromResponse(out, '/api/workspace/silent_deploy_from_url');
            }
            return out.body;
        });
    }

    function deployViaWriteBytesFallback(deployPackage) {
        var host = deployPackage.host_mount_path;
        var sid = deployPackage.session_id;
        var title = deployPackage.session_title || sid;
        var files = deployPackage.files || [];
        var changed = [];
        var writePath = '/api/tools/write_bytes';
        var folder = resolveDeploySessionFolder(deployPackage);

        return fetchSidecarCapabilities().then(function (caps) {
            if (caps.indexOf('silent_deploy') < 0 && caps.indexOf('write_bytes') < 0) {
                throw new Error(
                    '本机 Sidecar 不支持 silent_deploy / write_bytes（capabilities: '
                    + (caps.length ? caps.join(', ') : '未知') + '）。' + sidecarUpgradeHint()
                );
            }
            writePath = caps.indexOf('write_bytes') >= 0 ? '/api/tools/write_bytes' : '/api/tools/write_file';
            var chain = Promise.resolve();
            files.forEach(function (spec) {
                chain = chain.then(function () {
                    var destPath = buildResultFilePath(host, title, sid, spec.dest_name, deployPackage);
                    return postJson(writePath, {
                        file_path: destPath,
                        content: '',
                        content_b64: spec.content_b64,
                    }).then(function (out) {
                        if (out.res.status === 404 && writePath === '/api/tools/write_bytes') {
                            writePath = '/api/tools/write_file';
                            return postJson(writePath, {
                                file_path: destPath,
                                content: '',
                                content_b64: spec.content_b64,
                            });
                        }
                        if (!out.res.ok || !out.body || out.body.status !== 'success') {
                            throw sidecarErrorFromResponse(out, writePath);
                        }
                        changed.push(out.body.file_path || destPath);
                    });
                });
            });
            return chain.then(function () {
                return {
                    status: 'success',
                    changed_paths: changed,
                    copied_paths: changed,
                    session_folder: folder,
                    result_dir: buildResultFilePath(host, title, sid, '', deployPackage),
                    destination_dir: buildResultFilePath(host, title, sid, '', deployPackage),
                };
            });
        });
    }

    function deployViaInlineMode(deployPackage) {
        showRelayProgress('正在通过中继落盘小文件到本地…');
        return postJson('/api/workspace/silent_deploy', withSessionFolderPayload(deployPackage)).then(function (out) {
            if (out.res.status === 404) {
                _sidecarCapsCache = null;
                return deployViaWriteBytesFallback(deployPackage);
            }
            if (!out.res.ok || !out.body || out.body.status !== 'success') {
                if (out.res.status === 404 || /not found/i.test(formatSidecarApiError(out.body, out.res.status))) {
                    return deployViaWriteBytesFallback(deployPackage);
                }
                throw sidecarErrorFromResponse(out, '/api/workspace/silent_deploy');
            }
            return out.body;
        }).catch(function (err) {
            if (/not found|404/i.test(String(err && err.message || err))) {
                return deployViaWriteBytesFallback(deployPackage);
            }
            throw err;
        });
    }

    function usesUrlTransport(deployPackage) {
        var t = String(deployPackage.transport || '').toLowerCase();
        if (t === 'url' || t === 'mixed') return true;
        return !!(deployPackage.download_items && deployPackage.download_items.length);
    }

    function deployPackageViaLocalSidecar(deployPackage) {
        if (!hasDeployContent(deployPackage)) {
            return Promise.reject(new Error('无效的落盘包'));
        }
        var run;
        if (usesUrlTransport(deployPackage)) {
            run = deployViaUrlMode(deployPackage);
        } else {
            run = deployViaInlineMode(deployPackage);
        }
        return run.then(function (body) {
            finishRelayProgress(true, '本地落盘完成');
            return body;
        }).catch(function (err) {
            dismissRelayProgress();
            if (typeof showToast === 'function') {
                showToast('本机 Sidecar 落盘失败: ' + formatThrownError(err), 'danger');
            }
            throw err;
        });
    }

    function maybeDeployClientPackageViaSidecar(detail) {
        var pkg = (detail && (detail.client_deploy_package || detail.deploy_package)) || null;
        if (!hasDeployContent(pkg)) {
            return Promise.resolve(detail || null);
        }
        var key = deployKey(pkg);
        if (_deployInflight[key]) return _deployInflight[key];

        _deployInflight[key] = deployPackageViaLocalSidecar(pkg)
            .then(function (body) {
                return mergeDeployResult(
                    detail,
                    body,
                    pkg.host_mount_path,
                    pkg.session_title,
                    pkg.session_id,
                    pkg
                );
            })
            .catch(function (err) {
                var low = formatThrownError(err).toLowerCase();
                if (/not found|404/.test(low) && !/422/.test(low)) {
                    throw new Error('本机 Sidecar 版本过旧，缺少落盘接口。' + sidecarUpgradeHint());
                }
                throw err;
            })
            .finally(function () {
                delete _deployInflight[key];
            });

        return _deployInflight[key];
    }

    function dispatchSessionFilesChanged(detail) {
        if (!detail || typeof window.dispatchEvent !== 'function') return;
        window.dispatchEvent(new CustomEvent('omics:session-files-changed', { detail: detail }));
    }

    function handleFilesNotificationWithClientDeploy(detail) {
        return maybeDeployClientPackageViaSidecar(detail).then(function (merged) {
            if (merged) dispatchSessionFilesChanged(merged);
            return merged;
        }).catch(function (err) {
            if (typeof showToast === 'function') {
                showToast('本机 Sidecar 落盘失败: ' + formatThrownError(err), 'danger');
            }
            throw err;
        });
    }

    window.deployPackageViaLocalSidecar = deployPackageViaLocalSidecar;
    window.maybeDeployClientPackageViaSidecar = maybeDeployClientPackageViaSidecar;
    window.handleFilesNotificationWithClientDeploy = handleFilesNotificationWithClientDeploy;
    window.fetchSidecarCapabilities = fetchSidecarCapabilities;
    window.formatSidecarApiError = formatSidecarApiError;
    window.formatThrownError = formatThrownError;
})();
