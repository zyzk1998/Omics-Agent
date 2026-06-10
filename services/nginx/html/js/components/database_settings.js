/**
 * 挂载数据库设置（设置弹窗内表单）
 */
(function () {
    'use strict';

    function authHeadersMerge() {
        if (typeof window.mergeJsonAuthHeaders === 'function') return window.mergeJsonAuthHeaders();
        return Object.assign(
            { 'Content-Type': 'application/json' },
            typeof getAuthHeaders === 'function' ? getAuthHeaders() : {}
        );
    }

    function collectFormPayload() {
        var mountType = (document.getElementById('db-mount-type') || {}).value || 'local_volume';
        return {
            mount_type: mountType,
            is_auto_ingestion_enabled: !!(document.getElementById('db-auto-ingestion') || {}).checked,
            local_volume: {
                mount_path: ((document.getElementById('db-local-path') || {}).value || '').trim(),
            },
            hpc_slurm: {
                host: ((document.getElementById('db-hpc-host') || {}).value || '').trim(),
                port: parseInt((document.getElementById('db-hpc-port') || {}).value || '22', 10) || 22,
                username: ((document.getElementById('db-hpc-user') || {}).value || '').trim(),
                base_path: ((document.getElementById('db-hpc-base') || {}).value || '').trim(),
            },
            api_url: {
                endpoint: ((document.getElementById('db-api-endpoint') || {}).value || '').trim(),
                token: ((document.getElementById('db-api-token') || {}).value || '').trim(),
            },
        };
    }

    function applyPayloadToForm(cfg) {
        if (!cfg || typeof cfg !== 'object') return;
        var mt = document.getElementById('db-mount-type');
        if (mt) mt.value = cfg.mount_type || 'local_volume';
        var auto = document.getElementById('db-auto-ingestion');
        if (auto) auto.checked = !!cfg.is_auto_ingestion_enabled;
        var lv = cfg.local_volume || {};
        var hpc = cfg.hpc_slurm || {};
        var api = cfg.api_url || {};
        var set = function (id, v) { var el = document.getElementById(id); if (el) el.value = v != null ? v : ''; };
        set('db-local-path', lv.mount_path);
        set('db-hpc-host', hpc.host);
        set('db-hpc-port', hpc.port);
        set('db-hpc-user', hpc.username);
        set('db-hpc-base', hpc.base_path);
        set('db-api-endpoint', api.endpoint);
        set('db-api-token', api.token);
        toggleFieldGroups(cfg.mount_type || 'local_volume');
    }

    function toggleFieldGroups(mountType) {
        document.querySelectorAll('.settings-db-field-group').forEach(function (g) {
            g.classList.toggle('is-active', g.getAttribute('data-mount') === mountType);
        });
    }

    function loadDatabaseSettings() {
        return fetch('/api/settings/database', { headers: authHeadersMerge() })
            .then(function (r) { return r.json(); })
            .then(function (res) {
                if (res.config) applyPayloadToForm(res.config);
                window.__userDatabaseMountConfig = res.config || null;
                var lvPath = res.config && res.config.local_volume && res.config.local_volume.mount_path;
                var finish = function () {
                    if (typeof window.refreshWorkspaceIngestionBar === 'function') {
                        window.refreshWorkspaceIngestionBar();
                    }
                    return res.config;
                };
                if (!lvPath) {
                    return fetch('/api/ingestion/discover-mount', { headers: authHeadersMerge() })
                        .then(function (r) { return r.json(); })
                        .then(function (disc) {
                            if (disc && disc.default_path) {
                                var el = document.getElementById('db-local-path');
                                if (el && !String(el.value || '').trim()) {
                                    el.value = disc.default_path;
                                    el.placeholder = disc.default_path;
                                }
                            }
                            return finish();
                        })
                        .catch(function () { return finish(); });
                }
                return finish();
            })
            .catch(function () { return null; });
    }

    function saveDatabaseSettings() {
        var payload = collectFormPayload();
        return fetch('/api/settings/database', {
            method: 'PUT',
            headers: authHeadersMerge(),
            body: JSON.stringify(payload),
        })
            .then(function (r) { return r.json(); })
            .then(function (res) {
                if (res.status !== 'success') throw new Error(res.message || '保存失败');
                window.__userDatabaseMountConfig = res.config || payload;
                if (typeof showToast === 'function') showToast('数据库挂载配置已保存', 'success');
                if (typeof window.refreshWorkspaceIngestionBar === 'function') {
                    window.refreshWorkspaceIngestionBar();
                }
                return res;
            });
    }

    function testDatabaseConnection() {
        var payload = collectFormPayload();
        return fetch('/api/settings/database/test-connection', {
            method: 'POST',
            headers: authHeadersMerge(),
            body: JSON.stringify(payload),
        })
            .then(function (r) { return r.json(); })
            .then(function (res) {
                if (typeof showToast === 'function') {
                    showToast(res.message || (res.status === 'success' ? '关联测试通过' : '测试失败'), res.status === 'success' ? 'success' : 'danger');
                }
                return res;
            });
    }

    function bindDatabaseSettingsUi() {
        var mt = document.getElementById('db-mount-type');
        if (mt && !mt.__dbBound) {
            mt.__dbBound = true;
            mt.addEventListener('change', function () {
                toggleFieldGroups(mt.value);
            });
        }
        var testBtn = document.getElementById('db-settings-test');
        if (testBtn && !testBtn.__dbBound) {
            testBtn.__dbBound = true;
            testBtn.addEventListener('click', function () {
                testDatabaseConnection().catch(function (e) {
                    if (typeof showToast === 'function') showToast(e.message, 'danger');
                });
            });
        }
        var saveBtn = document.getElementById('db-settings-save');
        if (saveBtn && !saveBtn.__dbBound) {
            saveBtn.__dbBound = true;
            saveBtn.addEventListener('click', function () {
                saveDatabaseSettings().catch(function (e) {
                    if (typeof showToast === 'function') showToast(e.message, 'danger');
                });
            });
        }
    }

    function wrapOpenSettingsModal() {
        var orig = window.openSettingsModal;
        if (!orig || orig.__dbWrapped) return;
        window.openSettingsModal = function () {
            var ret = orig.apply(this, arguments);
            loadDatabaseSettings();
            bindDatabaseSettingsUi();
            return ret;
        };
        window.openSettingsModal.__dbWrapped = true;
    }

    document.addEventListener('DOMContentLoaded', function () {
        bindDatabaseSettingsUi();
        wrapOpenSettingsModal();
        loadDatabaseSettings();
    });

    window.loadDatabaseSettings = loadDatabaseSettings;
    window.openDatabaseSettingsPanel = function () {
        if (typeof openSettingsModal === 'function') openSettingsModal();
    };
})();
