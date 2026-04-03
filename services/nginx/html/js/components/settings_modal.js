/**
 * 系统设置弹窗 + MCP 插件中心（iOS 风格开关 + window 挂载）
 */
/** HPC-Open 与后端约定：enabled_mcps 仍使用 compute_scheduler */
var HPC_OPEN_MCP_KEY = 'compute_scheduler';

var MCP_PLUGIN_DEFS = [
    { key: 'web_search', label: '全网实时检索' },
    { key: 'authority_db', label: '权威数据库直连' },
    { key: 'private_data_engine', label: '私有化数据引擎' },
    { key: 'viz_oneclick', label: '科研成果一键可视化' },
];

function ensureMcpSwitchStyles() {
    if (document.getElementById('settings-mcp-ios-styles')) return;
    var st = document.createElement('style');
    st.id = 'settings-mcp-ios-styles';
    st.textContent =
        '.settings-mcp-section{margin-top:8px;padding-top:16px;border-top:1px solid var(--input-border,#e5e7eb);}' +
        '.settings-mcp-title{font-size:13px;font-weight:600;color:var(--text-color);opacity:.85;margin:0 0 12px;letter-spacing:.02em;}' +
        '.settings-mcp-list{display:flex;flex-direction:column;gap:14px;}' +
        '.settings-mcp-row{display:flex;align-items:center;justify-content:space-between;gap:16px;min-height:36px;}' +
        '.settings-mcp-label{font-size:14px;color:var(--text-color);flex:1;line-height:1.4;}' +
        '.mcp-ios-wrap{flex-shrink:0;cursor:pointer;display:inline-flex;align-items:center;}' +
        '.mcp-ios-wrap input.mcp-ios-input{position:absolute;opacity:0;width:0;height:0;pointer-events:none;}' +
        '.mcp-ios-track{position:relative;display:inline-block;width:51px;height:31px;background:#e5e5ea;border-radius:16px;transition:background .28s cubic-bezier(.4,0,.2,1);vertical-align:middle;box-shadow:inset 0 0 0 1px rgba(0,0,0,.06);}' +
        '.mcp-ios-thumb{position:absolute;top:2px;left:2px;width:27px;height:27px;border-radius:50%;background:#fff;box-shadow:0 2px 6px rgba(0,0,0,.18);transition:transform .28s cubic-bezier(.4,0,.2,1);}' +
        '.mcp-ios-wrap input.mcp-ios-input:checked + .mcp-ios-track{background:#34c759;box-shadow:inset 0 0 0 1px rgba(0,0,0,.04);}' +
        '.mcp-ios-wrap input.mcp-ios-input:checked + .mcp-ios-track .mcp-ios-thumb{transform:translateX(20px);}' +
        '.mcp-ios-wrap input.mcp-ios-input:focus-visible + .mcp-ios-track{outline:2px solid #4d6bfe;outline-offset:2px;}';
    document.head.appendChild(st);
}

function readEnabledMcpKeys() {
    return getEnabledMcpsFromStorage();
}

function writeEnabledMcpKeys(keys) {
    var uniq = [];
    var seen = {};
    (keys || []).forEach(function (k) {
        var s = String(k || '').trim();
        if (s && !seen[s]) {
            seen[s] = true;
            uniq.push(s);
        }
    });
    localStorage.setItem('enabled_mcps', JSON.stringify(uniq));
}

function syncMcpSwitchesFromStorage() {
    var set = {};
    readEnabledMcpKeys().forEach(function (k) {
        set[k] = true;
    });
    document.querySelectorAll('#settings-modal .mcp-ios-input').forEach(function (inp) {
        var key = inp.getAttribute('data-mcp-key');
        if (key) inp.checked = !!set[key];
    });
}

function closeSettingsModal() {
    var el = document.getElementById('settings-modal');
    if (el) {
        el.style.display = 'none';
        el.setAttribute('aria-hidden', 'true');
    }
}

function applyDefaultModelFromStorage() {
    var stored = localStorage.getItem('default_model') || 'deepseek-ai/DeepSeek-R1';
    var sel = document.getElementById('modelSelect');
    if (!sel) return;
    var opts = [].slice.call(sel.options || []);
    for (var i = 0; i < opts.length; i++) {
        if (opts[i].value === stored) {
            sel.value = stored;
            return;
        }
    }
    sel.value = opts[0] ? opts[0].value : stored;
}

function getEnabledMcpsFromStorage() {
    try {
        var raw = localStorage.getItem('enabled_mcps');
        if (!raw) return [];
        var parsed = JSON.parse(raw);
        if (!Array.isArray(parsed)) return [];
        return parsed
            .map(function (x) {
                return String(x).trim();
            })
            .filter(Boolean);
    } catch (e) {
        return [];
    }
}

/** 聊天输入区 HPC 标识：需 HPC-Open（compute_scheduler）开启且后端 MCP 已连接 */
function updateHpcContextBadge(optionalStatus) {
    var badge = document.getElementById('hpc-context-badge');
    if (!badge) return;
    var schedOn = getEnabledMcpsFromStorage().indexOf(HPC_OPEN_MCP_KEY) >= 0;
    if (!schedOn) {
        badge.classList.add('hidden');
        return;
    }
    function apply(data) {
        if (data && data.connected) {
            badge.classList.remove('hidden');
        } else {
            badge.classList.add('hidden');
        }
    }
    if (optionalStatus && typeof optionalStatus.connected === 'boolean') {
        apply(optionalStatus);
        return;
    }
    fetch('/api/config/mcp/status', {
        headers: typeof globalThis.getAuthHeaders === 'function' ? globalThis.getAuthHeaders() : {},
    })
        .then(function (r) {
            return r.ok ? r.json() : {};
        })
        .then(function (data) {
            apply(data);
        })
        .catch(function () {
            badge.classList.add('hidden');
        });
}

function ensureSettingsMcpUiMounted() {
    var form = document.querySelector('#settings-modal .settings-form');
    if (!form) return;
    var verOk = form.getAttribute('data-mcp-version') === '7';
    var v2 = form.querySelector('#settings-mcp-root[data-mcp-layout="v2"]');
    var hasScheduler = document.getElementById('mcp-switch-compute_scheduler');
    var hasHpcRow = document.getElementById('settings-hpc-mcp-row');
    if (verOk && v2 && hasScheduler && hasHpcRow) return;
    form.removeAttribute('data-mcp-version');
    initSettingsPanelMcpSwitches();
}

function openSettingsModal() {
    if (typeof window.closeAllPopovers === 'function') {
        window.closeAllPopovers();
    }
    var el = document.getElementById('settings-modal');
    if (!el) return;
    initSettingsModalBindings();
    ensureSettingsMcpUiMounted();
    var defaultModel = document.getElementById('setting-default-model');
    if (defaultModel) {
        defaultModel.value = localStorage.getItem('default_model') || 'deepseek-ai/DeepSeek-R1';
    }
    syncMcpSwitchesFromStorage();
    updateHpcContextBadge();
    el.style.display = 'flex';
    el.setAttribute('aria-hidden', 'false');
}

/** 去掉旧版脚本插入的重复 #settings-mcp-root，仅保留 index 中带 data-mcp-layout="v2" 的节点 */
function cleanupDuplicateMcpRoots(form) {
    var all = [].slice.call(form.querySelectorAll('#settings-mcp-root'));
    var v2 = null;
    for (var i = 0; i < all.length; i++) {
        if (all[i].getAttribute('data-mcp-layout') === 'v2') {
            v2 = all[i];
            break;
        }
    }
    all.forEach(function (el) {
        if (el !== v2) {
            el.remove();
        }
    });
    return v2;
}

/** 无静态 MCP 区块时的兜底挂载（与 index 中 v2 结构一致） */
function mountMcpSectionLegacy(form) {
    var localDataRow = null;
    var rows = form.querySelectorAll('.settings-row');
    for (var r = 0; r < rows.length; r++) {
        if (rows[r].querySelector('#btn-clear-cache')) {
            localDataRow = rows[r];
            break;
        }
    }
    var block = document.createElement('div');
    block.id = 'settings-mcp-root';
    block.className = 'settings-mcp-section';
    block.setAttribute('data-mcp-layout', 'v2');
    block.innerHTML =
        '<p class="settings-mcp-title">MCP 插件中心</p>' +
        '<div id="settings-mcp-list" class="settings-mcp-list">' +
        '<div class="settings-mcp-row" id="settings-hpc-mcp-row" data-mcp-row="1">' +
        '<span class="settings-mcp-label" id="mcp-switch-compute_scheduler-lbl">HPC-Open (超算中心)</span>' +
        '<label class="mcp-ios-wrap" for="mcp-switch-compute_scheduler">' +
        '<input type="checkbox" class="mcp-ios-input" id="mcp-switch-compute_scheduler" data-mcp-key="compute_scheduler" aria-labelledby="mcp-switch-compute_scheduler-lbl" />' +
        '<span class="mcp-ios-track"><span class="mcp-ios-thumb"></span></span>' +
        '</label></div></div>';
    if (localDataRow && localDataRow.parentNode === form) {
        form.insertBefore(block, localDataRow);
    } else {
        form.appendChild(block);
    }
}

function appendPluginMcpRows(list) {
    var hpcRow = document.getElementById('settings-hpc-mcp-row');
    if (!list || !hpcRow) return;
    list.querySelectorAll('[data-mcp-plugin-row="1"]').forEach(function (n) {
        n.remove();
    });
    MCP_PLUGIN_DEFS.forEach(function (def) {
        var row = document.createElement('div');
        row.className = 'settings-mcp-row';
        row.setAttribute('data-mcp-row', '1');
        row.setAttribute('data-mcp-plugin-row', '1');
        var lid = 'mcp-switch-' + def.key;
        row.innerHTML =
            '<span class="settings-mcp-label" id="' +
            lid +
            '-lbl">' +
            def.label +
            '</span>' +
            '<label class="mcp-ios-wrap" for="' +
            lid +
            '">' +
            '<input type="checkbox" class="mcp-ios-input" id="' +
            lid +
            '" data-mcp-key="' +
            def.key +
            '" aria-labelledby="' +
            lid +
            '-lbl" />' +
            '<span class="mcp-ios-track"><span class="mcp-ios-thumb"></span></span>' +
            '</label>';
        var inp = row.querySelector('.mcp-ios-input');
        inp.addEventListener('change', function () {
            var keys = readEnabledMcpKeys().filter(function (x) {
                return x !== def.key;
            });
            if (inp.checked) keys.push(def.key);
            writeEnabledMcpKeys(keys);
        });
        list.insertBefore(row, hpcRow);
    });
}

function bindHpcMcpSwitch() {
    var hpcRow = document.getElementById('settings-hpc-mcp-row');
    if (!hpcRow) return;
    var hpcLid = 'mcp-switch-' + HPC_OPEN_MCP_KEY;
    var hpcInp = document.getElementById(hpcLid);
    if (hpcInp && hpcRow.getAttribute('data-hpc-input-bound') !== '1') {
        hpcRow.setAttribute('data-hpc-input-bound', '1');
        hpcInp.addEventListener('change', function () {
            var keys = readEnabledMcpKeys().filter(function (x) {
                return x !== HPC_OPEN_MCP_KEY;
            });
            if (hpcInp.checked) keys.push(HPC_OPEN_MCP_KEY);
            writeEnabledMcpKeys(keys);
            updateHpcContextBadge();
        });
    }
}

function initSettingsPanelMcpSwitches() {
    var form = document.querySelector('#settings-modal .settings-form');
    if (!form) return;
    ensureMcpSwitchStyles();
    if (form.getAttribute('data-mcp-version') === '7') return;
    cleanupDuplicateMcpRoots(form);
    var list = document.getElementById('settings-mcp-list');
    var hpcRow = document.getElementById('settings-hpc-mcp-row');
    if (!list || !hpcRow) {
        mountMcpSectionLegacy(form);
        list = document.getElementById('settings-mcp-list');
        hpcRow = document.getElementById('settings-hpc-mcp-row');
    }
    if (!list || !hpcRow) return;
    appendPluginMcpRows(list);
    bindHpcMcpSwitch();
    syncMcpSwitchesFromStorage();
    form.setAttribute('data-mcp-version', '7');
}

/** Badge 关闭：隐藏标识，并强制 HPC 开关 OFF + 触发 change 以写 localStorage */
function disableHpcOpenAndSyncUi() {
    var badge = document.getElementById('hpc-context-badge');
    if (badge) badge.classList.add('hidden');
    var inp = document.getElementById('mcp-switch-compute_scheduler');
    if (inp && inp.checked) {
        inp.checked = false;
        inp.dispatchEvent(new Event('change', { bubbles: true }));
    } else {
        var keys = readEnabledMcpKeys().filter(function (x) {
            return x !== HPC_OPEN_MCP_KEY;
        });
        writeEnabledMcpKeys(keys);
        if (inp) inp.checked = false;
        updateHpcContextBadge();
    }
}

function initHpcContextBadgeCloseButton() {
    var badge = document.getElementById('hpc-context-badge');
    if (!badge || badge.getAttribute('data-badge-close-bound') === '1') return;
    badge.setAttribute('data-badge-close-bound', '1');
    var btn = badge.querySelector('.badge-close-btn');
    if (!btn) return;
    btn.addEventListener('click', function (e) {
        e.preventDefault();
        e.stopPropagation();
        disableHpcOpenAndSyncUi();
    });
}

function initSettingsModalBindings() {
    var root = document.getElementById('settings-modal');
    if (!root || root.getAttribute('data-js-bound') === '1') return;
    root.setAttribute('data-js-bound', '1');
    var settingsBackdrop = document.getElementById('settings-modal-backdrop');
    var settingsClose = document.getElementById('settings-modal-close');
    if (settingsBackdrop) settingsBackdrop.onclick = closeSettingsModal;
    if (settingsClose) settingsClose.onclick = closeSettingsModal;
    var settingDefaultModel = document.getElementById('setting-default-model');
    if (settingDefaultModel) {
        settingDefaultModel.onchange = function () {
            localStorage.setItem('default_model', settingDefaultModel.value);
            applyDefaultModelFromStorage();
        };
    }
    var btnClearCache = document.getElementById('btn-clear-cache');
    var btnSettingsCancel = document.getElementById('settings-btn-cancel');
    var btnSettingsConfirm = document.getElementById('settings-btn-confirm');
    if (btnSettingsCancel) btnSettingsCancel.onclick = closeSettingsModal;
    if (btnSettingsConfirm) btnSettingsConfirm.onclick = closeSettingsModal;
    if (btnClearCache) {
        btnClearCache.onclick = function () {
            if (!confirm('将清除所有本地缓存（不包含登录状态），确定吗？')) return;
            window.workflowCache = {};
            closeSettingsModal();
            window.location.reload();
        };
    }
    localStorage.setItem('default_model', 'deepseek-ai/DeepSeek-R1');
    applyDefaultModelFromStorage();
}

window.openSettingsModal = openSettingsModal;
window.closeSettingsModal = closeSettingsModal;
window.getEnabledMcpsFromStorage = getEnabledMcpsFromStorage;
window.applyDefaultModelFromStorage = applyDefaultModelFromStorage;
window.updateHpcContextBadge = updateHpcContextBadge;
window.disableHpcOpenAndSyncUi = disableHpcOpenAndSyncUi;

function bootstrapSettingsModalModule() {
    initSettingsModalBindings();
    initSettingsPanelMcpSwitches();
    initHpcContextBadgeCloseButton();
    updateHpcContextBadge();
}

if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', bootstrapSettingsModalModule);
} else {
    bootstrapSettingsModalModule();
}
