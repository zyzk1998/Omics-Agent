/**
 * 系统设置弹窗 + MCP 插件中心（iOS 风格开关 + window 挂载）
 */
var MCP_PLUGIN_DEFS = [
    { key: 'web_search', label: '全网实时检索' },
    { key: 'authority_db', label: '权威数据库直连' },
    { key: 'private_data_engine', label: '私有化数据引擎' },
    { key: 'compute_scheduler', label: '计算资源智能调度' },
    { key: 'viz_oneclick', label: '科研成果一键可视化' },
];

function ensureMcpSwitchStyles() {
    if (document.getElementById('settings-mcp-ios-styles')) return;
    var st = document.createElement('style');
    st.id = 'settings-mcp-ios-styles';
    st.textContent =
        '.settings-mcp-section{margin-top:8px;padding-top:16px;border-top:1px solid #e5e7eb;}' +
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

function openSettingsModal() {
    if (typeof window.closeAllPopovers === 'function') {
        window.closeAllPopovers();
    }
    var el = document.getElementById('settings-modal');
    if (!el) return;
    var defaultModel = document.getElementById('setting-default-model');
    if (defaultModel) {
        defaultModel.value = localStorage.getItem('default_model') || 'deepseek-ai/DeepSeek-R1';
    }
    syncMcpSwitchesFromStorage();
    el.style.display = 'flex';
    el.setAttribute('aria-hidden', 'false');
}

function initSettingsPanelMcpSwitches() {
    var form = document.querySelector('#settings-modal .settings-form');
    if (!form) return;
    ensureMcpSwitchStyles();
    if (form.getAttribute('data-mcp-version') === '3') return;
    form.querySelectorAll('[data-mcp-row], #settings-mcp-root').forEach(function (n) {
        n.remove();
    });
    form.removeAttribute('data-mcp-init');
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
    block.innerHTML =
        '<p class="settings-mcp-title">MCP 插件中心</p><div class="settings-mcp-list" id="settings-mcp-list"></div>';
    var list = block.querySelector('#settings-mcp-list');
    MCP_PLUGIN_DEFS.forEach(function (def) {
        var row = document.createElement('div');
        row.className = 'settings-mcp-row';
        row.setAttribute('data-mcp-row', '1');
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
        list.appendChild(row);
    });
    if (localDataRow && localDataRow.parentNode === form) {
        form.insertBefore(block, localDataRow);
    } else {
        form.appendChild(block);
    }
    syncMcpSwitchesFromStorage();
    form.setAttribute('data-mcp-version', '3');
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

document.addEventListener('DOMContentLoaded', () => {
    initSettingsModalBindings();
    if (typeof initSettingsPanelMcpSwitches === 'function') {
        initSettingsPanelMcpSwitches();
    }
});
