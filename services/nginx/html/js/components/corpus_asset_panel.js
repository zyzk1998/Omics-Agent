/**
 * SFT 语料 JSON 工作台面板：Live SSE corpus_asset + 历史时光机恢复。
 */
(function () {
    'use strict';

    function escapeHtml(text) {
        return String(text == null ? '' : text)
            .replace(/&/g, '&amp;')
            .replace(/</g, '&lt;')
            .replace(/>/g, '&gt;')
            .replace(/"/g, '&quot;');
    }

    function resolveCorpusJsonString(payload) {
        if (!payload || typeof payload !== 'object') return '';
        if (typeof payload.corpus_sft_json === 'string' && payload.corpus_sft_json.trim()) {
            return payload.corpus_sft_json.trim();
        }
        if (typeof payload.corpus_json === 'string' && payload.corpus_json.trim()) {
            return payload.corpus_json.trim();
        }
        if (Array.isArray(payload.corpus_sft_records) && payload.corpus_sft_records.length) {
            try {
                return JSON.stringify(payload.corpus_sft_records, null, 2);
            } catch (_e) {
                return '';
            }
        }
        if (Array.isArray(payload.records) && payload.records.length) {
            try {
                return JSON.stringify(payload.records, null, 2);
            } catch (_e2) {
                return '';
            }
        }
        return '';
    }

    function resolveCorpusStandardLabel(payload) {
        var standard = String((payload && (payload.corpus_standard || payload.standard)) || '').toLowerCase();
        if (standard === 'gold') return 'Gold Standard · 金标准';
        if (standard === 'silver') return 'Silver · 银标准（机器生成）';
        return 'SFT Ready';
    }

    function buildCorpusContainerHtml(payload, jsonText) {
        var label = resolveCorpusStandardLabel(payload);
        var modality = String((payload && (payload.corpus_modality || payload.modality)) || '').toUpperCase();
        var count = payload && (payload.corpus_count != null ? payload.corpus_count : payload.count);
        var title = 'SFT Ready JSON (VLM / NLP)';
        if (modality) title += ' · ' + modality;
        if (count != null && count !== '') title += ' · ' + count + ' 条';
        title += ' · ' + label;
        return (
            '<div class="omics-corpus-container" style="margin-top: 20px; border: 1px solid #e2e8f0; border-radius: 8px; background: #1e1e1e;">' +
            '<div class="corpus-header" style="padding: 8px 12px; background: #2d2d2d; color: #a3adc2; font-size: 12px; display: flex; justify-content: space-between; align-items: center;">' +
            '<span>' + escapeHtml(title) + '</span>' +
            '<button type="button" class="copy-json-btn" title="复制 JSON" aria-label="复制 JSON" ' +
            'style="background:none;border:none;color:#a3adc2;cursor:pointer;padding:4px 6px;line-height:1;display:inline-flex;align-items:center;">' +
            '<i class="bi bi-clipboard" aria-hidden="true"></i></button>' +
            '</div>' +
            '<pre style="margin: 0; padding: 12px; max-height: 400px; overflow-y: auto; color: #d4d4d4; font-family: monospace; font-size: 13px; white-space: pre-wrap; word-break: break-word;">' +
            '<code class="language-json">' + escapeHtml(jsonText) + '</code></pre></div>'
        );
    }

    function copyTextWithFallback(text) {
        if (navigator.clipboard && navigator.clipboard.writeText) {
            return navigator.clipboard.writeText(text);
        }
        return new Promise(function (resolve, reject) {
            try {
                var ta = document.createElement('textarea');
                ta.value = text;
                ta.setAttribute('readonly', '');
                ta.style.position = 'fixed';
                ta.style.left = '-9999px';
                ta.style.top = '0';
                document.body.appendChild(ta);
                ta.select();
                ta.setSelectionRange(0, text.length);
                var ok = document.execCommand('copy');
                document.body.removeChild(ta);
                if (ok) resolve();
                else reject(new Error('execCommand copy failed'));
            } catch (err) {
                reject(err);
            }
        });
    }

    function showCopyJsonSuccess(btn) {
        var icon = btn.querySelector('i');
        if (!icon) return;
        var oldClass = icon.className;
        icon.className = 'bi bi-check2 text-success';
        btn.setAttribute('title', '已复制');
        window.setTimeout(function () {
            icon.className = oldClass;
            btn.setAttribute('title', '复制 JSON');
        }, 2000);
    }

    function installCorpusCopyDelegation() {
        if (window.__corpusCopyDelegationInstalled) return;
        window.__corpusCopyDelegationInstalled = true;
        document.addEventListener('click', function (e) {
            var btn = e.target && e.target.closest ? e.target.closest('.copy-json-btn') : null;
            if (!btn) return;
            e.preventDefault();
            e.stopPropagation();
            var container = btn.closest('.omics-corpus-container');
            if (!container) return;
            var codeBlock = container.querySelector('code.language-json');
            var jsonText = codeBlock ? (codeBlock.textContent || '') : '';
            if (!jsonText.trim()) {
                if (typeof window.showToast === 'function') window.showToast('无可复制的 JSON 内容', 'warning');
                return;
            }
            copyTextWithFallback(jsonText)
                .then(function () {
                    showCopyJsonSuccess(btn);
                    if (typeof window.showToast === 'function') window.showToast('JSON 语料已复制', 'success');
                })
                .catch(function (err) {
                    console.error('[CorpusPanel] 复制失败', err);
                    if (typeof window.showToast === 'function') window.showToast('复制失败，请手动选择 JSON 文本', 'warning');
                });
        });
    }

    installCorpusCopyDelegation();

    function resolveExpertMountEl() {
        var workspaceRoot = document.getElementById('active-workspace-content');
        if (typeof window.ensureOmicsReportThreeSlots === 'function' && workspaceRoot) {
            var slots = window.ensureOmicsReportThreeSlots(workspaceRoot);
            if (slots && slots.expert) return slots.expert;
        }
        return document.getElementById('omics-expert-report-slot');
    }

    function removeAllCorpusContainers(scopeEl) {
        var root = scopeEl || document;
        root.querySelectorAll('.omics-corpus-container').forEach(function (el) {
            el.remove();
        });
    }

    function renderCorpusAssetPanel(targetEl, payload) {
        var jsonText = resolveCorpusJsonString(payload);
        if (!jsonText) {
            removeAllCorpusContainers(targetEl || document);
            return null;
        }
        var mountEl = targetEl || resolveExpertMountEl();
        if (!mountEl) return null;
        removeAllCorpusContainers(mountEl);
        var html = buildCorpusContainerHtml(payload, jsonText);
        mountEl.insertAdjacentHTML('beforeend', html);
        var container = mountEl.querySelector('.omics-corpus-container');
        mountEl.style.display = '';
        mountEl.classList.remove('d-none', 'hidden');
        return container;
    }

    function mountCorpusPanelAfterExpertReport(payload) {
        var panel = renderCorpusAssetPanel(null, payload);
        if (window.__omicsFrozenExecutionSnapshot && payload) {
            var jsonText = resolveCorpusJsonString(payload);
            window.__omicsFrozenExecutionSnapshot.corpus_sft_json = jsonText;
            window.__omicsFrozenExecutionSnapshot.corpus_standard = payload.corpus_standard || payload.standard || '';
            window.__omicsFrozenExecutionSnapshot.corpus_modality = payload.corpus_modality || payload.modality || '';
            window.__omicsFrozenExecutionSnapshot.corpus_count = payload.corpus_count != null ? payload.corpus_count : payload.count;
            if (Array.isArray(payload.records)) {
                window.__omicsFrozenExecutionSnapshot.corpus_sft_records = payload.records;
            }
        }
        return panel;
    }

    function restoreCorpusAssetFromTimeMachine(slots, payload) {
        if (!payload) return;
        var ex = (payload.execution_snapshot && typeof payload.execution_snapshot === 'object')
            ? payload.execution_snapshot
            : {};
        var merged = {
            corpus_sft_json: payload.corpus_sft_json || ex.corpus_sft_json,
            corpus_sft_records: payload.corpus_sft_records || ex.corpus_sft_records,
            corpus_standard: payload.corpus_standard || ex.corpus_standard,
            corpus_modality: payload.corpus_modality || ex.corpus_modality,
            corpus_count: payload.corpus_count != null ? payload.corpus_count : ex.corpus_count
        };
        if (!resolveCorpusJsonString(merged)) return;
        if (typeof window.setRightPanelState === 'function') {
            window.setRightPanelState('WORKBENCH');
        } else if (typeof window.setWorkspacePaneExpanded === 'function') {
            window.setWorkspacePaneExpanded(true);
        }
        var target = (slots && slots.expert) ? slots.expert : null;
        renderCorpusAssetPanel(target, merged);
    }

    window.renderCorpusAssetPanel = renderCorpusAssetPanel;
    window.mountCorpusPanelAfterExpertReport = mountCorpusPanelAfterExpertReport;
    window.restoreCorpusAssetFromTimeMachine = restoreCorpusAssetFromTimeMachine;
    window.resolveCorpusJsonString = resolveCorpusJsonString;
})();
