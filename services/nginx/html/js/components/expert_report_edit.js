/**
 * 专家分析报告 · 沉浸式 Markdown 编辑（overlay；UI 文案不强调「全屏」）
 */
(function () {
    'use strict';

    var _editingMarkdown = null;
    var _sourceCardEl = null;

    function authHeadersMerge() {
        if (typeof window.mergeJsonAuthHeaders === 'function') return window.mergeJsonAuthHeaders();
        return Object.assign(
            { 'Content-Type': 'application/json' },
            typeof getAuthHeaders === 'function' ? getAuthHeaders() : {}
        );
    }

    function icon(name) {
        return '<i class="bi bi-' + name + '" aria-hidden="true"></i>';
    }

    function getReportVersionOpts(markdown) {
        var opts = {};
        if (window.__expertReportVersion === 'final' || /最终版/.test(String(markdown || ''))) {
            opts.hitl_final = true;
        }
        return opts;
    }

    function ensureFullscreenEditor() {
        var host = document.getElementById('expert-report-fullscreen-editor');
        if (host) return host;
        host = document.createElement('div');
        host.id = 'expert-report-fullscreen-editor';
        host.className = 'expert-report-fullscreen-editor';
        host.setAttribute('aria-hidden', 'true');
        host.innerHTML =
            '<div class="expert-report-fullscreen-editor__panel" role="dialog" aria-labelledby="expert-report-edit-title">' +
            '<div class="expert-report-fullscreen-editor__header">' +
            '<h2 id="expert-report-edit-title" class="expert-report-fullscreen-editor__title">' +
            icon('pencil-square') + ' 编辑专家报告</h2>' +
            '<button type="button" class="expert-report-fullscreen-editor__close" data-expert-back aria-label="返回">' +
            icon('arrow-left') + ' 返回</button>' +
            '</div>' +
            '<textarea class="expert-report-fullscreen-editor__textarea" aria-label="编辑专家报告 Markdown"></textarea>' +
            '<div class="expert-report-fullscreen-editor__actions">' +
            '<button type="button" class="hitl-action-btn hitl-action-btn--primary" data-expert-save>' +
            icon('save') + ' 保存并退出</button>' +
            '<button type="button" class="hitl-action-btn hitl-action-btn--secondary" data-expert-cancel>' +
            icon('x-lg') + ' 取消</button>' +
            '</div></div>';
        document.body.appendChild(host);
        host.querySelector('[data-expert-back]').addEventListener('click', function () {
            exitEditMode(false);
        });
        host.querySelector('[data-expert-cancel]').addEventListener('click', function () {
            exitEditMode(false);
        });
        host.querySelector('[data-expert-save]').addEventListener('click', function () {
            var ta = host.querySelector('.expert-report-fullscreen-editor__textarea');
            saveExpertReport(ta ? ta.value : _editingMarkdown, _sourceCardEl);
        });
        host.addEventListener('click', function (e) {
            if (e.target === host) exitEditMode(false);
        });
        return host;
    }

    function exitEditMode(reRender) {
        document.body.classList.remove('expert-report-fullscreen-active');
        var host = document.getElementById('expert-report-fullscreen-editor');
        if (host) {
            host.classList.remove('is-open');
            host.setAttribute('aria-hidden', 'true');
        }
        if (_sourceCardEl) {
            var hdrBtn = _sourceCardEl.querySelector('.expert-report-edit-btn');
            if (hdrBtn) hdrBtn.innerHTML = icon('pencil') + ' 修改';
        }
        if (reRender && typeof window.renderExpertReportIntoSlot === 'function') {
            var slot = document.getElementById('omics-expert-report-slot');
            window.renderExpertReportIntoSlot(slot, _editingMarkdown, getReportVersionOpts(_editingMarkdown));
        }
    }

    function enterEditMode(cardEl, markdown) {
        _sourceCardEl = cardEl;
        _editingMarkdown = markdown || '';
        var host = ensureFullscreenEditor();
        var ta = host.querySelector('.expert-report-fullscreen-editor__textarea');
        if (ta) ta.value = _editingMarkdown;
        host.classList.add('is-open');
        host.setAttribute('aria-hidden', 'false');
        document.body.classList.add('expert-report-fullscreen-active');
        if (ta) {
            ta.focus();
            ta.setSelectionRange(ta.value.length, ta.value.length);
        }
        var hdrBtn = cardEl && cardEl.querySelector('.expert-report-edit-btn');
        if (hdrBtn) hdrBtn.innerHTML = icon('pencil-square') + ' 编辑中';
    }

    function injectEditButton(cardEl, markdown) {
        if (!cardEl) return;
        var header = cardEl.querySelector('.card-header');
        if (!header || header.querySelector('.expert-report-edit-btn')) return;

        var btn = document.createElement('button');
        btn.type = 'button';
        btn.className = 'expert-report-edit-btn';
        btn.innerHTML = icon('pencil') + ' 修改';
        btn.setAttribute('aria-label', '修改专家报告');
        btn.addEventListener('click', function () {
            var body = cardEl.querySelector('.expert-report-md-body') || cardEl.querySelector('.report-expert-body');
            var md = markdown || (body ? body.innerText : '');
            enterEditMode(cardEl, md);
        });
        header.appendChild(btn);
    }

    function saveExpertReport(markdown, cardEl) {
        var sid = typeof currentSessionId !== 'undefined' ? currentSessionId : null;
        if (!sid) {
            if (typeof showToast === 'function') showToast('无会话 ID，无法持久化', 'warning');
            return;
        }
        fetch('/api/sessions/' + encodeURIComponent(sid) + '/expert_report', {
            method: 'PUT',
            headers: authHeadersMerge(),
            body: JSON.stringify({ markdown: markdown }),
        })
            .then(function (r) { return r.json(); })
            .then(function (res) {
                if (res.status !== 'success') {
                    throw new Error(res.message || res.detail || '保存失败');
                }
                _editingMarkdown = markdown;
                exitEditMode(true);
                if (window.__omicsFrozenExecutionSnapshot) {
                    window.__omicsFrozenExecutionSnapshot.expert_report_markdown = markdown;
                }
                if (typeof showToast === 'function') showToast('专家知识已固化', 'success');
            })
            .catch(function (e) {
                if (typeof showToast === 'function') showToast('保存失败: ' + e.message, 'danger');
            });
    }

    function wrapRenderExpertReport() {
        var orig = window.renderExpertReportIntoSlot;
        if (!orig || orig.__editWrapped) return;
        window.renderExpertReportIntoSlot = function (slotEl, markdown, options) {
            orig(slotEl, markdown, options);
            if (!slotEl || !markdown) return;
            var card = slotEl.querySelector('.omics-expert-report-card');
            if (card) injectEditButton(card, markdown);
        };
        window.renderExpertReportIntoSlot.__editWrapped = true;
    }

    if (document.readyState === 'loading') {
        document.addEventListener('DOMContentLoaded', wrapRenderExpertReport);
    } else {
        wrapRenderExpertReport();
    }
    document.addEventListener('DOMContentLoaded', function () {
        setTimeout(wrapRenderExpertReport, 800);
    });
})();
