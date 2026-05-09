/**
 * 帮助与反馈：Modal + 管理员列表数据源（POST /api/feedbacks）
 */
(function () {
    'use strict';

    function escapeHtml(text) {
        if (text == null || text === '') return '';
        return String(text)
            .replace(/&/g, '&amp;')
            .replace(/</g, '&lt;')
            .replace(/>/g, '&gt;')
            .replace(/"/g, '&quot;');
    }

    function getAuthHeadersMerged() {
        if (typeof window.getAuthHeaders === 'function') {
            const h = window.getAuthHeaders();
            return Object.assign({ 'Content-Type': 'application/json' }, h || {});
        }
        return { 'Content-Type': 'application/json' };
    }

    function guestHeadersIfAny() {
        const g = (typeof window.GUEST_UUID === 'string' && window.GUEST_UUID)
            ? window.GUEST_UUID
            : (typeof localStorage !== 'undefined' && localStorage.getItem('guest_uuid')) || '';
        if (!g) return {};
        return { 'X-Guest-UUID': g };
    }

    function ensureFeedbackModal() {
        if (document.getElementById('feedback-modal')) return;
        const holder = document.createElement('div');
        holder.innerHTML =
            '<div class="modal fade" id="feedback-modal" tabindex="-1" aria-labelledby="feedback-modal-label" aria-hidden="true">' +
            '  <div class="modal-dialog modal-dialog-centered">' +
            '    <div class="modal-content feedback-modal-content">' +
            '      <div class="modal-header">' +
            '        <h5 class="modal-title" id="feedback-modal-label">💡 帮助与反馈</h5>' +
            '        <button type="button" class="btn-close" data-bs-dismiss="modal" aria-label="关闭"></button>' +
            '      </div>' +
            '      <div class="modal-body">' +
            '        <label class="form-label small text-muted" for="feedback-type-select">分类</label>' +
            '        <select id="feedback-type-select" class="form-select form-select-sm mb-3">' +
            '          <option value="issue">遇到问题</option>' +
            '          <option value="suggestion">功能建议</option>' +
            '          <option value="error">报错 / 异常</option>' +
            '          <option value="other">其他</option>' +
            '        </select>' +
            '        <label class="form-label small text-muted" for="feedback-content-text">详细描述</label>' +
            '        <textarea id="feedback-content-text" class="form-control feedback-textarea" rows="5" placeholder="请描述现象、复现步骤或建议…"></textarea>' +
            '        <input type="hidden" id="feedback-error-context" value="" />' +
            '      </div>' +
            '      <div class="modal-footer">' +
            '        <button type="button" class="btn btn-secondary btn-sm" data-bs-dismiss="modal">取消</button>' +
            '        <button type="button" class="btn btn-primary btn-sm" id="feedback-submit-btn">提交反馈</button>' +
            '      </div>' +
            '    </div>' +
            '  </div>' +
            '</div>';
        document.body.appendChild(holder.firstElementChild);

        document.getElementById('feedback-submit-btn').addEventListener('click', function () {
            void submitFeedbackForm();
        });
    }

    async function submitFeedbackForm() {
        const typeEl = document.getElementById('feedback-type-select');
        const contentEl = document.getElementById('feedback-content-text');
        const ctxEl = document.getElementById('feedback-error-context');
        const type = typeEl ? typeEl.value : 'other';
        const content = contentEl ? String(contentEl.value || '').trim() : '';
        const error_context = ctxEl ? String(ctxEl.value || '').trim() : '';
        if (!content && !error_context) {
            if (typeof window.showToast === 'function') window.showToast('请填写描述或上下文', 'warning');
            return;
        }
        const body = {
            type: type,
            content: content || '(空)',
            error_context: error_context || null,
            timestamp: new Date().toISOString(),
        };
        try {
            const res = await fetch('/api/feedbacks', {
                method: 'POST',
                headers: Object.assign({}, getAuthHeadersMerged(), guestHeadersIfAny()),
                body: JSON.stringify(body),
            });
            const data = await res.json().catch(function () {
                return {};
            });
            if (!res.ok) {
                throw new Error(data.detail || data.message || 'HTTP ' + res.status);
            }
            if (typeof window.showToast === 'function') window.showToast('反馈已提交，感谢支持！', 'success');
            const modalEl = document.getElementById('feedback-modal');
            if (modalEl && typeof window.bootstrap !== 'undefined' && window.bootstrap.Modal) {
                const inst = window.bootstrap.Modal.getInstance(modalEl);
                if (inst) inst.hide();
            }
            if (contentEl) contentEl.value = '';
            if (ctxEl) ctxEl.value = '';
        } catch (e) {
            const msg = e && e.message ? e.message : String(e);
            if (typeof window.showToast === 'function') window.showToast('提交失败：' + msg, 'danger');
            else alert(msg);
        }
    }

    function openFeedbackModal(prefillErrorText) {
        ensureFeedbackModal();
        const ctxEl = document.getElementById('feedback-error-context');
        const contentEl = document.getElementById('feedback-content-text');
        const typeEl = document.getElementById('feedback-type-select');
        if (prefillErrorText) {
            if (ctxEl) ctxEl.value = prefillErrorText;
            if (typeEl) typeEl.value = 'error';
            if (contentEl && !String(contentEl.value || '').trim()) {
                contentEl.value = '【右键上报】以下为我复制的报错文本：\n\n' + prefillErrorText;
            }
        } else {
            if (ctxEl) ctxEl.value = '';
        }
        const modalEl = document.getElementById('feedback-modal');
        if (!modalEl || typeof window.bootstrap === 'undefined' || !window.bootstrap.Modal) return;
        window.bootstrap.Modal.getOrCreateInstance(modalEl).show();
    }

    window.openFeedbackModal = openFeedbackModal;
    window.ensureFeedbackModal = ensureFeedbackModal;

    /** 紧凑右键菜单：上报错误 */
    let _errMenuEl = null;
    function removeErrorReportMenu() {
        if (_errMenuEl && _errMenuEl.parentNode) _errMenuEl.parentNode.removeChild(_errMenuEl);
        _errMenuEl = null;
    }

    function isErrorLikeTarget(el) {
        if (!el || !el.closest) return false;
        const n = el.closest(
            '.text-danger, .alert-danger, .omics-error-snippet, [data-omics-error], .toast-danger, .checklist-subtitle'
        );
        if (!n) return false;
        const t = (n.innerText || '').trim();
        return t.length >= 4;
    }

    document.addEventListener(
        'contextmenu',
        function (e) {
            const t = e.target;
            if (!isErrorLikeTarget(t)) return;
            const txt = (t.closest('.text-danger, .alert-danger, .omics-error-snippet, [data-omics-error], .toast-danger, .checklist-subtitle') || t).innerText || '';
            if (!String(txt).trim()) return;
            e.preventDefault();
            removeErrorReportMenu();
            const menu = document.createElement('div');
            menu.className = 'omics-error-report-menu popover-menu show';
            menu.style.position = 'fixed';
            menu.style.left = Math.min(e.clientX, window.innerWidth - 200) + 'px';
            menu.style.top = Math.min(e.clientY, window.innerHeight - 48) + 'px';
            menu.style.zIndex = '100100';
            menu.innerHTML =
                '<button type="button" class="popover-item" style="border:none;background:transparent;width:100%;text-align:left;">🚨 上报此错误</button>';
            menu.querySelector('button').onclick = function (ev) {
                ev.preventDefault();
                ev.stopPropagation();
                removeErrorReportMenu();
                openFeedbackModal(String(txt).trim().slice(0, 12000));
            };
            document.body.appendChild(menu);
            _errMenuEl = menu;
            setTimeout(function () {
                document.addEventListener(
                    'click',
                    function once() {
                        removeErrorReportMenu();
                        document.removeEventListener('click', once);
                    },
                    { once: true }
                );
            }, 0);
        },
        true
    );
})();
