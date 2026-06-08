/**
 * 统一 Auth 请求头（所有 components 必须通过此模块携带 Token / Guest UUID）
 */
(function () {
    'use strict';

    function resolveAuthHeaders() {
        if (typeof window.getAuthHeaders === 'function') {
            return window.getAuthHeaders();
        }
        if (typeof window.authHeaders === 'function') {
            return window.authHeaders();
        }
        return {};
    }

    window.resolveAuthHeaders = resolveAuthHeaders;

    window.mergeAuthHeaders = function mergeAuthHeaders(extra) {
        return Object.assign({}, resolveAuthHeaders(), extra || {});
    };

    window.mergeJsonAuthHeaders = function mergeJsonAuthHeaders() {
        return mergeAuthHeaders({ 'Content-Type': 'application/json' });
    };

    /** 401 / 登出 / 打开登录框前：拆除可能挡住登录 UI 的遮罩 */
    window.teardownAuthBlockingOverlays = function teardownAuthBlockingOverlays() {
        document.querySelectorAll('.ingestion-result-modal-overlay').forEach(function (el) {
            el.remove();
        });
        document.body.classList.remove('expert-report-fullscreen-active');
        var erEdit = document.getElementById('expert-report-fullscreen-editor');
        if (erEdit) {
            erEdit.classList.remove('is-open');
            erEdit.setAttribute('aria-hidden', 'true');
        }
    };

    if (document.readyState === 'loading') {
        document.addEventListener('DOMContentLoaded', teardownAuthBlockingOverlays);
    } else {
        teardownAuthBlockingOverlays();
    }
})();
