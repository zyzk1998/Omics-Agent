/**
 * IM 频道接入弹窗（集成中心占位）
 * 正式接入后在此挂载各平台 OAuth / Webhook 配置表单。
 */
(function () {
    var bindingsDone = false;

    function closeImChannelsModal() {
        var el = document.getElementById('im-channels-modal');
        if (!el) return;
        el.style.display = 'none';
        el.setAttribute('aria-hidden', 'true');
    }

    function initImChannelsModalBindings() {
        if (bindingsDone) return;
        bindingsDone = true;
        var backdrop = document.getElementById('im-channels-modal-backdrop');
        var closeBtn = document.getElementById('im-channels-modal-close');
        if (backdrop) backdrop.onclick = closeImChannelsModal;
        if (closeBtn) closeBtn.onclick = closeImChannelsModal;
        document.addEventListener('keydown', function (ev) {
            if (ev.key !== 'Escape') return;
            var el = document.getElementById('im-channels-modal');
            if (el && el.style.display !== 'none' && el.getAttribute('aria-hidden') === 'false') {
                closeImChannelsModal();
            }
        });
    }

    function openImChannelsModal() {
        if (typeof window.closeAllPopovers === 'function') {
            window.closeAllPopovers();
        }
        initImChannelsModalBindings();
        var el = document.getElementById('im-channels-modal');
        if (!el) return;
        el.style.display = 'flex';
        el.setAttribute('aria-hidden', 'false');
    }

    window.openImChannelsModal = openImChannelsModal;
    window.closeImChannelsModal = closeImChannelsModal;
})();
