/**
 * 用户站内通知：入口在用户头像下拉菜单；列表以居中 Modal 展示
 * API: GET /api/user/notifications/unread, POST .../read, POST .../read-all
 */

function _notifAuthHeaders() {
    if (typeof window.getAuthHeaders === 'function') {
        return window.getAuthHeaders();
    }
    return {};
}

function _notifEscapeHtml(text) {
    if (text == null || text === '') return '';
    return String(text)
        .replace(/&/g, '&amp;')
        .replace(/</g, '&lt;')
        .replace(/>/g, '&gt;')
        .replace(/"/g, '&quot;');
}

let _notifPollTimer = null;

function ensureNotificationsModal() {
    if (document.getElementById('user-notifications-modal')) return;
    const holder = document.createElement('div');
    holder.innerHTML =
        '<div class="modal fade" id="user-notifications-modal" tabindex="-1" aria-labelledby="user-notifications-modal-label" aria-hidden="true">' +
        '  <div class="modal-dialog modal-dialog-centered modal-dialog-scrollable">' +
        '    <div class="modal-content">' +
        '      <div class="modal-header">' +
        '        <h5 class="modal-title" id="user-notifications-modal-label">消息通知</h5>' +
        '        <button type="button" class="btn-close" data-bs-dismiss="modal" aria-label="关闭"></button>' +
        '      </div>' +
        '      <div class="modal-body p-0" id="user-notifications-modal-body"></div>' +
        '      <div class="modal-footer py-2">' +
        '        <button type="button" class="btn btn-sm btn-outline-secondary" data-notif-read-all="1">全部标为已读</button>' +
        '      </div>' +
        '    </div>' +
        '  </div>' +
        '</div>';
    document.body.appendChild(holder.firstElementChild);

    const modalEl = document.getElementById('user-notifications-modal');
    modalEl.addEventListener('click', function (ev) {
        const readAllBtn = ev.target.closest('[data-notif-read-all]');
        if (readAllBtn) {
            ev.preventDefault();
            void markAllNotificationsRead().then(function () {
                return refreshUserNotifications();
            });
            return;
        }
        const itemBtn = ev.target.closest('.user-notifications-item[data-notif-id]');
        if (itemBtn) {
            const nid = itemBtn.getAttribute('data-notif-id');
            if (nid) {
                void markNotificationRead(nid).then(function () {
                    itemBtn.classList.remove('user-notifications-item--unread');
                    return refreshUserNotifications();
                });
            }
        }
    });
}

async function fetchUnreadNotifications() {
    const res = await fetch('/api/user/notifications/unread', {
        method: 'GET',
        headers: _notifAuthHeaders(),
    });
    if (res.status === 401) return { unread_count: 0, items: [] };
    if (!res.ok) throw new Error('加载通知失败');
    return res.json();
}

function updateNotificationBadge(count) {
    const badge = document.getElementById('user-menu-unread-badge');
    const menuItem = document.getElementById('user-menu-notifications');
    const n = Number(count) || 0;
    if (!badge) return;
    if (n > 0) {
        badge.textContent = n > 99 ? '99+' : String(n);
        badge.style.display = 'inline-flex';
        badge.removeAttribute('aria-hidden');
        if (menuItem) menuItem.setAttribute('data-has-unread', '1');
    } else {
        badge.style.display = 'none';
        badge.setAttribute('aria-hidden', 'true');
        if (menuItem) menuItem.removeAttribute('data-has-unread');
    }
}

function renderNotificationsModalBody(items) {
    const body = document.getElementById('user-notifications-modal-body');
    if (!body) return;
    const list = Array.isArray(items) ? items : [];
    if (!list.length) {
        body.innerHTML = '<div class="user-notifications-panel__empty">暂无消息</div>';
        return;
    }
    const rows = list
        .map(function (n) {
            const unread = !n.is_read;
            const cls = unread ? ' user-notifications-item--unread' : '';
            return (
                '<button type="button" class="user-notifications-item' +
                cls +
                '" data-notif-id="' +
                String(n.id) +
                '">' +
                '<div class="user-notifications-item__title">' +
                _notifEscapeHtml(n.title || '通知') +
                '</div>' +
                '<div class="user-notifications-item__content">' +
                _notifEscapeHtml(n.content || '') +
                '</div>' +
                '<div class="user-notifications-item__time">' +
                _notifEscapeHtml(String(n.created_at || '').slice(0, 19)) +
                '</div>' +
                '</button>'
            );
        })
        .join('');
    body.innerHTML =
        '<div class="user-notifications-panel__list">' + rows + '</div>';
}

async function refreshUserNotifications() {
    if (typeof window.isLoggedIn === 'function' && !window.isLoggedIn()) {
        updateNotificationBadge(0);
        return;
    }
    try {
        const data = await fetchUnreadNotifications();
        updateNotificationBadge(data.unread_count);
        const modalEl = document.getElementById('user-notifications-modal');
        if (modalEl && modalEl.classList.contains('show')) {
            const res = await fetch('/api/user/notifications?limit=30', {
                method: 'GET',
                headers: _notifAuthHeaders(),
            });
            const full = res.ok ? await res.json() : { items: [] };
            renderNotificationsModalBody(full.items || []);
        }
    } catch (e) {
        console.warn('[notifications] poll failed', e);
    }
}

async function markNotificationRead(id) {
    await fetch('/api/user/notifications/' + encodeURIComponent(String(id)) + '/read', {
        method: 'POST',
        headers: Object.assign({ 'Content-Type': 'application/json' }, _notifAuthHeaders()),
    });
}

async function markAllNotificationsRead() {
    await fetch('/api/user/notifications/read-all', {
        method: 'POST',
        headers: Object.assign({ 'Content-Type': 'application/json' }, _notifAuthHeaders()),
    });
}

async function openUserNotificationsModal() {
    ensureNotificationsModal();
    const modalEl = document.getElementById('user-notifications-modal');
    if (!modalEl) return;
    if (typeof window.bootstrap === 'undefined' || !window.bootstrap.Modal) {
        console.error('[notifications] Bootstrap Modal 未加载');
        return;
    }
    try {
        const res = await fetch('/api/user/notifications?limit=30', {
            method: 'GET',
            headers: _notifAuthHeaders(),
        });
        const data = res.ok ? await res.json() : { items: [], unread_count: 0 };
        renderNotificationsModalBody(data.items || []);
        updateNotificationBadge(data.unread_count);
    } catch (_) {
        renderNotificationsModalBody([]);
    }
    window.bootstrap.Modal.getOrCreateInstance(modalEl).show();
}

function initUserNotifications() {
    ensureNotificationsModal();
    if (_notifPollTimer) clearInterval(_notifPollTimer);
    _notifPollTimer = setInterval(function () {
        void refreshUserNotifications();
    }, 60000);
    void refreshUserNotifications();
    window.refreshUserNotifications = refreshUserNotifications;
    window.openUserNotificationsModal = openUserNotificationsModal;
}

if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', initUserNotifications);
} else {
    initUserNotifications();
}

export { refreshUserNotifications, openUserNotificationsModal, initUserNotifications };
