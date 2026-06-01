/**
 * 用户站内通知铃铛：轮询未读数 + 下拉列表
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
let _notifPanelOpen = false;

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
    const badge = document.getElementById('user-notifications-badge');
    const bell = document.getElementById('user-notifications-bell');
    if (!badge || !bell) return;
    const n = Number(count) || 0;
    if (n > 0) {
        badge.textContent = n > 99 ? '99+' : String(n);
        badge.style.display = 'inline-flex';
        bell.setAttribute('data-has-unread', '1');
    } else {
        badge.style.display = 'none';
        bell.removeAttribute('data-has-unread');
    }
}

function renderNotificationsPanel(items) {
    const panel = document.getElementById('user-notifications-panel');
    if (!panel) return;
    const list = Array.isArray(items) ? items : [];
    if (!list.length) {
        panel.innerHTML =
            '<div class="user-notifications-panel__empty">暂无消息</div>' +
            '<div class="user-notifications-panel__footer">' +
            '<button type="button" class="user-notifications-panel__read-all" data-notif-read-all="1">全部标为已读</button>' +
            '</div>';
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
                '" role="menuitem">' +
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
    panel.innerHTML =
        '<div class="user-notifications-panel__header">消息通知</div>' +
        '<div class="user-notifications-panel__list">' +
        rows +
        '</div>' +
        '<div class="user-notifications-panel__footer">' +
        '<button type="button" class="user-notifications-panel__read-all" data-notif-read-all="1">全部标为已读</button>' +
        '</div>';
}

function positionNotificationsPanel() {
    const bell = document.getElementById('user-notifications-bell');
    const panel = document.getElementById('user-notifications-panel');
    if (!bell || !panel) return;
    const rect = bell.getBoundingClientRect();
    panel.style.position = 'fixed';
    panel.style.left = Math.max(8, rect.left) + 'px';
    panel.style.bottom = Math.max(8, window.innerHeight - rect.top + 8) + 'px';
    panel.style.top = 'auto';
    panel.style.right = 'auto';
    panel.style.zIndex = '10050';
    panel.style.minWidth = '320px';
    panel.style.maxWidth = 'min(420px, calc(100vw - 16px))';
}

async function refreshUserNotifications() {
    if (typeof window.isLoggedIn === 'function' && !window.isLoggedIn()) {
        const row = document.getElementById('sidebar-notifications-row');
        if (row) row.style.display = 'none';
        updateNotificationBadge(0);
        return;
    }
    const row = document.getElementById('sidebar-notifications-row');
    if (row) row.style.display = 'block';
    try {
        const data = await fetchUnreadNotifications();
        updateNotificationBadge(data.unread_count);
        if (_notifPanelOpen) {
            renderNotificationsPanel(data.items);
            positionNotificationsPanel();
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

function closeNotificationsPanel() {
    const panel = document.getElementById('user-notifications-panel');
    const bell = document.getElementById('user-notifications-bell');
    _notifPanelOpen = false;
    if (panel) {
        panel.style.display = 'none';
        panel.setAttribute('aria-hidden', 'true');
    }
    if (bell) bell.setAttribute('aria-expanded', 'false');
}

async function toggleNotificationsPanel() {
    const panel = document.getElementById('user-notifications-panel');
    const bell = document.getElementById('user-notifications-bell');
    if (!panel || !bell) return;
    if (_notifPanelOpen) {
        closeNotificationsPanel();
        return;
    }
    _notifPanelOpen = true;
    bell.setAttribute('aria-expanded', 'true');
    panel.style.display = 'block';
    panel.setAttribute('aria-hidden', 'false');
    try {
        const res = await fetch('/api/user/notifications?limit=30', {
            method: 'GET',
            headers: _notifAuthHeaders(),
        });
        const data = res.ok ? await res.json() : { items: [] };
        renderNotificationsPanel(data.items || []);
        updateNotificationBadge(data.unread_count);
    } catch (_) {
        renderNotificationsPanel([]);
    }
    positionNotificationsPanel();
}

function initUserNotificationsBell() {
    const bell = document.getElementById('user-notifications-bell');
    const panel = document.getElementById('user-notifications-panel');
    if (!bell || !panel) return;

    bell.addEventListener('click', function (ev) {
        ev.preventDefault();
        ev.stopPropagation();
        if (typeof window.isLoggedIn === 'function' && !window.isLoggedIn()) {
            if (typeof window.openAuthModal === 'function') window.openAuthModal();
            return;
        }
        void toggleNotificationsPanel();
    });

    panel.addEventListener('click', function (ev) {
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

    document.addEventListener('click', function (ev) {
        if (!_notifPanelOpen) return;
        if (ev.target.closest('#user-notifications-panel') || ev.target.closest('#user-notifications-bell')) {
            return;
        }
        closeNotificationsPanel();
    });

    window.addEventListener('resize', function () {
        if (_notifPanelOpen) positionNotificationsPanel();
    });

    if (_notifPollTimer) clearInterval(_notifPollTimer);
    _notifPollTimer = setInterval(function () {
        void refreshUserNotifications();
    }, 60000);

    void refreshUserNotifications();
    window.refreshUserNotifications = refreshUserNotifications;
}

if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', initUserNotificationsBell);
} else {
    initUserNotificationsBell();
}

export { refreshUserNotifications, initUserNotificationsBell };
