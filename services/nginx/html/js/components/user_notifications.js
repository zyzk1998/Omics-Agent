/**
 * 用户/管理员站内通知：头像菜单「消息通知」+ 右下角 Toast
 * API: GET /api/user/notifications/unread, POST .../read, POST .../read-all
 */

const _SEEN_NOTIF_STORAGE_KEY = 'omics_seen_notif_ids';
const _POLL_INTERVAL_MS = 45000;

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

function _hasNotificationIdentity() {
    try {
        if (localStorage.getItem('access_token')) return true;
        if (localStorage.getItem('guest_uuid')) return true;
    } catch (_e) {
        /* ignore */
    }
    return typeof window.getAuthHeaders === 'function';
}

function _loadSeenNotifIds() {
    try {
        const raw = sessionStorage.getItem(_SEEN_NOTIF_STORAGE_KEY);
        if (!raw) return new Set();
        const arr = JSON.parse(raw);
        if (!Array.isArray(arr)) return new Set();
        return new Set(arr.map(function (x) {
            return String(x);
        }));
    } catch (_e) {
        return new Set();
    }
}

function _saveSeenNotifIds(set) {
    try {
        const arr = Array.from(set).slice(-200);
        sessionStorage.setItem(_SEEN_NOTIF_STORAGE_KEY, JSON.stringify(arr));
    } catch (_e) {
        /* ignore */
    }
}

let _notifPollTimer = null;
let _seenNotifIds = _loadSeenNotifIds();
let _notifBootstrapped = false;

function showNotificationToast(title, content, meta) {
    const existing = document.getElementById('omics-notif-toast');
    if (existing) existing.remove();
    const el = document.createElement('div');
    el.id = 'omics-notif-toast';
    el.className = 'omics-notif-toast';
    el.setAttribute('role', 'status');
    const t = _notifEscapeHtml(title || '新消息');
    const c = _notifEscapeHtml(content || '');
    el.innerHTML =
        '<button type="button" class="omics-notif-toast__close" aria-label="关闭">&times;</button>' +
        '<div class="omics-notif-toast__title">' +
        t +
        '</div>' +
        (c ? '<div class="omics-notif-toast__content">' + c + '</div>' : '') +
        '<div class="omics-notif-toast__hint">点击查看全部</div>';
    el.querySelector('.omics-notif-toast__close').addEventListener('click', function (ev) {
        ev.stopPropagation();
        el.remove();
    });
    el.addEventListener('click', function () {
        el.remove();
        void openUserNotificationsModal();
        if (meta && meta.type && String(meta.type).indexOf('admin_') === 0) {
            if (typeof window.showAdminConsole === 'function') {
                window.showAdminConsole();
            }
        }
    });
    document.body.appendChild(el);
    setTimeout(function () {
        if (el.parentNode) {
            el.classList.add('omics-notif-toast--fade');
            setTimeout(function () {
                if (el.parentNode) el.remove();
            }, 320);
        }
    }, 8000);
}

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

function updateNotificationEntryVisibility() {
    const show = _hasNotificationIdentity();
    const menuItem = document.getElementById('user-menu-notifications');
    if (menuItem) menuItem.style.display = show ? 'flex' : 'none';
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
    const n = Number(count) || 0;
    const label = n > 99 ? '99+' : String(n);
    const badge = document.getElementById('user-menu-unread-badge');
    if (badge) {
        if (n > 0) {
            badge.textContent = label;
            badge.style.display = 'inline-flex';
            badge.removeAttribute('aria-hidden');
        } else {
            badge.style.display = 'none';
            badge.setAttribute('aria-hidden', 'true');
        }
    }
    const menuItem = document.getElementById('user-menu-notifications');
    if (menuItem) {
        if (n > 0) menuItem.setAttribute('data-has-unread', '1');
        else menuItem.removeAttribute('data-has-unread');
    }
}

function _toastNewNotifications(items) {
    if (!Array.isArray(items) || !items.length) return;
    let newest = null;
    items.forEach(function (n) {
        const id = String(n.id);
        if (_seenNotifIds.has(id)) return;
        _seenNotifIds.add(id);
        if (!newest) newest = n;
    });
    _saveSeenNotifIds(_seenNotifIds);
    if (!newest) return;
    if (!_notifBootstrapped) return;
    showNotificationToast(newest.title, newest.content, { type: newest.type, id: newest.id });
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
    body.innerHTML = '<div class="user-notifications-panel__list">' + rows + '</div>';
}

async function refreshUserNotifications() {
    updateNotificationEntryVisibility();
    if (!_hasNotificationIdentity()) {
        updateNotificationBadge(0);
        return;
    }
    try {
        const data = await fetchUnreadNotifications();
        updateNotificationBadge(data.unread_count);
        _toastNewNotifications(data.items || []);
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
    } finally {
        _notifBootstrapped = true;
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
        (data.items || []).forEach(function (n) {
            if (n && n.id != null) _seenNotifIds.add(String(n.id));
        });
        _saveSeenNotifIds(_seenNotifIds);
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
    }, _POLL_INTERVAL_MS);
    void refreshUserNotifications();
    window.refreshUserNotifications = refreshUserNotifications;
    window.openUserNotificationsModal = openUserNotificationsModal;
    window.showNotificationToast = showNotificationToast;
}

if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', initUserNotifications);
} else {
    initUserNotifications();
}

export { refreshUserNotifications, openUserNotificationsModal, initUserNotifications, showNotificationToast };
