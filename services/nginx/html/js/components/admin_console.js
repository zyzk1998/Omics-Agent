/**
 * 管理员控制台：Bootstrap Modal + Tab（技能审核 | 用户审核）
 * 技能：GET /api/admin/skills、POST .../review
 * 用户：GET /api/admin/users/pending、POST .../approve|reject
 */

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

function setInlineError(el, message) {
    if (!el) return;
    if (message) {
        el.textContent = message;
        el.classList.remove('d-none');
    } else {
        el.textContent = '';
        el.classList.add('d-none');
    }
}

let _adminActiveTab = 'skills';

function ensureAdminModal() {
    if (document.getElementById('admin-console-modal')) return;
    const holder = document.createElement('div');
    holder.innerHTML = `
<div class="modal fade" id="admin-console-modal" tabindex="-1" aria-labelledby="admin-console-modal-label" aria-hidden="true">
  <div class="modal-dialog modal-xl modal-dialog-scrollable">
    <div class="modal-content">
      <div class="modal-header">
        <h5 class="modal-title" id="admin-console-modal-label">管理员控制台</h5>
        <button type="button" class="btn-close" data-bs-dismiss="modal" aria-label="关闭"></button>
      </div>
      <div class="modal-body d-flex flex-column" style="min-height: 280px;">
        <div id="admin-console-error" class="alert alert-danger d-none py-2" role="alert"></div>
        <ul class="nav nav-tabs mb-2 flex-shrink-0" id="admin-console-tabs" role="tablist">
          <li class="nav-item" role="presentation">
            <button type="button" class="nav-link active" data-admin-tab="skills" id="admin-tab-btn-skills">技能审核</button>
          </li>
          <li class="nav-item" role="presentation">
            <button type="button" class="nav-link" data-admin-tab="users" id="admin-tab-btn-users">用户审核</button>
          </li>
        </ul>
        <div id="admin-panel-skills" class="admin-console-panel flex-fill" style="min-height: 0;">
          <div id="admin-console-table-wrap"><p class="text-muted mb-0">加载中…</p></div>
        </div>
        <div id="admin-panel-users" class="admin-console-panel flex-fill d-none" style="min-height: 0;">
          <div id="admin-users-table-wrap"><p class="text-muted mb-0">加载中…</p></div>
        </div>
      </div>
    </div>
  </div>
</div>`;
    document.body.appendChild(holder.firstElementChild);
    const wrap = document.getElementById('admin-console-table-wrap');
    wrap.addEventListener('click', (ev) => {
        const btn = ev.target.closest('[data-review-action]');
        if (!btn) return;
        const id = btn.getAttribute('data-skill-id');
        const action = btn.getAttribute('data-review-action');
        if (id && (action === 'approve' || action === 'reject')) {
            void reviewSkill(id, action);
        }
    });
    const uwrap = document.getElementById('admin-users-table-wrap');
    uwrap.addEventListener('click', (ev) => {
        const btn = ev.target.closest('[data-user-action]');
        if (!btn) return;
        const id = btn.getAttribute('data-user-id');
        const action = btn.getAttribute('data-user-action');
        if (id && (action === 'approve' || action === 'reject')) {
            void reviewUser(id, action);
        }
    });
    document.getElementById('admin-console-tabs').addEventListener('click', (ev) => {
        const b = ev.target.closest('[data-admin-tab]');
        if (!b) return;
        const tab = b.getAttribute('data-admin-tab');
        if (tab === 'skills' || tab === 'users') {
            setAdminConsoleTab(tab);
        }
    });
}

function setAdminConsoleTab(tab) {
    _adminActiveTab = tab;
    const skillsBtn = document.getElementById('admin-tab-btn-skills');
    const usersBtn = document.getElementById('admin-tab-btn-users');
    const panelSkills = document.getElementById('admin-panel-skills');
    const panelUsers = document.getElementById('admin-panel-users');
    if (skillsBtn) skillsBtn.classList.toggle('active', tab === 'skills');
    if (usersBtn) usersBtn.classList.toggle('active', tab === 'users');
    if (panelSkills) panelSkills.classList.toggle('d-none', tab !== 'skills');
    if (panelUsers) panelUsers.classList.toggle('d-none', tab !== 'users');
    if (tab === 'skills') void refreshAdminTable();
    else void refreshAdminUsersTable();
}

let _adminLoadSeq = 0;
let _adminUsersLoadSeq = 0;

async function fetchPendingSkills() {
    const res = await fetch('/api/admin/skills', { method: 'GET', headers: getAuthHeadersMerged() });
    if (res.status === 403) throw new Error('无管理员权限');
    if (!res.ok) {
        let detail = '加载失败';
        try {
            const j = await res.json();
            if (j && j.detail != null) {
                detail = typeof j.detail === 'string' ? j.detail : JSON.stringify(j.detail);
            }
        } catch (_) {
            /* ignore */
        }
        throw new Error(detail);
    }
    return res.json();
}

async function fetchPendingUsers() {
    const res = await fetch('/api/admin/users/pending', { method: 'GET', headers: getAuthHeadersMerged() });
    if (res.status === 403) throw new Error('无管理员权限');
    if (!res.ok) {
        let detail = '加载失败';
        try {
            const j = await res.json();
            if (j && j.detail != null) {
                detail = typeof j.detail === 'string' ? j.detail : JSON.stringify(j.detail);
            }
        } catch (_) {
            /* ignore */
        }
        throw new Error(detail);
    }
    return res.json();
}

async function reviewSkill(skillId, action) {
    const errEl = document.getElementById('admin-console-error');
    setInlineError(errEl, '');
    try {
        const res = await fetch(
            '/api/admin/skills/' + encodeURIComponent(String(skillId)) + '/review',
            {
                method: 'POST',
                headers: getAuthHeadersMerged(),
                body: JSON.stringify({ action }),
            }
        );
        if (res.status === 403) throw new Error('需要管理员权限');
        if (!res.ok) {
            let detail = '操作失败';
            try {
                const j = await res.json();
                if (j && j.detail != null) {
                    detail = typeof j.detail === 'string' ? j.detail : JSON.stringify(j.detail);
                }
            } catch (_) {
                /* ignore */
            }
            throw new Error(detail);
        }
        await refreshAdminTable();
    } catch (e) {
        const msg = e && e.message ? String(e.message) : '操作失败';
        setInlineError(errEl, msg);
    }
}

async function reviewUser(userId, action) {
    const errEl = document.getElementById('admin-console-error');
    setInlineError(errEl, '');
    const path = action === 'approve' ? 'approve' : 'reject';
    try {
        const res = await fetch(
            '/api/admin/users/' + encodeURIComponent(String(userId)) + '/' + path,
            {
                method: 'POST',
                headers: getAuthHeadersMerged(),
            }
        );
        if (res.status === 403) throw new Error('需要管理员权限');
        if (!res.ok) {
            let detail = '操作失败';
            try {
                const j = await res.json();
                if (j && j.detail != null) {
                    detail = typeof j.detail === 'string' ? j.detail : JSON.stringify(j.detail);
                }
            } catch (_) {
                /* ignore */
            }
            throw new Error(detail);
        }
        await refreshAdminUsersTable();
    } catch (e) {
        const msg = e && e.message ? String(e.message) : '操作失败';
        setInlineError(errEl, msg);
    }
}

async function refreshAdminTable() {
    const wrap = document.getElementById('admin-console-table-wrap');
    const errEl = document.getElementById('admin-console-error');
    if (!wrap) return;
    const seq = ++_adminLoadSeq;
    setInlineError(errEl, '');
    wrap.innerHTML = '<p class="text-muted mb-0">加载中…</p>';
    try {
        const list = await fetchPendingSkills();
        if (seq !== _adminLoadSeq) return;
        if (!Array.isArray(list) || list.length === 0) {
            wrap.innerHTML = '<p class="text-muted mb-0">暂无待审核技能</p>';
            return;
        }
        const rows = list
            .map((s) => {
                const id = s.id;
                const name = escapeHtml(s.name || '');
                const mc = escapeHtml(s.main_category || '—');
                const sc = escapeHtml(s.sub_category || '—');
                const author = escapeHtml(s.author_id || '—');
                const created = escapeHtml(String(s.created_at || '').slice(0, 19));
                const rawDesc = s.description || '';
                const descShort = escapeHtml(rawDesc.slice(0, 100));
                const ell = rawDesc.length > 100 ? '…' : '';
                return (
                    '<tr>' +
                    '<td>' +
                    name +
                    '</td><td>' +
                    mc +
                    '</td><td>' +
                    sc +
                    '</td><td>' +
                    author +
                    '</td><td><small class="text-secondary">' +
                    created +
                    '</small></td><td><small class="text-secondary">' +
                    descShort +
                    ell +
                    '</small></td>' +
                    '<td class="text-nowrap">' +
                    '<button type="button" class="btn btn-sm btn-success me-1" data-skill-id="' +
                    id +
                    '" data-review-action="approve">通过</button>' +
                    '<button type="button" class="btn btn-sm btn-outline-danger" data-skill-id="' +
                    id +
                    '" data-review-action="reject">驳回</button>' +
                    '</td></tr>'
                );
            })
            .join('');
        wrap.innerHTML =
            '<table class="table table-hover table-sm align-middle mb-0">' +
            '<thead><tr><th>名称</th><th>大类</th><th>标签</th><th>作者</th><th>提交时间</th><th>描述摘要</th><th>操作</th></tr></thead>' +
            '<tbody>' +
            rows +
            '</tbody></table>';
    } catch (e) {
        if (seq !== _adminLoadSeq) return;
        wrap.innerHTML = '';
        const msg = e && e.message ? String(e.message) : '加载失败';
        setInlineError(errEl, msg);
    }
}

async function refreshAdminUsersTable() {
    const wrap = document.getElementById('admin-users-table-wrap');
    const errEl = document.getElementById('admin-console-error');
    if (!wrap) return;
    const seq = ++_adminUsersLoadSeq;
    setInlineError(errEl, '');
    wrap.innerHTML = '<p class="text-muted mb-0">加载中…</p>';
    try {
        const list = await fetchPendingUsers();
        if (seq !== _adminUsersLoadSeq) return;
        if (!Array.isArray(list) || list.length === 0) {
            wrap.innerHTML = '<p class="text-muted mb-0">暂无待审核用户</p>';
            return;
        }
        const rows = list
            .map((u) => {
                const id = u.id;
                const un = escapeHtml(u.username || '');
                const em = escapeHtml(u.email != null && u.email !== '' ? u.email : '—');
                const created = escapeHtml(String(u.created_at || '').slice(0, 19));
                return (
                    '<tr>' +
                    '<td>' +
                    id +
                    '</td><td>' +
                    un +
                    '</td><td>' +
                    em +
                    '</td><td><small class="text-secondary">' +
                    created +
                    '</small></td>' +
                    '<td class="text-nowrap">' +
                    '<button type="button" class="btn btn-sm btn-success me-1" data-user-id="' +
                    id +
                    '" data-user-action="approve">通过</button>' +
                    '<button type="button" class="btn btn-sm btn-outline-danger" data-user-id="' +
                    id +
                    '" data-user-action="reject">拒绝</button>' +
                    '</td></tr>'
                );
            })
            .join('');
        wrap.innerHTML =
            '<table class="table table-hover table-sm align-middle mb-0">' +
            '<thead><tr><th>ID</th><th>用户名</th><th>邮箱</th><th>注册时间</th><th>操作</th></tr></thead>' +
            '<tbody>' +
            rows +
            '</tbody></table>';
    } catch (e) {
        if (seq !== _adminUsersLoadSeq) return;
        wrap.innerHTML = '';
        const msg = e && e.message ? String(e.message) : '加载失败';
        setInlineError(errEl, msg);
    }
}

function showAdminConsole() {
    try {
        ensureAdminModal();
        const modalEl = document.getElementById('admin-console-modal');
        if (!modalEl) return;
        if (typeof window.bootstrap === 'undefined' || !window.bootstrap.Modal) {
            console.error('[admin_console] Bootstrap Modal 未加载');
            return;
        }
        const modal = window.bootstrap.Modal.getOrCreateInstance(modalEl);
        modal.show();
        setAdminConsoleTab(_adminActiveTab || 'skills');
    } catch (e) {
        console.error('[admin_console]', e);
    }
}

window.showAdminConsole = showAdminConsole;
window.fetchPendingSkills = fetchPendingSkills;
window.reviewSkill = reviewSkill;
window.fetchPendingUsers = fetchPendingUsers;
window.reviewUser = reviewUser;
