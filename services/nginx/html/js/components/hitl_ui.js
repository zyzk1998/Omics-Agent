/**
 * HITL Phase 2/3：独立 LS 入口卡片 + iframe 模态 + resume_from_hitl 唤醒
 */
(function () {
    'use strict';

    if (!window.__LABEL_STUDIO_PUBLIC_URL) {
        window.__LABEL_STUDIO_PUBLIC_URL = '/label-studio';
    }

    var _lastHitlData = null;
    var _resumeInFlight = false;
    var _lsSubmitDetected = false;
    var _hitlCompletedMode = false;

    function authHeadersMerge() {
        if (typeof window.mergeJsonAuthHeaders === 'function') return window.mergeJsonAuthHeaders();
        return Object.assign(
            { 'Content-Type': 'application/json' },
            typeof getAuthHeaders === 'function' ? getAuthHeaders() : (typeof authHeaders === 'function' ? authHeaders() : {})
        );
    }

    function iconHtml(name) {
        return '<i class="bi bi-' + name + '" aria-hidden="true"></i>';
    }

    function normalizeLsIframeUrl(url) {
        var raw = String(url || '').trim();
        if (!raw) return '';
        try {
            var parsed = new URL(raw, window.location.origin);
            var match = parsed.pathname.match(/\/projects\/(\d+)\/data\/?$/);
            if (match && parsed.origin !== window.location.origin) {
                return '/label-studio/projects/' + match[1] + '/data';
            }
            if (parsed.origin === window.location.origin) {
                return parsed.pathname + parsed.search + parsed.hash;
            }
            return parsed.href;
        } catch (e) {
            if (raw.charAt(0) === '/') return raw;
            return raw;
        }
    }

    function resolveLsOpenUrl(payload) {
        if (!payload || typeof payload !== 'object') return '';
        var url = normalizeLsIframeUrl(payload.ls_url || payload.ls_project_url || '');
        if (url) return url;
        var pid = payload.project_id != null ? payload.project_id : payload.ls_project_id;
        if (pid == null) return '';
        var base = String(window.__LABEL_STUDIO_PUBLIC_URL || '/label-studio').replace(/\/$/, '');
        if (!base) return '';
        if (base.charAt(0) === '/') {
            return base + '/projects/' + encodeURIComponent(String(pid)) + '/data';
        }
        return base + '/projects/' + encodeURIComponent(String(pid)) + '/data';
    }

    function isHitlLaunchUnavailable(payload) {
        if (!payload) return true;
        if (payload.ls_unavailable) return true;
        if (String(payload.status || '').toLowerCase() === 'error') return true;
        var url = resolveLsOpenUrl(payload);
        if (!url) return true;
        return false;
    }

    function escapeHitlHtml(text) {
        return String(text == null ? '' : text)
            .replace(/&/g, '&amp;')
            .replace(/</g, '&lt;')
            .replace(/>/g, '&gt;')
            .replace(/"/g, '&quot;');
    }

    function resolveHitlErrorMessage(payload) {
        if (!payload || typeof payload !== 'object') {
            return 'Label Studio 服务未就绪，请确认 label-studio 容器已启动并通过健康检查。';
        }
        var raw = String(
            payload.error_detail
            || payload.message
            || payload.error
            || ''
        ).trim();
        if (!raw) {
            return 'Label Studio 服务未就绪（无 project_id / 无法建立连接）。';
        }
        if (/LABEL_STUDIO_API_KEY|API Key|Account & Settings|手动.*Token/i.test(raw)) {
            return 'Label Studio 服务未就绪或网络不可达，请检查 label-studio 容器状态与 gibh-network 连通性。';
        }
        return raw;
    }

    function resolveHitlMountAnchor() {
        var zones = window.currentMessageZones;
        if (zones && zones.logZone) return zones.logZone;
        if (zones && zones.cardZone) return zones.cardZone;
        var root = typeof getActiveChatContentEl === 'function' ? getActiveChatContentEl() : null;
        if (!root) root = document.getElementById('skillLiveContent') || document.getElementById('chatContent');
        if (!root) return null;
        var rows = root.querySelectorAll('.message-row.ai');
        for (var i = rows.length - 1; i >= 0; i--) {
            var row = rows[i];
            var logZone = row.querySelector('.process-log-zone');
            if (logZone) return logZone;
            var cardZone = row.querySelector('.cards-zone') || row.querySelector('.workflow-container');
            if (cardZone) return cardZone;
        }
        return null;
    }

    function ensureHitlArtifactRow(anchorEl) {
        if (!anchorEl) return null;
        if (typeof window.ensureHitlArtifactRow === 'function') {
            return window.ensureHitlArtifactRow(anchorEl);
        }
        var row = anchorEl.querySelector(':scope > .chat-action-cards-wrapper');
        if (!row) row = anchorEl.querySelector(':scope > .hitl-artifact-row');
        if (!row) {
            row = document.createElement('div');
            row.className = 'hitl-artifact-row chat-action-cards-wrapper chat-action-card-group';
            var left = document.createElement('div');
            left.className = 'hitl-artifact-row__left action-card-slot';
            var right = document.createElement('div');
            right.className = 'hitl-artifact-row__right action-card-slot';
            row.appendChild(left);
            row.appendChild(right);
            anchorEl.appendChild(row);
        } else {
            row.classList.add('chat-action-cards-wrapper', 'chat-action-card-group');
        }
        return {
            row: row,
            left: row.querySelector('.hitl-artifact-row__left'),
            right: row.querySelector('.hitl-artifact-row__right')
        };
    }

    function normalizeHitlPayload(data) {
        if (!data || typeof data !== 'object') return null;
        var status = String(data.status || '').toLowerCase();
        if (
            status === 'hitl_required'
            || status === 'hitl_completed'
            || data.ls_url
            || data.ls_project_url
            || data.project_id != null
            || data.ls_project_id != null
            || data.ls_unavailable
        ) {
            var merged = {
                status: 'hitl_required',
                ls_url: data.ls_url || data.ls_project_url || '',
                ls_project_url: data.ls_project_url || data.ls_url || '',
                project_id: data.project_id != null ? data.project_id : data.ls_project_id,
                ls_project_id: data.ls_project_id != null ? data.ls_project_id : data.project_id,
                scenario_type: data.scenario_type || 'scrna_cell_type_annotation',
                message: data.message || data.error_detail || '检测到需人工确认的步骤，请在 Label Studio 中完成专家复核。',
                error_detail: data.error_detail || data.message || '',
                ls_unavailable: !!data.ls_unavailable || String(data.status || '').toLowerCase() === 'error',
            };
            var resolved = resolveLsOpenUrl(merged);
            if (resolved) {
                merged.ls_url = resolved;
                merged.ls_project_url = resolved;
            }
            if (String(data.status || '').toLowerCase() === 'error') {
                merged.ls_unavailable = true;
            }
            return merged;
        }
        return null;
    }

    function extractHitlFromSteps(steps) {
        if (!Array.isArray(steps)) return null;
        for (var i = steps.length - 1; i >= 0; i--) {
            var step = steps[i];
            if (!step || String(step.status || '').toLowerCase() !== 'hitl_required') continue;
            var sr = step.step_result || {};
            var data = sr.data || sr.result || sr;
            var hitl = normalizeHitlPayload(data) || normalizeHitlPayload(sr);
            if (hitl) return hitl;
        }
        return null;
    }

    window.scanAndRenderHitlFromSteps = function scanAndRenderHitlFromSteps(steps) {
        var hitl = extractHitlFromSteps(steps);
        if (hitl) window.renderHitlActionCard(hitl);
        return hitl;
    };

    function setHitlEntryCardBusy(busy, message) {
        document.querySelectorAll('.hitl-entry-card').forEach(function (card) {
            if (busy) {
                card.classList.add('hitl-entry-card--busy');
                var actions = card.querySelector('.hitl-entry-card__actions');
                if (actions) {
                    actions.innerHTML =
                        '<span class="text-muted small">' +
                        iconHtml('arrow-repeat') + ' ' + (message || '正在生成最终版报告…') +
                        '</span>';
                }
            }
        });
    }

    function extractHitlFromSnapshot(snap) {
        if (!snap || typeof snap !== 'object') return null;
        var hitl = snap.hitl;
        if (hitl && typeof hitl === 'object') {
            return normalizeHitlPayload(Object.assign({ status: 'hitl_completed' }, hitl));
        }
        var ex = snap.execution_snapshot && typeof snap.execution_snapshot === 'object' ? snap.execution_snapshot : {};
        var fromSteps = extractHitlFromSteps(ex.steps_details || []);
        if (fromSteps) return fromSteps;
        if (window.__lastStateSnapshotHitl) {
            return normalizeHitlPayload(window.__lastStateSnapshotHitl);
        }
        return null;
    }

    function removeHitlEntryCards() {
        document.querySelectorAll('.hitl-entry-card').forEach(function (card) {
            card.remove();
        });
    }
    window.removeHitlEntryCards = removeHitlEntryCards;

    function mountHitlCompletedCard(payload) {
        if (!payload) return;
        removeHitlEntryCards();
        var anchor = resolveHitlMountAnchor();
        if (!anchor) {
            setTimeout(function () { mountHitlCompletedCard(payload); }, 400);
            return;
        }
        var slots = ensureHitlArtifactRow(anchor);
        if (!slots || !slots.right) return;

        var lsUrl = resolveLsOpenUrl(payload);
        var unavailable = isHitlLaunchUnavailable(payload);
        var errText = unavailable ? resolveHitlErrorMessage(payload) : '';
        var right = slots.right;
        var card = right.querySelector('.hitl-entry-card');
        if (!card) {
            card = document.createElement('div');
            card.className = 'hitl-entry-card action-card card shadow-sm mt-2';
            card.setAttribute('data-hitl-entry-card', '1');
            right.appendChild(card);
        }
        card.classList.remove('hitl-entry-card--disabled', 'hitl-entry-card--skipped', 'hitl-entry-card--busy');
        card.classList.add('hitl-entry-card--completed');
        card.innerHTML =
            '<div class="card-body p-3 chat-action-card__body d-flex flex-column h-100">' +
            '<div class="chat-action-card__meta hitl-entry-card__head">' +
            '<h6 class="hitl-entry-card__title mb-1 text-success">' + iconHtml('check-circle-fill') + ' 专家复核已完成</h6>' +
            '<p class="hitl-entry-card__desc text-muted mb-0 chat-action-card__subtitle">《专家分析报告（最终版）》已生成。如需修正标注，可重新打开 Label Studio 修改已提交数据，再次生成报告。</p>' +
            (unavailable && errText
                ? ('<div class="hitl-entry-card__error mt-2" role="alert">' +
                    iconHtml('exclamation-octagon-fill') + ' <strong>Label Studio 暂不可用</strong><br>' +
                    '<span class="hitl-entry-card__error-detail">' + escapeHitlHtml(errText) + '</span></div>')
                : '') +
            '</div>' +
            '<div class="chat-action-card__actions hitl-entry-card__actions mt-auto w-100 d-flex flex-wrap gap-2 justify-content-end">' +
            '<button type="button" class="btn btn-sm btn-outline-primary hitl-entry-reannotate-btn">' +
            iconHtml('arrow-repeat') + ' 重新标注</button>' +
            '</div></div>';

        _lastHitlData = payload;
        _hitlCompletedMode = true;
        window.__hitlAnnotateDisabled = false;
        window.__lastStateSnapshotHitl = payload;

        var reBtn = card.querySelector('.hitl-entry-reannotate-btn');
        if (reBtn) {
            if (unavailable || !lsUrl) {
                reBtn.disabled = true;
                reBtn.title = errText || 'Label Studio 暂不可用';
            } else {
                reBtn.disabled = false;
                reBtn.title = '重新打开 Label Studio 修改标注';
            }
            reBtn.addEventListener('click', function (ev) {
                ev.preventDefault();
                ev.stopPropagation();
                if (unavailable || !lsUrl) return;
                _lsSubmitDetected = false;
                openLabelStudioModal(lsUrl);
            });
        }
        slots.row.style.display = 'flex';
    }
    window.mountHitlCompletedCard = mountHitlCompletedCard;

    function isCorpusProcessingHitl(payload) {
        payload = payload || _lastHitlData || window.__lastStateSnapshotHitl || {};
        if (String(payload.scenario_type || '') === 'generic_corpus_processing') return true;
        if (String(payload.workflow_name || '') === '科学语料数据加工') return true;
        if (String(payload.skill_id || '') === 'skill_corpus_data_processing') return true;
        var skill = window._activeSkillPayload;
        if (skill && (skill.tool_id === 'skill_corpus_data_processing' || skill.id === 'skill_corpus_data_processing')) {
            return true;
        }
        return false;
    }

    function syncHitlModalCorpusSkipButton(payload) {
        var actions = document.querySelector('.hitl-ls-modal__footer-actions');
        if (!actions) return;
        var existing = document.getElementById('hitl-ls-modal-skip-btn');
        if (!isCorpusProcessingHitl(payload)) {
            if (existing) existing.remove();
            return;
        }
        if (!existing) {
            existing = document.createElement('button');
            existing.type = 'button';
            existing.id = 'hitl-ls-modal-skip-btn';
            existing.className = 'btn btn-outline-secondary btn-sm me-2 hitl-ls-modal__skip-btn';
            existing.textContent = '取消修改';
            existing.addEventListener('click', function (ev) {
                ev.preventDefault();
                closeLabelStudioModal();
                skipHitlReview();
            });
            actions.insertBefore(existing, actions.firstChild);
        }
    }

    function mountHitlEntryCard(payload) {
        if (!payload) return;
        removeHitlEntryCards();
        var anchor = resolveHitlMountAnchor();
        if (!anchor) {
            setTimeout(function () { mountHitlEntryCard(payload); }, 400);
            return;
        }
        var slots = ensureHitlArtifactRow(anchor);
        if (!slots || !slots.right) return;

        var lsUrl = resolveLsOpenUrl(payload);
        var unavailable = isHitlLaunchUnavailable(payload);
        var errText = unavailable ? resolveHitlErrorMessage(payload) : '';
        var isCorpusSkill = isCorpusProcessingHitl(payload);
        var cardTitle = isCorpusSkill
            ? iconHtml('tags') + ' 科学语料数据加工 / Label Studio'
            : iconHtml('bullseye') + ' 专家复核 / Label Studio';
        var cardSubtitle = isCorpusSkill
            ? '点击进入 Label Studio 进行可选人工标注'
            : '点击进入 Label Studio 完成专家复核';
        var skipLabel = isCorpusSkill ? '取消修改' : '取消修改';
        var right = slots.right;
        var card = right.querySelector('.hitl-entry-card');
        if (!card) {
            card = document.createElement('div');
            card.className = 'hitl-entry-card action-card card shadow-sm mt-2';
            card.setAttribute('data-hitl-entry-card', '1');
            right.appendChild(card);
        }
        card.classList.remove('hitl-entry-card--disabled', 'hitl-entry-card--skipped', 'hitl-entry-card--busy');
        card.innerHTML =
            '<div class="card-body p-3 chat-action-card__body d-flex flex-column h-100">' +
            '<div class="chat-action-card__meta hitl-entry-card__meta">' +
            '<h6 class="hitl-entry-card__title mb-1">' + cardTitle + '</h6>' +
            '<small class="text-muted chat-action-card__subtitle">' + cardSubtitle + '</small>' +
            (unavailable && errText
                ? ('<div class="hitl-entry-card__error mt-2" role="alert">' +
                    iconHtml('exclamation-octagon-fill') + ' <strong>Label Studio 未就绪</strong><br>' +
                    '<span class="hitl-entry-card__error-detail">' + escapeHitlHtml(errText) + '</span></div>')
                : '') +
            '</div>' +
            '<div class="chat-action-card__actions hitl-entry-card__actions mt-auto w-100 d-flex justify-content-end flex-wrap gap-2">' +
            '<button type="button" class="btn btn-sm btn-outline-secondary hitl-entry-skip-btn">' + skipLabel + '</button>' +
            '<button type="button" class="btn btn-sm btn-warning hitl-entry-open-btn">' +
            iconHtml('box-arrow-up-right') + ' 进入手动标注</button>' +
            '</div></div>';
        if (unavailable) {
            card.classList.add('hitl-entry-card--disabled');
        }

        var openBtn = card.querySelector('.hitl-entry-open-btn');
        var skipBtn = card.querySelector('.hitl-entry-skip-btn');
        if (openBtn) {
            if (unavailable) {
                openBtn.disabled = true;
                openBtn.classList.add('disabled');
                openBtn.title = errText || 'Label Studio 暂不可用';
            } else {
                openBtn.disabled = false;
                openBtn.classList.remove('disabled');
                openBtn.title = '在弹窗中打开 Label Studio';
            }
            openBtn.addEventListener('click', function (ev) {
                ev.preventDefault();
                ev.stopPropagation();
                if (unavailable || !lsUrl) {
                    console.error('[HITL] 进入手动标注被阻止:', errText || payload);
                    return;
                }
                _lsSubmitDetected = false;
                openLabelStudioModal(lsUrl);
            });
        }
        if (skipBtn) {
            skipBtn.addEventListener('click', function (ev) {
                ev.preventDefault();
                ev.stopPropagation();
                skipHitlReview();
            });
        }
        syncHitlModalCorpusSkipButton(payload);
        slots.row.style.display = 'flex';
    }

    function skipHitlReview() {
        closeLabelStudioModal();
        removeHitlEntryCards();
        window.__hitlAnnotateDisabled = true;
        if (typeof window.disableHitlActionCards === 'function') {
            window.disableHitlActionCards('已取消人工标注');
        }
        if (typeof showToast === 'function') {
            showToast('您已取消人工标注，保留当前 AI 自动生成结果', 'info');
        }
    }

    window.syncHitlEntryCard = function syncHitlEntryCard() {
        if (_hitlCompletedMode && _lastHitlData) {
            mountHitlCompletedCard(_lastHitlData);
            return;
        }
        if (!_lastHitlData || window.__hitlAnnotateDisabled) {
            removeHitlEntryCards();
            return;
        }
        mountHitlEntryCard(_lastHitlData);
    };

    function bridgeLabelStudioSession() {
        return fetch('/api/hitl/ls-session-bridge', {
            method: 'POST',
            headers: authHeadersMerge(),
            credentials: 'include',
        }).then(function (resp) {
            if (!resp.ok) {
                return resp.json().catch(function () { return {}; }).then(function (body) {
                    var msg = (body && (body.detail || body.message)) || ('HTTP ' + resp.status);
                    throw new Error(msg);
                });
            }
            return resp.json();
        });
    }

    function openLabelStudioModal(lsUrl) {
        var modal = document.getElementById('hitl-ls-modal');
        var iframe = document.getElementById('hitl-ls-iframe');
        if (!modal || !iframe) {
            if (typeof showToast === 'function') showToast('标注弹窗未加载，请硬刷新页面后重试', 'danger');
            return;
        }
        if (!lsUrl) {
            if (typeof showToast === 'function') showToast('缺少 Label Studio 项目链接', 'warning');
            return;
        }
        var targetUrl = normalizeLsIframeUrl(lsUrl);
        bridgeLabelStudioSession()
            .catch(function (err) {
                console.warn('[HITL] LS session bridge 失败，将尝试直接打开 iframe:', err);
            })
            .finally(function () {
                iframe.src = targetUrl;
                modal.classList.add('is-open');
                modal.setAttribute('aria-hidden', 'false');
                document.body.classList.add('hitl-ls-modal-open');
                var doneBtn = document.getElementById('hitl-ls-modal-done-btn');
                if (doneBtn) {
                    doneBtn.disabled = false;
                    doneBtn.classList.remove('hitl-ls-modal__done-btn--pulse');
                }
                syncHitlModalCorpusSkipButton(_lastHitlData);
            });
    }
    window.openLabelStudioModal = openLabelStudioModal;

    function closeLabelStudioModal() {
        var modal = document.getElementById('hitl-ls-modal');
        var iframe = document.getElementById('hitl-ls-iframe');
        if (!modal) return;
        modal.classList.remove('is-open');
        modal.setAttribute('aria-hidden', 'true');
        document.body.classList.remove('hitl-ls-modal-open');
        if (iframe) iframe.src = 'about:blank';
    }

    function highlightModalDoneButton() {
        var doneBtn = document.getElementById('hitl-ls-modal-done-btn');
        if (doneBtn) doneBtn.classList.add('hitl-ls-modal__done-btn--pulse');
        if (typeof showToast === 'function') {
            showToast('检测到 Label Studio 提交动作，请点击下方「继续生成报告」', 'info');
        }
    }

    function parseSseChunk(text, handler) {
        var parts = String(text || '').split('\n\n');
        parts.forEach(function (part) {
            if (!part.trim()) return;
            var lines = part.split('\n');
            var eventType = 'message';
            var dataLines = [];
            lines.forEach(function (line) {
                if (line.indexOf('event: ') === 0) eventType = line.substring(7).trim();
                else if (line.indexOf('data: ') === 0) dataLines.push(line.substring(6));
            });
            if (!dataLines.length) return;
            var raw = dataLines.join('\n');
            var data = raw;
            try { data = JSON.parse(raw); } catch (_e) { /* keep string */ }
            handler(eventType, data);
        });
    }

    function resumeFromHitl(projectId, trigger) {
        if (_resumeInFlight) return;
        var sid = typeof currentSessionId !== 'undefined' ? currentSessionId : null;
        if (!sid) {
            if (typeof showToast === 'function') showToast('无法识别当前会话', 'warning');
            return;
        }
        _resumeInFlight = true;
        setHitlEntryCardBusy(true, '正在唤醒 HITL 并生成最终版报告…');
        var doneBtn = document.getElementById('hitl-ls-modal-done-btn');
        if (doneBtn) {
            doneBtn.disabled = true;
            doneBtn.innerHTML = '<span class="spinner-border spinner-border-sm me-1"></span> 正在生成最终版报告…';
        }
        fetch('/api/sessions/' + encodeURIComponent(sid) + '/resume_from_hitl', {
            method: 'POST',
            headers: authHeadersMerge(),
            body: JSON.stringify({
                project_id: projectId != null ? projectId : (_lastHitlData && _lastHitlData.project_id),
                trigger: trigger || (_lsSubmitDetected ? 'postmessage_hint' : 'manual_confirm'),
            }),
        })
            .then(function (resp) {
                if (!resp.ok) {
                    return resp.text().then(function (t) {
                        throw new Error(t || ('HTTP ' + resp.status));
                    });
                }
                var reader = resp.body && resp.body.getReader ? resp.body.getReader() : null;
                if (!reader) return resp.text();
                var decoder = new TextDecoder();
                var buf = '';
                function pump() {
                    return reader.read().then(function (chunk) {
                        if (chunk.done) return;
                        buf += decoder.decode(chunk.value, { stream: true });
                        var idx;
                        while ((idx = buf.indexOf('\n\n')) >= 0) {
                            var slice = buf.slice(0, idx + 2);
                            buf = buf.slice(idx + 2);
                            parseSseChunk(slice, function (eventType, data) {
                                if (typeof window.handleServerEvent === 'function') {
                                    window.handleServerEvent(eventType, data);
                                }
                            });
                        }
                        return pump();
                    });
                }
                return pump();
            })
            .then(function () {
                closeLabelStudioModal();
                var completedPayload = _lastHitlData ? Object.assign({}, _lastHitlData) : null;
                if (completedPayload) {
                    _hitlCompletedMode = true;
                    mountHitlCompletedCard(completedPayload);
                }
                if (typeof showToast === 'function') {
                    showToast('专家知识已更新，《最终版报告》已覆盖生成', 'success');
                }
                if (typeof fetchSessions === 'function') fetchSessions();
            })
            .catch(function (e) {
                if (_hitlCompletedMode && _lastHitlData) {
                    mountHitlCompletedCard(_lastHitlData);
                } else {
                    setHitlEntryCardBusy(false);
                }
                if (typeof showToast === 'function') showToast('HITL 唤醒失败: ' + e.message, 'danger');
            })
            .finally(function () {
                _resumeInFlight = false;
                if (doneBtn) {
                    doneBtn.disabled = false;
                    doneBtn.innerHTML = '<i class="bi bi-check-circle"></i> 我已完成标注，继续生成报告';
                }
            });
    }
    window.resumeFromHitl = resumeFromHitl;

    window.disableHitlActionCards = function disableHitlActionCards(_reason) {
        window.__hitlAnnotateDisabled = true;
        _hitlCompletedMode = false;
        removeHitlEntryCards();
    };

    window.renderHitlActionCard = function renderHitlActionCard(data) {
        if (!data) return;
        var payload = normalizeHitlPayload(data);
        if (!payload) return;
        _lastHitlData = payload;
        _hitlCompletedMode = false;
        window.__hitlAnnotateDisabled = false;
        window.__lastStateSnapshotHitl = payload;
        mountHitlEntryCard(payload);
    };

    window.restoreHitlFromSnapshot = function restoreHitlFromSnapshot(snap, sessionStatus) {
        if (!snap || typeof snap !== 'object') return;
        var ex = snap.execution_snapshot && typeof snap.execution_snapshot === 'object' ? snap.execution_snapshot : {};
        if (snap.hitl_resumed || ex.hitl_final || snap.hitl_final || window.__expertReportVersion === 'final') {
            var completedPayload = extractHitlFromSnapshot(snap);
            if (completedPayload) {
                mountHitlCompletedCard(completedPayload);
            }
            return;
        }
        if (snap.hitl_skipped || ex.hitl_skipped) {
            return;
        }
        var waiting = String(sessionStatus || '').toLowerCase() === 'waiting_for_hitl';
        var hitl = snap.hitl;
        var optionalHitl = !!(snap.hitl_reannotation_enabled || ex.hitl_reannotation_enabled || (hitl && hitl.reannotation_enabled));
        if (!waiting && !snap.hitl_pending && !hitl && !optionalHitl) return;
        if (!hitl && optionalHitl) {
            hitl = { status: 'hitl_required', reannotation_enabled: true };
        }
        if (!waiting && !snap.hitl_pending && !optionalHitl && !hitl) return;
        var payload = hitl && typeof hitl === 'object' ? Object.assign({ status: 'hitl_required' }, hitl) : { status: 'hitl_required' };
        window.renderHitlActionCard(payload);
    };

    function installPostMessageListener() {
        if (window.__hitlPostMessageBound) return;
        window.__hitlPostMessageBound = true;
        window.addEventListener('message', function (ev) {
            if (!ev || !ev.data) return;
            var modal = document.getElementById('hitl-ls-modal');
            if (!modal || !modal.classList.contains('is-open')) return;
            var d = ev.data;
            if (typeof d === 'string') {
                try { d = JSON.parse(d); } catch (_e) { return; }
            }
            if (!d || typeof d !== 'object') return;
            var t = String(d.type || d.event || d.action || d.name || '').toLowerCase();
            if (
                t.indexOf('annotation') >= 0
                || t.indexOf('submit') >= 0
                || t.indexOf('complete') >= 0
                || d.hitl === true
            ) {
                _lsSubmitDetected = true;
                highlightModalDoneButton();
            }
        });
    }

    function installHitlModalHandlers() {
        var modal = document.getElementById('hitl-ls-modal');
        if (!modal) return;
        if (!modal.__hitlBound) {
            modal.__hitlBound = true;
            var closeBtn = document.getElementById('hitl-ls-modal-close');
            if (closeBtn) closeBtn.addEventListener('click', closeLabelStudioModal);
            modal.addEventListener('click', function (e) {
                if (e.target === modal) closeLabelStudioModal();
            });
        }
        var doneBtn = document.getElementById('hitl-ls-modal-done-btn');
        if (doneBtn && !doneBtn.__hitlResumeBound) {
            doneBtn.__hitlResumeBound = true;
            doneBtn.addEventListener('click', function (ev) {
                ev.preventDefault();
                ev.stopPropagation();
                resumeFromHitl(_lastHitlData && _lastHitlData.project_id, 'manual_confirm');
            });
        }
    }

    function wrapHandleServerEvent() {
        var orig = window.handleServerEvent;
        if (!orig || orig.__hitlWrapped) return;
        window.handleServerEvent = function (eventType, data) {
            if (eventType === 'hitl_action') {
                if (typeof data === 'string') {
                    try { data = JSON.parse(data); } catch (_e) { /* ignore */ }
                }
                window.renderHitlActionCard(data);
            }
            if (eventType === 'step_result' && data && data.report_data && data.report_data.steps_details) {
                window.scanAndRenderHitlFromSteps(data.report_data.steps_details);
            }
            if (eventType === 'done') {
                var st = String((data && data.status) || '').toLowerCase();
                if (data && data.hitl_resumed) {
                    /* 覆写完成：completed 卡由 diagnosis / resume 回调统一挂载，避免重复 entry 卡 */
                } else if ((st === 'waiting_for_hitl' || st === 'success') && window.__lastStateSnapshotHitl && !window.__hitlAnnotateDisabled) {
                    window.renderHitlActionCard(window.__lastStateSnapshotHitl);
                }
            }
            if (eventType === 'state_snapshot' && data && (data.hitl || data.hitl_pending)) {
                if (data.hitl) window.__lastStateSnapshotHitl = data.hitl;
                if (data.hitl_pending && data.hitl) {
                    window.renderHitlActionCard(data.hitl);
                }
            }
            if (eventType === 'diagnosis' && data && data.report_data && data.report_data.hitl_final) {
                var snapHitl = window.__lastStateSnapshotHitl || _lastHitlData;
                if (snapHitl && typeof window.mountHitlCompletedCard === 'function') {
                    _hitlCompletedMode = true;
                    window.mountHitlCompletedCard(snapHitl);
                }
            }
            return orig.apply(this, arguments);
        };
        window.handleServerEvent.__hitlWrapped = true;
    }

    function initHitlUi() {
        installPostMessageListener();
        installHitlModalHandlers();
        wrapHandleServerEvent();
    }

    if (document.readyState === 'loading') {
        document.addEventListener('DOMContentLoaded', initHitlUi);
    } else {
        initHitlUi();
    }
})();
