/**
 * 工作台宽度调度、智能媒体/表格容器、执行结果专注模式（全屏阅读）
 * 由 index.html 在 marked 初始化后使用 buildSmartContentShell 等 API
 */
(function () {
    'use strict';

    function escAttr(s) {
        return String(s == null ? '' : s)
            .replace(/&/g, '&amp;')
            .replace(/"/g, '&quot;')
            .replace(/</g, '&lt;')
            .replace(/>/g, '&gt;');
    }

    /**
     * 包裹表格 / 图片等「易撑破布局」的内容，配套右上角媒体工具条
     */
    function buildSmartContentShell(innerHtml) {
        return (
            '<div class="smart-content-wrapper" data-zoom="1">' +
            '<div class="smart-content-inner">' +
            innerHtml +
            '</div>' +
            '<div class="media-toolbar" role="toolbar" aria-label="内容缩放与全屏">' +
            '<button type="button" class="media-toolbar-btn" data-action="zoom-in" title="放大">+</button>' +
            '<button type="button" class="media-toolbar-btn" data-action="zoom-out" title="缩小">−</button>' +
            '<button type="button" class="media-toolbar-btn" data-action="fullscreen" title="全屏">⛶</button>' +
            '</div></div>'
        );
    }

    function workspaceContentNeedsReportWidth(html) {
        if (!html || typeof html !== 'string') return false;
        // 含 Markdown 围栏代码块时，渲染后会走 buildSmartContentShell，需与表格/图同宽
        return /markdown-table|md-table|class="table|&lt;table[\s>]|```\s*mermaid|```|!\[[^\]]*\]\(|<img[\s>]/i.test(
            html
        );
    }

    function setWorkspacePanelReportMode(on) {
        var wp = document.getElementById('workspace-pane');
        if (!wp) return;
        if (on) {
            wp.classList.add('panel-report-width');
        } else {
            wp.classList.remove('panel-report-width');
        }
    }

    function syncWorkspacePanelWidthFromContent(containerOrHtml) {
        var html = '';
        if (typeof containerOrHtml === 'string') {
            html = containerOrHtml;
        } else if (containerOrHtml && containerOrHtml.nodeType === 1) {
            if (containerOrHtml.querySelector) {
                if (
                    containerOrHtml.querySelector(
                        'table.markdown-table, table.md-table, .markdown-table, .smart-content-wrapper, .report-image, img.img-fluid, img.report-image, canvas, svg'
                    )
                ) {
                    setWorkspacePanelReportMode(true);
                    return;
                }
            }
            html = containerOrHtml.innerHTML || '';
        }
        setWorkspacePanelReportMode(workspaceContentNeedsReportWidth(html));
    }

    function getSmartZoom(wrap) {
        var z = parseFloat(wrap.getAttribute('data-zoom') || '1', 10);
        if (isNaN(z) || z <= 0) z = 1;
        return z;
    }

    function setSmartZoom(wrap, z) {
        z = Math.max(0.5, Math.min(2.5, z));
        wrap.setAttribute('data-zoom', String(z));
        var inner = wrap.querySelector('.smart-content-inner');
        if (!inner) return;
        if (typeof inner.style.zoom !== 'undefined' && inner.style.zoom !== null) {
            inner.style.zoom = (z * 100).toFixed(0) + '%';
        } else {
            inner.style.transform = 'scale(' + z + ')';
            inner.style.transformOrigin = '0 0';
        }
    }

    function ensureMediaFullscreenDialog() {
        var el = document.getElementById('gibh-smart-media-fullscreen');
        if (el) return el;
        el = document.createElement('div');
        el.id = 'gibh-smart-media-fullscreen';
        el.className = 'gibh-media-fullscreen-backdrop';
        el.setAttribute('role', 'dialog');
        el.setAttribute('aria-modal', 'true');
        el.setAttribute('aria-label', '全屏查看');
        el.style.display = 'none';
        el.innerHTML =
            '<div class="gibh-media-fullscreen-panel markdown-body">' +
            '<button type="button" class="gibh-media-fullscreen-close" aria-label="关闭">×</button>' +
            '<div class="gibh-media-fullscreen-body"></div>' +
            '</div>';
        document.body.appendChild(el);
        el.addEventListener('click', function (e) {
            if (e.target === el) closeMediaFullscreenDialog();
        });
        var closeBtn = el.querySelector('.gibh-media-fullscreen-close');
        if (closeBtn) {
            closeBtn.addEventListener('click', function (e) {
                e.stopPropagation();
                closeMediaFullscreenDialog();
            });
        }
        return el;
    }

    function openMediaFullscreenFromWrapper(wrap) {
        var inner = wrap && wrap.querySelector ? wrap.querySelector('.smart-content-inner') : null;
        if (!inner) return;
        var dlg = ensureMediaFullscreenDialog();
        var body = dlg.querySelector('.gibh-media-fullscreen-body');
        if (body) {
            body.innerHTML = inner.innerHTML;
        }
        dlg.style.display = 'flex';
    }

    function closeMediaFullscreenDialog() {
        var dlg = document.getElementById('gibh-smart-media-fullscreen');
        if (dlg) dlg.style.display = 'none';
    }

    function ensureFocusModeRoot() {
        var el = document.getElementById('gibh-focus-mode-root');
        if (el) return el;
        el = document.createElement('div');
        el.id = 'gibh-focus-mode-root';
        el.className = 'focus-mode-overlay';
        el.setAttribute('role', 'dialog');
        el.setAttribute('aria-modal', 'true');
        el.setAttribute('aria-label', '专注模式');
        el.style.display = 'none';
        el.innerHTML =
            '<div class="focus-mode-backdrop" data-focus-dismiss="1"></div>' +
            '<div class="focus-mode-card">' +
            '<button type="button" class="focus-mode-close" title="关闭" aria-label="关闭">×</button>' +
            '<div class="focus-mode-body markdown-body"></div>' +
            '</div>';
        document.body.appendChild(el);
        el.addEventListener('click', function (e) {
            if (e.target && e.target.getAttribute && e.target.getAttribute('data-focus-dismiss') === '1') {
                closeFocusMode();
            }
        });
        var xb = el.querySelector('.focus-mode-close');
        if (xb) {
            xb.addEventListener('click', function (e) {
                e.stopPropagation();
                closeFocusMode();
            });
        }
        return el;
    }

    function openFocusMode(contentHtml) {
        var root = ensureFocusModeRoot();
        var body = root.querySelector('.focus-mode-body');
        if (body) {
            body.innerHTML = contentHtml || '';
        }
        root.style.display = 'flex';
        document.body.classList.add('focus-mode-open');
    }

    function closeFocusMode() {
        var root = document.getElementById('gibh-focus-mode-root');
        if (root) {
            root.style.display = 'none';
        }
        document.body.classList.remove('focus-mode-open');
    }

    function openFocusModeFromAccordion(btn) {
        var h2 = btn && btn.closest ? btn.closest('h2.accordion-header') : null;
        if (!h2) return;
        var collapse = h2.nextElementSibling;
        if (!collapse || !collapse.classList || !collapse.classList.contains('accordion-collapse')) return;
        var accBody = collapse.querySelector('.accordion-body');
        if (!accBody) return;
        openFocusMode(accBody.innerHTML);
    }

    document.addEventListener(
        'click',
        function (e) {
            var t = e.target;
            if (!t || !t.closest) return;
            var tb = t.closest('.media-toolbar-btn');
            if (tb) {
                e.preventDefault();
                var wrap = tb.closest('.smart-content-wrapper');
                if (!wrap) return;
                var act = tb.getAttribute('data-action');
                var z = getSmartZoom(wrap);
                if (act === 'zoom-in') {
                    setSmartZoom(wrap, z * 1.1);
                } else if (act === 'zoom-out') {
                    setSmartZoom(wrap, z / 1.1);
                } else if (act === 'fullscreen') {
                    openMediaFullscreenFromWrapper(wrap);
                }
                return;
            }
            var fb = t.closest('.focus-mode-btn');
            if (fb) {
                e.stopPropagation();
                e.preventDefault();
                openFocusModeFromAccordion(fb);
            }
        },
        true
    );

    document.addEventListener('keydown', function (e) {
        if (e.key === 'Escape') {
            closeFocusMode();
            closeMediaFullscreenDialog();
        }
    });

    window.buildSmartContentShell = buildSmartContentShell;
    window.syncWorkspacePanelWidthFromContent = syncWorkspacePanelWidthFromContent;
    window.setWorkspacePanelReportMode = setWorkspacePanelReportMode;
    window.openFocusMode = openFocusMode;
    window.closeFocusMode = closeFocusMode;
    window.openFocusModeFromAccordion = openFocusModeFromAccordion;
    window.escAttr = escAttr;
})();
