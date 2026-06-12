/**
 * 文件树精准寻址 + VS Code 风格文本前景色状态（added/modified 仅染文件名与图标）
 */
(function () {
    'use strict';

    window.__fileTreeStatusRegistry = window.__fileTreeStatusRegistry || {};

    function normalizePathForCompare(path) {
        if (typeof window.normalizePathForCompare === 'function') {
            return window.normalizePathForCompare(path);
        }
        return String(path || '').trim().replace(/\\/g, '/').replace(/\/+$/, '').toLowerCase();
    }

    function normalizeMountPathString(path) {
        if (typeof window.normalizeMountPathString === 'function') {
            return window.normalizeMountPathString(path);
        }
        var p = String(path || '').trim();
        if (!p || p === '<PENDING_UPLOAD>' || p === '<待上传数据>' || p === '<user_input>') return '';
        if (p.length > 1 && p.charAt(0) === '<' && p.charAt(p.length - 1) === '>') return '';
        return p;
    }

    function pathBasename(path) {
        var p = String(path || '').replace(/\\/g, '/');
        var idx = p.lastIndexOf('/');
        return idx >= 0 ? p.slice(idx + 1) : p;
    }

    function extractCoreTokens(path) {
        var base = pathBasename(path);
        var stem = base.replace(/\.[^./\\]+$/, '');
        stem = stem.replace(/_\d{8}(_\d{6})?$/i, '');
        stem = stem.replace(/[-_]\d{4,}$/i, '');
        var tokens = [];
        if (stem) tokens.push(stem.toLowerCase());
        if (base) tokens.push(base.toLowerCase());
        stem.split(/[^a-zA-Z0-9\u4e00-\u9fff]+/).forEach(function (part) {
            if (part.length >= 4) tokens.push(part.toLowerCase());
        });
        var out = [];
        tokens.forEach(function (t) {
            if (out.indexOf(t) < 0) out.push(t);
        });
        return out;
    }

    function extractSessionIdFromPath(path) {
        var norm = normalizePathForCompare(path);
        var m = norm.match(/\/([0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12})(?:\/|$)/i);
        return m ? m[1] : '';
    }

    function scorePathMatch(targetPath, nodePath) {
        var normT = normalizePathForCompare(targetPath);
        var normN = normalizePathForCompare(nodePath);
        if (!normT || !normN) return 0;
        if (normT === normN) return 1000;

        var baseT = pathBasename(targetPath).toLowerCase();
        var baseN = pathBasename(nodePath).toLowerCase();
        if (baseT && baseT === baseN) return 900;

        if (normN.endsWith('/' + baseT) || normN.endsWith(baseT)) return 850;
        if (normT.endsWith('/' + baseN) || normT.endsWith(baseN)) return 840;

        var sid = extractSessionIdFromPath(targetPath);
        if (sid && normN.indexOf(sid) >= 0 && baseT && baseN === baseT) return 820;
        if (sid && normN.indexOf(sid) >= 0) return 500 + (baseT && normN.indexOf(baseT) >= 0 ? 100 : 0);

        var tokens = extractCoreTokens(targetPath);
        var tokenScore = 0;
        tokens.forEach(function (tok) {
            if (tok.length >= 4 && (normN.indexOf(tok) >= 0 || baseN.indexOf(tok) >= 0)) {
                tokenScore += 25;
            }
        });
        if (tokenScore > 0) return 400 + tokenScore;

        if (normT.length > 12 && (normN.indexOf(normT) >= 0 || normT.indexOf(normN) >= 0)) return 350;
        if (baseT.length >= 6 && baseN.indexOf(baseT) >= 0) return 300;
        if (baseT.length >= 6 && baseT.indexOf(baseN) >= 0) return 280;

        return 0;
    }

    function getFileTreeContainer() {
        return document.getElementById('file-tree-container');
    }

    function getWorkspaceRootPathNorm() {
        var sel = null;
        if (typeof window.getCachedWorkspaceSelection === 'function') {
            sel = window.getCachedWorkspaceSelection();
        }
        var root = sel && sel.workspace_path ? String(sel.workspace_path).trim() : '';
        return root ? normalizePathForCompare(root) : '';
    }

    function isWorkspaceRootTreeNode(el) {
        if (!el || !el.dataset) return false;
        var nodePath = normalizePathForCompare(el.dataset.path || '');
        var rootPath = getWorkspaceRootPathNorm();
        return !!(nodePath && rootPath && nodePath === rootPath);
    }

    function isMountCardElement(el) {
        if (!el || !el.closest) return false;
        return !!el.closest(
            '#sidebar-workspace-indicator, #right-sidebar-workspace-path, .omics-mount-card, .mount-path-card'
        );
    }

    function clearMountCardStatusClasses() {
        var selectors = [
            '#sidebar-workspace-indicator .omics-mount-card',
            '#right-sidebar-workspace-path .omics-mount-card',
            '#sidebar-workspace-indicator .mount-path-card',
            '#right-sidebar-workspace-path .mount-path-card',
        ];
        selectors.forEach(function (sel) {
            document.querySelectorAll(sel).forEach(function (el) {
                el.classList.remove(
                    'status-added',
                    'status-modified',
                    'file-tree-file--highlight',
                    'selected-node'
                );
            });
        });
    }

    function collectPathNodes(container) {
        if (!container) return [];
        return Array.prototype.slice.call(
            container.querySelectorAll('.file-tree-node[data-path], .file-tree-file[data-path]')
        );
    }

    function findFileTreeNodeByPath(targetPath, container, minScore) {
        container = container || getFileTreeContainer();
        var p = normalizeMountPathString(targetPath);
        if (!container || !p) return null;

        var nodes = collectPathNodes(container);
        var exact = null;
        var norm = normalizePathForCompare(p);
        nodes.forEach(function (el) {
            if (normalizePathForCompare(el.dataset.path || '') === norm) exact = el;
        });
        if (exact) {
            delete exact.dataset.fuzzyMatch;
            return exact;
        }

        var threshold = typeof minScore === 'number' ? minScore : 280;
        if (threshold >= 1000) return null;

        var best = null;
        var bestScore = -1;
        nodes.forEach(function (el) {
            var s = scorePathMatch(p, el.dataset.path || '');
            if (s > bestScore) {
                bestScore = s;
                best = el;
            }
        });
        if (best && bestScore >= threshold) {
            best.dataset.fuzzyMatch = '1';
            return best;
        }
        return null;
    }

    function clearAllFileTreeStatusClasses(container) {
        container = container || getFileTreeContainer();
        if (!container) return;
        container.querySelectorAll('.file-tree-node').forEach(function (node) {
            node.classList.remove('status-added', 'status-modified');
        });
        container.querySelectorAll('.file-tree-folder.file-tree-node').forEach(function (node) {
            if (isWorkspaceRootTreeNode(node)) {
                node.classList.remove('status-added', 'status-modified', 'file-tree-file--highlight');
            }
        });
        clearMountCardStatusClasses();
    }

    function expandAncestorsForNode(node) {
        if (!node) return;
        var p = node.parentElement;
        while (p) {
            if (p.tagName === 'DETAILS' || (p.classList && p.classList.contains('file-tree-folder'))) {
                p.open = true;
            }
            if (p.id === 'file-tree-container') break;
            p = p.parentElement;
        }
        var container = getFileTreeContainer();
        if (!container) return;
        var targetNorm = normalizePathForCompare(node.dataset.path || '');
        container.querySelectorAll('.file-tree-folder[data-path]').forEach(function (det) {
            var dp = normalizePathForCompare(det.dataset.path || '');
            if (dp && targetNorm && (targetNorm === dp || targetNorm.indexOf(dp + '/') === 0)) {
                det.open = true;
            }
        });
    }

    function clearSelectedFileTreeNodes(container) {
        container = container || getFileTreeContainer();
        if (!container) return;
        container.querySelectorAll('.file-tree-node.selected-node, .file-tree-file.selected-node').forEach(function (n) {
            n.classList.remove('selected-node');
        });
    }

    function focusFileTreeNode(node, options) {
        options = options || {};
        if (!node) return false;
        expandAncestorsForNode(node);
        var container = getFileTreeContainer();
        clearSelectedFileTreeNodes(container);
        node.classList.add('selected-node');
        var behavior = options.scrollBehavior || 'smooth';
        try {
            node.scrollIntoView({ behavior: behavior, block: 'center' });
        } catch (_e) {
            try {
                node.scrollIntoView(true);
            } catch (_e2) { /* ignore */ }
        }
        if (options.flashHighlight !== false) {
            node.classList.add('file-tree-file--highlight');
            window.setTimeout(function () {
                node.classList.remove('file-tree-file--highlight');
            }, options.flashMs || 3200);
        }
        return true;
    }

    function registerFileTreeStatus(filePath, status) {
        var p = normalizeMountPathString(filePath);
        if (!p) return;
        var norm = normalizePathForCompare(p);
        var st = status === 'modified' ? 'modified' : 'added';
        var prev = window.__fileTreeStatusRegistry[norm];
        if (prev && prev.status === 'modified') st = 'modified';
        window.__fileTreeStatusRegistry[norm] = { status: st, path: p };
    }

    function mergeFolderStatus(el, status) {
        if (!el || !el.classList) return;
        if (status === 'modified') {
            el.classList.remove('status-added');
            el.classList.add('status-modified');
            return;
        }
        if (!el.classList.contains('status-modified')) {
            el.classList.remove('status-added');
            el.classList.add('status-added');
        }
    }

    function bubbleStatusToFolderAncestors(node, status) {
        if (!node) return;
        var el = node.parentElement;
        while (el) {
            if (el.id === 'file-tree-container') break;
            if (isMountCardElement(el)) break;
            if (el.classList && el.classList.contains('file-tree-node')) {
                var t = el.dataset.type || '';
                if ((t === 'directory' || el.classList.contains('file-tree-folder')) && !isWorkspaceRootTreeNode(el)) {
                    mergeFolderStatus(el, status);
                }
            }
            el = el.parentElement;
        }
    }

    function applyIdeStyleHighlight(filePath, status, options) {
        options = options || {};
        var p = normalizeMountPathString(filePath);
        if (!p) return false;
        var st = status === 'modified' ? 'modified' : 'added';
        if (!options.skipRegistry) registerFileTreeStatus(p, st);

        var container = getFileTreeContainer();
        if (!container) return false;

        var node = findFileTreeNodeByPath(p, container, 1000);
        if (!node) return false;
        if (isMountCardElement(node) || isWorkspaceRootTreeNode(node)) return false;

        node.classList.remove('status-added', 'status-modified');
        node.classList.add('status-' + st);
        bubbleStatusToFolderAncestors(node, st);
        clearMountCardStatusClasses();
        return true;
    }

    function reapplyFileTreeStatusHighlights(container) {
        container = container || getFileTreeContainer();
        if (!container) return;
        clearAllFileTreeStatusClasses(container);
        clearMountCardStatusClasses();
        var registry = window.__fileTreeStatusRegistry || {};
        Object.keys(registry).forEach(function (key) {
            var entry = registry[key];
            if (entry && entry.path) {
                applyIdeStyleHighlight(entry.path, entry.status, { skipRegistry: true });
            }
        });
        clearMountCardStatusClasses();
    }

    function replaceFileTreeStatusFromBackend(detail) {
        detail = detail || {};
        window.__fileTreeStatusRegistry = {};
        var files = detail.changed_files;
        if (Array.isArray(files) && files.length) {
            files.forEach(function (item) {
                if (!item || !item.path) return;
                registerFileTreeStatus(item.path, item.status || 'added');
            });
        } else {
            (detail.changed_paths || []).forEach(function (p) {
                registerFileTreeStatus(p, 'added');
            });
        }
        (detail.sft_corpus_paths || []).forEach(function (p) {
            registerFileTreeStatus(p, 'added');
        });
        clearAllFileTreeStatusClasses();
        reapplyFileTreeStatusHighlights();
    }

    function locateAndFocusFileInTree(targetPath, options) {
        options = options || {};
        var p = normalizeMountPathString(targetPath);
        if (!p) return { found: false };
        var container = getFileTreeContainer();
        if (!container) return { found: false };

        var minScore = typeof options.minScore === 'number' ? options.minScore : 280;
        var node = findFileTreeNodeByPath(p, container, minScore);
        if (!node && options.relaxedFallback) {
            node = findFileTreeNodeByPath(p, container, 200);
        }
        if (!node) return { found: false, path: p };

        expandAncestorsForNode(node);
        focusFileTreeNode(node, options);
        return {
            found: true,
            path: p,
            node: node,
            fuzzy: node.dataset.fuzzyMatch === '1',
        };
    }

    function handleSessionFilesChangedEvent(detail) {
        replaceFileTreeStatusFromBackend(detail);
    }

    window.findFileTreeNodeByPath = findFileTreeNodeByPath;
    window.expandAncestorsForNode = expandAncestorsForNode;
    window.focusFileTreeNode = focusFileTreeNode;
    window.registerFileTreeStatus = registerFileTreeStatus;
    window.applyIdeStyleHighlight = applyIdeStyleHighlight;
    window.clearAllFileTreeStatusClasses = clearAllFileTreeStatusClasses;
    window.reapplyFileTreeStatusHighlights = reapplyFileTreeStatusHighlights;
    window.replaceFileTreeStatusFromBackend = replaceFileTreeStatusFromBackend;
    window.locateAndFocusFileInTree = locateAndFocusFileInTree;
    window.handleSessionFilesChangedHighlight = handleSessionFilesChangedEvent;
})();
