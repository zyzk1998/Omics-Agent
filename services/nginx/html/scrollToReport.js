window.scrollToReport = function(messageId) {
    try {
        // 优先通过 messageId 精确定位消息行
        let targetRow = null;
        if (messageId) {
            targetRow = document.querySelector('.message-row.ai[data-message-id="' + messageId + '"]');
        }
        if (!targetRow) {
            targetRow = document.querySelector('#chatContent .message-row.ai:last-of-type');
        }
        if (!targetRow) return;

        // 激活右侧工作台
        var wp = document.getElementById('workspace-pane');
        if (wp) wp.classList.add('workspace-active');
        if (typeof setWorkspaceTitle === 'function') setWorkspaceTitle('分析报告');

        // 滚动右侧报告容器到可见区域
        var reportRoot = document.getElementById('active-workspace-content') ||
                         document.querySelector('.unified-report-container:last-of-type');
        if (reportRoot && typeof reportRoot.scrollIntoView === 'function') {
            reportRoot.scrollIntoView({ behavior: 'smooth', block: 'start' });
            // 简单高亮闪烁效果，便于用户聚焦
            reportRoot.classList.add('highlight-report');
            setTimeout(function () { reportRoot.classList.remove('highlight-report'); }, 2000);
        }
    } catch (e) {
        console.warn('⚠️ [scrollToReport] 定位报告失败:', e);
    }
};

