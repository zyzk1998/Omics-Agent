#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
规划流程 + 执行流程 数据追踪 架构流程图
简约风格，避免遮挡，输出 VSDX (Visio) 格式
"""
import os
import sys

def main():
    try:
        import aspose.diagram as diagram
    except ImportError:
        print("请安装: pip install aspose-diagram-python")
        sys.exit(1)

    d = diagram.Diagram()
    page = d.pages[0]

    # 设置页面尺寸 14in x 10in
    try:
        page.page_sheet.page_width.value = 14.0
        page.page_sheet.page_height.value = 10.0
    except Exception:
        pass

    # 节点尺寸 (简约: 小矩形)
    W, H = 0.95, 0.45
    MARGIN = 0.03

    def box(x, y, text):
        """左下角 (x,y)，宽 W 高 H"""
        page.draw_rectangle(x, y, x + W, y + H)
        page.add_text(x + MARGIN, y + MARGIN, W - 2*MARGIN, H - 2*MARGIN, text)

    def arrow(x1, y1, x2, y2):
        page.draw_line(x1, y1, x2, y2)

    # ========== 规划流程 (Planning Flow) - 上半部分 ==========
    # Row 1: Y=9.0
    y1 = 9.0
    x = 0.2
    box(x, y1, "1.用户\n自然语言→HTTP POST"); x += W + 0.15
    box(x, y1, "2.编排器\nRule Check"); x += W + 0.15
    box(x, y1, "2.5.文件检查器\nFileMetadata"); x += W + 0.15
    box(x, y1, "3.LLM\nIntent→JSON"); x += W + 0.15
    box(x, y1, "4.编排器\nTask Object"); x += W + 0.15
    box(x, y1, "5.逻辑分支\nhas_file?"); x += W + 0.15
    # Row 2: Y=8.2
    y2 = 8.2
    x = 0.2
    box(x, y2, "6.规划器\nPlanning Context"); x += W + 0.15
    box(x, y2, "6.5.参数推荐\nParam Tuner"); x += W + 0.15
    box(x, y2, "7.检索器\nQuery Vector"); x += W + 0.15
    box(x, y2, "8.检索器\nTool Schemas"); x += W + 0.15
    box(x, y2, "9.LLM\nWorkflow JSON"); x += W + 0.15
    box(x, y2, "10.规划器\nComplete Workflow"); x += W + 0.15
    box(x, y2, "11.前端\nRendered UI"); x += W + 0.15

    # 规划流程连接线 (从左到右)
    def right_x(i): return 0.2 + i * (W + 0.15) + W
    def left_x(i): return 0.2 + i * (W + 0.15)
    cy = H / 2
    # Row1 内部 (0→1→2→3→4→5)
    for i in range(5):
        arrow(right_x(i), y1 + cy, left_x(i+1), y1 + cy)
    # Row1末→Row2首 (5→0)
    arrow(right_x(5), y1 + cy, left_x(0), y2 + cy)
    # Row2 内部 (0→1→2→3→4→5→6)
    for i in range(6):
        arrow(right_x(i), y2 + cy, left_x(i+1), y2 + cy)

    # ========== 执行流程 (Execution Flow) - 下半部分 ==========
    # Row 3: Y=6.5 (步骤 1-4)
    y3 = 6.5
    x = 0.2
    box(x, y3, "1.用户\n执行请求"); x += W + 0.15
    box(x, y3, "2.编排器\nAgent选择"); x += W + 0.15
    box(x, y3, "2.5.LLM\nIntent检测"); x += W + 0.15
    box(x, y3, "3.执行引擎\nStep Result"); x += W + 0.15
    box(x, y3, "4.工具层\nTool Output"); x += W + 0.15
    # Row 4: Y=5.7 (步骤 5-7)
    y4 = 5.7
    x = 0.2
    box(x, y4, "5.执行引擎\nAggregated"); x += W + 0.15
    box(x, y4, "5.5.上下文回滚\nEvidence"); x += W + 0.15
    box(x, y4, "6.领域智能体\nAI Report"); x += W + 0.15
    box(x, y4, "7.输出\nRendered UI"); x += W + 0.15

    # 执行流程连接线
    for i in range(4):
        arrow(right_x(i), y3 + cy, left_x(i+1), y3 + cy)
    arrow(right_x(4), y3 + cy, left_x(0), y4 + cy)  # Row3末→Row4首
    for i in range(3):
        arrow(right_x(i), y4 + cy, left_x(i+1), y4 + cy)

    # ========== 标题区 ==========
    page.add_text(0.2, 9.85, 5.0, 0.35, "Omics Agent — 项目架构流程图")
    page.add_text(0.2, 9.48, 4.0, 0.3, "规划流程数据追踪 (Planning Flow Data Trace)")
    page.add_text(0.2, 6.88, 4.0, 0.3, "执行流程数据追踪 (Execution Flow Data Trace)")

    base = os.path.join(os.path.dirname(__file__), "..")
    out_vsdx = os.path.join(base, "planning_execution_flowchart.vsdx")
    out_vsd = os.path.join(base, "planning_execution_flowchart.vsd")

    d.save(out_vsdx, diagram.SaveFileFormat.VSDX)
    d.save(out_vsd, diagram.SaveFileFormat.VSD)

    print(f"已生成: {out_vsdx}")
    print(f"已生成: {out_vsd}")
    print("\n若 WPS 导入 VSDX 出现「解析出错」，请尝试:")
    print("  1. 导入 .vsd 文件（兼容性通常更好）")
    print("  2. 或用 Microsoft Visio 打开 VSDX 后另存为 VSD，再导入 WPS")

if __name__ == "__main__":
    main()
