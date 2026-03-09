#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
将 Omics Agent 流程图转为 VSDX (Visio) 格式
使用 Aspose.Diagram 生成，可直接用 Visio/WPS 打开编辑
"""
import os
import sys

def px2in(px):
    return px / 100.0

def main():
    try:
        import aspose.diagram as diagram
    except ImportError:
        print("请安装: pip install aspose-diagram-python")
        sys.exit(1)

    d = diagram.Diagram()
    page = d.pages[0]

    # Visio 坐标: 10in x 6.2in, Y 向上
    # 模块: (名称, left, bottom, right, top, 文本)
    # 从 SVG 转换: left=px/100, bottom=(620-py-height)/100, right=(px+w)/100, top=(620-py)/100
    def rect(px, py, w, h, text):
        left = px2in(px)
        top = px2in(620 - py)
        right = px2in(px + w)
        bottom = px2in(620 - py - h)
        return left, bottom, right, top, text

    modules = [
        rect(30, 80, 940, 90, "用户层\n• 网页界面\n• 文件上传\n• 自然语言查询\n\n提供交互入口，支持自然语言查询与文件上传"),
        rect(30, 195, 200, 120, "编排器\n• 意图识别\n• 会话状态\n• 异步调度\n\n解析用户意图，调度后续流程"),
        rect(250, 195, 200, 120, "规划器\n• 推理链（逐步分解任务）\n• 动态规划\n• 错误分析\n\n基于推理链进行任务分解与动态规划"),
        rect(470, 195, 200, 120, "检索器\n• 语义检索\n• 工具匹配\n• 向量嵌入（可检索格式）\n\n语义匹配工具与文档，支撑决策"),
        rect(690, 195, 200, 120, "执行器\n• 动态调用\n• 数据校验\n• 自愈循环\n\n动态调用工具并校验数据，实现自愈"),
        rect(810, 295, 170, 140, "DeepSeek-R1 引擎\n• 混合专家架构\n• 上下文管理\n• 推理与逻辑\n• 大语言模型（通用AI推理）\n\n大语言模型核心，负责推理与上下文管理"),
        rect(810, 455, 170, 95, "知识库\n• 向量存储\n• 工具注册表\n\n存储向量与工具注册表，支持语义检索"),
        rect(30, 410, 760, 90, "输出层\n• 服务端推送\n• 可视化\n\n推送结果至前端，支持可视化展示"),
    ]

    for left, bottom, right, top, text in modules:
        page.draw_rectangle(left, bottom, right, top)
        # 在矩形内添加文本（左下角对齐，留边距）
        margin = 0.05
        tw = (right - left) - 2 * margin
        th = (top - bottom) - 2 * margin
        page.add_text(left + margin, bottom + margin, tw, th, text)

    # 连接线 (简化: 主要流程)
    # 用户层底 -> 编排器顶: (5, 1.7) to (1.3, 3.15)
    page.draw_line(5.0, 1.7, 1.3, 3.15)
    page.draw_line(2.3, 3.65, 2.5, 3.65)  # 编排器->规划器
    page.draw_line(4.5, 3.65, 4.7, 3.65)  # 规划器->检索器
    page.draw_line(6.7, 3.65, 6.9, 3.65)  # 检索器->执行器
    page.draw_line(7.9, 3.65, 7.9, 2.55)  # 执行器向下
    page.draw_line(7.9, 2.55, 4.1, 2.55)  # 到输出层
    page.draw_line(4.1, 2.55, 4.1, 2.1)   # 输出层

    # 标题
    page.add_text(3.5, 5.8, 2.4, 0.3, "Omics Agent — 系统拓扑蓝图")
    page.draw_rectangle(3.8, 5.5, 6.2, 5.82)
    page.add_text(5.0, 5.6, 2.2, 0.2, "国家重点研发计划")

    out = os.path.join(os.path.dirname(__file__), "..", "omics_agent_flowchart.vsdx")
    d.save(out, diagram.SaveFileFormat.VSDX)
    print(f"已生成: {out}")
    print("可用 Microsoft Visio 或 WPS 打开编辑")

if __name__ == "__main__":
    main()
