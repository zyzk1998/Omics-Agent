#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
生成 draw.io 格式流程图，用户可在 app.diagrams.net 打开后导出为 VSDX
draw.io 导出的 VSDX 通常与 WPS 兼容性更好
"""
import os
import xml.etree.ElementTree as ET
from xml.dom import minidom

def create_drawio():
    # draw.io mxfile 结构
    root = ET.Element("mxfile", {
        "host": "app.diagrams.net",
        "modified": "2025-01-01T00:00:00.000Z",
        "agent": "Omics Agent",
        "version": "22.1.0",
        "etag": "planning_execution",
        "type": "device"
    })
    diagram = ET.SubElement(root, "diagram", {"id": "flowchart", "name": "规划与执行流程"})
    mxgraph = ET.SubElement(diagram, "mxGraphModel", {
        "dx": "1422", "dy": "794", "grid": "1", "gridSize": "10",
        "guides": "1", "tooltips": "1", "connect": "1", "arrows": "1",
        "fold": "1", "page": "1", "pageScale": "1", "pageWidth": "1400",
        "pageHeight": "1000", "math": "0", "shadow": "0"
    })
    root_elem = ET.SubElement(mxgraph, "root")
    ET.SubElement(root_elem, "mxCell", {"id": "0"})
    ET.SubElement(root_elem, "mxCell", {"id": "1", "parent": "0"})

    cells = root_elem
    nodes = [
        (10, 50, 100, 100, 50, "1.用户\n自然语言→HTTP"),
        (11, 170, 100, 100, 50, "2.编排器\nRule Check"),
        (12, 290, 100, 100, 50, "2.5.文件检查器"),
        (13, 410, 100, 100, 50, "3.LLM\nIntent→JSON"),
        (14, 530, 100, 100, 50, "4.编排器\nTask Object"),
        (15, 650, 100, 100, 50, "5.逻辑分支\nhas_file?"),
        (20, 50, 250, 100, 50, "6.规划器\nPlanning Context"),
        (21, 170, 250, 100, 50, "6.5.参数推荐"),
        (22, 290, 250, 100, 50, "7.检索器\nQuery Vector"),
        (23, 410, 250, 100, 50, "8.检索器\nTool Schemas"),
        (24, 530, 250, 100, 50, "9.LLM\nWorkflow JSON"),
        (25, 650, 250, 100, 50, "10.规划器\nComplete"),
        (26, 770, 250, 100, 50, "11.前端\nRendered UI"),
        (30, 50, 450, 100, 50, "1.用户\n执行请求"),
        (31, 170, 450, 100, 50, "2.编排器\nAgent选择"),
        (32, 290, 450, 100, 50, "2.5.LLM\nIntent检测"),
        (33, 410, 450, 100, 50, "3.执行引擎\nStep Result"),
        (34, 530, 450, 100, 50, "4.工具层\nTool Output"),
        (40, 50, 600, 100, 50, "5.执行引擎\nAggregated"),
        (41, 170, 600, 100, 50, "5.5.上下文回滚\nEvidence"),
        (42, 290, 600, 100, 50, "6.领域智能体\nAI Report"),
        (43, 410, 600, 100, 50, "7.输出\nRendered UI"),
    ]

    for nid, x, y, w, h, text in nodes:
        cell = ET.SubElement(cells, "mxCell", {
            "id": str(nid), "value": text.replace("\n", "&#xa;"),
            "style": "rounded=1;whiteSpace=wrap;html=1;fillColor=#dae8fc;strokeColor=#6c8ebf;fontSize=11;",
            "vertex": "1", "parent": "1"
        })
        ET.SubElement(cell, "mxGeometry", {"x": str(x), "y": str(y), "width": str(w), "height": str(h), "as": "geometry"})

    # 连接线
    edges = [(10,11),(11,12),(12,13),(13,14),(14,15),(15,20),(20,21),(21,22),(22,23),(23,24),(24,25),(25,26),
             (30,31),(31,32),(32,33),(33,34),(34,40),(40,41),(41,42),(42,43)]
    for i, (a, b) in enumerate(edges):
        eid = 100 + i
        edge = ET.SubElement(cells, "mxCell", {
            "id": str(eid), "style": "edgeStyle=orthogonalEdgeStyle;rounded=0;html=1;exitX=1;exitY=0.5;entryX=0;entryY=0.5;",
            "edge": "1", "parent": "1", "source": str(a), "target": str(b)
        })
        ET.SubElement(edge, "mxGeometry", {"relative": "1", "as": "geometry"})

    # 标题
    title = ET.SubElement(cells, "mxCell", {
        "id": "200", "value": "Omics Agent — 规划与执行流程数据追踪",
        "style": "text;html=1;fontSize=16;fontStyle=1;align=center;verticalAlign=middle;",
        "vertex": "1", "parent": "1"
    })
    ET.SubElement(title, "mxGeometry", {"x": "400", "y": "20", "width": "400", "height": "40", "as": "geometry"})

    rough = ET.tostring(root, encoding="unicode", default_namespace="")
    reparsed = minidom.parseString(rough)
    pretty = reparsed.toprettyxml(indent="  ", encoding="utf-8").decode("utf-8")

    out = os.path.join(os.path.dirname(__file__), "..", "planning_execution_flowchart.drawio")
    with open(out, "w", encoding="utf-8") as f:
        f.write(pretty)
    print(f"已生成: {out}")
    print("请用 app.diagrams.net 打开，然后 文件→导出为→VSDX，导出的 VSDX 与 WPS 兼容性通常更好")

if __name__ == "__main__":
    create_drawio()
