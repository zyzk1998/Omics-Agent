#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
将原 Mermaid 占位图全部渲染为 PNG（Matplotlib），写入 述职PPT/assets/。
与 Mermaid 源等价的 docker 拓扑见 ../mermaid/docker_compose_topology.mmd。
用法: cd 述职PPT && python3 scripts/render_deck_raster_diagrams.py
"""
from __future__ import annotations

import os
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import patheffects

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "assets"


def _pick_font() -> str:
    for name in ("Noto Sans CJK SC", "Noto Sans CJK JP", "Source Han Sans SC", "WenQuanYi Zen Hei"):
        try:
            from matplotlib import font_manager

            if font_manager.findfont(name, fallback_to_default=False):
                return name
        except Exception:
            continue
    return "DejaVu Sans"


def setup() -> None:
    os.makedirs(OUT, exist_ok=True)
    plt.rcParams["font.sans-serif"] = [_pick_font(), "DejaVu Sans"]
    plt.rcParams["axes.unicode_minus"] = False


def box(ax, xy, w, h, text, fc, ec="#a78bfa", lw=1.5):
    r = mpatches.FancyBboxPatch(
        xy, w, h, boxstyle="round,pad=0.02,rounding_size=0.02",
        facecolor=fc, edgecolor=ec, linewidth=lw,
    )
    ax.add_patch(r)
    ax.text(xy[0] + w / 2, xy[1] + h / 2, text, ha="center", va="center", fontsize=9, color="#f8fafc", wrap=True)


def arrow(ax, p0, p1, color="#a78bfa"):
    ax.annotate(
        "", xy=p1, xytext=p0,
        arrowprops=dict(arrowstyle="-|>", color=color, lw=1.8, mutation_scale=12),
    )


def render_gantt() -> None:
    fig, ax = plt.subplots(figsize=(12, 4.2), dpi=150)
    ax.set_facecolor("#0f172a")
    fig.patch.set_facecolor("#020617")
    rows = [
        ("探索 9-10 月", 0, 61, "#7c3aed"),
        ("架构 11-12 月", 61, 61, "#6366f1"),
        ("契约 1-4 月", 122, 121, "#38bdf8"),
        ("持续优化 4 月起", 243, 180, "#22c55e"),
    ]
    y0 = 3
    for label, start, width, c in rows:
        ax.barh(y0, width, left=start, height=0.65, color=c, edgecolor="#e2e8f0", linewidth=0.6)
        ax.text(start + 2, y0, label, va="center", fontsize=10, color="#0f172a", fontweight="bold")
        y0 -= 1
    ax.set_xlim(0, 430)
    ax.set_ylim(-0.5, 3.8)
    ax.set_yticks([])
    ax.set_xlabel("相对天数（示意甘特，与 index 叙事月份对齐）", color="#94a3b8", fontsize=10)
    ax.set_title("架构演进时间轴（栅格图，替代 Mermaid gantt）", color="#f8fafc", fontsize=13, fontweight="bold", pad=10)
    for spine in ax.spines.values():
        spine.set_color("#334155")
    ax.tick_params(colors="#94a3b8")
    fig.tight_layout()
    fig.savefig(OUT / "img_slide_gantt.png", facecolor=fig.get_facecolor())
    plt.close(fig)


def render_p1_mini() -> None:
    fig, ax = plt.subplots(figsize=(11, 2.8), dpi=150)
    ax.set_facecolor("#020617")
    fig.patch.set_facecolor("#020617")
    ax.axis("off")
    xs = [0, 2.2, 4.4, 6.6, 8.8, 11]
    labels = ["用户", "bash/curl", "Ollama", "tools.json", "HTTP 工具", "回复"]
    for i, lb in enumerate(labels):
        box(ax, (xs[i], 0.35), 1.8, 0.9, lb, "#1e293b")
    for i in range(len(xs) - 1):
        arrow(ax, (xs[i] + 1.8, 0.8), (xs[i + 1], 0.8))
    ax.set_xlim(-0.2, 12.5)
    ax.set_ylim(0, 1.6)
    ax.set_title("agent-mini 闭环（替代 Mermaid flowchart LR）", color="#f8fafc", fontsize=12, fontweight="bold", loc="left", x=0)
    fig.tight_layout()
    fig.savefig(OUT / "img_slide_p1_flow.png", facecolor=fig.get_facecolor())
    plt.close(fig)


def render_p2_route() -> None:
    fig, ax = plt.subplots(figsize=(8.5, 5), dpi=150)
    ax.set_facecolor("#020617")
    fig.patch.set_facecolor("#020617")
    ax.axis("off")
    box(ax, (3.1, 3.0), 2.2, 0.75, "用户请求", "#334155")
    box(ax, (2.8, 1.85), 2.8, 0.85, "语义路由", "#4c1d95", ec="#c4b5fd", lw=2)
    for i, (t, y) in enumerate([("追问用户", 0.95), ("聊天通道", 0.55), ("快车道", 0.15), ("任务与算力", -0.35)]):
        box(ax, (0.2 + (i % 2) * 4.8, y), 2.1, 0.55, t, "#1e293b")
    arrow(ax, (4.2, 3.0), (4.2, 2.7))
    arrow(ax, (3.5, 1.85), (1.5, 1.2))
    arrow(ax, (4.2, 1.85), (4.2, 1.15))
    arrow(ax, (4.9, 1.85), (6.5, 1.15))
    arrow(ax, (4.2, 1.85), (4.2, 0.25))
    ax.set_title("关键词阶段之后：结构化路由分流（替代 Mermaid flowchart TD）", color="#f8fafc", fontsize=11, fontweight="bold", x=0.02, y=0.98, loc="left", transform=fig.transFigure)
    fig.tight_layout(rect=[0, 0.05, 1, 0.92])
    fig.savefig(OUT / "img_slide_p2_flow.png", facecolor=fig.get_facecolor())
    plt.close(fig)


def render_sequence() -> None:
    fig, ax = plt.subplots(figsize=(11, 3.6), dpi=150)
    ax.set_facecolor("#020617")
    fig.patch.set_facecolor("#020617")
    ax.axis("off")
    cols = [("用户", 1), ("编排器", 4), ("语义路由", 7)]
    y_lifeline = 0.35
    h_line = 2.0
    for name, x in cols:
        ax.text(x, 2.85, name, ha="center", fontsize=11, color="#e9d5ff", fontweight="bold")
        ax.plot([x, x], [y_lifeline, y_lifeline + h_line], color="#64748b", lw=1.2, linestyle="--")
    # messages
    ax.annotate("", xy=(4, 2.35), xytext=(1, 2.35), arrowprops=dict(arrowstyle="-|>", color="#a78bfa", lw=2))
    ax.text(2.5, 2.42, "问题与文件、算力状态", ha="center", fontsize=9, color="#cbd5e1")
    ax.annotate("", xy=(7, 1.85), xytext=(4, 1.85), arrowprops=dict(arrowstyle="-|>", color="#a78bfa", lw=2))
    ax.text(5.5, 1.92, "RouterInput", ha="center", fontsize=9, color="#cbd5e1")
    ax.annotate("", xy=(4, 1.35), xytext=(7, 1.35), arrowprops=dict(arrowstyle="-|>", color="#34d399", lw=2, connectionstyle="arc3,rad=0.1"))
    ax.text(5.5, 1.22, "RouterOutput", ha="center", fontsize=9, color="#cbd5e1")
    ax.text(4, 0.75, "任务链 / 轻量对话", ha="center", fontsize=9, color="#94a3b8")
    ax.set_title("语义路由时序（替代 Mermaid sequenceDiagram）", color="#f8fafc", fontsize=12, fontweight="bold", loc="left", x=0.02, y=0.98, transform=fig.transFigure)
    fig.tight_layout(rect=[0, 0, 1, 0.9])
    fig.savefig(OUT / "img_slide_route_sequence.png", facecolor=fig.get_facecolor())
    plt.close(fig)


def render_tools_chain() -> None:
    fig, ax = plt.subplots(figsize=(10, 2.4), dpi=150)
    ax.set_facecolor("#020617")
    fig.patch.set_facecolor("#020617")
    ax.axis("off")
    parts = ["gibh_agent/tools", "ToolRegistry", "JSON Schema", "digest", "tools 入模型"]
    x = 0.3
    for p in parts:
        box(ax, (x, 0.4), 1.55, 0.65, p, "#312e81")
        x += 1.75
    for i in range(len(parts) - 1):
        arrow(ax, (0.3 + (i + 1) * 1.75 - 0.2, 0.72), (0.3 + i * 1.75 + 1.55 + 0.05, 0.72))
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 1.5)
    ax.set_title("注册表到模型（替代 Mermaid flowchart LR）", color="#f8fafc", fontsize=11, fontweight="bold", loc="left", x=0)
    fig.tight_layout()
    fig.savefig(OUT / "img_slide_tools_flow.png", facecolor=fig.get_facecolor())
    plt.close(fig)


def render_arch_mesh() -> None:
    fig, ax = plt.subplots(figsize=(11, 6.2), dpi=150)
    ax.set_facecolor("#020617")
    fig.patch.set_facecolor("#020617")
    ax.axis("off")
    zones = [
        (0.15, 4.8, 1.1, 0.7, "Nginx", "#1e3a5f"),
        (0.15, 3.5, 2.4, 1.0, "api-server", "#4c1d95"),
        (2.85, 3.5, 2.3, 1.0, "Orchestrator", "#5b21b6"),
        (0.15, 2.0, 2.4, 1.1, "MySQL", "#14532d"),
        (2.85, 2.0, 2.3, 0.5, "Redis", "#9a3412"),
        (2.85, 2.55, 2.3, 0.55, "uploads/results", "#0f766e"),
        (5.45, 3.5, 2.5, 1.0, "worker-pyskills", "#7f1d1d"),
        (5.45, 2.0, 2.5, 1.1, "MCP Gateway", "#6d28d9"),
        (8.2, 3.5, 2.5, 1.0, "chem 子进程", "#831843"),
        (8.2, 2.0, 2.5, 1.1, "LLM API", "#0e7490"),
    ]
    for x, y, w, h, t, c in zones:
        box(ax, (x, y), w, h, t, c, ec="#94a3b8", lw=1)
    ax.set_title("服务网格示意（替代 Mermaid flowchart TB，与 compose 一致）", color="#f8fafc", fontsize=12, fontweight="bold", x=0.02, y=0.98, transform=fig.transFigure, loc="left")
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(OUT / "img_slide_arch_mesh.png", facecolor=fig.get_facecolor())
    plt.close(fig)


def render_docker_topology() -> None:
    """与 docker-compose.yml：nginx → api-server → mysql/redis/worker/mcp-gateway/worker-pyskills + 外部 LLM/HPC"""
    fig, ax = plt.subplots(figsize=(12, 7), dpi=150)
    ax.set_facecolor("#020617")
    fig.patch.set_facecolor("#020617")
    ax.axis("off")
    # dashed compose boundary
    boundary = mpatches.Rectangle((0.1, 0.35), 10.8, 5.5, fill=False, edgecolor="#7c3aed", linewidth=2, linestyle="--")
    ax.add_patch(boundary)
    ax.text(5.5, 6.05, "Docker Compose · gibh-network", ha="center", fontsize=12, color="#e9d5ff", fontweight="bold")

    box(ax, (4.2, 5.35), 2.0, 0.55, "用户 / 浏览器", "#334155")
    box(ax, (4.0, 4.55), 2.4, 0.6, "Nginx（可选 80/8018）", "#1e3a5f")
    box(ax, (3.3, 3.35), 3.6, 0.75, "api-server :8028\n(Gunicorn + GPU 可选)", "#4c1d95", lw=2)
    box(ax, (0.35, 1.85), 2.0, 0.65, "MySQL :3306", "#14532d")
    box(ax, (2.55, 1.85), 1.85, 0.65, "Redis", "#9a3412")
    box(ax, (4.55, 1.85), 2.0, 0.65, "Celery worker", "#5b21b6")
    box(ax, (6.75, 1.85), 1.85, 0.65, "worker-pyskills", "#7f1d1d")
    box(ax, (8.75, 1.85), 1.85, 0.65, "mcp-gateway :8002", "#6d28d9")
    box(ax, (3.5, 0.55), 4.5, 0.7, "uploads / results / data（宿主机挂载）", "#0f766e")

    # external
    box(ax, (0.25, 3.35), 2.4, 0.75, "LLM 云 API\n(DeepSeek 等)", "#0e7490", ec="#5eead4")
    box(ax, (8.95, 3.35), 1.85, 0.75, "HPC MCP\n宿主机", "#78350f", ec="#fdba74")

    arrow(ax, (5.2, 5.35), (5.2, 5.15))
    arrow(ax, (5.2, 4.55), (5.2, 4.1))
    arrow(ax, (5.2, 3.35), (2.6, 2.5))
    arrow(ax, (5.2, 3.35), (3.5, 2.5))
    arrow(ax, (5.2, 3.35), (5.5, 2.5))
    arrow(ax, (6.5, 3.35), (7.6, 2.5))
    arrow(ax, (7.8, 3.35), (9.6, 2.5))
    arrow(ax, (3.3, 3.72), (1.45, 4.1))
    arrow(ax, (6.9, 3.72), (9.4, 4.1))
    arrow(ax, (5.2, 1.85), (5.5, 1.25))

    ax.text(
        5.5, 0.15,
        "说明：api-server 经 MCP_GATEWAY_URL 调 mcp-gateway；网关转发 HPC_MCP_URL；Celery 与 API 共享代码卷与 Redis。",
        ha="center", fontsize=9, color="#94a3b8",
    )
    fig.tight_layout(rect=[0, 0.02, 1, 0.98])
    fig.savefig(ROOT / "deck_docker_topology.png", facecolor=fig.get_facecolor())
    plt.close(fig)


def main() -> None:
    setup()
    render_gantt()
    render_p1_mini()
    render_p2_route()
    render_sequence()
    render_tools_chain()
    render_arch_mesh()
    render_docker_topology()
    print("OK ->", OUT, "and", ROOT / "deck_docker_topology.png")


if __name__ == "__main__":
    main()
