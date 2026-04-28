#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
述职PPT 配图：从 data/deck_metrics.json 读数，输出到 assets/（PNG，高对比、大字轴标签）。
用法（仓库内）:
  cd 述职PPT && python3 scripts/build_deck_charts.py
依赖: matplotlib
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib import patheffects

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "deck_metrics.json"
OUT = ROOT / "assets"
FONT_CANDIDATES = [
    "Noto Sans CJK SC",
    "Noto Sans CJK JP",
    "Source Han Sans SC",
    "WenQuanYi Zen Hei",
    "AR PL UMing CN",
]


def _pick_font() -> str:
    for name in FONT_CANDIDATES:
        try:
            from matplotlib import font_manager

            if font_manager.findfont(name, fallback_to_default=False):
                return name
        except Exception:
            continue
    return "DejaVu Sans"


def _style_dark(ax: plt.Axes, title: str) -> None:
    ax.set_facecolor("#0f172a")
    ax.tick_params(colors="#e2e8f0", labelsize=12)
    ax.title.set_color("#f8fafc")
    ax.title.set_fontsize(15)
    ax.title.set_fontweight("bold")
    ax.title.set_text(title)
    for spine in ax.spines.values():
        spine.set_color("#475569")
    ax.yaxis.label.set_color("#cbd5e1")
    ax.xaxis.label.set_color("#cbd5e1")


def chart_tools_bar(data: dict) -> None:
    labels = ["本地注册\nToolRegistry", "MCP 固定清单\n(JSON 条数)", "对外能力合计\n(上两项相加)"]
    values = [
        data["tool_registry_count"],
        data["mcp_catalog_count"],
        data["tools_total"],
    ]
    fig, ax = plt.subplots(figsize=(11, 5.5), dpi=150)
    colors = ["#7c3aed", "#a78bfa", "#38bdf8"]
    bars = ax.bar(labels, values, color=colors, edgecolor="#c4b5fd", linewidth=1.2)
    _style_dark(ax, "目的：一眼看清「本地工具」和「超算侧清单」各多少、加起来多少")
    ax.set_ylabel("条数")
    ax.set_ylim(0, max(values) * 1.15)
    for b, v in zip(bars, values):
        ax.text(
            b.get_x() + b.get_width() / 2,
            b.get_height() + 2,
            str(v),
            ha="center",
            va="bottom",
            color="#f8fafc",
            fontsize=16,
            fontweight="bold",
        )
    fig.patch.set_facecolor("#020617")
    fig.tight_layout()
    fig.savefig(OUT / "chart_tools_stack.png", facecolor=fig.get_facecolor())
    plt.close(fig)


def chart_eval_distribution(data: dict) -> None:
    dist = data["eval_route_distribution"]
    labels = list(dist.keys())
    values = list(dist.values())
    fig, ax = plt.subplots(figsize=(11, 5.5), dpi=150)
    bars = ax.bar(labels, values, color="#6366f1", edgecolor="#a5b4fc", linewidth=1.0)
    _style_dark(
        ax,
        "目的：19 条脱机用例里，期望走哪条路由 —— 不是准确率，是「考卷长什么样」",
    )
    ax.set_ylabel("用例条数")
    ax.set_ylim(0, max(values) * 1.25)
    ax.tick_params(axis="x", labelrotation=15)
    for b, v in zip(bars, values):
        ax.text(
            b.get_x() + b.get_width() / 2,
            b.get_height() + 0.15,
            str(v),
            ha="center",
            va="bottom",
            color="#f8fafc",
            fontsize=13,
            fontweight="bold",
        )
    fig.patch.set_facecolor("#020617")
    fig.tight_layout()
    fig.savefig(OUT / "chart_eval_distribution.png", facecolor=fig.get_facecolor())
    plt.close(fig)


def chart_pass_rate(data: dict) -> None:
    pct = float(data["eval_semantic_pass_rate_pct"])
    fig, ax = plt.subplots(figsize=(9, 3.8), dpi=150)
    ax.set_facecolor("#0f172a")
    fig.patch.set_facecolor("#020617")
    ax.barh([0], [pct], height=0.45, color="#22c55e", edgecolor="#86efac", linewidth=1.2)
    ax.set_xlim(0, 100)
    ax.set_yticks([])
    ax.set_xlabel("百分比（%）")
    ax.axvline(95, color="#f97316", linestyle="--", linewidth=1.5, label="95% 参考线")
    ax.text(
        pct + 1.5,
        0,
        f"{pct:.1f}%",
        va="center",
        color="#f8fafc",
        fontsize=18,
        fontweight="bold",
    )
    ax.set_title(
        "语义路由脱机评测：最近一次批跑的 PASS 率（请用内网 eval 输出覆盖 JSON）",
        color="#f8fafc",
        fontsize=13,
        fontweight="bold",
        pad=12,
    )
    ax.tick_params(colors="#e2e8f0", labelsize=11)
    for spine in ax.spines.values():
        spine.set_color("#475569")
    ax.legend(loc="lower right", facecolor="#1e293b", edgecolor="#64748b", labelcolor="#e2e8f0")
    fig.tight_layout()
    fig.savefig(OUT / "chart_eval_pass_rate.png", facecolor=fig.get_facecolor())
    plt.close(fig)


def diagram_flywheel() -> None:
    fig, ax = plt.subplots(figsize=(7.5, 7.5), dpi=150)
    ax.set_facecolor("#020617")
    fig.patch.set_facecolor("#020617")
    ax.set_xlim(-1.2, 1.2)
    ax.set_ylim(-1.2, 1.2)
    ax.axis("off")
    title = "数据飞轮（示意）"
    ax.text(0, 1.05, title, ha="center", va="center", fontsize=16, fontweight="bold", color="#f8fafc")

    nodes = [
        (0, 0.55, "线上真实调用\n（日志可脱敏）"),
        (-0.65, -0.45, "难例回流\n（失败/歧义）"),
        (0.65, -0.45, "路由与 digest\n迭代"),
    ]
    for x, y, t in nodes:
        circ = plt.Circle((x, y), 0.38, color="#1e293b", ec="#7c3aed", lw=2)
        ax.add_patch(circ)
        ax.text(
            x,
            y,
            t,
            ha="center",
            va="center",
            fontsize=11,
            color="#e2e8f0",
            linespacing=1.25,
        )

    arr_kw = dict(arrowstyle="-|>", color="#a78bfa", mutation_scale=18, lw=2)
    ax.annotate("", xy=(-0.45, -0.25), xytext=(0.25, 0.35), arrowprops=arr_kw)
    ax.annotate("", xy=(0.45, -0.25), xytext=(-0.25, 0.35), arrowprops=arr_kw)
    ax.annotate("", xy=(0, 0.2), xytext=(-0.55, -0.35), arrowprops=arr_kw)
    ax.text(
        0,
        -1.05,
        "闭环：用得越多 → 越容易发现错路由 → 改 prompt/清单 → 下次更省算力、少返工",
        ha="center",
        fontsize=11,
        color="#94a3b8",
    )
    fig.tight_layout()
    fig.savefig(OUT / "diagram_flywheel.png", facecolor=fig.get_facecolor())
    plt.close(fig)


def diagram_llm_tools(data: dict) -> None:
    fig, ax = plt.subplots(figsize=(11, 4.2), dpi=150)
    ax.set_facecolor("#020617")
    fig.patch.set_facecolor("#020617")
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 3)
    ax.axis("off")

    def box(x, y, w, h, text, fc, ec):
        r = plt.Rectangle((x, y), w, h, facecolor=fc, edgecolor=ec, linewidth=2)
        ax.add_patch(r)
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=11, color="#f8fafc")

    box(0.3, 1.0, 1.8, 1.0, "用户 / 前端", "#334155", "#94a3b8")
    box(2.5, 1.0, 2.0, 1.0, "编排器\nAgentOrchestrator", "#4c1d95", "#c4b5fd")
    box(5.0, 1.85, 1.9, 0.75, "大模型\n函数调用", "#0f766e", "#5eead4")
    box(5.0, 0.35, 1.9, 0.75, f"ToolRegistry\n{data['tool_registry_count']} 条", "#7c3aed", "#ddd6fe")
    box(7.5, 1.0, 2.0, 1.0, "实际执行\n子进程 / Worker / MCP", "#991b1b", "#fecaca")

    ax.annotate(
        "",
        xy=(2.5, 1.5),
        xytext=(2.1, 1.5),
        arrowprops=dict(arrowstyle="-|>", color="#a78bfa", lw=2, mutation_scale=15),
    )
    ax.annotate(
        "",
        xy=(5.0, 1.5),
        xytext=(4.5, 1.5),
        arrowprops=dict(arrowstyle="-|>", color="#a78bfa", lw=2, mutation_scale=15),
    )
    ax.annotate(
        "",
        xy=(5.95, 1.85),
        xytext=(5.95, 1.75),
        arrowprops=dict(arrowstyle="-|>", color="#5eead4", lw=1.8, mutation_scale=12),
    )
    ax.annotate(
        "",
        xy=(5.95, 1.35),
        xytext=(5.95, 1.25),
        arrowprops=dict(arrowstyle="-|>", color="#5eead4", lw=1.8, mutation_scale=12),
    )
    ax.annotate(
        "",
        xy=(7.5, 1.5),
        xytext=(6.9, 1.5),
        arrowprops=dict(arrowstyle="-|>", color="#a78bfa", lw=2, mutation_scale=15),
    )
    ax.text(
        5.0,
        2.75,
        "大模型不直接摸磁盘：编排器把「工具 JSON Schema」从注册表取出来再喂给模型；执行时仍走注册表绑定的真实函数。",
        ha="center",
        fontsize=11,
        color="#cbd5e1",
        wrap=True,
    )
    fig.tight_layout()
    fig.savefig(OUT / "diagram_llm_toolregistry.png", facecolor=fig.get_facecolor())
    plt.close(fig)


def diagram_mcp_gateway() -> None:
    fig, ax = plt.subplots(figsize=(12, 5.0), dpi=150)
    ax.set_facecolor("#020617")
    fig.patch.set_facecolor("#020617")
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 4.2)
    ax.axis("off")

    def box(x, y, w, h, text, fc, ec):
        r = plt.Rectangle((x, y), w, h, facecolor=fc, edgecolor=ec, linewidth=2)
        ax.add_patch(r)
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=10, color="#f8fafc")

    box(0.2, 1.2, 1.3, 1.0, "浏览器\nSSE", "#334155", "#94a3b8")
    box(1.8, 1.2, 1.4, 1.0, "api-server", "#1e3a5f", "#7dd3fc")
    box(3.5, 1.2, 1.6, 1.0, "Orchestrator\n（已 MCP agent 化）", "#4c1d95", "#c4b5fd")
    box(5.5, 1.2, 1.5, 1.0, "MCP Client\n动态注册", "#5b21b6", "#ddd6fe")
    box(7.3, 1.2, 1.4, 1.0, "MCP Gateway\nHPC 入口", "#6d28d9", "#e9d5ff")
    box(9.0, 1.2, 2.6, 1.0, "Slurm / SSH\n作业与结果", "#7f1d1d", "#fecaca")

    arr = dict(arrowstyle="-|>", color="#a78bfa", lw=2, mutation_scale=14)
    for i, (x0, x1) in enumerate([(1.5, 1.8), (3.2, 3.5), (5.1, 5.5), (7.0, 7.3), (8.7, 9.0)]):
        ax.annotate("", xy=(x1, 1.7), xytext=(x0, 1.7), arrowprops=arr)

    box(3.5, 2.85, 5.4, 0.85, "现状：Job ID / 状态走编排内存与日志；禁止盲轮询，用事件或回调拉结果", "#0f172a", "#f97316")
    box(3.2, 0.15, 6.0, 0.85, "规划：作业元数据、配额、失败重试写入持久化层（如 DB），与网关会话对齐", "#0f172a", "#22c55e")

    ax.text(
        6.0,
        3.95,
        "MCP × HPC：自然语言 → 编排 → 网关 → 超算节点（异步 Job）",
        ha="center",
        fontsize=13,
        fontweight="bold",
        color="#f8fafc",
        path_effects=[patheffects.withStroke(linewidth=3, foreground="#0f172a")],
    )
    fig.tight_layout()
    fig.savefig(OUT / "diagram_mcp_gateway.png", facecolor=fig.get_facecolor())
    plt.close(fig)


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    plt.rcParams["font.sans-serif"] = [_pick_font()]
    plt.rcParams["axes.unicode_minus"] = False

    with open(DATA, encoding="utf-8") as f:
        raw = f.read()
        if raw.strip().startswith("{"):
            data = json.loads(raw)
        else:
            raise SystemExit("deck_metrics.json 格式错误")

    chart_tools_bar(data)
    chart_eval_distribution(data)
    chart_pass_rate(data)
    diagram_flywheel()
    diagram_llm_tools(data)
    diagram_mcp_gateway()
    print("已写入:", OUT)


if __name__ == "__main__":
    main()
