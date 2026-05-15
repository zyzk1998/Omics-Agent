#!/usr/bin/env python3
"""
CRISPR-Cas9 Genome Editing Simulator

模拟 CRISPR-Cas9 基因组编辑的完整流程：
  1. 向导RNA (gRNA) 验证
  2. 目标基因组序列中识别 PAM 位点 (NGG) 及 gRNA 匹配
  3. 评估递送效率（基于细胞类型）
  4. 模拟 DNA 双链断裂 (DSB) 和修复（NHEJ / HDR）
  5. 生成原始和编辑后的基因组序列文件

Usage:
    python3 crispr_cas9.py --guides GACGTCAGTCTAGCTAGCTA \
      --target "ATCGGACGTCAGTCTAGCTAGCTAGGCTAGCTAGCTAAGG" \
      --cell HEK293 --output /tmp/crispr_results
    
    python3 crispr_cas9.py \
      --guides "GACGTCAGTCTAGCTAGCTA,ATCGTAGCTAGCTAGCTGCA" \
      --target "ATCGGACGTCAGTCTAGCTAGCTAGGCTAGCTAGCTAAGG" \
      --cell K562 --format json --output /tmp/crispr_results
"""

import argparse
import json
import os
import random
import re
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Tuple

# --- 细胞类型递送效率基线 ---
CELL_DELIVERY_EFFICIENCY = {
    "HEK293": 0.85,
    "HeLa": 0.75,
    "K562": 0.65,
    "HepG2": 0.70,
    "U2OS": 0.80,
    "A549": 0.72,
    "MCF7": 0.68,
    "SH-SY5Y": 0.60,
    "Jurkat": 0.55,
    "mESC": 0.90,
    "hESC": 0.50,
    "iPSC": 0.55,
    "Primary T-cell": 0.40,
    "Primary fibroblast": 0.65,
    "Primary neuron": 0.30,
    "CHO": 0.78,
    "COS7": 0.82,
    "default": 0.70,
}

# PAM 序列正则（NGG 在靶链 3' 端）
PAM_PATTERN = re.compile(r'(?=([ACGT]GG))')

# 碱基互补配对
DNA_COMPLEMENT = {
    'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
    'a': 't', 't': 'a', 'c': 'g', 'g': 'c',
    'N': 'N', 'n': 'n',
}


def reverse_complement(seq: str) -> str:
    """DNA 反向互补"""
    return ''.join(DNA_COMPLEMENT.get(b, b) for b in reversed(seq))


def validate_guide_rna(guide: str) -> Tuple[bool, str]:
    """验证 gRNA 序列（应为 20nt，仅含 A/C/G/T/U）"""
    seq = guide.strip().upper().replace('U', 'T')
    if len(seq) != 20:
        return False, f"Length {len(seq)} != 20 nt"
    if not set(seq).issubset({'A', 'C', 'G', 'T'}):
        return False, f"Invalid characters in: {guide}"
    return True, seq


def find_pam_sites_both_strands(sequence: str) -> List[Dict]:
    """
    在正反向两条链上查找所有潜在 PAM (NGG) 位点。
    返回列表，每项包含：position, pam_seq, strand, target_20mer
    """
    seq = sequence.upper()
    rev = reverse_complement(seq)
    sites = []
    
    # 正向链 (+): 查找 "NGG"，gRNA 靶序列在 PAM 上游（左侧）
    for i in range(len(seq) - 2):
        if seq[i+1:i+3] == 'GG' and seq[i] in 'ACGT':
            # PAM 位置: i, i+1, i+2 = NGG
            # 靶序列在 PAM 上游 20bp
            start = i - 20
            if start >= 0:
                sites.append({
                    "position": i,
                    "pam": seq[i:i+3],
                    "strand": "+",
                    "target_20mer": seq[start:i],
                    "cut_site": i - 3  # PAM 上游 3bp
                })
    
    # 反向链 (-): 在反向互补序列上查找 "NGG"
    # 对应原始序列上为 "CCN"（因为反向互补后 NGG 变成 CCN）
    for i in range(len(seq) - 2):
        if seq[i:i+2] == 'CC' and seq[i+2] in 'ACGT':
            # 原始序列上: i, i+1, i+2 = CCN
            # 这是反向链上的 PAM 互补
            # gRNA 靶序列在 PAM 下游（右侧，即原始序列上 CCN 的右侧 20bp）
            end = i + 3 + 20
            if end <= len(seq):
                target_20 = seq[i+3:end]
                sites.append({
                    "position": i,
                    "pam": "CC" + seq[i+2] + " (rev)",
                    "strand": "-",
                    "target_20mer": target_20,
                    "cut_site": i + 3 + 3  # 对应反向链上 PAM 上游 3bp
                })
    
    return sites


def check_guide_match(guide_seq: str, site: Dict, max_mismatch: int = 3) -> Tuple[bool, int, str]:
    """
    检查 gRNA 是否与靶位点互补配对。
    
    正向链 (+): gRNA 与 target_20mer 互补
    反向链 (-): gRNA 与 reverse_complement(target_20mer) 互补
    """
    guide = guide_seq.upper()
    target_20 = site["target_20mer"].upper()
    
    if site["strand"] == "+":
        # gRNA 互补序列应与 target_20mer 相同（即 gRNA 的互补 = target_20mer）
        comp_guide = reverse_complement(guide)
        seq_to_match = target_20
    else:
        # 反向链：gRNA 应与 target_20mer 的反向互补相同
        # 即 gRNA = reverse_complement(target_20mer)
        # 等价于：reverse_complement(guide) = target_20mer
        comp_guide = reverse_complement(guide)
        seq_to_match = target_20
    
    mismatch = sum(1 for a, b in zip(comp_guide, seq_to_match) if a != b)
    
    return mismatch <= max_mismatch, mismatch, seq_to_match



def simulate_dsb_repair(sequence: str, cut_pos: int, guide_seq: str,
                          efficiency: float, seed: int = None) -> Tuple[str, str, int, int]:
    """
    模拟 Cas9 切割 + NHEJ 修复。
    
    Cas9 切割位置：PAM 上游 3bp 处（双链断裂产生平末端或 1-2nt overhang）
    NHEJ 修复：随机小片段插入/缺失（indel）。
    
    返回: (编辑后序列, 修复类型, 缺失长度, 插入长度)
    """
    if seed is not None:
        random.seed(seed + hash(guide_seq) % 10000)
    
    # 切割位点（PAM 上游 3bp 处）
    cut_point = cut_pos - 3
    if cut_point < 0 or cut_point >= len(sequence):
        return sequence, "no_cut", 0, 0
    
    # 根据递送效率和随机因素决定是否发生编辑
    if random.random() > efficiency:
        return sequence, "no_edit", 0, 0
    
    # NHEJ 修复模拟
    # 大多数 indel 是小的缺失（1-10bp），少数是插入（1-3bp）
    r = random.random()
    
    if r < 0.70:  # 70% 概率：缺失
        del_len = random.choices([1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                                  weights=[15, 18, 15, 12, 10, 8, 6, 5, 4, 7])[0]
        left = sequence[:cut_point]
        right = sequence[cut_point + del_len:]
        new_seq = left + right
        return new_seq, "deletion", del_len, 0
    
    elif r < 0.90:  # 20% 概率：插入
        ins_len = random.choices([1, 2, 3], weights=[50, 35, 15])[0]
        bases = ['A', 'C', 'G', 'T']
        insertion = ''.join(random.choices(bases, k=ins_len))
        new_seq = sequence[:cut_point] + insertion + sequence[cut_point:]
        return new_seq, "insertion", 0, ins_len
    
    else:  # 10% 概率：替换（substitution + small indel）
        sub_len = random.randint(1, 3)
        bases = ['A', 'C', 'G', 'T']
        replacement = ''.join(random.choices(bases, k=sub_len))
        left = sequence[:cut_point]
        right = sequence[cut_point + sub_len:]
        new_seq = left + replacement + right
        return new_seq, "substitution", sub_len, sub_len


def run_crispr_simulation(guides: List[str], target: str, cell_type: str,
                          output_dir: str = None, seed: int = None) -> dict:
    """
    运行完整的 CRISPR-Cas9 模拟。
    """
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    target = target.strip().upper()
    
    result = {
        "timestamp": timestamp,
        "cell_type": cell_type,
        "target_length": len(target),
        "guides": [],
        "total_editing_events": 0,
        "edited_sequence": target,
        "Status": "success",
        "Message": "CRISPR-Cas9 simulation completed"
    }
    
    # 递送效率
    delivery_eff = CELL_DELIVERY_EFFICIENCY.get(cell_type, CELL_DELIVERY_EFFICIENCY["default"])
    result["delivery_efficiency"] = delivery_eff
    
    log_lines = [
        "=" * 60,
        "CRISPR-Cas9 Genome Editing Simulation Log",
        f"Timestamp: {timestamp}",
        f"Cell Type: {cell_type}",
        f"Target Sequence Length: {len(target)} bp",
        f"Delivery Efficiency (baseline): {delivery_eff * 100:.1f}%",
        "=" * 60,
        "",
    ]
    
    # 查找所有 PAM 位点（双向）
    pam_sites = find_pam_sites_both_strands(target)
    log_lines.append(f"[Step 1] Target sequence scanned. {len(pam_sites)} potential PAM (NGG) sites found on both strands.")
    log_lines.append("")
    
    current_sequence = target
    total_edits = 0
    
    for guide_raw in guides:
        # 验证 gRNA
        valid, msg = validate_guide_rna(guide_raw)
        guide_entry = {
            "input": guide_raw,
            "validated": valid,
            "canonical": msg if valid else None,
            "error": None if valid else msg,
            "targets": []
        }
        
        if not valid:
            log_lines.append(f"[gRNA Validation] ❌ '{guide_raw}' — {msg}")
            guide_entry["error"] = msg
            result["guides"].append(guide_entry)
            continue
        
        guide_seq = msg  # 规范化后的序列
        log_lines.append(f"[gRNA Validation] ✅ '{guide_raw}' → canonical: {guide_seq}")
        
        # 在目标序列中搜索匹配位点
        matches_found = 0
        for site in pam_sites:
            matched, mismatch, matched_seq = check_guide_match(guide_seq, site, max_mismatch=3)
            if matched:
                matches_found += 1
                pam_pos = site["position"]
                pam = site["pam"]
                strand = site["strand"]
                cut_site = site["cut_site"]
                
                # 模拟编辑
                new_seq, repair_type, del_len, ins_len = simulate_dsb_repair(
                    current_sequence, cut_site, guide_seq,
                    efficiency=delivery_eff, seed=seed
                )
                
                # 记录靶点信息
                target_info = {
                    "pam_position": pam_pos,
                    "pam_sequence": pam,
                    "strand": strand,
                    "target_20mer": site["target_20mer"],
                    "matched_sequence": matched_seq,
                    "mismatches": mismatch,
                    "cut_position": cut_site,
                    "repair_type": repair_type,
                    "deletion_length": del_len,
                    "insertion_length": ins_len
                }
                
                if repair_type != "no_edit":
                    target_info["edit_successful"] = True
                    total_edits += 1
                    current_sequence = new_seq
                    log_lines.append(
                        f"  → Target #{matches_found}: strand {strand}, PAM@{pam_pos} ({pam}), "
                        f"cut@{cut_site}, {mismatch} mismatch(es), repair={repair_type} "
                        f"(-{del_len}bp / +{ins_len}bp)"
                    )
                else:
                    target_info["edit_successful"] = False
                    log_lines.append(
                        f"  → Target #{matches_found}: strand {strand}, PAM@{pam_pos} ({pam}), "
                        f"no edit (efficiency filtered)"
                    )
                
                guide_entry["targets"].append(target_info)
        
        if matches_found == 0:
            log_lines.append(f"  ⚠ No target sites found for this gRNA.")
        
        guide_entry["match_count"] = matches_found
        result["guides"].append(guide_entry)
        log_lines.append("")
    
    result["total_editing_events"] = total_edits
    result["edited_sequence"] = current_sequence
    
    # 最终总结
    log_lines.extend([
        "=" * 60,
        f"Summary: {total_edits} editing event(s) occurred.",
        f"Original length: {len(target)} bp",
        f"Edited length:   {len(current_sequence)} bp",
        f"Net change:      {len(current_sequence) - len(target):+d} bp",
        "=" * 60,
    ])
    
    result["research_log"] = "\n".join(log_lines)
    
    # 输出文件
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
        # 原始序列
        orig_path = os.path.join(output_dir, "original.txt")
        with open(orig_path, 'w') as fh:
            fh.write(f">original_target_{cell_type}\n")
            fh.write(target + "\n")
        result["original_file"] = orig_path
        
        # 编辑后序列
        mod_path = os.path.join(output_dir, "modified.txt")
        with open(mod_path, 'w') as fh:
            fh.write(f">edited_target_{cell_type}_events{total_edits}\n")
            fh.write(current_sequence + "\n")
        result["modified_file"] = mod_path
        
        log_lines.append(f"\n[Output] Original sequence written to: {orig_path}")
        log_lines.append(f"[Output] Modified sequence written to: {mod_path}")
    
    return result


def format_markdown(result: dict) -> str:
    """格式化为 Markdown"""
    if result["Status"] == "error":
        return f"❌ **Error:** {result['Message']}"
    
    lines = [
        f"**CRISPR-Cas9 Simulation Result**",
        "",
        f"* **Timestamp:** {result['timestamp']}",
        f"* **Cell Type:** {result['cell_type']}",
        f"* **Target Length:** {result['target_length']} bp",
        f"* **Delivery Efficiency:** {result['delivery_efficiency'] * 100:.1f}%",
        "",
        "**Guide RNA(s):**",
        "",
    ]
    
    for g in result["guides"]:
        status = "✅" if g["validated"] else "❌"
        lines.append(f"- {status} `{g['input']}`")
        if g["validated"]:
            lines.append(f"  - Canonical: `{g['canonical']}`")
            lines.append(f"  - Targets found: {g.get('match_count', 0)}")
            if g.get("targets"):
                for t in g["targets"]:
                    edit_mark = "✅" if t.get("edit_successful") else "⬜"
                    lines.append(
                        f"    - {edit_mark} PAM@{t['pam_position']} `{t['pam_sequence']}`, "
                        f"repair={t['repair_type']}"
                    )
        else:
            lines.append(f"  - Error: {g.get('error')}")
    
    lines.extend([
        "",
        "**Editing Summary:**",
        "",
        f"- Total editing events: **{result['total_editing_events']}**",
        f"- Original length: {result['target_length']} bp",
        f"- Edited length: {len(result['edited_sequence'])} bp",
        f"- Net change: {len(result['edited_sequence']) - result['target_length']:+d} bp",
        "",
    ])
    
    if result.get("original_file"):
        lines.append(f"📁 Original: `{result['original_file']}`")
        lines.append(f"📁 Modified: `{result['modified_file']}`")
    
    lines.append("")
    lines.append("**Research Log:**")
    lines.append("```")
    lines.append(result["research_log"])
    lines.append("```")
    
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(
        description='CRISPR-Cas9 Genome Editing Simulator'
    )
    parser.add_argument('--guides', '-g', required=True,
                        help='Comma-separated guide RNA sequences (20nt each)')
    parser.add_argument('--target', '-t', required=True,
                        help='Target genomic DNA sequence')
    parser.add_argument('--cell', '-c', default='HEK293',
                        help='Cell/tissue type (default: HEK293)')
    parser.add_argument('--output', '-o',
                        help='Output directory for result files')
    parser.add_argument('--format', '-f', choices=['json', 'markdown'],
                        default='json', help='Output format')
    parser.add_argument('--seed', type=int, default=None,
                        help='Random seed for reproducible simulation')
    
    args = parser.parse_args()
    
    guides = [g.strip() for g in args.guides.split(',')]
    
    result = run_crispr_simulation(
        guides=guides,
        target=args.target,
        cell_type=args.cell,
        output_dir=args.output,
        seed=args.seed
    )
    
    if args.format == 'json':
        output = json.dumps(result, indent=2, ensure_ascii=False)
    else:
        output = format_markdown(result)
    
    print(output)
    return 0 if result["Status"] == "success" else 1


if __name__ == '__main__':
    sys.exit(main())
