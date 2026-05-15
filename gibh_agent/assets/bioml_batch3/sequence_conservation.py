#!/usr/bin/env python3
"""
sequence_conservation.py

Protein Multiple Sequence Alignment (MSA), Phylogenetic Tree Construction,
and Conservation Analysis Tool.

Dependencies: biopython, scipy, numpy

Usage:
    python sequence_conservation.py --sequences "SEQ1" "SEQ2" --output-dir /tmp/results
    python sequence_conservation.py --json-input '{"protein_sequences":["SEQ1","SEQ2"]}' --output-dir /tmp/results
"""

import argparse
import json
import os
import sys
import tempfile
from typing import Dict, List, Optional, Tuple, Any

import numpy as np
from Bio import AlignIO, Phylo, SeqIO
from Bio.Align import MultipleSeqAlignment, PairwiseAligner
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from scipy.spatial.distance import squareform


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
AA_ALPHABET = "ACDEFGHIKLMNPQRSTVWY-"
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"  # without gap


# ---------------------------------------------------------------------------
# Input Handling
# ---------------------------------------------------------------------------
def parse_input(data: Dict[str, Any]) -> List[Tuple[str, str]]:
    """
    Parse protein_sequences from JSON input.
    Returns list of (seq_id, sequence) tuples.
    """
    sequences = data.get("protein_sequences", [])
    if not sequences:
        raise ValueError("No protein_sequences provided")
    
    result = []
    for i, seq in enumerate(sequences, 1):
        seq_str = str(seq).strip().upper()
        # Remove FASTA header if present
        if seq_str.startswith(">"):
            lines = seq_str.splitlines()
            header = lines[0][1:].strip().split()[0] or f"seq_{i}"
            seq_str = "".join(lines[1:]).replace(" ", "").upper()
            result.append((header, seq_str))
        else:
            result.append((f"seq_{i}", seq_str))
    
    # Validate: only standard amino acids + gaps
    for seq_id, seq in result:
        invalid = set(seq) - set(AA_ALPHABET)
        if invalid:
            raise ValueError(
                f"Sequence {seq_id} contains invalid characters: {invalid}"
            )
    
    return result


# ---------------------------------------------------------------------------
# Multiple Sequence Alignment (Progressive Pairwise)
# ---------------------------------------------------------------------------
def build_distance_matrix(sequences: List[Tuple[str, str]]) -> np.ndarray:
    """Build pairwise distance matrix using global alignment scores."""
    n = len(sequences)
    dist_matrix = np.zeros((n, n))
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 1
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -0.5
    
    max_len = max(len(seq) for _, seq in sequences)
    
    for i in range(n):
        for j in range(i + 1, n):
            _, seq1 = sequences[i]
            _, seq2 = sequences[j]
            # Normalize by max possible score (min length)
            alignment = aligner.align(seq1, seq2)[0]
            score = alignment.score
            norm = max(len(seq1), len(seq2))
            # Convert similarity to distance: d = 1 - score/max
            distance = max(0.0, 1.0 - score / norm)
            dist_matrix[i, j] = distance
            dist_matrix[j, i] = distance
    
    return dist_matrix


def progressive_alignment(sequences: List[Tuple[str, str]]) -> MultipleSeqAlignment:
    """
    Build MSA via progressive pairwise alignment.
    Uses guide tree from distance matrix (neighbor-joining order).
    """
    n = len(sequences)
    if n == 0:
        raise ValueError("No sequences")
    if n == 1:
        seq_id, seq = sequences[0]
        return MultipleSeqAlignment([SeqRecord(Seq(seq), id=seq_id)])
    
    # Distance matrix
    dist_matrix = build_distance_matrix(sequences)
    
    # Pairwise alignments
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 1
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -0.5
    
    # Progressive order: start with closest pair
    # Greedy: merge pair with smallest distance first
    profiles = []
    for seq_id, seq in sequences:
        profiles.append({
            "ids": [seq_id],
            "records": [SeqRecord(Seq(seq), id=seq_id)],
            "consensus": seq,
        })
    
    active = set(range(n))
    
    while len(active) > 1:
        # Find closest pair
        min_dist = float("inf")
        best_pair = (None, None)
        for i in active:
            for j in active:
                if i >= j:
                    continue
                d = dist_matrix[i, j]
                if d < min_dist:
                    min_dist = d
                    best_pair = (i, j)
        
        if best_pair[0] is None:
            break
        
        i, j = best_pair
        # Align the two consensus sequences
        ali = aligner.align(profiles[i]["consensus"], profiles[j]["consensus"])[0]
        
        # Extract aligned sequences with gaps
        aligned1 = str(ali[0])
        aligned2 = str(ali[1])
        
        # Map original records to new alignment
        new_records = []
        new_ids = []
        
        # Profile i records
        for rec in profiles[i]["records"]:
            new_seq = _map_sequence_to_alignment(rec.seq, profiles[i]["consensus"], aligned1)
            new_records.append(SeqRecord(Seq(new_seq), id=rec.id, description=rec.description))
            new_ids.append(rec.id)
        
        # Profile j records
        for rec in profiles[j]["records"]:
            new_seq = _map_sequence_to_alignment(rec.seq, profiles[j]["consensus"], aligned2)
            new_records.append(SeqRecord(Seq(new_seq), id=rec.id, description=rec.description))
            new_ids.append(rec.id)
        
        # Create merged profile
        merged = {
            "ids": new_ids,
            "records": new_records,
            "consensus": aligned1,  # Use aligned1 as consensus
        }
        
        # Replace i with merged, remove j
        profiles[i] = merged
        active.discard(j)
    
    # Build final MSA
    final_idx = list(active)[0]
    return MultipleSeqAlignment(profiles[final_idx]["records"])


def _map_sequence_to_alignment(original_seq: str, consensus_seq: str, aligned_consensus: str) -> str:
    """Map original sequence positions to aligned consensus positions."""
    # Build position map: consensus pos -> aligned pos
    pos_map = {}
    aligned_pos = 0
    consensus_pos = 0
    for char in aligned_consensus:
        if char != "-":
            pos_map[consensus_pos] = aligned_pos
            consensus_pos += 1
        aligned_pos += 1
    
    # Build result
    result = []
    orig_pos = 0
    for aligned_pos, char in enumerate(aligned_consensus):
        if char == "-":
            # Gap in consensus - check if original had gap here
            # For now, use gap
            result.append("-")
        else:
            if orig_pos < len(original_seq):
                result.append(original_seq[orig_pos])
                orig_pos += 1
            else:
                result.append("-")
    
    # Handle remaining original sequence chars
    while orig_pos < len(original_seq):
        result.append(original_seq[orig_pos])
        orig_pos += 1
    
    return "".join(result)


# ---------------------------------------------------------------------------
# Phylogenetic Tree Construction
# ---------------------------------------------------------------------------
def build_phylogenetic_tree(alignment: MultipleSeqAlignment, output_path: str) -> str:
    """Build NJ tree from MSA and save as Newick."""
    calculator = DistanceCalculator("blosum62")
    constructor = DistanceTreeConstructor(calculator, "nj")
    tree = constructor.build_tree(alignment)
    
    # Save as Newick
    Phylo.write(tree, output_path, "newick")
    return output_path


# ---------------------------------------------------------------------------
# Conservation Analysis
# ---------------------------------------------------------------------------
def calculate_shannon_entropy(column: List[str]) -> float:
    """Calculate Shannon entropy for an alignment column."""
    # Count amino acids (exclude gaps)
    counts = {}
    total = 0
    for aa in column:
        if aa != "-" and aa in AMINO_ACIDS:
            counts[aa] = counts.get(aa, 0) + 1
            total += 1
    
    if total == 0:
        return 0.0
    
    entropy = 0.0
    for count in counts.values():
        p = count / total
        if p > 0:
            entropy -= p * np.log2(p)
    
    return entropy


def calculate_conservation_score(column: List[str]) -> float:
    """
    Conservation score: 1 - (entropy / max_entropy).
    Higher = more conserved.
    Max entropy for 20 amino acids = log2(20) ≈ 4.32
    """
    entropy = calculate_shannon_entropy(column)
    max_entropy = np.log2(20)
    return max(0.0, 1.0 - entropy / max_entropy)


def get_consensus_residue(column: List[str]) -> str:
    """Get most frequent amino acid in column (ignoring gaps)."""
    counts = {}
    for aa in column:
        if aa != "-" and aa in AMINO_ACIDS:
            counts[aa] = counts.get(aa, 0) + 1
    
    if not counts:
        return "-"
    return max(counts, key=counts.get)


def analyze_conservation(alignment: MultipleSeqAlignment) -> Dict[str, Any]:
    """Perform conservation analysis on MSA."""
    n_seqs = len(alignment)
    n_cols = alignment.get_alignment_length()
    
    results = {
        "num_sequences": n_seqs,
        "alignment_length": n_cols,
        "columns": [],
        "summary": {
            "mean_entropy": 0.0,
            "mean_conservation": 0.0,
            "highly_conserved_sites": 0,  # conservation > 0.8
            "moderately_conserved_sites": 0,  # 0.5-0.8
            "variable_sites": 0,  # < 0.5
        }
    }
    
    total_entropy = 0.0
    total_conservation = 0.0
    
    for col_idx in range(n_cols):
        column = [str(alignment[i, col_idx]) for i in range(n_seqs)]
        entropy = calculate_shannon_entropy(column)
        conservation = calculate_conservation_score(column)
        consensus = get_consensus_residue(column)
        
        # Count amino acid frequencies
        freq = {}
        for aa in column:
            if aa != "-":
                freq[aa] = freq.get(aa, 0) + 1
        
        col_result = {
            "position": col_idx + 1,
            "entropy": round(entropy, 4),
            "conservation_score": round(conservation, 4),
            "consensus_residue": consensus,
            "amino_acid_frequencies": {k: v for k, v in freq.items()},
        }
        results["columns"].append(col_result)
        
        total_entropy += entropy
        total_conservation += conservation
        
        if conservation > 0.8:
            results["summary"]["highly_conserved_sites"] += 1
        elif conservation > 0.5:
            results["summary"]["moderately_conserved_sites"] += 1
        else:
            results["summary"]["variable_sites"] += 1
    
    results["summary"]["mean_entropy"] = round(total_entropy / n_cols, 4)
    results["summary"]["mean_conservation"] = round(total_conservation / n_cols, 4)
    
    return results


def write_conservation_report(conservation_data: Dict, output_path: str) -> str:
    """Write conservation analysis as a formatted text report."""
    lines = []
    lines.append("=" * 70)
    lines.append("PROTEIN SEQUENCE CONSERVATION ANALYSIS REPORT")
    lines.append("=" * 70)
    lines.append("")
    
    summary = conservation_data["summary"]
    lines.append(f"Number of sequences: {conservation_data['num_sequences']}")
    lines.append(f"Alignment length: {conservation_data['alignment_length']}")
    lines.append("")
    lines.append("-" * 70)
    lines.append("SUMMARY STATISTICS")
    lines.append("-" * 70)
    lines.append(f"Mean Shannon Entropy:      {summary['mean_entropy']:.4f}")
    lines.append(f"Mean Conservation Score: {summary['mean_conservation']:.4f}")
    lines.append("")
    lines.append(f"Highly conserved sites (score > 0.8):     {summary['highly_conserved_sites']}")
    lines.append(f"Moderately conserved sites (0.5-0.8):     {summary['moderately_conserved_sites']}")
    lines.append(f"Variable sites (score < 0.5):             {summary['variable_sites']}")
    lines.append("")
    
    lines.append("-" * 70)
    lines.append("PER-POSITION CONSERVATION DETAILS")
    lines.append("-" * 70)
    lines.append(f"{'Pos':>5} {'Cons':>8} {'Entropy':>10} {'Consensus':>12} {'Frequencies'}")
    lines.append("-" * 70)
    
    for col in conservation_data["columns"]:
        pos = col["position"]
        cons = col["conservation_score"]
        entropy = col["entropy"]
        consensus = col["consensus_residue"]
        freq_str = ", ".join([f"{k}:{v}" for k, v in col["amino_acid_frequencies"].items()])
        lines.append(f"{pos:>5} {cons:>8.4f} {entropy:>10.4f} {consensus:>12} {freq_str}")
    
    lines.append("")
    lines.append("=" * 70)
    lines.append("END OF REPORT")
    lines.append("=" * 70)
    
    with open(output_path, "w") as f:
        f.write("\n".join(lines))
    
    return output_path


def write_conservation_json(conservation_data: Dict, output_path: str) -> str:
    """Write conservation analysis as JSON."""
    with open(output_path, "w") as f:
        json.dump(conservation_data, f, indent=2)
    return output_path


# ---------------------------------------------------------------------------
# Main Pipeline
# ---------------------------------------------------------------------------
def run_analysis(input_data: Dict[str, Any], output_dir: str) -> Dict[str, Any]:
    """
    Main analysis pipeline.
    
    Returns dict with:
        - research_log: analysis log string
        - alignment_file_path: path to MSA FASTA
        - tree_file_path: path to Newick tree
        - conservation_file_path: path to conservation report
        - error: error message if any
    """
    log_lines = []
    
    try:
        # 1. Parse input
        log_lines.append("Step 1: Parsing input sequences...")
        sequences = parse_input(input_data)
        log_lines.append(f"  - Parsed {len(sequences)} sequences")
        for seq_id, seq in sequences:
            log_lines.append(f"  - {seq_id}: {len(seq)} residues")
        
        os.makedirs(output_dir, exist_ok=True)
        
        # 2. Multiple Sequence Alignment
        log_lines.append("")
        log_lines.append("Step 2: Performing multiple sequence alignment (progressive pairwise)...")
        alignment = progressive_alignment(sequences)
        align_length = alignment.get_alignment_length()
        log_lines.append(f"  - Alignment complete: {align_length} columns")
        
        # Save alignment
        alignment_path = os.path.join(output_dir, "alignment.fasta")
        AlignIO.write(alignment, alignment_path, "fasta")
        log_lines.append(f"  - Saved alignment to: {alignment_path}")
        
        # 3. Phylogenetic Tree
        log_lines.append("")
        log_lines.append("Step 3: Building phylogenetic tree (Neighbor-Joining)...")
        tree_path = os.path.join(output_dir, "phylogenetic_tree.nwk")
        build_phylogenetic_tree(alignment, tree_path)
        log_lines.append(f"  - Tree saved to: {tree_path}")
        
        # 4. Conservation Analysis
        log_lines.append("")
        log_lines.append("Step 4: Running conservation analysis...")
        conservation = analyze_conservation(alignment)
        
        summary = conservation["summary"]
        log_lines.append(f"  - Mean entropy: {summary['mean_entropy']:.4f}")
        log_lines.append(f"  - Mean conservation: {summary['mean_conservation']:.4f}")
        log_lines.append(f"  - Highly conserved sites (>0.8): {summary['highly_conserved_sites']}")
        log_lines.append(f"  - Moderately conserved (0.5-0.8): {summary['moderately_conserved_sites']}")
        log_lines.append(f"  - Variable sites (<0.5): {summary['variable_sites']}")
        
        conservation_path = os.path.join(output_dir, "conservation_report.txt")
        write_conservation_report(conservation, conservation_path)
        log_lines.append(f"  - Report saved to: {conservation_path}")
        
        # Also save JSON version
        conservation_json_path = os.path.join(output_dir, "conservation_data.json")
        write_conservation_json(conservation, conservation_json_path)
        log_lines.append(f"  - JSON data saved to: {conservation_json_path}")
        
        log_lines.append("")
        log_lines.append("Analysis complete!")
        
        return {
            "research_log": "\n".join(log_lines),
            "alignment_file_path": alignment_path,
            "tree_file_path": tree_path,
            "conservation_file_path": conservation_path,
            "conservation_json_path": conservation_json_path,
            "error": None,
        }
    
    except Exception as e:
        log_lines.append(f"ERROR: {str(e)}")
        return {
            "research_log": "\n".join(log_lines),
            "alignment_file_path": None,
            "tree_file_path": None,
            "conservation_file_path": None,
            "error": str(e),
        }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="Protein Sequence Conservation Analysis")
    parser.add_argument("--json-input", type=str, help="JSON string with protein_sequences array")
    parser.add_argument("--sequences", nargs="+", help="Protein sequences as arguments")
    parser.add_argument("--output-dir", type=str, default="/tmp/sequence_conservation", help="Output directory")
    parser.add_argument("--output-json", type=str, help="Write full result as JSON to this path")
    
    args = parser.parse_args()
    
    if args.json_input:
        input_data = json.loads(args.json_input)
    elif args.sequences:
        input_data = {"protein_sequences": args.sequences}
    else:
        print("Error: Provide --json-input or --sequences", file=sys.stderr)
        sys.exit(1)
    
    result = run_analysis(input_data, args.output_dir)
    
    # Print result as JSON
    print(json.dumps(result, indent=2))
    
    if args.output_json:
        with open(args.output_json, "w") as f:
            json.dump(result, f, indent=2)


if __name__ == "__main__":
    main()
