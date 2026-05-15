#!/usr/bin/env python3
"""
RNA Secondary Structure Analysis Tool

Analyzes RNA secondary structures represented in dot-bracket notation.
Computes base pairing, stems, loops, and estimates free energy when sequence is provided.

Usage:
    python rna_structure.py --input '{"structure": "..((..))", "sequence": "AAGGCCUU"}'
    python rna_structure.py --file structure.db
    python rna_structure.py --file structure.fa  (FASTA format)
"""

import argparse
import csv
import io
import json
import re
import sys
from dataclasses import asdict, dataclass, field
from typing import Dict, List, Optional, Tuple

# Simplified free energy parameters (kcal/mol) at 37°C
# Based on Turner nearest-neighbor model approximations
STACKING_ENERGY = {
    ('A', 'U', 'A', 'U'): -0.9,
    ('A', 'U', 'C', 'G'): -1.8,
    ('A', 'U', 'G', 'C'): -2.3,
    ('A', 'U', 'U', 'A'): -1.1,
    ('C', 'G', 'A', 'U'): -1.7,
    ('C', 'G', 'C', 'G'): -2.9,
    ('C', 'G', 'G', 'C'): -2.0,
    ('C', 'G', 'U', 'A'): -1.8,
    ('G', 'C', 'A', 'U'): -2.4,
    ('G', 'C', 'C', 'G'): -3.4,
    ('G', 'C', 'G', 'C'): -2.9,
    ('G', 'C', 'U', 'A'): -2.1,
    ('G', 'U', 'A', 'U'): -0.6,
    ('G', 'U', 'C', 'G'): -1.0,
    ('G', 'U', 'G', 'C'): -1.3,
    ('G', 'U', 'U', 'A'): -0.8,
    ('U', 'A', 'A', 'U'): -1.0,
    ('U', 'A', 'C', 'G'): -1.7,
    ('U', 'A', 'G', 'C'): -2.1,
    ('U', 'A', 'U', 'A'): -0.9,
    ('U', 'G', 'A', 'U'): -0.5,
    ('U', 'G', 'C', 'G'): -1.1,
    ('U', 'G', 'G', 'C'): -1.4,
    ('U', 'G', 'U', 'A'): -0.6,
}

# Terminal AU/GU penalty
TERMINAL_AU_PENALTY = 0.5

# Hairpin loop initiation penalties (simplified)
HAIRPIN_INITIATION = {
    3: 7.0,
    4: 5.5,
    5: 4.5,
    6: 4.2,
    7: 4.1,
    8: 4.1,
    9: 4.0,
    10: 4.0,
}

# Interior/bulge loop penalties (simplified per unpaired nucleotide)
INTERIOR_PER_NT = 0.5
BULGE_PER_NT = 3.0


@dataclass
class BasePair:
    i: int  # 0-based index of 5' base
    j: int  # 0-based index of 3' base
    pair_type: str  # e.g., "AU", "GC", "GU"


@dataclass
class Stem:
    id: int
    start_5p: int
    end_5p: int
    start_3p: int
    end_3p: int
    length: int
    pairs: List[BasePair] = field(default_factory=list)
    estimated_energy: float = 0.0
    stability_score: float = 0.0


@dataclass
class Loop:
    id: int
    loop_type: str  # hairpin, interior, bulge, multi
    start: int
    end: int
    size: int
    closing_stems: List[int] = field(default_factory=list)


@dataclass
class StructureResult:
    valid: bool = False
    error_message: str = ""
    structure: str = ""
    sequence: str = ""
    """与 structure 等长的规范序列（用于配对类型与能量）；可能与用户原始输入长度不同。"""
    sequence_original: str = ""
    """用户传入的原始序列（去空白、T→U），未补齐前；用于报告回显。"""
    length: int = 0
    alignment_note: str = ""
    """序列与点括号长度不一致时的对齐说明（如 3' 端 N 补齐 / 截断）。"""
    base_pairs: List[BasePair] = field(default_factory=list)
    stems: List[Stem] = field(default_factory=list)
    loops: List[Loop] = field(default_factory=list)
    unpaired_count: int = 0
    unpaired_percentage: float = 0.0
    estimated_free_energy: float = 0.0
    gc_content: float = 0.0


def validate_structure(structure: str) -> Tuple[bool, str]:
    """Validate dot-bracket notation."""
    if not structure:
        return False, "Empty structure"
    
    valid_chars = set('.()[]{}<><>')
    if not all(c in valid_chars for c in structure):
        invalid = set(structure) - valid_chars
        return False, f"Invalid characters in structure: {invalid}"
    
    # Use a stack to check pairing balance
    stack = []
    bracket_pairs = {'(': ')', '[': ']', '{': '}', '<': '>'}
    open_brackets = set(bracket_pairs.keys())
    close_brackets = set(bracket_pairs.values())
    
    for i, char in enumerate(structure):
        if char == '.':
            continue
        elif char in open_brackets:
            stack.append((char, i))
        elif char in close_brackets:
            if not stack:
                return False, f"Unmatched closing bracket '{char}' at position {i}"
            last_open, last_pos = stack.pop()
            if bracket_pairs[last_open] != char:
                return False, f"Mismatched brackets: '{last_open}' at {last_pos} and '{char}' at {i}"
    
    if stack:
        unmatched = ', '.join(f"'{b}' at {p}" for b, p in stack)
        return False, f"Unmatched opening brackets: {unmatched}"
    
    return True, ""


def parse_base_pairs(structure: str) -> List[BasePair]:
    """Extract all base pairs from structure."""
    pairs = []
    stack = []
    bracket_pairs = {'(': ')', '[': ']', '{': '}', '<': '>'}
    open_brackets = set(bracket_pairs.keys())
    
    for i, char in enumerate(structure):
        if char in open_brackets:
            stack.append((char, i))
        elif char in bracket_pairs.values():
            if stack:
                last_open, last_pos = stack.pop()
                pairs.append(BasePair(i=last_pos, j=i, pair_type=""))
    
    return sorted(pairs, key=lambda x: x.i)


def assign_pair_types(pairs: List[BasePair], sequence: str) -> List[BasePair]:
    """Assign pair types (AU, GC, GU) based on sequence."""
    result = []
    for bp in pairs:
        if bp.i < 0 or bp.j < 0 or bp.i >= len(sequence) or bp.j >= len(sequence):
            raise IndexError(f"base pair index out of range: i={bp.i}, j={bp.j}, len(seq)={len(sequence)}")
        b1 = sequence[bp.i].upper()
        b2 = sequence[bp.j].upper()
        
        # Canonical pairs: AU, UA, GC, CG, GU, UG
        pair_set = frozenset([b1, b2])
        if pair_set == {'A', 'U'}:
            bp.pair_type = "AU"
        elif pair_set == {'G', 'C'}:
            bp.pair_type = "GC"
        elif pair_set == {'G', 'U'}:
            bp.pair_type = "GU"
        else:
            bp.pair_type = f"{b1}-{b2}"
        result.append(bp)
    return result


def find_stems(pairs: List[BasePair]) -> List[Stem]:
    """Identify stems (continuous paired regions)."""
    if not pairs:
        return []
    
    stems = []
    current_stem_pairs = [pairs[0]]
    stem_id = 1
    
    for i in range(1, len(pairs)):
        prev = current_stem_pairs[-1]
        curr = pairs[i]
        
        # Check if consecutive (i+1, j-1)
        if curr.i == prev.i + 1 and curr.j == prev.j - 1:
            current_stem_pairs.append(curr)
        else:
            # Save current stem
            stem = Stem(
                id=stem_id,
                start_5p=current_stem_pairs[0].i,
                end_5p=current_stem_pairs[-1].i,
                start_3p=current_stem_pairs[-1].j,
                end_3p=current_stem_pairs[0].j,
                length=len(current_stem_pairs),
                pairs=current_stem_pairs
            )
            stems.append(stem)
            stem_id += 1
            current_stem_pairs = [curr]
    
    # Don't forget the last stem
    if current_stem_pairs:
        stem = Stem(
            id=stem_id,
            start_5p=current_stem_pairs[0].i,
            end_5p=current_stem_pairs[-1].i,
            start_3p=current_stem_pairs[-1].j,
            end_3p=current_stem_pairs[0].j,
            length=len(current_stem_pairs),
            pairs=current_stem_pairs
        )
        stems.append(stem)
    
    return stems


def find_loops(structure: str, stems: List[Stem]) -> List[Loop]:
    """Identify loop regions."""
    loops = []
    loop_id = 1
    
    if not stems:
        # Entirely unpaired - one big loop
        unpaired = [i for i, c in enumerate(structure) if c == '.']
        if unpaired:
            loops.append(Loop(
                id=loop_id,
                loop_type="unpaired",
                start=min(unpaired),
                end=max(unpaired),
                size=len(unpaired),
                closing_stems=[]
            ))
        return loops
    
    # Find hairpin loops (between ends of a single stem)
    for stem in stems:
        if stem.end_5p + 1 <= stem.start_3p - 1:
            loop_size = stem.start_3p - stem.end_5p - 1
            if loop_size >= 0:
                loop_type = "hairpin" if loop_size > 0 else "none"
                if loop_size == 0:
                    continue
                loops.append(Loop(
                    id=loop_id,
                    loop_type="hairpin",
                    start=stem.end_5p + 1,
                    end=stem.start_3p - 1,
                    size=loop_size,
                    closing_stems=[stem.id]
                ))
                loop_id += 1
    
    # Find interior/bulge/multi loops between stems
    # Simplified: look for unpaired regions between stems
    for i in range(len(stems) - 1):
        stem1 = stems[i]
        stem2 = stems[i + 1]
        
        # Check if there's a gap between stem1's 5' end and stem2's 5' end
        # or between their 3' ends
        if stem1.end_5p < stem2.start_5p and stem2.start_3p < stem1.start_3p:
            # Potential interior loop
            left_gap = stem2.start_5p - stem1.end_5p - 1
            right_gap = stem1.start_3p - stem2.end_3p - 1
            
            if left_gap > 0 or right_gap > 0:
                loop_size = left_gap + right_gap
                if left_gap > 0 and right_gap > 0:
                    loop_type = "interior"
                elif (left_gap > 0 and right_gap == 0) or (left_gap == 0 and right_gap > 0):
                    loop_type = "bulge"
                else:
                    loop_type = "multi"
                
                loops.append(Loop(
                    id=loop_id,
                    loop_type=loop_type,
                    start=stem1.end_5p + 1,
                    end=stem2.start_3p - 1 if stem2.start_3p > stem1.end_5p else stem1.start_3p - 1,
                    size=loop_size,
                    closing_stems=[stem1.id, stem2.id]
                ))
                loop_id += 1
    
    # Find unpaired regions at 5' and 3' ends
    # 5' end
    first_stem = min(stems, key=lambda s: s.start_5p)
    if first_stem.start_5p > 0:
        unpaired_5p = [c for c in structure[:first_stem.start_5p] if c == '.']
        if unpaired_5p:
            loops.append(Loop(
                id=loop_id,
                loop_type="terminal_5p",
                start=0,
                end=first_stem.start_5p - 1,
                size=len(unpaired_5p),
                closing_stems=[]
            ))
            loop_id += 1
    
    # 3' end
    last_stem = max(stems, key=lambda s: s.end_3p)
    if last_stem.end_3p < len(structure) - 1:
        unpaired_3p = [c for c in structure[last_stem.end_3p + 1:] if c == '.']
        if unpaired_3p:
            loops.append(Loop(
                id=loop_id,
                loop_type="terminal_3p",
                start=last_stem.end_3p + 1,
                end=len(structure) - 1,
                size=len(unpaired_3p),
                closing_stems=[]
            ))
            loop_id += 1
    
    return loops


def estimate_stem_energy(stem: Stem, sequence: str) -> float:
    """Estimate free energy contribution of a stem."""
    energy = 0.0
    
    for i, bp in enumerate(stem.pairs):
        b1 = sequence[bp.i].upper()
        b2 = sequence[bp.j].upper()
        
        # Terminal AU/GU penalty for first and last pair
        if i == 0 or i == len(stem.pairs) - 1:
            if bp.pair_type in ("AU", "GU"):
                energy += TERMINAL_AU_PENALTY
        
        # Stacking energy for consecutive pairs
        if i < len(stem.pairs) - 1:
            next_bp = stem.pairs[i + 1]
            b1_next = sequence[next_bp.i].upper()
            b2_next = sequence[next_bp.j].upper()
            
            key = (b1, b2, b1_next, b2_next)
            energy += STACKING_ENERGY.get(key, -1.5)  # Default stacking
    
    return energy


def estimate_loop_energy(loop: Loop) -> float:
    """Estimate free energy contribution of a loop."""
    if loop.loop_type == "hairpin":
        size = loop.size
        if size in HAIRPIN_INITIATION:
            return HAIRPIN_INITIATION[size]
        elif size < 3:
            return 10.0  # Very unstable
        else:
            return 4.0 + 0.1 * (size - 9)  # Linear extrapolation
    elif loop.loop_type == "interior":
        return loop.size * INTERIOR_PER_NT
    elif loop.loop_type == "bulge":
        return 3.5 + loop.size * BULGE_PER_NT
    elif loop.loop_type in ("terminal_5p", "terminal_3p", "unpaired"):
        return 0.0  # Terminal unpaired regions don't contribute
    else:
        return loop.size * 0.5


def calculate_stability_score(stem: Stem) -> float:
    """Calculate a stability score (0-100) for a stem."""
    score = 0.0
    
    # Length contribution (longer = more stable)
    score += min(stem.length * 5, 40)
    
    # GC content contribution
    gc_count = sum(1 for bp in stem.pairs if bp.pair_type == "GC")
    if stem.length > 0:
        gc_ratio = gc_count / stem.length
        score += gc_ratio * 40
    
    # AU/GU penalty
    au_gu_count = sum(1 for bp in stem.pairs if bp.pair_type in ("AU", "GU"))
    if stem.length > 0:
        au_gu_ratio = au_gu_count / stem.length
        score -= au_gu_ratio * 20
    
    return max(0, min(100, score))


def _normalize_sequence_input(raw: Optional[str]) -> str:
    """去空白，T→U，仅保留 IUPAC RNA 字母。"""
    if not raw:
        return ""
    return "".join(c for c in raw.upper().replace("T", "U").strip() if c in "AUCGN")


def _harmonize_sequence_to_structure(structure: str, seq_norm: str) -> Tuple[str, str]:
    """
    使序列与点括号等长以便索引配对。
    较短时在 3' 端用 N 补齐；较长时截断至结构长度（与结构权威对齐）。
    返回 (对齐后序列, 人类可读说明)。
    """
    n = len(structure)
    if not seq_norm:
        return "", ""
    if len(seq_norm) == n:
        return seq_norm, ""
    if len(seq_norm) < n:
        pad = "N" * (n - len(seq_norm))
        return seq_norm + pad, (
            f"输入序列长度为 {len(seq_norm)} nt，点括号为 {n} nt；已在 3' 端用 **N** 补齐至 {n} nt 以计算配对类型与能量。"
        )
    return seq_norm[:n], (
        f"输入序列长度为 {len(seq_norm)} nt，点括号为 {n} nt；已 **截断至 {n} nt**（与点括号对齐）以计算配对。"
    )


def analyze_structure(structure: str, sequence: Optional[str] = None) -> StructureResult:
    """Main analysis function."""
    structure = (structure or "").strip()
    result = StructureResult(structure=structure)

    # Validate
    valid, error = validate_structure(structure)
    result.valid = valid
    if not valid:
        result.error_message = error
        return result

    n = len(structure)
    result.length = n

    # Count unpaired
    result.unpaired_count = structure.count(".")
    result.unpaired_percentage = (result.unpaired_count / n * 100) if n else 0.0

    # Parse base pairs
    result.base_pairs = parse_base_pairs(structure)

    seq_raw = _normalize_sequence_input(sequence)
    result.sequence_original = seq_raw

    if seq_raw:
        seq_aligned, note = _harmonize_sequence_to_structure(structure, seq_raw)
        result.sequence = seq_aligned
        result.alignment_note = note
        if len(seq_aligned) == n:
            result.base_pairs = assign_pair_types(result.base_pairs, result.sequence)
            gc = sum(1 for c in result.sequence if c in "GC")
            result.gc_content = round(gc / n * 100, 2) if n else 0.0
    else:
        result.sequence = ""
        result.alignment_note = ""

    # Find stems
    result.stems = find_stems(result.base_pairs)
    
    # Find loops
    result.loops = find_loops(structure, result.stems)
    
    # Calculate energies if sequence provided (已与 structure 等长)
    if result.sequence and len(result.sequence) == n:
        total_energy = 0.0

        for stem in result.stems:
            stem.estimated_energy = estimate_stem_energy(stem, result.sequence)
            stem.stability_score = calculate_stability_score(stem)
            total_energy += stem.estimated_energy

        for loop in result.loops:
            total_energy += estimate_loop_energy(loop)

        result.estimated_free_energy = round(total_energy, 2)

    return result


def generate_txt_report(result: StructureResult) -> str:
    """Generate detailed TXT report."""
    lines = []
    lines.append("=" * 60)
    lines.append("RNA SECONDARY STRUCTURE ANALYSIS REPORT")
    lines.append("=" * 60)
    lines.append("")
    
    # Input Information
    lines.append("INPUT INFORMATION")
    lines.append("-" * 40)
    lines.append(f"Structure length: {result.length} nt")
    lines.append(f"Structure: {result.structure}")
    if result.sequence_original:
        lines.append(f"Sequence (input): {result.sequence_original}")
    if result.alignment_note:
        lines.append(f"Alignment note: {result.alignment_note}")
    if result.sequence:
        lines.append(f"Sequence (used): {result.sequence}")
        lines.append(f"GC content: {result.gc_content:.1f}%")
    lines.append("")
    
    # Validation
    lines.append("VALIDATION")
    lines.append("-" * 40)
    if result.valid:
        lines.append("Status: VALID")
        lines.append("The dot-bracket notation is well-formed and balanced.")
    else:
        lines.append("Status: INVALID")
        lines.append(f"Error: {result.error_message}")
    lines.append("")
    
    if not result.valid:
        return '\n'.join(lines)
    
    # Structure Features
    lines.append("STRUCTURE FEATURES")
    lines.append("-" * 40)
    lines.append(f"Total base pairs: {len(result.base_pairs)}")
    lines.append(f"Paired bases: {len(result.base_pairs) * 2} ({100 - result.unpaired_percentage:.1f}%)")
    lines.append(f"Unpaired bases: {result.unpaired_count} ({result.unpaired_percentage:.1f}%)")
    lines.append(f"Number of stems: {len(result.stems)}")
    lines.append(f"Number of loops: {len(result.loops)}")
    lines.append("")
    
    # Pair type distribution
    if result.sequence:
        lines.append("PAIR TYPE DISTRIBUTION")
        lines.append("-" * 40)
        pair_types = {}
        for bp in result.base_pairs:
            pair_types[bp.pair_type] = pair_types.get(bp.pair_type, 0) + 1
        for pt, count in sorted(pair_types.items()):
            lines.append(f"  {pt}: {count}")
        lines.append("")
    
    # Stem details
    lines.append("STEM DETAILS")
    lines.append("-" * 40)
    if result.stems:
        for stem in result.stems:
            lines.append(f"Stem #{stem.id}:")
            lines.append(f"  5' arm: [{stem.start_5p}-{stem.end_5p}] (length: {stem.length} bp)")
            lines.append(f"  3' arm: [{stem.start_3p}-{stem.end_3p}]")
            if result.sequence:
                lines.append(f"  Estimated energy: {stem.estimated_energy:.2f} kcal/mol")
                lines.append(f"  Stability score: {stem.stability_score:.1f}/100")
            lines.append("")
    else:
        lines.append("No stems detected.")
        lines.append("")
    
    # Loop details
    lines.append("LOOP DETAILS")
    lines.append("-" * 40)
    if result.loops:
        for loop in result.loops:
            lines.append(f"Loop #{loop.id}:")
            lines.append(f"  Type: {loop.loop_type}")
            lines.append(f"  Position: [{loop.start}-{loop.end}]")
            lines.append(f"  Size: {loop.size} nt")
            if loop.closing_stems:
                lines.append(f"  Closing stems: {', '.join(map(str, loop.closing_stems))}")
            lines.append("")
    else:
        lines.append("No loops detected.")
        lines.append("")
    
    # Energy
    if result.sequence:
        lines.append("FREE ENERGY ESTIMATION")
        lines.append("-" * 40)
        lines.append(f"Total estimated free energy: {result.estimated_free_energy:.2f} kcal/mol")
        lines.append("")
        lines.append("Note: Energy values are rough estimates based on simplified")
        lines.append("nearest-neighbor parameters. For accurate predictions, use")
        lines.append("specialized tools like ViennaRNA RNAfold or UNAFold.")
        lines.append("")
    
    lines.append("=" * 60)
    lines.append("END OF REPORT")
    lines.append("=" * 60)
    
    return '\n'.join(lines)


def generate_csv_data(result: StructureResult) -> str:
    """Generate CSV data for stems."""
    output = io.StringIO()
    writer = csv.writer(output)
    
    # Header
    headers = ["stem_id", "start_5p", "end_5p", "start_3p", "end_3p", 
               "length_bp", "pair_types", "estimated_energy_kcal_mol", "stability_score"]
    writer.writerow(headers)
    
    # Data
    for stem in result.stems:
        pair_types = ', '.join(bp.pair_type for bp in stem.pairs)
        writer.writerow([
            stem.id,
            stem.start_5p,
            stem.end_5p,
            stem.start_3p,
            stem.end_3p,
            stem.length,
            pair_types,
            f"{stem.estimated_energy:.2f}" if result.sequence else "N/A",
            f"{stem.stability_score:.1f}" if result.sequence else "N/A"
        ])
    
    return output.getvalue()


def parse_fasta(content: str) -> Tuple[Optional[str], Optional[str]]:
    """Parse FASTA format."""
    lines = content.strip().split('\n')
    sequence = []
    structure = None
    
    for line in lines:
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            continue
        elif line.startswith('#'):
            # Some FASTA variants include structure in comments
            continue
        else:
            # Check if it's a structure line (contains only .()[]{})
            if set(line).issubset({'.', '(', ')', '[', ']', '{', '}', '<', '>'}):
                structure = line
            else:
                sequence.append(line)
    
    seq_str = ''.join(sequence).upper() if sequence else None
    return seq_str, structure


def parse_dot_bracket_file(content: str) -> Tuple[Optional[str], Optional[str]]:
    """Parse dot-bracket file (sequence and structure on separate lines)."""
    lines = [l.strip() for l in content.strip().split('\n') if l.strip()]
    
    sequence = None
    structure = None
    
    for line in lines:
        if line.startswith('>') or line.startswith('#'):
            continue
        elif set(line).issubset({'.', '(', ')', '[', ']', '{', '}', '<', '>'}):
            structure = line
        elif set(line).issubset({'A', 'C', 'G', 'U', 'T', 'N', 'a', 'c', 'g', 'u', 't', 'n'}):
            sequence = line.upper().replace('T', 'U')
    
    return sequence, structure


def main():
    parser = argparse.ArgumentParser(description='RNA Secondary Structure Analysis')
    parser.add_argument('--input', type=str, help='JSON string with structure and optional sequence')
    parser.add_argument('--file', type=str, help='Input file (FASTA or dot-bracket format)')
    parser.add_argument('--output-txt', type=str, help='Output TXT report file path')
    parser.add_argument('--output-csv', type=str, help='Output CSV data file path')
    parser.add_argument('--format', type=str, choices=['json', 'txt', 'csv', 'all'], default='json',
                        help='Output format')
    
    args = parser.parse_args()
    
    structure = None
    sequence = None
    
    if args.input:
        try:
            data = json.loads(args.input)
            structure = (data.get("structure") or "").strip()
            sequence = data.get("sequence")
            if sequence is not None and not isinstance(sequence, str):
                sequence = str(sequence)
        except json.JSONDecodeError as e:
            print(json.dumps({"error": f"Invalid JSON input: {e}"}, indent=2))
            sys.exit(1)
    elif args.file:
        try:
            with open(args.file, 'r') as f:
                content = f.read()
            
            # Try to determine format
            if content.startswith('>'):
                sequence, structure = parse_fasta(content)
            else:
                sequence, structure = parse_dot_bracket_file(content)
            
            if not structure:
                print(json.dumps({"error": "Could not find structure in file"}, indent=2))
                sys.exit(1)
        except FileNotFoundError:
            print(json.dumps({"error": f"File not found: {args.file}"}, indent=2))
            sys.exit(1)
    else:
        # Read from stdin
        try:
            data = json.loads(sys.stdin.read())
            structure = (data.get("structure") or "").strip()
            sequence = data.get("sequence")
            if sequence is not None and not isinstance(sequence, str):
                sequence = str(sequence)
        except json.JSONDecodeError:
            print(json.dumps({"error": "No valid input provided"}, indent=2))
            sys.exit(1)
    
    # Run analysis
    result = analyze_structure(structure, sequence)
    
    # Output
    if args.format == 'json' or args.format == 'all':
        # Convert to dict for JSON serialization
        result_dict = {
            "valid": result.valid,
            "error_message": result.error_message,
            "structure": result.structure,
            "sequence": result.sequence,
            "sequence_original": result.sequence_original,
            "length": result.length,
            "alignment_note": result.alignment_note,
            "base_pairs": [{"i": bp.i, "j": bp.j, "pair_type": bp.pair_type} for bp in result.base_pairs],
            "stems": [{"id": s.id, "start_5p": s.start_5p, "end_5p": s.end_5p,
                       "start_3p": s.start_3p, "end_3p": s.end_3p, "length": s.length,
                       "estimated_energy": s.estimated_energy, "stability_score": s.stability_score}
                      for s in result.stems],
            "loops": [{"id": l.id, "loop_type": l.loop_type, "start": l.start,
                       "end": l.end, "size": l.size, "closing_stems": l.closing_stems}
                      for l in result.loops],
            "unpaired_count": result.unpaired_count,
            "unpaired_percentage": result.unpaired_percentage,
            "estimated_free_energy": result.estimated_free_energy,
            "gc_content": result.gc_content,
        }
        
        json_output = json.dumps(result_dict, indent=2)
        
        if args.output_txt or args.output_csv:
            print(json_output)
        else:
            print(json_output)
    
    if args.format in ('txt', 'all'):
        txt_report = generate_txt_report(result)
        if args.output_txt:
            with open(args.output_txt, 'w') as f:
                f.write(txt_report)
        elif args.format == 'txt':
            print(txt_report)
    
    if args.format in ('csv', 'all'):
        csv_data = generate_csv_data(result)
        if args.output_csv:
            with open(args.output_csv, 'w', newline='') as f:
                f.write(csv_data)
        elif args.format == 'csv':
            print(csv_data)


if __name__ == '__main__':
    main()
