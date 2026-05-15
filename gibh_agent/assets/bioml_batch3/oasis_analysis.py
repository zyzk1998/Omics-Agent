#!/usr/bin/env python3
"""
OASis Antibody Humanness Evaluation Script

Evaluates antibody humanness using peptide search in the Observed Antibody Space (OAS).
Supports multiple numbering schemes and CDR definitions.

Dependencies:
    - promb (pip install promb)
    - anarci (pip install anarci + hmmer)

Usage:
    python oasis_analysis.py --heavy-chain SEQ --light-chain SEQ \
        --name NAME --scheme kabat --cdr-definition kabat --threshold relaxed
"""

import argparse
import json
import sys
import os
from typing import Dict, List, Tuple, Optional

# Try to import promb and anarci
try:
    from promb import init_db
    PROMB_AVAILABLE = True
except ImportError:
    PROMB_AVAILABLE = False

try:
    from anarci import anarci
    ANARCI_AVAILABLE = True
except ImportError:
    ANARCI_AVAILABLE = False


# CDR definitions by numbering scheme
CDR_DEFINITIONS = {
    'kabat': {
        'heavy': {
            'CDR1': (31, 35),
            'CDR2': (50, 65),
            'CDR3': (95, 102),
        },
        'light': {
            'CDR1': (24, 34),
            'CDR2': (50, 56),
            'CDR3': (89, 97),
        }
    },
    'chothia': {
        'heavy': {
            'CDR1': (26, 32),
            'CDR2': (52, 56),
            'CDR3': (95, 102),
        },
        'light': {
            'CDR1': (24, 34),
            'CDR2': (50, 56),
            'CDR3': (89, 97),
        }
    },
    'imgt': {
        'heavy': {
            'CDR1': (27, 38),
            'CDR2': (56, 65),
            'CDR3': (105, 117),
        },
        'light': {
            'CDR1': (27, 38),
            'CDR2': (56, 65),
            'CDR3': (105, 117),
        }
    },
    'north': {
        'heavy': {
            'CDR1': (30, 35),
            'CDR2': (47, 58),
            'CDR3': (95, 102),
        },
        'light': {
            'CDR1': (24, 34),
            'CDR2': (50, 56),
            'CDR3': (89, 97),
        }
    }
}

THRESHOLDS = {
    'loose': '≥1% subjects',
    'relaxed': '≥10% subjects',
    'medium': '≥50% subjects',
    'strict': '≥90% subjects',
}


def number_antibody(sequence: str, scheme: str = 'kabat') -> Optional[List[Tuple]]:
    """Run ANARCI numbering on a single sequence."""
    if not ANARCI_AVAILABLE:
        return None
    
    results = anarci([('query', sequence)], scheme=scheme, output=False)
    if not results or len(results) < 3:
        return None
    
    # ANARCI returns: (numbering_list, chain_type_list, hit_table_list)
    # numbering_list = [ [ (numbered_seq, start, end) ], ... ]
    all_domains = results[0]
    if not all_domains or not all_domains[0]:
        return None
    
    # Get the first domain's results
    domain_results = all_domains[0]
    if not domain_results:
        return None
    
    # domain_results is a list of tuples: [(numbered_seq, start, end), ...]
    # Take the first one
    first_result = domain_results[0]
    if isinstance(first_result, (list, tuple)) and len(first_result) >= 1:
        numbered_seq = first_result[0]  # [((pos, insert), aa), ...]
        return numbered_seq
    
    return None


def extract_regions(numbered_seq: List[Tuple], chain_type: str, cdr_def: str) -> Dict[str, str]:
    """Extract CDR and framework regions from numbered sequence."""
    if cdr_def not in CDR_DEFINITIONS:
        cdr_def = 'kabat'
    
    definitions = CDR_DEFINITIONS[cdr_def]
    chain_key = 'heavy' if chain_type in ['H', 'heavy', 'Heavy'] else 'light'
    cdr_ranges = definitions.get(chain_key, {})
    
    regions = {
        'framework': [],
        'CDR1': [],
        'CDR2': [],
        'CDR3': [],
    }
    
    cdr_positions = {}
    for cdr_name, (start, end) in cdr_ranges.items():
        cdr_positions[cdr_name] = set(range(start, end + 1))
    
    for item in numbered_seq:
        # item format: ((position, insertion_code), amino_acid)
        if not isinstance(item, (list, tuple)) or len(item) < 2:
            continue
        
        numbering = item[0]
        aa = item[1]
        
        if isinstance(numbering, tuple) and len(numbering) >= 1:
            pos = numbering[0]  # position number
        elif isinstance(numbering, int):
            pos = numbering
        else:
            continue
        
        assigned = False
        for cdr_name, positions in cdr_positions.items():
            if pos in positions:
                regions[cdr_name].append(aa)
                assigned = True
                break
        if not assigned:
            regions['framework'].append(aa)
    
    return {
        k: ''.join(v) for k, v in regions.items()
    }


def compute_oasis_identity(sequence: str, db) -> float:
    """Compute OASis identity score for a sequence."""
    if not sequence or len(sequence) < 9:
        return 0.0
    return db.compute_peptide_content(sequence)


def run_analysis(
    heavy_chain: str,
    light_chain: str,
    name: str = 'Antibody',
    scheme: str = 'kabat',
    cdr_definition: str = 'kabat',
    threshold: str = 'relaxed',
    custom_threshold_pct: Optional[float] = None,
) -> Dict:
    """Run complete OASis humanness analysis."""
    
    if not PROMB_AVAILABLE:
        raise RuntimeError("promb is not installed. Install with: pip install promb")
    
    # Initialize promb database
    # Note: promb's built-in 'human-oas' corresponds to relaxed (>=10% subjects)
    # For other thresholds, we'd need the full OASis database (not available in promb)
    db = init_db('human-oas')
    
    results = {
        'antibody_name': name,
        'scheme': scheme,
        'cdr_definition': cdr_definition,
        'threshold': threshold,
        'threshold_description': THRESHOLDS.get(threshold, 'custom'),
        'warning': None,
        'heavy_chain': {},
        'light_chain': {},
        'overall': {},
    }
    
    # Warn if threshold other than relaxed is requested
    if threshold != 'relaxed' and custom_threshold_pct is None:
        results['warning'] = (
            f"Threshold '{threshold}' requires the full 22GB OASis database. "
            "Using promb's built-in human-oas database (relaxed, >=10% subjects). "
            "Install BioPhi with the full database for other thresholds."
        )
    
    # Process heavy chain
    heavy_numbered = number_antibody(heavy_chain, scheme) if ANARCI_AVAILABLE else None
    heavy_regions = None
    
    if heavy_numbered and ANARCI_AVAILABLE:
        heavy_regions = extract_regions(heavy_numbered, 'heavy', cdr_definition)
        
        results['heavy_chain']['full_sequence'] = heavy_chain
        results['heavy_chain']['length'] = len(heavy_chain)
        results['heavy_chain']['oasis_identity'] = compute_oasis_identity(heavy_chain, db)
        
        if heavy_regions:
            results['heavy_chain']['regions'] = {}
            for region_name, seq in heavy_regions.items():
                if seq:
                    score = compute_oasis_identity(seq, db)
                    results['heavy_chain']['regions'][region_name] = {
                        'sequence': seq,
                        'length': len(seq),
                        'oasis_identity': score,
                    }
    else:
        # Fallback without ANARCI
        results['heavy_chain']['full_sequence'] = heavy_chain
        results['heavy_chain']['length'] = len(heavy_chain)
        results['heavy_chain']['oasis_identity'] = compute_oasis_identity(heavy_chain, db)
        results['heavy_chain']['regions'] = None
        if not ANARCI_AVAILABLE:
            results['heavy_chain']['note'] = 'ANARCI not available. Install with: pip install anarci (requires hmmer)'
    
    # Process light chain
    light_numbered = number_antibody(light_chain, scheme) if ANARCI_AVAILABLE else None
    light_regions = None
    
    if light_numbered and ANARCI_AVAILABLE:
        light_regions = extract_regions(light_numbered, 'light', cdr_definition)
        
        results['light_chain']['full_sequence'] = light_chain
        results['light_chain']['length'] = len(light_chain)
        results['light_chain']['oasis_identity'] = compute_oasis_identity(light_chain, db)
        
        if light_regions:
            results['light_chain']['regions'] = {}
            for region_name, seq in light_regions.items():
                if seq:
                    score = compute_oasis_identity(seq, db)
                    results['light_chain']['regions'][region_name] = {
                        'sequence': seq,
                        'length': len(seq),
                        'oasis_identity': score,
                    }
    else:
        results['light_chain']['full_sequence'] = light_chain
        results['light_chain']['length'] = len(light_chain)
        results['light_chain']['oasis_identity'] = compute_oasis_identity(light_chain, db)
        results['light_chain']['regions'] = None
        if not ANARCI_AVAILABLE:
            results['light_chain']['note'] = 'ANARCI not available. Install with: pip install anarci (requires hmmer)'
    
    # Compute overall scores
    combined_seq = heavy_chain + light_chain
    results['overall']['combined_length'] = len(combined_seq)
    results['overall']['oasis_identity'] = compute_oasis_identity(combined_seq, db)
    
    # Average of heavy and light (more representative than combined)
    heavy_score = results['heavy_chain']['oasis_identity']
    light_score = results['light_chain']['oasis_identity']
    results['overall']['average_identity'] = (heavy_score + light_score) / 2.0
    
    # CDR vs Framework comparison
    if heavy_regions and light_regions and ANARCI_AVAILABLE:
        all_cdr_seqs = []
        all_fw_seqs = []
        
        for chain_regions in [heavy_regions, light_regions]:
            for region_name, seq in chain_regions.items():
                if region_name.startswith('CDR'):
                    all_cdr_seqs.append(seq)
                elif region_name == 'framework':
                    all_fw_seqs.append(seq)
        
        cdr_combined = ''.join(all_cdr_seqs)
        fw_combined = ''.join(all_fw_seqs)
        
        if cdr_combined:
            results['overall']['cdr_identity'] = compute_oasis_identity(cdr_combined, db)
            results['overall']['cdr_length'] = len(cdr_combined)
        if fw_combined:
            results['overall']['framework_identity'] = compute_oasis_identity(fw_combined, db)
            results['overall']['framework_length'] = len(fw_combined)
    
    return results


def format_text_report(data: Dict) -> str:
    """Format results as human-readable text."""
    lines = []
    
    lines.append("=" * 60)
    lines.append(f"OASis Humanness Evaluation Report")
    lines.append(f"Antibody: {data['antibody_name']}")
    lines.append(f"Scheme: {data['scheme']} | CDR Definition: {data['cdr_definition']}")
    lines.append(f"Threshold: {data['threshold']} ({data['threshold_description']})")
    lines.append("=" * 60)
    
    if data.get('warning'):
        lines.append(f"\n⚠️  WARNING: {data['warning']}\n")
    
    # Heavy chain
    lines.append("\n--- Heavy Chain ---")
    hc = data['heavy_chain']
    lines.append(f"Length: {hc['length']} residues")
    lines.append(f"OASis Identity: {hc['oasis_identity']:.4f} ({hc['oasis_identity']*100:.2f}%)")
    
    if hc.get('regions'):
        lines.append("\nRegion breakdown:")
        for region_name, region_data in hc['regions'].items():
            lines.append(f"  {region_name:12s}: {region_data['oasis_identity']:.4f} ({region_data['oasis_identity']*100:.2f}%) [{region_data['length']} aa]")
    elif hc.get('note'):
        lines.append(f"Note: {hc['note']}")
    
    # Light chain
    lines.append("\n--- Light Chain ---")
    lc = data['light_chain']
    lines.append(f"Length: {lc['length']} residues")
    lines.append(f"OASis Identity: {lc['oasis_identity']:.4f} ({lc['oasis_identity']*100:.2f}%)")
    
    if lc.get('regions'):
        lines.append("\nRegion breakdown:")
        for region_name, region_data in lc['regions'].items():
            lines.append(f"  {region_name:12s}: {region_data['oasis_identity']:.4f} ({region_data['oasis_identity']*100:.2f}%) [{region_data['length']} aa]")
    elif lc.get('note'):
        lines.append(f"Note: {lc['note']}")
    
    # Overall
    lines.append("\n--- Overall ---")
    ov = data['overall']
    lines.append(f"Combined Length: {ov['combined_length']} residues")
    lines.append(f"Overall Identity: {ov['oasis_identity']:.4f} ({ov['oasis_identity']*100:.2f}%)")
    lines.append(f"Average (H+L):    {ov['average_identity']:.4f} ({ov['average_identity']*100:.2f}%)")
    
    if 'cdr_identity' in ov:
        lines.append(f"\nCDR Regions:     {ov['cdr_identity']:.4f} ({ov['cdr_identity']*100:.2f}%) [{ov.get('cdr_length', 'N/A')} aa]")
    if 'framework_identity' in ov:
        lines.append(f"Framework:       {ov['framework_identity']:.4f} ({ov['framework_identity']*100:.2f}%) [{ov.get('framework_length', 'N/A')} aa]")
    
    lines.append("\n" + "=" * 60)
    lines.append("Interpretation:")
    score = ov['average_identity']
    if score >= 0.85:
        lines.append("  ✅ High humanness - Low predicted immunogenicity risk")
    elif score >= 0.65:
        lines.append("  ⚠️  Moderate humanness - Consider humanization")
    else:
        lines.append("  ❌ Low humanness - High predicted immunogenicity risk")
    lines.append("=" * 60)
    
    return '\n'.join(lines)


def main():
    parser = argparse.ArgumentParser(description='OASis Antibody Humanness Evaluation')
    parser.add_argument('--heavy-chain', required=True, help='Heavy chain amino acid sequence')
    parser.add_argument('--light-chain', required=True, help='Light chain amino acid sequence')
    parser.add_argument('--name', default='Antibody', help='Antibody name')
    parser.add_argument('--scheme', default='kabat', choices=['kabat', 'chothia', 'imgt', 'aho', 'martin'],
                        help='Numbering scheme')
    parser.add_argument('--cdr-definition', default='kabat', choices=['kabat', 'chothia', 'imgt', 'north'],
                        help='CDR definition')
    parser.add_argument('--threshold', default='relaxed', 
                        choices=['loose', 'relaxed', 'medium', 'strict'],
                        help='OASis threshold')
    parser.add_argument('--custom-threshold-pct', type=float, default=None,
                        help='Custom threshold percentage (0-100)')
    parser.add_argument('--output-format', default='json', choices=['json', 'text'],
                        help='Output format')
    
    args = parser.parse_args()
    
    try:
        results = run_analysis(
            heavy_chain=args.heavy_chain,
            light_chain=args.light_chain,
            name=args.name,
            scheme=args.scheme,
            cdr_definition=args.cdr_definition,
            threshold=args.threshold,
            custom_threshold_pct=args.custom_threshold_pct,
        )
        
        if args.output_format == 'json':
            print(json.dumps(results, indent=2))
        else:
            print(format_text_report(results))
            
    except Exception as e:
        import traceback
        traceback.print_exc()
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()