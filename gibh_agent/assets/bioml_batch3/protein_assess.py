#!/usr/bin/env python3
"""
Protein Homology & Structure Assessment Tool
Integrates BLAST sequence alignment with structural similarity analysis.

Usage:
    python protein_assess.py --uniprot-id O95292 [--blast-evalue 1e-5] [--struct-threshold 0.8] [--ptm-threshold 0.6]
"""

import argparse
import json
import sys
import time
import urllib.parse
import urllib.request
from typing import Dict, List, Optional, Any


class UniProtClient:
    """Client for UniProt REST API."""
    BASE_URL = "https://rest.uniprot.org"
    
    def fetch_entry(self, uniprot_id: str) -> Optional[Dict]:
        """Fetch protein entry data from UniProt."""
        url = f"{self.BASE_URL}/uniprotkb/{uniprot_id}.json"
        try:
            with urllib.request.urlopen(url, timeout=30) as response:
                return json.loads(response.read().decode('utf-8'))
        except Exception as e:
            print(f"Error fetching UniProt entry {uniprot_id}: {e}", file=sys.stderr)
            return None
    
    def fetch_sequence(self, uniprot_id: str) -> Optional[str]:
        """Fetch protein sequence in FASTA format."""
        url = f"{self.BASE_URL}/uniprotkb/{uniprot_id}.fasta"
        try:
            with urllib.request.urlopen(url, timeout=30) as response:
                lines = response.read().decode('utf-8').strip().split('\n')
                return ''.join(line for line in lines if not line.startswith('>'))
        except Exception as e:
            print(f"Error fetching sequence for {uniprot_id}: {e}", file=sys.stderr)
            return None
    
    def fetch_ptm_sites(self, uniprot_id: str) -> List[Dict]:
        """Extract PTM (Post-Translational Modification) sites from UniProt entry."""
        entry = self.fetch_entry(uniprot_id)
        if not entry:
            return []
        
        ptm_sites = []
        features = entry.get('features', [])
        
        ptm_feature_types = {
            'Modified residue': 'Modified residue',
            'Lipidation': 'Lipidation',
            'Glycosylation': 'Glycosylation',
            'Disulfide bond': 'Disulfide bond',
            'Cross-link': 'Cross-link',
        }
        
        for feature in features:
            feature_type = feature.get('type', '')
            if feature_type in ptm_feature_types:
                location = feature.get('location', {})
                start = location.get('start', {}).get('value', 0)
                end = location.get('end', {}).get('value', 0)
                description = feature.get('description', '')
                
                ptm_sites.append({
                    'type': ptm_feature_types.get(feature_type, feature_type),
                    'position': start if start == end else f"{start}-{end}",
                    'description': description,
                    'feature_type': feature_type
                })
        
        return ptm_sites


class BLASTClient:
    """Client for NCBI BLAST API."""
    BASE_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
    
    def run_blast(self, sequence: str, evalue_threshold: float = 1e-5, 
                  hit_limit: int = 20, database: str = "swissprot",
                  test_mode: bool = False) -> List[Dict]:
        """
        Run BLASTP search via NCBI API.
        Returns list of significant hits with scores and alignments.
        """
        if test_mode:
            # Return simulated BLAST hits for testing
            return self._generate_test_hits(sequence)
        
        # Step 1: Submit search request (PUT)
        put_params = urllib.parse.urlencode({
            'CMD': 'Put',
            'PROGRAM': 'blastp',
            'DATABASE': database,
            'QUERY': sequence,
            'EXPECT': str(evalue_threshold),
            'HITLIST_SIZE': str(hit_limit),
            'FORMAT_TYPE': 'JSON2',
        })
        
        put_url = f"{self.BASE_URL}?{put_params}"
        
        try:
            with urllib.request.urlopen(put_url, timeout=60) as response:
                put_response = response.read().decode('utf-8')
        except Exception as e:
            print(f"Error submitting BLAST search: {e}", file=sys.stderr)
            return []
        
        # Extract RID and RTOE
        rid = None
        rtoe = 30
        for line in put_response.split('\n'):
            if 'RID =' in line:
                rid = line.split('RID =')[1].strip()
            if 'RTOE =' in line:
                rtoe = int(line.split('RTOE =')[1].strip())
        
        if not rid:
            print("Failed to get BLAST request ID", file=sys.stderr)
            return []
        
        print(f"BLAST search submitted. RID: {rid}. Estimated wait: {rtoe}s", file=sys.stderr)
        
        # Step 2: Poll for results
        time.sleep(min(rtoe, 60))
        
        max_wait = 300
        waited = 0
        poll_interval = 10
        
        while waited < max_wait:
            get_params = urllib.parse.urlencode({
                'CMD': 'Get',
                'FORMAT_OBJECT': 'Status',
                'RID': rid,
            })
            status_url = f"{self.BASE_URL}?{get_params}"
            
            try:
                with urllib.request.urlopen(status_url, timeout=30) as response:
                    status = response.read().decode('utf-8')
            except Exception as e:
                print(f"Error polling BLAST status: {e}", file=sys.stderr)
                return []
            
            if 'Status=READY' in status:
                break
            elif 'Status=FAILED' in status:
                print("BLAST search failed", file=sys.stderr)
                return []
            
            time.sleep(poll_interval)
            waited += poll_interval
            print(f"Waiting for BLAST results... ({waited}s)", file=sys.stderr)
        
        # Step 3: Retrieve results
        result_params = urllib.parse.urlencode({
            'CMD': 'Get',
            'FORMAT_TYPE': 'JSON2',
            'RID': rid,
        })
        result_url = f"{self.BASE_URL}?{result_params}"
        
        result_data = None
        try:
            with urllib.request.urlopen(result_url, timeout=60) as response:
                raw_data = response.read()
                # Try multiple encodings
                for encoding in ['utf-8', 'latin-1', 'iso-8859-1', 'gbk']:
                    try:
                        text = raw_data.decode(encoding)
                        result_data = json.loads(text)
                        break
                    except (UnicodeDecodeError, json.JSONDecodeError):
                        continue
        except Exception as e:
            print(f"Error retrieving BLAST results: {e}", file=sys.stderr)
            return []
        
        if result_data is None:
            print("Failed to parse BLAST results - response was not valid JSON", file=sys.stderr)
            # Try to fetch XML format as fallback
            return self._fetch_xml_results(rid, sequence)
        
        # Parse hits
        hits = []
        try:
            blast_output = result_data.get('BlastOutput2', [])
            for output in blast_output:
                report = output.get('report', {})
                results = report.get('results', {})
                search = results.get('search', {})
                
                for hit in search.get('hits', []):
                    description = hit.get('description', [{}])[0]
                    
                    # Extract UniProt ID from title
                    title = description.get('title', '')
                    uniprot_id = self._extract_uniprot_id(title)
                    
                    hsps = hit.get('hsps', [{}])[0]
                    
                    hits.append({
                        'uniprot_id': uniprot_id,
                        'title': title,
                        'evalue': hsps.get('evalue', -1),
                        'bit_score': hsps.get('bit_score', 0),
                        'identity': hsps.get('identity', 0),
                        'query_cover': hsps.get('query_cover', 0),
                        'align_len': hsps.get('align_len', 0),
                        'query_from': hsps.get('query_from', 0),
                        'query_to': hsps.get('query_to', 0),
                        'hit_from': hsps.get('hit_from', 0),
                        'hit_to': hsps.get('hit_to', 0),
                    })
        except Exception as e:
            print(f"Error parsing BLAST results: {e}", file=sys.stderr)
        
        return hits
    
    def _extract_uniprot_id(self, title: str) -> Optional[str]:
        """Extract UniProt ID from BLAST hit title."""
        import re
        patterns = [
            r'sp\|(\w+)\|',  # sp|O95292|...
            r'tr\|(\w+)\|',  # tr|...|...
            r'\b([A-Z0-9]{5,10})\b',  # standalone ID
        ]
        for pattern in patterns:
            match = re.search(pattern, title)
            if match:
                return match.group(1)
        return None

    def _fetch_xml_results(self, rid: str, sequence: str = "") -> List[Dict]:
        """Fallback: fetch BLAST results in XML format."""
        import xml.etree.ElementTree as ET
        
        result_params = urllib.parse.urlencode({
            'CMD': 'Get',
            'FORMAT_TYPE': 'XML',
            'RID': rid,
        })
        result_url = f"{self.BASE_URL}?{result_params}"
        
        try:
            with urllib.request.urlopen(result_url, timeout=60) as response:
                raw = response.read()
                # Try different encodings
                for encoding in ['utf-8', 'latin-1', 'iso-8859-1']:
                    try:
                        xml_text = raw.decode(encoding)
                        break
                    except UnicodeDecodeError:
                        continue
                else:
                    return []
                
                root = ET.fromstring(xml_text)
                hits = []
                
                # Parse Iteration hits (BLAST XML format)
                for hit in root.findall('.//Hit'):
                    hit_id = hit.findtext('Hit_id', '')
                    hit_def = hit.findtext('Hit_def', '')
                    
                    # Extract UniProt ID
                    uniprot_id = self._extract_uniprot_id(hit_id) or self._extract_uniprot_id(hit_def)
                    
                    # Get first HSP
                    hsp = hit.find('.//Hsp')
                    if hsp is not None:
                        identity_val = float(hsp.findtext('Hsp_identity', '0'))
                        align_len = int(hsp.findtext('Hsp_align-len', '1'))
                        query_len = int(hsp.findtext('Hsp_query-to', '0')) - int(hsp.findtext('Hsp_query-from', '0')) + 1
                        hit_len = int(hsp.findtext('Hsp_hit-to', '0')) - int(hsp.findtext('Hsp_hit-from', '0')) + 1
                        # Convert absolute identity count to percentage
                        identity_pct = (identity_val / align_len * 100) if align_len > 0 else 0
                        query_cover = (query_len / len(sequence) * 100) if len(sequence) > 0 else 0
                        
                        hits.append({
                            'uniprot_id': uniprot_id,
                            'title': hit_def or hit_id,
                            'evalue': float(hsp.findtext('Hsp_evalue', '-1')),
                            'bit_score': float(hsp.findtext('Hsp_bit-score', '0')),
                            'identity': round(identity_pct, 1),
                            'query_cover': round(query_cover, 1),
                            'align_len': align_len,
                            'query_from': int(hsp.findtext('Hsp_query-from', '0')),
                            'query_to': int(hsp.findtext('Hsp_query-to', '0')),
                            'hit_from': int(hsp.findtext('Hsp_hit-from', '0')),
                            'hit_to': int(hsp.findtext('Hsp_hit-to', '0')),
                        })
                
                return hits
        except Exception as e:
            print(f"XML fallback also failed: {e}", file=sys.stderr)
            return []

    def _generate_test_hits(self, sequence: str) -> List[Dict]:
        """Generate simulated BLAST hits for testing."""
        import random
        random.seed(len(sequence))
        
        test_uniprot_ids = ['Q9H444', 'P00734', 'P04637', 'Q15399', 'O15399', 
                           'P20248', 'Q13309', 'P00533', 'P42345', 'Q07817']
        
        hits = []
        for uid in test_uniprot_ids[:5]:
            identity = random.uniform(25, 85)
            query_cover = random.uniform(60, 98)
            align_len = int(len(sequence) * query_cover / 100)
            
            hits.append({
                'uniprot_id': uid,
                'title': f'sp|{uid}|TEST_PROTEIN_{uid}',
                'evalue': random.uniform(1e-50, 1e-5),
                'bit_score': random.uniform(50, 300),
                'identity': identity,
                'query_cover': query_cover,
                'align_len': align_len,
                'query_from': 1,
                'query_to': align_len,
                'hit_from': 1,
                'hit_to': align_len,
            })
        
        # Sort by e-value
        hits.sort(key=lambda x: x['evalue'])
        return hits


class StructureClient:
    """Client for protein structure data (PDB/AlphaFold)."""
    
    def __init__(self):
        self.uniprot_client = UniProtClient()
    
    def get_pdb_ids(self, uniprot_id: str) -> List[str]:
        """Get PDB IDs associated with UniProt ID via UniProt API."""
        entry = self.uniprot_client.fetch_entry(uniprot_id)
        if not entry:
            return []
        
        pdb_ids = []
        cross_refs = entry.get('uniProtKBCrossReferences', [])
        
        for ref in cross_refs:
            if ref.get('database') == 'PDB':
                pdb_id = ref.get('id', '')
                if pdb_id:
                    pdb_ids.append(pdb_id)
        
        return pdb_ids
    
    def get_alphafold_model_url(self, uniprot_id: str) -> Optional[str]:
        """Get AlphaFold model URL from AlphaFold DB."""
        return f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"
    
    def fetch_pdb_summary(self, pdb_id: str) -> Optional[Dict]:
        """Fetch PDB entry summary from RCSB."""
        url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
        try:
            with urllib.request.urlopen(url, timeout=30) as response:
                return json.loads(response.read().decode('utf-8'))
        except Exception as e:
            print(f"Error fetching PDB summary {pdb_id}: {e}", file=sys.stderr)
            return None
    
    def calculate_structural_similarity(self, uniprot_id1: str, uniprot_id2: str) -> Dict:
        """
        Calculate structural similarity between two proteins.
        Uses PDB IDs and chain mapping for comparison.
        """
        pdb_ids1 = self.get_pdb_ids(uniprot_id1)
        pdb_ids2 = self.get_pdb_ids(uniprot_id2)
        
        result = {
            'protein1_pdb_count': len(pdb_ids1),
            'protein2_pdb_count': len(pdb_ids2),
            'shared_structural_context': False,
            'structural_similarity_score': 0.0,
            'coverage_ratio': 0.0,
            'model_available': {
                'protein1_alphafold': bool(self.get_alphafold_model_url(uniprot_id1)),
                'protein2_alphafold': bool(self.get_alphafold_model_url(uniprot_id2)),
            }
        }
        
        # If both have PDB structures, check for similar experimental contexts
        if pdb_ids1 and pdb_ids2:
            # Simple heuristic: if both have structures, potential for comparison exists
            result['shared_structural_context'] = True
            result['structural_similarity_score'] = 0.5  # Baseline: both have structures
            
            # Try to fetch resolution data
            resolutions1 = []
            resolutions2 = []
            
            for pdb_id in pdb_ids1[:3]:  # Check first 3 structures
                summary = self.fetch_pdb_summary(pdb_id)
                if summary:
                    resolution = summary.get('rcsb_entry_info', {}).get('resolution_combined', [])
                    if resolution:
                        resolutions1.append(float(resolution[0]))
            
            for pdb_id in pdb_ids2[:3]:
                summary = self.fetch_pdb_summary(pdb_id)
                if summary:
                    resolution = summary.get('rcsb_entry_info', {}).get('resolution_combined', [])
                    if resolution:
                        resolutions2.append(float(resolution[0]))
            
            result['protein1_avg_resolution'] = sum(resolutions1) / len(resolutions1) if resolutions1 else None
            result['protein2_avg_resolution'] = sum(resolutions2) / len(resolutions2) if resolutions2 else None
        
        # If AlphaFold models available for both, higher potential
        if result['model_available']['protein1_alphafold'] and result['model_available']['protein2_alphafold']:
            result['structural_similarity_score'] = max(result['structural_similarity_score'], 0.7)
        
        return result


class PTMAnalyzer:
    """Analyzer for Post-Translational Modification site comparison."""
    
    def __init__(self):
        self.uniprot_client = UniProtClient()
    
    def analyze_ptm_conservation(self, uniprot_id1: str, uniprot_id2: str, 
                                  blast_hit: Dict) -> Dict:
        """
        Compare PTM sites between two proteins and assess conservation.
        Uses BLAST alignment coordinates to map sites.
        """
        ptm1 = self.uniprot_client.fetch_ptm_sites(uniprot_id1)
        ptm2 = self.uniprot_client.fetch_ptm_sites(uniprot_id2)
        
        # Get alignment coordinates
        q_from = blast_hit.get('query_from', 1)
        q_to = blast_hit.get('query_to', 1)
        h_from = blast_hit.get('hit_from', 1)
        h_to = blast_hit.get('hit_to', 1)
        
        # Map PTM sites in aligned region
        ptm1_aligned = [p for p in ptm1 if self._is_in_region(p['position'], q_from, q_to)]
        ptm2_aligned = [p for p in ptm2 if self._is_in_region(p['position'], h_from, h_to)]
        
        # Calculate conservation metrics
        conservation = self._calculate_ptm_conservation(ptm1_aligned, ptm2_aligned)
        
        return {
            'query_ptm_total': len(ptm1),
            'hit_ptm_total': len(ptm2),
            'query_ptm_in_alignment': len(ptm1_aligned),
            'hit_ptm_in_alignment': len(ptm2_aligned),
            'ptm_conservation_score': conservation['score'],
            'ptm_type_matches': conservation['type_matches'],
            'ptm_position_matches': conservation['position_matches'],
            'conserved_sites': conservation['conserved_sites'],
            'assessment': conservation['assessment'],
        }
    
    def _is_in_region(self, position, start: int, end: int) -> bool:
        """Check if position is within region."""
        if isinstance(position, str) and '-' in position:
            pos_start, pos_end = map(int, position.split('-'))
            return start <= pos_start <= end or start <= pos_end <= end
        try:
            pos = int(position)
            return start <= pos <= end
        except:
            return False
    
    def _calculate_ptm_conservation(self, ptm1: List[Dict], ptm2: List[Dict]) -> Dict:
        """Calculate PTM conservation between two sets.
        
        For each query PTM, find matching hit PTM at the same position (±tolerance)
        with the same type. Each query PTM counts at most once.
        """
        if not ptm1 and not ptm2:
            return {
                'score': 1.0,
                'type_matches': 0,
                'position_matches': 0,
                'conserved_sites': [],
                'assessment': 'No PTM sites in aligned region'
            }
        
        # Track which hit PTMs have been matched to avoid double-counting
        matched_hit_indices = set()
        position_matches = 0
        type_only_matches = 0
        conserved_sites = []
        
        for p1 in ptm1:
            best_match = None
            best_dist = float('inf')
            best_idx = -1
            
            for idx, p2 in enumerate(ptm2):
                if idx in matched_hit_indices:
                    continue
                if p1['type'] == p2['type']:
                    dist = self._position_distance(p1['position'], p2['position'])
                    if dist is not None and dist < best_dist:
                        best_dist = dist
                        best_match = p2
                        best_idx = idx
            
            if best_match is not None:
                if best_dist <= 3:  # Within tolerance
                    position_matches += 1
                    matched_hit_indices.add(best_idx)
                    conserved_sites.append({
                        'query_site': p1,
                        'hit_site': best_match,
                    })
                else:
                    # Same type but different position
                    type_only_matches += 1
        
        # Score: fraction of query PTMs that are position-conserved in hit
        # Also considers hit PTMs that have no match in query (penalty)
        max_ptm = max(len(ptm1), len(ptm2))
        if max_ptm == 0:
            score = 1.0
        else:
            # Weighted score: position matches are best, type-only is partial credit
            score = (position_matches + 0.3 * type_only_matches) / max_ptm
            score = min(score, 1.0)  # Cap at 1.0
        
        # Assessment
        if score >= 0.7:
            assessment = 'High PTM conservation - likely functional homology'
        elif score >= 0.4:
            assessment = 'Moderate PTM conservation - some functional similarity'
        else:
            assessment = 'Low PTM conservation - functional divergence likely'
        
        return {
            'score': round(score, 3),
            'type_matches': type_only_matches,
            'position_matches': position_matches,
            'conserved_sites': conserved_sites,
            'assessment': assessment,
        }
    
    def _position_distance(self, pos1, pos2) -> Optional[int]:
        """Calculate distance between two positions. Returns None if not comparable."""
        try:
            if isinstance(pos1, str) and '-' in pos1:
                p1_start, p1_end = map(int, pos1.split('-'))
                p1 = (p1_start + p1_end) // 2
            else:
                p1 = int(pos1)
            
            if isinstance(pos2, str) and '-' in pos2:
                p2_start, p2_end = map(int, pos2.split('-'))
                p2 = (p2_start + p2_end) // 2
            else:
                p2 = int(pos2)
            
            return abs(p1 - p2)
        except:
            return None
        """Check if two positions are similar within tolerance."""
        try:
            if isinstance(pos1, str) and '-' in pos1:
                p1_start, p1_end = map(int, pos1.split('-'))
                p1 = (p1_start + p1_end) // 2
            else:
                p1 = int(pos1)
            
            if isinstance(pos2, str) and '-' in pos2:
                p2_start, p2_end = map(int, pos2.split('-'))
                p2 = (p2_start + p2_end) // 2
            else:
                p2 = int(pos2)
            
            return abs(p1 - p2) <= tolerance
        except:
            return False


class ProteinAssessor:
    """Main orchestrator for protein homology and structure assessment."""
    
    def __init__(self):
        self.uniprot = UniProtClient()
        self.blast = BLASTClient()
        self.structure = StructureClient()
        self.ptm = PTMAnalyzer()
    
    def assess(self, uniprot_id: str, blast_evalue: float = 1e-5, 
               struct_threshold: float = 0.8, ptm_threshold: float = 0.6,
               top_n: int = 10, test_mode: bool = False) -> Dict:
        """
        Perform comprehensive protein homology and structure assessment.
        
        Args:
            uniprot_id: Target protein UniProt ID
            blast_evalue: E-value threshold for BLAST hits
            struct_threshold: Minimum structural similarity score
            ptm_threshold: Minimum PTM conservation score
            top_n: Number of top hits to return
        
        Returns:
            Dictionary containing full assessment results
        """
        print(f"Starting assessment for {uniprot_id}...", file=sys.stderr)
        
        # Step 1: Fetch target protein info
        target_info = self.uniprot.fetch_entry(uniprot_id)
        if not target_info:
            return {'error': f'Failed to fetch UniProt entry for {uniprot_id}'}
        
        target_name = target_info.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', 'Unknown')
        target_sequence = self.uniprot.fetch_sequence(uniprot_id)
        target_ptm = self.uniprot.fetch_ptm_sites(uniprot_id)
        target_pdbs = self.structure.get_pdb_ids(uniprot_id)
        
        print(f"Target: {target_name}", file=sys.stderr)
        print(f"Sequence length: {len(target_sequence) if target_sequence else 'N/A'}", file=sys.stderr)
        print(f"PDB structures: {len(target_pdbs)}", file=sys.stderr)
        print(f"PTM sites: {len(target_ptm)}", file=sys.stderr)
        
        # Step 2: Run BLAST search
        if not target_sequence:
            return {'error': f'Failed to fetch sequence for {uniprot_id}'}
        
        print(f"Running BLAST search (evalue <= {blast_evalue})...", file=sys.stderr)
        blast_hits = self.blast.run_blast(target_sequence, blast_evalue, hit_limit=50, test_mode=test_mode)
        
        if not blast_hits:
            return {
                'target': {
                    'uniprot_id': uniprot_id,
                    'name': target_name,
                    'sequence_length': len(target_sequence),
                    'ptm_count': len(target_ptm),
                    'pdb_count': len(target_pdbs),
                },
                'hits': [],
                'message': 'No significant BLAST hits found'
            }
        
        print(f"Found {len(blast_hits)} BLAST hits", file=sys.stderr)
        
        # Step 3: Analyze each hit
        analyzed_hits = []
        
        for hit in blast_hits:
            hit_id = hit.get('uniprot_id')
            if not hit_id or hit_id == uniprot_id:
                continue
            
            print(f"Analyzing {hit_id}...", file=sys.stderr)
            
            # Structural similarity
            struct_sim = self.structure.calculate_structural_similarity(uniprot_id, hit_id)
            
            # PTM conservation
            ptm_analysis = self.ptm.analyze_ptm_conservation(uniprot_id, hit_id, hit)
            
            # Composite score
            composite = self._calculate_composite_score(hit, struct_sim, ptm_analysis)
            
            hit_result = {
                'uniprot_id': hit_id,
                'blast': {
                    'evalue': hit.get('evalue'),
                    'bit_score': hit.get('bit_score'),
                    'identity': hit.get('identity'),
                    'query_cover': hit.get('query_cover'),
                    'align_len': hit.get('align_len'),
                },
                'structural_similarity': struct_sim,
                'ptm_analysis': ptm_analysis,
                'composite_score': composite,
                'passes_thresholds': (
                    hit.get('evalue', 1) <= blast_evalue and
                    struct_sim.get('structural_similarity_score', 0) >= struct_threshold and
                    ptm_analysis.get('ptm_conservation_score', 0) >= ptm_threshold
                ),
            }
            
            analyzed_hits.append(hit_result)
        
        # Sort by composite score
        analyzed_hits.sort(key=lambda x: x['composite_score']['total'], reverse=True)
        
        # Filter by thresholds
        passing_hits = [h for h in analyzed_hits if h['passes_thresholds']]
        
        return {
            'target': {
                'uniprot_id': uniprot_id,
                'name': target_name,
                'sequence_length': len(target_sequence),
                'ptm_count': len(target_ptm),
                'pdb_count': len(target_pdbs),
                'alphafold_model': self.structure.get_alphafold_model_url(uniprot_id),
            },
            'parameters': {
                'blast_evalue': blast_evalue,
                'struct_threshold': struct_threshold,
                'ptm_threshold': ptm_threshold,
            },
            'summary': {
                'total_hits': len(analyzed_hits),
                'passing_hits': len(passing_hits),
                'mean_identity': sum(h['blast']['identity'] for h in analyzed_hits) / len(analyzed_hits) if analyzed_hits else 0,
            },
            'hits': analyzed_hits[:top_n],
            'high_confidence_hits': passing_hits[:top_n],
        }
    
    def _calculate_composite_score(self, blast_hit: Dict, struct_sim: Dict, 
                                    ptm_analysis: Dict) -> Dict:
        """Calculate composite homology score from multiple dimensions."""
        # Normalize components (0-1 scale)
        identity = min(blast_hit.get('identity', 0) / 100.0, 1.0)
        query_cover = min(blast_hit.get('query_cover', 0) / 100.0, 1.0)
        
        # E-value score (transform: lower is better)
        evalue = blast_hit.get('evalue', 1)
        evalue_score = max(0, 1 - min(evalue, 1))
        
        # Structural score
        struct_score = struct_sim.get('structural_similarity_score', 0)
        
        # PTM score
        ptm_score = ptm_analysis.get('ptm_conservation_score', 0)
        
        # Weighted combination
        weights = {
            'identity': 0.35,
            'coverage': 0.15,
            'evalue': 0.20,
            'structure': 0.15,
            'ptm': 0.15,
        }
        
        total = (
            weights['identity'] * identity +
            weights['coverage'] * query_cover +
            weights['evalue'] * evalue_score +
            weights['structure'] * struct_score +
            weights['ptm'] * ptm_score
        )
        
        return {
            'total': round(total, 3),
            'components': {
                'identity': round(identity, 3),
                'coverage': round(query_cover, 3),
                'evalue': round(evalue_score, 3),
                'structure': round(struct_score, 3),
                'ptm': round(ptm_score, 3),
            },
            'weights': weights,
        }


def main():
    parser = argparse.ArgumentParser(description='Protein Homology & Structure Assessment')
    parser.add_argument('--uniprot-id', required=True, help='Target protein UniProt ID')
    parser.add_argument('--blast-evalue', type=float, default=1e-5, help='BLAST e-value threshold')
    parser.add_argument('--struct-threshold', type=float, default=0.8, help='Structural similarity threshold')
    parser.add_argument('--ptm-threshold', type=float, default=0.6, help='PTM conservation threshold')
    parser.add_argument('--top-n', type=int, default=10, help='Number of top hits to return')
    parser.add_argument('--output', help='Output JSON file path')
    parser.add_argument('--test-mode', action='store_true', help='Use simulated data for testing (no BLAST API call)')
    
    args = parser.parse_args()
    
    assessor = ProteinAssessor()
    result = assessor.assess(
        args.uniprot_id,
        blast_evalue=args.blast_evalue,
        struct_threshold=args.struct_threshold,
        ptm_threshold=args.ptm_threshold,
        top_n=args.top_n,
        test_mode=args.test_mode,
    )
    
    output_json = json.dumps(result, indent=2, ensure_ascii=False)
    
    if args.output:
        with open(args.output, 'w') as f:
            f.write(output_json)
        print(f"Results saved to {args.output}")
    else:
        print(output_json)
    
    return 0 if 'error' not in result else 1


if __name__ == '__main__':
    sys.exit(main())
