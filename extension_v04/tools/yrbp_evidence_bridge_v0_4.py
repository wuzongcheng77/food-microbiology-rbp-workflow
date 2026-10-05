"""Sequence-derived E6 annotations and structure correspondence, without score imputation.

This module audits inputs; it does not identify a biological receptor or predict
binding. It never converts a residue count into an accessible surface handle.
"""
import hashlib
from pathlib import Path
from yrbp_input_contract_v0_3 import sequence_features

VERSION = '0.4.0'
AA3 = dict(zip(
    'ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR'.split(),
    'ACDEFGHIKLMNPQRSTVWY'))

def e6_annotations(candidate_id, sequence):
    f = sequence_features(sequence)
    return {
        'candidate_id': candidate_id, 'accession': candidate_id,
        'sequence_sha256': f['sequence_sha256'],
        'sequence_length_aa': len(sequence), 'cysteine_count': sequence.count('C'),
        'lysine_count': sequence.count('K'),
        'expression_annotation_source': 'computed_from_input_sequence_sha256',
        'expression_notes': 'Length and C count only; solubility/TM/expression not predicted',
        'immobilization_notes': 'K count only; surface accessibility and handles not established',
        'predicted_solubility_label': '', 'predicted_tm_helix_count': '',
        'signal_peptide_predicted': '', 'known_chaperone_required': '',
        'surface_lysine_proxy': '', 'cysteine_handle_count': '',
        'expression_feature_status': 'partial_raw_sequence_only',
        'chaperone_feature_status': 'not_annotated',
        'immobilization_feature_status': 'partial_raw_sequence_only',
        'engineering_scoring_ready': 'false',
        'normalized_score_status': 'not_generated',
    }

def structure_sequence_audit(sequence, spec, base):
    """Check ATOM residues in model 1; no structural quality/assembly inference.

    Only exact complete canonical chains pass automatically. Fragments, modified
    residues and unexpected chains require explicit review, never silent repair.
    """
    if not spec or not spec.get('protein_pdbqt'):
        return {'status': 'structure_required', 'chains': [], 'automatic_acceptance': False}
    path = (Path(base) / spec['protein_pdbqt']).resolve()
    if not path.is_file():
        return {'status': 'structure_file_missing', 'chains': [], 'automatic_acceptance': False}
    raw = path.read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    if digest != spec.get('protein_sha256'):
        return {'status': 'structure_hash_mismatch', 'chains': [], 'automatic_acceptance': False}
    chains, residues, conflicts = {}, {}, []
    for line in raw.decode('ascii', errors='replace').splitlines():
        if line.startswith('ENDMDL'): break
        if not line.startswith('ATOM  '): continue
        chain = line[21:22].strip() or '_'
        key = (chain, line[22:27])
        name = line[17:20].strip()
        if key in residues:
            if residues[key] != name: conflicts.append(list(key))
            continue
        residues[key] = name
        chains.setdefault(chain, []).append(AA3.get(name, '?'))
    rows = []
    for chain, values in chains.items():
        observed = ''.join(values)
        if '?' in observed: status = 'noncanonical_residue_requires_review'
        elif observed == sequence: status = 'exact_complete_chain'
        elif observed in sequence: status = 'contiguous_fragment_requires_mapping'
        else: status = 'sequence_mismatch_or_discontinuous_mapping'
        rows.append({'chain': chain, 'observed_residues': len(observed),
                     'input_residues': len(sequence), 'status': status,
                     'observed_sequence_sha256': hashlib.sha256(observed.encode()).hexdigest()})
    accepted = bool(rows) and not conflicts and all(r['status'] == 'exact_complete_chain' for r in rows)
    return {'status': 'exact_complete_chains' if accepted else 'structure_mapping_requires_review',
            'chains': rows, 'residue_conflicts': conflicts, 'protein_sha256': digest,
            'automatic_acceptance': accepted,
            'scope': 'Residue correspondence only; no fold, assembly, ligand or grid validation'}

def evidence_request(candidate_id, gate, structure):
    reasons = gate['reasons'][:]
    if not structure['automatic_acceptance']: reasons.append(structure['status'])
    return {'candidate_id': candidate_id,
            'decision': 'declarations_and_sequence_mapping_consistent' if not reasons else 'abstain_from_docking',
            'outstanding_reasons': ';'.join(reasons),
            'engineering_action': 'retain_raw_annotations; obtain_missing_predictor_evidence',
            'scope': 'Data-readiness decision; not a prediction of binding or failure'}
