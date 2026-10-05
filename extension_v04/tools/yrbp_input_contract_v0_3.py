"""RBP input provenance audit and raw sequence descriptors; no binding predictor.

The receptor gate checks declared evidence, not biological truth. An accepted
contract never establishes receptor identity or experimentally validates docking.
"""
import hashlib
import math
from collections import Counter

VERSION = '0.3.0'
AA = 'ACDEFGHIKLMNPQRSTVWY'

def sequence_features(sequence):
    if not sequence or set(sequence) - set(AA):
        raise ValueError('Expected a nonempty canonical AA20 sequence')
    counts = Counter(sequence)
    n = len(sequence)
    return {
        'sequence_sha256': hashlib.sha256(sequence.encode('ascii')).hexdigest(),
        'feature_schema': 'raw_composition_v1',
        'feature_execution': 'computed_from_input_sequence',
        'sequence_length': n,
        **{'count_' + a: counts[a] for a in AA},
        **{'fraction_' + a: counts[a] / n for a in AA},
        'composition_entropy_bits': -sum((c/n)*math.log2(c/n) for c in counts.values()),
        'interpretation': 'descriptive_only_not_binding_or_solubility_prediction',
    }

def receptor_contract(sequence_sha256, spec):
    """Fail closed for missing/mismatched declarations before any docking call."""
    c = (spec or {}).get('receptor_evidence', {})
    reasons = []
    required = ('sequence_sha256', 'target_taxon', 'target_context',
                'receptor_identity', 'receptor_class', 'evidence_scope',
                'evidence_source', 'reviewer', 'reviewed_on',
                'ligand_identity', 'ligand_class', 'method',
                'grid_rationale_source', 'structure_mapping_source')
    for key in required:
        if not c.get(key) or str(c[key]).strip().lower() in ('unknown', 'pending', 'not_done'):
            reasons.append('missing:' + key)
    if c.get('sequence_sha256') and c['sequence_sha256'] != sequence_sha256:
        reasons.append('sequence_identity_mismatch')
    if c.get('receptor_class') and c.get('ligand_class') and c['receptor_class'] != c['ligand_class']:
        reasons.append('receptor_ligand_class_mismatch')
    if c.get('method') != 'vina_prepared_ligand':
        reasons.append('unsupported_method_in_this_entry')
    if c.get('receptor_class') not in ('glycan', 'small_molecule'):
        reasons.append('receptor_class_outside_this_adapter')
    if c.get('evidence_scope') != 'candidate_specific':
        reasons.append('candidate_specific_evidence_not_established')
    if c.get('review_status') != 'reviewed':
        reasons.append('evidence_review_incomplete')
    return {'status': 'declaration_consistent' if not reasons else 'requires_review',
            'reasons': reasons,
            'scope': 'declaration_check_only_not_biological_validation'}
