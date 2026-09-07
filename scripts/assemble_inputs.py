"""One-time local extraction. Refuses to overwrite any destination file."""
import csv
import hashlib
import json
import io
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
PKG = Path(__file__).resolve().parents[1]
RUNS = ROOT / 'data/derived/YRBP_branch_runs'
PRE = RUNS / 'v5_4_8z2_G_pre1_target_engineering_composite_rank'
PRO = RUNS / 'v5_4_9g_G_pro1_C_verified_anchor_feedback_sequence_complete'
manifest = []

def write(path, content):
    path.parent.mkdir(parents=True, exist_ok=True)
    content = content.encode() if isinstance(content, str) else content
    if path.exists():
        if path.read_bytes() != content:
            raise FileExistsError(path)
        return
    with path.open('xb') as f:
        f.write(content)

def read(path):
    with path.open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))

def extract(source, destination, predicate=None, fields=None):
    target = PKG / destination
    if predicate is None and fields is None:
        write(target, source.read_bytes())
    else:
        rows = [r for r in read(source) if predicate is None or predicate(r)]
        fields = fields or list(rows[0])
        target.parent.mkdir(parents=True, exist_ok=True)
        f = io.StringIO(newline='')
        w = csv.DictWriter(f, fields, extrasaction='ignore', lineterminator='\n')
        w.writeheader(); w.writerows(rows)
        write(target, f.getvalue())
    manifest.append(dict(source=source.relative_to(ROOT).as_posix(), destination=destination,
                         source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                         sha256=hashlib.sha256(target.read_bytes()).hexdigest(),
                         extraction='byte copy' if predicate is None and fields is None else 'column/row subset'))

for name in ['configs','inputs','src','tests','outputs','figures','manifests','IP_REVIEW_PACKAGE']:
    (PKG/name).mkdir(exist_ok=True)
extract(next(PRE.glob('*full_rank_table.csv')), 'inputs/frozen_full.csv')
extract(next(PRE.glob('*top50.csv')), 'inputs/frozen_top50.csv')
extract(next(PRE.glob('*rulebook.csv')), 'configs/original_rulebook.csv')
rows = read(PKG/'inputs/frozen_full.csv')
ids = {r['source_accession'] for r in rows[:50]}
numcols = ['target_context_score','target_organism_score','phage_origin_score','docking_support_score',
           'solubility_no_chaperone_score','E6_5_solubility_aggregation_chaperone_score',
           'curated_no_external_chaperone_required_strength','curated_soluble_expression_strength',
           'global_engineering_residual_score']
extract(next(PRE.glob('*full_rank_table.csv')), 'inputs/pre_numeric.csv', fields=numcols)
extract(next(PRE.glob('*full_rank_table.csv')), 'inputs/pre_labels.csv', fields=['source_accession',
        'G_pre1_target_engineering_composite_rank','G_pre1_target_engineering_composite_score'])
extract(next(PRO.glob('*registry.csv')), 'inputs/feedback/registry.csv')
extract(next(PRO.glob('*post_feedback_top50.csv')), 'inputs/feedback/expected_top50.csv')
extract(next(PRO.glob('*score_components_top50.csv')), 'inputs/feedback/expected_components.csv')
extract(next(PRO.glob('*scoring_policy.yaml')), 'configs/original_feedback_policy.yaml')
extract(ROOT/'data/B0_standard_candidates/v5_2/B0_standard_candidate_records.csv',
        'inputs/feedback/sequences.csv', lambda r:r.get('source_accession') in ids,
        ['source_accession','candidate_record_id','sequence_normalized','sequence_sha256','identity_status'])
extract(RUNS/'v5_4_8u_E6_scoring_after_policy_repair/YRBP_v5_4_8u_E6_component_score_table.csv',
        'inputs/feedback/e6.csv',lambda r:r.get('source_accession') in ids)
extract(ROOT/'data/B3_interproscan/v5_2/B3_candidate_domain_evidence.csv',
        'inputs/feedback/domains.csv',lambda r:r.get('accession') in ids)
extract(RUNS/'v5_4_8r_pose_contact_QC_from_existing_Vina_outputs/YRBP_v5_4_8r_top_pose_contact_summary.csv',
        'inputs/feedback/contacts.csv',lambda r:r.get('source_accession') in ids)

# Extract pure helpers and calculation block, replacing candidate-specific shortcuts.
origin = RUNS/'v5_4_9g_G_pro1_C_verified_anchor_feedback_sequence_complete_generate.py'
s = origin.read_text(encoding='utf-8')
helpers = s[s.index('def read_csv('):s.index('if not TOP50.exists():')]
body = s[s.index('record_ids ='):s.index('evidence_defs =')]
body = body.replace('    if acc == ANCHOR:\n        return 1.0\n','')
body = body.replace('return 1.0 if acc == ANCHOR else jaccard', 'return jaccard')
loop = s[s.index('post, comps ='):s.index('top5 = post[:5]')]
loop = loop.replace('1.0 if acc == ANCHOR else 0.5', "status_scores.get(acc, 0.5)")
loop = loop.replace('"verified_positive" if acc == ANCHOR else "not_tested"',
                    '"verified_positive" if acc in status_scores else "not_tested"')
loop = loop.replace('"registered_7_C_verified_rows" if acc == ANCHOR else',
                    '"registered_direct_evidence" if acc in status_scores else')
body = (body+loop).replace('Q9T0Z9','verified anchor')
body = body.replace('B0.as_posix()', 'B0.relative_to(root).as_posix()').replace('B1.as_posix()', 'B1.relative_to(root).as_posix()')
setup = '''    TOP50 = root/'inputs/frozen_top50.csv'
    E6 = root/'inputs/feedback/e6.csv'
    DOMAIN = root/'inputs/feedback/domains.csv'
    CONTACT = root/'inputs/feedback/contacts.csv'
    B0 = root/'inputs/feedback/sequences.csv'
    B1 = root/'inputs/feedback/no_optional_b1.csv'
    top50 = read_csv(TOP50)
    accessions = [r['source_accession'] for r in top50]
    evidence = read_csv(root/'inputs/feedback/registry.csv')
    status_scores = register_verified(evidence)
    if len(status_scores) != 1:
        raise ValueError('This frozen single-anchor policy requires exactly one verified candidate')
    ANCHOR = next(iter(status_scores))
'''
registry = '''def register_verified(evidence):
    grouped = defaultdict(list)
    for row in evidence:
        if row['evidence_layer'] != 'C_verified':
            raise ValueError('Unexpected evidence layer')
        grouped[row['accession']].append(row)
    result = {}
    for candidate, records in grouped.items():
        if not all(r['evidence_status'] in ('pass','strong_pass') for r in records):
            raise ValueError('Unsupported evidence status; author review required')
        result[candidate] = 1.0
    return result

'''
write(PKG/'src/feedback.py', '"""Deterministic feedback; extracted from archived v5.4.9g, registry-driven."""\n'
      +'import csv\nimport hashlib\nimport math\nfrom collections import Counter, defaultdict\nfrom pathlib import Path\n\n'
      +helpers+registry+'def reproduce_feedback(root):\n'+setup
      +''.join('    '+line+'\n' for line in body.splitlines())+'    return post, comps\n')
manifest.append(dict(source=origin.relative_to(ROOT).as_posix(),destination='src/feedback.py',
    source_sha256=hashlib.sha256(origin.read_bytes()).hexdigest(),
    sha256=hashlib.sha256((PKG/'src/feedback.py').read_bytes()).hexdigest(),extraction='helpers and calculation block; generic registry status; no accession-specific self-similarity shortcuts'))
write(PKG/'manifests/input_provenance.json',json.dumps(manifest,indent=2)+'\n')
print('Extracted',len(manifest),'sources; all destinations created exclusively.')
