"""Read-only sequence identity checks for the frozen YRBP candidate set."""
import csv,hashlib
from pathlib import Path
REPO=Path(__file__).resolve().parents[1]
RANK=REPO/'data/derived/YRBP_branch_runs/v5_4_8z2_G_pre1_target_engineering_composite_rank/YRBP_v5_4_8z2_G_pre1_target_engineering_composite_full_rank_table.csv'
B0=REPO/'data/B0_standard_candidates/v5_2/B0_standard_candidate_records.csv'
AA=set('ACDEFGHIKLMNPQRSTVWY')
def load_registry(rank_path=RANK,b0_path=B0):
    with Path(rank_path).open(encoding='utf-8-sig') as f:ranks=list(csv.DictReader(f))
    with Path(b0_path).open(encoding='utf-8-sig') as f:base=list(csv.DictReader(f))
    by_id={}
    for row in base:
        key=row['candidate_record_id']
        if key in by_id and by_id[key]['sequence_normalized']!=row['sequence_normalized']:raise ValueError('Conflicting B0 record ID')
        by_id[key]=row
    result=[]
    for row in ranks:
        b=by_id.get(row['candidate_record_id'])
        if b is None:raise ValueError('Missing frozen sequence '+row['source_accession'])
        seq=''.join(b['sequence_normalized'].split()).upper()
        digest=hashlib.sha256(seq.encode()).hexdigest()
        if digest!=b['sequence_sha256']:raise ValueError('Stored sequence hash mismatch '+row['source_accession'])
        result.append({**row,'sequence':seq,'sequence_sha256':digest,'invalid_residues':''.join(sorted(set(seq)-AA)),'sequence_source':str(b0_path)})
    return result
def membership(accession,digest,registry):
    named=[r for r in registry if r['source_accession']==accession]
    matches=[r for r in registry if r['sequence_sha256']==digest]
    status=('accession_and_sequence_match' if named and any(r['sequence_sha256']==digest for r in named)
            else 'accession_sequence_mismatch' if named
            else 'sequence_match_under_other_id' if matches else 'outside_frozen_exact_sequence_set')
    return {'accession_in_historical_827':bool(named),'sequence_in_historical_827':bool(matches),
            'historical_identity_status':status,'historical_sequence_accessions':';'.join(r['source_accession'] for r in matches),
            'external_independence':'not_established_requires_homology_and_development_audit' if not matches else 'not_sequence_independent'}
