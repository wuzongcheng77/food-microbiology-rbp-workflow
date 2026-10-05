"""Portable sequence/audit reproduction. No D/E drive dependencies or docking."""
from pathlib import Path
import sys, subprocess, csv, json, hashlib
R=Path(__file__).resolve().parent
def call(args): subprocess.run([sys.executable,'-B','-X','utf8']+list(map(str,args)),check=True)
for name,fa,extra in [('full827','historical827.fasta',[]),('external','external.fasta',['--receptor-evidence',R/'data/external_receptor_evidence.json'])]:
    out=R/'reproduced'/name
    args=[R/'tools/yrbp_v0_2_rbp_repro_runner.py','--fasta',R/'data'/fa,'--out',out,'--rank-registry',R/'data/frozen_rank.csv','--sequence-registry',R/'data/sequence_registry.csv']+extra
    if out.exists(): args+=['--resume']
    call(args)
call([R/'analysis/audit_full827.py'])
with (R/'reproduced/full827_audit/ablation_summary.csv').open(encoding='utf-8-sig') as f: actual=list(csv.DictReader(f))
with (R/'reference/ablation_summary.csv').open(encoding='utf-8-sig') as f: expected=list(csv.DictReader(f))
assert actual==expected, 'Ablation summary differs'
manifest=json.loads((R/'reproduced/full827/fresh_feature_manifest.json').read_text())
assert (manifest['records'],manifest['unique_sequences'])==(827,744)
with (R/'reproduced/full827/fresh_e6_annotations.csv').open() as f: annotations=list(csv.DictReader(f))
assert all(r['engineering_scoring_ready']=='false' and not r['predicted_solubility_label'] for r in annotations)
gate=json.loads((R/'reproduced/external/receptor_contract_audit.json').read_text())['I7LEH0']
assert 'receptor_class_outside_this_adapter' in gate['reasons']
assert 'candidate_specific_evidence_not_established' in gate['reasons']
source_before=hashlib.sha256((R/'data/frozen_rank.csv').read_bytes()).hexdigest()
assert source_before=='85c3c3b9b8e5db3827fbdc7945426666226beb46b8461dd8605e187e50c94886'
if '--with-scientific-dependencies' in sys.argv:
    for name in ['full827','external']:
        out=R/'reproduced'/name; annotation=out/'fresh_e6_annotations.csv'
        call([R/'tools/yrbp_e6_annotate_engineering_features_v0_1.py','--e6_input_csv',annotation,'--annotation_csv',annotation,'--out_csv',out/'legacy_E6_consumed_annotations.csv'])
    call([R/'analysis/analyze_elisa.py'])
result={'status':'PASS','records':827,'unique_sequences':744,'ablation_scenarios':5,'external_case':'intake_and_abstention','new_candidate_predictive_validation':False}
(R/'reproduced/verification.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result))
