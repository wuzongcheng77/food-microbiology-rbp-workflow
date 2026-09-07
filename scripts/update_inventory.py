"""Refresh review-only inventories of this task's package, never original sources."""
import csv
import hashlib
import json
import re
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
REPO=ROOT.parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def writecsv(p,rows,fields):
    with p.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fields,lineterminator='\n');w.writeheader();w.writerows(rows)

sources=['README.md','ROOT_LAYOUT.md','PIPELINE_RUNBOOK.md','requirements.txt',
 'tools/yrbp_v5_4_8u_e6_scoring_after_policy_repair.py',
 'tools/yrbp_v5_4_8v_e7_unlocked_rank.py',
 'tools/yrbp_v5_4_8z2_g_pre1_target_engineering_composite_rank.py',
 'data/derived/YRBP_branch_runs/v5_4_9g_G_pro1_C_verified_anchor_feedback_sequence_complete_generate.py']
writecsv(ROOT/'manifests/original_code_and_constraints.csv',
 [dict(source=s,sha256=sha(REPO/s)) for s in sources],['source','sha256'])
for item in json.loads((ROOT/'manifests/input_provenance.json').read_text()):
    assert sha(REPO/item['source'])==item['source_sha256'],item['source']
for rel,h in json.loads((ROOT/'manifests/expected_output_hashes.json').read_text()).items():
    assert sha(ROOT/rel)==h,rel
withheld=ROOT/'IP_REVIEW_PACKAGE/WITHHELD_FILE_LIST.csv'
with withheld.open(encoding='utf-8') as f:rows=list(csv.DictReader(f))
for r in rows:
    if r['path']=='outputs/full827/':r['reason']='Initial task-generated development run removed; portable reference outputs retained'
writecsv(withheld,rows,['path','reason'])
patterns=[r'-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----',r'\bAKIA[0-9A-Z]{16}\b',
          r'\bgh[pousr]_[A-Za-z0-9]{30,}\b',
          r'(?i)(?:api_key|password|client_secret|access_token)\s*[:=]\s*[\"\'][^\"\']{12,}[\"\']']
findings=[];scanned=0
for p in ROOT.rglob('*'):
    if not p.is_file() or p.suffix.lower() in {'.png','.tiff','.pyc'}:continue
    try:text=p.read_text(encoding='utf-8-sig')
    except UnicodeError:continue
    scanned+=1
    if any(re.search(pattern,text) for pattern in patterns):findings.append(p.relative_to(ROOT).as_posix())
(ROOT/'IP_REVIEW_PACKAGE/SECRETS_SCAN_REPORT.md').write_text(
    '# Final local secrets scan\n\n'
    +f'Scanned {scanned} current text files including input CSVs. Binary images excluded. '
    +f'Pattern matches: {len(findings)}.\n\n'+json.dumps(findings,indent=2)+'\n\n'
    +'Patterns cover private-key headers, AWS access keys, GitHub tokens and credential-like assignments. '
    +'No credential values are disclosed. This is a bounded pattern scan, not proof of absence. '
    +'No accounts, connectors or deployment modules were extracted.\n',encoding='utf-8',newline='\n')
hits=[]
for p in sorted(ROOT.rglob('*')):
    if p.suffix not in {'.py','.md','.json','.yaml','.yml','.cff','.txt','.ps1'}:continue
    for number,line in enumerate(p.read_text(encoding='utf-8-sig').splitlines(),1):
        if 'Q9T0Z9' in line or 'Gp17' in line:
            hits.append(dict(path=p.relative_to(ROOT).as_posix(),line=number,
                classification='test/assertion' if 'tests' in p.parts else 'report/provenance/configuration'))
writecsv(ROOT/'manifests/identifier_occurrences.csv',hits,['path','line','classification'])
index=ROOT/'IP_REVIEW_PACKAGE/PUBLIC_FILE_LIST.csv'
rows=[]
for p in sorted(ROOT.rglob('*')):
    if not p.is_file() or '__pycache__' in p.parts:continue
    rows.append(dict(path=p.relative_to(ROOT).as_posix(),sha256='SELF_HASH_EXCLUDED' if p==index else sha(p),
       bytes='' if p==index else p.stat().st_size,status='PROPOSED_ONLY_PENDING_AUTHOR_IP_AND_LICENSE_REVIEW'))
writecsv(index,rows,['path','sha256','bytes','status'])
print('Original inputs intact; expected outputs verified; inventory:',len(rows),'files')
