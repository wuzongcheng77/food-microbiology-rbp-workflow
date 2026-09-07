"""Build local audit metadata; does not grant publication rights or change results."""
import csv
import hashlib
import importlib.metadata as md
import json
from pathlib import Path
import platform
import re
import sys

ROOT=Path(__file__).resolve().parents[1]
REPO=ROOT.parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(rel,text):
    p=ROOT/rel;p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('x',encoding='utf-8',newline='\n') as f:f.write(text)
def csvout(rel,rows,fields):
    p=ROOT/rel
    with p.open('x',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fields,lineterminator='\n',extrasaction='ignore');w.writeheader();w.writerows(rows)
def read(p):
    with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))

deps=['numpy','scipy','matplotlib','pillow','pytest','contourpy','cycler','fonttools','kiwisolver',
      'packaging','pyparsing','python-dateutil','six','iniconfig','pluggy','colorama']
write('requirements.lock.txt','# Exact versions observed locally; Python 3.12.3. No installation was performed.\n'
      +'\n'.join(f'{p}=={md.version(p)}' for p in deps)+'\n')
write('environment.yml','name: food-microbiology-rbp-v1\nchannels:\n  - conda-forge\ndependencies:\n  - python=3.12.3\n  - pip\n  - pip:\n    - -r requirements.lock.txt\n')
write('manifests/runtime.json',json.dumps(dict(python=sys.version,platform=platform.platform(),
      dependencies={p:md.version(p) for p in deps},portability='Tested on this Windows environment only; no clean-environment install performed'),indent=2)+'\n')
# Refresh only the manifest this task created, after final path normalization.
p=ROOT/'manifests/input_provenance.json';items=json.loads(p.read_text())
for item in items:item['sha256']=sha(ROOT/item['destination'])
p.write_text(json.dumps(items,indent=2)+'\n',encoding='utf-8',newline='\n')
expected={p.relative_to(ROOT).as_posix():sha(p) for folder in ['outputs/reference_full827','outputs/reference_top50','figures/Figure6']
          for p in sorted((ROOT/folder).iterdir()) if p.is_file()}
write('manifests/expected_output_hashes.json',json.dumps(expected,indent=2)+'\n')

compare=[dict(metric='G_pre1 candidate count',scope='full827',manuscript_or_frozen='827',reproduced='827',status='MATCH'),
         dict(metric='Q9T0Z9 G_pre1 rank',scope='full827',manuscript_or_frozen='1',reproduced='1',status='MATCH'),
         dict(metric='Q9T0Z9 G_pre1 score',scope='full827',manuscript_or_frozen='0.856763',reproduced='0.856763',status='MATCH'),
         dict(metric='G_pro1 scores and ranks',scope='top50',manuscript_or_frozen='archived v5.4.9g',reproduced='50/50',status='MATCH')]
for scope in ['full827','top50']:
    for r in read(ROOT/f'outputs/reference_{scope}/Supplementary_Table_S4.csv'):
        if r['analysis']=='full':continue
        for metric in ['target_rank','top5_overlap','top10_overlap','spearman_rho','kendall_tau']:
            compare.append(dict(metric=r['analysis']+'.'+metric,scope=scope,
                                manuscript_or_frozen='NOT_PROVIDED',reproduced=r[metric],status='AUTHOR_COMPARISON_PENDING'))
    s=json.loads((ROOT/f'outputs/reference_{scope}/summary.json').read_text())
    for metric in ['target_top1_frequency','target_top5_frequency','score_gap_mean','score_gap_quantiles_025_975']:
        compare.append(dict(metric=metric,scope=scope,manuscript_or_frozen='NOT_PROVIDED',
                            reproduced=str(s[metric]),status='AUTHOR_COMPARISON_PENDING'))
csvout('manifests/manuscript_B_comparison.csv',compare,['metric','scope','manuscript_or_frozen','reproduced','status'])

# Audit source occurrences without emitting credential material.
hits=[]
for p in sorted(ROOT.rglob('*')):
    if p.suffix not in {'.py','.md','.json','.yaml','.yml','.cff','.txt','.ps1'}:continue
    for number,line in enumerate(p.read_text(encoding='utf-8-sig').splitlines(),1):
        if 'Q9T0Z9' in line or 'Gp17' in line:
            hits.append(dict(path=p.relative_to(ROOT).as_posix(),line=number,
                classification='test/assertion' if 'tests' in p.parts else 'report/provenance/configuration'))
csvout('manifests/identifier_occurrences.csv',hits,['path','line','classification'])
secret_hits=[]
patterns=[('private_key',r'-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----'),
          ('aws_access_key',r'\bAKIA[0-9A-Z]{16}\b'),
          ('github_token',r'\bgh[pousr]_[A-Za-z0-9]{30,}\b'),
          ('credential_assignment',r'(?i)(?:api_key|password|client_secret|access_token)\s*[:=]\s*[\"\'][^\"\']{12,}[\"\']')]
scanned=0
for p in ROOT.rglob('*'):
    if not p.is_file() or p.suffix.lower() in {'.png','.tiff','.pyc'}:continue
    try:text=p.read_text(encoding='utf-8-sig')
    except UnicodeError:continue
    scanned+=1
    for kind,pattern in patterns:
        if re.search(pattern,text):secret_hits.append(dict(path=p.relative_to(ROOT).as_posix(),kind=kind))
write('IP_REVIEW_PACKAGE/SECRETS_SCAN_REPORT.md',f'''# Local secrets scan

Scanned {scanned} text files in this new package, including input CSVs. Binary images excluded.
Matched {len(secret_hits)} files/patterns. Values are never included in this report.

{json.dumps(secret_hits,indent=2)}

Patterns: private-key headers, AWS access-key identifiers, GitHub tokens and credential-like assignments.
This bounded pattern scan is not proof of absence. No credentials, accounts, network connectors or deployment modules were intentionally extracted.
''')
public=[]
for p in sorted(ROOT.rglob('*')):
    if not p.is_file() or '__pycache__' in p.parts:continue
    rel=p.relative_to(ROOT).as_posix()
    if rel.startswith('outputs/full827/') or rel.startswith('outputs/replay_'):continue
    public.append(dict(path=rel,sha256=sha(p),bytes=p.stat().st_size,status='PROPOSED_ONLY_PENDING_AUTHOR_IP_AND_LICENSE_REVIEW'))
csvout('IP_REVIEW_PACKAGE/PUBLIC_FILE_LIST.csv',public,['path','sha256','bytes','status'])
withheld=[dict(path='Repository modules outside this package',reason='Unrelated YAMPs/LAMPs/LYSIN/UI/accounts/connectors/deployment; excluded by scope'),
    dict(path='Full B0/B3 source tables',reason='Only Top50 sequence/domain subsets extracted; unrelated records withheld'),
    dict(path='Raw wet-lab measurements and manuscript source',reason='Not located/provided for this task; not invented or reconstructed'),
    dict(path='outputs/full827/',reason='Initial local development run; superseded by portable reference_full827; not proposed for public release'),
    dict(path='IP_REVIEW_PACKAGE/PUBLIC_FILE_LIST.csv',reason='Index itself has no recursive self-hash; review metadata only'),
    dict(path='Git metadata, environments, caches, secrets and upstream structures',reason='Excluded; only minimal frozen downstream inputs packaged')]
csvout('IP_REVIEW_PACKAGE/WITHHELD_FILE_LIST.csv',withheld,['path','reason'])
print('Audit manifests and dependency pins created.')
