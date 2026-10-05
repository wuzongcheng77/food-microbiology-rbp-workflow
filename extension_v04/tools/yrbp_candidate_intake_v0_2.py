"""Versioned protein intake and evidence audit in the original YRBP tool tree.

This is a partial execution entry, not a genome-to-RBP predictor. An absent
structure requests evidence; externally supplied feature values are not new
predictions. The frozen historical composite remains a separate version.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
from yrbp_candidate_intake_support_v0_2 import run
from yrbp_missing_evidence_v0_2 import summarize, rank_bounds

VERSION = '0.2.0'

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def audit(records, policy, evidence):
    if policy.get('version') is None or not policy.get('normalization_reference'):
        raise ValueError('Policy needs version and a declared normalization reference')
    weights = policy.get('weights', {})
    if not weights or any(isinstance(v,bool) or not isinstance(v,(int,float)) or not math.isfinite(v) or v<=0 for v in weights.values()):
        raise ValueError('Declare finite positive feature weights before comparing results')
    ids={r['candidate_id'] for r in records}
    if set(evidence)-ids:raise ValueError('Evidence contains IDs absent from input FASTA')
    result=[]
    for row in records:
        item=evidence.get(row['candidate_id'],{})
        if item and item.get('sequence_sha256')!=row['sequence_sha256']:
            raise ValueError('Evidence sequence hash mismatch: '+row['candidate_id'])
        features=item.get('features',{})
        if set(features)-set(weights):raise ValueError('Undeclared feature')
        values={}; provenance={}
        for key in weights:
            spec=features.get(key,{'status':'unknown'})
            status=spec.get('status','unknown')
            if status not in {'unknown','observed','computed','curated'}:raise ValueError('Unsupported evidence status')
            value=spec.get('value')
            if status=='unknown':
                if value is not None:raise ValueError('Unknown evidence cannot have a score')
                values[key]=None
            else:
                if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value) or not 0<=value<=1:
                    raise ValueError('Known feature requires a finite normalized value in [0,1]')
                if not spec.get('source'):raise ValueError('Known feature needs source provenance')
                values[key]=value
            provenance[key]=spec
        summary=summarize(values,weights)
        result.append({'candidate_id':row['candidate_id'],'sequence_sha256':row['sequence_sha256'],**summary,
                       'evidence':provenance,'literature':item.get('literature',{'status':'not_searched'}),
                       'wet_lab_status':item.get('wet_lab_status','not_tested')})
    return rank_bounds(result)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--fasta',required=True);p.add_argument('--manifest');p.add_argument('--vina')
    p.add_argument('--out',required=True);p.add_argument('--policy');p.add_argument('--evidence')
    p.add_argument('--resume',action='store_true',help='Continue an existing run dir; reuse fingerprint-matched completed poses as cache_reuse.')
    args=p.parse_args()
    if args.evidence and not args.policy:p.error('--evidence requires --policy')
    # Validate feature contracts before running expensive docking or creating output.
    from yrbp_candidate_intake_support_v0_2 import fasta
    records=[{'candidate_id':n,'sequence_sha256':hashlib.sha256(s.encode()).hexdigest()} for n,s in fasta(args.fasta)]
    policy=json.loads(Path(args.policy).read_text(encoding='utf-8-sig')) if args.policy else None
    evidence=json.loads(Path(args.evidence).read_text(encoding='utf-8-sig')) if args.evidence else {}
    audited=audit(records,policy,evidence) if policy else None
    output=Path(args.out).resolve()
    if output.exists() and not args.resume:p.error('Output must be a new run directory; existing results are never overwritten (pass --resume to continue a prior run)')
    rows=run(args.fasta,args.manifest,output,args.vina,resume=args.resume)
    if audited is not None:
        (output/'evidence_audit.json').write_text(json.dumps(audited,ensure_ascii=False,indent=2),encoding='utf-8')
    sources={'fasta':{'path':str(Path(args.fasta).resolve()),'sha256':sha(args.fasta)}}
    for name in ('manifest','policy','evidence'):
        value=getattr(args,name)
        if value:sources[name]={'path':str(Path(value).resolve()),'sha256':sha(value)}
    metadata={'version':VERSION,'sources':sources,'tool_sha256':sha(__file__),
              'candidate_count':len(rows),'evidence_audit_available':audited is not None,
              'stage_boundary':'Protein FASTA intake; optional prepared-structure docking; optional supplied-feature uncertainty audit.',
              'not_implemented':['genome ORF calling','automatic RBP identification','automatic folding in this entry','automatic full E6/E7 feature extraction'],
              'scope':'Computational execution and evidence bookkeeping only; no biological accuracy claim; rank bounds conditional on supplied normalized features and fixed policy.'}
    (output/'workflow_manifest.json').write_text(json.dumps(metadata,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps({'out':str(output),'candidates':len(rows),'version':VERSION}))

if __name__=='__main__':main()
