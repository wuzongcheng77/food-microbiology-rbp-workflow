"""Offline B-stage analysis. Existing output directories are never overwritten."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
from scipy.stats import spearmanr, kendalltau

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from scoring import load_numeric, components, scores, rounded, rank, frozen_ties, analysis_ties, WEIGHTS
from feedback import reproduce_feedback

def read(path):
    with path.open(encoding='utf-8-sig',newline='') as f:
        return list(csv.DictReader(f))

def csvout(path, rows, fields=None):
    fields = fields or list(rows[0])
    with path.open('x',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fields,lineterminator='\n',extrasaction='ignore')
        w.writeheader(); w.writerows(rows)

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def run(out,scope):
    out.mkdir(parents=True,exist_ok=False)
    cfg=json.loads((ROOT/'configs/analysis.json').read_text())
    n=827 if scope=='full827' else 50
    a=load_numeric(ROOT/'inputs/pre_numeric.csv')[:n]
    labels=read(ROOT/'inputs/pre_labels.csv')[:n]
    ids=[r['source_accession'] for r in labels]
    target=ids.index(cfg['report_target'])
    baseline=np.array([int(r['G_pre1_target_engineering_composite_rank']) for r in labels])
    expected=np.array([float(r['G_pre1_target_engineering_composite_score']) for r in labels])
    x=components(a)
    lex={v:i for i,v in enumerate(sorted(ids))}
    ties=frozen_ties(a,[lex[v] for v in ids])
    anon=analysis_ties(n,cfg['seed'])
    value=rounded(scores(x,WEIGHTS)); order,ranks=rank(value,ties)
    differences=[]
    rounded_only=rounded(scores(components(a,True),WEIGHTS))
    for i in range(n):
        if rounded_only[i]!=expected[i]:
            differences.append(dict(accession=ids[i],frozen=f'{expected[i]:.6f}',
                rounded_components=f'{rounded_only[i]:.6f}',reconstructed=f'{value[i]:.6f}',
                cause='intermediate_solubility_rounding',resolved=value[i]==expected[i]))
    csvout(out/'precision_differences.csv',differences,['accession','frozen','rounded_components','reconstructed','cause','resolved'])
    def ranking_file(name,v,o,r):
        csvout(out/(name+'.csv'),[dict(accession=ids[i],frozen_rank=int(baseline[i]),rank=int(r[i]),
                    score=f'{v[i]:.12f}') for i in o])
    ranking_file('baseline',value,order,ranks)
    if not np.array_equal(value,expected) or not np.array_equal(ranks,baseline):
        csvout(out/'BLOCKING_baseline_differences.csv',[dict(accession=ids[i],expected_score=expected[i],
            actual_score=value[i],expected_rank=baseline[i],actual_rank=ranks[i]) for i in range(n)
            if value[i]!=expected[i] or ranks[i]!=baseline[i]])
        raise ValueError('Frozen baseline discrepancy: author confirmation required')
    summaries=[]
    def summarize(name,v,t):
        o,r=rank(v,t); ranking_file(name,v,o,r)
        top=lambda k:set(np.flatnonzero(baseline<=k))
        row=dict(analysis=name,scope=scope,candidates=n,target_rank=int(r[target]),
            top5_overlap=len(set(o[:5])&top(5)),top10_overlap=len(set(o[:10])&top(10)),
            spearman_rho=f'{spearmanr(baseline,r).statistic:.12f}',
            kendall_tau=f'{kendalltau(baseline,r).statistic:.12f}')
        summaries.append(row)
        return row
    summarize('full',value,ties)
    names=cfg['components']
    # Sensitivity analyses use score then seeded anonymous ID, avoiding reintroduction
    # of omitted evidence through secondary tie-breaks. Scores remain unrounded.
    for j,name in enumerate(names):
        w=WEIGHTS.copy(); w[j]=0
        summarize('without_'+name,scores(x,w),anon)
        w=np.zeros(4); w[j]=1
        summarize('only_'+name,scores(x,w),anon)
    # Masking preserves frozen numerical columns; no reconstruction or metadata read.
    masked=components(a,True)
    summarize('identifier_masked',scores(masked,WEIGHTS),anon)
    w=WEIGHTS.copy(); w[0]=0
    summarize('strict_context_blind',scores(masked,w),anon)
    csvout(out/'Supplementary_Table_S4.csv',summaries)
    csvout(out/'Figure6_source_summary.csv',summaries)
    rng=np.random.default_rng(cfg['seed'])
    perturb=WEIGHTS*rng.uniform(.8,1.2,size=(cfg['iterations'],4))
    perturb/=perturb.sum(axis=1,keepdims=True)
    counts={k:np.zeros(n,dtype=int) for k in [1,5,10,50]}
    trials=[]
    for it,w in enumerate(perturb,1):
        v=scores(x,w); o,r=rank(v,anon)
        for k in counts: counts[k][o[:min(k,n)]]+=1
        rival=next(i for i in o if i!=target)
        trials.append(dict(iteration=it,**{name+'_weight':f'{w[j]:.17g}' for j,name in enumerate(names)},
                target_rank=int(r[target]),best_competitor=ids[rival],
                score_gap=f'{v[target]-v[rival]:.17g}'))
    csvout(out/'weight_perturbation_10000.csv',trials)
    freq=[dict(accession=ids[i],frozen_rank=int(baseline[i]),
        **{f'top{k}_frequency':f'{counts[k][i]/cfg["iterations"]:.6f}' for k in counts}) for i in range(n)]
    csvout(out/'candidate_inclusion_frequency.csv',freq)
    csvout(out/'Figure6_source_frequency.csv',freq)
    if n==827: csvout(out/'original_rank51_827_inclusion_frequency.csv',freq[50:])
    gaps=np.array([float(t['score_gap']) for t in trials]); tr=[t['target_rank'] for t in trials]
    summary=dict(scope=scope,candidates=n,seed=cfg['seed'],iterations=cfg['iterations'],
        baseline_score_matches=int(np.sum(value==expected)),baseline_rank_matches=int(np.sum(ranks==baseline)),
        target=cfg['report_target'],target_score=f'{value[target]:.6f}',target_rank=int(ranks[target]),
        target_top1_frequency=counts[1][target]/cfg['iterations'],target_top5_frequency=counts[5][target]/cfg['iterations'],
        target_rank_min=min(tr),target_rank_max=max(tr),
        target_rank_quantiles_025_975=np.quantile(tr,[.025,.975]).tolist(),
        score_gap_mean=float(gaps.mean()),score_gap_quantiles_025_975=np.quantile(gaps,[.025,.975]).tolist(),
        interval_interpretation='empirical central 95% perturbation interval; not a sampling confidence interval')
    post,comps=reproduce_feedback(ROOT)
    csvout(out/'G_pro1_top50.csv',post)
    csvout(out/'G_pro1_components.csv',comps)
    old={r['accession']:r for r in read(ROOT/'inputs/feedback/expected_top50.csv')}
    bad=[]
    for r in post:
        for field in ['G_pro1_rank','G_pro1_score','C_verified_score','C_verified_anchor_similarity_score','wet_lab_status']:
            if r[field]!=old[r['accession']][field]:
                bad.append(dict(accession=r['accession'],field=field,expected=old[r['accession']][field],actual=r[field],
                                cause='generic_registry_or_removed_self_similarity_shortcut_requires_review'))
    csvout(out/'G_pro1_differences.csv',bad,['accession','field','expected','actual','cause'])
    summary['G_pro1_matches']=50-len({r['accession'] for r in bad})
    with (out/'summary.json').open('x',encoding='utf-8',newline='\n') as f:
        f.write(json.dumps(summary,indent=2)+'\n')
    with (out/'output_hashes.json').open('x',encoding='utf-8',newline='\n') as f:
        f.write(json.dumps({p.name:digest(p) for p in sorted(out.iterdir()) if p.is_file() and p.name!='output_hashes.json'},indent=2)+'\n')
    print(json.dumps(summary,indent=2))
    if bad: raise ValueError('G_pro1 discrepancy: author confirmation required')

if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('--out',type=Path,required=True)
    p.add_argument('--scope',choices=['full827','top50'],default='full827'); args=p.parse_args()
    run(args.out,args.scope)
