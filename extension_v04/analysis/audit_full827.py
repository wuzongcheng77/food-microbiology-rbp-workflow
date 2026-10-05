import sys,csv,json,hashlib,collections
from pathlib import Path
root=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(root/'tools'))
from yrbp_sequence_registry_v0_2 import load_registry
rank=root/'data/frozen_rank.csv'
b0=root/'data/sequence_registry.csv'
out=root/'reproduced/full827_audit';out.mkdir(parents=True,exist_ok=True)
rows=load_registry(rank,b0)
def csvwrite(name,items):
    with (out/name).open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(items[0]));w.writeheader();w.writerows(items)
groups=collections.defaultdict(list)
for r in rows:groups[r['sequence_sha256']].append(r['source_accession'])
csvwrite('sequence_identity_827.csv',[{'accession':r['source_accession'],'candidate_record_id':r['candidate_record_id'],'sequence_sha256':r['sequence_sha256'],'length':len(r['sequence']),'invalid_residues':r['invalid_residues'],'exact_sequence_group_size':len(groups[r['sequence_sha256']]),'exact_sequence_group':';'.join(groups[r['sequence_sha256']])} for r in rows])
(out/'frozen827_sequences.fasta').write_text(''.join('>'+r['source_accession']+'\n'+r['sequence']+'\n' for r in rows),encoding='utf-8')
fields={'target_context_score':.55,'docking_support_score':.20,'solubility_no_chaperone_score':.20,'global_engineering_residual_score':.05}
baseline={r['source_accession']:sum(float(r[k])*w for k,w in fields.items()) for r in rows}
err=max(abs(baseline[r['source_accession']]-float(r['G_pre1_target_engineering_composite_score'])) for r in rows)
summaries=[];details=[]
for removed in [None,*fields]:
    intervals=[]
    for r in rows:
        lo=sum(float(r[k])*w for k,w in fields.items() if k!=removed)
        missing=fields.get(removed,0)
        intervals.append({'accession':r['source_accession'],'scenario':removed or 'complete','coverage':1-missing,'lower':lo,'upper':lo+missing,'historical_rank':int(r['G_pre1_target_engineering_composite_rank'])})
    for x in intervals:
        x['best_possible_rank']=1+sum(y['lower']>x['upper'] for y in intervals if y is not x)
        x['worst_possible_rank']=1+sum(y['upper']>=x['lower'] for y in intervals if y is not x)
    summaries.append({'scenario':removed or 'complete','coverage':intervals[0]['coverage'],'interval_width':intervals[0]['upper']-intervals[0]['lower'],'rank_ambiguous_candidates':sum(x['best_possible_rank']!=x['worst_possible_rank'] for x in intervals),'possible_top5_candidates':sum(x['best_possible_rank']<=5 for x in intervals),'guaranteed_top5_candidates':sum(x['worst_possible_rank']<=5 for x in intervals)})
    details+=intervals
csvwrite('full827_missingness_ablation.csv',details);csvwrite('ablation_summary.csv',summaries)
summary={'rows':len(rows),'unique_exact_sequences':len(groups),'duplicate_groups':{h:a for h,a in groups.items() if len(a)>1},'invalid_sequence_records':sum(bool(r['invalid_residues']) for r in rows),'max_score_reconstruction_error':err,'fixed_weights':fields,'rank_file_sha256':hashlib.sha256(rank.read_bytes()).hexdigest(),'sequence_source_sha256':hashlib.sha256(b0.read_bytes()).hexdigest(),'scope':'Retrospective full-pool arithmetic missingness ablation; not fresh docking, not biological accuracy evaluation. Exact deduplication is not homology filtering. Sequence source is current local B0 mapped by frozen candidate_record_id; historical temporal provenance requires source records.'}
(out/'full827_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({'rows':len(rows),'unique_sequences':len(groups),'scenarios':len(summaries)}))
