"""Missing-evidence interval policy, proposed retrospectively; not a trained predictor.
Unknown values span [0,1] instead of being classified as measured failures.
Weights must be declared before examining candidate outcomes. Literature is kept
as a separate evidence register and is not automatically added to this core.
"""
import math
def summarize(values,weights):
    if set(values)!=set(weights):raise ValueError('Feature contract mismatch')
    if any(not math.isfinite(w) or w<=0 for w in weights.values()):raise ValueError('Weights must be finite and positive')
    total=sum(weights.values());lower=0.;coverage=0.;missing=[]
    for k,v in values.items():
        w=weights[k]/total
        if v is None:missing.append(k);continue
        if not math.isfinite(v) or not 0<=v<=1:raise ValueError('Feature outside [0,1]')
        lower+=w*v;coverage+=w
    return {'lower':lower,'upper':lower+1-coverage,'coverage':coverage,'observed_weighted_mean':lower/coverage if coverage else None,'missing':missing,'status':'complete' if not missing else 'incomplete_requires_evidence'}
def rank_bounds(records):
    for i,r in enumerate(records):
        r['best_possible_rank']=1+sum(o['lower']>r['upper'] for j,o in enumerate(records) if i!=j)
        r['worst_possible_rank']=1+sum(o['upper']>=r['lower'] for j,o in enumerate(records) if i!=j)
    return records
