"""Numeric-only pre-feedback scoring; no identifiers, labels or outcome inputs."""
import csv
import numpy as np

FIELDS = ('target_context_score','target_organism_score','phage_origin_score',
          'docking_support_score','solubility_no_chaperone_score',
          'E6_5_solubility_aggregation_chaperone_score',
          'curated_no_external_chaperone_required_strength',
          'curated_soluble_expression_strength','global_engineering_residual_score')
WEIGHTS = np.array([.55,.20,.20,.05])

def load_numeric(path):
    with path.open(encoding='utf-8', newline='') as f:
        reader = csv.DictReader(f)
        if tuple(reader.fieldnames) != FIELDS:
            raise ValueError('Only the exact frozen numeric schema is accepted')
        a = np.array([[float(row[k]) for k in FIELDS] for row in reader])
    if not np.isfinite(a).all() or not ((a >= 0) & (a <= 1)).all():
        raise ValueError('Missing, nonfinite or out-of-range input')
    return a

def components(a, frozen_precision=False):
    # Archived composite was evaluated before the solubility subscore was rounded.
    sol = a[:,4] if frozen_precision else .60*a[:,5]+.25*a[:,6]+.15*a[:,7]
    return np.column_stack((a[:,0],a[:,3],sol,a[:,8]))

def scores(x, weights):
    w = np.asarray(weights, dtype=float)
    if w.shape != (4,) or (w < 0).any() or not np.isfinite(w).all() or w.sum() <= 0:
        raise ValueError('Invalid weights')
    w = w / w.sum()
    # Explicit left-to-right operations match archived Python float evaluation.
    return w[0]*x[:,0]+w[1]*x[:,1]+w[2]*x[:,2]+w[3]*x[:,3]

def rounded(values):
    return np.array([float(format(float(v), '.6f')) for v in values])

def anonymous_ties(n, seed):
    return np.random.default_rng(seed).permutation(n)

def rank(values, tie_keys):
    """Caller supplies only ordinal tie keys; never candidate-specific priorities."""
    order = sorted(range(len(values)), key=lambda i:(-float(values[i]),*tie_keys[i]))
    ranks = np.empty(len(order), dtype=int)
    ranks[order] = np.arange(1,len(order)+1)
    return np.array(order),ranks

def frozen_ties(a, lexical_ordinals):
    return [(-a[i,0],-a[i,1],-a[i,2],-a[i,3],-a[i,4],-a[i,8],int(lexical_ordinals[i]))
            for i in range(len(a))]

def analysis_ties(n, seed):
    return [(int(v),) for v in anonymous_ties(n,seed)]
