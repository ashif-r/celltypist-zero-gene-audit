"""Prospective query masks, deterministic partitions and paired failure accounting."""
import numpy as np
import pandas as pd
from scipy import sparse

def detection_counts(x):
    if x.ndim!=2:raise ValueError('Expected two-dimensional expression')
    values=x.data if sparse.issparse(x) else np.asarray(x)
    if not np.isfinite(values).all() or (values<0).any():raise ValueError('Nonfinite or negative expression')
    return np.asarray((x>0).sum(axis=0)).ravel()

def retained_mask(x,minimum_cells=1):
    if not isinstance(minimum_cells,int) or minimum_cells<1:raise ValueError('Invalid detection threshold')
    return detection_counts(x)>=minimum_cells

def partitions(cell_ids,size,seed):
    ids=np.asarray(cell_ids,dtype=str)
    if len(set(ids))!=len(ids):raise ValueError('Duplicate cell IDs')
    if not isinstance(size,int) or size<1:raise ValueError('Invalid size')
    ordered=np.sort(ids)
    shuffled=np.random.default_rng(seed).permutation(ordered)
    return [shuffled[start:start+size] for start in range(0,len(ids),size)]

def paired_summary(cell_ids,left,right):
    ids=list(cell_ids)
    if len(set(ids))!=len(ids):raise ValueError('Duplicate eligible IDs')
    def valid(v):return isinstance(v,(str,np.str_)) and bool(v) and not pd.isna(v)
    pairs=[(left.get(i),right.get(i)) for i in ids]
    usable=[(a,b) for a,b in pairs if valid(a) and valid(b)]
    n=len(usable);changed=sum(a!=b for a,b in usable)
    return {'eligible':len(ids),'valid':n,'changed':changed,'failed':len(ids)-n,'fraction':changed/n if n else None,'complete':n==len(ids)}

def stability_certificate(scores,shift):
    scores=np.asarray(scores,dtype=float);shift=np.asarray(shift,dtype=float)
    if scores.ndim!=2 or scores.shape[1]<2 or shift.shape!=(scores.shape[1],):raise ValueError('Score/shift shape mismatch')
    if not np.isfinite(scores).all() or not np.isfinite(shift).all():raise ValueError('Nonfinite scores/shift')
    top=np.partition(scores,-2,axis=1)[:,-2:]
    margin=top[:,1]-top[:,0]
    return margin>np.ptp(shift)

def material_effect(rows,threshold=.01):
    return any(r['valid']>0 and r['fraction'] is not None and r['fraction']>=threshold for r in rows)
