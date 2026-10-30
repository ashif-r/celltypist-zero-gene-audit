"""Independent float64 score reconstruction and fixed identifier policy."""
import numpy as np
import pandas as pd
from scipy import sparse

def identifier_keep(symbols):
    s=pd.Series(np.asarray(symbols,dtype=str))
    return (~s.duplicated(keep=False) & s.ne('') & s.ne('nan')).to_numpy()

def normalize_counts(x):
    x=sparse.csr_matrix(x,dtype=np.float64)
    x.eliminate_zeros()
    if not np.isfinite(x.data).all() or (x.data<0).any(): raise ValueError('Invalid counts')
    totals=np.asarray(x.sum(axis=1)).ravel()
    if (totals<=0).any(): raise ValueError('Empty library; exclude and report before normalization')
    y=x.multiply((10000/totals)[:,None]).tocsr(); y.data=np.log1p(y.data)
    return y,totals

def scores(x,genes,model):
    lookup={g:i for i,g in enumerate(model.features)}
    columns=np.array([j for j,g in enumerate(genes) if g in lookup],dtype=int)
    idx=np.array([lookup[genes[j]] for j in columns],dtype=int)
    if len(idx)==0: raise ValueError('No model features')
    z=x[:,columns].toarray() if sparse.issparse(x) else np.asarray(x)[:,columns].copy()
    z=z.astype(np.float64)
    if model.scaler.with_mean:z-=model.scaler.mean_[idx]
    z/=model.scaler.scale_[idx]; np.minimum(z,10,out=z)
    return z @ model.classifier.coef_[:,idx].T + model.classifier.intercept_

def zero_shift(removed,model):
    removed_set=set(removed)
    idx=np.asarray([i for i,g in enumerate(model.features) if g in removed_set],dtype=int)
    if not model.scaler.with_mean:return np.zeros(len(model.cell_types))
    zero=np.minimum(-model.scaler.mean_[idx]/model.scaler.scale_[idx],10)
    return -(model.classifier.coef_[:,idx]*zero).sum(axis=1)
