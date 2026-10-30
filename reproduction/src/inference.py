"""Native public API inference and independent score checks, with fixed query masks."""
import copy,logging
import numpy as np,pandas as pd
from scipy import sparse
from scipy.special import expit
from anndata import AnnData
from .gate_core import scores,zero_shift
from .design import retained_mask,stability_certificate

ATOL=RTOL=1e-6

def compare(expected,observed,name):
    if expected.shape!=observed.shape or not np.isfinite(observed).all():raise AssertionError(name+': invalid scores')
    residual=float(np.max(np.abs(expected-observed)))
    if not np.allclose(expected,observed,atol=ATOL,rtol=RTOL):raise AssertionError(name+': '+str(residual))
    return residual

def api_scores(x,genes,ids,model,chunk=512):
    from celltypist import annotate,logger
    logger.set_level(logging.ERROR)
    if len(set(genes))!=len(genes) or len(set(ids))!=len(ids):raise ValueError('Nonunique input IDs')
    # Validate the original matrix before removing features the native API ignores.
    values=x.data if sparse.issparse(x) else np.asarray(x)
    if not np.isfinite(values).all() or (values<0).any() or (values>9.22).any():
        raise ValueError('Invalid expression matrix before model-feature selection')
    model_genes=set(model.features)
    matched=np.array([g in model_genes for g in genes])
    if not matched.any():raise ValueError('No features overlap with the model')
    # No normalization or query-mask calculation occurs here. Preserve gene order.
    x=x[:,matched];genes=genes[matched]
    result=[];labels=[]
    for start in range(0,len(ids),chunk):
        a=AnnData(x[start:start+chunk].copy(),obs=pd.DataFrame(index=ids[start:start+chunk]),var=pd.DataFrame(index=genes))
        r=annotate(a,model=copy.deepcopy(model),majority_voting=False,mode='best match')
        assert a.raw is None
        np.testing.assert_array_equal(r.decision_matrix.index.to_numpy(),ids[start:start+chunk])
        np.testing.assert_array_equal(r.decision_matrix.columns.to_numpy(),model.cell_types)
        result.append(r.decision_matrix.to_numpy());labels.extend(r.predicted_labels.iloc[:,0].astype(str))
    out=np.vstack(result);native=np.asarray(labels,dtype=str)
    np.testing.assert_array_equal(native,np.asarray(model.cell_types,dtype=str)[out.argmax(axis=1)])
    return out

def direct_scores(x,genes,model,chunk=512):
    return np.vstack([scores(x[start:start+chunk],genes,model) for start in range(0,x.shape[0],chunk)])

def restore_zeros(filtered,mask,shape):
    c=filtered.tocoo();return sparse.csr_matrix((c.data,(c.row,np.flatnonzero(mask)[c.col])),shape=shape)

def query(x,genes,ids,model,A,path,do_D=False):
    """Predict and save B before independent native B; physically restore C."""
    path.mkdir(parents=True,exist_ok=False)
    mask=retained_mask(x);removed=genes[~mask];shift=zero_shift(removed,model);predicted=A+shift
    np.savez_compressed(path/'predicted_B.npz',scores=predicted,shift=shift,removed_genes=removed,cell_ids=ids,classes=np.asarray(model.cell_types,dtype=str))
    B=api_scores(x[:,mask],genes[mask],ids,model);b_residual=compare(predicted,B,'Analytic B')
    expected=predicted.argmax(axis=1);observed=B.argmax(axis=1)
    ordered=np.sort(predicted,axis=1);eps=ATOL+RTOL*np.abs(predicted).max(axis=1);near=(ordered[:,-1]-ordered[:,-2])<=2*eps
    mismatch=expected!=observed
    if np.any(mismatch&~near):raise AssertionError('Unexplained B label destination')
    restored=restore_zeros(x[:,mask],mask,x.shape)
    if (restored!=x).nnz:raise AssertionError('Restoration changed expression')
    C=api_scores(restored,genes,ids,model);c_residual=compare(A,C,'C recovery')
    # Invariant labels must agree away from the A tie band.
    atop=np.sort(A,axis=1);a_near=(atop[:,-1]-atop[:,-2])<=2*(ATOL+RTOL*np.abs(A).max(axis=1))
    if np.any((A.argmax(1)!=C.argmax(1))&~a_near):raise AssertionError('C label invariance failed')
    safe=stability_certificate(A,shift)
    if np.any(safe&(A.argmax(1)!=observed)):raise AssertionError('Stability certificate violated')
    np.savez_compressed(path/'observed.npz',A=A,B=B,C=C,cell_ids=ids,classes=np.asarray(model.cell_types,dtype=str),certified=safe,near_tie_B=near,label_prediction_mismatch=mismatch)
    dstatus='not_requested';d_residual=None
    if do_D:
        dm=retained_mask(x,3)
        try:
            D=api_scores(x[:,dm],genes[dm],ids,model)
            d_pred=direct_scores(x[:,dm],genes[dm],model)
            d_residual=compare(d_pred,D,'D independent scores')
            # Separately reconstruct A minus each removed feature's scaled/clipped term.
            lookup={g:i for i,g in enumerate(model.features)}
            j=np.asarray([j for j,g in enumerate(genes) if not dm[j] and g in lookup],dtype=int)
            mi=np.asarray([lookup[genes[jj]] for jj in j],dtype=int)
            for start in range(0,len(ids),512):
                z=x[start:start+512,j].toarray().astype(float)
                if model.scaler.with_mean:z-=model.scaler.mean_[mi]
                z/=model.scaler.scale_[mi];np.minimum(z,10,out=z)
                compare(A[start:start+512]-z@model.classifier.coef_[:,mi].T,D[start:start+512],'D cell-specific shift')
            np.savez_compressed(path/'secondary_D.npz',scores=D,predicted=d_pred,removed_genes=genes[~dm]);dstatus='passed'
        except ValueError as e:
            if not any(v in str(e) for v in ['Invalid expression matrix','No features overlap']):raise
            dstatus='input_contract_failure: '+str(e)
        # Assertion failures are numerical failures, not permissible input failures.
    model_genes=set(model.features)
    return {'B_max_score_residual':b_residual,'C_max_score_residual':c_residual,'B_destination_mismatches':int(mismatch.sum()),'B_near_ties':int(near.sum()),'removed_genes':len(removed),'removed_model_genes':sum(g in model_genes for g in removed),'certified_stable':int(safe.sum()),'shift_min':float(shift.min()),'shift_max':float(shift.max()),'shift_spread':float(np.ptp(shift)),'D_status':dstatus,'D_max_score_residual':d_residual}
