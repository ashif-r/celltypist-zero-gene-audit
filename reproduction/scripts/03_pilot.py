"""Exploratory PBMC3k implementation pilot; never computes main endpoints."""
from pathlib import Path
import sys,json,tarfile,hashlib,datetime,warnings,copy
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np,pandas as pd
from scipy import sparse,io
from scipy.special import expit
from anndata import AnnData
from celltypist.models import Model
from celltypist.classifier import Classifier
from src.gate_core import identifier_keep,normalize_counts,scores,zero_shift
R=Path(__file__).resolve().parents[1]; O=R/'results/gate/pilot';O.mkdir(exist_ok=True)
# The model-only screen is required and must predate expression predictions.
assert (R/'results/gate/checkpoint.json').is_file()
with tarfile.open(R/'data/pbmc3k.tar.gz') as t:
 def member(s):return next(m for m in t.getmembers() if m.name.endswith(s))
 f=pd.read_csv(t.extractfile(member('/genes.tsv')),sep='\t',header=None,names=['feature_id','symbol'])
 ids=pd.read_csv(t.extractfile(member('/barcodes.tsv')),header=None)[0].astype(str).to_numpy()
 x=io.mmread(t.extractfile(member('/matrix.mtx'))).T.tocsr().astype(np.float64)
assert x.shape==(len(ids),len(f)) and len(set(ids))==len(ids)
model=Model.load(str(R/'data/Immune_All_Low.pkl'))
assert model.scaler.with_mean and np.min(model.scaler.mean_)>=0
keep=identifier_keep(f.symbol);f['excluded']=~keep;f['in_model']=f.symbol.isin(model.features)
f.to_csv(O/'identifier_manifest.csv',index=False)
y,totals=normalize_counts(x)
removed_mass=np.asarray(x[:,~keep].sum(axis=1)).ravel()
np.savez_compressed(O/'normalization.npz',cell_ids=ids,original_library_counts=totals,excluded_counts=removed_mass)
y=y[:,keep].tocsr();genes=f.symbol.to_numpy()[keep]
assert len(set(genes))==len(genes)
controls=[]
def check(name,a,b):
 a=np.asarray(a);b=np.asarray(b)
 residual=float(np.max(np.abs(a-b))) if a.size else 0.
 ok=bool(np.allclose(a,b,atol=1e-6,rtol=1e-6))
 controls.append({'control':name,'passed':ok,'max_absolute_residual':residual})
 pd.DataFrame(controls).to_csv(O/'controls.csv',index=False)
 if not ok:raise AssertionError(name+': '+str(residual))
def api(z,gn,cellids,dense=False):
 out=[];labs=[];conf=[]
 for start in range(0,z.shape[0],512):
  v=z[start:start+512].toarray() if dense else z[start:start+512].copy()
  a=AnnData(v,obs=pd.DataFrame(index=cellids[start:start+512]),var=pd.DataFrame(index=gn))
  c=Classifier(a,model=copy.deepcopy(model))
  assert a.raw is None
  np.testing.assert_array_equal(c.indata_genes,gn)
  result=c.celltype();out.append(result.decision_matrix.to_numpy())
  labs.extend(result.predicted_labels.iloc[:,0].astype(str));conf.extend(result.probability_matrix.max(axis=1))
 return np.vstack(out),np.asarray(labs),np.asarray(conf)
A,la,ca=api(y,genes,ids)
np.savez_compressed(O/'reference.npz',scores=A,labels=la,confidence=ca,cell_ids=ids,classes=model.cell_types)
for start in range(0,len(ids),512):check('independent_A_'+str(start),A[start:start+512],scores(y[start:start+512],genes,model))
rows=[]
# Predetermined pilot queries: full PBMC pool and first 32/128/512 archive-order cells.
# These implementation probes are not natural samples or main-study endpoints.
for size in [len(ids),32,128,512]:
 q=y[:size]; qids=ids[:size];mask=np.asarray(q.getnnz(axis=0)>0).ravel()
 removed=genes[~mask];pred=A[:size]+zero_shift(removed,model)
 np.savez_compressed(O/f'predicted_B_{size}.npz',scores=pred,removed_genes=removed,cell_ids=qids)
 # Independent API inference only after saving analytic prediction.
 B,lb,cb=api(q[:,mask],genes[mask],qids)
 check('analytic_B_'+str(size),pred,B)
 expected=model.cell_types[pred.argmax(axis=1)]
 labels_agree=bool(np.array_equal(expected,lb))
 kept=q[:,mask].tocoo()
 restored=sparse.csr_matrix((kept.data,(kept.row,np.flatnonzero(mask)[kept.col])),shape=q.shape)
 assert (restored!=q).nnz==0
 C,lc,cc=api(restored,genes,qids);check('restore_C_'+str(size),A[:size],C)
 check('confidence_B_'+str(size),expit(pred.max(axis=1)),cb)
 np.savez_compressed(O/f'observed_B_{size}.npz',scores=B,labels=lb,confidence=cb)
 pd.DataFrame({'cell_id':qids,'A_label':la[:size],'B_label':lb,'predicted_B_label':expected,'A_confidence':ca[:size],'B_confidence':cb}).to_csv(O/f'paired_{size}.csv',index=False)
 rows.append({'query_cells':size,'removed_genes':len(removed),'removed_model_genes':int(np.isin(removed,model.features).sum()),'observed_flips':int((la[:size]!=lb).sum()),'predicted_flips':int((la[:size]!=expected).sum()),'destination_mismatches':int((lb!=expected).sum()),'score_max_residual':float(np.abs(pred-B).max())})
 if not labels_agree:raise AssertionError('Pilot label mismatch requires tie/alignment investigation')
q=y[:32];qids=ids[:32]
# Gene order, dense/sparse, fresh checkpoint/reversed order, fixed-schema query size.
perm=np.arange(len(genes))[::-1]
check('gene_order',A[:32],api(q[:,perm],genes[perm],qids)[0])
check('dense_sparse',A[:32],api(q,genes,qids,True)[0])
model=Model.load(str(R/'data/Immune_All_Low.pkl'))
reverse_mask=np.asarray(q.getnnz(axis=0)>0).ravel()
reverse_B=api(q[:,reverse_mask],genes[reverse_mask],qids)[0]
reverse_A=api(q,genes,qids)[0]
check('fresh_checkpoint_fixed_schema_reversed_order',A[:32],reverse_A)
check('reversed_B',reverse_B,np.load(O/'observed_B_32.npz')['scores'])
nonmodel=~np.isin(genes,model.features)&(np.asarray(q.getnnz(axis=0)).ravel()==0)
check('zero_nonmodel_deletion',A[:32],api(q[:,~nonmodel],genes[~nonmodel],qids)[0])
# D is secondary and may fail the API input contract independently.
dmask=np.asarray(q.getnnz(axis=0)>=3).ravel()
try:
 D,ld,cd=api(q[:,dmask],genes[dmask],qids)
 # Independent per-cell score calculation retaining exactly D's features.
 check('D_direct_score',D,scores(q[:,dmask],genes[dmask],model))
 lookup={g:i for i,g in enumerate(model.features)}
 j=np.array([j for j,g in enumerate(genes) if not dmask[j] and g in lookup]);mi=np.array([lookup[genes[k]] for k in j])
 scaled=np.minimum((q[:,j].toarray()-model.scaler.mean_[mi])/model.scaler.scale_[mi],10)
 check('D_shift',D,A[:32]-scaled @ model.classifier.coef_[:,mi].T)
 dstatus='passed'
except Exception as exc:dstatus=repr(exc)
pd.DataFrame(rows).to_csv(O/'pilot_summary.csv',index=False)
summary={'completed_at_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'cells':len(ids),'original_features':len(f),'excluded_feature_rows':int((~keep).sum()),'excluded_distinct_symbols':int(f.loc[~keep,'symbol'].nunique()),'scoring_features':len(genes),'model_overlap':int(np.isin(genes,model.features).sum()),'library_totals_preserved':True,'controls_passed':all(c['passed'] for c in controls),'primary_pilot_validated':True,'D_status':dstatus,'main_endpoints_computed':False,'queries':rows}
(O/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2),flush=True)
