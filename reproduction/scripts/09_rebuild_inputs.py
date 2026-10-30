"""Recover original counts for frozen scIB-selected cells; certify shared-count equality."""
from pathlib import Path
import sys,json,gzip,hashlib
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np,pandas as pd,h5py,anndata as ad
from anndata.utils import make_index_unique
from scipy import sparse,io
from src.gate_core import identifier_keep,normalize_counts
from src.design import detection_counts
R=Path(__file__).resolve().parents[1];D=R/'data/geo_original';O=R/'results/input_repair';obs=pd.read_csv(R/'results/gate/main_obs.csv.gz',index_col=0);shared=pd.read_csv(R/'results/gate/main_var.csv.gz').iloc[:,0].astype(str).to_numpy()
artifacts=json.loads((R/'manifests/geo_original.json').read_text())['artifacts'];lookup={(x['batch'],x['kind']):x for x in artifacts}
reconciliations=[];mapping=[];studies=[]
with h5py.File(R/'data/Immune_ALL_human.h5ad') as f:
 for study in ['Oetjen','Sun']:
  batches=list(obs.loc[obs.study.eq(study),'batch'].unique());frames={};ambiguous=set()
  for batch in batches:
   frames[batch]=pd.read_csv(D/lookup[batch,'genes']['name'],sep='\t',header=None,names=['feature_id','symbol']);v=frames[batch]
   ambiguous.update(v.symbol[~identifier_keep(v.symbol)])
  base=frames[batches[0]]
  for batch in batches:assert frames[batch].equals(base),'Original feature schemas differ: explicit alignment required'
  keep=~base.symbol.isin(ambiguous).to_numpy();assert base.symbol[keep].is_unique
  # Reproduce the authors' explicit unique-name conversion ONLY to reconcile the old artifact.
  # These aliases never enter scoring; every ambiguous original symbol is excluded below.
  aliases=make_index_unique(pd.Index(base.symbol));shared_idx=aliases.get_indexer(shared);assert (shared_idx>=0).all()
  for j in np.flatnonzero(keep):mapping.append({'study':study,'original_row':int(j),'feature_id':base.feature_id.iloc[j],'symbol':base.symbol.iloc[j]})
  ys=[];os=[];raw_available=[]
  for batch in batches:
   mi=lookup[batch,'matrix'];bi=lookup[batch,'barcodes'];raw=io.mmread(gzip.open(D/mi['name'],'rb')).T.tocsr();raw.eliminate_zeros()
   barcodes=pd.read_csv(D/bi['name'],header=None)[0].astype(str);assert barcodes.is_unique and raw.shape==(len(barcodes),len(base))
   if not np.isfinite(raw.data).all() or (raw.data<0).any() or not np.equal(raw.data,np.floor(raw.data)).all():raise ValueError('Invalid original counts')
   loc=np.flatnonzero(obs.batch.eq(batch));sub=obs.iloc[loc].copy();suffix='-'+batch
   assert all(str(i).endswith(suffix) for i in sub.index)
   stripped=pd.Index([str(i)[:-len(suffix)] for i in sub.index]);ri=pd.Index(barcodes).get_indexer(stripped);assert (ri>=0).all() and len(set(ri))==len(ri)
   x=raw[ri,:].tocsr();n_mismatch=0;maxdiff=0.
   for start in range(0,len(loc),256):
    expected=f['layers/counts'][loc[start:start+256],:];got=x[start:start+256,shared_idx].toarray();diff=got-expected
    n_mismatch+=int(np.count_nonzero(diff));maxdiff=max(maxdiff,float(np.abs(diff).max()))
   assert n_mismatch==0,(batch,n_mismatch,maxdiff)
   y,totals=normalize_counts(x);sub['original_barcode']=stripped.to_numpy();sub['original_library_total']=totals
   ys.append(y[:,keep]);os.append(sub);raw_available.append(x[:,keep])
   det=detection_counts(x[:,keep]);reconciliations.append({'study':study,'batch':batch,'cells':len(loc),'original_features':len(base),'fixed_excluded_rows':int((~keep).sum()),'scoring_features':int(keep.sum()),'verified_shared_genes':len(shared),'compared_count_entries':len(loc)*len(shared),'count_mismatches':n_mismatch,'maximum_count_difference':maxdiff,'all_cells_matched':True,'all_zero_scoring_genes':int((det==0).sum()),'matrix_sha256':mi['sha256'],'barcodes_sha256':bi['sha256']})
   print('Reconciled',batch,len(loc),'cells x',len(shared),'shared genes; zero mismatches',flush=True)
  var=base.loc[keep].copy();var.index=pd.Index(var.symbol,name='gene_symbol');merged_obs=pd.concat(os);combined=sparse.vstack(ys,format='csr')
  out=R/'data/prepared'/f'{study}_original.h5ad';adata=ad.AnnData(combined,obs=merged_obs,var=var);adata.layers['counts']=sparse.vstack(raw_available,format='csr');adata.write_h5ad(out,compression='gzip')
  pooldet=detection_counts(adata.layers['counts']);studies.append({'study':study,'cells':adata.n_obs,'features':adata.n_vars,'original_features':len(base),'pool_all_zero_genes':int((pooldet==0).sum()),'prepared_file':str(out.relative_to(R)),'sha256':hashlib.sha256(out.read_bytes()).hexdigest(),'normalization':'log1p CP10k over every original recorded feature, before fixed duplicate exclusions'})
pd.DataFrame(reconciliations).to_csv(O/'geo_scib_reconciliation.csv',index=False);pd.DataFrame(mapping).to_csv(R/'manifests/original_scoring_features.csv',index=False)
(R/'manifests/prepared_inputs.json').write_text(json.dumps({'main_predictions_computed':False,'studies':studies},indent=2)+'\n');print(json.dumps(studies,indent=2))
