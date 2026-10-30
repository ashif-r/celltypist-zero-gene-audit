"""Independent reconstruction audit directly from GEO files, without preprocessing helpers."""
from pathlib import Path
from collections import Counter
import json,gzip,hashlib
import numpy as np,pandas as pd,anndata as ad
from scipy import io,sparse
R=Path(__file__).resolve().parents[1];sources=json.loads((R/'manifests/geo_original.json').read_text())['artifacts'];sources={(v['batch'],v['kind']):v for v in sources};records=[]
def file(batch,kind):
 item=sources[batch,kind];p=R/'data/geo_original'/item['name'];assert hashlib.sha256(p.read_bytes()).hexdigest()==item['sha256'];return p
for study in ['Oetjen','Sun']:
 a=ad.read_h5ad(R/f'data/prepared/{study}_original.h5ad')
 for batch in a.obs.batch.unique():
  features=pd.read_csv(file(batch,'genes'),sep='\t',header=None,names=['id','symbol']);symbols=features.symbol.astype(str).tolist();counts=Counter(symbols);keep=np.array([counts[v]==1 and v not in ['', 'nan'] for v in symbols]);np.testing.assert_array_equal(np.array(symbols)[keep],a.var_names.to_numpy());np.testing.assert_array_equal(features.id.astype(str).to_numpy()[keep],a.var.feature_id.astype(str).to_numpy())
  raw=io.mmread(gzip.open(file(batch,'matrix'),'rb')).T.tocsr();raw.sum_duplicates();raw.eliminate_zeros();barcodes=pd.read_csv(file(batch,'barcodes'),header=None)[0].astype(str);assert barcodes.is_unique
  positions=np.flatnonzero(a.obs.batch.eq(batch));ids=a.obs.original_barcode.iloc[positions].astype(str);index=pd.Index(barcodes).get_indexer(ids);assert (index>=0).all() and len(set(index))==len(index);raw=raw[index];totals=np.asarray(raw.sum(1)).ravel();np.testing.assert_array_equal(totals,a.obs.original_library_total.iloc[positions].to_numpy());assert (totals>0).all()
  available=raw[:,keep];saved=a.layers['counts'][positions];assert (available!=saved).nnz==0
  # Independent elementwise calculation on sparse nonzeros; omitted entries are exactly zero.
  c=available.tocoo();expected=sparse.csr_matrix((np.log1p(c.data.astype(float)*10000/totals[c.row]),(c.row,c.col)),shape=available.shape);observed=a.X[positions];error=expected-observed;maxerr=float(np.abs(error.data).max()) if error.nnz else 0.;assert maxerr<1e-12
  records.append({'study':study,'sample':batch,'cells':len(positions),'all_original_rows_for_normalization':raw.shape[1],'scoring_columns':int(keep.sum()),'count_entries_checked_including_implicit_zeros':int(available.shape[0]*available.shape[1]),'count_mismatches':0,'maximum_normalized_value_difference':maxerr,'source_hashes_verified':True,'no_empty_libraries':True})
  print('Independently verified',batch,len(positions),'cells',flush=True)
report={'samples':records,'cells_verified':sum(v['cells'] for v in records),'scoring_count_entries_verified':sum(v['count_entries_checked_including_implicit_zeros'] for v in records),'all_original_normalization_totals_exact':True,'all_source_hashes_verified':True,'method':'Counter duplicate policy and direct raw-UMI elementwise log1p formula; no project preprocessing helper used'}
(R/'results/input_repair/independent_input_verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k!='samples'},indent=2))
