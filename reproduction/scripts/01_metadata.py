"""Read H5AD metadata only, without annotation or loading expression."""
from pathlib import Path
import json,hashlib
import h5py
from anndata.io import read_elem
R=Path(__file__).resolve().parents[1];p=R/'data/Immune_ALL_human.h5ad'
with h5py.File(p,'r') as f:
 obs=read_elem(f['obs']);var=read_elem(f['var'])
 obs.to_csv(R/'results/gate/main_obs.csv.gz')
 var.to_csv(R/'results/gate/main_var.csv.gz')
 summary={'cells':len(obs),'genes':len(var),'obs_columns':list(obs.columns),'var_columns':list(var.columns),'layers':list(f.get('layers',{})),'raw_present':'raw' in f,'batch_counts':obs['batch'].value_counts().to_dict() if 'batch' in obs else {},'identifier_provenance':'requires_original_source_feature_table_mapping','annotation_endpoints_computed':False}
 (R/'results/gate/main_metadata.json').write_text(json.dumps(summary,indent=2,default=str)+'\n')
 print(json.dumps(summary,indent=2,default=str))
 for c in ['batch','sample_ID','study','dataset','tissue']:
  if c in obs:print(c,obs[c].value_counts().to_dict())
