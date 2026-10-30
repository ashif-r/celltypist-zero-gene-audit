"""Labelled secondary input sensitivity; cannot substitute for original-count main data."""
from pathlib import Path
import sys,json,argparse
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from src.freeze import require_freeze
R=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser();p.add_argument('--prepare',action='store_true');p.add_argument('--reproduction',action='store_true');args=p.parse_args()
if not args.prepare:require_freeze(R)
import h5py,numpy as np,pandas as pd,anndata as ad
from scipy import sparse
from src.gate_core import normalize_counts
from src.inference import api_scores,direct_scores,compare,query
from src.design import retained_mask
import hashlib
obs=pd.read_csv(R/'results/gate/main_obs.csv.gz',index_col=0);genes=pd.read_csv(R/'results/gate/main_var.csv.gz').iloc[:,0].astype(str).to_numpy();mapping=pd.read_csv(R/'results/gate/main_original_feature_mapping.csv.gz');manifest=[];plan=[]
if args.prepare:
 with h5py.File(R/'data/Immune_ALL_human.h5ad') as f:
  for study in ['Oetjen','Sun']:
   loc=np.flatnonzero(obs.study.eq(study));excluded=set(mapping.loc[mapping.study.eq(study)&mapping.excluded.eq(True),'recorded_symbol']);keep=~np.isin(genes,list(excluded));chunks=[];totals=[]
   for start in range(0,len(loc),256):
    y,t=normalize_counts(f['layers/counts'][loc[start:start+256],:]);chunks.append(y[:,keep]);totals.extend(t)
   o=obs.iloc[loc].copy();o['available_feature_library_total']=totals;a=ad.AnnData(sparse.vstack(chunks,format='csr'),obs=o,var=pd.DataFrame(index=genes[keep]));out=R/'data/prepared'/f'{study}_scib.h5ad';a.write_h5ad(out,compression='gzip')
   manifest.append({'study':study,'prepared_file':str(out.relative_to(R)),'sha256':hashlib.sha256(out.read_bytes()).hexdigest(),'normalization':'available 12303 scIB features, before fixed identifier exclusions','role':'secondary_input_sensitivity_not_independent_dataset'})
   ids=a.obs_names.to_numpy(dtype=str);index={v:i for i,v in enumerate(ids)};units=[('pool',ids)]+[('sample',a.obs_names[a.obs.batch.eq(b)].to_numpy(dtype=str)) for b in sorted(a.obs.batch.unique())]
   for number,(context,members) in enumerate(units):
    pos=np.array([index[v] for v in members]);plan.append({'query_id':f'{study}_{number:05d}','study':study,'context':context,'cells':members.tolist(),'B_mask_sha256':hashlib.sha256(retained_mask(a.X[pos]).tobytes()).hexdigest()})
 (R/'manifests/secondary_inputs.json').write_text(json.dumps({'studies':manifest},indent=2)+'\n');(R/'manifests/secondary_query_plan.json').write_text(json.dumps(plan,separators=(',',':'))+'\n');print('Secondary inputs prepared without annotation.');sys.exit(0)
from celltypist.models import Model
model=Model.load(str(R/'data/Immune_All_Low.pkl'));items=json.loads((R/'manifests/secondary_inputs.json').read_text())['studies'];plan=json.loads((R/'manifests/secondary_query_plan.json').read_text());O=R/'results'/('scib_secondary_reproduction' if args.reproduction else 'scib_secondary');O.mkdir(exist_ok=False);rows=[]
for item in items:
 path=R/item['prepared_file'];assert hashlib.sha256(path.read_bytes()).hexdigest()==item['sha256'];a=ad.read_h5ad(path);ids=a.obs_names.to_numpy(dtype=str);genes=a.var_names.to_numpy(dtype=str);index={v:i for i,v in enumerate(ids)};A=api_scores(a.X,genes,ids,model);compare(direct_scores(a.X,genes,model),A,'Secondary A')
 for spec in [q for q in plan if q['study']==item['study']]:
  members=np.array(spec['cells']);ix=np.array([index[v] for v in members]);assert hashlib.sha256(retained_mask(a.X[ix]).tobytes()).hexdigest()==spec['B_mask_sha256'];qid=spec['query_id'];validation=query(a.X[ix],genes,members,model,A[ix],O/'queries'/qid,do_D=True);rows.append({'query_id':qid,'study':item['study'],'context':spec['context'],'cells':len(members),'sample':str(a.obs.loc[members[0],'batch']) if spec['context']=='sample' else '',**validation})
pd.DataFrame(rows).to_csv(O/'queries.csv',index=False);(O/'COMPLETED.json').write_text(json.dumps({'primary_controls_passed':True,'role':'secondary_input_sensitivity'},indent=2)+'\n')
