"""Compare optimized native calls with every preserved initial-run query."""
from pathlib import Path
import sys,json,time,copy
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np,pandas as pd,anndata as ad
from concurrent.futures import ThreadPoolExecutor,ProcessPoolExecutor
import multiprocessing
from src.main_worker import initialize,run_unit
from celltypist.models import Model
from src.inference import api_scores
from src.design import retained_mask
def main():
 R=Path(__file__).resolve().parents[1];old=R/'results/main_initial_performance_run';m=Model.load(str(R/'data/Immune_All_Low.pkl'));rows=pd.read_csv(old/'queries.csv');datasets={};lookups={}
 for study in rows.study.unique():
  a=ad.read_h5ad(R/f'data/prepared/{study}_original.h5ad');datasets[study]=a;lookups[study]={v:i for i,v in enumerate(a.obs_names)}
 def check(row):
  a=datasets[row.study];saved=np.load(old/'queries'/row.query_id/'observed.npz');ids=saved['cell_ids'];ix=np.array([lookups[row.study][v] for v in ids]);x=a.X[ix];genes=a.var_names.to_numpy(dtype=str);mask=retained_mask(x)
  B=api_scores(x[:,mask],genes[mask],ids,m);C=api_scores(x,genes,ids,m)
  np.testing.assert_array_equal(B,saved['B']);np.testing.assert_array_equal(C,saved['C'])
  return {'query_id':row.query_id,'cells':len(ids),'B_C_scores_exact':True}
 t=time.perf_counter()
 with ThreadPoolExecutor(max_workers=4) as pool:checks=list(pool.map(check,rows.itertuples()))
 report={'preserved_initial_run':'main_initial_performance_run','queries_compared':len(checks),'all_native_B_C_scores_exact':True,'normalization_and_masks_unchanged':True,'workers':4,'elapsed_seconds':time.perf_counter()-t,'queries':checks}
 (R/'results/gate/implementation_equivalence.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k!='queries'},indent=2))

 # Exercise isolated worker scheduling on eight predeclared saved queries.
 plan=json.loads((R/'manifests/query_plan.json').read_text());tasks=[]
 for number,q in enumerate(plan[:8]):
  tasks.append((number,(q['context'],np.array(q['cells']),q.get('size'),q.get('seed'))))
 out=R/'results/gate/process_worker_control';out.mkdir(exist_ok=False)
 with ProcessPoolExecutor(max_workers=4,mp_context=multiprocessing.get_context('spawn'),initializer=initialize,initargs=(str(R),'Oetjen',str(R/'data/prepared/Oetjen_original.h5ad'),str(old/'Oetjen_reference.npz'),str(out),0.0)) as pool:
  controls=list(pool.map(run_unit,tasks))
 for _,_,row in controls:
  for name in ['observed.npz','predicted_B.npz','secondary_D.npz']:
   left=out/'queries'/row['query_id']/name;right=old/'queries'/row['query_id']/name
   assert left.exists()==right.exists()
   if left.exists():
    l=np.load(left);r=np.load(right)
    for key in l.files:np.testing.assert_array_equal(l[key],r[key])
 report['process_worker_queries_exact']=len(controls)
 (R/'results/gate/implementation_equivalence.json').write_text(json.dumps(report,indent=2)+'\n')
 print('Isolated four-process worker control passed for',len(controls),'queries')

if __name__=='__main__':main()
