"""Frozen main experiment. Refuses inference without an authorized prospective record."""
from pathlib import Path
import sys,json,datetime,traceback,hashlib,argparse
from concurrent.futures import ProcessPoolExecutor
import multiprocessing
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from src.freeze import require_freeze,sha256
def main():
 R=Path(__file__).resolve().parents[1]
 parser=argparse.ArgumentParser();parser.add_argument('--run',default='main');args=parser.parse_args()
 if args.run not in ['main','reproduction']:raise ValueError('Unknown run name')
 record=require_freeze(R) # No checkpoint load, annotations or output directories before this guard.
 from src.main_worker import initialize,run_unit
 import numpy as np,pandas as pd,anndata as ad
 from scipy import sparse
 from celltypist.models import Model
 from src.inference import api_scores,direct_scores,compare,query
 from src.design import partitions,retained_mask
 config=json.loads((R/'config/analysis.json').read_text());prepared=json.loads((R/'manifests/prepared_inputs.json').read_text())
 O=R/'results'/args.run;O.mkdir(exist_ok=False)
 (O/'STARTED.json').write_text(json.dumps({'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'freeze':record,'run':args.run},indent=2)+'\n')
 model=Model.load(str(R/'data/Immune_All_Low.pkl'));validation=[];query_rows=[]
 prospective={q['query_id']:q for q in json.loads((R/'manifests/query_plan.json').read_text())}
 try:
  for item in prepared['studies']:
   p=R/item['prepared_file'];assert sha256(p)==item['sha256']
   a=ad.read_h5ad(p);genes=a.var_names.to_numpy(dtype=str);ids=a.obs_names.to_numpy(dtype=str);x=a.X.tocsr();idx={v:i for i,v in enumerate(ids)};study=item['study']
   print('Reference',study,len(ids),flush=True)
   A=api_scores(x,genes,ids,model)
   residual=compare(direct_scores(x,genes,model),A,'Independent A')
   np.savez_compressed(O/f'{study}_reference.npz',scores=A,cell_ids=ids,classes=np.asarray(model.cell_types,dtype=str))
   # Explicit fixed-schema controls on predetermined first 32 cells and model genes.
   n=min(32,len(ids));rev=np.arange(len(genes))[::-1]
   compare(A[:n],api_scores(x[:n,rev],genes[rev],ids[:n],model),'Main gene-order control')
   compare(A[:n],api_scores(sparse.csr_matrix(x[:n].toarray()),genes,ids[:n],model),'Main dense/sparse roundtrip')
   units=[('pool',ids,None,None)]+[('sample',a.obs_names[a.obs.batch.eq(b)].to_numpy(dtype=str),None,None) for b in sorted(a.obs.batch.unique())]
   for size in config['contexts']['stress_partition_sizes']:
    for seed in config['contexts']['stress_seeds']:
     for number,members in enumerate(partitions(ids,size,seed)):
      units.append(('stress',members,size,seed))
   # Bounded batches preserve planned output order and avoid queuing all sparse matrices.
   with ProcessPoolExecutor(max_workers=4,mp_context=multiprocessing.get_context('spawn'),initializer=initialize,initargs=(str(R),study,str(p),str(O/f'{study}_reference.npz'),str(O),residual)) as executor:
    for start in range(0,len(units),4):
     for number,validation_row,query_row in executor.map(run_unit,enumerate(units[start:start+4],start=start)):
      validation.append(validation_row);query_rows.append(query_row)
      pd.DataFrame(validation).to_csv(O/'validation.csv',index=False);pd.DataFrame(query_rows).to_csv(O/'queries.csv',index=False)
      if query_row['context']!='stress' or number%100==0:print(query_row['query_id'],query_row['context'],query_row['cells'],'validated',flush=True)
  (O/'COMPLETED.json').write_text(json.dumps({'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'primary_controls_passed':True,'queries':len(query_rows),'config_sha256':sha256(R/'config/analysis.json')},indent=2)+'\n')
 except Exception:
  (O/'FAILED.txt').write_text(traceback.format_exc());raise

if __name__=='__main__':main()
