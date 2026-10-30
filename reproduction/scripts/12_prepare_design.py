"""Freeze query memberships and deletion-mask digests without annotation."""
from pathlib import Path
import sys,json,hashlib
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np,pandas as pd,anndata as ad
from src.design import partitions,retained_mask
R=Path(__file__).resolve().parents[1];config=json.loads((R/'config/analysis.json').read_text());prepared=json.loads((R/'manifests/prepared_inputs.json').read_text());plan=[];cells=[];samples=[]
for item in prepared['studies']:
 a=ad.read_h5ad(R/item['prepared_file']);ids=a.obs_names.to_numpy(dtype=str);index={v:i for i,v in enumerate(ids)};study=item['study']
 units=[('pool',ids,None,None)]+[('sample',a.obs_names[a.obs.batch.eq(b)].to_numpy(dtype=str),None,None) for b in sorted(a.obs.batch.unique())]
 for size in config['contexts']['stress_partition_sizes']:
  for seed in config['contexts']['stress_seeds']:
   units.extend(('stress',members,size,seed) for members in partitions(ids,size,seed))
 for number,(context,members,size,seed) in enumerate(units):
  positions=np.array([index[v] for v in members]);mask=retained_mask(a.X[positions]);qid=f'{study}_{number:05d}'
  plan.append({'query_id':qid,'study':study,'context':context,'size':size,'seed':seed,'cells':members.tolist(),'B_mask_sha256':hashlib.sha256(mask.tobytes()).hexdigest()})
 for cid,row in a.obs.iterrows():cells.append({'cell_id':str(cid),'study':study,'sample':str(row.batch),'original_barcode':str(row.original_barcode),'normalization_total':float(row.original_library_total)})
 for sample,n in a.obs.batch.value_counts().items():samples.append({'study':study,'sample':str(sample),'cells':int(n),'prepared_file':item['prepared_file']})
(R/'manifests/query_plan.json').write_text(json.dumps(plan,separators=(',',':'))+'\n');pd.DataFrame(cells).to_csv(R/'manifests/cells.csv',index=False);pd.DataFrame(samples).to_csv(R/'manifests/samples.csv',index=False)
print('Prospective plan:',len(plan),'queries,',len(cells),'cells; no annotations computed.')
