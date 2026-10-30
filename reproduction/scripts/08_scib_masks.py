"""Input-only masks: no model prediction or annotation endpoints."""
from pathlib import Path
import h5py,numpy as np,pandas as pd,json
R=Path(__file__).resolve().parents[1];obs=pd.read_csv(R/'results/gate/main_obs.csv.gz',index_col=0);g=pd.read_csv(R/'results/gate/main_var.csv.gz').iloc[:,0].to_numpy();rows=[]
with h5py.File(R/'data/Immune_ALL_human.h5ad') as f:
 d=f['layers/counts']
 for study in ['Oetjen','Sun']:
  units=list(obs.loc[obs.study.eq(study),'batch'].unique())+[study+'_pool']
  for unit in units:
   idx=np.flatnonzero(obs.study.eq(study) if unit.endswith('_pool') else obs.batch.eq(unit));det=np.zeros(len(g),dtype=np.int64)
   for start in range(0,len(idx),256):det+=(d[idx[start:start+256],:]>0).sum(axis=0)
   rows.append({'study':study,'query':unit,'cells':len(idx),'recorded_features':len(g),'all_zero_features':int((det==0).sum()),'min_cells3_removed':int((det<3).sum())})
pd.DataFrame(rows).to_csv(R/'results/input_repair/scib_zero_masks.csv',index=False);print(pd.DataFrame(rows).to_string(index=False))
