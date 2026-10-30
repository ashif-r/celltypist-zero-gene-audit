"""Secondary reduced-schema endpoint tables and exhaustive archive reproduction check."""
from pathlib import Path
import json
import numpy as np,pandas as pd
from scipy.special import expit
R=Path(__file__).resolve().parents[1];O=R/'results/scib_secondary';S=R/'results/scib_secondary_reproduction'
assert (O/'COMPLETED.json').exists() and (S/'COMPLETED.json').exists()
q=pd.read_csv(O/'queries.csv');primary=[];drows=[];archives=0;arrays=0
for study,g in q.groupby('study'):
 poolrow=g.loc[g.context.eq('pool')].iloc[0];pool=np.load(O/'queries'/poolrow.query_id/'observed.npz');index={v:i for i,v in enumerate(pool['cell_ids'])}
 for row in g.itertuples():
  folder=O/'queries'/row.query_id;d=np.load(folder/'observed.npz');ids=d['cell_ids'];ix=np.array([index[v] for v in ids]);pr=np.load(folder/'predicted_B.npz')
  assert np.array_equal(d['A'],pool['A'][ix]);assert np.allclose(d['B'],pr['scores'],atol=1e-6,rtol=1e-6);assert np.array_equal(d['A'],d['C'])
  if row.context=='sample':
   for endpoint,left in [('mechanism',d['A']),('composition',pool['B'][ix])]:
    changed=sum(int(a!=b) for a,b in zip(left.argmax(1),d['B'].argmax(1)))
    primary.append({'study':study,'sample':row.sample,'endpoint':endpoint,'eligible':len(ids),'valid':len(ids),'failed':0,'changed':changed,'fraction':changed/len(ids),'role':'secondary_reduced_input'})
  dp=folder/'secondary_D.npz'
  if dp.exists():
   D=np.load(dp)['scores'];changed=int((D.argmax(1)!=d['A'].argmax(1)).sum());drows.append({'study':study,'sample':row.sample,'context':row.context,'valid':len(ids),'changed':changed,'fraction':changed/len(ids),'status':'passed'})
  else:drows.append({'study':study,'sample':row.sample,'context':row.context,'status':row.D_status})
  for name in ['observed.npz','predicted_B.npz','secondary_D.npz']:
   a=folder/name;b=S/'queries'/row.query_id/name;assert a.exists()==b.exists()
   if not a.exists():continue
   l=np.load(a);r=np.load(b);assert set(l.files)==set(r.files)
   for key in l.files:np.testing.assert_array_equal(l[key],r[key]);arrays+=1
   archives+=1
T=O/'tables';T.mkdir(exist_ok=True);pd.DataFrame(primary).to_csv(T/'primary.csv',index=False);pd.DataFrame(drows).to_csv(T/'secondary_D.csv',index=False)
report={'queries':len(q),'primary_comparisons':len(primary),'all_class_score_arrays_exact_on_reproduction':True,'archives_compared':archives,'arrays_compared':arrays,'max_B_score_residual':float(q.B_max_score_residual.max()),'max_C_score_residual':float(q.C_max_score_residual.max()),'B_destination_mismatches':int(q.B_destination_mismatches.sum()),'near_tie_B_count':int(q.B_near_ties.sum()),'interpretation':'Sensitivity of reduced input; available-feature normalization differs from original GEO input. No sole-cause attribution for cross-representation differences.'}
(O/'secondary_verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
