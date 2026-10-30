from pathlib import Path
import json,numpy as np,pandas as pd
R=Path(__file__).resolve().parents[1];O=R/'results/gate/pilot';rows=[]
ref=np.load(O/'reference.npz',allow_pickle=True)
for n in [2700,32,128,512]:
 p=pd.read_csv(O/f'paired_{n}.csv');a=np.load(O/f'predicted_B_{n}.npz')['scores'];b=np.load(O/f'observed_B_{n}.npz')['scores']
 assert len(p)==n and p.cell_id.is_unique
 assert (p.B_label==p.predicted_B_label).all()
 assert np.allclose(a,b,atol=1e-6,rtol=1e-6)
 flips=int(sum(x!=y for x,y in zip(p.A_label,p.B_label)))
 base_full=np.load(O/'observed_B_2700.npz')['labels'][:n]
 composition=int(sum(x!=y for x,y in zip(base_full,p.B_label)))
 epsilon=1e-6+1e-6*np.abs(b).max(axis=1);ordered=np.sort(b,axis=1)
 ties=int(((ordered[:,-1]-ordered[:,-2])<=2*epsilon).sum())
 row={'cells':n,'flips':flips,'fraction':flips/n,'pooled_vs_subquery_flips':composition,'B_near_ties':ties,'max_score_residual':float(np.abs(a-b).max())}
 for t in [.5,.9,.95]:row[f'confidence_crossings_{t}']=int(((p.A_confidence>=t)!=(p.B_confidence>=t)).sum())
 rows.append(row)
pd.DataFrame(rows).to_csv(O/'independent_verification.csv',index=False)
sample=pd.read_csv(R/'results/gate/sample_table.csv');eligible=sample.loc[sample.mechanism_eligible.eq(True)]
x=json.loads((R/'results/gate/input_audit.json').read_text());x['eligible_cells']=int(eligible.cells.sum());assert x['eligible_cells']==18410
x['eligible_natural_samples']=len(eligible);x['eligible_multisample_studies']=int((eligible.groupby('study').size()>=2).sum())
(R/'results/gate/input_audit.json').write_text(json.dumps(x,indent=2)+'\n')
print(pd.DataFrame(rows).to_string(index=False));print('Independent eligible count:',x['eligible_cells'])
