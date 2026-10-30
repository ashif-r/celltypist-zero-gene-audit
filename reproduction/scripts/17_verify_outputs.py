"""Independent endpoint and algebra audit; optional exhaustive clean-run comparison."""
from pathlib import Path
import argparse,json,pickle,warnings
import numpy as np,pandas as pd
R=Path(__file__).resolve().parents[1];pa=argparse.ArgumentParser();pa.add_argument('--run',default='main');pa.add_argument('--compare',default=None);args=pa.parse_args();O=R/'results'/args.run
if not (O/'COMPLETED.json').exists():raise RuntimeError('Run incomplete')
with warnings.catch_warnings():
 warnings.simplefilter('ignore');m=pickle.loads((R/'data/Immune_All_Low.pkl').read_bytes())
clf=m['Model'];scaler=m['Scaler_'];names=np.asarray(clf.features);contrib=clf.coef_*scaler.mean_/scaler.scale_;q=pd.read_csv(O/'queries.csv');audits=[];expected_primary={}
for study,g in q.groupby('study'):
 poolid=g.loc[g.context.eq('pool'),'query_id'].iloc[0];pool=np.load(O/'queries'/poolid/'observed.npz');pool_indices={str(cid):i for i,cid in enumerate(pool['cell_ids'])}
 for row in g.itertuples():
  folder=O/'queries'/row.query_id;a=np.load(folder/'observed.npz');pred=np.load(folder/'predicted_B.npz');zero=set(pred['removed_genes'].tolist());ix=np.array([i for i,name in enumerate(names) if name in zero],dtype=int);shift=contrib[:,ix].sum(1)
  assert np.allclose(shift,pred['shift'],atol=1e-10,rtol=1e-10)
  b_expected=a['A']+shift
  assert np.allclose(a['B'],b_expected,atol=1e-6,rtol=1e-6);assert np.allclose(a['A'],a['C'],atol=1e-6,rtol=1e-6)
  changes=int(sum(int(i!=j) for i,j in zip(a['A'].argmax(1),a['B'].argmax(1))))
  pi=np.array([pool_indices[str(v)] for v in a['cell_ids']]);composition=int(sum(int(i!=j) for i,j in zip(pool['B'][pi].argmax(1),a['B'].argmax(1))))
  if row.context=='sample':expected_primary[(study,row.sample,'mechanism')]=(len(a['cell_ids']),changes);expected_primary[(study,row.sample,'composition')]=(len(a['cell_ids']),composition)
  audits.append({'query_id':row.query_id,'independent_shift_residual':float(np.abs(shift-pred['shift']).max()),'B_score_residual':float(np.abs(a['B']-b_expected).max()),'mechanism_changes':changes,'composition_changes':composition})
  if args.compare:
   other=R/'results'/args.compare/'queries'/row.query_id
   for name in ['observed.npz','predicted_B.npz','secondary_D.npz']:
    p=folder/name;s=other/name
    assert p.exists()==s.exists(),(row.query_id,name)
    if not p.exists():continue
    x=np.load(p);y=np.load(s);assert set(x.files)==set(y.files)
    for key in x.files:
     if x[key].dtype.kind in 'USbiu':assert np.array_equal(x[key],y[key]),(row.query_id,name,key)
     else:assert np.allclose(x[key],y[key],atol=1e-6,rtol=1e-6,equal_nan=True),(row.query_id,name,key)
p=pd.read_csv(O/'tables/primary.csv')
for row in p.itertuples():
 n,k=expected_primary[(row.study,row.sample,row.endpoint)];assert row.valid==n and row.eligible==n and row.changed==k and row.failed==0;assert abs(row.fraction-k/n)<1e-12
assert len(expected_primary)==len(p)
if args.compare:
 second=pd.read_csv(R/'results'/args.compare/'tables/primary.csv');assert p[['study','sample','endpoint','eligible','valid','changed','failed']].equals(second[['study','sample','endpoint','eligible','valid','changed','failed']])
pd.DataFrame(audits).to_csv(O/'tables/independent_audit.csv',index=False)
summary={'all_queries_independently_audited':len(audits),'primary_comparisons_verified':len(p),'all_valid_denominators_verified':True,'model_algebra_recomputed_without_primary_helpers':True,'clean_reproduction_compared':args.compare,'score_tolerances':{'atol':1e-6,'rtol':1e-6},'labels_masks_and_counts_exact_when_compared':bool(args.compare)}
(O/('reproduction_verification.json' if args.compare else 'independent_verification.json')).write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
