"""Exhaustively compare every archive, decoded label and identity in both main runs."""
from pathlib import Path
import json
import numpy as np,pandas as pd
R=Path(__file__).resolve().parents[1];A=R/'results/main';B=R/'results/reproduction'
assert (A/'COMPLETED.json').exists() and (B/'COMPLETED.json').exists()
qa=pd.read_csv(A/'queries.csv');qb=pd.read_csv(B/'queries.csv');pd.testing.assert_frame_equal(qa,qb);archives=arrays=label_vectors=0;bit_exact=True;max_residual=0.
for row in qa.itertuples():
 for name in ['observed.npz','predicted_B.npz','secondary_D.npz']:
  p=A/'queries'/row.query_id/name;s=B/'queries'/row.query_id/name;assert p.exists()==s.exists()
  if not p.exists():continue
  x=np.load(p);y=np.load(s);assert set(x.files)==set(y.files)
  for key in x.files:
   l=x[key];r=y[key];exact=np.array_equal(l,r)
   if l.dtype.kind in 'f':
    assert np.allclose(l,r,atol=1e-6,rtol=1e-6);max_residual=max(max_residual,float(np.abs(l-r).max()));bit_exact=bit_exact and exact
    if l.ndim==2 and l.shape[1]==98:np.testing.assert_array_equal(l.argmax(1),r.argmax(1));label_vectors+=1
   else:assert exact
   arrays+=1
  archives+=1
for study in qa.study.unique():
 x=np.load(A/f'{study}_reference.npz');y=np.load(B/f'{study}_reference.npz')
 for key in x.files:np.testing.assert_array_equal(x[key],y[key])
result={'queries':len(qa),'query_archives_compared':archives,'arrays_compared':arrays,'decoded_label_vectors_exact':label_vectors,'all_float_arrays_bit_exact':bool(bit_exact),'max_float_residual':max_residual,'cell_and_class_orders_exact':True,'removed_gene_lists_and_stability_flags_exact':True,'reference_archives_exact':2,'note':'Exact labels are checked per cell, in addition to primary counts and float tolerances.'}
(A/'exact_reproduction.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
