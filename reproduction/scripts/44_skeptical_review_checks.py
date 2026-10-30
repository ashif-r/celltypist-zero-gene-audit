"""Post-analysis proof checks and predicted-count summaries; no new inference.

Reads immutable saved outputs. Reported abundances are model-assigned counts,
not validated biological abundances. Does not choose a new primary endpoint.
"""
from pathlib import Path
import json,sys,hashlib
import numpy as np,pandas as pd
R=Path(__file__).resolve().parents[1];sys.path.insert(0,str(R));from src.freeze import require_freeze
require_freeze(R);O=R/'results/reviewer_revision';O.mkdir(exist_ok=True)
c=pd.read_csv(R/'results/main/tables/per_cell_primary.csv.gz');q=pd.read_csv(R/'results/main/queries.csv')
rows=[];counts=[]
for ep,left in [('mechanism','A_label'),('composition','B_pool_label')]:
 transition=pd.read_csv(R/f'results/main/tables/{ep}_transitions.csv')
 for study,g in c.groupby('study'):
  labels=sorted(set(g[left])|set(g.B_sample_label));a=g[left].value_counts().reindex(labels,fill_value=0);b=g.B_sample_label.value_counts().reindex(labels,fill_value=0);delta=b-a;changed=int(g[left].ne(g.B_sample_label).sum());mass=int(delta.abs().sum()//2);assert delta.sum()==0 and mass<=changed
  tg=transition.loc[transition.study.eq(study)];ta=tg.groupby(left).cells.sum().reindex(labels,fill_value=0);tb=tg.groupby('B_sample_label').cells.sum().reindex(labels,fill_value=0);assert np.array_equal(a.to_numpy(),ta.to_numpy()) and np.array_equal(b.to_numpy(),tb.to_numpy())
  for label in labels:counts.append({'endpoint':ep,'study':study,'label':label,'reference_count':int(a[label]),'sample_B_count':int(b[label]),'count_delta':int(delta[label]),'denominator':len(g)})
  rows.append({'endpoint':ep,'study':study,'cells':len(g),'discordant_cells':changed,'discordance_percent':100*changed/len(g),'redistribution_mass_cells':mass,'total_variation_percent':100*mass/len(g),'largest_absolute_label_count_change':int(delta.abs().max()),'label_with_largest_absolute_change':str(delta.abs().idxmax())})
pd.DataFrame(rows).to_csv(O/'predicted_count_summary.csv',index=False);pd.DataFrame(counts).to_csv(O/'all_predicted_label_counts.csv',index=False)
# Independent recovery and exact all-competitor condition, using saved native scores.
# Strict inequalities express a unique winner; tied cases are flagged, not conflated.
reference_ties=0;native_B_ties=0;max_res=0.;changed_after=0;condition_mismatches=0;cert_violations=0;near_boundary=0;queries=0;cell_queries=0;checks=[]
for row in q.itertuples():
 folder=R/'results/main/queries'/row.query_id
 with np.load(folder/'observed.npz',allow_pickle=False) as ar,np.load(folder/'predicted_B.npz',allow_pickle=False) as pr:
  A=ar['A'];B=ar['B'];shift=pr['shift'];corrected=B-shift;res=float(np.max(np.abs(corrected-A)));max_res=max(max_res,res);ch=int(np.sum(corrected.argmax(1)!=A.argmax(1)));changed_after+=ch
  winner=A.argmax(1);lead=A[np.arange(len(A)),winner,None]-A;relative=shift[None,:]-shift[winner,None];gaps=lead-relative;gaps[np.arange(len(A)),winner]=np.inf
  # Check all-class inequalities against the native ranking. Numerical boundary cases have a saved tolerance-band check.
  minimum=gaps.min(1);eps=2*(1e-6+1e-6*np.abs(pr['scores']).max(1));near=np.abs(minimum)<=eps;near_boundary+=int(near.sum());unique_stays=minimum>0;native_stays=B.argmax(1)==winner;condition_mismatches+=int(np.sum((unique_stays!=native_stays)&~near))
  sorted_A=np.sort(A,axis=1);sorted_B=np.sort(B,axis=1);margin=sorted_A[:,-1]-sorted_A[:,-2];reference_ties+=int((margin==0).sum());native_B_ties+=int((sorted_B[:,-1]==sorted_B[:,-2]).sum());safe=margin>np.ptp(shift);cert_violations+=int(np.sum(safe&~native_stays));queries+=1;cell_queries+=len(A)
  checks.append({'query_id':row.query_id,'cells':len(A),'recovery_max_absolute_residual':res,'recovered_label_mismatches':ch,'exact_condition_mismatches_outside_tolerance':int(np.sum((unique_stays!=native_stays)&~near)),'condition_boundary_band_cells':int(near.sum())})
pd.DataFrame(checks).to_csv(O/'all_query_recovery_checks.csv',index=False)
assert queries==3794 and changed_after==0 and condition_mismatches==0 and cert_violations==0
assert max_res<1e-6 and reference_ties==0 and native_B_ties==0
assert sum(row['discordant_cells'] for row in rows if row['endpoint']=='mechanism')==351 and sum(row['discordant_cells'] for row in rows if row['endpoint']=='composition')==166
# Stress-test proof statements on deterministic synthetic arrays, independent of saved inference helpers.
rng=np.random.default_rng(270026);n=5000;k=9;z=rng.normal(size=(n,k));shift=rng.normal(size=k);winner=z.argmax(1);gaps=z[np.arange(n),winner,None]-z-(shift[None,:]-shift[winner,None]);gaps[np.arange(n),winner]=np.inf
assert np.array_equal(gaps.min(1)>0,(z+shift).argmax(1)==winner)
safe=(np.sort(z,axis=1)[:,-1]-np.sort(z,axis=1)[:,-2])>np.ptp(shift);assert np.all((z+shift).argmax(1)[safe]==winner[safe])
u=rng.normal(size=(23,k));assert np.ptp(u.sum(0))<=np.ptp(u,axis=1).sum()+1e-12
# Non-necessity: large common competitor shift creates a loose spread bound while the winner stays.
z0=np.array([10.,9.,-100.]);d0=np.array([0.,0.,2.]);assert z0.argmax()==(z0+d0).argmax() and 1<=np.ptp(d0)
# Strict boundary: equality can create a tie and is not certified.
assert np.ptp([0.,1.])==1 and np.array_equal(np.array([1.,0.])+[0.,1.],[1.,1.])
# Equal class shifts preserve ranks but shift native sigmoid confidence.
assert np.array_equal(z.argmax(1),(z+2).argmax(1))
report={'scope':'post-analysis skeptical-review checks; unchanged main inference and primary endpoints','queries':queries,'reference_exact_tie_occurrences':reference_ties,'native_B_exact_tie_occurrences':native_B_ties,'all_predicted_counts_independently_checked_against_transition_tables':True,'cell_query_occurrences':cell_queries,'score_offset_recovery_max_absolute_residual':max_res,'recovered_label_mismatches':changed_after,'all_competitor_condition_mismatches_outside_tolerance_band':condition_mismatches,'all_competitor_boundary_band_cell_query_occurrences':near_boundary,'spread_certificate_violations':cert_violations,'synthetic_vectors_checked':n,'synthetic_cases':['all-competitor equivalence','spread sufficiency','sum-of-gene-range bound','non-necessity counterexample','tie boundary','common-offset rank invariance'],'primary_mechanism_changes':351,'primary_composition_changes':166,'count_summary':rows,'no_biological_accuracy_claim':True,'public_release_authorized':False}
(O/'proof_and_impact_checks.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
