"""Independent aggregation from persisted all-class scores. No model inference."""
from pathlib import Path
import sys,argparse,json
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np,pandas as pd
from scipy.special import expit
from src.design import paired_summary,stability_certificate,material_effect
R=Path(__file__).resolve().parents[1];pa=argparse.ArgumentParser();pa.add_argument('--run',default='main');args=pa.parse_args();O=R/'results'/args.run
if not (O/'COMPLETED.json').exists():raise RuntimeError('Cannot summarize incomplete main run')
T=O/'tables';T.mkdir(exist_ok=True);qs=pd.read_csv(O/'queries.csv');primary=[];stress=[];class_rows=[];shift_rows=[];percell=[]
def contrast(ids,left,right,classes):
 li=left.argmax(1);ri=right.argmax(1);la=classes[li];lb=classes[ri]
 s=paired_summary(ids,dict(zip(ids,la)),dict(zip(ids,lb)))
 # A separate path uses integer class indices and explicit Python counts.
 assert s['changed']==sum(int(a!=b) for a,b in zip(li.tolist(),ri.tolist()))
 assert s['valid']==len(ids) and s['complete'] and np.isfinite(left).all() and np.isfinite(right).all()
 ca=expit(left.max(1));cb=expit(right.max(1));delta=cb-ca
 s.update(confidence_change_mean=float(delta.mean()),confidence_change_median=float(np.median(delta)),both_high_different=int(((li!=ri)&(ca>=.9)&(cb>=.9)).sum()))
 for q in [0,.25,.5,.75,.95,.99,1]:s['confidence_change_q'+str(q)]=float(np.quantile(delta,q))
 for t in [.5,.9,.95]:
  s[f'cross_up_{t}']=int(((ca<t)&(cb>=t)).sum());s[f'cross_down_{t}']=int(((ca>=t)&(cb<t)).sum())
 return s
for study in qs.study.unique():
 poolrow=qs.loc[qs.study.eq(study)&qs.context.eq('pool')].iloc[0];pool=np.load(O/'queries'/poolrow.query_id/'observed.npz');poolpred=np.load(O/'queries'/poolrow.query_id/'predicted_B.npz');pool_idx={v:i for i,v in enumerate(pool['cell_ids'])}
 for row in qs.loc[qs.study.eq(study)].itertuples():
  d=np.load(O/'queries'/row.query_id/'observed.npz');pred=np.load(O/'queries'/row.query_id/'predicted_B.npz');A=d['A'];B=d['B'];ids=d['cell_ids'];classes=d['classes'];ix=np.array([pool_idx[v] for v in ids]);BP=pool['B'][ix]
  assert np.array_equal(d['A'],pool['A'][ix])
  assert np.allclose(pred['scores'],B,atol=1e-6,rtol=1e-6)
  assert np.array_equal(pred['scores'].argmax(1),B.argmax(1)),'Destination mismatch needs explicit near-tie investigation'
  mech=contrast(ids,A,B,classes);comp=contrast(ids,BP,B,classes)
  csafe=stability_certificate(BP,pred['shift']-poolpred['shift']);assert not np.any(csafe&(BP.argmax(1)!=B.argmax(1)))
  if row.context=='sample':
   for endpoint,s in [('mechanism',mech),('composition',comp)]:primary.append(dict(query_id=row.query_id,study=study,sample=row.sample,endpoint=endpoint,certified_stable=int(d['certified'].sum()) if endpoint=='mechanism' else int(csafe.sum()),**s))
   aa=classes[A.argmax(1)];bb=classes[B.argmax(1)];pp=classes[BP.argmax(1)]
   for label in np.unique(aa):
    mask=aa==label;n=int(mask.sum())
    if n>=50:class_rows.append({'study':study,'sample':row.sample,'A_label':label,'A_label_cells':n,'mechanism_changes':int((aa[mask]!=bb[mask]).sum()),'composition_changes':int((pp[mask]!=bb[mask]).sum())})
   for i,cid in enumerate(ids):percell.append({'cell_id':cid,'study':study,'sample':row.sample,'A_label':aa[i],'B_sample_label':bb[i],'B_pool_label':pp[i],'A_confidence':float(expit(A[i].max())),'B_sample_confidence':float(expit(B[i].max())),'B_pool_confidence':float(expit(BP[i].max())),'mechanism_certified':bool(d['certified'][i]),'composition_certified':bool(csafe[i])})
  if row.context!='stress':
   for label,value in zip(classes,pred['shift']):shift_rows.append({'study':study,'query_id':row.query_id,'context':row.context,'sample':row.sample,'class':label,'score_shift':float(value)})
  else:
   stress.append({'study':study,'query_id':row.query_id,'size':int(row.size),'seed':int(row.seed),'actual_cells':len(ids),'mechanism_changed':mech['changed'],'composition_changed':comp['changed'],'mechanism_certified':int(d['certified'].sum()),'composition_certified':int(csafe.sum())})
p=pd.DataFrame(primary);p.to_csv(T/'primary.csv',index=False);pd.DataFrame(percell).to_csv(T/'per_cell_primary.csv.gz',index=False);pd.DataFrame(class_rows).to_csv(T/'per_A_label.csv',index=False);pd.DataFrame(shift_rows).to_csv(T/'class_shifts.csv',index=False);st=pd.DataFrame(stress);st.to_csv(T/'stress_queries.csv',index=False)
by=p.groupby(['study','endpoint'],as_index=False)[['eligible','valid','changed','failed','certified_stable']].sum();by['fraction']=by.changed/by.valid;by.to_csv(T/'study_summary.csv',index=False)
agg=st.groupby(['study','size','seed'],as_index=False)[['actual_cells','mechanism_changed','composition_changed','mechanism_certified','composition_certified']].sum();agg['mechanism_fraction']=agg.mechanism_changed/agg.actual_cells;agg['composition_fraction']=agg.composition_changed/agg.actual_cells;agg.to_csv(T/'stress_summary.csv',index=False)
for study,g in agg.groupby('study'):assert (g.actual_cells==sum(p.loc[p.study.eq(study)&p.endpoint.eq('mechanism'),'eligible'])).all()
result={'run':args.run,'eligible_cells':len(percell),'primary_comparisons':len(p),'all_primary_pairs_complete':bool(p.complete.all()),'material_effect_rule_met':material_effect(p.to_dict('records')),'material_threshold':.01,'equal_study_means':by.groupby('endpoint').fraction.mean().to_dict(),'publication_decision':'report_valid_results_regardless_of_effect_size_subject_to_release_approval','queries':len(qs)}
(O/'results.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
