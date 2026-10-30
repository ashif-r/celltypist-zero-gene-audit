"""Exhaustively decode saved predictions, confidence and composition diagnostics."""
from pathlib import Path
import json
import numpy as np,pandas as pd,anndata as ad
from scipy.special import expit
R=Path(__file__).resolve().parents[1];O=R/'results/main';T=O/'tables';assert (O/'independent_verification.json').exists()
q=pd.read_csv(O/'queries.csv');summaries=[];shiftrows=[];label_rows=[];normalization=[]
for study,g in q.groupby('study'):
 poolrow=g.loc[g.context.eq('pool')].iloc[0];pool=np.load(O/'queries'/poolrow.query_id/'observed.npz');pp=np.load(O/'queries'/poolrow.query_id/'predicted_B.npz');lookup={v:i for i,v in enumerate(pool['cell_ids'])}
 for row in g.itertuples():
  folder=O/'queries'/row.query_id;d=np.load(folder/'observed.npz');pr=np.load(folder/'predicted_B.npz');ix=np.array([lookup[v] for v in d['cell_ids']]);A=d['A'];B=d['B'];BP=pool['B'][ix];delta=pr['shift']-pp['shift'];comp_pred=BP+delta
  assert np.allclose(comp_pred,B,atol=1e-6,rtol=1e-6)
  for label,value in zip(d['classes'],pr['shift']):shiftrows.append({'query_id':row.query_id,'study':study,'context':row.context,'class':label,'score_shift':float(value)})
  for endpoint,left,pred in [('mechanism',A,pr['scores']),('composition',BP,comp_pred)]:
   lc=expit(left.max(1));rc=expit(B.max(1));pc=expit(pred.max(1));diff=rc-lc;li=left.argmax(1);ri=B.argmax(1);pi=pred.argmax(1)
   near=(np.sort(pred,axis=1)[:,-1]-np.sort(pred,axis=1)[:,-2])<=2*(1e-6+1e-6*np.abs(pred).max(1));mismatch=pi!=ri
   assert not (mismatch&~near).any()
   out={'query_id':row.query_id,'study':study,'context':row.context,'sample':row.sample,'endpoint':endpoint,'cells':len(B),'changed':int((li!=ri).sum()),'predicted_changed':int((li!=pi).sum()),'destination_mismatches':int(mismatch.sum()),'max_score_prediction_residual':float(np.abs(pred-B).max()),'max_confidence_prediction_residual':float(np.abs(pc-rc).max()),'mean_confidence_change':float(diff.mean()),'median_confidence_change':float(np.median(diff)),'both_high_different':int(((li!=ri)&(lc>=.9)&(rc>=.9)).sum())}
   for z in [0,.25,.5,.75,.95,.99,1]:out['confidence_change_q'+str(z)]=float(np.quantile(diff,z))
   for threshold in [.5,.9,.95]:
    up=(lc<threshold)&(rc>=threshold);down=(lc>=threshold)&(rc<threshold);decision_mismatch=(pc>=threshold)!=(rc>=threshold)
    assert not (decision_mismatch&(np.abs(pc-threshold)>1e-6)).any()
    out[f'cross_up_{threshold}']=int(up.sum());out[f'cross_down_{threshold}']=int(down.sum());out[f'predicted_threshold_mismatches_{threshold}']=int(decision_mismatch.sum())
   summaries.append(out)
  if row.context=='sample':
   labels=d['classes'][A.argmax(1)]
   for label in np.unique(labels):
    mask=labels==label;n=int(mask.sum())
    if n>=50:
     km=int((A.argmax(1)[mask]!=B.argmax(1)[mask]).sum());kc=int((BP.argmax(1)[mask]!=B.argmax(1)[mask]).sum());label_rows.append({'study':study,'sample':row.sample,'A_label':label,'A_label_cells':n,'mechanism_changes':km,'mechanism_fraction':km/n,'composition_changes':kc,'composition_fraction':kc/n})
 # Reconcile the inverse-log row-sum deficit created only by fixed identifier exclusions.
 a=ad.read_h5ad(R/f'data/prepared/{study}_original.h5ad');x=a.X.copy();x.data=np.expm1(x.data);scoring_sum=np.asarray(x.sum(1)).ravel();totals=a.obs.original_library_total.to_numpy();expected=10000*np.asarray(a.layers['counts'].sum(1)).ravel()/totals
 np.testing.assert_allclose(scoring_sum,expected,atol=1e-8,rtol=1e-12);assert np.all(totals>0)
 for i,cid in enumerate(a.obs_names):normalization.append({'cell_id':cid,'study':study,'sample':a.obs.batch.iloc[i],'original_library_total':float(totals[i]),'scoring_inverse_log_sum':float(scoring_sum[i]),'explained_deficit_from_excluded_rows':float(10000-scoring_sum[i])})
s=pd.DataFrame(summaries);s.to_csv(T/'all_query_confidence_and_composition.csv.gz',index=False);pd.DataFrame(shiftrows).to_csv(T/'all_query_class_shifts.csv.gz',index=False);pd.DataFrame(label_rows).to_csv(T/'per_A_label_rates.csv',index=False);pd.DataFrame(normalization).to_csv(T/'input_normalization.csv.gz',index=False)
report={'queries':len(q),'contrasts':len(s),'class_shift_rows':len(shiftrows),'all_class_composition_predictions_checked':True,'destination_mismatches':int(s.destination_mismatches.sum()),'max_composition_prediction_residual':float(s.loc[s.endpoint.eq('composition'),'max_score_prediction_residual'].max()),'max_confidence_prediction_residual':float(s.max_confidence_prediction_residual.max()),'threshold_decision_mismatches':int(s.filter(like='predicted_threshold_mismatches').sum().sum()),'normalization_rows_verified':len(normalization),'method_note':'Labels, sigmoids and thresholds decoded after inference from score predictions that were saved before each native B call; composition is additionally checked from pooled B plus the sample/partition-minus-pool shift.'}
(O/'complete_supplement_verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
