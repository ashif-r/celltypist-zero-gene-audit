"""Render prespecified summaries and descriptive supplements from saved outputs only."""
from pathlib import Path
import argparse,json
import numpy as np,pandas as pd
from scipy.special import expit
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parents[1];pa=argparse.ArgumentParser();pa.add_argument('--run',default='main');args=pa.parse_args();O=R/'results'/args.run;T=O/'tables';F=O/'figures';F.mkdir(exist_ok=True)
p=pd.read_csv(T/'primary.csv');st=pd.read_csv(T/'stress_summary.csv');validation=pd.read_csv(O/'validation.csv');queries=pd.read_csv(O/'queries.csv');plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':300})
colors={'Oetjen':'#176b87','Sun':'#b34a35'}
def save(fig,name):
 fig.savefig(F/(name+'.pdf'),bbox_inches='tight');fig.savefig(F/(name+'.png'),dpi=240,bbox_inches='tight');plt.close(fig)
fig,axes=plt.subplots(1,2,figsize=(10.2,4.2),layout='constrained')
for ax,endpoint,title in zip(axes,['mechanism','composition'],['A versus independently filtered sample B','Sample B versus within-study pooled B']):
 q=p.loc[p.endpoint.eq(endpoint)].sort_values(['study','sample']);labels=[v.replace('Oetjen_','Oetjen ').replace('Sun_sample','Sun ') for v in q['sample']]
 ax.bar(np.arange(len(q)),100*q.fraction,color=[colors[s] for s in q.study]);ax.set_xticks(np.arange(len(q)),labels,rotation=45,ha='right');ax.set_ylabel('Cells with different labels (%)');ax.set_title(title,fontsize=10);mx=float((100*q.fraction).max());ax.set_ylim(0,max(1,mx*1.30))
 for i,row in enumerate(q.itertuples()):ax.text(i,100*row.fraction+max(.03,mx*.025),f'{row.changed}/{row.valid}',ha='center',fontsize=7)
save(fig,'figure1_natural_samples')
# A deterministic display sample; full-query residuals are represented separately.
expected=[];observed=[]
for q in queries.loc[queries.context.eq('sample')].itertuples():
 d=np.load(O/'queries'/q.query_id/'observed.npz');pr=np.load(O/'queries'/q.query_id/'predicted_B.npz');n=min(32,len(d['cell_ids']));expected.extend(pr['scores'][:n].ravel());observed.extend(d['B'][:n].ravel())
fig,axes=plt.subplots(1,2,figsize=(10.2,4.0),layout='constrained');axes[0].scatter(expected,observed,s=2,alpha=.2,color='#176b87',rasterized=True);lo=min(expected);hi=max(expected);axes[0].plot([lo,hi],[lo,hi],color='black',lw=.8);axes[0].set(xlabel='Analytically predicted B score',ylabel='Native API B score',title='First 32 retained cells per natural sample')
for study,g in validation.groupby('study'):axes[1].scatter(g.removed_model_genes,g.B_max_score_residual,s=5,alpha=.35,label=study,color=colors[study],rasterized=True)
axes[1].set(xlabel='Removed model genes per query',ylabel='Maximum absolute score residual',title='All planned original-input queries');axes[1].ticklabel_format(axis='y',style='sci',scilimits=(0,0));axes[1].legend(frameon=False);save(fig,'figure2_mechanistic_validation')
fig,axes=plt.subplots(1,2,figsize=(10.2,4),layout='constrained')
for ax,field,title in zip(axes,['mechanism_fraction','composition_fraction'],['A versus partition B','Partition B versus study-pooled B']):
 for study,g in st.groupby('study'):
  for seed,z in g.groupby('seed'):z=z.sort_values('size');ax.plot(np.log2(z['size']),z[field]*100,color=colors[study],alpha=.2,lw=1)
  mean=g.groupby('size')[field].mean();ax.plot(np.log2(mean.index),mean.to_numpy()*100,'o-',color=colors[study],label=study,lw=2)
 ax.set_xticks(np.log2([32,128,512]),['32','128','512']);ax.set(xlabel='Target partition size (cells)',ylabel='Cells with different labels (%)',title=title);ax.set_ylim(bottom=0);ax.legend(frameon=False)
save(fig,'figure3_partition_sensitivity')
# Full transition counts; no biological lineage inference is assigned.
cells=pd.read_csv(T/'per_cell_primary.csv.gz')
for endpoint,left,right in [('mechanism','A_label','B_sample_label'),('composition','B_pool_label','B_sample_label')]:
 z=cells.groupby(['study','sample',left,right],dropna=False).size().reset_index(name='cells');z.to_csv(T/(endpoint+'_transitions.csv'),index=False)
# Secondary D outcomes for every requested natural sample/pool context.
drows=[]
for q in queries.loc[queries.context.ne('stress')].itertuples():
 folder=O/'queries'/q.query_id;dp=folder/'secondary_D.npz'
 if not dp.exists():
  drows.append({'query_id':q.query_id,'study':q.study,'context':q.context,'sample':q.sample,'status':'input_contract_failure; see validation.csv'});continue
 d=np.load(dp);a=np.load(folder/'observed.npz');A=a['A'];D=d['scores'];ca=expit(A.max(1));cd=expit(D.max(1))
 drows.append({'query_id':q.query_id,'study':q.study,'context':q.context,'sample':q.sample,'status':'passed','cells':len(A),'changed':int((A.argmax(1)!=D.argmax(1)).sum()),'fraction':float((A.argmax(1)!=D.argmax(1)).mean()),'mean_confidence_change':float((cd-ca).mean()),'removed_genes':len(d['removed_genes'])})
pd.DataFrame(drows).to_csv(T/'secondary_D.csv',index=False)
qc={'all_primary_pairs_complete':bool(p.complete.all()),'primary_comparisons':len(p),'max_B_score_residual':float(validation.B_max_score_residual.max()),'max_C_score_residual':float(validation.C_max_score_residual.max()),'B_destination_mismatches':int(validation.B_destination_mismatches.sum()),'near_tie_B_count_across_queries':int(validation.B_near_ties.sum()),'queries':len(queries),'note':'Query-level counts across repeated partitions are not unique-cell population counts.'}
(T/'validation_summary.json').write_text(json.dumps(qc,indent=2)+'\n')
print('Generated three main figures and complete transition/D supplements.')
