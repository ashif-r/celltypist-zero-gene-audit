"""Final figure styling, including the planned zero-restoration control."""
from pathlib import Path
import json
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parents[1];O=R/'results/main';T=O/'tables';F=O/'figures'
assert json.loads((T/'validation_summary.json').read_text())['max_C_score_residual']==0.0
p=pd.read_csv(T/'primary.csv');st=pd.read_csv(T/'stress_summary.csv');validation=pd.read_csv(O/'validation.csv');queries=pd.read_csv(O/'queries.csv');plt.rcParams.update({'font.family':'DejaVu Sans','font.size':12,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':300})
colors={'Oetjen':'#176b87','Sun':'#b34a35'}
def save(fig,name):
 fig.savefig(F/(name+'.pdf'),bbox_inches='tight');fig.savefig(F/(name+'.png'),dpi=240,bbox_inches='tight');plt.close(fig)
fig,axes=plt.subplots(1,2,figsize=(10.2,4.2),layout='constrained')
for ax,endpoint,title in zip(axes,['mechanism','composition'],['A versus independently filtered sample B','Sample B versus within-study pooled B']):
 q=p.loc[p.endpoint.eq(endpoint)].sort_values(['study','sample']);labels=[v.replace('Oetjen_','Oetjen ').replace('Sun_sample','Sun ') for v in q['sample']]
 ax.bar(np.arange(len(q)),100*q.fraction,color=[colors[s] for s in q.study]);ax.set_xticks(np.arange(len(q)),labels,rotation=45,ha='right');ax.set_ylabel('Cells with different labels (%)');ax.set_title(title,fontsize=13);mx=float((100*q.fraction).max());ax.set_ylim(0,max(1,mx*1.30))
 for i,row in enumerate(q.itertuples()):ax.text(i,100*row.fraction+max(.03,mx*.025)+(i%2)*max(.06,mx*.04),f'{row.changed}/{row.valid}',ha='center',fontsize=11)
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
 ax.set_xticks(np.log2([32,128,512]),['32','128','512']);ax.set(xlabel='Target partition size (cells)',ylabel='Cells with different labels (%)',title=title);ax.set_ylim(bottom=0);ax.axhline(0,color='#53616a',ls='--',lw=1.5,label='Restored C versus A: 0%');ax.legend(frameon=False,fontsize=11.5,loc='lower left')
save(fig,'figure3_partition_sensitivity')
print('Final figures use larger labels and include the fixed-schema control.')
