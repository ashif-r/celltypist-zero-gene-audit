"""Audit main inputs and recover original feature identity; no classifier inference."""
from pathlib import Path
import sys,json,hashlib
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import h5py,numpy as np,pandas as pd
from src.gate_core import identifier_keep
R=Path(__file__).resolve().parents[1];O=R/'results/gate'
obs=pd.read_csv(O/'main_obs.csv.gz',index_col=0);genes=pd.read_csv(O/'main_var.csv.gz').iloc[:,0].astype(str).to_numpy();assert obs.index.is_unique
sources={'Oetjen_A':'Oetjen_A_genes.tsv.gz','Oetjen_P':'GSM3396176_genes_P.tsv.gz','Oetjen_U':'GSM3396184_genes_U.tsv.gz','Sun_sample1_CS':'Sun_CS_genes.tsv.gz','Sun_sample2_KC':'GSM3665017_KC_PBL_genes.tsv.gz','Sun_sample3_TB':'GSM3665018_TB_PBL_genes.tsv.gz','Sun_sample4_TC':'GSM3665019_TC_PBL_genes.tsv.gz'}
frames={};meta=[]
for batch,name in sources.items():
 p=R/'data'/name;f=pd.read_csv(p,sep='\t',header=None,names=['feature_id','symbol']);frames[batch]=f
 assert set(genes)<=set(f.symbol),'Unknown original feature identity'
 meta.append({'batch':batch,'file':name,'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'original_features':len(f),'ambiguous_symbols':sorted(set(f.symbol[~identifier_keep(f.symbol)]))})
rows=[];feature_rows=[];cellrows=[]
with h5py.File(R/'data/Immune_ALL_human.h5ad') as h:
 counts=h['layers/counts']
 for study in ['Oetjen','Sun']:
  batches=obs.loc[obs.study.eq(study),'batch'].unique();amb=set()
  for b in batches:amb.update(frames[b].symbol[~identifier_keep(frames[b].symbol)])
  keep=~np.isin(genes,list(amb));assert len(set(genes[keep]))==keep.sum()
  for batch in batches:
   f=frames[batch];unique=f.loc[identifier_keep(f.symbol)].set_index('symbol')
   for symbol in genes:
    excluded=symbol in amb
    feature_rows.append({'study':study,'batch':batch,'recorded_symbol':symbol,'original_feature_id':unique.loc[symbol,'feature_id'] if not excluded else '|'.join(f.loc[f.symbol.eq(symbol),'feature_id']),'excluded':excluded,'reason':'symbol_ambiguous_in_study_source_tables' if excluded else 'exact_unique_original_symbol_match'})
   indices=np.flatnonzero(obs.batch.eq(batch));finite=True;nonnegative=True;integer=True;empty=0;mass=0.;excludedmass=0.;minlib=float('inf');maxlib=0.
   for start in range(0,len(indices),256):
    idx=indices[start:start+256];v=counts[idx,:].astype(np.float64)
    finite &= bool(np.isfinite(v).all());nonnegative &= bool((v>=0).all());integer &= bool(np.equal(v,np.floor(v)).all())
    totals=v.sum(axis=1);empty+=int((totals==0).sum());mass+=totals.sum();excludedmass+=v[:,~keep].sum();minlib=min(minlib,float(totals.min()));maxlib=max(maxlib,float(totals.max()))
    for i,total in zip(idx,totals):cellrows.append({'cell_id':str(obs.index[i]),'study':study,'batch':batch,'available_feature_library_total':float(total),'eligible_nonempty':bool(total>0)})
   assert finite and nonnegative and integer and empty==0
   rows.append({'study':study,'batch':batch,'natural_sample':batch,'donor':'three donors documented by scIB source' if study=='Oetjen' else 'donor identity unresolved; natural acquisition sample documented','cells':len(indices),'representation':'verified original counts layer; available-feature normalization planned','feature_count':len(genes),'fixed_excluded_symbols':int((~keep).sum()),'scoring_features':int(keep.sum()),'measured_support':'all retained symbols occur in every original sample feature table; unique original identity after fixed exclusions','finite':finite,'nonnegative':nonnegative,'integer_counts':integer,'empty_libraries':empty,'min_library':minlib,'max_library':maxlib,'excluded_count_fraction':float(excludedmass/mass),'mechanism_eligible':True,'composition_eligible':True,'status':'input_checks_passed; main remains frozen pending review/public record'})
for study in ['10X','Freytag','Villani']:
 sub=obs.loc[obs.study.eq(study)]
 rows.append({'study':study,'batch':study,'cells':len(sub),'representation':'TPM' if study=='Villani' else 'UMI','mechanism_eligible':False,'composition_eligible':False,'status':'TPM input not accepted as counts' if study=='Villani' else 'original sample-specific identifier provenance not yet recovered; single query unit'})
pd.DataFrame(rows).to_csv(O/'sample_table.csv',index=False)
pd.DataFrame(feature_rows).to_csv(O/'main_original_feature_mapping.csv.gz',index=False)
pd.DataFrame(cellrows).to_csv(O/'main_cell_eligibility.csv.gz',index=False)
(R/'manifests/feature_source_files.json').write_text(json.dumps(meta,indent=2)+'\n')
summary={'source_mapping':'exact original symbol matches after excluding all ambiguous symbols; no suffix stripping','measured_support_verified_studies':['Oetjen','Sun'],'eligible_natural_samples':7,'eligible_cells':sum(int(row['cells']) for row in rows if row.get('mechanism_eligible')),'eligible_multisample_studies':2,'available_feature_count':len(genes),'count_checks_passed':True,'main_predictions_computed':False,'limits':['Reduced cross-study intersection inherited from benchmark; no direction assumed for effect size','Original full transcriptome counts not reconstructed','10X/Freytag identifier provenance not resolved; Villani TPM excluded','Only 10x natural samples qualify, so no Smart-seq2 composition claim']}
(O/'input_audit.json').write_text(json.dumps(summary,indent=2)+'\n');print(pd.DataFrame(rows).to_string(index=False));print(json.dumps(summary,indent=2))
