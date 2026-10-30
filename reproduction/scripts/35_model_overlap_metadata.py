"""Report model-feature coverage using fixed schema metadata, without new inference."""
from pathlib import Path
import json
import pandas as pd,numpy as np
R=Path(__file__).resolve().parents[1];O=R/'results/main';genes=set(np.load(O/'tables/model_contributions.npz',allow_pickle=False)['genes']);features=pd.read_csv(R/'manifests/original_scoring_features.csv');v=pd.read_csv(O/'validation.csv');q=pd.read_csv(O/'queries.csv');rows=[]
for study,g in q.loc[q.context.ne('stress')].groupby('study'):
 schema=set(features.loc[features.study.eq(study),'symbol']);overlap=len(schema&genes);assert overlap>0
 for row in g.itertuples():
  check=v.loc[v.query_id.eq(row.query_id)].iloc[0];pr=np.load(O/'queries'/row.query_id/'predicted_B.npz');removed=sum(name in genes for name in pr['removed_genes']);assert removed==int(check.removed_model_genes)
  rows.append({'study':study,'context':row.context,'sample':row.sample,'query_id':row.query_id,'fixed_scoring_columns':len(schema),'checkpoint_model_genes':len(genes),'A_matched_model_genes':overlap,'model_genes_absent_or_identifier_excluded_in_A':len(genes)-overlap,'B_removed_model_genes':removed,'B_matched_model_genes':overlap-removed,'unmeasured_model_genes_zero_filled':False})
pd.DataFrame(rows).to_csv(O/'tables/model_feature_overlap.csv',index=False);print(pd.DataFrame(rows).to_string(index=False))
