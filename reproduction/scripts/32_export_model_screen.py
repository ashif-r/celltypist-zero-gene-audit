"""Export the model-only diagnostic using plain arrays, without object/pickle fields."""
from pathlib import Path
import json
import numpy as np,pandas as pd
R=Path(__file__).resolve().parents[1];T=R/'results/main/tables';T.mkdir(exist_ok=True)
# This is our trusted checkpoint-screen output, not an arbitrary uploaded archive.
old=np.load(R/'results/gate/model_only.npz',allow_pickle=True);a=old['a_kg'];genes=old['genes'].astype(str);classes=old['classes'].astype(str);assert a.shape==(98,6639)
p=T/'model_contributions.npz';np.savez_compressed(p,a_kg=a,genes=genes,classes=classes)
new=np.load(p,allow_pickle=False)
for key in new.files:assert new[key].dtype.kind!='O'
np.testing.assert_array_equal(new['a_kg'],a)
# Verify the per-gene diagnostic and the sum-of-ranges bound for all query masks.
rank=pd.read_csv(R/'results/gate/model_gene_screen.csv').set_index('gene').loc[genes]
np.testing.assert_allclose(rank.R_g.to_numpy(),np.ptp(a,axis=0),atol=1e-12,rtol=1e-12)
lookup={g:i for i,g in enumerate(genes)};checks=0;maxexcess=0.
for folder in (R/'results/main/queries').iterdir():
 if not folder.is_dir():continue
 pr=np.load(folder/'predicted_B.npz');ix=[lookup[g] for g in pr['removed_genes'] if g in lookup];bound=float(np.ptp(a[:,ix],axis=0).sum());spread=float(np.ptp(pr['shift']));maxexcess=max(maxexcess,spread-bound);assert spread<=bound+1e-10;checks+=1
report={'plain_array_export':True,'allow_pickle_false_loads_all_fields':True,'features':len(genes),'classes':len(classes),'gene_range_diagnostic_verified':True,'sum_gene_ranges_bound_verified_queries':checks,'positive_bound_violation':maxexcess}
(R/'results/main/model_screen_verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
