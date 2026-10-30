"""Isolated native-API workers; each query retains its predeclared identity and order."""
from pathlib import Path
import json,hashlib
import numpy as np,anndata as ad
from celltypist.models import Model
from .inference import query
from .design import retained_mask

def initialize(root,study,prepared_path,reference_path,output_path,residual):
 global R,STUDY,ADATA,GENES,IDS,INDEX,X,A,MODEL,OUTPUT,RESIDUAL,PLAN
 R=Path(root);STUDY=study;ADATA=ad.read_h5ad(prepared_path);GENES=ADATA.var_names.to_numpy(dtype=str);IDS=ADATA.obs_names.to_numpy(dtype=str);INDEX={v:i for i,v in enumerate(IDS)};X=ADATA.X.tocsr()
 ref=np.load(reference_path);np.testing.assert_array_equal(IDS,ref['cell_ids']);A=ref['scores'];MODEL=Model.load(str(R/'data/Immune_All_Low.pkl'));OUTPUT=Path(output_path);RESIDUAL=residual
 PLAN={q['query_id']:q for q in json.loads((R/'manifests/query_plan.json').read_text())}

def run_unit(task):
 number,(context,members,size,seed)=task
 positions=np.array([INDEX[v] for v in members]);qid=f'{STUDY}_{number:05d}';dest=OUTPUT/'queries'/qid;frozen=PLAN[qid]
 assert frozen['cells']==members.tolist() and frozen['context']==context
 assert hashlib.sha256(retained_mask(X[positions]).tobytes()).hexdigest()==frozen['B_mask_sha256']
 validation=query(X[positions],GENES,members,MODEL,A[positions],dest,do_D=context!='stress')
 validation.update(query_id=qid,study=STUDY,context=context,cells=len(members),A_direct_max_residual=RESIDUAL)
 sample=str(ADATA.obs.loc[members[0],'batch']) if context=='sample' else ''
 return number,validation,{'query_id':qid,'study':STUDY,'context':context,'sample':sample,'size':size,'seed':seed,'cells':len(members)}
