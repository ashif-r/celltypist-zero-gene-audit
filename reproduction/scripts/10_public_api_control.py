"""PBMC pilot control: public annotate API versus direct Classifier path."""
from pathlib import Path
import sys,tarfile,json
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np,pandas as pd
from scipy import io
from anndata import AnnData
import celltypist
from celltypist.models import Model
from celltypist.classifier import Classifier
from src.gate_core import identifier_keep,normalize_counts
R=Path(__file__).resolve().parents[1]
with tarfile.open(R/'data/pbmc3k.tar.gz') as t:
 def f(s):return t.extractfile(next(m for m in t.getmembers() if m.name.endswith(s)))
 genes=pd.read_csv(f('/genes.tsv'),sep='\t',header=None).iloc[:,1].to_numpy();ids=pd.read_csv(f('/barcodes.tsv'),header=None)[0].astype(str).to_numpy();counts=io.mmread(f('/matrix.mtx')).T.tocsr()[:32]
x,_=normalize_counts(counts);keep=identifier_keep(genes);a=AnnData(x[:,keep],obs=pd.DataFrame(index=ids[:32]),var=pd.DataFrame(index=genes[keep]))
p=str(R/'data/Immune_All_Low.pkl')
public=celltypist.annotate(a.copy(),model=Model.load(p),majority_voting=False,mode='best match')
direct=Classifier(a.copy(),model=Model.load(p)).celltype(mode='best match')
for attr in ['decision_matrix','probability_matrix']:
 np.testing.assert_array_equal(getattr(public,attr).to_numpy(),getattr(direct,attr).to_numpy())
np.testing.assert_array_equal(public.predicted_labels.to_numpy(),direct.predicted_labels.to_numpy())
np.savez_compressed(R/'results/gate/annotate_control.npz',scores=public.decision_matrix.to_numpy(),probabilities=public.probability_matrix.to_numpy(),labels=public.predicted_labels.iloc[:,0].astype(str).to_numpy())
(R/'results/gate/annotate_control.json').write_text(json.dumps({'dataset':'PBMC3k first 32 archive-order cells','cells':32,'public_API':'celltypist.annotate','mode':'best match','majority_voting':False,'scores_probabilities_labels_exactly_equal':True,'main_inference':False},indent=2)+'\n');print('Public annotate and Classifier paths exactly agree for pilot control.')
