from pathlib import Path
import sys,json,hashlib,warnings,platform,importlib.metadata,datetime,pickle
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np,pandas as pd
R=Path(__file__).resolve().parents[1]; O=R/'results/gate'; O.mkdir(parents=True,exist_ok=True)
with warnings.catch_warnings(record=True) as ws:
    raw=pickle.loads((R/'data/Immune_All_Low.pkl').read_bytes())
clf=raw['Model'];scaler=raw['Scaler_'];genes=np.asarray(clf.features)
assert clf.coef_.shape[1]==len(genes)==len(scaler.scale_)
assert len(set(genes))==len(genes)
assert np.isfinite(clf.coef_).all() and np.isfinite(scaler.mean_).all() and (scaler.scale_>0).all()
a=clf.coef_*scaler.mean_/scaler.scale_
np.savez_compressed(O/'model_only.npz',a_kg=a,genes=genes,classes=clf.classes_)
t=pd.DataFrame({'gene':genes,'R_g':np.ptp(a,axis=0),'M_g':np.abs(a).max(axis=0)})
t.sort_values(['R_g','gene'],ascending=[False,True]).to_csv(O/'model_gene_screen.csv',index=False)
settings={'timestamp_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'official_url':'https://celltypist.cog.sanger.ac.uk/models/Pan_Immune_CellTypist/v2/Immune_All_Low.pkl','sha256':hashlib.sha256((R/'data/Immune_All_Low.pkl').read_bytes()).hexdigest(),'description':raw['description'],'with_mean':bool(scaler.with_mean),'with_std':bool(scaler.with_std),'classes':len(clf.classes_),'features':len(genes),'python':platform.python_version(),'platform':platform.platform(),'warnings':[str(w.message) for w in ws],'packages':{p:importlib.metadata.version(p) for p in ['celltypist','numpy','scipy','scikit-learn','scanpy','anndata']},'class_order_sha256':hashlib.sha256('\n'.join(clf.classes_).encode()).hexdigest(),'feature_order_sha256':hashlib.sha256('\n'.join(genes).encode()).hexdigest(),'centred_route_supported':bool(scaler.with_mean)}
(O/'checkpoint.json').write_text(json.dumps(settings,indent=2,default=str)+'\n');print(json.dumps(settings,indent=2,default=str),flush=True)
