"""Retrieve authors' notebooks as evidence only; never execute downloaded code."""
from pathlib import Path
import urllib.request,json,concurrent.futures,hashlib
r=Path(__file__).resolve().parents[1];base='https://raw.githubusercontent.com/theislab/scib-reproducibility/3afbffd3674726e5146797be21cf6bd7470a2c5f/notebooks/data_preprocessing/immune_cells/'
paths=['Bone_Marrow_Human/BM_Oetjen_'+str(i)+'.ipynb' for i in [2,3]]+['PBMCs_Human/Sun_'+str(i)+'_PBMCs.ipynb' for i in [2,3,4]]
def one(p):
 b=urllib.request.urlopen(base+p,timeout=40).read();x=json.loads(b);text='\n'.join(''.join(c.get('source',[])) for c in x['cells']);out=r/'private/source_checks'/(p.split('/')[-1]+'.source.txt');out.write_text(text)
 return {'path':p,'url':base+p,'sha256':hashlib.sha256(b).hexdigest(),'relevant_lines':[s for s in text.splitlines() if any(t in s for t in ['GSM','gene_file','file_base','http'])]}
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:res=list(pool.map(one,paths))
(r/'private/source_checks/additional_notebooks.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res,indent=2))
