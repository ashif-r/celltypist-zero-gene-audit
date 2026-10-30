"""Retrieve original seven GEO samples; record checksums without running source code."""
from pathlib import Path
import json,urllib.request,hashlib,concurrent.futures,time
R=Path(__file__).resolve().parents[1]; D=R/'data/geo_original';D.mkdir(exist_ok=True)
specs=[('Oetjen_A','GSM3396161','A'),('Oetjen_P','GSM3396176','P'),('Oetjen_U','GSM3396184','U'),('Sun_sample1_CS','GSM3665016','CS'),('Sun_sample2_KC','GSM3665017','KC'),('Sun_sample3_TB','GSM3665018','TB'),('Sun_sample4_TC','GSM3665019','TC')]
jobs=[]
for batch,acc,code in specs:
 for kind in ['matrix','barcodes','genes']:
  name=f'{acc}_{kind}_{code}.'+('mtx.gz' if kind=='matrix' else 'tsv.gz') if batch.startswith('Oetjen') else f'{acc}_{code}_PBL_{kind}.'+('mtx.gz' if kind=='matrix' else 'tsv.gz')
  jobs.append({'batch':batch,'accession':acc,'kind':kind,'name':name,'url':f'https://ftp.ncbi.nlm.nih.gov/geo/samples/{acc[:-3]}nnn/{acc}/suppl/{name}'})
def fetch(item):
 p=D/item['name']
 if not p.exists():
  temp=p.with_suffix(p.suffix+'.part')
  for attempt in range(3):
   try:
    with urllib.request.urlopen(item['url'],timeout=120) as response,temp.open('wb') as out:
     while b:=response.read(1048576):out.write(b)
    temp.replace(p);break
   except Exception:
    if attempt==2:raise
    time.sleep(2)
 h=hashlib.sha256(p.read_bytes()).hexdigest();print(item['name'],p.stat().st_size,flush=True)
 return dict(item,sha256=h,bytes=p.stat().st_size)
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:items=list(pool.map(fetch,jobs))
(R/'manifests/geo_original.json').write_text(json.dumps({'artifacts':items},indent=2)+'\n')
