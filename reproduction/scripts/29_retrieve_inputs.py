"""Retrieve every required public file using preserved URLs and expected hashes.

Uses the GEO manifest for original-gene URLs as well as artifact URLs. The
historical audit writer omits URL fields on regeneration; never alter the frozen manifest.
"""
from pathlib import Path
import argparse,json,urllib.request,hashlib,shutil
R=Path(__file__).resolve().parents[1];pa=argparse.ArgumentParser();pa.add_argument('--cache',type=Path);args=pa.parse_args();D=R/'data';D.mkdir(exist_ok=True);(D/'geo_original').mkdir(exist_ok=True)
geo=json.loads((R/'manifests/geo_original.json').read_text())['artifacts'];genes={v['batch']:v for v in geo if v['kind']=='genes'}
jobs=[(v,D/v['name']) for v in json.loads((R/'manifests/artifacts.json').read_text())['artifacts']]
for v in json.loads((R/'manifests/feature_source_files.json').read_text()):
 source=genes[v['batch']];assert source['sha256']==v['sha256'];jobs.append((dict(source,name=v['file']),D/v['file']))
jobs.extend((v,D/'geo_original'/v['name']) for v in geo)
def digest(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  while b:=f.read(8388608):h.update(b)
 return h.hexdigest()
for item,path in jobs:
 if not path.exists():
  candidates=[]
  if args.cache:
   candidates=[args.cache/path.relative_to(D),args.cache/'geo_original'/item['name']]
  cached=next((p for p in candidates if p.is_file() and digest(p)==item['sha256']),None)
  if cached:shutil.copy2(cached,path)
  else:
   temp=path.with_suffix(path.suffix+'.part')
   with urllib.request.urlopen(item['url'],timeout=120) as response,temp.open('wb') as out:
    while b:=response.read(8388608):out.write(b)
   if digest(temp)!=item['sha256']:raise ValueError('Downloaded checksum mismatch: '+str(path))
   temp.replace(path)
 if digest(path)!=item['sha256']:raise ValueError('Checksum mismatch: '+str(path))
 print('Verified',path.relative_to(R),flush=True)
