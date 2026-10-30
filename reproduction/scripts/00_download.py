"""Download/check exact input artifacts from the saved manifests. No inference."""
from pathlib import Path
import json,urllib.request,hashlib
R=Path(__file__).resolve().parents[1];(R/'data').mkdir(exist_ok=True)
items=json.loads((R/'manifests/artifacts.json').read_text())['artifacts']
items+=json.loads((R/'manifests/feature_source_files.json').read_text())
for item in items:
 name=item.get('name',item.get('file'));target=R/'data'/name
 if not target.exists():
  print('Downloading',name,flush=True)
  with urllib.request.urlopen(item['url'],timeout=120) as response,target.open('wb') as out:
   while b:=response.read(8388608):out.write(b)
 h=hashlib.sha256()
 with target.open('rb') as stream:
  while b:=stream.read(8388608):h.update(b)
 if h.hexdigest()!=item['sha256']:raise ValueError('Checksum mismatch: '+name)
 print('Verified',name)
