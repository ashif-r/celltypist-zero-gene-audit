"""Prepare an allowlisted prospective bundle. Does not publish or run inference."""
from pathlib import Path
import json,hashlib,shutil,zipfile,datetime,argparse
parser=argparse.ArgumentParser();parser.add_argument('--implementation-revision',action='store_true');args=parser.parse_args()
R=Path(__file__).resolve().parents[1];D=R/'deliverables/clean_release';D.mkdir(exist_ok=True)
if any(D.iterdir()):raise RuntimeError('Clean release already exists; archive/review before replacement')
files=[R/n for n in ['protocol.md','amendments.md','requirements.lock.txt','run_all.sh','REPRODUCE.md','RELEASE_CONTENTS.md','.gitignore']]
for folder in ['config','src','scripts','tests','examples','manifests']:
 files.extend(p for p in (R/folder).rglob('*') if p.is_file() and '__pycache__' not in p.parts and not p.name.startswith('freeze_') and p.suffix in ['.py','.json','.csv'])
files=sorted(set(files));entries=[]
for p in files:
 rel=p.relative_to(R)
 assert not set(rel.parts)&{'private','data','logs','.venv'}
 q=D/rel;q.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,q);entries.append({'path':str(rel),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
m={'status':'implementation_revision_after_partial_main_outputs' if args.implementation_revision else 'prospective_main_analysis_bundle_after_disclosed_PBMC_pilot','protocol_version':'0.4-draft','created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'main_annotation_outcomes_seen':args.implementation_revision,'files':entries}
(R/'manifests/freeze_inputs.json').write_text(json.dumps(m,indent=2)+'\n');shutil.copy2(R/'manifests/freeze_inputs.json',D/'freeze_inputs.json')
z=R/('deliverables/implementation_revision_v0.4.1_private.zip' if args.implementation_revision else 'deliverables/protocol_freeze_v0.4_private.zip')
with zipfile.ZipFile(z,'w',zipfile.ZIP_DEFLATED) as out:
 for p in sorted(D.rglob('*')):
  if p.is_file():out.write(p,str(p.relative_to(D)))
with zipfile.ZipFile(z) as out:
 assert out.testzip() is None
 for name in out.namelist():assert not any(v in Path(name).parts for v in ['private','data','logs','.venv'])
info={'created_utc':m['created_utc'],'bundle_path':str(z.relative_to(R)),'bundle_sha256':hashlib.sha256(z.read_bytes()).hexdigest(),'bundle_bytes':z.stat().st_size,'manifest_path':'manifests/freeze_inputs.json','manifest_sha256':hashlib.sha256((R/'manifests/freeze_inputs.json').read_bytes()).hexdigest(),'file_count':len(entries),'public_authorization':False,'public_record':None}
(R/'deliverables/clean_release_summary.json').write_text(json.dumps(info,indent=2)+'\n');print(json.dumps(info,indent=2))
