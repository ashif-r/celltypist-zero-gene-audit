"""Rebuild in a new folder; never overwrite the reviewed project's outputs."""
from pathlib import Path
import argparse,sys,json,shutil,subprocess,os,datetime
R=Path(__file__).resolve().parents[1];pa=argparse.ArgumentParser();pa.add_argument('--output',type=Path,required=True);pa.add_argument('--reuse-downloads',type=Path);pa.add_argument('--prepare-only',action='store_true');args=pa.parse_args();O=args.output.expanduser().resolve()
if O.exists():raise RuntimeError('Output directory must not exist: '+str(O))
if O==R or R in O.parents and any(p in ['data','results','deliverables'] for p in O.relative_to(R).parts):raise ValueError('Choose a separate new run folder, e.g. runs/raw_rebuild')
O.mkdir(parents=True)
for folder in ['config','src','scripts','tests','examples','manifests']:
 for p in (R/folder).rglob('*'):
  if p.is_file() and '__pycache__' not in p.parts and p.suffix in ['.py','.json','.csv']:
   target=O/p.relative_to(R);target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,target)
for name in ['protocol.md','amendments.md','requirements.lock.txt','run_all.sh','REPRODUCE.md','RELEASE_CONTENTS.md','.gitignore']:
 shutil.copy2(R/name,O/name)
for folder in ['logs','results/gate','results/input_repair','data/prepared','deliverables']: (O/folder).mkdir(parents=True,exist_ok=True)
env=os.environ.copy();env.update(MPLCONFIGDIR=str(O/'.mpl'),NUMBA_CACHE_DIR=str(O/'.numba'),CELLTYPIST_FOLDER=str(O/'data/celltypist'),OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
def run(name,*arguments):
 print('Running',name,flush=True)
 with (O/'logs'/name).open('w') as out:subprocess.run([sys.executable,*map(str,arguments)],cwd=O,env=env,stdout=out,stderr=subprocess.STDOUT,check=True)
cache=['--cache',args.reuse_downloads.resolve()] if args.reuse_downloads else []
run('retrieve.log',O/'scripts/29_retrieve_inputs.py',*cache)
for number in ['01_metadata','02_checkpoint','03_pilot','04_input_audit','07_fetch_geo','08_scib_masks','09_rebuild_inputs','10_public_api_control','12_prepare_design']:
 run(number+'.log',O/'scripts'/f'{number}.py')
 if number=='04_input_audit':
  original=json.loads((R/'manifests/feature_source_files.json').read_text());rebuilt=json.loads((O/'manifests/feature_source_files.json').read_text())
  assert rebuilt==[{k:v for k,v in item.items() if k!='url'} for item in original], 'Rebuilt feature provenance differs'
  # The historical audit writer omits URL fields. Retain the verified original
  # retrieval metadata byte-for-byte, after verifying all regenerated fields.
  shutil.copy2(R/'manifests/feature_source_files.json',O/'manifests/feature_source_files.json')

run('secondary_prepare.log',O/'scripts/14_scib_secondary.py','--prepare')
# Frozen design/source/input hashes must still match; this is a later reproduction,
# not a new prospective freeze or registration. Do not rewrite the original record.
run('freeze_integrity.log','-c',"from src.freeze import require_freeze; from pathlib import Path; print(require_freeze(Path.cwd()))")
record={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'role':'later_raw_input_rebuild_not_new_preregistration','raw_input_rebuild_complete':True,'frozen_source_and_input_hashes_match':True,'mode':'prepare_only' if args.prepare_only else 'full'}
(O/'RAW_REBUILD_VERIFIED.json').write_text(json.dumps(record,indent=2)+'\n')
if args.prepare_only:print('Raw input rebuild and freeze hashes verified:',O);sys.exit(0)
for runname in ['main','reproduction']:
 run(runname+'.log',O/'scripts/11_main.py','--run',runname)
 run(runname+'_summary.log',O/'scripts/13_summarize.py','--run',runname)
run('figures.log',O/'scripts/16_figures_and_supplements.py','--run','main')
run('independent_audit.log',O/'scripts/17_verify_outputs.py','--run','main','--compare','reproduction')
run('exact_reproduction.log',O/'scripts/25_exact_reproduction.py')
run('secondary.log',O/'scripts/14_scib_secondary.py')
run('secondary_reproduction.log',O/'scripts/14_scib_secondary.py','--reproduction')
run('secondary_report.log',O/'scripts/19_secondary_report.py')
run('independent_input.log',O/'scripts/26_verify_original_inputs.py')
(O/'results/ANALYSIS_VERIFIED.json').write_text(json.dumps({'main_and_reproduction_complete':True,'independent_aggregation_and_algebra_audit':True,'secondary_reproduction_complete':True},indent=2)+'\n')
run('complete_supplements.log',O/'scripts/24_complete_supplements.py')
run('final_partition_figure.log',O/'scripts/27_final_partition_figure.py')
run('model_screen_export.log',O/'scripts/32_export_model_screen.py')
run('model_feature_overlap.log',O/'scripts/35_model_overlap_metadata.py')
print('Complete reproducible analysis:',O)
