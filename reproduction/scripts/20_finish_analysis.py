"""Continue the verified live main run into aggregation, audit and clean reproduction."""
from pathlib import Path
import sys,time,subprocess,json,datetime,argparse
import psutil
R=Path(__file__).resolve().parents[1];p=argparse.ArgumentParser();p.add_argument('--main-pid',type=int,required=True);args=p.parse_args();main=psutil.Process(args.main_pid)
assert any(str(v).endswith('11_main.py') for v in main.cmdline()),'Expected the owned scientific main worker'
while not (R/'results/main/COMPLETED.json').exists():
 if (R/'results/main/FAILED.txt').exists():raise RuntimeError('Main analysis failed; inspect preserved traceback')
 if not main.is_running() or main.status()==psutil.STATUS_ZOMBIE:raise RuntimeError('Main worker ended without COMPLETED')
 time.sleep(2)
def run(name,*argv):
 print('Running',name,flush=True)
 with (R/'logs'/name).open('w') as out:subprocess.run([sys.executable,*map(str,argv)],cwd=R,stdout=out,stderr=subprocess.STDOUT,check=True)
run('main_summary.log',R/'scripts/13_summarize.py','--run','main')
run('main_figures.log',R/'scripts/16_figures_and_supplements.py','--run','main')
run('independent_verification.log',R/'scripts/17_verify_outputs.py','--run','main')
run('reproduction.log',R/'scripts/11_main.py','--run','reproduction')
run('reproduction_summary.log',R/'scripts/13_summarize.py','--run','reproduction')
run('reproduction_verification.log',R/'scripts/17_verify_outputs.py','--run','main','--compare','reproduction')
run('secondary_report.log',R/'scripts/19_secondary_report.py')
(R/'results/ANALYSIS_VERIFIED.json').write_text(json.dumps({'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'main_and_reproduction_complete':True,'independent_aggregation_and_algebra_audit':True,'secondary_reproduction_complete':True},indent=2)+'\n');print('Analysis and reproduction verified.',flush=True)
