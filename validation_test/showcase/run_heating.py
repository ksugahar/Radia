"""Rerun the declared ESIM frequency cases with installed Radia."""
from pathlib import Path
import json
import subprocess
import sys
import platform
import hashlib
from datetime import datetime, timezone
import radia
import argparse
import math
from recording import runtime_metadata, scratch_directory

ROOT=Path(__file__).resolve().parents[2]
if __name__=='__main__':
    source=ROOT/'validation_test/esim_spatial'
    parser = argparse.ArgumentParser()
    parser.add_argument('--workdir', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args_cli = parser.parse_args()
    output = scratch_directory(args_cli.workdir)
    output.mkdir(parents=True,exist_ok=True)
    original=json.loads((source/'compute_record.json').read_text())
    solver=Path(radia.__file__).parent/'panels/calc_inductance.py'
    records=[]
    for case in original['cases']:
        args=case['command'][2:]
        target=output/f"f{case['frequency_hz']}.json"
        for i,arg in enumerate(args):
            if arg.startswith('validation_test/esim_spatial/inputs/'):
                args[i]=str(ROOT/arg)
        args[args.index('--output')+1]=str(target)
        log=output/f"f{case['frequency_hz']}.log"
        with log.open('w',encoding='utf-8') as stream:
            completed=subprocess.run([sys.executable,str(solver),*args],cwd=ROOT,stdout=stream,stderr=subprocess.STDOUT)
        if completed.returncode: raise RuntimeError(f'{case["frequency_hz"]}: solver failed; inspect {log}')
        result=json.loads(target.read_text())
        if result.get('status')!='ok' or not result.get('esim_converged') or not math.isfinite(result['P_wp_W']) or result['P_wp_W']<=0:
            raise AssertionError('ESIM convergence/power check failed')
        records.append(dict(frequency_hz=case['frequency_hz'],power_W=result['P_wp_W'],converged=result['esim_converged']))
        print(records[-1],flush=True)
    record=dict(schema='radia.selected_showcase.v1',case='heating',generated_at_utc=datetime.now(timezone.utc).isoformat(),platform_class=platform.system(),**runtime_metadata(__file__, threads=None),solver_sha256=hashlib.sha256(solver.read_bytes()).hexdigest(),result=dict(cases=records,checks={'positive_power_and_convergence':all(r['converged'] and math.isfinite(r['power_W']) and r['power_W']>0 for r in records)},scope='Weak-coupled ESIM integrated workpiece power; spatial-loss accuracy and transient temperature are separate checks.'))
    args_cli.output.parent.mkdir(parents=True, exist_ok=True)
    args_cli.output.write_text(json.dumps(record,indent=2,allow_nan=False),encoding='utf-8')
