"""Foreground, small ESIM showcase campaign; no installation mutation."""
import hashlib
import importlib.metadata as md
import json
import platform
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
import radia
from cubit_mesh_export.check import check_consistency

root = Path(__file__).resolve().parent
package = Path(radia.__file__).parent
calc = package / 'panels/calc_inductance.py'
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
report = check_consistency(str(root / 'ih_bem_sample_p1.vol'), threshold=10.0)
(root / 'vol_check.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
assert report['passed'], report
for family, name in (('materials','workpiece'),('boundaries','sibc')):
    row=next(r for r in report[family] if r['name']==name)
    assert abs(row['error_pct']) < 1.0, row
record = dict(generated_at_utc=datetime.now(timezone.utc).isoformat(),
              hostname=platform.node(), python=sys.version,
              versions={n: md.version(n) for n in ('radia','ngsolve','numpy','scipy','cubit-mesh-export')},
              solver_path=str(calc), inputs={p.name:sha(p) for p in root.iterdir() if p.suffix in ('.vol','.step','.txt')},
              sources={}, cases=[])
for rel in ('panels/calc_inductance.py','esim_cell_problem.py','bem_sibc_solver.py'):
    p=package/rel
    record['sources'][rel]={'sha256':sha(p), 'source':p.read_text(encoding='utf-8')}
record['sources']['run_cases.py']={'sha256':sha(Path(__file__)), 'source':Path(__file__).read_text(encoding='utf-8')}
try:
    for freq in (10000,50000,100000):
        stem=f'f{freq}'
        cmd=[sys.executable,str(calc),'--coil-step',str(root/'ih_fem_kelvin_demo_coil.step'),
             '--coil-solver','peec','--vol',str(root/'ih_bem_sample_p1.vol'),
             '--wp-label','sibc','--sigma','2e6','--mu-r','100','--half-thickness','0.005',
             '--coil-sigma','5.8e7','--impedance-model','esim','--bh-file',str(root/'em_sample_bh.txt'),
             '--esim-max-iter','30','--esim-tol','1e-3','--esim-relax','0.5',
             '--esim-anderson-m','5','--h1-order','1','--wp-bem-backend','intree-dense',
             '--frequency',str(freq),'--current','100','--esim-per-panel',
             '--output',str(root/(stem+'.json'))]
        start=time.perf_counter()
        print('RUN',stem,flush=True)
        with (root/(stem+'.log')).open('w',encoding='utf-8') as log:
            proc=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,timeout=600)
        item=dict(frequency_hz=freq,command=cmd,returncode=proc.returncode,elapsed_s=time.perf_counter()-start)
        record['cases'].append(item)
        assert proc.returncode == 0, item
        d=json.loads((root/(stem+'.json')).read_text(encoding='utf-8'))
        assert d.get('esim_converged') and d.get('esim_per_panel_H_t'), d.keys()
        q=Path(d['qsurf_sol'])
        (root/(stem+'_qsurf.sol')).write_bytes(q.read_bytes())
        item.update(power_W=d['P_wp_W'],iterations=d['esim_iterations'],result_sha256=sha(root/(stem+'.json')))
        print(item,flush=True)
finally:
    (root/'compute_record.json').write_text(json.dumps(record,indent=2),encoding='utf-8')
