import json,sys,hashlib
from pathlib import Path
import numpy as np
paths=[Path(p) for p in sys.argv[1:4]]
data=[json.loads(p.read_text()) for p in paths]
campaign=json.loads(Path(sys.argv[4]).read_text())
points=np.asarray(campaign['observation_points_m'])
mask=np.abs(points[:,0])<=campaign['gap_core_half_length_m']+1e-12
keys={tuple(np.round(p,12)):i for i,p in enumerate(points)}
indices=[keys[tuple(np.round(p*[1,1,-1],12))] for p in points]
def project(b): return (b+b[indices]*[-1,-1,1])*0.5
selected=[next(r for r in d['rows'] if r['rule']==2) for d in data]
fields=[np.asarray(r['field_T']) for r in selected]
core=[project(b)[mask] for b in fields]
inc=[float(np.linalg.norm(core[i+1]-core[i])/np.linalg.norm(core[-1])) for i in range(2)]
out=dict(status='BDM2 diagnostic, not three-engine pass or accuracy certificate',core_increments=inc,ratio=inc[1]/inc[0],primary_dof=[r['diagnostics']['ndof'] for r in selected],runtime_s=[r['diagnostics']['runtime_s'] for r in selected],native_residual=[r['solve_scalars']['last_solve_final_relative_residual'] for r in selected],source_hashes=[hashlib.sha256(p.read_bytes()).hexdigest() for p in paths],normalization='Fine core vector norm for both increments; campaign axial projection and core mask')
Path(sys.argv[5]).write_text(json.dumps(out,indent=2))
print(json.dumps(out,indent=2))
