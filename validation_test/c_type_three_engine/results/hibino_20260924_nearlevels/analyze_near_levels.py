import json,sys
from pathlib import Path
import numpy as np
data=[json.loads(Path(p).read_text()) for p in sys.argv[1:4]]
campaign=json.loads(Path(sys.argv[4]).read_text())
points=np.asarray(campaign['observation_points_m'])
mask=np.abs(points[:,0])<=campaign['gap_core_half_length_m']+1e-12
keys={tuple(np.round(p,12)):i for i,p in enumerate(points)}
indices=[keys[tuple(np.round(p*[1,1,-1],12))] for p in points]
def project(b): return (b+b[indices]*[-1,-1,1])*0.5
rows=[]
for j in range(3):
    fields=[np.asarray(d['rows'][j]['field_T']) for d in data]
    core=[project(b)[mask] for b in fields]
    inc=[float(np.linalg.norm(core[i+1]-core[i])/np.linalg.norm(core[-1])) for i in range(2)]
    rows.append(dict(intorder=data[0]['rows'][j]['intorder'],core_increments=inc,ratio=inc[1]/inc[0]))
out=dict(status='diagnostic_not_accuracy_certificate',rows=rows,normalization='Fine norm for both increments; median-plane axial projection and campaign core mask')
Path(sys.argv[5]).write_text(json.dumps(out,indent=2))
print(json.dumps(out,indent=2))
