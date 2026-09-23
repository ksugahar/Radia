"""Read-only field comparison to numerical BDM2 reference, never exact truth."""
import hashlib,json,sys
from pathlib import Path
import numpy as np
reference_path=Path(sys.argv[1]); campaign_path=Path(sys.argv[2])
ref=json.loads(reference_path.read_text()); c=json.loads(campaign_path.read_text())
p=np.asarray(c['observation_points_m']); keys={tuple(np.round(x,12)):i for i,x in enumerate(p)}
indices=[keys[tuple(np.round(x*[1,1,-1],12))] for x in p]
mask=np.abs(p[:,0])<=c['gap_core_half_length_m']+1e-12
b=np.asarray(ref['rows'][0]['field_T'])
assert b.shape==p.shape
b=((b+b[indices]*[-1,-1,1])*.5)[mask]
differences={name:float(np.linalg.norm(np.asarray(field)[mask]-b)/np.linalg.norm(b)) for name,field in c['median_plane_projected_fields_T'].items()}
result=dict(status='Cross-run diagnostic, not a new three-engine pass',reference='BDM2 finer numerical field, not analytic truth',relative_core_vector_differences=differences,reference_sha256=hashlib.sha256(reference_path.read_bytes()).hexdigest(),campaign_sha256=hashlib.sha256(campaign_path.read_bytes()).hexdigest(),warning='No speed ratio. Observation ordering relies on the shared retained driver and checked input contract; diagnostic JSON lacks its own point array.')
Path(sys.argv[3]).write_text(json.dumps(result,indent=2)); print(json.dumps(result,indent=2))
