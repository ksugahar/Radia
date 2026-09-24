"""Postprocess completed order-three evidence, without an accuracy certificate."""
import hashlib,json
from pathlib import Path
import numpy as np

base=Path('S:/Radia/validation_artifacts')
paths={
    'coarse_p3':base/'hdiv_order3_20260924/recovered/cached-three.json',
    'coarse_p2':base/'hdiv_mixed_20260924/cached-three.json',
    'finer_p2':base/'hdiv_fullfiner_20260924/recovered/cached-three.json',
}
data={key:json.loads(path.read_text()) for key,path in paths.items()}
current=data['coarse_p3']
for other in data.values():
    assert other['passed']
    assert other['observation_points_m']==current['observation_points_m']
    assert other['coil']==current['coil']
    assert other['gap_core_half_length_m']==current['gap_core_half_length_m']
    assert other['comparison_contract']['cad_authority_sha256']==current['comparison_contract']['cad_authority_sha256']
points=np.asarray(current['observation_points_m'])
mask=abs(points[:,0])<=current['gap_core_half_length_m']+1e-12
result={'kind':'p-refinement sensitivity, not absolute error or matched-error speed',
        'inputs':{key:{'path':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()} for key,path in paths.items()},
        'comparisons':{}}
for name in ('reduced_a','mixed_total_reduced_omega'):
    a=np.asarray(current['median_plane_projected_fields_T'][name])[mask]
    entry={}
    for label in ('coarse_p2','finer_p2'):
        b=np.asarray(data[label]['median_plane_projected_fields_T'][name])[mask]
        entry[label]=float(np.linalg.norm(a-b)/np.linalg.norm(b))
    result['comparisons'][name]=entry
out=base/'hdiv_order3_20260924/order_sensitivity.json'
if out.exists():
    raise FileExistsError(out)
out.write_text(json.dumps(result,indent=2))
print(json.dumps(result['comparisons'],indent=2))
