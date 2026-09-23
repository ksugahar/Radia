"""Inspect retained FEM refinement without mixing HDiv polynomial orders."""
import hashlib
import json
from pathlib import Path
import numpy as np

root = Path('S:/Radia/validation_artifacts')
paths = [root/'hdiv_medium_20260924/cached-three.json',
         root/'hdiv_fine_20260924/cached-three.json',
         root/'hdiv_fullfiner_20260924/recovered/cached-three.json']
data = [json.loads(p.read_text(encoding='utf-8-sig')) for p in paths]
for d in data:
    assert d['passed']
    assert d['observation_points_m'] == data[0]['observation_points_m']
    assert d['coil'] == data[0]['coil']
    assert d['gap_core_half_length_m'] == data[0]['gap_core_half_length_m']
points = np.asarray(data[0]['observation_points_m'])
mask = abs(points[:, 0]) <= data[0]['gap_core_half_length_m'] + 1e-12
rows = {}
for name in ['reduced_a', 'mixed_total_reduced_omega']:
    contracts = [d['engine_checkpoint_contracts'][name] for d in data]
    for key in ['mode', 'fem_order', 'implementation_sha256']:
        assert all(c[key] == contracts[0][key] for c in contracts), key
    b = [np.asarray(d['median_plane_projected_fields_T'][name])[mask] for d in data]
    denom = np.linalg.norm(b[-1])
    inc = [float(np.linalg.norm(b[i+1]-b[i])/denom) for i in range(2)]
    rows[name] = dict(core_increments=inc, contraction_ratio=inc[1]/inc[0],
                      primary_dofs=[d['engines'][name]['ndof'] for d in data])
out = dict(levels=['medium','fine','finer'], engines=rows,
           inputs=[dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in paths],
           interpretation='Mesh sensitivity only, common finer denominator. HDiv omitted because these campaigns change BDM1 to BDM2. No continuum error bound or matched-error speed claim.')
target = root/'hdiv_fullfiner_20260924/fem_refinement.json'
target.write_text(json.dumps(out,indent=2),encoding='utf-8')
print(json.dumps(out,indent=2))
