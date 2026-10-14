"""Independent projection checks and thin positive-area overlap control."""
import hashlib
import json
from pathlib import Path
import numpy as np
import shapely
from shapely.geometry import Polygon
from esrf6_surface_overlap_audit_v3 import triangle_contact

root=Path('C:/temp/esrf6-newton-4d85e72cc')
old=json.loads((root/'fem_overlap_n4.json').read_text())
triangles=np.array(old['audit']['contacts'][0]['triangle_coordinates_m'])
a,b=triangles
normal=np.cross(a[1]-a[0],a[2]-a[0])
normal/=np.linalg.norm(normal)

def projections(left,right):
    rows=[]
    for axis in range(3):
        pa,pb=Polygon(np.delete(left,axis,axis=1)),Polygon(np.delete(right,axis,axis=1))
        area=shapely.intersection(pa,pb,grid_size=1e-12).area
        usable=abs(normal[axis])>.1
        rows.append({'dropped_axis':axis,'normal_component':float(normal[axis]),
                     'usable_for_area':bool(usable),'projected_area_m2':float(area),
                     'physical_area_m2':float(area/abs(normal[axis])) if usable else None})
    return rows

positive=b+np.array([1e-6,0.,0.])
contact=triangle_contact(a,positive,1e-10,1e-16)
negative=triangle_contact(a,b,1e-10,1e-16)
assert negative is None
assert contact and contact['kind']=='coplanar_positive_area'
assert contact['overlap_area_m2']>100*1e-16
negative_axes=projections(a,b)
positive_axes=projections(a,positive)
areas=[x['physical_area_m2'] for x in positive_axes if x['usable_for_area']]
assert max(areas)/min(areas)-1 < 1e-4
assert all(x['physical_area_m2']<=1e-16 for x in negative_axes if x['usable_for_area'])
report={'assertions_passed':True,'length_tolerance_m':1e-10,'area_tolerance_m2':1e-16,
        'overlay_grid_m':1e-12,'positive_translation_m':[1e-6,0.,0.],
        'negative_contact':negative,'positive_contact':contact,
        'negative_three_projections':negative_axes,'positive_three_projections':positive_axes,
        'positive_projection_relative_spread':max(areas)/min(areas)-1,
        'original_triangle_coordinates_m':triangles.tolist(),
        'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'audit_script_sha256':hashlib.sha256((root/'esrf6_surface_overlap_audit_v3.py').read_bytes()).hexdigest(),
        'note':'Projection along a zero normal component collapses area and is explicitly unusable.'}
with (root/'surface_overlap_precision_controls.json').open('x',encoding='utf-8') as out:
    json.dump(report,out,indent=2,allow_nan=False)
print(json.dumps(report),flush=True)
