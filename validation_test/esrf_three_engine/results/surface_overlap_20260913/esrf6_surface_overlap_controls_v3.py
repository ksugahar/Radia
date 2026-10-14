"""Synthetic full-mesh controls for the ESRF6 read-only audit."""
import hashlib
import json
from pathlib import Path
from collections import Counter
from itertools import combinations
import ngsolve as ng
from netgen.meshing import Mesh, MeshPoint, Pnt, Element3D, Element2D, FaceDescriptor
from esrf6_surface_overlap_audit_v3 import audit, controls

root=Path('C:/temp/esrf6-newton-4d85e72cc')

def make(name, points, tets, include_shared):
    mesh=Mesh(dim=3)
    for point in points:
        mesh.Add(MeshPoint(Pnt(*point)))
    mesh.SetMaterial(1,'fixture')
    mesh.Add(FaceDescriptor(surfnr=1,domin=1,domout=0,bc=1))
    mesh.SetBCName(0,'fixture_boundary')
    faces=Counter()
    for tet in tets:
        mesh.Add(Element3D(1,tet))
        faces.update(tuple(sorted(face)) for face in combinations(tet,3))
    for face,count in faces.items():
        if count==1 or include_shared:
            mesh.Add(Element2D(1,list(face)))
    path=root/(name+'_v3.vol')
    if path.exists():
        raise FileExistsError(path)
    mesh.Save(str(path))
    return path

normal=make('control_two_owner',
            [(0,0,0),(1,0,0),(0,1,0),(0,0,1),(0,0,-1)],
            [(1,2,3,4),(1,3,2,5)],True)
overlap=make('control_partial_overlap',
             [(0,0,0),(1,0,0),(0,1,0),(0,0,1),
              (.2,.1,0),(1.2,.1,0),(.2,1.1,0),(.2,.1,-1)],
             [(1,2,3,4),(5,7,6,8)],False)
ng.SetNumThreads(2)
with ng.TaskManager():
    result={'narrow_phase':controls(),'normal_two_owner':audit(normal,2),
            'partial_overlap':audit(overlap,2)}
assert result['normal_two_owner']['detected_parent_pairs']==0
assert result['normal_two_owner']['topology']['boundary_face_owner_histogram'][2]==1
assert result['partial_overlap']['detected_parent_pairs']>0
assert any(c['kind']=='coplanar_positive_area' for c in result['partial_overlap']['contacts'])
result['assertions_passed']=True
result['control_script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
result['audit_script_sha256']=hashlib.sha256((root/'esrf6_surface_overlap_audit_v3.py').read_bytes()).hexdigest()
with (root/'surface_overlap_controls_v3.json').open('x',encoding='utf-8') as out:
    json.dump(result,out,indent=2,allow_nan=False)
print('Full-mesh controls PASS',flush=True)
