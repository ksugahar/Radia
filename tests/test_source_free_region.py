"""Independent elementary source-exclusion and topology contracts."""
from itertools import permutations
import importlib.util
from pathlib import Path
import numpy as np
import pytest
import ngsolve as ng
from netgen.meshing import Mesh,MeshPoint,Pnt,Element3D

spec=importlib.util.spec_from_file_location('source_free_region',
    Path(__file__).parents[1]/'src/radia/source_free_region.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
certify=module.certify_current_free_total_region


def blocks(hole=False):
    mesh=Mesh(dim=3);mesh.SetMaterial(1,'core');points={}
    for i in range(3 if hole else 1):
        for j in range(3 if hole else 1):
            if hole and (i,j)==(1,1):continue
            for perm in permutations(range(3)):
                p=np.array([i,j,0]);cell=[]
                for axis in (None,)+perm:
                    if axis is not None:p=p.copy();p[axis]+=1
                    key=tuple(p)
                    if key not in points:points[key]=mesh.Add(MeshPoint(Pnt(*map(float,key))))
                    cell.append(points[key])
                mesh.Add(Element3D(1,cell))
    result=ng.Mesh(mesh)
    with ng.TaskManager():
        result.Curve(1)
    return result


def test_disjoint_declared_source_and_contract():
    result=certify(blocks(),('core',),[[[2.,2.,2.],[3.,3.,3.]]])
    assert result['passed'] and result['betti_1_GF2']==0
    assert result['connected_components']==1
    assert result['source_support_completeness'].startswith('caller-owned')


@pytest.mark.parametrize('box',[
    [[[4.,4.,4.],[5.,5.,5.]]],
    [[[1.2,1.2,-1.],[1.8,1.8,2.]]],
])
def test_ring_rejects_source_free_claim(box):
    with pytest.raises(ValueError,match='first homology'):
        certify(blocks(True),('core',),box)


@pytest.mark.parametrize('box',[
    [[[.2,.2,.2],[.3,.3,.3]]],
    [[[1.,.2,.2],[2.,.3,.3]]],
    [[[-1.,.2,.2],[2.,.3,.3]]],
])
def test_containment_touching_and_crossing_reject(box):
    with pytest.raises(ValueError,match='intersects or touches'):
        certify(blocks(),('core',),box)


@pytest.mark.parametrize('boxes',[[],[[[0,0,0],[0,1,1]]],[[[0,0,0],[1,1,float('nan')]]]])
def test_incomplete_or_invalid_support_reject(boxes):
    with pytest.raises(ValueError,match='support boxes'):
        certify(blocks(),('core',),boxes)


def test_separating_axis_beyond_coordinate_boxes():
    tet=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
    assert module._box_separated(tet,np.array([.8,.8,.8]),np.array([.9,.9,.9]),1e-10)
    assert not module._box_separated(tet,np.array([.1,.1,.1]),np.array([.2,.2,.2]),1e-10)


def test_curved_geometry_is_not_certified_as_straight():
    mesh=blocks()
    with ng.TaskManager():
        mesh.Curve(2)
    with pytest.raises(ValueError,match='straight tetrahedra'):
        certify(mesh,('core',),[[[2.,2.,2.],[3.,3.,3.]]])


def test_vertex_touching_components_are_not_a_volume_manifold():
    mesh=Mesh(dim=3);mesh.SetMaterial(1,'core')
    points=[mesh.Add(MeshPoint(Pnt(*xyz))) for xyz in
        [(0,0,0),(1,0,0),(0,1,0),(0,0,1),(-1,0,0),(0,-1,0),(0,0,-1)]]
    mesh.Add(Element3D(1,[points[i] for i in (0,1,2,3)]))
    mesh.Add(Element3D(1,[points[i] for i in (0,4,5,6)]))
    mesh=ng.Mesh(mesh)
    with ng.TaskManager():
        mesh.Curve(1)
    with pytest.raises(ValueError,match='vertex link'):
        certify(mesh,('core',),[[[2.,2.,2.],[3.,3.,3.]]])
