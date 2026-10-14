"""Heavy dense-scatter contracts; self-authored mesh and summation bounds."""
import hashlib
import importlib.util
import json
from pathlib import Path
import ngsolve as ng
import numpy as np
import pytest
from radia import _radia_pybind as native


def nodes_for(points,triangles):
    corners=points[triangles]
    return np.ascontiguousarray(np.concatenate((corners,
        (corners[:,[0,1,2]]+corners[:,[1,2,0]])*.5),axis=1))


def test_large_p1_parallel_scatter_contribution_bound():
    reference=json.loads((Path(__file__).parent/'data/dense_scatter_ring_reference.json').read_text(encoding='utf-8'))
    path=Path(__file__).parents[1]/'induction_heating/run_loop_work_ring.py'
    spec=importlib.util.spec_from_file_location('scatter_ring',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    _,points,triangles=module.ring_mesh(reference['level'])
    nodes=nodes_for(points,triangles)
    assert len(points)>5000 and len(triangles)==reference['faces']
    assert hashlib.sha256(points.tobytes()+triangles.tobytes()+nodes.tobytes()).hexdigest()==reference['geometry_sha256']
    for threads in (1,4):
        ng.SetNumThreads(threads)
        with ng.TaskManager():
            matrices=native._AssembleSLDL_Galerkin(points,triangles,nodes,
                reference['regular_degree'],reference['singular_order'],threads)
        for entry in reference['entries']:
            count=2*entry['terms']+4
            gamma=count*np.finfo(float).eps/(1-count*np.finfo(float).eps)
            for name,matrix in zip(('SL','DL'),matrices):
                assert abs(matrix[entry['i'],entry['j']]-entry[name])<=gamma*entry['abs_'+name]
        del matrices


@pytest.mark.parametrize('ndof',[10,5001])
def test_curved_p2_small_and_large_storage_paths(ndof):
    points=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
    triangles=np.array([[0,2,1],[0,1,3],[0,3,2],[1,2,3]],dtype=np.int64)
    nodes=nodes_for(points,triangles)
    # A consistent curved displacement for each shared edge.
    edge_ids={edge:4+i for i,edge in enumerate(((0,1),(0,2),(0,3),(1,2),(1,3),(2,3)))}
    dofs=np.empty((4,6),dtype=np.int64);dofs[:,:3]=triangles
    for a,tri in enumerate(triangles):
        for local,(i,j) in enumerate(((0,1),(1,2),(2,0))):
            edge=tuple(sorted((int(tri[i]),int(tri[j]))))
            dofs[a,3+local]=edge_ids[edge]
            nodes[a,3+local,2]+=.01
    tables=[]
    for threads in (1,4):
        ng.SetNumThreads(threads)
        with ng.TaskManager():
            matrices=native._AssembleSLDL_Galerkin_P2(points,triangles,nodes,dofs,ndof,7,6,threads)
        # Padded inactive DOFs force the >5000 storage path without expensive
        # irrelevant pair integration. P1 above uses >5000 active vertices.
        tables.append([matrix[:10,:10].copy() for matrix in matrices])
        for matrix in matrices:
            assert np.count_nonzero(matrix[10:])==0
            assert np.count_nonzero(matrix[:,10:])==0
        del matrices
    # 9 incident triangle pairs maximum per corner, plus deterministic local
    # quadrature. Frobenius scaling avoids component-relative cancellation.
    # This small P2 gate complements the P1 per-entry contribution certificate.
    gamma=64*np.finfo(float).eps/(1-64*np.finfo(float).eps)
    for one,multi in zip(*tables):
        assert np.linalg.norm(one-multi)<=gamma*np.linalg.norm(one)
