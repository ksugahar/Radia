"""Native P1 entry ownership, orientation, cache and lifecycle contracts.

Numerical independence from the pre-refactor kernel is checked by
run_hacapk_entry_identity.py's separate baseline/compare workflow.
"""
import gc

import ngsolve as ng
import numpy as np
import pytest

from radia import _radia_pybind as native


def geometry(reverse=False, curved=False):
    tetra = np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]])
    faces = np.array([[0, 2, 1], [0, 1, 3], [1, 2, 3], [2, 0, 3]], dtype=np.int64)
    # Disconnected copies exercise all three regular-distance quadrature bins.
    vertices = np.concatenate([tetra+[distance, 0., 0.] for distance in (0., 3., 7., 15.)])
    triangles = np.concatenate([faces+4*k for k in range(4)])
    for k in range(len(triangles)):
        triangles[k] = np.roll(triangles[k], k % 3)
    if reverse:
        triangles = triangles[:, [0, 2, 1]].copy()
    corners = vertices[triangles]
    nodes = np.concatenate([corners, .5*(corners+np.roll(corners, -1, axis=1))], axis=1)
    if curved:
        nodes[:, 3:, 2] += .01
    return vertices, triangles, nodes


@pytest.mark.parametrize('reverse', [False, True])
@pytest.mark.parametrize('curved', [False, True])
def test_owned_entries_and_cache_independence(reverse, curved):
    vertices, triangles, nodes = geometry(reverse, curved)
    ng.SetNumThreads(1)
    with ng.TaskManager():
        expected = native._AssembleSLDL_Galerkin(vertices, triangles, nodes, 7, 6, 1)
        sl, dl, provider = native._CreateP1HACApKGeometry(vertices, triangles, nodes, 7, 6, 0)
        _, _, cached = native._CreateP1HACApKGeometry(vertices, triangles, nodes, 7, 6, 4096)
        # The managers/providers own snapshots, not borrowed NumPy buffers.
        vertices[:] = np.nan; triangles[:] = -1; nodes[:] = np.nan
        del vertices, triangles, nodes
        gc.collect()
        for i in range(16):
            for j in range(16):
                value = provider.Entry(i, j)
                assert value == cached.Entry(i, j)
                m = value[4]
                gamma = (2*m+4)*np.finfo(float).eps/(1-(2*m+4)*np.finfo(float).eps)
                assert abs(value[0]-expected[0][i, j]) <= gamma*value[2]
                assert abs(value[1]-expected[1][i, j]) <= gamma*value[3]
        assert provider.GetStats()['cache_bytes'] == 0
        assert 0 < cached.GetStats()['cache_bytes'] <= 4096
        for handle, matrix in zip((sl, dl), expected):
            with pytest.raises(RuntimeError, match='not been built'):
                handle.MatVec(np.ones(16))
            for _ in range(2):
                assert handle.BuildHMatrix(aca_eps=1e-10, leaf_size=64)
                for transpose in (False, True):
                    result = handle.MatVec(np.arange(16, dtype=float), transpose=transpose)
                    reference = (matrix.T if transpose else matrix) @ np.arange(16)
                    assert np.linalg.norm(result-reference) <= 1e-12*np.linalg.norm(reference)
            assert handle.GetStats()['source_dense_bytes'] == 0
            with pytest.raises(RuntimeError, match='no dense'):
                handle.ReleaseDenseEntries()


def test_bad_geometry_and_entry_indices_fail_loudly():
    vertices, triangles, nodes = geometry()
    with pytest.raises(TypeError):
        native._CreateP1HACApKGeometry(vertices, triangles, nodes)
    with ng.TaskManager():
        _, _, provider = native._CreateP1HACApKGeometry(vertices, triangles, nodes, 7, 6, 0)
        for i, j in ((-1, 0), (0, -1), (16, 0), (0, 16)):
            with pytest.raises(IndexError):
                provider.Entry(i, j)
        invalid = triangles.copy(); invalid[0, 0] = len(vertices)
        with pytest.raises(ValueError, match='outside'):
            native._CreateP1HACApKGeometry(vertices, invalid, nodes, 7, 6)
        invalid = nodes.copy(); invalid[0, 0, 0] += 1
        with pytest.raises(ValueError, match='corners'):
            native._CreateP1HACApKGeometry(vertices, triangles, invalid, 7, 6)
        with pytest.raises(ValueError, match='requires coordinates'):
            native._CreateP1HACApKGeometry(vertices.ravel(), triangles, nodes, 7, 6)


def test_solver_forwards_quadrature_without_dense_source(monkeypatch):
    import importlib.util
    from pathlib import Path
    from scipy.sparse import issparse
    from radia.bem_sibc_solver import ScalarBIESIBCSolver
    path = Path(__file__).resolve().parents[1] / 'induction_heating/run_loop_work_ring.py'
    spec = importlib.util.spec_from_file_location('on_demand_contract_ring', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    mesh, _, _ = module.ring_mesh(16)
    factory = native._CreateP1HACApKGeometry
    seen = []
    def capture(coords, triangles, nodes, regular, singular, cache):
        seen.append((regular, singular, cache))
        return factory(coords, triangles, nodes, regular, singular, cache)
    def forbid_dense(*args, **kwargs):
        raise AssertionError('On-demand route called the dense assembler')
    monkeypatch.setattr(native, '_CreateP1HACApKGeometry', capture)
    monkeypatch.setattr(native, '_AssembleSLDL_Galerkin', forbid_dense)
    with ng.TaskManager():
        solver = ScalarBIESIBCSolver(mesh, order=1, use_intree_bem=True,
            use_intree_hacapk=True, assemble_dense=False, intree_geom_order=1,
            intree_regular_quad_degree=5, intree_singular_n_q=4,
            hacapk_entry_cache_bytes=0, hacapk_construction="on-demand")
    assert seen == [(5, 4, 0)]
    assert solver.SL is None and solver.DL is None
    assert issparse(solver.M) and issparse(solver.K)
    assert solver._body_construction_route == 'p1-entry-on-demand'
    assert solver._entry_provider.GetStats()['cache_bytes'] == 0
    assert solver._SL_hacapk.GetStats()['source_dense_bytes'] == 0


def test_small_default_keeps_dense_entry_construction(monkeypatch):
    import importlib.util
    from pathlib import Path
    from radia.bem_sibc_solver import ScalarBIESIBCSolver
    path = Path(__file__).resolve().parents[1] / 'induction_heating/run_loop_work_ring.py'
    spec = importlib.util.spec_from_file_location('default_construction_ring', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    mesh, _, _ = module.ring_mesh(16)
    def forbid_on_demand(*args, **kwargs):
        raise AssertionError('Small default silently switched to the slower on-demand path')
    monkeypatch.setattr(native, '_CreateP1HACApKGeometry', forbid_on_demand)
    with ng.TaskManager():
        solver = ScalarBIESIBCSolver(mesh, order=1, use_intree_bem=True,
            use_intree_hacapk=True, assemble_dense=False, intree_geom_order=1,
            intree_regular_quad_degree=7, intree_singular_n_q=6)
    assert solver._body_construction_route == 'dense-entry'
    assert solver.SL is None and solver.DL is None
    assert not hasattr(solver, '_entry_provider')
    assert solver._SL_hacapk.GetStats()['construction_route'] == 'dense-entry'
    for bad in (-1, 1.5, True):
        with pytest.raises(ValueError, match='nonnegative integer'):
            ScalarBIESIBCSolver(mesh, order=1, use_intree_bem=True,
                use_intree_hacapk=True, assemble_dense=False, hacapk_entry_cache_bytes=bad)
    with pytest.raises(ValueError, match='requires non-dense in-tree P1'):
        ScalarBIESIBCSolver(mesh, assemble_dense=True, hacapk_construction='on-demand')


@pytest.mark.parametrize('panel',[False,True])
def test_explicit_genus0_on_demand_matches_retained_dense_body(panel):
    import netgen.meshing as nm
    from radia.bem_sibc_solver import ScalarBIESIBCSolver
    from radia.surface_impedance import PanelSurfaceImpedance
    points=.03*np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
    triangles=((0,2,1),(0,1,3),(0,3,2),(1,2,3))
    raw=nm.Mesh(dim=3)
    face=raw.Add(nm.FaceDescriptor(surfnr=1,domin=1,bc=1))
    vertices=[raw.Add(nm.MeshPoint(nm.Pnt(*point))) for point in points]
    for triangle in triangles:
        raw.Add(nm.Element2D(face,[vertices[i] for i in triangle]))
    mesh=ng.Mesh(raw)
    omega=2*np.pi*5e4
    z=(1+1j)*np.sqrt(omega*4e-7*np.pi/(2*5.8e7))
    impedance=PanelSurfaceImpedance(z*np.array([.8,1.,1.4,2.])) if panel else z
    outputs=[]
    with ng.TaskManager():
        for compressed in (False,True):
            solver=ScalarBIESIBCSolver(mesh,order=1,assemble_dense=not compressed,
                use_intree_bem=True,use_intree_hacapk=True,intree_geom_order=1,
                intree_regular_quad_degree=7,intree_singular_n_q=6,
                hacapk_construction='on-demand' if compressed else 'auto')
            outputs.append(solver.solve_hacapk(-ng.z,impedance,omega,
                tol=1e-12,maxiter=500,restart=100))
    old,new=outputs
    for key in ('H_t_rms','P_density','area'):
        assert abs(new[key]-old[key])<=1e-10*abs(old[key])
    assert np.linalg.norm(new['phi_vec']-old['phi_vec'])<=1e-10*np.linalg.norm(old['phi_vec'])
    assert new['linear_residual_rel']<=1e-6
