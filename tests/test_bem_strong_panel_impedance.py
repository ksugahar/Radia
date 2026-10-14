"""Constant reduction, local work and value-keyed preparation in strong coupling."""
import importlib.util
from pathlib import Path
import numpy as np
import pytest
import ngsolve as ng
import netgen.meshing as nm
from netgen.occ import WorkPlane, Axes, Axis, OCCGeometry, Pnt, X, Y, Z, Revolve
from radia.bem_coupled_solver import CoupledBEMSolver
from radia.peec_coupled_bem_solver import CoupledPEECBEMSolver
from radia.surface_impedance import PanelSurfaceImpedance


@pytest.fixture(scope='module')
def setup_cases():
    spec = importlib.util.spec_from_file_location('strong_panel_ring',
        Path(__file__).resolve().parents[1]/'validation_test/induction_heating/run_loop_work_ring.py')
    ring = importlib.util.module_from_spec(spec);spec.loader.exec_module(ring)
    torus, _, _ = ring.ring_mesh(16)
    points = .015*np.array([[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]],float)
    tri = np.array([[0,2,4],[2,1,4],[1,3,4],[3,0,4],[2,0,5],[1,2,5],[3,1,5],[0,3,5]])
    for _ in range(2):
        edges={};pp=list(points);tt=[]
        def mid(a,b):
            key=tuple(sorted((int(a),int(b))))
            if key not in edges:
                q=points[a]+points[b];q=.015*q/np.linalg.norm(q);edges[key]=len(pp);pp.append(q)
            return edges[key]
        for a,b,c in tri:
            ab,bc,ca=mid(a,b),mid(b,c),mid(c,a)
            tt.extend([(a,ab,ca),(ab,b,bc),(ca,bc,c),(ab,bc,ca)])
        points=np.array(pp);tri=np.array(tt)
    sphere=nm.Mesh(dim=3);fd=sphere.Add(nm.FaceDescriptor(surfnr=1,domin=1,bc=1))
    numbers=[sphere.Add(nm.MeshPoint(nm.Pnt(*p))) for p in points]
    for t in tri:sphere.Add(nm.Element2D(fd,[numbers[i] for i in t]))
    sphere=ng.Mesh(sphere)
    profile=WorkPlane(Axes(Pnt(.06,0,0),n=Y,h=X)).Circle(.005).Face()
    solid=Revolve(profile,Axis(Pnt(0,0,0),Z),355)
    for face in solid.faces:
        face.name='body' if face.mass>5*np.pi*.005**2 else ('source' if face.center[1]>=-1e-9 else 'sink')
    solid.name='coil'
    from radia.panels.surface_mesh_extract import _extract_surface_mesh_filtered
    coil=_extract_surface_mesh_filtered(ng.Mesh(OCCGeometry(solid).GenerateMesh(
        mp=nm.MeshingParameters(maxh=.03,curvaturesafety=.7,segmentsperedge=1))))
    angle=np.arange(96)*2*np.pi/96;paths=[]
    for radius,height in [(.06,.005),(.065,-.005)]:
        p=np.c_[radius*np.cos(angle),radius*np.sin(angle),np.full(96,height)]
        paths.append(list(zip(p,np.roll(p,-1,axis=0))))
    return [sphere,torus],coil,paths


def make_solver(cases, model, genus):
    bodies,coil,paths=cases
    if model=='surface':return CoupledBEMSolver(coil,bodies[genus], wp_loop_work_backend='fmm' if genus else 'dense')
    return CoupledPEECBEMSolver(paths,np.diag([.002,.003]),
                               np.array([[3e-7,1e-7],[1e-7,3.5e-7]]),bodies[genus],
                               wp_loop_work_backend='fmm' if genus else 'dense')


def solve(solver,z,genus):
    return solver.solve(z,2*np.pi*50000,max_iter=40,tol=1e-9,relax=.8,loop_dof=bool(genus))


def assert_surface_field_close(actual, expected, areas, *, label):
    """Compare complex vector fields in the surface L2 norm, at the same 1e-10.

    Componentwise relative errors are ill-defined at symmetry zeros. The
    area-weighted norm measures the whole physical field (including phase),
    is invariant to units, and needs no absolute tolerance in A/m.
    """
    actual, expected, areas = map(np.asarray, (actual, expected, areas))
    assert actual.shape == expected.shape == (len(areas), 3), label
    assert np.all(areas > 0) and np.all(np.isfinite(areas)), label
    assert np.all(np.isfinite(actual)) and np.all(np.isfinite(expected)), label
    norm = lambda field: np.sqrt(areas @ np.sum(abs(field)**2, axis=1))
    error, reference = norm(actual-expected), norm(expected)
    assert error <= 1e-10*reference, (
        f'{label}: surface L2 error {error:.6g} exceeds 1e-10 * {reference:.6g}')


@pytest.mark.parametrize('scale', [1e-6, 1., 1e6])
def test_surface_field_comparison_is_scale_invariant(scale):
    areas = np.array([1., 3.])
    reference = scale*np.array([[0., 1j, 2.], [2j, -1., 0.]])
    rounded = reference.copy()
    rounded[0, 0] += scale*2e-12  # a symmetry-zero component
    assert_surface_field_close(rounded, reference, areas, label='roundoff')
    changed = reference.copy()
    changed[0, 0] += scale*1e-8
    with pytest.raises(AssertionError, match='surface L2 error'):
        assert_surface_field_close(changed, reference, areas, label='changed field')
    # A zero reference is not an excuse for an arbitrary absolute error floor.
    with pytest.raises(AssertionError, match='surface L2 error'):
        assert_surface_field_close(scale*2e-12*np.ones_like(reference),
                                   np.zeros_like(reference), areas, label='zero')


def compare(a,b):
    for key in ('P_total','H_t_rms','L_total','Delta_L','R_total','Delta_R',
                'body_reaction_power_W','coil_loss_change_W','port_power_W'):
        np.testing.assert_allclose(a[key],b[key],rtol=1e-10,atol=1e-18,err_msg=key)
    np.testing.assert_array_equal(a['wp_a'], b['wp_a'])
    assert_surface_field_close(a['wp_J_re']+1j*a['wp_J_im'],
                               b['wp_J_re']+1j*b['wp_J_im'], b['wp_a'], label='wp_J')
    assert_surface_field_close(a['wp_H_t_tri'], b['wp_H_t_tri'],
                               b['wp_a'], label='wp_H_t_tri')
    for key in ('wp_q_tri','body_emf'):
        np.testing.assert_allclose(a[key],b[key],rtol=1e-10,atol=1e-12,err_msg=key)
    if 'wp_loop_alpha' in a:
        np.testing.assert_allclose(a['wp_loop_alpha'],b['wp_loop_alpha'],rtol=1e-10,atol=1e-12)


@pytest.mark.parametrize('model',['surface','filament'])
@pytest.mark.parametrize('genus',[0,1])
def test_uniform_and_rewritten_panel_strong_work(setup_cases,model,genus):
    omega=2*np.pi*50000;z=(1+1j)*np.sqrt(omega*4e-7*np.pi/(2*5.8e7))
    with ng.TaskManager():
        solver=make_solver(setup_cases,model,genus)
        scalar=solve(solver,z,genus)
        uniform=solve(solver,PanelSurfaceImpedance(np.full(len(solver.wp_tris),z)),genus)
        compare(scalar,uniform)
        np.testing.assert_array_equal(uniform['Z_s_per_panel'],np.full(len(solver.wp_tris),z))
        centers=solver.wp_nodes[solver.wp_tris].mean(axis=1)
        values=np.where(centers[:,2]>=0,3*z,.25*z)
        nonuniform=solve(solver,PanelSurfaceImpedance(values),genus)
        assert nonuniform['body_power_balance_relative_error']<.1
        assert nonuniform['coupling_residual']<=1e-9 and nonuniform['body_residual']<=1e-6
        np.testing.assert_allclose(nonuniform['wp_q_tri'],
            .5*values.real*np.sum(abs(nonuniform['wp_H_t_tri'])**2,axis=1),rtol=1e-12)
        assert nonuniform['P_total']==pytest.approx(solver._phi_poisson._areas@nonuniform['wp_q_tri'],rel=1e-12)
        key=('_prepared_loop_system' if genus else '_prepared_gauge_factor')
        first=getattr(solver.wp_solver,key)
        magnetic = getattr(solver.wp_solver, '_prepared_loop_magnetic_work', None)
        again=solve(solver,PanelSurfaceImpedance(values.copy()),genus)
        assert getattr(solver.wp_solver,key) is first
        compare(nonuniform,again)
        rewritten=solve(solver,PanelSurfaceImpedance(1.1*values),genus)
        assert getattr(solver.wp_solver,key) is not first
        if genus:
            assert solver.wp_solver._prepared_loop_magnetic_work is magnetic
        assert rewritten['body_impedance_value_hash']!=nonuniform['body_impedance_value_hash']
        assert abs(rewritten['P_total']/nonuniform['P_total']-1)>1e-4
        fresh=solve(make_solver(setup_cases,model,genus),PanelSurfaceImpedance(1.1*values),genus)
        compare(rewritten,fresh)


def test_panel_arrays_and_unbuilt_hacapk_handles_fail_loudly(setup_cases):
    with ng.TaskManager():
        solver=make_solver(setup_cases,'filament',0)
        with pytest.raises(ValueError,match='tagged'):
            solve(solver,np.full(len(solver.wp_tris),1e-4+1e-4j),0)
        with pytest.raises(ValueError,match='one value'):
            solve(solver,PanelSurfaceImpedance(np.full(len(solver.wp_tris)-1,1e-4+1e-4j)),0)
        from radia.bem_complete_reaction import solve_complete_body
        with pytest.raises(RuntimeError, match='HACApK handles not built'):
            solve_complete_body(solver.wp_solver, solver._phi_poisson,
                np.zeros(solver.wp_solver.mesh.nv, complex),
                PanelSurfaceImpedance(np.full(len(solver.wp_tris), 1e-4+1e-4j)),
                2*np.pi*50000, lambda x: np.zeros_like(x, dtype=complex),
                hacapk=True)



@pytest.mark.parametrize('genus',[0,1])
def test_zero_real_impedance_faces_keep_total_field(setup_cases,genus):
    z=1e-4+1e-4j
    with ng.TaskManager():
        solver=make_solver(setup_cases,'filament',genus)
        centers=solver.wp_nodes[solver.wp_tris].mean(axis=1)
        values=np.where(centers[:,2]>=0,2j*z.imag,z)
        result=solve(solver,PanelSurfaceImpedance(values),genus)
        assert np.all(result['wp_q_tri'][values.real==0]==0)
        assert np.all(np.linalg.norm(result['wp_H_t_tri'][values.real==0],axis=1)>0)
        assert result['body_power_balance_relative_error']<.1


def test_frequency_refreshes_factor_but_reuses_geometric_work(setup_cases):
    with ng.TaskManager():
        solver=make_solver(setup_cases,'filament',1)
        values=np.where(solver.wp_nodes[solver.wp_tris].mean(axis=1)[:,2]>=0,
                        3e-4+3e-4j,2e-5+2e-5j)
        z=PanelSurfaceImpedance(values)
        solve(solver,z,1)
        geometry=solver.wp_solver._prepared_loop_geometry
        magnetic=solver.wp_solver._prepared_loop_magnetic_work
        factor=solver.wp_solver._prepared_loop_system
        result=solver.solve(z,2*np.pi*55000,max_iter=40,tol=1e-9,relax=.8,loop_dof=True)
        assert solver.wp_solver._prepared_loop_geometry is geometry
        assert solver.wp_solver._prepared_loop_magnetic_work is magnetic
        assert solver.wp_solver._prepared_loop_system is not factor
        fresh=make_solver(setup_cases,'filament',1).solve(z,2*np.pi*55000,
            max_iter=40,tol=1e-9,relax=.8,loop_dof=True)
        compare(result,fresh)
