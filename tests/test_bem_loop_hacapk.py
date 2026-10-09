"""Compressed genus-1 body contracts on an independently authored ring."""
import importlib.util
from pathlib import Path
import weakref
import numpy as np
import ngsolve as ng
import pytest
from radia import _radia_pybind as native
from radia.bem_sibc_solver import ScalarBIESIBCSolver
from radia.bem_loop_extension import solve_loop_extended
from radia.surface_impedance import PanelSurfaceImpedance


def test_hacapk_transpose_and_released_entry_lifetime():
    rng = np.random.default_rng(7)
    coords = np.ascontiguousarray(rng.normal(size=(48,3)))
    matrix = np.ascontiguousarray(rng.normal(size=(48,48)))
    expected = matrix.copy()
    owner = weakref.ref(matrix)
    manager = native.HACApKBEMManager(coords, matrix)
    with pytest.raises(RuntimeError, match="before releasing"):
        manager.ReleaseDenseEntries()
    assert manager.BuildHMatrix(aca_eps=1e-12,leaf_size=8,eta=2.,max_rank=48,print_level=0)
    x = rng.normal(size=48)
    np.testing.assert_allclose(manager.MatVec(x), expected @ x,rtol=1e-11,atol=1e-12)
    np.testing.assert_allclose(manager.MatVec(x,transpose=True),expected.T @ x,rtol=1e-11,atol=1e-12)
    manager.ReleaseDenseEntries()
    del matrix
    assert owner() is None
    np.testing.assert_allclose(manager.MatVec(x,transpose=True),expected.T @ x,rtol=1e-11,atol=1e-12)
    with pytest.raises(RuntimeError,match="released"):
        manager.BuildHMatrix()


@pytest.fixture(scope="module")
def ring():
    path=Path(__file__).resolve().parents[1]/"validation_test/induction_heating/run_loop_work_ring.py"
    spec=importlib.util.spec_from_file_location("ring_hacapk",path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    mesh,points,tri=module.ring_mesh(16)
    solvers=[]
    with ng.TaskManager():
        for compressed in (False,True):
            solvers.append(ScalarBIESIBCSolver(mesh,order=1,assemble_dense=not compressed,
                use_intree_bem=True,intree_geom_order=1,intree_singular_n_q=6,
                intree_regular_quad_degree=7,use_intree_hacapk=compressed,hacapk_aca_eps=1e-10))
    return solvers,points,tri


@pytest.mark.parametrize("jump",[False,True])
@pytest.mark.parametrize("offset",[0.,.00015])
def test_hacapk_loop_matches_dense_without_dense_system(ring,jump,offset):
    solvers,points,tri=ring
    omega=2*np.pi*5e4;mu0=4e-7*np.pi
    z=(1+1j)*np.sqrt(omega*mu0/(2*5.8e7))
    impedance=PanelSurfaceImpedance(np.where(points[tri].mean(axis=1)[:,0]>=0,3*z,.25*z)) if jump else z
    angles=np.arange(96)*2*np.pi/96
    carrier=np.c_[.03*np.cos(angles),.03*np.sin(angles),np.full(96,offset)]
    a_inc=lambda p: .5*mu0*np.c_[-p[:,1],p[:,0],np.zeros(len(p))]
    outputs=[]
    with ng.TaskManager():
        for compressed,solver in zip((False,True),solvers):
            outputs.append(solve_loop_extended(solver,-points[:,2].astype(complex),impedance,omega,
                a_inc,carrier_ring=carrier,section_anchor=(.03,offset),hacapk=compressed,_reuse_prepared=True))
    dense,compressed=outputs
    for key in ("P_total","alpha","reaction_integral","P_frozen"):
        assert abs(compressed[key]-dense[key])/max(abs(dense[key]),1e-30)<1e-4
    assert compressed['linear_residual_rel']<=1e-6
    assert compressed['faraday_residual_rel']<=1e-6
    assert compressed['loop_frozen_residual_rel']<=1e-6
    assert compressed['loop_work_diagnostics']['magnetic_work_symmetry']<=1e-6
    assert compressed['loop_body_backend']=='hacapk'
    assert compressed['loop_work_diagnostics']['body_hmatrix_controls']==dict(aca_eps=1e-10,leaf_size=64,eta=2.,max_rank=-1)
    assert solvers[1].SL is None and solvers[1].DL is None
    assert not isinstance(solvers[1].M_inv,np.ndarray)
    assert 'A2' not in solvers[1]._prepared_loop_system
    assert abs(compressed['P_total']-compressed['P_reaction'])/compressed['P_total']<.1

from tests.test_bem_strong_panel_impedance import setup_cases, make_solver
from tests.test_panel_esim_complete_state import Law, evaluate
from radia.bem_coupled_solver import CoupledBEMSolver
from radia.peec_coupled_bem_solver import CoupledPEECBEMSolver
from radia.esim_panel_evaluator import solve_panel_esim


@pytest.mark.parametrize('model',['surface','filament'])
@pytest.mark.parametrize('genus',[0,1])
def test_released_hacapk_maps_support_strong_panel_esim_reuse(setup_cases,model,genus):
    bodies,coil,paths=setup_cases
    omega=2*np.pi*50000;z=(1+1j)*np.sqrt(omega*4e-7*np.pi/(2*5.8e7))
    with ng.TaskManager():
        solver=(CoupledBEMSolver(coil,bodies[genus],wp_hacapk=True,wp_gmres_tol=1e-10)
                if model=='surface' else CoupledPEECBEMSolver(paths,np.diag([.002,.003]),
                    np.array([[3e-7,1e-7],[1e-7,3.5e-7]]),bodies[genus],wp_hacapk=True,wp_gmres_tol=1e-10))
        handles=(solver.wp_solver._SL_hacapk,solver.wp_solver._DL_hacapk)
        mass=solver.wp_solver._mass_inverse
        def solve_em(values):
            return solver.solve(values,omega,max_iter=40,tol=1e-8,relax=.8,loop_dof=bool(genus))
        impedance,state,record=solve_panel_esim(evaluate(Law(z,True)),solve_em,
            np.full(len(solver.wp_tris),z),tolerance=1e-5,max_iter=30,
            relaxation=.8,inner_tolerance=1e-8)
    assert record['esim_panel_evaluation']['certification_cell_calls']==len(solver.wp_tris)
    assert record['esim_inner_history'][-1]['final_certification']
    np.testing.assert_array_equal(state['Z_s_per_panel'],impedance.values)
    assert (solver.wp_solver._SL_hacapk,solver.wp_solver._DL_hacapk)==handles
    assert solver.wp_solver._mass_inverse is mass
    if genus:
        assert state['wp_loop_body_backend']=='hacapk'
        assert state['wp_loop_faraday_residual']<=1e-6
    assert solver.wp_solver.SL is None and solver.wp_solver.DL is None

from tests.test_bem_loop_work_pairings import surface
from tests.test_bem_loop_work_pairings import test_complete_work_changes_covariantly_with_generic_lift as check_covariance


@pytest.mark.parametrize('zero_mean',[False,True])
def test_compressed_work_lift_transform_including_gauge(surface,zero_mean,monkeypatch):
    bem,points,*_=surface
    for name in ('SL','DL'):
        handle=native.HACApKBEMManager(np.ascontiguousarray(points),np.ascontiguousarray(getattr(bem,name)))
        assert handle.BuildHMatrix(aca_eps=1e-10,leaf_size=4,eta=2.,max_rank=4,print_level=0)
        handle.ReleaseDenseEntries()
        monkeypatch.setattr(bem,'_'+name+'_hacapk',handle)
    monkeypatch.setattr(bem,'_loop_body_backend','hacapk',raising=False)
    check_covariance(surface,zero_mean)

def test_hacapk_certifies_physical_residual_after_gmres(ring,monkeypatch):
    import scipy.sparse.linalg as sparse_linalg
    solvers,points,_=ring
    omega=2*np.pi*5e4;mu0=4e-7*np.pi
    z=(1+1j)*np.sqrt(omega*mu0/(2*5.8e7))
    monkeypatch.setattr(sparse_linalg,'gmres',lambda operator,rhs,**kwargs:(np.zeros_like(rhs),0))
    a_inc=lambda p:.5*mu0*np.c_[-p[:,1],p[:,0],np.zeros(len(p))]
    with ng.TaskManager(),pytest.raises(RuntimeError,match='true residual'):
        solve_loop_extended(solvers[1],-points[:,2].astype(complex),z,omega,a_inc,hacapk=True)


def test_hacapk_keeps_construction_guard_and_rejects_missing_handles():
    from radia.panels import calc_inductance as ci
    from types import SimpleNamespace
    args=SimpleNamespace(wp_bem_backend='hacapk',impedance_model='sibc',h1_order=1)
    with pytest.raises(ValueError,match='7000'):
        ci._resolve_workpiece_backend(args,1,7001)
    with pytest.raises(ValueError,match='multiple holes'):
        ci._resolve_workpiece_backend(args,2,100)
    assert args.wp_bem_backend=='hacapk'

@pytest.mark.parametrize('panel',[False,True])
def test_genus0_hacapk_released_storage_matches_retained_storage(setup_cases,panel):
    from radia.bem_sibc_solver import telegen_extract_coil_LR
    body=setup_cases[0][0]
    omega=2*np.pi*5e4;mu0=4e-7*np.pi
    z=(1+1j)*np.sqrt(omega*mu0/(2*5.8e7))
    solvers=[];outputs=[]
    with ng.TaskManager():
        for retain in (True,False):
            solver=ScalarBIESIBCSolver(body,order=1,assemble_dense=retain,
                use_intree_bem=True,use_intree_hacapk=True,hacapk_aca_eps=1e-10)
            solvers.append(solver)
            impedance=(PanelSurfaceImpedance(np.full(len(list(body.Elements(ng.BND))),z))
                       if panel else z)
            outputs.append(solver.solve_hacapk(-ng.z,impedance,omega,tol=1e-12,maxiter=500,restart=100))
    old,new=outputs
    for key in ('H_t_rms','P_density','area'):
        assert abs(new[key]-old[key])<=1e-10*abs(old[key])
    assert np.linalg.norm(new['phi_vec']-old['phi_vec'])<=1e-10*np.linalg.norm(old['phi_vec'])
    assert new['linear_residual_rel']<=1e-6
    assert solvers[1].SL is None and solvers[1].DL is None
    before,after=[telegen_extract_coil_LR(s,o['phi_vec'],1.,omega,z) for s,o in zip(solvers,outputs)]
    for key in ('R_port','L_port','P_diss'):
        assert abs(after[key]-before[key])<=1e-10*abs(before[key])
    print('genus0-storage',panel,'heat-before',old['P_density']*old['area'],
          'heat-after',new['P_density']*new['area'],'field-relative',
          np.linalg.norm(new['phi_vec']-old['phi_vec'])/np.linalg.norm(old['phi_vec']))
