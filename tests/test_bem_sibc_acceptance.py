"""BEM SIBC coefficient support and augmented-equation acceptance contracts."""
import numpy as np
import pytest

from radia import bem_sibc_solver as module


def _two_mode_surface():
    solver = module.ScalarBIESIBCSolver.__new__(module.ScalarBIESIBCSolver)
    solver._intree_lagrange_p2 = True
    solver.ndof = 2
    solver.M = np.eye(2)
    solver.M_inv = np.eye(2)
    solver.K = np.array([[1., -1.], [-1., 1.]])
    solver.SL = np.eye(2)
    solver.DL = np.zeros((2, 2))
    solver._c_gauge = np.ones(2)
    return solver


@pytest.mark.parametrize("impedance", [0j, 0.01+0.02j])
def test_uniform_surface_matches_independent_eigenmode_solution(impedance):
    solver = _two_mode_surface()
    omega = 1000.
    source = np.array([1., -1.])
    expected = source / (0.5 + 2*impedance/(1j*omega*module.MU_0))
    result = solver.solve(source, impedance, omega)
    repeated = solver.solve(source, np.full(2, impedance), omega)
    np.testing.assert_allclose(result['phi_vec'], expected, rtol=1e-13)
    np.testing.assert_array_equal(result['phi_vec'], repeated['phi_vec'])
    assert result['linear_residual_rel'] <= module.RELATIVE_LIMIT
    assert result['linear_residual_limit'] == module.RELATIVE_LIMIT


def test_variable_surface_is_rejected_before_wrong_observation_scaling():
    solver = _two_mode_surface()
    with pytest.raises(NotImplementedError, match='source-side coefficient-weighted'):
        solver.solve(np.array([1., -1.]), np.array([0.01, 0.02]), 1000.)


@pytest.mark.parametrize("info", [1, -1])
def test_gmres_failure_is_never_accepted_even_with_zero_residual(info):
    with pytest.raises(RuntimeError, match='GMRES failed'):
        module._check_bie_solution(np.eye(2), np.ones(2), np.ones(2), info=info)


@pytest.mark.parametrize("scale", [1., 1e-20])
def test_true_residual_uses_fixed_load_scale(scale):
    rhs = scale*np.array([1., 0.])
    with pytest.raises(RuntimeError, match='true relative residual'):
        module._check_bie_solution(np.eye(2), rhs*1.000002, rhs)
    assert module._check_bie_solution(np.eye(2), rhs, rhs) == 0.


def test_dense_bad_solution_is_rejected(monkeypatch):
    solver = _two_mode_surface()
    monkeypatch.setattr(module, 'scipy_solve', lambda matrix, rhs: np.zeros_like(rhs))
    with pytest.raises(RuntimeError, match='true relative residual'):
        solver._solve_with_gauge(np.eye(2), np.array([1., -1.]))

@pytest.mark.parametrize('order,symbol', [(1,'_AssembleSLDL_Galerkin'),(2,'_AssembleSLDL_Galerkin_P2')])
def test_selected_native_bem_rejects_missing_symbol(monkeypatch,order,symbol):
    import ngsolve as ng
    from netgen.occ import OCCGeometry, Sphere, Pnt
    from radia import _radia_pybind
    monkeypatch.delattr(_radia_pybind,symbol,raising=False)
    with ng.TaskManager():
        mesh=ng.Mesh(OCCGeometry(Sphere(Pnt(0,0,0),1)).GenerateMesh(maxh=1))
        with pytest.raises(RuntimeError,match='Rebuild/reinstall'):
            module.ScalarBIESIBCSolver(mesh,order=order,use_intree_bem=True)


def test_native_segment_field_missing_symbol_and_analytic_normal_path(monkeypatch):
    from radia import _radia_pybind
    segments=np.array([[[-1.,0.,0.],[1.,0.,0.]]])
    point=np.array([[0.,1.,0.]])
    field=module._h_segments_complex(segments,point,np.array([1.+2j]))
    # Finite straight wire at unit distance, endpoint angles +/- pi/4.
    expected=np.array([[0.,0.,np.sqrt(2)/(4*np.pi)*(1+2j)]])
    np.testing.assert_allclose(field,expected,rtol=1e-13,atol=1e-15)
    monkeypatch.delattr(_radia_pybind,'_HFromSegmentsComplex')
    with pytest.raises(RuntimeError,match='Rebuild/reinstall'):
        module._h_segments_complex(segments,point,np.array([1.+2j]))
