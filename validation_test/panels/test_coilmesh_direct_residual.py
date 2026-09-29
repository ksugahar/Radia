"""Complex A-V solves must pass the free-row residual before post-processing."""
import pytest
from netgen.occ import Box, Pnt, OCCGeometry, X

from radia.panels import calc_fem_coilmesh as solver


@pytest.mark.parametrize('model,fail_call', [('sibc', None), ('esim', None),
                                              ('sibc', 1), ('esim', 2)])
def test_coilmesh_checks_each_dirichlet_solve(tmp_path, monkeypatch, model, fail_call):
    body = Box(Pnt(0, 0, 0), Pnt(.02, .01, .01))
    body.mat('coil')
    body.faces.name = 'sibc'
    body.faces.Min(X).name = 'source'
    body.faces.Max(X).name = 'sink'
    mesh = OCCGeometry(body).GenerateMesh(maxh=.006)
    vol = tmp_path / 'coil.vol'
    mesh.Save(str(vol))
    # A constant surface impedance isolates the final ESIM re-solve contract.
    class ConstantImpedance:
        def solve(self, *args, **kwargs):
            return {'Z': .001 + .001j}
    monkeypatch.setattr(solver.EMMaterial, 'create_esim_solver',
                        lambda *args, **kwargs: ConstantImpedance())
    bh = tmp_path / 'bh.txt'
    bh.write_text('0 0\n100 0.1\n1000 0.5\n', encoding='utf-8')
    original = solver.apply_fe_inverse
    calls = []
    def checked(matrix, inverse, rhs, solution, free):
        calls.append(None)
        if len(calls) == fail_call:
            inverse = 0 * inverse
        return original(matrix, inverse, rhs, solution, free)
    monkeypatch.setattr(solver, 'apply_fe_inverse', checked)
    def solve():
        return solver.solve_fem_coilmesh(
            str(vol), frequency=1000., I_target=1., coil_sigma=1e5,
            wp_sigma=1e6, wp_mu_r=1., half_thickness=.01,
            impedance_model=model, bh_file=str(bh) if model == 'esim' else '',
            esim_max_iter=1)
    if fail_call is not None:
        with pytest.raises(RuntimeError, match='true residual check'):
            solve()
        assert len(calls) == fail_call
        return
    result = solve()
    residuals = result['linear_true_relative_residuals']
    assert len(residuals) == (2 if model == 'esim' else 1)
    assert all(0 <= value <= 1e-8 for value in residuals)
    assert result['P_total_W'] > 0
