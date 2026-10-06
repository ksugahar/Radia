"""Independent linear equivalence and falling-M retention for forward B(H)."""
import json
import numpy as np
import ngsolve as ng
from netgen.occ import Sphere, Pnt, OCCGeometry
from radia.vim import Solve

MU0 = 4e-7 * np.pi


def test_linear_equivalence_and_falling_m_retention(tmp_path):
    mesh = ng.Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), 1)).GenerateMesh(maxh=.7))
    with ng.TaskManager():
        linear = Solve(mesh, mu_r=3, H_ext=ng.CF((0, 0, 2)), order=1,
                       gram_eps=1e-10, tol=1e-10)
        forward = Solve(mesh, bh_table=[[0, 0], [20, MU0 * 60]],
                        H_ext=ng.CF((0, 0, 2)), order=1, gram_eps=1e-10,
                        nonlinear_solver="forward-newton", nl_tol=1e-9,
                        newton_inner_tol=1e-10, maxit=1000, nl_maxit=30)
        x = linear["gfM"].vec.FV().NumPy()
        y = forward["gfM"].vec.FV().NumPy()
        error = float(np.linalg.norm(x-y) / np.linalg.norm(x))
        assert error < 1e-7, error
        # M rises to 10 then falls to 5. Inverse-M clipping would never yield
        # an average above 5; this independent bounded-law check detects it.
        falling = Solve(mesh, bh_table=[[0, 0], [1, MU0 * 11], [10, MU0 * 15]],
                        H_ext=ng.CF((0, 0, 8)), order=1, gram_eps=1e-10,
                        nonlinear_solver="forward-newton", nl_tol=1e-9,
                        maxit=1000, nl_maxit=50)
    assert float(falling["M_avg"][2]) > 5.5, falling["M_avg"]
    for result in (forward, falling):
        assert result["nonlinear_final_relative_residual"] <= 1e-9
        assert result["constitutive_projection_relative_residual"] <= 1e-9
        assert result["material_magnetization_cap_applied"] is False
        assert all(a["true_relative_residual"] <= 1e-6 for a in result["linear_audit"])
    (tmp_path / "result.json").write_text(json.dumps(dict(
        passed=True, linear_coefficient_relative_error=error,
        falling_M_avg=falling["M_avg"].tolist(),
        linear_audit=falling["linear_audit"],
        nonlinear_residual=falling["nonlinear_final_relative_residual"])))
