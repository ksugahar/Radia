"""Golden: the PROMOTED planar (2D) HDiv-VIM solve layer (radia.vim._vim2d, dispatched by
vim.Solve on mesh.dim == 2).

Locks (the research-layer gates promoted to production, see memory hdiv-vim-tri-quad-motor):
  (1) LINEAR disk: Solve(mesh2d, mu_r) volume-average M matches the 2D
      Clausius-Mossotti relation M/H0 = chi / (1 + chi/2), demag factors (1/2, 1/2);
  (2) NONLINEAR disk via a saturating bh_table: deep-saturation sweep matches the analytic
      uniform fixed point M = Mof(H0 - M/2) (D = 1/2 exact) -- the scalar-chi Picard +
      safeguarded Anderson(1) converges from the linear regime through deep saturation;
  (3) ELLIPSE reluctance torque 3-way: closed form vs the volume torque mu0*A*(M x H0) vs the
      Maxwell-stress circle on the analytic charge field (maxwell_torque_circle + H_at) --
      the rotor-sweep machinery (N built once, per-angle solves) under golden;
  (4) fail-loud contract: 3D-only knobs raise on a 2D mesh; per-region dicts raise;
      non-convergence raises (no silent partial result).
Self-contained (OCC disk / ellipse); ~40 s.
"""
from __future__ import annotations

import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
pytest.importorskip("netgen.occ")
from netgen.occ import WorkPlane, OCCGeometry  # noqa: E402

from radia.vim import Solve, maxwell_torque_circle  # noqa: E402

MU0 = 4e-7 * np.pi
CHI0, MSAT = 1000.0, 1.2e6


def _disk_mesh(maxh=0.3, R=1.0):
    wp = WorkPlane().Circle(0, 0, R).Face()
    return ng.Mesh(OCCGeometry(wp, dim=2).GenerateMesh(maxh=maxh))


def _bh_table():
    """Sample the saturating law M(H) = chi0 H / (1 + chi0 H / Msat) as a [[H, B]] table dense
    enough that table interpolation is far below the gate tolerances."""
    H = np.logspace(0, 7.5, 400)
    M = CHI0 * H / (1.0 + CHI0 * H / MSAT)
    B = MU0 * (H + M)
    return np.stack([H, B], axis=1).tolist()


def _fixed_point(H0, D=0.5):
    """Analytic uniform-body root M = Mof(H0 - D*M) for the saturating law (odd law, monotone
    bracket [0, H0/D])."""
    Mof = lambda h: CHI0 * h / (1.0 + CHI0 * abs(h) / MSAT)
    lo, hi = 0.0, H0 / D
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if mid - Mof(H0 - D * mid) > 0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def test_2d_linear_disk_clausius_mossotti():
    mesh = _disk_mesh()
    with ng.TaskManager():
        res = Solve(mesh, 1000.0, ng.CoefficientFunction((1.0, 0.0)))
    assert res["linear_solver"] == "mass-riesz-cg-2d" and not res["nonlinear"]
    assert res["linear_iterations"] > 0
    chi = 999.0
    ref = chi / (1.0 + chi / 2.0)
    rel = abs(res["M_avg"][0] - ref) / ref
    assert rel < 3e-3, f"linear disk M/H0 {res['M_avg'][0]:.5f} vs CM2D {ref:.5f} (rel {rel:.1e})"
    Dx, Dy = res["demag_factors"]
    assert abs(Dx - 0.5) < 2e-3 and abs(Dy - 0.5) < 2e-3, f"disk demag ({Dx:.4f},{Dy:.4f}) off 1/2"


@pytest.mark.parametrize("H0,tol", [(1e3, 1e-3), (1e5, 1e-3), (1e6, 1e-3)])
def test_2d_nonlinear_disk_deep_saturation(H0, tol):
    """Linear regime -> knee -> deep saturation (M/Msat ~ 0.997) vs the analytic fixed point."""
    mesh = _disk_mesh()
    with ng.TaskManager():
        res = Solve(mesh, None, ng.CoefficientFunction((H0, 0.0)),
                               bh_table=_bh_table())
    assert res["nonlinear"] and res["iters"] >= 1
    Mref = _fixed_point(H0)
    rel = abs(res["M_avg"][0] - Mref) / Mref
    assert rel < tol, f"H0={H0:.0e}: M {res['M_avg'][0]:.1f} vs fixed point {Mref:.1f} (rel {rel:.1e})"


def test_2d_ellipse_reluctance_torque_three_way():
    """N built once; per-angle solve; torque agrees 3-way (closed form / volume / Maxwell circle)."""
    a_el, b_el = 0.2, 0.1
    Na, Nb = b_el / (a_el + b_el), a_el / (a_el + b_el)
    area = np.pi * a_el * b_el
    chi, H0 = 1000.0, 1e5
    wp = WorkPlane().Ellipse(a_el, b_el).Face()
    mesh = ng.Mesh(OCCGeometry(wp, dim=2).GenerateMesh(maxh=b_el / 3))
    body = None
    with ng.TaskManager():
        for th in (30.0, 60.0):
            Ha, Hb = H0 * np.cos(np.radians(th)), H0 * np.sin(np.radians(th))
            if body is None:
                res = Solve(mesh, chi + 1.0, ng.CoefficientFunction((Ha, Hb)))
                body = res["body"]                       # N built ONCE; reuse for the sweep
                m = res["m"]
            else:
                m = body.solve_linear(chi, body.project(ng.CoefficientFunction((Ha, Hb))))
            Ma = chi * Ha / (1 + chi * Na)
            Mb = chi * Hb / (1 + chi * Nb)
            T_ref = MU0 * area * (Ma * Hb - Mb * Ha)
            Mx, My = body.M_avg(m)
            T_vol = MU0 * area * (Mx * Hb - My * Ha)

            def H_tot(P):
                H = body.H_at(P, m)
                H[:, 0] += Ha
                H[:, 1] += Hb
                return H
            T_mx = maxwell_torque_circle(H_tot, 0.3)
            assert abs(T_vol - T_ref) / abs(T_ref) < 1e-2, \
                f"th={th}: volume torque {T_vol:.2f} vs closed form {T_ref:.2f}"
            assert abs(T_mx - T_ref) / abs(T_ref) < 2e-2, \
                f"th={th}: Maxwell torque {T_mx:.2f} vs closed form {T_ref:.2f}"


def test_2d_fail_loud_contract():
    mesh = _disk_mesh(maxh=0.6)
    H = ng.CoefficientFunction((1.0, 0.0))
    with ng.TaskManager():
        with pytest.raises(ValueError, match="linear_solver"):
            Solve(mesh, 1000.0, H, linear_solver="cpp-cg")
        with pytest.raises(ValueError, match="3D knob"):
            Solve(mesh, 1000.0, H, gram_eps=1e-6)
        with pytest.raises(NotImplementedError, match="per-region"):
            Solve(mesh, {"default": 1000.0}, H)
        with pytest.raises(ValueError, match="EXACTLY ONE"):
            Solve(mesh, 1000.0, H, bh_table=_bh_table())
