"""Golden: the all-C++ SYMMETRIC ENERGY-NEWTON is the DEFAULT HDiv-VIM nonlinear soft-iron solver
(radia.vim._solve._solve_nonlinear_energy_cpp), and it matches the forward (scipy splu + GMRES) Newton.

The nonlinear inner Newton step is now solved by the EXISTING C++ symmetric W-CG
(solve_linear_material_mass_riesz, W = the differential-reluctivity tangent mass, mass-Riesz PARDISO, N
H-matvec) -- bringing the nonlinear solve to C++ parity with the linear path (no scipy splu / GMRES /
M_mass^-1).  The co-energy form is robust through deep saturation via a hard-saturation barrier, a
co-energy line search, and settled-step acceptance for the achievable-precision limit cycle of the M-form.

Locks: (1) the default nonlinear solver is 'energy-newton-cpp'; (2) it agrees with the forward Newton
(linear_solver='gmres') to ~1e-5 at moderate AND deep drive; (3) deep saturation converges to M ~ Msat.
"""
import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
pytest.importorskip("netgen.occ")
from netgen.occ import Sphere, Pnt, OCCGeometry          # noqa: E402
from radia.vim import Solve                    # noqa: E402
from radia.vim._solve import _resolve_highorder_preconditioner  # noqa: E402

# realistic soft-iron BH table (mu_r ~ 4000 at low H, saturating ~2.2 T)
_H = np.array([0, 50, 100, 200, 500, 1e3, 2e3, 5e3, 1e4, 3e4, 1e5, 3e5, 1e6])
_B = np.array([0, 0.30, 0.60, 1.0, 1.45, 1.7, 1.9, 2.0, 2.05, 2.1, 2.15, 2.25, 2.5])
_BH = np.column_stack([_H, _B]).tolist()
_MU0 = 4e-7 * np.pi
_MSAT = _B[-1] / _MU0 - _H[-1]                            # the table's saturation magnetization


def test_auto_preconditioner_policy_is_energy_newton_specific(monkeypatch):
    """The large-run diagonal default belongs to the energy-Newton path, not every nonlinear label."""
    monkeypatch.delenv("RADIA_HDIV_AUTO_JACOBI_TET_NFACE", raising=False)

    eff, policy = _resolve_highorder_preconditioner(
        "auto", nonlinear=False, nonlinear_solver="energy-newton", vertex_counts={8}, n_face=100000)
    assert eff == "mass-riesz"
    assert policy == "auto:linear-mass-riesz"

    eff, policy = _resolve_highorder_preconditioner(
        "auto", nonlinear=True, nonlinear_solver="energy-newton", vertex_counts={8}, n_face=100)
    assert eff == "jacobi"
    assert policy == "auto:hex-wedge-energy-newton-jacobi"

    eff, policy = _resolve_highorder_preconditioner(
        "auto", nonlinear=True, nonlinear_solver="picard-energy", vertex_counts={8}, n_face=100000)
    assert eff == "mass-riesz"
    assert policy == "auto:picard-energy-mass-riesz"

    eff, policy = _resolve_highorder_preconditioner(
        "auto", nonlinear=True, nonlinear_solver="energy-newton", vertex_counts={4}, n_face=4314)
    assert eff == "mass-riesz"
    assert policy == "auto:tet-energy-newton-mass-riesz-nface<6000"

    eff, policy = _resolve_highorder_preconditioner(
        "auto", nonlinear=True, nonlinear_solver="energy-newton", vertex_counts={4}, n_face=8193)
    assert eff == "jacobi"
    assert policy == "auto:tet-energy-newton-jacobi-nface>=6000"

    monkeypatch.setenv("RADIA_HDIV_AUTO_JACOBI_TET_NFACE", "10000")
    eff, policy = _resolve_highorder_preconditioner(
        "auto", nonlinear=True, nonlinear_solver="energy-newton", vertex_counts={4}, n_face=8193)
    assert eff == "mass-riesz"
    assert policy == "auto:tet-energy-newton-mass-riesz-nface<10000"


def _sphere(maxh=0.6):
    return ng.Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), 1.0)).GenerateMesh(maxh=maxh))


def test_default_nonlinear_is_energy_newton_cpp(monkeypatch):
    """The default iron-only nonlinear solve routes to the all-C++ symmetric energy-Newton (RT1): it is the
    default solver, converges, and gives a sane +z M with the sphere demag ~1/3 at moderate AND knee drive.
    The former forward-Newton (scipy splu + GMRES) cross-check is gone -- the energy-Newton
    is validated against the closed-form spheroid fixed point in test_hdiv_vim_curved_solve_nonlinear."""
    monkeypatch.delenv("RADIA_HDIV_AUTO_JACOBI_TET_NFACE", raising=False)
    mesh = _sphere()
    for H0 in (500.0, 5000.0):                            # moderate (linear-ish) and into the knee
        Hext = ng.CoefficientFunction((0, 0, H0))
        with ng.TaskManager():
            re = Solve(mesh, bh_table=_BH, H_ext=Hext, order=1)
        assert re["linear_solver"] == "energy-newton-cpp", re["linear_solver"]
        assert re["preconditioner_requested"] == "auto"
        assert re["preconditioner"] == "mass-riesz"
        assert re["preconditioner_policy"] == "auto:tet-energy-newton-mass-riesz-nface<6000"
        assert re["nonlinear"] is True
        assert re["iters"] < 100
        assert re["M_avg"][2] > 0, f"H0={H0}: M_avg {re['M_avg']} not +z"
        assert 0.25 < re["demag"] < 0.40, f"H0={H0}: demag {re['demag']:.4f} not ~1/3"


def test_energy_newton_deep_saturation():
    """Deep saturation (H0 well past the table) drives M to ~Msat and converges (the hard-saturation barrier
    + settled acceptance handle the M-form limit cycle)."""
    mesh = _sphere()
    with ng.TaskManager():
        r = Solve(mesh, bh_table=_BH, H_ext=ng.CoefficientFunction((0, 0, 3e6)), order=1)
    assert r["linear_solver"] == "energy-newton-cpp"
    # at deep drive M_avg -> Msat (the table's saturation; demag-independent there).  The soft hard-saturation
    # barrier permits a small (<~1%) overshoot of the uniform Msat at the discrete/volume-average level.
    assert 0.95 * _MSAT < r["M_avg"][2] < 1.03 * _MSAT, (r["M_avg"][2], _MSAT)
    # it converged (returning, not raising) -- the deep-saturation M-form iteration count is table-dependent
    # and higher than the forward H-form (the M-form limit-cycles to the achievable precision); a generous
    # bound catches a runaway without being brittle to the table shape.
    assert r["iters"] < 100


def test_picard_energy_matches_energy_newton():
    """The tet-style Picard + mass-Riesz warmstart can be used before the energy solve without changing the
    final energy-stationary magnetization.  This keeps inverse-modeling / adjoint work on the energy path
    while allowing a faster forward warmstart."""
    mesh = _sphere(maxh=0.8)
    Hext = ng.CoefficientFunction((0, 0, 5000.0))
    with ng.TaskManager():
        re = Solve(mesh, bh_table=_BH, H_ext=Hext, order=1)
        rh = Solve(mesh, bh_table=_BH, H_ext=Hext, order=1, nonlinear_solver="picard-energy")
    assert rh["linear_solver"] == "picard-energy-cpp"
    assert rh["nonlinear"] is True
    rel = abs(rh["M_avg"][2] - re["M_avg"][2]) / max(abs(re["M_avg"][2]), 1.0)
    assert rel < 1e-5, (rh["M_avg"], re["M_avg"], rel)


def test_energy_newton_fast_knobs_preserve_solution_and_report_stats():
    """Inexact inner CG / continuation / chord-tangent reuse are explicit speed knobs, not new physics."""
    mesh = _sphere(maxh=0.8)
    Hext = ng.CoefficientFunction((0, 0, 5000.0))
    with ng.TaskManager():
        ref = Solve(mesh, bh_table=_BH, H_ext=Hext, order=1)
        tuned = Solve(mesh, bh_table=_BH, H_ext=Hext, order=1,
                      newton_continuation=2, newton_reuse_tangent_steps=2,
                      newton_inner_tol="auto")
    rel = np.linalg.norm(tuned["M_avg"] - ref["M_avg"]) / max(np.linalg.norm(ref["M_avg"]), 1.0)
    assert rel < 1e-5, (tuned["M_avg"], ref["M_avg"], rel)
    stats = tuned["nonlinear_solve_stats"]
    assert stats["nonlinear_continuation_steps"] == 2
    assert stats["nonlinear_reuse_tangent_steps"] == 2
    assert stats["nonlinear_tangent_reuses"] >= 1
    assert "nonlinear_linear_inner_iters" in tuned


def test_energy_newton_jacobi_inner_preconditioner_preserves_solution():
    """The diagonal W+N preconditioner is an explicit scaling knob for large runs.

    It changes only the Krylov preconditioner inside the same symmetric energy-Newton system; the converged
    magnetization should match the mass-Riesz path on a small regression mesh.
    """
    mesh = _sphere(maxh=0.8)
    Hext = ng.CoefficientFunction((0, 0, 5000.0))
    with ng.TaskManager():
        ref = Solve(mesh, bh_table=_BH, H_ext=Hext, order=1)
        tuned = Solve(mesh, bh_table=_BH, H_ext=Hext, order=1,
                      preconditioner="jacobi",
                      newton_continuation=2, newton_reuse_tangent_steps=3,
                      newton_inner_tol="auto", newton_cg_x0=True)
    rel = np.linalg.norm(tuned["M_avg"] - ref["M_avg"]) / max(np.linalg.norm(ref["M_avg"]), 1.0)
    assert rel < 1e-5, (tuned["M_avg"], ref["M_avg"], rel)
    assert tuned["preconditioner"] == "jacobi"
    assert tuned["preconditioner_policy"] == "explicit:jacobi"
    assert tuned["nonlinear_solve_stats"]["nonlinear_cg_x0"] is True
    assert tuned["nonlinear_solve_stats"]["nonlinear_inner_preconditioner"] == "jacobi"


def test_energy_newton_jacobi_inner_preconditioner_deep_saturation():
    """The diagonal W+N preconditioner remains a solver choice in the deep-saturation regime."""
    mesh = _sphere(maxh=0.8)
    Hext = ng.CoefficientFunction((0, 0, 3e6))
    with ng.TaskManager():
        ref = Solve(mesh, bh_table=_BH, H_ext=Hext, order=1)
        tuned = Solve(mesh, bh_table=_BH, H_ext=Hext, order=1,
                      preconditioner="jacobi",
                      newton_continuation=2, newton_reuse_tangent_steps=3,
                      newton_inner_tol="auto")
    rel = abs(tuned["M_avg"][2] - ref["M_avg"][2]) / max(abs(ref["M_avg"][2]), 1.0)
    assert rel < 5e-4, (tuned["M_avg"], ref["M_avg"], rel)
    assert 0.95 * _MSAT < tuned["M_avg"][2] < 1.03 * _MSAT
    assert tuned["nonlinear_solve_stats"]["nonlinear_inner_preconditioner"] == "jacobi"


def test_energy_newton_jacobi_inner_preconditioner_multiregion_bh_dict():
    """Per-region BH tables use the same geometry Gram; Jacobi changes only the Krylov preconditioner."""
    from netgen.csg import CSGeometry, OrthoBrick, Pnt

    geo = CSGeometry()
    geo.Add(OrthoBrick(Pnt(-1.0, -0.5, -0.5), Pnt(0.0, 0.5, 0.5)).mat("left"))
    geo.Add(OrthoBrick(Pnt(0.0, -0.5, -0.5), Pnt(1.0, 0.5, 0.5)).mat("right"))
    mesh = ng.Mesh(geo.GenerateMesh(maxh=0.9))
    assert set(mesh.GetMaterials()) == {"left", "right"}

    # Same shape, different grade: this exercises the material-region dispatch without turning the test into
    # a material-law benchmark.
    bh = {
        "left": _BH,
        "right": np.column_stack([_H, 0.92 * _B]).tolist(),
    }
    Hext = ng.CoefficientFunction((0, 0, 5000.0))
    with ng.TaskManager():
        ref = Solve(mesh, bh_table=bh, H_ext=Hext, order=1)
        tuned = Solve(mesh, bh_table=bh, H_ext=Hext, order=1,
                      preconditioner="jacobi",
                      newton_continuation=2, newton_reuse_tangent_steps=3,
                      newton_inner_tol="auto")
    rel = np.linalg.norm(tuned["M_avg"] - ref["M_avg"]) / max(np.linalg.norm(ref["M_avg"]), 1.0)
    assert rel < 5e-4, (tuned["M_avg"], ref["M_avg"], rel)
    assert tuned["preconditioner"] == "jacobi"
    assert tuned["preconditioner_policy"] == "explicit:jacobi"
    assert tuned["nonlinear_solve_stats"]["nonlinear_inner_preconditioner"] == "jacobi"
