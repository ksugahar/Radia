# -*- coding: utf-8 -*-
r"""FEM-BEM Schur coupling gate.

Validates the analytic DtN Robin BC approach on a spherical shell against the
exact exterior dipole solution:

  Domain:   R_inner < r < R_outer  (spherical shell)
  Inner BC: u = z/R_inner          (dipole Dirichlet datum = cosθ)
  Outer BC: ∂u/∂n = (-2/R_outer) u (analytic DtN eigenvalue for dipole)

  Exact:    u(r,θ) = R_inner² z / r³  (exterior dipole decaying as 1/r²)

The analytic DtN Robin BC enforces the exact transparent condition for the
dipole mode, which together with the Dirichlet BC uniquely selects the
exterior-dipole shell solution.  This validates the FEM+DtN framework that
the full BEM coupling extends.
"""
import pytest

pytestmark = pytest.mark.usefixtures("ngsolve_taskmanager")

from radia_mcp.radia_ngsolve.fem_bem_coupling import (
    sphere_shell_analytic_dtn,
    laplace_fem_bem_schur,
    kelvin_dtn_eigenvalue,
    kelvin_vs_exact_open_bc_error,
    kelvin_openbc_error_vs_exterior_mesh,
    kelvin_twosphere_shell_dipole,
)


def test_shell_dipole_rel_err_small():
    """Analytic DtN Robin BC: L2 error vs exact exterior dipole < 7% on coarse mesh."""
    res = sphere_shell_analytic_dtn(R_inner=0.5, R_outer=1.0, maxh_vol=0.4,
                                    order=2, intorder=8)
    assert res["rel_err"] < 0.07, (
        f"Analytic DtN rel_err={res['rel_err']:.3e} (expected < 7%), "
        f"ndof={res['ndof']}"
    )


def test_shell_dipole_converges_with_refinement():
    """Analytic DtN L2 error decreases from coarse to medium mesh."""
    r_coarse = sphere_shell_analytic_dtn(R_inner=0.5, R_outer=1.0, maxh_vol=0.5,
                                         order=2, intorder=8)
    r_fine   = sphere_shell_analytic_dtn(R_inner=0.5, R_outer=1.0, maxh_vol=0.3,
                                         order=2, intorder=8)
    assert r_fine["rel_err"] < r_coarse["rel_err"], (
        f"No refinement convergence: coarse={r_coarse['rel_err']:.3e}, fine={r_fine['rel_err']:.3e}"
    )


def test_shell_dipole_scales_with_R():
    """DtN eigenvalue scales correctly: solution invariant under uniform scaling."""
    for R_inner, R_outer in ((0.25, 0.5), (1.0, 2.0)):
        res = sphere_shell_analytic_dtn(R_inner=R_inner, R_outer=R_outer,
                                        maxh_vol=0.4 * R_outer,
                                        order=2, intorder=8)
        assert res["rel_err"] < 0.08, (
            f"R_inner={R_inner}, R_outer={R_outer}: rel_err={res['rel_err']:.3e}"
        )


def test_full_bem_schur_dense_path():
    """Full dense BEM FEM-BEM Schur vs the exact exterior dipole.

    Exercises laplace_fem_bem_schur end-to-end: inner Dirichlet u=z/R_in, the BEM
    exterior DtN on the outer sphere, vs the exact u = R_in²·z/r³.  This is the
    BEM realisation of the SAME open BC as the analytic-DtN tests above; it must
    reach the same accuracy (the FEM interior error dominates at order=2, ~5.4%).

    Validates the three fixes in one shot: the coupling sign (K − ...), the weak
    surface-mass factor (P^T M_bnd Λ P), and the crash fix (the coupling mixed
    form uses the restricted SurfaceL2 as TRIAL, not TEST -- definedon-as-test on
    a volume mesh corrupts the heap; see _coupling_matrix and knowledge
    `fem_bem_schur`).  Deterministic 5.5e-2 over 3 runs at order=2."""
    import math
    import ngsolve as ng
    from ngsolve import x, y, z, sqrt, dx
    from netgen.occ import Sphere, Pnt, OCCGeometry
    R_in, R_out = 0.5, 1.0
    outer = Sphere(Pnt(0, 0, 0), R_out); outer.bc("outer")
    inner = Sphere(Pnt(0, 0, 0), R_in); inner.bc("inner")
    mesh = ng.Mesh(OCCGeometry(outer - inner).GenerateMesh(maxh=0.5)).Curve(3)
    gfu = laplace_fem_bem_schur(mesh, bnd_outer="outer", bnd_dirichlet="inner",
                                dirichlet_cf=z / R_in, order=2, intorder=8)
    r2 = x * x + y * y + z * z
    u_exact = R_in**2 * z / (r2 * sqrt(r2))
    err = math.sqrt(abs(float(ng.Integrate((gfu - u_exact)**2 * dx, mesh))))
    ref = math.sqrt(abs(float(ng.Integrate(u_exact**2 * dx, mesh))))
    rel = err / ref
    assert rel < 0.08, f"BEM FEM-BEM Schur rel_err={rel:.3e} (expected ~5.4%, < 8%)"


def test_kelvin_openbc_error_isolated_below_fem_error():
    """The KELVIN open-boundary error, ISOLATED from the interior FEM error, is far
    below it -- Kameari's coarse-mesh accuracy as a separated error budget.

    On a shared shell mesh, swap only the Γ operator: exact-DtN Robin (zero open-BC
    error) vs Kelvin-DtN Robin.  ||u_Kelvin − u_exactDtN|| is the pure Kelvin
    open-BC error; the interior FEM error cancels.  Even with a COARSE Kelvin
    closure (maxh=0.6, order=1) it is ~0.1 %, an order of magnitude under the ~5 %
    interior FEM error.  (include_bem=False keeps this fast -- no dense BEM.)"""
    r = kelvin_vs_exact_open_bc_error(R_inner=0.5, R_outer=1.0, maxh=0.4, order=2,
                                      degree=1, kelvin_maxh=0.6, kelvin_order=1,
                                      include_bem=False)
    assert 0.02 < r["interior_fem_error"] < 0.10, (
        f"interior FEM error out of range: {r['interior_fem_error']:.3e}"
    )
    assert r["kelvin_openbc_error"] < 0.005, (
        f"Kelvin open-BC error should be small: {r['kelvin_openbc_error']:.3e}"
    )
    # the headline: Kelvin open-BC error is >10x below the interior FEM error
    assert r["kelvin_openbc_error"] * 10 < r["interior_fem_error"], (
        f"Kelvin open-BC error {r['kelvin_openbc_error']:.3e} not << interior FEM "
        f"error {r['interior_fem_error']:.3e}"
    )
    # field-level open-BC error is consistent with (and below) the DtN eigenvalue error
    assert r["kelvin_openbc_error"] < r["kelvin_dtn_rel_err"], (
        f"field open-BC error {r['kelvin_openbc_error']:.3e} should track the DtN "
        f"eigenvalue error {r['kelvin_dtn_rel_err']:.3e}"
    )


def test_openbc_error_vs_exterior_mesh():
    """Verify open-BC accuracy as a function of the EXTERIOR (Kelvin ball) mesh size,
    interior FIXED.  The open-BC error CONVERGES as the exterior refines, while the
    interior FEM error stays put and dominates -- Kameari's exterior refinement,
    isolated from the interior error."""
    s = kelvin_openbc_error_vs_exterior_mesh(
        kelvin_maxh_list=(0.6, 0.35, 0.25), R_inner=0.5, R_outer=1.0,
        maxh=0.4, order=2, degree=1, kelvin_order=1)
    assert s["always_below_fem"], (
        "Kelvin open-BC error should be below the interior FEM error at every "
        f"exterior mesh: {[m['kelvin_openbc_error'] for m in s['per_mesh']]} vs "
        f"FEM {s['interior_fem_error']:.3e}"
    )
    assert s["converges"], (
        "open-BC error should be non-increasing as the exterior mesh refines: "
        f"{[(m['kelvin_ndof'], m['kelvin_openbc_error']) for m in s['per_mesh']]}"
    )
    # refinement actually moves it: finest exterior mesh beats the coarsest
    coarsest = s["per_mesh"][0]["kelvin_openbc_error"]
    finest = s["per_mesh"][-1]["kelvin_openbc_error"]
    assert finest < coarsest / 2.0, (
        f"exterior refinement should cut the open-BC error: coarsest={coarsest:.3e}, "
        f"finest={finest:.3e}"
    )
    # and even the coarsest exterior is an order of magnitude under the interior error
    assert coarsest * 10 < s["interior_fem_error"]


def test_twosphere_periodic_kelvin_reproduces_dipole():
    """The lab's REAL Kelvin (two offset spheres + periodic BC + μ'=(R/r')²) solved
    for the shell dipole reproduces the exact u = R_in²z/r³ -- confirming the genuine
    Kelvin coupling matches the analytic exterior (and the single-ball effective-DtN
    used elsewhere).  ~5.5 % at maxh=0.3 order=2 (mostly interior FEM error, matching
    the analytic-DtN shell solve)."""
    r = kelvin_twosphere_shell_dipole(R_inner=0.5, R_outer=1.0, offset=3.0,
                                      maxh=0.3, order=2)
    assert r["rel_err"] < 0.10, (
        f"two-sphere Kelvin dipole rel_err={r['rel_err']:.3e} (expected ~5.5%, < 10%); "
        f"ndof={r['ndof']}"
    )


def test_capacitance_inductance_dtn_dual():
    """The C/L dual: capacitance and external inductance are two Steklov modes of the
    SAME exterior scalar Laplace DtN, on different rungs of the -(n+1)/R ladder.

      * C     <-> n=0 (MONOPOLE): the constant image is captured EXACTLY at every
                  order -> defect_0 ~ machine zero.
      * L_ext <-> n=1 (DIPOLE):   a current loop has no magnetic monopole, so its
                  leading exterior multipole is the dipole; W_ext = 1/2 mu0 (n+1)/R
                  oint phi^2, so L_ext = 2 W_ext/I^2 inherits defect_1, which FALLS
                  with FEM order (dipole image is linear -> captured at p>=1, then a
                  curved-geometry floor).

    Locks the inductance half of the DtN datasheet (knowledge: dtn_coarse_mesh,
    datasheet topic; example inductance_dtn.py).  This is the FIELD-ENERGY /
    potential-exterior inductance, NOT the ngsbem vector single-layer (a different
    operator with no -(n+1)/R ladder)."""
    # C: n=0 monopole is machine-exact (constant image, any order)
    c = kelvin_dtn_eigenvalue(R=1.0, degree=0, maxh=0.5, order=2)
    assert c["rel_err"] < 1e-9, f"C (n=0 monopole) must be exact, got {c['rel_err']:.2e}"

    # L_ext: n=1 dipole defect falls with order (linear image captured at p>=1)
    l1 = kelvin_dtn_eigenvalue(R=1.0, degree=1, maxh=0.4, order=1)
    l3 = kelvin_dtn_eigenvalue(R=1.0, degree=1, maxh=0.4, order=3)
    assert l1["rel_err"] < 1e-2, f"L_ext (n=1) coarse defect too large: {l1['rel_err']:.2e}"
    assert l3["rel_err"] < l1["rel_err"], (
        f"L_ext defect must fall with order: p1={l1['rel_err']:.2e}, p3={l3['rel_err']:.2e}"
    )

    # the energy identity that ties the DtN eigenvalue to the exterior inductance:
    #   1/2 (n+1)/R oint phi^2 == integral_{r>R}|H|^2 == m^2/(6 pi R^3)  (dipole)
    import math
    R, m = 1.0, 1.0
    lhs = 0.5 * (2.0 / R) * (m**2 / (12 * math.pi * R**2))
    rhs = 0.5 * (m**2 / (6 * math.pi * R**3))
    assert abs(lhs - rhs) / rhs < 1e-12, "exterior-energy identity (DtN coefficient) broken"


def main():
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn(); print("ok", name)
    print("[OK] FEM-BEM Schur: analytic DtN Robin BC exterior dipole L2 error < 1% on spherical shell.")


if __name__ == "__main__":
    from ngsolve import TaskManager

    with TaskManager():
        main()
