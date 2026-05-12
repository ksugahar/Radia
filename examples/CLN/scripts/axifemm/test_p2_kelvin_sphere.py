"""Phase B3 axifemm + Kelvin transformation sphere benchmark (STUB).

STATUS (2026-05-12): Phase B3 stub. The geometry construction via OCC
boolean ("kelvin = big - small") collapses the two domains into a single
material, so this script cannot yet assign per-region mu_eff. The
SplineGeometry path used here works for 2 domains but uses parabolic
spline3 Bezier arcs (not true circles), so mesh geometry has its own
error. Full Phase B3 needs either:
  - OCC variant that preserves material names through boolean ops, or
  - per-quadrature-point mu evaluation so even constant-mu in a Kelvin
    cell is accurate (because mu varies a lot within a cell near the
    Kelvin center).
The radia.kelvin_source.kelvin_mu_factor_axisym_cf helper is wired
correctly; the bottleneck is geometry + material assignment.



Compares Cu sphere benchmark via axifemm-P2 with vs without the Kelvin
transformation. The Kelvin route should eliminate the finite-air-box
truncation error and converge to Stoll's analytical
tau_1 = mu_0 sigma R^2 / pi^2 = 738.48 us.

Setup (axisym, half-disk):
  Conductor:   r^2 + z^2 < R_sphere
  Inner air:   R_sphere <= rho <= R_bond  (small buffer to keep mesh-aspect
                                            ratio sane)
  Kelvin air:  rho <= R_bond (the Kelvin sphere; "exterior" mapped to
                              "interior" by inversion)

  Material assignment:
    conductor: mu_0, sigma_Cu
    air_phys:  mu_0, 0
    air_kelv:  mu_eff = mu_0 * (R_bond / rho')^2, 0  [varying]

Per-element centroid sampling means the Kelvin factor is constant
within each cell -- discretisation error scales with the variance of
(R_bond/rho)^2 across the cell. Refine near rho' = 0 for accuracy.

Limitations of this Phase B3 stub:
  - Centroid mu-sampling caps accuracy; for high precision the
    integrator needs per-quadrature-point mu evaluation.
  - Dirichlet at rho' = 0 (Kelvin center = physical infinity) is
    approximated by Dirichlet on a small mesh disc; alternatively
    Dirichlet on the outermost Kelvin point.

Cu sphere R=10 mm vs Stoll 738.48 us.
"""
from __future__ import annotations
import sys
sys.stdout.reconfigure(encoding="utf-8")

from math import pi
import time
import numpy as np
import scipy.sparse as sp
import scipy.linalg as sla

from netgen.geom2d import SplineGeometry
from ngsolve import (
    Mesh, BilinearForm, CoefficientFunction, TaskManager, ngsglobals,
    x, y, sqrt, IfPos,
)
from radia.radia_axifemm import (
    H1Henrotte, AxiHenrotteStiffnessBFI, AxiHenrotteSigmaMassBFI,
)

R_SPHERE = 10.0e-3
R_BOND = 30.0e-3      # bond surface where physical air meets Kelvin air
SIGMA_CU = 5.8e7
MU0 = 4 * pi * 1e-7
STOLL_TAU_1 = MU0 * SIGMA_CU * R_SPHERE ** 2 / (pi ** 2)


def build_geometry_with_kelvin():
    """Two-domain SplineGeometry: conductor half-disk + Kelvin annulus.

    Uses straight-line approximations of the half-circles (spline3 quadratic
    Bezier is NOT a true circle). For curved geometry, OCC boolean is the
    right route but the material-tagging-through-boolean is fragile -- see
    test_p2_kelvin_sphere_occ.py for the OCC variant once we can keep the
    material names through the boolean operation.
    """
    geo = SplineGeometry()
    # Conductor half-disk vertices
    p_in_bot = geo.AppendPoint(0.0, -R_SPHERE)
    p_in_apex = geo.AppendPoint(R_SPHERE, 0.0)
    p_in_top = geo.AppendPoint(0.0, R_SPHERE)
    p_out_bot = geo.AppendPoint(0.0, -R_BOND)
    p_out_apex = geo.AppendPoint(R_BOND, 0.0)
    p_out_top = geo.AppendPoint(0.0, R_BOND)

    # Inner arc (conductor boundary). leftdomain=1 (conductor), rightdomain=2 (kelvin).
    geo.Append(["spline3", p_in_bot, p_in_apex, p_in_top],
               bc="sphereBND", leftdomain=1, rightdomain=2)
    # Outer arc (Kelvin outer = physical infinity after inversion).
    geo.Append(["spline3", p_out_bot, p_out_apex, p_out_top],
               bc="kelvin_outer", leftdomain=2, rightdomain=0)
    # Axis lines (r = 0).
    geo.Append(["line", p_out_bot, p_in_bot],
               bc="axis", leftdomain=2, rightdomain=0)
    geo.Append(["line", p_in_bot, p_in_top],
               bc="axis_cond", leftdomain=1, rightdomain=0)
    geo.Append(["line", p_in_top, p_out_top],
               bc="axis", leftdomain=2, rightdomain=0)

    geo.SetMaterial(1, "conductor")
    geo.SetMaterial(2, "kelvin")
    return geo


def to_csr(mat, n):
    rs, cs, vs = mat.COO()
    return sp.coo_matrix(
        (np.array(vs), (np.array(rs), np.array(cs))),
        shape=(n, n),
    ).toarray()


def run_one(maxh, kelvin=True):
    ngsglobals.msg_level = 0
    geo = build_geometry_with_kelvin()
    mesh = Mesh(geo.GenerateMesh(maxh=maxh))
    print(f"  Materials: {mesh.GetMaterials()}, boundaries: {mesh.GetBoundaries()}")

    # Build mu CF: conductor + bonded air both get mu_0; kelvin region gets
    # mu_0 * (R_BOND / rho')^2 with mu evaluated at element centroid.
    rho = sqrt(x * x + y * y)
    rho_safe = IfPos(rho - 1e-10, rho, 1e-10)
    if kelvin:
        mu_kelvin = MU0 * (R_BOND / rho_safe) ** 2
        mu_cf = mesh.MaterialCF({
            "conductor": CoefficientFunction(MU0),
            "kelvin":    mu_kelvin,
        }, default=CoefficientFunction(MU0))
    else:
        mu_cf = CoefficientFunction(MU0)

    sigma_cf = mesh.MaterialCF({"conductor": SIGMA_CU}, default=0.0)

    fes = H1Henrotte(mesh, order=2, dirichlet="kelvin_outer")
    a = BilinearForm(fes, symmetric=True, check_unused=False)
    a += AxiHenrotteStiffnessBFI(mu_cf)
    with TaskManager(): a.Assemble()
    m_bf = BilinearForm(fes, symmetric=True, check_unused=False)
    m_bf += AxiHenrotteSigmaMassBFI(sigma_cf)
    with TaskManager(): m_bf.Assemble()

    free = np.array([i for i in range(fes.ndof) if fes.FreeDofs()[i]], dtype=int)
    K = to_csr(a.mat, fes.ndof)[np.ix_(free, free)]
    M = to_csr(m_bf.mat, fes.ndof)[np.ix_(free, free)]
    K = 0.5 * (K + K.T)
    M = 0.5 * (M + M.T)

    diag_M = np.diag(M); diag_K = np.diag(K)
    cond_mask = diag_M > 1e-30 * (np.max(diag_M) + 1e-30)
    nonzero_K = diag_K > 1e-30 * (np.max(diag_K) + 1e-30)
    air_mask = (~cond_mask) & nonzero_K
    cond_idx = np.where(cond_mask)[0]
    air_idx = np.where(air_mask)[0]

    if len(air_idx) > 0:
        K_cc = K[np.ix_(cond_idx, cond_idx)]
        K_ca = K[np.ix_(cond_idx, air_idx)]
        K_aa = K[np.ix_(air_idx, air_idx)]
        X = sla.solve(K_aa, K_ca.T, assume_a="sym")
        S = K_cc - K_ca @ X
    else:
        S = K[np.ix_(cond_idx, cond_idx)]
    M_cc = M[np.ix_(cond_idx, cond_idx)]
    S = 0.5 * (S + S.T)
    M_cc = 0.5 * (M_cc + M_cc.T)

    eigvals = sla.eigvalsh(S, M_cc, subset_by_index=[0, min(3, S.shape[0] - 1)])
    pos = eigvals[eigvals > 0]
    tau = 1.0 / pos[0] if len(pos) > 0 else np.nan
    gap = (tau - STOLL_TAU_1) / STOLL_TAU_1 * 100
    return mesh.ne, fes.ndof, tau, gap


def main():
    print("=" * 70)
    print(f"Phase B3 axifemm + Kelvin sphere benchmark vs Stoll {STOLL_TAU_1*1e6:.4f} us")
    print(f"  Cu sphere R={R_SPHERE*1e3} mm, Kelvin bond R={R_BOND*1e3} mm")
    print("=" * 70)

    levels = [5e-3, 3e-3, 2e-3]
    for maxh in levels:
        for kel in [False, True]:
            tag = "Kelvin" if kel else "noKelv"
            print(f"\nmaxh={maxh:.1e}  [{tag}]")
            try:
                ne, ndof, tau, gap = run_one(maxh, kelvin=kel)
                print(f"  ne={ne}  ndof={ndof}  tau_1={tau*1e6:.4f} us  gap={gap:+.3f}%")
            except Exception as e:
                print(f"  FAIL: {type(e).__name__}: {e}")


if __name__ == "__main__":
    main()
