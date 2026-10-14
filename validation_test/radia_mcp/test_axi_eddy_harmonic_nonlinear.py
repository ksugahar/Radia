"""Nonlinear (saturating mu(|B|)) forced time-harmonic axisymmetric eddy currents
(solve.solve_axi_eddy_harmonic_nonlinear) -- the Picard layer unblocked by the
complex H1Henrotte field eval (src/ext/axifem complex CalcMatrix fix), since each
sweep samples mu_of_B(|B(A)|) at element centroids from the COMPLEX A_phi.

Two checks:
  (1) LINEAR LIMIT (definitive correctness): a constant mu(|B|) = mu0 must reduce
      EXACTLY to the linear solve_axi_eddy_harmonic (same K, M, RHS) -- P_eddy
      matches to machine precision.
  (2) SATURATION MONOTONICITY (physics): a magnetic conductor (mu_r=200, soft
      saturation at Bsat=1 T).  At a LOW applied amplitude |B| << Bsat so the
      volume-averaged effective <mu_r> stays near the unsaturated mu_r; at a HIGH
      amplitude |B| >~ Bsat so saturation drops <mu_r>.  <mu_r>_high < <mu_r>_low
      shows the nonlinear coupling (which needs the complex field eval) is live.
"""
import math
import os
import sys

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from ngsolve import (Mesh, CoefficientFunction, Conj, sqrt, Integrate, grad, dx,
                     TaskManager, x as r_cf)
from netgen.occ import OCCGeometry, MoveTo, Glue, X, Y
from radia_mcp.radia_ngsolve.solve import (solve_axi_eddy_harmonic,
                                           solve_axi_eddy_harmonic_nonlinear)

MU0 = 4e-7 * math.pi
R_DISK, T_DISK, SIGMA = 10e-3, 2e-3, 5.0e6      # steel-ish conductor
PI = math.pi


def build():
    air = MoveTo(0, -50e-3).Rectangle(50e-3, 100e-3).Face()
    disk = MoveTo(0, -T_DISK / 2).Rectangle(R_DISK, T_DISK).Face()
    disk.faces.name = "conductor"; air.faces.name = "air"; disk.maxh = 0.5e-3
    air.edges.Min(X).name = "axis"; air.edges.Max(X).name = "outer"
    air.edges.Min(Y).name = "outer"; air.edges.Max(Y).name = "outer"
    return Mesh(OCCGeometry(Glue([air, disk]), dim=2).GenerateMesh(maxh=5e-3))


def avg_mu_r(gfu, mesh, mu_of_B):
    """Volume(2*pi*r)-averaged effective relative permeability over the conductor,
    computed from the converged complex A_phi (uses the complex field eval)."""
    dA = grad(gfu)
    Br = -dA[1]
    Bz = dA[0] + gfu / r_cf
    Bmag = sqrt((Br * Conj(Br) + Bz * Conj(Bz)).real + 1e-30)
    mucf = mu_of_B(Bmag)
    cond = dx(mesh.Materials("conductor"))
    num = Integrate(mucf * 2 * PI * r_cf * cond, mesh)
    den = Integrate(2 * PI * r_cf * cond, mesh)
    return float((num / den).real) / MU0


def main():
    mesh = build()
    sigma_cf = mesh.MaterialCF({"conductor": SIGMA}, default=0.0)

    # ---- (1) LINEAR LIMIT: mu(|B|) = mu0 const  ==  linear solve_axi_eddy_harmonic
    omega1 = 2 * PI * 711.0
    applied1 = CoefficientFunction(0.5 * r_cf)             # uniform axial B0 = 1
    mu_const = lambda Bmag: CoefficientFunction(MU0)       # noqa: E731
    with TaskManager():
        _g_nl, P_nl = solve_axi_eddy_harmonic_nonlinear(
            mesh, mu_const, sigma_cf, omega1, applied1, order=1)
        _g_li, P_li = solve_axi_eddy_harmonic(
            mesh, CoefficientFunction(MU0), sigma_cf, omega1, applied1, order=1)
    rel = abs(P_nl - P_li) / abs(P_li)
    print(f"[1] linear limit: P_nl={P_nl:.6e}  P_lin={P_li:.6e}  rel={rel:.2e}")
    assert rel < 1e-6, f"nonlinear with const mu must equal linear solve (rel {rel:.2e})"

    # ---- (2) SATURATION: magnetic conductor, <mu_r> drops as amplitude rises ----
    MU_R, BSAT = 200.0, 1.0
    mu_r_cf = mesh.MaterialCF({"conductor": MU_R}, default=1.0)   # air stays mu_r=1

    def mu_sat(Bmag):
        # Frohlich/Kennelly soft saturation: mu_r -> 1 as |B| >> Bsat; air (mu_r_cf=1)
        # is unaffected because (mu_r_cf-1)=0 there regardless of |B|.
        mu_r_eff = 1.0 + (mu_r_cf - 1.0) / (1.0 + (Bmag / BSAT) ** 2)
        return MU0 * mu_r_eff

    # Low f so the skin depth exceeds the disk thickness and |B| penetrates the
    # bulk (at higher f the eddy skin effect shields the interior and saturation
    # stays confined to the surface, so the volume-averaged <mu_r> drops less).
    omega2 = 2 * PI * 5.0
    with TaskManager():
        g_lo, _ = solve_axi_eddy_harmonic_nonlinear(
            mesh, mu_sat, sigma_cf, omega2,
            CoefficientFunction(0.5 * 1e-3 * r_cf), order=1, relax=0.4, max_iter=80)    # B0 = 1e-3 T
        g_hi, _ = solve_axi_eddy_harmonic_nonlinear(
            mesh, mu_sat, sigma_cf, omega2,
            CoefficientFunction(0.5 * 15.0 * r_cf), order=1, relax=0.3, max_iter=150)    # B0 = 15 T
    mur_lo = avg_mu_r(g_lo, mesh, mu_sat)
    mur_hi = avg_mu_r(g_hi, mesh, mu_sat)
    print(f"[2] saturation: <mu_r> low-amp = {mur_lo:.2f}   high-amp = {mur_hi:.2f}"
          f"   (unsaturated mu_r = {MU_R:.0f})")

    assert mur_lo > 0.9 * MU_R, \
        f"low-amplitude <mu_r>={mur_lo:.1f} should stay near unsaturated {MU_R:.0f}"
    assert mur_hi < 0.85 * mur_lo, \
        f"high-amplitude <mu_r>={mur_hi:.1f} should be visibly saturated below {mur_lo:.1f}"

    print("\n[OK] solve_axi_eddy_harmonic_nonlinear validated: reduces EXACTLY to the\n"
          f"     linear harmonic solve for constant mu (rel {rel:.1e}), and a saturating\n"
          f"     mu(|B|) drops <mu_r> from {mur_lo:.0f} (|B|<<Bsat) to {mur_hi:.0f}"
          f" (|B|>~Bsat).\n     The Picard sweep samples mu_of_B(|B(A_phi)|) via the now-"
          "reliable complex\n     H1Henrotte field eval -- the C++ fix made this layer possible.")


if __name__ == "__main__":
    main()
