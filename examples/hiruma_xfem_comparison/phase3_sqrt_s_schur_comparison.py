"""Phase 3 — sqrt(s) Schur block CLN vs Hiruma XFEM vs CFEM on cylinder.

Builds the 4-way numerical comparison for paper 1 §V / paper 2 §VI:

    Method                                       DOF    Error r/d=15
    -----------------------------------------    ----   -----------
    CFEM fine                                    ~3.5k  < 0.2% (reference)
    CFEM coarse                                    44   -11%   (Phase 1)
    Hiruma XFEM                                    88   -2%    (Phase 2)
    Schur-complement augmented CLN                  5   ?      (this Phase)

Augmentation construction (paper 1 §III, "Schur-Complement Augmentation:
The sqrt(s) Block"):
    Y_R(s) = Y_CLN(s) + K_SIBC sqrt(s) / (s + d)
    where Y_CLN(s) = rational [N-1 / N] Pade at s=0
    extracted from the coarse-mesh matrices, and
    K_SIBC = perimeter * sqrt(sigma/mu) = 2 pi R sqrt(sigma/mu)
    for the 2D cylinder (per-unit-length admittance).

Y(s) here = total current per unit boundary E_z, i.e.
    Y(s) = integral_disk J_z(x,y; s) dA / E_z_boundary
For the unit-boundary-E_z setup (J0 = sigma, A_z=0 at r=R), V=1 so
Y(jw) = integral J_z dA.

For the CLN moment extraction, we use the coarse CFEM stiffness K0 and
the boundary-driven vector b, and compute Taylor moments
    a_n = sigma (-1)^n b^T K0^-1 (sigma K0^-1)^n b
following the formula in paper 1 §II.

Output:
    results_phase3.npz : Y_bessel, Y_cfem_coarse, Y_cfem_fine, Y_xfem,
                         Y_aug,  for each r/delta
    plus plot: |Y(jw)| log-log + r(f) asymptote
"""

import math
import numpy as np
import json
from pathlib import Path

from ngsolve import (
    Mesh, H1, BilinearForm, LinearForm, GridFunction,
    grad, dx, CoefficientFunction, Integrate, x, y, sqrt as ng_sqrt,
    Conj, exp as ng_exp,
)
from netgen.geom2d import SplineGeometry
from scipy.special import iv


# physics constants
R     = 0.01
SIGMA = 1.0e6
MU0   = 4.0e-7 * math.pi
NU    = 1.0 / MU0
J0    = SIGMA
PERIM = 2 * math.pi * R                       # cylinder perimeter
K_SIBC = PERIM * math.sqrt(SIGMA / MU0)        # surface-area-only scaling


def make_disk_mesh(maxh):
    geo = SplineGeometry()
    geo.AddCircle((0, 0), R, leftdomain=1, rightdomain=0, bc="outer")
    return Mesh(geo.GenerateMesh(maxh=maxh))


# =====================================================================
# Admittance extraction from each method
# =====================================================================

def bessel_Y(omega):
    """Closed-form Y(jw) for the cylinder driven by uniform J0=sigma.

    PDE:  (1/r) d/dr(r dA/dr) - kappa^2 A = -J0/nu
          where kappa = sqrt(j*omega*sigma*mu) = (1+j)/delta
    Solution (A=0 at r=R):
          A(r) = (J0/(j*omega*sigma)) * (1 - I_0(kappa r) / I_0(kappa R))
    Total current = integral of J_total = J0 - j*omega*sigma*A:
          int J_total dA = J0 * pi*R^2 - j*omega*sigma * int A dA
    int A dA = (J0/(j*omega*sigma)) * (pi*R^2 - 2*pi*R*I_1(kR)/(kappa*I_0(kR)))
    Substituting:
          Y = J0 * pi*R^2 - J0 * (pi*R^2 - 2*pi*R*I_1/(kappa*I_0))
            = J0 * 2*pi*R * I_1(kappa R) / (kappa * I_0(kappa R))
    With J0 = sigma:
          Y(jw) = sigma * 2*pi*R * I_1(kappa R) / (kappa * I_0(kappa R))

    This is the POSITIVE form (no catastrophic cancellation at low omega,
    unlike the algebraically equivalent  sigma*(pi*R^2 - 2*pi*R*I_1/(k*I_0))
    used in the buggy v1).

    Limits (correctness check):
      DC (omega -> 0, kR -> 0): I_1(kR) ~ kR/2, I_0(kR) ~ 1
        => Y -> sigma * 2*pi*R * (kR/2) / kappa = sigma * pi*R^2  (DC admittance)
      High freq (|kR| >> 1): I_1/I_0 -> 1
        => Y -> sigma * 2*pi*R / kappa = (2*pi*R) * sqrt(sigma/mu) / sqrt(j*omega)
             = K_SIBC / sqrt(j*omega)                              (SIBC asymptote)
    """
    import cmath
    kappa = cmath.sqrt(1j * omega * SIGMA * MU0)
    kR    = kappa * R
    Y     = SIGMA * 2 * math.pi * R * iv(1, kR) / (kappa * iv(0, kR))
    return Y


def cfem_Y(mesh, omega, return_diagnostics=False):
    """Plain CFEM solve, return Y(jw) = integral sigma * (-j*omega) * A_z dA."""
    fes = H1(mesh, order=1, complex=True, dirichlet="outer")
    u, v = fes.TnT()
    a = BilinearForm(fes, symmetric=True)
    a += NU * grad(u) * grad(v) * dx
    a += 1j * omega * SIGMA * u * v * dx
    a.Assemble()
    f = LinearForm(fes)
    f += J0 * v * dx
    f.Assemble()
    A = GridFunction(fes)
    A.vec.data = a.mat.Inverse(fes.FreeDofs(),
                                inverse="sparsecholesky") * f.vec
    # Y = integral of (-j*omega * sigma * A_z) dA (the eddy current piece) ...
    # Wait — need to be careful about sign convention.
    # Source J0 = sigma => drives E_z=1 at boundary, so total Y =
    #    int J_z dA = int(J0 - j w sigma A) dA = sigma*pi R^2 - jw sigma int A dA
    intA = Integrate(A, mesh)
    Y    = SIGMA * math.pi * R * R - 1j * omega * SIGMA * intA
    return Y


def xfem_Y(mesh, omega, integration_order=10):
    """Hiruma XFEM solve, return Y(jw)."""
    fes_std = H1(mesh, order=1, complex=True, dirichlet="outer")
    fes_enr = H1(mesh, order=1, complex=True, dirichlet="outer")
    fes     = fes_std * fes_enr
    (u_std, u_enr), (v_std, v_enr) = fes.TnT()

    gamma_c    = ng_sqrt(1j * omega * SIGMA * MU0)
    r_cart     = ng_sqrt(x * x + y * y + 1e-30)
    xi         = R - r_cart
    psi        = ng_exp(-gamma_c * xi)
    grad_psi_x = gamma_c * x / r_cart * psi
    grad_psi_y = gamma_c * y / r_cart * psi
    grad_psi   = CoefficientFunction((grad_psi_x, grad_psi_y))

    u_total = u_std + psi * u_enr
    v_total = v_std + psi * v_enr
    g_ut    = grad(u_std) + psi * grad(u_enr) + u_enr * grad_psi
    g_vt    = grad(v_std) + psi * grad(v_enr) + v_enr * grad_psi

    a = BilinearForm(fes, symmetric=True)
    a += (NU * g_ut * g_vt + 1j * omega * SIGMA * u_total * v_total) \
         * dx(bonus_intorder=integration_order)
    a.Assemble()
    f = LinearForm(fes)
    f += (J0 * v_total) * dx(bonus_intorder=integration_order)
    f.Assemble()
    A = GridFunction(fes)
    A.vec.data = a.mat.Inverse(fes.FreeDofs(),
                                inverse="sparsecholesky") * f.vec

    Az_total = A.components[0] + psi * A.components[1]
    intA = Integrate(Az_total, mesh, order=integration_order)
    Y    = SIGMA * math.pi * R * R - 1j * omega * SIGMA * intA
    return Y


# =====================================================================
# sqrt(s) Schur-block augmented CLN:
#   extract moments from coarse mesh, build Y_R(s) = Y_CLN(s) + K sqrt(s)/(s+d)
# =====================================================================

def augmented_cln_Y(mesh, basis_size=4):
    """Return omega -> Y_R(jw) for the Schur-complement augmented CLN.

    The augmentation is the single-DOF sqrt(s) Schur block of paper 1
    Theorem 1 ("Single-DOF asymptote-preserving augmentation").

    Volume-source-driven formulation (cylinder, uniform J0 = sigma):
        Y_FE(s) = sigma * area  -  s * sigma^2 * y_hat(s)
        y_hat(s) = b^T (K0 + s*sigma*M)^-1 b ,  b_i = integral of N_i over disk

    Taylor of y_hat(s) at s=0:
        y_hat(s) = sum_{n>=0} y_n s^n
        y_n = (-1)^n * b^T K0^-1 (sigma * M * K0^-1)^n b

    So Y_FE(s) Taylor coefficients are:
        a_0 = sigma * area
        a_n = - sigma^2 * y_{n-1}   for n >= 1

    Then Y_R(s) = Pade[N-1/N](Y_FE)(s) + K_SIBC * sqrt(s) / (s + d).

    Bug from v1: the previous version computed the port-driven formula
    a_n = sigma (-1)^n b_J0^T K0^-1 (sigma M K0^-1)^n b_J0  with
    b_J0 = J0 * b_const = sigma * b_const, which inserts an extra
    sigma^2 factor (vs the volume-source split above) and silently
    treats the DC offset as part of the rational piece, giving
    moments off by orders of magnitude and a Y_CLN with no DC offset
    at all.  Fix: split DC term out, compute moments of y_hat only.
    """
    fes = H1(mesh, order=1, dirichlet="outer")
    u, v = fes.TnT()
    a0 = BilinearForm(fes, symmetric=True)
    a0 += NU * grad(u) * grad(v) * dx
    a0.Assemble()
    M = BilinearForm(fes, symmetric=True)
    M += u * v * dx
    M.Assemble()
    # b_i = integral of N_i (no J0 multiplier — that becomes the sigma * area
    # DC term split out below)
    fL = LinearForm(fes)
    fL += 1.0 * v * dx
    fL.Assemble()

    K0_inv = a0.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")
    b = fL.vec.CreateVector(); b.data = fL.vec
    w = fL.vec.CreateVector(); w.data = K0_inv * b   # w_0 = K0^-1 b
    b_np = np.array(b.FV().NumPy())

    # y_hat moments
    y_hat = []
    y_hat.append(float(np.dot(b_np, np.array(w.FV().NumPy()))))      # y_0
    tmp = fL.vec.CreateVector()
    for n in range(1, 2 * basis_size):
        tmp.data = M.mat * w
        tmp.data *= SIGMA
        w.data = K0_inv * tmp
        y_hat.append(((-1) ** n) *
                      float(np.dot(b_np, np.array(w.FV().NumPy()))))
    y_hat = np.array(y_hat)

    # Y_FE(s) Taylor coefficients:  a_0 = sigma*area, a_n = -sigma^2 * y_{n-1}
    area = math.pi * R * R
    n_coeffs = 2 * basis_size
    a_coeffs = np.zeros(n_coeffs)
    a_coeffs[0] = SIGMA * area
    for n in range(1, n_coeffs):
        a_coeffs[n] = -SIGMA * SIGMA * y_hat[n - 1]

    # The bare a_coeffs span ~30+ decades (a_0 ~ 1e2, a_n decays geometrically),
    # so direct Pade in s is ill-conditioned (Hankel rcond ~ 1e-37).  Rescale
    # the s-axis by s_scale = |a_1/a_0| so the Taylor series of f(t) := Y_FE(t/s_scale)
    # has O(1) coefficients.  After Pade in t, undo the substitution.
    s_scale = abs(a_coeffs[1] / a_coeffs[0]) if a_coeffs[1] != 0 else 1.0
    # Taylor in t = s * s_scale:  f(t) = sum (a_n / s_scale^n) t^n
    scaled_coeffs = np.array([a_coeffs[n] / (s_scale ** n) for n in range(n_coeffs)])
    from scipy.interpolate import pade
    p_poly_t, q_poly_t = pade(scaled_coeffs, basis_size - 1, basis_size)

    # Crossover d: paper 1 §III calls d "a crossover frequency", treats it
    # as a free parameter (Theorem 1 uniqueness is independent of d).
    # The natural choice is the wall-band frequency where the rational
    # ROM transitions from bulk diffusion to surface skin behaviour.
    #
    # For the cylinder Y(s) = sigma * 2 pi R * I_1(kappa R) / (kappa I_0(kappa R))
    # the first pole is at s = -j_{0,1}^2 / (sigma mu R^2), where j_{0,1}
    # is the first zero of J_0 (Bessel function of the first kind).
    #
    # Extracting d from the Pade poles directly is unstable: the small
    # Pade poles are coarse-mesh discretisation artifacts, not the true
    # first pole.  (For basis_size=4 on this mesh, Pade finds poles at
    # -5.5e3, -1.5e4, -9.4e4; only the largest approximates the true
    # first pole at -4.6e4 within a factor ~2.)
    #
    # We use the analytical wall-band frequency for the cylinder.
    J01_SQ = 2.40483 ** 2
    d = J01_SQ / (SIGMA * MU0 * R * R)
    # Also report the Pade-extracted poles for diagnostic comparison.
    q_roots_t = np.roots(q_poly_t.coeffs[::-1])
    pade_poles_s = [float(-r.real / s_scale) for r in q_roots_t if r.real < 0]
    pade_poles_s.sort()

    def Y_R(omega):
        s = 1j * omega
        t = s * s_scale
        Y_cln = p_poly_t(t) / q_poly_t(t)
        # sqrt(s) Schur-block augmentation (paper 1 §III, Theorem 1)
        Y_schur = K_SIBC * np.sqrt(s) / (s + d)
        return Y_cln + Y_schur

    return Y_R, {
        "y_hat": y_hat.tolist(),
        "a_coeffs": a_coeffs.tolist(),
        "s_scale": float(s_scale),
        "Y_DC": float(a_coeffs[0]),
        "d_crossover": float(d),
        "d_source": "analytical wall band j_{0,1}^2 / (sigma mu R^2)",
        "pade_poles": pade_poles_s,
        "K_SIBC": float(K_SIBC),
        "basis_size": basis_size,
    }


# =====================================================================
# Main: 4-way comparison
# =====================================================================

def main():
    out_dir = Path(__file__).parent
    print(f"K_SIBC for cylinder = 2 pi R sqrt(sigma/mu) "
          f"= {K_SIBC:.6e} S/m^(1/2)", flush=True)
    print()

    mesh_coarse = make_disk_mesh(maxh=R / 3.5)
    mesh_fine   = make_disk_mesh(maxh=R / 30.0)
    print(f"coarse: nv={mesh_coarse.nv}, fine: nv={mesh_fine.nv}", flush=True)

    # Build sqrt(s)-Schur-block augmented CLN once
    print()
    print("=== Augmented CLN extraction (sqrt(s) Schur block) ===", flush=True)
    Y_aug, aug_diag = augmented_cln_Y(mesh_coarse, basis_size=4)
    print(f"  basis size = {aug_diag['basis_size']}", flush=True)
    print(f"  Y_DC (a_0) = {aug_diag['Y_DC']:.4e}  "
          f"(analytical sigma*pi*R^2 = {SIGMA*math.pi*R*R:.4e})", flush=True)
    print(f"  a_coeffs   = {aug_diag['a_coeffs']}", flush=True)
    print(f"  d crossover= {aug_diag['d_crossover']:.4e}", flush=True)
    print(f"  K_SIBC     = {aug_diag['K_SIBC']:.4e}", flush=True)

    p1 = np.load(out_dir / "results_phase1.npz")
    r_over_delta_list = p1["r_over_delta"]

    print()
    print("=== 4-way Y(jw) comparison ===", flush=True)
    print(f"{'r/d':>5} {'|Y_bess|':>12} {'|Y_fine|':>12} {'|Y_coarse|':>12} "
          f"{'|Y_xfem|':>12} {'|Y_aug|':>12}", flush=True)
    print(f"{'':>5} {'':>12} {'errFINE':>12} {'errCOARSE':>12} "
          f"{'errXFEM':>12} {'errAUG':>12}", flush=True)

    results = []
    for i, rd in enumerate(r_over_delta_list):
        omega = float(p1["omega"][i])
        Y_b = bessel_Y(omega)
        Y_c = cfem_Y(mesh_coarse, omega)
        Y_f = cfem_Y(mesh_fine,   omega)
        Y_x = xfem_Y(mesh_coarse, omega, integration_order=10)
        Y_s = Y_aug(omega)

        err_f = (abs(Y_f) - abs(Y_b)) / abs(Y_b) * 100
        err_c = (abs(Y_c) - abs(Y_b)) / abs(Y_b) * 100
        err_x = (abs(Y_x) - abs(Y_b)) / abs(Y_b) * 100
        err_s = (abs(Y_s) - abs(Y_b)) / abs(Y_b) * 100

        results.append({
            "r_over_delta": float(rd), "omega": omega,
            "Y_bessel": complex(Y_b), "Y_fine": complex(Y_f),
            "Y_coarse": complex(Y_c), "Y_xfem": complex(Y_x),
            "Y_aug": complex(Y_s),
            "err_fine_pct": float(err_f),
            "err_coarse_pct": float(err_c),
            "err_xfem_pct": float(err_x),
            "err_aug_pct": float(err_s),
        })
        print(f"{rd:5.2f} {abs(Y_b):12.4e} {abs(Y_f):12.4e} {abs(Y_c):12.4e} "
              f"{abs(Y_x):12.4e} {abs(Y_s):12.4e}", flush=True)
        print(f"{'err':>5} {'':>12} {err_f:+11.2f}% {err_c:+11.2f}% "
              f"{err_x:+11.2f}% {err_s:+11.2f}%", flush=True)

    # Save
    out_npz = out_dir / "results_phase3.npz"
    np.savez(out_npz,
             r_over_delta = np.array([r["r_over_delta"] for r in results]),
             omega        = np.array([r["omega"]        for r in results]),
             Y_bessel     = np.array([r["Y_bessel"]     for r in results]),
             Y_fine       = np.array([r["Y_fine"]       for r in results]),
             Y_coarse     = np.array([r["Y_coarse"]     for r in results]),
             Y_xfem       = np.array([r["Y_xfem"]       for r in results]),
             Y_aug     = np.array([r["Y_aug"]     for r in results]),
             err_fine_pct = np.array([r["err_fine_pct"] for r in results]),
             err_coarse_pct = np.array([r["err_coarse_pct"] for r in results]),
             err_xfem_pct = np.array([r["err_xfem_pct"] for r in results]),
             err_aug_pct = np.array([r["err_aug_pct"] for r in results]))
    print()
    print(f"saved: {out_npz}", flush=True)

    print()
    print("=== Phase 3 summary (DOF count + errors at r/d=15) ===", flush=True)
    last = results[-1]
    print(f"  Bessel:        DOF n/a   |Y| = {abs(last['Y_bessel']):.4e}",
          flush=True)
    print(f"  CFEM fine:     DOF ~3.5k err = {last['err_fine_pct']:+.3f}%",
          flush=True)
    print(f"  CFEM coarse:   DOF    44 err = {last['err_coarse_pct']:+.3f}%",
          flush=True)
    print(f"  Hiruma XFEM:   DOF    88 err = {last['err_xfem_pct']:+.3f}%",
          flush=True)
    print(f"  Augmented CLN (sqrt(s) Schur block): "
          f"DOF 5 err = {last['err_aug_pct']:+.3f}%", flush=True)


if __name__ == "__main__":
    main()
