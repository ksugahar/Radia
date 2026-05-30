"""Phase 4 — XFEM-CLN Hankel conditioning experiment.

Hypothesis (user 2026-05-25):
  XFEM enrichment with psi(x) = exp(-xi/delta_ref) absorbs the high-k
  skin modes into the enriched DOFs.  The Krylov-from-K_0 sequence at
  s=0 then only samples the "bulk" subspace whose tau_k spread is
  narrower than the full spectrum.  This SHOULD improve the Hankel
  matrix conditioning (kappa(H_N)) and EXTEND the float64 high-N
  reliability of canonical CLN extraction.

Hypothesis is FALSIFIABLE: compute Taylor moments
  mu_n = b^T (sigma K_0^{-1} M)^n K_0^{-1} b
for both CFEM (no enrichment) and XFEM (one ref-skin enrichment) on
the Hiruma 2023 cylinder.  Form Hankel matrices H_N, compute
condition numbers, predict float64 breakdown stage.

Expected outcomes:
  A. kappa(H_N) XFEM << kappa(H_N) CFEM at all N
     -> XFEM = MOR-stage stability enabler (NEW finding).
  B. kappa(H_N) XFEM ~ kappa(H_N) CFEM at all N
     -> Hypothesis rejected, XFEM utility limited to skin-resolution.
  C. Mixed (XFEM kappa slightly better)
     -> XFEM provides partial stage-stability extension.

Output:
  phase4_xfem_hankel_conditioning.{pdf,png}
  phase4_xfem_hankel_results.json
"""

import math
import json
import numpy as np
from pathlib import Path

from ngsolve import (
    Mesh, H1, BilinearForm, LinearForm, GridFunction, FESpace,
    grad, dx, CoefficientFunction, Integrate, x, y, sqrt as ng_sqrt,
    exp as ng_exp,
)
from ngsolve import TaskManager
from netgen.geom2d import SplineGeometry

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
import matplotlib.font_manager as fm
for _p in (r"C:\Windows\Fonts\times.ttf", r"C:\Windows\Fonts\timesbd.ttf"):
    try:
        fm.fontManager.addfont(_p)
    except Exception:
        pass
matplotlib.rcParams["font.family"] = "serif"
matplotlib.rcParams["font.serif"] = ["Times New Roman", "DejaVu Serif"]
import matplotlib.pyplot as plt


# Cylinder parameters (Hiruma 2023 §III.A)
R     = 0.01
SIGMA = 1.0e6
MU0   = 4.0e-7 * math.pi
NU    = 1.0 / MU0
J0    = SIGMA


def make_disk_mesh(maxh):
    geo = SplineGeometry()
    geo.AddCircle((0, 0), R, leftdomain=1, rightdomain=0, bc="outer")
    return Mesh(geo.GenerateMesh(maxh=maxh))


# =====================================================================
# Assemble static K_0 and sigmaM matrices (CFEM and XFEM variants)
# =====================================================================

def assemble_cfem(mesh, integration_order=4):
    """Build static (real-valued, frequency-independent) CFEM matrices
    for the volume-source diffusion problem on the Hiruma cylinder.

    Returns: (K0_mat, sigmaM_mat, b_vec, fes, ndof_free)
    """
    fes = H1(mesh, order=1, dirichlet="outer")  # REAL H1 for static moments
    u, v = fes.TnT()

    a_K0 = BilinearForm(fes, symmetric=True)
    a_K0 += NU * grad(u) * grad(v) * dx
    with TaskManager():
        a_K0.Assemble()

        a_M = BilinearForm(fes, symmetric=True)
        a_M += SIGMA * u * v * dx
        a_M.Assemble()

        f = LinearForm(fes)
        f += J0 * v * dx
        f.Assemble()

        return a_K0, a_M, f, fes


def assemble_xfem(mesh, delta_ref, integration_order=10):
    """Build static XFEM matrices with REAL enrichment
    psi(x) = exp(-xi/delta_ref) at fixed delta_ref (= skin depth at
    wall band center).  Frequency-independent.

    Returns: (K0_mat, sigmaM_mat, b_vec, fes, ndof_free)
    """
    fes_std = H1(mesh, order=1, dirichlet="outer")
    fes_enr = H1(mesh, order=1, dirichlet="outer")
    fes     = fes_std * fes_enr
    (u_std, u_enr), (v_std, v_enr) = fes.TnT()

    # Real-valued enrichment psi(x) = exp(-xi/delta_ref)
    r_cart     = ng_sqrt(x * x + y * y + 1e-30)
    xi         = R - r_cart
    psi        = ng_exp(-xi / delta_ref)
    # grad(psi) = -grad(xi)/delta_ref * psi = (1/delta_ref)(x,y)/r * psi
    grad_psi_x = x / (r_cart * delta_ref) * psi
    grad_psi_y = y / (r_cart * delta_ref) * psi
    grad_psi   = CoefficientFunction((grad_psi_x, grad_psi_y))

    g_us = grad(u_std)
    g_ue = psi * grad(u_enr) + u_enr * grad_psi
    g_vs = grad(v_std)
    g_ve = psi * grad(v_enr) + v_enr * grad_psi

    u_total = u_std + psi * u_enr
    v_total = v_std + psi * v_enr
    g_ut    = g_us + g_ue
    g_vt    = g_vs + g_ve

    a_K0 = BilinearForm(fes, symmetric=True)
    a_K0 += NU * g_ut * g_vt * dx(bonus_intorder=integration_order)
    a_K0.Assemble()

    a_M = BilinearForm(fes, symmetric=True)
    a_M += SIGMA * u_total * v_total * dx(bonus_intorder=integration_order)
    a_M.Assemble()

    f = LinearForm(fes)
    f += J0 * v_total * dx(bonus_intorder=integration_order)
    f.Assemble()

    return a_K0, a_M, f, fes


# =====================================================================
# Krylov moment extraction
# =====================================================================

def krylov_moments(a_K0, a_M, f, fes, N_moments=20):
    """Compute mu_n = b^T (K_0^{-1} sigmaM)^n K_0^{-1} b for n=0..N_moments-1.
    """
    K0_inv = a_K0.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")
    # v_0 = K_0^{-1} b
    v = GridFunction(fes)
    v.vec.data = K0_inv * f.vec
    # mu_0 = b^T v_0 = b^T K_0^{-1} b
    moments = []
    moments.append(float(f.vec.InnerProduct(v.vec)))
    # Recurrence: v_{n+1} = K_0^{-1} (sigmaM v_n)
    tmp = v.vec.CreateVector()
    for n in range(1, N_moments):
        tmp.data = a_M.mat * v.vec      # tmp = sigmaM v
        v.vec.data = K0_inv * tmp        # v   = K_0^{-1} tmp
        moments.append(float(f.vec.InnerProduct(v.vec)))
    return np.array(moments)


def hankel_matrices(moments):
    """Return list of (N, kappa) where Hankel H_N[i,j] = moments[i+j],
    i,j = 0..N-1, evaluated for N = 1, 2, ..., max possible."""
    n_mom = len(moments)
    results = []
    for N in range(1, n_mom // 2 + 1):
        H = np.zeros((N, N))
        for i in range(N):
            for j in range(N):
                H[i, j] = moments[i + j]
        # Condition number via SVD
        try:
            sv = np.linalg.svd(H, compute_uv=False)
            kappa = sv[0] / sv[-1] if sv[-1] > 0 else float('inf')
        except np.linalg.LinAlgError:
            kappa = float('inf')
        results.append((N, kappa, sv[0], sv[-1]))
    return results


# =====================================================================
# Main: run both CFEM and XFEM, compute kappa, plot
# =====================================================================

def main():
    out_dir = Path(__file__).parent
    print(f"output dir: {out_dir}", flush=True)

    # Use the coarse Hiruma mesh (44 DOF reference, matches Phase 1-3b)
    mesh = make_disk_mesh(maxh=R / 3.5)
    print(f"mesh vertices: {mesh.nv}", flush=True)

    # Wall band reference: omega_wall = 1/tau_1 = x_1^2 / (mu sigma R^2)
    # For first Bessel zero x_1 = 2.4048
    x1 = 2.4048255576957728
    tau_1 = MU0 * SIGMA * R**2 / x1**2
    omega_wall = 1.0 / tau_1
    delta_ref = math.sqrt(2.0 / (omega_wall * MU0 * SIGMA))
    print(f"  omega_wall = {omega_wall:.3e} rad/s (f = {omega_wall/(2*math.pi):.2e} Hz)")
    print(f"  delta_ref  = {delta_ref:.3e} m  (R/delta = {R/delta_ref:.3f})", flush=True)

    # CFEM
    print("\n=== CFEM (no enrichment) ===", flush=True)
    a_K0_c, a_M_c, f_c, fes_c = assemble_cfem(mesh)
    ndof_c = sum(fes_c.FreeDofs())
    print(f"  free DOFs: {ndof_c}", flush=True)
    mu_c = krylov_moments(a_K0_c, a_M_c, f_c, fes_c, N_moments=20)
    print(f"  first 10 moments magnitudes (log10):")
    for n in range(min(10, len(mu_c))):
        print(f"    mu_{n} = {mu_c[n]:+.4e}  (|.|: {abs(mu_c[n]):.2e})")
    hank_c = hankel_matrices(mu_c)

    # XFEM with real ref-skin enrichment
    print(f"\n=== XFEM (psi=exp(-xi/delta_ref), delta_ref={delta_ref:.2e}) ===",
          flush=True)
    a_K0_x, a_M_x, f_x, fes_x = assemble_xfem(mesh, delta_ref)
    ndof_x = sum(fes_x.FreeDofs())
    print(f"  free DOFs: {ndof_x} (CFEM: {ndof_c}, XFEM extra: {ndof_x - ndof_c})",
          flush=True)
    mu_x = krylov_moments(a_K0_x, a_M_x, f_x, fes_x, N_moments=20)
    print(f"  first 10 moments magnitudes (log10):")
    for n in range(min(10, len(mu_x))):
        print(f"    mu_{n} = {mu_x[n]:+.4e}  (|.|: {abs(mu_x[n]):.2e})")
    hank_x = hankel_matrices(mu_x)

    # Print Hankel condition number table
    print(f"\n=== Hankel condition number kappa(H_N) ===")
    print(f"{'N':>3}  {'CFEM kappa':>14}  {'XFEM kappa':>14}  {'ratio':>10}")
    for (N_c, k_c, _, _), (N_x, k_x, _, _) in zip(hank_c, hank_x):
        assert N_c == N_x
        ratio = k_c / k_x if k_x > 0 and k_x < float('inf') else float('nan')
        print(f"{N_c:>3}  {k_c:14.3e}  {k_x:14.3e}  {ratio:10.2f}")

    # Predict float64 breakdown stage
    FLOAT64_RCOND = 1e15   # rcond inverse threshold (sigma_max/sigma_min ~ 1/eps)
    N_break_c = next((N for (N, k, _, _) in hank_c if k > FLOAT64_RCOND), None)
    N_break_x = next((N for (N, k, _, _) in hank_x if k > FLOAT64_RCOND), None)
    print(f"\nFloat64 breakdown (kappa > {FLOAT64_RCOND:.0e}):")
    print(f"  CFEM: N_break = {N_break_c}")
    print(f"  XFEM: N_break = {N_break_x}")

    # Save JSON
    results = {
        "mesh_nv": mesh.nv,
        "ndof_cfem": ndof_c,
        "ndof_xfem": ndof_x,
        "omega_wall": omega_wall,
        "delta_ref": delta_ref,
        "Rdelta_ref": R / delta_ref,
        "moments_cfem": mu_c.tolist(),
        "moments_xfem": mu_x.tolist(),
        "hankel_cfem": [{"N": N, "kappa": k, "sigma_max": sm, "sigma_min": sn}
                        for (N, k, sm, sn) in hank_c],
        "hankel_xfem": [{"N": N, "kappa": k, "sigma_max": sm, "sigma_min": sn}
                        for (N, k, sm, sn) in hank_x],
        "N_break_cfem": N_break_c,
        "N_break_xfem": N_break_x,
        "float64_threshold": FLOAT64_RCOND,
    }
    json_path = out_dir / "phase4_xfem_hankel_results.json"
    with open(json_path, "w") as fp:
        json.dump(results, fp, indent=2)
    print(f"\nWrote {json_path}")

    # Plot kappa vs N
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.0, 2.8))

    Ns_c = [N for (N, _, _, _) in hank_c]
    ks_c = [k for (_, k, _, _) in hank_c]
    Ns_x = [N for (N, _, _, _) in hank_x]
    ks_x = [k for (_, k, _, _) in hank_x]

    ax1.semilogy(Ns_c, ks_c, color="royalblue", marker="o", lw=1.5,
                 markersize=5, label=f"CFEM (no enrich, {ndof_c} DOF)")
    ax1.semilogy(Ns_x, ks_x, color="red", marker="s", lw=1.5, markersize=5,
                 label=f"XFEM real-skin (psi=exp(-xi/d_ref), {ndof_x} DOF)")
    ax1.axhline(FLOAT64_RCOND, color="gray", linestyle="--", lw=0.8,
                label=f"float64 limit (κ ~ {FLOAT64_RCOND:.0e})")
    ax1.set_xlabel(r"Hankel size $N$")
    ax1.set_ylabel(r"Condition number $\kappa(H_N)$")
    ax1.set_title("(a) Hankel conditioning vs stage")
    ax1.grid(True, which="both", linestyle=":", alpha=0.35)
    ax1.legend(loc="lower right", fontsize=7, frameon=False)

    # Moment magnitudes
    ns = np.arange(len(mu_c))
    ax2.semilogy(ns, np.abs(mu_c), color="royalblue", marker="o", lw=1.0,
                 markersize=4, label="CFEM |moments|")
    ax2.semilogy(ns[:len(mu_x)], np.abs(mu_x), color="red", marker="s",
                 lw=1.0, markersize=4, label="XFEM |moments|")
    ax2.set_xlabel(r"Moment index $n$")
    ax2.set_ylabel(r"$|\mu_n| = |b^\top (K_0^{-1}\sigma M)^n K_0^{-1} b|$")
    ax2.set_title("(b) Krylov moment magnitudes")
    ax2.grid(True, which="both", linestyle=":", alpha=0.35)
    ax2.legend(loc="lower left", fontsize=7, frameon=False)

    fig.tight_layout(pad=0.4)
    out_base = out_dir / "phase4_xfem_hankel_conditioning"
    fig.savefig(str(out_base) + ".pdf")
    fig.savefig(str(out_base) + ".png", dpi=220)
    print(f"Wrote {out_base}.pdf, {out_base}.png")


if __name__ == "__main__":
    main()
