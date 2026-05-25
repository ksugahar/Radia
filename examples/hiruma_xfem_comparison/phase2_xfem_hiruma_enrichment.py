"""Phase 2 — Hiruma 2023 XFEM enrichment psi=exp(-gamma*xi) on cylinder.

Reference: Hiruma 2023 IEEE TMag 59(5), eq. (1)+(4):
    u_h(x) = sum_i  u_i N_i(x)              (standard part)
           + sum_{i in E} N_i(x) * psi(x) a_i     (enriched part)
    psi(x) = exp(-gamma * xi(x))
    gamma = sqrt(j*omega*sigma*mu),  xi(x) = signed distance from surface

For the cylinder of radius R: xi(x,y) = R - sqrt(x^2+y^2)  (positive inside).

This is a Partition-of-Unity FEM enrichment (NOT cut-element XFEM): the mesh
conforms to the boundary, and we enrich nodal DOFs near the surface with the
exponential skin-decay profile.

Strategy in NGSolve:
    - Compound FE space fes = fes_std x fes_enr
    - Both are H1, order=1, complex=True, Dirichlet on "outer"
    - The "selection" of enriched nodes is done by setting fes_enr DOFs in
      the bulk to zero (Dirichlet-like) — only a band near r=R is free.
      For simplicity, this prototype enriches ALL nodes; later refinement
      restricts to the surface band.
    - The trial / test function in the bilinear form:
         u_total(x) = u_std(x) + psi(x) * u_enr(x)
      with grad(u_total) = grad(u_std) + psi * grad(u_enr) + u_enr * grad(psi)
      grad(psi) = gamma * (x,y)/r * psi   (since xi = R - r => grad(xi)=-hat_r)

Output:
    results_phase2.npz : r_over_delta, P_xfem, E_xfem, ndof_xfem
                         (joined with Phase 1 data for comparison)
"""

import math
import numpy as np
import json
from pathlib import Path

from ngsolve import (
    Mesh, H1, BilinearForm, LinearForm, GridFunction, FESpace,
    grad, dx, CoefficientFunction, Integrate, x, y, sqrt as ng_sqrt,
    Conj, exp as ng_exp,
)
from netgen.geom2d import SplineGeometry


# physics (must match Phase 1)
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
# Hiruma XFEM solver on the cylinder
# =====================================================================

def xfem_solve(mesh, omega, integration_order=10):
    """Hiruma 2023 XFEM with psi=exp(-gamma*xi) for the cylinder problem.

    Solves the Galerkin form (no complex-conjugate on the test space,
    so the matrix is complex symmetric rather than skew-Hermitian as in
    Hiruma's eq. 6 — this is acceptable because we report norms, not
    matrix structure).

       a(u_t, v_t) = int_Omega [ nu grad(u_t) . grad(v_t)
                                 + j*omega*sigma * u_t * v_t ] dx
       f(v_t)      = int_Omega J0 * v_t dx

    with u_t = u_std + psi*u_enr and v_t = v_std + psi*v_enr.

    Returns: (P_avg, E_avg, ndof, total_time_diag).
    """
    # Build compound complex H1
    fes_std = H1(mesh, order=1, complex=True, dirichlet="outer")
    fes_enr = H1(mesh, order=1, complex=True, dirichlet="outer")
    fes     = fes_std * fes_enr
    (u_std, u_enr), (v_std, v_enr) = fes.TnT()

    # Enrichment function psi(x) = exp(-gamma*xi(x))
    gamma_c    = ng_sqrt(1j * omega * SIGMA * MU0)   # complex sqrt
    r_cart     = ng_sqrt(x * x + y * y + 1e-30)      # eps to avoid /0 at origin
    xi         = R - r_cart                           # distance from outer surface
    psi        = ng_exp(-gamma_c * xi)

    # grad(psi) = -gamma * grad(xi) * psi = -gamma * (-hat_r) * psi
    #           = gamma * (x,y)/r * psi
    grad_psi_x = gamma_c * x / r_cart * psi
    grad_psi_y = gamma_c * y / r_cart * psi
    grad_psi   = CoefficientFunction((grad_psi_x, grad_psi_y))

    # Build expanded grad operator:
    #   grad(u_std) standard
    #   grad(psi * u_enr) = psi * grad(u_enr) + u_enr * grad(psi)
    g_us  = grad(u_std)
    g_ue  = psi * grad(u_enr) + u_enr * grad_psi
    g_vs  = grad(v_std)
    g_ve  = psi * grad(v_enr) + v_enr * grad_psi

    # Total = standard + enriched contributions
    u_total = u_std + psi * u_enr
    v_total = v_std + psi * v_enr
    g_ut    = g_us + g_ue
    g_vt    = g_vs + g_ve

    # Bilinear form with HIGHER-order Gauss for the exponential
    a = BilinearForm(fes, symmetric=True)
    a += (NU * g_ut * g_vt + 1j * omega * SIGMA * u_total * v_total) \
         * dx(bonus_intorder=integration_order)
    a.Assemble()

    f = LinearForm(fes)
    f += (J0 * v_total) * dx(bonus_intorder=integration_order)
    f.Assemble()

    # Solve
    A = GridFunction(fes)
    A.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec

    # Reconstruct the total A_z field
    A_std_cmp, A_enr_cmp = A.components
    Az_total = A_std_cmp + psi * A_enr_cmp

    # Post-process: Joule loss + energy density area-averages
    E_z   = -1j * omega * Az_total
    P_cf  = 0.5 * SIGMA * (E_z * Conj(E_z))
    # B = curl(A_z e_z) = (dA_z/dy, -dA_z/dx, 0)
    grad_Az = grad(A_std_cmp) + psi * grad(A_enr_cmp) + A_enr_cmp * grad_psi
    B_x = grad_Az[1]
    B_y = -grad_Az[0]
    W_cf = (B_x * Conj(B_x) + B_y * Conj(B_y)) / (2 * MU0)

    area  = math.pi * R * R
    P_avg = Integrate(P_cf, mesh, order=integration_order).real / area
    W_avg = Integrate(W_cf, mesh, order=integration_order).real / area
    return P_avg, W_avg, fes.ndof, {}


# =====================================================================
# Main: re-run the sweep, this time with XFEM, and combine with Phase 1
# =====================================================================

def main():
    out_dir = Path(__file__).parent
    print(f"output dir: {out_dir}", flush=True)

    # Same coarse mesh as Phase 1
    mesh_coarse = make_disk_mesh(maxh=R / 3.5)
    print(f"coarse mesh: ne={mesh_coarse.ne}, nv={mesh_coarse.nv}", flush=True)

    # Load Phase 1 results to combine
    p1 = np.load(out_dir / "results_phase1.npz")
    r_over_delta_list = p1["r_over_delta"]

    print()
    print("=== XFEM sweep (Hiruma exp(-gamma*xi) enrichment) ===", flush=True)
    print(f"{'r/d':>6} {'omega':>12} {'P_bess':>12} {'P_xfem':>12} "
          f"{'err_xfem':>10} {'err_coarse':>10} (Phase 1)",
          flush=True)
    results = []
    for i, rd in enumerate(r_over_delta_list):
        omega = float(p1["omega"][i])
        P_bess = float(p1["P_bessel"][i])
        P_c    = float(p1["P_coarse"][i])
        try:
            P_x, E_x, n_x, _ = xfem_solve(mesh_coarse, omega,
                                           integration_order=10)
            err_x = (P_x - P_bess) / P_bess * 100
            err_c = (P_c - P_bess) / P_bess * 100
            results.append({
                "r_over_delta": float(rd), "omega": omega,
                "P_bessel": P_bess, "P_xfem": P_x, "E_xfem": E_x,
                "ndof_xfem": n_x, "err_xfem_pct": float(err_x),
                "err_coarse_pct": float(err_c),
            })
            print(f"{rd:6.2f} {omega:12.3e} {P_bess:12.4e} {P_x:12.4e} "
                  f"{err_x:+9.2f}% {err_c:+9.2f}%", flush=True)
        except Exception as e:
            print(f"{rd:6.2f} XFEM FAILED: {type(e).__name__}: {e}",
                  flush=True)
            import traceback
            traceback.print_exc()
            break

    # Save
    out_npz = out_dir / "results_phase2.npz"
    np.savez(out_npz,
             r_over_delta = np.array([r["r_over_delta"] for r in results]),
             omega        = np.array([r["omega"]        for r in results]),
             P_xfem       = np.array([r["P_xfem"]       for r in results]),
             E_xfem       = np.array([r["E_xfem"]       for r in results]),
             ndof_xfem    = np.array([r["ndof_xfem"]    for r in results]),
             err_xfem_pct = np.array([r["err_xfem_pct"] for r in results]))
    print()
    print(f"saved: {out_npz}", flush=True)
    if results:
        print(f"=== Phase 2 summary ===", flush=True)
        print(f"  XFEM DOF (all-node enrichment) = {results[0]['ndof_xfem']}",
              flush=True)
        print(f"  XFEM coarse-mesh error at r/d=15: "
              f"{results[-1]['err_xfem_pct']:+.2f}% (Hiruma claims < 3%)",
              flush=True)


if __name__ == "__main__":
    main()
