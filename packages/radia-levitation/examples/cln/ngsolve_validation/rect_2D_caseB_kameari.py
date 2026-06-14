"""2D rect Case B (uniform H_0 applied) — INDEPENDENT NGSolve implementation.

Different from rect_2D_caseB_via_duality.py which post-processes Case A:
this script computes Case B moments DIRECTLY via recursive Poisson solves
on H1 scalar.

Mathematical setup:
  Case A: -nabla^2 A_z + s mu sigma A_z = mu sigma E_z,  A_z|_d = 0
  Case B: -nabla^2 H_z + s mu sigma H_z = 0,            H_z|_d = H_0

  Substitute h = H_z - H_0:  -nabla^2 h + s mu sigma h = -s mu sigma H_0,  h|_d = 0
  Expand h(s) = sum_{k>=1} h_k s^k:
    s^1: -nabla^2 h_1 = -mu sigma H_0
    s^k: -nabla^2 h_k = -mu sigma h_{k-1}     (k >= 2)

  Set H_0 = 1 for normalization. Define phi_k via Case A's recursion
  (phi_0 satisfies -nabla^2 phi_0 = mu sigma, phi_k satisfies
   -nabla^2 phi_k = -mu sigma phi_{k-1}). Then h_{k+1} = -phi_k.

  Y_B(s) = (1 - g(s))/s where g(s) = average(H_z)/H_0.
  Y_B(s) = sum_k (-s)^k M_k^B

  Direct formula: M_k^B = (-1)^k <phi_k> / 1 (since H_0=1, ab cancels)
                       = (-1)^k * (1/(ab)) integral phi_k dxdy

Workflow:
  1. Solve recursive Poisson for phi_0, phi_1, ..., phi_K
  2. Compute M_k^B = (-1)^k * <phi_k> for k=0..K
  3. Extract Cauer-I from M_k^B sequence (mpmath dps=50)
  4. Compare to Mathematica reference (rect_2D_caseB_moments.wls)

This is ALGORITHMICALLY INDEPENDENT from the Case A Tanimoto-Kameari iteration
(rect_2D_kameari.py), serving as a cross-check at the FE level.
"""
import json
from math import pi
from pathlib import Path

from netgen.geom2d import SplineGeometry
from ngsolve import (
    Mesh, H1, BilinearForm, LinearForm, GridFunction,
    CoefficientFunction, grad, dx,
    Integrate, TaskManager, ngsglobals,
)
from mpmath import mp, mpf

mp.dps = 50

mu0 = 4 * pi * 1e-7
sigma = 5.8e7
ax, ay = 5e-3, 2e-3
H_MESH = 0.05e-3
ORDER = 3
N_MOMENTS = 14   # extract Cauer-I up to ~12 stages from 14 moments


def build_rect():
    geo = SplineGeometry()
    geo.AddRectangle((0, 0), (ax, ay), bcs=["bot", "right", "top", "left"])
    return geo


def caseB_moments_via_recursive_poisson(mesh, n_moments):
    """Solve recursive Poisson and compute M_k^B for k = 0..n_moments-1."""
    fes = H1(mesh, order=ORDER, dirichlet="bot|right|top|left")
    print(f"  H1 ndof = {fes.ndof}")

    u, v = fes.TnT()
    a = BilinearForm(fes)
    a += grad(u) * grad(v) * dx

    with TaskManager():
        a.Assemble()
        inv = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")

    inv_volume = 1.0 / (ax * ay)
    mu_sigma = mu0 * sigma

    # Stage 0: -nabla^2 phi_0 = mu sigma, phi_0|_d = 0
    f = LinearForm(fes)
    src_const = CoefficientFunction(mu_sigma)
    f += src_const * v * dx
    with TaskManager():
        f.Assemble()
    gfphi = GridFunction(fes, name="phi_0")
    with TaskManager():
        gfphi.vec.data = inv * f.vec

    M_B = []
    int_phi = float(Integrate(gfphi * dx, mesh))
    avg_phi = int_phi * inv_volume
    M_B_0 = avg_phi   # M_0^B = (-1)^0 * <phi_0>
    M_B.append(M_B_0)
    print(f"[k=0] <phi_0> = {avg_phi:.10e},  M_0^B = {M_B_0:.10e}")

    # Stage k = 1..n_moments-1: -nabla^2 phi_k = -mu sigma phi_{k-1}
    phi_prev = gfphi
    for k in range(1, n_moments):
        f_k = LinearForm(fes)
        # Source = -mu sigma * phi_{k-1}
        f_k += -mu_sigma * phi_prev * v * dx
        with TaskManager():
            f_k.Assemble()
        gfphi_k = GridFunction(fes, name=f"phi_{k}")
        with TaskManager():
            gfphi_k.vec.data = inv * f_k.vec
        int_phi = float(Integrate(gfphi_k * dx, mesh))
        avg_phi = int_phi * inv_volume
        M_B_k = ((-1) ** k) * avg_phi
        M_B.append(M_B_k)
        print(f"[k={k}] <phi_{k}> = {avg_phi:.10e},  M_{k}^B = {M_B_k:.10e}")
        phi_prev = gfphi_k

    return M_B


def moments_to_cauer(M_list, n_stages):
    """Cauer-I extraction via polynomial division (mpmath)."""
    K = len(M_list)
    Ypoly = [(-1)**k * mpf(M_list[k]) for k in range(K)]
    znum = [mpf(1)]
    zden = list(Ypoly)
    R_list = []
    L_list = []
    EPS = mpf(10)**(-mp.dps + 10)

    for kk in range(n_stages):
        while len(znum) > 1 and abs(znum[-1]) < EPS:
            znum.pop()
        while len(zden) > 1 and abs(zden[-1]) < EPS:
            zden.pop()
        if not zden or abs(zden[0]) < EPS:
            break
        R_k = znum[0] / zden[0]
        R_list.append(float(R_k))
        Lp = max(len(znum), len(zden))
        presid = []
        for i in range(Lp):
            zn = znum[i] if i < len(znum) else mpf(0)
            zd = zden[i] if i < len(zden) else mpf(0)
            presid.append(zn - R_k * zd)
        pas = presid[1:]
        if not pas or abs(pas[0]) < EPS:
            break
        L_k = pas[0] / zden[0]
        L_list.append(float(L_k))
        Ln = max(len(zden), len(pas))
        nresid = []
        for i in range(Ln):
            zd = zden[i] if i < len(zden) else mpf(0)
            ps = pas[i] if i < len(pas) else mpf(0)
            nresid.append(zd * L_k - ps)
        nresid_shifted = nresid[1:]
        znum_new = [L_k * c for c in pas]
        zden_new = nresid_shifted
        znum, zden = znum_new, zden_new

    return R_list, L_list


def main():
    ngsglobals.msg_level = 0
    print(f"2D rect bar {ax*1000}x{ay*1000} mm Cu — Case B INDEPENDENT")
    print(f"  Method: recursive Poisson moment computation -> Cauer-I")
    print(f"  h={H_MESH*1000} mm, order={ORDER}, n_moments={N_MOMENTS}")

    geo = build_rect()
    mesh = Mesh(geo.GenerateMesh(maxh=H_MESH))
    print(f"  ne={mesh.ne}, nv={mesh.nv}, nedge={mesh.nedge}")
    print()

    M_B = caseB_moments_via_recursive_poisson(mesh, N_MOMENTS)
    print()
    print(f"M_k^B (k=0..{N_MOMENTS-1}):")
    for k, m in enumerate(M_B):
        print(f"  M_{k}^B = {m:.16e}")
    print()

    R_B, L_B = moments_to_cauer(M_B, 12)
    print(f"Case B Cauer-I from independent NGSolve ({len(R_B)} stages):")
    for n in range(min(len(R_B), len(L_B))):
        print(f"  n={n}: R_B={R_B[n]:.10e}, L_B={L_B[n]:.10e}")

    out = {
        "method": "Case B independent NGSolve via recursive Poisson moments + Cauer-I",
        "geometry_mm": [ax*1000, ay*1000],
        "h_mm": H_MESH*1000, "order": ORDER, "ne": mesh.ne,
        "n_moments": N_MOMENTS,
        "M_B": M_B,
        "R_B": R_B, "L_B": L_B,
    }
    out_path = Path(__file__).parent / "rect_2D_caseB_kameari_results.json"
    out_path.write_text(json.dumps(out, indent=2))
    print(f"\nResults: {out_path}")


if __name__ == "__main__":
    main()
