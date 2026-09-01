"""
Phase 3b refined: relative-to-cube W(alpha) extraction.

The Phase 3a-attempt (04_*) hit a ~7% baseline offset at alpha=0 (cube) due
to: (i) FEM resolution at f=1e5, and (ii) imperfect c_2 + higher-order
Mellin truncation.  Either way, the BASELINE leakage is the same for all
alpha (it lives in the c_0/c_2 reconstruction), so it cancels when we look
at the DEVIATION from the cube reference:

   Delta_Sum_LW(alpha) = Sum_LW_FEM(alpha) - Sum_LW_FEM(0)

This deviation is purely the change in 4 non-right edges as their dihedrals
move away from pi/2.  The 8 right-angle edges and the 8 vertex contributions
should be identical between cube and parallelepiped (same vertex count, same
right-angle edges).  In fact for a sheared cuboid:
   - 4 horizontal y-edges (top+bot of front+back) stay at length L, dihedral pi/2
   - 4 vertical edges become length L_v = L sqrt(1 + alpha^2), dihedral pi/2
   - 4 horizontal x-edges (top+bot of left+right) stay at length L,
     but have dihedrals alpha_a = pi/2 - arctan(alpha) and alpha_o = pi/2 + arctan(alpha)
     (2 acute + 2 obtuse alternating around the body)

The 4 length-L right-angle edges + 4 length-L*sqrt(1+a^2) right-angle edges
contribute  (4 L + 4 L sqrt(1+a^2)) * W(pi/2) = (4L+4L sqrt(1+a^2)) * 4/pi

The 4 non-right edges contribute (note: all of length L, since shear is in x
which is the direction of these edges -- they shift in x but length stays L)
  2 * L * W(alpha_acute) + 2 * L * W(alpha_obtuse)

So:
  Sum_LW(alpha) = [4 L + 4 L sqrt(1+a^2) + 4 L] * 4/pi      ... wait, I need to be careful
                  + 2 L [W(alpha_a) + W(alpha_o)]

For cube (alpha=0): 12 L * 4/pi.

For parallelepiped:
  - 8 right-angle edges: 4 (y-direction, length L) + 4 (vertical, length L sqrt(1+a^2))
  - 4 non-right edges:   4 (x-direction at top/bot, length L)
                         [2 acute (pi/2-arctan a), 2 obtuse (pi/2+arctan a)]
  Sum_LW(alpha) = [4 L + 4 L sqrt(1+a^2)] * (4/pi)
                  + 2 L [W(pi/2 - arctan a) + W(pi/2 + arctan a)]

Define epsilon = arctan(alpha_shear).
W is locally symmetric around pi/2 (physical: acute and obtuse face are
boundary-layer mirror images), so:
  W(pi/2 - eps) + W(pi/2 + eps) = 2 W(pi/2) + W''(pi/2) eps^2 + O(eps^4)

Then:
  Sum_LW(alpha) - 12 L * (4/pi)
    = 4 L (sqrt(1+a^2) - 1) * (4/pi)              # vertical edge stretch
      + 2 L * W''(pi/2) * eps^2 + O(eps^4)        # dihedral perturbation

The first term is purely geometric (edge length change), known.
The second term gives W''(pi/2).

Subtraction strategy:
  Delta_geom(alpha) = 4 L (sqrt(1+a^2) - 1) * (4/pi)
  Delta_dihedral(alpha) = Sum_LW(alpha) - Sum_LW(0) - Delta_geom(alpha)
                        ≈ 2 L * W''(pi/2) * eps^2 + O(eps^4)

Solve for W''(pi/2): W'' = Delta_dihedral / (2 L * eps^2).
"""
from __future__ import annotations

import cmath
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from ngsolve import (
    Mesh, H1, BilinearForm, LinearForm, GridFunction,
    CoefficientFunction, grad, dx, Integrate, TaskManager,
    SetNumThreads
)
from netgen.occ import OCCGeometry, WorkPlane, Vec, Prism


SIGMA = 5.8e7
MU = 4 * math.pi * 1e-7
MS = MU * SIGMA
L = 5e-3


def build_parallelepiped(alpha_shear, maxh):
    bot = WorkPlane().Rectangle(L, L).Face()
    shape = Prism(bot, Vec(alpha_shear * L, 0, L))
    shape.mat("body")
    for f in shape.faces:
        f.name = "outer"
    geo = OCCGeometry(shape)
    return Mesh(geo.GenerateMesh(maxh=maxh)), shape


def solve_at(mesh, f, order, Y_DC):
    s = 1j * 2 * math.pi * f
    sMS = s * MS
    fes = H1(mesh, order=order, complex=True, dirichlet="outer")
    u, v = fes.TnT()
    a = BilinearForm(fes, symmetric=True)
    a += grad(u) * grad(v) * dx
    a += sMS * u * v * dx
    F = LinearForm(fes)
    F += -sMS * v * dx
    a.Assemble()
    F.Assemble()
    gfu = GridFunction(fes)
    gfu.vec.data = a.mat.Inverse(freedofs=fes.FreeDofs(), inverse="sparsecholesky") * F.vec
    u_avg = Integrate(gfu, mesh) / Integrate(CoefficientFunction(1), mesh)
    return Y_DC * (1 + u_avg)


def measure_geom(shape):
    V = float(shape.mass)
    S = 0.0
    for f in shape.faces:
        S += float(f.mass)
    edge_data = {}
    for e in shape.edges:
        c = (e.center[0], e.center[1], e.center[2])
        key = (round(c[0]*1e6), round(c[1]*1e6), round(c[2]*1e6))
        if key not in edge_data:
            edge_data[key] = float(e.mass)
    L_total = sum(edge_data.values())
    return V, S, L_total


def Y_baseline(s, S, n_verts=8):
    """c_0/sqrt(s) + c_2/s^1.5  (planar + vertex, NO edge correction)."""
    c_0 = S * math.sqrt(SIGMA / MU)
    c_2 = (48.0 / (math.pi * MU**1.5 * math.sqrt(SIGMA))) * (n_verts / 8.0)
    return c_0 / cmath.sqrt(s) + c_2 / s**1.5


def main():
    SetNumThreads(8)
    print("=== Phase 3b refined: W(alpha) relative to cube ===")
    print(f"L = {L*1e3} mm, sigma = {SIGMA:.2e}")
    print("Strategy: extract Sum_LW(alpha) - Sum_LW(0) - Delta_geom(alpha)")
    print("          = 2 L * W''(pi/2) * eps^2 + O(eps^4)")
    print()

    f = 1e5
    s = 1j * 2 * math.pi * f
    order = 3
    maxh = L / 8
    W_pi2 = 4.0 / math.pi

    alpha_sweep = [0.00, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.60, 0.80]
    results = []

    for alpha_shear in alpha_sweep:
        mesh, shape = build_parallelepiped(alpha_shear, maxh)
        V, S, L_total_geom = measure_geom(shape)
        Y_DC = V * SIGMA

        with TaskManager():
            Y_FEM = solve_at(mesh, f, order, Y_DC)

        Y_b = Y_baseline(s, S, n_verts=8)
        c_1_extracted = s * (Y_FEM - Y_b)
        Sum_LW = -MU * c_1_extracted.real

        eps = math.atan(alpha_shear)
        L_v = L * math.sqrt(1 + alpha_shear**2)  # vertical edge length
        Delta_geom = 4 * (L_v - L) * W_pi2  # 4 vertical edges stretched

        results.append({
            "alpha_shear": alpha_shear,
            "eps_rad": eps,
            "eps_deg": math.degrees(eps),
            "Sum_LW": Sum_LW,
            "Delta_geom": Delta_geom,
            "S": S,
            "L_total_geom": L_total_geom,
        })

    Sum_LW_cube = results[0]["Sum_LW"]
    print(f"  Baseline (alpha=0): Sum_LW_cube_FEM = {Sum_LW_cube*1e3:.4f} mm")
    print(f"  Theoretical cube Sum_LW = 12 L * 4/pi = {12*L*W_pi2*1e3:.4f} mm")
    print(f"  Baseline offset (FEM - theory) = {(Sum_LW_cube - 12*L*W_pi2)*1e3:.4f} mm")
    print()

    print(f"  {'a_shr':>6}  {'eps(deg)':>8}  {'L_v(mm)':>8}  {'D_geom(mm)':>10}  {'Sum_LW(mm)':>11}  {'Delta_dih(mm)':>13}  {'W_pp_local':>10}")
    for r in results:
        Delta_dih = r["Sum_LW"] - Sum_LW_cube - r["Delta_geom"]
        # W_pp = Delta_dih / (2 L * eps^2)
        if r["eps_rad"] > 1e-6:
            W_pp_local = Delta_dih / (2 * L * r["eps_rad"]**2)
        else:
            W_pp_local = float("nan")
        L_v_mm = L * math.sqrt(1 + r["alpha_shear"]**2) * 1e3
        print(f"  {r['alpha_shear']:6.2f}  {r['eps_deg']:8.2f}  {L_v_mm:8.4f}  "
              f"{r['Delta_geom']*1e3:10.4f}  {r['Sum_LW']*1e3:10.4f}  "
              f"{Delta_dih*1e3:12.4f}  {W_pp_local:10.4f}")
        r["Delta_dih"] = Delta_dih

    # Fit Delta_dih = K * eps^2 (through origin, weighted by small eps for cleaner fit)
    eps_sq = np.array([r["eps_rad"]**2 for r in results if r["eps_rad"] > 1e-6])
    dd = np.array([r["Delta_dih"] for r in results if r["eps_rad"] > 1e-6])
    if len(eps_sq) > 0 and np.sum(eps_sq**2) > 0:
        K = np.sum(eps_sq * dd) / np.sum(eps_sq**2)
        W_pp_fit = K / (2 * L)
        print()
        print(f"  Global fit:  Delta_dih = {K*1e3:.4f} mm * eps^2")
        print(f"  W''(pi/2)_global = {W_pp_fit:.4f}")
        print(f"  W''(pi/2) / W(pi/2) ratio = {W_pp_fit / W_pi2:.4f}")
        print()
        print("  Residual analysis:")
        for r in results:
            if r["eps_rad"] > 1e-6:
                pred = K * r["eps_rad"]**2
                actual = r["Delta_dih"]
                print(f"    eps = {r['eps_deg']:6.2f} deg ({r['eps_rad']:.4f} rad), "
                      f"actual = {actual*1e3:8.4f} mm, "
                      f"pred = {pred*1e3:8.4f} mm, "
                      f"residual = {(actual - pred)*1e3:8.4f} mm")

    # Optional: also fit cubic to check for asymmetric (eps^3) correction
    print()
    print("  Cubic-fit check:  Delta_dih = K2 eps^2 + K3 eps^3 + K4 eps^4 (asymmetry probe)")
    eps_arr = np.array([r["eps_rad"] for r in results if r["eps_rad"] > 1e-6])
    dd_arr = np.array([r["Delta_dih"] for r in results if r["eps_rad"] > 1e-6])
    A = np.column_stack([eps_arr**2, eps_arr**3, eps_arr**4])
    coeffs, *_ = np.linalg.lstsq(A, dd_arr, rcond=None)
    K2, K3, K4 = coeffs
    print(f"    K2 = {K2*1e3:8.4f} mm,  K3 = {K3*1e3:8.4f} mm,  K4 = {K4*1e3:8.4f} mm")
    print(f"    W''(pi/2) (from K2) = {K2/(2*L):.4f}")
    print(f"    K3/K2 ratio (should be ~0 if W is symmetric around pi/2) = {K3/K2:.4f}")


if __name__ == "__main__":
    main()
