"""
Phase 3b: numerical regression of W(alpha) from parallelepiped FEM sweep.

For a sheared cuboid with shear x' = x + alpha_shear * z:
   - 8 right-angle edges remain at dihedral pi/2
   - 4 non-right-angle edges have dihedrals
       alpha_acute  = pi/2 - arctan(alpha_shear)
       alpha_obtuse = pi/2 + arctan(alpha_shear)
   - 4 acute edges of length L_acute = L * sqrt(1 + alpha_shear^2)
   - similarly for obtuse (by reflection in x: same edge length L_obtuse = L_acute
     because the 4 affected edges are the 4 edges parallel to x at top and bottom)

Actually for the shear x' = x + alpha_shear * z prism construction:
   - 4 horizontal x-edges (top + bottom of left/right faces) have shifted x but same length L
   - 4 horizontal y-edges (front + back of top/bottom) stay at length L
   - 4 vertical edges have length L * sqrt(1 + alpha_shear^2)

Hmm, actually the 4 edges parallel to original x-axis at top/bottom (4 edges)
become the non-right-dihedral edges (2 acute + 2 obtuse).  Their length stays L
since shear preserves these edge directions.

Total Mellin c_1 contribution:
   c_1 = -(1/mu) [ 8 * L * W(pi/2)      (8 right-angle edges, length L)
                  + 2 * L * W(alpha_a)   (2 acute edges, length L)
                  + 2 * L * W(alpha_o) ] (2 obtuse edges, length L)

If we ASSUME W(alpha) has a symmetric structure around pi/2 to leading order
(physical symmetry: acute pinch vs obtuse open are mirror images of the
boundary layer geometry), then
   W(alpha_a) + W(alpha_o) = 2 W(pi/2) + W''(pi/2) * epsilon^2 + O(eps^4)
   where epsilon = arctan(alpha_shear)

Therefore c_1 deviation from cuboid prediction is purely O(epsilon^2):
   Delta c_1 / c_1_cube = (2/12) * (W''(pi/2)/W(pi/2)) * epsilon^2

Strategy:
   1. Run FEM at f = 1e5 (Mellin-dominated regime) for alpha_shear in 0.05..0.8
   2. Extract Y_FEM - c_0/sqrt(s) - c_2/s^1.5 = c_1(alpha) / s + residual
   3. Solve for c_1(alpha)
   4. Subtract cuboid prediction at right-angle edges
   5. Back out W(alpha_acute) + W(alpha_obtuse) sum as function of epsilon
   6. Fit quadratic in epsilon to extract W''(pi/2) -- the local curvature
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


def Y_planar_only(s, S):
    """Just c_0/sqrt(s) (planar SIBC, no edge correction)."""
    c_0 = S * math.sqrt(SIGMA / MU)
    return c_0 / cmath.sqrt(s)


def Y_with_vertex(s, S, n_verts=8):
    """c_0/sqrt(s) + c_2/s^1.5  (planar + vertex, NO edge)."""
    c_0 = S * math.sqrt(SIGMA / MU)
    c_2 = (48.0 / (math.pi * MU**1.5 * math.sqrt(SIGMA))) * (n_verts / 8.0)
    return c_0 / cmath.sqrt(s) + c_2 / s**1.5


def main():
    SetNumThreads(8)
    print("=== Phase 3b: W(alpha) numerical regression ===")
    print(f"L = {L*1e3} mm, sigma = {SIGMA:.2e}")
    print("Sweep alpha_shear = 0.05 to 0.8, extract c_1 = mu * s * (Y_FEM - c_0/sqrt(s) - c_2/s^1.5)")
    print()
    print(f"{'alpha_shear':>11}  {'eps_acute (deg)':>15}  {'Y_FEM_re':>11}  {'Y_FEM_im':>11}  {'Sum L_e W (m)':>13}  {'Avg W':>10}  {'rel.dev.from cube':>17}")

    # Use f = 1e5 where c_0/sqrt(s) dominates Mellin and c_1 is the residual we want
    f = 1e5
    s = 1j * 2 * math.pi * f
    order = 3

    # Reference: cube alpha=0 should give cuboid c_1 = -(16/(pi mu)) * (L+L+L) = -(48 L)/(pi mu)
    # which corresponds to sum_edges L_e W(alpha_e) = (12L) * W(pi/2) = (12L)(4/pi) = 48L/pi
    # So Sum_L_e_W = (4 L_total) / pi for cuboid. Cross-check.
    Sum_LW_cube_predicted = 12 * L * (4.0 / math.pi)
    print(f"  (Reference: cube Sum L_e * W(pi/2) = 12 L * 4/pi = {Sum_LW_cube_predicted*1e3:.4f} mm)")
    print()

    results = []
    alpha_sweep = [0.00, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.60, 0.80]

    for alpha_shear in alpha_sweep:
        maxh = L / 8
        try:
            mesh, shape = build_parallelepiped(alpha_shear, maxh)
        except Exception as e:
            print(f"  {alpha_shear:11.3f}  build failed: {e}")
            continue
        V, S, L_total_geom = measure_geom(shape)
        Y_DC = V * SIGMA

        with TaskManager():
            Y_FEM = solve_at(mesh, f, order, Y_DC)

        # Y_FEM ≈ c_0/sqrt(s) + c_1/s + c_2/s^1.5 + higher
        # → c_1 ≈ s * (Y_FEM - c_0/sqrt(s) - c_2/s^1.5)
        Y_subtract = Y_with_vertex(s, S, n_verts=8)
        c_1_extracted = s * (Y_FEM - Y_subtract)
        # c_1 should be approximately real negative for cuboid
        # c_1 = -(1/mu) * Sum_e L_e * W(alpha_e)
        # → Sum_e L_e * W(alpha_e) = -mu * Re(c_1)
        # (we take Re because c_1 by construction is real for the asymptote)
        Sum_LW = -MU * c_1_extracted.real
        avg_W = Sum_LW / L_total_geom  # average W per unit edge

        # Deviation from cube reference (cube has all 12 right-angle edges)
        eps_deg = math.degrees(math.atan(alpha_shear))
        rel_dev = (Sum_LW - Sum_LW_cube_predicted) / Sum_LW_cube_predicted * 100

        print(f"  {alpha_shear:11.3f}  {eps_deg:15.2f}  {Y_FEM.real:11.4e}  {Y_FEM.imag:11.4e}  {Sum_LW*1e3:12.4f}  {avg_W:10.4f}  {rel_dev:16.3f}%")

        results.append({
            "alpha_shear": alpha_shear,
            "eps_rad": math.atan(alpha_shear),
            "eps_deg": eps_deg,
            "S": S,
            "L_total_geom": L_total_geom,
            "Y_FEM": Y_FEM,
            "Sum_LW": Sum_LW,
            "avg_W": avg_W,
        })

    # Fit: Sum_LW(alpha) = Sum_LW_cube + 2 * L * (W(alpha_a) - W(pi/2)) + 2 * L * (W(alpha_o) - W(pi/2))
    # Symmetric expansion: W(pi/2 + eps) ≈ W(pi/2) + W' eps + W''/2 eps^2 + ...
    # acute = pi/2 - eps, obtuse = pi/2 + eps → W_a + W_o ≈ 2 W(pi/2) + W'' eps^2 + O(eps^4)
    # So Delta Sum_LW = 2 L * W'' * eps^2 + O(eps^4)
    print()
    print("=== Quadratic fit: Delta Sum_LW = 2 L * W''(pi/2) * eps^2 ===")
    eps_sq = np.array([r["eps_rad"]**2 for r in results])
    delta = np.array([r["Sum_LW"] - Sum_LW_cube_predicted for r in results])
    # Fit through origin: delta = K * eps^2
    K = np.sum(eps_sq * delta) / np.sum(eps_sq**2) if np.sum(eps_sq**2) > 0 else 0.0
    W_pp = K / (2 * L)
    W_at_pi_2 = 4.0 / math.pi
    print(f"  Fit:  Delta Sum_LW = {K*1e3:.4f} mm * eps^2  (eps in radians)")
    print(f"  W''(pi/2) ≈ {W_pp:.4f}")
    print(f"  W''(pi/2) / W(pi/2) ≈ {W_pp / W_at_pi_2:.4f}")
    print()
    print("  Residuals (mm):")
    for r in results:
        pred_delta = K * r["eps_rad"]**2
        actual_delta = r["Sum_LW"] - Sum_LW_cube_predicted
        print(f"    eps = {r['eps_deg']:6.2f} deg, "
              f"actual = {actual_delta*1e3:8.4f}, "
              f"pred (K*eps^2) = {pred_delta*1e3:8.4f}, "
              f"residual = {(actual_delta - pred_delta)*1e3:8.4f}")


if __name__ == "__main__":
    main()
