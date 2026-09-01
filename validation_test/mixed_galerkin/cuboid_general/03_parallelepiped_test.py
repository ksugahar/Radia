"""
Phase 2 experiment: true parallelepiped (sheared cuboid) where
dihedral angles depart from pi/2.

Setup: bottom face = L x L square, extrude in direction (alpha*L, 0, L).
Result: oblique prism with
   - 2 horizontal squares (top, bottom)
   - 2 vertical rectangles parallel to y-axis (front, back)
   - 2 tilted parallelograms (left, right)
   - 12 edges:
       4 parallel to y (right-angle dihedrals)
       4 vertical (right-angle dihedrals)
       4 parallel to x at top/bottom corners (NON-right dihedrals)

For alpha != 0, the 4 non-right edges have interior dihedrals:
   - 2 acute   (where body "pinches")
   - 2 obtuse  (where body "opens")

Naive Mellin: pretend all edges are right-angle:
   c_0 = S_total * sqrt(sigma/mu)              (S = total surface area, measured)
   c_1 = -(16/(pi mu)) * L_eff
        where L_eff = sum of edge lengths assuming right-angle W
   c_2 = +48 / (pi mu^1.5 sqrt(sigma))         (8 vertices, treating each as right octant)

If FEM matches naive Mellin within ~few % at high f, the framework
extrapolates to non-orthogonal hexahedra without requiring W(beta)
per edge.  If FEM disagrees by 10%+, we need wedge-specific corrections.

Test alpha values: 0.0 (cube control), 0.1 (mild), 0.3 (moderate),
0.5 (significant shear).
"""
from __future__ import annotations

import cmath
import math
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from ngsolve import (
    Mesh, H1, BilinearForm, LinearForm, GridFunction,
    CoefficientFunction, grad, dx, Integrate, TaskManager,
    SetNumThreads
)
from netgen.occ import (
    OCCGeometry, WorkPlane, Pnt, Vec, Prism
)


SIGMA = 5.8e7
MU = 4 * math.pi * 1e-7
MS = MU * SIGMA


# ---------------------------------------------------------------------------
# Parallelepiped construction
# ---------------------------------------------------------------------------


def build_parallelepiped(L, alpha, maxh):
    """Bottom L x L square, extruded in direction (alpha*L, 0, L).

    For alpha = 0: cube of side L.
    For alpha > 0: parallelepiped with 4 non-right dihedrals.
    """
    bot = WorkPlane().Rectangle(L, L).Face()
    shape = Prism(bot, Vec(alpha * L, 0, L))
    shape.mat("body")
    # All faces -> "outer" (Dirichlet on entire boundary)
    for f in shape.faces:
        f.name = "outer"
    geo = OCCGeometry(shape)
    return Mesh(geo.GenerateMesh(maxh=maxh)), shape


def measure_geometry(shape):
    """Compute V, S, edge lengths, vertex count from OCC shape."""
    V = float(shape.mass)
    S = 0.0
    for f in shape.faces:
        S += float(f.mass)  # face area
    # Unique edges (filter out duplicates)
    seen_edges = []
    for e in shape.edges:
        c = e.center
        key = (round(c[0]*1e6), round(c[1]*1e6), round(c[2]*1e6))
        if key not in seen_edges:
            seen_edges.append(key)
    # OCC may report edges per shell -- count unique by center
    edge_data = {}
    for e in shape.edges:
        c = (e.center[0], e.center[1], e.center[2])
        key = (round(c[0]*1e6), round(c[1]*1e6), round(c[2]*1e6))
        if key not in edge_data:
            edge_data[key] = float(e.mass)  # length
    L_total = sum(edge_data.values())
    n_edges = len(edge_data)

    # Unique vertices by position
    vert_set = set()
    for v in shape.vertices:
        p = v.p
        key = (round(p[0]*1e6), round(p[1]*1e6), round(p[2]*1e6))
        vert_set.add(key)
    n_verts = len(vert_set)

    return {
        "V": V,
        "S": S,
        "L_total_unique": L_total,
        "n_edges": n_edges,
        "n_verts": n_verts,
    }


def Y_mellin_naive(s, S, L_edge_sum, n_verts):
    """Naive Mellin: assume all edges are right-angle wedges.
        c_0 = S * sqrt(sigma/mu)
        c_1 = -(8/(pi mu)) * L_edge_sum
              [factor 8 because cube formula -16(a+b+c)/(pi mu) with 4 edges per axis
               = (16/4) * L_edge_sum_per_axis (a+b+c) = 16 * total_edge / 4 ... let me redo]

    Actually: cuboid total edge length = 4*(Lx+Ly+Lz), and Mellin c_1 = -16(Lx+Ly+Lz)/(pi mu)
    So c_1 = -4 * L_total / (pi mu)
            = -(4/pi) * (1/mu) * L_total
        where L_total = 4*(Lx+Ly+Lz).
    For cube L=5mm: 4*15mm = 60mm, c_1 = -4*60e-3/(pi mu) = -76.3e3.
    Cross-check: cube3d formula c_1 = -48L/(pi mu) = -48*5e-3/(pi mu) = -76.3e3. Matches.

    So generalized:
        c_1 = -(4/(pi mu)) * L_total_edges     (where L_total = sum of all edge lengths)
        c_2 = +(48/(pi mu^1.5 sqrt(sigma))) * (n_verts / 8)
              [scale by vertex count; cube has 8 vertices]
    """
    c_0 = S * math.sqrt(SIGMA / MU)
    c_1 = -(4.0 / (math.pi * MU)) * L_edge_sum
    c_2 = (48.0 / (math.pi * MU**1.5 * math.sqrt(SIGMA))) * (n_verts / 8.0)
    return c_0 / cmath.sqrt(s) + c_1 / s + c_2 / s**1.5


# ---------------------------------------------------------------------------
# FEM solve
# ---------------------------------------------------------------------------


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
    return Y_DC * (1 + u_avg), fes.ndof


def compute_dihedral_angles(L, alpha):
    """Closed-form: for the oblique prism shape,
       interior dihedrals at the 4 non-right edges:
           bottom-left:  acute   theta_acute  = 180 - arccos(-alpha/sqrt(1+alpha^2)) - 90
                                              = 90 - arctan(alpha)
           bottom-right: obtuse  theta_obtuse = 90 + arctan(alpha)
           top-left:     obtuse  90 + arctan(alpha)
           top-right:    acute   90 - arctan(alpha)
       (and 8 edges remain at exactly 90 degrees)
    """
    theta_acute = math.pi/2 - math.atan(alpha)
    theta_obtuse = math.pi/2 + math.atan(alpha)
    return {
        "n_right_angle": 8,
        "n_acute": 2,
        "theta_acute_deg": math.degrees(theta_acute),
        "n_obtuse": 2,
        "theta_obtuse_deg": math.degrees(theta_obtuse),
    }


# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------


def run_case(L, alpha, label, maxh_frac=8, order=3):
    print(f"\n=== Case {label}: parallelepiped L={L*1e3:.1f} mm, shear alpha = {alpha:.2f} ===")
    maxh = L / maxh_frac
    mesh, shape = build_parallelepiped(L, alpha, maxh)
    geom = measure_geometry(shape)
    dihedrals = compute_dihedral_angles(L, alpha)
    Y_DC = geom["V"] * SIGMA

    print(f"  V        = {geom['V']*1e9:.3f} mm^3 (cube ref: {L**3*1e9:.3f})")
    print(f"  S        = {geom['S']*1e6:.3f} mm^2 (cube ref: {6*L**2*1e6:.3f})")
    print(f"  L_total  = {geom['L_total_unique']*1e3:.3f} mm (cube ref: {12*L*1e3:.3f})")
    print(f"  n_edges  = {geom['n_edges']}, n_verts = {geom['n_verts']}")
    print(f"  Dihedrals: {dihedrals['n_right_angle']} right-angle, "
          f"{dihedrals['n_acute']} acute ({dihedrals['theta_acute_deg']:.1f} deg), "
          f"{dihedrals['n_obtuse']} obtuse ({dihedrals['theta_obtuse_deg']:.1f} deg)")
    print(f"  Y_DC     = {Y_DC:.4e}")
    print(f"  Mesh ne  = {mesh.ne}")
    print()
    print(f"  {'f (Hz)':>10}  {'|gL|':>8}  {'|Y_FEM|':>13}  {'|Y_naive_Mellin|':>16}  {'rel err':>10}")

    results = []
    for f in [3e4, 1e5, 3e5, 1e6]:
        s = 1j * 2 * math.pi * f
        with TaskManager():
            Y_FEM, ndof = solve_at(mesh, f, order, Y_DC)
        Y_M = Y_mellin_naive(s, geom["S"], geom["L_total_unique"], geom["n_verts"])
        gL = abs(cmath.sqrt(s * MS) * L)
        rel = abs(Y_FEM - Y_M) / abs(Y_M) * 100
        print(f"  {f:10.2e}  {gL:8.2f}  {abs(Y_FEM):13.6e}  {abs(Y_M):16.6e}  {rel:9.4f}%")
        results.append((f, abs(Y_FEM), abs(Y_M), rel))

    return results, geom, dihedrals


def main():
    SetNumThreads(8)
    print("=== Phase 2: parallelepiped (sheared) Mellin extrapolation test ===")
    print(f"L = 5 mm, sigma = {SIGMA:.2e} S/m, mu = {MU:.4e} H/m")
    print()
    print("Hypothesis: Mixed Galerkin framework's codim decomposition extends")
    print("to non-orthogonal hexahedra.  Naive Mellin uses measured geometry")
    print("(S, L_total, n_verts) but assumes right-angle wedge factor W(pi/2).")
    print("If FEM agrees within a few %, the framework is robust to dihedral")
    print("angle perturbations.")

    L = 5e-3
    all_results = []
    for alpha, label in [(0.0, "control (cube)"),
                         (0.1, "mild shear"),
                         (0.3, "moderate shear"),
                         (0.5, "significant shear")]:
        res, geom, dh = run_case(L, alpha, label=f"{label} alpha={alpha}")
        all_results.append((alpha, label, res, geom, dh))

    # Summary
    print()
    print("=" * 72)
    print("Summary: error vs naive Mellin (right-angle-wedge assumption)")
    print("=" * 72)
    print(f"{'alpha':>6}  {'dihedrals (acute/obtuse)':>26}  {'err@1e5':>10}  {'err@3e5':>10}")
    for alpha, label, res, geom, dh in all_results:
        err_1e5 = next(r[3] for r in res if abs(r[0] - 1e5) < 1)
        err_3e5 = next(r[3] for r in res if abs(r[0] - 3e5) < 1)
        dh_str = f"{dh['theta_acute_deg']:.1f}/{dh['theta_obtuse_deg']:.1f} deg"
        print(f"  {alpha:5.2f}  {dh_str:>26}  {err_1e5:9.4f}%  {err_3e5:9.4f}%")


if __name__ == "__main__":
    main()
