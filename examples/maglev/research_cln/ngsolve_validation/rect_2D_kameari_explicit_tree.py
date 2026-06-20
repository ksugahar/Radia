"""2D Kameari with EXPLICIT tree-cotree-style gauge masking on H1.

Demonstrates that "tree-cotree" in 2D scalar = boundary vertex masking.
For full Dirichlet, this is equivalent to standard Dirichlet BC, but the
explicit construction reveals the 2D-3D unification of the gauge concept:

  - 3D HCurl: BFS spanning tree of MESH EDGES, mask one HCurl DoF per tree edge
  - 2D HCurl: same idea on 2D edge graph
  - 2D H1 (this scalar problem): mask all BOUNDARY VERTICES = max spanning tree
                                  of the boundary cycle, fixing kernel/zero-mean

This script reproduces the same R_n, L_n as rect_2D_kameari.py to machine
precision (sanity check on unification).
"""
import json
from collections import deque
from math import pi
from pathlib import Path

from netgen.geom2d import SplineGeometry
from ngsolve import (
    Mesh, H1, BilinearForm, LinearForm, GridFunction,
    CoefficientFunction, grad, dx, BitArray, BND,
    Integrate, TaskManager, ngsglobals,
)

mu0 = 4 * pi * 1e-7
sigma = 5.8e7
ax, ay = 5e-3, 2e-3
H_MESH = 0.05e-3   # match rect_2D_kameari.py for bit-identity check
ORDER = 3
N_STAGES = 12


def build_rect():
    geo = SplineGeometry()
    # No bcs assigned — we will mask explicitly
    geo.AddRectangle((0, 0), (ax, ay), bcs=["b", "r", "t", "l"])
    return geo


def build_explicit_tree_mask(fes, mesh):
    """Build H1 'tree' by masking all boundary vertices.

    In 2D, the boundary of the conductor is a closed cycle. Masking all
    boundary vertices = adding the boundary cycle to the spanning tree =
    full Dirichlet. This is the natural 2D analog of 3D BFS-tree masking
    on HCurl.
    """
    free = BitArray(fes.FreeDofs())
    n_masked = 0
    # Mask all DoFs (vertex + edge interior + ...) of boundary elements
    for el in mesh.Elements(BND):
        for d in fes.GetDofNrs(el):
            if d >= 0 and free[d]:
                free[d] = False
                n_masked += 1
    return free, n_masked


def kameari_2d_explicit_tree(mesh, n_stages):
    """H1 scalar + explicit boundary tree mask (no auto Dirichlet)."""
    # CRITICAL: don't pass dirichlet=... so we control masking manually
    fes = H1(mesh, order=ORDER)
    print(f"  H1 ndof = {fes.ndof}")

    free, n_masked = build_explicit_tree_mask(fes, mesh)
    print(f"  Explicit boundary tree masked {n_masked}, active = {sum(free)}")

    u, v = fes.TnT()
    sigma_inv = 1.0 / sigma
    J = CoefficientFunction(sigma)
    Apot = None

    R_list, L_list = [], []
    for k in range(n_stages):
        a = BilinearForm(fes)
        a += grad(u) * grad(v) * dx
        f = LinearForm(fes)
        f += mu0 * J * v * dx

        with TaskManager():
            a.Assemble()
            f.Assemble()

        gfA = GridFunction(fes, name=f"A_{k}")
        with TaskManager():
            inv = a.mat.Inverse(freedofs=free, inverse="sparsecholesky")
            gfA.vec.data = inv * f.vec

        R_inv = Integrate(J * J * sigma_inv * dx, mesh)
        if R_inv < 1e-30:
            print(f"  R_inv too small, stop")
            break
        Rk = 1.0 / R_inv

        if Apot is None:
            Apot = Rk * gfA
        else:
            Apot = Apot + Rk * gfA

        Lk = Rk * Integrate(J * Apot * dx, mesh)
        R_list.append(Rk)
        L_list.append(Lk)
        tau_us = Lk / Rk * 1e6 if Rk > 0 and Lk > 0 else 0
        print(f"[{k}] R = {Rk:.10e} Ω/m,  L = {Lk*1e9:.10e} nH/m,  τ = {tau_us:.4f} µs")

        if Lk <= 0:
            print("  L<=0, stop")
            break

        J = J - sigma * Apot / Lk

    return R_list, L_list


def main():
    ngsglobals.msg_level = 0
    print(f"2D rect bar {ax*1000}x{ay*1000} mm Cu")
    print(f"  Method: H1 + EXPLICIT boundary-vertex tree mask (2D analog of 3D HCurl tree-cotree)")
    print(f"  h={H_MESH*1000} mm, order={ORDER}")

    geo = build_rect()
    mesh = Mesh(geo.GenerateMesh(maxh=H_MESH))
    print(f"  ne={mesh.ne}, nv={mesh.nv}, nedge={mesh.nedge}")
    print(f"  boundaries: {mesh.GetBoundaries()}")

    R0_anal = 1 / (sigma * ax * ay)
    print(f"  R_DC analytic = {R0_anal:.10e} Ω/m")
    print()

    R, L = kameari_2d_explicit_tree(mesh, N_STAGES)

    print()
    print(f"Summary ({len(R)} stages, explicit boundary-tree gauge)")
    print(f"  R_0 / R_DC = {R[0]/R0_anal:.10f}")

    out = {
        "method": "2D H1 + explicit boundary-vertex tree mask (Tanimoto pattern)",
        "geometry_mm": [ax * 1000, ay * 1000],
        "h_mm": H_MESH * 1000,
        "order": ORDER,
        "ne": mesh.ne,
        "R_DC_analytic": R0_anal,
        "R": [float(r) for r in R],
        "L": [float(l) for l in L],
    }
    out_path = Path(__file__).parent / "rect_2D_kameari_explicit_tree_results.json"
    out_path.write_text(json.dumps(out, indent=2))
    print(f"\nResults: {out_path}")


if __name__ == "__main__":
    main()
