"""2D rectangular bar (a x b) — Kameari CLN extraction (Tanimoto pattern).

Per-unit-length analysis. A_z(x,y) is scalar on 2D rectangular cross-section.
Full Dirichlet on perimeter. H1 elements (no gauge issue: -nabla^2 has no kernel).

Kameari iteration with ACCUMULATED A_potential (matching Tanimoto's formula):
  Stage k:
    Solve  -nabla^2 A_k = mu * J_k  with A_k=0 on perimeter (per unit length, scalar)
    R_k = 1 / <J_k^2 / sigma>           (scalar inner product)
    Apot_k = sum_{j<=k} R_j * A_j        (accumulated weighted potential)
    L_k = R_k * <J_k, Apot_k>            (= R_k^2 <J_k,A_k> for orthogonalized J)
    J_{k+1} = J_k - sigma * A_k / L_k

This pattern matches Tanimoto's W:/00_CAE/NGSolve/谷本/20241129_2次元CLN練習.ipynb
which validated for circular wire L_(2k+1) = mu/(8 pi (k+1)).
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

mu0 = 4 * pi * 1e-7
sigma = 5.8e7
ax, ay = 5e-3, 2e-3
H_MESH = 0.05e-3   # 50 µm — finer mesh for higher precision
ORDER = 3          # cubic elements
N_STAGES = 12


def build_rect():
    geo = SplineGeometry()
    geo.AddRectangle((0, 0), (ax, ay), bcs=["bot", "right", "top", "left"])
    return geo


def kameari_2d(mesh, n_stages):
    """Tanimoto pattern: accumulated A potential, L_k = R_k * <J_k, sum R_j A_j>."""
    fes = H1(mesh, order=ORDER, dirichlet="bot|right|top|left")
    print(f"  H1 ndof = {fes.ndof}")

    u, v = fes.TnT()
    sigma_inv = 1.0 / sigma

    J = CoefficientFunction(sigma)  # J_0 = sigma * E_0 with E_0 = 1
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
            inv = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")
            gfA.vec.data = inv * f.vec

        # R_k = 1 / <J_k * J_k / sigma> — same in 2D and 3D
        R_inv = Integrate(J * J * sigma_inv * dx, mesh)
        if R_inv < 1e-30:
            print(f"  R_inv too small ({R_inv:.2e}), stop")
            break
        Rk = 1.0 / R_inv

        # Accumulate Apotential
        if Apot is None:
            Apot = Rk * gfA
        else:
            Apot = Apot + Rk * gfA

        # L_k = R_k * <J_k, Apot_k>  (Tanimoto formula)
        Lk = Rk * Integrate(J * Apot * dx, mesh)

        R_list.append(Rk)
        L_list.append(Lk)
        tau_us = Lk / Rk * 1e6 if Rk > 0 and Lk > 0 else 0
        print(f"[{k}] R = {Rk:.10e} Ω/m,  L = {Lk*1e9:.10e} nH/m,  τ = {tau_us:.4f} µs")

        if Lk <= 0:
            print("  WARNING: L<=0, stop")
            break

        # Tanimoto J update: subtract weighted accumulated potential
        J = J - sigma * Apot / Lk

    return R_list, L_list


def main():
    ngsglobals.msg_level = 0
    print(f"2D rect bar {ax*1000}x{ay*1000} mm Cu (per unit length)")
    print(f"  h={H_MESH*1000} mm, order={ORDER}, sigma={sigma:.2e}")

    geo = build_rect()
    mesh = Mesh(geo.GenerateMesh(maxh=H_MESH))
    print(f"  ne={mesh.ne}, nv={mesh.nv}, nedge={mesh.nedge}")

    R0_anal = 1 / (sigma * ax * ay)
    print(f"  R_DC analytic = 1/(σ·a·b) = {R0_anal:.10e} Ω/m")
    print()

    R, L = kameari_2d(mesh, N_STAGES)

    print()
    print(f"Summary ({len(R)} stages)")
    print(f"  R_0 / R_DC = {R[0]/R0_anal:.10f}")

    out = {
        "method": "2D Kameari Tanimoto pattern (accumulated Apotential), H1 scalar",
        "geometry_mm": [ax * 1000, ay * 1000],
        "h_mm": H_MESH * 1000,
        "order": ORDER,
        "ne": mesh.ne,
        "R_DC_analytic": R0_anal,
        "R": [float(r) for r in R],
        "L": [float(l) for l in L],
    }
    out_path = Path(__file__).parent / "rect_2D_kameari_results.json"
    out_path.write_text(json.dumps(out, indent=2))
    print(f"\nResults: {out_path}")


if __name__ == "__main__":
    main()
