"""v24: Interior-PEC disk verification of Hiruma analytical formulas.

Matches the ∞-cylinder PEC setup (Dirichlet A=0 at conductor surface)
of E0_mode_analytical_Hiruma.m. NGSolve FEM via Hiruma 3-term recurrence
should reproduce analytical R_n, L_n for n=1, 2, 3.

Setup:
- 2D radial mesh: r ∈ [0, a], NO exterior, NO Kelvin
- Dirichlet A=0 at r=a (interior PEC)
- u = r × A_φ
- Hiruma 3-term recurrence
- Compare R_n, L_n with analytical
"""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace", line_buffering=True)
print("v24: Interior PEC disk - analytical Hiruma verification", flush=True)

from netgen.occ import (WorkPlane, Pnt, Axes, X, Z, MoveTo,
                         OCCGeometry)
from ngsolve import (Mesh, H1, BilinearForm, LinearForm,
                     GridFunction, CoefficientFunction, grad, dx,
                     x, y, IfPos, sqrt, Integrate, TaskManager, ngsglobals)
from math import pi
import json, time
from pathlib import Path

a = 10e-3        # disk radius
half_thick = 1e-3  # half-thickness (3D disk: 2mm full)
sigma = 5.8e7    # Cu
mu0 = 4 * pi * 1e-7
mur = 1.0
nu0 = 1.0 / (mu0 * mur)

maxh_disk = 0.05e-3  # very fine to match analytical
ORDER = 3
N_STAGES = 3


def R_analytical(n):
    coeffs = {1: 8, 2: 256, 3: 1944}
    return coeffs[n] / (a**4 * mu0 * mur * sigma**2 * pi)


def L_analytical(n):
    nums = {1: 4, 2: 36, 3: 144}
    dens = {1: 3, 2: 5, 3: 7}
    return nums[n] / (dens[n] * a**2 * sigma * pi)


def build_geo():
    """Half-rectangle 2D radial domain: r ∈ [0, a], z ∈ [-t, t]."""
    wp = WorkPlane(Axes((0, -half_thick, 0), n=Z, h=X))
    disk = wp.Rectangle(a, 2 * half_thick).Face()
    disk.name = "conductor"
    disk.maxh = maxh_disk
    # Edges
    for edge in disk.edges:
        cx = edge.center.x
        if cx < 0.001 * a:
            edge.name = "axis"
        elif abs(cx - a) < 0.01 * a:
            edge.name = "outer"
        else:
            edge.name = "topbot"
    return OCCGeometry(disk, dim=2)


def main():
    ngsglobals.msg_level = 0
    print(f"  Disk: a={a*1000}mm, t={2*half_thick*1000}mm, σ={sigma:.2e}, μᵣ={mur}",
          flush=True)
    print(f"  Order={ORDER}, maxh={maxh_disk*1e3}mm\n", flush=True)

    geo = build_geo()
    print("Mesh...", flush=True)
    with TaskManager():
        ngmesh = geo.GenerateMesh(maxh=maxh_disk)
    mesh = Mesh(ngmesh)
    print(f"  ne = {mesh.ne}", flush=True)
    print(f"  boundaries = {mesh.GetBoundaries()}", flush=True)

    # Interior PEC: Dirichlet u=0 at r=a (outer) AND axis (axis), AND topbot
    # (for axisymmetric pure radial problem we need topbot Dirichlet too,
    # so that solution has no z-dependence for thin disk)
    fes_u = H1(mesh, order=ORDER, dirichlet="axis|outer|topbot")
    print(f"  fes_u ndof = {fes_u.ndof}, FreeDofs = {sum(fes_u.FreeDofs())}",
          flush=True)
    u, v = fes_u.TnT()

    r_weight = IfPos(x - 1e-12, x, 1e-12)
    nu_cf = CoefficientFunction(nu0)
    sigma_cf = CoefficientFunction(sigma)

    # G ↔ curl-curl bilinear form (axisymmetric)
    a_G = BilinearForm(fes_u)
    a_G += nu_cf / r_weight * grad(u) * grad(v) * dx
    print("  Assembling G + factor...", flush=True)
    with TaskManager():
        a_G.Assemble()
        try:
            inv_G = a_G.mat.Inverse(fes_u.FreeDofs(), inverse="pardiso")
        except Exception:
            inv_G = a_G.mat.Inverse(fes_u.FreeDofs(),
                                     inverse="sparsecholesky")

    # C ↔ σ-mass with 1/r factor (axisymmetric)
    a_C = BilinearForm(fes_u)
    a_C += sigma_cf / r_weight * u * v * dx("conductor")
    with TaskManager():
        a_C.Assemble()

    twopi = 2 * pi

    # === Hiruma analytical-style iteration: J(1) = 1 (uniform) ===
    # Source vector b corresponds to "uniform J" excitation
    # weak form: ∫ J × v_test in conductor = ∫ 1 × v in conductor
    # but v here represents r × test_phi (axisymmetric u rep)
    # We need J source as "1" uniform in conductor (azimuthal current density unit)
    # → Linear functional: f(v) = ∫ J × v dr dz = ∫ 1 × v dr dz (axisym)
    f0 = LinearForm(fes_u)
    f0 += 1.0 * v * dx("conductor")
    with TaskManager():
        f0.Assemble()

    # w_1 = G^{-1} b
    w1 = GridFunction(fes_u, name="w_1")
    with TaskManager():
        w1.vec.data = inv_G * f0.vec

    # λ_1 = ∫ J × σ × J × 2π r dr dz = σ × ∫ 1 × 2π r in conductor
    # (per analytical formula: J = 1)
    lam_1_numerical = twopi * float(Integrate(
        sigma_cf / r_weight * 1.0 * 1.0 * r_weight * r_weight  # σ × J² × r
        * dx("conductor"), mesh)) / (2 * half_thick)
    # Hmm careful with axisym integrate vs analytical formula
    # Analytical: λ_1 = σ × ∫_0^a 1 × 2π r dr = σ × π × a² (per unit length)
    # Numerical (per unit thickness): integrate (σ × J²) × 2π × r over r,
    # divided by t to get per-thickness
    # Actually for direct equivalence, let me just compute:
    #   σ × π × a² (analytical, per unit length)
    lam_1_analytical_per_length = sigma * pi * a**2
    print(f"\n  λ_1 analytical (per unit length) = {lam_1_analytical_per_length:.6e}",
          flush=True)
    print(f"  λ_1 = σ × π × a² → R_anal(1) = {1/lam_1_analytical_per_length:.6e}",
          flush=True)
    print(f"\n  Hmm — but analytical R(1) = {R_analytical(1):.6e} (different units)",
          flush=True)

    # OK normalization issue. Let me just verify the per-stage L_n via B² form.
    # By the analytical iteration:
    # A(1) = σ × FE_sol(J=1) / λ(1)  → A is FE solution scaled by σ/λ_1
    # λ(2) = ∫ ν |∇A(1)|² × 2π r dr = "L_n_analytical formula"
    #
    # In our numerical: w_1 = G^{-1} b where b = ∫ J=1 × v in conductor.
    # By our weak form a(u, v) = ∫ ν/r × ∇u·∇v dr dz × 2π = ∫ J × v dr dz × 2π
    # we have a(w_1, w_1) = ∫ J × w_1 × 2π dr dz = ∫ ν/r |∇w_1|² × 2π
    # So a(w_1, w_1) = the "B² of w_1" integral.

    a_w1_w1 = twopi * float(Integrate(
        nu_cf / r_weight * grad(w1) * grad(w1) * dx, mesh))
    print(f"\n  a(w_1, w_1) = ∫_full ν/r |∇w_1|² × 2π = {a_w1_w1:.6e}",
          flush=True)

    # Compare: in analytical Hiruma, A(1) = σ × FE_sol / λ(1)
    # Our w_1 = FE_sol of (∫ J × v = ∫ 1 × v over conductor → ν/r ∇w_1 · ∇v test)
    # Our w_1 corresponds to "FE_sol" itself (no σ scaling, no λ_1 division)
    #
    # Analytical A(1) = σ × FE_sol / λ_1 = σ × w_1 / λ_1_anal
    # ∫ ν |∇A(1)|² = (σ/λ_1_anal)² × ∫ ν |∇w_1|² = (σ/λ_1_anal)² × a_w1_w1 / (2π)?
    #
    # Actually the analytical script integrates over (r) only (per unit length).
    # Our numerical integrates over (r, z). For uniform z (Dirichlet at top/bot
    # forces u → 0 there but in interior should be nearly z-uniform if t<<R).
    #
    # Hmm interior PEC at r=a forces A=0 at outer boundary AND at top/bot.
    # This makes A 2D: depends on both r AND z. Not purely 2D radial.
    #
    # For pure 2D radial comparison, would need different geometry (1D radial only).

    print("\n  NOTE: NGSolve 2D axisymmetric ≠ analytical 1D radial.", flush=True)
    print("  For exact comparison, need 1D radial mesh OR thick cylinder.",
          flush=True)
    print("\n  Per-thickness scaling: R_3D = R_anal/(2t), L_3D = L_anal × 2t",
          flush=True)
    print("  τ = L/R = (L_anal × 2t)/(R_anal/(2t)) = L_anal × R_anal × 4t² ≠ const",
          flush=True)
    print("  Hmm, scaling unclear without careful normalization.", flush=True)


if __name__ == "__main__":
    main()
