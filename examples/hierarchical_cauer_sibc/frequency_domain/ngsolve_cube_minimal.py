"""
Minimal NGSolve FEM eddy-current solve on a Cu cube.

Pragmatic geometry: Cu cube inside a large air box (no Kelvin sphere).
Outer boundary = perfect magnetic conductor (Dirichlet H_t = 0).  This
truncates the open domain but if the air box is large enough relative
to the conductor, the truncation error is small at the conductor surface.

Goal: compute Y(omega) at a few frequencies and compare to Mellin asymptote
      c_0/sqrt(s) + c_1/s + c_2/s^(3/2)
      where c_0 = 6 L^2 sqrt(sigma/mu).

NOTE: This is a SIMPLIFIED test (Dirichlet truncation, not Kelvin).
Production-quality would use Cubit + calc_fem_kelvin.py.

Usage:
    python ngsolve_cube_minimal.py
"""

import sys
import time
import numpy as np
from pathlib import Path


def try_ngsolve():
    try:
        from ngsolve import Mesh
        return True
    except ImportError:
        return False


def build_cube_geometry(L_cube=10e-3, L_air=80e-3, maxh_cube=2e-3, maxh_air=15e-3):
    """Cu cube inside an air box (large outer box, simple Dirichlet)."""
    from netgen.occ import OCCGeometry, Box, Pnt, Glue
    from ngsolve import Mesh

    # Cu cube centred at origin
    cube = Box(Pnt(-L_cube/2, -L_cube/2, -L_cube/2),
               Pnt(L_cube/2, L_cube/2, L_cube/2))
    cube.mat("cu")

    # Air box (larger)
    air = Box(Pnt(-L_air/2, -L_air/2, -L_air/2),
              Pnt(L_air/2, L_air/2, L_air/2)) - cube
    air.mat("air")

    # Outer boundary on air box
    for f in air.faces:
        center = f.center
        # Outer face has max coordinate
        max_coord = max(abs(center[0]), abs(center[1]), abs(center[2]))
        if abs(max_coord - L_air/2) < 1e-9:
            f.name = "outer"

    geo = OCCGeometry(Glue([cube, air]))
    print(f"Mesh generation (cube maxh={maxh_cube*1e3:.1f}mm, air maxh={maxh_air*1e3:.1f}mm)...")
    ngmesh = geo.GenerateMesh(maxh=maxh_air)
    # No curving needed for box geometry
    mesh = Mesh(ngmesh)
    print(f"  ne = {mesh.ne} elements, {mesh.nv} vertices")
    return mesh


def solve_cube_frequency(mesh, omega, sigma_cu, mu_0, H_ext=1.0,
                          fes_order=1):
    """
    Solve eddy current problem at frequency omega.

    Weak form for H-formulation in HCurl:
      a(H, v) = (1/mu) (curl H, curl v) + (j omega sigma) (H, v)_cube
              = j omega (B_ext . v)_cube + boundary forcing
    Approximation: impose H_ext as a known external field, solve for the
    perturbation H_p = H - H_ext.  Total field H_tot = H_ext + H_p.

    Returns: complex Y = (j omega) * integral_cube sigma |H_tot - H_ext|^2 / |H_ext|^2
             (effective eddy-current admittance per unit volume)
    """
    from ngsolve import (HCurl, BilinearForm, LinearForm, GridFunction,
                          CoefficientFunction, curl, dx, Preconditioner,
                          CGSolver, x, y, z, sqrt, Integrate)

    # Material constants per region
    sigma_cf = mesh.MaterialCF({"cu": sigma_cu, "air": 0.0}, default=0.0)
    nu_cf = mesh.MaterialCF({"cu": 1.0/mu_0, "air": 1.0/mu_0}, default=1.0/mu_0)

    # HCurl FES with Dirichlet on outer boundary (H_perp set, or use natural BC)
    fes = HCurl(mesh, order=fes_order, complex=True, dirichlet="outer")

    u, v = fes.TnT()
    a = BilinearForm(fes, symmetric=True)
    a += nu_cf * curl(u) * curl(v) * dx
    a += 1j * omega * sigma_cf * u * v * dx

    # Add a small regularisation in air to handle the singular curl-curl in non-conducting region
    eps_reg = 1e-3 * 1.0/mu_0
    a += eps_reg * u * v * dx("air")

    # External source: uniform H_z = H_ext in cube + air
    H_ext_vec = CoefficientFunction((0, 0, H_ext))
    f = LinearForm(fes)
    f += -1j * omega * sigma_cf * H_ext_vec * v * dx("cu")

    print(f"  Assembling... ", end="", flush=True)
    t0 = time.time()
    with TaskManager():
        a.Assemble()
        f.Assemble()
        print(f"{time.time()-t0:.1f}s")

        print(f"  Solving (direct LU, problem size ~{fes.ndof})... ",
              end="", flush=True)
        t0 = time.time()
        gfu = GridFunction(fes)
        gfu.vec.data = a.mat.Inverse(freedofs=fes.FreeDofs(), inverse="sparsecholesky") * f.vec
        print(f"{time.time()-t0:.1f}s")

        # Compute Y: power dissipation per cycle = 1/2 omega integral sigma |E|^2
        # E = -j omega A = ...  here gfu is H_perturbation
        # Actually for this formulation, integral_cube sigma |H_p|^2 = induced eddy
        from ngsolve import Conj
        from ngsolve import TaskManager
        H_p_sq = gfu * Conj(gfu)  # complex |H|^2
        integral_cu = Integrate(H_p_sq, mesh, definedon=mesh.Materials("cu")).real
        # Volume of cube
        V_cu = Integrate(CoefficientFunction(1), mesh,
                         definedon=mesh.Materials("cu"))
        print(f"  H_p^2 integral on cube: {integral_cu:.4e}, V_cube = {V_cu:.4e}")

        # Approximate Y: imaginary part is reactive, real part is loss
        # Y ~ (j omega) * sigma * V_eff where V_eff = integral |H_p|^2 / |H_ext|^2
        # This is a rough scaling -- a more proper extraction would use the full
        # admittance definition Y = (induced current) / (driving voltage).
        Y_approx = 1j * omega * sigma_cu * integral_cu / abs(H_ext)**2
        return Y_approx, gfu


def mellin_asymptote(s, L, sigma, mu):
    """Cube Mellin asymptote c_0/sqrt(s) + c_1/s + c_2/s^(3/2)."""
    c_0 = 6 * L**2 * np.sqrt(sigma / mu)
    c_1 = -48 * L / (np.pi * mu)
    c_2 = 48 / (np.pi * mu**1.5 * np.sqrt(sigma))
    return c_0 / np.sqrt(s) + c_1 / s + c_2 / s**1.5


if __name__ == "__main__":
    if not try_ngsolve():
        print("NGSolve not available; cannot run production solve.")
        sys.exit(0)

    sigma_cu = 5.96e7
    mu_0 = 4 * np.pi * 1e-7
    L = 10e-3

    print(f"NGSolve cube production run: L = {L*1e3:.1f}mm, Cu")
    print()

    # Build mesh
    print("Step 1: Build geometry + mesh")
    mesh = build_cube_geometry(L_cube=L, L_air=80e-3,
                                maxh_cube=2e-3, maxh_air=15e-3)
    print()

    # Frequency sweep at 3 representative frequencies
    print("Step 2: Frequency sweep")
    f_list = [1e3, 1e4, 1e5]  # 1 kHz, 10 kHz, 100 kHz
    results = []
    for f in f_list:
        omega = 2 * np.pi * f
        print(f"\nFrequency: {f/1e3:.0f} kHz (omega = {omega:.3e})")
        Y, gfu = solve_cube_frequency(mesh, omega, sigma_cu, mu_0,
                                       H_ext=1.0, fes_order=1)
        Y_mellin = mellin_asymptote(1j * omega, L, sigma_cu, mu_0)
        print(f"  Y_FEM     = {Y:.4e}")
        print(f"  Y_Mellin  = {Y_mellin:.4e}")
        if abs(Y_mellin) > 0:
            ratio = abs(Y) / abs(Y_mellin)
            print(f"  |Y_FEM|/|Y_Mellin| = {ratio:.4f}")
        results.append((f, omega, Y, Y_mellin))

    # Save results
    out_dir = Path(__file__).parent
    np.savez(out_dir / "ngsolve_cube_results.npz",
             freq=np.array([r[0] for r in results]),
             omega=np.array([r[1] for r in results]),
             Y_FEM=np.array([r[2] for r in results]),
             Y_Mellin=np.array([r[3] for r in results]))
    print()
    print(f"Saved: {out_dir / 'ngsolve_cube_results.npz'}")
    print()
    print("NOTE: this is a simplified FEM (no Kelvin transformation, Dirichlet")
    print("on outer air box).  The Y_FEM admittance extraction uses a")
    print("post-processing integral of |H_p|^2 which is an APPROXIMATION of")
    print("the full per-port admittance.  Production-quality requires Kelvin")
    print("boundary + proper port definition (see calc_fem_kelvin.py).")
