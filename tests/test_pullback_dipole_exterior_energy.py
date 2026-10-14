"""Reproduce §8.2 of pullback_derivation_3D.md.

Toroidal current loop:
- major radius a = 0.1 m
- wire radius b = 0.01 m
- magnetic moment m = I * pi * a^2 = 1 A m^2  (so I = 1 / (pi * 0.01) = 31.83 A)
- Kelvin sphere R = 1 m

Exterior magnetic energy (r > R):
    W_analytical = mu_0 * m^2 / (12 * pi * R^3) = 3.333e-8 J

Tests in isolation (no FE solve):
1. Build A_phys analytically as a magnetic dipole: A_phys = mu_0/(4 pi) m x r / |r|^3
2. Pull back to Omega' via Phase 2 (R/rho')^2 H A_phys(r_phys)
3. Compute W' = int_{Omega'} (1/2) nu' |curl(A_kext)|^2 dV'
   where nu' = (rho'/R)^2 nu_0
4. Compare to W_analytical

If this gives +0.33% (matching docs/kelvin/KELVIN_TRANSFORMATION.md §2.8),
the pullback + material-factor combo is correct in isolation.
"""
import sys
import math
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

import numpy as np

from netgen.occ import Sphere, Pnt, OCCGeometry, Vertex, Glue
from ngsolve import (
    Mesh, HCurl, GridFunction, CoefficientFunction, CF, Integrate,
    InnerProduct, curl, dx, x, y, z, sqrt, TaskManager,
)

MU_0 = 4e-7 * math.pi
NU_0 = 1.0 / MU_0


def main():
    print("=" * 60)
    print("Reproduce §8.2 dipole exterior energy validation")
    print("=" * 60)

    # Coil parameters (from §8.2)
    a_loop = 0.1       # major radius [m]
    R_K = 1.0          # Kelvin sphere radius [m]
    m_dipole = 1.0     # magnetic moment [A m^2]

    W_analytical = MU_0 * m_dipole**2 / (12.0 * math.pi * R_K**3)
    print(f"  m_dipole = {m_dipole} A m^2")
    print(f"  R_K = {R_K} m")
    print(f"  W_analytical (dipole exterior) = {W_analytical:.4e} J")
    print()

    # === Mesh: just the Kelvin exterior domain (kext alone) ===
    # In computational coords, kext is centered at offset, radius R_K.
    # We don't need the inner domain for this test — just the kext.
    print("[1/4] Building kext-only mesh...")
    offset = (5.0, 0.0, 0.0)   # arbitrary; physically this represents r > R
    kext_sphere = Sphere(Pnt(*offset), R_K)
    kext_sphere.mat("kelvin_air")
    for f in kext_sphere.faces:
        f.name = "kelvin_ext"
    kext_sphere.maxh = 0.15

    gnd = Vertex(Pnt(*offset))
    gnd.name = "GND"

    geo = Glue([kext_sphere, gnd])
    geo.solids[0].name = "kelvin_air"

    with TaskManager():
        ngmesh = OCCGeometry(geo).GenerateMesh(maxh=0.15, grading=0.5)
    mesh = Mesh(ngmesh)
    mesh.Curve(3)
    print(f"  mesh: {mesh.ne} tets, {mesh.nv} vertices")

    # === Phase 2 pullback A_kext from physical magnetic dipole at origin ===
    # Physical dipole: A_phys(r) = mu_0/(4pi) * m x r / |r|^3, m = m_dipole * z_hat
    # so A_phys = mu_0 m / (4 pi) * (-y, x, 0) / r^3.
    # Kelvin map: r_phys = c + (R^2 / |r' - c|^2) (r' - c)
    #            = c + (R/rho')^2 * (r' - c)  scaled vector.
    # Householder n = (r' - c) / rho'.
    print("[2/4] Building pullback A_kext via Phase 2 (analytic CF)...")

    ox, oy, oz = offset
    dxp = x - ox
    dyp = y - oy
    dzp = z - oz
    rho_p_sq = dxp * dxp + dyp * dyp + dzp * dzp + 1e-24
    rho_p = sqrt(rho_p_sq)

    # Mapped physical coordinate (in original physical frame, origin at 0)
    # Note: physical origin is at (0,0,0), Kelvin sphere is at offset.
    # Kelvin map is from offset (computational frame): r' is local in kext.
    # Physical exterior point is r_phys = (R^2 / rho'^2) * (r' - offset) RELATIVE to offset.
    # But the dipole sits at PHYSICAL origin (0,0,0), so we evaluate
    # A_phys at the OFFSET-RELATIVE physical point. In our mapping, the
    # physical exterior r > R corresponds to r' < R in offset-local coords,
    # WITH the inner sphere centered at origin and the kext at offset.
    # Specifically: r_phys = (R^2/rho_p^2) * (r' - offset) but interpreted
    # as a point (R/rho_p)*R away from origin in the physical setup.
    #
    # Simplification: r_phys.norm = R^2 / rho_p (purely radial scale-up).
    # For our axis-symmetric dipole at origin and kext offset along x-axis,
    # the proper mapping is:
    #   r' - offset has length rho_p in kext coords;
    #   the corresponding physical point lies at distance R^2/rho_p from origin,
    #   in the direction of (r' - offset).
    # But geometrically the kext is OUTSIDE the inner sphere — same direction
    # in physical space. So r_phys = (R/rho_p)^2 * (r' - offset).
    # Length: |r_phys| = R^2 / rho_p.

    r_phys_x = (R_K**2 / rho_p_sq) * dxp
    r_phys_y = (R_K**2 / rho_p_sq) * dyp
    r_phys_z = (R_K**2 / rho_p_sq) * dzp
    r_phys_norm = sqrt(r_phys_x**2 + r_phys_y**2 + r_phys_z**2 + 1e-24)
    r_phys3 = r_phys_norm**3

    # A_phys (dipole at origin): mu_0 m / (4 pi) * (m_hat x r_phys) / |r_phys|^3
    # m_hat = z_hat: m_hat x r_phys = (-r_phys_y, r_phys_x, 0)
    pref = MU_0 * m_dipole / (4.0 * math.pi)
    A_phys_x = pref * (-r_phys_y) / r_phys3
    A_phys_y = pref * (r_phys_x) / r_phys3
    A_phys_z = CoefficientFunction(0.0)

    # Phase 2 1-form pullback: A_kext = (R/rho')^2 (A_phys - 2 (A_phys . n) n)
    nx = dxp / rho_p
    ny = dyp / rho_p
    nz = dzp / rho_p
    A_dot_n = A_phys_x * nx + A_phys_y * ny + A_phys_z * nz
    refl_x = A_phys_x - 2.0 * A_dot_n * nx
    refl_y = A_phys_y - 2.0 * A_dot_n * ny
    refl_z = A_phys_z - 2.0 * A_dot_n * nz
    factor = (R_K / rho_p)**2

    A_kext_cf = CoefficientFunction((
        factor * refl_x,
        factor * refl_y,
        factor * refl_z,
    ))

    # === Project to HCurl (high order) and compute curl ===
    print("[3/4] Projecting A_kext to HCurl, computing curl...")
    fes = HCurl(mesh, order=3)
    gfA = GridFunction(fes)
    with TaskManager():
        gfA.Set(A_kext_cf, bonus_intorder=8)
    curlA = curl(gfA)
    print(f"  ndof = {fes.ndof}")

    # === Energy with nu' = (rho'/R)^2 nu_0 ===
    print("[4/4] Energy integral with nu' = (rho'/R)^2 nu_0...")
    nu_factor = rho_p_sq / R_K**2     # = (rho'/R)^2
    nu_kext = nu_factor * NU_0

    with TaskManager():
        W = float(Integrate(0.5 * nu_kext * InnerProduct(curlA, curlA),
                            mesh, order=10))

    print()
    print("=" * 60)
    print(f"  W (FEM kext, Phase 2 pullback) = {W:.4e} J")
    print(f"  W_analytical                  = {W_analytical:.4e} J")
    print(f"  rel err                       = "
          f"{(W / W_analytical - 1) * 100:+.2f}%")
    print("=" * 60)
    print()
    if abs(W / W_analytical - 1) < 0.05:
        print("  [OK] Pullback + nu' = (rho'/R)^2 nu_0 reproduces dipole")
        print("       exterior energy in isolation.")
        print("       => Phase 2 pullback formula is CORRECT.")
    else:
        print("  [NG] Pullback + nu' result deviates significantly.")
        print("       Either the pullback formula or nu' factor has a bug.")


if __name__ == "__main__":
    main()
