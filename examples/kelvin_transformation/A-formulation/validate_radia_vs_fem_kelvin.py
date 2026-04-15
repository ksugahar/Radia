"""validate_radia_vs_fem_kelvin.py

One-shot validation: SAME coil (square 2a x 2a torus) computed two ways
on the IDENTICAL physical setup:

  (A) FEM 3D HCurl + Kelvin, volume J source (Sugahara 2022 two-sphere
      with R_K inner + offset Kelvin-exterior sphere)
  (B) Radia rad.ObjArcCur analytical L = (1/I^2) int J . A dV
      over the coil volume (open domain, no BC truncation)

Both represent the IDENTICAL current distribution (uniform volume J
over a 2a x 2a square cross-section torus). For Kelvin FEM, the open
boundary at physical infinity is correctly modelled. So L_A should
match L_B up to the FEM discretization error (mesh + basis order).

Target: |L_A - L_B| / L_B < 1% (aiming for < 0.1% with refinement).
"""

from __future__ import annotations

import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.abspath(os.path.join(HERE, '..', '..', '..', 'src'))
SRC_RADIA = os.path.join(SRC, 'radia')
for p in (SRC, SRC_RADIA):
    if p not in sys.path:
        sys.path.insert(0, p)

import radia as rad   # noqa

from ngsolve import (HCurl, BilinearForm, LinearForm, GridFunction,
                      CoefficientFunction as CF, Integrate, InnerProduct,
                      Conj, curl, dx, x, y, z, sqrt,
                      Periodic, TaskManager, Mesh)
from netgen.occ import (WorkPlane, Axes, Pnt, Dir, Sphere, Axis, Vertex,
                         Glue, OCCGeometry, IdentificationType, MoveTo)


MU_0 = 4e-7 * math.pi
NU_0 = 1.0 / MU_0

R_coil = 0.030
a_coil = 0.003
gap_deg = 5.0
arc_deg = 360.0 - gap_deg
R_K = 0.06
offset = 0.15
maxh = 0.006
maxh_coil = 0.0012
order = 2
I_total = 1.0


def radia_L(nseg=80, phi_quad_per_seg=32, rho_quad=4, z_quad=4):
    """Radia self-inductance L = (1/I^2) int J . A dV for the square torus.

    Integrates J . A over the coil volume using Gauss-Legendre quadrature
    in (rho, phi, z) with phi stepped through Radia's arc pieces.
    """
    J0 = I_total / (2 * a_coil) ** 2   # square section

    rad.UtiDelAll()
    coil = rad.ObjArcCur(
        [0, 0, 0],
        [R_coil - a_coil, R_coil + a_coil],
        [0, math.radians(arc_deg)],
        2 * a_coil, nseg, 'man', 'z', J0)
    # Volume integrated J.A over the coil:
    # dV = rho d rho d phi d z
    t_r, w_r = np.polynomial.legendre.leggauss(rho_quad)
    t_z, w_z = np.polynomial.legendre.leggauss(z_quad)
    t_p, w_p = np.polynomial.legendre.leggauss(phi_quad_per_seg)
    # Map to [R-a, R+a], [-a, a], [0, arc_rad]
    rho = 0.5 * (t_r + 1) * (2 * a_coil) + (R_coil - a_coil)
    w_rho = 0.5 * (2 * a_coil) * w_r
    zz = 0.5 * (t_z + 1) * (2 * a_coil) - a_coil
    w_zz = 0.5 * (2 * a_coil) * w_z
    arc_rad = math.radians(arc_deg)
    phi = 0.5 * (t_p + 1) * arc_rad
    w_phi = 0.5 * arc_rad * w_p
    total = 0.0
    # For each quadrature point, compute A at (rho cos phi, rho sin phi, z)
    # and dot with J_phi = J0 * unit_phi = J0 * (-sin phi, cos phi, 0).
    for i_p, (ph, wp) in enumerate(zip(phi, w_phi)):
        cph, sph = math.cos(ph), math.sin(ph)
        Jphi_unit = np.array([-sph, cph, 0.0])
        J_vec = J0 * Jphi_unit
        for i_r, (r_, wr) in enumerate(zip(rho, w_rho)):
            for i_z, (zv, wz) in enumerate(zip(zz, w_zz)):
                pt = [r_ * cph, r_ * sph, zv]
                A = rad.Fld(coil, 'a', pt)
                JA = J_vec[0] * A[0] + J_vec[1] * A[1] + J_vec[2] * A[2]
                total += JA * r_ * wr * wz * wp   # dV = rho dr dphi dz
    L = total / I_total ** 2
    return L


def fem_L_kelvin_square():
    """FEM HCurl + Kelvin, square 2a x 2a cross-section torus, volume J."""
    J0 = I_total / (2 * a_coil) ** 2

    # Square cross-section torus (matches Radia ObjArcCur geometry)
    wp = WorkPlane(Axes(p=Pnt(R_coil, 0, 0), n=Dir(0, 1, 0), h=Dir(0, 0, 1)))
    wp.MoveTo(-a_coil, -a_coil)
    square = wp.Rectangle(2 * a_coil, 2 * a_coil).Face()
    torus = square.Revolve(Axis(Pnt(0, 0, 0), Dir(0, 0, 1)), arc_deg)
    torus.name = "coil"
    torus.maxh = maxh_coil

    inner_sphere = Sphere(Pnt(0, 0, 0), R_K)
    inner_sphere.maxh = maxh
    for f in inner_sphere.faces:
        f.name = "kelvin_int"
    inner_air = inner_sphere - torus
    inner_air.name = "air"

    outer_sphere = Sphere(Pnt(offset, 0, 0), R_K)
    outer_sphere.maxh = maxh * 2
    for f in outer_sphere.faces:
        f.name = "kelvin_ext"
    outer_sphere.name = "kelvin"

    gnd = Vertex(Pnt(offset, 0, 0)); gnd.name = "GND"

    for f in inner_air.faces:
        if f.name == "kelvin_int":
            k_int = f; break
    for f in outer_sphere.faces:
        if f.name == "kelvin_ext":
            k_ext = f; break
    k_int.Identify(k_ext, "periodic", IdentificationType.PERIODIC)

    shape = Glue([inner_air, torus, outer_sphere, gnd])
    geo = OCCGeometry(shape)
    with TaskManager():
        ngmesh = geo.GenerateMesh(maxh=maxh, grading=0.5)
    mesh = Mesh(ngmesh); mesh.Curve(order + 1)
    print(f"  mesh: {mesh.ne} tets, {mesh.nv} verts, "
          f"mats={mesh.GetMaterials()}")

    # Kelvin-modulated reluctivity
    r_prime_sq = (x - offset) ** 2 + y ** 2 + z ** 2
    kelvin_fac = R_K ** 2 / (r_prime_sq + 1e-20)
    nu_dict = {m: (NU_0 * kelvin_fac if "kelvin" in m.lower() else NU_0)
               for m in mesh.GetMaterials()}
    nu_cf = mesh.MaterialCF(nu_dict, default=NU_0)

    # J source in coil (φ direction around z)
    r_xy = sqrt(x * x + y * y)
    r_xy_s = r_xy + 1e-12
    J_src = J0 * CF((-y / r_xy_s, x / r_xy_s, 0))

    fes_base = HCurl(mesh, order=order, dirichlet_bbnd="GND")
    fes = Periodic(fes_base)
    u, v = fes.TnT()
    print(f"  HCurl Periodic ndof = {fes.ndof}")

    a_bf = BilinearForm(fes)
    a_bf += nu_cf * curl(u) * curl(v) * dx(bonus_intorder=4)
    a_bf += 1e-8 * NU_0 * u * v * dx

    f_lf = LinearForm(fes)
    f_lf += J_src * v * dx("coil")

    with TaskManager():
        a_bf.Assemble(); f_lf.Assemble()
    gfu = GridFunction(fes)
    t0 = time.perf_counter()
    with TaskManager():
        gfu.vec.data = a_bf.mat.Inverse(fes.FreeDofs(),
                                         inverse="pardiso") * f_lf.vec
    t_solve = time.perf_counter() - t0

    W = Integrate(0.5 * nu_cf * InnerProduct(curl(gfu), Conj(curl(gfu))),
                  mesh, order=10).real
    L = 2.0 * W / I_total ** 2
    print(f"  solve {t_solve:.1f}s, W = {W:.6e} J")
    return L


def main():
    print("=" * 60)
    print("Radia ObjArcCur (open-domain) vs FEM+Kelvin (square torus)")
    print("=" * 60)
    print(f"  R_coil = {R_coil*1e3:.0f} mm")
    print(f"  a_coil = {a_coil*1e3:.1f} mm (square 2a x 2a)")
    print(f"  gap    = {gap_deg} deg")
    print(f"  R_K    = {R_K*1e3:.0f} mm, offset = {offset*1e3:.0f} mm")
    print()

    print("[1/2] Radia ObjArcCur self-inductance L = int J.A/I^2 dV ...")
    t0 = time.perf_counter()
    L_rad = radia_L(nseg=200, phi_quad_per_seg=24, rho_quad=8, z_quad=8)
    print(f"  L_Radia = {L_rad*1e9:.4f} nH  "
          f"({time.perf_counter() - t0:.1f}s)")

    print("\n[2/2] FEM HCurl + Kelvin + square torus + volume J ...")
    t0 = time.perf_counter()
    L_fem = fem_L_kelvin_square()
    print(f"  L_FEM   = {L_fem*1e9:.4f} nH  "
          f"({time.perf_counter() - t0:.1f}s)")

    print()
    print("=" * 60)
    print(f"  L_Radia (open-domain)     = {L_rad*1e9:.4f} nH")
    print(f"  L_FEM   (Kelvin baseline) = {L_fem*1e9:.4f} nH")
    print(f"  relative diff  = {(L_fem - L_rad) / L_rad * 100:+.3f}%")
    print("=" * 60)


if __name__ == "__main__":
    main()
