"""test_nu_convention.py

A/B regression test of nu_kelvin conventions in the 3D full-A HCurl
FEM with Sugahara two-sphere Kelvin geometry (gapped torus, circular
cross-section, analytical L = 88.55 nH).

Canonical convention per Nagamine CEFC 2026 / Sugahara 2022 IEEE
TransMag (labeled "canonical" below):

    nu_kelvin = (rho'/R)^2 * nu_0

An "inverted_variant" with (R/rho')^2 is kept for A/B purposes: an
earlier (2026-04-15) session reported the inverted form winning with
+0.90% error on this benchmark and mislabeled it "sugahara", which
was incorrect. Nagamine CEFC 2026 derives (rho'/R)^2 rigorously
(pullback + bilinear energy functional) and validates numerically
on a toroidal loop to +0.33%. See
  examples/kelvin_transformation/CONVENTION.md
  examples/kelvin_transformation/docs/pullback_derivation_3D.md sec 8

Run this to verify the canonical convention reproduces L_analytical
on this benchmark. Any discrepancy from the old empirical result
(where the canonical gave +4.80% error) points to FEM setup issues
(GND / gauge regularization / mesh) and is an open investigation.

Labels:
  (a) "canonical"       : nu_kelvin = (rho'/R)^2 * nu_0   [Nagamine/Sugahara]
  (b) "inverted_variant": nu_kelvin = (R/rho')^2 * nu_0   [historical mislabel]
"""

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

import radia  # noqa

from ngsolve import (HCurl, BilinearForm, LinearForm, GridFunction,
                      CoefficientFunction as CF, Integrate, InnerProduct,
                      Conj, curl, dx, x, y, z, sqrt,
                      Periodic, TaskManager, Mesh, IfPos)
from netgen.occ import (WorkPlane, Axes, Pnt, Dir, Sphere, Axis, Vertex,
                         Glue, OCCGeometry, IdentificationType)


MU_0 = 4e-7 * math.pi
NU_0 = 1.0 / MU_0


def build_baseline_geometry():
    """Same geometry as Coil_3D_A_HCurl_with_Kelvin.py."""
    R_coil = 0.030
    a_coil = 0.003
    gap_deg = 5
    R_K = 0.06
    offset = 0.15
    maxh = 0.012
    maxh_coil = 0.005

    wp = WorkPlane(Axes(p=Pnt(R_coil, 0, 0),
                         n=Dir(0, 1, 0), h=Dir(0, 0, 1)))
    circle = wp.Circle(a_coil).Face()
    torus = circle.Revolve(Axis(Pnt(0, 0, 0), Dir(0, 0, 1)),
                            360 - gap_deg)
    torus.name = "coil"; torus.maxh = maxh_coil

    inner = Sphere(Pnt(0, 0, 0), R_K); inner.maxh = maxh
    for f in inner.faces:
        f.name = "kelvin_int"
    inner_air = inner - torus; inner_air.name = "air"

    outer = Sphere(Pnt(offset, 0, 0), R_K); outer.maxh = maxh * 2
    for f in outer.faces:
        f.name = "kelvin_ext"
    outer.name = "kelvin"

    gnd = Vertex(Pnt(offset, 0, 0)); gnd.name = "GND"

    for f in inner_air.faces:
        if f.name == "kelvin_int":
            ki = f; break
    for f in outer.faces:
        if f.name == "kelvin_ext":
            ke = f; break
    ki.Identify(ke, "periodic", IdentificationType.PERIODIC)

    shape = Glue([inner_air, torus, outer, gnd])
    geo = OCCGeometry(shape)
    with TaskManager():
        ngmesh = geo.GenerateMesh(maxh=maxh, grading=0.5)
    mesh = Mesh(ngmesh); mesh.Curve(2)
    return mesh, R_K, offset, R_coil, a_coil, gap_deg


def solve_with_convention(mesh, R_K, offset, R_coil, a_coil, gap_deg,
                           convention):
    """Solve full-A FEM with the specified nu_kelvin convention."""
    arc_deg = 360 - gap_deg
    I_total = 1.0
    J0 = I_total / (math.pi * a_coil ** 2)

    r_prime_sq = (x - offset) ** 2 + y ** 2 + z ** 2 + 1e-24
    if convention == "canonical":
        # Nagamine CEFC 2026 / Sugahara 2022 IEEE TransMag: nu' = (rho'/R)^2 * nu_0
        kelvin_fac = r_prime_sq / (R_K * R_K)
    elif convention == "inverted_variant":
        # Historical mislabel ("sugahara" in 2026-04-15 session); kept for A/B only
        kelvin_fac = R_K * R_K / r_prime_sq        # (R/rho')^2
    else:
        raise ValueError(convention)

    nu_dict = {m: (NU_0 * kelvin_fac if "kelvin" in m.lower() else NU_0)
               for m in mesh.GetMaterials()}
    nu_cf = mesh.MaterialCF(nu_dict, default=NU_0)

    r_xy = sqrt(x * x + y * y) + 1e-12
    J_src = J0 * CF((-y / r_xy, x / r_xy, 0))

    fes = Periodic(HCurl(mesh, order=1, dirichlet_bbnd="GND"))
    u, v = fes.TnT()

    a_bf = BilinearForm(fes)
    a_bf += nu_cf * curl(u) * curl(v) * dx(bonus_intorder=4)
    a_bf += 1e-8 * NU_0 * u * v * dx
    f_lf = LinearForm(fes)
    f_lf += J_src * v * dx("coil")
    with TaskManager():
        a_bf.Assemble(); f_lf.Assemble()
    gfu = GridFunction(fes)
    with TaskManager():
        gfu.vec.data = a_bf.mat.Inverse(fes.FreeDofs(),
                                         inverse="pardiso") * f_lf.vec

    W_inner = Integrate(0.5 * nu_cf * InnerProduct(curl(gfu),
                                                    Conj(curl(gfu))),
                        mesh,
                        definedon=mesh.Materials("air|coil"),
                        order=10).real
    W_kext = Integrate(0.5 * nu_cf * InnerProduct(curl(gfu),
                                                   Conj(curl(gfu))),
                       mesh,
                       definedon=mesh.Materials("kelvin"),
                       order=10).real
    L = 2.0 * (W_inner + W_kext) / I_total ** 2
    return L, W_inner, W_kext


def main():
    L_analytical = MU_0 * 0.030 * (math.log(8 * 0.030 / 0.003) - 2) \
        * (1.0 - 5.0 / 360.0)
    print("=" * 60)
    print("nu_kelvin convention test (full-A FEM, baseline geometry)")
    print("=" * 60)
    print(f"  Analytical (gapped torus, circular cross-section): "
          f"{L_analytical*1e9:.2f} nH")
    print()

    print("[1/3] Building baseline mesh...")
    t0 = time.perf_counter()
    mesh, R_K, offset, R_coil, a_coil, gap_deg = build_baseline_geometry()
    print(f"  {mesh.ne} tets, {mesh.nv} verts, "
          f"mats={mesh.GetMaterials()} ({time.perf_counter()-t0:.1f}s)")

    for conv in ["canonical", "inverted_variant"]:
        print(f"\n[{conv}] nu_kelvin = "
              f"{'(rho/R)^2' if conv=='canonical' else '(R/rho)^2'} * nu_0")
        t0 = time.perf_counter()
        L, W_in, W_k = solve_with_convention(
            mesh, R_K, offset, R_coil, a_coil, gap_deg, conv)
        L_in = 2 * W_in
        L_k = 2 * W_k
        print(f"  L_inner = {L_in*1e9:.3f} nH, "
              f"L_kelvin = {L_k*1e9:.3f} nH")
        print(f"  L_total = {L*1e9:.3f} nH  "
              f"({(L/L_analytical - 1)*100:+.2f}% vs analytical)  "
              f"({time.perf_counter()-t0:.1f}s)")


if __name__ == "__main__":
    main()
