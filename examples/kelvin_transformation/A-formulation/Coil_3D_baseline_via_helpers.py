"""Coil_3D_baseline_via_helpers.py

M1 smoke test: rebuild the same FEM problem as
``Coil_3D_A_HCurl_with_Kelvin.py`` but using the new helper API
(kelvin_geometry + kelvin_material). Goal: same mesh statistics and
same L within 0.1 %, in much less code.

Reference:
    Coil_3D_A_HCurl_with_Kelvin.py:    L = 89.44 nH (analytical 88.55,
                                                     err +1.0 %).
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

import radia  # noqa: F401   (MKL DLL paths)

from ngsolve import (HCurl, BilinearForm, LinearForm, GridFunction,
                      CoefficientFunction as CF, Integrate, InnerProduct,
                      Conj, curl, dx, x, y, z, sqrt,
                      Periodic, TaskManager, Mesh)
from netgen.occ import (WorkPlane, Axes, Pnt, Dir, Sphere, Axis,
                         OCCGeometry)


from kelvin_geometry import add_kelvin_exterior_domain
from kelvin_material import make_kelvin_nu_cf, NU_0


def main():
    R_coil = 0.030
    a_coil = 0.003
    gap_deg = 5.0
    R_K = 0.06
    offset = (0.15, 0.0, 0.0)
    maxh = 0.012
    maxh_coil = 0.005
    order = 1
    I_total = 1.0

    L_analytical = (4e-7 * math.pi) * R_coil \
        * (math.log(8 * R_coil / a_coil) - 2) * (1 - gap_deg / 360)

    print("=" * 60)
    print("M1 smoke: baseline rebuilt via Kelvin helper API")
    print("=" * 60)
    print(f"Coil  R={R_coil*1e3:.0f}mm a={a_coil*1e3:.1f}mm "
          f"gap={gap_deg}deg")
    print(f"Kelvin R_K={R_K*1e3:.0f}mm offset=(150, 0, 0)mm")
    print(f"L_analytical = {L_analytical*1e9:.2f} nH")

    # --- Inner physical part: torus (coil) inside a sphere of R_K ---
    wp = WorkPlane(Axes(p=Pnt(R_coil, 0, 0),
                         n=Dir(0, 1, 0), h=Dir(0, 0, 1)))
    circle = wp.Circle(a_coil).Face()
    torus = circle.Revolve(Axis(Pnt(0, 0, 0), Dir(0, 0, 1)),
                            360 - gap_deg)
    torus.name = "coil"; torus.maxh = maxh_coil

    inner_sphere = Sphere(Pnt(0, 0, 0), R_K); inner_sphere.maxh = maxh
    for f in inner_sphere.faces:
        f.name = "kelvin_int"
    inner_air = inner_sphere - torus; inner_air.name = "air"

    # The helper does periodic BC + GND + glue with the offset sphere.
    geometry, info = add_kelvin_exterior_domain(
        inner_air,                  # hmm but coil is part of the input...
        offset=offset, R_K=R_K, inner_maxh=maxh)
    # Above passes only the air; coil torus needs to be glued separately.
    # Workaround: glue manually here for now (M1 will refine API).
    from netgen.occ import Glue
    full_shape = Glue([inner_air, torus, info["gnd_vertex"]])
    # Re-do periodic identification on the glued shape's faces:
    for f in inner_air.faces:
        if f.name == "kelvin_int":
            inner_face = f; break
    # The offset Kelvin sphere must be added too:
    from netgen.occ import IdentificationType
    outer = Sphere(Pnt(*offset), R_K); outer.maxh = maxh * 2
    for f in outer.faces:
        f.name = "kelvin_ext"
    outer.name = "kelvin"
    for f in outer.faces:
        if f.name == "kelvin_ext":
            outer_face = f; break
    inner_face.Identify(outer_face, "periodic",
                         IdentificationType.PERIODIC)

    full_shape = Glue([inner_air, torus, outer, info["gnd_vertex"]])

    print("\n[1/3] Building mesh...")
    t0 = time.perf_counter()
    geo = OCCGeometry(full_shape)
    with TaskManager():
        ngmesh = geo.GenerateMesh(maxh=maxh, grading=0.5)
    mesh = Mesh(ngmesh); mesh.Curve(order + 1)
    print(f"  {mesh.ne} tets, {mesh.nv} verts, "
          f"mats={mesh.GetMaterials()} ({time.perf_counter()-t0:.1f}s)")

    print("\n[2/3] Material + source via helpers...")
    nu_cf = make_kelvin_nu_cf(mesh, R_K, offset, nu_0=NU_0)
    J0 = I_total / (math.pi * a_coil ** 2)
    r_xy = sqrt(x * x + y * y) + 1e-12
    J_src = J0 * CF((-y / r_xy, x / r_xy, 0))

    print("\n[3/3] FEM assemble + solve...")
    t0 = time.perf_counter()
    fes = Periodic(HCurl(mesh, order=order, dirichlet_bbnd="GND"))
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
    W = Integrate(0.5 * nu_cf * InnerProduct(curl(gfu), Conj(curl(gfu))),
                  mesh, order=10).real
    L = 2 * W / I_total ** 2
    print(f"  solve {time.perf_counter()-t0:.1f}s, ndof={fes.ndof}")

    print()
    print("=" * 60)
    print(f"  L (helper-built) = {L*1e9:.4f} nH")
    print(f"  L (analytical)   = {L_analytical*1e9:.4f} nH")
    print(f"  rel err          = {(L/L_analytical - 1)*100:+.3f}%")
    print(f"  Reference baseline (Coil_3D_A_HCurl_with_Kelvin.py): 89.44 nH")
    print("=" * 60)


if __name__ == "__main__":
    main()
