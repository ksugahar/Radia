"""
3D A-formulation (HCurl) with Periodic Kelvin transformation.

This script is the canonical baseline regression test for the Layer 1-3
Kelvin helper API (see ``docs/kelvin/api_plan.md``
M2 milestone). The whole workflow is now a handful of helper calls:

    inner_air+torus  --> add_kelvin_exterior_domain (L1)
    mesh             --> make_kelvin_nu_cf          (L2)
    J_src            --> solve_full_A_kelvin        (L3)
    gfu              --> inductance_from_energy

The reference from the hand-written predecessor (pre-2026-04-15) was
L = 89.44 nH (analytical torus 88.55, err +1.0%).

Usage:
    python Coil_3D_A_HCurl_with_Kelvin.py
"""
from __future__ import annotations

import math
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
for p in (os.path.abspath(os.path.join(HERE, '..', '..', '..', 'src')),
          os.path.abspath(os.path.join(HERE, '..', '..', '..', 'src',
                                        'radia'))):
    if p not in sys.path:
        sys.path.insert(0, p)

import radia  # noqa: F401 (load MKL DLLs)

from ngsolve import (CoefficientFunction as CF, Mesh, TaskManager, x, y,
                      sqrt)
from netgen.occ import (WorkPlane, Axes, Pnt, Dir, Sphere, Axis,
                         OCCGeometry)

from kelvin_geometry import add_kelvin_exterior_domain
from kelvin_solver import (solve_full_A_kelvin, inductance_from_energy,
                            NU_0)


R_coil, a_coil, gap_deg = 0.030, 0.003, 5.0
R_K, offset = 0.06, (0.15, 0.0, 0.0)
maxh, maxh_coil = 0.012, 0.005
order = 1
I_total = 1.0

MU_0 = 4e-7 * math.pi
J0 = I_total / (math.pi * a_coil ** 2)
L_analytical = MU_0 * R_coil * (math.log(8 * R_coil / a_coil) - 2) \
    * (1 - gap_deg / 360)


def build_inner_shapes():
    wp = WorkPlane(Axes(p=Pnt(R_coil, 0, 0), n=Dir(0, 1, 0), h=Dir(0, 0, 1)))
    torus = wp.Circle(a_coil).Face().Revolve(
        Axis(Pnt(0, 0, 0), Dir(0, 0, 1)), 360 - gap_deg)
    torus.name = "coil"
    torus.maxh = maxh_coil

    inner = Sphere(Pnt(0, 0, 0), R_K)
    inner.maxh = maxh
    for f in inner.faces:
        f.name = "kelvin_int"
    air = (inner - torus)
    air.name = "air"
    return [air, torus]


def main():
    print("=" * 60)
    print("3D HCurl A-formulation (driver-based)")
    print(f"Coil R={R_coil*1e3:.0f}mm a={a_coil*1e3:.1f}mm gap={gap_deg}deg")
    print(f"Kelvin R_K={R_K*1e3:.0f}mm offset={offset}")
    print(f"L_analytical = {L_analytical*1e9:.2f} nH")
    print("=" * 60)

    t0 = time.perf_counter()
    geometry, info = add_kelvin_exterior_domain(
        build_inner_shapes(), offset=offset, R_K=R_K, inner_maxh=maxh)
    with TaskManager():
        ngmesh = OCCGeometry(geometry).GenerateMesh(maxh=maxh, grading=0.5)
    mesh = Mesh(ngmesh)
    mesh.Curve(order + 1)
    print(f"  mesh: {mesh.ne} tets, {mesh.nv} verts, "
          f"mats={mesh.GetMaterials()}, {time.perf_counter()-t0:.1f}s")

    r_xy = sqrt(x * x + y * y) + 1e-12
    J_src = J0 * CF((-y / r_xy, x / r_xy, 0))
    t0 = time.perf_counter()
    res = solve_full_A_kelvin(mesh, J_src, R_K=R_K, offset=offset,
                                source_material="coil", order=order)
    L, _ = inductance_from_energy(res["gfu"], res["nu_cf"], mesh, I_total)
    print(f"  solve+energy: {time.perf_counter()-t0:.1f}s, ndof={res['fes'].ndof}")

    print("=" * 60)
    print(f"  L (driver) = {L*1e9:.2f} nH")
    print(f"  L (analyt) = {L_analytical*1e9:.2f} nH")
    print(f"  rel err    = {(L/L_analytical - 1)*100:+.2f}%")
    print("=" * 60)


if __name__ == "__main__":
    main()
