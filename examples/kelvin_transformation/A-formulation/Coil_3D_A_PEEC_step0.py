"""Coil_3D_A_PEEC_step0.py

Step 0 validation (pre-reduced-A, pre-Kelvin): verify that the
symbolic Biot-Savart A_s CF constructed by the closed-form
finite-segment formula is correct, using NGSolve purely as a
calculator for the magnetic energy integral.

Procedure
---------
1. Build a PEEC filament bundle matching the torus R=30mm, a=3mm.
2. Mesh a large bounded air sphere (radius >> R_coil).
3. Project A_s_cf onto an HCurl GridFunction via Set().
4. Compute W = int 0.5 nu_0 |curl(gfu_As)|^2 dV over the sphere.
5. L = 2 W / I_total^2.
6. Compare with analytical torus L = mu_0 R (ln(8R/a) - 2) (1 - gap/360).

Expected: L within a few percent of analytical. The mesh boundary
being at finite radius introduces a small truncation error (the
field still has finite energy outside), which is acceptable for
this sanity test.
"""

from __future__ import annotations

import argparse
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

import radia  # noqa: F401  -- MKL DLL paths

from ngsolve import (HCurl, GridFunction, CoefficientFunction, Integrate,
                      InnerProduct, curl, dx, x, y, z, sqrt, log,
                      TaskManager, Mesh)
from netgen.occ import Sphere, Pnt, OCCGeometry


MU_0 = 4e-7 * math.pi
NU_0 = 1.0 / MU_0


def biot_savart_A_cf(filament_paths, currents, x_cf, y_cf, z_cf):
    """Analytical finite-segment Biot-Savart as NGSolve CF."""
    Ax = CoefficientFunction(0.0)
    Ay = CoefficientFunction(0.0)
    Az = CoefficientFunction(0.0)
    for path, Ik in zip(filament_paths, currents):
        Ik_c = float(Ik)
        for p1, p2 in path:
            p1 = [float(c) for c in p1]
            p2 = [float(c) for c in p2]
            t = [p2[i] - p1[i] for i in range(3)]
            L_seg = math.sqrt(sum(ti * ti for ti in t))
            if L_seg < 1e-14:
                continue
            t_hat = [ti / L_seg for ti in t]
            dx1 = x_cf - p1[0]; dy1 = y_cf - p1[1]; dz1 = z_cf - p1[2]
            dx2 = x_cf - p2[0]; dy2 = y_cf - p2[1]; dz2 = z_cf - p2[2]
            r1 = sqrt(dx1 * dx1 + dy1 * dy1 + dz1 * dz1 + 1e-24)
            r2 = sqrt(dx2 * dx2 + dy2 * dy2 + dz2 * dz2 + 1e-24)
            a_sc = dx1 * t_hat[0] + dy1 * t_hat[1] + dz1 * t_hat[2]
            b_sc = dx2 * t_hat[0] + dy2 * t_hat[1] + dz2 * t_hat[2]
            coeff = (MU_0 * Ik_c / (4.0 * math.pi)) * \
                log((r1 + a_sc + 1e-24) / (r2 + b_sc + 1e-24))
            Ax = Ax + coeff * t_hat[0]
            Ay = Ay + coeff * t_hat[1]
            Az = Az + coeff * t_hat[2]
    return CoefficientFunction((Ax, Ay, Az))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nw", type=int, default=3)
    ap.add_argument("--nh", type=int, default=3)
    ap.add_argument("--n-arc", type=int, default=30)
    ap.add_argument("--R-sphere", type=float, default=0.30,
                    help="Radius of the air sphere containing the coil.")
    ap.add_argument("--maxh", type=float, default=0.015)
    ap.add_argument("--order", type=int, default=1)
    args = ap.parse_args()

    R_coil = 0.030
    a_coil = 0.003
    gap_deg = 5.0
    arc_deg = 360.0 - gap_deg
    I_total = 1.0

    L_analytical = MU_0 * R_coil * (math.log(8 * R_coil / a_coil) - 2) \
        * (1.0 - gap_deg / 360.0)

    print("=" * 60)
    print("Step 0: PEEC Biot-Savart A_s + NGSolve energy integral")
    print("=" * 60)
    print(f"Coil : R={R_coil*1e3:.0f} mm, a={a_coil*1e3:.1f} mm, "
          f"gap={gap_deg} deg, nw x nh = {args.nw} x {args.nh}")
    print(f"Sphere: R={args.R_sphere*1e3:.0f} mm, maxh={args.maxh*1e3:.1f} mm, "
          f"HCurl order={args.order}")
    print(f"L_analytical (torus) = {L_analytical*1e9:.2f} nH")
    print()

    # Build filament bundle
    print("[1/4] PEEC filament bundle...")
    t0 = time.perf_counter()
    from coil_builder import CoilBuilder
    cb = (CoilBuilder(current=I_total)
          .set_start([R_coil, 0, 0],
                     orientation=np.array([[1, 0, 0],
                                            [0, 1, 0],
                                            [0, 0, 1]]))
          .set_cross_section(width=2 * a_coil, height=2 * a_coil)
          .add_arc(radius=R_coil, arc_angle=arc_deg, tilt=0))
    paths, _ = cb.to_filaments(args.nw, args.nh, frequency=0.0,
                                sigma=5.8e7, n_arc=args.n_arc)
    currents = [I_total / len(paths)] * len(paths)
    print(f"  {len(paths)} filaments ({time.perf_counter() - t0:.1f}s)")

    # Air sphere mesh
    print("[2/4] Air sphere mesh...")
    t0 = time.perf_counter()
    sphere = Sphere(Pnt(0, 0, 0), args.R_sphere)
    sphere.name = "air"
    geo = OCCGeometry(sphere)
    with TaskManager():
        ngmesh = geo.GenerateMesh(maxh=args.maxh, grading=0.5)
    mesh = Mesh(ngmesh)
    mesh.Curve(args.order + 1)
    print(f"  {mesh.ne} tets, {mesh.nv} vertices "
          f"({time.perf_counter() - t0:.1f}s)")

    # A_s CF
    print("[3/4] Build A_s CF, project onto HCurl GridFunction...")
    t0 = time.perf_counter()
    A_s_cf = biot_savart_A_cf(paths, currents, x, y, z)
    fes = HCurl(mesh, order=args.order)
    gfu_As = GridFunction(fes)
    with TaskManager():
        gfu_As.Set(A_s_cf, dual=False)
    print(f"  project: {time.perf_counter() - t0:.1f}s, "
          f"ndof = {fes.ndof}")

    # Energy
    print("[4/4] Integrate 0.5 nu_0 |curl(A_s)|^2 ...")
    t0 = time.perf_counter()
    with TaskManager():
        W = Integrate(0.5 * NU_0 * InnerProduct(curl(gfu_As), curl(gfu_As)),
                      mesh, order=8).real
    print(f"  W = {W:.4e} J  ({time.perf_counter() - t0:.1f}s)")
    L = 2.0 * W / I_total ** 2

    print()
    print("=" * 60)
    print(f"  L (Step 0, FEM energy)     = {L * 1e9:.2f} nH")
    print(f"  L (analytical torus)       = {L_analytical * 1e9:.2f} nH")
    print(f"  rel err                    = "
          f"{(L / L_analytical - 1) * 100:+.2f}%")
    print("=" * 60)


if __name__ == "__main__":
    main()
