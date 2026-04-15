"""Coil_3D_A_HCurl_PEEC_source.py

Step 1 of Phase 2 validation (2026-04-15):

Port the 3D HCurl A-formulation + Kelvin baseline (
Coil_3D_A_HCurl_with_Kelvin.py) to use a PEEC filament bundle as the
source, instead of a volume J over a meshed coil region.

Target physics
--------------
In the reduced A-formulation, A_total = A_s + A_r where
  A_s  : Biot-Savart field from the filament bundle (analytical).
  A_r  : the FEM unknown, corrects for the material modulation
         (nu != nu_0 in the Kelvin shell).

Weak form
---------
  Find A_r in HCurl such that for all v in HCurl:

    a(A_r, v) = f(v)

    a(u, v) = int nu(r) curl(u) . curl(v) dV     (same as baseline)

    f(v)    = - int (nu(r) - nu_0) curl(A_s(r)) . curl(v) dV

The RHS f(v) is non-zero only where nu(r) != nu_0, i.e., the Kelvin
shell (Sugahara 2022: nu_kelvin = nu_0 (R_K / rho')^2). In the inner
"air" region, nu = nu_0, so the RHS vanishes and A_r = 0 there (exactly,
in the continuum limit). The air-only FEM contribution comes purely
from A_s, which is the analytical Biot-Savart.

Kelvin shell evaluation of A_s uses the Phase 2 1-form pullback
(kelvin_source.kelvin_pullback_vector):

    A_s_shell(r') = (R/rho')^2 [A_s(r_phys) - 2 (A_s.n') n']

with r_phys = Kelvin_inverse(r') = c + R^2 (r' - c) / |r' - c|^2.

Inductance
----------
  L = 2 W / I_total^2
  W = 0.5 * int nu |curl(A_total)|^2 dV

Comparison target
-----------------
The baseline (Coil_3D_A_HCurl_with_Kelvin.py) reports
L_analytical = mu_0 * R_coil * (ln(8 R_coil / a_coil) - 2) * (1 - gap/360)
and L_FEM ~ L_analytical within 1-2%.

For this script the coil cross-section is approximated by a
PEEC filament bundle with (nw, nh) filaments; single-filament
(nw=nh=1) gives zero-cross-section (diverges at self), so we
use nw=nh=3 or similar.

Run
---
    python Coil_3D_A_HCurl_PEEC_source.py [--nw 3 --nh 3]

At time of writing the goal is to verify that L_FEM from the PEEC
source agrees with L_FEM from the baseline volume-J source (both
modes use the same mesh).
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

import radia      # noqa: F401  -- MKL DLL paths

from ngsolve import (HCurl, BilinearForm, LinearForm, GridFunction,
                      CoefficientFunction, CF, Integrate, InnerProduct,
                      IfPos, Conj, curl, dx, x, y, z, sqrt, log,
                      Periodic, TaskManager, Mesh)
from netgen.occ import (WorkPlane, Axes, Pnt, Dir, Sphere, Axis, Vertex,
                         Glue, OCCGeometry, IdentificationType)


# =======================================================================
# Parameters (match the baseline)
# =======================================================================

R_coil = 0.030
a_coil = 0.003
gap_deg = 5.0
arc_deg = 360.0 - gap_deg
R_K = 0.06
offset = 0.15
maxh = 0.012
order = 1

MU_0 = 4e-7 * math.pi
NU_0 = 1.0 / MU_0
I_total = 1.0

L_analytical_torus = MU_0 * R_coil * (math.log(8 * R_coil / a_coil) - 2) \
    * (1.0 - gap_deg / 360.0)


def _biot_savart_A_cf(filament_paths, currents, x_cf, y_cf, z_cf):
    """Build an NGSolve CoefficientFunction for A_s(r) via finite-segment
    analytical Biot-Savart, evaluated at (x_cf, y_cf, z_cf).

    For a finite segment [p1, p2] with current I:
        A_seg(r) = (mu_0 I / 4 pi) * t_hat *
                   log( (r1 + a) / (r2 + b) )

    where
        t_hat = (p2 - p1) / |p2 - p1|
        r_i   = |r - p_i|
        a     = (r - p1) . t_hat
        b     = (r - p2) . t_hat

    This closed form comes from asinh(a/c) - asinh(b/c) simplified to
    the ratio form when the perpendicular distance c is kept implicit
    in the r_i norms. Numerically well-behaved away from the filament.
    """
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

            r1_sq = dx1 * dx1 + dy1 * dy1 + dz1 * dz1
            r2_sq = dx2 * dx2 + dy2 * dy2 + dz2 * dz2
            r1 = sqrt(r1_sq + 1e-24)
            r2 = sqrt(r2_sq + 1e-24)

            a_scal = (dx1 * t_hat[0] + dy1 * t_hat[1] + dz1 * t_hat[2])
            b_scal = (dx2 * t_hat[0] + dy2 * t_hat[1] + dz2 * t_hat[2])

            num = r1 + a_scal
            den = r2 + b_scal
            # Stabilize: shift both by a small absolute bias so ratio
            # never blows up on a segment endpoint-aligned observation.
            coeff = (MU_0 * Ik_c / (4.0 * math.pi)) * \
                log((num + 1e-24) / (den + 1e-24))
            Ax = Ax + coeff * t_hat[0]
            Ay = Ay + coeff * t_hat[1]
            Az = Az + coeff * t_hat[2]
    return CoefficientFunction((Ax, Ay, Az))


def kelvin_pullback_cf(A_inner_cf, A_kelvin_cf, mesh_materials,
                       kelvin_mats, offset_xyz, R):
    """Material-switched A_s CF: A_inner in non-Kelvin, pullback in Kelvin.

    A_inner_cf  : CF for A_s evaluated at physical coords (x, y, z).
    A_kelvin_cf : CF for A_s evaluated at Kelvin-mapped physical coords
                  (kelvin_x, kelvin_y, kelvin_z).
    kelvin_mats : list of material-name substrings to treat as Kelvin.
    offset_xyz  : (3,) Kelvin sphere center.
    R           : Kelvin sphere radius.

    Returns CF equal to:
        in Kelvin mat: (R/rho')^2 * (A_kelvin - 2 (A_kelvin . n') n')
        else         : A_inner
    """
    ox, oy, oz = offset_xyz
    dxp = x - ox
    dyp = y - oy
    dzp = z - oz
    rho_p_sq = dxp * dxp + dyp * dyp + dzp * dzp + 1e-24
    rho_p = sqrt(rho_p_sq)
    n_cf = CoefficientFunction((dxp / rho_p, dyp / rho_p, dzp / rho_p))
    A_dot_n = (A_kelvin_cf[0] * n_cf[0]
               + A_kelvin_cf[1] * n_cf[1]
               + A_kelvin_cf[2] * n_cf[2])
    refl_x = A_kelvin_cf[0] - 2.0 * A_dot_n * n_cf[0]
    refl_y = A_kelvin_cf[1] - 2.0 * A_dot_n * n_cf[1]
    refl_z = A_kelvin_cf[2] - 2.0 * A_dot_n * n_cf[2]
    factor = (R * R) / rho_p_sq
    A_shell = CoefficientFunction(
        (factor * refl_x, factor * refl_y, factor * refl_z))

    choice = []
    for m in mesh_materials:
        in_kelvin = any(kw in m.lower() for kw in kelvin_mats)
        choice.append(A_shell if in_kelvin else A_inner_cf)
    # NGSolve MaterialCF expects a CF per material; combine via list
    # order matching mesh.GetMaterials().
    return CoefficientFunction(choice)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nw", type=int, default=3)
    ap.add_argument("--nh", type=int, default=3)
    ap.add_argument("--n-arc", type=int, default=30)
    args = ap.parse_args()

    print("=" * 60)
    print("3D HCurl A-formulation + Kelvin, PEEC filament source")
    print("=" * 60)
    print(f"Coil : R={R_coil*1e3:.0f} mm, a={a_coil*1e3:.1f} mm, "
          f"gap={gap_deg} deg, nw x nh = {args.nw} x {args.nh}")
    print(f"Kelvin: R_K={R_K*1e3:.0f} mm, offset={offset*1e3:.0f} mm")
    print(f"L_analytical (torus)= {L_analytical_torus*1e9:.2f} nH")
    print()

    # === Build filament bundle via CoilBuilder ===
    print("[1/5] Building PEEC filament bundle...")
    t0 = time.perf_counter()
    from radia_coil_builder import CoilBuilder
    cb = (CoilBuilder(current=I_total)
          .set_start([R_coil, 0, 0],
                     orientation=np.array([[1, 0, 0],
                                            [0, 1, 0],
                                            [0, 0, 1]]))
          .set_cross_section(width=2 * a_coil, height=2 * a_coil)
          .add_arc(radius=R_coil, arc_angle=arc_deg, tilt=0))
    paths, _ = cb.to_filaments(args.nw, args.nh, frequency=0.0,
                                 sigma=5.8e7, n_arc=args.n_arc)
    dw = (2 * a_coil) / args.nw
    dh = (2 * a_coil) / args.nh
    # Use DC currents: each filament carries I_total / (nw*nh).
    currents = [I_total / (args.nw * args.nh)] * len(paths)
    print(f"  {len(paths)} filaments ({time.perf_counter() - t0:.1f}s)")

    # === Geometry: inner air sphere + Kelvin shell (no coil volume) ===
    print("[2/5] Building geometry (NO coil volume)...")
    t0 = time.perf_counter()
    inner_sphere = Sphere(Pnt(0, 0, 0), R_K)
    inner_sphere.maxh = maxh
    for f in inner_sphere.faces:
        f.name = "kelvin_int"
    inner_sphere.name = "air"

    outer_sphere = Sphere(Pnt(offset, 0, 0), R_K)
    outer_sphere.maxh = maxh * 2
    for f in outer_sphere.faces:
        f.name = "kelvin_ext"
    outer_sphere.name = "kelvin"

    gnd = Vertex(Pnt(offset, 0, 0))
    gnd.name = "GND"

    # Periodic identification
    k_int_face = None; k_ext_face = None
    for f in inner_sphere.faces:
        if f.name == "kelvin_int":
            k_int_face = f; break
    for f in outer_sphere.faces:
        if f.name == "kelvin_ext":
            k_ext_face = f; break
    if k_int_face and k_ext_face:
        k_int_face.Identify(k_ext_face, "periodic",
                             IdentificationType.PERIODIC)

    shape = Glue([inner_sphere, outer_sphere, gnd])
    geo = OCCGeometry(shape)
    with TaskManager():
        ngmesh = geo.GenerateMesh(maxh=maxh, grading=0.5)
    mesh = Mesh(ngmesh)
    mesh.Curve(order + 1)
    t_mesh = time.perf_counter() - t0
    mats = mesh.GetMaterials()
    print(f"  Mesh: {mesh.ne} tets, {mesh.nv} vertices ({t_mesh:.1f}s)")
    print(f"  Materials: {mats}")

    # === Build A_s CFs ===
    print("[3/5] Building A_s coefficient functions...")
    t0 = time.perf_counter()
    # Inner (physical) evaluation at (x, y, z):
    A_inner = _biot_savart_A_cf(paths, currents, x, y, z)
    # Kelvin-mapped physical coords:
    dxp = x - offset
    dyp = y
    dzp = z
    rho_p_sq = dxp * dxp + dyp * dyp + dzp * dzp + 1e-24
    kel_x = offset + (R_K * R_K / rho_p_sq) * dxp
    kel_y = (R_K * R_K / rho_p_sq) * dyp
    kel_z = (R_K * R_K / rho_p_sq) * dzp
    # Evaluate A_s at the Kelvin-mapped coord (physical exterior):
    A_kel_phys = _biot_savart_A_cf(paths, currents, kel_x, kel_y, kel_z)
    # Apply Phase 2 pullback to Kelvin shell, leave inner alone:
    A_s_cf = kelvin_pullback_cf(A_inner, A_kel_phys, mats,
                                 kelvin_mats=["kelvin"],
                                 offset_xyz=(offset, 0.0, 0.0), R=R_K)
    print(f"  A_s CF built ({time.perf_counter() - t0:.1f}s)")

    # === Reluctivity with Kelvin modulation ===
    kelvin_fac = R_K ** 2 / rho_p_sq
    nu_dict = {}
    for m in mats:
        if "kelvin" in m.lower():
            nu_dict[m] = NU_0 * kelvin_fac
        else:
            nu_dict[m] = NU_0
    nu_cf = mesh.MaterialCF(nu_dict, default=NU_0)

    # === FE space ===
    print("[4/5] FE space + assembly...")
    fes_base = HCurl(mesh, order=order, dirichlet_bbnd="GND")
    fes = Periodic(fes_base)
    u, v = fes.TnT()
    print(f"  DOFs: {fes.ndof}")

    a_bf = BilinearForm(fes)
    a_bf += nu_cf * curl(u) * curl(v) * dx(bonus_intorder=4)
    a_bf += 1e-8 * NU_0 * u * v * dx   # gauge regularization

    f_lf = LinearForm(fes)
    # Reduced-A RHS: - (nu - nu_0) * curl(A_s) . curl(v) over whole mesh
    # (non-zero only where nu != nu_0).
    f_lf += -(nu_cf - NU_0) * InnerProduct(curl(A_s_cf), curl(v)) \
        * dx(bonus_intorder=4)

    with TaskManager():
        a_bf.Assemble()
        f_lf.Assemble()
    print("  Assembly done.")

    # === Solve ===
    print("[5/5] Solve...")
    t0 = time.perf_counter()
    gfu = GridFunction(fes)
    with TaskManager():
        gfu.vec.data = a_bf.mat.Inverse(fes.FreeDofs(),
                                         inverse="pardiso") * f_lf.vec
    t_solve = time.perf_counter() - t0
    print(f"  solve: {t_solve:.1f}s")

    # === Inductance from total A = A_s + A_r ===
    A_total = A_s_cf + gfu
    W = Integrate(0.5 * nu_cf * InnerProduct(curl(A_total),
                                              Conj(curl(A_total))),
                  mesh, order=10).real
    L = 2.0 * W / I_total ** 2

    print()
    print("=" * 60)
    print(f"  L (PEEC source, reduced-A) = {L * 1e9:.2f} nH")
    print(f"  L (analytical torus)       = {L_analytical_torus * 1e9:.2f} nH")
    print(f"  rel err                    = "
          f"{(L / L_analytical_torus - 1) * 100:+.2f}%")
    print("=" * 60)


if __name__ == "__main__":
    main()
