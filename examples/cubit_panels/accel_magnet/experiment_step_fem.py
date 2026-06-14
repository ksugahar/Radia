#!/usr/bin/env python
"""Test STEP -> OCC -> Netgen -> FEM solve for accelerator dipole.

Compares FEM (Omega-reduced, linear mu_r=1000) against MSC reference (-976 mT).

Usage:
    python test_step_fem.py [--mu-r 1000] [--order 2] [--mesh-size 0.008]
"""

import sys
import os
import math
import time
import argparse

repo_root = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         '..', '..', '..')
sys.path.insert(0, os.path.join(repo_root, 'src'))
sys.path.insert(0, os.path.join(repo_root, 'src', 'radia'))

import numpy as np
import radia as rad


# ================================================================
# Constants
# ================================================================
MU_0 = 4e-7 * math.pi
NU_0 = 1.0 / MU_0
MSC_REFERENCE = -976e-3  # MSC Bz at origin [T]


# ================================================================
# Coil (same as test_occ_dipole.py)
# ================================================================
def create_coil(current_at=20000.0):
    """Create racetrack coil via CoilBuilder."""
    from coil_builder import CoilBuilder

    mm = 1e-3
    coil = (CoilBuilder(current=current_at)
        .set_start([47.5*mm, 100*mm, 0])
        .set_cross_section(width=35*mm, height=105*mm)
        .add_straight(62.5*mm)
        .add_arc(radius=22.5*mm, arc_angle=90)
        .add_straight(50*mm)
        .add_arc(radius=22.5*mm, arc_angle=90)
        .add_straight(62.5*mm)
        .add_arc(radius=22.5*mm, arc_angle=90)
        .add_straight(50*mm)
        .add_arc(radius=22.5*mm, arc_angle=90))

    rad.UtiDelAll()
    objs = coil.to_radia()
    container = rad.ObjCnt(objs)
    return container, coil


# ================================================================
# FEM solve (Omega-reduced, linear)
# ================================================================
def solve_fem(mesh, info, mu_r, current_at, fes_order=1):
    """Omega-reduced scalar potential solve with Kelvin."""
    from ngsolve import (Mesh, H1, L2, Periodic, GridFunction,
                         BilinearForm, LinearForm, grad, sqrt, dx,
                         InnerProduct, Integrate, CF, TaskManager, x, y, z)

    t0 = time.perf_counter()

    # Coil source field
    coil_container, _ = create_coil(current_at)
    H_s = rad.RadiaField(coil_container, 'h')

    materials = set(mesh.GetMaterials())
    boundaries = set(mesh.GetBoundaries())
    has_kelvin = any("kelvin" in m for m in materials)

    print(f"Materials: {materials}")
    print(f"Boundaries: {boundaries}")
    print(f"Has Kelvin: {has_kelvin}")

    # Kelvin weight CF: (a/r')^2 in kelvin domain, 1 elsewhere
    kc = info['kelvin_center']
    a_kelvin = info['kelvin_radius']
    dx_k = x - kc[0]
    dy_k = y - kc[1]
    dz_k = z - kc[2]
    rp_sq = dx_k*dx_k + dy_k*dy_k + dz_k*dz_k + 1e-20
    kelvin_fac = a_kelvin**2 / rp_sq

    kw_dict = {}
    for m in materials:
        if "kelvin" in m.lower():
            kw_dict[m] = kelvin_fac
        else:
            kw_dict[m] = CF(1.0)
    kelvin_weight = mesh.MaterialCF(kw_dict, default=CF(1.0))

    # Material CF: mu for omega-reduced
    mu_dict = {}
    for m in materials:
        if "yoke" in m.lower():
            mu_dict[m] = CF(MU_0 * mu_r)
        else:
            mu_dict[m] = CF(MU_0)
    mu_cf = mesh.MaterialCF(mu_dict, default=CF(MU_0))

    # Combined coefficient
    coeff = mu_cf * kelvin_weight

    # FE space
    dir_parts = []
    if "kelvin_ext" in boundaries:
        dir_parts.append("kelvin_ext")
    if "outer" in boundaries:
        dir_parts.append("outer")
    if "sym_normal" in boundaries:
        dir_parts.append("sym_normal")
    dirichlet_bnd = "|".join(dir_parts) if dir_parts else ""
    print(f"Dirichlet BC: {dirichlet_bnd}")

    base_fes = H1(mesh, order=fes_order, dirichlet=dirichlet_bnd)
    if has_kelvin:
        fes = Periodic(base_fes)
        print("FES: Periodic H1 (Omega-reduced)")
    else:
        fes = base_fes
        print("FES: H1 (Omega-reduced)")

    u, v = fes.TnT()
    ndof = fes.ndof
    print(f"DOFs: {ndof}")

    # Bilinear form
    a_bf = BilinearForm(fes)
    a_bf += coeff * grad(u) * grad(v) * dx(bonus_intorder=4)

    # Linear form
    f_lf = LinearForm(fes)
    f_lf += coeff * H_s * grad(v) * dx(bonus_intorder=4)

    # Solve
    print("Assembling...")
    with TaskManager():
        a_bf.Assemble()
        f_lf.Assemble()

    gfu = GridFunction(fes)
    print("Solving...")
    t_solve = time.perf_counter()

    with TaskManager():
        if ndof < 200000:
            gfu.vec.data = a_bf.mat.Inverse(
                fes.FreeDofs(), inverse="pardiso") * f_lf.vec
        else:
            from ngsolve import solvers, Preconditioner
            pre = Preconditioner(a_bf, "bddc")
            a_bf.Assemble()
            solvers.BVP(bf=a_bf, lf=f_lf, gf=gfu, pre=pre,
                        maxsteps=500, tol=1e-8)

    t_solve = time.perf_counter() - t_solve
    print(f"Solve time: {t_solve:.1f}s")

    # Post-process
    H_total = H_s - grad(gfu)
    B_field = mu_cf * H_total

    # B at origin
    try:
        B_origin = [float(v) for v in B_field(mesh(0, 0, 0))]
        Bz = B_origin[2]  # z-component = field/gap direction
        B_mag = math.sqrt(sum(b**2 for b in B_origin))
    except Exception as e:
        print(f"WARNING: B evaluation at origin failed: {e}")
        B_origin = [0, 0, 0]
        Bz = 0
        B_mag = 0

    t_total = time.perf_counter() - t0

    print(f"\n{'='*50}")
    print(f"RESULT:")
    print(f"  B at origin: [{B_origin[0]*1e3:.1f}, {B_origin[1]*1e3:.1f}, "
          f"{B_origin[2]*1e3:.1f}] mT")
    print(f"  |B|: {B_mag*1e3:.1f} mT")
    print(f"  Bz (gap): {Bz*1e3:.1f} mT")
    print(f"  MSC ref: {MSC_REFERENCE*1e3:.1f} mT")
    print(f"  Error: {(Bz - MSC_REFERENCE)/abs(MSC_REFERENCE)*100:.1f}%")
    print(f"  ndof: {ndof}, ne: {mesh.ne}")
    print(f"  Total time: {t_total:.1f}s")
    print(f"{'='*50}")

    return {
        'B_origin': B_origin,
        'B_mag': B_mag,
        'Bz': Bz,
        'error_pct': (Bz - MSC_REFERENCE) / abs(MSC_REFERENCE) * 100,
        'ndof': ndof,
        'ne': mesh.ne,
        't_solve': t_solve,
        't_total': t_total,
    }


# ================================================================
# Main
# ================================================================
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--mu-r', type=float, default=1000.0)
    parser.add_argument('--current', type=float, default=20000.0)
    parser.add_argument('--order', type=int, default=2,
                        help='Geometry curve order')
    parser.add_argument('--fes-order', type=int, default=1,
                        help='FE polynomial order')
    parser.add_argument('--mesh-size', type=float, default=None,
                        help='Yoke mesh size [m]')
    parser.add_argument('--kelvin-radius', type=float, default=None)
    parser.add_argument('--step', type=str, default=None,
                        help='STEP file (default: yoke.step)')
    args = parser.parse_args()

    step_file = args.step or os.path.join(
        os.path.dirname(os.path.abspath(__file__)), 'yoke.step')

    if not os.path.exists(step_file):
        print(f"ERROR: STEP file not found: {step_file}")
        sys.exit(1)

    from step_mesh_builder import build_mesh_from_step

    print(f"STEP: {step_file}")
    print(f"mu_r: {args.mu_r}")
    print(f"Current: {args.current} AT")
    print(f"Order: curve={args.order}, fes={args.fes_order}")

    mesh, info = build_mesh_from_step(
        step_file,
        symmetry="quarter_xz",
        curve_order=args.order,
        mesh_size=args.mesh_size,
        mesh_size_yoke=args.mesh_size,
        kelvin_radius=args.kelvin_radius,
    )

    result = solve_fem(mesh, info, args.mu_r, args.current,
                       fes_order=args.fes_order)
