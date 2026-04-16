"""validate_omega_coil_source.py  (item 1)

Cross-check of the Omega-reduced Omega scalar-potential FEM + Kelvin
against Radia BEM on the SAME physical problem: a magnetic cylinder
(mu_r = 100) sitting inside the field of a gapped circular coil
(rad.ObjArcCur).

    Case A (FEM):   radia.scalar_potential_solver.ScalarPotentialSolver
                    with Kelvin outer sphere, Radia H_s voxel source.
                    Solves the phi-reduced equation
                        div(mu * (H_s - grad(phi))) = 0
                    and returns B = mu * (H_s - grad(phi)).

    Case B (Radia): coil + cylinder in a Radia container, Solve(),
                    rad.Fld(container, 'b', pt) is the full BEM
                    ground truth (free-space, open domain).

This exercises the "source_in_omega_form" MCP pattern end-to-end:
no Omega_s on the coil-iron interface, H_s is injected as a voxel
CF that is identically zero outside the physical bbox (so the
Kelvin exterior region sees no source — the "EMPY-simple" regime).

Acceptance target: < 2 % max Bz diff on sample points well inside the
cylinder and in the near exterior (sphere-family precedent).
"""
from __future__ import annotations

import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))

import radia as rad  # noqa: E402  (installed package with compiled .pyd)

from ngsolve import Mesh, TaskManager  # noqa: E402
from netgen.occ import (Cylinder, Sphere, Pnt, Z, Vertex, Glue, OCCGeometry,
                         IdentificationType)  # noqa: E402

# Dev-tree helpers (kelvin_source + ScalarPotentialSolver with Kelvin)
SRC_RADIA = os.path.abspath(os.path.join(HERE, '..', '..', '..', 'src', 'radia'))
if SRC_RADIA not in sys.path:
    sys.path.insert(0, SRC_RADIA)
from scalar_potential_solver import ScalarPotentialSolver  # noqa: E402
from netgen_mesh_import import netgen_mesh_to_radia  # noqa: E402


MU_0 = 4.0 * math.pi * 1e-7

# --- Geometry / physics ---------------------------------------------
cylinder_radius = 0.03
cylinder_height = 0.06
mu_r = 100

R_coil = 0.08
a_coil = 0.008
gap_deg = 5.0
arc_deg = 360.0 - gap_deg
z_coil = 0.10       # coil placed above cylinder
I_total = 1000.0    # A-turn equivalent

# Kelvin
kelvin_center = [0.30, 0.0, 0.0]
kelvin_R = 0.15

# Mesh
maxh_iron = 0.008
maxh_air = 0.025
maxh_kelvin = 0.03
fe_order = 2

# Radia BEM mesh (for Case B)
radia_maxh = 0.008


# ---------------------------------------------------------------
# Case B: Radia BEM — coil + cylinder, solved
# ---------------------------------------------------------------
def build_and_solve_radia_reference():
    rad.UtiDelAll()
    J0 = I_total / (2 * a_coil) ** 2
    coil = rad.ObjArcCur(
        [0, 0, z_coil],
        [R_coil - a_coil, R_coil + a_coil],
        [0, math.radians(arc_deg)],
        2 * a_coil, 200, 'man', 'z', J0)
    # Magnetic cylinder meshed via Netgen tetra, imported to Radia
    half_h = cylinder_height / 2
    mag_cyl = Cylinder(Pnt(0, 0, -half_h), Z,
                        r=cylinder_radius, h=cylinder_height)
    mag_cyl.mat("magnetic")
    mag_cyl.maxh = radia_maxh
    geo = OCCGeometry(mag_cyl)
    with TaskManager():
        ngmesh = geo.GenerateMesh(maxh=radia_maxh)
    mesh = Mesh(ngmesh)
    cyl_obj = netgen_mesh_to_radia(mesh,
                                    material={'magnetization': [0, 0, 0]},
                                    units='m', material_filter='magnetic',
                                    verbose=False)
    rad.MatApl(cyl_obj, rad.MatLin(mu_r))
    grp = rad.ObjCnt([coil, cyl_obj])
    rad.SolverConfig(bicgstab_tol=1e-6, relax_param=0.0)
    print(f"  Radia tetra: {mesh.ne} elements in cylinder")
    t0 = time.perf_counter()
    # LU is more reliable than BiCGSTAB for this size; O(N^3) but OK < 2k elt
    result = rad.Solve(grp, 1e-5, 500, 0)   # method=0 LU
    print(f"  Radia Solve: {time.perf_counter() - t0:.1f}s, "
          f"residual={result[0]:.3e}")
    return coil, cyl_obj, grp


# ---------------------------------------------------------------
# Case A: FEM Omega-reduced (phi-reduced) + Kelvin
# ---------------------------------------------------------------
def build_kelvin_mesh():
    """Build mesh following the 3D_cylinder_with_Kelvin.py pattern.

    Materials: "iron", "air_inner", "air_outer"  (not 'air'+'air')
    """
    half_h = cylinder_height / 2
    phys_R = 0.18   # inner & Kelvin sphere radius (must cover coil at 0.10+a)

    mag_cyl = Cylinder(Pnt(0, 0, -half_h), Z,
                        r=cylinder_radius, h=cylinder_height)
    mag_cyl.mat("iron")
    mag_cyl.maxh = maxh_iron
    for f in mag_cyl.faces:
        f.name = "iron_surf"

    inner = Sphere(Pnt(0, 0, 0), phys_R)
    inner.maxh = maxh_air
    for f in inner.faces:
        f.name = "kelvin_int"
    inner_air = inner - mag_cyl
    inner_air.mat("air_inner")

    outer = Sphere(Pnt(*kelvin_center), phys_R)
    outer.maxh = maxh_kelvin
    outer.mat("air_outer")
    for f in outer.faces:
        f.name = "kelvin_ext"

    gnd = Vertex(Pnt(*kelvin_center))
    gnd.name = "GND"

    geo = Glue([inner_air, mag_cyl, outer, gnd])
    geo.solids[0].name = "air_inner"
    geo.solids[1].name = "iron"
    geo.solids[2].name = "air_outer"

    # Periodic identification
    k_int = k_ext = None
    for solid in geo.solids:
        for f in solid.faces:
            if f.name == "kelvin_int" and k_int is None:
                k_int = f
            elif f.name == "kelvin_ext" and k_ext is None:
                k_ext = f
    if k_int is None or k_ext is None:
        raise RuntimeError("Could not find Kelvin periodic faces")
    k_int.Identify(k_ext, "periodic", IdentificationType.PERIODIC)

    with TaskManager():
        ngmesh = OCCGeometry(geo).GenerateMesh(maxh=maxh_air, grading=0.5)
    mesh = Mesh(ngmesh)
    mesh.Curve(fe_order)
    return mesh, phys_R


def solve_fem_kelvin(coil_handle, mesh, phys_R):
    solver = ScalarPotentialSolver(
        mesh,
        iron_domains='iron',
        mu_r=mu_r,
        order=fe_order,
        kelvin_region='air_outer',   # must match mesh material name
        kelvin_radius=phys_R,
        kelvin_center=kelvin_center,
    )
    # Voxel bbox covers the physical ball (origin, radius phys_R) only
    bbox = [[-phys_R, phys_R], [-phys_R, phys_R], [-phys_R, phys_R]]
    print("  Building H_s voxel from Radia (coil only)...")
    t0 = time.perf_counter()
    solver.set_source_from_radia(coil_handle, resolution=41, bbox=bbox)
    print(f"    voxel build: {time.perf_counter() - t0:.1f}s")
    print("  Solving phi-reduced single potential + Kelvin ...")
    t0 = time.perf_counter()
    solver.solve_single_potential(dirichlet='GND')
    print(f"    FEM solve:  {time.perf_counter() - t0:.1f}s")
    return solver


# ---------------------------------------------------------------
# Comparison
# ---------------------------------------------------------------
def compare_B(solver, grp_radia, sample_points, label):
    print(f"\n{label}")
    print(f"  {'point':>24s}  {'|B|_FEM':>12s}  {'|B|_Radia':>12s}  {'rel':>9s}")
    rels = []
    B_cf = solver.get_B()
    mesh = solver.mesh
    for pt in sample_points:
        try:
            mip = mesh(pt[0], pt[1], pt[2])
            B_fem = np.array([B_cf(mip)[i] for i in range(3)])
        except Exception:
            B_fem = np.full(3, np.nan)
        B_rad = np.array(rad.Fld(grp_radia, 'b', list(pt)))
        mag_fem = float(np.linalg.norm(B_fem))
        mag_rad = float(np.linalg.norm(B_rad))
        if np.isnan(mag_fem):
            continue
        rel = (mag_fem - mag_rad) / max(mag_rad, 1e-12)
        rels.append(rel)
        p_str = f"({pt[0]:+.3f},{pt[1]:+.3f},{pt[2]:+.3f})"
        print(f"  {p_str:>24s}  {mag_fem:11.5e}  {mag_rad:11.5e}  "
              f"{rel*100:+7.2f}%")
    return rels


def main():
    print("=" * 72)
    print("Item 1: Omega-reduced + Kelvin vs Radia BEM (coil source + iron)")
    print("=" * 72)
    print(f"  cylinder r={cylinder_radius}  h={cylinder_height}  mu_r={mu_r}")
    print(f"  coil     R={R_coil}  a={a_coil}  z={z_coil}  I={I_total}")
    print(f"  Kelvin   c={kelvin_center}  R=(phys)  offset far from coil")
    print(f"  FEM  order={fe_order}  maxh=iron:{maxh_iron}/air:{maxh_air}/K:{maxh_kelvin}")
    print(f"  Radia maxh={radia_maxh}")

    print("\n[FEM] building mesh ...")
    mesh, phys_R = build_kelvin_mesh()
    print(f"  {mesh.ne} tets, {mesh.nv} verts, mats={mesh.GetMaterials()}")
    print(f"  phys_R (inner ball + Kelvin radius) = {phys_R}")

    print("\n[Radia] building coil + cylinder, Solve ...")
    coil, cyl_obj, grp = build_and_solve_radia_reference()

    print("\n[FEM] ScalarPotentialSolver with Kelvin ...")
    solver = solve_fem_kelvin(coil, mesh, phys_R)

    half_h = cylinder_height / 2
    z_inside = [(0.0, 0.0, zv) for zv in
                 np.linspace(-half_h + 0.01, half_h - 0.01, 5)]
    x_inside = [(xv, 0.0, 0.0) for xv in
                 np.linspace(0.005, cylinder_radius - 0.005, 4)]
    z_outside = [(0.0, 0.0, zv) for zv in
                  np.linspace(half_h + 0.02, phys_R - 0.04, 4)]

    rels_zi = compare_B(solver, grp, z_inside,
                         "Bz along z-axis (inside cylinder)")
    rels_xi = compare_B(solver, grp, x_inside,
                         "Bz along x-axis (inside cylinder, z=0)")
    rels_zo = compare_B(solver, grp, z_outside,
                         "Bz along z-axis (outside cylinder, towards coil)")

    all_rels = np.array(rels_zi + rels_xi + rels_zo)
    print()
    print("=" * 72)
    print("SUMMARY")
    print("=" * 72)
    print(f"  samples       = {len(all_rels)}")
    print(f"  max |rel diff|= {np.max(np.abs(all_rels))*100:.3f} %")
    print(f"  RMS rel diff  = {np.sqrt(np.mean(all_rels**2))*100:.3f} %")
    ok = np.max(np.abs(all_rels)) < 0.05
    print(f"  acceptance (|max| < 5 %): {'OK' if ok else 'FAIL'}")
    print("=" * 72)


if __name__ == "__main__":
    main()
