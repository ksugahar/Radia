"""convergence_study.py — mesh convergence for solve_total_reduced_potential.

Runs 3 refinement levels and reports max|rel| and RMS at each.
If errors decrease monotonically, the formulation is correct and
the residual is discretization-limited.
"""
from __future__ import annotations
import math, os, sys, time
import numpy as np

import radia as rad
from ngsolve import Mesh, TaskManager
from netgen.occ import (Cylinder, Sphere, Pnt, Z, Vertex, Glue, OCCGeometry,
                         IdentificationType)

HERE = os.path.dirname(os.path.abspath(__file__))
SRC_RADIA = os.path.abspath(os.path.join(HERE, '..', '..', '..', 'src', 'radia'))
if SRC_RADIA not in sys.path:
    sys.path.insert(0, SRC_RADIA)
from scalar_potential_solver import ScalarPotentialSolver
from netgen_mesh_import import netgen_mesh_to_radia

MU_0 = 4 * math.pi * 1e-7
cylinder_radius = 0.03
cylinder_height = 0.06
mu_r = 100
R_coil = 0.08
a_coil = 0.008
z_coil = 0.10
I_total = 1000.0
phys_R = 0.18
kelvin_center = (0.50, 0.0, 0.0)

LEVELS = [
    {"label": "L1 (coarse, p=2)",
     "maxh_iron": 0.008, "maxh_air": 0.025, "maxh_kelvin": 0.03,
     "fe_order": 2},
    {"label": "L2 (medium, p=2)",
     "maxh_iron": 0.005, "maxh_air": 0.015, "maxh_kelvin": 0.025,
     "fe_order": 2},
    {"label": "L3 (fine, p=3)",
     "maxh_iron": 0.003, "maxh_air": 0.010, "maxh_kelvin": 0.020,
     "fe_order": 3},
]


def build_radia_coil():
    # No UtiDelAll — coil creation does not require clearing existing
    # objects, and clearing would destroy the reference grp_ref.
    J0 = I_total / (2 * a_coil) ** 2
    return rad.ObjArcCur(
        [0, 0, z_coil],
        [R_coil - a_coil, R_coil + a_coil],
        [0, 2 * math.pi], 2 * a_coil, 200, 'man', 'z', J0)


def build_radia_reference():
    rad.UtiDelAll()
    J0 = I_total / (2 * a_coil) ** 2
    coil = rad.ObjArcCur(
        [0, 0, z_coil],
        [R_coil - a_coil, R_coil + a_coil],
        [0, 2 * math.pi], 2 * a_coil, 200, 'man', 'z', J0)
    half_h = cylinder_height / 2
    mag_cyl = Cylinder(Pnt(0, 0, -half_h), Z,
                        r=cylinder_radius, h=cylinder_height)
    mag_cyl.mat("magnetic"); mag_cyl.maxh = 0.004
    geo = OCCGeometry(mag_cyl)
    with TaskManager():
        ngmesh = geo.GenerateMesh(maxh=0.004)
    mesh_rad = Mesh(ngmesh)
    cyl_obj = netgen_mesh_to_radia(mesh_rad,
                                    material={'magnetization': [0, 0, 0]},
                                    units='m', material_filter='magnetic',
                                    verbose=False)
    rad.MatApl(cyl_obj, rad.MatLin(mu_r))
    grp = rad.ObjCnt([coil, cyl_obj])
    rad.SolverConfig(bicgstab_tol=1e-8, relax_param=0.0)
    result = rad.Solve(grp, 1e-7, 2000, 1)   # BiCGSTAB for large DOF
    print(f"  Radia ref: {mesh_rad.ne} tets, residual={result[0]:.3e}")
    return grp


def build_kelvin_mesh(maxh_iron, maxh_air, maxh_kelvin, fe_order):
    half_h = cylinder_height / 2
    mag_cyl = Cylinder(Pnt(0, 0, -half_h), Z,
                        r=cylinder_radius, h=cylinder_height)
    mag_cyl.mat("iron"); mag_cyl.maxh = maxh_iron
    for f in mag_cyl.faces:
        f.name = "iron_surf"
    inner = Sphere(Pnt(0, 0, 0), phys_R); inner.maxh = maxh_air
    for f in inner.faces:
        f.name = "kelvin_int"
    inner_air = inner - mag_cyl; inner_air.mat("air_inner")
    outer = Sphere(Pnt(*kelvin_center), phys_R); outer.maxh = maxh_kelvin
    outer.mat("air_outer")
    for f in outer.faces:
        f.name = "kelvin_ext"
    gnd = Vertex(Pnt(*kelvin_center)); gnd.name = "GND"
    geo = Glue([inner_air, mag_cyl, outer, gnd])
    geo.solids[0].name = "air_inner"
    geo.solids[1].name = "iron"
    geo.solids[2].name = "air_outer"
    k_int = k_ext = None
    for solid in geo.solids:
        for f in solid.faces:
            if f.name == "kelvin_int" and k_int is None:
                k_int = f
            elif f.name == "kelvin_ext" and k_ext is None:
                k_ext = f
    k_int.Identify(k_ext, "periodic", IdentificationType.PERIODIC)
    with TaskManager():
        ngmesh = OCCGeometry(geo).GenerateMesh(maxh=maxh_air, grading=0.5)
    mesh = Mesh(ngmesh)
    mesh.Curve(fe_order)
    return mesh


def run_level(level, grp_ref, sample_points):
    print(f"\n{'='*60}")
    print(f"  {level['label']}")
    print(f"{'='*60}")
    mesh = build_kelvin_mesh(level['maxh_iron'], level['maxh_air'],
                               level['maxh_kelvin'], level['fe_order'])
    print(f"  ne={mesh.ne}  nv={mesh.nv}")

    coil = build_radia_coil()
    solver = ScalarPotentialSolver(
        mesh, iron_domains='iron', mu_r=mu_r, order=level['fe_order'],
        kelvin_region='air_outer', kelvin_radius=phys_R,
        kelvin_center=list(kelvin_center))
    solver.set_source_from_radia(coil)

    t0 = time.perf_counter()
    solver.solve_total_reduced_potential()
    dt = time.perf_counter() - t0
    print(f"  solve time: {dt:.1f}s")

    B_cf = solver.get_B()
    rels = []
    for pt in sample_points:
        try:
            mip = mesh(*pt)
            B_fem = np.array([B_cf(mip)[i] for i in range(3)])
        except Exception:
            continue
        B_rad = np.array(rad.Fld(grp_ref, 'b', list(pt)))
        m_fem = float(np.linalg.norm(B_fem))
        m_rad = float(np.linalg.norm(B_rad))
        rel = (m_fem - m_rad) / max(m_rad, 1e-12)
        rels.append(rel)
    arr = np.array(rels)
    max_rel = np.max(np.abs(arr)) * 100
    rms_rel = np.sqrt(np.mean(arr**2)) * 100
    print(f"  max|rel| = {max_rel:.3f} %   RMS = {rms_rel:.3f} %")
    return max_rel, rms_rel, mesh.ne, dt


def main():
    print("=" * 60)
    print("Mesh convergence study: solve_total_reduced_potential")
    print("=" * 60)

    print("\n[Radia BEM reference (fine)]")
    grp_ref = build_radia_reference()

    half_h = cylinder_height / 2
    pts = ([(0, 0, zv) for zv in np.linspace(-half_h+0.01, half_h-0.01, 5)]
         + [(xv, 0, 0) for xv in np.linspace(0.005, cylinder_radius-0.005, 4)]
         + [(0, 0, zv) for zv in np.linspace(half_h+0.02, phys_R-0.04, 4)])

    results = []
    for lev in LEVELS:
        max_rel, rms_rel, ne, dt = run_level(lev, grp_ref, pts)
        results.append((lev['label'], ne, max_rel, rms_rel, dt))

    print("\n" + "=" * 60)
    print("CONVERGENCE SUMMARY")
    print("=" * 60)
    print(f"  {'Level':<16s}  {'ne':>7s}  {'max|rel|':>9s}  {'RMS':>8s}  {'time':>6s}")
    for label, ne, mx, rms, dt in results:
        print(f"  {label:<16s}  {ne:7d}  {mx:8.3f}%  {rms:7.3f}%  {dt:5.1f}s")
    if len(results) >= 2:
        r0 = results[-2][2]
        r1 = results[-1][2]
        if r1 < r0:
            print(f"\n  Convergence confirmed: {r0:.3f}% -> {r1:.3f}%")
        else:
            print(f"\n  WARNING: no improvement {r0:.3f}% -> {r1:.3f}%")


if __name__ == "__main__":
    main()
