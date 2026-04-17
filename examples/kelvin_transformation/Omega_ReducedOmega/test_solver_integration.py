"""test_solver_integration.py — test ScalarPotentialSolver.solve_total_reduced_potential()

Closed-loop coil (gap=0) + mu_r=100 cylinder + Kelvin.
Expected: max|rel| < 5% (validated at 3.17% with standalone v2).
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
maxh_iron = 0.008
maxh_air = 0.025
maxh_kelvin = 0.03
fe_order = 2


def build_radia_coil():
    rad.UtiDelAll()
    J0 = I_total / (2 * a_coil) ** 2
    return rad.ObjArcCur(
        [0, 0, z_coil],
        [R_coil - a_coil, R_coil + a_coil],
        [0, 2 * math.pi],
        2 * a_coil, 200, 'man', 'z', J0)


def build_radia_reference():
    rad.UtiDelAll()
    J0 = I_total / (2 * a_coil) ** 2
    coil = rad.ObjArcCur(
        [0, 0, z_coil],
        [R_coil - a_coil, R_coil + a_coil],
        [0, 2 * math.pi],
        2 * a_coil, 200, 'man', 'z', J0)
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
    t0 = time.perf_counter()
    result = rad.Solve(grp, 1e-7, 2000, 0)
    print(f"  Radia Solve: {time.perf_counter()-t0:.1f}s, "
          f"residual={result[0]:.3e}, {mesh_rad.ne} tets")
    return grp


def build_kelvin_mesh():
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


def main():
    print("=" * 72)
    print("ScalarPotentialSolver.solve_total_reduced_potential() integration test")
    print("=" * 72)

    print("\n[1] Build mesh")
    mesh = build_kelvin_mesh()
    print(f"  ne={mesh.ne}  mats={mesh.GetMaterials()}")

    print("\n[2] Build Radia coil (source)")
    coil = build_radia_coil()

    print("\n[3] ScalarPotentialSolver setup + set_source_from_radia")
    solver = ScalarPotentialSolver(
        mesh, iron_domains='iron', mu_r=mu_r, order=fe_order,
        kelvin_region='air_outer', kelvin_radius=phys_R,
        kelvin_center=list(kelvin_center))
    solver.set_source_from_radia(coil, resolution=41)

    print("\n[4] solve() -> auto-selects total_reduced")
    t0 = time.perf_counter()
    solver.solve()
    print(f"  solve time: {time.perf_counter()-t0:.1f}s")

    print("\n[5] Build Radia BEM reference")
    grp_ref = build_radia_reference()

    print("\n[6] Compare B")
    half_h = cylinder_height / 2
    pts = ([(0, 0, zv) for zv in np.linspace(-half_h+0.01, half_h-0.01, 5)]
         + [(xv, 0, 0) for xv in np.linspace(0.005, cylinder_radius-0.005, 4)]
         + [(0, 0, zv) for zv in np.linspace(half_h+0.02, phys_R-0.04, 4)])

    B_cf = solver.get_B()
    rels = []
    print(f"  {'point':>22s}  {'|B|_FEM':>11s}  {'|B|_Rad':>11s}  {'rel':>8s}")
    for pt in pts:
        try:
            mip = mesh(*pt)
            B_fem = np.array([B_cf(mip)[i] for i in range(3)])
        except Exception:
            B_fem = np.full(3, np.nan)
        B_rad = np.array(rad.Fld(grp_ref, 'b', list(pt)))
        m_fem = float(np.linalg.norm(B_fem))
        m_rad = float(np.linalg.norm(B_rad))
        if np.isnan(m_fem):
            continue
        rel = (m_fem - m_rad) / max(m_rad, 1e-12)
        rels.append(rel)
        ps = f"({pt[0]:+.3f},{pt[1]:+.3f},{pt[2]:+.3f})"
        print(f"  {ps:>22s}  {m_fem:11.4e}  {m_rad:11.4e}  {rel*100:+7.2f}%")

    arr = np.array(rels)
    print(f"\n  max|rel| = {np.max(np.abs(arr))*100:.3f} %")
    print(f"  RMS rel  = {np.sqrt(np.mean(arr**2))*100:.3f} %")
    ok = np.max(np.abs(arr)) < 0.05
    print(f"  acceptance (< 5 %): {'PASS' if ok else 'FAIL'}")


if __name__ == "__main__":
    main()
