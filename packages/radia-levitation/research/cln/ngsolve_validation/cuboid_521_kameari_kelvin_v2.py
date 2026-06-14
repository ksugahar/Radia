"""RETIRED field-recurrence experiment. Only independent setup remains executable.
The following is historical context, not current functionality:

Kameari + Kelvin v2: use solve_full_A_kelvin directly (production solver)
   with mesh.Curve(2), order=2, grading=0.5 — mirroring the canonical example
   Coil_3D_A_HCurl_with_Kelvin.py.

Key fixes vs v1:
  - mesh.Curve(2) for spherical face accuracy
  - order = 2 (1 was too low for curl-curl + conformal map)
  - grading = 0.5 in mesh generation
  - call solve_full_A_kelvin (= production code path) for each Kameari stage
"""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
sys.path.insert(0, r"W:/00_CAE/Radia/01_GitHub/src/radia")

from netgen.occ import Box, Sphere, Pnt, OCCGeometry
from ngsolve import (
    Mesh, CoefficientFunction, x, y, z, dx,
    Integrate, TaskManager, ngsglobals,
)
from kelvin_geometry import add_kelvin_exterior_domain
from kelvin_solver import solve_full_A_kelvin, NU_0
from math import pi
import json
import time
from pathlib import Path

mu0 = 4 * pi * 1e-7
sigma_Cu = 5.8e7

ax, ay, az = 5e-3, 2e-3, 1e-3
V_cond = ax * ay * az

R_K = 25e-3                     # 25 mm
OFFSET = (2.5 * R_K, 0, 0)       # = (62.5 mm, 0, 0); offset/R_K = 2.5 like the example

H_COND = 0.5e-3
H_AIR = 4.0e-3
ORDER = 2
N_STAGES = 8

R0_anal = 48 / (sigma_Cu * V_cond * (ax**2 + ay**2))


def build_geo():
    cuboid = Box(Pnt(-ax/2, -ay/2, -az/2), Pnt(ax/2, ay/2, az/2))
    cuboid.name = "conductor"
    cuboid.maxh = H_COND

    sphere_inner = Sphere(Pnt(0, 0, 0), R_K)
    for f in sphere_inner.faces:
        f.name = "kelvin_int"
    sphere_inner.maxh = H_AIR

    inner_air = sphere_inner - cuboid
    inner_air.name = "air"

    geo, info = add_kelvin_exterior_domain(
        [inner_air, cuboid],
        offset=OFFSET,
        R_K=R_K,
        inner_maxh=H_AIR,
    )
    return OCCGeometry(geo)


def kameari_stage_solve(mesh, J_cf):
    """Solve curl-curl A = J on the Kelvin geometry, return GridFunction."""
    res = solve_full_A_kelvin(
        mesh,
        J_source_cf=J_cf,
        R_K=R_K,
        offset=OFFSET,
        source_material="conductor",   # cuboid material name
        order=ORDER,
    )
    return res["gfu"], res["nu_cf"]


def main():
    ngsglobals.msg_level = 1
    print(f"=== Kameari + Kelvin v2 (production solver path) ===")
    print(f"  Conductor: {ax*1000}x{ay*1000}x{az*1000} mm Cu")
    print(f"  Kelvin: R_K = {R_K*1000} mm, offset = {[o*1000 for o in OFFSET]} mm")
    print(f"  Mesh: h_cond = {H_COND*1000} mm, h_air = {H_AIR*1000} mm, order = {ORDER}")
    print(f"  N_stages = {N_STAGES}\n")

    geo = build_geo()
    print("Generating mesh (grading=0.5)...")
    t0 = time.time()
    with TaskManager():
        ngmesh = geo.GenerateMesh(maxh=H_AIR, grading=0.5)
    mesh = Mesh(ngmesh)
    mesh.Curve(ORDER + 1)
    print(f"  ne = {mesh.ne}, nv = {mesh.nv}, "
          f"mats = {mesh.GetMaterials()}  ({time.time()-t0:.1f}s)\n")

    print(f"Reference R_0 (Parseval) = {R0_anal:.6e}")
    print(f"ELF Foster B_z tau_lead   = 11.51 us\n")

    sigma_cf = mesh.MaterialCF({"conductor": sigma_Cu, "air": 0.0, "kelvin": 0.0})
    sigma_inv_cf = mesh.MaterialCF({"conductor": 1.0/sigma_Cu, "air": 0.0, "kelvin": 0.0})

    A_ext = CoefficientFunction((-y, x, 0)) * 0.5  # B_0 = 1, z-direction

    diag = []
    J_cf = sigma_cf * A_ext
    Apot = None

    # Historical field-recurrence driver retired; assembly above is retained.
    return


if __name__ == "__main__":
    main()
