"""Kameari + Kelvin v8: GND as small sphere with surface Dirichlet (COMSOL style).

Insight from W:/30_CauerLadderNetwork/2020_11_04_線形のCLNの練習/CLN_H1_mode_Kelvin_NG.m:
   The COMSOL practice uses a small CIRCLE around the image-of-infinity, not
   a single vertex. The Dirichlet condition A=0 is applied on the SURFACE of
   this small circle, not at a single point.

In 3D, the analogue is a small SPHERE around the GND point in the Kelvin
material. Setting tangential A=0 on its surface gives a much stronger
constraint than a single-vertex Dirichlet.

v8 implementation:
  - Outer Kelvin sphere has a small inner sphere (radius ~ 0.02 * R_K) removed
  - On the small sphere's surface: Dirichlet A x n = 0
  - Periodic identification at outer Kelvin sphere face (as before)

Hypothesis: This stronger gauge fix at "infinity" makes the Kelvin pullback
properly enforce A → 0 at infinity, fixing the Kameari iteration's τ_0.
"""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace", line_buffering=True)
sys.path.insert(0, r"W:/00_CAE/Radia/01_GitHub/src/radia")
print("Starting v8...", flush=True)

from netgen.occ import Box, Sphere, Pnt, OCCGeometry, Glue, IdentificationType
from ngsolve import (
    Mesh, HCurl, Periodic, BilinearForm, LinearForm, GridFunction,
    CoefficientFunction, curl, dx, x, y, z,
    Integrate, TaskManager, ngsglobals,
)
from kelvin_material import make_kelvin_nu_cf, NU_0
from collections import deque
from math import pi
import json, time
from pathlib import Path

mu0 = 4 * pi * 1e-7
sigma_Cu = 5.8e7
ax, ay, az = 5e-3, 2e-3, 1e-3
V_cond = ax * ay * az

R_K = 25e-3
OFFSET = (2.5 * R_K, 0, 0)
R_GND_SPHERE = 0.02 * R_K   # 0.5 mm small sphere around GND (COMSOL convention)
H_COND = 0.5e-3
H_AIR = 2.0e-3
H_GND = 0.05e-3   # finer mesh near small GND sphere
ORDER = 2
N_STAGES = 6
BONUS_INT = 8
GAUGE_EPS = 1e-8

R0_anal = 48 / (sigma_Cu * V_cond * (ax**2 + ay**2))


def build_spanning_tree(mesh):
    nv = mesh.nv
    visited = [False] * nv
    tree_edges = []
    adj = [[] for _ in range(nv)]
    for ed in mesh.edges:
        v0, v1 = ed.vertices[0].nr, ed.vertices[1].nr
        adj[v0].append((v1, ed.nr))
        adj[v1].append((v0, ed.nr))
    visited[0] = True
    queue = deque([0])
    while queue:
        v = queue.popleft()
        for vn, edn in adj[v]:
            if not visited[vn]:
                visited[vn] = True
                tree_edges.append(edn)
                queue.append(vn)
    return tree_edges


def build_geo_v8():
    """Geometry with small sphere around GND for surface Dirichlet."""
    cuboid = Box(Pnt(-ax/2, -ay/2, -az/2), Pnt(ax/2, ay/2, az/2))
    cuboid.name = "conductor"
    cuboid.maxh = H_COND

    sphere_inner = Sphere(Pnt(0, 0, 0), R_K)
    for f in sphere_inner.faces:
        f.name = "kelvin_int"
    sphere_inner.maxh = H_AIR

    inner_air = sphere_inner - cuboid
    inner_air.name = "air"

    # Outer Kelvin sphere
    outer_kelvin = Sphere(Pnt(*OFFSET), R_K)
    for f in outer_kelvin.faces:
        f.name = "kelvin_ext"
    outer_kelvin.maxh = H_AIR * 2

    # Small sphere around GND (image of infinity) — to be removed from outer
    gnd_sphere = Sphere(Pnt(*OFFSET), R_GND_SPHERE)
    for f in gnd_sphere.faces:
        f.name = "GND_surface"  # this is what we Dirichlet-constrain
    gnd_sphere.maxh = H_GND

    # Kelvin material = outer_kelvin minus gnd_sphere
    kelvin_mat = outer_kelvin - gnd_sphere
    kelvin_mat.name = "kelvin"

    # Periodic identification: kelvin_int (inner) <-> kelvin_ext (outer)
    inner_face = None
    outer_face = None
    for f in inner_air.faces:
        if f.name == "kelvin_int":
            inner_face = f
            break
    for f in kelvin_mat.faces:
        if f.name == "kelvin_ext":
            outer_face = f
            break
    if inner_face is None or outer_face is None:
        raise RuntimeError(f"Failed to find faces: inner={inner_face}, outer={outer_face}")
    inner_face.Identify(outer_face, "kelvin_periodic", IdentificationType.PERIODIC)

    geometry = Glue([inner_air, cuboid, kelvin_mat])
    return OCCGeometry(geometry)


def main():
    ngsglobals.msg_level = 0
    print(f"=== Kameari + Kelvin v8: GND = small sphere with surface Dirichlet ===",
          flush=True)
    print(f"  R_K = {R_K*1000} mm, R_GND_SPHERE = {R_GND_SPHERE*1000} mm",
          flush=True)
    print(f"  H_COND = {H_COND*1000} mm, H_AIR = {H_AIR*1000} mm, "
          f"H_GND = {H_GND*1000} mm", flush=True)
    print(f"  Reference Parseval R_0 = {R0_anal:.6e}", flush=True)
    print(f"  Target tau_0 (mpmath) ~ 14.4 us\n", flush=True)

    geo = build_geo_v8()
    print("Generating mesh...", flush=True)
    t0 = time.time()
    with TaskManager():
        ngmesh = geo.GenerateMesh(maxh=H_AIR, grading=0.5)
    mesh = Mesh(ngmesh)
    mesh.Curve(ORDER + 1)
    print(f"  ne = {mesh.ne}, nv = {mesh.nv}", flush=True)
    print(f"  Materials: {mesh.GetMaterials()}", flush=True)
    print(f"  Boundaries: {set(mesh.GetBoundaries())}", flush=True)
    print(f"  ({time.time()-t0:.1f}s)\n", flush=True)

    sigma_cf = mesh.MaterialCF({"conductor": sigma_Cu, "air": 0.0, "kelvin": 0.0})
    sigma_inv_cf = mesh.MaterialCF({"conductor": 1.0/sigma_Cu, "air": 0.0,
                                     "kelvin": 0.0})
    A_ext = CoefficientFunction((-y, x, 0)) * 0.5

    # FE space: HCurl with Dirichlet on GND_surface (full surface, not vertex)
    # Periodic identification on inner-outer Kelvin faces
    bdry_set = set(mesh.GetBoundaries())
    print(f"  Available boundaries: {bdry_set}", flush=True)
    gnd_bnd = None
    for b in bdry_set:
        if b and "gnd" in b.lower():
            gnd_bnd = b
            break
    print(f"  Selected GND boundary: '{gnd_bnd}'", flush=True)

    fes = Periodic(HCurl(mesh, order=ORDER, dirichlet=gnd_bnd, nograds=True))
    print(f"  HCurl ndof = {fes.ndof}, active = {sum(fes.FreeDofs())}",
          flush=True)

    tree_edges = build_spanning_tree(mesh)
    fd = fes.FreeDofs()
    masked = 0
    for edge_nr in tree_edges:
        edge = mesh.edges[edge_nr]
        dofs = fes.GetDofNrs(edge)
        if dofs and fd[dofs[0]]:
            fd[dofs[0]] = False
            masked += 1
    print(f"  tree-cotree masked {masked} edges", flush=True)

    nu_cf = make_kelvin_nu_cf(mesh, R_K, OFFSET, nu_0=NU_0,
                              kelvin_mats=("kelvin",))

    u, v = fes.TnT()
    a = BilinearForm(fes)
    a += nu_cf * curl(u) * curl(v) * dx(bonus_intorder=BONUS_INT)
    a += GAUGE_EPS * NU_0 * u * v * dx(bonus_intorder=BONUS_INT)

    print("  Assembling+factor...", flush=True)
    t0 = time.time()
    with TaskManager():
        a.Assemble()
        try:
            inv = a.mat.Inverse(fes.FreeDofs(), inverse="pardiso")
        except Exception:
            inv = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")
    print(f"    {time.time()-t0:.1f}s\n", flush=True)

    diag = []
    J_cf = sigma_cf * A_ext
    gfApot = GridFunction(fes)
    gfApot.vec[:] = 0.0

    raise RuntimeError("The historical CLN field-reduction implementation was retired.")


if __name__ == "__main__":
    main()
