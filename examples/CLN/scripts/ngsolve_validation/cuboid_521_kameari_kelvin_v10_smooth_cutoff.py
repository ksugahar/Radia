"""Kameari + Kelvin v10: smooth cutoff for A_s to eliminate gauge divergence.

Root cause from v3-v9: A_ext = (B_0/2)(-y, x, 0) is unbounded at infinity,
making Kelvin pullback A_kelvin ~ 1/rho'^3 at GND (divergent).

v10 fix: replace A_ext with A_s = A_ext * cutoff(|r|/R_K) where cutoff is
a smooth bump function that's 1 inside conductor + boundary, → 0 at and
beyond the inner Kelvin sphere R_K.

Inside conductor (|r| << R_K): A_s ≈ A_ext exactly → uniform B_z applied
Outside R_K: A_s = 0 → no contribution to Kelvin pullback (bounded)
Transition zone: smooth interpolation, no singularity

The conductor's induced current J = curl T responds to LOCAL B (= B_0 z_hat
in conductor), and the iteration extracts the Cauer ladder from this
local response. The "wrong" B at infinity (= 0 instead of B_0) shouldn't
matter for the conductor's CLN — radiation BC is preserved for A_ind
(which decays as 1/r^3 dipole, Kelvin handles correctly).

Reference: 谷本 修論 CLN_AT formulation pattern + Sugahara Kelvin pullback.
"""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace", line_buffering=True)
sys.path.insert(0, r"W:/00_CAE/Radia/01_GitHub/src/radia")
print("Starting v10 (smooth cutoff A_s)...", flush=True)

from netgen.occ import Box, Sphere, Pnt, OCCGeometry
from ngsolve import (
    Mesh, HCurl, Periodic, BilinearForm, LinearForm, GridFunction,
    CoefficientFunction, curl, dx, x, y, z, sqrt, IfPos,
    Integrate, TaskManager, ngsglobals, exp,
)
from kelvin_geometry import add_kelvin_exterior_domain
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
H_COND = 0.5e-3
H_AIR = 2.0e-3
ORDER = 2
N_STAGES = 8
BONUS_INT = 8
GAUGE_EPS = 1e-8

# Cutoff: 1 inside r < r_inner, smooth transition r_inner < r < r_outer, 0 beyond
R_CUTOFF_INNER = 0.5 * R_K   # 12.5 mm — well past cuboid (5mm radius) but inside R_K
R_CUTOFF_OUTER = 0.9 * R_K   # 22.5 mm — well before R_K boundary

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
        [inner_air, cuboid], offset=OFFSET, R_K=R_K, inner_maxh=H_AIR)
    return OCCGeometry(geo)


def main():
    ngsglobals.msg_level = 0
    print(f"=== Kameari + Kelvin v10: smooth cutoff A_s ===", flush=True)
    print(f"  R_CUTOFF_INNER = {R_CUTOFF_INNER*1000} mm (=0.5 R_K)", flush=True)
    print(f"  R_CUTOFF_OUTER = {R_CUTOFF_OUTER*1000} mm (=0.9 R_K)", flush=True)
    print(f"  Inside cuboid (max ~3 mm radius): A_s = A_ext exactly", flush=True)
    print(f"  Inside R_CUTOFF_INNER: A_s = A_ext (cutoff = 1)", flush=True)
    print(f"  Beyond R_CUTOFF_OUTER: A_s = 0 (cutoff = 0)\n", flush=True)

    geo = build_geo()
    print("Generating mesh...", flush=True)
    t0 = time.time()
    with TaskManager():
        ngmesh = geo.GenerateMesh(maxh=H_AIR, grading=0.5)
    mesh = Mesh(ngmesh)
    mesh.Curve(ORDER + 1)
    print(f"  ne = {mesh.ne} ({time.time()-t0:.1f}s)\n", flush=True)

    sigma_cf = mesh.MaterialCF({"conductor": sigma_Cu, "air": 0.0, "kelvin": 0.0})
    sigma_inv_cf = mesh.MaterialCF({"conductor": 1.0/sigma_Cu, "air": 0.0,
                                     "kelvin": 0.0})

    # Smooth cutoff: cubic Hermite interpolation
    # cutoff(r) = 1                                  for r < R_inner
    # cutoff(r) = 1 - 3t² + 2t³, t = (r-R_in)/(R_out-R_in)  for R_in < r < R_out
    # cutoff(r) = 0                                  for r > R_outer
    r_phys = sqrt(x**2 + y**2 + z**2 + 1e-30)
    t = (r_phys - R_CUTOFF_INNER) / (R_CUTOFF_OUTER - R_CUTOFF_INNER)
    # Clamp t to [0, 1]
    t_clamped = IfPos(t, IfPos(1.0 - t, t, 1.0), 0.0)
    cutoff = 1.0 - t_clamped*t_clamped*(3.0 - 2.0*t_clamped)
    # Outside R_outer, force 0 (cutoff above already gives 1.0 - 1*(3-2)*1 = 0; good)
    # Inside R_inner, t < 0, t_clamped = 0, cutoff = 1.0 ✓

    # A_s with smooth cutoff
    A_s_localized = CoefficientFunction((-y, x, 0)) * 0.5 * cutoff

    # Initial impressed J source for Kameari
    # Inside conductor (well within R_CUTOFF_INNER): cutoff = 1, A_s = A_ext
    # So J_0 = sigma * A_s_localized = sigma * A_ext in conductor
    J_cf = sigma_cf * A_s_localized

    # FE space: full HCurl on Kelvin geometry, GND vertex Dirichlet
    fes = Periodic(HCurl(mesh, order=ORDER, dirichlet_bbnd="GND",
                         nograds=True))
    print(f"  HCurl ndof = {fes.ndof}", flush=True)
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

    print("  Assembling+factor (Kelvin)...", flush=True)
    t0 = time.time()
    with TaskManager():
        a.Assemble()
        try:
            inv = a.mat.Inverse(fes.FreeDofs(), inverse="pardiso")
        except Exception:
            inv = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")
    print(f"    {time.time()-t0:.1f}s\n", flush=True)

    # Verify R_0 (should match Parseval since cutoff = 1 in conductor)
    R_inv_check = float(Integrate(J_cf * J_cf * sigma_inv_cf
                                  * dx("conductor", bonus_intorder=BONUS_INT), mesh))
    R0_test = 1.0 / R_inv_check if R_inv_check > 1e-30 else float('inf')
    print(f"  R_0 verification: {R0_test:.6e} vs Parseval analytical "
          f"{R0_anal:.6e}  ratio = {R0_test/R0_anal:.6f}\n", flush=True)

    raise RuntimeError("The historical CLN field-reduction implementation was retired.")


if __name__ == "__main__":
    main()
