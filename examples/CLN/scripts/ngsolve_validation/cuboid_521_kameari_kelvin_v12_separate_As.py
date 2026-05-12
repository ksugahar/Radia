"""Kameari + Kelvin v12: A_s 別カウント版 (gfA_s 投影 + 内蔵 curl).

v11 反省:
  - A_s_phys = (B_0/2)(-y, x, 0) は無限遠で線形発散。
  - Kelvin pullback すると offset 点で 1/rho'^3 発散、curl は 1/rho'^4 発散。
  - (nu - nu_0) curl(A_s) の Kelvin 領域積分が破綻 (tau_0 = 1.95e9 us, L_1 < 0)。

v12 方針 (ご指導: 「A_s は別カウント」):
  - A_s は inner (conductor + air) のみで定義、Kelvin 領域では A_s = 0。
  - HCurl GridFunction gfA_s に Set() で投影 (内部領域は線形なので機械精度再現)。
  - 弱形式の curl 項は内蔵 curl(gfA_s) を使用 (CF の analytical curl を回避)。
  - Kelvin 領域は純粋に A_r のみが摂動を伝搬。
  - 2 段試験 (N_STAGES = 2) で formulation 健全性確認。

Reduced-A 弱形式 (A_s が inner 限定の場合):
    a(A_r, v) = (J_imp, v)_cond
              - integral over kelvin of (nu_kelvin - nu_0) curl(gfA_s) . curl(v) dV

  Kelvin 領域では gfA_s が界面近傍のみで非零 (HCurl tangential continuity)、
  curl(gfA_s) も bounded、特異性なし。

ご指導の「2 段試験」: L_0, L_1 共に正かつ tau 妥当なら formulation OK。
"""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace", line_buffering=True)
sys.path.insert(0, r"W:/00_CAE/Radia/01_GitHub/src/radia")
print("Starting v12 (A_s separated, gfA_s projection + builtin curl)...",
      flush=True)

from netgen.occ import Box, Sphere, Pnt, OCCGeometry
from ngsolve import (
    Mesh, HCurl, Periodic, BilinearForm, LinearForm, GridFunction,
    CoefficientFunction, curl, dx, x, y, z, sqrt,
    Integrate, TaskManager, ngsglobals,
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
N_STAGES = 2
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
    print(f"=== Kameari + Kelvin v12: A_s SEPARATED, 2-stage test ===",
          flush=True)
    print(f"  A_s defined only in inner (conductor + air), zero in Kelvin",
          flush=True)
    print(f"  curl(gfA_s) used via NGSolve builtin (no analytical CF curl)",
          flush=True)
    print(f"  Reference Parseval R_0 = {R0_anal:.6e}", flush=True)
    print(f"  Target tau_0 (mpmath) ~ 14.4 us\n", flush=True)

    geo = build_geo()
    print("Generating mesh...", flush=True)
    t0 = time.time()
    with TaskManager():
        ngmesh = geo.GenerateMesh(maxh=H_AIR, grading=0.5)
    mesh = Mesh(ngmesh)
    mesh.Curve(ORDER + 1)
    print(f"  ne = {mesh.ne}  ({time.time()-t0:.1f}s)", flush=True)
    print(f"  materials = {mesh.GetMaterials()}\n", flush=True)

    sigma_cf = mesh.MaterialCF({"conductor": sigma_Cu, "air": 0.0,
                                 "kelvin": 0.0})
    sigma_inv_cf = mesh.MaterialCF({"conductor": 1.0/sigma_Cu, "air": 0.0,
                                     "kelvin": 0.0})

    nu_cf = make_kelvin_nu_cf(mesh, R_K, OFFSET, nu_0=NU_0,
                              kelvin_mats=("kelvin",))

    # === A_s analytical CF: B_0 = 1 T uniform z, defined ONLY in inner ===
    A_s_phys_cf = CoefficientFunction((-y, x, 0)) * 0.5

    # === FE space ===
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

    u, v = fes.TnT()

    # === gfA_s: project A_s_phys_cf onto HCurl, only on inner domain ===
    # Set() with definedon restricts the canonical Nedelec interpolation to
    # those elements. Interface edges/faces shared with Kelvin elements get
    # values from the inner side; pure-Kelvin dofs stay zero.
    gfA_s = GridFunction(fes, name="A_s")
    gfA_s.vec[:] = 0.0
    inner_region = mesh.Materials("conductor|air")
    print("  Projecting A_s onto inner region...", flush=True)
    t0 = time.time()
    with TaskManager():
        gfA_s.Set(A_s_phys_cf, definedon=inner_region,
                  bonus_intorder=BONUS_INT)
    print(f"    {time.time()-t0:.1f}s", flush=True)

    # Diagnostic: in conductor, gfA_s should match A_s_phys exactly (linear)
    err_in_cond = float(Integrate(
        (gfA_s - A_s_phys_cf) * (gfA_s - A_s_phys_cf)
        * dx("conductor", bonus_intorder=BONUS_INT), mesh))
    norm_in_cond = float(Integrate(
        A_s_phys_cf * A_s_phys_cf
        * dx("conductor", bonus_intorder=BONUS_INT), mesh))
    rel_err_cond = (err_in_cond / norm_in_cond) ** 0.5 if norm_in_cond > 0 else 0
    print(f"  ||gfA_s - A_s_phys||/||A_s_phys|| in conductor = "
          f"{rel_err_cond:.2e}", flush=True)

    # Diagnostic: curl(gfA_s) should equal (0, 0, 1) in conductor
    curl_gfAs = curl(gfA_s)
    curl_target = CoefficientFunction((0, 0, 1))
    err_curl = float(Integrate(
        (curl_gfAs - curl_target) * (curl_gfAs - curl_target)
        * dx("conductor", bonus_intorder=BONUS_INT), mesh))
    rel_err_curl = (err_curl / V_cond) ** 0.5
    print(f"  ||curl(gfA_s) - z_hat|| / sqrt(V_cond) in conductor = "
          f"{rel_err_curl:.2e}\n", flush=True)

    # === Bilinear form ===
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

    # === Verify R_0 (impressed J^2/sigma in conductor) ===
    R_inv_check = float(Integrate(sigma_cf * gfA_s * gfA_s
                                  * dx("conductor", bonus_intorder=BONUS_INT),
                                  mesh))
    R0_test = 1.0 / R_inv_check if R_inv_check > 1e-30 else float('inf')
    print(f"  R_0 verification: {R0_test:.6e} vs Parseval {R0_anal:.6e}, "
          f"ratio {R0_test/R0_anal:.6f}\n", flush=True)

    # === Kameari iteration: reduced-A with A_s separated ===
    diag = []
    # Stage 0: J_imp = sigma * gfA_s in conductor (gfA_s = A_s_phys here)
    J_imp_cf = sigma_cf * gfA_s
    # gfApot accumulates A_r contributions across stages
    gfApot = GridFunction(fes)
    gfApot.vec[:] = 0.0

    raise RuntimeError("The historical CLN field-reduction implementation was retired.")


if __name__ == "__main__":
    main()
