"""Kameari + Kelvin v17: A-method gauge fix via Helmholtz-Hodge (2026-05-05).

Sugahara guidance (this session):
  - v15 sign-flips at stage 1 = FEM is NOT solving correctly (not a
    natural Kameari breakdown). A-method needs ∇·A = 0 gauge enforced.
  - Canonical recipe (memory project_3D_CLN_kameari_in_vacuum_hardness):
    nograds=True / type1=True, NO penalty, Helmholtz-Hodge OR tree-cotree.
  - For Kelvin BC (Periodic + GND vertex), tree-cotree alone is
    insufficient — Periodic identification creates equivalence classes
    that BFS spanning tree doesn't fully cover. Need explicit Helmholtz-
    Hodge projection.

v17 changes from v15:
  - REMOVE GAUGE_EPS penalty (canonical: no penalty)
  - ADD Helmholtz-Hodge projection after each FEM solve:
      A_div_free = A - grad(phi)
      where phi solves: ∫ ∇phi·∇psi dx = ∫ A·∇psi dx in matching H1 space
  - Keep nograds=True + tree-cotree gauge (canonical)
  - Keep Pardiso (next iteration: try shifted AMS preconditioner if v17
    insufficient, per Sugahara hint)

Expected if v17 works: clean stages 0..N (drift at machine precision,
all L_n positive), then natural Kameari breakdown onset at higher N
(this would be the paper motivation data — natural breakdown after a
healthy regime).
"""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace", line_buffering=True)
sys.path.insert(0, r"S:/Radia/01_GitHub/src/radia")
print("Starting v17 (no penalty + Helmholtz-Hodge projection)...", flush=True)

from netgen.occ import Box, Sphere, Pnt, OCCGeometry
from ngsolve import (
    Mesh, HCurl, H1, Periodic, BilinearForm, LinearForm, GridFunction,
    CoefficientFunction, curl, grad, dx, x, y, z,
    Integrate, TaskManager, ngsglobals,
)
from kelvin_geometry import add_kelvin_exterior_domain
from kelvin_material import (
    make_kelvin_nu_cf,
    make_reduced_potential_background_cf,
    NU_0,
)
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
# GAUGE_EPS removed — canonical recipe uses no penalty.
# H1 Poisson for Helmholtz-Hodge needs a tiny shift to be invertible
# (the constant function is in its kernel without Dirichlet on a face).
PHI_SHIFT = 1e-12

R0_anal = 48 / (sigma_Cu * V_cond * (ax**2 + ay**2))
ELF_TAU_LEAD = 11.51e-6


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
    print(f"=== v15: canonical CLN iteration with Kelvin ===", flush=True)
    print(f"  R_K = {R_K*1000:.1f} mm, OFFSET = {OFFSET}", flush=True)
    print(f"  Reference R_0 (Parseval) = {R0_anal:.6e}", flush=True)
    print(f"  ELF reference tau_lead = {ELF_TAU_LEAD*1e6:.2f} us", flush=True)
    print(f"  N_STAGES = {N_STAGES}\n", flush=True)

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

    A_s_cf = make_reduced_potential_background_cf(
        mesh,
        F_inner_factory=lambda xc, yc, zc: CoefficientFunction((-yc, xc, 0)) * 0.5,
        R_K=R_K, offset=OFFSET, kelvin_mats=("kelvin",), dim=3,
    )

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

    gfA_s = GridFunction(fes, name="A_s")
    gfA_s.vec[:] = 0.0
    print("  Projecting A_s onto HCurl (Periodic)...", flush=True)
    t0 = time.time()
    with TaskManager():
        gfA_s.Set(A_s_cf, bonus_intorder=BONUS_INT)
    print(f"    {time.time()-t0:.1f}s", flush=True)

    A_s_phys_cf = CoefficientFunction((-y, x, 0)) * 0.5
    err_in_cond = float(Integrate(
        (gfA_s - A_s_phys_cf) * (gfA_s - A_s_phys_cf)
        * dx("conductor", bonus_intorder=BONUS_INT), mesh))
    norm_in_cond = float(Integrate(
        A_s_phys_cf * A_s_phys_cf
        * dx("conductor", bonus_intorder=BONUS_INT), mesh))
    rel_err = (err_in_cond / norm_in_cond) ** 0.5 if norm_in_cond > 0 else 0
    print(f"  ||gfA_s - A_s_phys||/||A_s_phys|| in conductor = {rel_err:.2e}",
          flush=True)

    a = BilinearForm(fes)
    a += nu_cf * curl(u) * curl(v) * dx(bonus_intorder=BONUS_INT)
    # NO penalty term — canonical recipe relies on tree-cotree gauge alone
    # for the curl-curl operator. Helmholtz-Hodge projection is applied
    # post-solve to clean residual gradient components introduced through
    # the Periodic Kelvin boundary.

    print("  Assembling+factor (curl-curl)...", flush=True)
    t0 = time.time()
    with TaskManager():
        a.Assemble()
        try:
            inv = a.mat.Inverse(fes.FreeDofs(), inverse="pardiso")
        except Exception:
            inv = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")
    print(f"    {time.time()-t0:.1f}s\n", flush=True)

    # === Helmholtz-Hodge: H1 space + Poisson factorization ===
    # Same Periodic + GND vertex Dirichlet as HCurl, so H1 gradients live
    # in the same constraint structure.
    fes_phi = Periodic(H1(mesh, order=ORDER, dirichlet_bbnd="GND"))
    print(f"  H1 ndof (for Helmholtz-Hodge) = {fes_phi.ndof}", flush=True)
    phi, psi = fes_phi.TnT()
    a_phi = BilinearForm(fes_phi)
    a_phi += grad(phi) * grad(psi) * dx(bonus_intorder=BONUS_INT)
    a_phi += PHI_SHIFT * phi * psi * dx(bonus_intorder=BONUS_INT)
    print("  Assembling+factor (H1 Poisson for Helmholtz-Hodge)...",
          flush=True)
    t0 = time.time()
    with TaskManager():
        a_phi.Assemble()
        try:
            inv_phi = a_phi.mat.Inverse(fes_phi.FreeDofs(), inverse="pardiso")
        except Exception:
            inv_phi = a_phi.mat.Inverse(fes_phi.FreeDofs(),
                                        inverse="sparsecholesky")
    print(f"    {time.time()-t0:.1f}s\n", flush=True)

    def helmholtz_hodge_project(gfA_in):
        """Subtract grad(phi) such that ∫ (gfA - grad(phi)) · grad(psi) dx = 0
        for all psi in H1 (Periodic, GND-Dirichlet)."""
        f_phi = LinearForm(fes_phi)
        f_phi += gfA_in * grad(psi) * dx(bonus_intorder=BONUS_INT)
        with TaskManager():
            f_phi.Assemble()
        gf_phi = GridFunction(fes_phi)
        with TaskManager():
            gf_phi.vec.data = inv_phi * f_phi.vec
        gfA_proj = GridFunction(fes)
        with TaskManager():
            gfA_proj.Set(gfA_in - grad(gf_phi), bonus_intorder=BONUS_INT)
        return gfA_proj

    R_inv_check = float(Integrate(sigma_cf * gfA_s * gfA_s
                                  * dx("conductor", bonus_intorder=BONUS_INT),
                                  mesh))
    R0_test = 1.0 / R_inv_check if R_inv_check > 1e-30 else float('inf')
    print(f"  R_0 verification: {R0_test:.6e} vs Parseval {R0_anal:.6e}, "
          f"ratio {R0_test/R0_anal:.6f}\n", flush=True)

    # Breakdown demo: keep iterating past sign flip to record the full
    # Schmidt drift, energy norm growth, and L_n sign trajectory.
    diag = []
    J_history = []  # store gfJ_n for cross-product Schmidt drift diagnostic
    J_imp_cf = sigma_cf * gfA_s
    gfApot = GridFunction(fes)
    gfApot.vec[:] = 0.0
    sign_flip_first_n = None

    raise RuntimeError("The historical CLN field-reduction implementation was retired.")


if __name__ == "__main__":
    main()
