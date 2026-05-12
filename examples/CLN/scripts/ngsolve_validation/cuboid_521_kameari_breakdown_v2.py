"""RETIRED field-recurrence experiment. Only independent setup remains executable.
The following is historical context, not current functionality:

Robust sudden-breakdown hunt: project J_cf and Apot to GridFunction
at every stage so the CF tree doesn't grow unboundedly.

Periodic JSON snapshot every 5 stages so partial results are saved
even if a later stage crashes.
"""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace", line_buffering=True)
print("Starting...", flush=True)

from netgen.occ import Box, Pnt, OCCGeometry
from ngsolve import (
    Mesh, HCurl, BilinearForm, LinearForm, GridFunction,
    CoefficientFunction, curl, dx, x, y, z,
    Integrate, TaskManager, ngsglobals,
)
from collections import deque
from math import pi
import json, time
from pathlib import Path

mu0 = 4 * pi * 1e-7
sigma_Cu = 5.8e7
ax, ay, az = 5e-3, 2e-3, 1e-3
H_COND = 0.30e-3
N_STAGES = 30
ORDER = 3
BONUS_INT = 8
SNAPSHOT_EVERY = 3   # JSON dump every 3 stages

tau_TE_110_us = mu0 * sigma_Cu / (pi**2 * (1/ax**2 + 1/ay**2)) * 1e6

OUT_PATH = Path(__file__).parent / "cuboid_521_kameari_breakdown_v2.json"


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


def main():
    ngsglobals.msg_level = 0
    print(f"Reference tau_TE_110 = {tau_TE_110_us:.3f} us", flush=True)
    print(f"Setup: order={ORDER}, nograds=True, tree-cotree, bonus={BONUS_INT}\n",
          flush=True)

    cuboid = Box(Pnt(-ax/2, -ay/2, -az/2), Pnt(ax/2, ay/2, az/2))
    cuboid.mat("conductor").bc("conductor_surface")
    cuboid.maxh = H_COND
    print("Generating mesh...", flush=True)
    mesh = Mesh(OCCGeometry(cuboid).GenerateMesh(maxh=H_COND))
    print(f"  ne = {mesh.ne}, nv = {mesh.nv}", flush=True)

    fes = HCurl(mesh, order=ORDER, dirichlet="conductor_surface", nograds=True)
    print(f"  HCurl order={ORDER} ndof = {fes.ndof}", flush=True)

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

    sigma_inv_cf = 1.0 / sigma_Cu
    A_ext = CoefficientFunction((-y, x, 0)) * 0.5

    u, v = fes.TnT()
    a = BilinearForm(fes)
    a += (1.0/mu0) * curl(u) * curl(v) * dx(bonus_intorder=BONUS_INT)
    print("  Assembling+factorizing...", flush=True)
    t0 = time.time()
    with TaskManager():
        a.Assemble()
        inv = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")
    print(f"    {time.time()-t0:.1f}s\n", flush=True)

    diag = []
    J_history = []   # GridFunctions for cross-products

    # Initialize J as a GridFunction (Stage 0)
    gfJ = GridFunction(fes)
    gfJ.Set(sigma_Cu * A_ext)

    # Apot as accumulator GridFunction (initially zero)
    gfApot = GridFunction(fes)
    gfApot.vec[:] = 0.0

    breakdown_stage = None

    # Historical field-recurrence driver retired; assembly above is retained.
    return


if __name__ == "__main__":
    main()
