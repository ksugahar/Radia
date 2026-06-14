"""Diagnose why my order=3 Kameari shows drift jumping to 1e-5 at stage 2.

User's hypothesis: nograds=True alone should make high order work.
My v1: nograds=True + tree-cotree mask (had drift jump)

Test 4 configurations to isolate cause:
  (A) order=3, nograds=True, tree-cotree mask  -- my v1 setup (drift jump)
  (B) order=3, nograds=True, NO tree-cotree    -- pure nograds
  (C) order=3, nograds=False, tree-cotree mask -- pure tree-cotree (gradient kernel present)
  (D) order=3, bonus_intorder=8                -- higher integration order

If (B) gives clean machine-precision drift, the problem was tree-cotree
clashing with nograds at order >= 3.
If (B) still has drift, then need bonus_intorder or other fix.
"""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")

from netgen.occ import Box, Pnt, OCCGeometry
from ngsolve import (
    Mesh, HCurl, BilinearForm, LinearForm, GridFunction,
    CoefficientFunction, curl, dx, x, y, z,
    Integrate, TaskManager, ngsglobals,
)
from collections import deque
from math import pi
import numpy as np
import json
from pathlib import Path

mu0 = 4 * pi * 1e-7
sigma_Cu = 5.8e7
ax, ay, az = 5e-3, 2e-3, 1e-3
H_COND = 0.30e-3
N_STAGES = 12

tau_TE_110_us = mu0 * sigma_Cu / (pi**2 * (1/ax**2 + 1/ay**2)) * 1e6


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


def build_mesh():
    cuboid = Box(Pnt(-ax/2, -ay/2, -az/2), Pnt(ax/2, ay/2, az/2))
    cuboid.mat("conductor").bc("conductor_surface")
    cuboid.maxh = H_COND
    return Mesh(OCCGeometry(cuboid).GenerateMesh(maxh=H_COND))


def kameari_run(mesh, order, use_nograds, use_treecotree, bonus_intorder=0):
    fes = HCurl(mesh, order=order, dirichlet="conductor_surface",
                nograds=use_nograds)
    print(f"    HCurl ndof = {fes.ndof}, nograds={use_nograds}")
    if use_treecotree:
        tree_edges = build_spanning_tree(mesh)
        fd = fes.FreeDofs()
        masked = 0
        for edge_nr in tree_edges:
            edge = mesh.edges[edge_nr]
            dofs = fes.GetDofNrs(edge)
            if dofs and fd[dofs[0]]:
                fd[dofs[0]] = False
                masked += 1
        print(f"    tree-cotree masked {masked} edges")

    sigma_inv_cf = 1.0 / sigma_Cu
    A_ext = CoefficientFunction((-y, x, 0)) * 0.5

    u, v = fes.TnT()
    a = BilinearForm(fes)
    if bonus_intorder > 0:
        a += (1.0/mu0) * curl(u) * curl(v) * dx(bonus_intorder=bonus_intorder)
    else:
        a += (1.0/mu0) * curl(u) * curl(v) * dx
    with TaskManager():
        a.Assemble()
        try:
            inv = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")
        except Exception as e:
            print(f"    factor failed: {e}")
            return None

    fes_proj = HCurl(mesh, order=order, dirichlet="conductor_surface",
                     nograds=use_nograds)

    diag = []
    J_history = []
    J_cf = sigma_Cu * A_ext
    Apot = None

    for n in range(N_STAGES):
        if bonus_intorder > 0:
            f = LinearForm(fes)
            f += J_cf * v * dx(bonus_intorder=bonus_intorder)
        else:
            f = LinearForm(fes)
            f += J_cf * v * dx
        with TaskManager():
            f.Assemble()
        gfA = GridFunction(fes)
        with TaskManager():
            gfA.vec.data = inv * f.vec

        if bonus_intorder > 0:
            R_inv = float(Integrate(J_cf * J_cf * sigma_inv_cf
                                    * dx(bonus_intorder=bonus_intorder), mesh))
        else:
            R_inv = float(Integrate(J_cf * J_cf * sigma_inv_cf * dx, mesh))
        if R_inv < 1e-30 or R_inv < 0:
            break
        Rn = 1.0 / R_inv

        gfJ_n = GridFunction(fes_proj)
        gfJ_n.Set(J_cf)

        # Schmidt drift
        cross = []
        for gfJ_m in J_history:
            if bonus_intorder > 0:
                xnm = float(Integrate(J_cf * gfJ_m * sigma_inv_cf
                                      * dx(bonus_intorder=bonus_intorder), mesh))
            else:
                xnm = float(Integrate(J_cf * gfJ_m * sigma_inv_cf * dx, mesh))
            cross.append(xnm)
        max_cross = max((abs(c) for c in cross), default=0.0)
        rel_drift = max_cross / R_inv

        if Apot is None:
            Apot = Rn * gfA
        else:
            Apot = Apot + Rn * gfA

        if bonus_intorder > 0:
            Ln_int = float(Integrate(J_cf * Apot * dx(bonus_intorder=bonus_intorder), mesh))
        else:
            Ln_int = float(Integrate(J_cf * Apot * dx, mesh))
        Ln = Rn * Ln_int

        diag.append({"n": n, "R_n": Rn, "L_n": Ln,
                     "tau_us": Ln/Rn*1e6 if Rn > 0 and Ln > 0 else None,
                     "rel_drift": rel_drift, "energy_norm": R_inv})

        J_history.append(gfJ_n)

        if Ln <= 0:
            break
        J_cf = J_cf - sigma_Cu * Apot / Ln

    return diag


def main():
    ngsglobals.msg_level = 0
    print(f"Mesh: 5x2x1 mm closed PEC cuboid, h={H_COND*1000} mm")
    mesh = build_mesh()
    print(f"  ne = {mesh.ne}, nv = {mesh.nv}")
    print(f"Reference tau = {tau_TE_110_us:.3f} us\n")

    configs = [
        ("(A) order=3, nograds=T, tree-cotree", 3, True, True, 0),
        ("(B) order=3, nograds=T, NO tree-cotree", 3, True, False, 0),
        ("(C) order=3, nograds=F, tree-cotree only", 3, False, True, 0),
        ("(D) order=3, nograds=T+treecotree, bonus=8", 3, True, True, 8),
    ]

    all_results = []
    for name, order, ngs, tc, bonus in configs:
        print("=" * 78)
        print(name)
        print("=" * 78)
        try:
            d = kameari_run(mesh, order, ngs, tc, bonus_intorder=bonus)
        except Exception as e:
            print(f"    EXCEPTION: {e}")
            continue
        if d is None:
            continue
        # Pretty print drift trajectory
        print(f"    {'n':>3} {'rel_drift':>14} {'tau_n [us]':>14} {'energy':>14}")
        for r in d:
            tau_s = f"{r['tau_us']:.3f}" if r['tau_us'] else " N/A"
            print(f"    {r['n']:>3} {r['rel_drift']:>14.4e} {tau_s:>14} "
                  f"{r['energy_norm']:>14.4e}")
        bd = next((r['n'] for r in d if r['rel_drift'] > 1e-2), None)
        max_tau = max((r['tau_us'] for r in d if r['tau_us']), default=0)
        print(f"\n    Breakdown stage (drift > 1%): {bd}")
        print(f"    max tau_n = {max_tau:.3f} us "
              f"(analytical {tau_TE_110_us:.3f}, ratio {max_tau/tau_TE_110_us:.4f})")
        print()
        all_results.append({"config": name, "stages": d, "breakdown": bd,
                            "max_tau_us": max_tau})

    print("=" * 78)
    print("SUMMARY")
    print("=" * 78)
    print(f"{'config':>50} {'breakdown':>12} {'max tau':>12}")
    for r in all_results:
        bd = f"n={r['breakdown']}" if r['breakdown'] is not None else "(none in 12)"
        print(f"{r['config']:>50} {bd:>12} {r['max_tau_us']:>12.3f}")


if __name__ == "__main__":
    main()
