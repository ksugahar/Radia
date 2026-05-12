"""Robust sudden-breakdown hunt: project J_cf and Apot to GridFunction
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

    for n in range(N_STAGES):
        t_stage = time.time()

        f = LinearForm(fes)
        f += gfJ * v * dx(bonus_intorder=BONUS_INT)
        with TaskManager():
            f.Assemble()
        gfA = GridFunction(fes)
        with TaskManager():
            gfA.vec.data = inv * f.vec

        R_inv = float(Integrate(gfJ * gfJ * sigma_inv_cf
                                * dx(bonus_intorder=BONUS_INT), mesh))
        if R_inv < 0 or R_inv < 1e-30:
            print(f"  [{n}] energy_norm pathological ({R_inv:.2e}) — STOP",
                  flush=True)
            break
        Rn = 1.0 / R_inv

        # Cross-products with stored J_history
        cross = []
        for gfJ_m in J_history:
            xnm = float(Integrate(gfJ * gfJ_m * sigma_inv_cf
                                  * dx(bonus_intorder=BONUS_INT), mesh))
            cross.append(xnm)
        max_cross = max((abs(c) for c in cross), default=0.0)
        rel_drift = max_cross / R_inv

        # Update Apot: gfApot += Rn * gfA (in place, GridFunction arithmetic)
        gfApot.vec.data += Rn * gfA.vec

        Ln_int = float(Integrate(gfJ * gfApot * dx(bonus_intorder=BONUS_INT),
                                 mesh))
        Ln = Rn * Ln_int

        signL = "+" if Ln > 0 else ("-" if Ln < 0 else "0")
        signR = "+" if Rn > 0 else ("-" if Rn < 0 else "0")
        tau_us = Ln/Rn*1e6 if Rn > 0 and Ln > 0 else None

        elapsed = time.time() - t_stage
        tau_str = f"{tau_us:.3f}" if tau_us is not None else "N/A"
        print(f"  [{n:2d}] drift={rel_drift:.3e}  En={R_inv:.3e}  "
              f"R={Rn:.3e}  L={Ln:.3e}({signL})  tau={tau_str}us  "
              f"({elapsed:.1f}s)", flush=True)

        diag.append({"n": n, "R_n": Rn, "L_n": Ln, "tau_us": tau_us,
                     "rel_drift": rel_drift, "energy_norm": R_inv,
                     "sign_L": signL, "sign_R": signR,
                     "stage_time_s": elapsed})

        # Detect sudden breakdown
        if breakdown_stage is None and len(diag) >= 2:
            prev = diag[-2]['rel_drift']
            if prev > 0 and rel_drift / prev > 100:
                print(f"      *** SUDDEN: drift jumped {rel_drift/prev:.1e}x ***",
                      flush=True)
                breakdown_stage = n

        # Snapshot JSON every SNAPSHOT_EVERY stages
        if (n + 1) % SNAPSHOT_EVERY == 0:
            snap = {"reference_tau_us": tau_TE_110_us,
                    "h_cond_mm": H_COND*1000, "order": ORDER, "ne": mesh.ne,
                    "stages_completed": n+1, "breakdown_stage": breakdown_stage,
                    "stages": diag}
            OUT_PATH.write_text(json.dumps(snap, indent=2))
            print(f"      [snapshot at n={n}]", flush=True)

        # Save current J as GridFunction in history (deep copy via GridFunction)
        gfJ_save = GridFunction(fes)
        gfJ_save.vec.data = gfJ.vec.data
        J_history.append(gfJ_save)

        # Stop conditions
        if Ln <= 0:
            print(f"      *** SIGN FLIP at n={n} ***", flush=True)
            if breakdown_stage is None:
                breakdown_stage = n
            # Continue past breakdown for a few stages
        if abs(Ln) > 1e10 or abs(Rn) > 1e15:
            print(f"      *** EXPLOSION at n={n} *** STOP", flush=True)
            break

        if abs(Ln) < 1e-30:
            print(f"      L too small, STOP", flush=True)
            break

        # Update J: J_{n+1} = J_n - sigma * Apot / Ln, project to GridFunction
        # Express as: tmp = sigma * Apot / Ln, then gfJ -= tmp
        gfJ.vec.data -= (sigma_Cu / Ln) * gfApot.vec

    # Final summary
    max_tau = max((d['tau_us'] for d in diag if d['tau_us']), default=0)
    print(f"\nFinal: ran {len(diag)} stages, breakdown_stage = {breakdown_stage}",
          flush=True)
    print(f"max(tau_n) = {max_tau:.3f} us (analytical {tau_TE_110_us:.3f}, "
          f"ratio {max_tau/tau_TE_110_us:.4f})", flush=True)

    snap = {"reference_tau_us": tau_TE_110_us,
            "h_cond_mm": H_COND*1000, "order": ORDER, "ne": mesh.ne,
            "stages_completed": len(diag), "breakdown_stage": breakdown_stage,
            "max_tau_us": max_tau, "stages": diag}
    OUT_PATH.write_text(json.dumps(snap, indent=2))
    print(f"Saved: {OUT_PATH}", flush=True)


if __name__ == "__main__":
    main()
