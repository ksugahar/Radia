"""3D Kameari iteration on the VACUUM problem (uniform applied B_z, A-method).

Open-domain inner solver: NGSolve air-box truncation (AIR_SCALE = 5).
Same quality control as the closed-PEC isolation test:
  - HCurl order = 3, nograds = True
  - BFS spanning-tree gauge mask (interior + boundary)
  - bonus_intorder = 8
N = 30 stages, periodic snapshots every 3.

Track:
  - Schmidt drift |<J_n, J_m/sigma>|_cond / <J_n, J_n/sigma>_cond
  - tau_n = L_n/R_n
  - sign of L_n, R_n
  - exponential drift onset stage

Compare with:
  - air-box freq sweep (this same mesh, real Y(iw) from direct solve)
  - ELF + Vector Fitting tau_lead = 11.51 us
  - closed PEC analogue at N=11 healthy, N>=25 breakdown
"""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace", line_buffering=True)
print("Starting...", flush=True)

from netgen.occ import Box, Pnt, OCCGeometry, Glue
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
V_cond = ax * ay * az

AIR_SCALE = 5
H_COND = 0.30e-3
H_AIR = 1.5e-3
ORDER = 3
N_STAGES = 30
BONUS_INT = 8
SNAPSHOT_EVERY = 3

OUT_PATH = Path(__file__).parent / "cuboid_521_kameari_vacuum_v2.json"

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


def build_geometry():
    cuboid = Box(Pnt(-ax/2, -ay/2, -az/2), Pnt(ax/2, ay/2, az/2))
    cuboid.mat("conductor").bc("conductor_surface")
    cuboid.maxh = H_COND
    Bx, By, Bz = AIR_SCALE * ax / 2, AIR_SCALE * ay / 2, AIR_SCALE * az / 2
    air = Box(Pnt(-Bx, -By, -Bz), Pnt(Bx, By, Bz))
    air.mat("air").bc("outer_box")
    air.maxh = H_AIR
    return OCCGeometry(Glue([air - cuboid, cuboid]))


def main():
    ngsglobals.msg_level = 0
    print(f"=== 3D Kameari on VACUUM cuboid (A-method, air-box AIR_SCALE={AIR_SCALE}) ===",
          flush=True)
    print(f"Setup: order={ORDER}, nograds=True, tree-cotree, bonus={BONUS_INT}",
          flush=True)
    print(f"Reference Parseval R_0 = {R0_anal:.6e}", flush=True)
    print(f"ELF Foster B_z tau_lead = 11.51 us\n", flush=True)

    geo = build_geometry()
    print("Generating mesh...", flush=True)
    mesh = Mesh(geo.GenerateMesh(maxh=H_AIR))
    print(f"  ne = {mesh.ne}, nv = {mesh.nv}, mats = {mesh.GetMaterials()}",
          flush=True)

    fes = HCurl(mesh, order=ORDER, dirichlet="outer_box", nograds=True)
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

    sigma_cf = mesh.MaterialCF({"conductor": sigma_Cu, "air": 0.0})
    sigma_inv_cf = mesh.MaterialCF({"conductor": 1.0/sigma_Cu, "air": 0.0})
    A_ext = CoefficientFunction((-y, x, 0)) * 0.5  # B_0 = 1, z-direction

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
    J_history = []

    gfJ = GridFunction(fes)
    gfJ.Set(sigma_cf * A_ext)

    gfApot = GridFunction(fes)
    gfApot.vec[:] = 0.0

    breakdown_stage = None

    for n in range(N_STAGES):
        t_stage = time.time()

        f = LinearForm(fes)
        f += gfJ * v * dx("conductor", bonus_intorder=BONUS_INT)
        with TaskManager():
            f.Assemble()
        gfA = GridFunction(fes)
        with TaskManager():
            gfA.vec.data = inv * f.vec

        R_inv = float(Integrate(gfJ * gfJ * sigma_inv_cf
                                * dx("conductor", bonus_intorder=BONUS_INT), mesh))
        if R_inv < 0 or R_inv < 1e-30:
            print(f"  [{n}] energy_norm pathological ({R_inv:.2e}), STOP",
                  flush=True)
            break
        Rn = 1.0 / R_inv

        cross = []
        for gfJ_m in J_history:
            xnm = float(Integrate(gfJ * gfJ_m * sigma_inv_cf
                                  * dx("conductor", bonus_intorder=BONUS_INT), mesh))
            cross.append(xnm)
        max_cross = max((abs(c) for c in cross), default=0.0)
        rel_drift = max_cross / R_inv

        gfApot.vec.data += Rn * gfA.vec

        Ln_int = float(Integrate(gfJ * gfApot
                                  * dx("conductor", bonus_intorder=BONUS_INT), mesh))
        Ln = Rn * Ln_int

        signL = "+" if Ln > 0 else ("-" if Ln < 0 else "0")
        signR = "+" if Rn > 0 else "-"
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

        if breakdown_stage is None and rel_drift > 0.01:
            breakdown_stage = n
            print(f"      *** 1% threshold crossed at n={n} ***", flush=True)

        if (n + 1) % SNAPSHOT_EVERY == 0:
            snap = {"reference_R0_Parseval": R0_anal,
                    "elf_foster_tau_lead_us": 11.51,
                    "AIR_SCALE": AIR_SCALE,
                    "h_cond_mm": H_COND*1000, "h_air_mm": H_AIR*1000,
                    "order": ORDER, "ne": mesh.ne,
                    "stages_completed": n+1,
                    "breakdown_stage": breakdown_stage,
                    "stages": diag}
            OUT_PATH.write_text(json.dumps(snap, indent=2))
            print(f"      [snapshot at n={n}]", flush=True)

        gfJ_save = GridFunction(fes)
        gfJ_save.vec.data = gfJ.vec.data
        J_history.append(gfJ_save)

        if Ln <= 0:
            print(f"      *** SIGN FLIP at n={n} ***", flush=True)
            if breakdown_stage is None:
                breakdown_stage = n
        if abs(Ln) > 1e10 or abs(Rn) > 1e15:
            print(f"      *** EXPLOSION ***, STOP", flush=True)
            break
        if abs(Ln) < 1e-30:
            print(f"      L too small, STOP", flush=True)
            break

        gfJ.vec.data -= (sigma_Cu / Ln) * gfApot.vec

    max_tau = max((d['tau_us'] for d in diag if d['tau_us']), default=0)
    print(f"\nFinal: {len(diag)} stages, breakdown_stage = {breakdown_stage}",
          flush=True)
    print(f"max(tau_n) = {max_tau:.3f} us  (ELF tau_lead = 11.51 us, ratio "
          f"{max_tau/11.51:.4f})", flush=True)

    snap = {"reference_R0_Parseval": R0_anal,
            "elf_foster_tau_lead_us": 11.51,
            "AIR_SCALE": AIR_SCALE,
            "h_cond_mm": H_COND*1000, "h_air_mm": H_AIR*1000,
            "order": ORDER, "ne": mesh.ne,
            "stages_completed": len(diag),
            "breakdown_stage": breakdown_stage,
            "max_tau_us": max_tau,
            "stages": diag}
    OUT_PATH.write_text(json.dumps(snap, indent=2))
    print(f"Saved: {OUT_PATH}", flush=True)


if __name__ == "__main__":
    main()
