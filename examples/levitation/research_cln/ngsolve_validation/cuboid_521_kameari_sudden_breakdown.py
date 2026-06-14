"""Find the SUDDEN breakdown stage of 3D Kameari (closed PEC).

User's expert observation: "ある段数で途端に壊れる" — orthogonality
maintained for many stages then suddenly fails.

Push N to 30 with the optimal setup (order=3, nograds=True, tree-cotree,
bonus_intorder=8) and track:
  - rel_drift = max |<J_n, J_m/sigma>| / <J_n, J_n/sigma>
  - log(rel_drift) trajectory to find breakdown jump
  - sign of L_n, R_n
  - tau_n trajectory

Multiple mesh resolutions to see if breakdown stage scales with mesh:
  (A) h=0.30 mm (baseline, ne=1413)
  (B) h=0.20 mm (refined, ne=4385)
"""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace", line_buffering=True)

from netgen.occ import Box, Pnt, OCCGeometry
from ngsolve import (
    Mesh, HCurl, BilinearForm, LinearForm, GridFunction,
    CoefficientFunction, curl, dx, x, y, z,
    Integrate, TaskManager, ngsglobals,
)
from collections import deque
from math import pi
import numpy as np
import json, time
from pathlib import Path

mu0 = 4 * pi * 1e-7
sigma_Cu = 5.8e7
ax, ay, az = 5e-3, 2e-3, 1e-3
N_STAGES = 30
ORDER = 3
BONUS_INT = 8

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


def run_kameari_long(h_cond, n_stages):
    cuboid = Box(Pnt(-ax/2, -ay/2, -az/2), Pnt(ax/2, ay/2, az/2))
    cuboid.mat("conductor").bc("conductor_surface")
    cuboid.maxh = h_cond
    mesh = Mesh(OCCGeometry(cuboid).GenerateMesh(maxh=h_cond))
    print(f"  ne = {mesh.ne}, nv = {mesh.nv}")

    fes = HCurl(mesh, order=ORDER, dirichlet="conductor_surface", nograds=True)
    print(f"  HCurl order={ORDER} (nograds=True) ndof = {fes.ndof}")
    tree_edges = build_spanning_tree(mesh)
    fd = fes.FreeDofs()
    masked = 0
    for edge_nr in tree_edges:
        edge = mesh.edges[edge_nr]
        dofs = fes.GetDofNrs(edge)
        if dofs and fd[dofs[0]]:
            fd[dofs[0]] = False
            masked += 1
    print(f"  tree-cotree masked {masked} edges, active = {sum(fd)}")

    sigma_inv_cf = 1.0 / sigma_Cu
    A_ext = CoefficientFunction((-y, x, 0)) * 0.5

    u, v = fes.TnT()
    a = BilinearForm(fes)
    a += (1.0/mu0) * curl(u) * curl(v) * dx(bonus_intorder=BONUS_INT)
    print("  assembling + factor...")
    t0 = time.time()
    with TaskManager():
        a.Assemble()
        inv = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")
    print(f"    {time.time()-t0:.1f}s")

    fes_proj = HCurl(mesh, order=ORDER, dirichlet="conductor_surface", nograds=True)

    diag = []
    J_history = []
    J_cf = sigma_Cu * A_ext
    Apot = None
    breakdown_stage = None

    for n in range(n_stages):
        f = LinearForm(fes)
        f += J_cf * v * dx(bonus_intorder=BONUS_INT)
        with TaskManager():
            f.Assemble()
        gfA = GridFunction(fes)
        with TaskManager():
            gfA.vec.data = inv * f.vec

        R_inv = float(Integrate(J_cf * J_cf * sigma_inv_cf
                                * dx(bonus_intorder=BONUS_INT), mesh))
        if R_inv < 0 or R_inv < 1e-30:
            print(f"  [{n}] energy_norm pathological ({R_inv:.2e}), STOP")
            break
        Rn = 1.0 / R_inv

        gfJ_n = GridFunction(fes_proj)
        gfJ_n.Set(J_cf)

        cross = []
        for gfJ_m in J_history:
            xnm = float(Integrate(J_cf * gfJ_m * sigma_inv_cf
                                  * dx(bonus_intorder=BONUS_INT), mesh))
            cross.append(xnm)
        max_cross = max((abs(c) for c in cross), default=0.0)
        rel_drift = max_cross / R_inv

        if Apot is None:
            Apot = Rn * gfA
        else:
            Apot = Apot + Rn * gfA

        Ln_int = float(Integrate(J_cf * Apot * dx(bonus_intorder=BONUS_INT), mesh))
        Ln = Rn * Ln_int

        signL = "+" if Ln > 0 else ("-" if Ln < 0 else "0")
        signR = "+" if Rn > 0 else ("-" if Rn < 0 else "0")
        tau_us = Ln/Rn*1e6 if Rn > 0 and Ln > 0 else float('nan')

        diag.append({"n": n, "R_n": Rn, "L_n": Ln, "tau_us": tau_us,
                     "rel_drift": rel_drift, "energy_norm": R_inv,
                     "sign_L": signL, "sign_R": signR})

        # Print compactly
        print(f"  [{n:2d}] drift={rel_drift:.3e}  En={R_inv:.3e}  "
              f"R={Rn:.3e}  L={Ln:.3e}({signL})  tau={tau_us:.3f}us")

        # Detect sudden breakdown: drift jumps by 3+ orders of magnitude
        if breakdown_stage is None and len(diag) >= 2:
            prev_drift = diag[-2]['rel_drift']
            if prev_drift > 0 and rel_drift / prev_drift > 1000:
                print(f"      *** SUDDEN BREAKDOWN: drift jumped {rel_drift/prev_drift:.1e}x ***")
                # don't break, keep going to see post-breakdown behavior
                if breakdown_stage is None:
                    breakdown_stage = n
        if rel_drift > 0.01 and breakdown_stage is None:
            breakdown_stage = n

        J_history.append(gfJ_n)
        if Ln <= 0:
            print(f"      *** SIGN FLIP at n={n} ***")
            if breakdown_stage is None:
                breakdown_stage = n
            # keep going to see if it recovers
        if abs(Ln) > 1e10 or abs(Rn) > 1e15:
            print(f"      *** EXPLOSION ***, STOP")
            break

        # Update J for next stage (defensive: only if Ln nonzero)
        if abs(Ln) < 1e-30:
            print(f"      L too small, STOP")
            break
        J_cf = J_cf - sigma_Cu * Apot / Ln

    return diag, breakdown_stage, mesh.ne


def main():
    ngsglobals.msg_level = 0
    print(f"=== Sudden breakdown hunt: 3D Kameari closed PEC, N={N_STAGES} ===")
    print(f"Reference tau = {tau_TE_110_us:.3f} us")
    print(f"Setup: order={ORDER}, nograds=True, tree-cotree, bonus_intorder={BONUS_INT}\n")

    all_results = []
    for h_cond in [0.30e-3]:  # only baseline first to get quick result
        print("=" * 78)
        print(f"h_cond = {h_cond*1000} mm")
        print("=" * 78)
        try:
            diag, bd, ne = run_kameari_long(h_cond, N_STAGES)
            print(f"\n  Breakdown stage detected: {bd}")
            print(f"  Total stages run: {len(diag)}")
            max_tau = max((d['tau_us'] for d in diag if d['tau_us']), default=0)
            print(f"  Max tau over all stages: {max_tau:.3f} us "
                  f"(analytical {tau_TE_110_us:.3f}, ratio {max_tau/tau_TE_110_us:.4f})")
            all_results.append({"h_cond_mm": h_cond*1000, "ne": ne,
                                "breakdown_stage": bd, "stages": diag,
                                "max_tau_us": max_tau})
        except MemoryError:
            print("OOM, skip")
        print()

    # Detailed drift trajectories printed
    out_path = Path(__file__).parent / "cuboid_521_kameari_sudden_breakdown.json"
    out = {"reference_tau_us": tau_TE_110_us, "setup": {
        "order": ORDER, "nograds": True, "tree_cotree": True, "bonus_intorder": BONUS_INT
    }, "results": all_results}
    out_path.write_text(json.dumps(out, indent=2))
    print(f"Results: {out_path}")

    print("\n" + "=" * 78)
    print("SUMMARY: where does Kameari suddenly break?")
    print("=" * 78)
    for r in all_results:
        bd = r['breakdown_stage']
        print(f"  h={r['h_cond_mm']} mm, ne={r['ne']}: breakdown at n* = "
              f"{bd if bd is not None else '(none in '+str(N_STAGES)+')'}, "
              f"max tau = {r['max_tau_us']:.3f} us")


if __name__ == "__main__":
    main()
