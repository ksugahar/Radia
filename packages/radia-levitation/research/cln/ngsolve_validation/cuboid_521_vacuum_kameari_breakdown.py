"""Kameari iteration breakdown analysis on cuboid 5x2x1 vacuum problem.

Run Kameari recursion to N_STAGES = 12 and record per-stage diagnostics:
  - R_n, L_n, tau_n
  - ||J_n||^2_sigma  (energy norm, must monotonically decrease)
  - <J_n, J_m/sigma> for all m<n (Schmidt orthogonality drift)
  - Gram matrix condition number
  - sign of Ln, Rn (positivity of Cauer parameters)

Output JSON for downstream paper figure / breakdown stage detection.

Compare with:
  - ELF Foster ladder (cuboid_521_3dir_cauer.json): tau_lead = 11.51 us
  - PEC closed Kameari (existing _archive_closed_PEC_3D/cuboid_521_gauge_CLN_v3.py)
    which is reported to remain stable up to ~10 stages
"""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")

from netgen.occ import Box, Pnt, OCCGeometry, Glue
from ngsolve import (
    Mesh, HCurl, BilinearForm, LinearForm, GridFunction,
    CoefficientFunction, curl, dx, x, y, z,
    Integrate, TaskManager, ngsglobals, InnerProduct,
)
from collections import deque
import numpy as np
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
ORDER = 2
N_STAGES = 12

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


def kameari_breakdown(mesh, n_stages):
    fes = HCurl(mesh, order=ORDER, dirichlet="outer_box", nograds=True)
    print(f"  HCurl ndof = {fes.ndof}")

    tree_edges = build_spanning_tree(mesh)
    fd = fes.FreeDofs()
    masked = 0
    for edge_nr in tree_edges:
        edge = mesh.edges[edge_nr]
        dofs = fes.GetDofNrs(edge)
        if dofs and fd[dofs[0]]:
            fd[dofs[0]] = False
            masked += 1
    print(f"  Tree-cotree masked {masked} edges, active DoFs = {sum(fd)}")

    sigma_cf = mesh.MaterialCF({"conductor": sigma_Cu, "air": 0.0})
    sigma_inv_cf = mesh.MaterialCF({"conductor": 1.0/sigma_Cu, "air": 0.0})

    u, v = fes.TnT()
    A_ext = CoefficientFunction((-y, x, 0)) * 0.5  # B_0 = 1, normalized

    a_form = BilinearForm(fes)
    a_form += (1.0/mu0) * curl(u) * curl(v) * dx
    print("  Assembling curl-curl operator...")
    t0 = time.time()
    with TaskManager():
        a_form.Assemble()
    print(f"  Assembled in {time.time()-t0:.1f}s")

    print("  Factorizing (sparsecholesky)...")
    t0 = time.time()
    with TaskManager():
        inv = a_form.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")
    print(f"  Factored in {time.time()-t0:.1f}s")

    # Store J_n GridFunctions for cross-inner-product computation
    # (need these to evaluate <J_i, J_j/sigma> for Gram matrix)
    # In NGSolve, we can store J_n as CF expressions but for diagnostics
    # we need numerical evaluation of cross-products. Use GridFunction projection.
    fes_proj = HCurl(mesh, order=ORDER, dirichlet="outer_box", nograds=True)

    J_history = []   # list of GridFunctions representing J_n in conductor
    A_history = []   # list of A_n (raw, this stage solve)
    Apot_accum = None

    diag = []
    J_cf = sigma_cf * A_ext   # initial impressed source

    for n in range(n_stages):
        print(f"\n[Stage {n}] solve A_{n}...")
        t0 = time.time()
        f = LinearForm(fes)
        f += J_cf * v * dx("conductor")
        with TaskManager():
            f.Assemble()
        gfA_n = GridFunction(fes, name=f"A_{n}")
        with TaskManager():
            gfA_n.vec.data = inv * f.vec
        t_solve = time.time() - t0

        # Project current J_cf to a GridFunction so we can store and reuse
        gfJ_n = GridFunction(fes_proj)
        # Use L2 projection: solve (u, v) = (J_cf, v) for each new J
        # Cheap: GridFunction Set with CF
        try:
            gfJ_n.Set(J_cf, definedon=mesh.Materials("conductor"))
        except Exception:
            gfJ_n.Set(J_cf)

        # R_n = 1 / <J_n, J_n/sigma>_cond
        R_inv_int = float(Integrate(J_cf * J_cf * sigma_inv_cf * dx("conductor"), mesh))
        if abs(R_inv_int) < 1e-30:
            print(f"  R_inv too small ({R_inv_int:.2e}), stop")
            break
        Rn = 1.0 / R_inv_int

        # Tanimoto accumulated A_pot
        if Apot_accum is None:
            Apot_accum = Rn * gfA_n
        else:
            Apot_accum = Apot_accum + Rn * gfA_n

        Ln_int = float(Integrate(J_cf * Apot_accum * dx("conductor"), mesh))
        Ln = Rn * Ln_int

        tau_us = Ln / Rn * 1e6 if Rn > 0 and Ln > 0 else float('nan')

        # Schmidt orthogonality drift: <J_n, J_m/sigma> for m<n
        # In Kameari, J_n should be sigma-orthogonal to all J_m, m<n
        cross_inner = []
        for m, gfJ_m in enumerate(J_history):
            cross = float(Integrate(J_cf * gfJ_m * sigma_inv_cf * dx("conductor"), mesh))
            cross_inner.append(cross)

        # Sign flag
        signL = "+" if Ln > 0 else "-" if Ln < 0 else "0"
        signR = "+" if Rn > 0 else "-" if Rn < 0 else "0"

        print(f"  R_{n} = {Rn:.6e} ({signR})")
        print(f"  L_{n} = {Ln:.6e} ({signL}, {Ln*1e9:.4e} nH)")
        print(f"  tau_{n} = {tau_us:.4f} us")
        print(f"  <J_n, J_n/sigma> = {R_inv_int:.6e}  (energy norm, should decrease)")
        if cross_inner:
            max_cross = max(abs(c) for c in cross_inner)
            print(f"  max |<J_n, J_m/sigma>| (m<n) = {max_cross:.6e}  (Schmidt drift)")

        diag.append({
            "n": n,
            "R_n": Rn,
            "L_n": Ln,
            "tau_us": tau_us if not np.isnan(tau_us) else None,
            "energy_norm": R_inv_int,
            "cross_inner_with_prev": cross_inner,
            "sign_L": signL,
            "sign_R": signR,
            "solve_time_s": t_solve,
        })

        J_history.append(gfJ_n)

        if Ln <= 0 or Rn <= 0:
            print(f"  ** BREAKDOWN at stage {n}: L_n or R_n <= 0 **")
            # Continue one more stage to see how it evolves
            if n >= 2:
                break

        # Update J for next stage
        J_cf = J_cf - sigma_cf * Apot_accum / Ln

    # Build Gram matrix from cross_inner data
    # G[i,j] = <J_i, J_j/sigma> (only stored upper triangle from cross_inner_with_prev)
    nstages_done = len(diag)
    G = np.zeros((nstages_done, nstages_done))
    for i, d in enumerate(diag):
        G[i, i] = d["energy_norm"]
        for j, c in enumerate(d["cross_inner_with_prev"]):
            G[i, j] = c
            G[j, i] = c

    # Condition number of Gram matrix at each stage
    gram_cond = []
    for k in range(1, nstages_done + 1):
        Gk = G[:k, :k]
        try:
            evs = np.linalg.eigvalsh(Gk)
            cond = evs.max() / evs.min() if evs.min() > 0 else float('inf')
        except Exception:
            cond = float('nan')
        gram_cond.append(cond)
        diag[k-1]["gram_cond_through_n"] = cond

    return diag, G


def main():
    ngsglobals.msg_level = 1
    print(f"=== Vacuum cuboid Kameari BREAKDOWN analysis ===")
    print(f"  Conductor: {ax*1000}x{ay*1000}x{az*1000} mm Cu")
    print(f"  Air box: {AIR_SCALE}x conductor")
    print(f"  h_cond = {H_COND*1000} mm, h_air = {H_AIR*1000} mm, order = {ORDER}")
    print(f"  N_stages = {N_STAGES}")
    print()

    geo = build_geometry()
    print("Generating mesh...")
    mesh = Mesh(geo.GenerateMesh(maxh=H_AIR))
    print(f"  ne = {mesh.ne}, nv = {mesh.nv}")
    print(f"\nReference R_0 (Parseval) = {R0_anal:.6e}")
    print(f"ELF Foster B_z tau_lead   = 11.51 us\n")

    diag, G = kameari_breakdown(mesh, N_STAGES)

    print("\n" + "=" * 78)
    print("BREAKDOWN SUMMARY TABLE")
    print("=" * 78)
    print(f"{'n':>3} {'R_n':>12} {'L_n [nH]':>12} {'tau [us]':>10} "
          f"{'<J,J/σ>':>12} {'max |⟨J,J_m/σ⟩|':>16} {'gram cond':>12}")
    print("-" * 78)
    for d in diag:
        tau_s = f"{d['tau_us']:.3f}" if d['tau_us'] is not None else "  N/A"
        max_cross = (max(abs(c) for c in d['cross_inner_with_prev'])
                     if d['cross_inner_with_prev'] else 0.0)
        print(f"{d['n']:>3} {d['R_n']:>12.4e} {d['L_n']*1e9:>12.4e} {tau_s:>10} "
              f"{d['energy_norm']:>12.4e} {max_cross:>16.4e} "
              f"{d['gram_cond_through_n']:>12.4e}")

    # Detection: first stage where L_n < 0 OR Schmidt drift exceeds 1% of energy norm
    breakdown_stage = None
    for d in diag:
        max_cross = (max(abs(c) for c in d['cross_inner_with_prev'])
                     if d['cross_inner_with_prev'] else 0.0)
        if d['L_n'] <= 0 or d['R_n'] <= 0:
            breakdown_stage = d['n']
            print(f"\n** BREAKDOWN at stage n* = {breakdown_stage}: "
                  f"L_n or R_n changed sign **")
            break
        if max_cross > 0.01 * d['energy_norm']:
            breakdown_stage = d['n']
            print(f"\n** BREAKDOWN at stage n* = {breakdown_stage}: "
                  f"Schmidt drift exceeded 1% of energy norm **")
            break
    if breakdown_stage is None:
        print(f"\n** No breakdown detected within {N_STAGES} stages **")

    # Compare leading tau with ELF Foster
    if diag and diag[0]['tau_us']:
        tau0 = diag[0]['tau_us']
        elf_tau = 11.51
        ratio = tau0 / elf_tau
        print(f"\nLeading tau_0 (Kameari, vacuum) = {tau0:.3f} us")
        print(f"  vs ELF Foster tau_lead         = {elf_tau:.3f} us")
        print(f"  ratio Kameari/ELF              = {ratio:.4f}")

    out = {
        "method": "Kameari iteration on cuboid 5x2x1 vacuum (air box truncation)",
        "geometry_mm": [ax*1000, ay*1000, az*1000],
        "air_scale": AIR_SCALE,
        "h_cond_mm": H_COND*1000,
        "h_air_mm": H_AIR*1000,
        "order": ORDER,
        "ne": mesh.ne,
        "R0_analytic": R0_anal,
        "n_stages_attempted": N_STAGES,
        "n_stages_completed": len(diag),
        "breakdown_stage": breakdown_stage,
        "elf_foster_tau_lead_us": 11.51,
        "stages": diag,
        "gram_matrix": G.tolist(),
    }
    out_path = Path(__file__).parent / "cuboid_521_vacuum_kameari_breakdown.json"
    out_path.write_text(json.dumps(out, indent=2))
    print(f"\nResults: {out_path}")


if __name__ == "__main__":
    main()
