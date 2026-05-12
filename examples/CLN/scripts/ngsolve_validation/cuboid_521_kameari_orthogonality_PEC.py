"""Where does Schmidt orthogonality break in 3D Kameari?

Closed PEC cuboid (the BEST case — known healthy from analytical reference).
Run N=15 stages and track per-stage:
  - Schmidt drift |<J_n, J_m/sigma>| for ALL m < n
  - Energy norm <J_n, J_n/sigma> trajectory
  - Gram matrix condition number through stage n
  - Sign of R_n, L_n
  - tau_n = L_n/R_n trajectory

Compare three mesh/order configurations to map mesh-quality dependence:
  (1) h_cond = 0.30 mm, order = 2  (baseline)
  (2) h_cond = 0.20 mm, order = 2  (refined)
  (3) h_cond = 0.30 mm, order = 3  (higher order)

Tells us whether the breakdown is:
  - Mesh-resolution limited (refinement helps)
  - Iteration-fundamental (refinement doesn't help, breaks at same stage n*)
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
import json, time
from pathlib import Path

mu0 = 4 * pi * 1e-7
sigma_Cu = 5.8e7
ax, ay, az = 5e-3, 2e-3, 1e-3
V_cond = ax * ay * az
N_STAGES = 15

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


def build_geo(h_cond):
    cuboid = Box(Pnt(-ax/2, -ay/2, -az/2), Pnt(ax/2, ay/2, az/2))
    cuboid.mat("conductor").bc("conductor_surface")
    cuboid.maxh = h_cond
    return OCCGeometry(cuboid)


def kameari_with_orthogonality_tracking(mesh, n_stages, order):
    fes = HCurl(mesh, order=order, dirichlet="conductor_surface", nograds=True)
    print(f"  HCurl order={order} ndof = {fes.ndof}")

    tree_edges = build_spanning_tree(mesh)
    fd = fes.FreeDofs()
    masked = 0
    for edge_nr in tree_edges:
        edge = mesh.edges[edge_nr]
        dofs = fes.GetDofNrs(edge)
        if dofs and fd[dofs[0]]:
            fd[dofs[0]] = False
            masked += 1
    print(f"  tree-cotree masked {masked} edges")

    sigma_inv_cf = 1.0 / sigma_Cu
    A_ext = CoefficientFunction((-y, x, 0)) * 0.5

    u, v = fes.TnT()
    a = BilinearForm(fes)
    a += (1.0/mu0) * curl(u) * curl(v) * dx
    print("  assembling + factorizing...")
    t0 = time.time()
    with TaskManager():
        a.Assemble()
        inv = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")
    print(f"  factor: {time.time()-t0:.1f}s")

    # Project J_cf to a GridFunction at each stage so we can store it
    fes_proj = HCurl(mesh, order=order, dirichlet="conductor_surface", nograds=True)

    diag = []
    J_history = []  # list of GridFunctions for J_0, J_1, ...
    J_cf = sigma_Cu * A_ext
    Apot = None

    for n in range(n_stages):
        f = LinearForm(fes)
        f += J_cf * v * dx
        with TaskManager():
            f.Assemble()
        gfA = GridFunction(fes)
        with TaskManager():
            gfA.vec.data = inv * f.vec

        # R_n
        R_inv = float(Integrate(J_cf * J_cf * sigma_inv_cf * dx, mesh))
        if R_inv < 0:
            print(f"  [{n}] energy_norm < 0 ({R_inv:.2e}) — STOP")
            break
        if R_inv < 1e-30:
            print(f"  [{n}] energy_norm too small")
            break
        Rn = 1.0 / R_inv

        # Project current J_cf to a GridFunction (for cross-products)
        gfJ_n = GridFunction(fes_proj)
        try:
            gfJ_n.Set(J_cf)
        except Exception as e:
            print(f"  [{n}] J projection failed: {e}")
            break

        # Schmidt drift: cross products with all previous J_m
        cross = []
        for m, gfJ_m in enumerate(J_history):
            xnm = float(Integrate(J_cf * gfJ_m * sigma_inv_cf * dx, mesh))
            cross.append(xnm)
        max_cross_abs = max((abs(c) for c in cross), default=0.0)
        # Relative drift: max |<J_n, J_m/sigma>| / <J_n, J_n/sigma>
        rel_drift = max_cross_abs / R_inv if R_inv > 0 else float('nan')

        if Apot is None:
            Apot = Rn * gfA
        else:
            Apot = Apot + Rn * gfA

        Ln_int = float(Integrate(J_cf * Apot * dx, mesh))
        Ln = Rn * Ln_int
        tau_us = Ln/Rn * 1e6 if Rn > 0 and Ln > 0 else float('nan')

        signL = "+" if Ln > 0 else ("-" if Ln < 0 else "0")
        signR = "+" if Rn > 0 else ("-" if Rn < 0 else "0")

        print(f"  [{n:2d}] R={Rn:.3e}, L={Ln:.3e}, tau={tau_us:.3f}us, "
              f"En={R_inv:.3e}, max|⟨J,J⟩_rel|={rel_drift:.2e}  ({signL})")

        diag.append({
            "n": n, "R_n": Rn, "L_n": Ln,
            "tau_us": tau_us if not (tau_us != tau_us) else None,
            "energy_norm": R_inv,
            "rel_drift": rel_drift,
            "max_cross_abs": max_cross_abs,
            "cross_inner": cross,
            "sign_L": signL, "sign_R": signR,
        })

        J_history.append(gfJ_n)

        if Ln <= 0:
            print(f"  ** L_n <= 0, sign flip at stage {n} **")
            break

        J_cf = J_cf - sigma_Cu * Apot / Ln

    return diag


def detect_breakdown(diag):
    """Find first stage where orthogonality drift exceeds 1%."""
    drift_threshold = 1e-2  # 1% of energy norm
    for d in diag:
        if d['rel_drift'] > drift_threshold:
            return d['n']
    return None


def run_config(name, h_cond, order, n_stages):
    print("\n" + "=" * 78)
    print(f"CONFIG: {name}  (h_cond = {h_cond*1000} mm, order = {order})")
    print("=" * 78)
    geo = build_geo(h_cond)
    print("Generating mesh...")
    t0 = time.time()
    mesh = Mesh(geo.GenerateMesh(maxh=h_cond))
    print(f"  ne = {mesh.ne}, nv = {mesh.nv}  ({time.time()-t0:.1f}s)")
    diag = kameari_with_orthogonality_tracking(mesh, n_stages, order)
    bd = detect_breakdown(diag)
    print(f"\n  Breakdown stage (rel_drift > 1%): n* = {bd}")
    if diag:
        max_tau = max((d['tau_us'] for d in diag if d['tau_us']), default=0)
        max_n = next((d['n'] for d in diag if d['tau_us'] == max_tau), None)
        print(f"  max tau_n = {max_tau:.3f} us at stage {max_n}")
        print(f"  Analytical tau_TE_z(1,1,0) = {tau_TE_110_us:.3f} us")
        print(f"  ratio max(tau)/tau_anal = {max_tau/tau_TE_110_us:.4f}")
    return {"config": name, "h_cond_mm": h_cond*1000, "order": order,
            "ne": mesh.ne, "nv": mesh.nv,
            "breakdown_stage": bd, "stages": diag}


def main():
    ngsglobals.msg_level = 0
    print("=== Kameari orthogonality tracking on closed-PEC cuboid 5x2x1 ===")
    print(f"Reference: tau_TE_z(1,1,0) = {tau_TE_110_us:.3f} us")
    print(f"N_STAGES = {N_STAGES}")

    configs = [
        ("baseline",    0.30e-3, 2),
        ("refined",     0.20e-3, 2),
        ("higher_ord",  0.30e-3, 3),
    ]

    results = []
    for name, h_cond, order in configs:
        try:
            r = run_config(name, h_cond, order, N_STAGES)
            results.append(r)
        except MemoryError:
            print(f"  ** OOM for {name} **")
        except Exception as e:
            print(f"  ** {name} failed: {e} **")

    print("\n" + "=" * 78)
    print("MESH-DEPENDENCE OF BREAKDOWN STAGE")
    print("=" * 78)
    print(f"{'config':>15} {'ne':>8} {'order':>5} {'breakdown n*':>14} "
          f"{'max tau [us]':>14} {'tau_anal':>12}")
    for r in results:
        bd = r['breakdown_stage']
        bd_s = f"{bd}" if bd is not None else "(none in 15)"
        max_tau = max((d['tau_us'] for d in r['stages'] if d['tau_us']), default=0)
        print(f"{r['config']:>15} {r['ne']:>8} {r['order']:>5} {bd_s:>14} "
              f"{max_tau:>14.3f} {tau_TE_110_us:>12.3f}")

    out = {
        "reference": {"tau_TE_110_analytical_us": tau_TE_110_us},
        "configurations": results,
    }
    out_path = Path(__file__).parent / "cuboid_521_kameari_orthogonality_PEC.json"
    out_path.write_text(json.dumps(out, indent=2))
    print(f"\nResults: {out_path}")


if __name__ == "__main__":
    main()
