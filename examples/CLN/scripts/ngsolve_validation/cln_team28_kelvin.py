"""cln_team28_kelvin.py — Sugahara TEAM 28 axisymmetric CLN ported to
3D NGSolve with Kelvin two-sphere mesh.

Algorithm (from W:/30_CauerLadderNetwork/2020_11_04_線形のCLNの練習/
              2023_01_05_TEAM_28_軸対称/CLN_COMSOL.m):

  Initial: J_1 = -σ A_ext  (where A_ext is applied vector potential, here
                            (-y/2, x/2, 0)·B_0 for uniform B_z step)
  G_n = ∫_cond (J_n²/σ) dV    (Joule loss measure)
  R_n = 1/G_n
  Solve K_A · A_phi_n = J_n   on full mesh (cond + air + Kelvin)
  A_acc_n = A_acc_{n-1} + R_n × A_phi_n   (CLN vector potential, accumulated)
  L_n = R_n × ∫_cond J_n · A_acc_n dV
  τ_n = L_n / R_n = ∫_cond J_n · A_acc_n dV   (COMSOL/Sugahara CLN formula)
  Update: J_{n+1} = J_n - σ A_acc_n / L_n   (next current density)

NOTE: This is single-space A method (no T variable), matching COMSOL TEAM 28.
Helmholtz-Hodge correction on A_s is applied for non-axisymmetric geometries.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from collections import deque
from math import pi
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace", line_buffering=True)
sys.path.insert(0, r"S:/Radia/01_GitHub/src/radia")

import numpy as np
from netgen.occ import Cylinder, Box, Sphere, Pnt, X, Y, Z, OCCGeometry
from ngsolve import (
    Mesh, HCurl, H1, VectorH1, Periodic, BilinearForm, LinearForm,
    GridFunction, CoefficientFunction, curl, grad, dx, x, y, z,
    Integrate, InnerProduct, sqrt, TaskManager, ngsglobals,
)

from kelvin_geometry import add_kelvin_exterior_domain
from kelvin_material import (
    make_kelvin_nu_cf,
    make_reduced_potential_background_cf,
    NU_0,
)


# === Physical parameters ====================================================
mu0 = 4 * pi * 1e-7
sigma_Cu = 5.8e7
B0 = 1.0

# Kelvin parameters (match v5)
R_K_DEFAULT = 50e-3
OFFSET_DEFAULT = (5 * R_K_DEFAULT, 0.0, 0.0)

# FE parameters
ORDER = 2
H_COND = 1.0e-3
H_AIR = 12.0e-3
N_STAGES = 6
BONUS_INT = 8

GEOMETRIES = {
    "cylinder":   {"kind": "cyl",  "R": 10e-3, "t": 2e-3},
    "A1_square":  {"kind": "rect", "Lx": 17.72e-3, "Ly": 17.72e-3,
                    "Lz": 2e-3},
    "cuboid_521": {"kind": "rect", "Lx": 5e-3, "Ly": 2e-3, "Lz": 1e-3},
    "sphere":     {"kind": "sphere", "R": 10e-3},
}


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


def build_geo(geom_name, R_K, OFFSET):
    p = GEOMETRIES[geom_name]
    if p["kind"] == "cyl":
        cond = Cylinder(Pnt(0, 0, -p["t"] / 2), Z, r=p["R"], h=p["t"])
    elif p["kind"] == "sphere":
        cond = Sphere(Pnt(0, 0, 0), p["R"])
    else:
        Lx, Ly, Lz = p["Lx"], p["Ly"], p["Lz"]
        cond = Box(Pnt(-Lx/2, -Ly/2, -Lz/2),
                   Pnt( Lx/2,  Ly/2,  Lz/2))
    cond.mat("conductor")
    cond.maxh = H_COND
    inner_sphere = Sphere(Pnt(0, 0, 0), R_K)
    for f in inner_sphere.faces:
        f.name = "kelvin_int"
    inner_sphere.maxh = H_AIR
    inner_air = inner_sphere - cond
    inner_air.name = "air"
    inner_air.maxh = H_AIR
    geo, _ = add_kelvin_exterior_domain(
        [inner_air, cond], offset=OFFSET, R_K=R_K, inner_maxh=H_AIR,
    )
    return OCCGeometry(geo)


def main():
    global H_COND, H_AIR
    parser = argparse.ArgumentParser()
    parser.add_argument("--geom", choices=list(GEOMETRIES.keys()),
                         default="cylinder")
    parser.add_argument("--stages", type=int, default=N_STAGES)
    parser.add_argument("--order", type=int, default=ORDER)
    parser.add_argument("--h-cond", type=float, default=H_COND * 1000,
                         help="conductor mesh size in mm (default 1.0)")
    parser.add_argument("--h-air", type=float, default=H_AIR * 1000,
                         help="air mesh size in mm (default 12.0)")
    parser.add_argument("--r-k", type=float, default=R_K_DEFAULT * 1000,
                         help="Kelvin sphere radius in mm (default 50.0)")
    parser.add_argument("--out-tag", type=str, default=None)
    args = parser.parse_args()

    n_stages = args.stages
    ORDER_use = args.order
    h_cond_use = args.h_cond * 1e-3   # mm -> m
    h_air_use = args.h_air * 1e-3
    H_COND = h_cond_use
    H_AIR = h_air_use
    R_K = args.r_k * 1e-3
    OFFSET = (5 * R_K, 0.0, 0.0)

    ngsglobals.msg_level = 0
    print("=" * 78, flush=True)
    print(" CLN TEAM 28 (COMSOL/Sugahara) + Kelvin two-sphere", flush=True)
    print(f" geometry: {args.geom}", flush=True)
    print("=" * 78, flush=True)
    p = GEOMETRIES[args.geom]
    if p["kind"] == "cyl":
        print(f"  Cylinder R={p['R']*1000} mm, t={p['t']*1000} mm",
              flush=True)
    elif p["kind"] == "sphere":
        print(f"  Sphere R={p['R']*1000} mm", flush=True)
    else:
        print(f"  Rect prism {p['Lx']*1000}x{p['Ly']*1000}x{p['Lz']*1000} mm",
              flush=True)
    print(f"  σ = {sigma_Cu:.2e} S/m, B0 = {B0} T", flush=True)
    print(f"  Kelvin: R_K = {R_K*1000} mm, OFFSET = {OFFSET}", flush=True)
    print(f"  ORDER={ORDER_use}, H_COND={h_cond_use*1000}mm,"
          f" H_AIR={H_AIR*1000}mm, N_STAGES={n_stages}",
          flush=True)
    print()

    print(" Building geometry + mesh ...", flush=True)
    geo = build_geo(args.geom, R_K, OFFSET)
    t0 = time.time()
    with TaskManager():
        ngmesh = geo.GenerateMesh(maxh=H_AIR, grading=0.5)
    mesh = Mesh(ngmesh)
    mesh.Curve(ORDER_use)
    print(f"  GenerateMesh + Curve(p={ORDER_use}): {time.time()-t0:.1f} s",
          flush=True)
    print(f"  ne={mesh.ne}, nv={mesh.nv}", flush=True)
    print(f"  materials={mesh.GetMaterials()}", flush=True)

    sigma_cf = mesh.MaterialCF({"conductor": sigma_Cu, "air": 0.0,
                                 "kelvin": 0.0})
    nu_cf = make_kelvin_nu_cf(mesh, R_K, OFFSET, nu_0=NU_0,
                               kelvin_mats=("kelvin",))

    A_s_cf_raw = make_reduced_potential_background_cf(
        mesh,
        F_inner_factory=lambda xc, yc, zc: CoefficientFunction(
            (-yc, xc, 0)) * (B0 / 2.0),
        R_K=R_K, offset=OFFSET, kelvin_mats=("kelvin",), dim=3,
    )

    # fesA: HCurl on FULL mesh (cond + air + Kelvin)
    fesA = Periodic(HCurl(mesh, order=ORDER_use, dirichlet_bbnd="GND",
                           nograds=True))
    print(f"  fesA HCurl ndof = {fesA.ndof}", flush=True)

    tree_edges = build_spanning_tree(mesh)
    fdA = fesA.FreeDofs()
    masked = 0
    for edge_nr in tree_edges:
        edge = mesh.edges[edge_nr]
        dofs = fesA.GetDofNrs(edge)
        if dofs and fdA[dofs[0]]:
            fdA[dofs[0]] = False
            masked += 1
    print(f"  fesA tree-cotree masked {masked}, active = {sum(fdA)}",
          flush=True)

    A_tr, A_te = fesA.TnT()

    # σ-weighted H-H projection of A_s (only meaningful for non-axisymm)
    # CRITICAL: ORDER_PHI must MATCH HCurl ORDER (not ORDER+1) — Tier 1 fix
    # from project_3d_hcurl_hiruma_hodge_success.md memory:
    # higher-order H1 over-projects against grad-zero leak in HCurl.
    print(" σ-weighted Helmholtz-Hodge projection of A_s on cond"
          f" (ORDER_PHI={ORDER_use} matching HCurl) ...",
          flush=True)
    fes_phi = H1(mesh, order=ORDER_use,
                  definedon=mesh.Materials("conductor"))
    p_tr, p_te = fes_phi.TnT()
    a_phi = BilinearForm(fes_phi)
    a_phi += sigma_cf * grad(p_tr) * grad(p_te) * dx(
        "conductor", bonus_intorder=BONUS_INT)
    f_phi = LinearForm(fes_phi)
    f_phi += sigma_cf * A_s_cf_raw * grad(p_te) * dx(
        "conductor", bonus_intorder=BONUS_INT)
    with TaskManager():
        a_phi.Assemble()
        f_phi.Assemble()
        fd_phi = fes_phi.FreeDofs()
        pinned_phi = -1
        for ip in range(fes_phi.ndof):
            if fd_phi[ip]:
                fd_phi[ip] = False
                pinned_phi = ip
                break
        if pinned_phi >= 0:
            f_phi.vec[pinned_phi] = 0.0
        try:
            inv_phi = a_phi.mat.Inverse(fd_phi, inverse="pardiso")
        except Exception:
            inv_phi = a_phi.mat.Inverse(fd_phi, inverse="sparsecholesky")
        gf_phi = GridFunction(fes_phi)
        gf_phi.vec.data = inv_phi * f_phi.vec
    A_s_cf = A_s_cf_raw - grad(gf_phi)
    print(f"  pinned phi vertex {pinned_phi}", flush=True)

    # === DIAGNOSTICS (User questions: H-H ok? GND ok? solver ok?) ==========
    # 1. A_s H-H: |gf_phi| should be tiny if σA_s already div-free
    norm_phi = float(Integrate(gf_phi * gf_phi * dx("conductor"), mesh))
    norm_As_sq = float(Integrate(A_s_cf_raw * A_s_cf_raw * dx("conductor"), mesh))
    print(f"  [DIAG] A_s H-H: ∫|phi|²/∫|A_s|² = {norm_phi/norm_As_sq:.3e}"
          f" (should be ≪ 1 since σA_s = (-y,x,0) div-free)", flush=True)
    # 2. GND boundary effectiveness: count Dirichlet-eliminated DOFs
    dir_removed = sum(1 for f in fesA.FreeDofs() if not f) - masked
    print(f"  [DIAG] HCurl GND: dirichlet_bbnd='GND' removed {dir_removed} DOFs"
          f" (HCurl has NO vertex DOFs — GND vertex may be ineffective!)",
          flush=True)
    # 3. Galerkin identity check after stage 0 K_A solve
    print(f"  [DIAG] (Galerkin check ∫A·J = ∫ν|curl A|² will run at stage 0)",
          flush=True)

    # === ZZ error estimator setup ======================================
    # Recover B = curl(A_phi) onto VectorH1(order+1) via L²-projection.
    # Element-wise ||curl(A) - B_recovered||_L² indicates discretization
    # error & breakdown of the iteration (large at sharp corners).
    fes_B = VectorH1(mesh, order=ORDER_use + 1)
    u_B, v_B = fes_B.TnT()
    m_B_form = BilinearForm(fes_B, symmetric=True, check_unused=False)
    m_B_form += InnerProduct(u_B, v_B) * dx(bonus_intorder=BONUS_INT)
    print(" Assembling B-recovery mass matrix (VectorH1 order+1) ...",
          flush=True)
    with TaskManager():
        m_B_form.Assemble()
        m_B_inv = m_B_form.mat.Inverse(fes_B.FreeDofs(), inverse="pardiso")
    zz_per_stage = []   # store ZZ error per stage

    # Assemble + factor K_A (full mesh, with Kelvin)
    a_A = BilinearForm(fesA, check_unused=False)
    a_A += nu_cf * curl(A_tr) * curl(A_te) * dx(bonus_intorder=BONUS_INT)
    print(" Assembling K_A (full mesh + Kelvin nu) ...", flush=True)
    t0 = time.time()
    with TaskManager():
        a_A.Assemble()
        try:
            inv_A = a_A.mat.Inverse(fdA, inverse="pardiso")
        except Exception:
            inv_A = a_A.mat.Inverse(fdA, inverse="sparsecholesky")
    print(f"  K_A assemble + factor: {time.time()-t0:.1f} s", flush=True)

    # === Iteration ======================================================
    # Initial: J_1 = -σ A_s' (= initial eddy current pattern)
    # Note: sign convention follows COMSOL TEAM 28 (negative because step B
    # induces opposing current). We track magnitudes; sign conventions
    # cancel in τ = L/R = ∫J·A.
    J_n_cf = sigma_cf * A_s_cf  # CoefficientFunction (positive sign for our use)

    Apot_acc = GridFunction(fesA, name="Apot_acc")
    Apot_acc.vec[:] = 0.0  # Start with zero accumulator

    # Store all previous J_k for explicit Schmidt M⁻¹-orthogonalization in 3D.
    # Inner product: (J_a, J_b)_σ⁻¹ = ∫_cond J_a · J_b / σ dx (Joule-loss).
    # In axisym this is implicit; in 3D HCurl gauge errors break orthogonality.
    J_history_cf = []  # list of CoefficientFunction (frozen J_k expressions)
    G_history = []     # list of (J_k, J_k)_σ⁻¹ values

    Rn = []
    Ln = []
    diag = []

    for nStage in range(n_stages):
        t_stage = time.time()
        print(f"\n Stage {nStage+1}/{n_stages}", flush=True)

        # Schmidt disabled — snapshot fix should eliminate the need for it.

        # G_n = ∫_cond (J_n²/σ) dV  →  R_n = 1/G_n
        G_n = float(Integrate(
            J_n_cf * J_n_cf / sigma_Cu * dx(
                "conductor", bonus_intorder=BONUS_INT), mesh))
        R_n = 1.0 / G_n
        Rn.append(R_n)
        print(f"  G_{nStage} = ∫(J²/σ) = {G_n:.6e},  R_{nStage} = {R_n:.6e}",
              flush=True)

        # Save (orthogonalized) J_n to history for future Schmidt steps
        J_history_cf.append(J_n_cf)
        G_history.append(G_n)

        # Solve K_A · A_phi_n = J_n  on full mesh
        f_A = LinearForm(fesA)
        f_A += A_te * J_n_cf * dx(
            "conductor", bonus_intorder=BONUS_INT)
        with TaskManager():
            f_A.Assemble()

        gfAphi = GridFunction(fesA, name=f"Aphi_{nStage}")
        with TaskManager():
            gfAphi.vec.data = inv_A * f_A.vec

        # === ZZ error estimator for A_phi at this stage =================
        # Recover B = curl(A_phi) on VectorH1(order+1) via L²-projection.
        # eta = ||curl(A_phi) - B_recovered||_L² gives discretization error.
        gf_B_smooth = GridFunction(fes_B, name=f"B_smooth_{nStage}")
        f_B_form = LinearForm(fes_B)
        f_B_form += InnerProduct(curl(gfAphi), v_B) * dx(
            bonus_intorder=BONUS_INT)
        with TaskManager():
            f_B_form.Assemble()
            gf_B_smooth.vec.data = m_B_inv * f_B_form.vec
        # ZZ error: total + on cond + on air
        diff_cf = curl(gfAphi) - gf_B_smooth
        eta_total_sq = float(Integrate(
            InnerProduct(diff_cf, diff_cf) * dx(bonus_intorder=BONUS_INT),
            mesh))
        eta_cond_sq = float(Integrate(
            InnerProduct(diff_cf, diff_cf) * dx(
                "conductor", bonus_intorder=BONUS_INT), mesh))
        norm_B_total_sq = float(Integrate(
            InnerProduct(curl(gfAphi), curl(gfAphi)) * dx(
                bonus_intorder=BONUS_INT), mesh))
        eta_total = float(eta_total_sq) ** 0.5
        eta_cond = float(eta_cond_sq) ** 0.5
        norm_B = float(norm_B_total_sq) ** 0.5
        eta_rel = eta_total / norm_B if norm_B > 0 else float('nan')
        zz_per_stage.append({
            "stage": nStage, "eta_total": eta_total,
            "eta_cond": eta_cond, "norm_B": norm_B,
            "eta_rel": eta_rel,
        })
        print(f"  [ZZ] η_total = {eta_total:.4e}, η_cond = {eta_cond:.4e},"
              f" |B| = {norm_B:.4e}, η_rel = {eta_rel:.4%}",
              flush=True)

        # [DIAG] Galerkin: ∫ν|curl A_phi|² (full mesh) = ∫A_phi · J (cond) ?
        if nStage <= 1:
            E_curl = float(Integrate(
                nu_cf * curl(gfAphi) * curl(gfAphi) * dx(
                    bonus_intorder=BONUS_INT), mesh))
            E_pair = float(Integrate(
                gfAphi * J_n_cf * dx("conductor", bonus_intorder=BONUS_INT), mesh))
            ratio = E_curl / E_pair if abs(E_pair) > 1e-30 else float("nan")
            print(f"   [DIAG-{nStage}] ∫ν|curl A|² = {E_curl:.6e},"
                  f" ∫A·J on cond = {E_pair:.6e}, ratio = {ratio:.6f}"
                  f" (should = 1 if K_A solve is consistent)",
                  flush=True)
            # Also: check σ-divergence of A_phi on cond
            f_div_test = LinearForm(fes_phi)
            f_div_test += sigma_cf * gfAphi * grad(p_te) * dx(
                "conductor", bonus_intorder=BONUS_INT)
            with TaskManager():
                f_div_test.Assemble()
            div_norm = float(np.sqrt(np.sum(np.array(f_div_test.vec)**2)))
            print(f"   [DIAG-{nStage}] |∫σA·grad(test)| = {div_norm:.3e}"
                  f" (should = 0 if σA is div-free on cond)",
                  flush=True)

        # Accumulate: A_acc_n = A_acc_{n-1} + R_n × A_phi_n
        Apot_acc.vec.data += R_n * gfAphi.vec

        # === SNAPSHOT (CRITICAL FIX 2026-05-10) ============================
        # Apot_acc gets += updated next stage. CFs in J_n_cf must capture a
        # FROZEN copy, otherwise re-evaluating J_n_cf later sees the updated
        # Apot_acc — this is the GridFunction overwrite trap that silently
        # corrupted stage 1+ in all prior 3D Kameari attempts.
        snap_acc = GridFunction(fesA, name=f"snap_acc_{nStage}")
        snap_acc.vec.data = Apot_acc.vec  # explicit copy

        # === H-H projection of snap_acc on conductor =======================
        # Solve: ∇·(σ ∇φ) = ∇·(σ A_acc) on cond, Neumann BC
        # → (A_acc - ∇φ) is σ-divergence-free in conductor
        # In axisym this is automatic (A_θ only); in 3D it must be enforced
        # at every stage to keep J_{n+1} = J_n - σ A_acc/L_n div-free.
        f_phi_acc = LinearForm(fes_phi)
        f_phi_acc += sigma_cf * snap_acc * grad(p_te) * dx(
            "conductor", bonus_intorder=BONUS_INT)
        with TaskManager():
            f_phi_acc.Assemble()
        if pinned_phi >= 0:
            f_phi_acc.vec[pinned_phi] = 0.0
        gf_phi_acc = GridFunction(fes_phi, name=f"phi_acc_{nStage}")
        with TaskManager():
            gf_phi_acc.vec.data = inv_phi * f_phi_acc.vec
        # H-H projected accumulator (CoefficientFunction form, FROZEN refs)
        A_acc_proj_cf = snap_acc - grad(gf_phi_acc)

        # COMSOL formula:
        #   τ_n = ∫_cond J_n · A_acc_n dV  (use H-H projected A_acc)
        #   L_n = R_n × τ_n
        AJ_n = float(Integrate(
            A_acc_proj_cf * J_n_cf * dx(
                "conductor", bonus_intorder=BONUS_INT), mesh))
        tau_n_us = 1e6 * AJ_n
        L_n = R_n * AJ_n
        Ln.append(L_n)
        print(f"  ∫J·A_acc(proj) = {AJ_n:.6e}  →  τ_{nStage} = {tau_n_us:.4f} μs",
              flush=True)
        print(f"  L_{nStage} = R_n × τ_n = {L_n:.6e}", flush=True)

        # Update J: J_{n+1} = J_n - σ × A_acc_proj / L_n (div-free preserved)
        J_n_cf = J_n_cf - sigma_cf * A_acc_proj_cf / L_n

        diag.append({
            "stage": nStage,
            "G_n": G_n,
            "R_n": R_n,
            "AJ_integral": AJ_n,
            "tau_n_us": tau_n_us,
            "L_n": L_n,
            "stage_time_s": time.time() - t_stage,
        })
        print(f"  ({time.time()-t_stage:.1f} s)", flush=True)

    # Summary ------------------------------------------------------------
    print("\n" + "=" * 78, flush=True)
    print(" CLN TEAM 28 + Kelvin results", flush=True)
    print("=" * 78, flush=True)
    print(f"  {'stage':>5} {'R_n':>13} {'L_n':>13}"
          f" {'τ_n μs':>13}", flush=True)
    for s in diag:
        print(f"  {s['stage']:>5} {s['R_n']:>13.6e}"
              f" {s['L_n']:>13.6e} {s['tau_n_us']:>13.4f}",
              flush=True)

    # JSON ---------------------------------------------------------------
    out = {
        "method": "CLN TEAM 28 (COMSOL/Sugahara) + Kelvin",
        "geometry": args.geom,
        "params": {
            "sigma_S_per_m": sigma_Cu, "B0_T": B0,
            "ORDER": ORDER_use, "H_COND_mm": H_COND*1000,
            "H_AIR_mm": H_AIR*1000, "R_K_mm": R_K*1000,
            "OFFSET_mm": [v*1000 for v in OFFSET],
            "N_STAGES": n_stages, "BONUS_INT": BONUS_INT,
            "ne": mesh.ne, "fesA_ndof": fesA.ndof,
        },
        "Rn": Rn, "Ln": Ln, "stages": diag,
        "zz_per_stage": zz_per_stage,
    }
    if args.out_tag:
        tag = args.out_tag
    else:
        tag = f"{args.geom}_team28_kelvin"
    out_path = (Path(__file__).parent
                 / f"2026-05-10-cln_team28_{tag}_results.json")
    out_path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"\nResults: {out_path}", flush=True)


if __name__ == "__main__":
    main()
