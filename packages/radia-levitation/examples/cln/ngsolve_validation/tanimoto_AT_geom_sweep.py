"""tanimoto_AT_geom_sweep.py -- Tanimoto A-T alternating Poisson CLN
extraction for cylinder / square_prism (A1) / cuboid_521.

Algorithm (adapted from W:/00_CAE/NGSolve/谷本/修論/CLN_AT.ipynb):
  fesA = HCurl(mesh, ORDER, nograds=True, dirichlet="conductorBND") [in conductor only]
  fesT = HCurl(mesh, ORDER, nograds=True)                            [in conductor only]
  Stage 0 (T equation):
    a_T = ∫(1/σ) curl(T)·curl(W) dx
    f_T = SOURCE  (depends on physics — see below)
    R_0 = 1/∫σ E^2 dx  with E = J/σ, J = curl(T)
  Stage n ≥ 1:
    A equation:  a_A = ∫(1/μ) curl(A)·curl(N) dx,  f_A = ∫N·(R_{2n}·J) dx
                  → gfA, accumulate B = sum curl(gfA), L_{2n+1} = ∫B²/μ dx
    T equation:  a_T = ∫(1/σ) curl(T)·curl(W) dx,  f_T = ∫W·(-B/L_{2n+1}) dx
                  → gfT, accumulate J = sum curl(gfT), R_{2n+2} = 1/∫σ E² dx

Key structural property: J = curl(T), B = curl(A) — curl annihilates
HCurl gauge subspace (∇H1) on each propagation. NO M = ∫σuv matrix is
ever formed, so no K^-1 M step that triggers gauge collapse on rectangular.

Source choice (this script): for uniform applied B_z step at t=0+,
the impulse-equivalent E is E_s = -A_s where A_s = (-y, x, 0)·B0/2 is
the source vector potential giving curl(A_s) = B0 ẑ uniform.
Tanimoto's boundary integral becomes:
  f_T = -∫_∂cond (E_s × W·trace)·n ds   (with E_s = A_s used as drive)
We test this for cylinder first; if R_0 matches the v5/BEM reference
(R_0 ≈ 1.05e7 with our convention), we trust the source and apply to A1
and cuboid_521.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from math import pi
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace", line_buffering=True)

from netgen.occ import (
    Cylinder, Box, Pnt, X, Y, Z, OCCGeometry,
)
from ngsolve import (
    Mesh, HCurl, H1, BilinearForm, LinearForm, GridFunction,
    CoefficientFunction, curl, grad, dx, ds, x, y, z, specialcf,
    Cross, Integrate, TaskManager, ngsglobals, InnerProduct,
)


# === Physical parameters ====================================================
mu0 = 4 * pi * 1e-7
sigma_Cu = 5.8e7
B0 = 1.0

GEOM_TABLE = {
    "cylinder":   {"kind": "cyl",  "R": 10.0,  "t":  2.0},
    "square":     {"kind": "rect", "Lx": 17.72, "Ly": 17.72, "Lz": 2.0},
    "cuboid_521": {"kind": "rect", "Lx":  5.0,  "Ly":  2.0,  "Lz": 1.0},
}

ORDER_DEFAULT = 2
H_COND_DEFAULT_MM = 1.0
N_STAGES_DEFAULT = 4
BONUS_INT = 8


def build_geo(geom_name, h_cond):
    p = GEOM_TABLE[geom_name]
    if p["kind"] == "cyl":
        R = p["R"] * 1e-3
        t = p["t"] * 1e-3
        cond = Cylinder(Pnt(0, 0, -t/2), Z, r=R, h=t)
    else:
        Lx = p["Lx"] * 1e-3
        Ly = p["Ly"] * 1e-3
        Lz = p["Lz"] * 1e-3
        cond = Box(Pnt(-Lx/2, -Ly/2, -Lz/2),
                   Pnt( Lx/2,  Ly/2,  Lz/2))
    cond.faces.name = "conductorBND"
    cond.mat("conductor")
    cond.maxh = h_cond
    return OCCGeometry(cond)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--geom", choices=list(GEOM_TABLE.keys()),
                         default="cylinder")
    parser.add_argument("--order", type=int, default=ORDER_DEFAULT)
    parser.add_argument("--stages", type=int, default=N_STAGES_DEFAULT)
    parser.add_argument("--h-cond-mm", type=float, default=H_COND_DEFAULT_MM)
    parser.add_argument("--out-tag", type=str, default=None)
    args = parser.parse_args()

    ORDER = args.order
    n_stages = args.stages
    geom_name = args.geom
    H_COND = args.h_cond_mm * 1e-3

    p = GEOM_TABLE[geom_name]

    ngsglobals.msg_level = 0
    print("=" * 78, flush=True)
    print(f" Tanimoto A-T alternating Poisson CLN (geom={geom_name})",
          flush=True)
    print("=" * 78, flush=True)
    if p["kind"] == "cyl":
        V_mm3 = pi * p["R"]**2 * p["t"]
        print(f"  Cylinder R={p['R']}mm, t={p['t']}mm, V={V_mm3:.2f} mm^3",
              flush=True)
    else:
        V_mm3 = p["Lx"] * p["Ly"] * p["Lz"]
        print(f"  Rect prism {p['Lx']}x{p['Ly']}x{p['Lz']}mm,"
              f" V={V_mm3:.2f} mm^3", flush=True)
    print(f"  sigma={sigma_Cu:.2e} S/m, B0={B0} T,"
          f" mu0={mu0:.4e} H/m", flush=True)
    print(f"  ORDER={ORDER}, H_COND={H_COND*1000}mm, N_STAGES={n_stages},"
          f" bonus_intorder={BONUS_INT}",
          flush=True)
    print()

    # Build mesh ---------------------------------------------------------
    print(" Building geometry + mesh ...", flush=True)
    geo = build_geo(geom_name, H_COND)
    t0 = time.time()
    with TaskManager():
        ngmesh = geo.GenerateMesh(maxh=H_COND, grading=0.5)
    mesh = Mesh(ngmesh)
    mesh.Curve(ORDER)
    print(f"  GenerateMesh + Curve(p={ORDER}): {time.time()-t0:.1f} s",
          flush=True)
    print(f"  ne={mesh.ne}, nv={mesh.nv}, nedge={mesh.nedge}", flush=True)
    print(f"  materials={mesh.GetMaterials()},"
          f" boundaries={mesh.GetBoundaries()}",
          flush=True)

    # Spaces -------------------------------------------------------------
    fesA = HCurl(mesh, order=ORDER, nograds=True, dirichlet="conductorBND")
    fesT = HCurl(mesh, order=ORDER, nograds=True)
    print(f"  fesA ndof = {fesA.ndof} (active {sum(fesA.FreeDofs())})",
          flush=True)
    print(f"  fesT ndof = {fesT.ndof} (active {sum(fesT.FreeDofs())})",
          flush=True)

    A_tr, A_te = fesA.TnT()
    T_tr, T_te = fesT.TnT()

    # Helmholtz-Hodge projection of A_s to enforce A_s'·n = 0 on ∂cond
    # (so that div(σA_s') = 0 holds across boundary, source ⊥ kernel of curl·curl)
    # Solve: ∫∇φ·∇q dx = ∫A_s·∇q dx ∀q ∈ H1(cond), pin one DOF
    A_s_raw = CoefficientFunction((-y, x, 0)) * (B0 / 2.0)
    n_cf = specialcf.normal(mesh.dim)

    print(" Helmholtz-Hodge projection of A_s ...", flush=True)
    fes_phi_proj = H1(mesh, order=ORDER + 1)
    phi_tr, phi_te = fes_phi_proj.TnT()
    a_phi_proj = BilinearForm(fes_phi_proj)
    a_phi_proj += grad(phi_tr) * grad(phi_te) * dx(
        bonus_intorder=BONUS_INT)
    f_phi_proj = LinearForm(fes_phi_proj)
    f_phi_proj += A_s_raw * grad(phi_te) * dx(bonus_intorder=BONUS_INT)
    with TaskManager():
        a_phi_proj.Assemble()
        f_phi_proj.Assemble()
        fd_phi = fes_phi_proj.FreeDofs()
        # Pin first free DOF
        pinned_phi = -1
        for ip in range(fes_phi_proj.ndof):
            if fd_phi[ip]:
                fd_phi[ip] = False
                pinned_phi = ip
                break
        if pinned_phi >= 0:
            f_phi_proj.vec[pinned_phi] = 0.0
        try:
            inv_phi_proj = a_phi_proj.mat.Inverse(fd_phi, inverse="pardiso")
        except Exception:
            inv_phi_proj = a_phi_proj.mat.Inverse(
                fd_phi, inverse="sparsecholesky")
        gf_phi_proj = GridFunction(fes_phi_proj)
        gf_phi_proj.vec.data = inv_phi_proj * f_phi_proj.vec
    # Verify projection: ∫A_s'·∇q should be 0 for any q
    A_s_cf = A_s_raw - grad(gf_phi_proj)
    leak_norm = float(Integrate((A_s_cf * grad(gf_phi_proj))
                                  * dx(bonus_intorder=BONUS_INT), mesh))
    raw_norm = float(Integrate((A_s_raw * A_s_raw)
                                 * dx(bonus_intorder=BONUS_INT), mesh))
    print(f"  ||A_s||² = {raw_norm:.4e}, ⟨A_s', ∇φ⟩ = {leak_norm:.4e}",
          flush=True)
    print(f"  pinned vertex {pinned_phi}", flush=True)

    # Stage 0: T equation -----------------------------------------------
    print("\n Stage 0: T equation (R_0)", flush=True)
    a_T = BilinearForm(fesT, check_unused=False)
    a_T += (1.0 / sigma_Cu) * curl(T_tr) * curl(T_te) * dx(
        bonus_intorder=BONUS_INT)
    f_T = LinearForm(fesT)
    f_T += sigma_Cu * A_s_cf * T_te * dx(bonus_intorder=BONUS_INT)

    t0 = time.time()
    with TaskManager():
        a_T.Assemble()
        f_T.Assemble()
    print(f"  T form assemble: {time.time()-t0:.1f} s", flush=True)

    gfT = GridFunction(fesT, name="T_0")
    t0 = time.time()
    with TaskManager():
        try:
            inv_T = a_T.mat.Inverse(fesT.FreeDofs(), inverse="pardiso")
            gfT.vec.data = inv_T * f_T.vec
            T_solver = "pardiso"
        except Exception as e:
            print(f"  Pardiso failed ({e}), trying sparsecholesky",
                  flush=True)
            inv_T = a_T.mat.Inverse(fesT.FreeDofs(),
                                     inverse="sparsecholesky")
            gfT.vec.data = inv_T * f_T.vec
            T_solver = "sparsecholesky"
    print(f"  T solve ({T_solver}): {time.time()-t0:.1f} s", flush=True)

    # Compute R_0 from accumulated J --------------------------------------
    J_acc = curl(gfT)
    R_0 = 1.0 / float(Integrate(
        (1.0/sigma_Cu) * J_acc * J_acc * dx(bonus_intorder=BONUS_INT), mesh))
    print(f"  R_0 = {R_0:.6e}", flush=True)

    Rn = [R_0]
    Ln = []
    diag = []

    # Stages 1..N --------------------------------------------------------
    Apot_acc = None  # accumulated A
    B_acc = None     # accumulated B
    R_curr = R_0

    for nStage in range(n_stages):
        t_stage = time.time()
        print(f"\n Stage {nStage+1}/{n_stages}", flush=True)

        # A equation -----------------------------------------------------
        a_A = BilinearForm(fesA, check_unused=False)
        a_A += (1.0 / mu0) * curl(A_tr) * curl(A_te) * dx(
            bonus_intorder=BONUS_INT)
        f_A = LinearForm(fesA)
        # source = R_curr * J_acc (volume integral)
        f_A += A_te * (R_curr * J_acc) * dx(bonus_intorder=BONUS_INT)

        with TaskManager():
            a_A.Assemble()
            f_A.Assemble()

        gfA = GridFunction(fesA, name=f"A_{nStage}")
        with TaskManager():
            try:
                inv_A = a_A.mat.Inverse(fesA.FreeDofs(), inverse="pardiso")
            except Exception:
                inv_A = a_A.mat.Inverse(fesA.FreeDofs(),
                                         inverse="sparsecholesky")
            gfA.vec.data = inv_A * f_A.vec

        # Accumulate A and B
        if Apot_acc is None:
            Apot_acc = GridFunction(fesA)
            Apot_acc.vec.data = gfA.vec
            B_acc = curl(gfA)  # CoefficientFunction
        else:
            Apot_acc.vec.data += gfA.vec
            B_acc = B_acc + curl(gfA)

        L_curr = float(Integrate(B_acc * B_acc / mu0 * dx(
            bonus_intorder=BONUS_INT), mesh))
        Ln.append(L_curr)
        print(f"  L_{2*nStage+1} = {L_curr:.6e}", flush=True)

        # T equation -----------------------------------------------------
        a_T2 = BilinearForm(fesT, check_unused=False)
        a_T2 += (1.0 / sigma_Cu) * curl(T_tr) * curl(T_te) * dx(
            bonus_intorder=BONUS_INT)
        f_T2 = LinearForm(fesT)
        # source = -B_acc / L_curr  (volume integral)
        f_T2 += T_te * (-B_acc / L_curr) * dx(bonus_intorder=BONUS_INT)

        with TaskManager():
            a_T2.Assemble()
            f_T2.Assemble()

        gfT_new = GridFunction(fesT, name=f"T_{nStage+1}")
        with TaskManager():
            try:
                inv_T2 = a_T2.mat.Inverse(fesT.FreeDofs(),
                                            inverse="pardiso")
            except Exception:
                inv_T2 = a_T2.mat.Inverse(fesT.FreeDofs(),
                                            inverse="sparsecholesky")
            gfT_new.vec.data = inv_T2 * f_T2.vec

        # Accumulate J
        J_acc = J_acc + curl(gfT_new)

        R_curr = 1.0 / float(Integrate(
            (1.0/sigma_Cu) * J_acc * J_acc * dx(bonus_intorder=BONUS_INT),
            mesh))
        Rn.append(R_curr)
        print(f"  R_{2*(nStage+1)} = {R_curr:.6e}", flush=True)

        tau_pair_us = 1e6 * L_curr / Rn[nStage]  # tau = L_{2n+1} / R_{2n}
        diag.append({
            "stage": nStage + 1,
            "R_2n_idx": 2 * nStage,
            "R_2n": Rn[nStage],
            "L_2nplus1_idx": 2 * nStage + 1,
            "L_2nplus1": L_curr,
            "tau_pair_us": tau_pair_us,
            "R_2nplus2": R_curr,
            "stage_time_s": time.time() - t_stage,
        })
        print(f"  τ_pair[{nStage}] = L_{2*nStage+1} / R_{2*nStage}"
              f" = {tau_pair_us:.4f} μs"
              f"  ({time.time()-t_stage:.1f} s)",
              flush=True)

    # Summary table ------------------------------------------------------
    print("\n" + "=" * 78, flush=True)
    print(" Tanimoto A-T results", flush=True)
    print("=" * 78, flush=True)
    print(f"  {'stage':>5} {'R_2n':>13} {'L_2n+1':>13}"
          f" {'τ_pair us':>13}", flush=True)
    for s in diag:
        print(f"  {s['stage']:>5} {s['R_2n']:>13.6e}"
              f" {s['L_2nplus1']:>13.6e} {s['tau_pair_us']:>13.4f}",
              flush=True)

    # Save JSON ----------------------------------------------------------
    out = {
        "method": "Tanimoto A-T alternating Poisson CLN",
        "geometry": geom_name,
        "geom_params": p,
        "params": {
            "sigma_S_per_m": sigma_Cu, "B0_T": B0,
            "ORDER": ORDER, "H_COND_mm": H_COND * 1000,
            "N_STAGES": n_stages, "BONUS_INT": BONUS_INT,
            "ne": mesh.ne, "fesA_ndof": fesA.ndof, "fesT_ndof": fesT.ndof,
        },
        "Rn": Rn,
        "Ln": Ln,
        "stages": diag,
    }
    if args.out_tag:
        tag = args.out_tag
    else:
        tag = f"{geom_name}_p{ORDER}"
    out_path = (Path(__file__).parent
                 / f"2026-05-09-tanimoto_AT_{tag}_results.json")
    out_path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"\nResults: {out_path}", flush=True)


if __name__ == "__main__":
    main()
