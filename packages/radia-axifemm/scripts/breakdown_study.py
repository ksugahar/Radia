"""Failure-mode comparison: FE-Hiruma 3-term vs BEM-Foster classical Cauer.

Paper-(a) supporting evidence. Hypothesis:
  - FE-Hiruma 3-term (Lanczos-style recursion without reorthogonalization)
    accumulates orthogonality drift exponentially and **breaks
    catastrophically** at some stage k* (negative R, sign-flipped tau, NaN).
  - BEM-Foster classical Cauer extraction degrades **gracefully**: rung
    accuracy decreases monotonically, sign flips appear at a higher k, but
    no NaN/instability cliff.

This script pushes both methods to N_stages = 15 and dumps the per-stage
trace so the contrast is plotable for the manuscript.

Output: scripts/breakdown_study_results.json
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
TESTS = HERE.parent / "tests"
sys.path.insert(0, str(TESTS))

from test_hiruma_disk_q1 import (  # type: ignore
    make_structured_disk_quad_mesh, to_scipy_csr, hiruma_3term,
    R_DISK, T_DISK, SIGMA_CU, MU0, B0,
)

from math import pi

import scipy.sparse.linalg as spla

from ngsolve import (
    Mesh, BilinearForm, LinearForm, CoefficientFunction, TaskManager,
    x, dx, ngsglobals,
)
from radia.radia_axifemm import (
    H1Henrotte, AxiHenrotteStiffnessBFI, AxiHenrotteSigmaMassBFI,
)


def run_fe(order: int, NR_d: int, Nz_d: int, NR_a: int, Nz_a: int,
           R_a: float, Z_a: float, N_stages: int, label: str):
    print(f"\n--- FE Hiruma order={order}  {label}  N_stages={N_stages} ---")
    ngsglobals.msg_level = 0
    mesh = make_structured_disk_quad_mesh(NR_d, Nz_d, NR_a, Nz_a, R_a, Z_a)
    fes = H1Henrotte(mesh, order=order, dirichlet="axis|right|top|bot")
    n_free = sum(1 for f in fes.FreeDofs() if f)
    print(f"  ne={mesh.ne}  ndof={fes.ndof}  free={n_free}")

    mu_cf = CoefficientFunction(MU0)
    sigma_cf = mesh.MaterialCF({"conductor": SIGMA_CU}, default=0.0)
    A_imposed = B0 * x / 2

    a = BilinearForm(fes, symmetric=True)
    a += AxiHenrotteStiffnessBFI(mu_cf)
    with TaskManager():
        a.Assemble()

    m = BilinearForm(fes, symmetric=True)
    m += AxiHenrotteSigmaMassBFI(sigma_cf)
    with TaskManager():
        m.Assemble()

    b_form = LinearForm(fes)
    v = fes.TestFunction()
    b_form += sigma_cf * A_imposed * v * 2 * pi * x * dx
    with TaskManager():
        b_form.Assemble()
    b_vec = np.array(b_form.vec)

    K_csr = to_scipy_csr(a.mat, fes.ndof)
    M_csr = to_scipy_csr(m.mat, fes.ndof)
    free = np.array([i for i in range(fes.ndof) if fes.FreeDofs()[i]],
                    dtype=int)
    K_red = K_csr[free[:, None], free[None, :]]
    M_red = M_csr[free[:, None], free[None, :]]
    b_red = b_vec[free]

    res = hiruma_3term(K_red, M_red, b_red, N_stages, label=label)
    res["mesh"] = {"ne": mesh.ne, "ndof": fes.ndof, "free": int(len(free)),
                   "order": order}
    return res


def main():
    N_STAGES = 15

    print("=" * 78)
    print(f"Breakdown study: FE-Hiruma vs BEM-Foster, N_stages = {N_STAGES}")
    print("=" * 78)

    fe_q1 = run_fe(order=1, NR_d=80, Nz_d=16, NR_a=20, Nz_a=20,
                   R_a=500e-3, Z_a=500e-3, N_stages=N_STAGES,
                   label="Q1 very fine")
    fe_q2 = run_fe(order=2, NR_d=40, Nz_d=16, NR_a=15, Nz_a=15,
                   R_a=500e-3, Z_a=500e-3, N_stages=N_STAGES,
                   label="Q2 fine")

    # Prefer the high-moment BEM Cauer extraction (60 moments, 50-digit
    # mpmath) over the legacy 20-moment extraction; the high-moment version
    # extends the reliable Cauer ladder to k ~ 17 vs the legacy k ~ 5.
    bem_hp = HERE / "bem_cauer_high_moments_results.json"
    bem_lp = Path("W:/30_CauerLadderNetwork/2026_04_01_長方形CLN/"
                  "ngsolve_validation/bem_disk_axisym_cauer_python_results.json")
    if bem_hp.exists():
        bem = json.loads(bem_hp.read_text(encoding="utf-8"))
        bem_tau = bem["tau_pair_us_float"]
        print(f"\nLoaded BEM (high-moment 60): {bem_hp.name}")
        print(f"  BEM stages available: {len(bem_tau)}, "
              f"first sign flip at k = {bem['breakdown_k']['sign_flip']}")
    elif bem_lp.exists():
        bem = json.loads(bem_lp.read_text(encoding="utf-8"))
        bem_tau = [float(t) for t in bem.get("tau_pair_us",
                                             bem.get("bem_tau_rung_us"))]
        print(f"\nLoaded BEM (legacy 20-moment): {bem_lp.name}")
        print(f"  BEM stages available: {len(bem_tau)}")
    else:
        print(f"\nWARN: BEM reference not found")
        bem_tau = []

    # ---- Comparative breakdown table ----
    print("\n" + "=" * 78)
    print("Per-stage tau_pair (us) and orthogonality drift")
    print("=" * 78)
    print(f"{'k':>3}  {'BEM Cauer':>12}  "
          f"{'Q1 tau':>11}  {'Q1 drift':>10}  "
          f"{'Q2 tau':>11}  {'Q2 drift':>10}")
    rows = []
    for k in range(N_STAGES):
        b = bem_tau[k] if k < len(bem_tau) else None
        q1 = fe_q1["stages"][k]
        q2 = fe_q2["stages"][k]
        bs = f"{b:>+12.4f}" if b is not None else f"{'--':>12}"
        print(f"{k:>3}  {bs}  "
              f"{q1['tau_pair_us']:>+11.4f}  {q1['orth_drift']:>10.2e}  "
              f"{q2['tau_pair_us']:>+11.4f}  {q2['orth_drift']:>10.2e}")
        rows.append({
            "k": k,
            "bem_tau_us": b,
            "q1_tau_us": q1["tau_pair_us"],
            "q1_drift": q1["orth_drift"],
            "q1_R_2k": q1["R_2k"],
            "q1_L_2k_plus_1": q1["L_2k_plus_1"],
            "q2_tau_us": q2["tau_pair_us"],
            "q2_drift": q2["orth_drift"],
            "q2_R_2k": q2["R_2k"],
            "q2_L_2k_plus_1": q2["L_2k_plus_1"],
        })

    # ---- Breakdown detection ----
    def first_breakdown(stages, tau_key="tau_pair_us"):
        for s in stages:
            t = s[tau_key]
            if t is None or t < 0 or not np.isfinite(t):
                return s["k"]
        return None

    def first_neg_bem(taus):
        for k, t in enumerate(taus):
            if t < 0 or not np.isfinite(t):
                return k
        return None

    print("\n" + "=" * 78)
    print("Breakdown detection")
    print("=" * 78)
    k_q1 = first_breakdown(fe_q1["stages"])
    k_q2 = first_breakdown(fe_q2["stages"])
    k_bem = first_neg_bem(bem_tau) if bem_tau else None
    print(f"  FE Q1 (Hiruma 3-term):    first sign flip / NaN at k = {k_q1}")
    print(f"  FE Q2 (Hiruma 3-term):    first sign flip / NaN at k = {k_q2}")
    print(f"  BEM-Foster classical Cauer: first sign flip at k = {k_bem}")
    print()
    print("  -> Hiruma drift growth indicates Lanczos-style accuracy loss.")
    print("     A sign flip / NaN before k=BEM indicates catastrophic")
    print("     breakdown of the FE method, while BEM continues smoothly")
    print("     beyond its accuracy cliff.")

    out = {
        "config": {"N_stages": N_STAGES},
        "fe_q1": fe_q1,
        "fe_q2": fe_q2,
        "bem_tau_us": bem_tau,
        "comparison": rows,
        "breakdown_k": {"q1": k_q1, "q2": k_q2, "bem": k_bem},
    }
    out_path = HERE / "breakdown_study_results.json"
    out_path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"\nSaved: {out_path}")


if __name__ == "__main__":
    main()
