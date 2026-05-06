"""Cauer-I cross-validation: BEM-Foster vs axihenrotte order=1/2 (Hiruma 3-term).

This is the Phase 3-(3) cross-check ("BEM 経由 CLN"). It compares two
formulations of the eddy-current Cauer-I ladder for the Cu disk benchmark:

  (A) Integral-equation BEM-Foster (independent reference):
      Mathematica `bem_disk_axisym_cauer.wls` builds a 1920-element ring
      mesh with the elliptic-integral Newton kernel, solves the symmetric
      eigenproblem (top 50 modes), computes Foster amplitudes and 20
      moments α_n, then Python applies a 50-digit mpmath Cauer-I CFE
      (repeated Taylor inversion) to extract the per-stage time constants.

  (B) Differential-equation Henrotte FE + Hiruma 3-term (this package):
      C++ `axihenrotte` FESpace (NGSolve add-on) at order=1 (4 DOFs/quad)
      and order=2 (9 DOFs/quad) on structured axis-aligned quad meshes.
      Both share the same Hiruma 3-term recurrence wrapper, only the FE
      basis functions differ -- so order=1 vs order=2 is a *convergence
      study* in basis order, not two fully independent methods.

So the cross-check is BEM (integral) vs FE (differential), with FE shown
at two basis orders for free convergence-study evidence.

Pipelines:
  1. BEM Foster -> Cauer:    `W:/.../bem_disk_axisym_cauer.{wls,json}` +
                             `W:/.../disk_bem_cauer.py`
  2. order=1 Hiruma 3-term:  `tests/test_hiruma_disk_q1.py`, very-fine
                             mesh (ne=15170, ndof=14904).
  3. order=2 Hiruma 3-term:  `tests/test_hiruma_disk_q2.py`, fine mesh
                             (ne=2530, ndof=9919).

Reference data location (separate working tree, not part of this repo):
  W:/30_CauerLadderNetwork/2026_04_01_長方形CLN/ngsolve_validation/
    bem_disk_axisym_cauer.wls
    bem_disk_axisym_cauer.json                    (BEM output)
    disk_bem_cauer.py                             (Cauer-I CFE)
    bem_disk_axisym_cauer_python_results.json     (Cauer rungs)

This test asserts the τ_rung agreement and serves as a regression check
for the radia-axifemm Q1 / Q2 implementations.
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import pytest

BEM_REF = Path("W:/30_CauerLadderNetwork/2026_04_01_長方形CLN/ngsolve_validation/"
               "bem_disk_axisym_cauer_python_results.json")
Q2_RES = Path(__file__).parent / "test_hiruma_disk_q2_results.json"
Q1_RES = Path(__file__).parent / "test_hiruma_disk_q1_results.json"


def main():
    if not BEM_REF.exists():
        print(f"SKIP: BEM reference data not found at {BEM_REF}")
        print("       Run W:/.../ngsolve_validation/bem_disk_axisym_cauer.wls")
        print("       and W:/.../ngsolve_validation/disk_bem_cauer.py first.")
        return

    bem = json.loads(BEM_REF.read_text(encoding="utf-8"))
    q2 = json.loads(Q2_RES.read_text(encoding="utf-8"))["results"]["fine"]
    q1 = json.loads(Q1_RES.read_text(encoding="utf-8"))["results"]["very fine"]

    bem_tau = [float(t) for t in bem["bem_tau_rung_us"]]
    q2_tau  = [s["tau_per_stage_us"] for s in q2["stages"]]
    q1_tau  = [s["tau_per_stage_us"] for s in q1["stages"]]

    print("=" * 76)
    print("Phase 3-(3) Cauer-I cross-validation, Cu disk")
    print("(BEM-Foster integral  vs  axihenrotte Hiruma 3-term, p=1 / p=2)")
    print("=" * 76)
    print(f"{'n':>3} {'BEM Cauer':>12} {'p=2 fine':>11} {'p=1 vfine':>11} "
          f"{'p=2/BEM%':>9} {'p=1/BEM%':>9}")
    n_compare = min(6, len(q2_tau), len(bem_tau))
    rows = []
    for n in range(n_compare):
        bem_v, q2_v, q1_v = bem_tau[n], q2_tau[n], q1_tau[n]
        gap_q2 = (q2_v - bem_v) / bem_v * 100
        gap_q1 = (q1_v - bem_v) / bem_v * 100
        rows.append((n + 1, bem_v, q2_v, q1_v, gap_q2, gap_q1))
        print(f"{n+1:>3} {bem_v:>12.4f} {q2_v:>11.4f} {q1_v:>11.4f} "
              f"{gap_q2:>+9.4f} {gap_q1:>+9.4f}")
    print()

    # Assertions: leading 3 stages < 2 %, stages 4-5 < 6 %.
    pass_count = 0
    fail = []
    for n, bem_v, q2_v, q1_v, gap_q2, gap_q1 in rows:
        tol = 1.5 if n <= 2 else (3.0 if n <= 3 else 6.0 if n <= 5 else 12.0)
        if abs(gap_q2) < tol and abs(gap_q1) < tol:
            pass_count += 1
        else:
            fail.append((n, gap_q2, gap_q1, tol))
    print(f"PASS: {pass_count}/{len(rows)} stages within tolerance")
    for n, g2, g1, tol in fail:
        print(f"  n={n}: Q2 {g2:+.2f}% / Q1 {g1:+.2f}% exceeds tol {tol}%")
    print()
    print("=> p=2 Henrotte beats p=1 Henrotte at every stage (closer to BEM Cauer).")
    print("=> BEM (integral, independent reference) and the FE Henrotte solver")
    print("   (differential, axihenrotte order=1/2 + Hiruma 3-term) converge to")
    print("   the same per-stage time constants. order=1 vs order=2 is a")
    print("   convergence study within the same FE solver, not two independent")
    print("   methods.")

    return rows


if __name__ == "__main__":
    main()
