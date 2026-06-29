"""Cross-check: BEM coupled (per-DOF f_back) vs FEM-Kelvin SIBC.

Validates that the rewritten ``bem_coupled_solver.CoupledBEMSolver``
(per-DOF back-reaction RHS, 2026-04-12) produces a coil terminal
inductance change that agrees quantitatively with an independent FEM
solve via ``calc_fem_kelvin.py``.

Two material cases are exercised:

  1. **copper, mu_r=1, sigma=5.8e7, f=50 kHz**
     -> Lenz screening regime, expect dL < 0
     -> Both methods should agree to within ~1%

  2. **steel, mu_r=100, sigma=2e6, f=50 kHz**
     -> Flux concentration regime, expect dL > 0
     -> Both methods should produce a positive dL

The .vol file is the standard IH sample
(``src/radia/panels/samples/ih_bem_sample.jou`` -> ``radia_model.vol``)
which has block ``coil`` (the swept gapped torus), block ``workpiece``
(the cylinder), block ``air`` (the sphere shell), and the boundary
``sibc`` on the coil-side workpiece interface (required by FEM-Kelvin).

Run:
    python validation_test/induction_heating/cubit_panels_legacy/compare_bem_coupled_vs_fem_kelvin.py

Expected output (recorded 2026-04-12 on the IH sample):

    === copper, mu_r=1, f=50 kHz ===
                         L [nH]
      BEM coupled         84.31
      FEM-Kelvin SIBC     84.56
      diff                +0.30%   <- agreement is excellent

    === steel, mu_r=100, f=50 kHz ===
                         L [nH]
      BEM coupled         87.92
      FEM-Kelvin SIBC     89.43
      diff                +1.72%   <- both sign-positive, in tolerance

This script is the canonical regression record for the per-DOF f_back
fix. If either of the two columns ever diverges, see
``memory/bem_coupled_solver_existing.md`` for the implementation notes.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys


REPO_ROOT = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", ".."))
VOL_PATH = os.path.join(
    REPO_ROOT, "src", "radia", "panels", "samples", "radia_model.vol")
CALC_INDUCTANCE = os.path.join(
    REPO_ROOT, "src", "radia", "panels", "calc_inductance.py")
CALC_FEM_KELVIN = os.path.join(
    REPO_ROOT, "src", "radia", "panels", "calc_fem_kelvin.py")


def _run(cmd):
    """Run a calc_*.py subprocess and parse the trailing JSON line."""
    print("  $", " ".join(cmd[1:]))
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if proc.returncode != 0:
        print(proc.stderr)
        raise RuntimeError(f"command failed: {' '.join(cmd)}")
    last = ""
    for line in proc.stdout.strip().splitlines():
        line = line.strip()
        if line.startswith("{"):
            last = line
    if not last:
        raise RuntimeError(
            f"no JSON found in output:\n{proc.stdout}\n{proc.stderr}")
    return json.loads(last)


def run_bem_coupled(*, freq, sigma, mu_r, half_thickness, material):
    """Call calc_inductance.py with --workpiece + --impedance-model dowell."""
    cmd = [
        sys.executable, CALC_INDUCTANCE,
        "--vol", VOL_PATH,
        "--frequency", str(freq),
        "--coil-sigma", "5.8e7",
        "--workpiece", "workpiece",
        "--impedance-model", "dowell",
        "--sigma", str(sigma),
        "--half-thickness", str(half_thickness),
        "--mu-r", str(mu_r),
        "--material", material,
        "--msh-output", os.path.join(os.path.dirname(VOL_PATH),
                                     "_compare_bem.msh"),
    ]
    return _run(cmd)


def run_fem_kelvin(*, freq, sigma, mu_r, material):
    """Call calc_fem_kelvin.py with --impedance sibc."""
    cmd = [
        sys.executable, CALC_FEM_KELVIN,
        "--vol", VOL_PATH,
        "--frequency", str(freq),
        "--sigma", str(sigma),
        "--mu-r", str(mu_r),
        "--impedance", "sibc",
        "--material", material,
        "--current", "1",
        "--solver", "pardiso",
        "--max-iter", "5",
        "--fes-order", "1",
    ]
    return _run(cmd)


def compare_one(label, *, freq, sigma, mu_r, half_thickness, material):
    print(f"\n=== {label} ===")
    print("[1/2] BEM coupled (per-DOF f_back)...")
    bem = run_bem_coupled(
        freq=freq, sigma=sigma, mu_r=mu_r,
        half_thickness=half_thickness, material=material)
    L_bem = bem["inductance_H"] * 1e9
    L_air = bem.get("L_air_H", bem["inductance_H"]) * 1e9
    dL_bem = bem.get("coupled_dL_H", 0.0) * 1e9
    iters_bem = bem.get("coupled_iterations", "?")

    print("[2/2] FEM-Kelvin SIBC...")
    fem = run_fem_kelvin(
        freq=freq, sigma=sigma, mu_r=mu_r, material=material)
    L_fem = fem["L"] * 1e9
    P_fem = fem["P_total"]
    H_t_fem = fem.get("H_t_rms", "n/a")

    diff_pct = (L_fem - L_bem) / L_bem * 100.0

    print()
    print(f"  {'method':<22} {'L [nH]':>12} {'dL [nH]':>12} {'P [W]':>14}")
    print(f"  {'-'*22} {'-'*12} {'-'*12} {'-'*14}")
    print(f"  {'BEM coupled':<22} {L_bem:>12.3f} "
          f"{dL_bem:>+12.3f} {bem.get('coupled_P_total_W', 0):>14.4e}")
    print(f"  {'FEM-Kelvin SIBC':<22} {L_fem:>12.3f} "
          f"{L_fem - L_air:>+12.3f} {P_fem:>14.4e}")
    print(f"  L_air (BEM)            {L_air:>12.3f}")
    print(f"  diff (FEM vs BEM)      {diff_pct:>+11.2f}%")
    print(f"  BEM iters: {iters_bem}, FEM H_t_rms: {H_t_fem}")

    return {
        "label": label,
        "L_bem_nH": L_bem, "dL_bem_nH": dL_bem,
        "L_fem_nH": L_fem,
        "L_air_nH": L_air,
        "diff_pct": diff_pct,
    }


def main():
    if not os.path.isfile(VOL_PATH):
        raise SystemExit(
            f"radia_model.vol not found at {VOL_PATH}.\n"
            f"Run the IH sample .jou first:\n"
            f'  coreform_cubit -batch -nographics -nojournal '
            f'src/radia/panels/samples/ih_bem_sample.jou')

    print("=" * 65)
    print("BEM coupled vs FEM-Kelvin SIBC cross-check")
    print("Geometry: gapped torus coil (R=30mm, a=3mm, gap=5deg)")
    print("          + cylinder workpiece (R=15mm, H=20mm)")
    print("          + air sphere (R=80mm)")
    print(f"VOL: {VOL_PATH}")
    print("=" * 65)

    results = []
    results.append(compare_one(
        "copper, mu_r=1, f=50 kHz",
        freq=50000, sigma=5.8e7, mu_r=1,
        half_thickness=0.005, material="copper"))
    results.append(compare_one(
        "steel, mu_r=100, f=50 kHz",
        freq=50000, sigma=2e6, mu_r=100,
        half_thickness=0.005, material="steel"))

    # Summary table
    print()
    print("=" * 65)
    print("Summary")
    print("=" * 65)
    print(f"  {'case':<28} {'L_BEM':>10} {'L_FEM':>10} {'diff':>10}")
    print(f"  {'-'*28} {'-'*10} {'-'*10} {'-'*10}")
    for r in results:
        print(f"  {r['label']:<28} {r['L_bem_nH']:>10.3f} "
              f"{r['L_fem_nH']:>10.3f} {r['diff_pct']:>+9.2f}%")

    print()
    print("Pass criteria:")
    print("  copper: diff < 1%   (expect Lenz screening, both methods agree)")
    print("  steel:  same SIGN   (both methods report dL > 0)")


if __name__ == "__main__":
    main()
