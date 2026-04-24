"""Golden regression tests for calc_accel_magnet.py.

Locks the C-type dipole EM panel end-to-end path:

    em_sample.jou  (ELF CEFC 2020 C-Type geometry, half-z model)
       --> radia_export netgen   (via Cubit + Auto-Kelvin)
          --> em_sample.vol      (with cd3names "GND" from nodeset propagation)

    em_sample_coil.{step|py}
       --> _load_coil_script     (STEP via coil_from_step OR .py via build_coil)
          --> CoilBuilder        (current=1.0 for .step, NI=2000 for .py)

    calc_accel_magnet.py  (--formulation omega, --material linear, --mu-r 1000)
       --> FES:Periodic H1 + GND Dirichlet + Kelvin
       --> PARDISO linear solve
       --> {B_origin, W_mag, L, ndof, ne, converged}

The two test variants (coil-as-STEP and coil-as-Python) lock both the
linearity (B proportional to current) and the absolute numbers for
the given mesh discretization.  A regression in any of the following
will trip the test:

  - ExportNetgenCommand.cpp  (CD3 propagation / nodeset -> BBBND)
  - add_kelvin_cubit.py      (air + kelvin periodic pair)
  - add_periodic_kelvin      (graceful skip on Cubit-written ident)
  - coil_from_step.to_coil_builder  (tuple unpack + planar coil extract)
  - calc_accel_magnet._load_coil_script  (STEP dispatch)
  - calc_accel_magnet solve_accel  (BND+BBBND union, Dirichlet GND)

Run::

    pytest tests/panels/test_em_golden.py -v -m slow

Marked slow because each run is ~12 min (PARDISO on 456k-element mesh).
"""
from __future__ import annotations

import json
import os
import subprocess
import sys

import pytest


HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
CALC = os.path.join(ROOT, "src", "radia", "panels", "calc_accel_magnet.py")
GOLDEN = os.path.join(HERE, "golden", "em_sample_mu1000.json")


def _load_golden() -> dict:
    with open(GOLDEN, "r", encoding="utf-8") as f:
        return json.load(f)


def _resolve(path: str) -> str:
    """Resolve a path as absolute or repo-relative."""
    if os.path.isabs(path):
        return path
    return os.path.join(ROOT, path)


def _require_sample(g: dict) -> tuple[str, str, str]:
    """Resolve sample file paths; skip test with a helpful message if
    any of them is missing.  The .vol is regenerated at deploy time
    (not tracked in git), so a fresh repo clone starts with a skipped
    test until deploy has been run on this machine.
    """
    vol = _resolve(g["sample"]["vol"])
    step = _resolve(g["sample"]["step"])
    coil_py = _resolve(g["sample"]["coil_script"])
    jou = _resolve(g["sample"]["generator_jou"])
    missing = [p for p in (vol, step, coil_py) if not os.path.isfile(p)]
    if missing:
        pytest.skip(
            "EM golden sample files missing: "
            + ", ".join(os.path.relpath(p, ROOT) for p in missing)
            + f"\nRegenerate the .vol via Cubit batch using {jou}"
              f" + auto_add_kelvin_from_current_model();"
              f" see tests/panels/golden/em_sample_mu1000.json"
              f" for the exact command template.")
    return vol, step, coil_py


def _run(coil_path: str, vol: str, phys: dict, timeout_s: int = 1500) -> dict:
    cmd = [sys.executable, CALC,
           "--coil-script", coil_path,
           "--vol", vol,
           "--formulation", phys["formulation"],
           "--material", phys["material"],
           "--mu-r", str(phys["mu_r"]),
           "--fes-order", str(phys["fes_order"]),
           "--max-iter", str(phys["max_iter"]),
           "--tol", str(phys["tol"]),
           "--solver", phys["solver"]]
    proc = subprocess.run(cmd, capture_output=True, text=True,
                          timeout=timeout_s)
    assert proc.returncode == 0, (
        f"calc_accel_magnet failed (rc={proc.returncode}):\n"
        f"STDOUT tail:\n{proc.stdout[-2000:]}\n"
        f"STDERR tail:\n{proc.stderr[-2000:]}")
    lines = [ln for ln in proc.stdout.splitlines()
             if ln.strip().startswith("{")]
    assert lines, f"no JSON in stdout:\n{proc.stdout[-2000:]}"
    return json.loads(lines[-1])


def _assert_close(got: float, expected: float, tol_pct: float, label: str):
    assert expected != 0, f"{label}: expected is zero, cannot compute pct"
    dev_pct = abs(got - expected) / abs(expected) * 100
    assert dev_pct <= tol_pct, (
        f"{label}: got {got:.4e}, expected {expected:.4e} "
        f"({dev_pct:.2f}% dev, tol {tol_pct:.2f}%)")


@pytest.mark.slow
def test_em_sample_coil_step_current_1A():
    """STEP coil input (current=1.0 baked by coil_from_step.to_coil_builder)."""
    g = _load_golden()
    vol, step, _coil_py = _require_sample(g)

    result = _run(step, vol, g["physics"])

    exp = g["expected_linear"]["I_1A_coil_step"]
    tol = g["tolerance"]

    # Exact checks (mesh-topology invariants)
    if tol["ndof_exact"]:
        assert result["ndof"] == exp["ndof"], \
            f"ndof drift: {result['ndof']} vs {exp['ndof']}"
    if tol["ne_exact"]:
        assert result["ne"] == exp["ne"], \
            f"ne drift: {result['ne']} vs {exp['ne']}"

    assert result["converged"] == exp["converged"]

    # Field quantities
    _assert_close(result["B_origin_mag"], exp["B_origin_mag_T"],
                  tol["B_origin_pct"], "B_origin_mag")
    _assert_close(result["W_mag"], exp["W_mag_J"],
                  tol["W_mag_pct"], "W_mag")
    _assert_close(result["L"], exp["L_H"], tol["L_pct"], "L")


@pytest.mark.slow
def test_em_sample_coil_py_linearity_NI_2000():
    """Python coil build_coil() bakes NI=2000; B should scale linearly,
    L and normalised energy should be invariant."""
    g = _load_golden()
    vol, _step, coil_py = _require_sample(g)

    # Skip gracefully if the NI=2000 golden value hasn't been populated
    # yet (initial check-in placeholder).
    exp_2k = g["expected_linear"]["I_2000A_coil_py"]
    if isinstance(exp_2k.get("B_origin_mag_T_approx"), str):
        pytest.skip("NI=2000 golden value not yet populated in em_sample_mu1000.json")

    result = _run(coil_py, vol, g["physics"])
    ref_1A = g["expected_linear"]["I_1A_coil_step"]
    tol = g["tolerance"]

    # Scaling check: B at NI=2000 should be 2000x B at NI=1.
    expected_mag = ref_1A["B_origin_mag_T"] * 2000.0
    _assert_close(result["B_origin_mag"], expected_mag,
                  tol["B_origin_pct"], "B_origin_mag (NI=2000 scaling)")

    # L is geometry-only, independent of current.
    _assert_close(result["L"], ref_1A["L_H"], tol["L_pct"],
                  "L (invariant with I)")

    # W ∝ I², 4e6x from NI=1 to NI=2000.
    expected_W = ref_1A["W_mag_J"] * (2000.0 ** 2)
    _assert_close(result["W_mag"], expected_W,
                  tol["W_mag_pct"], "W_mag (NI=2000 scaling)")

    assert result["converged"]
