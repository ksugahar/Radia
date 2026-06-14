"""Golden lock for radia-streamfunction Volume-3D mode.

Runs calc_streamfunction_volume.py as a subprocess (no panel) on the frozen
hollow-cylinder fixture and asserts the foliated-Clebsch volume stream-function
designer reproduces the validated band:

  - the two-codebase Biot-Savart invariant (my straight-segment law vs Radia
    rad.ObjFlmCur+rad.Fld) agree to < 1e-3  -- the golden cannot pass on a wrong
    Biot-Savart;
  - the equal-current wires reproduce the target axial Bz (field residual band);
  - the helicity gate is ~0 (Clebsch-feasible) and the inverse fit is tight;
  - the standard housekeeping keys are present (ne / ndof / t_*_s / gmsh_file).

Mirrors examples/vim/foliated_solenoid_wires.py (Stage B), promoted to the
panel tier.
"""
import json
import os
import subprocess
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
CALC = os.path.join(HERE, "..", "..", "src", "radia", "panels",
                    "calc_streamfunction_volume.py")
FIXTURE = os.path.join(HERE, "fixtures", "streamfunction_volume_tube.vol")
GEN = os.path.join(HERE, "fixtures", "generate_streamfunction_volume_tube.py")
GOLDEN = os.path.join(HERE, "golden", "streamfunction_volume.json")


def _load_golden():
    with open(GOLDEN, encoding="utf-8") as f:
        return json.load(f)


def _ensure_fixture():
    """Generate the tube .vol on demand (fixtures are NOT tracked in git;
    regenerated from generate_streamfunction_volume_tube.py)."""
    if os.path.isfile(FIXTURE):
        return
    proc = subprocess.run([sys.executable, GEN],
                          capture_output=True, text=True, timeout=120)
    if proc.returncode != 0 or not os.path.isfile(FIXTURE):
        pytest.skip(f"could not generate tube fixture:\n"
                    f"STDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}")


@pytest.mark.slow
def test_streamfunction_volume_golden(tmp_path):
    _ensure_fixture()
    out_json = tmp_path / "result.json"
    out_msh = tmp_path / "result.msh"
    cmd = [sys.executable, CALC,
           "--vol", FIXTURE,
           "--output", str(out_json),
           "--msh-output", str(out_msh)]
    subprocess.run(cmd, check=True, timeout=600)

    got = json.loads(out_json.read_text())
    g = _load_golden()

    # --- the load-bearing invariant: two independent Biot-Savart codebases ---
    assert got["two_codebase"] < 1e-3, (
        f"two-codebase Biot-Savart disagree: {got['two_codebase']:.2e} "
        "(straight-segment vs rad.ObjFlmCur+Fld)")

    # --- helicity gate (radial foliation is Clebsch -> ~0) ---
    assert got["h_rel"] < 5e-2, f"helicity gate failed: H_rel={got['h_rel']:.2e}"

    # --- inverse design fits the target ---
    assert got["fit_res"] < 1e-3, f"ACA+TSVD fit drift: {got['fit_res']:.2e}"

    # --- equal-current wires reproduce the target axial Bz ---
    assert got["field_res"] < 0.06, f"field residual too high: {got['field_res']:.3f}"
    assert abs(got["field_res"] - g["field_res"]) < 0.01, (
        f"field_res drift: got {got['field_res']:.4f}, golden {g['field_res']:.4f}")
    assert got["equal_current_spread"] < 1e-6, (
        f"wires not equal-current: spread={got['equal_current_spread']:.2e}")

    # --- wire count locked to the frozen mesh (small band for contour edges) ---
    assert abs(got["n_wires"] - g["n_wires"]) <= 2, (
        f"n_wires drift: got {got['n_wires']}, golden {g['n_wires']}")
    assert got["n_wires"] >= 3 * 3, "too few wires (n_leaves=3 expects >= 9)"

    # --- mesh / solve invariants ---
    assert got["ne"] == g["ne"], f"mesh ne changed: {got['ne']} vs {g['ne']}"
    assert got["k_aca"] == g["k_aca"], "ACA rank changed (gauge truncation)"

    # --- standard housekeeping keys (Result Output Policy) ---
    for k in ("ne", "ndof", "t_total_s", "n_wires", "field_res",
              "two_codebase", "gmsh_file"):
        assert k in got, f"calc output missing required key '{k}': {sorted(got)}"
    assert os.path.exists(out_msh), "wire .msh overlay not written"
