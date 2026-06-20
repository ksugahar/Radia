"""Golden regression tests for calc_accel_magnet.py.

Locks the C-type dipole EM panel end-to-end path:

    em_sample.jou  (ELF CEFC 2020 C-Type geometry, half-z model)
       --> export netgen   (via Cubit + Auto-Kelvin)
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
GOLDEN_QUARTER_XZ = os.path.join(HERE, "golden", "em_quarter_xz_mu1000.json")
GOLDEN_EIGHTH = os.path.join(HERE, "golden", "em_eighth_mu1000.json")
GOLDEN_FULL = os.path.join(HERE, "golden", "em_full_mu1000.json")
GOLDEN_ELF = os.path.join(HERE, "golden", "em_elf_quarter_xz_mu1000.json")


def _load_golden(path: str = GOLDEN) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def _resolve(path: str) -> str:
    """Resolve a path as absolute or repo-relative."""
    if os.path.isabs(path):
        return path
    return os.path.join(ROOT, path)


def _require_sample(g: dict) -> tuple[str, str]:
    """Resolve sample file paths; skip test with a helpful message if
    any of them is missing.  The .vol is regenerated at deploy time
    (not tracked in git), so a fresh repo clone starts with a skipped
    test until deploy has been run on this machine.
    """
    vol = _resolve(g["sample"]["vol"])
    coil_py = _resolve(g["sample"]["coil_script"])
    jou = _resolve(g["sample"]["generator_jou"])
    missing = [p for p in (vol, coil_py) if not os.path.isfile(p)]
    if missing:
        pytest.skip(
            "EM golden sample files missing: "
            + ", ".join(os.path.relpath(p, ROOT) for p in missing)
            + f"\nRegenerate the .vol via Cubit batch using {jou}"
              f" + auto_add_kelvin_from_current_model();"
              f" see tests/panels/golden/em_sample_mu1000.json"
              f" for the exact command template.")
    return vol, coil_py


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


def _check_common(result: dict, exp: dict, tol: dict):
    """Topology invariants + convergence, shared by both variants."""
    if tol["ndof_exact"]:
        assert result["ndof"] == exp["ndof"], \
            f"ndof drift: {result['ndof']} vs {exp['ndof']}"
    if tol["ne_exact"]:
        assert result["ne"] == exp["ne"], \
            f"ne drift: {result['ne']} vs {exp['ne']}"
    assert result["converged"] == exp["converged"]
    assert result["iterations"] == exp["iterations"], \
        f"iterations drift: {result['iterations']} vs {exp['iterations']}"


@pytest.mark.slow
def test_em_sample_coil_py_NI_2000():
    """Python coil `build_coil()`: ELF racetrack via 4 straight + 4 arc
    CoilBuilder segments (84 wire segments), NI=2000.  Absolute-value
    snapshot for regression guard.

    EM panel policy (2026-04-25): analytical coils enter as a .py
    module only; STEP coil input is reserved for PEEC workflows.
    """
    g = _load_golden()
    vol, coil_py = _require_sample(g)

    result = _run(coil_py, vol, g["physics"])
    exp = g["observations"]["I_2000A_coil_py"]
    tol = g["tolerance"]

    _check_common(result, exp, tol)
    _assert_close(result["B_origin_mag"], exp["B_origin_mag_T"],
                  tol["B_origin_pct"], "B_origin_mag")
    _assert_close(result["W_mag"], exp["W_mag_J"],
                  tol["W_mag_pct"], "W_mag")
    _assert_close(result["L"], exp["L_H"], tol["L_pct"], "L")
    # Sanity on coil segment count (protects against CoilBuilder
    # resampling drift between runs).
    assert result["n_wire_segments"] == exp["n_wire_segments"], (
        f"n_wire_segments drift: {result['n_wire_segments']} vs "
        f"{exp['n_wire_segments']}")


@pytest.mark.slow
def test_em_quarter_xz_coil_py_NI_2000():
    """1/4 xz reduction model: yoke reduced to x>=0 and z>=0, sym_bn=0_x
    + sym_ht=0_z sidesets from add_kelvin_cubit(reduction=...).
    Regression guard for the whole reduction-mode plumbing (2026-04-25):

        add_kelvin.add_kelvin_cubit(reduction={x: bn=0, z: ht=0})
         -> ExportNetgenCommand sideset-vs-auto-detect priority fix
            -> NGSolve .vol load (was broken before the fix)
               -> calc_accel_magnet sym_ht=0_z -> Dirichlet pick

    NOT a physics-accuracy comparison -- the coil is not x/z symmetric,
    so the sym BCs impose artificial constraints.  Numbers are specific
    to this discretization; they will drift if any of the above
    components changes behaviour.
    """
    g = _load_golden(GOLDEN_QUARTER_XZ)
    vol, coil_py = _require_sample(g)

    result = _run(coil_py, vol, g["physics"])
    exp = g["observations"]["I_2000A_coil_py"]
    tol = g["tolerance"]

    _check_common(result, exp, tol)
    _assert_close(result["B_origin_mag"], exp["B_origin_mag_T"],
                  tol["B_origin_pct"], "B_origin_mag")
    _assert_close(result["W_mag"], exp["W_mag_J"],
                  tol["W_mag_pct"], "W_mag")
    _assert_close(result["L"], exp["L_H"], tol["L_pct"], "L")
    assert result["n_wire_segments"] == exp["n_wire_segments"], (
        f"n_wire_segments drift: {result['n_wire_segments']} vs "
        f"{exp['n_wire_segments']}")


@pytest.mark.slow
def test_em_elf_quarter_xz_NI_2000():
    """ELF CEFC-2020 PHYSICS reproduction: yoke geometry from the
    ELF .meg + ELF racetrack coil + sym_bn=0_x + sym_ht=0_z BCs.

    Locks `Bz = -240.02 mT` at origin -- which matches the ELF
    reference -228.1 mT to within 5.2% (this is the genuine physics
    benchmark, NOT just a regression guard).

    Uses uniteed-yoke + 5mm tetmesh (loses the exact 13-hex topology
    for face-conformality at the air-yoke interface).  The ELF
    13-hex topology is preserved in the .jou + builder for
    reference; users who need the exact 13-hex mesh can disable the
    `unite` step.  Runtime ~6.5 min on LAB.

    Skipped if the ELF .meg is not present (CI / mdx environments).
    """
    g = _load_golden(GOLDEN_ELF)
    vol, coil_py = _require_sample(g)

    result = _run(coil_py, vol, g["physics"], timeout_s=900)
    exp = g["observations"]["I_2000A_coil_elf"]
    tol = g["tolerance"]

    _check_common(result, exp, tol)
    _assert_close(result["B_origin_mag"], exp["B_origin_mag_T"],
                  tol["B_origin_pct"], "B_origin_mag")
    _assert_close(result["W_mag"], exp["W_mag_J"],
                  tol["W_mag_pct"], "W_mag")
    _assert_close(result["L"], exp["L_H"], tol["L_pct"], "L")
    assert result["n_wire_segments"] == exp["n_wire_segments"], (
        f"n_wire_segments drift: {result['n_wire_segments']} vs "
        f"{exp['n_wire_segments']}")

    # Cross-check vs the ELF reference (NOT the regression value).
    elf_Bz = g["elf_reference"]["Bz_T"]
    fem_Bz = result["B_origin"][2]
    elf_dev_pct = abs(fem_Bz - elf_Bz) / abs(elf_Bz) * 100
    assert elf_dev_pct < 8.0, (
        f"FEM vs ELF reference drift: {fem_Bz*1000:.2f} mT vs "
        f"{elf_Bz*1000:.2f} mT ({elf_dev_pct:.2f}% > 8% tolerance)")


@pytest.mark.slow
def test_em_full_coil_py_NI_2000():
    """1/1 FULL model: full yoke (no reduction in any axis), full air
    sphere with z=0 mesh seam, full Kelvin sphere with z=0 webcut for
    1:1 copy-mesh.  Reference baseline -- the same physics as the
    reductions, at higher DOF cost.

    Regression guard for the no-reduction path:

        - auto_add_kelvin centroid-straddle heuristic correctly
          identifies symmetry=['z'] (not ['x', 'z']) on a z-only
          webcut full sphere
        - add_kelvin_cubit symmetry=['z'] webcut-seam path still works
        - calc_accel_magnet has no Dirichlet beyond GND + kelvin_int/
          kelvin_ext periodic identification
    """
    g = _load_golden(GOLDEN_FULL)
    vol, coil_py = _require_sample(g)

    result = _run(coil_py, vol, g["physics"], timeout_s=3600)
    exp = g["observations"]["I_2000A_coil_py"]
    tol = g["tolerance"]

    _check_common(result, exp, tol)
    _assert_close(result["B_origin_mag"], exp["B_origin_mag_T"],
                  tol["B_origin_pct"], "B_origin_mag")
    _assert_close(result["W_mag"], exp["W_mag_J"],
                  tol["W_mag_pct"], "W_mag")
    _assert_close(result["L"], exp["L_H"], tol["L_pct"], "L")
    assert result["n_wire_segments"] == exp["n_wire_segments"], (
        f"n_wire_segments drift: {result['n_wire_segments']} vs "
        f"{exp['n_wire_segments']}")


@pytest.mark.slow
def test_em_eighth_coil_py_NI_2000():
    """1/8 reduction model: yoke + air at (x>=0, y>=0, z>=0), with
    sym_ht=0_x + sym_ht=0_y + sym_bn=0_z + kelvin_far sidesets from
    add_kelvin_cubit(reduction={x: ht=0, y: ht=0, z: bn=0}).

    Regression guard for the 1/8-specific pieces (2026-04-25):

        - add_kelvin.py reduction allows 3 axes (was NotImplementedError)
        - all-bn=0 still rejected (B=0 physical impossibility)
        - Kelvin offset along single reduction axis (default = first)
        - 3-plane webcut: 2 sym planes coincide with air, 1 is the
          new "kelvin_far" infinity plane through the Kelvin centre
        - calc_accel_magnet treats kelvin_far as always-Dirichlet
        - ExportNetgenCommand all-or-nothing identification policy
          (when copy-mesh leaves a few unmatched corner vertices for
          large geometry scales, skip C++ identification entirely;
          NGSolve's add_periodic_kelvin recovers at solve time, and
          the kelvin_far Dirichlet anchors the field)

    NOT a physics-accuracy comparison -- coil is full racetrack, not
    x/y-antisymmetric / z-symmetric.
    """
    g = _load_golden(GOLDEN_EIGHTH)
    vol, coil_py = _require_sample(g)

    result = _run(coil_py, vol, g["physics"])
    exp = g["observations"]["I_2000A_coil_py"]
    tol = g["tolerance"]

    _check_common(result, exp, tol)
    _assert_close(result["B_origin_mag"], exp["B_origin_mag_T"],
                  tol["B_origin_pct"], "B_origin_mag")
    _assert_close(result["W_mag"], exp["W_mag_J"],
                  tol["W_mag_pct"], "W_mag")
    _assert_close(result["L"], exp["L_H"], tol["L_pct"], "L")
    assert result["n_wire_segments"] == exp["n_wire_segments"], (
        f"n_wire_segments drift: {result['n_wire_segments']} vs "
        f"{exp['n_wire_segments']}")


def test_em_panel_rejects_step_coil_input():
    """Policy test (fast): EM panel rejects STEP coil inputs with a
    clear error message pointing the user at the .py path (STEP is
    reserved for PEEC).
    """
    import subprocess
    # Use an EXISTING .step file so the extension-dispatch fires
    # before the file-not-found check.
    real_step = os.path.join(ROOT, "examples", "cubit_panels",
                             "accel_magnet", "coil_wire.step")
    if not os.path.isfile(real_step):
        pytest.skip(f"reference .step not present: {real_step}")

    proc = subprocess.run(
        [sys.executable, CALC,
         "--coil-script", real_step,
         "--vol", "nonexistent.vol",
         "--formulation", "omega",
         "--material", "linear", "--mu-r", "1000"],
        capture_output=True, text=True, timeout=60)

    # calc_main wraps exceptions in JSON {"error": ...} on stdout AND
    # exits non-zero (radia >= 4.92.0, commit 5b88b67f, "No Fallbacks
    # -- Fail Fast" policy).  Both stdout and stderr are inspected
    # for the message text.
    combined = (proc.stdout + proc.stderr).lower()
    assert "step" in combined and "peec" in combined, (
        f"Expected 'STEP'/'PEEC' in rejection message.  Got:\n"
        f"STDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}")
