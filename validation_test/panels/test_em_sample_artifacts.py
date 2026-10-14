"""Smoke tests for the bundled EM panel sample artifacts.

The EM panel ships three canonical artifacts in `panels/samples/`,
matching the panel's three input slots (Coil script / .vol / BH file):

    em_sample_coil.py   -> CoilBuilder source for `--coil-script`
    em_sample.jou + em_sample.vol -> mesh for `--vol`
    em_sample_bh.txt    -> BH curve for `--bh-file` (optional, falls
                            back to the built-in STEEL_BH table)

These are sub-second smoke tests.  The full panel chain
(coil + mesh + solver) is regression-locked by
`validation_test/panels/test_em_golden.py::test_em_sample_coil_py_NI_2000`.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
SAMPLES_DIR = REPO / "src" / "radia" / "panels" / "samples"
SAMPLE_BH = SAMPLES_DIR / "em_sample_bh.txt"
SAMPLE_COIL = SAMPLES_DIR / "em_sample_coil.py"
SAMPLE_JOU = SAMPLES_DIR / "em_sample.jou"


# -----------------------------------------------------------------
# BH curve
# -----------------------------------------------------------------

def test_bh_sample_ships():
    """The bundled BH file must exist in the source tree."""
    assert SAMPLE_BH.exists(), (
        f"em_sample_bh.txt missing at {SAMPLE_BH}; check the .gitignore "
        f"exception + pyproject.toml package-data inclusion.")


def test_bh_sample_loads():
    """The file parses via numpy.loadtxt with default whitespace delim."""
    sys.path.insert(0, str(REPO / "src" / "radia"))
    from em_material import _load_bh_file

    rows = _load_bh_file(str(SAMPLE_BH))
    assert isinstance(rows, list)
    assert len(rows) >= 50, f"too few BH points: {len(rows)}"
    Hs = [row[0] for row in rows]
    Bs = [row[1] for row in rows]
    assert Hs == sorted(Hs), "H values must be monotonically increasing"
    assert Bs == sorted(Bs), "B values must be monotonically increasing"
    assert Hs[0] == 0.0 and Bs[0] == 0.0, "curve must start at (0, 0)"
    assert Hs[-1] >= 100_000, f"saturation H too low: {Hs[-1]}"
    assert Bs[-1] >= 2.0, f"saturation B too low: {Bs[-1]}"


def test_bh_sample_matches_builtin():
    """Shipped file ~ built-in STEEL_BH (file is full-precision source)."""
    sys.path.insert(0, str(REPO / "src" / "radia"))
    from em_material import STEEL_BH, _load_bh_file

    loaded = _load_bh_file(str(SAMPLE_BH))
    assert len(loaded) == len(STEEL_BH), (
        f"point count mismatch: file={len(loaded)} builtin={len(STEEL_BH)}")
    max_rel = 0.0
    for (h_a, b_a), (h_b, b_b) in zip(loaded, STEEL_BH):
        if abs(h_b) > 1e-9:
            max_rel = max(max_rel, abs(h_a - h_b) / abs(h_b))
        if abs(b_b) > 1e-9:
            max_rel = max(max_rel, abs(b_a - b_b) / abs(b_b))
    assert max_rel < 1e-3, (
        f"shipped file deviates from built-in STEEL_BH by {max_rel:.2e} "
        f"(expected < 1e-3 since file is the full-precision source)")


def test_emmaterial_from_name_bh_file_path():
    """from_name('steel', bh_file=...) accepts the shipped path."""
    sys.path.insert(0, str(REPO / "src" / "radia"))
    from em_material import EMMaterial

    mat = EMMaterial.from_name("steel", bh_file=str(SAMPLE_BH))
    assert mat.is_nonlinear
    assert mat.bh_curve is not None
    assert len(mat.bh_curve) >= 50
    assert mat.sigma == 2e6
    assert mat.mu_r == 100.0


# -----------------------------------------------------------------
# Coil script
# -----------------------------------------------------------------

def test_coil_sample_ships():
    """The bundled coil .py must exist in the source tree."""
    assert SAMPLE_COIL.exists(), (
        f"em_sample_coil.py missing at {SAMPLE_COIL}")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_coil_sample_build_coil():
    """`build_coil()` returns a CoilBuilder with non-trivial segments.

    Suppress UserWarning: CoilBuilder's scipy `Rotation.as_euler` hits
    the canonical gimbal-lock case for axis-aligned racetracks, which
    emits a benign warning that the project's
    `filterwarnings = error` policy would otherwise turn into a
    failure.  The CoilBuilder still produces correct geometry.
    """
    sys.path.insert(0, str(REPO / "src" / "radia"))
    sys.path.insert(0, str(SAMPLES_DIR))
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "em_sample_coil", str(SAMPLE_COIL))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    coil = mod.build_coil()
    assert len(coil.segments) >= 4, (
        f"too few coil segments: {len(coil.segments)}")
    assert hasattr(coil, "write_step")


# -----------------------------------------------------------------
# yoke + air + Kelvin .jou recipe
# -----------------------------------------------------------------

def test_jou_sample_ships():
    """The bundled .jou recipe must exist."""
    assert SAMPLE_JOU.exists(), (
        f"em_sample.jou missing at {SAMPLE_JOU}")


def test_jou_sample_declares_blocks():
    """The .jou must build the canonical block triple (yoke / air / kelvin
    set up via Auto-Kelvin entry, not the .jou itself).
    """
    text = SAMPLE_JOU.read_text(encoding="utf-8")
    # The .jou itself only authors yoke + air; the Auto-Kelvin entry adds
    # the kelvin block.  Verify the .jou's piece.
    for label in ("block", "yoke", "air"):
        assert label in text.lower(), (
            f"em_sample.jou missing expected '{label}' content")


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
