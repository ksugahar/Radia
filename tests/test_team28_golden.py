"""Golden test: the axisymmetric full-FEM TEAM 28 levitation force.

Locks the full-FEM (mixed phi-B, infinite-element shell) solve of
docs/maglev/demos/team28/team28_axisym_fem.py at dZ=0:
  - the legacy lab integral matches the lab ground truth -2.1928 N (band 0.5%);
  - it is an upward lift of order 2 N;
  - the legacy integral Re[B_r J_t] is 2x the physical time-averaged force, and
    the physical lift at dZ=0 is close to the disk weight, consistent with the
    published 11.3 mm stationary height.

Runs a real axisymmetric NGSolve eddy-current solve; skipped cleanly if ngsolve /
netgen are not importable in the active env.
"""
import json
import os
import sys
from pathlib import Path

import pytest

pytest.importorskip("ngsolve")
pytest.importorskip("netgen.occ")

_HERE = os.path.dirname(os.path.abspath(__file__))
_TEAM28 = os.path.join(_HERE, "..", "docs", "maglev", "demos", "team28")
sys.path.insert(0, _TEAM28)

LAB_REF = -2.1928   # N, lab full-FEM ground truth at dZ=0 (legacy integral)
REFERENCE = json.loads(
    (Path(_HERE).parent / "validation_test/maglev/team28_reference.json").read_text(encoding="utf-8"))
PUB_LEVITATION_MM = REFERENCE["stationary_height_mm"]
DISK_BOTTOM_DZ0_MM = 10.8
DISK_WEIGHT_N = 1.055


@pytest.fixture(scope="module")
def forces():
    from team28_axisym_fem import solve_force_pair  # noqa: E402
    return solve_force_pair()


def test_full_fem_matches_lab_ground_truth(forces):
    f_legacy, _ = forces
    rel = abs(f_legacy - LAB_REF) / abs(LAB_REF)
    assert rel < 0.005, (
        f"full-FEM force {f_legacy:.4f} N deviates {rel*100:.2f}% from the "
        f"lab ground truth {LAB_REF} N (band 0.5%)")


def test_force_is_upward_lift_of_order_2N(forces):
    f_legacy, _ = forces
    # sign convention: negative = upward lift; magnitude ~2.19 N at dZ=0
    assert f_legacy < 0.0
    assert 2.0 < abs(f_legacy) < 2.4


def test_model_parameters_match_the_published_definition():
    import team28_axisym_fem as model

    assert PUB_LEVITATION_MM == 11.3  # official Model A, section II
    for name, key in (("N1", "inner_coil_turns"), ("N2", "outer_coil_turns"),
                      ("I1", "current_peak_A"), ("I2", "current_peak_A"),
                      ("FREQ", "frequency_Hz"), ("SIGMA_AL", "conductivity_S_per_m"),
                      ("aluminium_w", "disk_radius_m"), ("aluminium_h", "disk_thickness_m"),
                      ("coil_1_r", "inner_coil_center_radius_m"),
                      ("coil_2_r", "outer_coil_center_radius_m"),
                      ("coil_1_w", "inner_coil_width_m"),
                      ("coil_2_w", "outer_coil_width_m"),
                      ("coil_1_h", "coil_height_m"), ("coil_2_h", "coil_height_m")):
        assert getattr(model, name) == pytest.approx(REFERENCE[key]), name


def test_force_convention_2x_and_published_height(forces):
    """The legacy integral is 2x the physical force; physical lift ~ weight at dZ=0.

    Regression guard for the 2026-06-20 convention fix: balancing the 2x integral
    against the 1x weight gave a spurious 14.9 mm height.
    """
    f_legacy, f_phys = forces
    assert abs(abs(f_legacy / f_phys) - 2.0) < 0.02, \
        f"legacy/physical = {f_legacy/f_phys:.3f} (expect ~2.0)"
    assert 0.9 < abs(f_phys) / DISK_WEIGHT_N < 1.2, \
        f"physical lift {abs(f_phys):.3f} N vs weight {DISK_WEIGHT_N} N"
    assert abs(DISK_BOTTOM_DZ0_MM - PUB_LEVITATION_MM) / PUB_LEVITATION_MM < 0.10
