"""Golden test: 3D (axisymmetric) CLN reproduces the TEAM 28 levitation force.

Locks the verified result of examples/maglev/team28/:
  - the full-FEM (split K+sN) levitation force at dZ=0 == lab ground truth
    -2.1928 N (within a hard band),
  - the N-stage CLN/Cauer reduced force CONVERGES to it (stage 3 < 1%,
    stage 5 <= 0.05%).

Runs a real axisymmetric NGSolve eddy-current solve (~20-40 s); skipped
cleanly if ngsolve / netgen are not importable in the active env.
"""
import os
import sys

import pytest

pytest.importorskip("ngsolve")
pytest.importorskip("netgen.occ")
pytest.importorskip("scipy")

_HERE = os.path.dirname(os.path.abspath(__file__))
_TEAM28 = os.path.join(_HERE, "..", "examples", "maglev", "team28")
sys.path.insert(0, _TEAM28)

LAB_REF = -2.1928   # N, lab full-FEM ground truth at dZ=0 (50Hz_可動 .mat)


@pytest.fixture(scope="module")
def forces():
    from team28_cln_force import cln_forces  # noqa: E402
    fz_full, stage_forces = cln_forces(max_stage=6)
    return fz_full, stage_forces


def test_full_fem_matches_lab_ground_truth(forces):
    fz_full, _ = forces
    rel = abs(fz_full - LAB_REF) / abs(LAB_REF)
    assert rel < 0.005, (
        f"full-FEM force {fz_full:.4f} N deviates {rel*100:.2f}% from the "
        f"lab ground truth {LAB_REF} N (band 0.5%)")


def test_force_is_upward_lift_of_order_2N(forces):
    fz_full, _ = forces
    # sign convention: negative = upward lift; magnitude ~2.19 N at dZ=0
    assert fz_full < 0.0
    assert 2.0 < abs(fz_full) < 2.4


def test_cln_converges_to_full(forces):
    fz_full, sf = forces
    assert len(sf) >= 5, "expected at least 5 CLN stages"
    err = [abs(f - fz_full) / abs(fz_full) for f in sf]
    # stage 1 is the eddy-free DC response -> large error
    assert err[0] > 0.5
    # convergence: by stage 3 within 1%, by stage 5 within 0.05%
    assert err[2] < 0.01, f"stage 3 rel err {err[2]*100:.3f}% (expect <1%)"
    assert err[4] < 5e-4, f"stage 5 rel err {err[4]*100:.4f}% (expect <0.05%)"
    # monotone-ish: stage 5 is at least as good as stage 3
    assert err[4] <= err[2]
