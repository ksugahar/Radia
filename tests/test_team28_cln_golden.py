"""CLN-stage convergence of the TEAM 28 levitation force.

Scheduled for removal with the repository's Cauer-ladder code.  The full-FEM
ground truth, the force convention and the published height are locked by
tests/test_team28_golden.py on the CLN-free solver team28_axisym_fem.py.

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
_TEAM28 = os.path.join(_HERE, "..", "docs", "maglev", "demos", "team28")
sys.path.insert(0, _TEAM28)


@pytest.fixture(scope="module")
def forces():
    from team28_cln_force import cln_forces  # noqa: E402
    fz_full, stage_forces = cln_forces(max_stage=6)
    return fz_full, stage_forces


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
