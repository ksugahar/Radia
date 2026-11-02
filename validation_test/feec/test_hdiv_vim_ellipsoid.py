"""Golden test: HDiv-VIM demag factor on a SMOOTH NON-SPHERICAL body (production #3).

The sphere/cube tests fix the demag factor at the isotropic 1/3.  A prolate spheroid has a non-trivial
ANALYTIC demag factor != 1/3 (e.g. a 2:1 prolate has N_z ~ 0.1736 along the long axis), so it checks
that the production C++ analytic charge-Gram path gets a non-1/3 demag right -- broadening the operator's
validated domain to smooth bodies with anisotropic demag.

The C++ analytic charge Gram lands within ~0.3% of the analytic N_z, confirming the analytic surface
integral is correct for a general (non-equilateral, curved-surface) triangulation.

Input: the checked named 1338-TET fixture (maxh 1.0) shared with the regression lock, not a mesh
generated inside the test; the 1.5% criterion is unchanged.
"""
import hashlib
import os
import sys
from math import log, sqrt
from pathlib import Path

import pytest

pytest.importorskip("ngsolve")

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "vim_legacy"))
from radia.vim import Solve  # noqa: E402
import ngsolve as ng  # noqa: E402


def _prolate_Nz_analytic(c, a=1.0):
    """Demag factor along the long (z) axis of a prolate spheroid (semi-axes a=b<c)."""
    e = sqrt(1.0 - (a / c) ** 2)
    return (1.0 - e * e) / e ** 2 * ((1.0 / (2.0 * e)) * log((1.0 + e) / (1.0 - e)) - 1.0)


PROLATE_LABELS = (("magnetic_body",), ("exterior",))
PROLATE_FIXTURE = Path(__file__).resolve().parents[2] / "tests" / "feec" / "data" / "prolate_named"
# Checked VOL handoff: mesh sha256 -> sha256 of its recorded strict label-check report
# against prolate-label-contract.json. No other mesh reaches Solve.
CHECKED_PROLATE_VOLS = {
    "27c3778663b1124e37ee5268e8fbd998d62d3550de5572378cd3eaf3bd6c109d":
        "a85188adc458c6e0faf3e79f4a7bbbffba745c2a1e37e5fd60a7ee3cd2de36c1",
}


def run_prolate_case(mesh_path):
    """Run one production case on a checked named prolate VOL; return the production output.

    Only a VOL whose sha256 is in CHECKED_PROLATE_VOLS is loaded; any other file raises
    before Solve. There is no generation, remeshing or save path. The loaded labels must
    still be exactly PROLATE_LABELS. The analytic oracle belongs only to the heavy
    validation entry below.
    """
    digest = hashlib.sha256(Path(mesh_path).read_bytes()).hexdigest()
    if digest not in CHECKED_PROLATE_VOLS:
        raise ValueError(f"prolate mesh sha256: actual={digest}; expected one of "
                         f"{sorted(CHECKED_PROLATE_VOLS)} (checked VOL handoff required)")
    with ng.TaskManager():
        mesh = ng.Mesh(str(mesh_path))
        labels = (tuple(mesh.GetMaterials()), tuple(sorted(set(mesh.GetBoundaries()))))
        if labels != PROLATE_LABELS:
            raise ValueError(f"prolate mesh labels: actual={labels}; expected={PROLATE_LABELS}")
        return Solve(mesh, mu_r=1000.0, H_ext=ng.CoefficientFunction((0, 0, 1.0)))


def assert_prolate_validation(Dw):
    """Existing independent heavy criterion; never used by a regression lock."""
    Na = _prolate_Nz_analytic(2.0)
    assert 0.15 < Na < 0.20, f"analytic prolate N_z sanity: {Na:.4f}"
    assert abs(Dw - Na) < 1.5e-2 * Na, f"prolate N_z {Dw:.4f} not within 1.5% of analytic {Na:.4f}"


def test_prolate_spheroid_demag_factor():
    assert_prolate_validation(run_prolate_case(PROLATE_FIXTURE / "prolate_named.vol")["demag"])