"""Golden: the T-Omega cohomology CUT solve (gmsh-free) reproduces the wire field.

Exercises CohomologyCutSolver.solve() -- not just setup_from_mesh -- the
PRIMARY use of radia.cohomology (the T-Omega total-scalar-potential cut, of
which the Clebsch current-linking demo is one instance).  A straight wire
(current I along z) through an annular air region (b1 = 1):
  * Ampere's law  oint H.dl = NI  (exact, by the unit-circulation cut);
  * the wire field  H_phi = I / (2 pi r)  to FEM accuracy.
"""
import sys
from pathlib import Path

import pytest

pytest.importorskip("ngsolve")
pytest.importorskip("netgen.occ")

EXDIR = Path(__file__).resolve().parents[2] / "examples" / "cohomology"
sys.path.insert(0, str(EXDIR))


def test_tomega_wire_cohomology_cut():
    import tomega_wire as tw
    r = tw.solve(maxh=0.010)                       # slightly coarser for CI
    assert r["b1"] == 1, r                         # one current loop (annulus)
    assert r["ampere_rel_err"] < 1e-3, r           # cut carries NI exactly
    assert r["wire_field_err"] < 0.05, r           # H_phi ~ I/(2 pi r)
