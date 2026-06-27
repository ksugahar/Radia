"""Golden: reduced-potential Clebsch hodograph, FEM (path B increment 1).

Runs validation_test/clebsch_legendre/reduced_clebsch_hodograph_fem.py and asserts:
  - the reduced scalar potential is SOLVED by FEM and matches the exact 2-D cylinder;
  - the air field B = H_s - grad(phi) matches the exact 2-D dipole;
  - the hodograph's two coordinates are consistent: B(flux psi) == B(potential Phi),
    which requires the DUAL boundary conditions (psi Dirichlet on the far boundary
    only, iron interface Neumann/natural -- not Dirichlet).
"""
import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
EXDIR = os.path.join(HERE, "..", "clebsch_legendre")


@pytest.mark.slow
def test_reduced_clebsch_hodograph_fem_golden():
    sys.path.insert(0, EXDIR)
    try:
        import reduced_clebsch_hodograph_fem as M
    except Exception as e:  # pragma: no cover - env-dependent (ngsolve)
        pytest.skip(f"reduced-potential Clebsch hodograph (FEM) not importable: {e}")

    res, ok = M.main()
    assert ok, f"solver self-check failed: {res}"
    assert res["err_phi"] < 1e-2, f"reduced phi not solved: {res['err_phi']:.2e}"
    assert res["err_field"] < 1e-2, f"air field off exact dipole: {res['err_field']:.2e}"
    assert res["agree"] < 1e-2, f"hodograph flux/potential disagree: {res['agree']:.2e}"
