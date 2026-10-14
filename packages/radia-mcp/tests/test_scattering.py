r"""High-frequency physical-optics RCS of canonical PEC shapes -- closed-form regression test.

Pure textbook physical optics (Knott / Shaeffer / Tuley *Radar Cross Section* Ch. 5-6,
Balanis, Jackson): the flat-plate sinc-squared pattern and the dihedral / trihedral corner
reflector peaks.  Analytic helpers only (no field solve), so this is a closed-form check:
exact value at a concrete point, one non-tautological independent cross-check, and the
ValueError guard for each function.
"""
import math
import os
import sys

import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.scattering import (
    physical_optics_flat_plate_rcs,
    dihedral_corner_reflector_peak_rcs,
)


def test_flat_plate_rcs():
    # exact broadside value: a 1 m x 1 m plate at lambda = 1 m -> sigma_max = 4 pi exactly.
    out = physical_optics_flat_plate_rcs(1.0, 1.0, 1.0, theta_deg=0.0)
    assert math.isclose(out["rcs_max"], 4.0 * math.pi, rel_tol=1e-12)
    assert math.isclose(out["rcs"], 4.0 * math.pi, rel_tol=1e-12)        # sinc(0)=1, cos0=1
    # dBsm is exactly 10*log10(rcs) (defining identity, here vs an independent recompute).
    assert math.isclose(out["rcs_dbsm"], 10.0 * math.log10(out["rcs"]), rel_tol=1e-12)

    # NON-TAUTOLOGICAL independent gate: the first PO pattern null is set by the sinc zero
    # k*a*sin(theta) = pi  =>  theta = arcsin(lambda / (2 a)).  Build that angle from the
    # sinc root (no reference to the sigma formula) and demand the pattern vanishes there.
    a, lam = 2.0, 1.0
    theta_null = math.degrees(math.asin(lam / (2.0 * a)))
    null = physical_optics_flat_plate_rcs(a, 3.0, lam, theta_deg=theta_null)
    assert null["rcs"] < 1e-12 * null["rcs_max"]                          # deep null
    assert null["rcs_dbsm"] < -100.0                                      # >100 dB below peak
    # just inside the null the pattern is still strongly down from broadside (a real lobe).
    near = physical_optics_flat_plate_rcs(a, 3.0, lam, theta_deg=theta_null * 0.9)
    assert 0.0 < near["rcs"] < 0.1 * near["rcs_max"]

    # bad input -> ValueError
    with pytest.raises(ValueError):
        physical_optics_flat_plate_rcs(0.0, 1.0, 1.0)
    with pytest.raises(ValueError):
        physical_optics_flat_plate_rcs(1.0, 1.0, -0.5)


def test_dihedral_corner_reflector_rcs():
    a, h, lam = 0.30, 0.30, 0.03
    out = dihedral_corner_reflector_peak_rcs(a, h, lam)
    # exact closed forms
    assert math.isclose(out["rcs_dihedral"], 8.0 * math.pi * a**2 * h**2 / lam**2, rel_tol=1e-12)
    assert math.isclose(out["rcs_trihedral"], 12.0 * math.pi * a**4 / lam**2, rel_tol=1e-12)
    assert math.isclose(out["rcs_dihedral_dbsm"], 10.0 * math.log10(out["rcs_dihedral"]),
                        rel_tol=1e-12)

    # NON-TAUTOLOGICAL independent gate: the dihedral peak is exactly TWICE the broadside RCS
    # of a single flat plate of the SAME lit area a x h (8 pi vs 4 pi prefactor) -- the
    # double-bounce signature.  Cross-check against the *other* function, not its own formula.
    plate = physical_optics_flat_plate_rcs(a, h, lam, theta_deg=0.0)
    assert math.isclose(out["rcs_dihedral"], 2.0 * plate["rcs_max"], rel_tol=1e-12)

    # scaling law (independent): both peaks ~ 1/lambda^2 -> halving lambda quadruples RCS.
    out_half = dihedral_corner_reflector_peak_rcs(a, h, lam / 2.0)
    assert math.isclose(out_half["rcs_dihedral"], 4.0 * out["rcs_dihedral"], rel_tol=1e-12)

    # bad input -> ValueError
    with pytest.raises(ValueError):
        dihedral_corner_reflector_peak_rcs(-1.0, h, lam)
    with pytest.raises(ValueError):
        dihedral_corner_reflector_peak_rcs(a, h, 0.0)


if __name__ == "__main__":
    test_flat_plate_rcs()
    test_dihedral_corner_reflector_rcs()
    print("[OK] physical-optics RCS: flat-plate sinc pattern + dihedral/trihedral peaks.")
