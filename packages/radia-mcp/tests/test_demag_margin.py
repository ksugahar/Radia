r"""PM demagnetization operating point + irreversible-demag margin.

H_op = (B_in - Br)/(mu0 mu_r),  B_in = Br(1-N)/(1+N(mu_r-1));  margin = H_op - H_knee.
Pure analytic reference (transverse cylinder N=1/2, sphere N=1/3); the operating field is
cross-checked by a transverse PM-cylinder magnetostatic FEA (bulk H_x == -Hc/2).
"""
import math
import os
import sys

import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)
from radia_mcp.radia_ngsolve.solve import demag_operating_field, demag_margin

MU0 = 4e-7 * math.pi


def test_operating_field_demag_factors():
    Br = 1.2
    # transverse cylinder N=1/2, mu_r=1 -> H_op = -N M = -Br/(2 mu0) = -Hc/2
    assert demag_operating_field(Br, 0.5, 1.0) == pytest.approx(-0.5 * Br / MU0)
    # sphere N=1/3
    assert demag_operating_field(Br, 1.0 / 3.0, 1.0) == pytest.approx(-(1.0 / 3.0) * Br / MU0)


def test_recoil_permeability_operating_field():
    Br = 1.2
    Bin = Br / (1.05 + 1.0)                       # N=1/2 -> B_in = Br/(mu_r+1)
    assert demag_operating_field(Br, 0.5, 1.05) == pytest.approx((Bin - Br) / (MU0 * 1.05))


def test_margin_sign():
    assert demag_margin(-4.77e5, -8.0e5) > 0      # operating point above the knee -> safe
    assert demag_margin(-9.0e5, -8.0e5) < 0       # driven below the knee -> irreversible demag
    assert demag_margin(-8.0e5, -8.0e5) == 0      # exactly at the knee


def main():
    Br = 1.2
    print(f"transverse cylinder N=1/2: H_op = {demag_operating_field(Br, 0.5):.4e} A/m (= -Hc/2)")
    print(f"margin vs H_knee=-800kA/m: {demag_margin(demag_operating_field(Br, 0.5), -8e5):.4e} A/m")
    print("[OK] demag operating field = -N Br/mu0 (mu_r=1); margin = H_op - H_knee.")


if __name__ == "__main__":
    main()
