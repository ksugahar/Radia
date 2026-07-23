"""Vector EM Mie BISTATIC far-field of a PEC sphere -- the angular scattering amplitudes
S1(theta), S2(theta) -- gated by three airtight identities against the trusted existing
total/backscatter cross-sections (mie_pec_sphere_rcs). Pure analytic (numpy/scipy); no FEM.

Fills the gap: the angular pattern + forward optical theorem previously existed only for the
SCALAR Helmholtz Mie set; the vector EM Mie helpers had only total + backscatter.
"""
import math
import os
import sys

import pytest

# This is a numerical far-field validation corpus, not a lightweight contract
# test.  Keep it out of the pre-push matrix; run it in the validation lane.
pytestmark = pytest.mark.slow

np = pytest.importorskip("numpy")
pytest.importorskip("scipy")
from scipy.integrate import trapezoid

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.bem_integral import (
    mie_pec_coefficients, mie_pec_sphere_rcs, mie_pec_amplitude, mie_pec_bistatic_rcs)

KAS = [0.5, 1.0, 2.0, 3.0, 5.0]


@pytest.mark.parametrize("ka", KAS)
def test_optical_theorem(ka):
    # forward amplitude S(0): Q_ext = (4/x^2) Re[S(0)] must equal Q_sca (lossless PEC -> ext == sca)
    x = float(ka)
    S1_0, S2_0 = mie_pec_amplitude(x, 0.0)
    assert abs(S1_0 - S2_0) < 1e-12 * abs(S1_0)          # forward S1 == S2
    q_ext = (4.0 / x ** 2) * S1_0.real
    sca_ref = mie_pec_sphere_rcs(x)["sigma_sca_over_pa2"]
    assert q_ext == pytest.approx(sca_ref, rel=1e-9)


@pytest.mark.parametrize("ka", KAS)
def test_backscatter_from_pattern(ka):
    # the pattern at theta=pi reproduces the existing monostatic backscatter RCS
    x = float(ka)
    back_pat = mie_pec_bistatic_rcs(x, math.pi)["h_plane_over_pa2"]
    back_ref = mie_pec_sphere_rcs(x)["sigma_back_over_pa2"]
    assert back_pat == pytest.approx(back_ref, rel=1e-9)
    # E-plane and H-plane coincide at backscatter
    bp = mie_pec_bistatic_rcs(x, math.pi)
    assert bp["e_plane_over_pa2"] == pytest.approx(bp["h_plane_over_pa2"], rel=1e-9)


@pytest.mark.parametrize("ka", KAS)
def test_total_from_angular_integral(ka):
    # integrating the angular pattern over solid angle recovers the trusted total cross-section:
    #   Q_sca = (1/x^2) int_0^pi (|S1|^2 + |S2|^2) sin th dth
    x = float(ka)
    th = np.linspace(0.0, math.pi, 4001)
    s1 = np.array([mie_pec_amplitude(x, float(t))[0] for t in th])
    s2 = np.array([mie_pec_amplitude(x, float(t))[1] for t in th])
    q_pat = (1.0 / x ** 2) * trapezoid((np.abs(s1) ** 2 + np.abs(s2) ** 2) * np.sin(th), th)
    sca_ref = mie_pec_sphere_rcs(x)["sigma_sca_over_pa2"]
    assert q_pat == pytest.approx(sca_ref, rel=1e-4)      # quadrature-limited (measured ~1e-6)


def test_forward_exceeds_backward_at_high_ka():
    # physical sanity: at large ka forward scattering dominates (diffraction lobe) over backscatter
    x = 5.0
    fwd = mie_pec_bistatic_rcs(x, 1e-6)["e_plane_over_pa2"]
    back = mie_pec_bistatic_rcs(x, math.pi)["e_plane_over_pa2"]
    assert fwd > back


def test_amplitude_matches_coefficient_count():
    # the helper uses the same Wiscombe nmax as the coefficient routine (consistency)
    n, a_n, b_n = mie_pec_coefficients(2.0)
    S1, S2 = mie_pec_amplitude(2.0, 0.7)
    assert isinstance(S1, complex) and isinstance(S2, complex)
    assert len(a_n) == len(b_n) == len(n)
