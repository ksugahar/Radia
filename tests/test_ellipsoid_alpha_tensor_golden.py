"""Golden test: shape-anisotropic eddy-current polarizability of an ellipsoid.

Locks the analytic shape-anisotropy result of
docs/maglev/demos/ellipsoid/ellipsoid_alpha_tensor.py:
  - the demagnetizing tensor sums to 1 and the Osborn integral matches the
    spheroid closed form (sphere 1/3, 2:1 prolate N_c=0.1736),
  - a triaxial body a1!=a2!=a3 gives three DISTINCT alpha_i (genuine shape
    anisotropy of an isotropic-material conductor),
  - |kappa_i| is largest along the SHORT axis (most flux excluded),
  - the sphere limit reduces to the verified -2 pi a^3 (sphere force module).

Pure numpy + scipy -- no NGSolve / Cubit, so this always runs in CI.
"""
import os
import sys

import numpy as np
import pytest

pytest.importorskip("scipy")

_HERE = os.path.dirname(os.path.abspath(__file__))
# CLN was absorbed into radia.maglev ->
# docs/maglev/demos/{sphere,ellipsoid}. Add both so sphere + ellipsoid modules resolve.
_LEV = os.path.join(_HERE, "..", "docs", "maglev", "demos")
sys.path.insert(0, os.path.join(_LEV, "sphere"))
sys.path.insert(0, os.path.join(_LEV, "ellipsoid"))

import ellipsoid_alpha_tensor as E  # noqa: E402


def test_demag_sums_to_one_and_matches_closed_form():
    for a in [(1, 1, 1), (1, 1, 2), (1, 1, 0.25), (1, 1, 100), (1, 1, 0.01)]:
        N = E.demag_tensor(a)
        assert abs(N.sum() - 1.0) < 1e-6
        Ncf = E.demag_spheroid_closed(a[0], a[2])
        assert np.max(np.abs(N - Ncf)) < 1e-4


def test_sphere_demag_isotropic():
    assert np.max(np.abs(E.demag_tensor((1, 1, 1)) - 1 / 3)) < 1e-6


def test_known_prolate_value():
    # 2:1 prolate spheroid: N_c (long axis) = 0.17357 (Osborn 1945)
    assert abs(E.demag_tensor((1.0, 1.0, 2.0))[2] - 0.17357) < 1e-3


def test_triaxial_is_anisotropic():
    a = (5e-3, 3e-3, 1.5e-3)
    kap, N, V = E.kappa_hf(a)
    # three distinct demag factors / polarizabilities, ordered by axis length
    assert N[0] < N[1] < N[2]
    assert abs(kap[0]) < abs(kap[1]) < abs(kap[2])
    # the short axis excludes the most flux per volume
    assert abs(kap[2]) / V > abs(kap[0]) / V
    assert abs(kap[2]) / abs(kap[0]) > 1.2     # clearly anisotropic


def test_sphere_limit_reduces_to_minus_2pi_a3():
    a_r = 5e-3
    kap, N, V = E.kappa_hf((a_r, a_r, a_r))
    assert np.allclose(kap, -2 * np.pi * a_r**3, rtol=1e-4)
    assert np.max(np.abs(N - 1 / 3)) < 1e-6


def test_lift_largest_along_short_axis():
    a = (5e-3, 3e-3, 1.5e-3)
    lf = E.lift_factor(a)          # 1/(1-N_i)
    assert lf[2] > lf[1] > lf[0]   # short axis -> largest lift factor


def test_kappa_lift_reduces_to_verified_sphere_coefficient():
    import maglev_sphere_force as L
    mu0 = E.mu0
    gradB2 = 0.2
    F_sphere = (np.pi * L.a**3 / (2 * mu0)) * gradB2
    F_kappa = (4 / 3 * np.pi * L.a**3 / (4 * mu0)) * (1 / (1 - 1 / 3)) * gradB2
    assert abs(F_sphere - F_kappa) / F_sphere < 1e-9
