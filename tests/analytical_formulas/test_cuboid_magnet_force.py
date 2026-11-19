"""Guards: closed-form cuboid-magnet force (Akoun-Yonnet) in radia.analytical_formulas.

Assumes: the surface-charge model (rigid uniform magnetisation).  The oracle is
independent of the closed-form algebra: direct Gauss-Legendre integration of
Coulomb forces between the z-face charge sheets, plus the far-field dipole limit.
"""
import math

import numpy as np
import pytest

from radia.analytical_formulas.cuboid_magnet_force import MU0, cuboid_magnet_force


def _charge_quadrature(J1, s1, J2, s2, d, n=40):
    g, w = np.polynomial.legendre.leggauss(n)

    def face(cx, cy, sx, sy):
        X, Y = np.meshgrid(cx + 0.5 * sx * g, cy + 0.5 * sy * g, indexing="ij")
        return X.ravel(), Y.ravel(), (np.outer(w, w) * 0.25 * sx * sy).ravel()

    F = np.zeros(3)
    x1, y1, w1 = face(0.0, 0.0, s1[0], s1[1])
    x2, y2, w2 = face(d[0], d[1], s2[0], s2[1])
    for z1, q1 in ((0.5 * s1[2], J1), (-0.5 * s1[2], -J1)):
        for z2, q2 in ((d[2] + 0.5 * s2[2], J2), (d[2] - 0.5 * s2[2], -J2)):
            R = np.stack([x2[None, :] - x1[:, None], y2[None, :] - y1[:, None],
                          np.full((len(x1), len(x2)), z2 - z1)])
            r3 = np.sqrt((R ** 2).sum(0)) ** 3
            F += q1 * q2 * (R / r3 * (w1[:, None] * w2[None, :])).sum((1, 2))
    return F / (4 * math.pi * MU0)


CASES = [
    (1.2, (0.01, 0.01, 0.01), 1.2, (0.01, 0.02, 0.005), (0.004, 0.003, 0.012)),
    (1.0, (0.02, 0.01, 0.004), -0.8, (0.006, 0.006, 0.01), (-0.01, 0.002, 0.009)),
    (1.3, (0.01, 0.01, 0.01), 1.3, (0.01, 0.01, 0.01), (0.0, 0.0, 0.0125)),
    (1.1, (0.01, 0.01, 0.01), 1.1, (0.01, 0.01, 0.01), (0.012, 0.0, -0.003)),   # side by side, below
    (1.0, (0.008, 0.012, 0.006), 1.0, (0.004, 0.004, 0.004), (0.003, -0.004, -0.009)),
    (1.0, (0.01, 0.01, 0.01), 1.0, (0.01, 0.01, 0.01), (0.02, 0.03, 0.09)),   # near the accepted far edge
]


@pytest.mark.parametrize("J1,s1,J2,s2,d", CASES)
def test_matches_independent_charge_quadrature(J1, s1, J2, s2, d):
    got = cuboid_magnet_force(J1, s1, J2, s2, d)
    ref = _charge_quadrature(J1, s1, J2, s2, d)
    err = np.abs(got - ref).max() / np.linalg.norm(ref)
    assert err < 1e-6, f"closed form {got} vs quadrature {ref} (rel {err:.2e})"


def test_far_field_dipole_limit_and_newton_third_law():
    s1, s2 = (0.01, 0.01, 0.01), (0.01, 0.02, 0.005)
    dz = 0.1
    m1 = 1.0 / MU0 * np.prod(s1); m2 = 1.0 / MU0 * np.prod(s2)
    got = cuboid_magnet_force(1.0, s1, 1.0, s2, (0.0, 0.0, dz))[2]
    ref = -3 * MU0 * m1 * m2 / (2 * math.pi * dz ** 4)
    assert abs(got / ref - 1) < 3e-2, f"F_z {got} vs dipole {ref}"  # finite-size term ~2 % at d = 10 s
    d = (0.004, 0.003, 0.012)
    f12 = cuboid_magnet_force(1.0, s1, 1.0, s2, d)
    f21 = cuboid_magnet_force(1.0, s2, 1.0, s1, tuple(-x for x in d))
    assert np.allclose(f12, -f21, rtol=1e-10, atol=0), f"{f12} vs -{f21}"


def test_far_separation_is_rejected_not_sign_flipped():
    # Two 10 mm cubes 10 m apart: the corner sum cancels completely in double precision.
    # The result must never be returned with the wrong sign; the formula refuses instead.
    with pytest.raises(ValueError, match="separation too large"):
        cuboid_magnet_force(1.0, (0.01,) * 3, 1.0, (0.01,) * 3, (0.0, 0.0, 10.0))
    for dz in (0.03, 0.05, 0.1):
        f = cuboid_magnet_force(1.0, (0.01,) * 3, 1.0, (0.01,) * 3, (0.0, 0.0, dz))[2]
        m = 1e-6 / MU0
        assert f < 0 and abs(f / (-3 * MU0 * m * m / (2 * math.pi * dz ** 4)) - 1) < 2e-2, (dz, f)


def test_rejects_overlap_and_bad_sizes():
    with pytest.raises(ValueError, match="overlap"):
        cuboid_magnet_force(1, (0.01,) * 3, 1, (0.01,) * 3, (0.005, 0, 0))
    with pytest.raises(ValueError, match="positive"):
        cuboid_magnet_force(1, (0.0, 0.01, 0.01), 1, (0.01,) * 3, (0, 0, 0.1))
