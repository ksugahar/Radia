"""Opt-in oriented-box Galerkin kernel of the PEEC builder and the energy gate."""
from __future__ import annotations

import math

import numpy as np
import pytest

MU0_4PI = 1e-7


def _matrix(bars, **kw):
    from radia.peec_matrices import PEECBuilder

    b = PEECBuilder(**kw)
    for p0, p1, w, h, wd in bars:
        n0 = b.add_node_at(*map(float, p0))
        n1 = b.add_node_at(*map(float, p1))
        b.add_connected_segment(n0, n1, float(w), float(h))
        if wd is not None and kw.get("oriented_box_order"):
            b.set_last_width_dir(*map(float, wd))
    return np.asarray(b.build_topology()["L"])


def _rot(axis, angle):
    axis = np.asarray(axis, float) / np.linalg.norm(axis)
    k = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
    return np.eye(3) + math.sin(angle) * k + (1 - math.cos(angle)) * k @ k


def _pair(offset, *, rot=np.eye(3)):
    """Two 10 mm z-bars, 1 mm wide along y, 0.2 mm high along x, offset in x."""
    out = []
    for x in (0.0, offset):
        p0, p1 = rot @ [x, 0, 0], rot @ [x, 0, 10e-3]
        out.append((p0, p1, 1e-3, 0.2e-3, rot @ [0, 1, 0]))
    return out


def test_default_builder_is_unchanged():
    """Guards: the kernel is opt-in; the default and order 0 give
    bit-identical L.
    Assumes: two parallel bars 0.2 mm apart."""
    bars = _pair(0.4e-3)
    assert np.array_equal(_matrix(bars), _matrix(bars, oriented_box_order=0))


@pytest.mark.parametrize("gap", [0.5, 0.2e-3, 0.0, -0.1e-3, -0.2e-3 + 1e-6])
def test_aligned_pairs_match_the_cartesian_volume_rule(gap):
    """Guards: self and mutual entries for separated, touching and
    overlapping parallel bars equal the independent Laplace-transform
    Cartesian volume rule within 5e-4, also after an arbitrary rotation
    (frame handling), and the 2x2 matrix is positive semidefinite.
    Assumes: Cartesian rule order 64 is converged to 1e-5 for these bars."""
    ref = _matrix(_pair(0.2e-3 + gap), volume_order=64)
    for rot in (np.eye(3), _rot([1, 2, 3], 0.7)):
        got = _matrix(_pair(0.2e-3 + gap, rot=rot), oriented_box_order=8)
        assert np.allclose(got, ref, rtol=5e-4, atol=0), (gap, got, ref)
        ev = np.linalg.eigvalsh(got)
        assert ev[0] >= -4 * np.finfo(float).eps * ev[-1], (gap, ev)


def _gauss_box(c, axes, half, n, sub):
    x, w = np.polynomial.legendre.leggauss(n)
    edges = np.linspace(-1, 1, sub + 1)
    pts, wts = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        pts.append(0.5 * (hi - lo) * x + 0.5 * (hi + lo))
        wts.append(0.5 * (hi - lo) * w)
    t, tw = np.concatenate(pts), np.concatenate(wts)
    u, v, s = np.meshgrid(t, t, t, indexing="ij")
    wt = (tw[:, None, None] * tw[None, :, None] * tw[None, None, :]).ravel()
    loc = np.stack([u.ravel() * half[0], v.ravel() * half[1], s.ravel() * half[2]], 1)
    return c + loc @ np.asarray(axes), wt * np.prod(half)


def test_skew_plates_match_an_independent_six_dimensional_integral():
    """Guards: a skew pair of short, wide plates (the transverse bars of
    the fin lattice: 15 degrees apart, 10 um apart at the near corner)
    matches a plain 6-D Gauss-Legendre integral of 1/r within 2e-3, with
    the d_i . d_j orientation factor and its sign reversal.
    Assumes: 12-point composite Gauss per axis resolves the 10 um gap."""
    length, width, height = 0.5e-3, 1e-3, 0.1e-3
    th = math.radians(15)
    da = np.array([1.0, 0, 0])
    db = np.array([math.cos(th), math.sin(th), 0])
    ez = np.array([0, 0, 1.0])
    ca = 0.5 * length * da
    hb = np.cross(db, ez)
    cb = np.array([0, height / 2 + 10e-6 + height / 2, 0]) + 0.5 * length * db
    bars = [(ca - 0.5 * length * da, ca + 0.5 * length * da, width, height, ez),
            (cb - 0.5 * length * db, cb + 0.5 * length * db, width, height, ez)]
    got = _matrix(bars, oriented_box_order=8)

    pa, wa = _gauss_box(ca, [da, ez, np.cross(da, ez)], (length / 2, width / 2, height / 2), 6, 2)
    pb, wb = _gauss_box(cb, [db, ez, hb], (length / 2, width / 2, height / 2), 6, 2)
    integral = 0.0
    for chunk in range(0, len(pa), 512):
        d = np.linalg.norm(pa[chunk:chunk + 512, None] - pb[None], axis=2)
        integral += float(wa[chunk:chunk + 512] @ (1 / d) @ wb)
    area = width * height
    ref = MU0_4PI * float(da @ db) * integral / area ** 2
    assert abs(got[0, 1] - ref) <= 2e-3 * abs(ref), (got[0, 1], ref)

    flipped = _matrix([bars[0], (bars[1][1], bars[1][0], width, height, ez)], oriented_box_order=8)
    assert flipped[0, 1] == -got[0, 1]
    assert np.linalg.eigvalsh(got)[0] > 0


def test_overlapping_support_keeps_positive_energy():
    """Guards: the overlap case that made the filament average indefinite
    (two 10 mm bars, centres 1 um apart) gives |L_ij| <= L_ii and no
    negative eigenvalue beyond roundoff.
    Assumes: identical bars, width axis x."""
    bars = []
    for x in (0.0, 1e-6):
        bars.append(((x, 0, 0), (x, 0, 10e-3), 1e-3, 0.2e-3, (1, 0, 0)))
    m = _matrix(bars, oriented_box_order=8)
    assert abs(m[0, 1]) <= m[0, 0] * (1 + 1e-12), m
    ev = np.linalg.eigvalsh(m)
    assert ev[0] >= -4 * np.finfo(float).eps * ev[-1], ev


def test_energy_gate_separates_null_modes_from_negative_energy():
    """Guards: the fin-volume admission accepts positive definite and
    roundoff-null (duplicate basis) matrices, records the null modes, and
    raises with a witness on negative energy; it never loads the diagonal.
    Assumes: synthetic 3x3 matrices."""
    from radia.fin_volume import require_energy_admissible

    assert require_energy_admissible(np.diag([1.0, 2.0, 3.0]))["cholesky"] is True
    dup = np.array([[1.0, 1.0, 0.0], [1.0, 1.0, 0.0], [0.0, 0.0, 2.0]])
    rec = require_energy_admissible(dup)
    assert rec["cholesky"] is False and rec["n_null"] == 1, rec
    bad = np.array([[1.0, 2.0, 0.0], [2.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    with pytest.raises(ValueError, match="negative energy.*witness branches"):
        require_energy_admissible(bad)
    with pytest.raises(ValueError, match="not symmetric"):
        require_energy_admissible(np.array([[1.0, 0.1], [0.0, 1.0]]))


def test_oriented_box_kernel_rejects_bad_input():
    """Guards: fail-loud contract of the opt-in kernel.
    Assumes: none."""
    from radia.peec_matrices import PEECBuilder

    for order in (-1, 33):
        with pytest.raises(Exception, match="order"):
            PEECBuilder(oriented_box_order=order)
    with pytest.raises(Exception, match="exclusive"):
        PEECBuilder(volume_order=8, oriented_box_order=8)
    b = PEECBuilder(oriented_box_order=8)
    n0, n1 = b.add_node_at(0.0, 0.0, 0.0), b.add_node_at(0.0, 0.0, 1e-3)
    b.add_connected_segment(n0, n1, 1e-4, 1e-4)
    with pytest.raises(Exception, match="parallel"):
        b.set_last_width_dir(0.0, 0.0, 1.0)
