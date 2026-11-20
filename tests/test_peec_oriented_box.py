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
    # Equal volumes: the outer/inner choice must not follow insertion order.
    swapped = _matrix(bars[::-1], oriented_box_order=8)
    assert swapped[0, 1] == got[0, 1], (swapped[0, 1], got[0, 1])


def test_near_coincident_centres_keep_the_full_order():
    """Guards: the quadrature count is clamped in floating point, so a
    centre distance far below the bar size (1e-150 m) keeps the maximum
    order instead of overflowing the integer conversion, and the entry
    equals the coincident-pair self term.
    Assumes: two identical 10 mm bars."""
    def pair(dx):
        return [((x, 0, 0), (x, 0, 10e-3), 1e-3, 0.2e-3, (1, 0, 0)) for x in (0.0, dx)]

    near = _matrix(pair(1e-150), oriented_box_order=24)
    assert np.isfinite(near).all(), near
    assert abs(near[0, 1] - near[0, 0]) <= 1e-9 * near[0, 0], near


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


def test_energy_gate_rejects_malformed_input_and_handles_zero_scale():
    """Guards: the gate is fail-closed on non-finite, empty, non-square or
    complex input, and a zero matrix reports finite diagnostics (all modes
    null) instead of NaN.
    Assumes: synthetic matrices."""
    from radia.fin_volume import require_energy_admissible

    for bad, match in ((np.diag([1.0, np.nan]), "non-finite"), (np.diag([1.0, np.inf]), "non-finite"),
                       (np.zeros((0, 0)), "non-empty"), (np.ones((2, 3)), "square"),
                       (np.eye(2) * (1 + 1j), "real")):
        with pytest.raises(ValueError, match=match):
            require_energy_admissible(bad)
    rec = require_energy_admissible(np.zeros((2, 2)))
    assert rec["n_null"] == 2 and rec["min_eig_over_max"] is None, rec


def _junction(gap, angle_deg):
    """Fin-like junction: a long longitudinal bar along z and a short, wide
    transverse plate leaving its side at angle_deg to the bar (90 = the
    straight-lattice orthogonal case), width axis shared along y."""
    th = math.radians(angle_deg)
    d_plate = np.array([math.sin(th), 0.0, math.cos(th)])
    bar = (np.array([0.0, 0, -1e-3]), np.array([0.0, 0, 1e-3]), 0.3e-3, 0.2e-3, (0, 1, 0))
    start = np.array([0.1e-3 + gap, 0.0, 0.0])
    plate = (start, start + 0.8e-3 * d_plate, 0.6e-3, 0.15e-3, (0, 1, 0))
    return [bar, plate]


def test_junction_entries_are_order_independent_of_insertion_and_converge():
    """Guards: for a skew transverse/longitudinal junction (80 degrees, the
    curved-lattice case) the entry matches a plain 6-D Gauss integral at
    near contact (10 um gap) within 3e-3, is bitwise independent of the
    insertion order, and for touching and overlapping support it converges
    with the quadrature order (orders 24 vs 32 within 1e-4) and stays PSD.
    Assumes: a 12-point composite Gauss per axis resolves the 10 um gap."""
    bars = _junction(10e-6, 80.0)
    got = _matrix(bars, oriented_box_order=24)
    swapped = _matrix(bars[::-1], oriented_box_order=24)
    assert got[0, 1] == swapped[0, 1] and got[0, 0] == swapped[1, 1]

    frames = []
    for p0, p1, w, h, wd in bars:
        d = (p1 - p0) / np.linalg.norm(p1 - p0)
        ew = np.asarray(wd, float) - np.dot(wd, d) * d
        ew /= np.linalg.norm(ew)
        frames.append((0.5 * (p0 + p1), [d, ew, np.cross(d, ew)],
                       (0.5 * np.linalg.norm(p1 - p0), w / 2, h / 2), d, w * h))
    pa, wa = _gauss_box(frames[0][0], frames[0][1], frames[0][2], 6, 2)
    pb, wb = _gauss_box(frames[1][0], frames[1][1], frames[1][2], 6, 2)
    integral = 0.0
    for chunk in range(0, len(pa), 512):
        dist = np.linalg.norm(pa[chunk:chunk + 512, None] - pb[None], axis=2)
        integral += float(wa[chunk:chunk + 512] @ (1 / dist) @ wb)
    ref = MU0_4PI * float(frames[0][3] @ frames[1][3]) * integral / (frames[0][4] * frames[1][4])
    assert abs(got[0, 1] - ref) <= 3e-3 * abs(ref), (got[0, 1], ref)

    for gap in (0.0, -0.05e-3):
        m24 = _matrix(_junction(gap, 80.0), oriented_box_order=24)
        m32 = _matrix(_junction(gap, 80.0), oriented_box_order=32)
        assert abs(m24[0, 1] - m32[0, 1]) <= 1e-4 * abs(m32[0, 1]), (gap, m24[0, 1], m32[0, 1])
        assert np.linalg.eigvalsh(m24)[0] > 0, (gap, m24)


def test_oblique_crossing_uses_the_composite_rule():
    """Guards: when the outer box crosses the inner box surface obliquely (a
    long fin bar under a wide transverse plate of the same layer, geometry of
    a real fin pair, rounded) the entry matches an independent reference
    within 2e-4 at order 24 (the single tensor rule was 3.5e-3 off), is
    bitwise independent of insertion order and flips sign exactly under
    current reversal.
    Assumes: reference -2.6883734243e-10 H from a NumPy integral of the
    closed-form sheet potential over the plate thickness, outer composite
    Gauss 4 points x 16 sub-intervals per axis (4 -> 16: 4.5e-4, 8 -> 16:
    1.9e-5)."""
    mm = 1e-3
    d_b = np.array([-0.8011, 0.5986, 0.0])
    d_b /= np.linalg.norm(d_b)
    w_b = np.array([-d_b[1], d_b[0], 0.0])
    c_b = np.array([-1.53, -2.88, 0.0]) * mm
    bar = (np.array([0.0, 4.037, 0]) * mm, np.array([0.0, -4.037, 0]) * mm, 1.206 * mm, 0.233 * mm, (1, 0, 0))
    plate = (c_b - 0.6435 * mm * d_b, c_b + 0.6435 * mm * d_b, 7.812 * mm, 0.426 * mm, w_b)
    ref = -2.6883734243e-10
    got = _matrix([bar, plate], oriented_box_order=24)
    assert abs(got[0, 1] / ref - 1) <= 2e-4, (got[0, 1], ref)
    swapped = _matrix([plate, bar], oriented_box_order=24)
    assert swapped[0, 1] == got[0, 1] and swapped[0, 0] == got[1, 1]
    reversed_plate = (plate[1], plate[0], *plate[2:])
    flipped = _matrix([bar, reversed_plate], oriented_box_order=24)
    assert flipped[0, 1] == -got[0, 1] and flipped[1, 1] == got[1, 1]


def _plates(thick, gap, tilt=0.0):
    """Two parallel 2 mm x 1 mm plates of thickness thick stacked along x
    with surface gap gap (negative = overlap); the second is tilted by tilt
    rad about y (near-parallel cross axes)."""
    out = []
    for x, r in ((0.0, np.eye(3)), (thick + gap, _rot([0, 1, 0], tilt))):
        c = np.array([x, 0, 0])
        out.append((c + r @ [0, 0, -1e-3], c + r @ [0, 0, 1e-3], 1e-3, thick, (0, 1, 0)))
    return out


@pytest.mark.parametrize("thick", [20e-6, 2e-6])
def test_thin_plates_through_contact_match_the_cartesian_rule(thick):
    """Guards: for thin, nearly parallel plates (tilt 1e-7 rad, so the SAT
    cross axes are near-degenerate) the mutual entry through separation,
    tolerance-scale gaps, contact and overlap matches the independent
    aligned Cartesian volume rule in relative (2e-3) and absolute
    (1e-6 L_ii) terms, varies continuously across contact and stays PSD.
    Assumes: tilt 1e-7 changes the exact entry by far less than 1e-6;
    Cartesian order 64 is converged for these plates."""
    gaps = [5 * thick, thick, 1e-12, 0.0, -1e-12, -0.5 * thick]
    prev = None
    for gap in gaps:
        ref = _matrix(_plates(thick, gap), volume_order=64)
        got = _matrix(_plates(thick, gap, tilt=1e-7), oriented_box_order=24)
        err = abs(got[0, 1] - ref[0, 1])
        assert err <= 2e-3 * abs(ref[0, 1]), (thick, gap, got[0, 1], ref[0, 1])
        assert err <= 1e-6 * ref[0, 0] or err <= 2e-3 * abs(ref[0, 1]), (thick, gap, err)
        ev = np.linalg.eigvalsh(got)
        assert ev[0] >= -4 * np.finfo(float).eps * ev[-1], (thick, gap, ev)
        if prev is not None and abs(gap) <= 1e-12:
            assert abs(got[0, 1] - prev) <= 1e-5 * abs(prev), (thick, gap, got[0, 1], prev)
        if abs(gap) <= 1e-12:
            prev = got[0, 1]


def test_nested_child_cell_matches_the_cartesian_rule():
    """Guards: a small bar wholly inside a larger one (the containment
    shortcut of the composite rule) matches the independent Cartesian
    volume rule within 2e-3, also after a rigid rotation, is bitwise
    independent of insertion order and keeps |L_ij| <= sqrt(L_ii L_jj).
    Assumes: Cartesian order 64 is converged for the nested pair."""
    big = ((0.0, 0, -1e-3), (0.0, 0, 1e-3), 1e-3, 0.4e-3, (0, 1, 0))
    small = ((0.05e-3, 0.1e-3, -0.3e-3), (0.05e-3, 0.1e-3, 0.4e-3), 0.3e-3, 0.1e-3, (0, 1, 0))
    ref = _matrix([big, small], volume_order=64)
    rot = _rot([1, 2, 3], 0.7)
    for r in (np.eye(3), rot):
        bars = [(r @ p0, r @ p1, w, h, r @ np.asarray(wd, float)) for p0, p1, w, h, wd in (big, small)]
        got = _matrix(bars, oriented_box_order=24)
        assert np.allclose(got, ref, rtol=2e-3, atol=0), (got, ref)
        swapped = _matrix(bars[::-1], oriented_box_order=24)
        assert swapped[0, 1] == got[0, 1] and swapped[0, 0] == got[1, 1]
        assert got[0, 1] ** 2 <= got[0, 0] * got[1, 1], got
