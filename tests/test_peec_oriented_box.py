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
    rad about y."""
    out = []
    for x, r in ((0.0, np.eye(3)), (thick + gap, _rot([0, 1, 0], tilt))):
        c = np.array([x, 0, 0])
        out.append((c + r @ [0, 0, -1e-3], c + r @ [0, 0, 1e-3], 1e-3, thick, (0, 1, 0)))
    return out


@pytest.mark.parametrize("thick", [20e-6, 2e-6])
def test_thin_plates_agree_with_the_cartesian_rule_and_are_continuous(thick):
    """Guards: thin aligned plates from separation through contact to
    overlap match the independent Cartesian rule (order 128) within 1e-5;
    inside the SAT contact tolerance (about 2.2e-15 m here) the entry is
    constant to 1e-12 (only the geometry moves), and just beyond it (-1e-14 m, composite path) it moves
    by less than 1e-6; tilts bracketing the 1e-12 cross-axis cutoff (1e-13,
    1e-11 rad) give the same entry within 1e-12.
    Assumes: measured order-24 error at most 3.2e-7 on frozen 0024997
    (controls/diag_amend2). Not discriminating: both rules pass, and with the
    shared width axis the tilted cross axes coincide with face normals, so it
    does not show whether the SAT skip fires; the skip only selects the
    quadrature path, which no entry value here can reveal."""
    for gap in (-0.5 * thick, -0.1 * thick, 0.1 * thick):
        ref = _matrix(_plates(thick, gap), volume_order=128)[0, 1]
        got = _matrix(_plates(thick, gap), oriented_box_order=24)[0, 1]
        assert abs(got / ref - 1) <= 1e-5, (thick, gap, got, ref)
    at = [_matrix(_plates(thick, g), oriented_box_order=24)[0, 1] for g in (-2e-16, 0.0, 2e-16)]
    assert max(abs(a / at[1] - 1) for a in at) <= 1e-12, at
    beyond = _matrix(_plates(thick, -1e-14), oriented_box_order=24)[0, 1]
    assert abs(beyond / at[1] - 1) <= 1e-6, (beyond, at[1])
    t13, t11 = (_matrix(_plates(thick, -0.5 * thick, a), oriented_box_order=24)[0, 1] for a in (1e-13, 1e-11))
    assert abs(t13 / t11 - 1) <= 1e-12, (t13, t11)
    m = _matrix(_plates(thick, -0.5 * thick, 1e-7), oriented_box_order=24)
    ev = np.linalg.eigvalsh(m)
    assert ev[0] >= -4 * np.finfo(float).eps * ev[-1], ev



def test_thin_overlap_absolute_error_discriminates_the_composite_rule():
    """Guards: aligned 2 um plates overlapping by half their thickness give a
    mutual entry within 1e-6 L_ii (absolute) of the Cartesian rule at order
    256. Discriminating: composite rule 1.22e-7 L_ii, single-tensor rule
    (frozen 42f3032) 1.18e-6 L_ii; both pass the 1e-5 relative guard above.
    Assumes: Cartesian 128 -> 256 shifts this entry by 1.1e-7 relative."""
    bars = _plates(2e-6, -1e-6)
    ref = _matrix(bars, volume_order=256)
    got = _matrix(bars, oriented_box_order=24)
    assert abs(got[0, 1] - ref[0, 1]) <= 1e-6 * ref[0, 0], (got[0, 1], ref[0, 1], ref[0, 0])


_BIG = ((0.0, 0, -1e-3), (0.0, 0, 1e-3), 1e-3, 0.4e-3, (0, 1, 0))


def _small(zc):
    return ((0.05e-3, 0.1e-3, zc - 0.35e-3), (0.05e-3, 0.1e-3, zc + 0.35e-3), 0.3e-3, 0.1e-3, (0, 1, 0))


def test_partially_penetrating_bar_uses_contained_descendant_cells():
    """Guards: a small bar half inside a larger one (end face crossing; the
    composite rule's descendant cells inside the big bar take the contained
    path) matches the Cartesian rule (order 128) within 1e-6, in the
    identity and a rotated frame, bitwise independent of insertion order.
    Discriminating: the single-tensor rule (frozen 42f3032) is 3.0e-5 off,
    the composite rule 3.5e-8 (controls/diag_amend2).
    Assumes: Cartesian order 128 converged to well below 1e-6 here. Covers
    the composite route as a whole; it does not isolate the descendant
    containment shortcut, which only selects the quadrature path."""
    ref = _matrix([_BIG, _small(0.9e-3)], volume_order=128)
    for r in (np.eye(3), _rot([1, 2, 3], 0.7)):
        bars = [(r @ p0, r @ p1, w, h, r @ np.asarray(wd, float))
                for p0, p1, w, h, wd in (_BIG, _small(0.9e-3))]
        got = _matrix(bars, oriented_box_order=24)
        assert abs(got[0, 1] / ref[0, 1] - 1) <= 1e-6, (got[0, 1], ref[0, 1])
        swapped = _matrix(bars[::-1], oriented_box_order=24)
        assert swapped[0, 1] == got[0, 1] and swapped[0, 0] == got[1, 1]


def test_contained_top_level_box_matches_the_cartesian_rule():
    """Guards: a small bar wholly inside a larger one (top-level containment
    shortcut, no subdivision) matches the Cartesian rule within 2e-3 in the
    identity and a rotated frame and keeps |L_ij|^2 <= L_ii L_jj.
    Assumes: Cartesian order 64 converged for the nested pair; coarse guard,
    not a composite-rule discriminator."""
    small = ((0.05e-3, 0.1e-3, -0.3e-3), (0.05e-3, 0.1e-3, 0.4e-3), 0.3e-3, 0.1e-3, (0, 1, 0))
    ref = _matrix([_BIG, small], volume_order=64)
    for r in (np.eye(3), _rot([1, 2, 3], 0.7)):
        bars = [(r @ p0, r @ p1, w, h, r @ np.asarray(wd, float)) for p0, p1, w, h, wd in (_BIG, small)]
        got = _matrix(bars, oriented_box_order=24)
        assert np.allclose(got, ref, rtol=2e-3, atol=0), (got, ref)
        assert got[0, 1] ** 2 <= got[0, 0] * got[1, 1], got


def _builder(bars, order=24):
    from radia.peec_matrices import PEECBuilder

    b = PEECBuilder(oriented_box_order=order)
    for p0, p1, w, h, wd in bars:
        n0 = b.add_node_at(*map(float, p0))
        n1 = b.add_node_at(*map(float, p1))
        b.add_connected_segment(n0, n1, float(w), float(h))
        b.set_last_width_dir(*map(float, wd))
    return b


def _trace(bars, i=0, j=1, order=24):
    return _builder(bars, order).oriented_box_trace(i, j)


def test_sat_cross_axis_cutoff_is_bracketed_by_the_trace():
    """Guards: near-parallel half-overlapping plates reach the SAT edge-cross
    axes; a 1e-13 rad tilt (cross norm below the 1e-12 cutoff) skips the
    tilted pairs that a 1e-11 rad tilt evaluates, and the entry is the same
    within 1e-12 on both sides of the cutoff.
    Assumes: the diagnostic trace counts each skipped and evaluated cross
    axis of every CellStraddles call."""
    lo, hi = (_trace(_plates(2e-6, -1e-6, a)) for a in (1e-13, 1e-11))
    assert lo["composite"] and hi["composite"], (lo, hi)
    assert lo["cross_axes_skipped"] > hi["cross_axes_skipped"], (lo, hi)
    assert hi["cross_axes_evaluated"] > lo["cross_axes_evaluated"], (lo, hi)
    g = [_matrix(_plates(2e-6, -1e-6, a), oriented_box_order=24)[0, 1] for a in (1e-13, 1e-11)]
    assert abs(g[0] / g[1] - 1) <= 1e-12, g


def test_trace_reports_contained_descendants_after_subdivision():
    """Guards: the partially penetrating bar takes the composite route, and
    after subdivision descendant cells inside the big bar take the corner
    containment shortcut (contained_descendants > 0, none at the top level);
    the wholly contained bar takes the top-level shortcut without splitting;
    a separated pair never enters the composite route; and tracing on the
    same builder leaves its entries and assembled L bitwise unchanged and
    is symmetric in (i, j).
    Assumes: penetration geometry as in the composite-rule witness above."""
    pen = _trace([_BIG, _small(0.9e-3)])
    assert pen["composite"] and pen["splits"] > 0, pen
    assert pen["contained_descendants"] > 0 and pen["contained_top"] == 0, pen
    inside = ((0.05e-3, 0.1e-3, -0.3e-3), (0.05e-3, 0.1e-3, 0.4e-3), 0.3e-3, 0.1e-3, (0, 1, 0))
    top = _trace([_BIG, inside])
    assert not top["composite"] and top["contained_top"] == 1 and top["splits"] == 0, top
    far = _trace(_pair(0.6e-3))
    assert not far["composite"] and far["leaves"] == 0, far
    b = _builder([_BIG, _small(0.9e-3)])
    entries = [b.inductance_entry(i, j) for i in (0, 1) for j in (0, 1)]
    t01, t10 = b.oriented_box_trace(0, 1), b.oriented_box_trace(1, 0)
    assert t01 == t10 == pen, (t01, t10, pen)
    assert [b.inductance_entry(i, j) for i in (0, 1) for j in (0, 1)] == entries
    after = np.asarray(b.build_topology()["L"])
    assert np.array_equal(after, _matrix([_BIG, _small(0.9e-3)], oriented_box_order=24)), after
    assert [after[i, j] for i in (0, 1) for j in (0, 1)] == entries


def test_tilted_plate_contact_switches_the_route():
    """Guards: for 2 um plates tilted 1e-7 rad about y the composite route is
    taken just inside the analytic first contact (gap g_c - 1e-14 m) and the
    plain route at and just beyond it (g_c, g_c + 1e-14 m).
    Assumes: g_c = 1 mm sin(t) - (t_p / 2)(1 - cos(t)), the gap at which the
    tilted plate's nearest edge touches the first plate. Routing only: the
    order-24 entry moves by 1.7e-6 relative across this bracket (reviewer
    diagnostic) and no accuracy tolerance is adopted for it."""
    thick, tilt = 2e-6, 1e-7
    gc = 1e-3 * math.sin(tilt) - 0.5 * thick * (1 - math.cos(tilt))
    routes = [_trace(_plates(thick, gc + d, tilt))["composite"] for d in (-1e-14, 0.0, 1e-14)]
    assert routes == [True, False, False], (gc, routes)


def test_trace_rejects_bad_input():
    """Guards: fail-loud trace contract (kernel off, index out of range).
    Assumes: none."""
    from radia.peec_matrices import PEECBuilder

    b = PEECBuilder()
    n0, n1 = b.add_node_at(0.0, 0.0, 0.0), b.add_node_at(0.0, 0.0, 1e-3)
    b.add_connected_segment(n0, n1, 1e-4, 1e-4)
    with pytest.raises(Exception, match="oriented-box"):
        b.oriented_box_trace(0, 0)
    with pytest.raises(Exception, match="index"):
        _trace(_pair(0.6e-3), 0, 2)
