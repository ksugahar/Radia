"""Native axisymmetric volume PEEC kernel (rad_peec_axisym via radia.peec_axisym).

Guards: the ring partial inductance and loop-voltage solve of
``peec_matrices.axisym_ring_inductance`` / ``axisym_ring_loop_solve``: the far filament
mutual against Maxwell's formula with SciPy's complete elliptic integrals, the DC
limit of the circuit, a thin square ring against the geometric-mean-distance
inductance, symmetry and power balance, and loud refusal of invalid sections.
Assumes: structured meridian triangles; no FEM or meshing library.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

pm = pytest.importorskip("radia.peec_matrices")
if not hasattr(pm, "axisym_ring_inductance"):
    pytest.fail("radia.peec_matrices lacks axisym_ring_inductance: rebuild the native module")

from radia.peec_axisym import solve_axisym_ring_peec  # noqa: E402

MU0 = 4e-7*math.pi


def square_section(r0, z0, a, n):
    """n x n squares of side a/n, each split into two triangles."""
    xs = np.linspace(r0 - a/2, r0 + a/2, n + 1)
    zs = np.linspace(z0 - a/2, z0 + a/2, n + 1)
    pts = np.array([(x, z) for z in zs for x in xs])
    idx = lambda i, j: j*(n + 1) + i
    tris = []
    for j in range(n):
        for i in range(n):
            tris += [(idx(i, j), idx(i + 1, j), idx(i + 1, j + 1)), (idx(i, j), idx(i + 1, j + 1), idx(i, j + 1))]
    return pts, np.array(tris, dtype=np.int32)


def maxwell_mutual(r1, z1, r2, z2):
    from scipy.special import ellipe, ellipk
    m = 4*r1*r2/((r1 + r2)**2 + (z1 - z2)**2)
    k = math.sqrt(m)
    return MU0*math.sqrt(r1*r2)*((2/k - k)*ellipk(m) - (2/k)*ellipe(m))


def test_far_pair_is_the_filament_mutual_of_the_centroids():
    pts = np.array([[0.10, 0.0], [0.1001, 0.0], [0.10, 0.0001],
                    [0.13, 0.02], [0.1301, 0.02], [0.13, 0.0201]])
    M, area, c = pm.axisym_ring_inductance(pts, np.array([[0, 1, 2], [3, 4, 5]], dtype=np.int32))
    expected = maxwell_mutual(c[0, 0], c[0, 1], c[1, 0], c[1, 1])
    assert M[0, 1] == M[1, 0]
    assert abs(M[0, 1]/expected - 1) < 1e-12, (M[0, 1], expected)


def test_dc_limit_is_the_parallel_ring_resistance():
    pts, tris = square_section(0.1, 0.0, 2e-3, 4)
    out = solve_axisym_ring_peec(pts, tris, 1e-3, 5.8e7)
    r_cells = 2*math.pi*out["centroid_m"][:, 0]/(5.8e7*out["area_m2"])
    expected = 1.0/np.sum(1.0/r_cells)
    assert abs(out["R_ohm"]/expected - 1) < 1e-9, (out["R_ohm"], expected)
    assert out["power_balance"] < 1e-9, out["power_balance"]
    assert out["true_relative_residual"] < 1e-12, out["true_relative_residual"]


def test_thin_square_ring_matches_the_gmd_inductance():
    R, a = 0.1, 1e-3
    pts, tris = square_section(R, 0.0, a, 8)
    out = solve_axisym_ring_peec(pts, tris, 1.0, 5.8e7)       # skin depth 66 mm >> a: uniform current
    gmd = 0.447049*a                                          # square cross-section
    expected = MU0*R*(math.log(8*R/gmd) - 2.0)
    assert abs(out["L_H"]/expected - 1) < 1e-4, (out["L_H"], expected)


def test_matrix_is_symmetric_with_positive_self_terms():
    pts, tris = square_section(0.05, 0.01, 3e-3, 5)
    M, area, c = pm.axisym_ring_inductance(pts, tris)
    assert np.array_equal(M, M.T)
    assert np.all(np.diag(M) > np.abs(M - np.diag(np.diag(M))).max(axis=1)), "self term not dominant"


@pytest.mark.parametrize("name, mutate, match", [
    ("axis", lambda p, t: (p - np.array([p[:, 0].max(), 0.0]), t), "off the axis"),
    ("degenerate", lambda p, t: (p, np.vstack([t, [[0, 0, 1]]]).astype(np.int32)), "degenerate"),
    ("index", lambda p, t: (p, np.vstack([t, [[0, 1, len(p)]]]).astype(np.int32)), "references point"),
    ("nan", lambda p, t: (np.where(np.arange(len(p))[:, None] == 0, np.nan, p), t), "off the axis"),
])
def test_invalid_sections_are_refused(name, mutate, match):
    pts, tris = square_section(0.1, 0.0, 1e-3, 2)
    p, t = mutate(pts, tris)
    with pytest.raises(ValueError, match=match):
        pm.axisym_ring_inductance(p, t)


def test_adapter_refuses_bad_shapes_and_labels():
    pts, tris = square_section(0.1, 0.0, 1e-3, 2)
    with pytest.raises(ValueError, match="points"):
        solve_axisym_ring_peec(pts[:, :1], tris, 1e3, 5.8e7)
    with pytest.raises(ValueError, match="integer"):
        solve_axisym_ring_peec(pts, tris.astype(float), 1e3, 5.8e7)
    with pytest.raises(ValueError, match="region"):
        solve_axisym_ring_peec(pts, tris, 1e3, 5.8e7, region=["a"])


def test_reversed_winding_gives_the_same_matrix():
    pts, tris = square_section(0.2, 0.0, 1e-3, 3)
    flipped = tris[:, [0, 2, 1]].copy()
    mixed = tris.copy()
    mixed[::2] = flipped[::2]
    M, _, _ = pm.axisym_ring_inductance(pts, tris)
    for other in (flipped, mixed):
        Mo, _, _ = pm.axisym_ring_inductance(pts, other)
        assert np.abs(Mo - M).max() <= 1e-15*np.abs(M).max(), np.abs(Mo - M).max()
    assert np.all(np.diag(M) > 0)


def test_far_mutual_is_the_coaxial_dipole_limit():
    tri = np.array([[0.2, 0.0], [0.2001, 0.0], [0.2, 0.0001]])
    for dz in (10.0, 1e3, 1e4):
        M, _, c = pm.axisym_ring_inductance(np.vstack([tri, tri + [0, dz]]),
                                            np.array([[0, 1, 2], [3, 4, 5]], dtype=np.int32))
        r, d = c[0, 0], dz
        dipole = MU0*math.pi*r**4/(2*d**3)*(1 - 3*(r/d)**2)   # equal coaxial loops, to second order
        assert M[0, 1] > 0 and abs(M[0, 1]/dipole - 1) < 20*(r/d)**4 + 1e-12, (dz, M[0, 1], dipole)


def test_series_and_closed_form_agree_at_the_switch():
    # m = 4 r1 r2 / ((r1 + r2)^2 + dz^2) = 0.5 at dz = sqrt(4 r1 r2 / 0.5 - (r1 + r2)^2)
    r = 0.2
    dz0 = math.sqrt(8*r*r - 4*r*r)
    vals = []
    for dz in (dz0*(1 - 1e-9), dz0*(1 + 1e-9)):
        tri = np.array([[r, 0.0], [r + 1e-7, 0.0], [r, 1e-7]])
        M, _, _ = pm.axisym_ring_inductance(np.vstack([tri, tri + [0, dz]]),
                                            np.array([[0, 1, 2], [3, 4, 5]], dtype=np.int32))
        vals.append(M[0, 1])
    assert abs(vals[0]/vals[1] - 1) < 1e-8, vals


def test_loop_solve_refuses_a_nonsymmetric_matrix():
    M = np.array([[2e-6, 1e-6], [0.0, 3e-6]])
    with pytest.raises(ValueError, match="not symmetric"):
        pm.axisym_ring_loop_solve(M, np.array([1e-6, 2e-6]), np.array([[0.2, 0], [0.21, 0]]), 150e3, 5.8e7)
