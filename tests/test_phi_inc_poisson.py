"""Golden tests: compute_phi_inc_surface_poisson (surface-Poisson phi_inc).

Pure numpy/scipy (no NGSolve): an icosphere in a uniform incident field
has the closed-form potential psi = -H0 . r (H = -grad psi), which the
Laplace-Beltrami projection must recover; a rotational (Killing) field
has NO gradient part, so the fail-loud gate must fire.
"""
from __future__ import annotations

import numpy as np
import pytest

from radia.bem_sibc_solver import compute_phi_inc_surface_poisson


def _icosphere(n_subdiv=2):
    """Unit icosphere: 12-vertex icosahedron + midpoint subdivision."""
    t = (1.0 + np.sqrt(5.0)) / 2.0
    verts = np.array([
        [-1, t, 0], [1, t, 0], [-1, -t, 0], [1, -t, 0],
        [0, -1, t], [0, 1, t], [0, -1, -t], [0, 1, -t],
        [t, 0, -1], [t, 0, 1], [-t, 0, -1], [-t, 0, 1],
    ], dtype=float)
    verts /= np.linalg.norm(verts, axis=1)[:, None]
    faces = [
        (0, 11, 5), (0, 5, 1), (0, 1, 7), (0, 7, 10), (0, 10, 11),
        (1, 5, 9), (5, 11, 4), (11, 10, 2), (10, 7, 6), (7, 1, 8),
        (3, 9, 4), (3, 4, 2), (3, 2, 6), (3, 6, 8), (3, 8, 9),
        (4, 9, 5), (2, 4, 11), (6, 2, 10), (8, 6, 7), (9, 8, 1),
    ]
    verts = list(map(tuple, verts))
    for _ in range(n_subdiv):
        cache = {}

        def midpoint(a, b):
            key = (min(a, b), max(a, b))
            if key not in cache:
                p = np.asarray(verts[a]) + np.asarray(verts[b])
                p /= np.linalg.norm(p)
                verts.append(tuple(p))
                cache[key] = len(verts) - 1
            return cache[key]

        new_faces = []
        for (a, b, c) in faces:
            ab, bc, ca = midpoint(a, b), midpoint(b, c), midpoint(c, a)
            new_faces += [(a, ab, ca), (b, bc, ab), (c, ca, bc),
                          (ab, bc, ca)]
        faces = new_faces
    return np.asarray(verts), np.asarray(faces, dtype=np.int64)


def test_uniform_field_recovers_linear_potential():
    """H = H0 z_hat  ->  psi = -H0 z (mean-zero on the sphere)."""
    pts, tris = _icosphere(2)          # 162 verts / 320 tris
    H0 = 123.0
    H = np.zeros((len(pts), 3))
    H[:, 2] = H0
    psi, resid = compute_phi_inc_surface_poisson(pts, tris, H)

    psi_exact = -H0 * pts[:, 2]
    rel = (np.linalg.norm(psi.real - psi_exact)
           / np.linalg.norm(psi_exact))
    assert rel < 0.02, f"psi vs -H0*z rel error {rel:.3e}"
    assert abs(np.mean(psi)) < 1e-9 * H0, "mean-zero gauge violated"
    assert np.linalg.norm(psi.imag) < 1e-9 * H0
    # faceting-limited consistency on a 2-subdiv icosphere
    assert resid < 0.10, f"grad-consistency residual {resid:.3e}"


def test_complex_field_scales_linearly():
    """Complex H phasor -> the same complex factor on psi (linearity)."""
    pts, tris = _icosphere(1)
    fac = 2.0 - 3.0j
    H = np.zeros((len(pts), 3), dtype=complex)
    H[:, 2] = fac
    psi, _ = compute_phi_inc_surface_poisson(pts, tris, H)
    H1 = np.zeros((len(pts), 3))
    H1[:, 2] = 1.0
    psi1, _ = compute_phi_inc_surface_poisson(pts, tris, H1)
    assert np.allclose(psi, fac * psi1, atol=1e-12)


def test_rotational_field_fails_loud():
    """A Killing (rotation) field has zero gradient part: the gate fires."""
    pts, tris = _icosphere(2)
    H = np.stack([-pts[:, 1], pts[:, 0], np.zeros(len(pts))], axis=1)
    with pytest.raises(ValueError, match="not a surface gradient"):
        compute_phi_inc_surface_poisson(pts, tris, H,
                                        max_grad_residual=0.10)
    # without the gate it returns the (near-total) residual for diagnosis
    psi, resid = compute_phi_inc_surface_poisson(pts, tris, H)
    assert resid > 0.9, f"rotational field residual {resid:.3f} (expect ~1)"
    assert np.linalg.norm(psi) < 0.05 * np.linalg.norm(pts[:, 2])


def test_winding_invariance():
    """Flipping triangle winding must not change psi (unsigned projector)."""
    pts, tris = _icosphere(1)
    H = np.zeros((len(pts), 3))
    H[:, 2] = 1.0
    psi_a, _ = compute_phi_inc_surface_poisson(pts, tris, H)
    tris_flipped = tris[:, [0, 2, 1]].copy()
    psi_b, _ = compute_phi_inc_surface_poisson(pts, tris_flipped, H)
    assert np.allclose(psi_a, psi_b, atol=1e-12)


def test_h_shape_mismatch_raises():
    pts, tris = _icosphere(1)
    with pytest.raises(ValueError, match="every surface vertex"):
        compute_phi_inc_surface_poisson(pts, tris, np.zeros((5, 3)))


@pytest.mark.parametrize("bad_value", [np.nan, np.inf, -np.inf])
def test_nonfinite_h_fails_loud(bad_value):
    pts, tris = _icosphere(1)
    H = np.zeros((len(pts), 3), dtype=complex)
    H[0, 0] = bad_value
    with pytest.raises(ValueError, match="only finite values"):
        compute_phi_inc_surface_poisson(pts, tris, H,
                                        max_grad_residual=0.10)
