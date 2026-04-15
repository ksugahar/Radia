"""kelvin_source.py

Kelvin-transformed source field evaluation for reduced A-formulation FEM
with PEEC filament coils.

Scope (2026-04-15)
------------------
* Source (PEEC filament coil) is assumed to live in the INNER PHYSICAL
  domain (not geometrically inside the Kelvin exterior domain). This covers the
  standard IH setup where R_coil + a_coil <= R_air.
* Phase 1: scalar Kelvin factor (uniform-H_s analogue).
* Phase 2: correct 1-form pullback Jacobian for the Kelvin coordinate
  map. Phase 1 is a ``factor_mode='scalar'`` approximation; Phase 2 is
  ``factor_mode='pullback'`` (exact).
* Phase 3: source filaments that extend into the Kelvin exterior domain -- NOT
  supported here; document_user_discussed_skip_2026-04-15.

Geometry conventions
--------------------
Kelvin map (3D sphere) centered at ``c`` with radius ``R``:

    phi: r -> r' = c + (R**2 / |r - c|**2) * (r - c)

The map is involutive: phi(phi(r)) = r. Points with |r-c| > R are the
physical exterior and are mapped to the Kelvin exterior domain (|r'-c| < R).

A-formulation source (Biot-Savart vector potential)
---------------------------------------------------
For a filament bundle with currents I_k carried by polyline path_k:

    A_s(r) = (mu_0 / 4 pi) * sum_k I_k * integral_{path_k} dl' / |r - r'|

Evaluated via Gauss-Legendre quadrature on each polyline segment.
"""

from __future__ import annotations

import math

import numpy as np

MU_0 = 4e-7 * math.pi


# ---------- Kelvin coordinate map ----------------------------------------


def kelvin_map_3d(points, center, R):
    """Apply 3D sphere Kelvin inversion ``r' = c + R^2 (r-c) / |r-c|^2``.

    Args:
        points: (N, 3) or (3,) array of world coordinates.
        center: (3,) Kelvin sphere center.
        R: Kelvin sphere radius.

    Returns:
        Same-shape array of transformed coordinates.
    """
    p = np.asarray(points, dtype=float)
    c = np.asarray(center, dtype=float)
    single = p.ndim == 1
    if single:
        p = p.reshape(1, 3)
    d = p - c
    r2 = np.sum(d * d, axis=1)
    # Protect against singularity at r=0 (center of Kelvin sphere maps to
    # physical infinity, undefined for a finite-distance source).
    r2 = np.where(r2 < 1e-30, 1e-30, r2)
    out = c + (R * R / r2)[:, None] * d
    return out[0] if single else out


def is_in_kelvin_exterior_domain(points, center, R):
    """True where |p - center| < R (within the Kelvin exterior domain)."""
    p = np.asarray(points, dtype=float)
    c = np.asarray(center, dtype=float)
    return np.sum((p - c) ** 2, axis=-1) < R * R


# ---------- Kelvin factor (Phase 1 scalar, Phase 2 pullback) --------------


def kelvin_pullback_vector(A_phys, r_prime, center, R):
    """Phase 2: exact 1-form pullback of a vector A under 3D sphere Kelvin.

    Kelvin map in d-coords (d = r - center):
        phi: d -> d' = (R^2 / rho^2) * d,    rho = |d|, rho' = R^2/rho

    Jacobian at d:
        J(d) = (R^2 / rho^2) * (I - 2 n n^T),    n = d / rho    (Householder)

    Since phi is involutive, J(d)^{-1} = (rho^2/R^2)(I - 2 n n^T), and in
    exterior-domain coordinates where rho = R^2/rho' and n = n' (radial direction
    preserved):

        A'(r') = J^{-1}(r) . A(r)
               = (R / rho')^2 * (I - 2 n' n'^T) . A(r_phys)
               = (R / rho')^2 * [A_phys - 2 (A_phys . n') n']

    So the exact 1-form transformation is the scalar Phase-1 factor
    (R/rho')^2 composed with a Householder reflection that flips the
    RADIAL component of A. For A fields with no radial component (e.g.
    the azimuthal A from a uniform B_0 z_hat background) this reduces
    to the scalar factor; for localized Biot-Savart sources it does
    not.

    Args:
        A_phys: (N, 3) physical-frame A evaluated at r_phys = Kelvin_inv(r').
        r_prime: (N, 3) exterior-frame coordinates.
        center:  (3,) Kelvin sphere center.
        R: Kelvin sphere radius.

    Returns:
        (N, 3) pulled-back A', same dtype as A_phys.
    """
    r_prime = np.asarray(r_prime, dtype=float)
    A_phys = np.asarray(A_phys)
    c = np.asarray(center, dtype=float)
    single = r_prime.ndim == 1
    if single:
        r_prime = r_prime.reshape(1, 3)
        A_phys = A_phys.reshape(1, 3)
    d_prime = r_prime - c
    rho_prime = np.sqrt(np.sum(d_prime * d_prime, axis=1))
    rho_prime = np.where(rho_prime < 1e-30, 1e-30, rho_prime)
    n = d_prime / rho_prime[:, None]
    A_dot_n = np.sum(A_phys * n, axis=1)
    A_refl = A_phys - 2.0 * A_dot_n[:, None] * n
    factor = (R / rho_prime) ** 2
    out = factor[:, None] * A_refl
    return out[0] if single else out


def kelvin_factor_scalar(r_prime, center, R):
    """Phase 1: uniform-H_s-analogue scalar factor (R/rho')^2.

    Extracted from the H-formulation docs section 4.3:
        H'_s = -H_0 * (R / rho')^2 * z_hat'

    For Biot-Savart A_s with non-uniform spatial dependence, this is an
    APPROXIMATION. See ``kelvin_factor_pullback`` for the exact Jacobian.

    Args:
        r_prime: (N, 3) Kelvin exterior domain coordinates.
        center: Kelvin sphere center.
        R: Kelvin sphere radius.

    Returns:
        (N,) array of scalar factors. Magnitude only; sign (-1 for
        H_s direction flip) is applied by the caller when converting
        physical H_s to Kelvin exterior domain H'_s.
    """
    p = np.asarray(r_prime, dtype=float)
    c = np.asarray(center, dtype=float)
    rho_prime = np.sqrt(np.sum((p - c) ** 2, axis=-1))
    rho_prime = np.where(rho_prime < 1e-30, 1e-30, rho_prime)
    return (R / rho_prime) ** 2


# ---------- Biot-Savart vector potential from filaments -------------------


def biot_savart_A_at_points(obs_points, filament_paths, currents,
                             n_quad=8, chunk_size=2048):
    """Vector potential A(r) at obs points from a filament bundle.

        A(r) = (mu_0 / 4 pi) * sum_k I_k * integral_{path_k} dl' / |r - r'|

    Each filament is a polyline given as a list of ``(p1, p2)`` segment
    endpoint tuples. Segments are integrated via Gauss-Legendre
    quadrature. Per-filament currents may be real or complex.

    Args:
        obs_points: (N, 3) observation coordinates.
        filament_paths: list of K filaments, each a list of segment
            ``(p1, p2)`` tuples or a (n_seg, 2, 3) array.
        currents: length-K complex (or real) currents.
        n_quad: Gauss-Legendre points per segment (default 8).
        chunk_size: max rows of the (N_quad, N_obs) distance matrix
            per chunk (memory cap).

    Returns:
        (N, 3) complex (or real) array.
    """
    obs = np.asarray(obs_points, dtype=float)
    if obs.ndim == 1:
        obs = obs.reshape(1, 3)
    N_obs = obs.shape[0]

    # Flatten (segment, current) pairs.
    segs = []
    curr = []
    for path, Ik in zip(filament_paths, currents):
        Ikc = complex(Ik)
        for p1, p2 in path:
            segs.append((np.asarray(p1, float), np.asarray(p2, float)))
            curr.append(Ikc)
    if not segs:
        return np.zeros((N_obs, 3), dtype=complex)

    segs = np.array(segs)                   # (M, 2, 3)
    curr = np.array(curr, dtype=complex)    # (M,)
    M = segs.shape[0]

    t_gl, w_gl = np.polynomial.legendre.leggauss(n_quad)
    t01 = 0.5 * (t_gl + 1.0)                # [0, 1]
    w01 = 0.5 * w_gl

    # Quadrature points along each segment + tangent vectors.
    # q_pts[i, q] = segs[i, 0] + t01[q] * (segs[i, 1] - segs[i, 0])
    dl = segs[:, 1, :] - segs[:, 0, :]      # (M, 3)
    seg_len = np.linalg.norm(dl, axis=1)    # (M,)
    dl_hat = dl / np.where(seg_len[:, None] > 1e-30, seg_len[:, None], 1.0)
    q_pts = segs[:, 0:1, :] + t01[None, :, None] * dl[:, None, :]   # (M, Q, 3)

    # Accumulate A in chunks over M (quadrature source side).
    is_complex = np.iscomplexobj(curr)
    A = np.zeros((N_obs, 3), dtype=complex if is_complex else float)

    for s in range(0, M, chunk_size):
        e = min(s + chunk_size, M)
        mc = e - s
        qs = q_pts[s:e]               # (mc, Q, 3)
        # Distances (N_obs, mc, Q)
        dx = obs[:, None, None, 0] - qs[None, :, :, 0]
        dy = obs[:, None, None, 1] - qs[None, :, :, 1]
        dz = obs[:, None, None, 2] - qs[None, :, :, 2]
        dist = np.sqrt(dx * dx + dy * dy + dz * dz)
        dist = np.where(dist < 1e-14, 1e-14, dist)
        kernel = 1.0 / dist           # (N_obs, mc, Q)
        # Segment contribution to A:
        #   (mu_0/4pi) * I_k * seg_len_k * sum_q w_q * (dl_hat_k / |r - q|)
        weight = (MU_0 / (4 * math.pi)) * w01[None, None, :]     # (1, 1, Q)
        # sum over Q -> scalar per (obs, seg); then multiply by dl_hat * I * seg_len
        scalar = np.sum(weight * kernel, axis=-1)                # (N_obs, mc)
        I_len = (curr[s:e] * seg_len[s:e])                       # (mc,)
        # A_contrib[n, c] along direction dl_hat[c] is scalar[n, c] * I_len[c] * dl_hat[c]
        contrib = scalar[:, :, None] * (I_len[:, None] * dl_hat[s:e])  # (N_obs, mc, 3)
        A += np.sum(contrib, axis=1)

    if not is_complex:
        A = A.real
    return A


# ---------- Convenience: A_s at obs with Kelvin handling ------------------


def A_s_at_obs_with_kelvin(obs_points, filament_paths, currents,
                            kelvin_center=None, R_kelvin=None,
                            factor_mode='scalar', n_quad=8):
    """Evaluate A_s at observation points, handling Kelvin exterior domain.

    For obs points in the INNER physical domain: Biot-Savart directly.
    For obs points in the Kelvin exterior domain (|r'-c| < R): map r' -> physical
    r = c + R^2 (r'-c)/|r'-c|^2, evaluate Biot-Savart at r, then apply
    Kelvin factor.

    Args:
        obs_points: (N, 3)
        filament_paths, currents: as in biot_savart_A_at_points.
        kelvin_center, R_kelvin: if None, Kelvin handling is disabled
            (Biot-Savart everywhere). If both given, Kelvin exterior domain quad
            points are remapped.
        factor_mode: 'scalar' (Phase 1 approximation, (R/rho')**2 on
            the vector magnitude) or 'pullback' (Phase 2, exact
            1-form Jacobian -- to be implemented).
        n_quad: quadrature per segment.

    Returns:
        (N, 3) array of A values matching input dtype.
    """
    obs = np.asarray(obs_points, dtype=float)
    if obs.ndim == 1:
        obs = obs.reshape(1, 3)

    if kelvin_center is None or R_kelvin is None:
        return biot_savart_A_at_points(obs, filament_paths, currents,
                                       n_quad=n_quad)

    c = np.asarray(kelvin_center, dtype=float)
    in_kelvin_ext = is_in_kelvin_exterior_domain(obs, c, R_kelvin)

    A = np.zeros((obs.shape[0], 3), dtype=complex)
    # Inner points: direct evaluation.
    if np.any(~in_kelvin_ext):
        A_inner = biot_savart_A_at_points(
            obs[~in_kelvin_ext], filament_paths, currents, n_quad=n_quad)
        A[~in_kelvin_ext] = A_inner

    # Kelvin exterior domain points: map + evaluate + factor.
    if np.any(in_kelvin_ext):
        r_prime = obs[in_kelvin_ext]
        r_phys = kelvin_map_3d(r_prime, c, R_kelvin)
        A_phys = biot_savart_A_at_points(
            r_phys, filament_paths, currents, n_quad=n_quad)
        if factor_mode == 'scalar':
            f = kelvin_factor_scalar(r_prime, c, R_kelvin)   # (N_ext,)
            A[in_kelvin_ext] = A_phys * f[:, None]
        elif factor_mode == 'pullback':
            A[in_kelvin_ext] = kelvin_pullback_vector(
                A_phys, r_prime, c, R_kelvin)
        else:
            raise ValueError(f"Unknown factor_mode: {factor_mode!r}")

    if not np.iscomplexobj(np.asarray(list(currents))):
        A = A.real
    return A
