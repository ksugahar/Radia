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

Two kinds of Kelvin factor (don't confuse them!)
------------------------------------------------
This module exposes helpers for two DIFFERENT concepts:

1. **Solution pullback** (PEEC source evaluation): transform a vector
   potential A_phys defined on the physical exterior into the
   computational frame:
     A_comp(r') = (R/rho')^2 H A_phys(r_phys)        (1-form pullback)
     B_comp(r') = -(R/rho')^4 H B_phys(r_phys)       (2-form pullback)
   Functions: kelvin_pullback_vector, kelvin_pullback_B_pseudovector,
   kelvin_factor_scalar, A_s_at_obs_with_kelvin.

2. **Material modulation** (FEM bilinear form coefficient): the nu or
   mu that, applied in the transformed domain Omega', makes the FEM
   bilinear form equal the physical-domain energy. For 3D spherical
   (conformal) Kelvin (Nagamine CEFC 2026 eq. 9):
     nu_ext = (rho'/R)^2 * nu_0          [HCurl A-formulation]
     mu_ext = (R/rho')^2 * mu_0          [H1 Omega/H-formulation; reciprocal]
   Functions: kelvin_nu_factor_{3d,axisym,2d}_cf,
     kelvin_mu_factor_{3d,axisym,2d}_cf, build_material_cf.

Reference: H. Nagamine, T. Yamaguchi, K. Sugahara, "A Pullback-Based
Formulation of Kelvin Transformation in Electromagnetic Field Analysis,"
CEFC 2026 (Thessaloniki) id 350 (with Sugahara as co-author); see also
Sugahara 2022 IEEE TransMag 58(9) [ref [3] in Nagamine]. Canonical
declaration: examples/kelvin_transformation/CONVENTION.md.
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
    """Exact 1-form pullback of A under 3D sphere Kelvin inversion.

    See examples/kelvin_transformation/docs/pullback_derivation_3D.md
    for the full derivation. Summary:

    Kelvin map (origin at `center`, radius `R`):
        phi:  d -> d' = (R^2 / |d|^2) d,   d = r - c,  d' = r' - c.

    Forward Jacobian:
        J^j_i = dd'^j / dd^i  =  (rho'/R)^2 * H^j_i,    H = I - 2 n n^T,
                                                       n = d' / rho'
    Inverse Jacobian (involutive):
        (J^-1)^j_i = dd^j / dd'^i  =  (R/rho')^2 * H^j_i.

    For a 1-form omega, the pullback (line-integral preserving) is

        omega_comp_i (r') = (dd^j / dd'^i) * omega_phys_j (r_phys)
                         = (R/rho')^2 * H^j_i * omega_phys_j

    In matrix form (H symmetric):

        A_comp(r') = (R / rho')^2 * H * A_phys(r_phys)
                  = (R / rho')^2 * [A_phys - 2 (A_phys . n) n]

    The factor is (R/rho')^2 (NOT (rho'/R)^2): the inverse Jacobian has
    determinant (R/rho')^6, and a 1-form picks up one factor of the
    inverse Jacobian per index. The Householder reflection flips the
    radial component of A while preserving tangential components.

    Internal consistency: applying curl in the computational frame to
    this A_comp recovers the 2-form pullback of B (see
    kelvin_pullback_B_pseudovector). Verified numerically in
    test_kelvin_source.py and in pullback_derivation_3D.md sec 6.

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
    factor = (R / rho_prime) ** 2   # 1-form A pullback factor
    out = factor[:, None] * A_refl
    return out[0] if single else out


def kelvin_pullback_B_pseudovector(B_phys, r_prime, center, R):
    """3D Kelvin pullback of B as a 2-form pseudovector.

        B_comp(r') = -(R/rho')^4 * H * B_phys(r_phys)

    Derived in pullback_derivation_3D.md section 5: the 2-form
    transforms with two inverse-Jacobian indices, picking up
    (R/rho')^4 in magnitude. Hodge-dual to a vector adds a sign of
    det(H) = -1 from the Levi-Civita identity. So B is a pseudovector
    that flips its sign through the Householder reflection AND picks
    up an overall minus from the Hodge identity.

    Args:
        B_phys: (N, 3) physical-frame B at r_phys = Kelvin_inv(r').
        r_prime: (N, 3) computational coordinates.
        center, R: Kelvin sphere center, radius.

    Returns:
        (N, 3) pulled-back B'.
    """
    r_prime = np.asarray(r_prime, dtype=float)
    B_phys = np.asarray(B_phys)
    c = np.asarray(center, dtype=float)
    single = r_prime.ndim == 1
    if single:
        r_prime = r_prime.reshape(1, 3)
        B_phys = B_phys.reshape(1, 3)
    d_prime = r_prime - c
    rho_prime = np.sqrt(np.sum(d_prime * d_prime, axis=1))
    rho_prime = np.where(rho_prime < 1e-30, 1e-30, rho_prime)
    n = d_prime / rho_prime[:, None]
    B_dot_n = np.sum(B_phys * n, axis=1)
    B_refl = B_phys - 2.0 * B_dot_n[:, None] * n
    factor = -(R / rho_prime) ** 4   # 2-form pseudovector factor
    out = factor[:, None] * B_refl
    return out[0] if single else out


def kelvin_factor_scalar(r_prime, center, R):
    """Scalar magnitude factor for 1-form A pullback: (R/rho')^2.

    Drops the Householder reflection from the full 1-form pullback
    (kelvin_pullback_vector). Exact when A_phys is purely tangential
    at the evaluation point (n . A_phys = 0); in particular for the
    azimuthal A_phys of a uniform B_0 z_hat background, Householder
    is the identity and only this factor remains.

    See pullback_derivation_3D.md section 4 for the derivation: the
    inverse Kelvin Jacobian (dr/dr') has determinant (R/rho')^6 in 3D
    and a 1-form picks up (R/rho')^2 per index.

    NOTE on history: an earlier (2026-04-15) "fix" replaced this with
    (rho'/R)^2 based on an empirical observation that the (R/rho')^2
    formula appeared to over-count energy in the Kelvin exterior. That
    over-count was real but came from the missing Householder, NOT
    from a wrong scalar factor. The proper resolution is to use the
    full 1-form pullback (kelvin_pullback_vector); the scalar form is
    only correct as the tangential limit.
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


# ---------- FEM (NGSolve) material CoefficientFunction factory --------
#
# Helpers for the canonical Kelvin material modulation per
#   H. Nagamine, T. Yamaguchi, K. Sugahara,
#   "A Pullback-Based Formulation of Kelvin Transformation in
#    Electromagnetic Field Analysis," CEFC 2026 id 350
# (with Sugahara as co-author); see also Sugahara 2022 IEEE TransMag
# 58(9) [ref [3] in Nagamine]. Full declaration in
#   examples/kelvin_transformation/CONVENTION.md
# and derivation in
#   examples/kelvin_transformation/docs/pullback_derivation_3D.md sec 8.
#
# 3D spherical (conformal) Kelvin:
#     nu' = (rho'/R)^2 * nu_0          [HCurl A-formulation]
#     mu' = (R/rho')^2 * mu_0          [H1 Omega / H-formulation]
# (pointwise reciprocals; mu * nu = 1)
#
# These factors are the MATERIAL modulation entering the FEM bilinear
# form. They are DIFFERENT from the solution pullback factors (R/rho')^2
# on A and -(R/rho')^4 on B used by kelvin_pullback_vector and
# kelvin_pullback_B_pseudovector above.


def kelvin_nu_factor_3d_cf(center, R, coords=None):
    """3D spherical Kelvin nu factor (rho'/R)^2 as NGSolve CF.

    Multiply by nu_0 to produce the Kelvin exterior-domain reluctivity:
        nu_ext = nu_0 * kelvin_nu_factor_3d_cf(...)

    Vanishes at rho'=0 (image of physical infinity), equals nu_0 at
    rho'=R (continuous with air domain). Reference: Nagamine CEFC 2026
    eq. (9).

    Args:
        center: (cx, cy, cz) -- Kelvin sphere center in world coords.
        R: Kelvin sphere radius.
        coords: optional 3-tuple of NGSolve CFs (default: (x, y, z)).

    Returns:
        NGSolve CoefficientFunction.
    """
    from ngsolve import x, y, z, sqrt, IfPos
    if coords is None:
        coords = (x, y, z)
    cx, cy, cz = center
    rho2 = (coords[0] - cx) ** 2 + (coords[1] - cy) ** 2 + (coords[2] - cz) ** 2
    rho_prime = sqrt(rho2)
    rho_safe = IfPos(rho_prime - 1e-10, rho_prime, 1e-10)
    return (rho_safe / R) ** 2


def kelvin_mu_factor_3d_cf(center, R, coords=None):
    """3D spherical Kelvin mu factor (R/rho')^2 as NGSolve CF.

    Multiply by mu_0 to produce the Kelvin exterior-domain permeability:
        mu_ext = mu_0 * kelvin_mu_factor_3d_cf(...)

    Diverges at rho'=0 (mu_r -> infinity at image of infinity), equals
    mu_0 at rho'=R. Reciprocal of kelvin_nu_factor_3d_cf.
    """
    from ngsolve import x, y, z, sqrt, IfPos
    if coords is None:
        coords = (x, y, z)
    cx, cy, cz = center
    rho2 = (coords[0] - cx) ** 2 + (coords[1] - cy) ** 2 + (coords[2] - cz) ** 2
    rho_prime = sqrt(rho2)
    rho_safe = IfPos(rho_prime - 1e-10, rho_prime, 1e-10)
    return (R / rho_safe) ** 2


def kelvin_nu_factor_axisym_cf(z_offset, R, r_coord=None, z_coord=None):
    """Axisym (r,z) Kelvin nu factor (rho'/R)^2 with Z-offset.

    rho' = sqrt(r^2 + (z - z_offset)^2)  (3D spherical distance in the
    meridional plane; this is 3D sphere Kelvin viewed axisymmetrically).

    NGSolve axisym convention: r = x, z = y. Pass explicit r_coord /
    z_coord when your mesh uses a different embedding.

    Multiply by nu_0 for A-formulation nu_ext.
    """
    from ngsolve import x, y, sqrt, IfPos
    if r_coord is None:
        r_coord = x
    if z_coord is None:
        z_coord = y
    z_local = z_coord - z_offset
    rho_prime = sqrt(r_coord ** 2 + z_local ** 2)
    rho_safe = IfPos(rho_prime - 1e-10, rho_prime, 1e-10)
    return (rho_safe / R) ** 2


def kelvin_mu_factor_axisym_cf(z_offset, R, r_coord=None, z_coord=None):
    """Axisym (r,z) Kelvin mu factor (R/rho')^2. Reciprocal of nu factor."""
    from ngsolve import x, y, sqrt, IfPos
    if r_coord is None:
        r_coord = x
    if z_coord is None:
        z_coord = y
    z_local = z_coord - z_offset
    rho_prime = sqrt(r_coord ** 2 + z_local ** 2)
    rho_safe = IfPos(rho_prime - 1e-10, rho_prime, 1e-10)
    return (R / rho_safe) ** 2


def kelvin_nu_factor_2d_cf(offset, R, x_coord=None, y_coord=None):
    """2D Cartesian Kelvin nu factor (rho'/R)^2 with (x, y) offset.

    rho' = sqrt((x - x_off)^2 + (y - y_off)^2)

    For 2D CARTESIAN conformal Kelvin. For 2D CYLINDRICAL (z-axis)
    Kelvin, the modulation is anisotropic (Nagamine eq. 12):
    nu' = diag(1, 1, (rho'/R)^4) nu -- use a tensor CF rather than
    this scalar helper.
    """
    from ngsolve import x, y, sqrt, IfPos
    if x_coord is None:
        x_coord = x
    if y_coord is None:
        y_coord = y
    xo, yo = offset
    rho_prime = sqrt((x_coord - xo) ** 2 + (y_coord - yo) ** 2)
    rho_safe = IfPos(rho_prime - 1e-10, rho_prime, 1e-10)
    return (rho_safe / R) ** 2


def kelvin_mu_factor_2d_cf(offset, R, x_coord=None, y_coord=None):
    """2D Cartesian Kelvin mu factor (R/rho')^2. Reciprocal of nu factor."""
    from ngsolve import x, y, sqrt, IfPos
    if x_coord is None:
        x_coord = x
    if y_coord is None:
        y_coord = y
    xo, yo = offset
    rho_prime = sqrt((x_coord - xo) ** 2 + (y_coord - yo) ** 2)
    rho_safe = IfPos(rho_prime - 1e-10, rho_prime, 1e-10)
    return (R / rho_safe) ** 2


def build_material_cf(mesh, base_value, kelvin_factor_cf, *,
                      overrides=None,
                      outer_keyword="outer"):
    """Build a mesh-material-indexed CoefficientFunction.

    For each material in ``mesh.GetMaterials()``:
      - if present in ``overrides`` dict, use overrides[material]
      - elif name contains ``outer_keyword`` (case-insensitive),
        use ``base_value * kelvin_factor_cf``
      - otherwise use ``base_value``

    Typical usage:

        # A-formulation axisym with Z-offset:
        nu_cf = build_material_cf(
            mesh, nu0,
            kelvin_nu_factor_axisym_cf(z_offset, a),
            overrides={"magnetic": nu0 / mu_r},
        )

        # Omega-formulation 3D sphere:
        mu_cf = build_material_cf(
            mesh, mu0,
            kelvin_mu_factor_3d_cf(center, R),
            overrides={"magnetic": mu_r * mu0},
        )

    Args:
        mesh: NGSolve Mesh.
        base_value: float -- baseline material (nu_0 or mu_0).
        kelvin_factor_cf: CF -- pick the matching helper:
            nu-direction: kelvin_nu_factor_{3d,axisym,2d}_cf
            mu-direction: kelvin_mu_factor_{3d,axisym,2d}_cf
        overrides: {material_name: value_or_cf} explicit overrides.
        outer_keyword: substring to match Kelvin outer-domain name.

    Returns:
        NGSolve CoefficientFunction indexed by mesh.GetMaterials().
    """
    from ngsolve import CoefficientFunction
    overrides = overrides or {}
    values = []
    for mat in mesh.GetMaterials():
        if mat in overrides:
            values.append(overrides[mat])
        elif outer_keyword in mat.lower():
            values.append(base_value * kelvin_factor_cf)
        else:
            values.append(base_value)
    return CoefficientFunction(values)
