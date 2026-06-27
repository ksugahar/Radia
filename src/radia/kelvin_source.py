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
   bilinear form equal the physical-domain energy.

   3D spherical (conformal) Kelvin (Nagamine CEFC 2026 eq. 9). The
   axisymmetric (r, z) case is the SAME basic formula applied in the
   meridional plane (the 3D Kelvin sphere viewed axisymmetrically):
     nu_ext = (rho'/R)^2 * nu_0          [HCurl A-formulation]
     mu_ext = (R/rho')^2 * mu_0          [H1 Omega/H-formulation; reciprocal]
   Helpers: kelvin_nu_factor_{3d,axisym}_cf,
     kelvin_mu_factor_{3d,axisym}_cf, build_material_cf.

   2D cylindrical (non-conformal; Nagamine eq. 12) -- ANISOTROPIC:
     nu' = diag(1, 1, (rho'/R)^4) * nu       (only nu_zz modulated)
     mu' = diag(1, 1, (R/rho')^4) * mu       (only mu_zz modulated)
   Which slot enters the bilinear form depends on which B / H
   components the problem actually uses -- not on whether it's
   "A-form" vs "Omega-form":
   - In-plane (Hx, Hy) or in-plane B (e.g. 2D A-form with
     A = A_z z_hat, so curl(A) is in-plane): the z-slot is NEVER
     contracted; the effective Kelvin factor is **1 (identity)**.
     This is the COMMON case. Helper: kelvin_factor_2d_inplane_cf().
   - Axial B_z (e.g. 2D A-form with A in-plane, so curl(A) = B_z z_hat)
     or axial H_z: nu_zz / mu_zz enter; factor is (rho'/R)^4 or
     (R/rho')^4 respectively. Helper: kelvin_{nu,mu}_factor_2d_axial_cf.

Reference: H. Nagamine, T. Yamaguchi, K. Sugahara, "A Pullback-Based
Formulation of Kelvin Transformation in Electromagnetic Field Analysis,"
CEFC 2026 (Thessaloniki) id 350 (with Sugahara as co-author); see also
Sugahara 2022 IEEE TransMag 58(9) [ref [3] in Nagamine]. Canonical
declaration: docs/kelvin/legacy_assets/kelvin_transformation/CONVENTION.md.
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

    See docs/kelvin/KELVIN_TRANSFORMATION.md §2
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


# ---------- GridFunction evaluation in physical coordinates ---------------


def eval_Omega_physical_from_gf(gf_Omega_comp, mesh, r_phys, center, R):
    """Physical Omega(r_phys) from a GridFunction storing Omega_comp on
    the Kelvin outer-sphere FEM mesh.

    Omega is a SCALAR potential (0-form). Under the Kelvin inversion the
    0-form pullback is trivial function composition — no Jacobian factor
    and no Householder reflection:

        Omega_phys(r_phys) = Omega_comp(r'),   r' = kelvin_map_3d(r_phys, c, R)

    Provided as the 0-form member of the eval_{Omega, A, B}_physical_from_gf
    family so that the differential-form hierarchy is complete:

        0-form (Omega)  : factor 1                          (scalar)
        1-form (A)      : (R/rho')^2 H                      (covector)
        2-form (B)      : -(R/rho')^4 H  (pseudovector)     (2-form)

    Args:
        gf_Omega_comp: scalar NGSolve GridFunction or CoefficientFunction,
            callable as ``gf_Omega_comp(mip)``.
        mesh: NGSolve Mesh.
        r_phys: (3,) physical-frame coordinate.
        center, R: Kelvin sphere center and radius in computational coords.

    Returns:
        float Omega_phys(r_phys).
    """
    r_phys = np.asarray(r_phys, dtype=float).reshape(3)
    r_prime = kelvin_map_3d(r_phys, center, R)
    mip = mesh(float(r_prime[0]), float(r_prime[1]), float(r_prime[2]))
    return float(gf_Omega_comp(mip))


def eval_A_physical_from_gf(gf_A_comp, mesh, r_phys, center, R):
    """Physical A(r_phys) from a GridFunction storing A_comp on the Kelvin
    outer-sphere FEM mesh.

    For evaluation of FEM results at physical points outside the inner
    domain: the Kelvin exterior stores A_comp in the COMPUTATIONAL frame,
    which is not directly the physical A at r_phys. Because the Kelvin
    map is involutive (``phi(phi(r)) = r``) and the 1-form pullback
    factor ``(R/rho)^2 H`` is self-inverse up to (R/rho')^2 -> (rho'/R)^-2,
    the inverse pullback is obtained by calling the same formula with
    the query point as the "r_prime" argument:

        A_phys(r_phys) = (R / |r_phys - c|)**2 * H(n) * A_comp(r')
                       = kelvin_pullback_vector(A_comp, r_phys, c, R)

    where r' = kelvin_map_3d(r_phys, c, R) is the computational coordinate
    and n = (r_phys - c) / |r_phys - c| = (r' - c) / |r' - c|.

    Args:
        gf_A_comp: NGSolve 3-component GridFunction (HCurl or VectorH1)
            containing the computed A in the computational frame. Must
            be callable as ``gf_A_comp[i](mip)`` at a mesh point.
        mesh: NGSolve Mesh (the one gf_A_comp lives on).
        r_phys: (3,) physical-frame coordinate where A is wanted.
        center: (3,) Kelvin outer-sphere center in computational coords.
        R: Kelvin radius.

    Returns:
        (3,) numpy array A_phys(r_phys).
    """
    r_phys = np.asarray(r_phys, dtype=float).reshape(3)
    r_prime = kelvin_map_3d(r_phys, center, R)
    mip = mesh(float(r_prime[0]), float(r_prime[1]), float(r_prime[2]))
    A_comp = np.array([gf_A_comp[0](mip), gf_A_comp[1](mip),
                        gf_A_comp[2](mip)], dtype=float)
    return kelvin_pullback_vector(A_comp, r_phys, center, R)


def eval_B_physical_from_gf(gf_B_comp, mesh, r_phys, center, R):
    """Physical B(r_phys) from a GridFunction storing B_comp on the Kelvin
    outer-sphere FEM mesh.

    Counterpart of eval_A_physical_from_gf for B (2-form pseudovector).
    By the same involutivity argument:

        B_phys(r_phys) = -(R / |r_phys - c|)**4 * H(n) * B_comp(r')
                       = kelvin_pullback_B_pseudovector(B_comp, r_phys, c, R)

    Args:
        gf_B_comp: NGSolve 3-component GridFunction (or CoefficientFunction)
            containing the computed B in the computational frame.
            Callable as ``gf_B_comp[i](mip)`` at a mesh point, or
            ``gf_B_comp(mip)`` returning a 3-tuple.
        mesh: NGSolve Mesh.
        r_phys: (3,) physical-frame coordinate.
        center, R: Kelvin sphere center and radius in computational coords.

    Returns:
        (3,) numpy array B_phys(r_phys).
    """
    r_phys = np.asarray(r_phys, dtype=float).reshape(3)
    r_prime = kelvin_map_3d(r_phys, center, R)
    mip = mesh(float(r_prime[0]), float(r_prime[1]), float(r_prime[2]))
    try:
        B_comp = np.array([gf_B_comp[0](mip), gf_B_comp[1](mip),
                            gf_B_comp[2](mip)], dtype=float)
    except TypeError:
        # CoefficientFunction returning a 3-tuple
        val = gf_B_comp(mip)
        B_comp = np.array([val[0], val[1], val[2]], dtype=float)
    return kelvin_pullback_B_pseudovector(B_comp, r_phys, center, R)


# ---------- Radia source evaluation with Kelvin pullback -----------------


def eval_H_from_radia_at_physical(radia_obj, r_phys):
    """Evaluate H_s from a Radia object at a physical-frame point.

    Simply calls ``rad.Fld(radia_obj, 'h', r_phys)``.  Provided for
    symmetry with the Kelvin-aware variants below.

    Args:
        radia_obj: Radia object handle (coil, magnet, etc.).
        r_phys: (3,) or (N, 3) physical-frame coordinates.

    Returns:
        (3,) or (N, 3) H field in A/m.
    """
    import radia as rad
    return np.asarray(rad.Fld(radia_obj, 'h', np.asarray(r_phys).tolist()),
                      dtype=float)


def eval_B_from_radia_at_physical(radia_obj, r_phys):
    """Evaluate B_s from a Radia object at a physical-frame point.

    Simply calls ``rad.Fld(radia_obj, 'b', r_phys)``.

    Args:
        radia_obj: Radia object handle.
        r_phys: (3,) or (N, 3) physical-frame coordinates.

    Returns:
        (3,) or (N, 3) B field in Tesla.
    """
    import radia as rad
    return np.asarray(rad.Fld(radia_obj, 'b', np.asarray(r_phys).tolist()),
                      dtype=float)


def eval_H_from_radia_in_kelvin(radia_obj, r_kelvin, center, R):
    """Evaluate physical H_s at a Kelvin-mesh coordinate point.

    Given a point ``r_kelvin`` in the Kelvin exterior mesh region,
    compute the PHYSICAL H field of the Radia source at the
    corresponding physical exterior point, with proper 1-form pullback.

    Steps:
        1. r_phys = kelvin_map_3d(r_kelvin, center, R)  # inverse map
        2. H_phys = rad.Fld(radia_obj, 'h', r_phys)
        3. H_comp = kelvin_pullback_vector(H_phys, r_kelvin, center, R)

    The returned H_comp is what the Kelvin mesh "sees" — use this when
    injecting an H-source into the Kelvin exterior FE space.

    Args:
        radia_obj: Radia object handle (coil, magnet, etc.).
        r_kelvin: (3,) or (N, 3) coordinate in the Kelvin mesh region.
        center: (3,) Kelvin sphere center.
        R: Kelvin sphere radius.

    Returns:
        (3,) or (N, 3) H field in Kelvin computational frame.

    Example::

        # Evaluate coil H field "seen" inside the Kelvin mesh at (0.4, 0, 0)
        H_kelvin = eval_H_from_radia_in_kelvin(coil, [0.4, 0, 0],
                                                center=[0.5, 0, 0], R=0.18)
    """
    import radia as rad
    r_kelvin = np.asarray(r_kelvin, dtype=float)
    r_phys = kelvin_map_3d(r_kelvin, center, R)
    H_phys = np.asarray(rad.Fld(radia_obj, 'h',
                                 r_phys.tolist()), dtype=float)
    return kelvin_pullback_vector(H_phys, r_kelvin, center, R)


def eval_B_from_radia_in_kelvin(radia_obj, r_kelvin, center, R):
    """Evaluate physical B_s at a Kelvin-mesh coordinate point.

    Counterpart of ``eval_H_from_radia_in_kelvin`` for B (2-form
    pseudovector). Uses ``kelvin_pullback_B_pseudovector``.

    Steps:
        1. r_phys = kelvin_map_3d(r_kelvin, center, R)
        2. B_phys = rad.Fld(radia_obj, 'b', r_phys)
        3. B_comp = kelvin_pullback_B_pseudovector(B_phys, r_kelvin, center, R)

    Args:
        radia_obj: Radia object handle.
        r_kelvin: (3,) or (N, 3) coordinate in Kelvin mesh region.
        center: (3,) Kelvin sphere center.
        R: Kelvin sphere radius.

    Returns:
        (3,) or (N, 3) B field in Kelvin computational frame (Tesla).
    """
    import radia as rad
    r_kelvin = np.asarray(r_kelvin, dtype=float)
    r_phys = kelvin_map_3d(r_kelvin, center, R)
    B_phys = np.asarray(rad.Fld(radia_obj, 'b',
                                 r_phys.tolist()), dtype=float)
    return kelvin_pullback_B_pseudovector(B_phys, r_kelvin, center, R)


def eval_field_from_radia_with_kelvin(radia_obj, points, center, R,
                                       field_type='b'):
    """Evaluate Radia source field at arbitrary points, with automatic
    Kelvin handling for points inside the exterior domain.

    For each point:
    - If ``|p - center| < R`` (inside Kelvin mesh region): evaluate
      at the physical exterior point (Kelvin inverse map) and apply
      the appropriate pullback.
    - Otherwise (physical region): evaluate directly via rad.Fld.

    This is the student-friendly "just give me the field" function.

    Args:
        radia_obj: Radia object handle.
        points: (3,) or (N, 3) observation points.
        center: (3,) Kelvin sphere center.
        R: Kelvin sphere radius.
        field_type: 'h' for H [A/m] or 'b' for B [T].

    Returns:
        (3,) or (N, 3) field values.

    Example::

        # Evaluate coil B everywhere (auto Kelvin for exterior points)
        pts = np.array([[0, 0, 0], [0.4, 0, 0], [0, 0, 0.5]])
        B = eval_field_from_radia_with_kelvin(coil, pts, [0.5,0,0], 0.18)
    """
    import radia as rad
    points = np.asarray(points, dtype=float)
    single = points.ndim == 1
    if single:
        points = points.reshape(1, 3)

    in_kelvin = is_in_kelvin_exterior_domain(points, center, R)
    result = np.zeros_like(points)

    # Physical-region points: direct evaluation
    phys_mask = ~in_kelvin
    if np.any(phys_mask):
        pts_phys = points[phys_mask]
        result[phys_mask] = np.asarray(
            rad.Fld(radia_obj, field_type, pts_phys.tolist()),
            dtype=float).reshape(-1, 3)

    # Kelvin-region points: inverse map + evaluate + pullback
    if np.any(in_kelvin):
        pts_kelvin = points[in_kelvin]
        pts_phys_mapped = kelvin_map_3d(pts_kelvin, center, R)
        field_phys = np.asarray(
            rad.Fld(radia_obj, field_type, pts_phys_mapped.tolist()),
            dtype=float).reshape(-1, 3)
        if field_type == 'h':
            result[in_kelvin] = kelvin_pullback_vector(
                field_phys, pts_kelvin, center, R)
        elif field_type == 'b':
            result[in_kelvin] = kelvin_pullback_B_pseudovector(
                field_phys, pts_kelvin, center, R)

    return result[0] if single else result


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
            the vector magnitude -- exact only for tangential A) or
            'pullback' (Phase 2, full 1-form Jacobian with Householder
            reflection via kelvin_pullback_vector).
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


def B_s_at_obs_with_kelvin(obs_points, filament_paths, currents,
                            kelvin_center=None, R_kelvin=None):
    """Evaluate B_s = mu_0 * H_s at observation points with Kelvin handling.

    For the reduced-A RHS, we need curl(A_s) = B_s. Instead of projecting
    A_s into HCurl and taking curl (which gives piecewise-constant B and
    large errors in the Kelvin domain), this function evaluates B_s
    DIRECTLY at each point for subsequent H1 projection.

    Inner domain: B = mu_0 * Biot-Savart H.
    Kelvin exterior: map to physical, compute B_phys, apply 2-form pullback.

    Args:
        obs_points: (N, 3) coordinates.
        filament_paths: list of K filaments, each [(p1,p2), ...].
        currents: length-K currents.
        kelvin_center, R_kelvin: Kelvin sphere (None = no Kelvin).

    Returns:
        (N, 3) real array of B values [T].
    """
    from biot_savart import h_segments_batch

    obs = np.asarray(obs_points, dtype=float)
    if obs.ndim == 1:
        obs = obs.reshape(1, 3)
    N = obs.shape[0]

    # Flatten all filament segments with per-segment current
    all_segs = []
    for path, Ik in zip(filament_paths, currents):
        for p1, p2 in path:
            all_segs.append((p1, p2, float(np.real(Ik))))

    def _eval_B(pts):
        """B = mu_0 * sum_seg H_seg at points."""
        B = np.zeros((len(pts), 3), dtype=float)
        for p1, p2, I_seg in all_segs:
            H = h_segments_batch([(p1, p2)], pts, current=I_seg)
            B += MU_0 * H
        return B

    if kelvin_center is None or R_kelvin is None:
        return _eval_B(obs)

    c = np.asarray(kelvin_center, dtype=float)
    in_kelvin = is_in_kelvin_exterior_domain(obs, c, R_kelvin)

    B = np.zeros((N, 3), dtype=float)

    phys_idx = np.where(~in_kelvin)[0]
    if len(phys_idx) > 0:
        B[phys_idx] = _eval_B(obs[phys_idx])

    kelv_idx = np.where(in_kelvin)[0]
    if len(kelv_idx) > 0:
        r_phys = kelvin_map_3d(obs[kelv_idx], c, R_kelvin)
        B_phys = _eval_B(r_phys)
        B[kelv_idx] = kelvin_pullback_B_pseudovector(
            B_phys, obs[kelv_idx], c, R_kelvin)

    return B


# ---------- FEM (NGSolve) material CoefficientFunction factory --------
#
# Helpers for the canonical Kelvin material modulation per
#   H. Nagamine, T. Yamaguchi, K. Sugahara,
#   "A Pullback-Based Formulation of Kelvin Transformation in
#    Electromagnetic Field Analysis," CEFC 2026 id 350
# (with Sugahara as co-author); see also Sugahara 2022 IEEE TransMag
# 58(9) [ref [3] in Nagamine]. Full declaration in
#   docs/kelvin/legacy_assets/kelvin_transformation/CONVENTION.md
# and derivation in
#   docs/kelvin/KELVIN_TRANSFORMATION.md §2 sec 8.
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


def kelvin_nu_factor_3d_cf(center, R, coords=None, rho_min=None):
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
        rho_min: minimum rho' for regularization.  If None, uses 1e-10
            (legacy behaviour). Set to the Kelvin-region mesh element
            size (e.g. maxh_kelvin) to cap the metric and stabilize
            higher-order (order >= 3) FEM.

    Returns:
        NGSolve CoefficientFunction.
    """
    from ngsolve import x, y, z, sqrt, IfPos
    if coords is None:
        coords = (x, y, z)
    if rho_min is None:
        rho_min = 1e-10
    cx, cy, cz = center
    rho2 = (coords[0] - cx) ** 2 + (coords[1] - cy) ** 2 + (coords[2] - cz) ** 2
    rho_prime = sqrt(rho2)
    rho_safe = IfPos(rho_prime - rho_min, rho_prime, rho_min)
    return (rho_safe / R) ** 2


def kelvin_mu_factor_3d_cf(center, R, coords=None, rho_min=None):
    """3D spherical Kelvin mu factor (R/rho')^2 as NGSolve CF.

    Multiply by mu_0 to produce the Kelvin exterior-domain permeability:
        mu_ext = mu_0 * kelvin_mu_factor_3d_cf(...)

    Diverges at rho'=0 (mu_r -> infinity at image of infinity), equals
    mu_0 at rho'=R. Reciprocal of kelvin_nu_factor_3d_cf.

    Args:
        center, R, coords: see kelvin_nu_factor_3d_cf.
        rho_min: minimum rho' for regularization.  Caps the mu factor
            at (R/rho_min)^2.  Set to the Kelvin-region mesh element
            size to stabilize higher-order FEM (order >= 3).
            Default None -> 1e-10 (legacy).
    """
    from ngsolve import x, y, z, sqrt, IfPos
    if coords is None:
        coords = (x, y, z)
    if rho_min is None:
        rho_min = 1e-10
    cx, cy, cz = center
    rho2 = (coords[0] - cx) ** 2 + (coords[1] - cy) ** 2 + (coords[2] - cz) ** 2
    rho_prime = sqrt(rho2)
    rho_safe = IfPos(rho_prime - rho_min, rho_prime, rho_min)
    return (R / rho_safe) ** 2


def kelvin_nu_factor_axisym_cf(z_offset, R, r_coord=None, z_coord=None):
    """Axisym (r,z) Kelvin nu factor (rho'/R)^2 with Z-offset.

    rho' = sqrt(r^2 + (z - z_offset)^2)  (3D spherical distance in the
    meridional plane; this is 3D sphere Kelvin viewed axisymmetrically).

    NGSolve axisym convention: r = x, z = y. Pass explicit r_coord /
    z_coord when your mesh uses a different embedding.

    Multiply by nu_0 for A-formulation nu_ext.

    IMPORTANT -- NGSolve has no built-in axisymmetric mode. The r-weight
    from the axisymmetric volume element (2 pi r dr dz) must be written
    explicitly by the caller:

        a += nu_cf * grad(u) * grad(v) * r_coord * dx

    The Kelvin factor returned here ONLY handles the sphere-inversion
    pullback; it does NOT absorb the r-weight. A common bug is to drop
    the ``* r_coord`` factor and get results that are off by O(r) near
    the axis.
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
    """Axisym (r,z) Kelvin mu factor (R/rho')^2. Reciprocal of nu factor.

    Same r-weight caveat as kelvin_nu_factor_axisym_cf: the caller must
    write ``* r_coord * dx`` explicitly (NGSolve has no built-in axisym
    mode). This helper only returns the sphere-inversion pullback.
    """
    from ngsolve import x, y, sqrt, IfPos
    if r_coord is None:
        r_coord = x
    if z_coord is None:
        z_coord = y
    z_local = z_coord - z_offset
    rho_prime = sqrt(r_coord ** 2 + z_local ** 2)
    rho_safe = IfPos(rho_prime - 1e-10, rho_prime, 1e-10)
    return (R / rho_safe) ** 2


def kelvin_factor_2d_inplane_cf():
    """2D Kelvin factor for IN-PLANE (Bx, By) / (Hx, Hy) problems: 1.

    Under 2D cylindrical Kelvin inversion k(rho,phi,z) = (R^2/rho, phi, z)
    the bilinear form ``integral(coeff * grad(u) . grad(v)) dA`` is
    INVARIANT when grad(u) is in-plane: the gradient factor (rho'/R)^2
    and the volume factor (R/rho')^4 cancel with the in-plane metric,
    leaving the coefficient unmodulated.

    Nagamine CEFC 2026 eq. 12 gives the full anisotropic tensor
    ``nu' = diag(1, 1, (rho'/R)^4) nu`` (and reciprocal for mu'). The
    z-slot -- (rho'/R)^4 for nu, (R/rho')^4 for mu -- only enters the
    bilinear form when the problem involves the Z-component of B or H.
    For the **common** 2D magnetostatic setups below, the z-slot is
    never contracted and this helper returns literally **1**:

    - Omega-form with scalar potential phi(x, y): grad(phi) = H in-plane
    - A-form with A = A_z(x, y) z_hat: curl(A) = B in-plane

    This helper is the uniform-API counterpart of the 3D/axisym factor
    helpers, so ``build_material_cf(mesh, mu0, kelvin_factor_2d_inplane_cf())``
    composes without special-casing. Kelvin factor IS 1 here, not absent.

    The deprecated names kelvin_{mu,nu}_factor_2d_Hx_Hy_cf and
    kelvin_nu_factor_2d_Az_cf are aliases of this helper (the A_z name
    was misleading: A_z scalar formulations ALSO land in the in-plane
    case above, so the factor is 1, not (rho'/R)^4).
    """
    from ngsolve import CoefficientFunction
    return CoefficientFunction(1.0)


def kelvin_nu_factor_2d_axial_cf(offset, R, x_coord=None, y_coord=None):
    """2D Kelvin nu factor for AXIAL (B_z) problems: (rho'/R)^4.

    Nagamine CEFC 2026 eq. 12, z-slot of ``nu' = diag(1, 1, (rho'/R)^4) nu``.
    Applies when B has a z-component -- e.g. 2D A-formulation with an
    in-plane vector potential A = (A_x(x,y), A_y(x,y), 0), for which
    curl(A) = B_z z_hat and nu_zz contracts into the bilinear form.

    For the common A = A_z z_hat case (B in-plane, B_z=0), use
    kelvin_factor_2d_inplane_cf() instead -- the axial slot does NOT
    enter the bilinear form there.

    rho' = sqrt((x - x_off)^2 + (y - y_off)^2)

    Args:
        offset: (x_off, y_off) -- Kelvin circle center in world coords.
        R: Kelvin radius.

    Returns:
        NGSolve CoefficientFunction (rho'/R)^4.
    """
    from ngsolve import x, y, sqrt, IfPos
    if x_coord is None:
        x_coord = x
    if y_coord is None:
        y_coord = y
    xo, yo = offset
    rho_prime = sqrt((x_coord - xo) ** 2 + (y_coord - yo) ** 2)
    rho_safe = IfPos(rho_prime - 1e-10, rho_prime, 1e-10)
    return (rho_safe / R) ** 4


def kelvin_mu_factor_2d_axial_cf(offset, R, x_coord=None, y_coord=None):
    """2D Kelvin mu factor for AXIAL (H_z) problems: (R/rho')^4.

    Reciprocal of kelvin_nu_factor_2d_axial_cf. Applies when H has a
    z-component -- unusual in standard 2D magnetostatics (where the
    common in-plane-H case uses factor = 1).
    """
    from ngsolve import x, y, sqrt, IfPos
    if x_coord is None:
        x_coord = x
    if y_coord is None:
        y_coord = y
    xo, yo = offset
    rho_prime = sqrt((x_coord - xo) ** 2 + (y_coord - yo) ** 2)
    rho_safe = IfPos(rho_prime - 1e-10, rho_prime, 1e-10)
    return (R / rho_safe) ** 4


# Deprecated aliases -- keep callers working while the ecosystem migrates.
# Previous misleading signatures (especially kelvin_nu_factor_2d_Az_cf taking
# (offset, R) and returning (rho'/R)^4) are corrected to the in-plane factor
# of 1, matching the actual physics of the A_z scalar formulation.
def kelvin_mu_factor_2d_Hx_Hy_cf():
    """DEPRECATED alias. Use kelvin_factor_2d_inplane_cf()."""
    return kelvin_factor_2d_inplane_cf()


def kelvin_nu_factor_2d_Hx_Hy_cf():
    """DEPRECATED alias. Use kelvin_factor_2d_inplane_cf()."""
    return kelvin_factor_2d_inplane_cf()


def kelvin_nu_factor_2d_Az_cf(offset=None, R=None, x_coord=None, y_coord=None):
    """DEPRECATED. A_z scalar formulation has in-plane B, so factor = 1.

    The original implementation returned (rho'/R)^4, but that is the
    nu_zz slot (for axial B_z, e.g. in-plane A) and does NOT enter the
    bilinear form when A = A_z z_hat. Use kelvin_factor_2d_inplane_cf()
    (returns 1) for the standard A_z scalar case, or
    kelvin_nu_factor_2d_axial_cf(offset, R) if you actually have B_z.
    """
    return kelvin_factor_2d_inplane_cf()


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


# ---------- HCurl projection of the Biot-Savart source ---------------------


def curl_of_vector_cf(A_cf):
    """Build curl(A) as a symbolic CF via component Diff.

    NGSolve does not overload ``curl()`` for an arbitrary
    ``VectorialCoefficientFunction``, so we assemble it manually:

        curl(A) = (dA_z/dy - dA_y/dz, dA_x/dz - dA_z/dx, dA_y/dx - dA_x/dy)

    This is the natural formula when A is a smooth vector field built
    from symbolic expressions in (x, y, z). NGSolve's symbolic
    differentiation handles sqrt / log terms via the chain rule.
    """
    from ngsolve import x, y, z, CF
    cx = A_cf[2].Diff(y) - A_cf[1].Diff(z)
    cy = A_cf[0].Diff(z) - A_cf[2].Diff(x)
    cz = A_cf[1].Diff(x) - A_cf[0].Diff(y)
    return CF((cx, cy, cz))


def biot_savart_A_cf(filament_paths, fil_currents):
    """Symbolic NGSolve CoefficientFunction for Biot-Savart A from a filament bundle.

    Uses the closed-form A for a straight current segment (avoids any
    L2 projection / quadrature-truncation error):

        A_seg(r) = (mu_0 I / 4 pi) * dl_hat * ln((|r-p2| + sB)/(|r-p1| + sA))
          where sA = -dl_hat . (r - p1), sB = -dl_hat . (r - p2)
                dl_hat = (p2 - p1) / |p2 - p1|

    Sums over all segments of all filaments. No Kelvin pullback here —
    valid only for evaluating A_s at points in the INNER (physical)
    domain. Do NOT use this CF to evaluate in a Kelvin-transformed
    region; the PEEC+FEM scattered workflow only needs A_s on the
    workpiece sibc boundary (inner) and on the filament line (inner),
    so the restriction is not a practical limit there.

    Args:
        filament_paths: list of K filaments, each a list of
            ``(p1, p2)`` segment endpoint tuples.
        fil_currents: length-K currents (real or complex).

    Returns:
        NGSolve 3-component CoefficientFunction (complex if any
        current is complex).
    """
    from ngsolve import x, y, z, sqrt, log, CF

    # Detect complex currents
    is_complex = any(np.iscomplexobj(I) or (isinstance(I, complex) and I.imag != 0)
                     for I in fil_currents)

    A_x = CF(0.0j if is_complex else 0.0)
    A_y = CF(0.0j if is_complex else 0.0)
    A_z = CF(0.0j if is_complex else 0.0)

    four_pi = 4 * math.pi
    for path, I_k in zip(filament_paths, fil_currents):
        I_val = complex(I_k) if is_complex else float(I_k)
        pref = MU_0 * I_val / four_pi
        for (p1, p2) in path:
            p1 = np.asarray(p1, dtype=float)
            p2 = np.asarray(p2, dtype=float)
            dl = p2 - p1
            L = np.linalg.norm(dl)
            if L < 1e-15:
                continue
            dlh = dl / L
            R1x = x - p1[0]; R1y = y - p1[1]; R1z = z - p1[2]
            R2x = x - p2[0]; R2y = y - p2[1]; R2z = z - p2[2]
            # Small epsilon guards for points ON or very close to the filament
            eps = 1e-18
            L1 = sqrt(R1x * R1x + R1y * R1y + R1z * R1z + eps)
            L2 = sqrt(R2x * R2x + R2y * R2y + R2z * R2z + eps)
            sA = -(dlh[0] * R1x + dlh[1] * R1y + dlh[2] * R1z)
            sB = -(dlh[0] * R2x + dlh[1] * R2y + dlh[2] * R2z)
            # log argument; both numerator and denominator are positive
            # for observation off the segment
            logterm = log((L2 + sB + eps) / (L1 + sA + eps))
            contrib = pref * logterm
            A_x = A_x + contrib * dlh[0]
            A_y = A_y + contrib * dlh[1]
            A_z = A_z + contrib * dlh[2]

    return CF((A_x, A_y, A_z))


def biot_savart_B_cf(filament_paths, fil_currents):
    """Symbolic NGSolve CoefficientFunction for Biot-Savart B from a filament bundle.

    Closed-form B from a finite straight current segment (the curl of the
    closed-form A in biot_savart_A_cf):

        B_seg(r) = (mu_0 I / 4 pi) * (dl_hat x R1) / d_perp^2
                                   * (l1/L1 - l2/L2)
          where  R1 = r - p1
                 dl_hat = (p2 - p1) / |p2 - p1|
                 l1 = dl_hat . R1, l2 = dl_hat . R2
                 L1 = |R1|, L2 = |R2|
                 d_perp^2 = L1^2 - l1^2

    Equivalent to ``curl(biot_savart_A_cf(...))`` but evaluated as a
    direct CF (no NGSolve AD on the symbolic A, no L2 projection).
    Used by the scattered FEM-SIBC formulation to supply the surface
    term ``-nu * (n x B_s) . v ds`` on the workpiece SIBC boundary
    (the term that integration-by-parts of curl-curl in the scattered
    decomposition demands).

    Args:
        filament_paths: list of K filaments, each a list of
            ``(p1, p2)`` segment endpoint tuples.
        fil_currents: length-K currents (real or complex).

    Returns:
        NGSolve 3-component CoefficientFunction (complex if any
        current is complex).
    """
    from ngsolve import x, y, z, sqrt, CF

    is_complex = any(np.iscomplexobj(I) or (isinstance(I, complex) and I.imag != 0)
                     for I in fil_currents)

    Bx = CF(0.0j if is_complex else 0.0)
    By = CF(0.0j if is_complex else 0.0)
    Bz = CF(0.0j if is_complex else 0.0)

    four_pi = 4 * math.pi
    eps = 1e-18
    for path, I_k in zip(filament_paths, fil_currents):
        I_val = complex(I_k) if is_complex else float(I_k)
        pref = MU_0 * I_val / four_pi
        for (p1, p2) in path:
            p1 = np.asarray(p1, dtype=float)
            p2 = np.asarray(p2, dtype=float)
            dl = p2 - p1
            L = np.linalg.norm(dl)
            if L < 1e-15:
                continue
            dlh = dl / L
            R1x = x - p1[0]; R1y = y - p1[1]; R1z = z - p1[2]
            R2x = x - p2[0]; R2y = y - p2[1]; R2z = z - p2[2]
            L1 = sqrt(R1x * R1x + R1y * R1y + R1z * R1z + eps)
            L2 = sqrt(R2x * R2x + R2y * R2y + R2z * R2z + eps)
            # sA = -l1, sB = -l2 (matching biot_savart_A_cf convention)
            sA = -(dlh[0] * R1x + dlh[1] * R1y + dlh[2] * R1z)
            sB = -(dlh[0] * R2x + dlh[1] * R2y + dlh[2] * R2z)
            # angle factor (l1/L1 - l2/L2) = (-sA/L1) - (-sB/L2) = sB/L2 - sA/L1
            angle = sB / L2 - sA / L1
            # d_perp^2 = R1.R1 - l1^2 = R1.R1 - sA^2  (sA^2 = l1^2)
            d_perp_sq = (R1x * R1x + R1y * R1y + R1z * R1z) - sA * sA + eps
            # dl_hat x R1
            cx = dlh[1] * R1z - dlh[2] * R1y
            cy = dlh[2] * R1x - dlh[0] * R1z
            cz = dlh[0] * R1y - dlh[1] * R1x
            factor = pref * angle / d_perp_sq
            Bx = Bx + factor * cx
            By = By + factor * cy
            Bz = Bz + factor * cz

    return CF((Bx, By, Bz))


def project_A_s_to_hcurl(mesh, filament_paths, fil_currents,
                         kelvin_center=None, R_kelvin=None,
                         factor_mode='pullback', is_complex=None,
                         order=1, n_quad=8):
    """Project the Biot-Savart A_s (filament bundle) onto an HCurl GridFunction.

    Handles both real and complex per-filament currents, with optional
    Kelvin-exterior pullback for meshes that include a Kelvin-transformed
    outer domain.

    Workaround: ``HCurl(complex=True).Set(H1(complex=True))`` currently
    projects to zero in NGSolve 6.2.2603.  This helper projects the real
    and imaginary parts separately through a real HCurl space, then
    combines the dof vectors into the complex target GridFunction.

    Args:
        mesh: NGSolve Mesh.
        filament_paths: list of K filaments, each a list of ``(p1, p2)``
            segment endpoint tuples.
        fil_currents: length-K currents.  If any is complex, the output
            is complex; ``is_complex=`` can override.
        kelvin_center, R_kelvin: exterior-Kelvin sphere parameters.
            Pass both to enable pullback; ``None`` disables.
        factor_mode: ``'pullback'`` (full 1-form pullback, default) or
            ``'scalar'`` (Phase 1 approximation).
        is_complex: override for complex HCurl target.  Default: auto-
            detect from ``fil_currents``.
        order: HCurl (and intermediate H1) polynomial order.
        n_quad: quadrature per segment inside ``A_s_at_obs_with_kelvin``.

    Returns:
        ``GridFunction`` on ``HCurl(mesh, order=order, complex=...)``
        holding the projected A_s.
    """
    from ngsolve import H1, HCurl, GridFunction

    Ik_arr = np.asarray(list(fil_currents))
    if is_complex is None:
        is_complex = bool(np.iscomplexobj(Ik_arr)
                          and np.any(Ik_arr.imag != 0))

    nv = mesh.nv
    obs = np.array([mesh.vertices[v].point for v in range(nv)])

    if kelvin_center is not None and R_kelvin is not None:
        A_vals = A_s_at_obs_with_kelvin(
            obs, filament_paths, fil_currents,
            kelvin_center=kelvin_center, R_kelvin=R_kelvin,
            factor_mode=factor_mode, n_quad=n_quad)
    else:
        A_vals = biot_savart_A_at_points(
            obs, filament_paths, fil_currents, n_quad=n_quad)

    def _project_real(A_real):
        fes_h1 = H1(mesh, order=order, dim=3, complex=False)
        gf_h1 = GridFunction(fes_h1)
        fv = gf_h1.vec.FV()
        for vi in range(nv):
            for k in range(3):
                fv[vi * 3 + k] = float(A_real[vi, k])
        fes_c = HCurl(mesh, order=order, complex=False)
        gf_c = GridFunction(fes_c)
        gf_c.Set(gf_h1)
        return gf_c

    if not is_complex:
        A_real = A_vals.real if np.iscomplexobj(A_vals) else np.asarray(A_vals, float)
        return _project_real(A_real)

    # Complex path: project Re/Im through the real HCurl space separately
    # (HCurl.Set on a complex H1 currently returns the zero vector).
    gf_As_re = _project_real(A_vals.real)
    gf_As_im = _project_real(A_vals.imag)

    fes_c = HCurl(mesh, order=order, complex=True)
    gf_c = GridFunction(fes_c)
    v_re = np.asarray(gf_As_re.vec.FV())
    v_im = np.asarray(gf_As_im.vec.FV())
    fv_c = gf_c.vec.FV()
    for i in range(len(v_re)):
        fv_c[i] = complex(float(v_re[i]), float(v_im[i]))
    return gf_c


# ---------- Back-reaction: line integral of A along filaments ---------------


def line_integral_A_filaments(gfu, mesh, filament_paths, n_quad=4):
    """Line integral of vector potential A along each filament path.

    Computes  Delta_Phi_k = integral_{path_k} A . dl  for each filament,
    using Gauss-Legendre quadrature with *n_quad* points per segment.

    This is the core primitive for PEEC+FEM back-reaction: the scattered
    field A_r (from FEM-SIBC solve) links through each PEEC filament,
    producing a flux-linkage change that modifies the effective inductance.

    Args:
        gfu: NGSolve GridFunction (HCurl, possibly complex).  This is
            typically A_r (the scattered/reduced field from the FEM solve).
        mesh: NGSolve Mesh on which *gfu* is defined.
        filament_paths: list of K filaments, each a list of ``(p1, p2)``
            segment endpoint tuples (coordinates in meters).
        n_quad: Gauss-Legendre quadrature points per segment (default 4).

    Returns:
        (K,) complex ndarray of flux linkages Delta_Phi_k [Wb].
    """
    t_gl, w_gl = np.polynomial.legendre.leggauss(n_quad)
    t01 = 0.5 * (t_gl + 1.0)   # map [-1,1] -> [0,1]
    w01 = 0.5 * w_gl

    n_fil = len(filament_paths)
    delta_phi = np.zeros(n_fil, dtype=complex)

    for k, path in enumerate(filament_paths):
        for p1, p2 in path:
            p1a = np.asarray(p1, float)
            p2a = np.asarray(p2, float)
            dl = p2a - p1a
            if np.linalg.norm(dl) < 1e-15:
                continue
            for q in range(n_quad):
                pt = p1a + t01[q] * dl
                try:
                    A_val = gfu(mesh(*pt))
                except Exception:
                    continue  # skip quadrature points outside mesh
                dot_A_dl = sum(complex(A_val[i]) * dl[i] for i in range(3))
                delta_phi[k] += w01[q] * dot_A_dl

    return delta_phi


def compute_back_reaction(fil_currents, delta_phi, I_total):
    """Compute Delta_L and Delta_R from PEEC filament flux linkages.

    Physics:  The scattered field A_r (from FEM-SIBC on the workpiece)
    threads through each PEEC filament.  The resulting impedance change
    at the coil port is

        Delta_Z = j * omega * sum_k conj(I_k) * Delta_Phi_k / |I_total|^2

    Decomposed into inductance and resistance changes:

        Delta_L =  Re(sum_k conj(I_k) * Delta_Phi_k) / |I_total|^2
        Delta_R = -omega * Im(sum_k conj(I_k) * Delta_Phi_k) / |I_total|^2

    Sign conventions (validated by physical expectation):
        Cu  workpiece  (Lenz law):           Delta_L < 0
        Steel workpiece (flux concentration): Delta_L > 0

    Args:
        fil_currents: length-K complex array of per-filament currents [A].
        delta_phi: length-K complex array from line_integral_A_filaments [Wb].
        I_total: scalar total port current [A].

    Returns:
        dict with keys:
            Delta_L           -- inductance change [H] (real)
            Delta_R_over_omega -- resistance change / omega [Ohm*s/rad]
            flux_linkage_sum  -- raw complex sum (diagnostic)
    """
    I_k = np.asarray(fil_currents, dtype=complex)
    dphi = np.asarray(delta_phi, dtype=complex)

    flux_sum = np.sum(np.conj(I_k) * dphi)
    I2 = abs(float(I_total)) ** 2

    return {
        'Delta_L': float(flux_sum.real / I2),
        'Delta_R_over_omega': float(-flux_sum.imag / I2),
        'flux_linkage_sum': complex(flux_sum),
    }
