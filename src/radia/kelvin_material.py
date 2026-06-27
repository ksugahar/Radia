"""kelvin_material.py

Layer 2 of the Kelvin helper API: build NGSolve CoefficientFunctions
for the Kelvin-modulated material parameter (nu) and for external
source / background fields applied in the Kelvin exterior domain.

Canonical convention (Nagamine, Yamaguchi, Sugahara,
"A Pullback-Based Formulation of Kelvin Transformation in
Electromagnetic Field Analysis," CEFC 2026 id 350; see also
Sugahara 2022 IEEE TransMag 58(9) [ref [3] in Nagamine]):

    3D spherical (conformal): nu' = (rho'/R)^2 * nu_0
                              mu' = (R/rho')^2 * mu_0   (reciprocal)

Derivation via pullback of orthonormal 1-form basis + bilinear energy
functional (Nagamine eq. 9). Numerical validation on toroidal current
loop: analytical dipole exterior energy 3.333e-8 J vs FEM on Omega'
= 3.344e-8 J (+0.33%).

TWO DIFFERENT BACKGROUND-FIELD CONVENTIONS (DON'T CONFUSE):

(A) Solution / 1-form pullback for source localized in physical region
    (e.g. PEEC filament coil in inner air):
        A_comp(r') = (R/rho')^2 * H * A_phys(r_phys)        (1-form)
        B_comp(r') = -(R/rho')^4 * H * B_phys(r_phys)       (2-form)
    Helper: make_kelvin_aware_A_s_cf (full pullback, evaluates A_phys
    at the Kelvin-mapped physical point r_phys = T(r')).
    Has 1/rho'^3 singularity at offset for unbounded A_phys.

(B) Reduced-potential background field for source defined globally
    (e.g. uniform B_0 z_hat applied at infinity, dipole / quadrupole
    background fields). Per Sugahara-Nagamine-Kameari internal note
    (2026, reflected in docs/kelvin/KELVIN_TRANSFORMATION.md §7), the
    metric-tensor scaling gives:
        H_s'(r') = -(rho'/R)^2 * H_s(r')   (3D, evaluated at comp coords)
        A_s'(r') = -(rho'/R)^2 * A_s(r')   (3D, A is also a 1-form)
    Helper: make_reduced_potential_background_cf (vanishes at offset =
    no singularity, sign-flipped for opposite-normal periodic BC matching
    at the Kelvin sphere boundary).
    See docs/Reduced_Potential_Kelvin.md for the derivation.

These two conventions are INCOMPATIBLE: (A) is a covariant 1-form
pullback that preserves line integrals exactly but is singular at
infinity; (B) is a metric-tensor-derived "engineering" formula that is
finite, bounded, and chosen to make the discrete weak form give the
right physical answer for reduced-potential formulations.

See:
    docs/kelvin/CONVENTION.md
        (material modulation)
    docs/kelvin/KELVIN_TRANSFORMATION.md §2  (1-form)
    docs/kelvin/KELVIN_TRANSFORMATION.md §7  (red-pot)
"""

from __future__ import annotations

import math

import numpy as np

MU_0 = 4e-7 * math.pi
NU_0 = 1.0 / MU_0


def make_kelvin_nu_cf(mesh, R_K, offset, nu_0=NU_0,
                       kelvin_mats=("kelvin",)):
    """NGSolve CoefficientFunction for nu(r) with Kelvin modulation.

    nu(r) = nu_0 in non-Kelvin materials,
    nu(r) = nu_0 * (|r' - offset| / R_K)^2 in Kelvin materials.

    Canonical Nagamine CEFC 2026 / Sugahara 2022 convention for 3D
    spherical (conformal) Kelvin: nu' = (rho'/R)^2 * nu_0. Vanishes at
    rho' = 0 (image of infinity), equals nu_0 at rho' = R (continuous).

    Args:
        mesh: NGSolve Mesh.
        R_K: Kelvin sphere radius.
        offset: 3-tuple, Kelvin sphere center.
        nu_0: vacuum reluctivity (default 1/mu_0).
        kelvin_mats: substring(s) used to detect Kelvin materials.
            A material whose name (lowercased) contains any of these
            substrings receives the Kelvin modulation.

    Returns:
        Scalar CoefficientFunction nu(x, y, z).
    """
    from ngsolve import x, y, z

    ox, oy, oz = offset
    rho_prime_sq = (x - ox) ** 2 + (y - oy) ** 2 + (z - oz) ** 2 + 1e-24
    kelvin_fac = rho_prime_sq / (R_K * R_K)    # (rho'/R)^2

    nu_dict = {}
    for m in mesh.GetMaterials():
        ml = m.lower()
        is_kelvin = any(kw in ml for kw in kelvin_mats)
        nu_dict[m] = (nu_0 * kelvin_fac if is_kelvin else nu_0)
    return mesh.MaterialCF(nu_dict, default=nu_0)


def make_kelvin_aware_A_s_cf(mesh, A_phys_factory, R_K, offset,
                              kelvin_mats=("kelvin",)):
    """Build an A_s vector CF using the FULL 1-form pullback (Convention A).

    *** USE THIS FOR PEEC FILAMENT COIL SOURCES ***
    *** localized in the inner physical domain.   ***

    For a source field defined globally (e.g. uniform B_0 z_hat applied
    at infinity, dipole / quadrupole background fields), use
    ``make_reduced_potential_background_cf`` instead — the full pullback
    is singular at offset for unbounded physical fields.

    In non-Kelvin materials: A_s_comp(x, y, z) = A_phys_factory(x, y, z).
    In Kelvin materials: evaluate A_phys at the Kelvin-mapped physical
    point, then apply the Phase 2 1-form pullback (Householder + scalar
    factor (R/rho')^2).

    The Kelvin-mapped position from a computational point r' (within
    the Kelvin sphere centered at offset, radius R_K) is
        r_phys = offset + (R_K^2 / |r' - offset|^2) * (r' - offset)

    Per-component MaterialCF switching is built in so the result is a
    well-formed VectorCF on the entire mesh.

    Args:
        mesh: NGSolve Mesh.
        A_phys_factory: callable ``(x_cf, y_cf, z_cf) -> VectorCF``
            returning a 3-vector CF for A_phys evaluated at the
            given (x, y, z) CFs. Lets the caller plug in a symbolic
            Biot-Savart sum, a wrapper around Radia.Fld, or any other
            analytical formula.
        R_K, offset: Kelvin sphere parameters.
        kelvin_mats: substring(s) used to detect Kelvin materials.

    Returns:
        VectorCoefficientFunction (length 3) suitable for use in
        LinearForm or BilinearForm.
    """
    from ngsolve import x, y, z, sqrt, CoefficientFunction as CF

    A_inner = A_phys_factory(x, y, z)

    ox, oy, oz = offset
    dxp = x - ox
    dyp = y - oy
    dzp = z - oz
    rho_p_sq = dxp * dxp + dyp * dyp + dzp * dzp + 1e-24
    inv = R_K * R_K / rho_p_sq                 # (R/rho')^2 scalar
    kel_x = ox + inv * dxp
    kel_y = oy + inv * dyp
    kel_z = oz + inv * dzp

    A_at_phys = A_phys_factory(kel_x, kel_y, kel_z)

    # 1-form pullback (Phase 2): A_comp = (R/rho')^2 * Householder(A_at_phys)
    # n = (r' - offset) / rho'.
    rho_p = sqrt(rho_p_sq)
    nx = dxp / rho_p
    ny = dyp / rho_p
    nz = dzp / rho_p
    A_dot_n = (A_at_phys[0] * nx + A_at_phys[1] * ny + A_at_phys[2] * nz)
    refl_x = A_at_phys[0] - 2.0 * A_dot_n * nx
    refl_y = A_at_phys[1] - 2.0 * A_dot_n * ny
    refl_z = A_at_phys[2] - 2.0 * A_dot_n * nz
    A_kelvin_x = inv * refl_x
    A_kelvin_y = inv * refl_y
    A_kelvin_z = inv * refl_z

    def _switch(kelvin_comp, inner_comp):
        d = {}
        for m in mesh.GetMaterials():
            ml = m.lower()
            is_kelvin = any(kw in ml for kw in kelvin_mats)
            d[m] = kelvin_comp if is_kelvin else inner_comp
        return mesh.MaterialCF(d, default=inner_comp)

    Ax = _switch(A_kelvin_x, A_inner[0])
    Ay = _switch(A_kelvin_y, A_inner[1])
    Az = _switch(A_kelvin_z, A_inner[2])
    return CF((Ax, Ay, Az))


def make_reduced_potential_background_cf(mesh, F_inner_factory, R_K, offset,
                                          kelvin_mats=("kelvin",),
                                          dim=3):
    """Build a reduced-potential background field CF (Convention B).

    *** USE THIS FOR GLOBALLY-DEFINED BACKGROUND FIELDS ***
    *** (uniform B, dipole, quadrupole at infinity)      ***

    Per Sugahara-Nagamine-Kameari internal note (2026, see
    docs/kelvin/KELVIN_TRANSFORMATION.md §7): a 1-form background
    field (H_s for H-formulation, A_s for
    A-formulation) transforms in the Kelvin exterior with a
    metric-tensor scaling factor:

        F_s'(r') = -(rho'/R)^2 * F_s(r' - offset)   in 3D
        F_s'(r') = -F_s(r' - offset)                in 2D

    where F_s on the RHS is evaluated at coordinates measured from the
    KELVIN SPHERE CENTER (the offset point), NOT at the global origin
    or at the Kelvin-mapped physical point r_phys = T(r'). The local
    evaluation is essential for position-dependent A_s = (B_0/2)(-y,x,0):
    using global coords would introduce a spurious offset-dependent term.

    This differs from the proper 1-form pullback (make_kelvin_aware_A_s_cf)
    in three crucial ways:

    1. NO Householder reflection — preserves the source's direction.
    2. Evaluated at LOCAL coords (offset-relative), NOT at r_phys.
       The physical-frame functional form is reused with offset as the
       new origin, multiplied by the scalar factor.
    3. Vanishes at offset (rho' -> 0) — no singularity, even when the
       physical field is unbounded at infinity (e.g. uniform B_z).

    The sign flip ensures matching at the periodic Kelvin boundary
    (rho' = R) where the inner and exterior normals are opposite.

    Args:
        mesh: NGSolve Mesh.
        F_inner_factory: callable ``(x_cf, y_cf, z_cf) -> VectorCF``
            returning the background field (3-vector) evaluated at the
            given coordinates. The factory receives GLOBAL coords for
            inner-region values and LOCAL (offset-relative) coords for
            Kelvin-region values. For uniform B_0 z_hat, A-formulation:
            ``lambda x, y, z: CF((-y, x, 0)) * (B_0 / 2)``.
        R_K: Kelvin sphere radius.
        offset: 3-tuple, Kelvin sphere center.
        kelvin_mats: substring(s) used to detect Kelvin materials.
        dim: 2 or 3 (default 3). Selects between 2D and 3D scaling.

    Returns:
        VectorCoefficientFunction (length 3) with:
            inner materials:  F_s(x, y, z)                       (global)
            kelvin materials: -(rho'/R)^2 * F_s(x-ox, y-oy, z-oz) (3D, local)
                              -F_s(x-ox, y-oy, z-oz)              (2D, local)

    Example::

        # Uniform B_0 z_hat applied at infinity, A-formulation:
        A_s_cf = make_reduced_potential_background_cf(
            mesh,
            lambda xc, yc, zc: CF((-yc, xc, 0)) * 0.5,   # B_0 = 1
            R_K=R_K, offset=offset, kelvin_mats=("kelvin",))

        # Uniform H_0 z_hat applied at infinity, H-formulation:
        # (constant field; local vs global eval gives the same result)
        H_s_cf = make_reduced_potential_background_cf(
            mesh,
            lambda xc, yc, zc: CF((0, 0, 1)),
            R_K=R_K, offset=offset, kelvin_mats=("kelvin",))

    See:
        docs/kelvin/KELVIN_TRANSFORMATION.md §7
    """
    from ngsolve import x, y, z, CoefficientFunction as CF

    if dim not in (2, 3):
        raise ValueError(f"dim must be 2 or 3, got {dim}")

    ox, oy, oz = offset
    # Inner region: evaluate at global coordinates
    F_inner = F_inner_factory(x, y, z)
    # Kelvin region: evaluate at local (offset-relative) coordinates
    F_local = F_inner_factory(x - ox, y - oy, z - oz)

    if dim == 3:
        rho2 = ((x - ox) ** 2 + (y - oy) ** 2 + (z - oz) ** 2 + 1e-24)
        kelvin_factor = -rho2 / (R_K * R_K)            # -(rho'/R)^2
    else:  # dim == 2
        kelvin_factor = -1.0                            # 2D: sign flip only

    F_kelvin_x = kelvin_factor * F_local[0]
    F_kelvin_y = kelvin_factor * F_local[1]
    F_kelvin_z = kelvin_factor * F_local[2]

    def _switch(kelvin_comp, inner_comp):
        d = {}
        for m in mesh.GetMaterials():
            ml = m.lower()
            is_kelvin = any(kw in ml for kw in kelvin_mats)
            d[m] = kelvin_comp if is_kelvin else inner_comp
        return mesh.MaterialCF(d, default=inner_comp)

    Fx = _switch(F_kelvin_x, F_inner[0])
    Fy = _switch(F_kelvin_y, F_inner[1])
    Fz = _switch(F_kelvin_z, F_inner[2])
    return CF((Fx, Fy, Fz))
