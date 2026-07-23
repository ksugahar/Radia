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
    The right-hand sides are evaluated at LOCAL (offset-relative)
    coordinates r' - offset, NOT at the Kelvin-mapped physical point.

    The 0-form (scalar potential) counterpart, needed by T-Omega and by
    any Omega-reduced route that carries a background POTENTIAL rather
    than a background FIELD in the Kelvin exterior, is
        Omega_s'(r') = -(rho'/R)^2 * Omega_s(r' - offset)   (3D)
    Helper: make_reduced_potential_scalar_cf.

    WARNING (verified symbolically, locked by
    tests/test_reduced_potential_background.py): the 0-form and 1-form
    Convention B rules are NOT gradient-consistent.  For a uniform
    background H_s = H_0 z_hat, Omega_s = -H_0 z,

        curl(H_s' 1-form B) = (2 H_0 / R^2) (-y', x', 0)  != 0

    so the 1-form Convention B field admits NO scalar potential at all,
    and
        -grad(Omega_s' 0-form B) - H_s' (1-form B)
            = -(2 H_0 / R^2) (x' z', y' z', z'^2)  != 0 .

    Pick ONE per formulation and state which: a field-driven weak form
    (H-formulation linear form int(mu H_s . grad v)) uses the 1-form
    helper; a potential-driven weak form (T-Omega, which needs Omega_s
    itself) uses the 0-form helper.  Do not mix them in one model and do
    not assume one can be differentiated into the other.

    See docs/kelvin/KELVIN_TRANSFORMATION.md 7.4 / 7.6 for the derivation
    and the Convention A vs B disambiguation table.

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


def _kelvin_mapped_coords(R_K, offset, phys_center):
    """(kx, ky, kz) = phys_center + (R/rho')^2 (r' - offset), and rho'^2.

    ``phys_center`` is the centre of the PHYSICAL inversion sphere.  It is an
    explicit argument because the repository contains both conventions:
    ``make_kelvin_aware_A_s_cf`` hard-codes ``phys_center = offset`` (the
    mapped point stays near the Kelvin ball), whereas the two-sphere geometry
    of docs/kelvin/KELVIN_TRANSFORMATION.md 3 puts the physical domain at the
    ORIGIN and uses the offset purely as a meshing translation.  Passing it
    explicitly keeps the choice visible instead of silently inherited.
    """
    from ngsolve import x, y, z

    ox, oy, oz = offset
    cx, cy, cz = phys_center
    dxp, dyp, dzp = x - ox, y - oy, z - oz
    rho2 = dxp * dxp + dyp * dyp + dzp * dzp + 1e-24
    scale = (R_K * R_K) / rho2                    # (R/rho')^2
    return (cx + scale * dxp, cy + scale * dyp, cz + scale * dzp), rho2


def make_kelvin_aware_Omega_s_cf(mesh, Omega_phys_factory, R_K, offset,
                                 phys_center=(0.0, 0.0, 0.0),
                                 kelvin_mats=("kelvin",)):
    """Twisted 0-form pullback of a background SCALAR potential (Convention A).

    *** USE THIS FOR T-OMEGA / Omega-reduced WITH A DECAYING SOURCE ***
    *** (real coils, dipoles -- anything that falls off at infinity) ***

    The magnetic scalar potential is a **twisted 0-form**, and the Kelvin
    inversion is orientation-reversing (`det Dk = -R^6/rho'^6 < 0`), so its
    right-hand representative picks up `s_k = sgn(det Dk) = -1` and NO metric
    factor at all (a 0-form pullback has exponent 0):

        Omega_s'(r') = - Omega_s( k(r') ) ,   k(r') = offset + (R/rho')^2 (r'-offset)

    This is the 0-form entry of the same Convention A family as
    :func:`make_kelvin_aware_A_s_cf`, and it is the rule that makes the
    potential route WORK, because pullback commutes with the exterior
    derivative (`g*(d w) = d(g* w)`).  Verified to machine zero
    (``tests/test_reduced_potential_background.py``): with the matching
    twisted 1-form rule

        H_s'(r') = -(R/rho')^2 * Householder * H_s(k(r'))
                   (radial component +, tangential components -)

    one has EXACTLY

        H_s'  ==  -grad'( Omega_s' ) .

    Contrast :func:`make_reduced_potential_scalar_cf`, which applies the
    1-form Convention B factor `-(rho'/R)^2` to a scalar.  That is not a
    0-form pullback and is not gradient-consistent; measured cost on the
    magnetic-sphere golden: a factor 4/3
    (``validation_test/kelvin_source/test_kelvin_exterior_source_routes.py``).

    REGULARITY.  For a source that DECAYS at infinity the pullback is regular
    at the offset: a dipole `|H_s| ~ 1/r^3` maps to `|H_s'| ~ rho'/R^4`, and
    the potential vanishes quadratically (`Omega_s'(0,0,t) = O(t^2)`).  For a
    background that does NOT decay -- a uniform field applied at infinity --
    `Omega_s'` diverges like `R^2/rho'^2`, because the uniform-field potential
    is genuinely unbounded at infinity.  That is physics, not a defect of the
    rule: there is no bounded 0-form representative in that case, so use the
    1-form route (:func:`make_reduced_potential_background_cf`) instead.

    Args:
        mesh: NGSolve Mesh.
        Omega_phys_factory: callable ``(x_cf, y_cf, z_cf) -> scalar CF``
            giving the PHYSICAL background potential at the given point.
        R_K: Kelvin sphere radius.
        offset: 3-tuple, Kelvin sphere center.
        kelvin_mats: substring(s) used to detect Kelvin materials.

    Returns:
        Scalar CoefficientFunction with
            inner materials:  Omega_s(x, y, z)
            kelvin materials: -Omega_s(k(r'))

    Reference: https://www.ele.kindai.ac.jp/laboratory/sugahara/elemag/geometry09.php
    (twisted-form sign table, `Phi_m = -k* Phi'_m`); docs/kelvin/
    KELVIN_TRANSFORMATION.md 2.3 (0-form exponent 0) and 7.4.
    """
    from ngsolve import x, y, z

    (kx, ky, kz), _ = _kelvin_mapped_coords(R_K, offset, phys_center)

    Omega_inner = Omega_phys_factory(x, y, z)
    Omega_kelvin = -Omega_phys_factory(kx, ky, kz)   # twisted: s_k = -1

    d = {}
    for m in mesh.GetMaterials():
        ml = m.lower()
        is_kelvin = any(kw in ml for kw in kelvin_mats)
        d[m] = Omega_kelvin if is_kelvin else Omega_inner
    return mesh.MaterialCF(d, default=Omega_inner)


def make_kelvin_aware_H_s_cf(mesh, H_phys_factory, R_K, offset,
                             phys_center=(0.0, 0.0, 0.0),
                             kelvin_mats=("kelvin",)):
    """Twisted 1-form pullback of a background FIELD (Convention A).

    Partner of :func:`make_kelvin_aware_Omega_s_cf`.  `H` is a twisted 1-form,
    so on top of the straight 1-form pullback `(R/rho')^2 * Householder` it
    carries `s_k = -1`:

        H_s'(r') = -(R/rho')^2 * (I - 2 n n^T) * H_s(k(r'))

    i.e. the RADIAL component keeps its sign and the TANGENTIAL components
    flip -- the sign table of the geometry09 note.  Use this instead of
    :func:`make_kelvin_aware_A_s_cf` whenever the quantity being pulled back
    is twisted (`H`, `J`, `D`); `A` and `B` are straight and do NOT take the
    extra minus.

    Together with `make_kelvin_aware_Omega_s_cf` this satisfies
    `H_s' = -grad'(Omega_s')` exactly, which is the whole point of using a
    genuine pullback rather than the Convention B engineering formula.

    Args / Returns: as :func:`make_kelvin_aware_Omega_s_cf`, but the factory
    returns a 3-vector CF and the result is a VectorCF.
    """
    from ngsolve import x, y, z, sqrt, CoefficientFunction as CF

    ox, oy, oz = offset
    (kx, ky, kz), rho2 = _kelvin_mapped_coords(R_K, offset, phys_center)
    scale = (R_K * R_K) / rho2                   # (R/rho')^2

    H_inner = H_phys_factory(x, y, z)
    H_at_k = H_phys_factory(kx, ky, kz)

    # Householder about n = (r' - offset)/rho', then the twisted sign s_k = -1.
    rho = sqrt(rho2)
    nx, ny, nz = (x - ox) / rho, (y - oy) / rho, (z - oz) / rho
    h_dot_n = H_at_k[0] * nx + H_at_k[1] * ny + H_at_k[2] * nz
    refl = (H_at_k[0] - 2 * h_dot_n * nx,
            H_at_k[1] - 2 * h_dot_n * ny,
            H_at_k[2] - 2 * h_dot_n * nz)
    factor = -scale                               # -(R/rho')^2
    H_kelvin = tuple(factor * c for c in refl)

    def _switch(kelvin_comp, inner_comp):
        d = {}
        for m in mesh.GetMaterials():
            ml = m.lower()
            is_kelvin = any(kw in ml for kw in kelvin_mats)
            d[m] = kelvin_comp if is_kelvin else inner_comp
        return mesh.MaterialCF(d, default=inner_comp)

    return CF(tuple(_switch(H_kelvin[i], H_inner[i]) for i in range(3)))


def make_reduced_potential_scalar_cf(mesh, Phi_inner_factory, R_K, offset,
                                     kelvin_mats=("kelvin",),
                                     dim=3):
    """Convention B applied to a scalar -- NOT a 0-form pullback.

    *** PREFER make_kelvin_aware_Omega_s_cf FOR T-OMEGA / Omega-reduced. ***

    This applies the 1-form Convention B factor `-(rho'/R)^2` to a scalar.
    It is bounded at the offset, but it is NOT the pullback of a 0-form (a
    0-form pullback carries NO metric factor) and it is NOT gradient-
    consistent with any field rule.  Differentiating it into a field
    overshoots the magnetic-sphere golden by exactly 4/3
    (``validation_test/kelvin_source/test_kelvin_exterior_source_routes.py``,
    route B-0).  It is kept because the T-Omega design note proposes it and
    because its behaviour is contract-locked, but it should not be used to
    drive a weak form.  Use :func:`make_kelvin_aware_Omega_s_cf`.

    The formula, for the record:

        Omega_s'(r') = -(rho'/R)^2 * Omega_s(r' - offset)    (3D)
        Omega_s'(r') = -Omega_s(r' - offset)                 (2D)

    evaluated at LOCAL (offset-relative) coordinates, bounded at the
    offset, sign-flipped at the periodic Kelvin boundary.

    Why it cannot drive a weak form.  Locked by
    ``tests/test_reduced_potential_background.py``: for a uniform
    background ``H_s = H_0 z_hat`` with ``Omega_s = -H_0 z``,

        curl(H_s' from the 1-form Convention B helper)
            = (2 H_0 / R^2) (-y', x', 0)   !=  0

    so the 1-form Convention B exterior field admits no scalar potential
    at all, and

        -grad(Omega_s' from this helper) - H_s' (1-form Convention B)
            = -(2 H_0 / R^2) (x' z', y' z', z'^2)  !=  0 .

    Neither of those is a defect of the helpers -- Convention B is an
    engineering formula, not a pullback, so its 0-form and 1-form flavours
    are simply unrelated.  The genuine pullback family
    (:func:`make_kelvin_aware_Omega_s_cf` / :func:`make_kelvin_aware_H_s_cf`)
    IS gradient-consistent, to machine zero, because pullback commutes with
    the exterior derivative.

    Summary of which helper drives which weak form:

    ==============================  ===================================
    situation                       helper
    ==============================  ===================================
    background FIELD, non-decaying  make_reduced_potential_background_cf
    background POTENTIAL, decaying  make_kelvin_aware_Omega_s_cf
    background FIELD, decaying      make_kelvin_aware_H_s_cf
    background POTENTIAL, uniform   (none exists -- unbounded at infinity;
                                     use the 1-form route)
    ==============================  ===================================

    Args:
        mesh: NGSolve Mesh.
        Phi_inner_factory: callable ``(x_cf, y_cf, z_cf) -> scalar CF``
            returning the background potential at the given coordinates.
            Receives GLOBAL coords for inner-region values and LOCAL
            (offset-relative) coords for Kelvin-region values.  For a
            uniform ``H_0 z_hat`` applied at infinity:
            ``lambda xc, yc, zc: -H_0 * zc``.
        R_K: Kelvin sphere radius.
        offset: 3-tuple, Kelvin sphere center.
        kelvin_mats: substring(s) used to detect Kelvin materials.
        dim: 2 or 3 (default 3).  Selects between 2D and 3D scaling.

    Returns:
        Scalar CoefficientFunction with:
            inner materials:  Omega_s(x, y, z)                       (global)
            kelvin materials: -(rho'/R)^2 * Omega_s(x-ox, y-oy, z-oz) (3D, local)
                              -Omega_s(x-ox, y-oy, z-oz)             (2D, local)

    Example::

        # Uniform H_0 z_hat applied at infinity, T-Omega background:
        Omega_s_cf = make_reduced_potential_scalar_cf(
            mesh,
            lambda xc, yc, zc: -H_0 * zc,
            R_K=R_K, offset=offset, kelvin_mats=("kelvin",))

    See:
        docs/kelvin/KELVIN_TRANSFORMATION.md 7.4 (rule), 7.6 (A vs B)
        validation_test/maglev/research_cln/ngsolve_validation/
            cuboid_521_T_Omega_Kelvin_design.md (the T-Omega route that
            requires this 0-form variant)
    """
    from ngsolve import x, y, z

    if dim not in (2, 3):
        raise ValueError(f"dim must be 2 or 3, got {dim}")

    ox, oy, oz = offset
    # Inner region: evaluate at global coordinates
    Phi_inner = Phi_inner_factory(x, y, z)
    # Kelvin region: evaluate at local (offset-relative) coordinates
    Phi_local = Phi_inner_factory(x - ox, y - oy, z - oz)

    if dim == 3:
        rho2 = ((x - ox) ** 2 + (y - oy) ** 2 + (z - oz) ** 2 + 1e-24)
        kelvin_factor = -rho2 / (R_K * R_K)            # -(rho'/R)^2
    else:  # dim == 2
        kelvin_factor = -1.0                            # 2D: sign flip only

    Phi_kelvin = kelvin_factor * Phi_local

    d = {}
    for m in mesh.GetMaterials():
        ml = m.lower()
        is_kelvin = any(kw in ml for kw in kelvin_mats)
        d[m] = Phi_kelvin if is_kelvin else Phi_inner
    return mesh.MaterialCF(d, default=Phi_inner)
