"""kelvin_material.py

Layer 2 of the Kelvin helper API: build NGSolve CoefficientFunctions
for the Kelvin-modulated material parameter (nu) and for an external
A_s source field with the Kelvin pullback automatically applied in the
Kelvin exterior domain.

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

Solution pullback (A 1-form, DIFFERENT from material factor):
    A_comp(r') = (R/rho')^2 * H * A_phys(r_phys)
    B_comp(r') = -(R/rho')^4 * H * B_phys(r_phys)    (2-form pseudovector)

See examples/kelvin_transformation/CONVENTION.md for the declaration
and examples/kelvin_transformation/docs/pullback_derivation_3D.md sec 8
for the derivation.
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
    """Build an A_s vector CF that handles inner / Kelvin domains.

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
