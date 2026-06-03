"""
Analytical magnetic field computation for cylindrical permanent magnets.

This module provides analytical formulas for computing the B-field of
uniformly magnetized cylindrical magnets using elliptic integrals.

The implementation is based on:
    [1] Derby, N., Olbert, S., "Cylindrical Magnets and Ideal Solenoids",
        American Journal of Physics, Vol. 78(3), pp. 229-235, 2010.
    [2] Caciagli, A., et al., "Exact expression for the magnetic field of
        a finite cylinder with arbitrary uniform magnetization",
        Journal of Magnetism and Magnetic Materials, 456, 423-432, 2018.
    [3] Bulirsch, R., "Numerical Calculation of Elliptic Integrals and
        Elliptic Functions. III", Numerische Mathematik 13, 305-315, 1969.

Usage:
    # Create a cylinder field source
    from radia.cylindrical_magnet import CylindricalMagnet

    # Cylinder with radius 10mm, height 20mm, centered at origin
    # Magnetization: 1.2 T equivalent in z-direction
    cyl = CylindricalMagnet(
        center=[0, 0, 0],
        radius=10.0,     # mm
        height=20.0,     # mm
        magnetization=[0, 0, 955000]  # A/m (Mz = 1.2T / mu0)
    )

    # Get field at point
    B = cyl.get_B([5, 0, 15])  # Returns [Bx, By, Bz] in Tesla

    # Use as background field with Radia
    import radia as rad
    bckg = rad.ObjBckg(cyl.get_B)
"""

import numpy as np
from typing import List, Tuple

# Physical constants
PI = np.pi
MU0 = 4.0 * np.pi * 1.0e-7  # H/m


def _cel(kc: float, p: float, c: float, s: float) -> float:
    """
    Bulirsch's generalized complete elliptic integral cel(kc, p, c, s).

    This function computes the generalized complete elliptic integral
    as defined by Bulirsch (1969). The standard elliptic integrals
    can be expressed in terms of cel:

        K(m) = cel(sqrt(1-m), 1, 1, 1)
        E(m) = cel(sqrt(1-m), 1, 1, 1-m)
        Pi(n, m) = cel(sqrt(1-m), 1-n, 1, 1)

    Parameters:
        kc: complementary modulus, kc = sqrt(1 - k^2), kc != 0
        p:  parameter
        c:  parameter
        s:  parameter

    Reference: Kirby2009, based on Bulirsch1969
    """
    if kc == 0.0:
        raise ValueError("cel: kc=0 not allowed")

    errtol = 1.0e-8
    k = abs(kc)
    pp = p
    cc = c
    ss = s
    em = 1.0

    if p > 0.0:
        pp = np.sqrt(p)
        ss = s / pp
    else:
        f = kc * kc
        q = 1.0 - f
        g = 1.0 - pp
        f = f - pp
        q = q * (ss - c * pp)
        pp = np.sqrt(f / g)
        cc = (c - ss) / g
        ss = -q / (g * g * pp) + cc * pp

    f = cc
    cc = cc + ss / pp
    g = k / pp
    ss = 2.0 * (ss + f * g)
    pp = g + pp
    g = em
    em = k + em
    kk = k

    while abs(g - k) > g * errtol:
        k = 2.0 * np.sqrt(kk)
        kk = k * em
        f = cc
        cc = cc + ss / pp
        g = kk / pp
        ss = 2.0 * (ss + f * g)
        pp = g + pp
        g = em
        em = k + em

    return (PI / 2.0) * (ss + cc * em) / (em * (em + pp))


def _axial_cylinder_bfield(rho: float, z: float, R: float, L: float, Mz: float
                           ) -> Tuple[float, float]:
    """
    B-field of axially magnetized cylinder at point (rho, z).

    Uses Derby-Olbert (2010) formulation.

    Parameters:
        rho: radial distance from axis [mm]
        z:   axial position [mm]
        R:   cylinder radius [mm]
        L:   cylinder half-height [mm]
        Mz:  axial magnetization [A/m]

    Returns:
        (Brho, Bz): field components in Tesla
    """
    if R <= 0.0 or L <= 0.0:
        return 0.0, 0.0

    # Polarization J = mu0 * M [T]
    J = MU0 * Mz

    # Make dimensionless by dividing by R
    r = rho / R
    zn = z / R
    z0 = L / R

    # Important quantities
    zph = zn + z0
    zmh = zn - z0
    dpr = 1.0 + r
    dmr = 1.0 - r

    eps = 1.0e-15
    if dpr < eps:
        dpr = eps

    sq0 = np.sqrt(zmh * zmh + dpr * dpr)
    sq1 = np.sqrt(zph * zph + dpr * dpr)

    k0_sq = (zmh * zmh + dmr * dmr) / (zmh * zmh + dpr * dpr)
    k1_sq = (zph * zph + dmr * dmr) / (zph * zph + dpr * dpr)
    k0 = np.sqrt(k0_sq)
    k1 = np.sqrt(k1_sq)
    gamma = dmr / dpr
    gamma2 = gamma * gamma

    # On-axis case
    if r < 1.0e-10:
        Brho = 0.0
        Bz = J / 2.0 * (zph / sq1 - zmh / sq0)
        return Brho, Bz

    # Radial field component
    Br_norm = (_cel(k1, 1.0, 1.0, -1.0) / sq1 - _cel(k0, 1.0, 1.0, -1.0) / sq0) / PI

    # Axial field component
    Bz_norm = (1.0 / dpr) * (
        zph * _cel(k1, gamma2, 1.0, gamma) / sq1 -
        zmh * _cel(k0, gamma2, 1.0, gamma) / sq0
    ) / PI

    Brho = J * Br_norm
    Bz = J * Bz_norm
    return Brho, Bz


def _diametric_cylinder_hfield(x: float, y: float, z: float,
                                R: float, L: float,
                                Mx: float, My: float
                                ) -> Tuple[float, float, float]:
    """
    H-field of diametrically magnetized cylinder at point (x, y, z).

    Uses Caciagli (2018) formulation.

    Parameters:
        x, y, z: observation point [mm]
        R:       cylinder radius [mm]
        L:       cylinder half-height [mm]
        Mx, My:  transverse magnetization [A/m]

    Returns:
        (Hx, Hy, Hz): field components in A/m
    """
    if R <= 0.0 or L <= 0.0:
        return 0.0, 0.0, 0.0

    # Convert to cylindrical coordinates
    r = np.sqrt(x * x + y * y)
    phi = np.arctan2(y, x)

    # Magnetization direction
    M_mag = np.sqrt(Mx * Mx + My * My)
    if M_mag < 1.0e-15:
        return 0.0, 0.0, 0.0
    tetta = np.arctan2(My, Mx)
    phi_rel = phi - tetta

    # Make dimensionless
    rn = r / R
    zn = z / R
    z0 = L / R

    zp = zn + z0
    zm = zn - z0
    zp2 = zp * zp
    zm2 = zm * zm
    r2 = rn * rn

    # Small r case: Taylor series
    if rn < 0.05:
        zpp = zp2 + 1.0
        zmm = zm2 + 1.0
        sqrt_p = np.sqrt(zpp)
        sqrt_m = np.sqrt(zmm)

        frac1 = zp / sqrt_p
        frac2 = zm / sqrt_m

        r3 = r2 * rn
        r4 = r3 * rn
        r5 = r4 * rn

        term1 = frac1 - frac2
        term2 = (frac1 / (zpp ** 2) - frac2 / (zmm ** 2)) * r2 / 8.0
        term3 = (
            (3.0 - 4.0 * zp2) * frac1 / (zpp ** 4) -
            (3.0 - 4.0 * zm2) * frac2 / (zmm ** 4)
        ) / 64.0 * r4

        Hr = -np.cos(phi_rel) / 4.0 * (term1 + 9.0 * term2 + 25.0 * term3)
        Hphi = np.sin(phi_rel) / 4.0 * (term1 + 3.0 * term2 + 5.0 * term3)
        Hz_cyl = -np.cos(phi_rel) / 4.0 * (
            rn * (1.0 / (zpp * sqrt_p) - 1.0 / (zmm * sqrt_m)) +
            3.0 / 8.0 * r3 * (
                (1.0 - 4.0 * zp2) / (zpp ** 3 * sqrt_p) -
                (1.0 - 4.0 * zm2) / (zmm ** 3 * sqrt_m)
            ) +
            15.0 / 64.0 * r5 * (
                (1.0 - 12.0 * zp2 + 8.0 * zp2 ** 2) / (zpp ** 5 * sqrt_p) -
                (1.0 - 12.0 * zm2 + 8.0 * zm2 ** 2) / (zmm ** 5 * sqrt_m)
            )
        )
    else:
        # General case
        from scipy.special import ellipe, ellipk

        rp = rn + 1.0
        rm = rn - 1.0
        rp2 = rp * rp
        rm2 = rm * rm

        ap2 = zp2 + rm2
        am2 = zm2 + rm2
        ap = np.sqrt(ap2)
        am = np.sqrt(am2)

        argp = -4.0 * rn / ap2
        argm = -4.0 * rn / am2

        eps = 1.0e-15
        if abs(rm) < eps:
            argc = -1.0e16
            one_over_rm = 0.0
        else:
            argc = -4.0 * rn / rm2
            one_over_rm = 1.0 / rm

        elle_p = ellipe(argp)
        elle_m = ellipe(argm)
        ellk_p = ellipk(argp)
        ellk_m = ellipk(argm)

        ellpi_p = _cel(np.sqrt(1.0 - argp), 1.0 - argc, 1.0, 1.0)
        ellpi_m = _cel(np.sqrt(1.0 - argm), 1.0 - argc, 1.0, 1.0)

        Hr = -np.cos(phi_rel) / (4.0 * PI * r2) * (
            -zm * am * elle_m +
            zp * ap * elle_p +
            zm / am * (2.0 + zm2) * ellk_m -
            zp / ap * (2.0 + zp2) * ellk_p +
            (zm / am * ellpi_m - zp / ap * ellpi_p) * rp * (r2 + 1.0) * one_over_rm
        )

        Hphi = np.sin(phi_rel) / (4.0 * PI * r2) * (
            +zm * am * elle_m -
            zp * ap * elle_p -
            zm / am * (2.0 + zm2 + 2.0 * r2) * ellk_m +
            zp / ap * (2.0 + zp2 + 2.0 * r2) * ellk_p +
            zm / am * rp2 * ellpi_m -
            zp / ap * rp2 * ellpi_p
        )

        Hz_cyl = -np.cos(phi_rel) / (2.0 * PI * rn) * (
            +am * elle_m -
            ap * elle_p -
            (1.0 + zm2 + r2) / am * ellk_m +
            (1.0 + zp2 + r2) / ap * ellk_p
        )

    # Convert cylindrical to Cartesian
    cos_phi = np.cos(phi)
    sin_phi = np.sin(phi)

    Hx_unit = Hr * cos_phi - Hphi * sin_phi
    Hy_unit = Hr * sin_phi + Hphi * cos_phi

    Hx = M_mag * Hx_unit
    Hy = M_mag * Hy_unit
    Hz = M_mag * Hz_cyl
    return Hx, Hy, Hz


def _current_loop_vector_potential(rho: float, z: float, a: float, I: float) -> float:
    """
    Vector potential A_phi of a circular current loop.

    Uses elliptic integral formulation from Ortner et al. (2025).

    A_phi = (mu0*I / 2*pi) * sqrt((a+rho)^2 + z^2) / rho
            * [ (a^2 + rho^2 + z^2) / ((a+rho)^2 + z^2) * K(kappa) - E(kappa) ]

    where kappa = sqrt(4*a*rho / ((a+rho)^2 + z^2))

    Parameters:
        rho: radial distance from axis [mm]
        z:   axial position relative to loop [mm]
        a:   loop radius [mm]
        I:   current [A]

    Returns:
        A_phi: azimuthal vector potential [T*mm]
    """
    from scipy.special import ellipe, ellipk

    # On-axis singularity: A_phi = 0 by symmetry
    if rho < 1.0e-15:
        return 0.0

    if a <= 0.0:
        return 0.0

    # Distance factors
    apr = a + rho
    apr2 = apr * apr
    z2 = z * z

    denom = apr2 + z2
    if denom < 1.0e-30:
        return 0.0

    # Elliptic integral argument kappa^2
    kappa_sq = 4.0 * a * rho / denom
    kappa = np.sqrt(kappa_sq)

    # Avoid singularity when kappa -> 1 (on the loop)
    if kappa >= 0.9999:
        kappa = 0.9999
        kappa_sq = kappa * kappa

    # Elliptic integrals (scipy uses m = kappa^2 as argument)
    K_val = ellipk(kappa_sq)
    E_val = ellipe(kappa_sq)

    # Formula components
    sqrt_denom = np.sqrt(denom)
    bracket = (a * a + rho * rho + z2) / denom * K_val - E_val

    # A_phi = (mu0 * I / 2*pi) * sqrt(denom) / rho * bracket
    A_phi = (MU0 * I / (2.0 * PI)) * sqrt_denom / rho * bracket

    return A_phi


def _axial_cylinder_vector_potential(rho: float, z: float, R: float, L: float,
                                      Mz: float, n_gauss: int = 20,
                                      length_to_m: float = 0.001) -> float:
    """
    Vector potential A_phi of an axially magnetized cylinder.

    For an axially magnetized cylinder, the equivalent surface current
    density on the curved surface is K_phi = Mz [A/m].

    The total A_phi is obtained by integrating over the height:
    A_phi = integral from -L to +L of A_phi_loop(z') dz'

    where A_phi_loop is the vector potential of a current loop at height z'
    with current dI = Mz * dz' (per unit length in m).

    Parameters:
        rho: radial distance from axis [mm]
        z:   axial position [mm]
        R:   cylinder radius [mm]
        L:   cylinder half-height [mm]
        Mz:  axial magnetization [A/m]
        n_gauss: number of Gauss quadrature points

    Returns:
        A_phi: azimuthal vector potential in T*user_length_unit

    Unit Analysis:
        - Surface current density: K = Mz [A/m]
        - Current element: dI = K * dz' [m] = Mz * dz * length_to_m [A]
        - _current_loop_vector_potential returns A in T*m (SI)
        - Result converted to user units: A_phi / length_to_m
    """
    from numpy.polynomial.legendre import leggauss

    # On-axis: A_phi = 0 by symmetry
    if rho < 1.0e-15:
        return 0.0

    if R <= 0.0 or L <= 0.0:
        return 0.0

    # Gauss-Legendre quadrature
    xi, wi = leggauss(n_gauss)

    # Transform from [-1, 1] to [-L, +L]
    # z' = L * xi, dz' = L * d(xi)
    A_phi = 0.0
    for i in range(n_gauss):
        z_prime = L * xi[i]
        dz = L * wi[i]  # dz is in mm

        # Surface current density K = Mz [A/m]
        # Current element dI = K * dz' [m] = Mz * dz * length_to_m [A]
        dI = Mz * dz * length_to_m  # [A]

        # Vector potential from this loop element
        A_loop = _current_loop_vector_potential(rho, z - z_prime, R, dI)
        A_phi += A_loop

    # _current_loop_vector_potential returns A in T*m (SI).
    # Convert to user length units: T*m / length_to_m = T*user_unit
    return A_phi / length_to_m


class CylindricalMagnet:
    """
    Analytical cylindrical permanent magnet field source.

    This class computes the magnetic field of a uniformly magnetized
    cylinder using analytical formulas based on elliptic integrals.

    Attributes:
        center: [x, y, z] center position
        radius: cylinder radius
        height: cylinder height (half-height = height/2)
        magnetization: [Mx, My, Mz] magnetization in A/m
        axis: 'x', 'y', or 'z' - cylinder axis direction

    Note: Magnetization in A/m. For NdFeB with Br=1.2T:
          M = Br / mu0 = 1.2 / (4*pi*1e-7) ~ 955000 A/m
    """

    def __init__(self, center: List[float], radius: float, height: float,
                 magnetization: List[float], axis: str = 'z',
                 units: str = 'mm'):
        """
        Initialize cylindrical magnet.

        Parameters:
            center: [x, y, z] center position
            radius: cylinder radius
            height: total height
            magnetization: [Mx, My, Mz] in A/m
            axis: cylinder axis direction ('x', 'y', or 'z')
            units: length unit for geometry ('mm' or 'm')
        """
        self._length_to_m = 0.001 if units == 'mm' else 1.0
        self.center = np.array(center, dtype=float)
        self.radius = float(radius)
        self.height = float(height)
        self.half_height = height / 2.0
        self.magnetization = np.array(magnetization, dtype=float)
        self.axis = axis.lower()

        if self.axis not in ['x', 'y', 'z']:
            raise ValueError("axis must be 'x', 'y', or 'z'")

    def get_B(self, point: List[float]) -> List[float]:
        """
        Get magnetic field B at observation point.

        Parameters:
            point: [x, y, z] observation point (in constructor length units)

        Returns:
            [Bx, By, Bz] magnetic field in Tesla
        """
        # Transform to local coordinates (cylinder at origin, axis along z)
        p = np.array(point, dtype=float) - self.center

        # Rotate to align cylinder axis with z-axis
        if self.axis == 'z':
            p_local = p.copy()
            M_local = self.magnetization.copy()
        elif self.axis == 'x':
            # Rotate: x->z, y->x, z->y
            p_local = np.array([p[1], p[2], p[0]])
            M_local = np.array([self.magnetization[1], self.magnetization[2], self.magnetization[0]])
        else:  # axis == 'y'
            # Rotate: y->z, z->x, x->y
            p_local = np.array([p[2], p[0], p[1]])
            M_local = np.array([self.magnetization[2], self.magnetization[0], self.magnetization[1]])

        # Compute field in local coordinates
        rho = np.sqrt(p_local[0]**2 + p_local[1]**2)
        z_local = p_local[2]

        Bx_local, By_local, Bz_local = 0.0, 0.0, 0.0

        # Axial magnetization contribution (Mz)
        Mz_local = M_local[2]
        if abs(Mz_local) > 1.0e-15:
            Brho, Bz = _axial_cylinder_bfield(rho, z_local, self.radius, self.half_height, Mz_local)
            if rho > 1.0e-15:
                Bx_local += Brho * p_local[0] / rho
                By_local += Brho * p_local[1] / rho
            Bz_local += Bz

        # Transverse magnetization contribution (Mx, My)
        Mx_local = M_local[0]
        My_local = M_local[1]
        if abs(Mx_local) > 1.0e-15 or abs(My_local) > 1.0e-15:
            Hx, Hy, Hz = _diametric_cylinder_hfield(
                p_local[0], p_local[1], p_local[2],
                self.radius, self.half_height,
                Mx_local, My_local
            )
            Bx_local += MU0 * Hx
            By_local += MU0 * Hy
            Bz_local += MU0 * Hz

        # Rotate back to global coordinates
        if self.axis == 'z':
            Bx, By, Bz = Bx_local, By_local, Bz_local
        elif self.axis == 'x':
            # Rotate back: z->x, x->y, y->z
            Bx, By, Bz = Bz_local, Bx_local, By_local
        else:  # axis == 'y'
            # Rotate back: z->y, x->z, y->x
            Bx, By, Bz = By_local, Bz_local, Bx_local

        return [Bx, By, Bz]

    def get_H(self, point: List[float]) -> List[float]:
        """
        Get H-field at observation point.

        For points outside the magnet: H = B / mu0
        For points inside: this is only approximate.

        Parameters:
            point: [x, y, z] observation point (in constructor length units)

        Returns:
            [Hx, Hy, Hz] H-field in A/m
        """
        B = self.get_B(point)
        return [B[0] / MU0, B[1] / MU0, B[2] / MU0]

    def get_A(self, point: List[float], n_gauss: int = 20) -> List[float]:
        """
        Get vector potential A at observation point.

        The vector potential has only an azimuthal component A_phi for
        axially magnetized cylinders. For diametric (transverse)
        magnetization, A is not implemented and get_A RAISES (No-Fallback:
        it does not silently return a wrong / zero A).

        Parameters:
            point: [x, y, z] observation point (in constructor length units)
            n_gauss: number of Gauss quadrature points for integration

        Returns:
            [Ax, Ay, Az] vector potential in T*length_unit
        """
        # Transform to local coordinates (cylinder at origin, axis along z)
        p = np.array(point, dtype=float) - self.center

        # Rotate to align cylinder axis with z-axis
        if self.axis == 'z':
            p_local = p.copy()
            M_local = self.magnetization.copy()
        elif self.axis == 'x':
            p_local = np.array([p[1], p[2], p[0]])
            M_local = np.array([self.magnetization[1], self.magnetization[2], self.magnetization[0]])
        else:  # axis == 'y'
            p_local = np.array([p[2], p[0], p[1]])
            M_local = np.array([self.magnetization[2], self.magnetization[0], self.magnetization[1]])

        # Diametric (transverse) magnetization contributes to A through a volume
        # integral that is NOT implemented here.  Raise rather than SILENTLY
        # return an axial-only (wrong) A or zero -- a plausible-but-wrong A is
        # worse than a clear error (No-Fallback).  get_B()/get_H() are
        # implemented for all magnetization directions.
        if abs(M_local[0]) > 1.0e-15 or abs(M_local[1]) > 1.0e-15:
            raise NotImplementedError(
                "CylindricalMagnet.get_A() is implemented for AXIAL "
                "magnetization only; this magnet has a diametric (transverse "
                "Mx/My) component whose vector-potential contribution would be "
                "silently dropped. Use get_B() / get_H() (implemented for all "
                "directions), or implement the diametric A volume integral.")

        # Compute A in local coordinates
        rho = np.sqrt(p_local[0]**2 + p_local[1]**2)
        z_local = p_local[2]
        phi = np.arctan2(p_local[1], p_local[0])

        Ax_local, Ay_local, Az_local = 0.0, 0.0, 0.0

        # Axial magnetization contribution (Mz)
        Mz_local = M_local[2]
        if abs(Mz_local) > 1.0e-15:
            A_phi = _axial_cylinder_vector_potential(
                rho, z_local, self.radius, self.half_height, Mz_local, n_gauss,
                self._length_to_m
            )
            # Convert A_phi to Cartesian: Ax = -A_phi * sin(phi), Ay = A_phi * cos(phi)
            Ax_local = -A_phi * np.sin(phi)
            Ay_local = A_phi * np.cos(phi)
            # Az = 0 for axial magnetization

        # (Diametric magnetization was rejected above with a clear raise --
        #  its A would require a volume integral, not implemented.)

        # Rotate back to global coordinates
        if self.axis == 'z':
            Ax, Ay, Az = Ax_local, Ay_local, Az_local
        elif self.axis == 'x':
            Ax, Ay, Az = Az_local, Ax_local, Ay_local
        else:  # axis == 'y'
            Ax, Ay, Az = Ay_local, Az_local, Ax_local

        return [Ax, Ay, Az]

    def __call__(self, point: List[float]) -> List[float]:
        """
        Callable interface for use with rad.ObjBckg().

        Parameters:
            point: [x, y, z] observation point (in constructor length units)

        Returns:
            [Bx, By, Bz] magnetic field in Tesla
        """
        return self.get_B(point)


class RingMagnet(CylindricalMagnet):
    """
    Analytical ring magnet (hollow cylinder) field source.

    This class computes the magnetic field of a uniformly magnetized
    hollow cylinder (ring magnet) by superposition of outer and inner
    solid cylinders.
    """

    def __init__(self, center: List[float], inner_radius: float, outer_radius: float,
                 height: float, magnetization: List[float], axis: str = 'z',
                 units: str = 'mm'):
        """
        Initialize ring magnet.

        Parameters:
            center: [x, y, z] center position
            inner_radius: inner radius (must be positive)
            outer_radius: outer radius
            height: total height
            magnetization: [Mx, My, Mz] in A/m
            axis: cylinder axis direction ('x', 'y', or 'z')
            units: length unit for geometry ('mm' or 'm')
        """
        if inner_radius <= 0:
            raise ValueError("inner_radius must be positive")
        if inner_radius >= outer_radius:
            raise ValueError("inner_radius must be less than outer_radius")
        super().__init__(center, outer_radius, height, magnetization, axis, units)
        self.inner_radius = float(inner_radius)
        self.outer_radius = float(outer_radius)

    def get_B(self, point: List[float]) -> List[float]:
        """
        Get magnetic field B at observation point.

        Computed as B_outer - B_inner (superposition).

        Parameters:
            point: [x, y, z] observation point (in constructor length units)

        Returns:
            [Bx, By, Bz] magnetic field in Tesla
        """
        saved_radius = self.radius
        try:
            # Field from outer cylinder
            self.radius = self.outer_radius
            B_outer = super().get_B(point)

            # Field from inner cylinder
            self.radius = self.inner_radius
            B_inner = super().get_B(point)
        finally:
            self.radius = saved_radius

        return [
            B_outer[0] - B_inner[0],
            B_outer[1] - B_inner[1],
            B_outer[2] - B_inner[2]
        ]
