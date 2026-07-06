"""
Analytical magnetic field computation for permanent magnets.

This module provides analytical formulas for computing the magnetic field of
uniformly magnetized permanent magnets of various shapes:
- SphericalMagnet: Homogeneously magnetized sphere (dipole field outside)
- CuboidMagnet: Rectangular block magnet (Yang/Camacho formulation)
- CurrentLoop: Circular current loop (Ortner formulation)

References:
    [1] Derby, N., Olbert, S., "Cylindrical Magnets and Ideal Solenoids",
        American Journal of Physics, Vol. 78(3), pp. 229-235, 2010.
    [2] Caciagli, A., et al., "Exact expression for the magnetic field of
        a finite cylinder with arbitrary uniform magnetization",
        Journal of Magnetism and Magnetic Materials, 456, 423-432, 2018.
    [3] Yang, Z.J., et al., "Potential and force between a magnet and
        a bulk Y1Ba2Cu3O7 superconductor studied by a mechanical pendulum",
        Supercond. Sci. Technol. 3(12):591, 1990.
    [4] Engel-Herbert, R., Hesjedal, T., "Calculation of the magnetic stray
        field of a uniaxial magnetic domain",
        J. Appl. Phys. 97(7):074504-4, 2005.
    [5] Camacho, J.M., Sosa, V., "Alternative method to calculate the magnetic
        field of permanent magnets with azimuthal symmetry",
        Rev. Mex. Fis. E 59 (2013) 8-17.
    [6] Ortner, M., et al., "Numerically stable and computationally efficient
        expression for the magnetic field of a current loop",
        Magnetism 2023, 3(1), 11-31.

Usage:
    from radia.analytical_magnet import SphericalMagnet, CuboidMagnet, CurrentLoop

    # Spherical magnet
    sphere = SphericalMagnet(center=[0,0,0], diameter=20.0, magnetization=[0,0,955000])
    B = sphere.get_B([10, 0, 0])

    # Cuboid magnet
    cuboid = CuboidMagnet(center=[0,0,0], dimensions=[20,20,10], magnetization=[0,0,955000])
    B = cuboid.get_B([15, 0, 0])

    # Current loop
    loop = CurrentLoop(center=[0,0,0], diameter=50.0, current=100.0)
    B = loop.get_B([0, 0, 25])
"""

import numpy as np
from typing import List, Tuple

# Physical constants
PI = np.pi
MU0 = 4.0 * np.pi * 1.0e-7  # H/m


# =============================================================================
# Spherical Magnet
# =============================================================================

class SphericalMagnet:
    """
    Analytical spherical permanent magnet field source.

    A uniformly magnetized sphere produces:
    - Outside: A magnetic dipole field
    - Inside: A uniform field B = (2/3) * mu0 * M

    Attributes:
        center: [x, y, z] center position
        diameter: sphere diameter
        magnetization: [Mx, My, Mz] magnetization in A/m

    Note: Magnetization in A/m. For NdFeB with Br=1.2T:
          M = Br / mu0 = 1.2 / (4*pi*1e-7) ~ 955000 A/m
    """

    def __init__(self, center: List[float], diameter: float,
                 magnetization: List[float], units: str = 'mm'):
        """
        Initialize spherical magnet.

        Parameters:
            center: [x, y, z] center position
            diameter: sphere diameter
            magnetization: [Mx, My, Mz] in A/m
            units: length unit for geometry ('mm' or 'm')
        """
        self._length_to_m = 0.001 if units == 'mm' else 1.0
        self.center = np.array(center, dtype=float)
        self.diameter = float(diameter)
        self.radius = diameter / 2.0
        self.magnetization = np.array(magnetization, dtype=float)
        # Polarization J = mu0 * M [T]
        self.polarization = MU0 * self.magnetization

    def get_B(self, point: List[float]) -> List[float]:
        """
        Get magnetic field B at observation point.

        Parameters:
            point: [x, y, z] observation point (in constructor length units)

        Returns:
            [Bx, By, Bz] magnetic field in Tesla
        """
        # Position vector from center to observation point
        r_vec = np.array(point, dtype=float) - self.center
        r = np.linalg.norm(r_vec)

        # Check if inside sphere
        if r < self.radius:
            # Inside: B = (2/3) * J = (2/3) * mu0 * M
            return list(self.polarization * 2.0 / 3.0)

        # Outside: dipole field
        # B = (mu0 / 4*pi) * (3*(m.r)*r/r^5 - m/r^3)
        # where m = (4/3)*pi*R^3 * M is the magnetic moment
        # This simplifies to:
        # B = (J * R^3 / 3) * (3*(M_hat.r_hat)*r_hat - M_hat) / r^3
        # Or using polarization directly:
        # B = (R^3 / (3*r^3)) * (3*(J.r_hat)*r_hat - J)

        r_hat = r_vec / r
        J_dot_r = np.dot(self.polarization, r_hat)
        factor = (self.radius ** 3) / (3.0 * r ** 3)

        B = factor * (3.0 * J_dot_r * r_hat - self.polarization)
        return list(B)

    def get_H(self, point: List[float]) -> List[float]:
        """
        Get H-field at observation point.

        Parameters:
            point: [x, y, z] observation point (in constructor length units)

        Returns:
            [Hx, Hy, Hz] H-field in A/m
        """
        B = np.array(self.get_B(point))
        r_vec = np.array(point, dtype=float) - self.center
        r = np.linalg.norm(r_vec)

        if r < self.radius:
            # Inside: H = B/mu0 - M
            return list(B / MU0 - self.magnetization)
        else:
            # Outside: H = B/mu0
            return list(B / MU0)

    def get_A(self, point: List[float]) -> List[float]:
        """
        Get vector potential A at observation point.

        For a uniformly magnetized sphere:
        - Outside: A = (mu0/4*pi) * (m x r) / r^3, where m = (4/3)*pi*R^3 * M
        - Inside: A = (mu0/3) * (M x r)

        Parameters:
            point: [x, y, z] observation point (in constructor length units)

        Returns:
            [Ax, Ay, Az] vector potential in T*m (SI units)
        """
        # Convert to meters for SI calculation
        r_vec_m = (np.array(point, dtype=float) - self.center) * self._length_to_m  # [m]
        r_m = np.linalg.norm(r_vec_m)  # [m]
        radius_m = self.radius * self._length_to_m  # [m]

        # Magnetic moment: m = (4/3)*pi*R^3 * M [A*m^2]
        m = (4.0 / 3.0) * PI * (radius_m ** 3) * self.magnetization

        if r_m < radius_m:
            # Inside sphere: A = (mu0/3) * (M x r)
            A = (MU0 / 3.0) * np.cross(self.magnetization, r_vec_m)
        else:
            # Outside sphere: dipole field A = (mu0/4*pi) * (m x r) / r^3
            if r_m < 1e-15:
                A = np.array([0.0, 0.0, 0.0])
            else:
                A = (MU0 / (4.0 * PI)) * np.cross(m, r_vec_m) / (r_m ** 3)

        return list(A)

    def __call__(self, point: List[float]) -> List[float]:
        """Callable interface for use with rad.ObjBckg()."""
        return self.get_B(point)


# =============================================================================
# Cuboid Magnet
# =============================================================================

def _cuboid_bfield_component(x: float, y: float, z: float,
                              a: float, b: float, c: float,
                              Jx: float, Jy: float, Jz: float) -> Tuple[float, float, float]:
    """
    B-field of a cuboid magnet at point (x, y, z).

    Uses the Yang/Engel-Herbert/Camacho formulation based on magnetic surface charges.
    Following magpylib's implementation for numerical stability.

    Parameters:
        x, y, z: observation point relative to cuboid center [mm]
        a, b, c: half-dimensions of cuboid [mm]
        Jx, Jy, Jz: polarization components [T]

    Returns:
        (Bx, By, Bz) in Tesla
    """
    # Map to "bottomQ4" quadrant for numerical stability (Cichon 2019)
    # bottomQ4 has: x > 0, y < 0, z < 0
    maskx = x < 0
    masky = y > 0
    maskz = z > 0

    # Create sign matrix for field transformation
    qs = np.ones((3, 3))
    qs_flipx = np.array([[1, -1, -1], [-1, 1, 1], [-1, 1, 1]])
    qs_flipy = np.array([[1, -1, 1], [-1, 1, -1], [1, -1, 1]])
    qs_flipz = np.array([[1, 1, -1], [1, 1, -1], [-1, -1, 1]])

    if maskx:
        x = -x
        qs = qs * qs_flipx
    if masky:
        y = -y
        qs = qs * qs_flipy
    if maskz:
        z = -z
        qs = qs * qs_flipz

    # Vertex distances
    xma, xpa = x - a, x + a
    ymb, ypb = y - b, y + b
    zmc, zpc = z - c, z + c

    xma2, xpa2 = xma ** 2, xpa ** 2
    ymb2, ypb2 = ymb ** 2, ypb ** 2
    zmc2, zpc2 = zmc ** 2, zpc ** 2

    # Distance to 8 vertices
    r000 = np.sqrt(xma2 + ymb2 + zmc2)
    pmp = np.sqrt(xpa2 + ymb2 + zpc2)
    pmm = np.sqrt(xpa2 + ymb2 + zmc2)
    mmp = np.sqrt(xma2 + ymb2 + zpc2)
    mpm = np.sqrt(xma2 + ypb2 + zmc2)
    ppp = np.sqrt(xpa2 + ypb2 + zpc2)
    ppm = np.sqrt(xpa2 + ypb2 + zmc2)
    mpp = np.sqrt(xma2 + ypb2 + zpc2)

    # Avoid division by zero
    eps = 1e-15

    # Log terms for field computation
    with np.errstate(divide='ignore', invalid='ignore'):
        # Term for Bx from Jx, By from Jy contributions
        ff2x = np.log((xma + r000 + eps) * (xpa + ppm + eps) *
                      (xpa + pmp + eps) * (xma + mpp + eps)) - \
               np.log((xpa + pmm + eps) * (xma + mpm + eps) *
                      (xma + mmp + eps) * (xpa + ppp + eps))

        ff2y = np.log((-ymb + r000 + eps) * (-ypb + ppm + eps) *
                      (-ymb + pmp + eps) * (-ypb + mpp + eps)) - \
               np.log((-ymb + pmm + eps) * (-ypb + mpm + eps) *
                      (ymb - mmp + eps) * (ypb - ppp + eps))

        ff2z = np.log((-zmc + r000 + eps) * (-zmc + ppm + eps) *
                      (-zpc + pmp + eps) * (-zpc + mpp + eps)) - \
               np.log((-zmc + pmm + eps) * (zmc - mpm + eps) *
                      (-zpc + mmp + eps) * (zpc - ppp + eps))

    # Handle NaN from log
    if np.isnan(ff2x):
        ff2x = 0.0
    if np.isnan(ff2y):
        ff2y = 0.0
    if np.isnan(ff2z):
        ff2z = 0.0

    # Arctan terms for diagonal field contributions
    ff1x = (np.arctan2(ymb * zmc, xma * r000 + eps) -
            np.arctan2(ymb * zmc, xpa * pmm + eps) -
            np.arctan2(ypb * zmc, xma * mpm + eps) +
            np.arctan2(ypb * zmc, xpa * ppm + eps) -
            np.arctan2(ymb * zpc, xma * mmp + eps) +
            np.arctan2(ymb * zpc, xpa * pmp + eps) +
            np.arctan2(ypb * zpc, xma * mpp + eps) -
            np.arctan2(ypb * zpc, xpa * ppp + eps))

    ff1y = (np.arctan2(xma * zmc, ymb * r000 + eps) -
            np.arctan2(xpa * zmc, ymb * pmm + eps) -
            np.arctan2(xma * zmc, ypb * mpm + eps) +
            np.arctan2(xpa * zmc, ypb * ppm + eps) -
            np.arctan2(xma * zpc, ymb * mmp + eps) +
            np.arctan2(xpa * zpc, ymb * pmp + eps) +
            np.arctan2(xma * zpc, ypb * mpp + eps) -
            np.arctan2(xpa * zpc, ypb * ppp + eps))

    ff1z = (np.arctan2(xma * ymb, zmc * r000 + eps) -
            np.arctan2(xpa * ymb, zmc * pmm + eps) -
            np.arctan2(xma * ypb, zmc * mpm + eps) +
            np.arctan2(xpa * ypb, zmc * ppm + eps) -
            np.arctan2(xma * ymb, zpc * mmp + eps) +
            np.arctan2(xpa * ymb, zpc * pmp + eps) +
            np.arctan2(xma * ypb, zpc * mpp + eps) -
            np.arctan2(xpa * ypb, zpc * ppp + eps))

    # Combine contributions from each polarization component
    # Bx contributions
    bx_from_Jx = Jx * ff1x * qs[0, 0]
    bx_from_Jy = Jy * ff2z * qs[1, 0]
    bx_from_Jz = Jz * ff2y * qs[2, 0]

    # By contributions
    by_from_Jx = Jx * ff2z * qs[0, 1]
    by_from_Jy = Jy * ff1y * qs[1, 1]
    by_from_Jz = -Jz * ff2x * qs[2, 1]

    # Bz contributions
    bz_from_Jx = Jx * ff2y * qs[0, 2]
    bz_from_Jy = -Jy * ff2x * qs[1, 2]
    bz_from_Jz = Jz * ff1z * qs[2, 2]

    Bx = (bx_from_Jx + bx_from_Jy + bx_from_Jz) / (4.0 * PI)
    By = (by_from_Jx + by_from_Jy + by_from_Jz) / (4.0 * PI)
    Bz = (bz_from_Jx + bz_from_Jy + bz_from_Jz) / (4.0 * PI)

    return Bx, By, Bz


class CuboidMagnet:
    """
    Analytical cuboid (rectangular block) permanent magnet field source.

    Uses the Yang/Engel-Herbert/Camacho formulation based on magnetic surface charges.

    Attributes:
        center: [x, y, z] center position
        dimensions: [a, b, c] side lengths
        magnetization: [Mx, My, Mz] magnetization in A/m
    """

    def __init__(self, center: List[float], dimensions: List[float],
                 magnetization: List[float], units: str = 'mm'):
        """
        Initialize cuboid magnet.

        Parameters:
            center: [x, y, z] center position
            dimensions: [a, b, c] side lengths (NOT half-dimensions)
            magnetization: [Mx, My, Mz] in A/m
            units: length unit for geometry ('mm' or 'm')
        """
        self._length_to_m = 0.001 if units == 'mm' else 1.0
        self.center = np.array(center, dtype=float)
        self.dimensions = np.array(dimensions, dtype=float)
        self.half_dim = self.dimensions / 2.0
        self.magnetization = np.array(magnetization, dtype=float)
        # Polarization J = mu0 * M [T]
        self.polarization = MU0 * self.magnetization

    def _is_inside(self, point: np.ndarray) -> bool:
        """Check if point is inside the cuboid."""
        rel = np.abs(point - self.center)
        return np.all(rel <= self.half_dim * (1.0 + 1e-10))

    def _is_on_edge(self, point: np.ndarray) -> bool:
        """Check if point is on an edge (numerical instability region)."""
        rel = np.abs(point - self.center)
        tol = 1e-12 * self.half_dim
        on_surface = np.abs(rel - self.half_dim) < tol
        inside = rel < self.half_dim + tol
        # On edge if on 2+ surfaces and inside along remaining dimension
        return np.sum(on_surface) >= 2 and np.all(inside)

    def get_B(self, point: List[float]) -> List[float]:
        """
        Get magnetic field B at observation point.

        Parameters:
            point: [x, y, z] observation point (in constructor length units)

        Returns:
            [Bx, By, Bz] magnetic field in Tesla
        """
        p = np.array(point, dtype=float)

        # Check for edge singularity
        if self._is_on_edge(p):
            return [0.0, 0.0, 0.0]

        # Relative position from center
        rel = p - self.center

        Bx, By, Bz = _cuboid_bfield_component(
            rel[0], rel[1], rel[2],
            self.half_dim[0], self.half_dim[1], self.half_dim[2],
            self.polarization[0], self.polarization[1], self.polarization[2]
        )

        return [Bx, By, Bz]

    def get_H(self, point: List[float]) -> List[float]:
        """
        Get H-field at observation point.

        Parameters:
            point: [x, y, z] observation point (in constructor length units)

        Returns:
            [Hx, Hy, Hz] H-field in A/m
        """
        B = np.array(self.get_B(point))
        p = np.array(point, dtype=float)

        if self._is_inside(p):
            # Inside: H = B/mu0 - M
            return list(B / MU0 - self.magnetization)
        else:
            # Outside: H = B/mu0
            return list(B / MU0)

    def get_A(self, point: List[float]) -> List[float]:
        """
        Get vector potential A at observation point.

        Uses the analytical closed-form solution based on the equivalent surface
        current model. A uniformly magnetized cuboid is equivalent to surface
        currents K = M x n flowing on its faces.

        For z-magnetization, currents flow around the four side faces, equivalent
        to two rectangular current loops at z = +c and z = -c.

        The vector potential of a rectangular current loop has closed-form
        expressions involving logarithms (NIST formula).

        Parameters:
            point: [x, y, z] observation point (in constructor length units)

        Returns:
            [Ax, Ay, Az] vector potential in T*m (SI units)

        Reference:
            NIST J. Res. 105(4), 2000 - Calculation of Mutual Inductance
        """
        return self._get_A_analytical(point)

    def _get_A_analytical(self, point: List[float]) -> List[float]:
        """
        Analytical vector potential using equivalent surface current model.

        For a uniformly magnetized cuboid, the equivalent surface current density
        is K = M x n on each face, where n is the outward normal.

        For z-magnetization M = (0, 0, Mz):
        - Right face (x = +a, n = +x): K = Mz * (+z x +x) = Mz * (-y) -> K flows in -y
        - Left face (x = -a, n = -x):  K = Mz * (+z x -x) = Mz * (+y) -> K flows in +y
        - Front face (y = +b, n = +y): K = Mz * (+z x +y) = Mz * (-x) -> K flows in -x
        - Back face (y = -b, n = -y):  K = Mz * (+z x -y) = Mz * (+x) -> K flows in +x
        - Top/bottom faces: K = 0 (M parallel to n)

        The vector potential is computed by integrating:
        A(r) = (mu0 / 4*pi) * integral_S [K / |r - r'|] dS'

        This uses the analytical formula for integrating 1/|r - r'| over a rectangle.
        """
        # Convert to SI units (meters)
        p = (np.array(point, dtype=float) - self.center) * self._length_to_m  # [m]
        a = self.half_dim[0] * self._length_to_m  # half-width in x [m]
        b = self.half_dim[1] * self._length_to_m  # half-width in y [m]
        c = self.half_dim[2] * self._length_to_m  # half-width in z [m]
        M = self.magnetization  # [A/m]

        Ax, Ay, Az = 0.0, 0.0, 0.0
        pf = MU0 / (4.0 * PI)

        # Mz contribution: surface currents on 4 side faces
        # K = M x n = (0,0,Mz) x n
        if abs(M[2]) > 1e-15:
            Mz = M[2]
            # Right face (x = +a, n = +x): K = (0,0,Mz) x (1,0,0) = (0, Mz, 0)
            # Integrate over y' in [-b, b], z' in [-c, c]
            Ay += pf * Mz * self._rect_integral(p[0] - a, p[1], p[2], b, c)

            # Left face (x = -a, n = -x): K = (0,0,Mz) x (-1,0,0) = (0, -Mz, 0)
            Ay += pf * (-Mz) * self._rect_integral(p[0] + a, p[1], p[2], b, c)

            # Front face (y = +b, n = +y): K = (0,0,Mz) x (0,1,0) = (-Mz, 0, 0)
            Ax += pf * (-Mz) * self._rect_integral(p[1] - b, p[0], p[2], a, c)

            # Back face (y = -b, n = -y): K = (0,0,Mz) x (0,-1,0) = (Mz, 0, 0)
            Ax += pf * Mz * self._rect_integral(p[1] + b, p[0], p[2], a, c)

        # Mx contribution: surface currents on 4 side faces
        # K = M x n = (Mx,0,0) x n
        # x-faces have K=0 (M parallel to n)
        if abs(M[0]) > 1e-15:
            Mx = M[0]
            # Front face (y = +b, n = +y): K = (Mx,0,0) x (0,1,0) = (0, 0, Mx)
            Az += pf * Mx * self._rect_integral(p[1] - b, p[2], p[0], c, a)

            # Back face (y = -b, n = -y): K = (Mx,0,0) x (0,-1,0) = (0, 0, -Mx)
            Az += pf * (-Mx) * self._rect_integral(p[1] + b, p[2], p[0], c, a)

            # Top face (z = +c, n = +z): K = (Mx,0,0) x (0,0,1) = (0, -Mx, 0)
            Ay += pf * (-Mx) * self._rect_integral(p[2] - c, p[1], p[0], b, a)

            # Bottom face (z = -c, n = -z): K = (Mx,0,0) x (0,0,-1) = (0, Mx, 0)
            Ay += pf * Mx * self._rect_integral(p[2] + c, p[1], p[0], b, a)

        # My contribution: surface currents on 4 side faces
        # K = M x n = (0,My,0) x n
        # y-faces have K=0 (M parallel to n)
        if abs(M[1]) > 1e-15:
            My = M[1]
            # Right face (x = +a, n = +x): K = (0,My,0) x (1,0,0) = (0, 0, -My)
            Az += pf * (-My) * self._rect_integral(p[0] - a, p[2], p[1], c, b)

            # Left face (x = -a, n = -x): K = (0,My,0) x (-1,0,0) = (0, 0, My)
            Az += pf * My * self._rect_integral(p[0] + a, p[2], p[1], c, b)

            # Top face (z = +c, n = +z): K = (0,My,0) x (0,0,1) = (My, 0, 0)
            Ax += pf * My * self._rect_integral(p[2] - c, p[0], p[1], a, b)

            # Bottom face (z = -c, n = -z): K = (0,My,0) x (0,0,-1) = (-My, 0, 0)
            Ax += pf * (-My) * self._rect_integral(p[2] + c, p[0], p[1], a, b)

        return [Ax, Ay, Az]

    def _rect_integral(self, d: float, u: float, v: float,
                       hu: float, hv: float) -> float:
        """
        Analytical integral of 1/|r - r'| over a rectangle.

        The rectangle lies in the plane at distance d from the observation point,
        with the observation point projected onto the rectangle at (u, v).
        Rectangle spans u' in [-hu, hu], v' in [-hv, hv].

        Uses the formula from Urankar (1980) and Ravaud (2009):
        integral_S [1/r] dS = sum of F(corner terms)

        where F involves arcsinh and arctan functions.
        """
        result = 0.0
        signs = [(1, 1), (-1, 1), (-1, -1), (1, -1)]

        for su, sv in signs:
            uu = su * hu - u  # u' - u
            vv = sv * hv - v  # v' - v
            rr = np.sqrt(d*d + uu*uu + vv*vv)

            # Avoid division by zero
            eps = 1e-30

            # F = uu * arcsinh(vv / sqrt(d^2 + uu^2)) + vv * arcsinh(uu / sqrt(d^2 + vv^2))
            #     - d * arctan(uu * vv / (d * rr))
            term1 = uu * np.arcsinh(vv / (np.sqrt(d*d + uu*uu) + eps)) if abs(d*d + uu*uu) > eps*eps else 0.0
            term2 = vv * np.arcsinh(uu / (np.sqrt(d*d + vv*vv) + eps)) if abs(d*d + vv*vv) > eps*eps else 0.0
            term3 = -d * np.arctan2(uu * vv, d * rr + eps) if abs(d) > eps else 0.0

            result += su * sv * (term1 + term2 + term3)

        return result

    def _rect_loop_A_xy(self, x: float, y: float, z: float,
                        a: float, b: float, I: float) -> Tuple[float, float]:
        """
        Vector potential (Ax, Ay) of a rectangular current loop in the xy-plane.

        Loop has corners at (±a, ±b, 0) with current I flowing counterclockwise
        when viewed from +z. The observation point is at (x, y, z).

        Based on NIST formula:
        Ax = (mu0*I/4*pi) * ln[(r1+a+x)(r2-a+x) / (r3-a+x)(r4+a+x)]
        Ay = (mu0*I/4*pi) * ln[(r2+b+y)(r3-b+y) / (r4-b+y)(r1+b+y)]

        where r1, r2, r3, r4 are distances from corners to observation point.
        """
        # Small epsilon to avoid log(0)
        eps = 1e-30

        # Distances from 4 corners to observation point
        # Corner 1: (+a, +b, 0)
        r1 = np.sqrt((x - a)**2 + (y - b)**2 + z**2)
        # Corner 2: (-a, +b, 0)
        r2 = np.sqrt((x + a)**2 + (y - b)**2 + z**2)
        # Corner 3: (-a, -b, 0)
        r3 = np.sqrt((x + a)**2 + (y + b)**2 + z**2)
        # Corner 4: (+a, -b, 0)
        r4 = np.sqrt((x - a)**2 + (y + b)**2 + z**2)

        # Prefactor
        pf = MU0 * I / (4.0 * PI)

        # Ax contribution (from y-directed segments at y = +b and y = -b)
        # Segment at y = +b: from (-a, +b) to (+a, +b) -> contributes to Ax
        # Segment at y = -b: from (+a, -b) to (-a, -b) -> contributes to -Ax
        num_x = (r1 + a - x + eps) * (r2 + a + x + eps)
        den_x = (r3 + a + x + eps) * (r4 + a - x + eps)
        if num_x > 0 and den_x > 0:
            Ax = pf * np.log(num_x / den_x)
        else:
            Ax = 0.0

        # Ay contribution (from x-directed segments at x = +a and x = -a)
        num_y = (r2 + b - y + eps) * (r3 + b + y + eps)
        den_y = (r4 + b + y + eps) * (r1 + b - y + eps)
        if num_y > 0 and den_y > 0:
            Ay = pf * np.log(num_y / den_y)
        else:
            Ay = 0.0

        return Ax, Ay

    def _rect_loop_A_yz(self, y: float, z: float, x: float,
                        b: float, c: float, I: float) -> Tuple[float, float]:
        """
        Vector potential (Ay, Az) of a rectangular loop in the yz-plane at x=const.
        Loop corners at (0, ±b, ±c).
        """
        eps = 1e-30

        r1 = np.sqrt(x**2 + (y - b)**2 + (z - c)**2)
        r2 = np.sqrt(x**2 + (y + b)**2 + (z - c)**2)
        r3 = np.sqrt(x**2 + (y + b)**2 + (z + c)**2)
        r4 = np.sqrt(x**2 + (y - b)**2 + (z + c)**2)

        pf = MU0 * I / (4.0 * PI)

        num_y = (r1 + b - y + eps) * (r2 + b + y + eps)
        den_y = (r3 + b + y + eps) * (r4 + b - y + eps)
        if num_y > 0 and den_y > 0:
            Ay = pf * np.log(num_y / den_y)
        else:
            Ay = 0.0

        num_z = (r2 + c - z + eps) * (r3 + c + z + eps)
        den_z = (r4 + c + z + eps) * (r1 + c - z + eps)
        if num_z > 0 and den_z > 0:
            Az = pf * np.log(num_z / den_z)
        else:
            Az = 0.0

        return Ay, Az

    def _rect_loop_A_xz(self, z: float, x: float, y: float,
                        c: float, a: float, I: float) -> Tuple[float, float]:
        """
        Vector potential (Az, Ax) of a rectangular loop in the xz-plane at y=const.
        Loop corners at (±a, 0, ±c).
        """
        eps = 1e-30

        r1 = np.sqrt((x - a)**2 + y**2 + (z - c)**2)
        r2 = np.sqrt((x + a)**2 + y**2 + (z - c)**2)
        r3 = np.sqrt((x + a)**2 + y**2 + (z + c)**2)
        r4 = np.sqrt((x - a)**2 + y**2 + (z + c)**2)

        pf = MU0 * I / (4.0 * PI)

        num_z = (r1 + c - z + eps) * (r2 + c - z + eps)
        den_z = (r3 + c + z + eps) * (r4 + c + z + eps)
        if num_z > 0 and den_z > 0:
            Az = pf * np.log(num_z / den_z)
        else:
            Az = 0.0

        num_x = (r2 + a + x + eps) * (r3 + a + x + eps)
        den_x = (r4 + a - x + eps) * (r1 + a - x + eps)
        if num_x > 0 and den_x > 0:
            Ax = pf * np.log(num_x / den_x)
        else:
            Ax = 0.0

        return Az, Ax

    def __call__(self, point: List[float]) -> List[float]:
        """Callable interface for use with rad.ObjBckg()."""
        return self.get_B(point)


# =============================================================================
# Current Loop
# =============================================================================

def _cel_iter(qc: float, p: float, g: float, cc: float, ss: float,
               em: float, kk: float) -> float:
    """
    Iterative part of Bulirsch cel algorithm.

    Reference: Implementation from Ortner et al. (2022), Magnetism 2023, 3(1), 11-31.
    Based on magpylib implementation.
    """
    import math
    while math.fabs(g - qc) >= qc * 1e-8:
        qc = 2.0 * math.sqrt(kk)
        kk = qc * em
        f = cc
        cc = cc + ss / p
        g = kk / p
        ss = 2.0 * (ss + f * g)
        p = p + g
        g = em
        em = em + qc
    return 1.5707963267948966 * (ss + cc * em) / (em * (em + p))


class CurrentLoop:
    """
    Analytical circular current loop field source.

    Uses the Ortner et al. (2022) numerically stable formulation based on
    generalized complete elliptic integrals.

    Attributes:
        center: [x, y, z] center position
        diameter: loop diameter
        current: current in Amperes
        axis: loop axis direction ('x', 'y', or 'z')
    """

    def __init__(self, center: List[float], diameter: float, current: float,
                 axis: str = 'z', units: str = 'mm'):
        """
        Initialize current loop.

        Parameters:
            center: [x, y, z] center position
            diameter: loop diameter
            current: current in Amperes (positive = CCW when viewed from +axis)
            axis: loop axis direction ('x', 'y', or 'z')
            units: length unit for geometry ('mm' or 'm')
        """
        self._length_to_m = 0.001 if units == 'mm' else 1.0
        self.center = np.array(center, dtype=float)
        self.diameter = float(diameter)
        self.radius = diameter / 2.0
        self.current = float(current)
        self.axis = axis.lower()

        if self.axis not in ['x', 'y', 'z']:
            raise ValueError("axis must be 'x', 'y', or 'z'")

    def _field_cylindrical(self, rho: float, z: float) -> Tuple[float, float]:
        """
        H-field in cylindrical coordinates (Hr, Hz) in A/m.

        This uses the Ortner et al. (2022) formulation.
        Input: rho, z in mm (same unit as self.radius)
        Output: Hr, Hz in A/m

        Parameters:
            rho: radial distance from axis [mm]
            z: axial distance from loop plane [mm]

        Returns:
            (Hr, Hz) in A/m
        """
        # Convert to meters for SI calculations
        r0 = self.radius * self._length_to_m  # loop radius [m]
        rho_m = rho * self._length_to_m  # [m]
        z_m = z * self._length_to_m  # [m]
        I = self.current  # current [A]

        if r0 == 0:
            return 0.0, 0.0

        # On-axis case: H_z = I * r0^2 / (2 * (z^2 + r0^2)^(3/2))
        # With SI units (meters): H has units A/m
        if rho_m < 1e-15 * r0:
            Hz = I * r0 ** 2 / (2.0 * (z_m ** 2 + r0 ** 2) ** 1.5)
            return 0.0, Hz

        # Singularity check (on the loop)
        if abs(rho_m - r0) < 1e-15 * r0 and abs(z_m) < 1e-15 * r0:
            return 0.0, 0.0

        # Dimensionless ratios (for Ortner formulation)
        r = rho_m / r0
        zn = z_m / r0

        # Field computation from Ortner et al. (2022)
        z2 = zn ** 2
        x0 = z2 + (r + 1.0) ** 2
        k2 = 4.0 * r / x0
        q2 = (z2 + (r - 1.0) ** 2) / x0

        k = np.sqrt(k2)
        q = np.sqrt(q2)
        p = 1.0 + q
        # pf factor with 1/r0 in meters for correct A/m output
        pf = k / np.sqrt(r) / q2 / 20.0 / r0 * 1e-6 * I

        # cel* for Hr (using _cel_iter)
        cc = k2 * k2
        ss = 2.0 * cc * q / p
        Hr = pf * zn / r * _cel_iter(q, p, 1.0, cc, ss, p, q)

        # cel** for Hz
        cc = k2 * (k2 - (q2 + 1.0) / r)
        ss = 2.0 * k2 * q * (k2 / p - p / r)
        Hz = -pf * _cel_iter(q, p, 1.0, cc, ss, p, q)

        # Scale to A/m (factor from Ortner paper: 1e7 / (4*pi))
        scale = 795774.7154594767
        return Hr * scale, Hz * scale

    def get_B(self, point: List[float]) -> List[float]:
        """
        Get magnetic field B at observation point.

        Parameters:
            point: [x, y, z] observation point (in constructor length units)

        Returns:
            [Bx, By, Bz] magnetic field in Tesla
        """
        # Transform to local coordinates
        p = np.array(point, dtype=float) - self.center

        # Rotate to align loop axis with z-axis
        if self.axis == 'z':
            p_local = p.copy()
        elif self.axis == 'x':
            p_local = np.array([p[1], p[2], p[0]])
        else:  # axis == 'y'
            p_local = np.array([p[2], p[0], p[1]])

        # Cylindrical coordinates
        rho = np.sqrt(p_local[0] ** 2 + p_local[1] ** 2)
        phi = np.arctan2(p_local[1], p_local[0])
        z = p_local[2]

        # Get H-field in cylindrical coordinates
        Hr, Hz = self._field_cylindrical(rho, z)

        # Convert H to B: B = mu0 * H
        Br = MU0 * Hr
        Bz = MU0 * Hz

        # Convert to Cartesian (Bphi = 0 by symmetry)
        Bx_local = Br * np.cos(phi)
        By_local = Br * np.sin(phi)
        Bz_local = Bz

        # Rotate back to global coordinates
        if self.axis == 'z':
            Bx, By, Bz = Bx_local, By_local, Bz_local
        elif self.axis == 'x':
            Bx, By, Bz = Bz_local, Bx_local, By_local
        else:  # axis == 'y'
            Bx, By, Bz = By_local, Bz_local, Bx_local

        return [Bx, By, Bz]

    def get_H(self, point: List[float]) -> List[float]:
        """
        Get H-field at observation point.

        Parameters:
            point: [x, y, z] observation point (in constructor length units)

        Returns:
            [Hx, Hy, Hz] H-field in A/m
        """
        B = self.get_B(point)
        return [B[0] / MU0, B[1] / MU0, B[2] / MU0]

    def _vector_potential_cylindrical(self, rho: float, z: float) -> float:
        """
        Vector potential A_phi of a circular current loop.

        Uses the standard elliptic integral formula:
        A_phi = (mu0*I / pi) * sqrt(a/rho) * [(1 - k^2/2)*K(k) - E(k)] / k

        where k^2 = 4*a*rho / ((a+rho)^2 + z^2)

        Parameters:
            rho: radial distance from axis [mm]
            z: axial distance from loop plane [mm]

        Returns:
            A_phi in T*m (compatible with curl(A) = B in Tesla)
        """
        from scipy.special import ellipe, ellipk

        # Convert to meters for SI calculation
        a = self.radius * self._length_to_m  # loop radius [m]
        rho_m = rho * self._length_to_m  # [m]
        z_m = z * self._length_to_m  # [m]
        I = self.current  # [A]

        # On-axis singularity: A_phi = 0 by symmetry
        if rho_m < 1e-15:
            return 0.0

        if a <= 0.0:
            return 0.0

        # Distance factors
        apr = a + rho_m
        apr2 = apr * apr
        z2 = z_m * z_m

        denom = apr2 + z2
        if denom < 1e-30:
            return 0.0

        # Elliptic integral argument k^2
        k_sq = 4.0 * a * rho_m / denom
        k = np.sqrt(k_sq)

        # Avoid singularity when k -> 1 (on the loop)
        if k >= 0.9999:
            k = 0.9999
            k_sq = k * k

        # Elliptic integrals K(m) and E(m) where m = k^2
        K_val = ellipk(k_sq)
        E_val = ellipe(k_sq)

        # Standard formula: A_phi = (mu0*I / pi) * sqrt(a/rho) * [(1 - k^2/2)*K - E] / k
        factor = (MU0 * I / PI) * np.sqrt(a / rho_m)
        bracket = ((1.0 - k_sq / 2.0) * K_val - E_val) / k

        A_phi = factor * bracket

        return A_phi

    def get_A(self, point: List[float]) -> List[float]:
        """
        Get vector potential A at observation point.

        The vector potential of a current loop has only an azimuthal component A_phi.

        Parameters:
            point: [x, y, z] observation point (in constructor length units)

        Returns:
            [Ax, Ay, Az] vector potential in T*m (SI units for curl(A) = B)
        """
        # Transform to local coordinates
        p = np.array(point, dtype=float) - self.center

        # Rotate to align loop axis with z-axis
        if self.axis == 'z':
            p_local = p.copy()
        elif self.axis == 'x':
            p_local = np.array([p[1], p[2], p[0]])
        else:  # axis == 'y'
            p_local = np.array([p[2], p[0], p[1]])

        # Cylindrical coordinates
        rho = np.sqrt(p_local[0] ** 2 + p_local[1] ** 2)
        phi = np.arctan2(p_local[1], p_local[0])
        z = p_local[2]

        # Get A_phi
        A_phi = self._vector_potential_cylindrical(rho, z)

        # Convert to Cartesian: Ax = -A_phi * sin(phi), Ay = A_phi * cos(phi)
        Ax_local = -A_phi * np.sin(phi)
        Ay_local = A_phi * np.cos(phi)
        Az_local = 0.0

        # Rotate back to global coordinates
        if self.axis == 'z':
            Ax, Ay, Az = Ax_local, Ay_local, Az_local
        elif self.axis == 'x':
            Ax, Ay, Az = Az_local, Ax_local, Ay_local
        else:  # axis == 'y'
            Ax, Ay, Az = Ay_local, Az_local, Ax_local

        return [Ax, Ay, Az]

    def __call__(self, point: List[float]) -> List[float]:
        """Callable interface for use with rad.ObjBckg()."""
        return self.get_B(point)


# =============================================================================
# Re-export cylindrical magnet classes for convenience
# =============================================================================

try:
    from radia.cylindrical_magnet import CylindricalMagnet, RingMagnet  # noqa: F401
except ImportError:
    from cylindrical_magnet import CylindricalMagnet, RingMagnet  # noqa: F401

__all__ = [
    'SphericalMagnet',
    'CuboidMagnet',
    'CurrentLoop',
    'CylindricalMagnet',
    'RingMagnet',
]
