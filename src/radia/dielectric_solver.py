"""
Dielectric Solver using Surface Charge Method

This module implements the surface charge method for dielectric analysis,
designed to integrate with the Loop-Star PEEC framework.

The surface charge density σ_e on dielectric boundaries satisfies:
    ε_1 * E_1n - ε_2 * E_2n = σ_e

The potential is computed using the P-matrix (potential coefficient matrix):
    φ_i = Σ_j P_ij * σ_j

where P_ij = (1/4πε_0) * ∫∫ 1/|r_i - r'_j| dS'

For integration with Loop-Star:
    - Dielectric charges become Star DOFs
    - P_SS includes both conductor and dielectric contributions
    - Coupling terms P_cd link conductor and dielectric surfaces

References:
    [1] Rao, Wilton, Glisson, "Electromagnetic Scattering by Surfaces
        of Arbitrary Shape", IEEE Trans. AP, 1982.
    [2] Harrington, "Field Computation by Moment Methods", 1968.

Author: Radia Development Team
License: BSD-style AND MIT
"""

import numpy as np
from typing import Optional, Tuple, List, Dict
from dataclasses import dataclass

# Physical constants
EPS_0 = 8.854187817e-12  # F/m
INV_FOUR_PI_EPS0 = 1.0 / (4.0 * np.pi * EPS_0)


@dataclass
class DielectricMesh:
    """
    Surface mesh for dielectric boundary.

    Attributes:
        vertices: Vertex coordinates [n_vertices, 3]
        triangles: Triangle connectivity [n_triangles, 3] (0-indexed)
        eps_r_inside: Relative permittivity inside the surface
        eps_r_outside: Relative permittivity outside the surface
        centroids: Triangle centroids [n_triangles, 3]
        areas: Triangle areas [n_triangles]
        normals: Outward unit normals [n_triangles, 3]
    """
    vertices: np.ndarray
    triangles: np.ndarray
    eps_r_inside: float = 1.0
    eps_r_outside: float = 1.0
    centroids: Optional[np.ndarray] = None
    areas: Optional[np.ndarray] = None
    normals: Optional[np.ndarray] = None

    def __post_init__(self):
        """Compute centroids, areas, and normals if not provided."""
        if self.centroids is None or self.areas is None or self.normals is None:
            self._compute_geometry()

    def _compute_geometry(self):
        """Compute triangle centroids, areas, and normals."""
        n_tri = len(self.triangles)
        self.centroids = np.zeros((n_tri, 3))
        self.areas = np.zeros(n_tri)
        self.normals = np.zeros((n_tri, 3))

        for i, tri in enumerate(self.triangles):
            v0 = self.vertices[tri[0]]
            v1 = self.vertices[tri[1]]
            v2 = self.vertices[tri[2]]

            # Centroid
            self.centroids[i] = (v0 + v1 + v2) / 3.0

            # Area and normal via cross product
            edge1 = v1 - v0
            edge2 = v2 - v0
            cross = np.cross(edge1, edge2)
            area = 0.5 * np.linalg.norm(cross)
            self.areas[i] = area

            if area > 1e-30:
                self.normals[i] = cross / (2.0 * area)
            else:
                self.normals[i] = np.array([0.0, 0.0, 1.0])


class DielectricSolver:
    """
    Surface charge method solver for dielectric bodies.

    Computes the P-matrix (potential coefficient matrix) for dielectric
    surfaces and solves for surface charge distribution given boundary
    conditions.

    The solver can handle:
    - Single dielectric body in free space
    - Multiple dielectric bodies with different permittivities
    - Coupling with conductor surfaces (for PEEC integration)

    Example:
        >>> # Create mesh for dielectric sphere
        >>> mesh = DielectricMesh(vertices, triangles, eps_r_inside=4.0)
        >>>
        >>> # Create solver
        >>> solver = DielectricSolver()
        >>> solver.add_dielectric(mesh, name='sphere')
        >>>
        >>> # Build P-matrix
        >>> solver.build()
        >>>
        >>> # Get P-matrix for Loop-Star integration
        >>> P_dd = solver.get_P_matrix()
    """

    def __init__(self):
        self._dielectrics: Dict[str, DielectricMesh] = {}
        self._P = None  # Full P-matrix
        self._P_blocks = {}  # Block P-matrices for multi-body
        self._is_built = False
        self._n_total = 0  # Total number of DOFs

    def add_dielectric(self, mesh: DielectricMesh, name: str = 'default'):
        """
        Add a dielectric body to the solver.

        Parameters:
            mesh: DielectricMesh object with surface mesh and permittivity
            name: Unique name for this dielectric body
        """
        self._dielectrics[name] = mesh
        self._is_built = False

    def build(self, use_gauss_quadrature: bool = True, n_quad: int = 4):
        """
        Build the P-matrix (potential coefficient matrix).

        The P-matrix relates surface charge to potential:
            φ = P * σ

        For a dielectric interface with ε_1 inside and ε_2 outside:
            P_ij = (1/4πε_eff) * ∫∫ 1/|r_i - r'| dS'

        where ε_eff = (ε_1 + ε_2) / 2 for the interface.

        Parameters:
            use_gauss_quadrature: Use Gaussian quadrature for integration
            n_quad: Number of quadrature points per triangle
        """
        # Collect all triangles
        all_centroids = []
        all_areas = []
        all_vertices = []
        all_triangles = []
        all_eps_eff = []
        vertex_offset = 0

        for name, mesh in self._dielectrics.items():
            n_tri = len(mesh.triangles)
            all_centroids.append(mesh.centroids)
            all_areas.append(mesh.areas)

            # Offset triangle indices
            triangles_offset = mesh.triangles + vertex_offset
            all_triangles.append(triangles_offset)
            all_vertices.append(mesh.vertices)
            vertex_offset += len(mesh.vertices)

            # Effective permittivity at interface
            eps_eff = 0.5 * (mesh.eps_r_inside + mesh.eps_r_outside) * EPS_0
            all_eps_eff.extend([eps_eff] * n_tri)

        centroids = np.vstack(all_centroids)
        areas = np.concatenate(all_areas)
        vertices = np.vstack(all_vertices)
        triangles = np.vstack(all_triangles)
        eps_eff = np.array(all_eps_eff)

        self._n_total = len(centroids)
        n = self._n_total

        # Build P-matrix
        self._P = np.zeros((n, n))

        if use_gauss_quadrature:
            # Gauss quadrature points for triangle (barycentric)
            if n_quad == 1:
                # 1-point rule (centroid)
                xi = np.array([[1/3, 1/3, 1/3]])
                wi = np.array([1.0])
            elif n_quad == 3:
                # 3-point rule
                xi = np.array([
                    [2/3, 1/6, 1/6],
                    [1/6, 2/3, 1/6],
                    [1/6, 1/6, 2/3]
                ])
                wi = np.array([1/3, 1/3, 1/3])
            elif n_quad == 4:
                # 4-point rule
                xi = np.array([
                    [1/3, 1/3, 1/3],
                    [3/5, 1/5, 1/5],
                    [1/5, 3/5, 1/5],
                    [1/5, 1/5, 3/5]
                ])
                wi = np.array([-27/48, 25/48, 25/48, 25/48])
            else:
                # Default to centroid
                xi = np.array([[1/3, 1/3, 1/3]])
                wi = np.array([1.0])

            for i in range(n):
                r_obs = centroids[i]

                for j in range(n):
                    tri_j = triangles[j]
                    v0 = vertices[tri_j[0]]
                    v1 = vertices[tri_j[1]]
                    v2 = vertices[tri_j[2]]
                    A_j = areas[j]
                    eps_j = eps_eff[j]

                    if i == j:
                        # Self-term: use analytical formula
                        # For flat triangle: P_ii ≈ (1/2πε) * sqrt(A/π)
                        # This is approximate; exact formula is more complex
                        r_eff = np.sqrt(A_j / np.pi)
                        self._P[i, j] = 1.0 / (2.0 * np.pi * eps_j) * np.log(
                            (1.0 + np.sqrt(2.0)) * r_eff / (r_eff * 0.1 + 1e-10)
                        )
                        # Simpler approximation
                        self._P[i, j] = 1.0 / (4.0 * np.pi * eps_j * r_eff) * 1.5
                    else:
                        # Quadrature integration
                        integral = 0.0
                        for k, (xik, wk) in enumerate(zip(xi, wi)):
                            r_src = xik[0] * v0 + xik[1] * v1 + xik[2] * v2
                            dist = np.linalg.norm(r_obs - r_src)
                            if dist > 1e-30:
                                integral += wk / dist

                        self._P[i, j] = A_j * integral / (4.0 * np.pi * eps_j)
        else:
            # Simple centroid-to-centroid
            for i in range(n):
                for j in range(n):
                    eps_j = eps_eff[j]
                    if i == j:
                        r_eff = np.sqrt(areas[j] / np.pi)
                        self._P[i, j] = 1.0 / (4.0 * np.pi * eps_j * r_eff) * 1.5
                    else:
                        dist = np.linalg.norm(centroids[i] - centroids[j])
                        if dist > 1e-30:
                            self._P[i, j] = areas[j] / (4.0 * np.pi * eps_j * dist)

        self._is_built = True

    def get_P_matrix(self) -> np.ndarray:
        """
        Get the P-matrix for Loop-Star integration.

        The P-matrix can be used directly as the Star-Star block
        in the Loop-Star impedance system:
            Z_SS = P / s = P / (jω)

        Returns:
            P-matrix [n_dof, n_dof] in units of [1/F]
        """
        if not self._is_built:
            raise RuntimeError("Solver not built. Call build() first.")
        return self._P

    def get_capacitance_matrix(self) -> np.ndarray:
        """
        Get the capacitance matrix C = inv(P).

        The capacitance matrix relates charge to voltage:
            Q = C * V

        Returns:
            Capacitance matrix [n_dof, n_dof] in units of [F]
        """
        if not self._is_built:
            raise RuntimeError("Solver not built. Call build() first.")
        return np.linalg.inv(self._P)

    def solve_static(self, V_boundary: np.ndarray) -> Tuple[np.ndarray, float]:
        """
        Solve for surface charge given boundary potentials.

        Parameters:
            V_boundary: Potential at each triangle centroid [V]

        Returns:
            sigma: Surface charge density [C/m^2]
            total_charge: Total charge on surface [C]
        """
        if not self._is_built:
            raise RuntimeError("Solver not built. Call build() first.")

        # σ = P^{-1} * φ
        sigma = np.linalg.solve(self._P, V_boundary)

        # Total charge
        areas = []
        for mesh in self._dielectrics.values():
            areas.extend(mesh.areas)
        areas = np.array(areas)
        total_charge = np.sum(sigma * areas)

        return sigma, total_charge

    def compute_coupling_matrix(self, conductor_centroids: np.ndarray,
                                 conductor_areas: np.ndarray) -> np.ndarray:
        """
        Compute coupling P-matrix between conductor and dielectric surfaces.

        This matrix represents the potential at dielectric points due to
        conductor surface charges, and vice versa.

        Parameters:
            conductor_centroids: Conductor triangle centroids [n_cond, 3]
            conductor_areas: Conductor triangle areas [n_cond]

        Returns:
            P_cd: Coupling matrix [n_dielectric, n_conductor]
        """
        if not self._is_built:
            raise RuntimeError("Solver not built. Call build() first.")

        # Collect dielectric centroids
        dielectric_centroids = []
        for mesh in self._dielectrics.values():
            dielectric_centroids.append(mesh.centroids)
        dielectric_centroids = np.vstack(dielectric_centroids)

        n_d = len(dielectric_centroids)
        n_c = len(conductor_centroids)

        P_cd = np.zeros((n_d, n_c))

        for i in range(n_d):
            r_obs = dielectric_centroids[i]
            for j in range(n_c):
                r_src = conductor_centroids[j]
                dist = np.linalg.norm(r_obs - r_src)
                if dist > 1e-30:
                    # Use free space permittivity for conductor coupling
                    P_cd[i, j] = conductor_areas[j] / (4.0 * np.pi * EPS_0 * dist)

        return P_cd


def create_sphere_mesh(center: np.ndarray, radius: float,
                       n_theta: int = 10, n_phi: int = 20) -> DielectricMesh:
    """
    Create a triangular mesh for a sphere surface.

    Parameters:
        center: Sphere center [x, y, z]
        radius: Sphere radius
        n_theta: Number of divisions in polar angle
        n_phi: Number of divisions in azimuthal angle

    Returns:
        DielectricMesh object
    """
    vertices = []

    # Generate vertices
    for i in range(n_theta + 1):
        theta = np.pi * i / n_theta
        for j in range(n_phi):
            phi = 2.0 * np.pi * j / n_phi
            x = center[0] + radius * np.sin(theta) * np.cos(phi)
            y = center[1] + radius * np.sin(theta) * np.sin(phi)
            z = center[2] + radius * np.cos(theta)
            vertices.append([x, y, z])

    vertices = np.array(vertices)

    # Generate triangles
    triangles = []
    for i in range(n_theta):
        for j in range(n_phi):
            # Current row indices
            i0 = i * n_phi + j
            i1 = i * n_phi + (j + 1) % n_phi
            # Next row indices
            i2 = (i + 1) * n_phi + j
            i3 = (i + 1) * n_phi + (j + 1) % n_phi

            # Two triangles per quad (skip degenerate at poles)
            if i > 0:
                triangles.append([i0, i2, i1])
            if i < n_theta - 1:
                triangles.append([i1, i2, i3])

    triangles = np.array(triangles)

    return DielectricMesh(vertices, triangles)


def create_box_mesh(center: np.ndarray, size: np.ndarray,
                    n_div: int = 5) -> DielectricMesh:
    """
    Create a triangular mesh for a box surface.

    Parameters:
        center: Box center [x, y, z]
        size: Box dimensions [dx, dy, dz]
        n_div: Number of divisions per edge

    Returns:
        DielectricMesh object
    """
    hx, hy, hz = size / 2.0
    cx, cy, cz = center

    vertices = []
    triangles = []

    # Generate 6 faces
    faces = [
        # (normal direction, two tangent axes, offset)
        ([1, 0, 0], [0, 1, 0], [0, 0, 1], hx),   # +x
        ([-1, 0, 0], [0, -1, 0], [0, 0, 1], hx),  # -x
        ([0, 1, 0], [-1, 0, 0], [0, 0, 1], hy),   # +y
        ([0, -1, 0], [1, 0, 0], [0, 0, 1], hy),   # -y
        ([0, 0, 1], [1, 0, 0], [0, 1, 0], hz),    # +z
        ([0, 0, -1], [-1, 0, 0], [0, 1, 0], hz),  # -z
    ]

    vertex_offset = 0

    for normal, t1, t2, offset in faces:
        normal = np.array(normal)
        t1 = np.array(t1)
        t2 = np.array(t2)

        # Half-sizes in tangent directions
        h1 = np.abs(t1 @ np.array([hx, hy, hz]))
        h2 = np.abs(t2 @ np.array([hx, hy, hz]))

        # Generate grid vertices on this face
        face_vertices = []
        for i in range(n_div + 1):
            for j in range(n_div + 1):
                u = -h1 + 2 * h1 * i / n_div
                v = -h2 + 2 * h2 * j / n_div
                p = center + normal * offset + t1 * u + t2 * v
                face_vertices.append(p)

        vertices.extend(face_vertices)

        # Generate triangles
        for i in range(n_div):
            for j in range(n_div):
                i0 = vertex_offset + i * (n_div + 1) + j
                i1 = i0 + 1
                i2 = i0 + (n_div + 1)
                i3 = i2 + 1

                triangles.append([i0, i2, i1])
                triangles.append([i1, i2, i3])

        vertex_offset += len(face_vertices)

    return DielectricMesh(np.array(vertices), np.array(triangles))


# =============================================================================
# Integration with Loop-Star: Conductor + Dielectric + Magnetic
# =============================================================================

class UnifiedSurfaceSolver:
    """
    Unified surface solver for conductors, dielectrics, and magnetic materials.

    Surface quantities:
    - Conductor: J_s (surface current) + σ_c (surface charge)
    - Dielectric: σ_d (surface charge)
    - Magnetic: σ_m (surface magnetic charge) or M (magnetization)

    System equation (Loop-Star formulation):

        [Z_LL    Z_LS_c   Z_LS_d   Z_LM  ] [I_L  ]   [V_L  ]
        [Z_SL_c  Z_SS_c   Z_cd     Z_cM  ] [Q_c  ] = [V_S_c]
        [Z_SL_d  Z_dc     Z_SS_d   Z_dM  ] [Q_d  ]   [V_S_d]
        [Z_ML    Z_Mc     Z_Md     Z_MM  ] [M    ]   [H_ext]

    where:
        Z_LL = R + sL (conductor inductance from SIBC/EFIE)
        Z_SS_c = P_cc / s (conductor-conductor potential)
        Z_SS_d = P_dd / s (dielectric-dielectric potential)
        Z_MM = N / (s*mu) (magnetic reluctance, mu can be complex)
        Z_LS = Loop-Star coupling (geometry dependent)
        Z_LM = current-magnetization coupling (Biot-Savart)

    Key features:
    1. Conductor current creates B-field that magnetizes magnetic material
    2. Magnetization creates H-field that affects conductor (back-EMF)
    3. Surface charges couple electrostatically
    4. All couplings are handled through unified matrix system

    Physical interpretation:
    - Loop currents flow on conductor surfaces (skin effect via SIBC)
    - Surface charges accumulate at conductor and dielectric boundaries
    - Magnetic materials are magnetized by external and induced fields
    """

    def __init__(self):
        # Conductor Loop (inductive)
        self._L_LL = None  # Inductance matrix
        self._R_LL = None  # DC Resistance matrix
        self._n_loops = 0

        # Conductor Star (capacitive)
        self._P_cc = None  # Conductor potential coefficient
        self._L_Lc = None  # Loop-Star_c coupling
        self._n_stars_c = 0

        # Dielectric Star (capacitive)
        self._P_dd = None  # Dielectric potential coefficient
        self._P_cd = None  # Conductor-dielectric coupling
        self._n_stars_d = 0

        # Magnetic material block
        self._N_mm = None  # Demagnetization matrix
        self._mu_r = 1.0   # Relative permeability (can be complex)
        self._L_Lm = None  # Loop-magnetic coupling (Biot-Savart)
        self._n_magnetic = 0

        # Cross-coupling matrices
        self._P_cm = None  # Conductor-magnetic electric coupling
        self._P_dm = None  # Dielectric-magnetic electric coupling

        # Skin effect parameters
        self._use_skin_effect = False
        self._skin_effect_model = 'dowell'  # 'dowell' or 'sibc'
        self._sigma = 5.8e7  # Conductivity [S/m] (copper default)
        self._conductor_height = None  # Conductor height for Dowell [m]
        self._conductor_width = None   # Conductor width for Dowell [m]
        self._mu_conductor = 4.0 * np.pi * 1e-7  # Permeability of conductor

    def set_conductor_loop(self, L: np.ndarray, R: np.ndarray):
        """
        Set conductor Loop matrices.

        The inductance L comes from:
        - Partial inductance (PEEC)
        - Or SIBC integral equation

        The resistance R includes:
        - DC resistance
        - Skin effect (SIBC surface impedance)

        Parameters:
            L: Inductance matrix [H]
            R: Resistance matrix [Ohm] (can be frequency-dependent)
        """
        self._L_LL = np.asarray(L)
        self._R_LL = np.asarray(R)
        self._n_loops = L.shape[0]

    def set_conductor_star(self, P_cc: np.ndarray,
                           L_Lc: Optional[np.ndarray] = None):
        """
        Set conductor Star (charge) matrices.

        The potential coefficient P_cc relates surface charge to potential:
            φ = P_cc * σ

        The Loop-Star coupling L_Lc represents the inductive coupling
        between current loops and charge distributions.

        Parameters:
            P_cc: Conductor potential coefficient [1/F]
            L_Lc: Loop-Star inductive coupling [H] (optional)
        """
        self._P_cc = np.asarray(P_cc)
        self._n_stars_c = P_cc.shape[0]
        self._L_Lc = np.asarray(L_Lc) if L_Lc is not None else None

    def set_dielectric_star(self, P_dd: np.ndarray,
                            P_cd: Optional[np.ndarray] = None):
        """
        Set dielectric Star matrices.

        Parameters:
            P_dd: Dielectric potential coefficient [1/F]
            P_cd: Conductor-dielectric coupling [1/F]
        """
        self._P_dd = np.asarray(P_dd)
        self._n_stars_d = P_dd.shape[0]
        self._P_cd = np.asarray(P_cd) if P_cd is not None else None

    def set_magnetic_material(self, N: np.ndarray, mu_r: complex = 1000.0,
                               L_Lm: Optional[np.ndarray] = None):
        """
        Set magnetic material matrices.

        The demagnetization matrix N relates magnetization to internal field:
            H_demag = -N * M

        The susceptibility χ = μ_r - 1 relates magnetization to applied field:
            M = χ * H_applied

        Parameters:
            N: Demagnetization matrix [dimensionless]
            mu_r: Complex relative permeability (can include loss: μ' - jμ")
            L_Lm: Loop-magnetic coupling [H/A] (Biot-Savart coefficient)
        """
        self._N_mm = np.asarray(N)
        self._mu_r = mu_r
        self._n_magnetic = N.shape[0]
        self._L_Lm = np.asarray(L_Lm) if L_Lm is not None else None

    def set_skin_effect(self, sigma: float = 5.8e7,
                        model: str = 'dowell',
                        conductor_height: Optional[float] = None,
                        conductor_width: Optional[float] = None,
                        mu_r_conductor: float = 1.0):
        """
        Configure skin effect model for conductor losses.

        For power electronics and WPT, skin effect significantly increases
        AC resistance at high frequencies.

        Two models are available:
        1. Dowell model: For rectangular conductors (foil windings, PCB traces)
           - F_R(xi) = xi * (sinh(2*xi) + sin(2*xi)) / (cosh(2*xi) - cos(2*xi))
           - F_L(xi) = (3/(2*xi)) * (sinh(2*xi) - sin(2*xi)) / (cosh(2*xi) - cos(2*xi))
           - xi = h / delta, delta = sqrt(2 / (omega*mu*sigma))

        2. SIBC model: Surface impedance boundary condition
           - Z_s = (1 + j) / (sigma * delta)
           - Applies to general conductor shapes

        Parameters:
            sigma: Conductivity [S/m] (default: 5.8e7 for copper)
            model: 'dowell' or 'sibc'
            conductor_height: Conductor height [m] (required for Dowell)
            conductor_width: Conductor width [m] (optional, for aspect ratio)
            mu_r_conductor: Relative permeability of conductor (default: 1.0)
        """
        self._use_skin_effect = True
        self._skin_effect_model = model.lower()
        self._sigma = sigma
        self._conductor_height = conductor_height
        self._conductor_width = conductor_width
        self._mu_conductor = 4.0 * np.pi * 1e-7 * mu_r_conductor

    def _compute_skin_depth(self, omega: float) -> float:
        """Compute skin depth delta = sqrt(2 / (omega * mu * sigma))."""
        if omega < 1e-10:
            return 1e10  # Very large at DC
        return np.sqrt(2.0 / (omega * self._mu_conductor * self._sigma))

    def _apply_dowell_correction(self, R_dc: float, L_dc: float,
                                   omega: float) -> Tuple[float, float]:
        """
        Apply Dowell skin effect correction.

        F_R(xi) = xi * (sinh(2*xi) + sin(2*xi)) / (cosh(2*xi) - cos(2*xi))
        F_L(xi) = (3/(2*xi)) * (sinh(2*xi) - sin(2*xi)) / (cosh(2*xi) - cos(2*xi))

        where xi = h / delta, delta = sqrt(2 / (omega*mu*sigma))

        Returns:
            (R_ac, L_ac): AC resistance and inductance
        """
        if self._conductor_height is None or omega < 1e-10:
            return R_dc, L_dc

        delta = self._compute_skin_depth(omega)
        xi = self._conductor_height / delta

        if xi < 0.01:
            # Small argument: F_R -> 1, F_L -> 1
            return R_dc, L_dc
        elif xi > 50:
            # Large argument: F_R -> xi, F_L -> 3/(2*xi)
            F_R = xi
            F_L = 1.5 / xi
        else:
            # General case
            sinh_2xi = np.sinh(2.0 * xi)
            sin_2xi = np.sin(2.0 * xi)
            cosh_2xi = np.cosh(2.0 * xi)
            cos_2xi = np.cos(2.0 * xi)

            denom = cosh_2xi - cos_2xi
            if abs(denom) < 1e-30:
                denom = 1e-30

            F_R = xi * (sinh_2xi + sin_2xi) / denom
            F_L = (1.5 / xi) * (sinh_2xi - sin_2xi) / denom

        R_ac = R_dc * F_R
        L_ac = L_dc * F_L

        return R_ac, L_ac

    def _apply_sibc_correction(self, R_dc: float, omega: float) -> complex:
        """
        Apply SIBC (Surface Impedance Boundary Condition) correction.

        Z_s = (1 + j) / (sigma * delta) = (1 + j) * sqrt(omega * mu / (2 * sigma))

        The AC resistance increase is:
            R_ac / R_dc = (h / delta) for h >> delta (thick conductor)

        Returns:
            Complex surface impedance contribution
        """
        if omega < 1e-10:
            return R_dc

        delta = self._compute_skin_depth(omega)

        # Surface impedance per unit area
        Z_s = (1 + 1j) * np.sqrt(omega * self._mu_conductor / (2.0 * self._sigma))

        # For a conductor with height h:
        # If h >> delta: current flows in skin layer, R_ac ~ R_dc * (h / delta)
        # If h << delta: current uniform, R_ac ~ R_dc
        if self._conductor_height is not None:
            xi = self._conductor_height / delta
            if xi > 1.0:
                # Skin effect regime
                return R_dc * xi * (1 + 1j) / np.sqrt(2)
            else:
                # Transition regime
                return R_dc * (1 + 0.5 * xi ** 2 * 1j)
        else:
            # No height specified, return complex surface impedance
            return R_dc + 1j * omega * self._mu_conductor * delta / 2

    def set_conductor_magnetic_coupling(self, L_Lm: np.ndarray):
        """
        Set the coupling between conductor currents and magnetic material.

        This matrix represents the Biot-Savart law:
            H_at_magnetic = L_Lm * I_conductor

        And reciprocally, the back-EMF from magnetization:
            V_induced = -d/dt(L_Lm^T * M)

        Parameters:
            L_Lm: Current-to-field coupling [A/m per A = 1/m]
        """
        self._L_Lm = np.asarray(L_Lm)

    def build_system_matrix(self, s: complex) -> np.ndarray:
        """
        Build the full system impedance matrix at complex frequency s.

        Returns the block matrix:
            [Z_LL    Z_LS_c   Z_LS_d   Z_LM  ]
            [Z_SL_c  Z_SS_c   Z_cd     Z_cM  ]
            [Z_SL_d  Z_dc     Z_SS_d   Z_dM  ]
            [Z_ML    Z_Mc     Z_Md     Z_MM  ]
        """
        n_L = self._n_loops
        n_Sc = self._n_stars_c
        n_Sd = self._n_stars_d
        n_M = self._n_magnetic
        n_total = n_L + n_Sc + n_Sd + n_M

        Z = np.zeros((n_total, n_total), dtype=complex)

        # Block indices
        i_L = slice(0, n_L)
        i_Sc = slice(n_L, n_L + n_Sc)
        i_Sd = slice(n_L + n_Sc, n_L + n_Sc + n_Sd)
        i_M = slice(n_L + n_Sc + n_Sd, n_total)

        # Z_LL: Conductor Loop impedance (with skin effect if enabled)
        if self._L_LL is not None:
            if self._R_LL.ndim == 1:
                R_dc = np.diag(self._R_LL)
            else:
                R_dc = self._R_LL

            # Apply skin effect correction if enabled
            if self._use_skin_effect and np.abs(s) > 1e-10:
                omega = np.abs(np.imag(s))

                if self._skin_effect_model == 'dowell':
                    # Apply Dowell correction to each diagonal element
                    R_ac = np.zeros_like(R_dc, dtype=complex)
                    L_ac = np.zeros_like(self._L_LL, dtype=complex)

                    for i in range(n_L):
                        R_dc_i = np.real(R_dc[i, i])
                        L_dc_i = np.real(self._L_LL[i, i])
                        R_i, L_i = self._apply_dowell_correction(R_dc_i, L_dc_i, omega)
                        R_ac[i, i] = R_i
                        L_ac[i, i] = L_i

                        # Off-diagonal elements: no skin effect correction
                        for j in range(n_L):
                            if i != j:
                                R_ac[i, j] = R_dc[i, j]
                                L_ac[i, j] = self._L_LL[i, j]

                    Z[i_L, i_L] = R_ac + s * L_ac

                elif self._skin_effect_model == 'sibc':
                    # Apply SIBC correction
                    Z_sibc = np.zeros_like(R_dc, dtype=complex)
                    for i in range(n_L):
                        R_dc_i = np.real(R_dc[i, i])
                        Z_sibc[i, i] = self._apply_sibc_correction(R_dc_i, omega)
                        for j in range(n_L):
                            if i != j:
                                Z_sibc[i, j] = R_dc[i, j]

                    Z[i_L, i_L] = Z_sibc + s * self._L_LL
                else:
                    # No correction
                    Z[i_L, i_L] = R_dc + s * self._L_LL
            else:
                # No skin effect
                Z[i_L, i_L] = R_dc + s * self._L_LL

        # Z_SS_c: Conductor Star impedance (capacitive)
        if self._P_cc is not None:
            Z[i_Sc, i_Sc] = self._P_cc / s

        # Z_SS_d: Dielectric Star impedance (capacitive)
        if self._P_dd is not None:
            Z[i_Sd, i_Sd] = self._P_dd / s

        # Z_MM: Magnetic material impedance (reluctance formulation)
        #
        # magnetic-material equation: (I + χ*N) * M = χ * H_applied
        # where χ = μ_r - 1 is the susceptibility
        #
        # Rewrite: M = χ * (I + χ*N)^{-1} * H_applied
        #
        # For impedance analogy:
        # - "Current" = M (magnetization)
        # - "Voltage" = H (magnetic field)
        # - Z_MM = (I + χ*N) / χ (reluctance-like)
        #
        # Note: This is NOT multiplied by 1/s because the magnetic response is quasi-static
        # (magnetization responds instantaneously to H in MQS regime)
        if self._N_mm is not None:
            chi = self._mu_r - 1.0
            I_mat = np.eye(n_M)
            if abs(chi) > 1e-10:
                # Standard case: Z_MM = (I + χ*N) / χ
                Z[i_M, i_M] = (I_mat + chi * self._N_mm) / chi
            else:
                # Non-magnetic (χ → 0): M → 0, use large impedance
                Z[i_M, i_M] = I_mat * 1e10

        # Z_LS_c: Loop to conductor Star coupling
        if self._L_Lc is not None:
            Z[i_L, i_Sc] = s * self._L_Lc
            Z[i_Sc, i_L] = s * self._L_Lc.T

        # Z_cd: Conductor Star to Dielectric Star coupling
        if self._P_cd is not None:
            Z[i_Sc, i_Sd] = self._P_cd / s
            Z[i_Sd, i_Sc] = self._P_cd.T / s

        # Z_LM: Loop to Magnetic coupling (Biot-Savart + back-EMF)
        if self._L_Lm is not None:
            # Current creates H-field at magnetic: H = L_Lm * I
            # Magnetization creates flux linkage: ψ = L_Lm^T * M
            # Induced voltage: V = -s * ψ = -s * L_Lm^T * M
            Z[i_L, i_M] = -s * self._L_Lm
            Z[i_M, i_L] = self._L_Lm.T  # H from current

        return Z

    def solve_frequency(self, frequencies: np.ndarray,
                         H_external: Optional[np.ndarray] = None) -> dict:
        """
        Solve for impedance and field quantities at given frequencies.

        The port impedance is computed at the Loop DOFs (conductor ports).
        For a single-port system at Loop[0]:
            Z_port = V_1 / I_1

        The magnetic material block affects impedance through mutual coupling
        (Z_LM block), NOT as a parallel circuit element. The coupling enhances
        inductance via: L_eff = L_self + L_mutual(magnetization).

        For multiple parallel branches (Loop + Star_c + Star_d):
            Y_port = sum of admittances at port nodes only

        Parameters:
            frequencies: Array of frequencies [Hz]
            H_external: External magnetic field [A/m] (optional)

        Returns:
            dict with frequency response data
        """
        n_L = self._n_loops
        n_Sc = self._n_stars_c
        n_Sd = self._n_stars_d
        n_M = self._n_magnetic
        n_total = n_L + n_Sc + n_Sd + n_M

        n_freq = len(frequencies)
        Z_result = np.zeros(n_freq, dtype=complex)
        I_L_result = np.zeros((n_freq, max(n_L, 1)), dtype=complex)
        M_result = np.zeros((n_freq, max(n_M, 1)), dtype=complex)

        for idx, f in enumerate(frequencies):
            omega = 2.0 * np.pi * f
            s = 1j * omega if f > 1e-10 else 1e-10

            Z = self.build_system_matrix(s)

            # Compute admittance
            try:
                Y = np.linalg.inv(Z)
            except np.linalg.LinAlgError:
                Z_result[idx] = Z[0, 0]
                continue

            # Port impedance computation
            #
            # The system has 4 blocks: [Loop, Star_c, Star_d, Magnetic]
            # The electrical port is at the Loop DOFs only.
            #
            # For port impedance, we need Z_port = V_1 / I_1 where
            # the excitation is V_1 = 1V at Loop[0] and V = 0 elsewhere.
            #
            # From Z * X = V, we get X = Z^{-1} * V = Y * V
            # Since V = [1, 0, 0, ...], the current at port is I_1 = Y[0,0]
            # And Z_port = V_1 / I_1 = 1 / Y[0,0]
            #
            # However, for PARALLEL connection of Loop with Star (capacitor),
            # we need to consider that Loop and Star share the same port voltage.
            # The correct formula depends on the circuit topology.

            if n_L >= 1:
                # Check if we have a simple parallel LC circuit
                # (1 Loop + Star branches connected in parallel)
                if n_L == 1 and (n_Sc + n_Sd) > 0 and n_M == 0:
                    # Pure electrical system: parallel connection
                    # Sum admittances from Loop and Star DOFs only
                    n_electric = n_L + n_Sc + n_Sd
                    Y_port = np.sum(Y[:n_electric, :n_electric])
                    Z_result[idx] = 1.0 / Y_port if abs(Y_port) > 1e-30 else np.inf

                elif n_L == 1 and n_M > 0 and (n_Sc + n_Sd) == 0:
                    # Pure magnetic coupling (inductor + magnetic core)
                    # The magnetic material is NOT in parallel - it couples
                    # through Z_LM which modifies effective inductance
                    #
                    # Solve for I_1 given V_1 = 1:
                    # From row 0: Z_LL * I_L + Z_LM * M = V_L = 1
                    # From row M: Z_ML * I_L + Z_MM * M = 0 (no external H)
                    #
                    # Solve for M: M = -Z_MM^{-1} * Z_ML * I_L
                    # Substitute: (Z_LL - Z_LM * Z_MM^{-1} * Z_ML) * I_L = 1
                    # Z_eff = Z_LL - Z_LM * Z_MM^{-1} * Z_ML (Schur complement)
                    #
                    # For single DOF: Z_eff = Z_LL[0,0] - Z_LM @ inv(Z_MM) @ Z_ML
                    i_L = slice(0, n_L)
                    i_M = slice(n_L, n_total)

                    Z_LL = Z[i_L, i_L]
                    Z_LM = Z[i_L, i_M]
                    Z_ML = Z[i_M, i_L]
                    Z_MM = Z[i_M, i_M]

                    try:
                        Z_MM_inv = np.linalg.inv(Z_MM)
                        Z_eff = Z_LL - Z_LM @ Z_MM_inv @ Z_ML
                        Z_result[idx] = Z_eff[0, 0]
                    except np.linalg.LinAlgError:
                        Z_result[idx] = Z_LL[0, 0]

                elif n_L == 1 and (n_Sc + n_Sd) > 0 and n_M > 0:
                    # Full system: Loop + Star + Magnetic
                    # Use Schur complement to eliminate magnetic DOFs first,
                    # then compute parallel admittance of electrical DOFs
                    n_electric = n_L + n_Sc + n_Sd
                    i_E = slice(0, n_electric)
                    i_M = slice(n_electric, n_total)

                    Z_EE = Z[i_E, i_E]
                    Z_EM = Z[i_E, i_M]
                    Z_ME = Z[i_M, i_E]
                    Z_MM = Z[i_M, i_M]

                    try:
                        Z_MM_inv = np.linalg.inv(Z_MM)
                        # Effective electrical impedance after eliminating M
                        Z_eff_electric = Z_EE - Z_EM @ Z_MM_inv @ Z_ME

                        # Now compute port admittance (parallel connection)
                        Y_eff = np.linalg.inv(Z_eff_electric)
                        Y_port = np.sum(Y_eff)
                        Z_result[idx] = 1.0 / Y_port if abs(Y_port) > 1e-30 else np.inf
                    except np.linalg.LinAlgError:
                        Z_result[idx] = Z[0, 0]

                else:
                    # General case: use Y[0,0]
                    Z_result[idx] = 1.0 / Y[0, 0] if abs(Y[0, 0]) > 1e-30 else np.inf

            # Solve for currents and magnetization with unit excitation
            V = np.zeros(n_total, dtype=complex)
            if n_L > 0:
                V[0] = 1.0  # Unit voltage excitation

            try:
                X = np.linalg.solve(Z, V)
                if n_L > 0:
                    I_L_result[idx, :n_L] = X[:n_L]
                if n_M > 0:
                    M_result[idx, :n_M] = X[n_L + n_Sc + n_Sd:]
            except np.linalg.LinAlgError:
                pass

        # Compute derived quantities
        R_result = np.real(Z_result)
        X_result = np.imag(Z_result)
        L_result = np.zeros(n_freq)
        for i, f in enumerate(frequencies):
            if f > 1e-10:
                L_result[i] = X_result[i] / (2.0 * np.pi * f)

        return {
            'frequencies': frequencies,
            'impedance': Z_result,
            'resistance': R_result,
            'reactance': X_result,
            'inductance': L_result,
            'loop_currents': I_L_result,
            'magnetization': M_result
        }


class CombinedLoopStarDielectric:
    """
    Combined Loop-Star-Dielectric solver.

    Integrates:
    - Loop (conductor currents, inductive)
    - Star_c (conductor charges, capacitive)
    - Star_d (dielectric charges, capacitive)

    System equation:
        [Z_LL    Z_LS_c   Z_LS_d ] [I_L  ]   [V_L  ]
        [Z_SL_c  Z_SS_c   Z_cd   ] [Q_c  ] = [V_S_c]
        [Z_SL_d  Z_dc     Z_SS_d ] [Q_d  ]   [V_S_d]

    where all Star blocks are P/s (potential coefficient / frequency).
    """

    def __init__(self):
        # Loop matrices (from PEEC)
        self._L_LL = None
        self._R_LL = None
        self._n_loops = 0

        # Conductor Star matrices
        self._P_cc = None  # Conductor-conductor P
        self._L_Lc = None  # Loop-conductor Star coupling
        self._n_stars_c = 0

        # Dielectric Star matrices
        self._P_dd = None  # Dielectric-dielectric P
        self._P_cd = None  # Conductor-dielectric coupling
        self._L_Ld = None  # Loop-dielectric coupling (usually zero)
        self._n_stars_d = 0

    def set_loop_matrices(self, L_LL: np.ndarray, R_LL: np.ndarray,
                          L_Lc: Optional[np.ndarray] = None):
        """Set Loop (inductive) matrices from PEEC."""
        self._L_LL = np.asarray(L_LL)
        self._R_LL = np.asarray(R_LL)
        self._n_loops = L_LL.shape[0]
        self._L_Lc = np.asarray(L_Lc) if L_Lc is not None else None

    def set_conductor_star(self, P_cc: np.ndarray):
        """Set conductor Star (capacitive) P-matrix."""
        self._P_cc = np.asarray(P_cc)
        self._n_stars_c = P_cc.shape[0]

    def set_dielectric_star(self, P_dd: np.ndarray,
                            P_cd: Optional[np.ndarray] = None):
        """
        Set dielectric Star P-matrices.

        Parameters:
            P_dd: Dielectric-dielectric P-matrix
            P_cd: Conductor-dielectric coupling P-matrix
        """
        self._P_dd = np.asarray(P_dd)
        self._n_stars_d = P_dd.shape[0]
        self._P_cd = np.asarray(P_cd) if P_cd is not None else None

    def solve_frequency(self, frequencies: np.ndarray) -> dict:
        """
        Solve for impedance at given frequencies.

        Returns:
            dict with 'frequencies', 'impedance', 'resistance', 'inductance'
        """
        n_L = self._n_loops
        n_Sc = self._n_stars_c
        n_Sd = self._n_stars_d
        n_total = n_L + n_Sc + n_Sd

        n_freq = len(frequencies)
        Z_result = np.zeros(n_freq, dtype=complex)

        for idx, f in enumerate(frequencies):
            omega = 2.0 * np.pi * f
            s = 1j * omega if f > 1e-10 else 1e-10

            # Build full impedance matrix
            Z_full = np.zeros((n_total, n_total), dtype=complex)

            # Z_LL block
            if self._R_LL.ndim == 1:
                R_LL = np.diag(self._R_LL)
            else:
                R_LL = self._R_LL
            Z_full[:n_L, :n_L] = R_LL + s * self._L_LL

            # Z_SS_c block (conductor Star)
            if self._P_cc is not None:
                Z_full[n_L:n_L+n_Sc, n_L:n_L+n_Sc] = self._P_cc / s

            # Z_SS_d block (dielectric Star)
            if self._P_dd is not None:
                Z_full[n_L+n_Sc:, n_L+n_Sc:] = self._P_dd / s

            # Z_cd coupling (conductor-dielectric)
            if self._P_cd is not None:
                Z_full[n_L:n_L+n_Sc, n_L+n_Sc:] = self._P_cd / s
                Z_full[n_L+n_Sc:, n_L:n_L+n_Sc] = self._P_cd.T / s

            # Loop-Star coupling
            if self._L_Lc is not None:
                Z_full[:n_L, n_L:n_L+n_Sc] = s * self._L_Lc
                Z_full[n_L:n_L+n_Sc, :n_L] = s * self._L_Lc.T

            # Compute admittance and port impedance
            try:
                Y_full = np.linalg.inv(Z_full)

                # For single port at Loop[0]
                if n_L == 1 and n_Sc + n_Sd > 0:
                    # Parallel connection: sum all admittances
                    Y_total = np.sum(Y_full)
                    Z_result[idx] = 1.0 / Y_total if abs(Y_total) > 1e-30 else np.inf
                else:
                    Z_result[idx] = 1.0 / Y_full[0, 0] if abs(Y_full[0, 0]) > 1e-30 else np.inf

            except np.linalg.LinAlgError:
                Z_result[idx] = Z_full[0, 0]

        # Compute derived quantities
        R_result = np.real(Z_result)
        X_result = np.imag(Z_result)
        L_result = np.zeros(n_freq)
        for i, f in enumerate(frequencies):
            if f > 1e-10:
                L_result[i] = X_result[i] / (2.0 * np.pi * f)

        return {
            'frequencies': frequencies,
            'impedance': Z_result,
            'resistance': R_result,
            'reactance': X_result,
            'inductance': L_result
        }
