"""
RWG-EFIE Solver for Induction Heating Analysis

This module provides a full 3D surface element solver using RWG (Rao-Wilton-Glisson)
basis functions and the Electric Field Integral Equation (EFIE).

Unlike the PEEC-based approach which uses linear segments, this solver:
- Supports arbitrary 3D surface meshes
- Handles curved coils and complex workpiece geometries
- Accurately captures eddy current distribution on surfaces

Key Features:
- RWG basis functions for surface current representation
- MQS (Magneto-Quasi-Static) formulation
- ESIM surface impedance for skin effect modeling
- Coil-workpiece coupling via mutual impedance

Physical Formulation:
    EFIE: n x (E_inc + E_scat) = Z_s * J_s  on conductor surface

    where:
        E_inc = incident electric field (from external sources)
        E_scat = scattered field from surface currents
        Z_s = surface impedance (ESIM)
        J_s = surface current density

    In MQS formulation:
        Z * I = V
        Z_mn = j*omega*L_mn + R_mn

        L_mn = (mu_0/4*pi) * integral{ f_m . f_n / |r-r'| } dA dA'
        R_mn = (1/sigma*delta) * integral{ f_m . f_n } dA  (skin effect)

Reference:
    S.M. Rao, D.R. Wilton, A.W. Glisson, "Electromagnetic scattering by surfaces
    of arbitrary shape," IEEE TAP, 1982.

Author: Radia Development Team
Date: 2026-01-08
"""

import numpy as np
from scipy.constants import mu_0
from scipy.sparse.linalg import gmres, bicgstab
import warnings


class RWGTriangle:
    """Triangle element for RWG mesh."""

    def __init__(self, vertices, vertex_indices, normal=None):
        """
        Initialize a triangle.

        Parameters:
            vertices: 3x3 array of vertex coordinates [[x1,y1,z1], [x2,y2,z2], [x3,y3,z3]]
            vertex_indices: [i1, i2, i3] global vertex indices
            normal: Optional normal vector (computed if not provided)
        """
        self.vertices = np.array(vertices)
        self.vertex_indices = list(vertex_indices)

        # Compute area and normal
        v0, v1, v2 = self.vertices
        e1 = v1 - v0
        e2 = v2 - v0
        n = np.cross(e1, e2)
        self.area = np.linalg.norm(n) / 2.0

        if normal is not None:
            self.normal = np.array(normal) / np.linalg.norm(normal)
        elif self.area > 1e-30:
            self.normal = n / (2.0 * self.area)
        else:
            self.normal = np.array([0, 0, 1])

        # Centroid
        self.center = np.mean(self.vertices, axis=0)

        # Edge indices (to be set by mesh)
        self.edges = [-1, -1, -1]
        self.edge_signs = [1, 1, 1]


class RWGEdge:
    """Edge element for RWG basis function."""

    def __init__(self, v1_idx, v2_idx, v1, v2):
        """
        Initialize an edge.

        Parameters:
            v1_idx, v2_idx: Global vertex indices (v1_idx < v2_idx)
            v1, v2: Vertex coordinates
        """
        self.vertex_indices = (min(v1_idx, v2_idx), max(v1_idx, v2_idx))
        self.vertices = (np.array(v1), np.array(v2))
        self.length = np.linalg.norm(v2 - v1)
        self.midpoint = (v1 + v2) / 2.0
        self.direction = (v2 - v1) / self.length if self.length > 1e-15 else np.array([1, 0, 0])

        # Adjacent triangles
        self.tri_plus = -1   # T+ triangle index
        self.tri_minus = -1  # T- triangle index
        self.local_idx_plus = -1   # Local edge index in T+
        self.local_idx_minus = -1  # Local edge index in T-

        # Boundary flag
        self.is_boundary = True  # Set to False when second triangle is found


class RWGMesh:
    """
    RWG mesh for EFIE discretization.

    Manages triangular surface mesh with RWG basis function support.
    """

    def __init__(self):
        self.vertices = []
        self.triangles = []
        self.edges = []
        self.num_interior_edges = 0

        # Material properties per region (for multi-material support)
        self.conductivity = 5.8e7  # Default: copper

    def create_rectangular_plate(self, center, Lx, Ly, normal, nx=10, ny=10):
        """
        Create a rectangular plate mesh (open surface for workpiece).

        Parameters:
            center: [x, y, z] plate center
            Lx, Ly: Plate dimensions
            normal: Plate normal direction
            nx, ny: Number of divisions
        """
        self.vertices = []
        self.triangles = []

        center = np.array(center)
        normal = np.array(normal)
        normal = normal / np.linalg.norm(normal)

        # Create local coordinate system
        if abs(normal[2]) < 0.9:
            u = np.cross(normal, [0, 0, 1])
        else:
            u = np.cross(normal, [1, 0, 0])
        u = u / np.linalg.norm(u)
        v = np.cross(normal, u)

        # Generate vertices
        dx = Lx / nx
        dy = Ly / ny

        for j in range(ny + 1):
            for i in range(nx + 1):
                xi = -Lx / 2.0 + i * dx
                yi = -Ly / 2.0 + j * dy
                p = center + xi * u + yi * v
                self.vertices.append(p)

        # Generate triangles
        nVx = nx + 1
        for j in range(ny):
            for i in range(nx):
                v00 = j * nVx + i
                v10 = j * nVx + (i + 1)
                v01 = (j + 1) * nVx + i
                v11 = (j + 1) * nVx + (i + 1)

                # Two triangles per quad
                verts1 = [self.vertices[v00], self.vertices[v10], self.vertices[v11]]
                tri1 = RWGTriangle(verts1, [v00, v10, v11])
                self.triangles.append(tri1)

                verts2 = [self.vertices[v00], self.vertices[v11], self.vertices[v01]]
                tri2 = RWGTriangle(verts2, [v00, v11, v01])
                self.triangles.append(tri2)

        self._build_edge_connectivity()

    def create_circular_disk(self, center, radius, normal, nr=5, ntheta=16):
        """
        Create a circular disk mesh (open surface).

        Parameters:
            center: Disk center
            radius: Disk radius
            normal: Disk normal direction
            nr: Number of radial divisions
            ntheta: Number of angular divisions
        """
        self.vertices = []
        self.triangles = []

        center = np.array(center)
        normal = np.array(normal)
        normal = normal / np.linalg.norm(normal)

        # Create local coordinate system
        if abs(normal[2]) < 0.9:
            u = np.cross(normal, [0, 0, 1])
        else:
            u = np.cross(normal, [1, 0, 0])
        u = u / np.linalg.norm(u)
        v = np.cross(normal, u)

        # Center vertex
        self.vertices.append(center.copy())

        # Ring vertices
        dr = radius / nr
        dtheta = 2.0 * np.pi / ntheta

        for r in range(1, nr + 1):
            ri = r * dr
            for t in range(ntheta):
                theta = t * dtheta
                p = center + ri * (np.cos(theta) * u + np.sin(theta) * v)
                self.vertices.append(p)

        # Inner triangles (center to first ring)
        for t in range(ntheta):
            t_next = (t + 1) % ntheta
            verts = [self.vertices[0], self.vertices[1 + t], self.vertices[1 + t_next]]
            tri = RWGTriangle(verts, [0, 1 + t, 1 + t_next])
            self.triangles.append(tri)

        # Ring triangles
        for r in range(1, nr):
            ring0_start = 1 + (r - 1) * ntheta
            ring1_start = 1 + r * ntheta

            for t in range(ntheta):
                t_next = (t + 1) % ntheta

                # Triangle 1
                v0, v1, v2 = ring0_start + t, ring1_start + t, ring1_start + t_next
                verts1 = [self.vertices[v0], self.vertices[v1], self.vertices[v2]]
                tri1 = RWGTriangle(verts1, [v0, v1, v2])
                self.triangles.append(tri1)

                # Triangle 2
                v0, v1, v2 = ring0_start + t, ring1_start + t_next, ring0_start + t_next
                verts2 = [self.vertices[v0], self.vertices[v1], self.vertices[v2]]
                tri2 = RWGTriangle(verts2, [v0, v1, v2])
                self.triangles.append(tri2)

        self._build_edge_connectivity()

    def create_spiral_coil(self, center, inner_radius, outer_radius, pitch, num_turns,
                           axis, wire_radius, num_around=8):
        """
        Create a spiral coil mesh (closed surface).

        Parameters:
            center: Coil center
            inner_radius, outer_radius: Spiral radii
            pitch: Axial advance per turn
            num_turns: Number of turns
            axis: Coil axis direction
            wire_radius: Wire cross-section radius
            num_around: Panels around wire circumference
        """
        self.vertices = []
        self.triangles = []

        center = np.array(center)
        axis = np.array(axis)
        axis = axis / np.linalg.norm(axis)

        # Create local coordinate system
        if abs(axis[2]) < 0.9:
            u = np.cross(axis, [0, 0, 1])
        else:
            u = np.cross(axis, [1, 0, 0])
        u = u / np.linalg.norm(u)
        v = np.cross(axis, u)

        # Generate wire path
        total_angle = 2 * np.pi * num_turns
        n_path_points = int(num_turns * 30)  # 30 points per turn
        path = []

        for i in range(n_path_points + 1):
            t = i / n_path_points
            phi = t * total_angle
            R = inner_radius + t * (outer_radius - inner_radius)
            z = t * num_turns * pitch

            p = center + R * (np.cos(phi) * u + np.sin(phi) * v) + z * axis
            path.append(p)

        # Generate wire mesh along path
        is_closed = False  # Spiral is open
        n_rings = len(path)

        for ring, p in enumerate(path):
            # Compute tangent
            if ring == 0:
                tangent = path[1] - path[0]
            elif ring == n_rings - 1:
                tangent = path[-1] - path[-2]
            else:
                tangent = path[ring + 1] - path[ring - 1]

            tangent = tangent / np.linalg.norm(tangent)

            # Find perpendicular directions
            if abs(tangent[2]) < 0.9:
                perp1 = np.cross(tangent, [0, 0, 1])
            else:
                perp1 = np.cross(tangent, [1, 0, 0])
            perp1 = perp1 / np.linalg.norm(perp1)
            perp2 = np.cross(tangent, perp1)

            # Generate ring vertices
            for i in range(num_around):
                theta = 2.0 * np.pi * i / num_around
                vert = p + wire_radius * (np.cos(theta) * perp1 + np.sin(theta) * perp2)
                self.vertices.append(vert)

        # Generate triangles
        for seg in range(n_rings - 1):
            ring0 = seg * num_around
            ring1 = (seg + 1) * num_around

            for i in range(num_around):
                i_next = (i + 1) % num_around

                # Triangle 1
                v0, v1, v2 = ring0 + i, ring1 + i, ring1 + i_next
                verts1 = [self.vertices[v0], self.vertices[v1], self.vertices[v2]]
                tri1 = RWGTriangle(verts1, [v0, v1, v2])
                self.triangles.append(tri1)

                # Triangle 2
                v0, v1, v2 = ring0 + i, ring1 + i_next, ring0 + i_next
                verts2 = [self.vertices[v0], self.vertices[v1], self.vertices[v2]]
                tri2 = RWGTriangle(verts2, [v0, v1, v2])
                self.triangles.append(tri2)

        self._build_edge_connectivity()

    def create_loop_coil(self, center, radius, normal, wire_radius, num_around=8, num_along=30):
        """
        Create a circular loop coil mesh (closed surface).

        Parameters:
            center: Loop center
            radius: Loop radius
            normal: Loop normal direction
            wire_radius: Wire cross-section radius
            num_around: Panels around wire circumference
            num_along: Panels along loop
        """
        self.vertices = []
        self.triangles = []

        center = np.array(center)
        normal = np.array(normal)
        normal = normal / np.linalg.norm(normal)

        # Create local coordinate system
        if abs(normal[2]) < 0.9:
            u = np.cross(normal, [0, 0, 1])
        else:
            u = np.cross(normal, [1, 0, 0])
        u = u / np.linalg.norm(u)
        v = np.cross(normal, u)

        # Generate wire path (closed loop)
        path = []
        for i in range(num_along):
            phi = 2 * np.pi * i / num_along
            p = center + radius * (np.cos(phi) * u + np.sin(phi) * v)
            path.append(p)

        # Generate wire mesh (closed loop)
        n_rings = num_along

        for ring, p in enumerate(path):
            # Compute tangent for closed loop
            ring_next = (ring + 1) % n_rings
            ring_prev = (ring - 1) % n_rings
            tangent = path[ring_next] - path[ring_prev]
            tangent = tangent / np.linalg.norm(tangent)

            # Find perpendicular directions
            if abs(tangent[2]) < 0.9:
                perp1 = np.cross(tangent, [0, 0, 1])
            else:
                perp1 = np.cross(tangent, [1, 0, 0])
            perp1 = perp1 / np.linalg.norm(perp1)
            perp2 = np.cross(tangent, perp1)

            # Generate ring vertices
            for i in range(num_around):
                theta = 2.0 * np.pi * i / num_around
                vert = p + wire_radius * (np.cos(theta) * perp1 + np.sin(theta) * perp2)
                self.vertices.append(vert)

        # Generate triangles (closed loop)
        for seg in range(n_rings):
            ring0 = seg * num_around
            ring1 = ((seg + 1) % n_rings) * num_around

            for i in range(num_around):
                i_next = (i + 1) % num_around

                # Triangle 1
                v0, v1, v2 = ring0 + i, ring1 + i, ring1 + i_next
                verts1 = [self.vertices[v0], self.vertices[v1], self.vertices[v2]]
                tri1 = RWGTriangle(verts1, [v0, v1, v2])
                self.triangles.append(tri1)

                # Triangle 2
                v0, v1, v2 = ring0 + i, ring1 + i_next, ring0 + i_next
                verts2 = [self.vertices[v0], self.vertices[v1], self.vertices[v2]]
                tri2 = RWGTriangle(verts2, [v0, v1, v2])
                self.triangles.append(tri2)

        self._build_edge_connectivity()

    def _build_edge_connectivity(self):
        """Build edge connectivity and RWG basis function support."""
        self.edges = []
        edge_map = {}  # (v1, v2) -> edge index

        for tri_idx, tri in enumerate(self.triangles):
            v_idx = tri.vertex_indices

            # Process 3 edges (edge i is opposite to vertex i)
            for local_edge in range(3):
                v1_idx = v_idx[(local_edge + 1) % 3]
                v2_idx = v_idx[(local_edge + 2) % 3]

                # Canonical edge key (smaller index first)
                key = (min(v1_idx, v2_idx), max(v1_idx, v2_idx))

                if key not in edge_map:
                    # New edge
                    v1 = self.vertices[v1_idx]
                    v2 = self.vertices[v2_idx]
                    edge = RWGEdge(v1_idx, v2_idx, v1, v2)
                    edge.tri_plus = tri_idx
                    edge.local_idx_plus = local_edge

                    edge_idx = len(self.edges)
                    self.edges.append(edge)
                    edge_map[key] = edge_idx

                    tri.edges[local_edge] = edge_idx
                    tri.edge_signs[local_edge] = 1 if v1_idx < v2_idx else -1
                else:
                    # Existing edge - this triangle is T-
                    edge_idx = edge_map[key]
                    edge = self.edges[edge_idx]
                    edge.tri_minus = tri_idx
                    edge.local_idx_minus = local_edge
                    edge.is_boundary = False

                    tri.edges[local_edge] = edge_idx
                    tri.edge_signs[local_edge] = -1 if v1_idx < v2_idx else 1

        # Count interior edges
        self.num_interior_edges = sum(1 for e in self.edges if not e.is_boundary)

    def get_boundary_edges(self):
        """Get indices of boundary edges."""
        return [i for i, e in enumerate(self.edges) if e.is_boundary]

    def is_closed(self):
        """Check if mesh is a closed surface."""
        return self.num_interior_edges == len(self.edges)

    def merge_mesh(self, other):
        """
        Merge another mesh into this mesh.

        Parameters:
            other: RWGMesh to merge
        """
        vertex_offset = len(self.vertices)

        # Add vertices
        self.vertices.extend(other.vertices)

        # Add triangles with offset indices
        for tri in other.triangles:
            new_v_idx = [v + vertex_offset for v in tri.vertex_indices]
            new_verts = [self.vertices[i] for i in new_v_idx]
            new_tri = RWGTriangle(new_verts, new_v_idx, tri.normal)
            self.triangles.append(new_tri)

        # Rebuild connectivity
        self._build_edge_connectivity()


class RWGEFIESolver:
    """
    EFIE Solver using RWG basis functions.

    Solves the MQS (Magneto-Quasi-Static) EFIE:
        Z * I = V

    where:
        Z = j*omega*L + R  (impedance matrix)
        I = current coefficients on interior edges
        V = voltage excitation
    """

    def __init__(self, mesh):
        """
        Initialize solver.

        Parameters:
            mesh: RWGMesh object
        """
        self.mesh = mesh
        self.frequency = 1e6  # Default 1 MHz
        self.omega = 2 * np.pi * self.frequency
        self.conductivity = 5.8e7  # Default copper

        # Matrices
        self.L_matrix = None  # Inductance matrix
        self.R_matrix = None  # Resistance matrix
        self.Z_matrix = None  # Combined impedance
        self.rhs = None       # Right-hand side
        self.solution = None  # Current coefficients

        # Results
        self.impedance = 0j
        self.total_current = 0j

        # Port definition
        self.port_edges = []
        self.voltage_excitation = 1.0 + 0j

    def set_frequency(self, frequency):
        """Set operating frequency [Hz]."""
        self.frequency = frequency
        self.omega = 2 * np.pi * frequency

    def set_conductivity(self, sigma):
        """Set material conductivity [S/m]."""
        self.conductivity = sigma

    def define_port(self, edge_indices):
        """
        Define port for impedance calculation.

        Parameters:
            edge_indices: List of edge indices forming the port
        """
        self.port_edges = list(edge_indices)

    def set_voltage_excitation(self, V):
        """Set voltage excitation at port [V]."""
        self.voltage_excitation = complex(V)

    def _get_interior_edge_mapping(self):
        """Get mapping from interior edges to unknown indices."""
        edge_to_unknown = {}
        unknown_idx = 0
        for i, edge in enumerate(self.mesh.edges):
            if not edge.is_boundary:
                edge_to_unknown[i] = unknown_idx
                unknown_idx += 1
        return edge_to_unknown

    def _compute_L_element(self, edge_m, edge_n):
        """
        Compute inductance matrix element L_mn.

        L_mn = (mu_0/4*pi) * integral{ f_m . f_n / |r-r'| } dA dA'

        Uses numerical quadrature over triangles.
        """
        # 7-point Gauss quadrature for triangle
        quad_pts = np.array([
            [0.333333333333333, 0.333333333333333],
            [0.059715871789770, 0.470142064105115],
            [0.470142064105115, 0.059715871789770],
            [0.470142064105115, 0.470142064105115],
            [0.797426985353087, 0.101286507323456],
            [0.101286507323456, 0.797426985353087],
            [0.101286507323456, 0.101286507323456]
        ])
        quad_wts = np.array([
            0.225000000000000,
            0.132394152788506, 0.132394152788506, 0.132394152788506,
            0.125939180544827, 0.125939180544827, 0.125939180544827
        ])

        L_mn = 0.0

        # Integrate over triangles adjacent to edge_m and edge_n
        tris_m = []
        if edge_m.tri_plus >= 0:
            tris_m.append((edge_m.tri_plus, +1))
        if edge_m.tri_minus >= 0:
            tris_m.append((edge_m.tri_minus, -1))

        tris_n = []
        if edge_n.tri_plus >= 0:
            tris_n.append((edge_n.tri_plus, +1))
        if edge_n.tri_minus >= 0:
            tris_n.append((edge_n.tri_minus, -1))

        for tri_m_idx, sign_m in tris_m:
            tri_m = self.mesh.triangles[tri_m_idx]
            v_m = np.array([self.mesh.vertices[i] for i in tri_m.vertex_indices])

            # Free vertex for RWG in T_m
            if sign_m > 0:
                free_m = self.mesh.vertices[tri_m.vertex_indices[edge_m.local_idx_plus]]
            else:
                free_m = self.mesh.vertices[tri_m.vertex_indices[edge_m.local_idx_minus]]

            coeff_m = edge_m.length / (2.0 * tri_m.area)

            for tri_n_idx, sign_n in tris_n:
                tri_n = self.mesh.triangles[tri_n_idx]
                v_n = np.array([self.mesh.vertices[i] for i in tri_n.vertex_indices])

                # Free vertex for RWG in T_n
                if sign_n > 0:
                    free_n = self.mesh.vertices[tri_n.vertex_indices[edge_n.local_idx_plus]]
                else:
                    free_n = self.mesh.vertices[tri_n.vertex_indices[edge_n.local_idx_minus]]

                coeff_n = edge_n.length / (2.0 * tri_n.area)

                # Double quadrature
                for q_m in range(len(quad_wts)):
                    L1_m, L2_m = quad_pts[q_m]
                    L3_m = 1 - L1_m - L2_m
                    r_m = L1_m * v_m[0] + L2_m * v_m[1] + L3_m * v_m[2]

                    # RWG value at r_m
                    if sign_m > 0:
                        f_m = coeff_m * (r_m - free_m)
                    else:
                        f_m = coeff_m * (free_m - r_m)

                    for q_n in range(len(quad_wts)):
                        L1_n, L2_n = quad_pts[q_n]
                        L3_n = 1 - L1_n - L2_n
                        r_n = L1_n * v_n[0] + L2_n * v_n[1] + L3_n * v_n[2]

                        # RWG value at r_n
                        if sign_n > 0:
                            f_n = coeff_n * (r_n - free_n)
                        else:
                            f_n = coeff_n * (free_n - r_n)

                        # Green's function 1/R
                        R = np.linalg.norm(r_m - r_n)
                        if R < 1e-12:
                            # Self-term: use regularized Green's function
                            R = 0.5 * np.sqrt(tri_m.area / np.pi)

                        G = 1.0 / R

                        # Integrand: f_m . f_n * G
                        integrand = np.dot(f_m, f_n) * G

                        # Weighted contribution
                        w = quad_wts[q_m] * quad_wts[q_n] * tri_m.area * tri_n.area
                        L_mn += w * integrand

        # Scale by mu_0 / (4*pi)
        L_mn *= mu_0 / (4 * np.pi)

        return L_mn

    def _compute_R_element(self, edge_m, edge_n):
        """
        Compute resistance matrix element R_mn.

        For skin effect: R_mn = R_s * integral{ f_m . f_n } dA
        where R_s = 1/(sigma * delta) is the surface resistance.
        """
        # Skin depth
        delta = np.sqrt(2.0 / (self.omega * mu_0 * self.conductivity))
        R_s = 1.0 / (self.conductivity * delta)

        # Only non-zero if edges share a triangle
        common_tris = []
        if edge_m.tri_plus >= 0:
            if edge_m.tri_plus == edge_n.tri_plus or edge_m.tri_plus == edge_n.tri_minus:
                common_tris.append(edge_m.tri_plus)
        if edge_m.tri_minus >= 0:
            if edge_m.tri_minus == edge_n.tri_plus or edge_m.tri_minus == edge_n.tri_minus:
                common_tris.append(edge_m.tri_minus)

        if not common_tris:
            return 0.0

        # Quadrature
        quad_pts = np.array([
            [0.333333333333333, 0.333333333333333],
            [0.059715871789770, 0.470142064105115],
            [0.470142064105115, 0.059715871789770],
            [0.470142064105115, 0.470142064105115],
        ])
        quad_wts = np.array([0.225, 0.132394, 0.132394, 0.132394])

        R_mn = 0.0

        for tri_idx in common_tris:
            tri = self.mesh.triangles[tri_idx]
            v = np.array([self.mesh.vertices[i] for i in tri.vertex_indices])

            # Determine signs for this triangle
            sign_m = 1 if tri_idx == edge_m.tri_plus else -1
            sign_n = 1 if tri_idx == edge_n.tri_plus else -1

            # Free vertices
            if sign_m > 0:
                free_m = self.mesh.vertices[tri.vertex_indices[edge_m.local_idx_plus]]
            else:
                free_m = self.mesh.vertices[tri.vertex_indices[edge_m.local_idx_minus]]

            if sign_n > 0:
                free_n = self.mesh.vertices[tri.vertex_indices[edge_n.local_idx_plus]]
            else:
                free_n = self.mesh.vertices[tri.vertex_indices[edge_n.local_idx_minus]]

            coeff_m = edge_m.length / (2.0 * tri.area)
            coeff_n = edge_n.length / (2.0 * tri.area)

            for q in range(len(quad_wts)):
                L1, L2 = quad_pts[q]
                L3 = 1 - L1 - L2
                r = L1 * v[0] + L2 * v[1] + L3 * v[2]

                # RWG values
                if sign_m > 0:
                    f_m = coeff_m * (r - free_m)
                else:
                    f_m = coeff_m * (free_m - r)

                if sign_n > 0:
                    f_n = coeff_n * (r - free_n)
                else:
                    f_n = coeff_n * (free_n - r)

                integrand = np.dot(f_m, f_n)
                R_mn += quad_wts[q] * tri.area * integrand

        return R_s * R_mn

    def assemble_matrices(self):
        """Assemble impedance matrices."""
        edge_to_unknown = self._get_interior_edge_mapping()
        N = self.mesh.num_interior_edges

        if N == 0:
            warnings.warn("No interior edges found. Cannot solve.")
            return

        self.L_matrix = np.zeros((N, N), dtype=complex)
        self.R_matrix = np.zeros((N, N), dtype=complex)

        # Get interior edges
        interior_edges = [(i, e) for i, e in enumerate(self.mesh.edges) if not e.is_boundary]

        for idx_m, (edge_m_idx, edge_m) in enumerate(interior_edges):
            for idx_n, (edge_n_idx, edge_n) in enumerate(interior_edges):
                if idx_n >= idx_m:  # Exploit symmetry
                    L_mn = self._compute_L_element(edge_m, edge_n)
                    R_mn = self._compute_R_element(edge_m, edge_n)

                    self.L_matrix[idx_m, idx_n] = L_mn
                    self.L_matrix[idx_n, idx_m] = L_mn

                    self.R_matrix[idx_m, idx_n] = R_mn
                    self.R_matrix[idx_n, idx_m] = R_mn

        # Combined impedance matrix: Z = j*omega*L + R
        j_omega = 1j * self.omega
        self.Z_matrix = j_omega * self.L_matrix + self.R_matrix

    def assemble_rhs(self):
        """Assemble right-hand side (voltage excitation)."""
        edge_to_unknown = self._get_interior_edge_mapping()
        N = self.mesh.num_interior_edges

        self.rhs = np.zeros(N, dtype=complex)

        # Apply voltage across port edges
        if self.port_edges:
            # Simple model: distribute voltage across port edges
            V_per_edge = self.voltage_excitation / len(self.port_edges)

            for edge_idx in self.port_edges:
                if edge_idx in edge_to_unknown:
                    unknown_idx = edge_to_unknown[edge_idx]
                    edge = self.mesh.edges[edge_idx]
                    self.rhs[unknown_idx] = V_per_edge * edge.length

    def solve(self):
        """Solve the EFIE system."""
        self.assemble_matrices()
        self.assemble_rhs()

        N = self.mesh.num_interior_edges
        if N == 0:
            return

        # Solve Z * I = V
        try:
            self.solution = np.linalg.solve(self.Z_matrix, self.rhs)
        except np.linalg.LinAlgError:
            # Fallback to GMRES for singular systems
            self.solution, _ = gmres(self.Z_matrix, self.rhs, tol=1e-6)

        # Compute port impedance
        if self.port_edges:
            edge_to_unknown = self._get_interior_edge_mapping()

            # Total current through port
            I_total = 0j
            for edge_idx in self.port_edges:
                if edge_idx in edge_to_unknown:
                    unknown_idx = edge_to_unknown[edge_idx]
                    I_total += self.solution[unknown_idx]

            self.total_current = I_total
            if abs(I_total) > 1e-30:
                self.impedance = self.voltage_excitation / I_total
            else:
                self.impedance = np.inf + 0j

    def get_impedance(self):
        """Get port impedance [Ohm]."""
        return self.impedance

    def get_current(self):
        """Get total port current [A]."""
        return self.total_current

    def compute_B_field(self, point):
        """
        Compute magnetic field B at a point.

        B = curl(A) = (mu_0/4*pi) * integral{ J x R / R^3 } dA

        Parameters:
            point: [x, y, z] observation point

        Returns:
            B: [Bx, By, Bz] magnetic field [T]
        """
        if self.solution is None:
            return np.zeros(3, dtype=complex)

        point = np.array(point)
        B = np.zeros(3, dtype=complex)

        edge_to_unknown = self._get_interior_edge_mapping()

        # Quadrature points
        quad_pts = np.array([
            [0.333333333333333, 0.333333333333333],
            [0.059715871789770, 0.470142064105115],
            [0.470142064105115, 0.059715871789770],
        ])
        quad_wts = np.array([0.225, 0.132394, 0.132394])

        for edge_idx, edge in enumerate(self.mesh.edges):
            if edge.is_boundary:
                continue

            unknown_idx = edge_to_unknown[edge_idx]
            I_n = self.solution[unknown_idx]

            # Integrate over adjacent triangles
            tris = []
            if edge.tri_plus >= 0:
                tris.append((edge.tri_plus, +1))
            if edge.tri_minus >= 0:
                tris.append((edge.tri_minus, -1))

            for tri_idx, sign in tris:
                tri = self.mesh.triangles[tri_idx]
                v = np.array([self.mesh.vertices[i] for i in tri.vertex_indices])

                # Free vertex
                if sign > 0:
                    free = self.mesh.vertices[tri.vertex_indices[edge.local_idx_plus]]
                else:
                    free = self.mesh.vertices[tri.vertex_indices[edge.local_idx_minus]]

                coeff = edge.length / (2.0 * tri.area)

                for q in range(len(quad_wts)):
                    L1, L2 = quad_pts[q]
                    L3 = 1 - L1 - L2
                    r_src = L1 * v[0] + L2 * v[1] + L3 * v[2]

                    # RWG value (surface current direction)
                    if sign > 0:
                        J = coeff * (r_src - free)
                    else:
                        J = coeff * (free - r_src)

                    J *= I_n  # Scale by current coefficient

                    # Biot-Savart kernel
                    R_vec = point - r_src
                    R = np.linalg.norm(R_vec)

                    if R > 1e-12:
                        dB = np.cross(J, R_vec) / (R ** 3)
                        B += quad_wts[q] * tri.area * dB

        # Scale by mu_0 / (4*pi)
        B *= mu_0 / (4 * np.pi)

        return B


class CoupledEFIESolver:
    """
    Coupled coil-workpiece EFIE solver.

    Solves the coupled system:
        [Z_cc  Z_cw] [I_c]   [V_c]
        [Z_wc  Z_ww] [I_w] = [0  ]

    where:
        Z_cc = coil self-impedance
        Z_ww = workpiece self-impedance
        Z_cw, Z_wc = mutual impedance (coil-workpiece coupling)
        I_c = coil current
        I_w = workpiece eddy current
        V_c = applied voltage
    """

    def __init__(self, coil_mesh, workpiece_mesh):
        """
        Initialize coupled solver.

        Parameters:
            coil_mesh: RWGMesh for coil
            workpiece_mesh: RWGMesh for workpiece
        """
        self.coil_mesh = coil_mesh
        self.workpiece_mesh = workpiece_mesh

        self.frequency = 1e6
        self.omega = 2 * np.pi * self.frequency
        self.coil_conductivity = 5.8e7   # Copper
        self.workpiece_conductivity = 5e6  # Steel
        self.workpiece_mu_r = 1.0  # Relative permeability

        # Matrices
        self.Z_cc = None  # Coil-coil
        self.Z_ww = None  # Workpiece-workpiece
        self.Z_cw = None  # Coil-workpiece
        self.Z_wc = None  # Workpiece-coil

        # Solution
        self.I_coil = None
        self.I_workpiece = None
        self.impedance = 0j
        self.power_workpiece = 0.0

        # Port
        self.port_edges = []
        self.voltage_excitation = 1.0 + 0j

    def set_frequency(self, frequency):
        """Set operating frequency [Hz]."""
        self.frequency = frequency
        self.omega = 2 * np.pi * frequency

    def set_coil_conductivity(self, sigma):
        """Set coil conductivity [S/m]."""
        self.coil_conductivity = sigma

    def set_workpiece_conductivity(self, sigma):
        """Set workpiece conductivity [S/m]."""
        self.workpiece_conductivity = sigma

    def set_workpiece_permeability(self, mu_r):
        """Set workpiece relative permeability."""
        self.workpiece_mu_r = mu_r

    def define_port(self, edge_indices):
        """Define port edges for coil excitation."""
        self.port_edges = list(edge_indices)

    def set_voltage_excitation(self, V):
        """Set voltage excitation [V]."""
        self.voltage_excitation = complex(V)

    def _get_skin_depth(self, sigma, mu_r=1.0):
        """Calculate skin depth."""
        if self.omega > 0 and sigma > 0:
            return np.sqrt(2.0 / (self.omega * mu_0 * mu_r * sigma))
        return 1e-3

    def _compute_mutual_L(self, mesh1, mesh2, edge1_idx, edge2_idx):
        """
        Compute mutual inductance element between edges in different meshes.

        L_mn = (mu_0/4*pi) * integral{ f_m . f_n / |r-r'| } dA dA'
        """
        edge1 = mesh1.edges[edge1_idx]
        edge2 = mesh2.edges[edge2_idx]

        if edge1.is_boundary or edge2.is_boundary:
            return 0j

        # Quadrature points
        quad_pts = np.array([
            [0.333333333333333, 0.333333333333333],
            [0.059715871789770, 0.470142064105115],
            [0.470142064105115, 0.059715871789770],
            [0.470142064105115, 0.470142064105115],
        ])
        quad_wts = np.array([0.225, 0.132394, 0.132394, 0.132394])

        L_mn = 0.0

        # Triangles for edge1
        tris1 = []
        if edge1.tri_plus >= 0:
            tris1.append((edge1.tri_plus, +1, edge1.local_idx_plus))
        if edge1.tri_minus >= 0:
            tris1.append((edge1.tri_minus, -1, edge1.local_idx_minus))

        # Triangles for edge2
        tris2 = []
        if edge2.tri_plus >= 0:
            tris2.append((edge2.tri_plus, +1, edge2.local_idx_plus))
        if edge2.tri_minus >= 0:
            tris2.append((edge2.tri_minus, -1, edge2.local_idx_minus))

        for tri1_idx, sign1, local1 in tris1:
            tri1 = mesh1.triangles[tri1_idx]
            v1 = np.array([mesh1.vertices[i] for i in tri1.vertex_indices])
            free1 = mesh1.vertices[tri1.vertex_indices[local1]]
            coeff1 = sign1 * edge1.length / (2.0 * tri1.area)

            for tri2_idx, sign2, local2 in tris2:
                tri2 = mesh2.triangles[tri2_idx]
                v2 = np.array([mesh2.vertices[i] for i in tri2.vertex_indices])
                free2 = mesh2.vertices[tri2.vertex_indices[local2]]
                coeff2 = sign2 * edge2.length / (2.0 * tri2.area)

                # Double integration
                for q1 in range(len(quad_wts)):
                    L1_1, L2_1 = quad_pts[q1]
                    L3_1 = 1 - L1_1 - L2_1
                    r1 = L1_1 * v1[0] + L2_1 * v1[1] + L3_1 * v1[2]

                    # RWG basis at r1
                    if sign1 > 0:
                        f1 = coeff1 * (r1 - free1)
                    else:
                        f1 = -coeff1 * (free1 - r1)

                    for q2 in range(len(quad_wts)):
                        L1_2, L2_2 = quad_pts[q2]
                        L3_2 = 1 - L1_2 - L2_2
                        r2 = L1_2 * v2[0] + L2_2 * v2[1] + L3_2 * v2[2]

                        # RWG basis at r2
                        if sign2 > 0:
                            f2 = coeff2 * (r2 - free2)
                        else:
                            f2 = -coeff2 * (free2 - r2)

                        # Green's function
                        R = np.linalg.norm(r1 - r2)
                        if R > 1e-12:
                            G = 1.0 / R
                        else:
                            G = 0.0

                        # Integrand
                        integrand = np.dot(f1, f2) * G
                        L_mn += quad_wts[q1] * quad_wts[q2] * tri1.area * tri2.area * integrand

        # Scale by mu_0 / (4*pi)
        L_mn *= mu_0 / (4 * np.pi)

        return complex(L_mn, 0)

    def _assemble_self_impedance(self, mesh, sigma, mu_r=1.0):
        """Assemble self-impedance matrix for a mesh."""
        # Get interior edges
        interior_edges = [i for i, e in enumerate(mesh.edges) if not e.is_boundary]
        N = len(interior_edges)

        if N == 0:
            return np.zeros((0, 0), dtype=complex), []

        # Create temporary solver to use existing methods
        solver = RWGEFIESolver(mesh)
        solver.set_frequency(self.frequency)
        solver.set_conductivity(sigma)
        solver.assemble_matrices()

        return solver.Z_matrix, interior_edges

    def assemble_system(self):
        """Assemble the coupled EFIE system."""
        # Get interior edges for both meshes
        coil_interior = [i for i, e in enumerate(self.coil_mesh.edges) if not e.is_boundary]
        work_interior = [i for i, e in enumerate(self.workpiece_mesh.edges) if not e.is_boundary]

        N_c = len(coil_interior)
        N_w = len(work_interior)

        print(f"  Coil interior edges: {N_c}")
        print(f"  Workpiece interior edges: {N_w}")
        print(f"  Total DOF: {N_c + N_w}")

        if N_c == 0:
            raise ValueError("Coil has no interior edges")

        # Assemble coil self-impedance
        print("  Assembling coil self-impedance...")
        self.Z_cc, _ = self._assemble_self_impedance(
            self.coil_mesh, self.coil_conductivity
        )

        # Assemble workpiece self-impedance (if has interior edges)
        if N_w > 0:
            print("  Assembling workpiece self-impedance...")
            self.Z_ww, _ = self._assemble_self_impedance(
                self.workpiece_mesh, self.workpiece_conductivity, self.workpiece_mu_r
            )

            # Assemble mutual impedance blocks
            print("  Assembling mutual impedance (coil-workpiece)...")
            self.Z_cw = np.zeros((N_c, N_w), dtype=complex)
            self.Z_wc = np.zeros((N_w, N_c), dtype=complex)

            j_omega = 1j * self.omega

            for i, edge_c_idx in enumerate(coil_interior):
                for j, edge_w_idx in enumerate(work_interior):
                    L_cw = self._compute_mutual_L(
                        self.coil_mesh, self.workpiece_mesh,
                        edge_c_idx, edge_w_idx
                    )
                    self.Z_cw[i, j] = j_omega * L_cw
                    self.Z_wc[j, i] = j_omega * L_cw  # Reciprocity

            # Build full system matrix
            self.Z_full = np.block([
                [self.Z_cc, self.Z_cw],
                [self.Z_wc, self.Z_ww]
            ])
        else:
            # Coil only
            self.Z_full = self.Z_cc
            self.Z_cw = None
            self.Z_wc = None
            self.Z_ww = None

        self._N_c = N_c
        self._N_w = N_w
        self._coil_interior = coil_interior
        self._work_interior = work_interior

    def assemble_rhs(self):
        """Assemble right-hand side vector."""
        N_c = self._N_c
        N_w = self._N_w

        self.rhs = np.zeros(N_c + N_w, dtype=complex)

        # Apply voltage at coil port
        if self.port_edges:
            # Map edge index to unknown index
            edge_to_unknown = {e: i for i, e in enumerate(self._coil_interior)}

            V_per_edge = self.voltage_excitation / len(self.port_edges)

            for edge_idx in self.port_edges:
                if edge_idx in edge_to_unknown:
                    unknown_idx = edge_to_unknown[edge_idx]
                    edge = self.coil_mesh.edges[edge_idx]
                    self.rhs[unknown_idx] = V_per_edge * edge.length

    def solve(self, verbose=True):
        """
        Solve the coupled system.

        Returns:
            result: Dictionary with solution data
        """
        if verbose:
            print("\nCoupled EFIE Solver")
            print("=" * 50)
            print(f"Frequency: {self.frequency/1e3:.1f} kHz")

        # Assemble system
        self.assemble_system()
        self.assemble_rhs()

        # Solve
        if verbose:
            print("  Solving linear system...")

        try:
            solution = np.linalg.solve(self.Z_full, self.rhs)
        except np.linalg.LinAlgError:
            solution, _ = gmres(self.Z_full, self.rhs, tol=1e-6)

        # Extract coil and workpiece currents
        self.I_coil = solution[:self._N_c]
        if self._N_w > 0:
            self.I_workpiece = solution[self._N_c:]
        else:
            self.I_workpiece = np.array([])

        # Compute port impedance
        if self.port_edges:
            edge_to_unknown = {e: i for i, e in enumerate(self._coil_interior)}

            I_total = 0j
            for edge_idx in self.port_edges:
                if edge_idx in edge_to_unknown:
                    I_total += self.I_coil[edge_to_unknown[edge_idx]]

            if abs(I_total) > 1e-30:
                self.impedance = self.voltage_excitation / I_total
            else:
                self.impedance = np.inf + 0j

        # Compute power in workpiece
        if self._N_w > 0 and self.Z_ww is not None:
            # P = 0.5 * Re(I* . Z . I)
            self.power_workpiece = 0.5 * np.real(
                np.conj(self.I_workpiece) @ self.Z_ww @ self.I_workpiece
            )
        else:
            self.power_workpiece = 0.0

        result = {
            'impedance': self.impedance,
            'resistance': self.impedance.real,
            'reactance': self.impedance.imag,
            'inductance': self.impedance.imag / (2 * np.pi * self.frequency) if self.frequency > 0 else 0,
            'power_workpiece': self.power_workpiece,
            'coil_dof': self._N_c,
            'workpiece_dof': self._N_w,
        }

        if verbose:
            print()
            print("Results:")
            print(f"  Coil impedance: {result['resistance']:.6f} + j{result['reactance']:.6f} Ohm")
            print(f"  Inductance: {result['inductance']*1e6:.3f} uH")
            if self._N_w > 0:
                print(f"  Power in workpiece: {result['power_workpiece']:.3f} W")
            print()

        return result

    def compute_B_field(self, point):
        """
        Compute B field at a point from both coil and workpiece currents.

        Parameters:
            point: [x, y, z] observation point

        Returns:
            B: [Bx, By, Bz] complex magnetic field [T]
        """
        point = np.array(point)
        B = np.zeros(3, dtype=complex)

        # Quadrature points
        quad_pts = np.array([
            [0.333333333333333, 0.333333333333333],
            [0.059715871789770, 0.470142064105115],
            [0.470142064105115, 0.059715871789770],
        ])
        quad_wts = np.array([0.225, 0.132394, 0.132394])

        # Contribution from coil
        if self.I_coil is not None:
            B += self._compute_B_from_mesh(
                self.coil_mesh, self._coil_interior, self.I_coil,
                point, quad_pts, quad_wts
            )

        # Contribution from workpiece
        if self.I_workpiece is not None and len(self.I_workpiece) > 0:
            B += self._compute_B_from_mesh(
                self.workpiece_mesh, self._work_interior, self.I_workpiece,
                point, quad_pts, quad_wts
            )

        return B

    def _compute_B_from_mesh(self, mesh, interior_edges, currents, point, quad_pts, quad_wts):
        """Compute B field contribution from a mesh."""
        B = np.zeros(3, dtype=complex)

        for unknown_idx, edge_idx in enumerate(interior_edges):
            edge = mesh.edges[edge_idx]
            I_n = currents[unknown_idx]

            # Integrate over adjacent triangles
            tris = []
            if edge.tri_plus >= 0:
                tris.append((edge.tri_plus, +1, edge.local_idx_plus))
            if edge.tri_minus >= 0:
                tris.append((edge.tri_minus, -1, edge.local_idx_minus))

            for tri_idx, sign, local_idx in tris:
                tri = mesh.triangles[tri_idx]
                v = np.array([mesh.vertices[i] for i in tri.vertex_indices])
                free = mesh.vertices[tri.vertex_indices[local_idx]]
                coeff = edge.length / (2.0 * tri.area)

                for q in range(len(quad_wts)):
                    L1, L2 = quad_pts[q]
                    L3 = 1 - L1 - L2
                    r_src = L1 * v[0] + L2 * v[1] + L3 * v[2]

                    # RWG value
                    if sign > 0:
                        J = coeff * (r_src - free)
                    else:
                        J = coeff * (free - r_src)

                    J *= I_n

                    # Biot-Savart kernel
                    R_vec = point - r_src
                    R = np.linalg.norm(R_vec)

                    if R > 1e-12:
                        dB = np.cross(J, R_vec) / (R ** 3)
                        B += quad_wts[q] * tri.area * dB

        B *= mu_0 / (4 * np.pi)
        return B


class InductionHeatingSolver:
    """
    Full induction heating solver using RWG-EFIE.

    Combines:
    - Coil mesh (closed surface)
    - Workpiece mesh (open or closed surface)
    - Coupled EFIE formulation
    """

    def __init__(self):
        self.coil_mesh = None
        self.workpiece_mesh = None
        self.combined_mesh = None
        self.solver = None
        self.coupled_solver = None

        self.frequency = 1e6
        self.coil_conductivity = 5.8e7   # Copper
        self.workpiece_conductivity = 5e6  # Steel
        self.workpiece_mu_r = 1.0  # Relative permeability

    def create_spiral_coil(self, center, inner_radius, outer_radius, pitch, num_turns,
                           axis, wire_radius, num_around=8):
        """Create spiral coil mesh."""
        self.coil_mesh = RWGMesh()
        self.coil_mesh.create_spiral_coil(
            center, inner_radius, outer_radius, pitch, num_turns,
            axis, wire_radius, num_around
        )
        self.coil_mesh.conductivity = self.coil_conductivity

    def create_loop_coil(self, center, radius, normal, wire_radius, num_around=8, num_along=30):
        """Create loop coil mesh."""
        self.coil_mesh = RWGMesh()
        self.coil_mesh.create_loop_coil(
            center, radius, normal, wire_radius, num_around, num_along
        )
        self.coil_mesh.conductivity = self.coil_conductivity

    def create_plate_workpiece(self, center, Lx, Ly, normal, nx=10, ny=10):
        """Create rectangular plate workpiece mesh."""
        self.workpiece_mesh = RWGMesh()
        self.workpiece_mesh.create_rectangular_plate(center, Lx, Ly, normal, nx, ny)
        self.workpiece_mesh.conductivity = self.workpiece_conductivity

    def create_disk_workpiece(self, center, radius, normal, nr=5, ntheta=16):
        """Create circular disk workpiece mesh."""
        self.workpiece_mesh = RWGMesh()
        self.workpiece_mesh.create_circular_disk(center, radius, normal, nr, ntheta)
        self.workpiece_mesh.conductivity = self.workpiece_conductivity

    def set_frequency(self, frequency):
        """Set operating frequency [Hz]."""
        self.frequency = frequency

    def set_workpiece_permeability(self, mu_r):
        """Set workpiece relative permeability."""
        self.workpiece_mu_r = mu_r

    def solve(self, verbose=True, coupled=True):
        """
        Solve the induction heating problem.

        Parameters:
            verbose: Print progress and results
            coupled: Use coupled coil-workpiece solver if workpiece exists

        Returns:
            result: Dictionary with solution data
        """
        if self.coil_mesh is None:
            raise ValueError("Coil mesh not created")

        # Determine port edges
        boundary_edges = self.coil_mesh.get_boundary_edges()
        if boundary_edges:
            port_edges = boundary_edges[:2]
        else:
            interior_edges = [i for i, e in enumerate(self.coil_mesh.edges) if not e.is_boundary]
            port_edges = [interior_edges[0]] if interior_edges else []

        # Use coupled solver if workpiece exists
        if coupled and self.workpiece_mesh is not None:
            return self._solve_coupled(port_edges, verbose)
        else:
            return self._solve_coil_only(port_edges, verbose)

    def _solve_coil_only(self, port_edges, verbose):
        """Solve coil-only problem."""
        self.solver = RWGEFIESolver(self.coil_mesh)
        self.solver.set_frequency(self.frequency)
        self.solver.set_conductivity(self.coil_conductivity)

        if port_edges:
            self.solver.define_port(port_edges)

        self.solver.set_voltage_excitation(1.0)
        self.solver.solve()

        Z = self.solver.get_impedance()
        I = self.solver.get_current()

        result = {
            'impedance': Z,
            'impedance_real': Z.real,
            'impedance_imag': Z.imag,
            'resistance': Z.real,
            'inductance': Z.imag / (2 * np.pi * self.frequency) if self.frequency > 0 else 0,
            'current_magnitude': abs(I),
            'current_phase_deg': np.angle(I) * 180 / np.pi,
            'num_triangles_coil': len(self.coil_mesh.triangles),
            'num_edges_coil': len(self.coil_mesh.edges),
            'num_interior_edges': self.coil_mesh.num_interior_edges,
            'power_workpiece': 0.0,
            'coupled': False,
        }

        if verbose:
            print("\nRWG-EFIE Induction Heating Analysis (Coil Only)")
            print("=" * 50)
            print(f"Frequency: {self.frequency/1e3:.1f} kHz")
            print(f"Coil triangles: {result['num_triangles_coil']}")
            print(f"Interior edges (DOF): {result['num_interior_edges']}")
            print()
            print("Results:")
            print(f"  Impedance: {result['resistance']:.6f} + j{result['impedance_imag']:.6f} Ohm")
            print(f"  Inductance: {result['inductance']*1e6:.3f} uH")
            print(f"  Current: {result['current_magnitude']:.6f} A")
            print()

        return result

    def _solve_coupled(self, port_edges, verbose):
        """Solve coupled coil-workpiece problem."""
        self.coupled_solver = CoupledEFIESolver(self.coil_mesh, self.workpiece_mesh)
        self.coupled_solver.set_frequency(self.frequency)
        self.coupled_solver.set_coil_conductivity(self.coil_conductivity)
        self.coupled_solver.set_workpiece_conductivity(self.workpiece_conductivity)
        self.coupled_solver.set_workpiece_permeability(self.workpiece_mu_r)

        if port_edges:
            self.coupled_solver.define_port(port_edges)

        self.coupled_solver.set_voltage_excitation(1.0)

        if verbose:
            print("\nRWG-EFIE Induction Heating Analysis (Coupled)")
            print("=" * 50)
            print(f"Frequency: {self.frequency/1e3:.1f} kHz")
            print(f"Coil conductivity: {self.coil_conductivity:.2e} S/m")
            print(f"Workpiece conductivity: {self.workpiece_conductivity:.2e} S/m")
            print(f"Workpiece mu_r: {self.workpiece_mu_r}")

        result_coupled = self.coupled_solver.solve(verbose=verbose)

        # Get current through port for result
        I_total = 0j
        if port_edges and self.coupled_solver.I_coil is not None:
            edge_to_unknown = {e: i for i, e in enumerate(self.coupled_solver._coil_interior)}
            for edge_idx in port_edges:
                if edge_idx in edge_to_unknown:
                    I_total += self.coupled_solver.I_coil[edge_to_unknown[edge_idx]]

        result = {
            'impedance': result_coupled['impedance'],
            'impedance_real': result_coupled['resistance'],
            'impedance_imag': result_coupled['reactance'],
            'resistance': result_coupled['resistance'],
            'inductance': result_coupled['inductance'],
            'current_magnitude': abs(I_total),
            'current_phase_deg': np.angle(I_total) * 180 / np.pi,
            'num_triangles_coil': len(self.coil_mesh.triangles),
            'num_edges_coil': len(self.coil_mesh.edges),
            'num_interior_edges': self.coil_mesh.num_interior_edges,
            'num_triangles_workpiece': len(self.workpiece_mesh.triangles) if self.workpiece_mesh else 0,
            'workpiece_dof': result_coupled['workpiece_dof'],
            'power_workpiece': result_coupled['power_workpiece'],
            'coupled': True,
        }

        return result

    def compute_B_field(self, points):
        """
        Compute B field at multiple points.

        Parameters:
            points: Nx3 array of observation points

        Returns:
            B: Nx3 array of B field values [T]
        """
        # Check which solver was used
        if self.coupled_solver is not None:
            points = np.atleast_2d(points)
            B = np.zeros((len(points), 3), dtype=complex)
            for i, pt in enumerate(points):
                B[i] = self.coupled_solver.compute_B_field(pt)
            return B
        elif self.solver is not None and self.solver.solution is not None:
            points = np.atleast_2d(points)
            B = np.zeros((len(points), 3), dtype=complex)
            for i, pt in enumerate(points):
                B[i] = self.solver.compute_B_field(pt)
            return B
        else:
            raise ValueError("Must solve() first")


# Convenience functions

def create_induction_heating_model(coil_type='spiral', workpiece_type='plate', **kwargs):
    """
    Create an induction heating model.

    Parameters:
        coil_type: 'spiral' or 'loop'
        workpiece_type: 'plate', 'disk', or None
        **kwargs: Geometry parameters

    Returns:
        solver: InductionHeatingSolver object
    """
    solver = InductionHeatingSolver()

    # Default parameters
    defaults = {
        # Coil
        'coil_center': [0, 0, 0.02],
        'coil_inner_radius': 0.02,
        'coil_outer_radius': 0.05,
        'coil_pitch': 0.005,
        'coil_num_turns': 3,
        'coil_axis': [0, 0, 1],
        'wire_radius': 0.002,
        'loop_radius': 0.05,

        # Workpiece
        'workpiece_center': [0, 0, 0],
        'workpiece_Lx': 0.1,
        'workpiece_Ly': 0.1,
        'workpiece_radius': 0.05,
        'workpiece_normal': [0, 0, 1],

        # Mesh
        'num_around': 8,
        'nx': 10,
        'ny': 10,
        'nr': 5,
        'ntheta': 16,

        # Physics
        'frequency': 50e3,
        'coil_conductivity': 5.8e7,
        'workpiece_conductivity': 5e6,
    }
    defaults.update(kwargs)

    solver.frequency = defaults['frequency']
    solver.coil_conductivity = defaults['coil_conductivity']
    solver.workpiece_conductivity = defaults['workpiece_conductivity']

    # Create coil
    if coil_type == 'spiral':
        solver.create_spiral_coil(
            defaults['coil_center'],
            defaults['coil_inner_radius'],
            defaults['coil_outer_radius'],
            defaults['coil_pitch'],
            defaults['coil_num_turns'],
            defaults['coil_axis'],
            defaults['wire_radius'],
            defaults['num_around']
        )
    elif coil_type == 'loop':
        solver.create_loop_coil(
            defaults['coil_center'],
            defaults['loop_radius'],
            defaults['coil_axis'],
            defaults['wire_radius'],
            defaults['num_around']
        )

    # Create workpiece
    if workpiece_type == 'plate':
        solver.create_plate_workpiece(
            defaults['workpiece_center'],
            defaults['workpiece_Lx'],
            defaults['workpiece_Ly'],
            defaults['workpiece_normal'],
            defaults['nx'],
            defaults['ny']
        )
    elif workpiece_type == 'disk':
        solver.create_disk_workpiece(
            defaults['workpiece_center'],
            defaults['workpiece_radius'],
            defaults['workpiece_normal'],
            defaults['nr'],
            defaults['ntheta']
        )

    return solver


class NonlinearCoupledEFIESolver:
    """
    Nonlinear coupled EFIE solver with ESIM surface impedance.

    This solver handles nonlinear ferromagnetic workpieces by using
    the ESIM (Effective Surface Impedance Method) to compute field-dependent
    surface impedance Z(H).

    The nonlinear problem is solved iteratively:
        1. Solve linear EFIE with current Z estimates
        2. Compute H field at workpiece surface
        3. Update Z from ESI table: Z = Z(|H|)
        4. Repeat until convergence
    """

    def __init__(self, coil_mesh, workpiece_mesh, esi_table=None):
        """
        Initialize nonlinear coupled solver.

        Parameters:
            coil_mesh: RWGMesh for coil
            workpiece_mesh: RWGMesh for workpiece
            esi_table: ESITable object with Z(H) data (optional)
        """
        self.coil_mesh = coil_mesh
        self.workpiece_mesh = workpiece_mesh
        self.esi_table = esi_table

        self.frequency = 1e6
        self.omega = 2 * np.pi * self.frequency
        self.coil_conductivity = 5.8e7
        self.workpiece_conductivity = 5e6
        self.workpiece_mu_r = 1.0

        # For linear skin effect (when no ESI table)
        self._use_esim = esi_table is not None

        # Solution
        self.I_coil = None
        self.I_workpiece = None
        self.impedance = 0j
        self.power_workpiece = 0.0
        self.converged = False
        self.iterations = 0

        # Port
        self.port_edges = []
        self.voltage_excitation = 1.0 + 0j

        # Convergence settings
        self.max_iterations = 20
        self.tolerance = 1e-3

    def set_frequency(self, frequency):
        """Set operating frequency [Hz]."""
        self.frequency = frequency
        self.omega = 2 * np.pi * frequency

    def set_esi_table(self, esi_table):
        """Set ESI table for nonlinear material."""
        self.esi_table = esi_table
        self._use_esim = esi_table is not None

    def define_port(self, edge_indices):
        """Define port edges for coil excitation."""
        self.port_edges = list(edge_indices)

    def set_voltage_excitation(self, V):
        """Set voltage excitation [V]."""
        self.voltage_excitation = complex(V)

    def _get_surface_impedance(self, H_abs, edge_idx=None):
        """
        Get surface impedance for given |H|.

        If ESI table is available, use Z(|H|) from table.
        Otherwise, use linear skin effect model.
        """
        if self._use_esim and self.esi_table is not None:
            # Get from ESI table
            Z = self.esi_table.get_impedance(H_abs)
            return Z
        else:
            # Linear skin effect: Z_s = (1+j) / (sigma * delta)
            sigma = self.workpiece_conductivity
            mu_r = self.workpiece_mu_r
            if self.omega > 0 and sigma > 0:
                delta = np.sqrt(2.0 / (self.omega * mu_0 * mu_r * sigma))
                Z_s = (1 + 1j) / (sigma * delta)
                return Z_s
            return 0j

    def _compute_H_at_workpiece(self, I_coil, I_workpiece):
        """Compute |H| at each workpiece interior edge from coil and workpiece currents."""
        work_interior = [i for i, e in enumerate(self.workpiece_mesh.edges) if not e.is_boundary]

        H_abs = np.zeros(len(work_interior))

        # Quadrature points
        quad_pts = np.array([
            [0.333333333333333, 0.333333333333333],
        ])
        quad_wts = np.array([1.0])

        for unknown_idx, edge_idx in enumerate(work_interior):
            edge = self.workpiece_mesh.edges[edge_idx]

            # Compute H at edge midpoint
            point = edge.midpoint

            # Contribution from coil
            H_coil = self._compute_H_from_I(
                self.coil_mesh, self._coil_interior, I_coil, point
            )

            # Contribution from workpiece
            H_work = self._compute_H_from_I(
                self.workpiece_mesh, work_interior, I_workpiece, point
            )

            H_total = H_coil + H_work
            H_abs[unknown_idx] = np.abs(np.linalg.norm(H_total))

        return H_abs

    def _compute_H_from_I(self, mesh, interior_edges, currents, point):
        """Compute H field at a point from edge currents."""
        # H = B / mu_0 where B is from Biot-Savart
        H = np.zeros(3, dtype=complex)

        for unknown_idx, edge_idx in enumerate(interior_edges):
            edge = mesh.edges[edge_idx]
            I_n = currents[unknown_idx]

            # Simplified: use edge midpoint as source
            tris = []
            if edge.tri_plus >= 0:
                tris.append((edge.tri_plus, +1, edge.local_idx_plus))
            if edge.tri_minus >= 0:
                tris.append((edge.tri_minus, -1, edge.local_idx_minus))

            for tri_idx, sign, local_idx in tris:
                tri = mesh.triangles[tri_idx]

                # Use triangle center as source
                r_src = tri.center
                free = mesh.vertices[tri.vertex_indices[local_idx]]

                # RWG coefficient
                coeff = edge.length / (2.0 * tri.area)
                if sign > 0:
                    J = coeff * (r_src - free) * I_n
                else:
                    J = coeff * (free - r_src) * I_n

                # Biot-Savart for H
                R_vec = point - r_src
                R = np.linalg.norm(R_vec)

                if R > 1e-12:
                    # H = (1/4*pi) * J x R / R^3
                    dH = np.cross(J, R_vec) / (R ** 3) / (4 * np.pi)
                    H += tri.area * dH

        return H

    def solve(self, verbose=True):
        """
        Solve the nonlinear coupled system iteratively.

        Returns:
            result: Dictionary with solution data
        """
        if verbose:
            print("\nNonlinear Coupled EFIE Solver" + (" (ESIM)" if self._use_esim else " (Linear)"))
            print("=" * 50)
            print(f"Frequency: {self.frequency/1e3:.1f} kHz")

        # Get interior edges
        coil_interior = [i for i, e in enumerate(self.coil_mesh.edges) if not e.is_boundary]
        work_interior = [i for i, e in enumerate(self.workpiece_mesh.edges) if not e.is_boundary]

        self._coil_interior = coil_interior
        self._work_interior = work_interior

        N_c = len(coil_interior)
        N_w = len(work_interior)

        if verbose:
            print(f"  Coil DOF: {N_c}")
            print(f"  Workpiece DOF: {N_w}")

        if N_c == 0:
            raise ValueError("Coil has no interior edges")

        # Initialize Z for workpiece (initial guess)
        Z_work = np.full(N_w, self._get_surface_impedance(1000.0))  # Initial guess at H=1000 A/m

        # Assemble coil self-impedance (constant)
        solver_coil = RWGEFIESolver(self.coil_mesh)
        solver_coil.set_frequency(self.frequency)
        solver_coil.set_conductivity(self.coil_conductivity)
        solver_coil.assemble_matrices()
        Z_cc = solver_coil.Z_matrix

        # Iterative solve
        self.converged = False
        I_prev = None

        for iteration in range(self.max_iterations):
            if verbose:
                print(f"\n  Iteration {iteration + 1}:")

            # Assemble workpiece impedance with current Z
            Z_ww = self._assemble_workpiece_impedance(Z_work)

            # Assemble mutual impedance
            Z_cw, Z_wc = self._assemble_mutual_impedance()

            # Build full system
            if N_w > 0:
                Z_full = np.block([
                    [Z_cc, Z_cw],
                    [Z_wc, Z_ww]
                ])
            else:
                Z_full = Z_cc

            # Assemble RHS
            rhs = np.zeros(N_c + N_w, dtype=complex)
            if self.port_edges:
                edge_to_unknown = {e: i for i, e in enumerate(coil_interior)}
                V_per_edge = self.voltage_excitation / len(self.port_edges)
                for edge_idx in self.port_edges:
                    if edge_idx in edge_to_unknown:
                        unknown_idx = edge_to_unknown[edge_idx]
                        edge = self.coil_mesh.edges[edge_idx]
                        rhs[unknown_idx] = V_per_edge * edge.length

            # Solve
            try:
                solution = np.linalg.solve(Z_full, rhs)
            except np.linalg.LinAlgError:
                solution, _ = gmres(Z_full, rhs, tol=1e-6)

            I_coil = solution[:N_c]
            I_work = solution[N_c:] if N_w > 0 else np.array([])

            # Check convergence
            if I_prev is not None:
                diff = np.linalg.norm(solution - I_prev) / (np.linalg.norm(solution) + 1e-30)
                if verbose:
                    print(f"    Relative change: {diff:.6e}")
                if diff < self.tolerance:
                    self.converged = True
                    if verbose:
                        print(f"    Converged!")
                    break

            I_prev = solution.copy()

            # Update Z for workpiece based on H field (if ESIM)
            if self._use_esim and N_w > 0:
                H_abs = self._compute_H_at_workpiece(I_coil, I_work)
                for i in range(N_w):
                    Z_work[i] = self._get_surface_impedance(H_abs[i])
                if verbose:
                    print(f"    H_max at workpiece: {np.max(H_abs):.1f} A/m")

        self.iterations = iteration + 1
        self.I_coil = I_coil
        self.I_workpiece = I_work

        # Compute port impedance
        if self.port_edges:
            edge_to_unknown = {e: i for i, e in enumerate(coil_interior)}
            I_total = sum(I_coil[edge_to_unknown[e]] for e in self.port_edges if e in edge_to_unknown)
            if abs(I_total) > 1e-30:
                self.impedance = self.voltage_excitation / I_total
            else:
                self.impedance = np.inf + 0j

        # Compute power in workpiece
        if N_w > 0:
            self.power_workpiece = 0.5 * np.real(np.conj(I_work) @ Z_ww @ I_work)
        else:
            self.power_workpiece = 0.0

        result = {
            'impedance': self.impedance,
            'resistance': self.impedance.real,
            'reactance': self.impedance.imag,
            'inductance': self.impedance.imag / (2 * np.pi * self.frequency) if self.frequency > 0 else 0,
            'power_workpiece': self.power_workpiece,
            'coil_dof': N_c,
            'workpiece_dof': N_w,
            'converged': self.converged,
            'iterations': self.iterations,
        }

        if verbose:
            print()
            print("Results:")
            print(f"  Impedance: {result['resistance']:.6f} + j{result['reactance']:.6f} Ohm")
            print(f"  Inductance: {result['inductance']*1e6:.3f} uH")
            print(f"  Power in workpiece: {result['power_workpiece']:.3f} W")
            print(f"  Converged: {result['converged']} in {result['iterations']} iterations")
            print()

        return result

    def _assemble_workpiece_impedance(self, Z_surface):
        """Assemble workpiece self-impedance with given surface impedance."""
        work_interior = [i for i, e in enumerate(self.workpiece_mesh.edges) if not e.is_boundary]
        N_w = len(work_interior)

        if N_w == 0:
            return np.zeros((0, 0), dtype=complex)

        # Use CoupledEFIESolver's method to get L matrix
        solver_work = RWGEFIESolver(self.workpiece_mesh)
        solver_work.set_frequency(self.frequency)
        solver_work.set_conductivity(self.workpiece_conductivity)
        solver_work.assemble_matrices()

        # Replace R matrix with ESIM surface impedance
        j_omega = 1j * self.omega
        Z_ww = j_omega * solver_work.L_matrix

        # Add surface impedance contribution
        for i, (edge_idx, edge) in enumerate([(idx, self.workpiece_mesh.edges[idx]) for idx in work_interior]):
            # Diagonal term: Z_s * effective length / width
            Z_ww[i, i] += Z_surface[i]

        return Z_ww

    def _assemble_mutual_impedance(self):
        """Assemble coil-workpiece mutual impedance."""
        coil_interior = self._coil_interior
        work_interior = self._work_interior

        N_c = len(coil_interior)
        N_w = len(work_interior)

        if N_w == 0:
            return np.zeros((N_c, 0), dtype=complex), np.zeros((0, N_c), dtype=complex)

        Z_cw = np.zeros((N_c, N_w), dtype=complex)
        Z_wc = np.zeros((N_w, N_c), dtype=complex)

        j_omega = 1j * self.omega

        # Quadrature points
        quad_pts = np.array([
            [0.333333333333333, 0.333333333333333],
            [0.059715871789770, 0.470142064105115],
            [0.470142064105115, 0.059715871789770],
            [0.470142064105115, 0.470142064105115],
        ])
        quad_wts = np.array([0.225, 0.132394, 0.132394, 0.132394])

        for i, edge_c_idx in enumerate(coil_interior):
            for j, edge_w_idx in enumerate(work_interior):
                L_cw = self._compute_mutual_L(edge_c_idx, edge_w_idx, quad_pts, quad_wts)
                Z_cw[i, j] = j_omega * L_cw
                Z_wc[j, i] = j_omega * L_cw

        return Z_cw, Z_wc

    def _compute_mutual_L(self, edge_c_idx, edge_w_idx, quad_pts, quad_wts):
        """Compute mutual inductance element."""
        edge_c = self.coil_mesh.edges[edge_c_idx]
        edge_w = self.workpiece_mesh.edges[edge_w_idx]

        if edge_c.is_boundary or edge_w.is_boundary:
            return 0.0

        L_mn = 0.0

        tris_c = []
        if edge_c.tri_plus >= 0:
            tris_c.append((edge_c.tri_plus, +1, edge_c.local_idx_plus))
        if edge_c.tri_minus >= 0:
            tris_c.append((edge_c.tri_minus, -1, edge_c.local_idx_minus))

        tris_w = []
        if edge_w.tri_plus >= 0:
            tris_w.append((edge_w.tri_plus, +1, edge_w.local_idx_plus))
        if edge_w.tri_minus >= 0:
            tris_w.append((edge_w.tri_minus, -1, edge_w.local_idx_minus))

        for tri_c_idx, sign_c, local_c in tris_c:
            tri_c = self.coil_mesh.triangles[tri_c_idx]
            v_c = np.array([self.coil_mesh.vertices[k] for k in tri_c.vertex_indices])
            free_c = self.coil_mesh.vertices[tri_c.vertex_indices[local_c]]
            coeff_c = sign_c * edge_c.length / (2.0 * tri_c.area)

            for tri_w_idx, sign_w, local_w in tris_w:
                tri_w = self.workpiece_mesh.triangles[tri_w_idx]
                v_w = np.array([self.workpiece_mesh.vertices[k] for k in tri_w.vertex_indices])
                free_w = self.workpiece_mesh.vertices[tri_w.vertex_indices[local_w]]
                coeff_w = sign_w * edge_w.length / (2.0 * tri_w.area)

                for q1 in range(len(quad_wts)):
                    L1_1, L2_1 = quad_pts[q1]
                    L3_1 = 1 - L1_1 - L2_1
                    r1 = L1_1 * v_c[0] + L2_1 * v_c[1] + L3_1 * v_c[2]

                    if sign_c > 0:
                        f1 = coeff_c * (r1 - free_c)
                    else:
                        f1 = -coeff_c * (free_c - r1)

                    for q2 in range(len(quad_wts)):
                        L1_2, L2_2 = quad_pts[q2]
                        L3_2 = 1 - L1_2 - L2_2
                        r2 = L1_2 * v_w[0] + L2_2 * v_w[1] + L3_2 * v_w[2]

                        if sign_w > 0:
                            f2 = coeff_w * (r2 - free_w)
                        else:
                            f2 = -coeff_w * (free_w - r2)

                        R = np.linalg.norm(r1 - r2)
                        if R > 1e-12:
                            G = 1.0 / R
                        else:
                            G = 0.0

                        L_mn += quad_wts[q1] * quad_wts[q2] * tri_c.area * tri_w.area * np.dot(f1, f2) * G

        L_mn *= mu_0 / (4 * np.pi)
        return L_mn


class pFFTAccelerator:
    """
    Pre-corrected FFT (pFFT) accelerator for fast matrix-vector products.

    The pFFT method decomposes the matrix-vector product into:
        Z * I = Z_near * I + Z_far * I

    The far-field part uses FFT-based convolution on a uniform grid:
        Z_far * I = P^T * G_fft * P * I

    where:
        P = projection operator (triangles -> grid)
        G_fft = Green's function convolution via FFT
        P^T = interpolation operator (grid -> triangles)

    Complexity: O(N log N) vs O(N^2) for dense matrix-vector product

    Reference:
        J.R. Phillips, J.K. White, "A precorrected-FFT method for electrostatic
        analysis of complicated 3-D structures," IEEE TCAD, 1997.
    """

    def __init__(self, mesh, frequency):
        """
        Initialize pFFT accelerator.

        Parameters:
            mesh: RWGMesh object
            frequency: Operating frequency [Hz]
        """
        self.mesh = mesh
        self.frequency = frequency
        self.omega = 2 * np.pi * frequency

        # Grid parameters
        self.grid_size = None  # [Nx, Ny, Nz]
        self.grid_h = None     # Grid spacing
        self.grid_origin = None  # Grid origin

        # Near-field parameters
        self.near_radius = None  # Near-field cutoff radius

        # Precomputed data
        self.G_fft = None       # FFT of Green's function
        self.P_matrices = None  # Projection matrices for each basis function
        self.near_indices = None  # Near-field interaction indices
        self.Z_near = None      # Near-field correction matrix (sparse)

        # Edge mapping
        self.interior_edges = None
        self.edge_to_unknown = None

        self._initialized = False

    def setup(self, grid_order=3, near_layers=2):
        """
        Set up the pFFT accelerator.

        Parameters:
            grid_order: Projection order (number of grid cells for each basis function)
            near_layers: Number of grid layers for near-field
        """
        print("  Setting up pFFT accelerator...")

        # Get interior edges
        self.interior_edges = [i for i, e in enumerate(self.mesh.edges) if not e.is_boundary]
        self.edge_to_unknown = {e: i for i, e in enumerate(self.interior_edges)}
        N = len(self.interior_edges)

        if N == 0:
            warnings.warn("No interior edges found")
            return

        # Compute bounding box
        all_vertices = np.array(self.mesh.vertices)
        bbox_min = np.min(all_vertices, axis=0)
        bbox_max = np.max(all_vertices, axis=0)

        # Estimate characteristic mesh size
        mesh_h = self._estimate_mesh_size()
        print(f"    Characteristic mesh size: {mesh_h*1000:.3f} mm")

        # Set grid spacing (should be comparable to mesh size)
        self.grid_h = mesh_h * grid_order
        self.near_radius = near_layers * self.grid_h

        # Determine grid dimensions
        extent = bbox_max - bbox_min
        # Pad the grid for FFT convolution
        pad = 2 * self.near_radius
        self.grid_origin = bbox_min - pad

        self.grid_size = np.ceil((extent + 2 * pad) / self.grid_h).astype(int)
        # Make grid sizes suitable for FFT (powers of 2 or smooth numbers)
        self.grid_size = self._next_fft_size(self.grid_size)

        Nx, Ny, Nz = self.grid_size
        print(f"    Grid size: {Nx} x {Ny} x {Nz}")
        print(f"    Grid spacing: {self.grid_h*1000:.3f} mm")
        print(f"    Near-field radius: {self.near_radius*1000:.3f} mm")

        # Precompute FFT of Green's function
        print("    Computing FFT of Green's function...")
        self._compute_green_fft()

        # Build projection matrices
        print("    Building projection matrices...")
        self._build_projection_matrices()

        # Identify near-field pairs and compute corrections
        print("    Computing near-field corrections...")
        self._compute_near_field()

        self._initialized = True
        print("    pFFT setup complete.")

    def _estimate_mesh_size(self):
        """Estimate characteristic mesh element size."""
        areas = [tri.area for tri in self.mesh.triangles]
        avg_area = np.mean(areas)
        return np.sqrt(avg_area)  # Characteristic length

    def _next_fft_size(self, sizes):
        """Find next suitable FFT size (power of 2 or smooth number)."""
        result = []
        for n in sizes:
            # Find next power of 2
            p2 = 2 ** int(np.ceil(np.log2(n)))
            # Or use smooth number (product of small primes)
            result.append(max(p2, 8))  # Minimum 8
        return np.array(result)

    def _compute_green_fft(self):
        """Compute FFT of the free-space Green's function."""
        Nx, Ny, Nz = self.grid_size
        h = self.grid_h

        # Green's function kernel: G(r) = 1 / (4*pi*|r|)
        # Build on a doubled grid for proper FFT convolution (zero-padding)
        G = np.zeros((2*Nx, 2*Ny, 2*Nz))

        # Grid coordinates (centered at origin for convolution)
        for i in range(2*Nx):
            ix = i if i < Nx else i - 2*Nx
            for j in range(2*Ny):
                iy = j if j < Ny else j - 2*Ny
                for k in range(2*Nz):
                    iz = k if k < Nz else k - 2*Nz

                    r = np.sqrt((ix*h)**2 + (iy*h)**2 + (iz*h)**2)

                    if r > 1e-15:
                        G[i, j, k] = 1.0 / (4 * np.pi * r)
                    else:
                        # Regularization for self-term
                        G[i, j, k] = 1.0 / (4 * np.pi * h * 0.5)

        # FFT of Green's function (stored for reuse)
        self.G_fft = np.fft.fftn(G)

    def _build_projection_matrices(self):
        """Build projection matrices for each RWG basis function."""
        Nx, Ny, Nz = self.grid_size
        h = self.grid_h
        N = len(self.interior_edges)

        # For each interior edge, compute projection weights to nearby grid points
        # Using trilinear interpolation (8 grid points per source point)

        # We use centroid of each adjacent triangle as representative source point
        self.P_matrices = []  # List of (grid_indices, weights) for each edge

        for unknown_idx, edge_idx in enumerate(self.interior_edges):
            edge = self.mesh.edges[edge_idx]

            # Collect contributions from adjacent triangles
            proj_data = {'indices': [], 'weights': []}

            tris = []
            if edge.tri_plus >= 0:
                tris.append((edge.tri_plus, +1))
            if edge.tri_minus >= 0:
                tris.append((edge.tri_minus, -1))

            for tri_idx, sign in tris:
                tri = self.mesh.triangles[tri_idx]

                # Use triangle centroid
                centroid = tri.center

                # Find grid cell
                local = (centroid - self.grid_origin) / h
                i0, j0, k0 = int(local[0]), int(local[1]), int(local[2])

                # Trilinear interpolation weights
                fx, fy, fz = local - np.array([i0, j0, k0])

                # 8 corner contributions
                for di in [0, 1]:
                    for dj in [0, 1]:
                        for dk in [0, 1]:
                            ii = i0 + di
                            jj = j0 + dj
                            kk = k0 + dk

                            if 0 <= ii < Nx and 0 <= jj < Ny and 0 <= kk < Nz:
                                # Trilinear weight
                                wx = (1 - fx) if di == 0 else fx
                                wy = (1 - fy) if dj == 0 else fy
                                wz = (1 - fz) if dk == 0 else fz
                                w = wx * wy * wz

                                # Scale by triangle contribution and basis function
                                coeff = edge.length / (2.0 * tri.area)
                                w *= coeff * tri.area * sign

                                if abs(w) > 1e-15:
                                    proj_data['indices'].append((ii, jj, kk))
                                    proj_data['weights'].append(w)

            self.P_matrices.append(proj_data)

    def _compute_near_field(self):
        """Compute near-field interaction indices and correction matrix."""
        N = len(self.interior_edges)

        # Build spatial hash for fast neighbor lookup
        edge_centers = []
        for edge_idx in self.interior_edges:
            edge = self.mesh.edges[edge_idx]
            edge_centers.append(edge.midpoint)
        edge_centers = np.array(edge_centers)

        # Find near-field pairs (O(N) with spatial hashing)
        near_pairs = []

        # Simple O(N^2) check for now (can be optimized with k-d tree)
        for m in range(N):
            for n in range(m, N):
                dist = np.linalg.norm(edge_centers[m] - edge_centers[n])
                if dist < self.near_radius:
                    near_pairs.append((m, n))

        self.near_indices = near_pairs
        print(f"    Near-field pairs: {len(near_pairs)} (out of {N*(N+1)//2} total)")

        # Compute near-field correction: Z_near_correct = Z_near_direct - Z_near_fft
        # For efficiency, store as sparse matrix

        # Get edges for computation
        interior_edge_objects = [self.mesh.edges[idx] for idx in self.interior_edges]

        # Compute near-field matrix elements
        self.Z_near = np.zeros((N, N), dtype=complex)

        solver_temp = RWGEFIESolver(self.mesh)
        solver_temp.set_frequency(self.frequency)
        solver_temp.set_conductivity(self.mesh.conductivity)

        for m, n in near_pairs:
            edge_m = interior_edge_objects[m]
            edge_n = interior_edge_objects[n]

            # Direct computation
            L_direct = solver_temp._compute_L_element(edge_m, edge_n)
            R_direct = solver_temp._compute_R_element(edge_m, edge_n)
            Z_direct = 1j * self.omega * L_direct + R_direct

            # FFT-based estimate for this pair (simplified)
            # For near-field, FFT gives an approximation which we correct
            Z_fft = self._compute_fft_element(m, n)

            # Correction
            correction = Z_direct - Z_fft
            self.Z_near[m, n] = correction
            if m != n:
                self.Z_near[n, m] = correction

    def _compute_fft_element(self, m, n):
        """Compute single matrix element using FFT-based method."""
        # This is the approximate value from FFT convolution
        # For proper correction, we compute what the FFT would give

        Nx, Ny, Nz = self.grid_size

        # Project edge m to grid
        proj_m = self.P_matrices[m]
        proj_n = self.P_matrices[n]

        # Simplified: use centroid-based estimate
        # In full implementation, this would use the actual FFT convolution

        # Sum over grid interactions
        Z_approx = 0j

        for (im, jm, km), wm in zip(proj_m['indices'], proj_m['weights']):
            for (in_, jn, kn), wn in zip(proj_n['indices'], proj_n['weights']):
                # Distance in grid indices
                di = im - in_
                dj = jm - jn
                dk = km - kn

                # Green's function from grid
                r = self.grid_h * np.sqrt(di**2 + dj**2 + dk**2)
                if r > 1e-15:
                    G = 1.0 / (4 * np.pi * r)
                else:
                    G = 1.0 / (4 * np.pi * self.grid_h * 0.5)

                Z_approx += wm * wn * G

        # Scale by mu_0 for inductance
        Z_approx *= mu_0 * 1j * self.omega

        return Z_approx

    def matvec(self, I):
        """
        Compute matrix-vector product Z * I using pFFT.

        Parameters:
            I: Input current vector (length N)

        Returns:
            V: Output voltage vector (length N)
        """
        if not self._initialized:
            raise RuntimeError("pFFT not initialized. Call setup() first.")

        N = len(self.interior_edges)
        Nx, Ny, Nz = self.grid_size

        # Step 1: Project currents to grid
        grid_q = np.zeros((Nx, Ny, Nz), dtype=complex)

        for m in range(N):
            proj = self.P_matrices[m]
            for (i, j, k), w in zip(proj['indices'], proj['weights']):
                grid_q[i, j, k] += w * I[m]

        # Step 2: Zero-pad for FFT convolution
        grid_q_pad = np.zeros((2*Nx, 2*Ny, 2*Nz), dtype=complex)
        grid_q_pad[:Nx, :Ny, :Nz] = grid_q

        # Step 3: FFT convolution
        Q_fft = np.fft.fftn(grid_q_pad)
        V_fft = self.G_fft * Q_fft
        grid_v = np.fft.ifftn(V_fft)[:Nx, :Ny, :Nz]

        # Step 4: Interpolate from grid to edges
        V_far = np.zeros(N, dtype=complex)

        for m in range(N):
            proj = self.P_matrices[m]
            for (i, j, k), w in zip(proj['indices'], proj['weights']):
                V_far[m] += w * grid_v[i, j, k]

        # Scale by j*omega*mu_0 for inductance term
        V_far *= 1j * self.omega * mu_0

        # Step 5: Add near-field correction
        V_near = self.Z_near @ I

        return V_far + V_near


class pFFTEFIESolver:
    """
    EFIE Solver with pFFT acceleration.

    Uses pFFT for O(N log N) matrix-vector products within an iterative solver.
    Suitable for large-scale problems.
    """

    def __init__(self, mesh):
        """
        Initialize pFFT-accelerated solver.

        Parameters:
            mesh: RWGMesh object
        """
        self.mesh = mesh
        self.frequency = 1e6
        self.omega = 2 * np.pi * self.frequency
        self.conductivity = 5.8e7

        self.pfft = None
        self.solution = None
        self.impedance = 0j
        self.total_current = 0j

        self.port_edges = []
        self.voltage_excitation = 1.0 + 0j

        # Iterative solver settings
        self.tol = 1e-6
        self.maxiter = 500

    def set_frequency(self, frequency):
        """Set operating frequency [Hz]."""
        self.frequency = frequency
        self.omega = 2 * np.pi * frequency

    def set_conductivity(self, sigma):
        """Set conductivity [S/m]."""
        self.conductivity = sigma
        self.mesh.conductivity = sigma

    def define_port(self, edge_indices):
        """Define port edges."""
        self.port_edges = list(edge_indices)

    def set_voltage_excitation(self, V):
        """Set voltage excitation [V]."""
        self.voltage_excitation = complex(V)

    def setup(self, grid_order=3, near_layers=2):
        """
        Set up pFFT accelerator.

        Parameters:
            grid_order: Grid resolution relative to mesh
            near_layers: Near-field cutoff in grid layers
        """
        self.pfft = pFFTAccelerator(self.mesh, self.frequency)
        self.pfft.setup(grid_order, near_layers)

    def _get_interior_edge_mapping(self):
        """Get mapping from interior edges to unknown indices."""
        edge_to_unknown = {}
        unknown_idx = 0
        for i, edge in enumerate(self.mesh.edges):
            if not edge.is_boundary:
                edge_to_unknown[i] = unknown_idx
                unknown_idx += 1
        return edge_to_unknown

    def solve(self, verbose=True):
        """
        Solve using pFFT-accelerated iterative solver.

        Uses GMRES or BiCGSTAB with pFFT for matrix-vector products.

        Parameters:
            verbose: Print progress

        Returns:
            result: Solution dictionary
        """
        if self.pfft is None:
            if verbose:
                print("Setting up pFFT accelerator...")
            self.setup()

        edge_to_unknown = self._get_interior_edge_mapping()
        N = self.mesh.num_interior_edges

        if N == 0:
            warnings.warn("No interior edges found")
            return {}

        if verbose:
            print(f"\npFFT-Accelerated EFIE Solver")
            print("=" * 50)
            print(f"Frequency: {self.frequency/1e3:.1f} kHz")
            print(f"Total DOF: {N}")

        # Assemble RHS
        rhs = np.zeros(N, dtype=complex)

        if self.port_edges:
            V_per_edge = self.voltage_excitation / len(self.port_edges)
            for edge_idx in self.port_edges:
                if edge_idx in edge_to_unknown:
                    unknown_idx = edge_to_unknown[edge_idx]
                    edge = self.mesh.edges[edge_idx]
                    rhs[unknown_idx] = V_per_edge * edge.length

        # Define linear operator for iterative solver
        from scipy.sparse.linalg import LinearOperator

        def matvec(x):
            return self.pfft.matvec(x)

        A = LinearOperator((N, N), matvec=matvec, dtype=complex)

        # Solve with GMRES
        if verbose:
            print("  Solving with GMRES...")

        # Callback for progress
        iteration_count = [0]
        def callback(xk):
            iteration_count[0] += 1
            if verbose and iteration_count[0] % 50 == 0:
                print(f"    Iteration {iteration_count[0]}")

        self.solution, info = gmres(A, rhs, tol=self.tol, maxiter=self.maxiter,
                                     callback=callback, callback_type='x')

        if info != 0:
            if verbose:
                print(f"  Warning: GMRES did not converge (info={info})")
        elif verbose:
            print(f"  Converged in {iteration_count[0]} iterations")

        # Compute port impedance
        if self.port_edges:
            I_total = 0j
            for edge_idx in self.port_edges:
                if edge_idx in edge_to_unknown:
                    unknown_idx = edge_to_unknown[edge_idx]
                    I_total += self.solution[unknown_idx]

            self.total_current = I_total
            if abs(I_total) > 1e-30:
                self.impedance = self.voltage_excitation / I_total
            else:
                self.impedance = np.inf + 0j

        result = {
            'impedance': self.impedance,
            'resistance': self.impedance.real,
            'reactance': self.impedance.imag,
            'inductance': self.impedance.imag / (2 * np.pi * self.frequency) if self.frequency > 0 else 0,
            'current': self.total_current,
            'iterations': iteration_count[0],
            'dof': N,
        }

        if verbose:
            print()
            print("Results:")
            print(f"  Impedance: {result['resistance']:.6f} + j{result['reactance']:.6f} Ohm")
            print(f"  Inductance: {result['inductance']*1e6:.3f} uH")
            print()

        return result

    def compute_B_field(self, point):
        """Compute B field at a point."""
        if self.solution is None:
            return np.zeros(3, dtype=complex)

        # Use same method as RWGEFIESolver
        solver_temp = RWGEFIESolver(self.mesh)
        solver_temp.solution = self.solution
        return solver_temp.compute_B_field(point)


class HybridPEECRWGSolver:
    """
    Hybrid PEEC + RWG-EFIE solver.

    Combines:
    - PEEC (Partial Element Equivalent Circuit) for coils
    - RWG-EFIE for workpieces

    This approach leverages the strengths of both methods:
    - PEEC: Efficient for wire/coil structures with simple cross-sections
    - RWG-EFIE: Accurate for arbitrary 3D surfaces (workpieces)

    The coupling between PEEC segments and RWG triangles is computed using
    mutual inductance integrals.

    System matrix:
        [Z_peec   Z_pw  ] [I_p]   [V_p]
        [Z_wp    Z_rwg  ] [I_w] = [0  ]

    where:
        Z_peec = PEEC impedance matrix (wire segments)
        Z_rwg  = RWG-EFIE impedance matrix (surface triangles)
        Z_pw, Z_wp = coupling between PEEC and RWG
    """

    def __init__(self):
        """Initialize hybrid solver."""
        # Coil (PEEC) model
        self.coil_segments = []  # List of wire segments
        self.coil_conductivity = 5.8e7
        self.coil_wire_radius = 0.002

        # Workpiece (RWG) model
        self.workpiece_mesh = None
        self.workpiece_conductivity = 5e6
        self.workpiece_mu_r = 1.0

        # Frequency
        self.frequency = 1e6
        self.omega = 2 * np.pi * self.frequency

        # Matrices
        self.Z_peec = None
        self.Z_rwg = None
        self.Z_pw = None
        self.Z_wp = None

        # Solution
        self.I_peec = None
        self.I_rwg = None
        self.impedance = 0j
        self.power_workpiece = 0.0

        # Port
        self.port_segments = []
        self.voltage_excitation = 1.0 + 0j

    def set_frequency(self, frequency):
        """Set operating frequency [Hz]."""
        self.frequency = frequency
        self.omega = 2 * np.pi * frequency

    def add_coil_segment(self, start, end, radius=None):
        """
        Add a PEEC wire segment to the coil.

        Parameters:
            start: [x, y, z] start point
            end: [x, y, z] end point
            radius: Wire radius (optional, uses default)
        """
        r = radius if radius is not None else self.coil_wire_radius
        self.coil_segments.append({
            'start': np.array(start),
            'end': np.array(end),
            'radius': r,
            'length': np.linalg.norm(np.array(end) - np.array(start)),
        })

    def create_circular_coil_peec(self, center, radius, normal, n_segments=36):
        """
        Create circular coil using PEEC segments.

        Parameters:
            center: Coil center
            radius: Coil radius
            normal: Coil normal direction
            n_segments: Number of segments
        """
        center = np.array(center)
        normal = np.array(normal)
        normal = normal / np.linalg.norm(normal)

        # Create local coordinate system
        if abs(normal[2]) < 0.9:
            u = np.cross(normal, [0, 0, 1])
        else:
            u = np.cross(normal, [1, 0, 0])
        u = u / np.linalg.norm(u)
        v = np.cross(normal, u)

        # Generate segment endpoints
        dtheta = 2 * np.pi / n_segments
        for i in range(n_segments):
            theta1 = i * dtheta
            theta2 = (i + 1) * dtheta

            p1 = center + radius * (np.cos(theta1) * u + np.sin(theta1) * v)
            p2 = center + radius * (np.cos(theta2) * u + np.sin(theta2) * v)

            self.add_coil_segment(p1, p2)

    def create_spiral_coil_peec(self, center, inner_radius, outer_radius,
                                 pitch, num_turns, axis, n_segments_per_turn=36):
        """
        Create spiral coil using PEEC segments.

        Parameters:
            center: Coil center
            inner_radius, outer_radius: Spiral radii
            pitch: Axial advance per turn
            num_turns: Number of turns
            axis: Coil axis
            n_segments_per_turn: Segments per turn
        """
        center = np.array(center)
        axis = np.array(axis)
        axis = axis / np.linalg.norm(axis)

        # Create local coordinate system
        if abs(axis[2]) < 0.9:
            u = np.cross(axis, [0, 0, 1])
        else:
            u = np.cross(axis, [1, 0, 0])
        u = u / np.linalg.norm(u)
        v = np.cross(axis, u)

        # Total segments
        n_total = int(num_turns * n_segments_per_turn)
        total_angle = 2 * np.pi * num_turns

        for i in range(n_total):
            t1 = i / n_total
            t2 = (i + 1) / n_total

            theta1 = t1 * total_angle
            theta2 = t2 * total_angle
            R1 = inner_radius + t1 * (outer_radius - inner_radius)
            R2 = inner_radius + t2 * (outer_radius - inner_radius)
            z1 = t1 * num_turns * pitch
            z2 = t2 * num_turns * pitch

            p1 = center + R1 * (np.cos(theta1) * u + np.sin(theta1) * v) + z1 * axis
            p2 = center + R2 * (np.cos(theta2) * u + np.sin(theta2) * v) + z2 * axis

            self.add_coil_segment(p1, p2)

    def set_workpiece_mesh(self, mesh):
        """Set workpiece mesh (RWGMesh)."""
        self.workpiece_mesh = mesh

    def create_workpiece_plate(self, center, Lx, Ly, normal, nx=10, ny=10):
        """Create rectangular plate workpiece."""
        self.workpiece_mesh = RWGMesh()
        self.workpiece_mesh.create_rectangular_plate(center, Lx, Ly, normal, nx, ny)

    def create_workpiece_disk(self, center, radius, normal, nr=5, ntheta=16):
        """Create circular disk workpiece."""
        self.workpiece_mesh = RWGMesh()
        self.workpiece_mesh.create_circular_disk(center, radius, normal, nr, ntheta)

    def define_port(self, segment_indices):
        """Define port segments for coil excitation."""
        self.port_segments = list(segment_indices)

    def _compute_segment_self_inductance(self, segment):
        """
        Compute self-inductance of a wire segment.

        L_self = (mu_0 / 2*pi) * L * [ln(2*L/r) - 1]

        where L is segment length and r is wire radius.
        """
        L = segment['length']
        r = segment['radius']

        if L < 1e-15:
            return 0.0

        return (mu_0 / (2 * np.pi)) * L * (np.log(2 * L / r) - 1)

    def _compute_segment_mutual_inductance(self, seg1, seg2):
        """
        Compute mutual inductance between two wire segments.

        Uses Neumann formula:
            M = (mu_0 / 4*pi) * integral{ dl1 . dl2 / |r1 - r2| }
        """
        # Simple midpoint approximation for well-separated segments
        c1 = (seg1['start'] + seg1['end']) / 2
        c2 = (seg2['start'] + seg2['end']) / 2

        R = np.linalg.norm(c1 - c2)

        if R < 1e-12:
            return 0.0

        # Direction vectors
        d1 = seg1['end'] - seg1['start']
        d2 = seg2['end'] - seg2['start']

        L1 = seg1['length']
        L2 = seg2['length']

        if L1 < 1e-15 or L2 < 1e-15:
            return 0.0

        # Neumann formula (midpoint approximation)
        cos_angle = np.dot(d1, d2) / (L1 * L2)
        M = (mu_0 / (4 * np.pi)) * L1 * L2 * cos_angle / R

        return M

    def _compute_segment_resistance(self, segment):
        """
        Compute AC resistance of a wire segment with skin effect.

        R_ac = R_dc * (r / delta) for r >> delta
        """
        L = segment['length']
        r = segment['radius']

        # Skin depth
        delta = np.sqrt(2.0 / (self.omega * mu_0 * self.coil_conductivity))

        # DC resistance
        R_dc = L / (np.pi * r**2 * self.coil_conductivity)

        # AC factor (simplified for thick wire)
        if r > delta:
            ac_factor = r / (2 * delta)
        else:
            ac_factor = 1.0

        return R_dc * ac_factor

    def _compute_mutual_peec_rwg(self, segment_idx, edge_idx):
        """
        Compute mutual inductance between PEEC segment and RWG edge.

        M = (mu_0 / 4*pi) * integral{ J_seg . f_rwg / |r1 - r2| } dA dL
        """
        segment = self.coil_segments[segment_idx]
        edge = self.workpiece_mesh.edges[edge_idx]

        if edge.is_boundary:
            return 0.0

        # Quadrature points for triangle
        quad_pts = np.array([
            [0.333333333333333, 0.333333333333333],
            [0.059715871789770, 0.470142064105115],
            [0.470142064105115, 0.059715871789770],
        ])
        quad_wts = np.array([0.225, 0.132394, 0.132394])

        # Segment points (3-point quadrature along segment)
        seg_pts = np.array([0.1127, 0.5, 0.8873])
        seg_wts = np.array([0.278, 0.444, 0.278])

        M = 0.0

        # Segment direction (current direction)
        d_seg = segment['end'] - segment['start']
        L_seg = segment['length']
        if L_seg < 1e-15:
            return 0.0
        J_dir = d_seg / L_seg

        # RWG triangles
        tris = []
        if edge.tri_plus >= 0:
            tris.append((edge.tri_plus, +1, edge.local_idx_plus))
        if edge.tri_minus >= 0:
            tris.append((edge.tri_minus, -1, edge.local_idx_minus))

        for tri_idx, sign, local_idx in tris:
            tri = self.workpiece_mesh.triangles[tri_idx]
            v = np.array([self.workpiece_mesh.vertices[i] for i in tri.vertex_indices])
            free = self.workpiece_mesh.vertices[tri.vertex_indices[local_idx]]
            coeff = edge.length / (2.0 * tri.area)

            for q_seg in range(len(seg_wts)):
                # Point on segment
                t = seg_pts[q_seg]
                r_seg = segment['start'] + t * d_seg

                for q_tri in range(len(quad_wts)):
                    L1, L2 = quad_pts[q_tri]
                    L3 = 1 - L1 - L2
                    r_tri = L1 * v[0] + L2 * v[1] + L3 * v[2]

                    # RWG basis value
                    if sign > 0:
                        f_rwg = coeff * (r_tri - free)
                    else:
                        f_rwg = coeff * (free - r_tri)

                    # Distance
                    R = np.linalg.norm(r_seg - r_tri)
                    if R > 1e-12:
                        # Integrand: J_seg . f_rwg / R
                        integrand = np.dot(J_dir, f_rwg) / R
                        w = seg_wts[q_seg] * quad_wts[q_tri] * L_seg * tri.area
                        M += w * integrand

        M *= mu_0 / (4 * np.pi)
        return M

    def assemble_system(self, verbose=True):
        """Assemble the hybrid PEEC-RWG system."""
        N_p = len(self.coil_segments)

        work_interior = []
        if self.workpiece_mesh is not None:
            work_interior = [i for i, e in enumerate(self.workpiece_mesh.edges) if not e.is_boundary]
        N_w = len(work_interior)

        if verbose:
            print(f"  PEEC segments: {N_p}")
            print(f"  RWG interior edges: {N_w}")
            print(f"  Total DOF: {N_p + N_w}")

        j_omega = 1j * self.omega

        # PEEC matrix
        if verbose:
            print("  Assembling PEEC matrix...")

        self.Z_peec = np.zeros((N_p, N_p), dtype=complex)

        for i in range(N_p):
            seg_i = self.coil_segments[i]
            # Self terms
            L_self = self._compute_segment_self_inductance(seg_i)
            R_self = self._compute_segment_resistance(seg_i)
            self.Z_peec[i, i] = j_omega * L_self + R_self

            # Mutual terms
            for j in range(i + 1, N_p):
                seg_j = self.coil_segments[j]
                M = self._compute_segment_mutual_inductance(seg_i, seg_j)
                self.Z_peec[i, j] = j_omega * M
                self.Z_peec[j, i] = j_omega * M

        # RWG matrix
        if N_w > 0:
            if verbose:
                print("  Assembling RWG matrix...")

            solver_rwg = RWGEFIESolver(self.workpiece_mesh)
            solver_rwg.set_frequency(self.frequency)
            solver_rwg.set_conductivity(self.workpiece_conductivity)
            solver_rwg.assemble_matrices()
            self.Z_rwg = solver_rwg.Z_matrix

            # Coupling matrices
            if verbose:
                print("  Assembling coupling matrices...")

            self.Z_pw = np.zeros((N_p, N_w), dtype=complex)
            self.Z_wp = np.zeros((N_w, N_p), dtype=complex)

            for i in range(N_p):
                for j, edge_idx in enumerate(work_interior):
                    M = self._compute_mutual_peec_rwg(i, edge_idx)
                    self.Z_pw[i, j] = j_omega * M
                    self.Z_wp[j, i] = j_omega * M

            # Build full system
            self.Z_full = np.block([
                [self.Z_peec, self.Z_pw],
                [self.Z_wp, self.Z_rwg]
            ])
        else:
            self.Z_full = self.Z_peec
            self.Z_pw = None
            self.Z_wp = None
            self.Z_rwg = None

        self._N_p = N_p
        self._N_w = N_w
        self._work_interior = work_interior

    def solve(self, verbose=True):
        """
        Solve the hybrid PEEC-RWG system.

        Returns:
            result: Dictionary with solution data
        """
        if verbose:
            print("\nHybrid PEEC + RWG Solver")
            print("=" * 50)
            print(f"Frequency: {self.frequency/1e3:.1f} kHz")

        # Assemble
        self.assemble_system(verbose)

        N_p = self._N_p
        N_w = self._N_w

        # RHS
        rhs = np.zeros(N_p + N_w, dtype=complex)

        if self.port_segments:
            V_per_seg = self.voltage_excitation / len(self.port_segments)
            for seg_idx in self.port_segments:
                if 0 <= seg_idx < N_p:
                    rhs[seg_idx] = V_per_seg
        else:
            # Default: excite first segment
            rhs[0] = self.voltage_excitation

        # Solve
        if verbose:
            print("  Solving linear system...")

        try:
            solution = np.linalg.solve(self.Z_full, rhs)
        except np.linalg.LinAlgError:
            solution, _ = gmres(self.Z_full, rhs, tol=1e-6)

        self.I_peec = solution[:N_p]
        self.I_rwg = solution[N_p:] if N_w > 0 else np.array([])

        # Compute impedance
        if self.port_segments:
            I_total = sum(self.I_peec[i] for i in self.port_segments if i < N_p)
        else:
            I_total = self.I_peec[0]

        if abs(I_total) > 1e-30:
            self.impedance = self.voltage_excitation / I_total
        else:
            self.impedance = np.inf + 0j

        # Compute power in workpiece
        if N_w > 0 and self.Z_rwg is not None:
            self.power_workpiece = 0.5 * np.real(
                np.conj(self.I_rwg) @ self.Z_rwg @ self.I_rwg
            )
        else:
            self.power_workpiece = 0.0

        result = {
            'impedance': self.impedance,
            'resistance': self.impedance.real,
            'reactance': self.impedance.imag,
            'inductance': self.impedance.imag / (2 * np.pi * self.frequency) if self.frequency > 0 else 0,
            'power_workpiece': self.power_workpiece,
            'peec_dof': N_p,
            'rwg_dof': N_w,
        }

        if verbose:
            print()
            print("Results:")
            print(f"  Impedance: {result['resistance']:.6f} + j{result['reactance']:.6f} Ohm")
            print(f"  Inductance: {result['inductance']*1e6:.3f} uH")
            if N_w > 0:
                print(f"  Power in workpiece: {result['power_workpiece']:.3f} W")
            print()

        return result

    def compute_B_field(self, point):
        """
        Compute B field at a point from both PEEC and RWG currents.

        Parameters:
            point: [x, y, z] observation point

        Returns:
            B: [Bx, By, Bz] complex magnetic field [T]
        """
        point = np.array(point)
        B = np.zeros(3, dtype=complex)

        # Contribution from PEEC segments (Biot-Savart)
        if self.I_peec is not None:
            for i, segment in enumerate(self.coil_segments):
                I = self.I_peec[i]

                # Segment midpoint
                c = (segment['start'] + segment['end']) / 2
                d = segment['end'] - segment['start']
                L = segment['length']

                if L < 1e-15:
                    continue

                # Biot-Savart for segment
                R_vec = point - c
                R = np.linalg.norm(R_vec)

                if R > 1e-12:
                    dB = (mu_0 / (4 * np.pi)) * I * np.cross(d, R_vec) / (R ** 3)
                    B += dB

        # Contribution from RWG currents
        if self.I_rwg is not None and len(self.I_rwg) > 0 and self.workpiece_mesh is not None:
            # Use RWG field computation
            quad_pts = np.array([
                [0.333333333333333, 0.333333333333333],
            ])
            quad_wts = np.array([1.0])

            for unknown_idx, edge_idx in enumerate(self._work_interior):
                edge = self.workpiece_mesh.edges[edge_idx]
                I_n = self.I_rwg[unknown_idx]

                tris = []
                if edge.tri_plus >= 0:
                    tris.append((edge.tri_plus, +1, edge.local_idx_plus))
                if edge.tri_minus >= 0:
                    tris.append((edge.tri_minus, -1, edge.local_idx_minus))

                for tri_idx, sign, local_idx in tris:
                    tri = self.workpiece_mesh.triangles[tri_idx]
                    v = np.array([self.workpiece_mesh.vertices[k] for k in tri.vertex_indices])
                    free = self.workpiece_mesh.vertices[tri.vertex_indices[local_idx]]
                    coeff = edge.length / (2.0 * tri.area)

                    for q in range(len(quad_wts)):
                        L1, L2 = quad_pts[q]
                        L3 = 1 - L1 - L2
                        r_src = L1 * v[0] + L2 * v[1] + L3 * v[2]

                        if sign > 0:
                            J = coeff * (r_src - free)
                        else:
                            J = coeff * (free - r_src)

                        J *= I_n

                        R_vec = point - r_src
                        R = np.linalg.norm(R_vec)

                        if R > 1e-12:
                            dB = np.cross(J, R_vec) / (R ** 3)
                            B += quad_wts[q] * tri.area * dB

            B *= mu_0 / (4 * np.pi)

        return B


def create_nonlinear_solver(coil_mesh, workpiece_mesh, bh_curve, sigma, frequency, **kwargs):
    """
    Create a nonlinear coupled EFIE solver with ESIM.

    This is a convenience function to create a NonlinearCoupledEFIESolver
    with automatic ESI table generation from a B-H curve.

    Parameters:
        coil_mesh: RWGMesh for coil
        workpiece_mesh: RWGMesh for workpiece
        bh_curve: B-H curve data [[H1, B1], [H2, B2], ...] (H in A/m, B in T)
        sigma: Workpiece conductivity [S/m]
        frequency: Operating frequency [Hz]
        **kwargs: Additional parameters for solver

    Returns:
        solver: NonlinearCoupledEFIESolver object

    Example:
        >>> # Steel B-H curve
        >>> bh_curve = [
        ...     [0, 0], [100, 0.5], [500, 1.0], [1000, 1.3],
        ...     [5000, 1.6], [10000, 1.8], [50000, 2.0]
        ... ]
        >>> solver = create_nonlinear_solver(
        ...     coil_mesh, workpiece_mesh,
        ...     bh_curve, sigma=5e6, frequency=50e3
        ... )
        >>> result = solver.solve()
    """
    try:
        from .esim_cell_problem import ESITable
    except ImportError:
        from esim_cell_problem import ESITable

    # Generate ESI table from B-H curve
    esi_table = ESITable(bh_curve=bh_curve, sigma=sigma, frequency=frequency)

    # Create solver
    solver = NonlinearCoupledEFIESolver(coil_mesh, workpiece_mesh, esi_table)
    solver.set_frequency(frequency)
    solver.workpiece_conductivity = sigma

    # Apply additional parameters
    if 'coil_conductivity' in kwargs:
        solver.coil_conductivity = kwargs['coil_conductivity']
    if 'max_iterations' in kwargs:
        solver.max_iterations = kwargs['max_iterations']
    if 'tolerance' in kwargs:
        solver.tolerance = kwargs['tolerance']
    if 'port_edges' in kwargs:
        solver.define_port(kwargs['port_edges'])
    if 'voltage' in kwargs:
        solver.set_voltage_excitation(kwargs['voltage'])

    return solver


def create_induction_heating_model_with_esim(
    coil_type='loop',
    workpiece_type='plate',
    bh_curve=None,
    **kwargs
):
    """
    Create an induction heating model with ESIM nonlinear material.

    This function creates coil and workpiece meshes, generates the ESI table,
    and sets up the NonlinearCoupledEFIESolver.

    Parameters:
        coil_type: 'loop' or 'spiral'
        workpiece_type: 'plate', 'disk', or None
        bh_curve: B-H curve data [[H1, B1], [H2, B2], ...] (required)
        **kwargs: Geometry and physics parameters

    Geometry Parameters (all in meters):
        coil_center: Coil center position
        loop_radius: Loop radius (for 'loop' type)
        inner_radius, outer_radius, pitch, num_turns: For 'spiral' type
        wire_radius: Wire radius
        workpiece_center: Workpiece center
        workpiece_Lx, workpiece_Ly: Plate dimensions
        workpiece_radius: Disk radius

    Physics Parameters:
        frequency: Operating frequency [Hz]
        coil_conductivity: Coil conductivity [S/m] (default: 5.8e7 for copper)
        workpiece_conductivity: Workpiece conductivity [S/m] (default: 5e6)

    Mesh Parameters:
        num_around: Panels around wire
        num_along: Panels along coil
        nx, ny: Workpiece mesh divisions

    Returns:
        solver: NonlinearCoupledEFIESolver object

    Example:
        >>> # Steel B-H curve
        >>> bh_curve = [
        ...     [0, 0], [100, 0.5], [500, 1.0], [1000, 1.3],
        ...     [5000, 1.6], [10000, 1.8]
        ... ]
        >>> solver = create_induction_heating_model_with_esim(
        ...     coil_type='loop',
        ...     workpiece_type='plate',
        ...     bh_curve=bh_curve,
        ...     coil_center=[0, 0, 0.02],
        ...     loop_radius=0.04,
        ...     frequency=50e3,
        ...     workpiece_conductivity=5e6,
        ... )
        >>> result = solver.solve()
    """
    if bh_curve is None:
        raise ValueError("bh_curve is required for ESIM model")

    # Default parameters
    defaults = {
        # Coil
        'coil_center': [0, 0, 0.02],
        'loop_radius': 0.05,
        'inner_radius': 0.02,
        'outer_radius': 0.05,
        'pitch': 0.005,
        'num_turns': 3,
        'wire_radius': 0.002,
        'num_around': 8,
        'num_along': 24,

        # Workpiece
        'workpiece_center': [0, 0, 0],
        'workpiece_Lx': 0.1,
        'workpiece_Ly': 0.1,
        'workpiece_radius': 0.05,
        'nx': 10,
        'ny': 10,
        'nr': 5,
        'ntheta': 16,

        # Physics
        'frequency': 50e3,
        'coil_conductivity': 5.8e7,
        'workpiece_conductivity': 5e6,

        # Solver
        'max_iterations': 20,
        'tolerance': 1e-3,
    }
    defaults.update(kwargs)

    # Create coil mesh
    coil_mesh = RWGMesh()
    if coil_type == 'loop':
        coil_mesh.create_loop_coil(
            center=defaults['coil_center'],
            radius=defaults['loop_radius'],
            normal=[0, 0, 1],
            wire_radius=defaults['wire_radius'],
            num_around=defaults['num_around'],
            num_along=defaults['num_along']
        )
    elif coil_type == 'spiral':
        coil_mesh.create_spiral_coil(
            center=defaults['coil_center'],
            inner_radius=defaults['inner_radius'],
            outer_radius=defaults['outer_radius'],
            pitch=defaults['pitch'],
            num_turns=defaults['num_turns'],
            axis=[0, 0, 1],
            wire_radius=defaults['wire_radius'],
            num_around=defaults['num_around']
        )

    # Create workpiece mesh
    workpiece_mesh = None
    if workpiece_type == 'plate':
        workpiece_mesh = RWGMesh()
        workpiece_mesh.create_rectangular_plate(
            defaults['workpiece_center'],
            defaults['workpiece_Lx'],
            defaults['workpiece_Ly'],
            [0, 0, 1],  # normal
            defaults['nx'],
            defaults['ny']
        )
    elif workpiece_type == 'disk':
        workpiece_mesh = RWGMesh()
        workpiece_mesh.create_circular_disk(
            defaults['workpiece_center'],
            defaults['workpiece_radius'],
            [0, 0, 1],  # normal
            defaults['nr'],
            defaults['ntheta']
        )

    # Create solver with ESIM
    solver = create_nonlinear_solver(
        coil_mesh,
        workpiece_mesh,
        bh_curve,
        defaults['workpiece_conductivity'],
        defaults['frequency'],
        coil_conductivity=defaults['coil_conductivity'],
        max_iterations=defaults['max_iterations'],
        tolerance=defaults['tolerance'],
    )

    # Find and set port
    interior_edges = [i for i, e in enumerate(coil_mesh.edges) if not e.is_boundary]
    if interior_edges:
        solver.define_port([interior_edges[0]])

    return solver


def create_hybrid_solver(coil_type='loop', workpiece_type='plate', **kwargs):
    """
    Create a hybrid PEEC + RWG solver for induction heating.

    Parameters:
        coil_type: 'loop' or 'spiral'
        workpiece_type: 'plate', 'disk', or None
        **kwargs: Geometry parameters

    Returns:
        solver: HybridPEECRWGSolver object
    """
    defaults = {
        # Coil
        'coil_center': [0, 0, 0.02],
        'coil_radius': 0.05,
        'coil_inner_radius': 0.02,
        'coil_outer_radius': 0.05,
        'coil_pitch': 0.005,
        'coil_num_turns': 3,
        'coil_axis': [0, 0, 1],
        'coil_n_segments': 36,
        'wire_radius': 0.002,

        # Workpiece
        'workpiece_center': [0, 0, 0],
        'workpiece_Lx': 0.1,
        'workpiece_Ly': 0.1,
        'workpiece_radius': 0.05,
        'workpiece_normal': [0, 0, 1],
        'nx': 8,
        'ny': 8,
        'nr': 4,
        'ntheta': 12,

        # Physics
        'frequency': 50e3,
        'coil_conductivity': 5.8e7,
        'workpiece_conductivity': 5e6,
    }
    defaults.update(kwargs)

    solver = HybridPEECRWGSolver()
    solver.set_frequency(defaults['frequency'])
    solver.coil_conductivity = defaults['coil_conductivity']
    solver.workpiece_conductivity = defaults['workpiece_conductivity']
    solver.coil_wire_radius = defaults['wire_radius']

    # Create coil
    if coil_type == 'loop':
        solver.create_circular_coil_peec(
            defaults['coil_center'],
            defaults['coil_radius'],
            defaults['coil_axis'],
            defaults['coil_n_segments']
        )
    elif coil_type == 'spiral':
        solver.create_spiral_coil_peec(
            defaults['coil_center'],
            defaults['coil_inner_radius'],
            defaults['coil_outer_radius'],
            defaults['coil_pitch'],
            defaults['coil_num_turns'],
            defaults['coil_axis'],
            defaults['coil_n_segments']
        )

    # Create workpiece
    if workpiece_type == 'plate':
        solver.create_workpiece_plate(
            defaults['workpiece_center'],
            defaults['workpiece_Lx'],
            defaults['workpiece_Ly'],
            defaults['workpiece_normal'],
            defaults['nx'],
            defaults['ny']
        )
    elif workpiece_type == 'disk':
        solver.create_workpiece_disk(
            defaults['workpiece_center'],
            defaults['workpiece_radius'],
            defaults['workpiece_normal'],
            defaults['nr'],
            defaults['ntheta']
        )

    # Set port (first segment)
    solver.define_port([0])

    return solver


# ============================================================================
# LOOP-STAR Decomposition for Low-Frequency Stability
# ============================================================================

class LoopStarDecomposition:
    """
    Loop-Star Decomposition for RWG Basis Functions.

    This class implements the Loop-Star (or Helmholtz) decomposition of RWG
    basis functions to eliminate the low-frequency breakdown in EFIE.

    The decomposition separates the RWG current into:
    - Loop (solenoidal) component: divergence-free, represents magnetic effects
    - Star (irrotational) component: curl-free, represents electric charge effects

    For MQS (Magneto-Quasi-Static) problems, only the Loop basis is needed
    since there is no charge accumulation (div J = 0).

    Mathematical basis:
        J_RWG = T_LS * [J_loop; J_star]

    where T_LS is the Loop-Star transformation matrix.

    Benefits:
    - Eliminates low-frequency breakdown (omega -> 0 instability)
    - Separates inductive and capacitive effects
    - More stable conditioning at low frequencies

    Reference:
        G. Vecchi, "Loop-Star Decomposition of Basis Functions in the
        Discretization of the EFIE," IEEE TAP, 1999.

    Author: Radia Development Team
    Date: 2026-01-08
    """

    def __init__(self, mesh):
        """
        Initialize Loop-Star decomposition.

        Parameters:
            mesh: RWGMesh object
        """
        self.mesh = mesh

        # Counts
        self.num_loops = 0      # = num_triangles for closed surface
        self.num_stars = 0      # = num_vertices - 1 (for connected surface)
        self.num_rwg = 0        # = num_interior_edges

        # Transformation matrices (sparse)
        self.T_loop = None      # RWG-to-Loop transformation (N_rwg x N_loop)
        self.T_star = None      # RWG-to-Star transformation (N_rwg x N_star)
        self.T_ls = None        # Combined [T_loop | T_star]

        # Vertex-triangle connectivity
        self.vertex_triangles = None  # vertex_idx -> list of triangle indices
        self.vertex_edges = None      # vertex_idx -> list of (edge_idx, sign)

        # Spanning tree for Star basis
        self.spanning_tree_edges = None
        self.non_tree_edges = None

        self._build_decomposition()

    def _build_decomposition(self):
        """Build the Loop-Star decomposition matrices."""
        # Get interior edges only
        interior_edges = [i for i, e in enumerate(self.mesh.edges) if not e.is_boundary]
        self.num_rwg = len(interior_edges)

        if self.num_rwg == 0:
            warnings.warn("No interior edges found - mesh may be open")
            return

        # Build vertex connectivity
        self._build_vertex_connectivity()

        # Build spanning tree for Star basis
        self._build_spanning_tree()

        # Build Loop basis (triangle-based)
        self._build_loop_basis()

        # Build Star basis (vertex-based)
        self._build_star_basis()

        # Combine into T_ls
        if self.T_loop is not None and self.T_star is not None:
            from scipy import sparse
            self.T_ls = sparse.hstack([self.T_loop, self.T_star], format='csr')

    def _build_vertex_connectivity(self):
        """Build vertex-to-triangle and vertex-to-edge connectivity."""
        num_vertices = len(self.mesh.vertices)

        self.vertex_triangles = [[] for _ in range(num_vertices)]
        self.vertex_edges = [[] for _ in range(num_vertices)]

        # Vertex-triangle connectivity
        for tri_idx, tri in enumerate(self.mesh.triangles):
            for v_idx in tri.vertex_indices:
                self.vertex_triangles[v_idx].append(tri_idx)

        # Vertex-edge connectivity (only interior edges)
        for edge_idx, edge in enumerate(self.mesh.edges):
            if edge.is_boundary:
                continue

            v1, v2 = edge.vertex_indices

            # Edge points from v1 to v2, so:
            # - At v1: current flows "out" (+1)
            # - At v2: current flows "in" (-1)
            self.vertex_edges[v1].append((edge_idx, +1))
            self.vertex_edges[v2].append((edge_idx, -1))

    def _build_spanning_tree(self):
        """Build spanning tree for Star basis using BFS."""
        num_vertices = len(self.mesh.vertices)

        # Find vertices that are connected to interior edges
        connected_vertices = set()
        for edge_idx, edge in enumerate(self.mesh.edges):
            if not edge.is_boundary:
                connected_vertices.add(edge.vertex_indices[0])
                connected_vertices.add(edge.vertex_indices[1])

        if not connected_vertices:
            return

        # BFS to build spanning tree
        from collections import deque

        visited = set()
        tree_edges = []
        non_tree_edges = []

        start_vertex = next(iter(connected_vertices))
        queue = deque([start_vertex])
        visited.add(start_vertex)

        while queue:
            v = queue.popleft()

            for edge_idx, sign in self.vertex_edges[v]:
                edge = self.mesh.edges[edge_idx]
                v1, v2 = edge.vertex_indices

                # Find the other vertex
                other_v = v2 if v == v1 else v1

                if other_v not in visited:
                    visited.add(other_v)
                    queue.append(other_v)
                    tree_edges.append(edge_idx)
                elif edge_idx not in tree_edges:
                    # This edge closes a loop
                    if edge_idx not in non_tree_edges:
                        non_tree_edges.append(edge_idx)

        self.spanning_tree_edges = tree_edges
        self.non_tree_edges = non_tree_edges

        # Number of loops = number of non-tree edges (for closed surfaces)
        # Actually for EFIE, num_loops = num_triangles (use triangle loops)

    def _build_loop_basis(self):
        """
        Build Loop basis functions.

        Each triangle defines a loop: the circulation of RWG functions around
        the triangle's edges.

        For triangle t with edges e1, e2, e3:
            Loop_t = s1*f_e1 + s2*f_e2 + s3*f_e3

        where s_i = +1/-1 depending on edge orientation.

        This gives div(Loop_t) = 0 (solenoidal).
        """
        from scipy import sparse

        num_triangles = len(self.mesh.triangles)
        self.num_loops = num_triangles

        # Map edge index to column index in RWG matrix
        edge_to_col = {}
        col = 0
        for i, e in enumerate(self.mesh.edges):
            if not e.is_boundary:
                edge_to_col[i] = col
                col += 1

        # Build Loop transformation matrix (N_rwg x N_loops)
        # T_loop[i, j] = coefficient of RWG_i in Loop_j
        rows = []
        cols = []
        data = []

        for tri_idx, tri in enumerate(self.mesh.triangles):
            # Get the three edges of this triangle
            for local_edge in range(3):
                edge_idx = tri.edges[local_edge]
                edge = self.mesh.edges[edge_idx]

                if edge.is_boundary:
                    continue  # Skip boundary edges

                # Sign: +1 if edge direction matches triangle circulation
                # The sign is determined by whether this triangle is T+ or T-
                sign = tri.edge_signs[local_edge]

                if edge_idx in edge_to_col:
                    rows.append(edge_to_col[edge_idx])
                    cols.append(tri_idx)
                    data.append(float(sign))

        if rows:
            self.T_loop = sparse.csr_matrix(
                (data, (rows, cols)),
                shape=(self.num_rwg, self.num_loops)
            )
        else:
            self.T_loop = sparse.csr_matrix((self.num_rwg, 0))

    def _build_star_basis(self):
        """
        Build Star basis functions.

        Each interior vertex (except one reference) defines a star:
        the sum of RWG functions radiating from that vertex.

        For vertex v with adjacent edges e1, e2, ..., en:
            Star_v = sum_i s_i * f_ei

        where s_i = +1 if edge points away from v, -1 otherwise.

        This gives curl(Star_v) = 0 (irrotational).
        """
        from scipy import sparse

        # Find interior vertices (connected to interior edges)
        interior_vertices = set()
        for edge_idx, edge in enumerate(self.mesh.edges):
            if not edge.is_boundary:
                interior_vertices.add(edge.vertex_indices[0])
                interior_vertices.add(edge.vertex_indices[1])

        interior_vertices = sorted(interior_vertices)

        # Remove one vertex (reference) to avoid singular matrix
        # Use the last vertex as reference
        if len(interior_vertices) > 0:
            star_vertices = interior_vertices[:-1]
        else:
            star_vertices = []

        self.num_stars = len(star_vertices)

        # Map edge index to column index in RWG matrix
        edge_to_col = {}
        col = 0
        for i, e in enumerate(self.mesh.edges):
            if not e.is_boundary:
                edge_to_col[i] = col
                col += 1

        # Build Star transformation matrix (N_rwg x N_stars)
        rows = []
        cols = []
        data = []

        for star_idx, v_idx in enumerate(star_vertices):
            # Get edges adjacent to this vertex
            for edge_idx, sign in self.vertex_edges[v_idx]:
                if edge_idx in edge_to_col:
                    rows.append(edge_to_col[edge_idx])
                    cols.append(star_idx)
                    data.append(float(sign))

        if rows:
            self.T_star = sparse.csr_matrix(
                (data, (rows, cols)),
                shape=(self.num_rwg, self.num_stars)
            )
        else:
            self.T_star = sparse.csr_matrix((self.num_rwg, 0))

    def transform_to_loop_star(self, Z_rwg):
        """
        Transform RWG impedance matrix to Loop-Star basis.

        Z_ls = T_ls^T * Z_rwg * T_ls

        Parameters:
            Z_rwg: RWG impedance matrix (N_rwg x N_rwg)

        Returns:
            Z_ls: Loop-Star impedance matrix
        """
        if self.T_ls is None:
            return Z_rwg

        from scipy import sparse

        # Z_ls = T_ls^T * Z_rwg * T_ls
        if sparse.issparse(Z_rwg):
            Z_ls = self.T_ls.T @ Z_rwg @ self.T_ls
        else:
            Z_ls = self.T_ls.T @ Z_rwg @ self.T_ls.toarray()

        return Z_ls

    def transform_to_rwg(self, I_ls):
        """
        Transform Loop-Star currents back to RWG basis.

        I_rwg = T_ls * I_ls

        Parameters:
            I_ls: Loop-Star current vector

        Returns:
            I_rwg: RWG current vector
        """
        if self.T_ls is None:
            return I_ls

        return self.T_ls @ I_ls

    def transform_rhs_to_loop_star(self, V_rwg):
        """
        Transform RWG RHS to Loop-Star basis.

        V_ls = T_ls^T * V_rwg

        Parameters:
            V_rwg: RWG voltage vector

        Returns:
            V_ls: Loop-Star voltage vector
        """
        if self.T_ls is None:
            return V_rwg

        return self.T_ls.T @ V_rwg

    def get_loop_matrix(self, Z_rwg):
        """
        Get the Loop-Loop block of the impedance matrix.

        Z_LL = T_loop^T * Z_rwg * T_loop

        This is useful for MQS where only loops contribute.
        """
        if self.T_loop is None:
            return None

        from scipy import sparse

        if sparse.issparse(Z_rwg):
            return self.T_loop.T @ Z_rwg @ self.T_loop
        else:
            return self.T_loop.T @ Z_rwg @ self.T_loop.toarray()

    def solve_mqs_loop_only(self, Z_rwg, V_rwg):
        """
        Solve MQS problem using only Loop basis.

        For MQS: div(J) = 0, so only solenoidal currents exist.
        This reduces the problem size from N_rwg to N_loops.

        Parameters:
            Z_rwg: RWG impedance matrix
            V_rwg: RWG voltage vector

        Returns:
            I_rwg: RWG current solution
        """
        if self.T_loop is None or self.T_loop.shape[1] == 0:
            warnings.warn("Loop basis not available")
            return None

        # Transform to Loop basis
        Z_loop = self.get_loop_matrix(Z_rwg)
        V_loop = self.T_loop.T @ V_rwg

        # Solve in Loop space
        from scipy.linalg import solve
        import numpy as np

        Z_loop_dense = Z_loop.toarray() if hasattr(Z_loop, 'toarray') else Z_loop
        I_loop = solve(Z_loop_dense, V_loop)

        # Transform back to RWG
        I_rwg = self.T_loop @ I_loop

        return I_rwg


class LoopStarEFIESolver:
    """
    EFIE Solver with Loop-Star decomposition for low-frequency stability.

    This solver uses Loop-Star decomposition to avoid the low-frequency
    breakdown that occurs in standard EFIE when omega -> 0.

    The decomposition separates the problem into:
    - Inductive (Loop) part: j*omega*L
    - Capacitive (Star) part: 1/(j*omega*C)

    At low frequencies, the standard formulation has:
        Z = j*omega*L + 1/(j*omega)*P

    where L is bounded but 1/(j*omega)*P -> infinity as omega -> 0.

    With Loop-Star decomposition:
    - Loop-Loop block: j*omega*L_LL (inductive)
    - Star-Star block: 1/(j*omega)*P_SS (capacitive)
    - Cross terms: bounded

    For MQS problems (div J = 0), only the Loop equations are needed.

    Reference:
        F.P. Andriulli et al., "A Multiplicative Calderon Preconditioner
        for the EFIE," IEEE TAP, 2008.
    """

    def __init__(self, mesh):
        """
        Initialize Loop-Star EFIE solver.

        Parameters:
            mesh: RWGMesh object
        """
        self.mesh = mesh
        self.decomposition = LoopStarDecomposition(mesh)

        self.frequency = 1e6
        self.omega = 2 * np.pi * self.frequency
        self.conductivity = 5.8e7

        # Matrices in RWG basis
        self.L_matrix = None
        self.R_matrix = None
        self.P_matrix = None  # Scalar potential (capacitive) matrix

        # Matrices in Loop-Star basis
        self.Z_LL = None  # Loop-Loop
        self.Z_LS = None  # Loop-Star
        self.Z_SL = None  # Star-Loop
        self.Z_SS = None  # Star-Star

        # Solution
        self.solution = None
        self.I_loop = None
        self.I_star = None
        self.impedance = 0j
        self.total_current = 0j

        # Port
        self.port_edges = []
        self.voltage_excitation = 1.0 + 0j

        # Solver settings
        self.use_mqs = True  # Use MQS approximation (Loop-only)
        self.use_full_wave = False  # Use full EFIE (Loop + Star)

    def set_frequency(self, frequency):
        """Set operating frequency [Hz]."""
        self.frequency = frequency
        self.omega = 2 * np.pi * frequency

    def set_conductivity(self, sigma):
        """Set material conductivity [S/m]."""
        self.conductivity = sigma

    def define_port(self, edge_indices):
        """Define port for impedance calculation."""
        self.port_edges = list(edge_indices)

    def set_voltage_excitation(self, V):
        """Set voltage excitation at port [V]."""
        self.voltage_excitation = complex(V)

    def assemble_matrices(self, verbose=False):
        """Assemble impedance matrices in Loop-Star basis."""
        if verbose:
            print("Assembling Loop-Star matrices...")

        # First, assemble RWG matrices using base solver
        rwg_solver = RWGEFIESolver(self.mesh)
        rwg_solver.set_frequency(self.frequency)
        rwg_solver.set_conductivity(self.conductivity)

        # Get interior edges
        interior_edges = [i for i, e in enumerate(self.mesh.edges) if not e.is_boundary]
        N = len(interior_edges)

        if N == 0:
            warnings.warn("No interior edges")
            return

        # Assemble L matrix (inductance)
        if verbose:
            print("  Assembling L matrix...")

        self.L_matrix = np.zeros((N, N), dtype=complex)
        edge_to_idx = {e: i for i, e in enumerate(interior_edges)}

        for m, edge_m_idx in enumerate(interior_edges):
            edge_m = self.mesh.edges[edge_m_idx]
            for n, edge_n_idx in enumerate(interior_edges):
                edge_n = self.mesh.edges[edge_n_idx]
                self.L_matrix[m, n] = rwg_solver._compute_L_element(edge_m, edge_n)

        # Assemble R matrix (resistance)
        if verbose:
            print("  Assembling R matrix...")

        self.R_matrix = np.zeros((N, N), dtype=complex)
        for m, edge_m_idx in enumerate(interior_edges):
            edge_m = self.mesh.edges[edge_m_idx]
            for n, edge_n_idx in enumerate(interior_edges):
                edge_n = self.mesh.edges[edge_n_idx]
                self.R_matrix[m, n] = rwg_solver._compute_R_element(edge_m, edge_n)

        # Assemble P matrix (scalar potential)
        # P_mn = (1/4*pi*eps_0) * integral{ div(f_m) * div(f_n) / |r-r'| } dA dA'
        if verbose:
            print("  Assembling P matrix...")

        self.P_matrix = self._assemble_P_matrix(interior_edges)

        # Transform to Loop-Star basis
        if verbose:
            print("  Transforming to Loop-Star basis...")

        self._transform_to_loop_star()

        if verbose:
            print(f"  Loop DOF: {self.decomposition.num_loops}")
            print(f"  Star DOF: {self.decomposition.num_stars}")
            print(f"  Total DOF: {self.decomposition.num_loops + self.decomposition.num_stars}")

    def _assemble_P_matrix(self, interior_edges):
        """Assemble scalar potential matrix P."""
        N = len(interior_edges)
        P = np.zeros((N, N), dtype=complex)

        eps_0 = 8.854e-12  # Permittivity of free space

        # Quadrature points
        quad_pts = np.array([
            [0.333333333333333, 0.333333333333333],
            [0.059715871789770, 0.470142064105115],
            [0.470142064105115, 0.059715871789770],
            [0.470142064105115, 0.470142064105115],
        ])
        quad_wts = np.array([0.225, 0.132394, 0.132394, 0.132394])

        for m, edge_m_idx in enumerate(interior_edges):
            edge_m = self.mesh.edges[edge_m_idx]

            # Divergence of RWG = l_n / A_+ in T+, -l_n / A_- in T-
            div_m_plus = edge_m.length / self.mesh.triangles[edge_m.tri_plus].area if edge_m.tri_plus >= 0 else 0
            div_m_minus = -edge_m.length / self.mesh.triangles[edge_m.tri_minus].area if edge_m.tri_minus >= 0 else 0

            for n, edge_n_idx in enumerate(interior_edges):
                edge_n = self.mesh.edges[edge_n_idx]

                div_n_plus = edge_n.length / self.mesh.triangles[edge_n.tri_plus].area if edge_n.tri_plus >= 0 else 0
                div_n_minus = -edge_n.length / self.mesh.triangles[edge_n.tri_minus].area if edge_n.tri_minus >= 0 else 0

                P_mn = 0.0

                # Integrate over all triangle pairs
                tris_m = []
                if edge_m.tri_plus >= 0:
                    tris_m.append((edge_m.tri_plus, div_m_plus))
                if edge_m.tri_minus >= 0:
                    tris_m.append((edge_m.tri_minus, div_m_minus))

                tris_n = []
                if edge_n.tri_plus >= 0:
                    tris_n.append((edge_n.tri_plus, div_n_plus))
                if edge_n.tri_minus >= 0:
                    tris_n.append((edge_n.tri_minus, div_n_minus))

                for tri_m_idx, div_m in tris_m:
                    tri_m = self.mesh.triangles[tri_m_idx]
                    v_m = np.array([self.mesh.vertices[i] for i in tri_m.vertex_indices])

                    for tri_n_idx, div_n in tris_n:
                        tri_n = self.mesh.triangles[tri_n_idx]
                        v_n = np.array([self.mesh.vertices[i] for i in tri_n.vertex_indices])

                        # Double quadrature
                        for q_m in range(len(quad_wts)):
                            L1_m, L2_m = quad_pts[q_m]
                            L3_m = 1 - L1_m - L2_m
                            r_m = L1_m * v_m[0] + L2_m * v_m[1] + L3_m * v_m[2]

                            for q_n in range(len(quad_wts)):
                                L1_n, L2_n = quad_pts[q_n]
                                L3_n = 1 - L1_n - L2_n
                                r_n = L1_n * v_n[0] + L2_n * v_n[1] + L3_n * v_n[2]

                                # Green's function
                                R = np.linalg.norm(r_m - r_n)
                                if R < 1e-12:
                                    R = 0.5 * np.sqrt(tri_m.area / np.pi)
                                G = 1.0 / R

                                # Integrand: div_m * div_n * G
                                integrand = div_m * div_n * G

                                w = quad_wts[q_m] * quad_wts[q_n] * tri_m.area * tri_n.area
                                P_mn += w * integrand

                P[m, n] = P_mn / (4 * np.pi * eps_0)

        return P

    def _transform_to_loop_star(self):
        """Transform matrices to Loop-Star basis."""
        T_loop = self.decomposition.T_loop
        T_star = self.decomposition.T_star

        if T_loop is None:
            return

        # Z_rwg = j*omega*L + R + P/(j*omega)
        # For low-frequency stability, we work with L, R, P separately

        # Loop-Loop block: T_loop^T * X * T_loop
        self.L_LL = T_loop.T @ self.L_matrix @ T_loop
        self.R_LL = T_loop.T @ self.R_matrix @ T_loop
        self.P_LL = T_loop.T @ self.P_matrix @ T_loop

        if T_star is not None and T_star.shape[1] > 0:
            # Star-Star block
            self.L_SS = T_star.T @ self.L_matrix @ T_star
            self.R_SS = T_star.T @ self.R_matrix @ T_star
            self.P_SS = T_star.T @ self.P_matrix @ T_star

            # Cross blocks
            self.L_LS = T_loop.T @ self.L_matrix @ T_star
            self.R_LS = T_loop.T @ self.R_matrix @ T_star
            self.P_LS = T_loop.T @ self.P_matrix @ T_star

            self.L_SL = T_star.T @ self.L_matrix @ T_loop
            self.R_SL = T_star.T @ self.R_matrix @ T_loop
            self.P_SL = T_star.T @ self.P_matrix @ T_loop

    def solve(self, verbose=False):
        """
        Solve the EFIE using Loop-Star decomposition.

        For MQS (self.use_mqs=True): Only solve Loop equations.
        For full-wave (self.use_full_wave=True): Solve complete system.
        """
        if self.L_matrix is None:
            self.assemble_matrices(verbose)

        if self.use_mqs:
            return self._solve_mqs(verbose)
        else:
            return self._solve_full_wave(verbose)

    def _solve_mqs(self, verbose=False):
        """Solve using MQS approximation (Loop-only)."""
        if verbose:
            print("Solving MQS (Loop-only) system...")

        N_loop = self.decomposition.num_loops

        if N_loop == 0:
            warnings.warn("No Loop DOF available")
            return {'impedance': 0j}

        # MQS impedance: Z = j*omega*L + R (no P term)
        jOmega = 1j * self.omega
        Z_LL = jOmega * np.asarray(self.L_LL.todense()) + np.asarray(self.R_LL.todense())

        # Build RHS in Loop basis
        V_rwg = self._build_rhs_rwg()
        V_loop = self.decomposition.T_loop.T @ V_rwg

        # Solve
        from scipy.linalg import solve
        self.I_loop = solve(Z_LL, V_loop)

        # Transform back to RWG
        self.solution = self.decomposition.T_loop @ self.I_loop

        # Compute impedance
        self._compute_impedance()

        result = {
            'impedance': self.impedance,
            'resistance': self.impedance.real,
            'reactance': self.impedance.imag,
            'inductance': self.impedance.imag / self.omega if self.omega > 0 else 0,
            'current': self.total_current,
            'dof_loop': N_loop,
            'dof_star': 0,
        }

        if verbose:
            print(f"  Impedance: {result['resistance']:.6f} + j{result['reactance']:.6f} Ohm")
            print(f"  Inductance: {result['inductance']*1e6:.3f} uH")

        return result

    def _solve_full_wave(self, verbose=False):
        """Solve full Loop-Star system."""
        if verbose:
            print("Solving full Loop-Star system...")

        N_loop = self.decomposition.num_loops
        N_star = self.decomposition.num_stars
        N_total = N_loop + N_star

        if N_total == 0:
            return {'impedance': 0j}

        jOmega = 1j * self.omega
        inv_jOmega = 1.0 / jOmega if abs(jOmega) > 1e-30 else 0

        # Build full Loop-Star system
        Z_ls = np.zeros((N_total, N_total), dtype=complex)

        # Loop-Loop block
        Z_ls[:N_loop, :N_loop] = (
            jOmega * np.asarray(self.L_LL.todense()) +
            np.asarray(self.R_LL.todense()) +
            inv_jOmega * np.asarray(self.P_LL.todense())
        )

        if N_star > 0:
            # Star-Star block
            Z_ls[N_loop:, N_loop:] = (
                jOmega * np.asarray(self.L_SS.todense()) +
                np.asarray(self.R_SS.todense()) +
                inv_jOmega * np.asarray(self.P_SS.todense())
            )

            # Cross blocks
            Z_ls[:N_loop, N_loop:] = (
                jOmega * np.asarray(self.L_LS.todense()) +
                np.asarray(self.R_LS.todense()) +
                inv_jOmega * np.asarray(self.P_LS.todense())
            )
            Z_ls[N_loop:, :N_loop] = (
                jOmega * np.asarray(self.L_SL.todense()) +
                np.asarray(self.R_SL.todense()) +
                inv_jOmega * np.asarray(self.P_SL.todense())
            )

        # Build RHS
        V_rwg = self._build_rhs_rwg()
        V_loop = self.decomposition.T_loop.T @ V_rwg
        V_star = self.decomposition.T_star.T @ V_rwg if N_star > 0 else np.array([])

        V_ls = np.concatenate([V_loop, V_star])

        # Solve
        from scipy.linalg import solve
        I_ls = solve(Z_ls, V_ls)

        self.I_loop = I_ls[:N_loop]
        self.I_star = I_ls[N_loop:] if N_star > 0 else np.array([])

        # Transform back to RWG
        self.solution = self.decomposition.transform_to_rwg(I_ls)

        # Compute impedance
        self._compute_impedance()

        result = {
            'impedance': self.impedance,
            'resistance': self.impedance.real,
            'reactance': self.impedance.imag,
            'inductance': self.impedance.imag / self.omega if self.omega > 0 else 0,
            'current': self.total_current,
            'dof_loop': N_loop,
            'dof_star': N_star,
        }

        if verbose:
            print(f"  Impedance: {result['resistance']:.6f} + j{result['reactance']:.6f} Ohm")
            print(f"  Inductance: {result['inductance']*1e6:.3f} uH")

        return result

    def _build_rhs_rwg(self):
        """Build RHS vector in RWG basis."""
        interior_edges = [i for i, e in enumerate(self.mesh.edges) if not e.is_boundary]
        N = len(interior_edges)

        V_rwg = np.zeros(N, dtype=complex)

        # Map edge index to unknown index
        edge_to_idx = {e: i for i, e in enumerate(interior_edges)}

        # Set voltage at port edges
        for edge_idx in self.port_edges:
            if edge_idx in edge_to_idx:
                V_rwg[edge_to_idx[edge_idx]] = self.voltage_excitation

        return V_rwg

    def _compute_impedance(self):
        """Compute port impedance from solution."""
        if self.solution is None:
            self.impedance = 0j
            self.total_current = 0j
            return

        interior_edges = [i for i, e in enumerate(self.mesh.edges) if not e.is_boundary]
        edge_to_idx = {e: i for i, e in enumerate(interior_edges)}

        I_total = 0j
        for edge_idx in self.port_edges:
            if edge_idx in edge_to_idx:
                I_total += self.solution[edge_to_idx[edge_idx]]

        self.total_current = I_total
        if abs(I_total) > 1e-30:
            self.impedance = self.voltage_excitation / I_total
        else:
            self.impedance = np.inf + 0j


# Example usage
if __name__ == "__main__":
    print("RWG-EFIE Solver Test")
    print("=" * 60)

    # Create a simple loop coil
    solver = create_induction_heating_model(
        coil_type='loop',
        workpiece_type=None,
        loop_radius=0.05,
        wire_radius=0.002,
        frequency=50e3,
        num_around=8
    )

    result = solver.solve(verbose=True)

    # Compute B field on axis
    print("B field on axis (z > coil):")
    for z in [0.03, 0.05, 0.07, 0.10]:
        B = solver.compute_B_field([0, 0, z])
        B_mag = np.abs(np.linalg.norm(B))
        print(f"  z={z*100:.0f}cm: |B| = {B_mag*1e6:.3f} uT")

    print()
    print("Test completed!")
