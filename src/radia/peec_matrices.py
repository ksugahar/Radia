"""
PEEC Matrix Construction Python Wrapper

Re-exports from C++ pybind11 module (peec_matrices.pyd).

This module provides Loop-Star PEEC matrix construction for quasi-static
electromagnetic analysis using Darwin approximation.

Usage:
    from radia.peec_matrices import PEECBuilder, create_wire_peec, create_loop_peec

    # Method 1: Using PEECBuilder class
    builder = PEECBuilder()
    builder.create_wire([0,0,0], [0.1,0,0], 1e-3, 1e-3, 10)
    L, R, P, M_LS = builder.build()

    # Port impedance at 1 MHz
    Z = builder.compute_impedance(1e6, [1.0] + [0.0]*9)

    # Method 2: Using convenience functions
    L, R, P, M_LS = create_wire_peec([0,0,0], [0.1,0,0], 1e-3, 1e-3, 10)
"""

import numpy as np

# Try to import from pybind11 module (C++ implementation)
try:
    from peec_matrices import PEECBuilder, create_wire_peec, create_loop_peec
    _USE_CPP = True
except ImportError:
    # Fall back to pure Python implementation if C++ module not available
    _USE_CPP = False


if not _USE_CPP:
    # Pure Python fallback implementation
    # This is used when the C++ module is not available

    class PEECBuilder:
        """
        PEEC Matrix Builder (Pure Python Fallback)

        Constructs L, P, R, M_LS matrices from conductor segments.
        Uses Darwin approximation: G(r) = 1/(4*pi*r)
        """

        # Physical constants
        MU_0 = 4 * np.pi * 1e-7  # H/m
        EPS_0 = 8.854187817e-12  # F/m
        INV_FOUR_PI = 1.0 / (4 * np.pi)

        def __init__(self):
            self.segments = []
            self.nodes = []
            self._matrices_built = False
            self._L = None
            self._R = None
            self._P = None
            self._M_LS = None
            self._Zs_diag = None

        def add_segment(self, center, direction, length, width, height, sigma=5.8e7):
            """Add a conductor segment."""
            self.segments.append({
                'center': np.array(center),
                'direction': np.array(direction) / np.linalg.norm(direction),
                'length': length,
                'width': width,
                'height': height,
                'sigma': sigma
            })
            self._matrices_built = False

        def add_node(self, position, area):
            """Add a charge node."""
            self.nodes.append({
                'position': np.array(position),
                'area': area
            })
            self._matrices_built = False

        def create_wire(self, start, end, width, height, n_segments, sigma=5.8e7):
            """Create wire segments along a straight path."""
            start = np.array(start)
            end = np.array(end)
            direction = end - start
            total_length = np.linalg.norm(direction)
            direction = direction / total_length
            seg_length = total_length / n_segments

            for i in range(n_segments):
                center = start + direction * (i + 0.5) * seg_length
                self.add_segment(center, direction, seg_length, width, height, sigma)

            self._matrices_built = False
            return n_segments

        def create_loop(self, center, radius, normal, width, height, n_segments, sigma=5.8e7):
            """Create circular loop segments."""
            center = np.array(center)
            normal = np.array(normal)
            normal = normal / np.linalg.norm(normal)

            # Create orthogonal basis
            if abs(normal[2]) < 0.9:
                u = np.cross(normal, [0, 0, 1])
            else:
                u = np.cross(normal, [1, 0, 0])
            u = u / np.linalg.norm(u)
            v = np.cross(normal, u)

            dtheta = 2 * np.pi / n_segments
            for i in range(n_segments):
                theta = i * dtheta + dtheta / 2
                seg_center = center + radius * (np.cos(theta) * u + np.sin(theta) * v)
                # Tangent direction
                direction = -np.sin(theta) * u + np.cos(theta) * v
                seg_length = radius * dtheta
                self.add_segment(seg_center, direction, seg_length, width, height, sigma)

            self._matrices_built = False
            return n_segments

        def _auto_generate_nodes(self):
            """Auto-generate nodes from segment endpoints."""
            if not self.segments:
                return

            node_positions = []
            tol = 1e-10

            for seg in self.segments:
                p1 = seg['center'] - seg['direction'] * seg['length'] / 2
                p2 = seg['center'] + seg['direction'] * seg['length'] / 2
                avg_area = seg['width'] * seg['height']

                for p in [p1, p2]:
                    is_new = True
                    for existing in node_positions:
                        if np.linalg.norm(p - existing['position']) < tol:
                            is_new = False
                            break
                    if is_new:
                        node_positions.append({'position': p, 'area': avg_area})

            self.nodes = node_positions

        def _compute_self_inductance(self, seg):
            """Compute self-inductance using GMD approximation."""
            l = seg['length']
            a = np.sqrt(seg['width'] * seg['height'] / np.pi)  # Equivalent radius
            if a < 1e-15:
                return 0.0
            return self.MU_0 * l * (np.log(2 * l / a) - 1) * self.INV_FOUR_PI

        def _compute_mutual_inductance(self, seg_i, seg_j):
            """Compute mutual inductance using Neumann formula."""
            ci = seg_i['center']
            cj = seg_j['center']
            di = seg_i['direction']
            dj = seg_j['direction']

            r = np.linalg.norm(cj - ci)
            if r < 1e-15:
                return 0.0

            cos_angle = np.dot(di, dj)
            li = seg_i['length']
            lj = seg_j['length']

            return self.MU_0 * self.INV_FOUR_PI * li * lj * cos_angle / r

        def _compute_self_potential(self, node):
            """Compute self-potential coefficient using disk approximation."""
            a = np.sqrt(node['area'] / np.pi)  # Equivalent radius
            if a < 1e-15:
                return 0.0
            return 1.0 / (4 * self.EPS_0 * a)

        def _compute_mutual_potential(self, node_i, node_j):
            """Compute mutual potential coefficient."""
            r = np.linalg.norm(node_j['position'] - node_i['position'])
            if r < 1e-15:
                return 0.0
            return self.INV_FOUR_PI / (self.EPS_0 * r)

        def build(self, include_star=True):
            """Build PEEC matrices."""
            if not self.nodes:
                self._auto_generate_nodes()

            n_loop = len(self.segments)
            n_star = len(self.nodes)

            # L matrix (inductance)
            L = np.zeros((n_loop, n_loop))
            for i in range(n_loop):
                L[i, i] = self._compute_self_inductance(self.segments[i])
                for j in range(i + 1, n_loop):
                    M = self._compute_mutual_inductance(self.segments[i], self.segments[j])
                    L[i, j] = M
                    L[j, i] = M

            # R vector (resistance)
            R = np.array([
                seg['length'] / (seg['sigma'] * seg['width'] * seg['height'])
                for seg in self.segments
            ])

            self._L = L
            self._R = R

            if include_star and n_star > 0:
                # P matrix (potential coefficient)
                P = np.zeros((n_star, n_star))
                for i in range(n_star):
                    P[i, i] = self._compute_self_potential(self.nodes[i])
                    for j in range(i + 1, n_star):
                        Pij = self._compute_mutual_potential(self.nodes[i], self.nodes[j])
                        P[i, j] = Pij
                        P[j, i] = Pij

                # M_LS matrix (Loop-Star coupling)
                M_LS = np.zeros((n_loop, n_star))
                # Simple incidence matrix
                tol = 1e-10
                for i, seg in enumerate(self.segments):
                    p1 = seg['center'] - seg['direction'] * seg['length'] / 2
                    p2 = seg['center'] + seg['direction'] * seg['length'] / 2
                    for j, node in enumerate(self.nodes):
                        if np.linalg.norm(p1 - node['position']) < tol:
                            M_LS[i, j] = -1.0
                        elif np.linalg.norm(p2 - node['position']) < tol:
                            M_LS[i, j] = 1.0

                self._P = P
                self._M_LS = M_LS
                self._matrices_built = True
                return L, R, P, M_LS
            else:
                self._matrices_built = True
                return L, R, None, None

        def compute_impedance(self, freq_hz, port_vector):
            """Compute port impedance at given frequency."""
            if not self._matrices_built:
                self.build(include_star=True)

            port_vector = np.array(port_vector)
            omega = 2 * np.pi * freq_hz

            n_loop = len(self.segments)
            n_star = len(self.nodes)

            # Build Z matrix
            Z_LL = np.diag(self._R) + 1j * omega * self._L
            if self._Zs_diag is not None:
                Z_LL += np.diag(self._Zs_diag)

            if n_star > 0 and self._P is not None:
                Z_SS = self._P / (1j * omega)
                Z_LS = 1j * omega * self._M_LS

                n_total = n_loop + n_star
                Z = np.zeros((n_total, n_total), dtype=complex)
                Z[:n_loop, :n_loop] = Z_LL
                Z[:n_loop, n_loop:] = Z_LS
                Z[n_loop:, :n_loop] = Z_LS.T
                Z[n_loop:, n_loop:] = Z_SS

                # Full port vector
                full_port = np.zeros(n_total, dtype=complex)
                full_port[:n_loop] = port_vector

                # Solve
                I = np.linalg.solve(Z, full_port)
                return np.dot(full_port.conj(), I)
            else:
                I = np.linalg.solve(Z_LL, port_vector)
                return np.dot(port_vector.conj(), I)

        def set_surface_impedance(self, Zs_diag):
            """Set surface impedance for skin effect."""
            self._Zs_diag = np.array(Zs_diag)

        def compute_impedance_with_Zs(self, freq_hz, port_vector):
            """Compute impedance including surface impedance."""
            return self.compute_impedance(freq_hz, port_vector)

        def clear(self):
            """Clear all geometry."""
            self.segments = []
            self.nodes = []
            self._matrices_built = False
            self._L = None
            self._R = None
            self._P = None
            self._M_LS = None
            self._Zs_diag = None

        @property
        def n_loop(self):
            return len(self.segments)

        @property
        def n_star(self):
            return len(self.nodes)

        @property
        def num_segments(self):
            return len(self.segments)

        @property
        def num_nodes(self):
            return len(self.nodes)

    def create_wire_peec(start, end, width, height, n_segments, sigma=5.8e7, include_star=True):
        """Create PEEC matrices for a straight wire."""
        builder = PEECBuilder()
        builder.create_wire(start, end, width, height, n_segments, sigma)
        return builder.build(include_star)

    def create_loop_peec(center, radius, normal, width, height, n_segments, sigma=5.8e7, include_star=True):
        """Create PEEC matrices for a circular loop."""
        builder = PEECBuilder()
        builder.create_loop(center, radius, normal, width, height, n_segments, sigma)
        return builder.build(include_star)


# Export public API
__all__ = ['PEECBuilder', 'create_wire_peec', 'create_loop_peec']
