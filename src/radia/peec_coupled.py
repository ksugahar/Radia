"""
peec_coupled.py

Coupled PEEC + MMM Solver

Computes port impedance Z(f) for conductors near magnetic materials.
Uses Biot-Savart for H-field from conductors, Radia Solve + A-field
for magnetic material response, yielding a coupling matrix Delta_L.

Effective impedance:
    Z_eff(f) = diag(R + Zs(f)) + jw * (L_air + Delta_L)

For linear materials, Delta_L is computed ONCE and reused for all frequencies.

Usage:
    import radia as rad
    from peec_matrices import PEECBuilder
    from peec_coupled import CoupledPEECSolver

    # Build conductor topology
    builder = PEECBuilder()
    n1 = builder.add_node_at(0, 0, 0)
    n2 = builder.add_node_at(0.1, 0, 0)
    builder.add_connected_segment(n1, n2, 1e-3, 1e-3)
    builder.add_port(n1, n2)
    topo = builder.build_topology()

    # Create magnetic object
    rad.FldUnits('m')
    verts = [[-0.02,-0.02,-0.02],[0.02,-0.02,-0.02],[0.02,0.02,-0.02],[-0.02,0.02,-0.02],
             [-0.02,-0.02,0.02],[0.02,-0.02,0.02],[0.02,0.02,0.02],[-0.02,0.02,0.02]]
    core = rad.ObjHexahedron(verts, [0, 0, 0])
    mat = rad.MatLin(999)
    rad.MatApl(core, mat)

    # Coupled solver
    solver = CoupledPEECSolver(topo, [core])
    solver.compute_coupling_matrix()
    Z = solver.compute_port_impedance(1e6)

Part of Radia project
"""

import numpy as np
from peec_topology import PEECCircuitSolver

# Physical constants
MU_0 = 4.0 * np.pi * 1e-7  # H/m


def biot_savart_finite_filament(p1, p2, obs_point):
    """
    Compute H-field at obs_point from a finite current filament (p1 -> p2)
    carrying unit current (I = 1 A).

    Formula:
        H = (1 / 4*pi*d) * (cos(alpha1) - cos(alpha2)) * e_perp

    where d is perpendicular distance from obs_point to filament axis,
    alpha1, alpha2 are angles at p1 and p2.

    Args:
        p1: Start point of filament [x, y, z] (meters)
        p2: End point of filament [x, y, z] (meters)
        obs_point: Observation point [x, y, z] (meters)

    Returns:
        H-field vector [Hx, Hy, Hz] in A/m (for I = 1 A)
    """
    p1 = np.asarray(p1, dtype=float)
    p2 = np.asarray(p2, dtype=float)
    r = np.asarray(obs_point, dtype=float)

    # Filament direction vector
    dl = p2 - p1
    L = np.linalg.norm(dl)
    if L < 1e-20:
        return np.zeros(3)

    e_l = dl / L  # unit direction

    # Vector from p1 to observation point
    r1 = r - p1
    r2 = r - p2

    # Perpendicular distance
    cross = np.cross(e_l, r1)
    d = np.linalg.norm(cross)

    if d < 1e-20:
        # Point on filament axis - H is zero (or singular)
        return np.zeros(3)

    # Unit perpendicular direction (H direction)
    e_perp = cross / d  # This gives e_l x (r - p1) / |e_l x (r - p1)|

    # Cosines of angles
    r1_mag = np.linalg.norm(r1)
    r2_mag = np.linalg.norm(r2)

    if r1_mag < 1e-20 or r2_mag < 1e-20:
        return np.zeros(3)

    cos_alpha1 = np.dot(e_l, r1) / r1_mag
    cos_alpha2 = np.dot(e_l, r2) / r2_mag

    # H = (1 / 4*pi*d) * (cos_alpha1 - cos_alpha2) * e_perp
    # Note: The correct Biot-Savart for finite wire is:
    #   H = I/(4*pi*d) * (cos(alpha1) - cos(alpha2)) * e_phi
    # where alpha1 and alpha2 are angles FROM the wire axis
    H = (1.0 / (4.0 * np.pi * d)) * (cos_alpha1 - cos_alpha2) * e_perp

    return H


class CoupledPEECSolver(PEECCircuitSolver):
    """
    Coupled PEEC + MMM solver for conductors near magnetic materials.

    Extends PEECCircuitSolver with magnetic material coupling via Delta_L matrix.

    The coupling matrix Delta_L[i][j] represents the additional inductance
    on segment i due to unit current in segment j, caused by the magnetized
    magnetic material.

    Physics:
        1. Unit current in segment j produces H-field via Biot-Savart
        2. H-field magnetizes the material (Radia Solve with ObjBckg)
        3. Magnetized material produces vector potential A at other segments
        4. Delta_L[i][j] = dot(A(center_i), dir_i) * length_i

    For linear materials, Delta_L is frequency-independent and computed once.
    """

    def __init__(self, topology_dict, magnetic_objects=None):
        """
        Initialize coupled solver.

        Args:
            topology_dict: Dict from PEECBuilder.build_topology() containing:
                L, R, segment_nodes, n_nodes, ports, n_loop,
                segment_centers, segment_directions, segment_lengths
            magnetic_objects: List of Radia object handles (hex/tet with material)
                             or None for uncoupled mode
        """
        super().__init__(topology_dict)

        # Store segment geometry for Biot-Savart and coupling
        self.seg_centers = np.array(topology_dict['segment_centers'])
        self.seg_directions = np.array(topology_dict['segment_directions'])
        self.seg_lengths = np.array(topology_dict['segment_lengths'])

        # Node positions for Biot-Savart endpoints
        nodes = topology_dict.get('node_positions', None)
        if nodes is not None:
            self.node_positions = np.array(nodes)
        else:
            # Reconstruct from segment center/direction/length
            self._reconstruct_endpoints()

        # Magnetic objects
        self.magnetic_objects = magnetic_objects or []
        self.Delta_L = None  # Coupling matrix (computed on demand)

    def _reconstruct_endpoints(self):
        """Reconstruct segment endpoints from center/direction/length."""
        n_seg = self.n_loop
        self.seg_p1 = np.zeros((n_seg, 3))
        self.seg_p2 = np.zeros((n_seg, 3))

        for i in range(n_seg):
            half = self.seg_directions[i] * self.seg_lengths[i] * 0.5
            self.seg_p1[i] = self.seg_centers[i] - half
            self.seg_p2[i] = self.seg_centers[i] + half

    def compute_coupling_matrix(self, solver_method=0, solver_prec=0.0001,
                                solver_maxiter=1000):
        """
        Compute Delta_L matrix (n_seg x n_seg) from magnetic material response.

        For each column j:
        1. Compute H-field from unit current in segment j (Biot-Savart)
        2. Set as background field via rad.ObjBckg()
        3. Solve for induced magnetization M via rad.Solve()
        4. Compute A-field at each segment center
        5. Delta_L[i][j] = dot(A(center_i), dir_i) * length_i

        Args:
            solver_method: Radia solver method (0=LU, 1=BiCGSTAB, 2=HACApK)
            solver_prec: Solver precision
            solver_maxiter: Maximum iterations

        Returns:
            Delta_L matrix (n_seg x n_seg) in Henries
        """
        try:
            import sys
            import os
            sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
            from _radia_pybind import (ObjBckg, ObjCnt, Solve, Fld, UtiDel,
                                       FldUnits)
        except ImportError:
            import radia as rad
            ObjBckg = rad.ObjBckg
            ObjCnt = rad.ObjCnt
            Solve = rad.Solve
            Fld = rad.Fld
            UtiDel = rad.UtiDel

        if len(self.magnetic_objects) == 0:
            self.Delta_L = np.zeros((self.n_loop, self.n_loop))
            return self.Delta_L

        n_seg = self.n_loop
        Delta_L = np.zeros((n_seg, n_seg))

        # Create container of magnetic objects
        mag_container = ObjCnt(self.magnetic_objects)

        for j in range(n_seg):
            p1_j = self.seg_p1[j]
            p2_j = self.seg_p2[j]

            # H-field function from unit current in segment j
            def h_field_from_seg_j(point, _p1=p1_j, _p2=p2_j):
                H = biot_savart_finite_filament(_p1, _p2, point)
                # Convert H to B for ObjBckg (which expects B in Tesla)
                B = MU_0 * H
                return [B[0], B[1], B[2]]

            # Create background field
            bkg = ObjBckg(h_field_from_seg_j)

            # Create container with magnetic objects + background
            container = ObjCnt([mag_container, bkg])

            # Solve for induced magnetization
            Solve(container, solver_prec, solver_maxiter, solver_method)

            # Compute A-field at each segment center from magnetized material
            for i in range(n_seg):
                center_i = self.seg_centers[i].tolist()
                A_vec = Fld(mag_container, 'a', center_i)
                A_vec = np.array(A_vec)

                # Delta_L[i][j] = dot(A, dir_i) * length_i
                Delta_L[i, j] = np.dot(A_vec, self.seg_directions[i]) * self.seg_lengths[i]

            # Clean up background field
            UtiDel(bkg)

        self.Delta_L = Delta_L
        return Delta_L

    def compute_port_impedance(self, freq, Zs=None):
        """
        Compute port impedance at a given frequency with magnetic coupling.

        Uses L_total = L_air + Delta_L for the inductance matrix.

        Args:
            freq: Frequency in Hz
            Zs: Surface impedance array (n_loop,) complex, or None for DC

        Returns:
            Complex port impedance [Ohm]
        """
        # Temporarily replace L with L_total = L_air + Delta_L
        L_original = self.L
        if self.Delta_L is not None:
            self.L = L_original + self.Delta_L

        Z = super().compute_port_impedance(freq, Zs)

        # Restore original L
        self.L = L_original
        return Z

    def frequency_sweep(self, freqs, Zs_func=None):
        """
        Compute port impedance over a frequency range with magnetic coupling.

        For linear materials, Delta_L is frequency-independent. Call
        compute_coupling_matrix() once before frequency_sweep().

        Args:
            freqs: Array of frequencies [Hz]
            Zs_func: Optional callable(freq) returning Zs array (n_loop,)

        Returns:
            Z_port: Complex array of port impedances, shape (len(freqs),)
        """
        freqs = np.asarray(freqs)
        Z_port = np.zeros(len(freqs), dtype=complex)

        for i, f in enumerate(freqs):
            Zs = Zs_func(f) if Zs_func is not None else None
            Z_port[i] = self.compute_port_impedance(f, Zs)

        return Z_port

    def get_effective_inductance(self):
        """
        Get effective inductance matrix L_total = L_air + Delta_L.

        Returns:
            L_total: (n_seg x n_seg) inductance matrix [H]
        """
        if self.Delta_L is not None:
            return self.L + self.Delta_L
        return self.L.copy()
