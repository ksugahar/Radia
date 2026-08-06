"""
peec_shield.py

Shield conductor modeling for PEEC analysis.

Adds aluminum (or other conductor) shield as PEEC filament grid to
PEECBuilder. Shield currents are induced by magnetic coupling with
coil segments - no external port excitation.

The shield effect on coil impedance appears as reflected impedance:
    Delta_Z_shield(f) = Z_port(with_shield) - Z_port(without_shield)
    Delta_R_shield = Re(Delta_Z_shield)
    Delta_L_shield = Im(Delta_Z_shield) / omega

Physics:
    The shield is modeled as a grid of PEEC filaments forming closed
    current loops. Each filament has:
    - R = rho * length / (dual_width * thickness)
    - L = self + mutual inductance (Neumann integral with all segments)
    - Mutual coupling with coil segments (also via Neumann integral)

    The MNA solver treats the shield as a separate connected component
    with no port excitation. Shield currents are determined by the
    coupling through the full mutual inductance matrix.

Usage:
    from radia.peec_matrices import PEECBuilder
    from radia.peec_shield import add_shield_mesh, add_coil_from_segments
    from radia.peec_topology import PEECCircuitSolver

    builder = PEECBuilder()

    # Add coil
    n_start, n_end = add_coil_from_segments(
        builder, p1_m, p2_m, wire_w, wire_h, sigma_cu)
    builder.add_port(n_start, n_end)

    # Add shield (no port - floating conductor)
    add_shield_mesh(builder, [0, 0, -0.005], 0.05, 0.05,
                    thickness=1e-3, sigma=3.5e7, nx=8, ny=8)

    # Solve
    topo = builder.build_topology()
    solver = PEECCircuitSolver(topo)
    Z = solver.compute_port_impedance(freq)

Part of Radia project
"""

import numpy as np
from collections import Counter


# Aluminum conductivity [S/m]
SIGMA_AL = 3.5e7


def add_shield_mesh(builder, center, size_x, size_y, thickness, sigma,
                    nx, ny):
    """Add a rectangular conducting shield as PEEC filament grid.

    Creates a grid of (nx+1)*(ny+1) nodes connected by
    nx*(ny+1) horizontal + (nx+1)*ny vertical segments.
    The shield forms closed loops - no port is added.

    Cross-section of each segment:
        - height = shield thickness
        - width = dual mesh width (full for interior, half for boundary)

    Args:
        builder: PEECBuilder instance (peec_matrices.PEECBuilder)
        center: [x, y, z] center of shield [m]
        size_x: Shield extent in x-direction [m]
        size_y: Shield extent in y-direction [m]
        thickness: Shield thickness [m]
        sigma: Electrical conductivity [S/m] (e.g. 3.5e7 for aluminum)
        nx: Number of divisions in x
        ny: Number of divisions in y

    Returns:
        dict with:
            n_nodes: number of shield nodes added
            n_segments: number of shield segments added
            node_ids: 2D array (ny+1, nx+1) of PEECBuilder node IDs
    """
    cx, cy, cz = center[0], center[1], center[2]
    dx = size_x / nx
    dy = size_y / ny

    # Create grid of nodes
    node_ids = np.zeros((ny + 1, nx + 1), dtype=int)
    for iy in range(ny + 1):
        for ix in range(nx + 1):
            x = cx - size_x / 2 + ix * dx
            y = cy - size_y / 2 + iy * dy
            z = cz
            node_ids[iy, ix] = builder.add_node_at(x, y, z)

    n_segments = 0

    # Horizontal segments (along x-direction)
    for iy in range(ny + 1):
        # Dual mesh width for this row
        if iy == 0 or iy == ny:
            w_dual = dy / 2  # boundary: half width
        else:
            w_dual = dy  # interior: full width

        for ix in range(nx):
            n_from = int(node_ids[iy, ix])
            n_to = int(node_ids[iy, ix + 1])
            builder.add_connected_segment(n_from, n_to, w_dual, thickness,
                                          sigma)
            n_segments += 1

    # Vertical segments (along y-direction)
    for ix in range(nx + 1):
        # Dual mesh width for this column
        if ix == 0 or ix == nx:
            w_dual = dx / 2  # boundary: half width
        else:
            w_dual = dx  # interior: full width

        for iy in range(ny):
            n_from = int(node_ids[iy, ix])
            n_to = int(node_ids[iy + 1, ix])
            builder.add_connected_segment(n_from, n_to, w_dual, thickness,
                                          sigma)
            n_segments += 1

    n_nodes = (nx + 1) * (ny + 1)
    return {
        'n_nodes': n_nodes,
        'n_segments': n_segments,
        'node_ids': node_ids,
    }


def add_coil_from_segments(builder, p1_m, p2_m, wire_width, wire_height,
                           sigma=5.8e7):
    """Add coil segments to PEECBuilder, maintaining connectivity.

    Matches segment endpoints by proximity to find unique nodes,
    then adds nodes and segments to the builder.

    Terminal nodes (appearing in only one segment) are identified
    as port candidates.

    Args:
        builder: PEECBuilder instance
        p1_m: (N, 3) array of segment start points [m]
        p2_m: (N, 3) array of segment end points [m]
        wire_width: Wire cross-section width [m]
        wire_height: Wire cross-section height [m]
        sigma: Conductivity [S/m] (default: 5.8e7 for copper)

    Returns:
        (node_start, node_end): PEECBuilder node IDs for the coil terminals
    """
    p1_m = np.asarray(p1_m)
    p2_m = np.asarray(p2_m)
    n_seg = len(p1_m)

    # Find unique node positions by proximity matching
    tol = 1e-8  # 10 nm tolerance
    unique_points = []

    def find_or_create(pt):
        for idx, up in enumerate(unique_points):
            if np.linalg.norm(pt - up) < tol:
                return idx
        unique_points.append(pt.copy())
        return len(unique_points) - 1

    seg_nf = []
    seg_nt = []
    for i in range(n_seg):
        nf = find_or_create(p1_m[i])
        nt = find_or_create(p2_m[i])
        seg_nf.append(nf)
        seg_nt.append(nt)

    # Find terminal nodes (degree 1 = endpoint of the path)
    degree = Counter()
    for nf, nt in zip(seg_nf, seg_nt):
        degree[nf] += 1
        degree[nt] += 1

    terminals = [n for n, d in degree.items() if d == 1]

    is_closed_loop = len(terminals) < 2

    if is_closed_loop:
        # Closed loop: break it by splitting the last segment's endpoint.
        # Create a duplicate node at the same position but with different ID.
        # The last segment connects to the new node instead of the original.
        # Port spans from original node to new node (across the break).
        break_node = seg_nt[-1]  # node where loop closes
        break_pt = unique_points[break_node].copy()
        new_node_idx = len(unique_points)
        unique_points.append(break_pt)
        seg_nt[-1] = new_node_idx  # redirect last segment to new node
        terminals = [break_node, new_node_idx]

    # Add all unique nodes to PEECBuilder
    local_to_builder = {}
    for idx, pt in enumerate(unique_points):
        local_to_builder[idx] = builder.add_node_at(pt[0], pt[1], pt[2])

    # Add all segments
    for i in range(n_seg):
        nf = local_to_builder[seg_nf[i]]
        nt = local_to_builder[seg_nt[i]]
        builder.add_connected_segment(nf, nt, wire_width, wire_height, sigma)

    return (local_to_builder[terminals[0]], local_to_builder[terminals[1]])


def compute_shield_effect(coil_p1_m, coil_p2_m, wire_width, wire_height,
                          shield_center, shield_size_x, shield_size_y,
                          shield_thickness, shield_sigma, shield_nx, shield_ny,
                          freqs, coil_sigma=5.8e7):
    """Compute shield effect on coil impedance over a frequency range.

    Builds two PEEC models (with and without shield), solves both,
    and returns the impedance difference.

    Args:
        coil_p1_m, coil_p2_m: (N, 3) coil segment endpoints [m]
        wire_width, wire_height: Coil wire cross-section [m]
        shield_center: [x, y, z] center of shield [m]
        shield_size_x, shield_size_y: Shield dimensions [m]
        shield_thickness: Shield thickness [m]
        shield_sigma: Shield conductivity [S/m]
        shield_nx, shield_ny: Shield mesh divisions
        freqs: Array of frequencies [Hz]
        coil_sigma: Coil conductivity [S/m]

    Returns:
        dict with:
            freqs: frequency array [Hz]
            Z_with: Z_port with shield (complex array)
            Z_without: Z_port without shield (complex array)
            Delta_R: Re(Delta_Z) [Ohm]
            Delta_L: Im(Delta_Z)/omega [H]
            L_with: Inductance with shield [H]
            L_without: Inductance without shield [H]
            n_coil_segments: number of coil filaments
            n_shield_segments: number of shield filaments
    """
    from radia.peec_matrices import PEECBuilder
    from radia.peec_topology import PEECCircuitSolver

    freqs = np.asarray(freqs, dtype=np.float64)

    # --- Model 1: Coil only ---
    builder1 = PEECBuilder()
    port1 = add_coil_from_segments(builder1, coil_p1_m, coil_p2_m,
                                   wire_width, wire_height, coil_sigma)
    builder1.add_port(port1[0], port1[1])
    topo1 = builder1.build_topology()
    n_coil = topo1['n_loop']

    solver1 = PEECCircuitSolver(topo1)
    Z_without = solver1.frequency_sweep(freqs)

    # --- Model 2: Coil + Shield ---
    builder2 = PEECBuilder()
    port2 = add_coil_from_segments(builder2, coil_p1_m, coil_p2_m,
                                   wire_width, wire_height, coil_sigma)
    builder2.add_port(port2[0], port2[1])
    shield_info = add_shield_mesh(builder2, shield_center,
                                  shield_size_x, shield_size_y,
                                  shield_thickness, shield_sigma,
                                  shield_nx, shield_ny)
    topo2 = builder2.build_topology()

    solver2 = PEECCircuitSolver(topo2)
    Z_with = solver2.frequency_sweep(freqs)

    # Compute differences
    Delta_Z = Z_with - Z_without
    omega = 2 * np.pi * freqs

    L_without = np.where(omega > 0, Z_without.imag / omega, 0.0)
    L_with = np.where(omega > 0, Z_with.imag / omega, 0.0)
    Delta_L = L_with - L_without

    return {
        'freqs': freqs,
        'Z_with': Z_with,
        'Z_without': Z_without,
        'Delta_R': Delta_Z.real,
        'Delta_L': Delta_L,
        'L_with': L_with,
        'L_without': L_without,
        'n_coil_segments': n_coil,
        'n_shield_segments': shield_info['n_segments'],
    }


def compute_shield_effect_2port(
        tx_p1_m, tx_p2_m, tx_wire_w, tx_wire_h,
        rx_p1_m, rx_p2_m, rx_wire_w, rx_wire_h,
        shield_center, shield_size_x, shield_size_y,
        shield_thickness, shield_sigma, shield_nx, shield_ny,
        freq, coil_sigma=5.8e7):
    """Compute shield effect on a 2-port (TX/RX) coil system.

    Builds two PEEC models (with and without shield) for the full
    TX+RX system, solves both, and returns Z-matrix, L, k differences.

    Args:
        tx_p1_m, tx_p2_m: TX coil segment endpoints [m]
        tx_wire_w, tx_wire_h: TX wire cross-section [m]
        rx_p1_m, rx_p2_m: RX coil segment endpoints [m]
        rx_wire_w, rx_wire_h: RX wire cross-section [m]
        shield_center: [x, y, z] center of shield [m]
        shield_size_x, shield_size_y: Shield dimensions [m]
        shield_thickness: Shield thickness [m]
        shield_sigma: Shield conductivity [S/m]
        shield_nx, shield_ny: Shield mesh divisions
        freq: Frequency [Hz]
        coil_sigma: Coil conductivity [S/m]

    Returns:
        dict with Z_matrix, L, k for both with and without shield cases,
        plus Delta values.
    """
    from radia.peec_matrices import PEECBuilder
    from radia.peec_topology import PEECCircuitSolver

    omega = 2 * np.pi * freq

    # --- Model 1: TX + RX only ---
    b1 = PEECBuilder()
    tx1 = add_coil_from_segments(b1, tx_p1_m, tx_p2_m,
                                 tx_wire_w, tx_wire_h, coil_sigma)
    rx1 = add_coil_from_segments(b1, rx_p1_m, rx_p2_m,
                                 rx_wire_w, rx_wire_h, coil_sigma)
    b1.add_port(tx1[0], tx1[1])
    b1.add_port(rx1[0], rx1[1])
    topo1 = b1.build_topology()
    s1 = PEECCircuitSolver(topo1)
    r1 = s1.compute_coupling_coefficient(freq)

    # --- Model 2: TX + RX + Shield ---
    b2 = PEECBuilder()
    tx2 = add_coil_from_segments(b2, tx_p1_m, tx_p2_m,
                                 tx_wire_w, tx_wire_h, coil_sigma)
    rx2 = add_coil_from_segments(b2, rx_p1_m, rx_p2_m,
                                 rx_wire_w, rx_wire_h, coil_sigma)
    b2.add_port(tx2[0], tx2[1])
    b2.add_port(rx2[0], rx2[1])
    shield_info = add_shield_mesh(b2, shield_center,
                                  shield_size_x, shield_size_y,
                                  shield_thickness, shield_sigma,
                                  shield_nx, shield_ny)
    topo2 = b2.build_topology()
    s2 = PEECCircuitSolver(topo2)
    r2 = s2.compute_coupling_coefficient(freq)

    return {
        # Without shield
        'Z_no_shield': r1['Z_matrix'],
        'L1_no': r1['L1'], 'L2_no': r1['L2'], 'M_no': r1['M'],
        'k_no': r1['k'],
        'R1_no': r1['Z_matrix'][0, 0].real,
        'R2_no': r1['Z_matrix'][1, 1].real,
        # With shield
        'Z_with_shield': r2['Z_matrix'],
        'L1_with': r2['L1'], 'L2_with': r2['L2'], 'M_with': r2['M'],
        'k_with': r2['k'],
        'R1_with': r2['Z_matrix'][0, 0].real,
        'R2_with': r2['Z_matrix'][1, 1].real,
        # Differences
        'Delta_k': r2['k'] - r1['k'],
        'Delta_L1': r2['L1'] - r1['L1'],
        'Delta_L2': r2['L2'] - r1['L2'],
        'Delta_M': r2['M'] - r1['M'],
        'Delta_R1': r2['Z_matrix'][0, 0].real - r1['Z_matrix'][0, 0].real,
        'Delta_R2': r2['Z_matrix'][1, 1].real - r1['Z_matrix'][1, 1].real,
        # Info
        'n_shield_segments': shield_info['n_segments'],
        'freq': freq,
    }
