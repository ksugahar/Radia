"""peec_bundle.py

Filament-bundle PEEC helper for the "N parallel filaments between
a single port+ and port-" workflow.

Given a list of N filament polylines (e.g. from
CoilBuilder.to_filaments()), build a PEECCircuitSolver where every
filament is electrically in parallel between two logical nodes.
Calling solver.compute_branch_currents(freq, [I_total]) then returns
per-filament complex AC currents including full (R + jwM) coupling.

Usage:
    from radia.peec_bundle import build_bundle_solver
    solver, seg_of_fil, port_p, port_m = build_bundle_solver(
        filament_paths, dw, dh, sigma)
    I_branch = solver.compute_branch_currents(freq, [I_total])
    # Per-filament current (series chain, so all segments equal):
    I_k = np.array([np.mean(I_branch[seg_of_fil[k]]) for k in range(N)])

For non-uniform cross-sections (tapers / lofts), pass ``cell_wh``
(a list of per-piece (w, h) tuples matching ``filament_paths``)
instead of scalar ``dw, dh``.

This helper does NOT touch the C++ PEEC layout rules: each filament
becomes its own series chain of PEEC segments (correct physical
geometry for L_kk, R_k, M_kj), and the electrical "bundle" topology
is imposed AFTER build_topology() by remapping segment_nodes so that
all filament starts map to logical node 0 and all filament ends to
logical node 1.
"""

import numpy as np

from peec_matrices import PEECBuilder
from peec_topology import PEECCircuitSolver


def build_bundle_solver(filament_paths, dw, dh, sigma, cell_wh=None):
    """Build a PEECCircuitSolver for N parallel filaments sharing one port.

    Args:
        filament_paths: list of N filaments. filament_paths[k] is a list
            of ((p0, p1), (p0, p1), ...) polyline tuples; consecutive
            tuples share endpoints. Points are 3-tuples in meters.
        dw, dh: scalar sub-filament cross-section (m). Used when cell_wh
            is None (uniform cross-section case).
        sigma: conductor conductivity (S/m).
        cell_wh: optional per-piece cross-section for non-uniform cases
            (lofts, tapers). List of N lists; cell_wh[k][i] is the
            (w_sub, h_sub) tuple for the i-th polyline piece of
            filament k. Takes precedence over scalar (dw, dh). All
            filaments at the same piece index typically share (w_sub,
            h_sub) in a rect loft, but nothing requires that.

    Returns:
        solver: PEECCircuitSolver with a single port, port_id=0,
            between logical nodes (port_plus=0) and (port_minus=1).
        seg_of_filament: list of N lists. seg_of_filament[k] is the
            list of PEEC segment indices belonging to filament k
            (series chain -- they all carry the same current).
        port_plus, port_minus: logical node IDs (always 0, 1).
    """
    builder = PEECBuilder()
    N = len(filament_paths)
    use_cell_wh = cell_wh is not None

    phys_start_ids = []
    phys_end_ids = []
    phys_internal_ids = []
    seg_indices = [[] for _ in range(N)]
    seg_counter = 0

    for k, path in enumerate(filament_paths):
        # Polyline waypoints from consecutive tuples.
        waypoints = [path[0][0]]
        for (_p1, p2) in path:
            waypoints.append(p2)

        node_ids = []
        for (x, y, z) in waypoints:
            nid = builder.add_node_at(float(x), float(y), float(z))
            node_ids.append(nid)

        phys_start_ids.append(node_ids[0])
        phys_end_ids.append(node_ids[-1])
        phys_internal_ids.append(node_ids[1:-1])

        for i in range(len(node_ids) - 1):
            if use_cell_wh:
                w_i, h_i = cell_wh[k][i]
            else:
                w_i, h_i = dw, dh
            builder.add_connected_segment(node_ids[i], node_ids[i + 1],
                                          float(w_i), float(h_i),
                                          sigma=sigma)
            seg_indices[k].append(seg_counter)
            seg_counter += 1

    topo = builder.build_topology()
    seg_nodes = np.array(topo['segment_nodes'], dtype=np.int64)

    start_set = set(phys_start_ids)
    end_set = set(phys_end_ids)

    all_internal = []
    for lst in phys_internal_ids:
        all_internal.extend(lst)
    all_internal = sorted(set(all_internal))

    logical_map = {}
    for pid in start_set:
        logical_map[pid] = 0
    for pid in end_set:
        logical_map[pid] = 1
    for i, pid in enumerate(all_internal):
        logical_map[pid] = 2 + i

    seg_nodes_new = np.empty_like(seg_nodes)
    for f in range(seg_nodes.shape[0]):
        seg_nodes_new[f, 0] = logical_map[int(seg_nodes[f, 0])]
        seg_nodes_new[f, 1] = logical_map[int(seg_nodes[f, 1])]

    topo['segment_nodes'] = seg_nodes_new
    topo['n_nodes'] = 2 + len(all_internal)
    topo['ports'] = [(0, 1, 0)]

    solver = PEECCircuitSolver(topo)
    return solver, seg_indices, 0, 1


def filament_currents(I_branch, seg_of_filament):
    """Aggregate per-filament current from per-segment branch current array.

    In a bundle built by build_bundle_solver, each filament is a series
    chain of PEEC segments with no interior taps, so all segments of one
    filament carry the same current. We take the mean to suppress any
    round-off noise.

    Args:
        I_branch: (n_loop,) complex array from
            solver.compute_branch_currents(freq, port_currents).
        seg_of_filament: second return value of build_bundle_solver.

    Returns:
        I: (N,) complex array. sum(I) equals port_currents[0] up to
           numerical precision.
    """
    N = len(seg_of_filament)
    out = np.zeros(N, dtype=complex)
    for k in range(N):
        out[k] = np.mean(I_branch[seg_of_filament[k]])
    return out
