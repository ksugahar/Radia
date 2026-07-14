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

from radia.peec_matrices import PEECBuilder
from radia.peec_topology import PEECCircuitSolver


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
    _assert_solver_L_finite(solver, seg_indices, filament_paths)
    return solver, seg_indices, 0, 1


def _assert_solver_L_finite(solver, seg_indices, filament_paths):
    """Fail fast (CLAUDE.md "No Fallbacks - Fail Fast, Fail Loud") when
    the assembled mutual-inductance matrix contains non-finite entries.

    Background (2026-05-16, keiko 1turn_coil_loft_outsideline.step):
    when the open-spine path samples a coil with a sharp corner at
    the lead-cap junction, parallel transport of the cross-section
    perimeter across the corner places adjacent filament segments
    near-coincident.  The Ruehli mutual-inductance kernel is singular
    on coincident segment pairs, and the dense kernel returned NaN /
    Inf silently rather than raising.  Downstream
    ``compute_port_impedance`` then propagated the NaN into a
    ``L_coil_nH = NaN`` JSON output that exited 0 -- the user saw
    a "successful" run with garbage numbers.

    This guard converts the silent NaN into a hard ValueError that
    points at the offending filament/segment pairs and (per CLAUDE.md
    "STEP-Only Centerline: Auto-Detect or Fail") at the documented
    fix path: regenerate the STEP with vertex-aligned cross-sections
    so the lateral surface is a single dominant BSPLINE / TORUS,
    routing the centerline through Predicate 1 (UV-map sampling)
    which does NOT use parallel-transport.
    """
    L = getattr(solver, "L", None)
    if L is None:
        return  # HACApK path: no dense L assembled
    L = np.asarray(L)
    bad = ~np.isfinite(L)
    if not bad.any():
        return  # clean L, normal path

    # Map the N_seg-indexed L matrix back to (filament, segment_in_fil)
    # so the diagnostic names something the user can act on.
    n_seg_total = L.shape[0]
    fil_of_seg = np.full(n_seg_total, -1, dtype=np.int64)
    pos_in_fil = np.full(n_seg_total, -1, dtype=np.int64)
    for k, segs in enumerate(seg_indices):
        for j, s in enumerate(segs):
            if 0 <= s < n_seg_total:
                fil_of_seg[s] = k
                pos_in_fil[s] = j

    rows, cols = np.where(bad)
    n_bad = int(bad.sum())
    diag_bad = int((~np.isfinite(np.diag(L))).sum())
    # Sample first 8 pairs for the diagnostic (avoid spamming on
    # massively degenerate geometries).
    samples = []
    for r, c in zip(rows[:8], cols[:8]):
        fr, sr = int(fil_of_seg[r]), int(pos_in_fil[r])
        fc, sc = int(fil_of_seg[c]), int(pos_in_fil[c])
        samples.append(
            f"  L[{int(r)},{int(c)}] = fil{fr}.seg{sr} <-> fil{fc}.seg{sc}")

    n_fil = len(filament_paths)
    n_per_fil = n_seg_total // max(n_fil, 1)
    raise ValueError(
        "PEEC mutual-inductance matrix L contains "
        f"{n_bad} non-finite entries ({diag_bad} on the diagonal).  "
        "This means the open-spine cross-section sampler placed "
        "filament segments near-coincident across a sharp spine corner "
        "(typically a lead-to-arc junction) and the Ruehli kernel "
        "returned NaN / Inf silently.\n"
        f"Geometry: {n_fil} filaments x {n_per_fil} segments each.\n"
        "Sample degenerate pairs (within-filament, near a corner):\n"
        + "\n".join(samples) +
        "\nFIX: regenerate the STEP with vertex-aligned cross-sections "
        "(e.g. gen_1turn_coil_loft.py vertex_side='inner' or 'outer', "
        "or build123d sweep(circle_profile, smooth_spline_path)) so "
        "the lateral surface is a single dominant BSPLINE / TORUS.  "
        "That routes the centerline through Predicate 1 (UV-map "
        "sampling) which samples directly from the lateral surface "
        "without parallel transport, eliminating the singular pair.")


def build_loop_bundle_impedance(solver, seg_of_filament):
    """Reduce a parallel-filament PEEC bundle to filament-level (R_f, L_f).

    For a bundle of K filaments where each filament is a SERIES CHAIN
    of segments with shared port start/end nodes (the topology produced
    by build_bundle_solver), every segment in filament k carries the same
    current I_f[k]. Under this constraint the nodal MNA system collapses
    exactly to a K x K Loop-form problem:

        R_f[k, k]  = sum_{i in seg_of_filament[k]} R_dc[i, i]
        L_f[k, l]  = sum_{i in seg_of_filament[k]}
                     sum_{j in seg_of_filament[l]} L[i, j]

    This is a lossless reduction (not approximate) for the series-chain
    bundle topology and scales the system from O(N_seg^3) to O(K^3).

    Returns (R_f, L_f) as K x K real arrays.
    """
    K = len(seg_of_filament)
    R_full = solver.R_dc
    L_full = solver.L
    N_seg = R_full.shape[0]
    A = np.zeros((K, N_seg), dtype=np.float64)
    for k, segs in enumerate(seg_of_filament):
        for i in segs:
            A[k, int(i)] = 1.0
    # R_full may be a 1D vector (per-segment R) or 2D diagonal matrix.
    R_diag_vec = R_full if R_full.ndim == 1 else np.diag(R_full)
    R_f = np.diag(A @ R_diag_vec)            # K x K diagonal
    L_f = A @ L_full @ A.T                    # K x K dense symmetric
    return R_f, L_f


def solve_loop_bundle(R_f, L_f, frequency, I_port=1.0, Zs_fil=None, emf=None):
    """Solve the K x K Loop-form bundle system at one frequency.

    Extended system (K+1 x K+1) with port voltage V_port as an unknown:

        [ Z_fil  -1 ] [I_f   ]   [-emf  ]
        [ 1^T     0 ] [V_port] = [I_port]

    Z_fil = R_f + jw * L_f + diag(Zs_fil) (optional per-filament SIBC
    surface impedance, e.g. Dowell).

    ``emf`` (optional, length K complex) is a per-filament externally
    induced EMF (voltage) added to the filament KVL row, so the filament
    equation becomes ``Z_fil I_f = V_port 1 - emf``.  This is the entry
    point for a workpiece back-reaction: emf[k] = jw * flux_wp,k where
    flux_wp,k = integral over filament k of A_wp . dl is the flux linked
    by filament k from the workpiece reaction vector potential (see
    ``peec_coupled_bem_solver.CoupledPEECBEMSolver``).  Default None
    reproduces the uncoupled air solve exactly.

    Returns (I_f, V_port) where I_f is (K,) complex.
    """
    import math
    omega = 2.0 * math.pi * frequency
    K = R_f.shape[0]
    Z_fil = R_f.astype(complex) + 1j * omega * L_f.astype(complex)
    if Zs_fil is not None:
        Z_fil = Z_fil + np.diag(np.asarray(Zs_fil, dtype=complex))
    M = np.zeros((K + 1, K + 1), dtype=complex)
    M[:K, :K] = Z_fil
    M[:K, K] = -1.0
    M[K, :K] = 1.0
    b = np.zeros(K + 1, dtype=complex)
    if emf is not None:
        b[:K] = -np.asarray(emf, dtype=complex)
    b[K] = I_port
    x = np.linalg.solve(M, b)
    return x[:K], complex(x[K])


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
