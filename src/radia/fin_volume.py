"""Unreduced volume-lattice PEEC of a swept fin conductor (experimental).

Each spine station of :func:`radia.fin_sweep.fin_graph_from_step` is filled
with a lattice of ``n_layers`` concentric rings (the surface lanes scaled
toward the section centroid by ``rho_j = 1 - (j / n_layers)**exponent``,
graded to the surface) plus one centre node.  Every station shares the
same (layer, lane) topology, so branches are

* longitudinal: (s, j, k) -> (s + 1, j, k), and centre -> centre,
* circumferential: (s, j, k) -> (s, j, k + 1) at interior stations,
* radial: (s, j, k) -> (s, j + 1, k), innermost ring -> centre.

All nodes of the first and last station are one terminal node each.
Longitudinal bars tile the section by the lattice dual areas; transverse
bars take the axial dual as one side.  The coil carries DC bar resistance
plus partial inductance only: skin and proximity effect are left to the
inductive coupling, so this route is independent of the surface-impedance
sheet of :mod:`radia.fin_surface_solver`.

The rectangular partial-element kernel orients a bar's width itself; that
approximation is shared with the surface route and is not validated for
non-parallel bars.  No reduction is applied here.
"""
from __future__ import annotations

import math
import numbers
import time
from types import SimpleNamespace
from dataclasses import dataclass, field

import numpy as np

LONGITUDINAL, CIRCUMFERENTIAL, RADIAL = 0, 1, 2
BOX_ORDER = 16  # fin_beak 24x6 s25 vs order 24: R -0.014 %, L +0.043 %, P_wp -0.001 %, 4.2x faster L; order 12 fails the energy gate


@dataclass(frozen=True)
class VolumeLattice:
    nodes: np.ndarray        # (N, 3) m; terminal nodes hold the station centroid
    branches: np.ndarray     # (M, 2) node indices, current from [:, 0] to [:, 1]
    branch_xyz: np.ndarray   # (M, 2, 3) m, the bar's own end points
    kind: np.ndarray         # (M,) LONGITUDINAL / CIRCUMFERENTIAL / RADIAL
    station: np.ndarray      # (M,) owning station (longitudinal: lower station)
    width: np.ndarray        # (M,) m
    height: np.ndarray       # (M,) m
    width_dir: np.ndarray    # (M, 3) unit cross-section width axis; zero = kernel frame
    port_plus: int
    port_minus: int
    rho: np.ndarray          # (n_layers,) ring scale factors, rho[0] = 1
    gates: dict = field(default_factory=dict)


def lattice_branch_count(n_stations, n_lanes, n_layers):
    ring = int(n_lanes) * int(n_layers)
    return (int(n_stations) - 1) * (ring + 1) + 2 * (int(n_stations) - 2) * ring


def _polygon_area(uv):
    x, y = uv[..., 0], uv[..., 1]
    return 0.5 * (np.sum(x * np.roll(y, -1, axis=-1), axis=-1)
                  - np.sum(np.roll(x, -1, axis=-1) * y, axis=-1))


def _inside(poly, pts):
    """Even-odd test of ``pts`` (n, 2) against the closed polygon (m, 2)."""
    x, y = pts[:, 0][:, None], pts[:, 1][:, None]
    x0, y0 = poly[:, 0][None, :], poly[:, 1][None, :]
    x1, y1 = np.roll(poly[:, 0], -1)[None, :], np.roll(poly[:, 1], -1)[None, :]
    crosses = (y0 > y) != (y1 > y)
    with np.errstate(divide="ignore", invalid="ignore"):
        x_at = x0 + (y - y0) * (x1 - x0) / (y1 - y0)
    return np.count_nonzero(crosses & (x < x_at), axis=1) % 2 == 1


def build_volume_lattice(sweep, n_layers, *, exponent=2.0):
    """Fill every station of ``sweep`` (a FinSweepGraph) with the lattice.

    Raises when a lattice cell is inverted or a point leaves its section
    (the radial map needs a section that is star-shaped about its centroid).
    """
    n_layers = int(n_layers)
    if n_layers < 1:
        raise ValueError("n_layers must be >= 1")
    stations, feats = sweep.stations, sweep.features
    S, K, J = len(stations), int(sweep.graph.n_lanes), n_layers
    if S < 3:
        raise ValueError("the volume lattice needs >= 3 stations")
    rho = 1.0 - (np.arange(J) / J) ** float(exponent)

    # uv lattice per station: (S, J, K, 2), centre at the uv origin.
    lane_uv = np.array([f.lane_uv for f in feats], dtype=float)
    if lane_uv.shape != (S, K, 2):
        raise ValueError("station lanes do not match the sweep graph")
    uv = rho[None, :, None, None] * lane_uv[:, None, :, :]
    xyz = np.empty((S, J, K, 3))
    for s, st in enumerate(stations):
        xyz[s] = st.to_xyz(uv[s].reshape(-1, 2)).reshape(J, K, 3)
    centre = np.array([st.centroid for st in stations], dtype=float)

    # Signed cells and dual areas in each section plane.
    dual = np.zeros((S, J, K))
    dual_c = np.zeros(S)
    min_cell = math.inf
    for s in range(S):
        q = uv[s]
        for j in range(J - 1):
            for k in range(K):
                k1 = (k + 1) % K
                a = float(_polygon_area(np.array([q[j, k], q[j, k1], q[j + 1, k1], q[j + 1, k]])))
                min_cell = min(min_cell, a)
                for jj, kk in ((j, k), (j, k1), (j + 1, k1), (j + 1, k)):
                    dual[s, jj, kk] += a / 4
        for k in range(K):
            k1 = (k + 1) % K
            a = float(_polygon_area(np.array([q[J - 1, k], q[J - 1, k1], [0.0, 0.0]])))
            min_cell = min(min_cell, a)
            dual[s, J - 1, k] += a / 3
            dual[s, J - 1, k1] += a / 3
            dual_c[s] += a / 3
    if not min_cell > 0:
        raise ValueError(f"inverted lattice cell (signed area {min_cell:.3e} m^2); "
                         "the section is not star-shaped about its centroid")
    lattice_area = dual.sum(axis=(1, 2)) + dual_c
    outline_area = np.array([st.area for st in stations])

    outside = 0
    for s, st in enumerate(stations):
        if J > 1:
            outside += int(np.count_nonzero(~_inside(st.uv, uv[s, 1:].reshape(-1, 2))))
        outside += int(not _inside(st.uv, np.zeros((1, 2)))[0])
    if outside:
        raise ValueError(f"{outside} interior lattice points lie outside their section")

    # Node numbering: terminals first, then interior stations.
    port_plus, port_minus = 0, 1
    nid = np.full((S, J * K + 1), -1, dtype=np.int64)
    nid[0, :], nid[-1, :] = port_plus, port_minus
    nid[1:-1] = 2 + np.arange((S - 2) * (J * K + 1)).reshape(S - 2, J * K + 1)
    n_nodes = 2 + (S - 2) * (J * K + 1)
    nodes = np.empty((n_nodes, 3))
    nodes[port_plus], nodes[port_minus] = centre[0], centre[-1]
    for s in range(1, S - 1):
        nodes[nid[s, :J * K]] = xyz[s].reshape(-1, 3)
        nodes[nid[s, J * K]] = centre[s]

    def pt(s, j, k):
        return centre[s] if j == J else xyz[s, j, k]

    def col(j, k):
        return J * K if j == J else j * K + k

    # Duals: circumferential half-spans, radial half-spans, axial half-spans.
    circ = np.zeros((S, J, K))
    for j in range(J):
        nxt = np.linalg.norm(np.roll(xyz[:, j], -1, axis=1) - xyz[:, j], axis=2)
        circ[:, j] = 0.5 * (nxt + np.roll(nxt, 1, axis=1))
    ring_pts = np.concatenate([xyz, np.broadcast_to(centre[:, None, None, :], (S, 1, K, 3))], axis=1)
    step = np.linalg.norm(np.diff(ring_pts, axis=1), axis=3)  # (S, J, K) j -> j+1
    rad = np.zeros((S, J, K))
    rad[:, 0] = 0.5 * step[:, 0]
    for j in range(1, J):
        rad[:, j] = 0.5 * (step[:, j - 1] + step[:, j])
    seg = np.linalg.norm(np.diff(xyz, axis=0), axis=3)          # (S-1, J, K)
    ax = np.zeros((S, J, K))
    ax[1:-1] = 0.5 * (seg[:-1] + seg[1:])
    seg_c = np.linalg.norm(np.diff(centre, axis=0), axis=1)
    ax_c = np.zeros(S)
    ax_c[1:-1] = 0.5 * (seg_c[:-1] + seg_c[1:])

    # In-plane ring tangent at every lattice point (longitudinal bar width axis).
    ring_t = np.roll(xyz, -1, axis=2) - np.roll(xyz, 1, axis=2)
    ring_t /= np.linalg.norm(ring_t, axis=3, keepdims=True)
    axial = np.array([st.tangent for st in stations], dtype=float)
    axial /= np.linalg.norm(axial, axis=1, keepdims=True)
    none = np.zeros(3)

    rows = []  # (n0, n1, p0, p1, kind, station, width, height, width_dir)
    for s in range(S - 1):
        for j in range(J):
            for k in range(K):
                area = 0.5 * (dual[s, j, k] + dual[s + 1, j, k])
                w = 0.5 * (circ[s, j, k] + circ[s + 1, j, k])
                rows.append((nid[s, col(j, k)], nid[s + 1, col(j, k)], pt(s, j, k),
                             pt(s + 1, j, k), LONGITUDINAL, s, w, area / w, ring_t[s, j, k]))
        area = 0.5 * (dual_c[s] + dual_c[s + 1])
        rows.append((nid[s, J * K], nid[s + 1, J * K], centre[s], centre[s + 1],
                     LONGITUDINAL, s, math.sqrt(area), math.sqrt(area), none))
    for s in range(1, S - 1):
        for j in range(J):
            for k in range(K):
                k1 = (k + 1) % K
                rows.append((nid[s, col(j, k)], nid[s, col(j, k1)], pt(s, j, k),
                             pt(s, j, k1), CIRCUMFERENTIAL, s,
                             0.5 * (ax[s, j, k] + ax[s, j, k1]),
                             0.5 * (rad[s, j, k] + rad[s, j, k1]), axial[s]))
                c_out = circ[s, j, k]
                c_in = circ[s, j + 1, k] if j + 1 < J else c_out
                a_in = ax[s, j + 1, k] if j + 1 < J else ax_c[s]
                rows.append((nid[s, col(j, k)], nid[s, col(j + 1, k)], pt(s, j, k),
                             pt(s, j + 1, k), RADIAL, s,
                             0.5 * (ax[s, j, k] + a_in), 0.5 * (c_out + c_in), axial[s]))
    branches = np.array([[r[0], r[1]] for r in rows], dtype=np.int64)
    branch_xyz = np.array([[r[2], r[3]] for r in rows], dtype=float)
    kind = np.array([r[4] for r in rows], dtype=np.int8)
    owner = np.array([r[5] for r in rows], dtype=np.int64)
    width = np.array([r[6] for r in rows], dtype=float)
    height = np.array([r[7] for r in rows], dtype=float)
    width_dir = np.array([r[8] for r in rows], dtype=float)
    if len(rows) != lattice_branch_count(S, K, J):
        raise RuntimeError("volume lattice branch count does not match its formula")

    lengths = np.linalg.norm(branch_xyz[:, 1] - branch_xyz[:, 0], axis=1)
    if not (np.min(lengths) > 0 and np.min(width) > 0 and np.min(height) > 0):
        raise ValueError("zero-length or zero-section volume bar")
    if np.any(branches[:, 0] == branches[:, 1]):
        raise RuntimeError("volume lattice has a self-loop branch")

    import scipy.sparse as sp
    from scipy.sparse.csgraph import connected_components

    m = len(branches)
    adj = sp.coo_matrix((np.ones(m), (branches[:, 0], branches[:, 1])),
                        shape=(n_nodes, n_nodes))
    n_comp = int(connected_components(adj, directed=False)[0])
    if n_comp != 1:
        raise ValueError(f"volume lattice has {n_comp} connected components")
    gates = {
        "n_nodes": int(n_nodes), "n_branches": int(m),
        "branches_by_kind": {name: int(np.count_nonzero(kind == code)) for name, code in
                             (("longitudinal", LONGITUDINAL),
                              ("circumferential", CIRCUMFERENTIAL), ("radial", RADIAL))},
        "connected_components": n_comp,
        "incidence_rank": int(n_nodes - n_comp),
        "loop_dimension": int(m - n_nodes + n_comp),
        "min_signed_cell_area_m2": float(min_cell),
        "points_outside_section": 0,
        "segment_length_m": {"min": float(np.min(lengths)), "median": float(np.median(lengths)),
                             "max": float(np.max(lengths))},
        "lattice_over_outline_area": {"min": float(np.min(lattice_area / outline_area)),
                                      "max": float(np.max(lattice_area / outline_area))},
        "first_layer_depth_m": {"min": float(np.min(step[:, 0])), "max": float(np.max(step[:, 0]))},
        "layer_rho": rho.tolist(),
        "each_branch_owned_once": True,
    }
    return VolumeLattice(nodes, branches, branch_xyz, kind, owner, width, height,
                         width_dir, port_plus, port_minus, rho, gates)



def require_energy_admissible(L):
    """Fail loudly unless the inductance matrix is positive semidefinite.

    The input must be a real, finite, non-empty, square and exactly
    symmetric matrix.  A Cholesky factorisation that succeeds proves
    positive definiteness.  When it fails, the symmetric spectrum decides:
    eigenvalues within the roundoff scale n * eps * max|lambda| are null or
    roundoff modes (admissible, counted); anything below is negative
    magnetic energy and raises with the witness eigenvalue and its dominant
    branches.  No diagonal loading or eigenvalue clipping.
    """
    L = np.asarray(L)
    if L.ndim != 2 or L.shape[0] != L.shape[1] or L.shape[0] == 0:
        raise ValueError(f"inductance matrix must be non-empty and square, got shape {L.shape}")
    if not np.isrealobj(L) or not np.issubdtype(L.dtype, np.number):
        raise ValueError(f"inductance matrix must be real, got dtype {L.dtype}")
    L = L.astype(float, copy=False)
    bad = ~np.isfinite(L)
    if bad.any():
        i, j = (int(k) for k in np.argwhere(bad)[0])
        raise ValueError(f"inductance matrix has {int(bad.sum())} non-finite entries, "
                         f"first at ({i}, {j}) = {L[i, j]}")
    asym = float(np.abs(L - L.T).max())
    if asym > 0.0:
        raise ValueError(f"inductance matrix is not symmetric (max |L - L^T| = {asym:.3e} H)")
    n = len(L)
    try:
        np.linalg.cholesky(L)
        return {"cholesky": True, "n": n, "min_eig_over_max": None, "n_null": 0}
    except np.linalg.LinAlgError:
        pass
    ev, vec = np.linalg.eigh(L)
    scale = float(max(abs(ev[0]), abs(ev[-1])))
    tol = n * np.finfo(float).eps * scale
    if ev[0] < -tol:
        top = np.argsort(-np.abs(vec[:, 0]))[:6].tolist()
        raise ValueError(
            f"inductance matrix has negative energy: min eigenvalue {ev[0]:.6e} H "
            f"({ev[0] / scale:.3e} of the spectral scale; roundoff scale {tol:.3e} H), "
            f"{int(np.sum(ev < -tol))} negative modes; witness branches {top}")
    return {"cholesky": False, "n": n,
            "min_eig_over_max": float(ev[0] / scale) if scale > 0.0 else None,
            "n_null": int(np.sum(np.abs(ev) <= tol))}


def energy_receipt(L, branches, current):
    """Source-bound record of the admitted L and the solved current's energy.

    Hashes bind the receipt to the exact matrix and incidence; the
    magnetic energies of the real and imaginary parts of the peak phasor
    (1/2 x^T L x each) are those of the actual fin current mode, and their
    sum is the peak-phasor stored energy.
    """
    import hashlib

    L = np.ascontiguousarray(np.asarray(L, dtype=float))
    br = np.ascontiguousarray(np.asarray(branches, dtype=np.int64))
    cur = np.asarray(current, dtype=complex)
    e_re = 0.5 * float(cur.real @ L @ cur.real)
    e_im = 0.5 * float(cur.imag @ L @ cur.imag)
    return {"L_sha256": hashlib.sha256(L.tobytes()).hexdigest(),
            "L_shape": list(L.shape),
            "incidence_sha256": hashlib.sha256(br.tobytes()).hexdigest(),
            "energy_real_part_J": e_re, "energy_imag_part_J": e_im}


def fin_volume_lattice(step_path, *, n_layers, n_lanes, n_stations, n_outline=None,
                       lane_grading=None, tip_lanes=None, exponent=2.0, route="auto",
                       cad_units_per_meter="auto", max_branches=None):
    """STEP -> (sweep, volume lattice) within the dense-kernel branch budget."""
    from radia.fin_sweep import DEFAULT_MAX_BRANCHES, fin_graph_from_step

    if n_lanes is None or n_stations is None:
        raise ValueError("fin-volume needs explicit n_lanes and n_stations")
    budget = DEFAULT_MAX_BRANCHES if max_branches is None else int(max_branches)
    m = lattice_branch_count(n_stations, n_lanes, n_layers)
    if m > budget:
        raise ValueError(f"the volume lattice needs {m} branches ({n_lanes} lanes x "
                         f"{n_layers} layers x {n_stations} stations), over the "
                         f"{budget} the dense partial-element kernel is budgeted for")
    sweep = fin_graph_from_step(
        step_path, n_lanes=n_lanes, n_stations=n_stations, n_outline=n_outline,
        lane_grading=lane_grading, tip_lanes=tip_lanes,
        cad_units_per_meter=cad_units_per_meter, route=route, max_branches=budget)
    actual = int(sweep.meta["n_stations"])
    m = lattice_branch_count(actual, n_lanes, n_layers)
    if m > budget:  # feature-aligned stations can add stations beyond n_stations
        raise ValueError(f"the volume lattice needs {m} branches ({n_lanes} lanes x "
                         f"{n_layers} layers x {actual} stations after feature-aligned "
                         f"placement, {n_stations} requested), over the {budget} the "
                         f"dense partial-element kernel is budgeted for")
    lat = build_volume_lattice(sweep, n_layers, exponent=exponent)
    if len(lat.branch_xyz) != m:
        raise RuntimeError(f"volume lattice has {len(lat.branch_xyz)} branches, "
                           f"lattice_branch_count says {m}")
    return sweep, lat


def assemble_fin_volume(step_path, *, sigma, n_layers, n_lanes, n_stations,
                        n_outline=None, lane_grading=None, tip_lanes=None,
                        exponent=2.0, route="auto", cad_units_per_meter="auto",
                        max_branches=None, box_order=BOX_ORDER):
    """STEP -> volume lattice and its built PEEC branch model (CPU).

    Returns ``(sweep, lattice, builder, topology, energy)``; ``topology``
    carries the lattice incidence (``segment_nodes``, ``n_nodes``) and the
    one terminal port, and ``energy`` is the
    :func:`require_energy_admissible` record of its L.
    """
    from radia.peec_matrices import PEECBuilder

    if sigma <= 0:
        raise ValueError("fin-volume needs sigma > 0")
    sweep, lat = fin_volume_lattice(
        step_path, n_layers=n_layers, n_lanes=n_lanes, n_stations=n_stations,
        n_outline=n_outline, lane_grading=lane_grading, tip_lanes=tip_lanes,
        exponent=exponent, route=route, cad_units_per_meter=cad_units_per_meter,
        max_branches=max_branches)

    builder = PEECBuilder(oriented_box_order=int(box_order))
    for (p0, p1), w, h, wd in zip(lat.branch_xyz, lat.width, lat.height, lat.width_dir):
        n0 = builder.add_node_at(*map(float, p0))
        n1 = builder.add_node_at(*map(float, p1))
        builder.add_connected_segment(n0, n1, float(w), float(h), sigma=float(sigma))
        if np.any(wd):
            builder.set_last_width_dir(*map(float, wd))
    topo = builder.build_topology()
    if int(topo["n_loop"]) != len(lat.branches):
        raise RuntimeError("PEEC kernel dropped a volume bar")
    topo["segment_nodes"] = lat.branches.copy()
    topo["n_nodes"] = len(lat.nodes)
    topo["ports"] = [(lat.port_plus, lat.port_minus, 0)]
    energy = require_energy_admissible(np.asarray(topo["L"]))
    return sweep, lat, builder, topo, energy


def solve_fin_volume(step_path, *, frequency, sigma, current_A, n_layers,
                     n_lanes, n_stations, n_outline=None, lane_grading=None,
                     tip_lanes=None, exponent=2.0, route="auto",
                     cad_units_per_meter="auto", max_branches=None,
                     box_order=BOX_ORDER, backend="cpu", prima_order=None):
    """STEP -> volume-lattice PEEC branch currents (peak phasors).

    ``backend`` selects the route explicitly: ``"cpu"`` builds L with the
    C++ builder and solves on the host; ``"gpu"`` assembles, gates and
    solves on the device (:func:`solve_fin_volume_gpu`). No fallback.

    Every inductance entry, self included, is the Galerkin integral of
    uniform currents over the bars' oriented boxes (kernel option
    ``oriented_box_order`` = ``box_order``; each bar's width axis comes from
    the lattice).  Filament mutuals make the lattice L indefinite: the
    transverse bars are short, wide plates.  Before the solve, L must pass
    :func:`require_energy_admissible`.

    ``prima_order`` (CPU only, opt-in) replaces the dense solve by a PRIMA
    model of that order expanded at ``s0 = 2 pi frequency``
    (:func:`radia.peec_matrices.prima_reduce`).  The branch currents are
    recovered from the projection basis and must pass the same KCL/KVL true
    residual gate as the unreduced solve, or the call raises.  The order is
    the caller's choice: it is accepted only through that gate, never
    selected or raised here.
    """
    from radia.fin_surface_solver import coil_true_residuals
    from radia.peec_topology import PEECCircuitSolver

    if frequency <= 0 or sigma <= 0:
        raise ValueError("fin-volume needs frequency > 0 and sigma > 0")
    if backend == "gpu":
        if prima_order is not None:
            raise ValueError("prima_order is CPU-only; the GPU route has no PRIMA reduction")
        sweep, lat = fin_volume_lattice(
            step_path, n_layers=n_layers, n_lanes=n_lanes, n_stations=n_stations,
            n_outline=n_outline, lane_grading=lane_grading, tip_lanes=tip_lanes,
            exponent=exponent, route=route, cad_units_per_meter=cad_units_per_meter,
            max_branches=max_branches)
        out = solve_fin_volume_gpu(lat, frequency=frequency, sigma=sigma,
                                   current_A=current_A, box_order=box_order)
        out["sweep"] = sweep
        return out
    if backend != "cpu":
        raise ValueError(f"backend must be 'cpu' or 'gpu', got {backend!r}")
    if prima_order is not None and (isinstance(prima_order, bool)
                                    or not isinstance(prima_order, numbers.Integral)
                                    or prima_order < 1):
        raise ValueError(f"prima_order must be a positive integer, got {prima_order!r}")
    sweep, lat, _builder, topo, energy = assemble_fin_volume(
        step_path, sigma=sigma, n_layers=n_layers, n_lanes=n_lanes,
        n_stations=n_stations, n_outline=n_outline, lane_grading=lane_grading,
        tip_lanes=tip_lanes, exponent=exponent, route=route,
        cad_units_per_meter=cad_units_per_meter, max_branches=max_branches,
        box_order=box_order)
    prima = None
    if prima_order is None:
        peec = PEECCircuitSolver(topo)
        r_branch = np.asarray(np.diag(peec.R_dc) if peec.R_dc.ndim == 2 else peec.R_dc,
                              dtype=float)
        current = np.asarray(peec.compute_branch_currents(frequency, [float(current_A)]),
                             dtype=complex)
    else:
        current, prima = _prima_branch_currents(topo, int(prima_order), frequency,
                                                float(current_A))
        r_branch = np.asarray(topo["R"], dtype=float)
        peec = SimpleNamespace(segment_nodes=np.asarray(lat.branches), n_nodes=len(lat.nodes),
                               ports=topo["ports"], L=np.asarray(topo["L"]))
    residuals = coil_true_residuals(peec, r_branch.astype(complex), current,
                                    frequency, float(current_A))
    i2 = abs(current_A) ** 2
    r_total = float(np.real(np.vdot(current, r_branch * current)) / i2)
    l_total = float(np.real(np.vdot(current, np.asarray(peec.L) @ current)) / i2)
    delta = math.sqrt(2 / (2 * math.pi * frequency * 4e-7 * math.pi * sigma))
    paths = [[(tuple(map(float, p0)), tuple(map(float, p1)))] for p0, p1 in lat.branch_xyz]
    return {"sweep": sweep, "lattice": lat, "current": current,
            "R_ohm": r_total, "L_H": l_total, "coil_true_residual": residuals,
            "box_order": int(box_order), "inductance_energy": energy,
            "energy_receipt": energy_receipt(topo["L"], lat.branches, current),
            "coil_power_W": 0.5 * r_total * i2, "skin_depth_m": delta,
            "n_branches": int(len(lat.branches)), "n_nodes": int(len(lat.nodes)),
            "paths": paths, "I_fil": current, "backend": "cpu", "prima": prima}


def _prima_branch_currents(topo, order, frequency, current_A):
    """Branch currents of the PRIMA model expanded at the operating frequency."""
    from radia.peec_matrices import prima_reduce

    s0 = 2 * math.pi * float(frequency)
    plus, minus = (int(k) for k in topo["ports"][0][:2])
    branches = np.ascontiguousarray(np.asarray(topo["segment_nodes"], dtype=np.int32))
    t0 = time.perf_counter()
    red = prima_reduce(np.asarray(topo["L"], dtype=float), np.asarray(topo["R"], dtype=float),
                       branches, int(topo["n_nodes"]), [plus], [minus], order, s0)
    t_reduce = time.perf_counter() - t0
    G, E, B, V = (np.asarray(red[k], dtype=float) for k in ("G", "E", "B", "V"))
    t0 = time.perf_counter()
    xq = np.linalg.solve(G + 2j * math.pi * float(frequency) * E, B[:, 0] * current_A)
    current = V[:len(branches)] @ xq
    t_eval = time.perf_counter() - t0
    return current, {"order": order, "q": int(red["q"]), "s0": s0,
                     "deflated": list(red["deflated"]),
                     "max_solve_residual": float(red["max_solve_residual"]),
                     "min_E_eig_over_max": float(red["min_E_eig_over_max"]),
                     "t_reduce_s": t_reduce, "t_eval_s": t_eval}


def bar_geometry(lattice):
    """Bar centre, unit direction, length and width axis exactly as the C++
    builder forms them from the two end points (AddConnectedSegment and
    SetLastSegmentWidthDir), so CPU and GPU routes see identical bars."""
    p = np.asarray(lattice.branch_xyz, dtype=float)
    p1, p2 = p[:, 0], p[:, 1]
    center = 0.5 * (p1 + p2)
    d = p2 - p1
    length = np.sqrt(d[:, 0] * d[:, 0] + d[:, 1] * d[:, 1] + d[:, 2] * d[:, 2])
    if np.any(length < 1e-15):
        raise ValueError("volume lattice has a zero-length bar")
    d = d / length[:, None]
    wd = np.asarray(lattice.width_dir, dtype=float)
    has = np.any(wd != 0.0, axis=1)
    dot = wd[:, 0] * d[:, 0] + wd[:, 1] * d[:, 1] + wd[:, 2] * d[:, 2]
    w = wd - dot[:, None] * d
    n = np.sqrt(w[:, 0] * w[:, 0] + w[:, 1] * w[:, 1] + w[:, 2] * w[:, 2])
    if np.any(has & ~(n > 1e-9)):
        raise ValueError("width axis is parallel to a bar")
    w = np.where(has[:, None], w / np.where(n > 0, n, 1.0)[:, None], 0.0)
    return center, d, length, w, has


def solve_fin_volume_gpu(lattice, *, frequency, sigma, current_A, box_order=BOX_ORDER):
    """Device-resident solve of a built volume lattice.

    L (oriented-box Galerkin, the same shared kernel as the CPU builder) is
    assembled on the GPU, gated for energy there (Cholesky, else the
    spectrum), Z = R + j omega L is factorised there and the port solve runs
    there; only bar geometry, R and incidence go up, and only branch
    currents, node voltages and scalars come back. Raises when CUDA, the
    module or the solve fails; there is no CPU fallback.
    """
    from radia.gpu import _native

    if frequency <= 0 or sigma <= 0:
        raise ValueError("fin-volume needs frequency > 0 and sigma > 0")
    mod = _native("_peec_fin_cuda")
    center, d, length, wd, has = bar_geometry(lattice)
    width = np.asarray(lattice.width, dtype=float)
    height = np.asarray(lattice.height, dtype=float)
    resistance = length / (float(sigma) * (width * height))
    nodes = np.ascontiguousarray(lattice.branches, dtype=np.int32)
    n_nodes = len(lattice.nodes)
    plus, minus = int(lattice.port_plus), int(lattice.port_minus)
    if plus == minus:
        raise ValueError("port terminals coincide")
    parent = np.arange(n_nodes)

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a
    for a, b in nodes:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb
    roots = {find(k) for k in range(n_nodes)}
    if len(roots) != 1:
        raise ValueError(f"volume lattice has {len(roots)} disconnected parts")
    node_map = np.full(n_nodes, -1, dtype=np.int32)
    keep = np.flatnonzero(np.arange(n_nodes) != minus)  # the CPU solver grounds the negative terminal
    node_map[keep] = np.arange(len(keep), dtype=np.int32)
    injection = np.zeros(len(keep), dtype=complex)
    injection[node_map[plus]] = complex(current_A)
    current, voltage, meta = mod.solve(
        np.ascontiguousarray(center.ravel()), np.ascontiguousarray(d.ravel()),
        np.ascontiguousarray(wd.ravel()), np.ascontiguousarray(has, dtype=np.uint8),
        np.ascontiguousarray(length), np.ascontiguousarray(width), np.ascontiguousarray(height),
        np.ascontiguousarray(resistance), int(box_order), float(frequency),
        np.ascontiguousarray(nodes.ravel()), node_map, int(len(keep)), injection)
    i2 = abs(current_A) ** 2
    r_total = float(meta["IRI"]) / i2
    l_total = float(meta["ILI"]) / i2
    delta = math.sqrt(2 / (2 * math.pi * frequency * 4e-7 * math.pi * sigma))
    paths = [[(tuple(map(float, p0)), tuple(map(float, p1)))] for p0, p1 in lattice.branch_xyz]
    return {"lattice": lattice, "current": np.asarray(current), "R_ohm": r_total, "L_H": l_total,
            "port_voltage_V": complex(voltage[node_map[plus]]), "gpu": dict(meta),
            "box_order": int(box_order), "coil_power_W": 0.5 * r_total * i2,
            "true_residual": {"kvl_relative": float(meta["true_relative_residual"]),
                              "kcl_relative": float(meta["kcl_relative_residual"]),
                              "limit": 1e-6},
            "skin_depth_m": delta, "n_branches": int(len(lattice.branches)),
            "n_nodes": int(n_nodes), "paths": paths, "I_fil": np.asarray(current),
            "backend": "gpu"}
