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
from dataclasses import dataclass, field

import numpy as np

LONGITUDINAL, CIRCUMFERENTIAL, RADIAL = 0, 1, 2
BOX_ORDER = 8


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

    A Cholesky factorisation that succeeds proves positive definiteness.
    When it fails, the symmetric spectrum decides: eigenvalues within the
    roundoff scale n * eps * lambda_max are null or roundoff modes
    (admissible, counted); anything below is negative magnetic energy and
    raises with the witness eigenvalue and its dominant branches.  No
    diagonal loading or eigenvalue clipping.
    """
    L = np.asarray(L, dtype=float)
    asym = float(np.abs(L - L.T).max())
    if asym > 0.0:
        raise ValueError(f"inductance matrix is not symmetric (max |L - L^T| = {asym:.3e} H)")
    try:
        np.linalg.cholesky(L)
        return {"cholesky": True, "min_eig_over_max": None, "n_null": 0}
    except np.linalg.LinAlgError:
        pass
    ev, vec = np.linalg.eigh(L)
    tol = len(L) * np.finfo(float).eps * ev[-1]
    if ev[0] < -tol:
        top = np.argsort(-np.abs(vec[:, 0]))[:6].tolist()
        raise ValueError(
            f"inductance matrix has negative energy: min eigenvalue {ev[0]:.6e} H "
            f"({ev[0] / ev[-1]:.3e} of the largest; roundoff scale {tol:.3e} H), "
            f"{int(np.sum(ev < -tol))} negative modes; witness branches {top}")
    return {"cholesky": False, "min_eig_over_max": float(ev[0] / ev[-1]),
            "n_null": int(np.sum(np.abs(ev) <= tol))}

def solve_fin_volume(step_path, *, frequency, sigma, current_A, n_layers,
                     n_lanes, n_stations, n_outline=None, lane_grading=None,
                     tip_lanes=None, exponent=2.0, route="auto",
                     cad_units_per_meter="auto", max_branches=None,
                     box_order=BOX_ORDER):
    """STEP -> unreduced volume-lattice PEEC branch currents (peak phasors).

    Every inductance entry, self included, is the Galerkin integral of
    uniform currents over the bars' oriented boxes (kernel option
    ``oriented_box_order`` = ``box_order``; each bar's width axis comes from
    the lattice).  Filament mutuals make the lattice L indefinite: the
    transverse bars are short, wide plates.  Before the solve, L must pass
    :func:`require_energy_admissible`.
    """
    from radia.fin_surface_solver import coil_true_residuals
    from radia.fin_sweep import DEFAULT_MAX_BRANCHES, fin_graph_from_step
    from radia.peec_matrices import PEECBuilder
    from radia.peec_topology import PEECCircuitSolver

    if frequency <= 0 or sigma <= 0:
        raise ValueError("fin-volume needs frequency > 0 and sigma > 0")
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
    lat = build_volume_lattice(sweep, n_layers, exponent=exponent)

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
    peec = PEECCircuitSolver(topo)

    r_branch = np.asarray(np.diag(peec.R_dc) if peec.R_dc.ndim == 2 else peec.R_dc, dtype=float)
    current = np.asarray(peec.compute_branch_currents(frequency, [float(current_A)]),
                         dtype=complex)
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
            "coil_power_W": 0.5 * r_total * i2, "skin_depth_m": delta,
            "n_branches": int(len(lat.branches)), "n_nodes": int(len(lat.nodes)),
            "paths": paths, "I_fil": current}
