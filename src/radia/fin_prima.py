"""PRIMA port models of the fin volume-lattice PEEC (CPU C++ reduction).

The full port model is the modified nodal analysis of the lattice

    (G + s E) x = B u,  y = B^T x,  x = [I; V]
    E = [[L, 0], [0, 0]],  G = [[R, A^T], [-A, 0]],  B = [0; P]

with the port current as input and the port voltage as output.  The
reduction runs in C++ (:func:`radia.peec_matrices.prima_reduce`), which takes
L, R and the incidence together from one assembled topology and requires L
itself to be positive definite:
block Arnoldi on ``(G + s0 E)^-1 E`` from ``(G + s0 E)^-1 B`` and the
congruence ``V^T E V``, ``V^T G V``, ``V^T B``, which keeps passivity.

The reduced order is chosen by the measured port-impedance error against the
unreduced model at the operating frequency; the tolerance is the caller's
(user-owned acceptance), never a default here.
"""
from __future__ import annotations

import math

import numpy as np

from radia.prima_export import DescriptorPortModel


def _uniform_field_emf_columns(lattice):
    """Return the branch RHS coefficients of s H for exp(+j omega t).

    Each bar path includes short legs to its circuit nodes. At collapsed
    terminals these nodes are the station centroids. The legs make a change
    of coordinate origin a discrete voltage gradient at the circuit nodes.
    Columns correspond to Hx, Hy and Hz in A/m. Coordinates are in metres.
    """
    nodes = np.asarray(lattice.nodes, dtype=float)
    branches = np.asarray(lattice.branches)
    bars = np.asarray(lattice.branch_xyz, dtype=float)
    start, end = nodes[branches[:, 0]], nodes[branches[:, 1]]
    p0, p1 = bars[:, 0], bars[:, 1]
    area = np.cross(start, p0) + np.cross(p0, p1) + np.cross(p1, end)
    # A = mu0 (H cross r) / 2; the impressed voltage is -s integral(A dl).
    return -(2e-7 * math.pi) * area



def _reduce_uniform_field(topo, lattice, order, s0):
    """Build a split current/potential basis for one port and three field axes."""
    from scipy import linalg
    from scipy.sparse import csr_matrix

    L, R = np.asarray(topo["L"]), np.asarray(topo["R"])
    branches = np.asarray(lattice.branches)
    count, nodes = len(branches), len(lattice.nodes)
    plus, minus = topo["ports"][0][:2]
    indices = np.arange(count)
    incidence = csr_matrix((np.r_[-np.ones(count), np.ones(count)],
                            (branches.T.ravel(), np.tile(indices, 2))),
                           shape=(nodes, count))
    keep = np.arange(nodes) != minus
    A = incidence[keep]
    source = np.zeros(nodes)
    source[plus], source[minus] = 1., -1.
    flux = _uniform_field_emf_columns(lattice)

    # Solve the real shifted saddle system by eliminating branch current.
    if order >= count + nodes - 1:
        raise ValueError("uniform-field PRIMA order must be below the full state dimension")
    linalg.cholesky(L, check_finite=True)
    K = np.array(L, order="F", copy=True)
    K *= s0
    K[indices, indices] += R
    factor = linalg.cho_factor(K, overwrite_a=True)
    lift = linalg.cho_solve(factor, A.T.toarray())
    schur = linalg.cho_factor(A @ lift, overwrite_a=True)

    def shifted_solve(branch_rhs, node_rhs):
        current = linalg.cho_solve(factor, branch_rhs)
        potential = linalg.cho_solve(schur, node_rhs - A @ current)
        return np.r_[current + lift @ potential, potential]

    basis = np.zeros((count + nodes - 1, order))
    rank = 0
    deflated = 0

    def append(vector):
        nonlocal rank, deflated
        norm = linalg.norm(vector)
        active = basis[:, :rank]
        for _ in range(2):
            vector -= active @ (active.T @ vector)
        remainder = linalg.norm(vector)
        if remainder <= 1e-10 * norm:
            deflated += 1
        elif rank == order:
            raise ValueError("PRIMA order is below the independent field start-block rank")
        else:
            basis[:, rank] = vector / remainder
            rank += 1

    zero_node = np.zeros(nodes - 1)
    append(shifted_solve(np.zeros(count), source[keep]))
    for column in flux.T:
        append(shifted_solve(column, zero_node))
    index = 0
    while index < rank < order:
        append(shifted_solve(L @ basis[:count, index], zero_node))
        index += 1

    def range_basis(matrix):
        vectors, singular, _ = linalg.svd(matrix, full_matrices=False)
        return vectors[:, singular > singular[0] * 1e-10]

    potential = range_basis(basis[count:, :rank])
    current = range_basis(np.column_stack((range_basis(basis[:count, :rank]),
                                           A.T @ potential)))
    coupling = potential.T @ (A @ current)
    m, n = current.shape[1], potential.shape[1]
    G = np.block([[current.T @ (R[:, None] * current), -coupling.T],
                  [coupling, np.zeros((n, n))]])
    E = np.zeros_like(G)
    E[:m, :m] = current.T @ (L @ current)
    B = np.zeros((m + n, 4))
    B[m:, 0] = potential.T @ source[keep]
    B[:m, 1:] = current.T @ flux
    return (G, E, B, current, potential, incidence, source, flux,
            {"q": m + n, "arnoldi_rank": rank, "deflated": deflated, "s0": s0})


def _evaluate_uniform_field(model, topo, frequency, drive):
    """Lift and gate each input and the combined drive."""
    G, E, B, current, potential, incidence, source, flux, info = model
    s = 2j * math.pi * frequency
    rhs = B.astype(complex)
    rhs[:, 1:] *= s
    # Evaluate all four inputs, including zero columns, then the physical drive.
    scales = np.where(np.abs(drive) > 0, np.abs(drive), 1.)
    rhs *= scales
    responses = np.linalg.solve(G + s * E, rhs)
    combined = responses @ (drive / scales)
    lifted = np.column_stack((responses, combined))
    m = current.shape[1]
    currents = current @ lifted[:m]
    voltages = np.zeros((len(source), 5), dtype=complex)
    minus = topo["ports"][0][1]
    voltages[np.arange(len(source)) != minus] = potential @ lifted[m:]
    drives = np.column_stack((np.diag(scales), drive))
    injection = source[:, None] * drives[0]
    voltage = incidence.T @ voltages + s * flux @ drives[1:]
    kcl_error = np.linalg.norm(incidence @ currents - injection, axis=0)
    current_scale = np.where(drives[0] != 0, np.linalg.norm(injection, axis=0),
                             np.linalg.norm(currents, axis=0))
    kvl_error = np.linalg.norm(np.asarray(topo["R"])[:, None] * currents
                               + s * (np.asarray(topo["L"]) @ currents) - voltage, axis=0)
    tiny = np.finfo(float).tiny
    kcl = kcl_error / np.maximum(current_scale, tiny)
    kvl = kvl_error / np.maximum(np.linalg.norm(voltage, axis=0), tiny)
    if not np.all((kcl <= 1e-6) & (kvl <= 1e-6)):
        raise RuntimeError(f"{frequency:g} Hz uniform-field PRIMA failed true KCL/KVL: {kcl=}, {kvl=}")
    return currents[:, -1], dict(info, kcl_relative=kcl, kvl_relative=kvl)

def _ports(topology):
    ports = list(topology["ports"])
    if not ports:
        raise ValueError("the topology has no port")
    return [int(p[0]) for p in ports], [int(p[1]) for p in ports]


def _port_model(reduced, port_names, metadata):
    G, E, B = (np.asarray(reduced[k], dtype=float) for k in ("G", "E", "B"))
    p = B.shape[1]
    names = tuple(port_names) if port_names else tuple(f"P{k + 1}" for k in range(p))
    return DescriptorPortModel(E, -G, B, B.T.copy(), np.zeros((p, p)), "current", names,
                               dict(metadata))


def reduce_topology(topology, order, *, s0=0.0, port_names=None):
    """PRIMA model of one assembled topology.

    ``topology`` holds ``L``, ``R``, ``segment_nodes`` (branch node pairs),
    ``n_nodes`` and ``ports`` of the same model, so energy and incidence
    cannot come from different models.  Returns ``(DescriptorPortModel, info)``.
    """
    from radia.peec_matrices import prima_reduce

    plus, minus = _ports(topology)
    branches = np.ascontiguousarray(np.asarray(topology["segment_nodes"], dtype=np.int32))
    reduced = prima_reduce(np.asarray(topology["L"], dtype=float),
                           np.asarray(topology["R"], dtype=float), branches,
                           int(topology["n_nodes"]), plus, minus, int(order), float(s0))
    info = {k: reduced[k] for k in ("q", "n_full", "s0", "deflated",
                                    "max_solve_residual", "min_E_eig_over_max")}
    model = _port_model(reduced, port_names, {"source": "radia.fin_prima", **info,
                                              "order_requested": int(order)})
    return model, info


def full_port_impedance(L, R, branches, n_nodes, ports, frequency):
    """Unreduced port impedance matrix from a dense complex MNA solve (numpy).

    Independent of the C++ reduction; used as the order-selection reference.
    """
    L = np.asarray(L, dtype=float)
    R = np.asarray(R, dtype=float)
    branches = np.asarray(branches, dtype=np.int64)
    n = len(R)
    plus = [int(p[0]) for p in ports]
    minus = [int(p[1]) for p in ports]
    ground = minus[0]
    keep = np.array([k for k in range(int(n_nodes)) if k != ground])
    pos = -np.ones(int(n_nodes), dtype=np.int64)
    pos[keep] = np.arange(len(keep))
    A = np.zeros((len(keep), n))
    for b, (u, v) in enumerate(branches):
        if pos[u] >= 0:
            A[pos[u], b] -= 1.0
        if pos[v] >= 0:
            A[pos[v], b] += 1.0
    P = np.zeros((len(keep), len(plus)))
    for k, (a, c) in enumerate(zip(plus, minus)):
        if pos[a] >= 0:
            P[pos[a], k] += 1.0
        if pos[c] >= 0:
            P[pos[c], k] -= 1.0
    s = 2j * math.pi * float(frequency)
    M = np.block([[np.diag(R) + s * L, A.T], [-A, np.zeros((len(keep), len(keep)))]])
    rhs = np.vstack([np.zeros((n, len(plus))), P])
    x = np.linalg.solve(M, rhs)
    residual = np.linalg.norm(M @ x - rhs) / np.linalg.norm(rhs)
    if not residual <= 1e-6:
        raise RuntimeError(f"unreduced MNA true residual {residual:.3e} exceeds 1e-6")
    return rhs.T @ x


def order_error_table(topology, orders, frequency, *, s0=0.0):
    """Relative port-impedance error of each PRIMA order at ``frequency``.

    The error is ``||Z_q - Z_full||_F / ||Z_full||_F`` with ``Z_full`` from
    :func:`full_port_impedance` on the same L, R and incidence.
    """
    plus, minus = _ports(topology)
    ports = list(zip(plus, minus))
    z_full = full_port_impedance(topology["L"], topology["R"], topology["segment_nodes"],
                                 topology["n_nodes"], ports, frequency)
    s = 2j * math.pi * float(frequency)
    rows = []
    for order in orders:
        model, info = reduce_topology(topology, order, s0=s0)
        z_q = model.response(s)
        rows.append({"order": int(order), "q": int(info["q"]),
                     "Z_ohm": z_q.tolist() if z_q.size > 1 else complex(z_q[0, 0]),
                     "rel_error": float(np.linalg.norm(z_q - z_full) / np.linalg.norm(z_full)),
                     "max_solve_residual": info["max_solve_residual"]})
    return {"frequency_Hz": float(frequency), "s0": float(s0),
            "Z_full_ohm": z_full.tolist() if z_full.size > 1 else complex(z_full[0, 0]),
            "rows": rows}


def select_order(table, tolerance):
    """Smallest order of ``table`` whose measured error is <= ``tolerance``.

    ``tolerance`` is user-owned and required.  Raises when no order meets it.
    """
    if tolerance is None or not tolerance > 0:
        raise ValueError("select_order needs an explicit tolerance > 0")
    passing = [row for row in table["rows"] if row["rel_error"] <= tolerance]
    if not passing:
        best = min(row["rel_error"] for row in table["rows"])
        raise ValueError(f"no PRIMA order meets {tolerance:g} at "
                         f"{table['frequency_Hz']:g} Hz (best {best:.3e})")
    return min(passing, key=lambda row: row["q"])


__all__ = ["full_port_impedance", "order_error_table", "reduce_topology", "select_order"]
