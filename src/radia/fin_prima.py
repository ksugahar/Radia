"""PRIMA port models of the fin volume-lattice PEEC (CPU C++ reduction).

The full port model is the modified nodal analysis of the lattice

    (G + s E) x = B u,  y = B^T x,  x = [I; V]
    E = [[L, 0], [0, 0]],  G = [[R, A^T], [-A, 0]],  B = [0; P]

with the port current as input and the port voltage as output.  The
reduction runs in C++ (:meth:`radia.peec_matrices.PEECBuilder.prima_reduce`):
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


def reduce_topology(builder, topology, order, *, s0=0.0, port_names=None):
    """PRIMA model of a built ``PEECBuilder`` with the incidence in ``topology``.

    ``topology`` holds ``segment_nodes`` (branch node pairs), ``n_nodes`` and
    ``ports``; ``builder`` holds the L and R of those branches
    (``build_topology`` called).  Returns ``(DescriptorPortModel, info)``.
    """
    plus, minus = _ports(topology)
    branches = np.ascontiguousarray(np.asarray(topology["segment_nodes"], dtype=np.int32))
    reduced = builder.prima_reduce(branches, int(topology["n_nodes"]), plus, minus,
                                   int(order), float(s0))
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


def order_error_table(builder, topology, orders, frequency, *, s0=0.0):
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
        model, info = reduce_topology(builder, topology, order, s0=s0)
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
