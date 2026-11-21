"""Fin-capable surface PEEC solve of one STEP conductor (terminal + fin metrics).

Shared by ``validation_test/induction_heating/fin_peec_from_step.py`` and the
``--coil-solver fin-surface`` backend of ``panels/calc_inductance.py`` so both
report the same numbers from the same code.  Exploratory kernel: the
surface-branch partial elements are the existing rectangular Ruehli
kernel with Leontovich sheet impedance; no workpiece back-reaction here.
"""

from __future__ import annotations

import math

import numpy as np


def coil_true_residuals(peec, z_diag, current, frequency, port_current):
    """True residuals of the solved branch currents against the MNA equations.

    The kernel solves ``Z i = -A^T V`` with ``A i`` injecting the port
    current, where ``Z = diag(z_diag) + j w L``.  KCL is measured against the
    fixed port current.  KVL is measured as the part of ``Z i`` outside the
    range of ``A^T`` (the best node potentials, from the reduced graph
    Laplacian), relative to ``||Z i||``.  Raises above the shared limit.
    """
    import scipy.sparse as sp
    from scipy.sparse.linalg import splu

    from radia._residual_gate import RELATIVE_LIMIT

    m = len(current)
    src = np.asarray(peec.segment_nodes[:, 0], dtype=np.int64)
    dst = np.asarray(peec.segment_nodes[:, 1], dtype=np.int64)
    n_nodes = int(peec.n_nodes)
    a = sp.csr_matrix((np.r_[-np.ones(m), np.ones(m)],
                       (np.r_[src, dst], np.r_[np.arange(m), np.arange(m)])),
                      shape=(n_nodes, m))
    pos, neg = (int(k) for k in peec.ports[0][:2])
    inj = np.zeros(n_nodes, dtype=complex)
    inj[pos], inj[neg] = -complex(port_current), complex(port_current)
    kcl = float(np.linalg.norm(a @ current - inj)) / max(abs(port_current), 1e-300)

    omega = 2 * math.pi * float(frequency)
    z_i = np.asarray(z_diag) * current + 1j * omega * (np.asarray(peec.L) @ current)
    keep = np.arange(n_nodes) != neg
    ak = a[keep]
    lu = splu((ak @ ak.T).tocsc().astype(complex))
    v = lu.solve(-(ak @ z_i))
    kvl = float(np.linalg.norm(z_i + ak.T @ v)) / max(float(np.linalg.norm(z_i)), 1e-300)
    for name, value in (("KCL", kcl), ("KVL", kvl)):
        if not math.isfinite(value) or value > RELATIVE_LIMIT:
            raise RuntimeError(
                f"fin-surface coil {name} true relative residual {value:.3e} "
                f"exceeds {RELATIVE_LIMIT:g}")
    return {"kcl_relative": kcl, "kvl_relative": kvl,
            "limit": RELATIVE_LIMIT}


def solve_fin_surface(step_path, *, frequency, sigma, n_lanes=None,
                      n_stations=None, n_outline=None, lane_grading=None,
                      tip_lanes=None, cad_units_per_meter="auto",
                      route="auto", probe_points=41, current_A=1.0):
    """STEP -> fin-graded graph -> branch currents -> per-station fin metrics.

    Discretisation arguments left as ``None`` are measured from the
    geometry (:func:`radia.fin_sweep.auto_fin_resolution`), and
    ``cad_units_per_meter="auto"`` is resolved from the conductor extent
    when the STEP settles it -- so a beak-fin STEP needs nothing but the
    file, the frequency and the conductivity.  The values actually used
    and the measurements behind them are in ``sweep.meta``.

    Returns a dict with ``sweep`` (FinSweepGraph), ``current`` (branch
    currents at ``current_A``), ``widths``, ``R_ohm``, ``L_external_H``,
    ``stations`` (JSON-ready per-station metrics) and ``paths`` /
    ``I_fil`` in the ``calc_inductance`` coil_data convention (one
    two-point polyline per branch).
    """
    from radia.fin_section import feature_panel_weights
    from radia.fin_sweep import fin_graph_from_step, probe_points_3d
    from radia.peec_fin_topology import assemble_experimental_fin_peec
    from radia.peec_proximity import _biot_savart_H

    if frequency <= 0 or sigma <= 0:
        raise ValueError("fin-surface needs frequency > 0 and sigma > 0 "
                         "(thin-skin sheet model)")
    mu0 = 4e-7 * math.pi
    omega = 2 * math.pi * frequency
    delta = math.sqrt(2 / (omega * mu0 * sigma))
    rs = 1.0 / (sigma * delta)

    sweep = fin_graph_from_step(
        step_path, n_lanes=n_lanes, n_stations=n_stations,
        n_outline=n_outline, lane_grading=lane_grading, tip_lanes=tip_lanes,
        cad_units_per_meter=cad_units_per_meter, route=route)
    # fin_graph_from_step measured whatever was left as None.
    tip_lanes = int(sweep.meta["resolution"]["tip_lanes"])
    graph = sweep.graph
    peec, widths = assemble_experimental_fin_peec(
        graph, sweep.rings, sigma=sigma, sheet_depth=delta)
    r_branch = np.diag(peec.R_dc) if peec.R_dc.ndim == 2 else peec.R_dc
    r_branch = np.asarray(r_branch, dtype=float)
    current = np.asarray(peec.compute_branch_currents(
        frequency, [float(current_A)], Zs=1j * r_branch), dtype=complex)
    if not np.all(np.isfinite(current)):
        raise RuntimeError("fin-surface PEEC returned non-finite branch currents")
    residuals = coil_true_residuals(peec, (1 + 1j) * r_branch, current,
                                    frequency, float(current_A))
    i2 = abs(current_A) ** 2
    r_total = float(np.real(np.vdot(current, r_branch * current)) / i2)
    l_external = float(np.real(np.vdot(current, peec.L @ current)) / i2)
    widths = np.asarray(widths, dtype=float)

    starts = graph.branch_xyz[:, 0, :]
    ends = graph.branch_xyz[:, 1, :]
    per_station = []
    for feat, st in zip(sweep.features, sweep.stations):
        s = feat.index
        if s >= graph.n_stations - 1:
            break
        lo, hi = s * graph.n_lanes, (s + 1) * graph.n_lanes
        i_lane = current[lo:hi]
        w_lane = widths[lo:hi]
        k_lane = i_lane / w_lane
        loss = rs * np.abs(i_lane) ** 2 / w_lane            # W/m per lane
        entry = {"station": int(s), "fin": feat.fin is not None,
                 "graded": bool(feat.graded),
                 "lane_current_sum_A": float(abs(np.sum(i_lane))),
                 "lateral_loss_W_per_m": float(np.sum(loss)),
                 "max_absK_over_mean": float(
                     np.max(np.abs(k_lane)) / (abs(np.sum(i_lane)) / np.sum(w_lane)))}
        fin = feat.fin
        if fin is not None:
            wb = feature_panel_weights(feat.analysis, fin, feat.lane_uv)
            wt = feature_panel_weights(feat.analysis, fin, feat.lane_uv, tip=True)
            lane_s = feat.analysis.arclength_of(feat.lane_uv)
            axial = fin.project(feat.lane_uv)
            offset = max(5.0 * delta, 4.0 * fin.tip_radius)
            probe = probe_points_3d(feat, st, offset, n=probe_points)
            h = _biot_savart_H(probe, starts, ends, current)
            h_abs = np.sqrt(np.sum(np.abs(h) ** 2, axis=1))
            entry.update({
                "beak_current_fraction": float(abs(np.sum(i_lane * wb)) /
                                               abs(np.sum(i_lane))),
                "beak_loss_fraction": float(np.sum(loss * wb) / np.sum(loss)),
                "tip_loss_fraction": float(np.sum(loss * wt) / np.sum(loss)),
                "beak_current_centroid_axial_m": float(
                    np.sum(np.abs(i_lane) * wb * axial) /
                    np.sum(np.abs(i_lane) * wb)),
                "tip_lanes": int(tip_lanes),
                "tip_arc_lane_centers": int(np.count_nonzero(
                    fin.on_tip_arc(lane_s, feat.analysis.perimeter))),
                "tip_overlap_majority_panels": int(np.count_nonzero(wt > 0.5)),
                "fin_length_m": float(fin.length),
                "tip_radius_m": float(fin.tip_radius),
                "tip_radius_over_delta": float(fin.tip_radius / delta),
                "probe_offset_m": float(offset),
                "probe_center_H_abs_A_per_m": float(h_abs[len(h_abs) // 2]),
                "probe_H_abs_A_per_m": h_abs.tolist(),
                "probe_xyz_m": probe.tolist(),
            })
        per_station.append(entry)

    paths = [[(tuple(map(float, p0)), tuple(map(float, p1)))]
             for p0, p1 in zip(starts, ends)]
    return {
        "sweep": sweep, "current": current, "widths": widths,
        "R_ohm": r_total, "L_external_H": l_external,
        "coil_true_residual": residuals,
        "skin_depth_m": delta, "sheet_resistance_ohm": rs,
        "resolution": sweep.meta["resolution"],
        "auto_resolution": sweep.meta["auto_resolution"],
        "n_branches": int(len(graph.branches)),
        "n_nodes": int(len(graph.nodes)),
        "stations": per_station,
        "paths": paths, "I_fil": current,
    }
