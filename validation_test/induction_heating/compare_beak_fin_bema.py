"""Exploratory beak-fin STEP comparison; not a fin-PEEC acceptance gate."""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "src" / "radia" / "panels"))


def _fill_periodic(values, valid):
    """Periodically interpolate unsampled perimeter lanes."""
    values = np.asarray(values, dtype=complex).copy()
    valid = np.asarray(valid, dtype=bool)
    if np.count_nonzero(valid) < 3:
        raise ValueError("need at least three sampled perimeter lanes")
    idx = np.flatnonzero(valid)
    query = np.arange(len(values))
    xp = np.r_[idx - len(values), idx, idx + len(values)]
    for component in ("real", "imag"):
        fp0 = getattr(values[idx], component)
        fp = np.tile(fp0, 3)
        part = np.interp(query, xp, fp)
        if component == "real":
            values.real = part
        else:
            values.imag = part
    return values


def _probe_offset(fin, delta):
    """Probe distance beyond the tip: outside the SIBC breakdown zone."""
    return max(5.0 * delta, 4.0 * fin.tip_radius)


def _distribution_metrics(k_surface, ds, xy, sheet_resistance, analysis,
                          fin, delta):
    from radia.fin_section import feature_panel_weights

    panel_current = np.asarray(k_surface, complex) * np.asarray(ds, float)
    panel_current /= np.sum(panel_current)
    loss = sheet_resistance * np.abs(panel_current) ** 2 / ds
    # Regions and probe come from the CAD outline (radia.fin_section), not
    # from fixture constants.
    beak_weight = feature_panel_weights(analysis, fin, xy)
    tip_weight = feature_panel_weights(analysis, fin, xy, tip=True)
    probe_xy = fin.probe_points(offset=_probe_offset(fin, delta), n=41)
    probe_axis = (probe_xy - probe_xy[len(probe_xy) // 2]) @ fin.normal
    offset = probe_xy[:, None, :] - xy[None, :, :]
    radius2 = np.sum(offset**2, axis=2)
    hx = np.sum(-panel_current[None, :] * offset[:, :, 1] /
                (2 * math.pi * radius2), axis=1)
    hy = np.sum(panel_current[None, :] * offset[:, :, 0] /
                (2 * math.pi * radius2), axis=1)
    h_abs = np.sqrt(np.abs(hx)**2 + np.abs(hy)**2)
    return {
        "beak_current_fraction": float(abs(np.sum(
            panel_current * beak_weight))),
        "beak_loss_fraction": float(np.sum(loss * beak_weight) /
                                    np.sum(loss)),
        "tip_loss_fraction": float(np.sum(loss * tip_weight) / np.sum(loss)),
        # centroid measured along the fin axis from the root chord
        "beak_current_centroid_x_m": float(
            np.sum(np.abs(panel_current) * beak_weight * fin.project(xy)) /
            np.sum(np.abs(panel_current) * beak_weight)
            + fin.root_mid @ fin.axis),
        "beak_current_centroid_axial_m": float(
            np.sum(np.abs(panel_current) * beak_weight * fin.project(xy)) /
            np.sum(np.abs(panel_current) * beak_weight)),
        "probe_center_H_abs_A_per_m": float(h_abs[len(h_abs) // 2]),
        "probe_center_xy_m": probe_xy[len(probe_xy) // 2].tolist(),
        "probe_axis_m": probe_axis.tolist(),
        "probe_H_abs_A_per_m": h_abs.tolist(),
    }


def solve_bema(step_path: Path, *, frequency: float, maxh: float,
               bema_solver: str = "lu", curve_order: int = 2,
               fes_order: int = 0, mesh_only: bool = False):
    """Solve the BEM-A reference once and return everything the comparison needs.

    Split out of ``run`` because the reference depends on the fixture, the
    mesh and the basis, and on nothing the PEEC side varies.  A refinement
    sweep over ``n_peri`` / ``n_stations`` would otherwise re-solve the same
    dense system at every point, which on the 48 mm fixture is the whole cost
    of the sweep.

    ``curve_order`` curves the BEM-A volume mesh before its boundary is used,
    so the 0.25 mm rounded tip is an arc rather than a polygon; ``fes_order``
    is the HDivSurface order of the BEM-A current.  Both were fixed at 1 and 0
    before 2026-09-21, when the tip deficit turned out to be a p-convergence
    problem: one order step moved R by +1.69% where curving moved it by
    +0.067%.  The same solve on the same boundary; only the representation
    changes.
    """
    from netgen.occ import OCCGeometry, Pnt
    from ngsolve import BND, Mesh, TaskManager

    from radia.bem.coil_inductance_ngsolve import (
        compute_centroids_areas_J,
        compute_inductance_source_sink,
    )

    sigma = 5.8e7
    mu0 = 4e-7 * math.pi
    omega = 2 * math.pi * frequency
    delta = math.sqrt(2 / (omega * mu0 * sigma))
    zs = (1 + 1j) / (sigma * delta)

    # OCCGeometry's STEP units are mm; scale to SI before NGSolve assembly.
    solid = OCCGeometry(str(step_path)).shape.Scale(Pnt(0, 0, 0), 1e-3)
    z_end = max(float(face.center[2]) for face in solid.faces)
    for face in solid.faces:
        z = face.center[2]
        face.name = ("source" if abs(z) < 1e-8 else
                     "sink" if abs(z - z_end) < 1e-8 else "body")
    with TaskManager():
        mesh = Mesh(OCCGeometry(solid).GenerateMesh(maxh=maxh))
        # The boundary of the curved volume mesh IS the BEM-A surface:
        # HDivSurface lives on it directly (the solver compresses away the
        # interior-edge DOFs), and Curve() has the CAD to project onto.  The
        # earlier flat extraction discarded the curvature and could not be
        # curved afterwards.
        mesh.Curve(int(curve_order))
        n_surface_faces = int(mesh.GetNE(BND))
        n_surface_vertices = len({int(v.nr) for el in mesh.Elements(BND)
                                  for v in el.vertices})
        print(f"BEM-A mesh: nface={n_surface_faces} nv={n_surface_vertices} "
              f"curve_order={curve_order} fes_order={fes_order}",
              file=sys.stderr, flush=True)
        if mesh_only:
            return {"mesh_only": True, "n_surface_faces": n_surface_faces,
                    "n_vertices": n_surface_vertices}
        bem = compute_inductance_source_sink(
            mesh, "source", "sink", fes_order=int(fes_order), omega=omega,
            Z_s_complex=zs, solver=bema_solver)
    cen, area, j_re = compute_centroids_areas_J(mesh, bem["gf_J"])
    _cen_im, _area_im, j_im = compute_centroids_areas_J(
        mesh, bem["gf_J_im"])
    j_complex = j_re + 1j * j_im
    j_abs2 = np.sum(np.abs(j_complex) ** 2, axis=1)
    end_tol = max(1e-12, z_end * 1e-8)
    source_cap = np.abs(cen[:, 2]) <= end_tol
    sink_cap = np.abs(cen[:, 2] - z_end) <= end_tol
    caps = source_cap | sink_cap
    # For a unit peak terminal current, Pavg=R/2. Hence the equivalent
    # resistance of a region is Rs*integral(|K|^2 dS).
    def region_resistance(mask):
        return float(zs.real * np.sum(j_abs2[mask] * area[mask]))

    cap_r_source = region_resistance(source_cap)
    cap_r_sink = region_resistance(sink_cap)
    cap_r = cap_r_source + cap_r_sink
    lateral_r = region_resistance(~caps)
    return {
        "mesh_only": False,
        "fixture": str(step_path), "frequency_hz": frequency,
        "sigma_S_per_m": sigma, "skin_depth_m": delta,
        "sigma": sigma, "omega": omega, "delta": delta, "zs": zs,
        "z_end": z_end, "cen": cen, "area": area, "j_complex": j_complex,
        "bema": {"R_ohm": float(bem["R"]), "L_H": float(bem["L"]),
                 "residual": float(bem["residual"]),
                 "n_J": int(bem["n_J"]), "n_f": int(bem["n_f"]),
                 "n_surface_faces": n_surface_faces,
                 "curve_order": int(curve_order),
                 "fes_order": int(fes_order),
                 "solver": bema_solver,
                 "R_source_cap_ohm": cap_r_source,
                 "R_sink_cap_ohm": cap_r_sink,
                 "R_caps_ohm": cap_r,
                 "R_lateral_ohm": lateral_r,
                 "R_partition_relative_closure": float(
                     (cap_r + lateral_r) / bem["R"] - 1)},
        "_R": float(bem["R"]), "_L": float(bem["L"]),
    }


def compare_peec(state: dict, step_path: Path, *, n_peri: int,
                 n_stations: int, experimental_peec: bool = False,
                 profile_csv: Path | None = None,
                 lane_grading: str = "uniform"):
    """Build the PEEC topology at one refinement and compare it to ``state``."""
    from radia.peec_fin_topology import (
        assemble_experimental_fin_peec,
        build_hybrid_surface_topology_from_straight_prism_step,
    )

    sigma = state["sigma"]
    delta = state["delta"]
    zs = state["zs"]
    frequency = state["frequency_hz"]
    z_end = state["z_end"]
    cen = state["cen"]
    area = state["area"]
    j_complex = state["j_complex"]

    graph, cad = build_hybrid_surface_topology_from_straight_prism_step(
        step_path, n_peri=n_peri, n_stations=n_stations,
        lane_grading=lane_grading)
    cad = dict(cad)
    analysis = cad.pop("section_analysis")
    fin = analysis.primary
    if fin is None:
        raise ValueError("no fin detected on the STEP section; the beak "
                         "comparison needs a protruding fin")

    result = {
        "fixture": state["fixture"], "frequency_hz": frequency,
        "sigma_S_per_m": sigma, "skin_depth_m": delta,
        "bema": state["bema"],
        "peec_topology": {**cad, "n_branches": len(graph.branches),
                          "n_nodes": len(graph.nodes)},
        "status": "BEM-A baseline only; PEEC physical R/L and current/loss comparison pending",
    }
    bem = {"R": state["_R"], "L": state["_L"]}

    # Gate 1 asks for station alignment, branch geometry, and the rejection of
    # missing tips and zero-area cells.  Report them as data rather than as an
    # exception so a refinement sweep can record where a level fails.
    station_z = np.array([graph.branch_xyz[s * graph.n_lanes, 0, 2]
                          for s in range(graph.n_stations - 1)]
                         + [graph.branch_xyz[
                             (graph.n_stations - 2) * graph.n_lanes, 1, 2]])
    spacing = np.diff(station_z)
    checks = {
        "n_stations": int(graph.n_stations),
        "n_lanes": int(graph.n_lanes),
        "stations_monotonic": bool(np.all(spacing > 0)),
        "station_spacing_relative_spread": float(
            (spacing.max() - spacing.min()) / spacing.mean())
        if spacing.size else 0.0,
        "axial_step_m": float(spacing.mean()) if spacing.size else 0.0,
    }
    result["geometry_checks"] = checks
    if experimental_peec:
        # The existing rectangular PEEC kernel is *not* a validated SIBC
        # discretization on a folded fin. This is a discrepancy probe only.
        rings = np.empty((graph.n_stations, graph.n_lanes, 3))
        for station in range(graph.n_stations - 1):
            lo = station * graph.n_lanes
            hi = lo + graph.n_lanes
            rings[station] = graph.branch_xyz[lo:hi, 0]
        rings[-1] = graph.branch_xyz[
            (graph.n_stations - 2) * graph.n_lanes:
            (graph.n_stations - 1) * graph.n_lanes, 1]
        peec, widths = assemble_experimental_fin_peec(
            graph, rings, sigma=sigma, sheet_depth=delta)
        # Leontovich sheet reactance equals sheet resistance at this omega.
        r_branch = (np.diag(peec.R_dc) if peec.R_dc.ndim == 2
                    else peec.R_dc)
        z_peec = peec.compute_port_impedance(
            frequency, Zs=1j * r_branch)
        current = peec.compute_branch_currents(
            frequency, [1.0], Zs=1j * r_branch)
        r_power = float(np.real(np.vdot(current, r_branch * current)))
        l_external = float(np.real(np.vdot(current, peec.L @ current)))
        result["experimental_peec"] = {
            "R_ohm": r_power,
            "L_external_H": l_external,
            "Z_port_ohm": [float(z_peec.real), float(z_peec.imag)],
            "relative_R_error": float(r_power / bem["R"] - 1),
            "relative_L_error": float(l_external / bem["L"] - 1),
        }
        jz = j_complex[:, 2]
        zlo, zhi = 0.35 * z_end, 0.65 * z_end
        body = (cen[:, 2] > zlo) & (cen[:, 2] < zhi) & (area > 0)

        # The middle PEEC axial segment supplies one K sample per CAD lane.
        mid_segment = (graph.n_stations - 1) // 2
        first = mid_segment * graph.n_lanes
        last = first + graph.n_lanes
        xy = rings[mid_segment, :, :2]
        k_peec = current[first:last] / np.asarray(widths)[first:last]

        # Assign each BEM triangle to its nearest sampled perimeter lane.
        # Area weighting integrates over the selected axial strip and avoids
        # bias from locally refined triangles at the rounded nose.
        dist2 = np.sum((cen[body, None, :2] - xy[None, :, :]) ** 2, axis=2)
        lane = np.argmin(dist2, axis=1)
        k_bem = np.zeros(graph.n_lanes, dtype=complex)
        lane_area = np.zeros(graph.n_lanes)
        for k in range(graph.n_lanes):
            take = lane == k
            if np.any(take):
                lane_area[k] = np.sum(area[body][take])
                k_bem[k] = np.sum(jz[body][take] * area[body][take]) / lane_area[k]
        valid = lane_area > 0
        k_bem = _fill_periodic(k_bem, valid)
        # Remove the arbitrary global phasor before a complex-field norm.
        phase = np.angle(np.vdot(k_peec[valid], k_bem[valid]))
        k_peec_aligned = k_peec * np.exp(1j * phase)
        scale = (np.vdot(k_peec_aligned[valid], k_bem[valid]) /
                 np.vdot(k_peec_aligned[valid], k_peec_aligned[valid])).real
        k_peec_aligned *= scale
        rel_l2 = float(np.linalg.norm(k_peec_aligned[valid] - k_bem[valid]) /
                       np.linalg.norm(k_bem[valid]))
        beak = fin.fin_mask(xy)
        tip = fin.tip_mask(xy)
        mid_widths = np.asarray(widths)[first:last]
        checks.update({
            "min_cell_width_m": float(mid_widths.min()),
            "no_zero_area_cells": bool(mid_widths.min() > 0.0),
            "beak_lanes": int(np.count_nonzero(beak)),
            "tip_lanes": int(np.count_nonzero(tip)),
            "tip_present": bool(np.count_nonzero(tip) > 0),
            "every_lane_sampled_by_bema": bool(np.all(valid)),
        })
        checks["passed"] = bool(
            checks["stations_monotonic"] and checks["no_zero_area_cells"]
            and checks["tip_present"] and checks["beak_lanes"] > 0
            and checks["station_spacing_relative_spread"] <= 1e-9)
        rs = zs.real
        loss_bem = 0.5 * rs * np.abs(k_bem) ** 2
        loss_peec = 0.5 * rs * np.abs(k_peec_aligned) ** 2
        metrics_bem = _distribution_metrics(
            k_bem, widths[first:last], xy, rs, analysis, fin, delta)
        metrics_peec = _distribution_metrics(
            k_peec, widths[first:last], xy, rs, analysis, fin, delta)
        h_bem = np.asarray(metrics_bem["probe_H_abs_A_per_m"])
        h_peec = np.asarray(metrics_peec["probe_H_abs_A_per_m"])
        result["surface_current_profile"] = {
            "axial_window_m": [zlo, zhi],
            "fin_regions": {
                "beak_lanes": int(np.count_nonzero(beak)),
                "tip_lanes": int(np.count_nonzero(tip)),
                "tip_arc_lanes": int(np.count_nonzero(fin.on_tip_arc(
                    analysis.arclength_of(xy), analysis.perimeter))),
                "probe_offset_m": _probe_offset(fin, delta),
            },
            "valid_lanes": int(np.count_nonzero(valid)),
            "complex_relative_L2": rel_l2,
            "beak_mean_absK_ratio_peec_over_bema": float(
                np.mean(np.abs(k_peec_aligned[beak])) /
                np.mean(np.abs(k_bem[beak]))),
            "tip_mean_absK_ratio_peec_over_bema": float(
                np.mean(np.abs(k_peec_aligned[tip])) /
                np.mean(np.abs(k_bem[tip]))),
            "beak_mean_loss_ratio_peec_over_bema": float(
                np.mean(loss_peec[beak]) / np.mean(loss_bem[beak])),
            "tip_mean_loss_ratio_peec_over_bema": float(
                np.mean(loss_peec[tip]) / np.mean(loss_bem[tip])),
            "integral_metrics": {"bema": metrics_bem, "peec": metrics_peec},
            "probe_H_relative_L2": float(
                np.linalg.norm(h_peec - h_bem) / np.linalg.norm(h_bem)),
            "probe_center_H_ratio_peec_over_bema": float(
                metrics_peec["probe_center_H_abs_A_per_m"] /
                metrics_bem["probe_center_H_abs_A_per_m"]),
        }
        profile = result["surface_current_profile"]
        profile["acceptance_limits"] = {
            "beak_current_and_loss_fraction_relative": 0.03,
            "centroid_absolute_m": float(np.mean(widths[first:last])),
            "probe_H_relative_L2": 0.02,
        }
        current_share_error = abs(
            metrics_peec["beak_current_fraction"] /
            metrics_bem["beak_current_fraction"] - 1)
        loss_share_error = abs(
            metrics_peec["beak_loss_fraction"] /
            metrics_bem["beak_loss_fraction"] - 1)
        centroid_error = abs(
            metrics_peec["beak_current_centroid_x_m"] -
            metrics_bem["beak_current_centroid_x_m"])
        profile["integral_metric_errors"] = {
            "beak_current_fraction_relative": float(current_share_error),
            "beak_loss_fraction_relative": float(loss_share_error),
            "beak_current_centroid_absolute_m": float(centroid_error),
        }
        profile["accepted"] = bool(
            current_share_error <= 0.03
            and loss_share_error <= 0.03
            and centroid_error <= np.mean(widths[first:last])
            and profile["probe_H_relative_L2"] <= 0.02)
        result["status"] = ("accepted" if profile["accepted"] else
                            "rejected: integral fin-delivery metrics mismatch")
        if profile_csv is not None:
            profile_csv.parent.mkdir(parents=True, exist_ok=True)
            with profile_csv.open("w", newline="", encoding="utf-8") as f:
                writer = csv.writer(f)
                writer.writerow(["lane", "x_m", "y_m", "region",
                                 "K_bema_abs_A_per_m", "K_peec_abs_A_per_m",
                                 "loss_bema_W_per_m2", "loss_peec_W_per_m2"])
                for k in range(graph.n_lanes):
                    region = "tip" if tip[k] else "beak" if beak[k] else "body"
                    writer.writerow([k, xy[k, 0], xy[k, 1], region,
                                     abs(k_bem[k]), abs(k_peec_aligned[k]),
                                     loss_bem[k], loss_peec[k]])
    return result


def run(step_path: Path, *, frequency: float, maxh: float,
        n_peri: int, n_stations: int, mesh_only: bool = False,
        experimental_peec: bool = False, profile_csv: Path | None = None,
        bema_solver: str = "lu", lane_grading: str = "uniform",
        curve_order: int = 2, fes_order: int = 0):
    """One BEM-A reference and one PEEC comparison against it."""
    state = solve_bema(step_path, frequency=frequency, maxh=maxh,
                       bema_solver=bema_solver, curve_order=curve_order,
                       fes_order=fes_order, mesh_only=mesh_only)
    if mesh_only:
        cad = build_section_topology_cad(step_path, n_peri=n_peri,
                                         n_stations=n_stations,
                                         lane_grading=lane_grading)
        return {"n_surface_faces": state["n_surface_faces"],
                "n_vertices": state["n_vertices"], "cad": cad}
    return compare_peec(state, step_path, n_peri=n_peri,
                        n_stations=n_stations,
                        experimental_peec=experimental_peec,
                        profile_csv=profile_csv, lane_grading=lane_grading)


def build_section_topology_cad(step_path: Path, *, n_peri: int,
                               n_stations: int, lane_grading: str):
    """The CAD summary ``--mesh-only`` reports, without the section analysis."""
    from radia.peec_fin_topology import (
        build_hybrid_surface_topology_from_straight_prism_step,
    )

    _graph, cad = build_hybrid_surface_topology_from_straight_prism_step(
        step_path, n_peri=n_peri, n_stations=n_stations,
        lane_grading=lane_grading)
    cad = dict(cad)
    cad.pop("section_analysis", None)
    return cad


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--step", type=Path, default=(
        ROOT / "tests" / "coil_from_cad" / "fixtures" /
        "beak_fin_short.step"))
    parser.add_argument("--frequency", type=float, default=150_000.0)
    parser.add_argument("--maxh", type=float, default=0.003)
    parser.add_argument("--n-peri", type=int, default=64)
    parser.add_argument("--n-stations", type=int, default=5)
    parser.add_argument("--mesh-only", action="store_true")
    parser.add_argument("--experimental-peec", action="store_true")
    parser.add_argument("--profile-csv", type=Path)
    parser.add_argument("--bema-solver", choices=("lu", "cocr", "hacapk_cocr"),
                        default="lu")
    parser.add_argument("--lane-grading", choices=("uniform", "auto"),
                        default="uniform",
                        help="auto: grade PEEC lanes toward the detected fin tip")
    parser.add_argument("--curve-order", type=int, default=2,
                        help="geometry order of the BEM-A boundary (2 makes "
                             "the rounded tip an arc; 1 reproduces the "
                             "pre-2026-09-21 flat surface)")
    parser.add_argument("--fes-order", type=int, default=0,
                        help="HDivSurface order of the BEM-A current; 0 is "
                             "RT0, the pre-2026-09-21 default, and does not "
                             "p-converge the tip")
    args = parser.parse_args()
    print(json.dumps(run(args.step, frequency=args.frequency,
                         maxh=args.maxh, n_peri=args.n_peri,
                         n_stations=args.n_stations,
                         mesh_only=args.mesh_only,
                         experimental_peec=args.experimental_peec,
                         profile_csv=args.profile_csv,
                         bema_solver=args.bema_solver,
                         lane_grading=args.lane_grading,
                         curve_order=args.curve_order,
                         fes_order=args.fes_order), indent=2))


if __name__ == "__main__":
    main()
