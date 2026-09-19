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


def run(step_path: Path, *, frequency: float, maxh: float,
        n_peri: int, n_stations: int, mesh_only: bool = False,
        experimental_peec: bool = False, profile_csv: Path | None = None,
        bema_solver: str = "lu"):
    from netgen.occ import OCCGeometry, Pnt
    from ngsolve import Mesh, TaskManager
    from surface_mesh_extract import _extract_surface_mesh_filtered

    from radia.bem.coil_inductance_ngsolve import (
        compute_centroids_areas_J,
        compute_inductance_source_sink,
    )
    from radia.peec_fin_topology import (
        assemble_experimental_fin_peec,
        build_hybrid_surface_topology_from_straight_prism_step,
    )

    sigma = 5.8e7
    mu0 = 4e-7 * math.pi
    omega = 2 * math.pi * frequency
    delta = math.sqrt(2 / (omega * mu0 * sigma))
    zs = (1 + 1j) / (sigma * delta)

    graph, cad = build_hybrid_surface_topology_from_straight_prism_step(
        step_path, n_peri=n_peri, n_stations=n_stations)

    # OCCGeometry's STEP units are mm; scale to SI before NGSolve assembly.
    solid = OCCGeometry(str(step_path)).shape.Scale(Pnt(0, 0, 0), 1e-3)
    z_end = max(float(face.center[2]) for face in solid.faces)
    for face in solid.faces:
        z = face.center[2]
        face.name = ("source" if abs(z) < 1e-8 else
                     "sink" if abs(z - z_end) < 1e-8 else "body")
    with TaskManager():
        mesh = Mesh(OCCGeometry(solid).GenerateMesh(maxh=maxh))
        surface = _extract_surface_mesh_filtered(mesh, keep_label="")
        print(f"BEM-A mesh: nface={surface.nface} nv={surface.nv}",
              file=sys.stderr, flush=True)
        if mesh_only:
            return {"n_surface_faces": int(surface.nface),
                    "n_vertices": int(surface.nv), "cad": cad}
        bem = compute_inductance_source_sink(
            surface, "source", "sink", omega=omega, Z_s_complex=zs,
            solver=bema_solver)
    cen, area, j_re = compute_centroids_areas_J(surface, bem["gf_J"])
    _cen_im, _area_im, j_im = compute_centroids_areas_J(
        surface, bem["gf_J_im"])
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
    result = {
        "fixture": str(step_path), "frequency_hz": frequency,
        "sigma_S_per_m": sigma, "skin_depth_m": delta,
        "bema": {"R_ohm": float(bem["R"]), "L_H": float(bem["L"]),
                 "residual": float(bem["residual"]),
                 "n_J": int(bem["n_J"]), "n_f": int(bem["n_f"]),
                 "n_surface_faces": int(surface.nface),
                 "solver": bema_solver,
                 "R_source_cap_ohm": cap_r_source,
                 "R_sink_cap_ohm": cap_r_sink,
                 "R_caps_ohm": cap_r,
                 "R_lateral_ohm": lateral_r,
                 "R_partition_relative_closure": float(
                     (cap_r + lateral_r) / bem["R"] - 1)},
        "peec_topology": {**cad, "n_branches": len(graph.branches),
                          "n_nodes": len(graph.nodes)},
        "status": "BEM-A baseline only; PEEC physical R/L and current/loss comparison pending",
    }
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
        # Remove the arbitrary global phasor before a complex-field norm.
        phase = np.angle(np.vdot(k_peec[valid], k_bem[valid]))
        k_peec_aligned = k_peec * np.exp(1j * phase)
        scale = (np.vdot(k_peec_aligned[valid], k_bem[valid]) /
                 np.vdot(k_peec_aligned[valid], k_peec_aligned[valid])).real
        k_peec_aligned *= scale
        rel_l2 = float(np.linalg.norm(k_peec_aligned[valid] - k_bem[valid]) /
                       np.linalg.norm(k_bem[valid]))
        beak = xy[:, 0] >= 1.6e-3
        tip = xy[:, 0] >= 3.5e-3
        rs = zs.real
        loss_bem = 0.5 * rs * np.abs(k_bem) ** 2
        loss_peec = 0.5 * rs * np.abs(k_peec_aligned) ** 2
        result["surface_current_profile"] = {
            "axial_window_m": [zlo, zhi],
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
        }
        profile = result["surface_current_profile"]
        profile["acceptance_limits"] = {
            "complex_relative_L2_max": 0.05,
            "absK_mean_ratio_range": [0.95, 1.05],
            "loss_mean_ratio_range": [0.90, 1.10],
        }
        profile["accepted"] = bool(
            profile["complex_relative_L2"] <= 0.05
            and all(0.95 <= profile[name] <= 1.05 for name in (
                "beak_mean_absK_ratio_peec_over_bema",
                "tip_mean_absK_ratio_peec_over_bema"))
            and all(0.90 <= profile[name] <= 1.10 for name in (
                "beak_mean_loss_ratio_peec_over_bema",
                "tip_mean_loss_ratio_peec_over_bema")))
        result["status"] = ("accepted" if profile["accepted"] else
                            "rejected: beak/tip current distribution mismatch")
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
    args = parser.parse_args()
    print(json.dumps(run(args.step, frequency=args.frequency,
                         maxh=args.maxh, n_peri=args.n_peri,
                         n_stations=args.n_stations,
                         mesh_only=args.mesh_only,
                         experimental_peec=args.experimental_peec,
                         profile_csv=args.profile_csv,
                         bema_solver=args.bema_solver), indent=2))


if __name__ == "__main__":
    main()
