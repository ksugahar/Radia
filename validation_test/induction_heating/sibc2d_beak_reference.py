"""Independent 2-D SIBC reference for the straight beak-fin cross-section."""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))


def solve_beak_sibc_2d(step_path: Path, *, frequency=150_000.0,
                       sigma=5.8e7, n_peri=256):
    """Solve the constant-E 2-D Leontovich integral equation.

    Unknowns are panel currents ``I_k=K_k ds_k`` with sum(I_k)=1 A.
    The logarithmic vector-potential kernel is Galerkin midpoint for
    distinct panels and uses the exact straight-panel self average
    ``-log(ds)+3/2``. The arbitrary logarithm reference cancels through
    the terminal-voltage Lagrange multiplier.
    """
    from radia.peec_fin_topology import (
        build_hybrid_surface_topology_from_straight_prism_step,
    )

    graph, _ = build_hybrid_surface_topology_from_straight_prism_step(
        step_path, n_peri=n_peri, n_stations=3)
    xy = graph.branch_xyz[:n_peri, 0, :2]
    ds = 0.5 * (np.linalg.norm(xy - np.roll(xy, 1, axis=0), axis=1)
                + np.linalg.norm(np.roll(xy, -1, axis=0) - xy, axis=1))
    if np.min(ds) <= 0:
        raise ValueError("degenerate perimeter panel")
    distance = np.linalg.norm(xy[:, None, :] - xy[None, :, :], axis=2)
    np.fill_diagonal(distance, 1.0)
    mu0 = 4e-7 * math.pi
    kernel = mu0 / (2 * math.pi) * np.log(1.0 / distance)
    np.fill_diagonal(kernel, mu0 / (2 * math.pi) *
                     (np.log(1.0 / ds) + 1.5))
    omega = 2 * math.pi * frequency
    delta = math.sqrt(2 / (omega * mu0 * sigma))
    zs = (1 + 1j) / (sigma * delta)
    impedance = np.diag(zs / ds) + 1j * omega * kernel
    system = np.block([
        [impedance, -np.ones((n_peri, 1))],
        [np.ones((1, n_peri)), np.zeros((1, 1))],
    ])
    solution = np.linalg.solve(system, np.r_[np.zeros(n_peri), 1.0])
    panel_current = solution[:n_peri]
    k_surface = panel_current / ds
    resistance = float(zs.real * np.sum(np.abs(panel_current) ** 2 / ds))
    beak = xy[:, 0] >= 1.6e-3
    tip = xy[:, 0] >= 3.5e-3
    mean_k = 1.0 / np.sum(ds)
    loss = zs.real * np.abs(panel_current) ** 2 / ds
    probe = np.array([5.0e-3, 0.0])
    offset = probe - xy
    radius2 = np.sum(offset**2, axis=1)
    hx = np.sum(-panel_current * offset[:, 1] / (2 * math.pi * radius2))
    hy = np.sum(panel_current * offset[:, 0] / (2 * math.pi * radius2))
    return {
        "n_peri": n_peri,
        "perimeter_m": float(np.sum(ds)),
        "skin_depth_m": delta,
        "resistance_ohm_per_m": resistance,
        "beak_mean_absK_over_mean": float(np.mean(np.abs(k_surface[beak])) /
                                           mean_k),
        "tip_mean_absK_over_mean": float(np.mean(np.abs(k_surface[tip])) /
                                          mean_k),
        "max_absK_over_mean": float(np.max(np.abs(k_surface)) / mean_k),
        "beak_current_fraction": float(abs(np.sum(panel_current[beak]))),
        "beak_loss_fraction": float(np.sum(loss[beak]) / np.sum(loss)),
        "tip_loss_fraction": float(np.sum(loss[tip]) / np.sum(loss)),
        "beak_current_centroid_x_m": float(
            np.sum(np.abs(panel_current[beak]) * xy[beak, 0]) /
            np.sum(np.abs(panel_current[beak]))),
        "probe_xy_m": probe.tolist(),
        "probe_H_abs_A_per_m": float(np.sqrt(abs(hx)**2 + abs(hy)**2)),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--step", type=Path, default=(
        ROOT / "tests" / "coil_from_cad" / "fixtures" /
        "beak_fin_straight.step"))
    parser.add_argument("--n-peri", type=int, default=256)
    args = parser.parse_args()
    print(json.dumps(solve_beak_sibc_2d(args.step, n_peri=args.n_peri),
                     indent=2))


def test_beak_sibc_2d_golden():
    """Moderate mesh must resolve the converged R and rounded-tip K."""
    result = solve_beak_sibc_2d(
        ROOT / "tests" / "coil_from_cad" / "fixtures" /
        "beak_fin_straight.step", n_peri=256)
    assert abs(result["resistance_ohm_per_m"] / 0.005242 - 1) < 0.01
    assert abs(result["tip_mean_absK_over_mean"] / 2.0 - 1) < 0.05
    assert abs(result["beak_current_fraction"] / 0.31 - 1) < 0.03
    assert abs(result["beak_loss_fraction"] / 0.34 - 1) < 0.03
    assert abs(result["beak_current_centroid_x_m"] - 2.86e-3) < 0.05e-3
    assert abs(result["probe_H_abs_A_per_m"] / 41.78 - 1) < 0.01


if __name__ == "__main__":
    main()
