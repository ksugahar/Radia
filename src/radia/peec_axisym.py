"""Axisymmetric volume PEEC of a revolved conductor section.

Under an axisymmetric loop-voltage drive every meridian triangle of the section
is one azimuthal ring (P0 current); the circuit ``(R + j omega M) I = V`` uses
the Galerkin ring partial inductance of ``peec_matrices.axisym_ring_inductance``
and is solved by ``peec_matrices.axisym_ring_loop_solve`` (C++, LAPACK).
This module is the thin Python adapter: it checks inputs, runs the two native
steps and sums the losses over named regions.

Normalization: per 1 A peak terminal current, ``exp(+j omega t)``.
Assembly rows run in parallel only under a caller-owned ``ngsolve.TaskManager``.
"""

from __future__ import annotations

import time

import numpy as np

from .peec_matrices import axisym_ring_inductance, axisym_ring_loop_solve


def solve_axisym_ring_peec(points, triangles, frequency_hz, sigma, region=None,
                           near_factor=3.0, moderate_factor=10.0):
    """Impedance and ring-current losses of a revolved section.

    Args:
        points: (n, 2) meridian points (r, z) [m], every r > 0.
        triangles: (m, 3) zero-based point indices.
        frequency_hz, sigma: drive frequency [Hz] and conductivity [S/m].
        region: optional (m,) labels; losses are summed per label.

    Returns:
        dict with R_ohm, L_H, impedance_ohm, current_A (at V = 1), cell_loss_W_per_A2,
        region_loss_W_per_A2, power_balance (2 * total loss vs Re Z), true_relative_residual,
        cells, assembly_s, solve_s.
    """
    points = np.ascontiguousarray(points, dtype=float)
    triangles = np.ascontiguousarray(triangles)
    if points.ndim != 2 or points.shape[1] != 2:
        raise ValueError(f"points must have shape (n, 2), got {points.shape}")
    if triangles.ndim != 2 or triangles.shape[1] != 3 or not np.issubdtype(triangles.dtype, np.integer):
        raise ValueError(f"triangles must be an integer (m, 3) array, got {triangles.dtype} {triangles.shape}")
    if triangles.min() < 0 or triangles.max() >= len(points) or triangles.max() > np.iinfo(np.int32).max:
        raise ValueError("triangle indices outside the point array")
    if region is not None:
        region = np.asarray(region)
        if region.shape != (len(triangles),):
            raise ValueError(f"region must have one label per triangle, got {region.shape}")
    t0 = time.perf_counter()
    M, area, centroid = axisym_ring_inductance(points, triangles.astype(np.int32),
                                               near_factor=near_factor, moderate_factor=moderate_factor)
    t1 = time.perf_counter()
    sol = axisym_ring_loop_solve(M, area, centroid, float(frequency_hz), float(sigma))
    t2 = time.perf_counter()
    omega = 2*np.pi*float(frequency_hz)
    z = complex(sol["impedance"])
    loss = np.asarray(sol["cell_loss"])
    total = float(loss.sum())
    regional = None if region is None else {str(k): float(loss[region == k].sum()) for k in np.unique(region)}
    return dict(R_ohm=z.real, L_H=z.imag/omega, impedance_ohm=z, current_A=np.asarray(sol["current"]),
                cell_loss_W_per_A2=loss, region_loss_W_per_A2=regional,
                power_balance=abs(2*total - z.real)/z.real,
                true_relative_residual=float(sol["true_relative_residual"]),
                cells=int(len(triangles)), area_m2=area, centroid_m=centroid,
                assembly_s=t1 - t0, solve_s=t2 - t1)
