"""Verify the Loop-form bundle reduction matches full nodal MNA.

build_loop_bundle_impedance collapses a series-chain filament bundle's
nodal MNA (n_nodes x n_nodes) to a K x K (K = N_filaments) Loop-form
problem. For series-chain bundles this is exact; the test confirms that
solve_loop_bundle gives the same per-filament currents as the full MNA.

Also times both to show the scaling benefit as K grows.
"""

import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.abspath(os.path.join(HERE, '..', '..', '..', 'src'))
SRC_RADIA = os.path.join(SRC, 'radia')
for p in (SRC, SRC_RADIA):
    if p not in sys.path:
        sys.path.insert(0, p)

import radia  # noqa: F401  -- sets up MKL DLL paths


def build_test_coil(R, a, gap_deg, nw, nh, n_arc, freq, sigma):
    from coil_builder import CoilBuilder
    cb = (CoilBuilder(current=1.0)
          .set_start([R, 0, 0],
                     orientation=np.array([[1, 0, 0],
                                            [0, 1, 0],
                                            [0, 0, 1]]))
          .set_cross_section(width=2 * a, height=2 * a)
          .add_arc(radius=R, arc_angle=360.0 - gap_deg, tilt=0))
    paths, _ = cb.to_filaments(nw, nh, frequency=freq, sigma=sigma, n_arc=n_arc)
    return paths


def run(nw, nh, n_arc):
    R, a = 0.030, 0.003
    gap_deg = 5.0
    sigma_cu = 5.8e7
    freq = 7000.0

    print(f"\n=== nw x nh = {nw} x {nh}, n_arc = {n_arc} "
          f"(K = {nw*nh}, N_seg = {nw*nh*n_arc}) ===")

    from peec_bundle import (build_bundle_solver, filament_currents,
                              build_loop_bundle_impedance, solve_loop_bundle)
    paths = build_test_coil(R, a, gap_deg, nw, nh, n_arc, freq, sigma_cu)
    dw = (2 * a) / nw
    dh = (2 * a) / nh

    t0 = time.perf_counter()
    solver, seg_of_fil, _, _ = build_bundle_solver(paths, dw, dh, sigma_cu)
    t_build = time.perf_counter() - t0

    # Full MNA solve (existing path)
    t0 = time.perf_counter()
    I_branch = solver.compute_branch_currents(freq, [1.0])
    I_fil_mna = filament_currents(I_branch, seg_of_fil)
    t_mna = time.perf_counter() - t0

    # Loop reduction + solve
    t0 = time.perf_counter()
    R_f, L_f = build_loop_bundle_impedance(solver, seg_of_fil)
    t_reduce = time.perf_counter() - t0
    t0 = time.perf_counter()
    I_fil_loop, V_port = solve_loop_bundle(R_f, L_f, freq, I_port=1.0)
    t_loop = time.perf_counter() - t0

    err = np.max(np.abs(I_fil_mna - I_fil_loop)) / np.max(np.abs(I_fil_mna))
    print(f"  build_bundle_solver : {t_build*1e3:7.1f} ms")
    print(f"  full MNA solve      : {t_mna*1e3:7.1f} ms")
    print(f"  Loop reduction      : {t_reduce*1e3:7.1f} ms")
    print(f"  Loop K x K solve    : {t_loop*1e3:7.1f} ms")
    print(f"  max rel err |I_fil| : {err:.2e}")
    if err > 1e-6:
        print(f"  WARNING: Loop result differs from MNA by {err*100:.3f}%")


if __name__ == "__main__":
    for nw, nh in [(3, 3), (5, 5), (7, 7)]:
        run(nw, nh, n_arc=30)
