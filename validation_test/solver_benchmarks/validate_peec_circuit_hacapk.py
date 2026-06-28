"""validate_peec_circuit_hacapk.py -- Step 4 stage 4 validation

Compare PEECCircuitSolver(use_hacapk=True) against the existing dense
C++ MNA path on a helical coil topology (filaments in series, single
port).  Both paths must agree on the scalar port impedance Z_11 within
the combined nested BiCGSTAB tolerance.

Run:
    python validation_test/solver_benchmarks/validate_peec_circuit_hacapk.py
"""

from __future__ import annotations

import math
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))

import radia  # triggers MKL DLL directory registration

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src", "radia"))

from peec_matrices import PEECBuilder
from peec_topology import PEECCircuitSolver


def build_helical_topology(turns=6, seg_per_turn=12,
                           radius=0.050, pitch=0.006,
                           w=0.003, h=0.003, sigma=5.8e7):
    builder = PEECBuilder()
    nodes = []
    for i in range(turns * seg_per_turn + 1):
        t = 2.0 * math.pi * i / seg_per_turn
        z = pitch * (i / seg_per_turn)
        nodes.append(builder.add_node_at(radius * math.cos(t),
                                         radius * math.sin(t), z))
    for i in range(turns * seg_per_turn):
        builder.add_connected_segment(nodes[i], nodes[i + 1],
                                      w, h, sigma=sigma)
    builder.add_port(nodes[0], nodes[-1])
    return builder.build_topology(), len(nodes), turns * seg_per_turn


def main():
    turns, spt = 6, 12
    topo, n_nodes, n_loop = build_helical_topology(turns=turns,
                                                   seg_per_turn=spt)
    print(f"[helical] turns={turns} seg/turn={spt}  "
          f"n_loop={n_loop}  n_nodes={n_nodes}")

    freqs = np.array([0.0, 1e3, 10e3, 100e3, 1e6, 10e6])

    # ---- Dense reference (C++ MNASolver) ----
    solver_dense = PEECCircuitSolver(topo)
    t0 = time.perf_counter()
    Z_dense = np.array([solver_dense.compute_port_impedance(f)
                        for f in freqs])
    t_dense = time.perf_counter() - t0

    # ---- HACApK Python-MNA ----
    solver_hmat = PEECCircuitSolver(
        topo,
        use_hacapk=True,
        hacapk_params={
            "aca_eps": 1e-8,
            "leaf_size": 32,
            "eta": 2.0,
            "max_rank": 200,
            "print_level": 0,
        },
        bicgstab_outer_tol=1e-8,
        bicgstab_inner_tol=1e-10,
        bicgstab_maxiter=2000,
    )
    t0 = time.perf_counter()
    Z_hmat = np.array([solver_hmat.compute_port_impedance(f)
                       for f in freqs])
    t_hmat = time.perf_counter() - t0

    print(f"\n{'f [Hz]':>10s}  {'Re(Z_dense)':>12s} {'Im(Z_dense)':>12s}  "
          f"{'Re(Z_hmat)':>12s} {'Im(Z_hmat)':>12s}  {'rel_err':>10s}")
    max_err = 0.0
    for f, zd, zh in zip(freqs, Z_dense, Z_hmat):
        rel = abs(zh - zd) / max(abs(zd), 1e-30)
        max_err = max(max_err, rel)
        print(f"{f:10.3e}  {zd.real:12.4e} {zd.imag:12.4e}  "
              f"{zh.real:12.4e} {zh.imag:12.4e}  {rel:10.2e}")

    print(f"\ndense path:   {t_dense*1e3:.1f} ms total, "
          f"{t_dense*1e3/len(freqs):.2f} ms/freq")
    print(f"HACApK path:  {t_hmat*1e3:.1f} ms total, "
          f"{t_hmat*1e3/len(freqs):.2f} ms/freq")
    print(f"max rel err over {len(freqs)} freqs: {max_err:.2e}")

    # Correctness threshold: BiCGSTAB nested tol allows ~1e-6 drift
    tol = 1e-6
    if max_err > tol:
        print(f"[FAIL] max rel err {max_err:.2e} > tol {tol:.2e}")
        sys.exit(1)
    print(f"[PASS] max rel err {max_err:.2e} <= tol {tol:.2e}")


if __name__ == "__main__":
    main()
