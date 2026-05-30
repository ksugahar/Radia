"""
Test: ShieldBEMSIBC with mu_r > 1
===================================

Verify that ShieldBEMSIBC correctly handles magnetic conductors (mu_r > 1).

Physical expectations:
  - Higher mu_r -> smaller skin depth: delta = sqrt(2/(omega*mu*sigma))
  - Higher mu_r -> larger |Zs|: Zs = (1+j)/(sigma*delta)
  - Loss should change predictably with mu_r
  - The solver should converge without errors

Test cases:
  1. mu_r=1   (aluminum baseline)
  2. mu_r=100  (magnetic steel)
  3. mu_r=1000 (soft iron)
"""

import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))

MU_0 = 4e-7 * np.pi


def main():
    try:
        from netgen.occ import Box, Pnt, OCCGeometry
        from ngsolve import Mesh, BND, TaskManager
    except ImportError:
        print("SKIP: NGSolve/Netgen not available")
        return

    try:
        from ngsolve.bem import MaxwellSingleLayerPotentialOperator
    except ImportError:
        print("SKIP: ngsolve.bem not available")
        return

    from ngbem_eddy import ShieldBEMSIBC, _biot_savart_A

    print("=" * 70)
    print("Test: ShieldBEMSIBC with mu_r > 1")
    print("=" * 70)

    # --- Problem setup ---
    width, height, depth = 0.1, 0.1, 0.002  # 100x100x2 mm plate
    sigma = 1e6  # S/m (mild steel conductivity)
    freq = 10e3  # 10 kHz
    wire_height = 0.01

    maxh = 0.025
    block = Box(Pnt(-width/2, -height/2, -depth/2),
                Pnt(width/2, height/2, depth/2))
    block.solids.name = "conductor"
    block.faces.name = "surface"
    geo = OCCGeometry(block)
    with TaskManager():
        ngmesh = geo.GenerateMesh(maxh=maxh)
        mesh = Mesh(ngmesh)

    n_bnd = sum(1 for _ in mesh.Elements(BND))
    print(f"\nMesh: {n_bnd} boundary faces")

    # Wire excitation
    wire_center = np.array([0, 0, wire_height + depth/2])
    wire_dir = np.array([1.0, 0, 0])
    wire_length = 0.2

    def A_inc(points):
        return _biot_savart_A(points, wire_center, wire_dir, wire_length)

    # --- Test with multiple mu_r values ---
    mu_r_values = [1.0, 100.0, 1000.0]
    results = {}

    print(f"\nFrequency: {freq/1e3:.0f} kHz")
    print(f"sigma = {sigma:.1e} S/m")
    print()
    print(f"{'mu_r':>8s}  {'delta(mm)':>10s}  {'|Zs|':>12s}  "
          f"{'P_loss(W)':>12s}  {'Status':>8s}")
    print(f"{'-'*8}  {'-'*10}  {'-'*12}  {'-'*12}  {'-'*8}")

    for mu_r in mu_r_values:
        omega = 2 * np.pi * freq
        mu = mu_r * MU_0
        delta = np.sqrt(2.0 / (omega * mu * sigma))

        solver = ShieldBEMSIBC(mesh, sigma, mu_r=mu_r)
        solver.assemble(intorder=4)
        result = solver.solve(freq, A_inc_func=A_inc)
        P_loss = solver.compute_loss()

        results[mu_r] = {
            'delta': delta,
            'Zs': result['Zs'],
            'P_loss': P_loss,
            'J_max': result['J_loop_max'],
        }

        status = "OK" if P_loss > 0 else "FAIL"
        print(f"{mu_r:8.0f}  {delta*1e3:10.4f}  {abs(result['Zs']):12.4e}  "
              f"{P_loss:12.4e}  {status:>8s}")

    # --- Validation checks ---
    print("\n--- Validation ---")
    all_pass = True

    # Check 1: Skin depth decreases with mu_r
    deltas = [results[mu_r]['delta'] for mu_r in mu_r_values]
    test_delta_decrease = all(deltas[i] > deltas[i+1]
                              for i in range(len(deltas)-1))
    print(f"  [{'PASS' if test_delta_decrease else 'FAIL'}] "
          f"Skin depth decreases with mu_r")
    if not test_delta_decrease:
        all_pass = False

    # Check 2: |Zs| increases with mu_r
    zs_abs = [abs(results[mu_r]['Zs']) for mu_r in mu_r_values]
    test_zs_increase = all(zs_abs[i] < zs_abs[i+1]
                            for i in range(len(zs_abs)-1))
    print(f"  [{'PASS' if test_zs_increase else 'FAIL'}] "
          f"|Zs| increases with mu_r")
    if not test_zs_increase:
        all_pass = False

    # Check 3: All losses positive
    test_loss_positive = all(results[mu_r]['P_loss'] > 0
                              for mu_r in mu_r_values)
    print(f"  [{'PASS' if test_loss_positive else 'FAIL'}] "
          f"All losses positive")
    if not test_loss_positive:
        all_pass = False

    # Check 4: delta ratio matches theory
    # delta(mu_r) / delta(1) = 1/sqrt(mu_r)
    for mu_r in [100.0, 1000.0]:
        ratio_actual = results[mu_r]['delta'] / results[1.0]['delta']
        ratio_theory = 1.0 / np.sqrt(mu_r)
        rel_err = abs(ratio_actual - ratio_theory) / ratio_theory
        test_ratio = rel_err < 0.01  # 1% tolerance
        print(f"  [{'PASS' if test_ratio else 'FAIL'}] "
              f"delta ratio mu_r={mu_r:.0f}: "
              f"actual={ratio_actual:.6f}, "
              f"theory={ratio_theory:.6f}, err={rel_err*100:.2f}%")
        if not test_ratio:
            all_pass = False

    # Check 5: |Zs| ratio matches theory
    # |Zs| = sqrt(omega*mu/(2*sigma)) * sqrt(2)
    # |Zs(mu_r)| / |Zs(1)| = sqrt(mu_r)
    for mu_r in [100.0, 1000.0]:
        ratio_actual = abs(results[mu_r]['Zs']) / abs(results[1.0]['Zs'])
        ratio_theory = np.sqrt(mu_r)
        rel_err = abs(ratio_actual - ratio_theory) / ratio_theory
        test_ratio = rel_err < 0.01
        print(f"  [{'PASS' if test_ratio else 'FAIL'}] "
              f"|Zs| ratio mu_r={mu_r:.0f}: "
              f"actual={ratio_actual:.4f}, "
              f"theory={ratio_theory:.4f}, err={rel_err*100:.2f}%")
        if not test_ratio:
            all_pass = False

    # --- Frequency sweep for mu_r=100 ---
    print(f"\n--- Frequency sweep for mu_r=100 ---")
    freqs = [1e3, 10e3, 100e3]
    solver_100 = ShieldBEMSIBC(mesh, sigma, mu_r=100.0)
    solver_100.assemble(intorder=4)

    losses_100 = []
    for f in freqs:
        delta_f = np.sqrt(2.0 / (2*np.pi*f * 100.0*MU_0 * sigma))
        solver_100.solve(f, A_inc_func=A_inc)
        P_f = solver_100.compute_loss()
        losses_100.append(P_f)
        print(f"  f={f/1e3:8.1f} kHz, delta={delta_f*1e3:.4f} mm, "
              f"P={P_f:.4e} W")

    test_monotonic = all(losses_100[i] < losses_100[i+1]
                          for i in range(len(losses_100)-1))
    print(f"  [{'PASS' if test_monotonic else 'FAIL'}] "
          f"Loss monotonic with freq (mu_r=100)")
    if not test_monotonic:
        all_pass = False

    # --- Summary ---
    print("\n" + "=" * 70)
    print(f"Overall: {'ALL PASSED' if all_pass else 'SOME FAILED'}")
    print("=" * 70)


if __name__ == '__main__':
    main()
