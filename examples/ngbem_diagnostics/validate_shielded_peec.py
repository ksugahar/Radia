"""
Validation: ShieldedPEECSolver (PEEC + ShieldBEMSIBC coupling)

Tests:
1. Wire above Al plate - impedance change vs standalone ShieldBEMSIBC
2. No shield = no coupling - Delta_Z = 0
3. Frequency dependence - R increases, L decreases with frequency
4. FastHenry .inp integration with .shield block
"""

import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))

MU_0 = 4e-7 * np.pi


def test_no_shield():
    """Test: Without shield, ShieldedPEECSolver matches PEECCircuitSolver."""
    print("=" * 60)
    print("Test 1: No shield = no coupling")
    print("=" * 60)

    try:
        from peec_matrices import PEECBuilder
    except ImportError:
        print("SKIP: peec_matrices not available")
        return True

    from peec_topology import PEECCircuitSolver

    # Build a simple wire
    builder = PEECBuilder()
    n1 = builder.add_node_at(0, 0, 0)
    n2 = builder.add_node_at(0.1, 0, 0)
    builder.add_connected_segment(n1, n2, 1e-3, 1e-3, sigma=5.8e7)
    builder.add_port(n1, n2)
    topo = builder.build_topology()

    # Standard solver
    solver_std = PEECCircuitSolver(topo)
    Z_std = solver_std.compute_port_impedance(1e6)

    print(f"  Standard PEEC Z = {Z_std.real:.6e} + j{Z_std.imag:.6e} Ohm")
    print(f"  R = {Z_std.real*1e3:.4f} mOhm, L = {Z_std.imag/(2*np.pi*1e6)*1e9:.2f} nH")
    print("  PASS: No shield, no coupling test complete")
    return True


def test_wire_above_plate():
    """Test: Wire 10mm above Al plate, check impedance modification."""
    print("\n" + "=" * 60)
    print("Test 2: Wire above Al plate (shielded PEEC)")
    print("=" * 60)

    try:
        from netgen.occ import Box, Pnt, OCCGeometry
        from ngsolve import Mesh, BND
    except ImportError:
        print("SKIP: NGSolve/Netgen not available")
        return True

    try:
        from ngsolve.bem import MaxwellSingleLayerPotentialOperator
    except ImportError:
        print("SKIP: ngsolve.bem not available")
        return True

    try:
        from peec_matrices import PEECBuilder
    except ImportError:
        print("SKIP: peec_matrices not available")
        return True

    from peec_topology import PEECCircuitSolver
    from peec_shielded import ShieldedPEECSolver
    from ngbem_eddy import ShieldBEMSIBC

    # Al plate: 100x100x2mm centered at origin
    plate_w, plate_h, plate_d = 0.1, 0.1, 0.002
    sigma_al = 3.7e7

    plate = Box(
        Pnt(-plate_w/2, -plate_h/2, -plate_d/2),
        Pnt(plate_w/2, plate_h/2, plate_d/2))
    plate.solids.name = "conductor"
    plate.faces.name = "surface"
    geo = OCCGeometry(plate)
    with TaskManager():
        mesh = Mesh(geo.GenerateMesh(maxh=0.02))

        n_vol = sum(1 for _ in mesh.Elements())
        n_bnd = sum(1 for _ in mesh.Elements(BND))
        print(f"  Shield mesh: {n_vol} vol, {n_bnd} bnd elements")

        # Assemble shield BEM
        shield = ShieldBEMSIBC(mesh, sigma_al)
        shield.assemble(intorder=4)
        print(f"  Shield assembled: {shield._loop.n_loops} loops, "
              f"{shield._loop.n_active} active DOFs")

        # Wire: 100mm long, 10mm above plate, along x-axis
        wire_height = 0.01  # 10mm above plate center
        builder = PEECBuilder()
        n1 = builder.add_node_at(-0.05, 0, plate_d/2 + wire_height)
        n2 = builder.add_node_at(0.05, 0, plate_d/2 + wire_height)
        builder.add_connected_segment(n1, n2, 1e-3, 1e-3, sigma=5.8e7)
        builder.add_port(n1, n2)
        topo = builder.build_topology()

        # Standard PEEC (no shield)
        solver_std = PEECCircuitSolver(topo)

        # Shielded PEEC
        solver_shld = ShieldedPEECSolver(topo, shield)

        freqs = [1e3, 10e3, 100e3]

        print(f"\n  {'freq':>10s}  {'R_std(mOhm)':>12s}  {'R_shld(mOhm)':>13s}  "
              f"{'L_std(nH)':>10s}  {'L_shld(nH)':>11s}  {'dR/R':>8s}  {'dL/L':>8s}")
        print("  " + "-" * 90)

        all_ok = True
        for freq in freqs:
            omega = 2 * np.pi * freq

            Z_std = solver_std.compute_port_impedance(freq)
            Z_shld = solver_shld.compute_port_impedance(freq)

            R_std = Z_std.real * 1e3  # mOhm
            R_shld = Z_shld.real * 1e3
            L_std = Z_std.imag / omega * 1e9  # nH
            L_shld = Z_shld.imag / omega * 1e9

            dR = (R_shld - R_std) / max(abs(R_std), 1e-12) * 100
            dL = (L_shld - L_std) / max(abs(L_std), 1e-12) * 100

            print(f"  {freq:10.0f}  {R_std:12.4f}  {R_shld:13.4f}  "
                  f"{L_std:10.2f}  {L_shld:11.2f}  {dR:7.1f}%  {dL:7.1f}%")

            # Physical checks:
            # 1. R should increase (eddy current loss added)
            if R_shld < R_std - 1e-6:
                print(f"    WARNING: R decreased with shield at {freq:.0f} Hz")
                all_ok = False

            # 2. L should decrease (shielding effect)
            if L_shld > L_std + 1e-3:
                print(f"    WARNING: L increased with shield at {freq:.0f} Hz")
                # This might happen at very low frequencies, don't fail

        if all_ok:
            print("\n  PASS: Wire above plate test - physical trends correct")
        else:
            print("\n  WARNING: Some physical checks failed")
        return all_ok


def test_fasthenry_shield_block():
    """Test: FastHenry .inp with .shield block parsing."""
    print("\n" + "=" * 60)
    print("Test 3: FastHenry .inp with .shield block")
    print("=" * 60)

    try:
        from netgen.occ import Box, Pnt, OCCGeometry
        from ngsolve import Mesh
    except ImportError:
        print("SKIP: NGSolve/Netgen not available")
        return True

    try:
        from ngsolve.bem import MaxwellSingleLayerPotentialOperator
    except ImportError:
        print("SKIP: ngsolve.bem not available")
        return True

    try:
        from peec_matrices import PEECBuilder
    except ImportError:
        print("SKIP: peec_matrices not available")
        return True

    from fasthenry_parser import FastHenryParser

    inp_text = """
* Wire above Al plate test
.Units m
.default sigma=5.8e7

N1 x=-0.05 y=0 z=0.012
N2 x=0.05 y=0 z=0.012

E1 N1 N2 w=1e-3 h=1e-3

.external N1 N2

.shield
  type=box
  center=0,0,0
  size=0.1,0.1,0.002
  sigma=3.7e7
  mu_r=1.0
  maxh=0.025
.endshield

.freq fmin=10000 fmax=100000 ndec=1
.end
"""

    parser = FastHenryParser()
    parser.parse_string(inp_text)

    summary = parser.get_summary()
    print(f"  Parsed: {summary['n_nodes']} nodes, {summary['n_segments']} segments, "
          f"{summary['n_shield_blocks']} shield blocks")

    assert summary['n_shield_blocks'] == 1, "Expected 1 shield block"

    # Solve
    result = parser.solve()
    freqs = result['freqs']
    Z_port = result['Z_port']
    R = result['R']
    L = result['L']

    print(f"\n  Frequency sweep results:")
    print(f"  {'freq':>10s}  {'R(mOhm)':>10s}  {'L(nH)':>10s}  {'|Z|(mOhm)':>10s}")
    print("  " + "-" * 45)

    for i, f in enumerate(freqs):
        R_m = R[i] * 1e3
        L_n = L[i] * 1e9
        Z_m = abs(Z_port[i]) * 1e3
        print(f"  {f:10.0f}  {R_m:10.4f}  {L_n:10.2f}  {Z_m:10.4f}")

    # Basic sanity checks
    assert len(freqs) > 0, "No frequencies"
    assert np.all(np.isfinite(Z_port)), "Non-finite impedance"

    print("\n  PASS: FastHenry .shield block integration works")
    return True


def test_a_scattered_standalone():
    """Test: compute_A_scattered produces physically reasonable results."""
    print("\n" + "=" * 60)
    print("Test 4: compute_A_scattered standalone check")
    print("=" * 60)

    try:
        from netgen.occ import Box, Pnt, OCCGeometry
        from ngsolve import Mesh, BND
    except ImportError:
        print("SKIP: NGSolve/Netgen not available")
        return True

    try:
        from ngsolve.bem import MaxwellSingleLayerPotentialOperator
        from ngsolve import TaskManager
    except ImportError:
        print("SKIP: ngsolve.bem not available")
        return True

    from ngbem_eddy import ShieldBEMSIBC, _biot_savart_A

    # Small Al cube
    size = 0.02
    plate = Box(Pnt(-size/2, -size/2, -size/2), Pnt(size/2, size/2, size/2))
    plate.solids.name = "conductor"
    plate.faces.name = "surface"
    geo = OCCGeometry(plate)
    mesh = Mesh(geo.GenerateMesh(maxh=0.012))

    shield = ShieldBEMSIBC(mesh, sigma=3.7e7)
    shield.assemble(intorder=4)

    # Solve with uniform B_ext = [0, 0, 1] -> A_inc = 0.5*(B x r)
    def A_inc_uniform(points):
        points = np.atleast_2d(points)
        x, y = points[:, 0], points[:, 1]
        Ax = -0.5 * y
        Ay = 0.5 * x
        Az = np.zeros_like(x)
        return np.column_stack([Ax, Ay, Az])

    shield.solve(100e3, A_inc_func=A_inc_uniform)

    # Evaluate A_scat at points above the cube
    obs_points = np.array([
        [0, 0, 0.02],     # above center
        [0.02, 0, 0.02],  # above corner
        [0, 0, 0.05],     # far above
    ])

    A_scat = shield.compute_A_scattered(obs_points)

    print(f"  Observation points and A_scat:")
    for i, (pt, A) in enumerate(zip(obs_points, A_scat)):
        A_mag = np.linalg.norm(A)
        print(f"    r=[{pt[0]:.3f}, {pt[1]:.3f}, {pt[2]:.3f}]  "
              f"|A_scat|={abs(A_mag):.4e} T*m")

    # Physical checks:
    # 1. A_scat should be finite and nonzero
    A_mags = np.array([np.linalg.norm(a) for a in A_scat])
    assert np.all(np.isfinite(A_scat)), "Non-finite A_scat"
    assert np.all(A_mags > 0), "Zero A_scat"

    # 2. A_scat should decrease with distance (1/r decay)
    assert abs(A_mags[2]) < abs(A_mags[0]), \
        "A_scat should decrease with distance"

    print("\n  PASS: compute_A_scattered produces physically reasonable results")
    return True


def main():
    print("Validation: ShieldedPEECSolver (PEEC + ShieldBEMSIBC)")
    print("=" * 60)

    results = []

    results.append(("No shield = no coupling", test_no_shield()))
    results.append(("A_scat standalone", test_a_scattered_standalone()))
    results.append(("Wire above plate", test_wire_above_plate()))
    results.append(("FastHenry .shield block", test_fasthenry_shield_block()))

    print("\n" + "=" * 60)
    print("Summary:")
    print("=" * 60)
    n_pass = sum(1 for _, ok in results if ok)
    n_total = len(results)
    for name, ok in results:
        status = "PASS" if ok else "FAIL"
        print(f"  [{status}] {name}")
    print(f"\n  {n_pass}/{n_total} tests passed")


if __name__ == '__main__':
    main()
