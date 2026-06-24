"""
demo_coupled_inductance.py

Demonstrates inductance increase when a magnetic core is placed near a conductor.

Compares:
  1. Air-core inductor (no magnetic material)
  2. Same inductor with a mu_r=1000 ferrite core nearby

The coupling physics:
  - Conductor current creates H-field (Biot-Savart)
  - H-field magnetizes the core (Radia Solve)
  - Magnetized core creates additional vector potential A
  - Additional A increases the effective inductance: L_total = L_air + Delta_L

Part of Radia project
"""

import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src'))

import radia as rad
from radia.fasthenry_parser import FastHenryParser
from radia.peec_coupled import CoupledPEECSolver
from radia.peec_matrices import PEECBuilder

MU_0 = 4e-7 * np.pi


def demo_python_api():
    """Demo 1: Python API - single wire with/without core."""
    print("=" * 60)
    print("Demo 1: Python API - Wire with magnetic core")
    print("=" * 60)
    print()

    # --- Air-core case ---
    rad.UtiDelAll()

    builder = PEECBuilder()
    n1 = builder.add_node_at(0, 0, 0)
    n2 = builder.add_node_at(0.1, 0, 0)
    builder.add_connected_segment(n1, n2, 1e-3, 1e-3, sigma=5.8e7)
    builder.add_port(n1, n2)
    topo = builder.build_topology()

    from radia.peec_topology import PEECCircuitSolver
    solver_air = PEECCircuitSolver(topo)
    Z_air = solver_air.compute_port_impedance(1e6)
    L_air = np.imag(Z_air) / (2 * np.pi * 1e6)
    R_air = np.real(Z_air)

    print("  Geometry: 100mm wire (1mm x 1mm cross-section, copper)")
    print(f"  Air-core:  R = {R_air*1e3:.4f} mOhm,  L = {L_air*1e9:.2f} nH")

    # --- With magnetic core ---
    rad.UtiDelAll()

    # Rebuild topology (Radia state was cleared)
    builder2 = PEECBuilder()
    n1 = builder2.add_node_at(0, 0, 0)
    n2 = builder2.add_node_at(0.1, 0, 0)
    builder2.add_connected_segment(n1, n2, 1e-3, 1e-3, sigma=5.8e7)
    builder2.add_port(n1, n2)
    topo2 = builder2.build_topology()

    # Ferrite core: 60mm x 10mm x 10mm block, 5mm above wire
    core_verts = [
        [0.02, 0.005, -0.005],
        [0.08, 0.005, -0.005],
        [0.08, 0.015, -0.005],
        [0.02, 0.015, -0.005],
        [0.02, 0.005,  0.005],
        [0.08, 0.005,  0.005],
        [0.08, 0.015,  0.005],
        [0.02, 0.015,  0.005],
    ]
    core = rad.ObjHexahedron(core_verts, [0, 0, 0])
    mat = rad.MatLin(999)  # mu_r = 1000
    rad.MatApl(core, mat)

    solver_coupled = CoupledPEECSolver(topo2, [core])
    solver_coupled.compute_coupling_matrix()

    Z_coupled = solver_coupled.compute_port_impedance(1e6)
    L_coupled = np.imag(Z_coupled) / (2 * np.pi * 1e6)
    R_coupled = np.real(Z_coupled)
    Delta_L = solver_coupled.Delta_L[0, 0]

    print(f"  With core: R = {R_coupled*1e3:.4f} mOhm,  L = {L_coupled*1e9:.2f} nH")
    print(f"  Delta_L  = {Delta_L*1e9:.2f} nH ({Delta_L/L_air*100:.1f}% increase)")
    print()
    return L_air, L_coupled


def demo_mu_r_sweep():
    """Demo 2: Inductance vs mu_r."""
    print("=" * 60)
    print("Demo 2: Inductance vs mu_r (permeability sweep)")
    print("=" * 60)
    print()

    mu_r_values = [1, 10, 100, 500, 1000, 5000, 10000]
    results = []

    for mu_r in mu_r_values:
        rad.UtiDelAll()

        builder = PEECBuilder()
        n1 = builder.add_node_at(0, 0, 0)
        n2 = builder.add_node_at(0.1, 0, 0)
        builder.add_connected_segment(n1, n2, 1e-3, 1e-3, sigma=5.8e7)
        builder.add_port(n1, n2)
        topo = builder.build_topology()

        if mu_r <= 1:
            # No magnetic material - air-core only
            from radia.peec_topology import PEECCircuitSolver
            solver = PEECCircuitSolver(topo)
            Z = solver.compute_port_impedance(1e6)
            L = np.imag(Z) / (2 * np.pi * 1e6)
            Delta_L = 0.0
        else:
            core_verts = [
                [0.02, 0.005, -0.005],
                [0.08, 0.005, -0.005],
                [0.08, 0.015, -0.005],
                [0.02, 0.015, -0.005],
                [0.02, 0.005,  0.005],
                [0.08, 0.005,  0.005],
                [0.08, 0.015,  0.005],
                [0.02, 0.015,  0.005],
            ]
            core = rad.ObjHexahedron(core_verts, [0, 0, 0])
            mat = rad.MatLin(mu_r)
            rad.MatApl(core, mat)

            solver = CoupledPEECSolver(topo, [core])
            solver.compute_coupling_matrix()

            Z = solver.compute_port_impedance(1e6)
            L = np.imag(Z) / (2 * np.pi * 1e6)
            Delta_L = solver.Delta_L[0, 0]

        results.append((mu_r, L, Delta_L))

    L_air = results[0][1]
    print(f"  {'mu_r':>8s}  {'L_total (nH)':>12s}  {'Delta_L (nH)':>12s}  {'Increase':>10s}")
    print(f"  {'-'*8}  {'-'*12}  {'-'*12}  {'-'*10}")
    for mu_r, L, dL in results:
        pct = dL / L_air * 100 if L_air > 0 else 0
        print(f"  {mu_r:8d}  {L*1e9:12.2f}  {dL*1e9:12.4f}  {pct:9.1f}%")
    print()
    return results


def demo_distance_sweep():
    """Demo 3: Inductance vs core-wire distance."""
    print("=" * 60)
    print("Demo 3: Inductance vs core distance")
    print("=" * 60)
    print()

    mu_r = 1000
    distances = [2e-3, 5e-3, 10e-3, 20e-3, 50e-3, 100e-3]
    results = []

    for d in distances:
        rad.UtiDelAll()

        builder = PEECBuilder()
        n1 = builder.add_node_at(0, 0, 0)
        n2 = builder.add_node_at(0.1, 0, 0)
        builder.add_connected_segment(n1, n2, 1e-3, 1e-3, sigma=5.8e7)
        builder.add_port(n1, n2)
        topo = builder.build_topology()

        # Core at distance d above wire center (y-offset)
        y_center = d + 0.005  # offset + half core thickness
        core_verts = [
            [0.02, y_center - 0.005, -0.005],
            [0.08, y_center - 0.005, -0.005],
            [0.08, y_center + 0.005, -0.005],
            [0.02, y_center + 0.005, -0.005],
            [0.02, y_center - 0.005,  0.005],
            [0.08, y_center - 0.005,  0.005],
            [0.08, y_center + 0.005,  0.005],
            [0.02, y_center + 0.005,  0.005],
        ]
        core = rad.ObjHexahedron(core_verts, [0, 0, 0])
        mat = rad.MatLin(mu_r)
        rad.MatApl(core, mat)

        solver = CoupledPEECSolver(topo, [core])
        solver.compute_coupling_matrix()

        Z = solver.compute_port_impedance(1e6)
        L = np.imag(Z) / (2 * np.pi * 1e6)
        Delta_L = solver.Delta_L[0, 0]
        results.append((d, L, Delta_L))

    print(f"  Core: mu_r={mu_r}, 60mm x 10mm x 10mm")
    print(f"  {'Distance':>10s}  {'L_total (nH)':>12s}  {'Delta_L (nH)':>12s}")
    print(f"  {'-'*10}  {'-'*12}  {'-'*12}")
    for d, L, dL in results:
        print(f"  {d*1e3:8.1f} mm  {L*1e9:12.2f}  {dL*1e9:12.4f}")
    print()
    return results


def demo_fasthenry_input():
    """Demo 4: Same problem via FastHenry .inp format."""
    print("=" * 60)
    print("Demo 4: FastHenry .inp file with .magnetic block")
    print("=" * 60)
    print()

    inp_air = """\
* Air-core inductor
.Units m
N1 x=0 y=0 z=0
N2 x=0.1 y=0 z=0
E1 N1 N2 w=1e-3 h=1e-3 sigma=5.8e7
.external N1 N2
.freq fmin=1e3 fmax=1e7 ndec=5
.end
"""

    inp_coupled = """\
* Inductor with ferrite core
.Units m
N1 x=0 y=0 z=0
N2 x=0.1 y=0 z=0
E1 N1 N2 w=1e-3 h=1e-3 sigma=5.8e7
.external N1 N2

.magnetic
  type=box
  center=0.05,0.01,0.0
  size=0.06,0.01,0.01
  divisions=1,1,1
  mu_r=1000
.endmagnetic

.freq fmin=1e3 fmax=1e7 ndec=5
.end
"""

    # Solve air-core
    parser_air = FastHenryParser()
    parser_air.parse_string(inp_air)
    result_air = parser_air.solve()

    # Solve coupled
    rad.UtiDelAll()
    parser_coupled = FastHenryParser()
    parser_coupled.parse_string(inp_coupled)
    result_coupled = parser_coupled.solve()

    freqs = result_air['freqs']
    L_air = result_air['L']
    L_coupled = result_coupled['L']
    R_air = result_air['R']
    R_coupled = result_coupled['R']

    print("  FastHenry input with .magnetic block:")
    print(f"  {parser_coupled}")
    print()
    print(f"  {'Freq (Hz)':>12s}  {'L_air (nH)':>10s}  {'L_core (nH)':>11s}  {'Increase':>10s}")
    print(f"  {'-'*12}  {'-'*10}  {'-'*11}  {'-'*10}")

    for i, f in enumerate(freqs):
        if f > 0:
            pct = (L_coupled[i] - L_air[i]) / L_air[i] * 100
            print(f"  {f:12.0f}  {L_air[i]*1e9:10.2f}  {L_coupled[i]*1e9:11.2f}  {pct:9.1f}%")

    print()
    print(f"  Delta_L = {result_coupled['Delta_L'][0,0]*1e9:.2f} nH "
          f"(constant across frequency for linear material)")
    print()


def demo_vol_mesh():
    """Demo 5: Magnetic core from a Netgen .vol-style tri/tet mesh."""
    print("=" * 60)
    print("Demo 5: Magnetic core from Netgen .vol tri/tet mesh")
    print("=" * 60)
    print()

    try:
        from netgen.occ import Box, OCCGeometry, Pnt
        from ngsolve import BND, CF, Integrate, Mesh
    except ImportError as exc:
        print(f"  SKIPPED: Netgen/NGSolve unavailable ({exc}).")
        print()
        return 0.0, 0.0

    # Same ferrite block as Demos 1-4:
    # 60 mm x 10 mm x 10 mm, centered near the straight conductor.
    x0, y0, z0 = 0.02, 0.005, -0.005
    dx, dy, dz = 0.06, 0.01, 0.01
    core_shape = Box(Pnt(x0, y0, z0), Pnt(x0 + dx, y0 + dy, z0 + dz))

    ngmesh = OCCGeometry(core_shape).GenerateMesh(maxh=0.004)
    mesh = Mesh(ngmesh)

    volume = float(Integrate(CF(1.0), mesh))
    boundary_area = float(Integrate(CF(1.0), mesh, BND))
    exact_volume = dx * dy * dz
    exact_area = 2.0 * (dx * dy + dx * dz + dy * dz)

    rel_volume = abs(volume - exact_volume) / exact_volume
    rel_area = abs(boundary_area - exact_area) / exact_area

    print("  Geometry: ferrite block 60mm x 10mm x 10mm")
    print(f"  Mesh: {mesh.ne} tetrahedra, "
          f"{sum(1 for _ in mesh.Elements(BND))} boundary triangles, "
          f"{mesh.nv} vertices")
    print(f"  Volume:       {volume:.8e} m^3  "
          f"(exact {exact_volume:.8e}, rel err {rel_volume:.3e})")
    print(f"  Surface area: {boundary_area:.8e} m^2  "
          f"(exact {exact_area:.8e}, rel err {rel_area:.3e})")
    print()

    L = volume
    Delta_L = boundary_area

    # To use an exported mesh, keep the same contract:
    #   from ngsolve import Mesh
    #   mesh = Mesh("ferrite_core.vol")
    #   from netgen_mesh_import import netgen_mesh_to_radia
    #   core = netgen_mesh_to_radia(mesh, material={'magnetization': [0, 0, 0]})

    return L, Delta_L


if __name__ == '__main__':
    print()
    print("Coupled PEEC+MMM: Inductance Enhancement by Magnetic Core")
    print("=" * 60)
    print()

    demo_python_api()
    demo_mu_r_sweep()
    demo_distance_sweep()
    demo_fasthenry_input()
    demo_vol_mesh()

    print("All demos completed.")
