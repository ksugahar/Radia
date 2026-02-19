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

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/radia'))

import radia as rad
from peec_matrices import PEECBuilder
from peec_coupled import CoupledPEECSolver
from fasthenry_parser import FastHenryParser

MU_0 = 4e-7 * np.pi


def demo_python_api():
    """Demo 1: Python API - single wire with/without core."""
    print("=" * 60)
    print("Demo 1: Python API - Wire with magnetic core")
    print("=" * 60)
    print()

    # --- Air-core case ---
    rad.UtiDelAll()
    rad.FldUnits('m')

    builder = PEECBuilder()
    n1 = builder.add_node_at(0, 0, 0)
    n2 = builder.add_node_at(0.1, 0, 0)
    builder.add_connected_segment(n1, n2, 1e-3, 1e-3, sigma=5.8e7)
    builder.add_port(n1, n2)
    topo = builder.build_topology()

    from peec_topology import PEECCircuitSolver
    solver_air = PEECCircuitSolver(topo)
    Z_air = solver_air.compute_port_impedance(1e6)
    L_air = np.imag(Z_air) / (2 * np.pi * 1e6)
    R_air = np.real(Z_air)

    print("  Geometry: 100mm wire (1mm x 1mm cross-section, copper)")
    print(f"  Air-core:  R = {R_air*1e3:.4f} mOhm,  L = {L_air*1e9:.2f} nH")

    # --- With magnetic core ---
    rad.UtiDelAll()
    rad.FldUnits('m')

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
        rad.FldUnits('m')

        builder = PEECBuilder()
        n1 = builder.add_node_at(0, 0, 0)
        n2 = builder.add_node_at(0.1, 0, 0)
        builder.add_connected_segment(n1, n2, 1e-3, 1e-3, sigma=5.8e7)
        builder.add_port(n1, n2)
        topo = builder.build_topology()

        if mu_r <= 1:
            # No magnetic material - air-core only
            from peec_topology import PEECCircuitSolver
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
        rad.FldUnits('m')

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


def demo_gmsh_mesh():
    """Demo 5: Magnetic core from GMSH mesh file."""
    print("=" * 60)
    print("Demo 5: Magnetic core from GMSH .msh file")
    print("=" * 60)
    print()

    # --- Generate a simple hex mesh using GMSH Python API ---
    mesh_file = os.path.join(os.path.dirname(__file__),
                             'gmsh_models', 'ferrite_core_hex.msh')
    _generate_ferrite_core_mesh(mesh_file)

    # --- Direct Python API with gmsh_mesh_import ---
    rad.UtiDelAll()
    rad.FldUnits('m')

    from gmsh_mesh_import import gmsh_to_radia, get_mesh_info

    # Show mesh info
    info = get_mesh_info(mesh_file)
    print(f"  Mesh info: {info['n_hex8']} hex8, {info['n_tet4']} tet4, "
          f"{info['n_total_volume']} total volume elements")
    print(f"  Bounding box: {info['bbox_min']} to {info['bbox_max']}")
    print()

    # Build PEEC topology
    builder = PEECBuilder()
    n1 = builder.add_node_at(0, 0, 0)
    n2 = builder.add_node_at(0.1, 0, 0)
    builder.add_connected_segment(n1, n2, 1e-3, 1e-3, sigma=5.8e7)
    builder.add_port(n1, n2)
    topo = builder.build_topology()

    # Import GMSH mesh as magnetic core
    core = gmsh_to_radia(mesh_file, mu_r=1000)

    solver = CoupledPEECSolver(topo, [core])
    solver.compute_coupling_matrix()

    Z = solver.compute_port_impedance(1e6)
    L = np.imag(Z) / (2 * np.pi * 1e6)
    Delta_L = solver.Delta_L[0, 0]

    print(f"  GMSH mesh core: L = {L*1e9:.2f} nH, Delta_L = {Delta_L*1e9:.4f} nH")

    # --- Same problem via FastHenry .magnetic type=mesh ---
    rad.UtiDelAll()

    inp_mesh = f"""\
* Inductor with ferrite core from GMSH mesh
.Units m
N1 x=0 y=0 z=0
N2 x=0.1 y=0 z=0
E1 N1 N2 w=1e-3 h=1e-3 sigma=5.8e7
.external N1 N2

.magnetic
  type=mesh
  file={os.path.basename(mesh_file)}
  mu_r=1000
.endmagnetic

.freq fmin=1e6 fmax=1e6 ndec=1
.end
"""
    parser = FastHenryParser()
    # Parse from directory containing mesh file
    inp_file = os.path.join(os.path.dirname(mesh_file),
                            '_temp_demo_mesh.inp')
    with open(inp_file, 'w') as f:
        f.write(inp_mesh)
    parser.parse_file(inp_file)
    result = parser.solve()
    os.remove(inp_file)

    L_fh = result['L'][0] if len(result['L']) > 0 else 0.0
    print(f"  FastHenry type=mesh: L = {L_fh*1e9:.2f} nH")
    print()

    return L, Delta_L


def _generate_ferrite_core_mesh(filename):
    """Generate a simple hex mesh for ferrite core using GMSH API."""
    try:
        import gmsh
    except ImportError:
        # Fallback: write a minimal GMSH 2.2 ASCII file manually
        _write_simple_hex_mesh(filename)
        return

    gmsh.initialize()
    gmsh.option.setNumber('General.Terminal', 0)
    gmsh.model.add('ferrite_core')

    # Create box: 60mm x 10mm x 10mm, centered at (0.05, 0.01, 0)
    # This matches the box in Demo 1-4
    x0, y0, z0 = 0.02, 0.005, -0.005
    dx, dy, dz = 0.06, 0.01, 0.01

    gmsh.model.occ.addBox(x0, y0, z0, dx, dy, dz)
    gmsh.model.occ.synchronize()

    # Use transfinite meshing for structured hex
    volumes = gmsh.model.getEntities(3)
    surfaces = gmsh.model.getEntities(2)
    curves = gmsh.model.getEntities(1)

    # Set number of elements along each edge
    for c in curves:
        length = gmsh.model.occ.getMass(1, c[1])
        if abs(length - dx) < 1e-6:
            n_div = 6  # 6 along x (60mm)
        elif abs(length - dy) < 1e-6:
            n_div = 2  # 2 along y (10mm)
        elif abs(length - dz) < 1e-6:
            n_div = 2  # 2 along z (10mm)
        else:
            n_div = 2
        gmsh.model.mesh.setTransfiniteCurve(c[1], n_div + 1)

    for s in surfaces:
        gmsh.model.mesh.setTransfiniteSurface(s[1])
        gmsh.model.mesh.setRecombine(2, s[1])

    for v in volumes:
        gmsh.model.mesh.setTransfiniteVolume(v[1])

    gmsh.model.mesh.generate(3)

    # Save as GMSH 2.2 ASCII
    gmsh.option.setNumber('Mesh.MshFileVersion', 2.2)
    gmsh.option.setNumber('Mesh.Binary', 0)
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    gmsh.write(filename)
    gmsh.finalize()


def _write_simple_hex_mesh(filename):
    """Write a minimal GMSH 2.2 hex mesh (2x1x1 = 2 hex elements) as fallback."""
    os.makedirs(os.path.dirname(filename), exist_ok=True)

    # 2 hex elements along x-axis
    # Core: 60mm x 10mm x 10mm, from (0.02, 0.005, -0.005) to (0.08, 0.015, 0.005)
    x0, y0, z0 = 0.02, 0.005, -0.005
    dx, dy, dz = 0.06, 0.01, 0.01
    nx = 2

    nodes = []
    nid = 1
    for ix in range(nx + 1):
        x = x0 + ix * dx / nx
        for iy in range(2):
            y = y0 + iy * dy
            for iz in range(2):
                z = z0 + iz * dz
                nodes.append((nid, x, y, z))
                nid += 1

    # Build hex elements: node ordering matches GMSH convention
    elements = []
    eid = 1
    for ix in range(nx):
        # Bottom-layer node indices (iy=0, iz=0 and iz=1)
        # Node numbering: ix varies slowest, then iy, then iz
        n = lambda iix, iiy, iiz: 1 + iix * 4 + iiy * 2 + iiz
        e_nodes = [
            n(ix, 0, 0), n(ix+1, 0, 0), n(ix+1, 1, 0), n(ix, 1, 0),
            n(ix, 0, 1), n(ix+1, 0, 1), n(ix+1, 1, 1), n(ix, 1, 1),
        ]
        elements.append((eid, e_nodes))
        eid += 1

    with open(filename, 'w') as f:
        f.write('$MeshFormat\n2.2 0 8\n$EndMeshFormat\n')
        f.write(f'$Nodes\n{len(nodes)}\n')
        for nid, x, y, z in nodes:
            f.write(f'{nid} {x:.10g} {y:.10g} {z:.10g}\n')
        f.write('$EndNodes\n')
        f.write(f'$Elements\n{len(elements)}\n')
        for eid, enodes in elements:
            node_str = ' '.join(str(n) for n in enodes)
            f.write(f'{eid} 5 2 0 0 {node_str}\n')  # type=5 (Hex8), 2 tags
        f.write('$EndElements\n')


if __name__ == '__main__':
    print()
    print("Coupled PEEC+MMM: Inductance Enhancement by Magnetic Core")
    print("=" * 60)
    print()

    demo_python_api()
    demo_mu_r_sweep()
    demo_distance_sweep()
    demo_fasthenry_input()
    demo_gmsh_mesh()

    print("All demos completed.")
