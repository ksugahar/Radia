"""
ngbem PEEC: Material coupling (ferrite + shield) on a flat-plate conductor.

Demonstrates:
  1. Air-only baseline (MQS loop impedance)
  2. Ferrite core coupling (analytical image method: Delta_L = L / (mu_r + 1))
  3. Al shield coupling (BEM+SIBC eddy current solve: Delta_Z)
  4. Combined ferrite + shield

Requires: ngsolve>=6.2.2601, ngsolve-ngsbem, matplotlib
"""
import sys, os
import numpy as np
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from ngbem_peec import NGBEMPEECSolver, create_plate_mesh, MU_0
from ngbem_eddy import ShieldBEMSIBC


def extract_edge_geometry(mesh):
    """Extract edge midpoints, directions, lengths from NGSolve surface mesh."""
    verts = {}
    for v in mesh.vertices:
        verts[v.nr] = np.array(v.point)
    centers, directions, lengths = [], [], []
    for edge in mesh.edges:
        v0, v1 = [v.nr for v in edge.vertices]
        p0, p1 = verts[v0], verts[v1]
        diff = p1 - p0
        length = np.linalg.norm(diff)
        centers.append(0.5 * (p0 + p1))
        directions.append(diff / length)
        lengths.append(length)
    return {
        'centers': np.array(centers),
        'directions': np.array(directions),
        'lengths': np.array(lengths),
    }


def compute_delta_L(L_air, mu_r):
    """Image method coupling: Delta_L = L_air / (mu_r + 1)."""
    return L_air / (mu_r + 1.0)


def main():
    # --- Geometry ---
    width = 0.01       # 10 mm
    height = 0.01      # 10 mm
    maxh = 0.003       # ~3 mm element size
    thickness = 35e-6   # 35 um copper
    sigma = 5.8e7       # Cu conductivity [S/m]

    print("=== ngbem PEEC: Material coupling demo ===")
    print(f"Conductor: {width*1e3:.0f} x {height*1e3:.0f} mm, "
          f"t = {thickness*1e6:.0f} um Cu")

    # --- Mesh and assembly ---
    mesh = create_plate_mesh(width, height, maxh, label="conductor")
    solver = NGBEMPEECSolver(mesh, conductor_label="conductor",
                              sigma=sigma, thickness=thickness,
                              order=0, intorder=5)
    solver.assemble()
    print(f"Conductor: {solver.n_loop} loop DOFs, {mesh.ne} elements")

    # --- Case 1: Air only ---
    freqs = np.array([100.0, 1e3, 10e3, 100e3])
    Z_air = solver.solve_frequency(freqs, mode='mqs')

    # --- Case 2: +Ferrite (analytical) ---
    mu_r = 1000.0
    Delta_L = compute_delta_L(solver.L, mu_r)
    L_with_core = solver.L + Delta_L * (mu_r - 1.0)

    Z_core = np.zeros(len(freqs), dtype=complex)
    for k, f in enumerate(freqs):
        omega = 2 * np.pi * f
        Z_branch = np.diag(solver.R_loop.astype(complex)) + 1j * omega * L_with_core
        Y = np.linalg.inv(Z_branch)
        e = np.ones(solver.n_loop) / solver.n_loop
        Z_core[k] = 1.0 / (e @ Y @ e)

    # --- Case 3: +Shield (BEM) ---
    from netgen.occ import Box, Pnt, OCCGeometry
    from ngsolve import Mesh as NGMesh
    from ngsolve import TaskManager

    print("\nAssembling shield BEM solver...")
    t0 = time.perf_counter()
    shield_plate = Box(Pnt(-0.006, -0.006, 0.005), Pnt(0.006, 0.006, 0.0055))
    shield_plate.solids.name = "conductor"
    shield_plate.faces.name = "surface"
    with TaskManager():
        shield_mesh = NGMesh(OCCGeometry(shield_plate).GenerateMesh(maxh=0.003))

        shield = ShieldBEMSIBC(shield_mesh, sigma=3.7e7)
        shield.assemble(intorder=4)
        t_shield = time.perf_counter() - t0
        print(f"Shield: {shield._loop.n_loops} loops, {shield._loop.n_active} active DOFs")
        print(f"Shield assembly: {t_shield*1e3:.0f} ms")

        edge_geom = extract_edge_geometry(mesh)
        topo_dict = {
            'segment_centers': edge_geom['centers'],
            'segment_directions': edge_geom['directions'],
            'segment_lengths': edge_geom['lengths'],
        }

        Z_shield = np.zeros(len(freqs), dtype=complex)
        for k, f in enumerate(freqs):
            omega = 2 * np.pi * f
            Delta_Z = shield.compute_impedance_matrix(f, topo_dict)
            Z_branch = np.diag(solver.R_loop.astype(complex)) + 1j * omega * solver.L + Delta_Z
            Y = np.linalg.inv(Z_branch)
            e = np.ones(solver.n_loop) / solver.n_loop
            Z_shield[k] = 1.0 / (e @ Y @ e)

        # --- Case 4: +Both ---
        Z_both = np.zeros(len(freqs), dtype=complex)
        for k, f in enumerate(freqs):
            omega = 2 * np.pi * f
            Delta_Z = shield.compute_impedance_matrix(f, topo_dict)
            Z_branch = np.diag(solver.R_loop.astype(complex)) + 1j * omega * L_with_core + Delta_Z
            Y = np.linalg.inv(Z_branch)
            e = np.ones(solver.n_loop) / solver.n_loop
            Z_both[k] = 1.0 / (e @ Y @ e)

        # --- Results ---
        print(f"\n{'Case':>20s}  {'f [Hz]':>10s}  {'L [nH]':>10s}  "
              f"{'R [mOhm]':>10s}  {'dL [nH]':>10s}")
        print("=" * 70)
        for k, f in enumerate(freqs):
            omega = 2 * np.pi * f
            L_a = np.imag(Z_air[k]) / omega * 1e9
            L_c = np.imag(Z_core[k]) / omega * 1e9
            L_s = np.imag(Z_shield[k]) / omega * 1e9
            L_b = np.imag(Z_both[k]) / omega * 1e9
            R_a = np.real(Z_air[k]) * 1e3
            R_c = np.real(Z_core[k]) * 1e3
            R_s = np.real(Z_shield[k]) * 1e3
            R_b = np.real(Z_both[k]) * 1e3
            print(f"{'Air':>20s}  {f:10.0f}  {L_a:10.2f}  {R_a:10.4f}  {'--':>10s}")
            print(f"{'+ Ferrite':>20s}  {f:10.0f}  {L_c:10.2f}  {R_c:10.4f}  {L_c-L_a:+10.2f}")
            print(f"{'+ Shield':>20s}  {f:10.0f}  {L_s:10.2f}  {R_s:10.4f}  {L_s-L_a:+10.2f}")
            print(f"{'+ Both':>20s}  {f:10.0f}  {L_b:10.2f}  {R_b:10.4f}  {L_b-L_a:+10.2f}")
            if k < len(freqs) - 1:
                print("-" * 70)

        # --- Physics checks ---
        L_air_vals = np.imag(Z_air) / (2 * np.pi * freqs) * 1e9
        L_core_vals = np.imag(Z_core) / (2 * np.pi * freqs) * 1e9
        L_shield_vals = np.imag(Z_shield) / (2 * np.pi * freqs) * 1e9
        R_shield_vals = np.real(Z_shield) * 1e3
        R_air_vals = np.real(Z_air) * 1e3

        dL_core = L_core_vals - L_air_vals
        dL_shield = L_shield_vals - L_air_vals
        dR_shield = R_shield_vals - R_air_vals

        print(f"\n--- Physics checks ---")
        print(f"Ferrite dL > 0:  {np.all(dL_core > 0)}  "
              f"(min dL = {np.min(dL_core):+.2f} nH)")
        print(f"Shield  dL < 0:  {np.all(dL_shield < 0)}  "
              f"(min dL = {np.min(dL_shield):+.2f} nH)")
        print(f"Shield  dR > 0:  {np.all(dR_shield > 0)}  "
              f"(min dR = {np.min(dR_shield):+.4f} mOhm)")


if __name__ == '__main__':
    main()
