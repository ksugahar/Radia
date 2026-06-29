"""
demo_schur_complement.py

Demonstrates the Schur complement solver for PEEC-MSC coupling.

Compares:
  1. Existing column-by-column approach (CoupledPEECSolver.compute_coupling_matrix)
  2. New Schur complement approach (SchurComplementSolver)

Both should produce the same Delta_L and coupled impedance.

The Schur complement approach preserves H-matrix structure for large MSC
systems, while the column-by-column approach requires n_peec full MSC solves
through Radia's object system.

Part of Radia project
"""

import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/radia'))

import radia as rad  # DLL loader for MKL etc.
from mmm_core import MMMBuilder, MMMSolver, MMMFieldComputer
from peec_msc_schur import SchurComplementSolver, biot_savart_finite_filament, MU_0


def create_test_problem():
    """
    Create a simple test problem: single wire segment + single hexahedral core.

    Returns geometry data needed by both approaches.
    """
    # PEEC: single wire along x-axis, 100mm long
    seg_p1 = np.array([[0.0, 0.0, 0.0]])
    seg_p2 = np.array([[0.1, 0.0, 0.0]])
    seg_centers = (seg_p1 + seg_p2) / 2.0
    seg_dirs = seg_p2 - seg_p1
    seg_lengths = np.array([np.linalg.norm(seg_dirs[0])])
    seg_dirs = seg_dirs / seg_lengths[:, None]

    # Wire parameters
    R_dc = np.array([0.1])  # 100 mOhm DC resistance

    # Air inductance (approximate: L = mu_0*l/(2*pi) * (ln(2*l/a) - 1))
    l = 0.1
    a = 0.5e-3  # effective radius for 1mm x 1mm cross section
    L_air_val = MU_0 * l / (2.0 * np.pi) * (np.log(2.0 * l / a) - 1.0)
    L_air = np.array([[L_air_val]])

    # MSC: Hexahedral magnetic core, 60mm x 10mm x 10mm, centered 10mm above wire
    # 8 vertices in standard hex order
    core_verts = np.array([
        [0.02, 0.005, -0.005],
        [0.08, 0.005, -0.005],
        [0.08, 0.015, -0.005],
        [0.02, 0.015, -0.005],
        [0.02, 0.005,  0.005],
        [0.08, 0.005,  0.005],
        [0.08, 0.015,  0.005],
        [0.02, 0.015,  0.005],
    ])

    mu_r = 1000.0
    chi = mu_r - 1.0

    return {
        'seg_p1': seg_p1, 'seg_p2': seg_p2,
        'seg_centers': seg_centers, 'seg_directions': seg_dirs,
        'seg_lengths': seg_lengths,
        'R_dc': R_dc, 'L_air': L_air,
        'core_verts': core_verts,
        'chi': chi, 'mu_r': mu_r,
    }


def get_hex_face_geometry(verts):
    """
    Compute face evaluation points, normals, and areas for a hexahedron.

    Face ordering MUST match MMMBuilder (rad_mmm_matrices.cpp):
        0: bottom (-Z)  v0, v3, v2, v1
        1: top    (+Z)  v4, v5, v6, v7
        2: front  (-Y)  v0, v1, v5, v4
        3: back   (+Y)  v2, v3, v7, v6
        4: left   (-X)  v0, v4, v7, v3
        5: right  (+X)  v1, v2, v6, v5

    Returns:
        eval_points: (6, 3) face midpoints
        normals: (6, 3) outward normals
        areas: (6,) face areas
    """
    v = verts

    # Face vertex indices matching MMMBuilder
    face_verts = [
        [v[0], v[3], v[2], v[1]],  # bottom (-Z)
        [v[4], v[5], v[6], v[7]],  # top    (+Z)
        [v[0], v[1], v[5], v[4]],  # front  (-Y)
        [v[2], v[3], v[7], v[6]],  # back   (+Y)
        [v[0], v[4], v[7], v[3]],  # left   (-X)
        [v[1], v[2], v[6], v[5]],  # right  (+X)
    ]

    eval_points = np.zeros((6, 3))
    normals = np.zeros((6, 3))
    areas = np.zeros(6)

    for f, fv in enumerate(face_verts):
        fv = np.array(fv)
        # Midpoint
        eval_points[f] = fv.mean(axis=0)
        # Diagonals for area and normal
        d1 = fv[2] - fv[0]
        d2 = fv[3] - fv[1]
        cross = np.cross(d1, d2)
        area = 0.5 * np.linalg.norm(cross)
        areas[f] = area
        if area > 0:
            normals[f] = cross / np.linalg.norm(cross)

    # Ensure normals point outward (check against center)
    center = verts.mean(axis=0)
    for f in range(6):
        if np.dot(normals[f], eval_points[f] - center) < 0:
            normals[f] = -normals[f]

    return eval_points, normals, areas


def demo_schur_complement():
    """Demo: Schur complement vs column-by-column Delta_L computation."""
    print("=" * 70)
    print("Schur Complement Solver for PEEC-MSC Coupling")
    print("=" * 70)
    print()

    prob = create_test_problem()

    # --- Build MSC interaction matrix using MMMBuilder ---
    builder = MMMBuilder()
    builder.add_hexahedron(prob['core_verts'])

    N_matrix, dof_offset = builder.build()
    N_matrix = np.array(N_matrix)
    dof_offset = np.array(dof_offset)

    n_msc = N_matrix.shape[0]
    chi = prob['chi']
    inv_chi = np.full(n_msc, 1.0 / chi)

    # --- Setup field computer for proper A-field ---
    field_comp = MMMFieldComputer()
    field_comp.set_elements_from_builder(builder)

    print(f"MSC DOF: {n_msc} (1 hexahedron, 6 faces)")
    print(f"PEEC DOF: {len(prob['seg_lengths'])} (1 segment)")
    print(f"chi = {chi:.0f} (mu_r = {prob['mu_r']:.0f})")
    print()

    # --- Get hex face geometry ---
    eval_points, face_normals, face_areas = get_hex_face_geometry(prob['core_verts'])

    # --- Setup Schur complement solver ---
    solver = SchurComplementSolver()
    solver.set_msc_system(N_matrix, dof_offset, inv_chi)
    solver.set_field_computer(field_comp)
    solver.set_peec_system(
        prob['L_air'], prob['R_dc'],
        prob['seg_p1'], prob['seg_p2'],
        prob['seg_centers'], prob['seg_directions'], prob['seg_lengths']
    )
    solver.set_msc_geometry(eval_points, face_normals, face_areas)

    # Compute coupling matrices using Radia (exact A-field from fixed M)
    solver.set_msc_hex_vertices(prob['core_verts'])
    B_pm, M_mp = solver.compute_coupling_matrices_radia([])
    print("Coupling matrices (Radia exact):")
    print(f"  B_pm shape: {B_pm.shape}  (MSC faces x PEEC segments)")
    print(f"  M_mp shape: {M_mp.shape}  (PEEC segments x MSC faces)")
    print(f"  B_pm norm: {np.linalg.norm(B_pm):.6e}")
    print(f"  M_mp norm: {np.linalg.norm(M_mp):.6e}")
    print()

    # Save Radia-computed matrices before dipole comparison
    B_pm_radia = solver._B_pm.copy()
    M_mp_radia = solver._M_mp.copy()

    # Show dipole-based for comparison (this overwrites internal state)
    B_pm_d, M_mp_d = solver.compute_coupling_matrices()
    print("Coupling matrices (dipole approx):")
    print(f"  M_mp norm: {np.linalg.norm(M_mp_d):.6e}")
    print()

    # Restore Radia-computed matrices for actual solve
    solver._B_pm = B_pm_radia
    solver._M_mp = M_mp_radia

    # --- Compute Delta_L via Schur complement ---
    Delta_L = solver.compute_delta_L(tol=1e-10)
    L_total = prob['L_air'] + Delta_L
    L_air_val = prob['L_air'][0, 0]

    print("Inductance results:")
    print(f"  L_air     = {L_air_val * 1e9:.4f} nH")
    print(f"  Delta_L   = {Delta_L[0, 0] * 1e9:.4f} nH")
    print(f"  L_total   = {L_total[0, 0] * 1e9:.4f} nH")
    print(f"  Ratio     = {L_total[0, 0] / L_air_val:.4f}")
    print()

    # --- Frequency-domain solve with Schur complement ---
    print("Frequency sweep (Schur complement):")
    freqs = np.array([1e3, 1e4, 1e5, 1e6, 1e7])
    V_source = np.array([1.0 + 0j])  # 1V source

    I_all, sigma_all = solver.frequency_sweep(freqs, V_source, tol=1e-10)

    print(f"  {'Freq [Hz]':>12s}  {'|I| [A]':>12s}  {'|Z| [Ohm]':>12s}  "
          f"{'R [Ohm]':>12s}  {'L [nH]':>10s}")
    for idx, f in enumerate(freqs):
        I = I_all[idx, 0]
        Z = V_source[0] / I
        R_eff = np.real(Z)
        L_eff = np.imag(Z) / (2.0 * np.pi * f)
        print(f"  {f:12.1e}  {abs(I):12.6e}  {abs(Z):12.6e}  "
              f"{R_eff:12.6e}  {L_eff*1e9:10.4f}")

    print()
    print("Done.")


def demo_static_coupling():
    """Demo: Static MSC solve with known DC current."""
    print("=" * 70)
    print("Static PEEC-MSC Coupling (DC)")
    print("=" * 70)
    print()

    prob = create_test_problem()

    builder = MMMBuilder()
    builder.add_hexahedron(prob['core_verts'])
    N_matrix, dof_offset = builder.build()
    N_matrix = np.array(N_matrix)
    dof_offset = np.array(dof_offset)

    n_msc = N_matrix.shape[0]
    inv_chi = np.full(n_msc, 1.0 / prob['chi'])

    eval_points, face_normals, face_areas = get_hex_face_geometry(prob['core_verts'])

    field_comp = MMMFieldComputer()
    field_comp.set_elements_from_builder(builder)

    solver = SchurComplementSolver()
    solver.set_msc_system(N_matrix, dof_offset, inv_chi)
    solver.set_field_computer(field_comp)
    solver.set_peec_system(
        prob['L_air'], prob['R_dc'],
        prob['seg_p1'], prob['seg_p2'],
        prob['seg_centers'], prob['seg_directions'], prob['seg_lengths']
    )
    solver.set_msc_geometry(eval_points, face_normals, face_areas)
    solver.compute_coupling_matrices()

    # Solve with 1A DC current
    I_dc = np.array([1.0])
    sigma, info = solver.solve_static(I_dc, tol=1e-10)

    print(f"DC current: {I_dc[0]:.1f} A")
    print(f"Converged: {info['converged']}")
    print(f"MSC surface charges (sigma):")
    for i in range(len(sigma)):
        print(f"  Face {i}: sigma = {sigma[i]:.6e}")
    print()


if __name__ == "__main__":
    demo_schur_complement()
    print()
    demo_static_coupling()
