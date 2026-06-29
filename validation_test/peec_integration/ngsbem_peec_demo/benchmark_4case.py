"""
4-case material coupling benchmark: ngbem surface BEM vs PEEC filament.

Cases:
  1. Air only         - ngbem HelmholtzSL / PEEC Neumann
  2. +Ferrite core    - Image method: L_total = L_air * 2*mu_r/(mu_r+1)
  3. +Al Shield       - ShieldBEMSIBC eddy current (PEEC filaments)
  4. +Both            - Combined

Frame geometry: 10mm x 10mm, 1mm trace width, 35um Cu thickness.
"""
import sys
import os
import numpy as np
from scipy.linalg import null_space

# Ensure MKL DLLs are found
os.add_dll_directory(os.path.join(sys.prefix, 'Library', 'bin'))

# Add module path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/radia'))

MU_0 = 4.0 * np.pi * 1e-7
WIDTH = 0.01       # 10 mm
TRACE_W = 1e-3     # 1 mm
hw = TRACE_W / 2.0
SIGMA_CU = 5.8e7   # Cu conductivity
THICKNESS = 35e-6   # Cu thickness
SIGMA_AL = 3.7e7   # Al shield conductivity
MU_R_FERRITE = 1000

# Shield geometry
SHIELD_Z = 0.005     # 5mm above frame
SHIELD_T = 0.5e-3    # 0.5mm thick
SHIELD_SIZE = 0.014   # 14mm x 14mm


# ---- Helper: extract ngsolve matrix ----
def extract_dense(mat, n):
    ei = mat.CreateColVector()
    col = mat.CreateColVector()
    M = np.zeros((n, n), dtype=complex)
    for i in range(n):
        ei[:] = 0; ei[i] = 1.0
        mat.Mult(ei, col)
        for j in range(n):
            M[j, i] = col[j]
    return M


def extract_rect(mat, nrow, ncol):
    ei = mat.CreateRowVector()
    col = mat.CreateColVector()
    M = np.zeros((nrow, ncol))
    for i in range(ncol):
        ei[:] = 0; ei[i] = 1.0
        mat.Mult(ei, col)
        for j in range(nrow):
            M[j, i] = col[j]
    return M


# ---- ngbem surface BEM approach ----
def compute_ngbem_L_air():
    """Compute L_air using ngbem RT0 harmonic + HelmholtzSL."""
    from netgen.occ import WorkPlane, OCCGeometry, Axes, Pnt, Dir
    from ngsolve import (Mesh, HDivSurface, SurfaceL2, BilinearForm,
                         ds, BND, TaskManager, InnerProduct)
    from ngsolve import div as ng_div
    from ngsolve.bem import HelmholtzSL

    # Build frame mesh
    wp_o = WorkPlane(Axes(Pnt(-hw, -hw, 0), Dir(0, 0, 1), Dir(1, 0, 0)))
    outer = wp_o.Rectangle(WIDTH + TRACE_W, WIDTH + TRACE_W).Face()
    wp_i = WorkPlane(Axes(Pnt(hw, hw, 0), Dir(0, 0, 1), Dir(1, 0, 0)))
    inner = wp_i.Rectangle(WIDTH - TRACE_W, WIDTH - TRACE_W).Face()
    frame = outer - inner
    frame.faces.name = "conductor"
    geo = OCCGeometry(frame)
    mesh = Mesh(geo.GenerateMesh(maxh=0.001))
    n_el = sum(1 for _ in mesh.Elements(BND))

    fes_J = HDivSurface(mesh, order=0)
    n_J = fes_J.ndof
    n_v = mesh.nv
    fes_L2 = SurfaceL2(mesh, order=0)
    n_f = fes_L2.ndof

    # D matrix (divergence)
    u_J = fes_J.TrialFunction()
    q_L2 = fes_L2.TestFunction()
    bf_D = BilinearForm(trialspace=fes_J, testspace=fes_L2)
    bf_D += ng_div(u_J.Trace()) * q_L2 * ds
    bf_D.Assemble()
    D = extract_rect(bf_D.mat, n_f, n_J)

    # C matrix (surface curl = signed incidence)
    C = np.zeros((n_J, n_v))
    for e_idx, edge in enumerate(mesh.edges):
        verts = list(edge.vertices)
        C[e_idx, verts[0].nr] = -1
        C[e_idx, verts[1].nr] = +1

    # M_J (mass matrix)
    u2, v2 = fes_J.TnT()
    bf_M = BilinearForm(fes_J)
    bf_M += InnerProduct(u2.Trace(), v2.Trace()) * ds
    bf_M.Assemble()
    M_J = np.real(extract_dense(bf_M.mat, n_J))

    # Harmonic (M_J-orthogonal Hodge decomposition)
    constraint = np.vstack([D, C.T @ M_J])
    null_h = null_space(constraint, rcond=1e-10)
    c_h = null_h[:, 0]
    energy = c_h @ M_J @ c_h

    # V_A via HelmholtzSL
    kappa = 0.01
    u3, v3 = fes_J.TnT()
    with TaskManager():
        V1_op = HelmholtzSL(
            u3.Trace() * ds("conductor", bonus_intorder=6),
            kappa
        ) * v3.Trace() * ds("conductor", bonus_intorder=6)

    V1_mat = extract_dense(V1_op.mat, n_J)
    V_A = np.real(c_h @ V1_mat @ c_h)

    # L from V_A
    R_sheet = 1.0 / (SIGMA_CU * THICKNESS)
    perimeter = 4 * WIDTH
    R_loop = R_sheet * perimeter / TRACE_W
    LR_ratio = MU_0 * V_A / (R_sheet * energy)
    L = LR_ratio * R_loop

    return {
        'L_air': L,
        'R_loop': R_loop,
        'n_el': n_el,
        'n_J': n_J,
        'V_A': V_A,
        'energy': energy,
    }


# ---- PEEC filament approach ----
def compute_peec_4cases(freq_test=10000):
    """Compute all 4 cases using PEEC filaments + ShieldBEMSIBC."""
    from peec_matrices import PEECBuilder
    from peec_topology import PEECCircuitSolver

    half = WIDTH / 2.0
    w = TRACE_W
    h = THICKNESS
    sigma = SIGMA_CU

    builder = PEECBuilder()
    n1 = builder.add_node_at(-half, -half, 0)
    n2 = builder.add_node_at(half, -half, 0)
    n3 = builder.add_node_at(half, half, 0)
    n4 = builder.add_node_at(-half, half, 0)
    n5 = builder.add_node_at(-half, -half, 0)  # port- (same location)

    builder.add_connected_segment(n1, n2, w, h, sigma)
    builder.add_connected_segment(n2, n3, w, h, sigma)
    builder.add_connected_segment(n3, n4, w, h, sigma)
    builder.add_connected_segment(n4, n5, w, h, sigma)
    builder.add_port(n1, n5)

    topo = builder.build_topology()
    solver = PEECCircuitSolver(topo)

    n_seg = topo['n_loop']
    L_mat = topo['L']
    R_vec = topo['R']

    # Case 1: Air only
    Z_air = solver.compute_port_impedance(freq=1.0)
    L_air = np.imag(Z_air) / (2 * np.pi * 1.0)
    R_dc = np.real(Z_air)

    # Case 2: Ferrite (image method)
    factor_ferrite = 2.0 * MU_R_FERRITE / (MU_R_FERRITE + 1.0)
    L_ferrite = L_air * factor_ferrite

    # Case 3: Shield
    omega_test = 2.0 * np.pi * freq_test
    delta_shield = np.sqrt(2.0 / (omega_test * MU_0 * SIGMA_AL))

    # Build shield mesh
    from netgen.occ import Box, Pnt, OCCGeometry
    from ngsolve import Mesh as NGMesh

    shield_box = Box(
        Pnt(-SHIELD_SIZE / 2, -SHIELD_SIZE / 2, SHIELD_Z),
        Pnt(SHIELD_SIZE / 2, SHIELD_SIZE / 2, SHIELD_Z + SHIELD_T))
    shield_box.solids.name = "conductor"
    shield_box.faces.name = "surface"
    shield_geo = OCCGeometry(shield_box)
    shield_mesh = NGMesh(shield_geo.GenerateMesh(maxh=0.003))

    from ngbem_eddy import ShieldBEMSIBC

    shield = ShieldBEMSIBC(shield_mesh, sigma=SIGMA_AL)
    print("  Assembling shield BEM...")
    shield.assemble(intorder=4)
    print(f"  Assembly: {shield.t_assemble:.1f} s")
    print(f"  Skin depth at {freq_test / 1e3:.0f} kHz: {delta_shield * 1e3:.2f} mm")

    # Segment geometry from PEEC topology
    seg_centers = np.array(topo['segment_centers'])
    seg_dirs = np.array(topo['segment_directions'])
    seg_lens = np.array(topo['segment_lengths'])

    peec_topo = {
        'segment_centers': seg_centers,
        'segment_directions': seg_dirs,
        'segment_lengths': seg_lens,
    }

    print(f"  Computing Delta_Z ({n_seg}x{n_seg}) at {freq_test / 1e3:.0f} kHz...")
    Delta_Z = shield.compute_impedance_matrix(freq_test, peec_topo)

    # For single-loop (series segments): Z_loop = ones^T @ Z_branch @ ones
    ones = np.ones(n_seg)

    # Air branch impedance
    Z_branch_air = np.diag(R_vec.astype(complex)) + 1j * omega_test * L_mat
    Z_loop_air = ones @ Z_branch_air @ ones

    # Shield branch impedance
    Z_branch_shield = Z_branch_air + Delta_Z
    Z_loop_shield = ones @ Z_branch_shield @ ones

    Delta_Z_loop = ones @ Delta_Z @ ones
    Delta_L_shield = np.imag(Delta_Z_loop) / omega_test
    Delta_R_shield = np.real(Delta_Z_loop)

    L_shield = L_air + Delta_L_shield
    R_shield = R_dc + Delta_R_shield

    # Case 4: Both
    L_both = L_ferrite + Delta_L_shield
    R_both = R_dc + Delta_R_shield

    return {
        'L_air': L_air,
        'R_dc': R_dc,
        'L_ferrite': L_ferrite,
        'factor_ferrite': factor_ferrite,
        'L_shield': L_shield,
        'R_shield': R_shield,
        'Delta_L_shield': Delta_L_shield,
        'Delta_R_shield': Delta_R_shield,
        'L_both': L_both,
        'R_both': R_both,
        'delta_shield': delta_shield,
        'freq_test': freq_test,
        'n_seg': n_seg,
        'L_mat': L_mat,
        'Delta_Z': Delta_Z,
    }


def main():
    print("=" * 70)
    print("4-Case Material Coupling Benchmark")
    print("Frame: 10mm x 10mm, 1mm trace, 35um Cu")
    print("=" * 70)

    # ---- ngbem surface BEM (air only) ----
    print("\n--- ngbem Surface BEM (Air only) ---")
    ngbem_result = compute_ngbem_L_air()
    L_ngbem = ngbem_result['L_air']
    R_ngbem = ngbem_result['R_loop']
    print(f"  {ngbem_result['n_el']} elements, {ngbem_result['n_J']} edge DOFs")
    print(f"  L_air = {L_ngbem * 1e9:.2f} nH")
    print(f"  R_loop = {R_ngbem * 1e3:.2f} mOhm")

    # ---- PEEC filament model (all 4 cases) ----
    print("\n--- PEEC Filament Model ---")
    peec = compute_peec_4cases(freq_test=10000)
    L_peec = peec['L_air']
    R_peec = peec['R_dc']
    print(f"  {peec['n_seg']} filaments (4 sides)")

    # Summary
    print(f"\n{'=' * 70}")
    print(f"{'':>30s}  {'PEEC':>12s}  {'ngbem':>12s}  {'Ratio':>8s}")
    print(f"{'':>30s}  {'(filament)':>12s}  {'(surface)':>12s}")
    print(f"{'-' * 70}")
    print(f"{'L_air [nH]':>30s}  {L_peec * 1e9:12.2f}  "
          f"{L_ngbem * 1e9:12.2f}  {L_ngbem / L_peec:8.3f}")
    print(f"{'R_dc [mOhm]':>30s}  {R_peec * 1e3:12.2f}  "
          f"{R_ngbem * 1e3:12.2f}  {R_ngbem / R_peec:8.3f}")

    # PEEC 4-case table
    print(f"\n{'=' * 70}")
    print(f"PEEC 4-Case Results (f={peec['freq_test'] / 1e3:.0f} kHz, "
          f"delta={peec['delta_shield'] * 1e3:.2f} mm)")
    print(f"{'-' * 70}")
    print(f"{'Case':>25s}  {'L [nH]':>10s}  {'R [mOhm]':>10s}  "
          f"{'L/L_air':>8s}  {'Check':>6s}")
    print(f"{'-' * 70}")

    cases = [
        ('1. Air only', peec['L_air'], peec['R_dc'], 1.0, 'OK'),
        ('2. +Ferrite (mu_r=1000)',
         peec['L_ferrite'], peec['R_dc'],
         peec['L_ferrite'] / peec['L_air'],
         'OK' if peec['L_ferrite'] > peec['L_air'] else 'FAIL'),
        ('3. +Al Shield',
         peec['L_shield'], peec['R_shield'],
         peec['L_shield'] / peec['L_air'],
         'OK' if peec['Delta_L_shield'] < 0 and peec['Delta_R_shield'] > 0 else 'FAIL'),
        ('4. +Ferrite + Shield',
         peec['L_both'], peec['R_both'],
         peec['L_both'] / peec['L_air'], ''),
    ]

    for name, L, R, ratio, check in cases:
        print(f"{name:>25s}  {L * 1e9:10.2f}  {R * 1e3:10.2f}  "
              f"{ratio:8.4f}  {check:>6s}")

    print(f"{'-' * 70}")

    # Physics checks
    print(f"\n  Shield effects:")
    print(f"    Delta_L = {peec['Delta_L_shield'] * 1e9:.3f} nH "
          f"({peec['Delta_L_shield'] / peec['L_air'] * 100:.2f}%)")
    print(f"    Delta_R = {peec['Delta_R_shield'] * 1e3:.4f} mOhm "
          f"({peec['Delta_R_shield'] / peec['R_dc'] * 100:.2f}%)")
    print(f"    Lenz (Delta_L < 0): "
          f"{'OK' if peec['Delta_L_shield'] < 0 else 'FAIL'}")
    print(f"    Loss (Delta_R > 0): "
          f"{'OK' if peec['Delta_R_shield'] > 0 else 'FAIL'}")

    # Delta_Z matrix inspection
    print(f"\n  Delta_Z matrix ({peec['n_seg']}x{peec['n_seg']}):")
    DZ = peec['Delta_Z']
    print(f"    |Re| range: [{np.min(np.abs(np.real(DZ))):.3e}, "
          f"{np.max(np.abs(np.real(DZ))):.3e}]")
    print(f"    |Im| range: [{np.min(np.abs(np.imag(DZ))):.3e}, "
          f"{np.max(np.abs(np.imag(DZ))):.3e}]")

    # Relative effects (model-independent)
    print(f"\n  Ferrite factor: {peec['factor_ferrite']:.4f}")
    print(f"  (same for any model: 2*mu_r/(mu_r+1))")

    print("\nDone.")


if __name__ == '__main__':
    main()
