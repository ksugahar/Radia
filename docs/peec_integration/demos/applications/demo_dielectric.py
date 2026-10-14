"""
Demonstration of Dielectric Solver with Loop-Star Integration

This script demonstrates the surface charge method for dielectric analysis
and its integration with the Loop-Star PEEC framework.

Test cases:
1. Parallel plate capacitor (analytical validation)
2. LC circuit with dielectric capacitor
3. Coil near dielectric body

Author: Radia Development Team
Date: 2026-01-15
"""

import numpy as np
import matplotlib.pyplot as plt

from radia.dielectric_solver import (
    DielectricSolver, DielectricMesh,
    create_box_mesh, create_sphere_mesh,
    CombinedLoopStarDielectric, EPS_0
)


def test_parallel_plate_capacitor():
    """
    Test: Parallel plate capacitor

    Analytical capacitance:
        C = ε_0 * ε_r * A / d

    where A = plate area, d = plate separation
    """
    print("\n" + "=" * 70)
    print("Test 1: Parallel Plate Capacitor")
    print("=" * 70)

    # Geometry
    plate_size = 0.1  # 10 cm x 10 cm plates
    separation = 0.01  # 1 cm gap
    eps_r = 4.0  # Dielectric between plates

    A = plate_size ** 2
    d = separation

    # Analytical capacitance
    C_analytical = EPS_0 * eps_r * A / d
    print(f"Plate size: {plate_size*100:.1f} cm x {plate_size*100:.1f} cm")
    print(f"Separation: {separation*100:.1f} cm")
    print(f"Dielectric eps_r: {eps_r}")
    print(f"Analytical C = {C_analytical*1e12:.3f} pF")

    # Create mesh for two plates
    n_div = 5

    # Top plate mesh
    vertices_top = []
    for i in range(n_div + 1):
        for j in range(n_div + 1):
            x = -plate_size/2 + plate_size * i / n_div
            y = -plate_size/2 + plate_size * j / n_div
            z = separation / 2
            vertices_top.append([x, y, z])
    vertices_top = np.array(vertices_top)

    # Bottom plate mesh
    vertices_bottom = vertices_top.copy()
    vertices_bottom[:, 2] = -separation / 2

    # Triangles for each plate
    triangles = []
    for i in range(n_div):
        for j in range(n_div):
            i0 = i * (n_div + 1) + j
            i1 = i0 + 1
            i2 = i0 + (n_div + 1)
            i3 = i2 + 1
            triangles.append([i0, i2, i1])
            triangles.append([i1, i2, i3])
    triangles = np.array(triangles)

    # Create dielectric mesh (the dielectric slab between plates)
    # For this test, we treat the plate surfaces as conductors
    # and compute coupling between them

    # Create solver for top plate
    mesh_top = DielectricMesh(
        vertices=vertices_top,
        triangles=triangles,
        eps_r_inside=eps_r,
        eps_r_outside=1.0
    )

    # Create solver for bottom plate
    mesh_bottom = DielectricMesh(
        vertices=vertices_bottom,
        triangles=triangles,
        eps_r_inside=1.0,
        eps_r_outside=eps_r
    )

    # Build P-matrix for conductor surfaces
    solver = DielectricSolver()
    solver.add_dielectric(mesh_top, name='top_plate')
    solver.add_dielectric(mesh_bottom, name='bottom_plate')
    solver.build()

    P = solver.get_P_matrix()
    n = len(P)

    print(f"\nP-matrix size: {n} x {n}")
    print(f"P[0,0] = {P[0,0]:.3e}")
    print(f"P[0,{n//2}] = {P[0, n//2]:.3e}")

    # Compute capacitance from P-matrix
    # For two-terminal capacitor: C = 1 / (P_11 + P_22 - P_12 - P_21)
    # where 1 = all top plate DOFs, 2 = all bottom plate DOFs
    n_per_plate = n // 2

    # Sum P blocks
    P_11 = np.sum(P[:n_per_plate, :n_per_plate])
    P_22 = np.sum(P[n_per_plate:, n_per_plate:])
    P_12 = np.sum(P[:n_per_plate, n_per_plate:])
    P_21 = np.sum(P[n_per_plate:, :n_per_plate])

    P_total = P_11 + P_22 - P_12 - P_21
    C_numerical = n_per_plate ** 2 / P_total if abs(P_total) > 1e-30 else 0.0

    print(f"\nP_11 (top-top) = {P_11:.3e}")
    print(f"P_22 (bottom-bottom) = {P_22:.3e}")
    print(f"P_12 (top-bottom) = {P_12:.3e}")
    print(f"P_21 (bottom-top) = {P_21:.3e}")

    # Alternative: use capacitance matrix
    C_matrix = solver.get_capacitance_matrix()

    # Sum all elements for equivalent capacitance
    # C_total = sum(C_ij) for two-terminal device
    C_sum = np.sum(C_matrix)
    print(f"\nSum of C_matrix: {C_sum*1e12:.3f} pF")

    # For two-terminal: C = C_11 + C_22 + 2*C_12
    C_11 = np.sum(C_matrix[:n_per_plate, :n_per_plate])
    C_22 = np.sum(C_matrix[n_per_plate:, n_per_plate:])
    C_12 = np.sum(C_matrix[:n_per_plate, n_per_plate:])

    print(f"C_11 = {C_11*1e12:.3f} pF")
    print(f"C_22 = {C_22*1e12:.3f} pF")
    print(f"C_12 = {C_12*1e12:.3f} pF")

    # Two-terminal capacitance: apply +V to top, -V to bottom
    # Q_top = C_11*V + C_12*(-V) = (C_11 - C_12)*V
    # C_two_terminal = Q_top / (2V) = (C_11 - C_12) / 2
    C_two_terminal = (C_11 - C_12)

    print(f"\nNumerical C (two-terminal) = {C_two_terminal*1e12:.3f} pF")
    print(f"Analytical C = {C_analytical*1e12:.3f} pF")

    error = abs(C_two_terminal - C_analytical) / C_analytical * 100
    print(f"Error: {error:.1f}%")

    return C_analytical, C_two_terminal


def test_lc_circuit_with_dielectric():
    """
    Test: LC circuit with dielectric capacitor

    Parallel LC resonance with:
    - Inductor: L = 10 uH
    - Dielectric capacitor: C = ε_0 * ε_r * A / d

    SRF = 1 / (2*pi*sqrt(L*C))
    """
    print("\n" + "=" * 70)
    print("Test 2: LC Circuit with Dielectric Capacitor")
    print("=" * 70)

    # Inductor parameters
    L = 10e-6  # 10 uH
    R = 0.1    # 0.1 Ohm

    # Capacitor parameters (designed for ~16 MHz resonance)
    # C = 1 / ((2*pi*f)^2 * L)
    f_target = 16e6  # 16 MHz
    C_target = 1.0 / ((2.0 * np.pi * f_target) ** 2 * L)
    print(f"Target frequency: {f_target/1e6:.1f} MHz")
    print(f"Required C: {C_target*1e12:.2f} pF")

    # Design capacitor geometry
    eps_r = 10.0
    d = 0.001  # 1 mm separation
    A = C_target * d / (EPS_0 * eps_r)
    plate_size = np.sqrt(A)
    print(f"Capacitor design:")
    print(f"  eps_r = {eps_r}")
    print(f"  separation = {d*1000:.1f} mm")
    print(f"  plate size = {plate_size*1000:.1f} mm x {plate_size*1000:.1f} mm")

    # Analytical capacitance
    C_analytical = EPS_0 * eps_r * A / d
    f_analytical = 1.0 / (2.0 * np.pi * np.sqrt(L * C_analytical))
    print(f"Analytical C: {C_analytical*1e12:.3f} pF")
    print(f"Analytical SRF: {f_analytical/1e6:.2f} MHz")

    # Create combined Loop-Star-Dielectric solver
    solver = CombinedLoopStarDielectric()

    # Set Loop matrices (inductor)
    L_LL = np.array([[L]])
    R_LL = np.array([[R]])
    solver.set_loop_matrices(L_LL, R_LL)

    # Set dielectric Star (capacitor as P-matrix)
    # For ideal parallel plate: P = d / (ε_0 * ε_r * A)
    P_dd = np.array([[d / (EPS_0 * eps_r * A)]])
    solver.set_dielectric_star(P_dd)

    # Frequency sweep
    frequencies = np.logspace(6, 8, 300)  # 1 MHz to 100 MHz
    result = solver.solve_frequency(frequencies)

    # Find resonance
    idx_max = np.argmax(np.abs(result['impedance']))
    f_resonance = frequencies[idx_max]
    Z_max = np.abs(result['impedance'][idx_max])

    print(f"\nNumerical SRF: {f_resonance/1e6:.2f} MHz")
    print(f"|Z| at SRF: {Z_max:.0f} Ohm")

    error = abs(f_resonance - f_analytical) / f_analytical * 100
    print(f"Error: {error:.1f}%")

    # Plot
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))

    ax = axes[0, 0]
    ax.loglog(frequencies / 1e6, np.abs(result['impedance']))
    ax.axvline(x=f_analytical / 1e6, color='r', linestyle='--', label=f'Analytical SRF = {f_analytical/1e6:.1f} MHz')
    ax.axvline(x=f_resonance / 1e6, color='g', linestyle=':', label=f'Numerical SRF = {f_resonance/1e6:.1f} MHz')
    ax.set_xlabel('Frequency [MHz]')
    ax.set_ylabel('|Z| [Ohm]')
    ax.set_title('LC Circuit Impedance (Dielectric Capacitor)')
    ax.legend()
    ax.grid(True, which='both', alpha=0.3)

    ax = axes[0, 1]
    phase = np.angle(result['impedance']) * 180 / np.pi
    ax.semilogx(frequencies / 1e6, phase)
    ax.axhline(y=0, color='k', linestyle=':', alpha=0.3)
    ax.axvline(x=f_analytical / 1e6, color='r', linestyle='--')
    ax.set_xlabel('Frequency [MHz]')
    ax.set_ylabel('Phase [deg]')
    ax.set_title('Impedance Phase')
    ax.grid(True, which='both', alpha=0.3)

    ax = axes[1, 0]
    ax.semilogx(frequencies / 1e6, result['resistance'])
    ax.set_xlabel('Frequency [MHz]')
    ax.set_ylabel('R [Ohm]')
    ax.set_title('Resistance')
    ax.grid(True, which='both', alpha=0.3)

    ax = axes[1, 1]
    ax.semilogx(frequencies[1:] / 1e6, result['inductance'][1:] * 1e6)
    ax.set_xlabel('Frequency [MHz]')
    ax.set_ylabel('L [uH]')
    ax.set_title('Effective Inductance')
    ax.grid(True, which='both', alpha=0.3)

    plt.tight_layout()
    plt.savefig('demo_dielectric_lc.png', dpi=150)
    print("\nSaved: demo_dielectric_lc.png")

    return f_analytical, f_resonance


def test_sphere_mesh():
    """
    Test: Dielectric sphere mesh generation
    """
    print("\n" + "=" * 70)
    print("Test 3: Dielectric Sphere Mesh")
    print("=" * 70)

    # Create sphere mesh
    center = np.array([0.0, 0.0, 0.0])
    radius = 0.05  # 5 cm
    mesh = create_sphere_mesh(center, radius, n_theta=10, n_phi=20)

    print(f"Sphere radius: {radius*100:.1f} cm")
    print(f"Vertices: {len(mesh.vertices)}")
    print(f"Triangles: {len(mesh.triangles)}")
    print(f"Total area: {np.sum(mesh.areas):.6f} m^2")
    print(f"Analytical area (4*pi*r^2): {4*np.pi*radius**2:.6f} m^2")

    # Create solver
    mesh.eps_r_inside = 4.0
    mesh.eps_r_outside = 1.0
    solver = DielectricSolver()
    solver.add_dielectric(mesh, name='sphere')
    solver.build()

    P = solver.get_P_matrix()
    print(f"\nP-matrix size: {P.shape}")
    print(f"P diagonal mean: {np.mean(np.diag(P)):.3e}")
    print(f"P off-diagonal mean: {np.mean(P - np.diag(np.diag(P))):.3e}")

    # Capacitance of dielectric sphere in free space
    # C = 4*pi*eps_0*eps_r*a (isolated sphere)
    C_analytical = 4 * np.pi * EPS_0 * mesh.eps_r_inside * radius
    print(f"\nAnalytical capacitance (isolated sphere): {C_analytical*1e12:.3f} pF")

    return mesh


def test_inductor_with_magnetic_core():
    """
    Test: Inductor with magnetic core (UnifiedSurfaceSolver)

    Validates that magnetic material coupling enhances inductance.

    For an inductor with magnetic core:
        L_eff = L_air * mu_eff

    where mu_eff depends on geometry and demagnetization factor.
    """
    print("\n" + "=" * 70)
    print("Test 4: Inductor with Magnetic Core")
    print("=" * 70)

    from radia.dielectric_solver import UnifiedSurfaceSolver

    # Inductor parameters (air core)
    L_air = 10e-6  # 10 uH
    R = 0.1  # 0.1 Ohm

    # Magnetic core parameters
    mu_r = 100.0  # Relative permeability

    # Demagnetization factor (approximation for cylindrical core)
    # N ~ 0.2 for length/diameter ~ 2
    N = 0.2

    # Expected inductance enhancement
    # For a core inside coil: L_eff = L_air * mu_eff
    # where mu_eff = mu_r / (1 + (mu_r - 1) * N)
    mu_eff = mu_r / (1.0 + (mu_r - 1.0) * N)
    L_expected = L_air * mu_eff

    print(f"Air core inductance: {L_air*1e6:.2f} uH")
    print(f"Core mu_r: {mu_r}")
    print(f"Demagnetization factor N: {N}")
    print(f"Effective mu: {mu_eff:.2f}")
    print(f"Expected L_eff: {L_expected*1e6:.2f} uH")

    # Create UnifiedSurfaceSolver
    solver = UnifiedSurfaceSolver()

    # Set Loop (conductor) matrices
    L_LL = np.array([[L_air]])
    R_LL = np.array([[R]])
    solver.set_conductor_loop(L_LL, R_LL)

    # Set magnetic material
    # N-matrix is the demagnetization tensor (scalar for isotropic)
    N_mm = np.array([[N]])

    # Coupling: current creates H-field at core
    # For unit current, H_at_core = n_turns / length ~ some coupling coefficient
    # We use L_Lm as the coupling factor
    # The induced flux from magnetization: psi = L_Lm^T * M
    # Back-EMF: V = -s * psi = -s * L_Lm^T * M
    #
    # For correct inductance enhancement, the coupling must satisfy:
    # L_eff = L_air + (L_Lm^T * chi * inv(I + chi*N) * L_Lm)
    #
    # With chi = mu_r - 1, for single DOF:
    # L_eff = L_air + L_Lm^2 * chi / (1 + chi*N)
    #
    # To get L_expected = L_air * mu_eff:
    # L_Lm^2 * chi / (1 + chi*N) = L_air * (mu_eff - 1)
    # L_Lm = sqrt(L_air * (mu_eff - 1) * (1 + chi*N) / chi)

    chi = mu_r - 1.0
    coupling_factor = np.sqrt(L_air * (mu_eff - 1.0) * (1.0 + chi * N) / chi)
    L_Lm = np.array([[coupling_factor]])

    print(f"\nCoupling coefficient L_Lm: {coupling_factor*1e6:.4f} uH/A")

    solver.set_magnetic_material(N_mm, mu_r, L_Lm)

    # Frequency sweep
    frequencies = np.logspace(3, 7, 100)  # 1 kHz to 10 MHz
    result = solver.solve_frequency(frequencies)

    # Check effective inductance at low frequency
    idx_low = 10  # ~1.6 kHz
    L_numerical = result['inductance'][idx_low]
    f_test = frequencies[idx_low]

    print(f"\nAt f = {f_test/1e3:.1f} kHz:")
    print(f"  Numerical L_eff: {L_numerical*1e6:.2f} uH")
    print(f"  Expected L_eff: {L_expected*1e6:.2f} uH")

    error = abs(L_numerical - L_expected) / L_expected * 100
    print(f"  Error: {error:.1f}%")

    # Plot
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))

    ax = axes[0, 0]
    ax.loglog(frequencies / 1e3, np.abs(result['impedance']))
    ax.set_xlabel('Frequency [kHz]')
    ax.set_ylabel('|Z| [Ohm]')
    ax.set_title('Inductor with Magnetic Core: Impedance')
    ax.grid(True, which='both', alpha=0.3)

    ax = axes[0, 1]
    phase = np.angle(result['impedance']) * 180 / np.pi
    ax.semilogx(frequencies / 1e3, phase)
    ax.axhline(y=90, color='r', linestyle='--', alpha=0.5, label='Pure inductance')
    ax.set_xlabel('Frequency [kHz]')
    ax.set_ylabel('Phase [deg]')
    ax.set_title('Impedance Phase')
    ax.legend()
    ax.grid(True, which='both', alpha=0.3)

    ax = axes[1, 0]
    ax.semilogx(frequencies / 1e3, result['inductance'] * 1e6)
    ax.axhline(y=L_expected * 1e6, color='r', linestyle='--',
               label=f'Expected: {L_expected*1e6:.1f} uH')
    ax.axhline(y=L_air * 1e6, color='g', linestyle=':',
               label=f'Air core: {L_air*1e6:.1f} uH')
    ax.set_xlabel('Frequency [kHz]')
    ax.set_ylabel('L [uH]')
    ax.set_title('Effective Inductance')
    ax.legend()
    ax.grid(True, which='both', alpha=0.3)

    ax = axes[1, 1]
    ax.semilogx(frequencies / 1e3, result['resistance'])
    ax.set_xlabel('Frequency [kHz]')
    ax.set_ylabel('R [Ohm]')
    ax.set_title('Resistance')
    ax.grid(True, which='both', alpha=0.3)

    plt.tight_layout()
    plt.savefig('demo_dielectric_magnetic_core.png', dpi=150)
    print("\nSaved: demo_dielectric_magnetic_core.png")

    return L_expected, L_numerical


def test_lc_with_magnetic_core():
    """
    Test: LC circuit with dielectric capacitor AND magnetic core

    This tests the full UnifiedSurfaceSolver with:
    - Loop (conductor inductance)
    - Star (dielectric capacitor)
    - Magnetic (core material)

    Expected: Resonance at lower frequency due to enhanced inductance.
    """
    print("\n" + "=" * 70)
    print("Test 5: LC Circuit with Magnetic Core")
    print("=" * 70)

    from radia.dielectric_solver import UnifiedSurfaceSolver

    # Inductor parameters
    L_air = 10e-6  # 10 uH
    R = 0.1  # 0.1 Ohm

    # Magnetic core parameters
    mu_r = 100.0
    N = 0.2  # Demagnetization factor
    mu_eff = mu_r / (1.0 + (mu_r - 1.0) * N)
    L_eff = L_air * mu_eff

    # Capacitor (dielectric)
    C = 10e-12  # 10 pF

    # Expected resonance frequencies
    f_air = 1.0 / (2.0 * np.pi * np.sqrt(L_air * C))
    f_core = 1.0 / (2.0 * np.pi * np.sqrt(L_eff * C))

    print(f"Air core L: {L_air*1e6:.2f} uH")
    print(f"Enhanced L: {L_eff*1e6:.2f} uH")
    print(f"Capacitance: {C*1e12:.2f} pF")
    print(f"Resonance (air core): {f_air/1e6:.2f} MHz")
    print(f"Resonance (magnetic core): {f_core/1e6:.2f} MHz")

    # Create solver
    solver = UnifiedSurfaceSolver()

    # Loop matrices
    solver.set_conductor_loop(np.array([[L_air]]), np.array([[R]]))

    # Dielectric Star (capacitor as P = 1/C)
    P_dd = np.array([[1.0 / C]])
    solver.set_dielectric_star(P_dd)

    # Magnetic material
    N_mm = np.array([[N]])
    chi = mu_r - 1.0
    coupling = np.sqrt(L_air * (mu_eff - 1.0) * (1.0 + chi * N) / chi)
    L_Lm = np.array([[coupling]])
    solver.set_magnetic_material(N_mm, mu_r, L_Lm)

    # Frequency sweep
    frequencies = np.logspace(5, 8, 300)  # 100 kHz to 100 MHz
    result = solver.solve_frequency(frequencies)

    # Find resonance (peak impedance)
    idx_max = np.argmax(np.abs(result['impedance']))
    f_resonance = frequencies[idx_max]
    Z_max = np.abs(result['impedance'][idx_max])

    print(f"\nNumerical resonance: {f_resonance/1e6:.2f} MHz")
    print(f"|Z| at resonance: {Z_max:.0f} Ohm")

    error = abs(f_resonance - f_core) / f_core * 100
    print(f"Error vs expected: {error:.1f}%")

    # Also test effective inductance at low frequency
    idx_low = 20
    L_numerical = result['inductance'][idx_low]
    print(f"\nEffective L at {frequencies[idx_low]/1e3:.1f} kHz: {L_numerical*1e6:.2f} uH")

    # Plot
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))

    ax = axes[0, 0]
    ax.loglog(frequencies / 1e6, np.abs(result['impedance']))
    ax.axvline(x=f_air / 1e6, color='g', linestyle=':', label=f'Air SRF = {f_air/1e6:.1f} MHz')
    ax.axvline(x=f_core / 1e6, color='r', linestyle='--', label=f'Core SRF = {f_core/1e6:.1f} MHz')
    ax.axvline(x=f_resonance / 1e6, color='b', linestyle='-', alpha=0.5,
               label=f'Numerical = {f_resonance/1e6:.1f} MHz')
    ax.set_xlabel('Frequency [MHz]')
    ax.set_ylabel('|Z| [Ohm]')
    ax.set_title('LC with Magnetic Core: Impedance')
    ax.legend()
    ax.grid(True, which='both', alpha=0.3)

    ax = axes[0, 1]
    phase = np.angle(result['impedance']) * 180 / np.pi
    ax.semilogx(frequencies / 1e6, phase)
    ax.axhline(y=0, color='k', linestyle=':', alpha=0.3)
    ax.axvline(x=f_core / 1e6, color='r', linestyle='--')
    ax.set_xlabel('Frequency [MHz]')
    ax.set_ylabel('Phase [deg]')
    ax.set_title('Impedance Phase')
    ax.grid(True, which='both', alpha=0.3)

    ax = axes[1, 0]
    ax.semilogx(frequencies / 1e6, result['inductance'] * 1e6)
    ax.axhline(y=L_eff * 1e6, color='r', linestyle='--',
               label=f'Expected: {L_eff*1e6:.1f} uH')
    ax.axhline(y=L_air * 1e6, color='g', linestyle=':',
               label=f'Air core: {L_air*1e6:.1f} uH')
    ax.set_xlabel('Frequency [MHz]')
    ax.set_ylabel('L [uH]')
    ax.set_title('Effective Inductance')
    ax.legend()
    ax.grid(True, which='both', alpha=0.3)

    ax = axes[1, 1]
    ax.semilogx(frequencies / 1e6, result['resistance'])
    ax.set_xlabel('Frequency [MHz]')
    ax.set_ylabel('R [Ohm]')
    ax.set_title('Resistance')
    ax.grid(True, which='both', alpha=0.3)

    plt.tight_layout()
    plt.savefig('demo_dielectric_lc_magnetic.png', dpi=150)
    print("\nSaved: demo_dielectric_lc_magnetic.png")

    return f_core, f_resonance


def main():
    """Main demonstration."""
    print("=" * 70)
    print("Dielectric Solver Demonstration")
    print("=" * 70)

    # Test 1: Parallel plate capacitor
    test_parallel_plate_capacitor()

    # Test 2: LC circuit with dielectric
    test_lc_circuit_with_dielectric()

    # Test 3: Sphere mesh
    test_sphere_mesh()

    # Test 4: Inductor with magnetic core
    test_inductor_with_magnetic_core()

    # Test 5: LC with magnetic core
    test_lc_with_magnetic_core()

    print("\n" + "=" * 70)
    print("Demo completed!")
    print("=" * 70)

    plt.show()


if __name__ == "__main__":
    main()
