"""Quick validation: updated ngsbem_eddy.py with separate-space approach."""
import sys
import os
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src'))

MU_0 = 4.0 * np.pi * 1e-7


def main():
    from radia.ngsbem_eddy import EddyCurrentFEMBEM, create_conductor_mesh
    from ngsolve import Integrate, CF
    from ngsolve import TaskManager

    mu_r = 1.0
    sigma = 3.5e7
    depth = 0.002
    width = 0.01
    height = 0.01
    freq = 5000.0
    omega = 2.0 * np.pi * freq
    Hz_inc = 1.0 / MU_0

    mesh = create_conductor_mesh(width, height, depth, maxh=0.004,
                                  conductor_label="conductor",
                                  surface_label="surface")
    with TaskManager():
        volume = Integrate(CF(1), mesh).real
        print(f"Volume: {volume:.6e} m^3")

        # === FEM reference ===
        print("\n--- FEM (Dirichlet BC) ---")
        eddy_fem = EddyCurrentFEMBEM(mesh, sigma=sigma, mu_r=mu_r, order=2,
                                       conductor_label="conductor",
                                       surface_label="surface")
        eddy_fem.assemble_fem(freq)
        eddy_fem.solve(B_ext=[0, 0, 1.0], mode='fem')
        Hz_avg_fem = Integrate(eddy_fem.get_total_field(), mesh) / volume
        mu_eff_fem = mu_r * Hz_avg_fem / Hz_inc
        eddy_fem.print_summary()
        print(f"  mu_eff = {mu_eff_fem.real:.6f}{mu_eff_fem.imag:+.6f}j")

        # === FEM-BEM ===
        print("\n--- FEM-BEM (Costabel, separate spaces) ---")
        eddy_bem = EddyCurrentFEMBEM(mesh, sigma=sigma, mu_r=mu_r, order=2,
                                       conductor_label="conductor",
                                       surface_label="surface")
        eddy_bem.assemble_fembem(freq)
        eddy_bem.solve(B_ext=[0, 0, 1.0], mode='fembem', printrates=False)
        Hz_avg_bem = Integrate(eddy_bem.get_total_field(), mesh) / volume
        mu_eff_bem = mu_r * Hz_avg_bem / Hz_inc
        eddy_bem.print_summary()
        print(f"  mu_eff = {mu_eff_bem.real:.6f}{mu_eff_bem.imag:+.6f}j")

        # === Analytical ===
        mu = mu_r * MU_0
        gamma = np.sqrt(1j * omega * mu * sigma)
        gd2 = gamma * depth / 2.0
        mu_eff_ana = mu_r * np.tanh(gd2) / gd2

        # === Summary ===
        print(f"\n=== mu_eff comparison ===")
        print(f"  FEM:        {mu_eff_fem.real:.6f}{mu_eff_fem.imag:+.6f}j")
        print(f"  FEM-BEM:    {mu_eff_bem.real:.6f}{mu_eff_bem.imag:+.6f}j")
        print(f"  Analytical: {mu_eff_ana.real:.6f}{mu_eff_ana.imag:+.6f}j")
        print(f"\n  |BEM - FEM|: {abs(mu_eff_bem - mu_eff_fem):.6f}")
        print(f"  |BEM - ana|: {abs(mu_eff_bem - mu_eff_ana):.6f}")
        print(f"  |FEM - ana|: {abs(mu_eff_fem - mu_eff_ana):.6f}")


if __name__ == '__main__':
    main()
