"""
TEAM Problem 7: A-method Solver

Solves the time-harmonic eddy current problem using the A-phi formulation:
  curl(1/mu * curl(A)) + j*omega*sigma*A = J0

For the conducting region:
  curl(1/mu * curl(A)) + j*omega*sigma*A = 0

For the non-conducting region (air, coil):
  curl(1/mu * curl(A)) = J0

Reference:
- https://www.comsol.com/blogs/solving-team-problem-7-acdc-module/
- https://modules.freefem.org/modules/team7/
"""
import os
import sys

os.environ['NETGENDIR'] = r'S:\NGSolve\01_GitHub\install_ngsolve\bin'
os.environ['PATH'] = r'S:\NGSolve\01_GitHub\install_ngsolve\bin;' + os.environ.get('PATH', '')
sys.path.insert(0, r'S:\NGSolve\01_GitHub\install_ngsolve\Lib\site-packages')

from numpy import pi, sqrt
import numpy as np
from ngsolve import *
from ngsolve.krylovspace import CGSolver

from team7_geometry import (
    create_team7_geometry,
    get_coil_current_density,
    get_material_properties,
    measurement_line_A1B1,
    measurement_line_A2B2,
    freq_50, freq_200,
    mu0, sigma_al
)


def solve_time_harmonic(mesh, freq, order=1, verbose=True):
    """
    Solve time-harmonic eddy current problem at given frequency.

    Args:
        mesh: NGSolve mesh
        freq: Frequency [Hz]
        order: Polynomial order for HCurl elements
        verbose: Print progress

    Returns:
        gfA: Complex GridFunction for vector potential A
        B: Magnetic flux density B = curl(A)
        J_eddy: Eddy current density J = -j*omega*sigma*A
    """
    omega = 2 * pi * freq

    if verbose:
        print(f"\nSolving time-harmonic problem at {freq} Hz...")
        print(f"  omega = {omega:.2f} rad/s")

    # Important: Curve order should match FES polynomial order for good convergence
    with TaskManager():
        mesh.Curve(order)

        # Get material properties
        mu, sigma, inv_mu = get_material_properties(mesh)

        # Get source current density
        J0 = get_coil_current_density(mesh)

        # Complex HCurl space
        fes = HCurl(mesh, order=order, dirichlet="outer", complex=True)
        A, v = fes.TnT()

        if verbose:
            print(f"  DOFs: {fes.ndof}")

        # Bilinear form:
        # curl(1/mu * curl(A)) + j*omega*sigma*A = J0
        # Weak form: <1/mu * curl(A), curl(v)> + <j*omega*sigma*A, v> = <J0, v>
        a = BilinearForm(fes)
        a += inv_mu * InnerProduct(curl(A), curl(v)) * dx
        a += 1j * omega * sigma * InnerProduct(A, v) * dx
        # Small regularization for gauge (only needed in non-conducting regions)
        a += 1e-12 * InnerProduct(A, v) * dx
        a.Assemble()

        # Linear form
        f = LinearForm(fes)
        f += InnerProduct(J0, v) * dx
        f.Assemble()

        # Solve
        gfA = GridFunction(fes)

        if verbose:
            print("  Assembling and solving...")

        # Try direct solver for small problems, iterative for large
        ndof = fes.ndof
        if ndof < 150000:
            if verbose:
                print(f"  Using direct solver (ndof={ndof})")
            try:
                gfA.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="umfpack") * f.vec
            except:
                if verbose:
                    print("  UMFPACK failed, trying pardiso...")
                gfA.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec
        else:
            if verbose:
                print(f"  Using iterative solver (ndof={ndof})")
            from ngsolve.krylovspace import GMResSolver
            from ngsolve import TaskManager
            pre = a.mat.CreateSmoother(fes.FreeDofs())
            solver = GMResSolver(a.mat, pre, maxiter=5000, tol=1e-8, printrates=verbose)
            gfA.vec.data = solver * f.vec

        if verbose:
            print("  Solution completed.")

        # Compute derived quantities
        B = curl(gfA)
        J_eddy = -1j * omega * sigma * gfA

        return gfA, B, J_eddy


def compute_measurements(mesh, B, J_eddy, line_params, verbose=True):
    """
    Compute field values along a measurement line.

    Args:
        mesh: NGSolve mesh
        B: Magnetic flux density CoefficientFunction
        J_eddy: Eddy current density CoefficientFunction
        line_params: Dictionary with 'x_values', 'y', 'z'
        verbose: Print results

    Returns:
        results: Dictionary with x, Bz_real, Bz_imag, Jy_real, Jy_imag
    """
    x_values = line_params['x_values']
    y_val = line_params['y']
    z_val = line_params['z']

    Bz_real = []
    Bz_imag = []
    Jy_real = []
    Jy_imag = []

    for x_mm in x_values:
        x_val = x_mm * 1e-3
        try:
            pnt = mesh(x_val, y_val, z_val)
            B_val = B(pnt)
            J_val = J_eddy(pnt)

            Bz_real.append(B_val[2].real)
            Bz_imag.append(B_val[2].imag)
            Jy_real.append(J_val[1].real)
            Jy_imag.append(J_val[1].imag)
        except:
            Bz_real.append(float('nan'))
            Bz_imag.append(float('nan'))
            Jy_real.append(float('nan'))
            Jy_imag.append(float('nan'))

    results = {
        'x': x_values,
        'Bz_real': Bz_real,
        'Bz_imag': Bz_imag,
        'Jy_real': Jy_real,
        'Jy_imag': Jy_imag
    }

    if verbose:
        print(f"\nMeasurement line: y={y_val*1000:.0f}mm, z={z_val*1000:.0f}mm")
        print(f"{'x[mm]':>8s} {'Bz_re[mT]':>12s} {'Bz_im[mT]':>12s} {'Jy_re[MA/m2]':>14s} {'Jy_im[MA/m2]':>14s}")
        print("-" * 65)
        for i, x_mm in enumerate(x_values):
            Bz_re = Bz_real[i] * 1000 if not np.isnan(Bz_real[i]) else float('nan')
            Bz_im = Bz_imag[i] * 1000 if not np.isnan(Bz_imag[i]) else float('nan')
            Jy_re = Jy_real[i] / 1e6 if not np.isnan(Jy_real[i]) else float('nan')
            Jy_im = Jy_imag[i] / 1e6 if not np.isnan(Jy_imag[i]) else float('nan')
            print(f"{x_mm:8.0f} {Bz_re:12.4f} {Bz_im:12.4f} {Jy_re:14.4f} {Jy_im:14.4f}")

    return results


def compute_power_loss(mesh, J_eddy, sigma):
    """
    Compute time-averaged Joule power loss in conductor.

    P = 0.5 * Re(∫ J* · J / sigma dV)

    Args:
        mesh: NGSolve mesh
        J_eddy: Complex eddy current density
        sigma: Conductivity CoefficientFunction

    Returns:
        P_loss: Power loss [W]
    """
    # J_eddy is complex, J* is conjugate
    # P = 0.5 * ∫ |J|^2 / sigma dV
    integrand = InnerProduct(J_eddy, Conj(J_eddy)) / (sigma + 1e-10)
    P_loss = 0.5 * Integrate(integrand, mesh, definedon=mesh.Materials("plate")).real

    return P_loss


def compute_magnetic_energy(mesh, B, inv_mu):
    """
    Compute time-averaged magnetic energy.

    W = 0.25 * Re(∫ B* · (1/mu) · B dV)

    Args:
        mesh: NGSolve mesh
        B: Complex magnetic flux density
        inv_mu: 1/mu CoefficientFunction

    Returns:
        W_mag: Magnetic energy [J]
    """
    integrand = inv_mu * InnerProduct(B, Conj(B))
    W_mag = 0.25 * Integrate(integrand, mesh).real

    return W_mag


if __name__ == "__main__":
    print("=" * 70)
    print("TEAM Problem 7: A-method Solver")
    print("=" * 70)

    # Create geometry (use coarser mesh to avoid memory issues)
    mesh = create_team7_geometry(maxh=0.08)  # 80 mm mesh for faster testing

    # Get material properties
    mu, sigma, inv_mu = get_material_properties(mesh)

    # ============================================================
    # Solve at 50 Hz
    # ============================================================
    print("\n" + "=" * 70)
    print("Analysis at 50 Hz")
    print("=" * 70)

    gfA_50, B_50, J_eddy_50 = solve_time_harmonic(mesh, freq_50, order=1)

    # Compute global quantities
    P_loss_50 = compute_power_loss(mesh, J_eddy_50, sigma)
    W_mag_50 = compute_magnetic_energy(mesh, B_50, inv_mu)

    print(f"\nGlobal results at 50 Hz:")
    print(f"  Joule power loss: {P_loss_50:.4f} W")
    print(f"  Magnetic energy:  {W_mag_50:.6e} J")

    # Measurements along line A1-B1
    results_A1B1_50 = compute_measurements(mesh, B_50, J_eddy_50, measurement_line_A1B1)

    # ============================================================
    # Solve at 200 Hz
    # ============================================================
    print("\n" + "=" * 70)
    print("Analysis at 200 Hz")
    print("=" * 70)

    gfA_200, B_200, J_eddy_200 = solve_time_harmonic(mesh, freq_200, order=1)

    # Compute global quantities
    P_loss_200 = compute_power_loss(mesh, J_eddy_200, sigma)
    W_mag_200 = compute_magnetic_energy(mesh, B_200, inv_mu)

    print(f"\nGlobal results at 200 Hz:")
    print(f"  Joule power loss: {P_loss_200:.4f} W")
    print(f"  Magnetic energy:  {W_mag_200:.6e} J")

    # Measurements along line A1-B1
    results_A1B1_200 = compute_measurements(mesh, B_200, J_eddy_200, measurement_line_A1B1)

    # ============================================================
    # Summary
    # ============================================================
    print("\n" + "=" * 70)
    print("SUMMARY: TEAM Problem 7 Results")
    print("=" * 70)
    print(f"""
Mesh elements: {mesh.ne}
Polynomial order: 2

Frequency    Joule Loss [W]    Magnetic Energy [J]
--------------------------------------------------
  50 Hz      {P_loss_50:12.4f}      {W_mag_50:.6e}
 200 Hz      {P_loss_200:12.4f}      {W_mag_200:.6e}

Note: Compare Bz values along line A1-B1 with experimental data
from TEAM Problem 7 benchmark.
""")
