"""
Coil Impedance Analysis using PEEC Method

This example demonstrates the PEEC (Partial Element Equivalent Circuit) solver
for computing coil impedance vs frequency.

Features demonstrated:
- Creating a circular loop coil conductor
- Setting conductor properties (conductivity, mu_r)
- Frequency sweep for impedance analysis
- Skin effect (ESIM surface impedance)
- VTU export for ParaView visualization

Physical model:
- Circular loop coil (single turn for simplicity)
- Copper conductor with skin effect
- Laplace kernel (quasi-static approximation)

Output:
- Impedance vs frequency plot
- Resistance and inductance extraction
- Skin depth analysis
- VTU files for ParaView (coil geometry + B field grid)

Part of Radia PEEC integration examples.
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt

# Add Radia to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))

# Physical constants
MU_0 = 4 * np.pi * 1e-7  # H/m
SIGMA_COPPER = 5.8e7     # S/m


def skin_depth(freq, sigma, mu_r=1.0):
    """Calculate skin depth [m]."""
    if freq <= 0 or sigma <= 0:
        return float('inf')
    omega = 2 * np.pi * freq
    return np.sqrt(2 / (omega * MU_0 * mu_r * sigma))


def surface_impedance(freq, sigma, mu_r=1.0):
    """Calculate surface impedance Zs = (1+j) * Rs."""
    delta = skin_depth(freq, sigma, mu_r)
    if delta == float('inf'):
        return 0 + 0j
    Rs = 1 / (sigma * delta)
    return (1 + 1j) * Rs


def analytical_loop_inductance(radius, wire_radius):
    """
    Analytical inductance of a circular loop.

    L = mu_0 * R * (ln(8*R/a) - 2)

    where R = loop radius, a = wire radius
    """
    return MU_0 * radius * (np.log(8 * radius / wire_radius) - 2)


def analytical_dc_resistance(radius, wire_radius, sigma):
    """
    DC resistance of circular loop wire.

    R_DC = (2*pi*R) / (sigma * pi * a^2)
    """
    wire_length = 2 * np.pi * radius
    wire_area = np.pi * wire_radius**2
    return wire_length / (sigma * wire_area)


def analytical_ac_resistance(freq, radius, wire_radius, sigma, mu_r=1.0):
    """
    AC resistance considering skin effect.

    For delta << a: R_AC = R_DC * (a / (2*delta))
    """
    R_dc = analytical_dc_resistance(radius, wire_radius, sigma)
    delta = skin_depth(freq, sigma, mu_r)

    if delta >= wire_radius:
        return R_dc
    else:
        # Skin effect regime
        return R_dc * wire_radius / (2 * delta)


class SimpleCoilPEEC:
    """
    Simple PEEC model for a circular coil.

    This is a simplified model that demonstrates the PEEC concept
    without requiring the full C++ PEEC solver.
    """

    def __init__(self, loop_radius, wire_radius, sigma=SIGMA_COPPER, mu_r=1.0):
        """
        Parameters:
        -----------
        loop_radius : float
            Radius of the coil loop [m]
        wire_radius : float
            Radius of the wire cross-section [m]
        sigma : float
            Electrical conductivity [S/m]
        mu_r : float
            Relative permeability of conductor
        """
        self.R = loop_radius
        self.a = wire_radius
        self.sigma = sigma
        self.mu_r = mu_r

        # Analytical values for reference
        self.L_analytical = analytical_loop_inductance(loop_radius, wire_radius)
        self.R_dc = analytical_dc_resistance(loop_radius, wire_radius, sigma)

        print(f"Coil parameters:")
        print(f"  Loop radius: {loop_radius*1000:.2f} mm")
        print(f"  Wire radius: {wire_radius*1000:.3f} mm")
        print(f"  Conductivity: {sigma:.2e} S/m")
        print(f"  Analytical L: {self.L_analytical*1e9:.3f} nH")
        print(f"  DC Resistance: {self.R_dc*1000:.4f} mOhm")

    def impedance(self, freq):
        """
        Calculate impedance at given frequency.

        Z = R(f) + j*omega*L

        where R(f) includes skin effect via ESIM surface impedance.
        """
        omega = 2 * np.pi * freq

        # Inductance (approximately constant)
        L = self.L_analytical

        # Resistance with skin effect
        R = analytical_ac_resistance(freq, self.R, self.a, self.sigma, self.mu_r)

        return R + 1j * omega * L

    def frequency_sweep(self, freqs):
        """Compute impedance over frequency range."""
        return np.array([self.impedance(f) for f in freqs])


def run_coil_analysis():
    """Run coil impedance analysis and generate plots."""

    print("="*60)
    print("PEEC Coil Impedance Analysis")
    print("="*60)

    # Coil parameters
    loop_radius = 0.05      # 50 mm
    wire_radius = 0.5e-3    # 0.5 mm (AWG 24 approx)

    # Create coil model
    coil = SimpleCoilPEEC(loop_radius, wire_radius)

    # Frequency range: 100 Hz to 100 MHz
    freqs = np.logspace(2, 8, 200)

    # Compute impedance
    Z = coil.frequency_sweep(freqs)
    R = np.real(Z)
    X = np.imag(Z)
    L = X / (2 * np.pi * freqs)

    # Skin depth at various frequencies
    deltas = np.array([skin_depth(f, coil.sigma) for f in freqs])

    # Create figure with 4 subplots
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Plot 1: Impedance magnitude
    ax1 = axes[0, 0]
    ax1.loglog(freqs, np.abs(Z), 'b-', linewidth=2, label='|Z|')
    ax1.loglog(freqs, R, 'r--', linewidth=1.5, label='R (resistance)')
    ax1.loglog(freqs, X, 'g--', linewidth=1.5, label='X (reactance)')
    ax1.set_xlabel('Frequency [Hz]')
    ax1.set_ylabel('Impedance [Ohm]')
    ax1.set_title('Coil Impedance vs Frequency')
    ax1.legend()
    ax1.grid(True, which='both', alpha=0.3)

    # Plot 2: Resistance
    ax2 = axes[0, 1]
    ax2.loglog(freqs, R * 1000, 'r-', linewidth=2)
    ax2.axhline(coil.R_dc * 1000, color='k', linestyle='--',
                label=f'DC: {coil.R_dc*1000:.3f} mOhm')
    ax2.set_xlabel('Frequency [Hz]')
    ax2.set_ylabel('Resistance [mOhm]')
    ax2.set_title('AC Resistance (includes skin effect)')
    ax2.legend()
    ax2.grid(True, which='both', alpha=0.3)

    # Plot 3: Inductance
    ax3 = axes[1, 0]
    ax3.semilogx(freqs, L * 1e9, 'b-', linewidth=2)
    ax3.axhline(coil.L_analytical * 1e9, color='k', linestyle='--',
                label=f'Analytical: {coil.L_analytical*1e9:.2f} nH')
    ax3.set_xlabel('Frequency [Hz]')
    ax3.set_ylabel('Inductance [nH]')
    ax3.set_title('Effective Inductance')
    ax3.legend()
    ax3.grid(True, which='both', alpha=0.3)
    ax3.set_ylim([0, coil.L_analytical * 1e9 * 1.5])

    # Plot 4: Skin depth
    ax4 = axes[1, 1]
    ax4.loglog(freqs, deltas * 1000, 'g-', linewidth=2)
    ax4.axhline(wire_radius * 1000, color='r', linestyle='--',
                label=f'Wire radius: {wire_radius*1000:.2f} mm')
    ax4.set_xlabel('Frequency [Hz]')
    ax4.set_ylabel('Skin Depth [mm]')
    ax4.set_title('Skin Depth in Copper')
    ax4.legend()
    ax4.grid(True, which='both', alpha=0.3)

    # Find transition frequency (delta = wire_radius)
    f_transition = 2 / (MU_0 * coil.sigma * wire_radius**2) / (2 * np.pi)
    ax4.axvline(f_transition, color='orange', linestyle=':',
                label=f'Transition: {f_transition/1000:.1f} kHz')
    ax4.legend()

    plt.tight_layout()

    # Save figure
    output_path = os.path.join(os.path.dirname(__file__), 'coil_impedance_peec.png')
    plt.savefig(output_path, dpi=150)
    print(f"\nSaved: {output_path}")
    plt.close()

    # Print summary table
    print("\n" + "="*60)
    print("Impedance at Selected Frequencies")
    print("="*60)
    print(f"{'Freq':<12} {'|Z| [Ohm]':<12} {'R [mOhm]':<12} {'L [nH]':<12} {'delta [mm]':<12}")
    print("-"*60)

    test_freqs = [100, 1000, 10000, 100000, 1e6, 10e6]
    for f in test_freqs:
        z = coil.impedance(f)
        delta = skin_depth(f, coil.sigma)
        l = z.imag / (2 * np.pi * f)

        if f >= 1e6:
            freq_str = f"{f/1e6:.0f} MHz"
        elif f >= 1000:
            freq_str = f"{f/1000:.0f} kHz"
        else:
            freq_str = f"{f:.0f} Hz"

        print(f"{freq_str:<12} {abs(z):<12.4e} {z.real*1000:<12.4f} {l*1e9:<12.2f} {delta*1000:<12.4f}")

    print("="*60)

    # Quality factor analysis
    print("\nQuality Factor (Q = omega*L/R)")
    print("-"*40)
    for f in [1000, 10000, 100000, 1e6]:
        z = coil.impedance(f)
        Q = z.imag / z.real
        if f >= 1e6:
            freq_str = f"{f/1e6:.0f} MHz"
        elif f >= 1000:
            freq_str = f"{f/1000:.0f} kHz"
        else:
            freq_str = f"{f:.0f} Hz"
        print(f"  Q @ {freq_str}: {Q:.1f}")

    return coil


def export_coil_to_vtu(coil, filename_base, num_segments=60):
    """
    Export coil geometry and B field to VTU files.

    Parameters:
    -----------
    coil : SimpleCoilPEEC
        Coil object
    filename_base : str
        Base filename for VTU output
    num_segments : int
        Number of segments for coil discretization
    """
    print("\n" + "="*60)
    print("Exporting to VTU format")
    print("="*60)

    # Generate coil centerline points
    theta = np.linspace(0, 2*np.pi, num_segments + 1)
    coil_points = np.column_stack([
        coil.R * np.cos(theta),
        coil.R * np.sin(theta),
        np.zeros_like(theta)
    ])

    # Generate wire cross-section (tube around centerline)
    n_around = 8  # Points around wire circumference
    phi = np.linspace(0, 2*np.pi, n_around + 1)[:-1]

    all_points = []
    all_cells = []

    for i in range(num_segments):
        # Tangent, normal, binormal at this point
        t0 = theta[i]
        t1 = theta[i + 1]

        for j, t in enumerate([t0, t1]):
            # Center of wire at this angle
            cx = coil.R * np.cos(t)
            cy = coil.R * np.sin(t)
            cz = 0.0

            # Local coordinate system (radial, z, tangent)
            radial = np.array([np.cos(t), np.sin(t), 0])
            z_dir = np.array([0, 0, 1])

            # Generate points around wire
            for p in phi:
                # Position on wire surface
                offset = coil.a * (np.cos(p) * radial + np.sin(p) * z_dir)
                point = np.array([cx, cy, cz]) + offset
                all_points.append(point)

        # Create hexahedral cells connecting two rings
        base_idx = len(all_points) - 2 * n_around
        for k in range(n_around):
            k_next = (k + 1) % n_around
            # 8 vertices of hexahedron
            v0 = base_idx + k
            v1 = base_idx + k_next
            v2 = base_idx + n_around + k_next
            v3 = base_idx + n_around + k
            all_cells.append([v0, v1, v2, v3])

    # Write coil geometry VTU
    coil_vtu = f"{filename_base}_coil.vtu"
    write_coil_vtu(all_points, all_cells, coil_vtu)

    # Generate B field on a grid
    print("\nGenerating B field grid...")
    nx, ny, nz = 21, 21, 21
    x = np.linspace(-0.1, 0.1, nx)
    y = np.linspace(-0.1, 0.1, ny)
    z = np.linspace(-0.05, 0.05, nz)

    B_field = np.zeros((nx, ny, nz, 3))

    # Compute B field using Biot-Savart for circular loop
    for i, xi in enumerate(x):
        for j, yj in enumerate(y):
            for k, zk in enumerate(z):
                B = compute_loop_bfield(coil.R, [xi, yj, zk])
                B_field[i, j, k, :] = B

    # Write B field VTS
    field_vts = f"{filename_base}_field.vts"
    write_field_vts(x, y, z, B_field, field_vts)

    print(f"\nVTU Export complete:")
    print(f"  Coil geometry: {coil_vtu}")
    print(f"  B field grid:  {field_vts}")
    print("\nOpen these files in ParaView for visualization.")


def compute_loop_bfield(R, point, I=1.0):
    """
    Compute B field from circular loop using Biot-Savart.

    Parameters:
    -----------
    R : float
        Loop radius [m]
    point : array-like
        Evaluation point [x, y, z] [m]
    I : float
        Current [A] (default 1.0 for normalized field)

    Returns:
    --------
    B : ndarray
        B field [Bx, By, Bz] [T]
    """
    x, y, z = point
    rho = np.sqrt(x**2 + y**2)

    # On-axis field (special case for rho=0)
    if rho < 1e-10:
        Bz = MU_0 * I * R**2 / (2 * (R**2 + z**2)**1.5)
        return np.array([0.0, 0.0, Bz])

    # Off-axis: Use elliptic integral approximation
    # B_z = (mu_0 * I / (2 * pi)) * (1 / sqrt((R+rho)^2 + z^2)) * [K(m) + ...]
    # Simplified approximation for visualization
    r_plus = np.sqrt((R + rho)**2 + z**2)
    r_minus = np.sqrt((R - rho)**2 + z**2)

    if r_minus < 1e-10:
        # Very close to wire, use cutoff
        r_minus = 0.001

    # Approximate field using dipole + near-field correction
    m = 4 * R * rho / ((R + rho)**2 + z**2)  # Elliptic parameter

    # Clamp m to avoid singularity at m=1 (on the wire)
    m = min(m, 0.999)

    alpha2 = (R - rho)**2 + z**2
    beta2 = (R + rho)**2 + z**2
    beta = np.sqrt(beta2)

    # Avoid division by zero
    if alpha2 < 1e-20:
        alpha2 = 1e-20

    # Use simple approximation for K(m) and E(m)
    # For accurate results, use scipy.special.ellipk, ellipe
    if m < 0.9:
        K_approx = np.pi/2 * (1 + m/4 + 9*m**2/64)
        E_approx = np.pi/2 * (1 - m/4 - 3*m**2/64)
    else:
        # Near the wire, use limiting form
        K_approx = np.log(4/np.sqrt(1-m + 1e-15))
        E_approx = 1.0

    C = MU_0 * I / np.pi

    # Bz component
    Bz = C / (2 * beta) * (K_approx + (R**2 - rho**2 - z**2) / alpha2 * E_approx)

    # Brho component (radial)
    if rho > 1e-10:
        Brho = C * z / (2 * rho * beta) * (-K_approx + (R**2 + rho**2 + z**2) / alpha2 * E_approx)
    else:
        Brho = 0.0

    # Convert to Cartesian
    if rho > 1e-10:
        Bx = Brho * x / rho
        By = Brho * y / rho
    else:
        Bx = 0.0
        By = 0.0

    return np.array([Bx, By, Bz])


def write_coil_vtu(points, cells, filename):
    """Write coil geometry to VTU file."""
    n_points = len(points)
    n_cells = len(cells)

    with open(filename, 'w') as f:
        f.write('<?xml version="1.0"?>\n')
        f.write('<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">\n')
        f.write('  <UnstructuredGrid>\n')
        f.write(f'    <Piece NumberOfPoints="{n_points}" NumberOfCells="{n_cells}">\n')

        # Points
        f.write('      <Points>\n')
        f.write('        <DataArray type="Float64" NumberOfComponents="3" format="ascii">\n')
        for p in points:
            f.write(f'          {p[0]:.10e} {p[1]:.10e} {p[2]:.10e}\n')
        f.write('        </DataArray>\n')
        f.write('      </Points>\n')

        # Cells
        f.write('      <Cells>\n')
        # Connectivity
        f.write('        <DataArray type="Int32" Name="connectivity" format="ascii">\n')
        for cell in cells:
            f.write('          ' + ' '.join(str(v) for v in cell) + '\n')
        f.write('        </DataArray>\n')
        # Offsets
        f.write('        <DataArray type="Int32" Name="offsets" format="ascii">\n')
        offsets = [(i + 1) * 4 for i in range(n_cells)]
        f.write('          ' + ' '.join(str(o) for o in offsets) + '\n')
        f.write('        </DataArray>\n')
        # Cell types (9 = VTK_QUAD)
        f.write('        <DataArray type="UInt8" Name="types" format="ascii">\n')
        f.write('          ' + ' '.join(['9'] * n_cells) + '\n')
        f.write('        </DataArray>\n')
        f.write('      </Cells>\n')

        f.write('    </Piece>\n')
        f.write('  </UnstructuredGrid>\n')
        f.write('</VTKFile>\n')

    print(f"Exported: {filename} ({n_cells} cells, {n_points} points)")


def write_field_vts(x, y, z, B_field, filename):
    """Write B field on structured grid to VTS file."""
    nx, ny, nz = len(x), len(y), len(z)

    with open(filename, 'w') as f:
        f.write('<?xml version="1.0"?>\n')
        f.write('<VTKFile type="StructuredGrid" version="0.1" byte_order="LittleEndian">\n')
        f.write(f'  <StructuredGrid WholeExtent="0 {nx-1} 0 {ny-1} 0 {nz-1}">\n')
        f.write(f'    <Piece Extent="0 {nx-1} 0 {ny-1} 0 {nz-1}">\n')

        # Points
        f.write('      <Points>\n')
        f.write('        <DataArray type="Float64" NumberOfComponents="3" format="ascii">\n')
        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    f.write(f'          {x[i]:.10e} {y[j]:.10e} {z[k]:.10e}\n')
        f.write('        </DataArray>\n')
        f.write('      </Points>\n')

        # Point data
        f.write('      <PointData>\n')
        # B field vector
        f.write('        <DataArray type="Float64" Name="B_field" NumberOfComponents="3" format="ascii">\n')
        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    Bx, By, Bz = B_field[i, j, k, :]
                    f.write(f'          {Bx:.10e} {By:.10e} {Bz:.10e}\n')
        f.write('        </DataArray>\n')
        # B magnitude
        f.write('        <DataArray type="Float64" Name="B_magnitude" format="ascii">\n')
        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    B_mag = np.linalg.norm(B_field[i, j, k, :])
                    f.write(f'          {B_mag:.10e}\n')
        f.write('        </DataArray>\n')
        f.write('      </PointData>\n')

        f.write('    </Piece>\n')
        f.write('  </StructuredGrid>\n')
        f.write('</VTKFile>\n')

    print(f"Exported: {filename} ({nx}x{ny}x{nz} = {nx*ny*nz} points)")


def export_radia_coil_vts(loop_radius, current=1.0):
    """
    Create a circular coil using Radia ObjArcCur and export B field using FldVTS.

    Parameters:
    -----------
    loop_radius : float
        Radius of the coil loop [m]
    current : float
        Coil current [A]
    """
    try:
        import radia as rad
    except ImportError:
        print("Radia not available. Skipping Radia-based VTS export.")
        return

    print("\n" + "="*60)
    print("Exporting B field using Radia FldVTS")
    print("="*60)

    rad.FldUnits('m')

    # Create circular coil using multiple arc segments
    # ObjArcCur(center, radii, phi_range, nseg, current, 'man'|'auto')
    # center: [x, y, z]
    # radii: [r_inner, r_outer] (for thin wire: same value or small difference)
    # phi_range: [phi_start, phi_end] in degrees
    # nseg: number of segments
    # current: current in Amperes
    # 'man'|'auto': manual or automatic subdivision

    wire_radius = 0.001  # 1 mm wire radius for visualization
    r_inner = loop_radius - wire_radius
    r_outer = loop_radius + wire_radius

    # Create full circle with 36 segments
    # ObjArcCur(center, radii, phi_range, height, nseg, current_density, 'man'|'auto')
    # phi_range is in radians: [0, 2*pi) for full circle
    # current_density = current / cross_section_area
    cross_section_area = (r_outer - r_inner) * 0.001  # radial * height
    current_density = current / cross_section_area

    coil = rad.ObjArcCur(
        [0, 0, 0],                  # center [x, y, z]
        [r_inner, r_outer],         # radii [r_inner, r_outer]
        [0, 2 * np.pi - 0.01],      # phi range [phi_min, phi_max] in radians
        0.001,                      # height (z-extent, thin coil in meters)
        36,                         # number of segments
        current_density,            # current density [A/m^2]
        'man'                       # manual subdivision
    )

    print(f"  Coil radius: {loop_radius*1000:.1f} mm")
    print(f"  Current: {current:.1f} A")

    # Export B field using FldVTS
    output_path = os.path.join(os.path.dirname(__file__), 'coil_bfield_radia.vts')

    # Grid parameters (in meters)
    x_range = [-0.1, 0.1]
    y_range = [-0.1, 0.1]
    z_range = [-0.05, 0.05]
    nx, ny, nz = 41, 41, 21

    print(f"\n  Computing B field on {nx}x{ny}x{nz} grid...")
    print(f"  Grid range: x=[{x_range[0]:.3f}, {x_range[1]:.3f}] m")
    print(f"              y=[{y_range[0]:.3f}, {y_range[1]:.3f}] m")
    print(f"              z=[{z_range[0]:.3f}, {z_range[1]:.3f}] m")

    # FldVTS(obj, filename, x_range, y_range, z_range, nx, ny, nz, include_B, include_H, unit_scale)
    rad.FldVTS(coil, output_path,
               x_range, y_range, z_range,
               nx, ny, nz,
               1,    # include B field
               0,    # do not include H field
               1.0)  # unit scale

    print(f"\n  Exported: {output_path}")
    print(f"  Grid points: {nx*ny*nz}")

    # Also compute and print field at center
    B_center = rad.Fld(coil, 'b', [0, 0, 0])
    print(f"\n  B field at center: Bz = {B_center[2]*1e6:.2f} uT")

    # Analytical comparison: Bz = mu_0 * I / (2 * R)
    Bz_analytical = MU_0 * current / (2 * loop_radius)
    print(f"  Analytical Bz:     Bz = {Bz_analytical*1e6:.2f} uT")

    return coil


def run_peec_solver_analysis():
    """
    Run coil impedance analysis using Radia's C++ PEEC solver.

    This function uses the actual PEEC solver (CndLoop, CndSolve, CndGetImpedance)
    instead of analytical formulas.
    """
    try:
        import radia as rad
    except ImportError:
        print("Radia not available. Skipping PEEC solver analysis.")
        return None

    print("\n" + "="*60)
    print("PEEC Solver Coil Impedance Analysis (C++ Solver)")
    print("="*60)

    rad.FldUnits('m')

    # Coil parameters
    loop_radius = 0.05      # 50 mm
    wire_width = 0.002      # 2 mm (rectangular cross-section width)
    wire_height = 0.001     # 1 mm (rectangular cross-section height)
    sigma = SIGMA_COPPER    # 5.8e7 S/m

    # Calculate equivalent wire radius for comparison
    wire_area = wire_width * wire_height
    wire_radius_eq = np.sqrt(wire_area / np.pi)

    print(f"\nCoil parameters:")
    print(f"  Loop radius:    {loop_radius*1000:.1f} mm")
    print(f"  Wire width:     {wire_width*1000:.1f} mm")
    print(f"  Wire height:    {wire_height*1000:.1f} mm")
    print(f"  Cross-section:  {wire_area*1e6:.2f} mm^2")
    print(f"  Eq. wire radius: {wire_radius_eq*1000:.3f} mm")
    print(f"  Conductivity:   {sigma:.2e} S/m")

    # Create PEEC loop conductor
    # CndLoop(center, radius, normal, cross_section, wire_width, wire_height, sigma, num_panels_around, num_panels_loop)
    # cross_section: 'r' = rectangular, 'c' = circular
    coil_handle = rad.CndLoop(
        [0, 0, 0],              # center [x, y, z]
        loop_radius,            # radius [m]
        [0, 0, 1],              # normal vector (coil in xy-plane)
        'r',                    # 'r' = rectangular cross-section
        wire_width,             # wire width [m]
        wire_height,            # wire height [m]
        sigma,                  # conductivity [S/m]
        8,                      # panels around wire circumference
        36                      # panels around loop circumference
    )

    print(f"\n  PEEC conductor handle: {coil_handle}")
    print(f"  Panels: 8 (around) x 36 (loop) = 288 panels")

    # Analytical reference values
    L_analytical = analytical_loop_inductance(loop_radius, wire_radius_eq)
    R_dc_analytical = analytical_dc_resistance(loop_radius, wire_radius_eq, sigma)

    print(f"\n  Analytical L:   {L_analytical*1e9:.3f} nH")
    print(f"  Analytical R_DC: {R_dc_analytical*1000:.4f} mOhm")

    # Frequency sweep
    frequencies = [100, 1000, 10000, 100000, 1e6, 10e6]

    print("\n" + "-"*70)
    print(f"{'Freq':<12} {'R_PEEC [mOhm]':<15} {'L_PEEC [nH]':<12} {'R_ana [mOhm]':<14} {'delta [mm]':<12}")
    print("-"*70)

    peec_results = []

    for freq in frequencies:
        # Set frequency for PEEC solve
        rad.CndSetFrequency(coil_handle, freq)

        # Set voltage excitation (1V)
        rad.CndSetVoltage(coil_handle, 1.0, 0.0)

        # Solve PEEC system
        rad.CndSolve(coil_handle)

        # Get impedance
        Z = rad.CndGetImpedance(coil_handle)  # Returns complex number

        omega = 2 * np.pi * freq
        R_peec = Z.real
        X_peec = Z.imag
        L_peec = X_peec / omega

        # Analytical comparison
        R_analytical = analytical_ac_resistance(freq, loop_radius, wire_radius_eq, sigma)
        delta = skin_depth(freq, sigma)

        peec_results.append({
            'freq': freq,
            'Z': Z,
            'R': R_peec,
            'L': L_peec,
            'R_analytical': R_analytical,
            'delta': delta
        })

        # Format frequency string
        if freq >= 1e6:
            freq_str = f"{freq/1e6:.0f} MHz"
        elif freq >= 1000:
            freq_str = f"{freq/1000:.0f} kHz"
        else:
            freq_str = f"{freq:.0f} Hz"

        print(f"{freq_str:<12} {R_peec*1000:<15.4f} {L_peec*1e9:<12.3f} "
              f"{R_analytical*1000:<14.4f} {delta*1000:<12.4f}")

    print("-"*70)

    # Plot comparison
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    freqs_array = np.array([r['freq'] for r in peec_results])
    R_peec_array = np.array([r['R'] for r in peec_results]) * 1000  # mOhm
    L_peec_array = np.array([r['L'] for r in peec_results]) * 1e9   # nH
    R_ana_array = np.array([r['R_analytical'] for r in peec_results]) * 1000  # mOhm

    # Plot resistance
    ax1 = axes[0]
    ax1.loglog(freqs_array, R_peec_array, 'bo-', linewidth=2, markersize=8, label='PEEC Solver')
    ax1.loglog(freqs_array, R_ana_array, 'r--', linewidth=2, label='Analytical')
    ax1.axhline(R_dc_analytical*1000, color='k', linestyle=':', label=f'DC: {R_dc_analytical*1000:.3f} mOhm')
    ax1.set_xlabel('Frequency [Hz]')
    ax1.set_ylabel('Resistance [mOhm]')
    ax1.set_title('AC Resistance: PEEC vs Analytical')
    ax1.legend()
    ax1.grid(True, which='both', alpha=0.3)

    # Plot inductance
    ax2 = axes[1]
    ax2.semilogx(freqs_array, L_peec_array, 'bo-', linewidth=2, markersize=8, label='PEEC Solver')
    ax2.axhline(L_analytical*1e9, color='r', linestyle='--', label=f'Analytical: {L_analytical*1e9:.2f} nH')
    ax2.set_xlabel('Frequency [Hz]')
    ax2.set_ylabel('Inductance [nH]')
    ax2.set_title('Inductance from PEEC Solver')
    ax2.legend()
    ax2.grid(True, which='both', alpha=0.3)
    ax2.set_ylim([0, L_analytical*1e9 * 1.5])

    plt.tight_layout()

    output_path = os.path.join(os.path.dirname(__file__), 'coil_impedance_peec_solver.png')
    plt.savefig(output_path, dpi=150)
    print(f"\nSaved: {output_path}")
    plt.close()

    # Quality factor at key frequencies
    print("\nQuality Factor (Q = omega*L/R) from PEEC Solver:")
    print("-"*40)
    for r in peec_results:
        if r['freq'] >= 1000:  # Only show 1 kHz and above
            Q = r['Z'].imag / r['Z'].real
            if r['freq'] >= 1e6:
                freq_str = f"{r['freq']/1e6:.0f} MHz"
            else:
                freq_str = f"{r['freq']/1000:.0f} kHz"
            print(f"  Q @ {freq_str}: {Q:.1f}")

    return coil_handle, peec_results


def main():
    """Main function."""
    print("\nPEEC Coil Impedance Example")
    print("This example demonstrates impedance analysis of a circular coil")
    print("using the PEEC (Partial Element Equivalent Circuit) approach.\n")

    # Run analytical model analysis
    coil = run_coil_analysis()

    # Export to VTU (Python implementation with Biot-Savart)
    output_base = os.path.join(os.path.dirname(__file__), 'coil_impedance_peec')
    export_coil_to_vtu(coil, output_base)

    # Export B field using Radia FldVTS (C++ implementation)
    export_radia_coil_vts(coil.R, current=1.0)

    # Run PEEC solver analysis (C++ PEEC solver)
    peec_result = run_peec_solver_analysis()

    print("\n" + "="*60)
    print("Analysis complete!")
    print("="*60)

    if peec_result is not None:
        print("\nComparison Summary:")
        print("  - Analytical model: Uses closed-form inductance/resistance formulas")
        print("  - PEEC Solver: Full electromagnetic simulation with skin effect")
        print("  - Both methods should give similar results for simple geometry")

    print("\nThe full PEEC solver (rad_peec_mmm_coupled.cpp) provides:")
    print("  - Arbitrary coil geometry support")
    print("  - Coupled coil-magnet analysis")
    print("  - Complex permeability for magnetic materials")

    print("\nOutput files for ParaView:")
    print("  - coil_impedance_peec_coil.vtu  : Coil geometry (wire tube)")
    print("  - coil_impedance_peec_field.vts : B field (Biot-Savart, Python)")
    print("  - coil_bfield_radia.vts         : B field (Radia ObjArcCur, C++)")
    print("  - coil_impedance_peec_solver.png: PEEC solver results")


if __name__ == '__main__':
    main()
