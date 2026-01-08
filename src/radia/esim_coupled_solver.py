"""
ESIM Coupled Solver for Induction Heating Analysis

This module implements the coupled solver that combines:
1. FastImp coil (conductor model with eddy currents)
2. ESIM workpiece (nonlinear ferromagnetic material)

The solver uses fixed-point iteration to handle the nonlinear material behavior
where the surface impedance Z depends on the tangential field magnitude |H_t|.

Reference:
    K. Hollaus, M. Kaltenbacher, J. Schoberl, "A Nonlinear Effective Surface
    Impedance in a Magnetic Scalar Potential Formulation," IEEE Trans. Magnetics,
    2025, DOI: 10.1109/TMAG.2025.3613932

Author: Radia Development Team
Date: 2026-01-08
"""

import numpy as np
from scipy.constants import mu_0

try:
    from .esim_workpiece import ESIMWorkpiece, create_esim_block, create_esim_cylinder
    from .esim_cell_problem import ESITable, generate_esi_table_from_bh_curve
except ImportError:
    from esim_workpiece import ESIMWorkpiece, create_esim_block, create_esim_cylinder
    from esim_cell_problem import ESITable, generate_esi_table_from_bh_curve


class InductionHeatingCoil:
    """
    Wrapper class for induction heating coil using Radia's FastImp conductor API.

    This class provides a Python interface to create and analyze spiral/loop coils
    for induction heating applications.
    """

    def __init__(self, coil_type='spiral', **kwargs):
        """
        Initialize an induction heating coil.

        Parameters:
            coil_type: 'spiral', 'loop', or 'custom'

        For 'spiral' coil:
            center: [x, y, z] center coordinates [m]
            inner_radius: Inner radius [m]
            outer_radius: Outer radius [m]
            pitch: Height per turn [m]
            num_turns: Number of turns
            axis: [ax, ay, az] coil axis direction
            wire_width: Wire width [m]
            wire_height: Wire height [m] (optional, for rectangular)
            cross_section: 'r' (rectangular) or 'c' (circular)
            conductivity: Conductivity [S/m] (default: 5.8e7 for copper)
            num_panels_around: Panels around wire cross-section

        For 'loop' coil:
            center: [x, y, z] center coordinates [m]
            radius: Loop radius [m]
            normal: [nx, ny, nz] loop normal direction
            wire_width: Wire width [m]
            wire_height: Wire height [m] (optional)
            cross_section: 'r' or 'c'
            conductivity: [S/m]
            num_panels_around: Panels around wire
            num_panels_loop: Panels around loop circumference
        """
        self.coil_type = coil_type
        self.params = kwargs
        self.handle = None
        self.frequency = None
        self.current = 1.0  # Default 1 A

        # Use analytical model for now (FastImp integration requires separate solver)
        # FastImp conductor API computes impedance, not DC/low-freq field from current
        # For induction heating coil field, we use Biot-Savart analytical model
        self._rad = None
        self._create_analytical_model()

    def _create_coil(self):
        """Create the coil using Radia's FastImp API."""
        if self._rad is None:
            return

        rad = self._rad

        if self.coil_type == 'spiral':
            self.handle = rad.CndSpiral(
                self.params.get('center', [0, 0, 0]),
                self.params.get('inner_radius', 0.02),
                self.params.get('outer_radius', 0.05),
                self.params.get('pitch', 0.005),
                self.params.get('num_turns', 5),
                self.params.get('axis', [0, 0, 1]),
                self.params.get('cross_section', 'r'),
                self.params.get('wire_width', 0.003),
                self.params.get('wire_height', 0.002),
                self.params.get('conductivity', 5.8e7),
                self.params.get('num_panels_around', 8)
            )
        elif self.coil_type == 'loop':
            self.handle = rad.CndLoop(
                self.params.get('center', [0, 0, 0]),
                self.params.get('radius', 0.05),
                self.params.get('normal', [0, 0, 1]),
                self.params.get('cross_section', 'r'),
                self.params.get('wire_width', 0.003),
                self.params.get('wire_height', 0.002),
                self.params.get('conductivity', 5.8e7),
                self.params.get('num_panels_around', 8),
                self.params.get('num_panels_loop', 36)
            )

    def _create_analytical_model(self):
        """Create analytical coil model for testing without radia."""
        # Store coil geometry for analytical field calculation
        if self.coil_type == 'spiral':
            self.center = np.array(self.params.get('center', [0, 0, 0]))
            self.inner_radius = self.params.get('inner_radius', 0.02)
            self.outer_radius = self.params.get('outer_radius', 0.05)
            self.pitch = self.params.get('pitch', 0.005)
            self.num_turns = self.params.get('num_turns', 5)
            self.axis = np.array(self.params.get('axis', [0, 0, 1]))
            self.axis = self.axis / np.linalg.norm(self.axis)
        elif self.coil_type == 'loop':
            self.center = np.array(self.params.get('center', [0, 0, 0]))
            self.radius = self.params.get('radius', 0.05)
            self.normal = np.array(self.params.get('normal', [0, 0, 1]))
            self.normal = self.normal / np.linalg.norm(self.normal)

    def set_frequency(self, frequency):
        """Set the operating frequency [Hz]."""
        self.frequency = frequency
        if self._rad is not None and self.handle is not None:
            self._rad.CndSetFrequency(self.handle, frequency)

    def set_current(self, current):
        """Set the coil current [A]."""
        self.current = current

    def compute_field_at_point(self, point):
        """
        Compute B field at a single point.

        Parameters:
            point: [x, y, z] coordinates [m]

        Returns:
            B: Complex B field [Bx, By, Bz] in Tesla
        """
        if self._rad is not None and self.handle is not None:
            # Use FastImp for field computation
            B_complex = self._rad.CndFld(self.handle, 'b', point)
            # CndFld returns [Bx_re, By_re, Bz_re, Bx_im, By_im, Bz_im]
            if len(B_complex) == 6:
                return np.array([
                    B_complex[0] + 1j * B_complex[3],
                    B_complex[1] + 1j * B_complex[4],
                    B_complex[2] + 1j * B_complex[5]
                ])
            else:
                return np.array(B_complex[:3])
        else:
            # Use analytical model
            return self._compute_field_analytical(point)

    def _compute_field_analytical(self, point):
        """
        Compute B field using analytical formulas (Biot-Savart for circular loops).

        For spiral coil, approximates as stack of circular loops.
        """
        point = np.array(point)

        if self.coil_type == 'loop':
            return self._biot_savart_loop(point, self.center, self.radius,
                                          self.normal, self.current)
        elif self.coil_type == 'spiral':
            # Sum contribution from each turn
            B_total = np.zeros(3, dtype=complex)

            for i in range(self.num_turns):
                # Position along spiral
                t = i / max(self.num_turns - 1, 1)
                R = self.inner_radius + t * (self.outer_radius - self.inner_radius)
                z_offset = i * self.pitch

                # Turn center
                turn_center = self.center + z_offset * self.axis

                # Add contribution from this turn
                B_turn = self._biot_savart_loop(point, turn_center, R,
                                                self.axis, self.current)
                B_total += B_turn

            return B_total

    def _biot_savart_loop(self, point, center, radius, normal, current):
        """
        Compute B field from a circular current loop using Biot-Savart law.

        Uses the analytical formula for on-axis field and approximation for off-axis.
        """
        # Vector from loop center to point
        r = point - center

        # Component along loop axis
        z = np.dot(r, normal)

        # Perpendicular distance from axis
        r_perp = r - z * normal
        rho = np.linalg.norm(r_perp)

        # Simple on-axis formula (accurate for rho << radius)
        denom = (radius**2 + z**2)**(3/2)
        if denom < 1e-20:
            return np.zeros(3, dtype=complex)

        # Bz on axis
        Bz = mu_0 * current * radius**2 / (2 * denom)

        # Off-axis correction (first order)
        if rho > 1e-10 and radius > 1e-10:
            # Radial component (approximate)
            Br = 3 * mu_0 * current * radius**2 * z * rho / (4 * denom * (radius**2 + z**2))

            # Unit radial direction
            if rho > 1e-10:
                r_hat = r_perp / rho
            else:
                r_hat = np.array([1, 0, 0])

            B = Bz * normal + Br * r_hat
        else:
            B = Bz * normal

        return B.astype(complex)

    def compute_field_batch(self, points):
        """
        Compute B field at multiple points.

        Parameters:
            points: List of [x, y, z] coordinates [m]

        Returns:
            B_list: List of complex B field vectors
        """
        return [self.compute_field_at_point(p) for p in points]

    def compute_tangential_field(self, point, normal):
        """
        Compute tangential magnetic field H_t at a point on a surface.

        Parameters:
            point: [x, y, z] coordinates [m]
            normal: [nx, ny, nz] surface normal (outward)

        Returns:
            H_t: Complex tangential H field magnitude [A/m]
        """
        B = self.compute_field_at_point(point)
        H = B / mu_0  # H = B/mu_0 in air

        # Tangential component: H_t = H - (H . n) * n
        normal = np.array(normal)
        normal = normal / np.linalg.norm(normal)

        H_normal = np.dot(H, normal) * normal
        H_tangential = H - H_normal

        # Return magnitude of tangential field
        return np.linalg.norm(H_tangential)

    @property
    def num_panels(self):
        """Get number of surface panels in the coil model."""
        if self._rad is not None and self.handle is not None:
            return self._rad.CndNumPanels(self.handle)
        return 0


class ESIMCoupledSolver:
    """
    Coupled solver for induction heating with ESIM workpiece.

    This solver combines:
    1. FastImp-based coil model (or analytical model)
    2. ESIM workpiece with nonlinear surface impedance

    The coupling is achieved through fixed-point iteration:
    1. Compute B field from coil at workpiece surface
    2. Extract tangential H field
    3. Look up Z(|H_t|) from ESI table
    4. Update workpiece impedance
    5. Repeat until convergence

    Impedance Calculation:
    The total coil impedance seen from the power source includes:
    - Z_coil = R_coil + j*omega*L_coil  (coil self-impedance)
    - Z_reflected = k^2 * Z_workpiece    (reflected from workpiece)

    where:
    - R_coil: Coil AC resistance (including skin effect)
    - L_coil: Coil self-inductance
    - k: Coupling coefficient between coil and workpiece
    - Z_workpiece: Effective workpiece impedance from ESIM
    """

    def __init__(self, coil, workpiece, frequency):
        """
        Initialize the coupled solver.

        Parameters:
            coil: InductionHeatingCoil object
            workpiece: ESIMWorkpiece object
            frequency: Operating frequency [Hz]
        """
        self.coil = coil
        self.workpiece = workpiece
        self.frequency = frequency
        self.omega = 2 * np.pi * frequency

        # Set frequency on coil
        self.coil.set_frequency(frequency)

        # Solver state
        self.converged = False
        self.iterations = 0
        self.residual_history = []

        # Impedance calculation results
        self.Z_coil_self = None       # Coil self-impedance [Ohm]
        self.Z_reflected = None       # Reflected impedance from workpiece [Ohm]
        self.Z_total = None           # Total impedance [Ohm]
        self.coupling_factor = None   # Coupling coefficient k
        self.L_coil = None            # Coil self-inductance [H]
        self.R_coil = None            # Coil AC resistance [Ohm]
        self.M_mutual = None          # Mutual inductance [H]

    def compute_coil_field_on_workpiece(self):
        """
        Compute the B field from coil at all workpiece panel centers.

        Returns:
            B_fields: Dict {panel_id: complex B vector}
            H_tangential: Dict {panel_id: complex H_t magnitude}
        """
        B_fields = {}
        H_tangential = {}

        for panel in self.workpiece.panels:
            center = panel.center
            normal = panel.normal

            # Compute B at panel center
            B = self.coil.compute_field_at_point(center.tolist())
            B_fields[panel.panel_id] = B

            # Extract tangential H
            H = B / mu_0
            H_n = np.dot(H, normal) * normal
            H_t = H - H_n
            H_t_mag = np.linalg.norm(H_t)

            H_tangential[panel.panel_id] = H_t_mag

        return B_fields, H_tangential

    def solve(self, tol=1e-4, max_iter=50, relaxation=0.5, verbose=True):
        """
        Solve the coupled induction heating problem with fixed-point iteration.

        Parameters:
            tol: Convergence tolerance (relative change in Z)
            max_iter: Maximum number of iterations
            relaxation: Under-relaxation parameter (0 < alpha <= 1)
            verbose: Print iteration progress

        Returns:
            result: Dict with solution data
        """
        if verbose:
            print(f"ESIM Coupled Solver")
            print(f"  Frequency: {self.frequency/1000:.1f} kHz")
            print(f"  Workpiece panels: {self.workpiece.num_panels}")
            print(f"  Tolerance: {tol}")
            print()

        # Initialize: compute field from coil
        B_fields, H_tangential = self.compute_coil_field_on_workpiece()

        # Set initial tangential field on workpiece
        for panel_id, H_t in H_tangential.items():
            self.workpiece.set_tangential_field(panel_id, H_t)

        # Store previous Z values for convergence check
        Z_prev = {p.panel_id: p.Z_surface for p in self.workpiece.panels}

        self.residual_history = []

        for iteration in range(max_iter):
            # Update impedances based on current H field
            self.workpiece.update_all_impedances()

            # Get new Z values
            Z_new = {p.panel_id: p.Z_surface for p in self.workpiece.panels}

            # Check convergence (relative change in Z)
            max_rel_change = 0.0
            for panel_id in Z_new:
                if abs(Z_prev[panel_id]) > 1e-20:
                    rel_change = abs(Z_new[panel_id] - Z_prev[panel_id]) / abs(Z_prev[panel_id])
                    max_rel_change = max(max_rel_change, rel_change)

            self.residual_history.append(max_rel_change)

            if verbose:
                P_total, Q_total = self.workpiece.compute_power_losses()
                print(f"  Iter {iteration+1:3d}: max_rel_change = {max_rel_change:.2e}, "
                      f"P = {P_total:.1f} W, Q = {Q_total:.1f} var")

            if max_rel_change < tol:
                self.converged = True
                self.iterations = iteration + 1
                break

            # Under-relaxation
            for panel_id in Z_new:
                Z_relaxed = (1 - relaxation) * Z_prev[panel_id] + relaxation * Z_new[panel_id]
                # Apply relaxed Z to panel
                self.workpiece.panels[panel_id].Z_surface = Z_relaxed

            Z_prev = {p.panel_id: p.Z_surface for p in self.workpiece.panels}
        else:
            self.converged = False
            self.iterations = max_iter

        # Final power computation
        P_total, Q_total = self.workpiece.compute_power_losses()

        # Compute impedances
        self.compute_coil_self_impedance()
        self.compute_reflected_impedance(P_total, Q_total)
        self.compute_total_impedance()
        impedance_summary = self.get_impedance_summary()

        # Get summary
        summary = self.workpiece.get_summary()

        result = {
            'converged': self.converged,
            'iterations': self.iterations,
            'P_total': P_total,
            'Q_total': Q_total,
            'S_total': np.sqrt(P_total**2 + Q_total**2),
            'power_factor': P_total / np.sqrt(P_total**2 + Q_total**2) if P_total > 0 else 0,
            'max_P_density': summary['max_P_density'],
            'residual_history': self.residual_history,
            'H_tangential': H_tangential,
            'B_fields': B_fields,
            # Impedance results
            'impedance': impedance_summary,
            'Z_coil_self': self.Z_coil_self,
            'Z_reflected': self.Z_reflected,
            'Z_total': self.Z_total,
        }

        if verbose:
            print()
            print(f"Solution {'converged' if self.converged else 'did NOT converge'} "
                  f"in {self.iterations} iterations")
            print(f"  Total power: P = {P_total:.1f} W, Q = {Q_total:.1f} var")
            print(f"  Power factor: {result['power_factor']:.3f}")
            print(f"  Max power density: {summary['max_P_density']/1e3:.2f} kW/m^2")
            print()
            print("Impedance Analysis:")
            print(f"  Coil: L = {impedance_summary['L_coil_uH']:.3f} uH, "
                  f"R = {impedance_summary['R_coil_mOhm']:.3f} mOhm")
            print(f"  Z_coil_self = {self.Z_coil_self.real*1e3:.3f} + j{self.Z_coil_self.imag*1e3:.3f} mOhm")
            print(f"  Z_reflected = {self.Z_reflected.real*1e3:.3f} + j{self.Z_reflected.imag*1e3:.3f} mOhm")
            print(f"  Z_total     = {self.Z_total.real*1e3:.3f} + j{self.Z_total.imag*1e3:.3f} mOhm")
            print(f"  |Z_total|   = {abs(self.Z_total)*1e3:.3f} mOhm, phase = {np.angle(self.Z_total, deg=True):.1f} deg")
            print(f"  Efficiency (P_wp/P_total): {impedance_summary['efficiency']*100:.1f}%")

        return result

    def get_power_distribution(self):
        """
        Get the power distribution over the workpiece surface.

        Returns:
            power_data: List of panel power data
        """
        return self.workpiece.get_power_distribution()

    def get_field_distribution(self):
        """
        Get the field distribution over the workpiece surface.

        Returns:
            field_data: List of dicts with panel field information
        """
        field_data = []
        for panel in self.workpiece.panels:
            field_data.append({
                'panel_id': panel.panel_id,
                'center': panel.center.tolist(),
                'normal': panel.normal.tolist(),
                'H_tangential': float(abs(panel.H_tangential)),
                'Z_surface': complex(panel.Z_surface),
            })
        return field_data

    def compute_coil_self_impedance(self):
        """
        Compute the coil self-impedance Z_coil = R_coil + j*omega*L_coil.

        For spiral coil:
            L = mu_0 * N^2 * R_avg^2 / (2 * R_avg)  (simplified formula)
            L = mu_0 * N^2 * R_avg * (ln(8*R_avg/a) - 2)  (more accurate)

        For loop coil:
            L = mu_0 * R * (ln(8*R/a) - 2)  (Neumann formula for thin ring)

        where:
            N = number of turns
            R_avg = average radius
            a = effective wire radius

        Returns:
            Z_coil: Complex coil self-impedance [Ohm]
        """
        if self.coil.coil_type == 'loop':
            R = self.coil.radius
            wire_w = self.coil.params.get('wire_width', 0.003)
            wire_h = self.coil.params.get('wire_height', wire_w)

            # Effective wire radius (geometric mean for rectangular)
            a_eff = np.sqrt(wire_w * wire_h) / 2

            # Neumann formula for self-inductance of circular loop
            if R > a_eff:
                L = mu_0 * R * (np.log(8 * R / a_eff) - 2)
            else:
                L = mu_0 * R  # Fallback for very thick wire

            N = 1

        elif self.coil.coil_type == 'spiral':
            N = self.coil.num_turns
            R_inner = self.coil.inner_radius
            R_outer = self.coil.outer_radius
            R_avg = (R_inner + R_outer) / 2
            wire_w = self.coil.params.get('wire_width', 0.003)
            wire_h = self.coil.params.get('wire_height', wire_w)

            # Effective wire radius
            a_eff = np.sqrt(wire_w * wire_h) / 2

            # Wheeler's formula for planar spiral (approximate)
            # L = mu_0 * N^2 * R_avg / 2 * (ln(8*R_avg/a) - 2)
            # More accurate: use modified Wheeler formula
            c = (R_outer - R_inner) / 2  # Coil width
            rho = (R_outer - R_inner) / (R_outer + R_inner)  # Fill factor

            if rho < 0.9:
                # Wheeler's formula for flat spiral
                # L = 31.33 * mu_0 * N^2 * R_avg^2 / (8*R_avg + 11*c)  [in SI]
                # Simplified: L = K * mu_0 * N^2 * R_avg
                K = 1.0  # Geometry factor
                L = mu_0 * N**2 * R_avg * K * (np.log(8 * R_avg / a_eff) - 2)
            else:
                # Thin solenoid approximation
                L = mu_0 * N**2 * np.pi * R_avg**2 / (N * self.coil.pitch)

        else:
            # Unknown type, use simple estimate
            R_avg = 0.05
            N = 1
            L = mu_0 * N**2 * R_avg

        # Store inductance
        self.L_coil = L

        # Compute AC resistance (includes skin effect)
        sigma_coil = self.coil.params.get('conductivity', 5.8e7)
        wire_w = self.coil.params.get('wire_width', 0.003)
        wire_h = self.coil.params.get('wire_height', wire_w)

        # Skin depth in copper
        delta = np.sqrt(2 / (self.omega * mu_0 * sigma_coil))

        # Wire length
        if self.coil.coil_type == 'loop':
            wire_length = 2 * np.pi * self.coil.radius
        elif self.coil.coil_type == 'spiral':
            # Approximate spiral length
            R_avg = (self.coil.inner_radius + self.coil.outer_radius) / 2
            wire_length = 2 * np.pi * R_avg * self.coil.num_turns
        else:
            wire_length = 1.0

        # DC resistance
        A_wire = wire_w * wire_h
        R_dc = wire_length / (sigma_coil * A_wire)

        # AC resistance factor (skin effect)
        # For rectangular conductor, Rac/Rdc ~ (w + h) / (4 * delta) for delta << w, h
        if delta < min(wire_w, wire_h) / 2:
            # High frequency: current flows in skin layer
            perimeter = 2 * (wire_w + wire_h)
            A_eff = perimeter * delta  # Effective area for skin current
            R_ac = wire_length / (sigma_coil * A_eff)
        else:
            # Low frequency: uniform current
            R_ac = R_dc

        self.R_coil = R_ac

        # Total self-impedance
        self.Z_coil_self = R_ac + 1j * self.omega * L

        return self.Z_coil_self

    def compute_reflected_impedance(self, P_workpiece, Q_workpiece):
        """
        Compute the reflected impedance from workpiece to coil.

        The reflected impedance represents the loading effect of the workpiece
        on the coil. It is computed from the power balance:

            Z_reflected = (P + jQ) / I^2

        where P and Q are the real and reactive power absorbed by the workpiece.

        Alternative formulation using coupling coefficient:
            Z_reflected = omega^2 * M^2 / Z_workpiece
                        = k^2 * omega * L_coil * omega * L_workpiece / Z_workpiece

        Parameters:
            P_workpiece: Real power absorbed by workpiece [W]
            Q_workpiece: Reactive power (inductive) [var]

        Returns:
            Z_reflected: Complex reflected impedance [Ohm]
        """
        I = self.coil.current

        if abs(I) < 1e-20:
            self.Z_reflected = 0j
            return self.Z_reflected

        # Impedance from power balance
        # P = Re(Z) * I^2, Q = Im(Z) * I^2
        R_reflected = P_workpiece / (I**2)
        X_reflected = Q_workpiece / (I**2)

        self.Z_reflected = R_reflected + 1j * X_reflected

        # Estimate mutual inductance and coupling factor
        if self.L_coil is not None and self.L_coil > 0:
            # From Q_reflected = omega * M^2 / L_workpiece_eff
            # Approximate L_workpiece_eff from workpiece geometry
            # For a slab: L_eff ~ mu_0 * A / delta where A is area, delta is skin depth

            # Total workpiece area
            A_workpiece = sum(p.area for p in self.workpiece.panels)

            # Average skin depth (from first panel's Z_surface)
            if self.workpiece.panels:
                Z_avg = np.mean([abs(p.Z_surface) for p in self.workpiece.panels])
                sigma = self.workpiece.esi_table.sigma
                # Z_s = (1+j) * rho / delta = (1+j) / (sigma * delta)
                # |Z_s| = sqrt(2) / (sigma * delta)
                # delta = sqrt(2) / (sigma * |Z_s|)
                if Z_avg > 1e-20:
                    delta_est = np.sqrt(2) / (sigma * Z_avg)
                else:
                    delta_est = 0.001  # Default 1mm
            else:
                delta_est = 0.001

            # Effective workpiece inductance (very rough estimate)
            # L_workpiece_eff ~ mu_0 * A / (pi * delta) for induced currents
            L_workpiece_eff = mu_0 * A_workpiece / (np.pi * delta_est)

            # Mutual inductance from reflected reactance
            # X_reflected = omega * M^2 / L_workpiece_eff
            if abs(X_reflected) > 1e-20 and L_workpiece_eff > 1e-20:
                M_squared = X_reflected * L_workpiece_eff / self.omega
                if M_squared > 0:
                    self.M_mutual = np.sqrt(M_squared)
                    self.coupling_factor = self.M_mutual / np.sqrt(self.L_coil * L_workpiece_eff)
                else:
                    self.M_mutual = 0
                    self.coupling_factor = 0
            else:
                self.M_mutual = 0
                self.coupling_factor = 0

        return self.Z_reflected

    def compute_total_impedance(self):
        """
        Compute total system impedance seen from the power source.

        Z_total = Z_coil + Z_reflected
                = (R_coil + R_reflected) + j*(omega*L_coil + X_reflected)

        The real part represents total power dissipation (coil + workpiece).
        The imaginary part represents total reactive power (coil inductance + workpiece).

        Returns:
            Z_total: Complex total impedance [Ohm]
        """
        if self.Z_coil_self is None:
            self.compute_coil_self_impedance()

        if self.Z_reflected is None:
            # Use current power values
            P_total, Q_total = self.workpiece.compute_power_losses()
            self.compute_reflected_impedance(P_total, Q_total)

        self.Z_total = self.Z_coil_self + self.Z_reflected

        return self.Z_total

    def get_impedance_summary(self):
        """
        Get a summary of all impedance values.

        Returns:
            summary: Dict with impedance information
        """
        # Ensure impedances are computed
        if self.Z_total is None:
            self.compute_total_impedance()

        summary = {
            # Coil self-impedance
            'L_coil_uH': self.L_coil * 1e6 if self.L_coil else 0,
            'R_coil_mOhm': self.R_coil * 1e3 if self.R_coil else 0,
            'Z_coil_self': self.Z_coil_self,

            # Reflected impedance
            'Z_reflected': self.Z_reflected,
            'R_reflected_mOhm': self.Z_reflected.real * 1e3 if self.Z_reflected else 0,
            'X_reflected_mOhm': self.Z_reflected.imag * 1e3 if self.Z_reflected else 0,

            # Total impedance
            'Z_total': self.Z_total,
            'R_total_mOhm': self.Z_total.real * 1e3 if self.Z_total else 0,
            'X_total_mOhm': self.Z_total.imag * 1e3 if self.Z_total else 0,
            'Z_total_magnitude_mOhm': abs(self.Z_total) * 1e3 if self.Z_total else 0,
            'phase_deg': np.angle(self.Z_total, deg=True) if self.Z_total else 0,

            # Coupling
            'M_mutual_uH': self.M_mutual * 1e6 if self.M_mutual else 0,
            'coupling_factor': self.coupling_factor if self.coupling_factor else 0,

            # Efficiency estimate (P_workpiece / P_total)
            'efficiency': (self.Z_reflected.real / self.Z_total.real
                          if self.Z_total and self.Z_total.real > 0 else 0),
        }

        return summary


def solve_induction_heating(coil_params, workpiece_params, frequency,
                            bh_curve, sigma, tol=1e-4, max_iter=50, verbose=True):
    """
    High-level function to solve induction heating problem.

    Parameters:
        coil_params: Dict with coil parameters (see InductionHeatingCoil)
        workpiece_params: Dict with workpiece parameters (geometry, panels)
        frequency: Operating frequency [Hz]
        bh_curve: BH curve data [[H, B], ...]
        sigma: Workpiece conductivity [S/m]
        tol: Convergence tolerance
        max_iter: Maximum iterations
        verbose: Print progress

    Returns:
        result: Solution result dict
    """
    # Create coil
    coil_type = coil_params.pop('type', 'spiral')
    coil = InductionHeatingCoil(coil_type=coil_type, **coil_params)

    # Create workpiece
    geometry = workpiece_params.get('geometry', 'block')

    if geometry == 'block':
        workpiece = create_esim_block(
            center=workpiece_params.get('center', [0, 0, -0.01]),
            dimensions=workpiece_params.get('dimensions', [0.1, 0.1, 0.02]),
            bh_curve=bh_curve,
            sigma=sigma,
            frequency=frequency,
            panels_per_side=workpiece_params.get('panels_per_side', 5)
        )
    elif geometry == 'cylinder':
        workpiece = create_esim_cylinder(
            center=workpiece_params.get('center', [0, 0, 0]),
            radius=workpiece_params.get('radius', 0.05),
            height=workpiece_params.get('height', 0.02),
            bh_curve=bh_curve,
            sigma=sigma,
            frequency=frequency,
            panels_radial=workpiece_params.get('panels_radial', 8),
            panels_axial=workpiece_params.get('panels_axial', 5)
        )
    else:
        raise ValueError(f"Unknown workpiece geometry: {geometry}")

    # Create and run solver
    solver = ESIMCoupledSolver(coil, workpiece, frequency)
    result = solver.solve(tol=tol, max_iter=max_iter, verbose=verbose)

    # Add additional data
    result['power_distribution'] = solver.get_power_distribution()
    result['field_distribution'] = solver.get_field_distribution()

    return result


# Example usage and test
if __name__ == "__main__":
    print("ESIM Coupled Solver Test")
    print("=" * 60)

    # Steel BH curve
    bh_curve_steel = [
        [0, 0],
        [100, 0.2],
        [250, 0.5],
        [500, 0.9],
        [1000, 1.3],
        [2500, 1.6],
        [5000, 1.8],
        [10000, 1.95],
        [50000, 2.1],
    ]

    sigma_steel = 2e6  # S/m (hot steel)
    freq = 50000  # 50 kHz

    # Create coil
    coil = InductionHeatingCoil(
        coil_type='spiral',
        center=[0, 0, 0.02],  # 20mm above workpiece
        inner_radius=0.03,
        outer_radius=0.06,
        pitch=0.005,
        num_turns=3,
        axis=[0, 0, 1],
        wire_width=0.003,
        wire_height=0.002,
        cross_section='r',
        conductivity=5.8e7
    )
    coil.set_current(100)  # 100 A

    # Create workpiece
    workpiece = create_esim_block(
        center=[0, 0, -0.01],  # 10mm below origin
        dimensions=[0.1, 0.1, 0.02],  # 100mm x 100mm x 20mm
        bh_curve=bh_curve_steel,
        sigma=sigma_steel,
        frequency=freq,
        panels_per_side=6
    )

    print(f"Coil: 3-turn spiral, R_in=30mm, R_out=60mm, I=100A")
    print(f"Workpiece: 100mm x 100mm x 20mm steel block")
    print(f"Frequency: {freq/1000} kHz")
    print(f"Conductivity: {sigma_steel/1e6} MS/m")
    print()

    # Create and run solver
    solver = ESIMCoupledSolver(coil, workpiece, freq)
    result = solver.solve(tol=1e-4, max_iter=20, verbose=True)

    print()
    print("Power Distribution on Top Face:")
    print("-" * 50)
    power_data = solver.get_power_distribution()

    # Show top face panels (first n^2 panels where n=panels_per_side)
    n_top = 36  # 6x6 for panels_per_side=6
    top_panels = power_data[:n_top]

    print(f"{'Panel':>6} {'P_loss [W]':>12} {'P_density [kW/m^2]':>18}")
    for pd in top_panels[:9]:  # Show first 9
        print(f"{pd['panel_id']:>6} {pd['P_loss']:>12.3f} {pd['P_density']/1e3:>18.2f}")
    print("...")

    print()
    print("Test completed!")
