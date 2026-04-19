"""
ESIM Coupled Solver for Induction Heating and WPT Analysis

This module implements coupled solvers that combine:
1. FastImp coil (conductor model with eddy currents)
2. ESIM workpiece (nonlinear ferromagnetic material)
3. WPT (Wireless Power Transfer) multi-coil coupling analysis

Key Features:
- Fixed-point iteration for nonlinear material (ESIM)
- Neumann integral for mutual inductance calculation
- Coupling coefficient k = M / sqrt(L1*L2)
- Two-coil WPT impedance matrix analysis

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
        """Legacy stub — CndSpiral/CndLoop APIs have been removed.

        Coil field is computed via the analytical model (_create_analytical_model).
        This method is never called (self._rad is always None).
        """
        pass

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
            # Use unified rad.Fld() for field computation
            # For conductor objects, Fld() returns 6 values:
            # [Bx_re, By_re, Bz_re, Bx_im, By_im, Bz_im]
            B_result = self._rad.Fld(self.handle, 'b', point)
            if hasattr(B_result, '__len__') and len(B_result) == 6:
                return np.array([
                    B_result[0] + 1j * B_result[3],
                    B_result[1] + 1j * B_result[4],
                    B_result[2] + 1j * B_result[5]
                ])
            elif hasattr(B_result, '__len__') and len(B_result) == 3:
                # Static field (real only)
                return np.array(B_result[:3], dtype=complex)
            else:
                return np.array([0, 0, 0], dtype=complex)
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
        """Compute B field from a circular current loop."""
        from radia.biot_savart import h_filament, MU0
        center = np.asarray(center); normal = np.asarray(normal)
        normal = normal / np.linalg.norm(normal)
        if abs(normal[2]) < 0.9:
            t1 = np.cross(normal, [0, 0, 1])
        else:
            t1 = np.cross(normal, [1, 0, 0])
        t1 = t1 / np.linalg.norm(t1)
        t2 = np.cross(normal, t1)
        n_seg = 100
        H = np.zeros(3, dtype=complex)
        for i in range(n_seg):
            a1 = 2 * np.pi * i / n_seg
            a2 = 2 * np.pi * (i + 1) / n_seg
            p1 = center + radius * (np.cos(a1) * t1 + np.sin(a1) * t2)
            p2 = center + radius * (np.cos(a2) * t1 + np.sin(a2) * t2)
            H += h_filament(p1, p2, point, current=float(np.real(current)))
        return (MU0 * H).astype(complex)

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

    def get_wire_segments(self, n_segments=100):
        """
        Get discretized wire segments for Neumann integral calculation.

        Parameters:
            n_segments: Number of segments to discretize each turn

        Returns:
            segments: List of [start_point, end_point] for each segment
            dl_vectors: List of dl vectors (end - start)
            midpoints: List of segment midpoints
        """
        segments = []
        dl_vectors = []
        midpoints = []

        if self.coil_type == 'loop':
            # Single circular loop
            R = self.radius
            center = self.center
            normal = self.normal

            # Create local coordinate system
            if abs(normal[2]) < 0.9:
                u = np.cross(normal, [0, 0, 1])
            else:
                u = np.cross(normal, [1, 0, 0])
            u = u / np.linalg.norm(u)
            v = np.cross(normal, u)

            # Discretize loop
            dphi = 2 * np.pi / n_segments
            for i in range(n_segments):
                phi1 = i * dphi
                phi2 = (i + 1) * dphi

                p1 = center + R * (np.cos(phi1) * u + np.sin(phi1) * v)
                p2 = center + R * (np.cos(phi2) * u + np.sin(phi2) * v)

                segments.append([p1, p2])
                dl_vectors.append(p2 - p1)
                midpoints.append((p1 + p2) / 2)

        elif self.coil_type == 'spiral':
            # Multi-turn spiral coil
            N = self.num_turns
            R_inner = self.inner_radius
            R_outer = self.outer_radius
            pitch = self.pitch
            center = self.center
            axis = self.axis

            # Create local coordinate system
            if abs(axis[2]) < 0.9:
                u = np.cross(axis, [0, 0, 1])
            else:
                u = np.cross(axis, [1, 0, 0])
            u = u / np.linalg.norm(u)
            v = np.cross(axis, u)

            # Total angle for all turns
            total_angle = 2 * np.pi * N
            n_total_segments = n_segments * N
            dphi = total_angle / n_total_segments

            for i in range(n_total_segments):
                # Parameter along spiral (0 to N turns)
                t1 = i / n_total_segments
                t2 = (i + 1) / n_total_segments

                phi1 = t1 * total_angle
                phi2 = t2 * total_angle

                # Radius varies linearly from inner to outer
                R1 = R_inner + t1 * (R_outer - R_inner)
                R2 = R_inner + t2 * (R_outer - R_inner)

                # Height along axis
                z1 = t1 * N * pitch
                z2 = t2 * N * pitch

                # 3D positions
                p1 = center + R1 * (np.cos(phi1) * u + np.sin(phi1) * v) + z1 * axis
                p2 = center + R2 * (np.cos(phi2) * u + np.sin(phi2) * v) + z2 * axis

                segments.append([p1, p2])
                dl_vectors.append(p2 - p1)
                midpoints.append((p1 + p2) / 2)

        return segments, dl_vectors, midpoints

    def compute_self_inductance_neumann(self, n_segments=100):
        """
        Compute self-inductance using Neumann integral with GMD correction.

        For self-inductance, the Neumann integral diverges when segments overlap.
        We use the Geometric Mean Distance (GMD) correction for the diagonal terms.

        L = (mu_0 / 4*pi) * SUM_i SUM_j (dl_i . dl_j) / R_ij

        where R_ij = |r_i - r_j| for i != j, and R_ii = GMD for i == j.

        Parameters:
            n_segments: Number of segments for discretization

        Returns:
            L: Self-inductance [H]
        """
        segments, dl_vectors, midpoints = self.get_wire_segments(n_segments)
        n = len(segments)

        if n == 0:
            return 0.0

        # Wire dimensions for GMD calculation
        wire_w = self.params.get('wire_width', 0.003)
        wire_h = self.params.get('wire_height', wire_w)

        # GMD for rectangular cross-section (approximate)
        # GMD ~ 0.2235 * (w + h) for square
        # GMD ~ exp(-1/4) * sqrt(w*h) for rectangle
        a_eff = np.exp(-0.25) * np.sqrt(wire_w * wire_h)

        L = 0.0
        for i in range(n):
            for j in range(n):
                dl_i = dl_vectors[i]
                dl_j = dl_vectors[j]
                dot_product = np.dot(dl_i, dl_j)

                if i == j:
                    # Self term: use GMD
                    # For a segment of length l, self-inductance contribution
                    # L_self = (mu_0 / 4*pi) * l * (ln(2*l/a) - 1)
                    l_seg = np.linalg.norm(dl_i)
                    if l_seg > a_eff:
                        L += (mu_0 / (4 * np.pi)) * l_seg * (np.log(2 * l_seg / a_eff) - 1)
                else:
                    # Mutual term: standard Neumann
                    r_ij = np.linalg.norm(midpoints[i] - midpoints[j])
                    if r_ij > 1e-15:
                        L += (mu_0 / (4 * np.pi)) * dot_product / r_ij

        return L


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


def compute_mutual_inductance(coil1, coil2, n_segments=100):
    """
    Compute mutual inductance between two coils using Neumann integral.

    M = (mu_0 / 4*pi) * oint oint (dl_1 . dl_2) / |r_12|

    This is the fundamental formula for mutual inductance calculation
    that accounts for arbitrary coil geometries.

    Parameters:
        coil1: First InductionHeatingCoil object
        coil2: Second InductionHeatingCoil object
        n_segments: Number of segments per turn for discretization

    Returns:
        M: Mutual inductance [H]
    """
    # Get wire segments for both coils
    _, dl1_vectors, midpoints1 = coil1.get_wire_segments(n_segments)
    _, dl2_vectors, midpoints2 = coil2.get_wire_segments(n_segments)

    n1 = len(midpoints1)
    n2 = len(midpoints2)

    if n1 == 0 or n2 == 0:
        return 0.0

    # Compute Neumann integral
    M = 0.0
    for i in range(n1):
        for j in range(n2):
            dl1 = dl1_vectors[i]
            dl2 = dl2_vectors[j]
            r12 = np.linalg.norm(midpoints1[i] - midpoints2[j])

            if r12 > 1e-15:
                dot_product = np.dot(dl1, dl2)
                M += dot_product / r12

    M *= mu_0 / (4 * np.pi)

    return M


def compute_coupling_coefficient(coil1, coil2, n_segments=100,
                                  L1=None, L2=None, M=None):
    """
    Compute coupling coefficient k between two coils.

    k = M / sqrt(L1 * L2)

    where:
        M = mutual inductance
        L1, L2 = self-inductances of coil1 and coil2

    Parameters:
        coil1: First InductionHeatingCoil object
        coil2: Second InductionHeatingCoil object
        n_segments: Number of segments for discretization
        L1: Pre-computed self-inductance of coil1 (optional)
        L2: Pre-computed self-inductance of coil2 (optional)
        M: Pre-computed mutual inductance (optional)

    Returns:
        k: Coupling coefficient (0 <= k <= 1)
        L1: Self-inductance of coil1 [H]
        L2: Self-inductance of coil2 [H]
        M: Mutual inductance [H]
    """
    # Compute self-inductances if not provided
    if L1 is None:
        L1 = coil1.compute_self_inductance_neumann(n_segments)

    if L2 is None:
        L2 = coil2.compute_self_inductance_neumann(n_segments)

    # Compute mutual inductance if not provided
    if M is None:
        M = compute_mutual_inductance(coil1, coil2, n_segments)

    # Compute coupling coefficient
    if L1 > 0 and L2 > 0:
        k = abs(M) / np.sqrt(L1 * L2)
    else:
        k = 0.0

    # Clip to physical range [0, 1]
    k = min(k, 1.0)

    return k, L1, L2, M


class WPTCoupledSolver:
    """
    Wireless Power Transfer (WPT) coupled solver for two-coil systems.

    This solver computes the electrical characteristics of a two-coil WPT system:
    - Mutual inductance M via Neumann integral
    - Coupling coefficient k = M / sqrt(L1 * L2)
    - Impedance matrix [Z11, Z12; Z21, Z22]
    - Power transfer efficiency and optimal load

    Impedance Matrix:
        V1 = Z11*I1 + Z12*I2
        V2 = Z21*I1 + Z22*I2

        where:
            Z11 = R1 + j*omega*L1  (primary self-impedance)
            Z22 = R2 + j*omega*L2  (secondary self-impedance)
            Z12 = Z21 = j*omega*M  (mutual impedance)

    Resonant Topologies:
        - S-S: Series-Series (both sides series capacitors)
        - S-P: Series-Primary, Parallel-Secondary
        - P-S: Parallel-Primary, Series-Secondary
        - P-P: Parallel-Parallel
    """

    def __init__(self, coil_primary, coil_secondary, frequency):
        """
        Initialize the WPT coupled solver.

        Parameters:
            coil_primary: Primary (transmitter) InductionHeatingCoil
            coil_secondary: Secondary (receiver) InductionHeatingCoil
            frequency: Operating frequency [Hz]
        """
        self.coil1 = coil_primary
        self.coil2 = coil_secondary
        self.frequency = frequency
        self.omega = 2 * np.pi * frequency

        # Set frequency on coils
        self.coil1.set_frequency(frequency)
        self.coil2.set_frequency(frequency)

        # Inductance values
        self.L1 = None  # Primary self-inductance [H]
        self.L2 = None  # Secondary self-inductance [H]
        self.M = None   # Mutual inductance [H]
        self.k = None   # Coupling coefficient

        # Resistance values
        self.R1 = None  # Primary resistance [Ohm]
        self.R2 = None  # Secondary resistance [Ohm]
        self.Rm = None  # Mutual resistance (proximity effect) [Ohm]

        # Impedance matrix
        self.Z_matrix = None  # 2x2 complex impedance matrix

    def compute_inductances(self, n_segments=100):
        """
        Compute all inductance values using Neumann integral.

        Parameters:
            n_segments: Segments per turn for discretization

        Returns:
            L1, L2, M, k: Inductances [H] and coupling coefficient
        """
        # Self-inductances
        self.L1 = self.coil1.compute_self_inductance_neumann(n_segments)
        self.L2 = self.coil2.compute_self_inductance_neumann(n_segments)

        # Mutual inductance
        self.M = compute_mutual_inductance(self.coil1, self.coil2, n_segments)

        # Coupling coefficient
        if self.L1 > 0 and self.L2 > 0:
            self.k = abs(self.M) / np.sqrt(self.L1 * self.L2)
        else:
            self.k = 0.0

        return self.L1, self.L2, self.M, self.k

    def compute_ac_resistance(self, coil, wire_length=None):
        """
        Compute AC resistance of a coil including skin effect.

        Parameters:
            coil: InductionHeatingCoil object
            wire_length: Wire length [m] (computed if not provided)

        Returns:
            R_ac: AC resistance [Ohm]
        """
        sigma = coil.params.get('conductivity', 5.8e7)
        wire_w = coil.params.get('wire_width', 0.003)
        wire_h = coil.params.get('wire_height', wire_w)

        # Skin depth
        delta = np.sqrt(2 / (self.omega * mu_0 * sigma))

        # Wire length (approximate)
        if wire_length is None:
            if coil.coil_type == 'loop':
                wire_length = 2 * np.pi * coil.radius
            elif coil.coil_type == 'spiral':
                R_avg = (coil.inner_radius + coil.outer_radius) / 2
                wire_length = 2 * np.pi * R_avg * coil.num_turns
            else:
                wire_length = 1.0

        # DC resistance
        A_wire = wire_w * wire_h
        R_dc = wire_length / (sigma * A_wire)

        # AC resistance with skin effect
        if delta < min(wire_w, wire_h) / 2:
            perimeter = 2 * (wire_w + wire_h)
            A_eff = perimeter * delta
            R_ac = wire_length / (sigma * A_eff)
        else:
            R_ac = R_dc

        return R_ac

    def compute_mutual_resistance(self, n_segments=100):
        """
        Compute mutual resistance (Rm) due to proximity effect between two coils.

        The mutual resistance accounts for eddy current losses induced in one coil
        by the alternating magnetic field of the other coil.

        Physical mechanism:
        - The magnetic field from coil 1 penetrates coil 2's conductor
        - This induces eddy currents in coil 2's wire
        - These eddy currents dissipate power (I^2*Rm loss)
        - The effect is reciprocal (Rm12 = Rm21)

        Calculation Method:
        For I1 = 1A in coil 1, the external magnetic field H_ext at coil 2's
        wire location is computed. The power loss in coil 2 due to this external
        field is:

            P_prox = Re(Z_s) * |H_ext|^2 * (wire surface area)

        where Z_s is the surface impedance from ESIM cell problem.

        The mutual resistance is defined such that:
            P_prox = Rm * |I1|^2

        Therefore:
            Rm = Re(Z_s) * |H_ext/I1|^2 * (wire surface area)

        Parameters:
            n_segments: Number of segments for field calculation

        Returns:
            Rm: Mutual resistance [Ohm]
        """
        # Import ESIM cell problem solver
        from .esim_cell_problem import ESIMCellProblemSolver

        # Get coil parameters
        sigma1 = self.coil1.params.get('conductivity', 5.8e7)
        sigma2 = self.coil2.params.get('conductivity', 5.8e7)
        wire_w1 = self.coil1.params.get('wire_width', 0.003)
        wire_h1 = self.coil1.params.get('wire_height', wire_w1)
        wire_w2 = self.coil2.params.get('wire_width', 0.003)
        wire_h2 = self.coil2.params.get('wire_height', wire_w2)

        # Compute mutual inductance if not done
        if self.M is None:
            self.compute_inductances(n_segments)

        # Compute AC resistances if not already done
        if self.R1 is None:
            self.R1 = self.compute_ac_resistance(self.coil1)
        if self.R2 is None:
            self.R2 = self.compute_ac_resistance(self.coil2)

        # Wire lengths
        if self.coil1.coil_type == 'loop':
            L_wire1 = 2 * np.pi * self.coil1.radius
        elif self.coil1.coil_type == 'spiral':
            R_avg1 = (self.coil1.inner_radius + self.coil1.outer_radius) / 2
            L_wire1 = 2 * np.pi * R_avg1 * self.coil1.num_turns
        else:
            L_wire1 = 1.0

        if self.coil2.coil_type == 'loop':
            L_wire2 = 2 * np.pi * self.coil2.radius
        elif self.coil2.coil_type == 'spiral':
            R_avg2 = (self.coil2.inner_radius + self.coil2.outer_radius) / 2
            L_wire2 = 2 * np.pi * R_avg2 * self.coil2.num_turns
        else:
            L_wire2 = 1.0

        # Wire perimeters
        perimeter1 = 2 * (wire_w1 + wire_h1)
        perimeter2 = 2 * (wire_w2 + wire_h2)

        # Wire surface areas
        S_wire1 = L_wire1 * perimeter1
        S_wire2 = L_wire2 * perimeter2

        # ==== Proximity effect calculation ====
        #
        # Physical model:
        # When I1 = 1A flows in coil 1, it creates a magnetic field at coil 2's
        # wire location. This field is approximately:
        #
        #   H_ext = B_ext / mu_0 = M * I1 / (mu_0 * A_coil2)
        #
        # where A_coil2 is the effective area of coil 2.
        #
        # More accurately, H_ext at the wire surface can be estimated from:
        #   - The mutual flux linkage: Phi_12 = M * I1
        #   - This flux passes through coil 2's area
        #   - H at wire ~ B / mu_0 ~ (M * I1) / (mu_0 * A_eff)
        #
        # For WPT coils, a simpler approximation using Ampere's law:
        #   H_ext ~ I1 / (2 * pi * gap)  for loosely coupled coils
        #   H_ext ~ k * I1 * N / (2 * R)  for closely coupled coils
        #
        # We use a combined formula based on the mutual flux:
        #   H_ext = k * sqrt(omega * L1 / mu_0) / sqrt(R_avg)
        #
        # This gives H field per unit current that accounts for coupling.

        Rm = 0.0

        if self.k is not None and self.k > 0 and self.L1 > 0 and self.L2 > 0:
            # Effective coil radii
            if self.coil1.coil_type == 'loop':
                R1_eff = self.coil1.radius
            else:
                R1_eff = (self.coil1.inner_radius + self.coil1.outer_radius) / 2

            if self.coil2.coil_type == 'loop':
                R2_eff = self.coil2.radius
            else:
                R2_eff = (self.coil2.inner_radius + self.coil2.outer_radius) / 2

            # Number of turns
            N1 = getattr(self.coil1, 'num_turns', 1)
            N2 = getattr(self.coil2, 'num_turns', 1)

            # Gap between coils (z-direction distance)
            z1 = self.coil1.center[2] if hasattr(self.coil1, 'center') else 0
            z2 = self.coil2.center[2] if hasattr(self.coil2, 'center') else 0
            gap = abs(z2 - z1)
            if gap < 0.001:  # Minimum gap 1mm
                gap = 0.001

            # ==== H_ext at coil 2 due to I1 = 1A in coil 1 ====
            # Using Biot-Savart approximation for multi-turn coil:
            #   B_center ~ mu_0 * N * I / (2 * R)  (at center of loop)
            # At distance z on axis:
            #   B_z ~ mu_0 * N * I * R^2 / (2 * (R^2 + z^2)^1.5)
            #
            # For proximity effect, we need H at the wire surface, which is
            # roughly the field at the gap location.

            # Field from coil 1 at coil 2's location (I1 = 1A)
            # Using simplified on-axis formula
            B_at_coil2 = mu_0 * N1 * R1_eff**2 / (2 * (R1_eff**2 + gap**2)**1.5)
            H_ext2 = B_at_coil2 / mu_0  # H = B/mu_0 in air

            # Field from coil 2 at coil 1's location (I2 = 1A)
            B_at_coil1 = mu_0 * N2 * R2_eff**2 / (2 * (R2_eff**2 + gap**2)**1.5)
            H_ext1 = B_at_coil1 / mu_0

            # Clamp to reasonable range
            H_ext1 = max(H_ext1, 0.1)
            H_ext2 = max(H_ext2, 0.1)
            H_ext1 = min(H_ext1, 1e6)
            H_ext2 = min(H_ext2, 1e6)

            # ==== Proximity effect power loss calculation ====
            #
            # Compute the external H field at wire locations using Biot-Savart
            # integration, then calculate eddy current loss.
            #
            # For mutual proximity between separate WPT coils (external field source),
            # we use a hybrid approach:
            #
            # 1. ESIM (surface) model: P = Re(Z_s) * H^2 * S
            #    - Appropriate when field is uniform across wire cross-section
            #    - Surface impedance Z_s = (1+j)/(sigma*delta)
            #
            # 2. Dowell (volume) model: P ~ omega^2 * H^2 * h^3 * L / delta
            #    - Appropriate when field varies across wire cross-section
            #    - More relevant for internal proximity (adjacent turns)
            #
            # For external mutual coupling, the ESIM model is more appropriate,
            # but with an enhancement factor for field gradient effects.
            #
            # Enhancement factor accounts for:
            # - Non-uniform field distribution along wire
            # - Field gradient across wire cross-section
            # - Typically 2-4x for closely coupled WPT coils

            # Skin depths
            delta1 = np.sqrt(2 / (self.omega * mu_0 * sigma1))
            delta2 = np.sqrt(2 / (self.omega * mu_0 * sigma2))

            # ===== Biot-Savart field integration =====
            # Compute actual H field at wire locations using full Biot-Savart
            try:
                _, dl_vec1, mid1 = self.coil1.get_wire_segments(n_segments)
                _, dl_vec2, mid2 = self.coil2.get_wire_segments(n_segments)

                # Compute H at coil2 wire from coil1 (I1 = 1A)
                H2_values = []
                for pt in mid2:
                    B = np.array([0.0, 0.0, 0.0])
                    for i, dl in enumerate(dl_vec1):
                        r_vec = pt - mid1[i]
                        r_mag = np.linalg.norm(r_vec)
                        if r_mag > 1e-10:
                            dB = (mu_0 / (4 * np.pi)) * np.cross(dl, r_vec) / (r_mag**3)
                            B += dB
                    H2_values.append(np.linalg.norm(B) / mu_0)
                H2_values = np.array(H2_values)
                H_ext2 = np.sqrt(np.mean(H2_values**2))  # RMS

                # Compute H at coil1 wire from coil2 (I2 = 1A)
                H1_values = []
                for pt in mid1:
                    B = np.array([0.0, 0.0, 0.0])
                    for i, dl in enumerate(dl_vec2):
                        r_vec = pt - mid2[i]
                        r_mag = np.linalg.norm(r_vec)
                        if r_mag > 1e-10:
                            dB = (mu_0 / (4 * np.pi)) * np.cross(dl, r_vec) / (r_mag**3)
                            B += dB
                    H1_values.append(np.linalg.norm(B) / mu_0)
                H1_values = np.array(H1_values)
                H_ext1 = np.sqrt(np.mean(H1_values**2))  # RMS

                # Field non-uniformity factor: max/avg ratio indicates gradient
                # Higher gradient means more eddy current loss
                grad_factor2 = H2_values.max() / (H2_values.mean() + 1e-10)
                grad_factor1 = H1_values.max() / (H1_values.mean() + 1e-10)

            except Exception:
                # Fallback to simple on-axis formula if Biot-Savart fails
                grad_factor1 = grad_factor2 = 1.5

            # Clamp to reasonable range
            H_ext1 = max(H_ext1, 0.1)
            H_ext2 = max(H_ext2, 0.1)

            # Surface impedance (classical formula)
            Z_s1_real = 1 / (sigma1 * delta1)
            Z_s2_real = 1 / (sigma2 * delta2)

            # Enhancement factor for proximity effect
            # Based on coupling coefficient and field gradient
            # Higher k -> more field concentration -> higher losses
            k_factor = 1 + 2 * self.k  # Ranges from 1 (k=0) to 3 (k=1)
            enhance1 = k_factor * grad_factor1
            enhance2 = k_factor * grad_factor2

            # Power loss in coil 2 due to H_ext from coil 1 (I1 = 1A)
            # P_2 = enhance * Re(Z_s) * H_ext^2 * S_wire
            Rm_from_coil1 = enhance2 * Z_s2_real * (H_ext2**2) * S_wire2

            # Power loss in coil 1 due to H_ext from coil 2 (I2 = 1A)
            Rm_from_coil2 = enhance1 * Z_s1_real * (H_ext1**2) * S_wire1

            # Total mutual resistance (average for symmetry)
            Rm = 0.5 * (Rm_from_coil1 + Rm_from_coil2)

        # Physical upper bound: Rm should not exceed self-resistance
        # For WPT systems with high coupling (k > 0.3), Rm/R can be 10% to 100%
        R_avg = np.sqrt(self.R1 * self.R2)
        Rm = min(Rm, R_avg)

        self.Rm = Rm
        return Rm

    def compute_impedance_matrix(self, include_mutual_resistance=True):
        """
        Compute the 2x2 impedance matrix of the coupled coil system.

        Z = [Z11, Z12]
            [Z21, Z22]

        where:
            Z11 = R1 + j*omega*L1
            Z22 = R2 + j*omega*L2
            Z12 = Z21 = Rm + j*omega*M  (includes mutual resistance from proximity effect)

        The mutual resistance Rm accounts for eddy current losses induced in one
        coil by the magnetic field of the other coil (proximity effect).

        Parameters:
            include_mutual_resistance: If True, include Rm in off-diagonal terms

        Returns:
            Z_matrix: 2x2 complex numpy array
        """
        # Compute inductances if not done
        if self.L1 is None or self.L2 is None or self.M is None:
            self.compute_inductances()

        # Compute resistances
        if self.R1 is None:
            self.R1 = self.compute_ac_resistance(self.coil1)
        if self.R2 is None:
            self.R2 = self.compute_ac_resistance(self.coil2)

        # Compute mutual resistance
        if include_mutual_resistance and self.Rm is None:
            self.compute_mutual_resistance()

        # Build impedance matrix
        Z11 = self.R1 + 1j * self.omega * self.L1
        Z22 = self.R2 + 1j * self.omega * self.L2

        # Off-diagonal: Rm + j*omega*M (mutual impedance with proximity loss)
        if include_mutual_resistance and self.Rm is not None:
            Z12 = self.Rm + 1j * self.omega * self.M
        else:
            Z12 = 1j * self.omega * self.M

        Z21 = Z12  # Reciprocity

        self.Z_matrix = np.array([[Z11, Z12],
                                   [Z21, Z22]], dtype=complex)

        return self.Z_matrix

    def compute_resonant_capacitors(self, topology='SS'):
        """
        Compute resonant capacitors for the given topology.

        Parameters:
            topology: 'SS', 'SP', 'PS', or 'PP'

        Returns:
            C1, C2: Capacitances [F] for primary and secondary
        """
        if self.L1 is None or self.L2 is None:
            self.compute_inductances()

        omega = self.omega

        if topology == 'SS':
            # Series-Series: C = 1 / (omega^2 * L)
            C1 = 1 / (omega**2 * self.L1)
            C2 = 1 / (omega**2 * self.L2)
        elif topology == 'SP':
            # Series-Primary, Parallel-Secondary
            C1 = 1 / (omega**2 * self.L1)
            # For parallel: C2 = 1/(omega^2 * L2) * (1 - k^2)
            C2 = 1 / (omega**2 * self.L2 * (1 - self.k**2)) if self.k < 1 else np.inf
        elif topology == 'PS':
            # Parallel-Primary, Series-Secondary
            C1 = 1 / (omega**2 * self.L1 * (1 - self.k**2)) if self.k < 1 else np.inf
            C2 = 1 / (omega**2 * self.L2)
        elif topology == 'PP':
            # Parallel-Parallel
            C1 = 1 / (omega**2 * self.L1 * (1 - self.k**2)) if self.k < 1 else np.inf
            C2 = 1 / (omega**2 * self.L2 * (1 - self.k**2)) if self.k < 1 else np.inf
        else:
            raise ValueError(f"Unknown topology: {topology}")

        return C1, C2

    def compute_transfer_efficiency(self, R_load):
        """
        Compute power transfer efficiency for a given load resistance.

        For S-S topology at resonance:
            eta = k^2 * Q1 * Q2 / (1 + k^2 * Q1 * Q2) * R_load / (R2 + R_load)

        Parameters:
            R_load: Load resistance [Ohm]

        Returns:
            eta: Transfer efficiency (0 to 1)
            P_load: Power delivered to load [W] (for I1 = 1A)
            P_total: Total input power [W]
        """
        if self.L1 is None or self.R1 is None:
            self.compute_impedance_matrix()

        # Quality factors
        Q1 = self.omega * self.L1 / self.R1
        Q2 = self.omega * self.L2 / self.R2

        # For S-S resonant at omega_0:
        # Input impedance: Z_in = R1 + (omega*M)^2 / (R2 + R_load)
        # Reflected resistance: R_ref = (omega*M)^2 / (R2 + R_load)

        omega_M = self.omega * self.M
        R_ref = (omega_M)**2 / (self.R2 + R_load)
        Z_in = self.R1 + R_ref

        # For primary current I1 = 1A:
        P_in = abs(Z_in)  # Input power (real part * I1^2)

        # Secondary current at resonance:
        I2_mag = omega_M / (self.R2 + R_load)  # |I2/I1|

        # Power to load
        P_load = I2_mag**2 * R_load

        # Efficiency
        eta = P_load / (self.R1 + R_ref) if (self.R1 + R_ref) > 0 else 0

        return eta, P_load, P_in

    def compute_optimal_load(self):
        """
        Compute the optimal load resistance for maximum efficiency.

        For S-S topology:
            R_load_opt = R2 * sqrt(1 + k^2 * Q1 * Q2)

        Returns:
            R_load_opt: Optimal load resistance [Ohm]
            eta_max: Maximum efficiency
        """
        if self.L1 is None or self.R1 is None:
            self.compute_impedance_matrix()

        Q1 = self.omega * self.L1 / self.R1
        Q2 = self.omega * self.L2 / self.R2

        # Optimal load for maximum efficiency
        R_load_opt = self.R2 * np.sqrt(1 + self.k**2 * Q1 * Q2)

        # Maximum efficiency
        kQ = self.k * np.sqrt(Q1 * Q2)
        eta_max = kQ**2 / (1 + np.sqrt(1 + kQ**2))**2

        return R_load_opt, eta_max

    def analyze(self, n_segments=100, R_load=None, topology='SS', verbose=True):
        """
        Perform complete WPT system analysis.

        Parameters:
            n_segments: Segments for inductance calculation
            R_load: Load resistance [Ohm] (uses optimal if None)
            topology: Resonant topology ('SS', 'SP', 'PS', 'PP')
            verbose: Print results

        Returns:
            result: Dict with all analysis results
        """
        # Compute inductances
        self.compute_inductances(n_segments)

        # Compute impedance matrix
        self.compute_impedance_matrix()

        # Compute resonant capacitors
        C1, C2 = self.compute_resonant_capacitors(topology)

        # Compute optimal load
        R_load_opt, eta_max = self.compute_optimal_load()

        # Use specified or optimal load
        if R_load is None:
            R_load = R_load_opt

        # Compute efficiency at operating point
        eta, P_load, P_in = self.compute_transfer_efficiency(R_load)

        # Quality factors
        Q1 = self.omega * self.L1 / self.R1
        Q2 = self.omega * self.L2 / self.R2

        result = {
            # Inductances
            'L1_uH': self.L1 * 1e6,
            'L2_uH': self.L2 * 1e6,
            'M_uH': self.M * 1e6,
            'k': self.k,

            # Resistances
            'R1_mOhm': self.R1 * 1e3,
            'R2_mOhm': self.R2 * 1e3,
            'Rm_mOhm': self.Rm * 1e3 if self.Rm is not None else 0.0,

            # Quality factors
            'Q1': Q1,
            'Q2': Q2,

            # Resonant capacitors
            'C1_nF': C1 * 1e9,
            'C2_nF': C2 * 1e9,

            # Efficiency
            'R_load_opt_Ohm': R_load_opt,
            'R_load_Ohm': R_load,
            'eta_max': eta_max,
            'eta': eta,

            # Power (for I1 = 1A)
            'P_in_W': P_in,
            'P_load_W': P_load,

            # Impedance matrix
            'Z_matrix': self.Z_matrix,
        }

        if verbose:
            print()
            print("=" * 60)
            print("WPT System Analysis Results")
            print("=" * 60)
            print()
            print(f"Frequency: {self.frequency/1000:.1f} kHz")
            print(f"Topology: {topology}")
            print()
            print("Inductance Parameters:")
            print(f"  L1 (primary)  = {result['L1_uH']:.3f} uH")
            print(f"  L2 (secondary)= {result['L2_uH']:.3f} uH")
            print(f"  M (mutual)    = {result['M_uH']:.3f} uH")
            print(f"  k (coupling)  = {result['k']:.4f}")
            print()
            print("Resistance Parameters:")
            print(f"  R1 (primary)  = {result['R1_mOhm']:.4f} mOhm")
            print(f"  R2 (secondary)= {result['R2_mOhm']:.4f} mOhm")
            print(f"  Rm (mutual)   = {result['Rm_mOhm']:.4f} mOhm  (proximity effect)")
            print()
            print("Quality Factors:")
            print(f"  Q1 = {result['Q1']:.1f}")
            print(f"  Q2 = {result['Q2']:.1f}")
            print()
            print("Resonant Capacitors:")
            print(f"  C1 = {result['C1_nF']:.2f} nF")
            print(f"  C2 = {result['C2_nF']:.2f} nF")
            print()
            print("Efficiency Analysis:")
            print(f"  Optimal load:  R_load = {result['R_load_opt_Ohm']:.2f} Ohm")
            print(f"  Max efficiency: eta_max = {result['eta_max']*100:.1f}%")
            print(f"  At R_load = {result['R_load_Ohm']:.2f} Ohm:")
            print(f"    Efficiency: eta = {result['eta']*100:.1f}%")
            print(f"    Power (I1=1A): P_load = {result['P_load_W']:.2f} W")
            print()
            print("Impedance Matrix (Z):")
            print(f"  Z11 = {self.Z_matrix[0,0].real:.4f} + j{self.Z_matrix[0,0].imag:.4f} Ohm")
            print(f"  Z12 = {self.Z_matrix[0,1].real:.4f} + j{self.Z_matrix[0,1].imag:.4f} Ohm  (Rm + j*omega*M)")
            print(f"  Z22 = {self.Z_matrix[1,1].real:.4f} + j{self.Z_matrix[1,1].imag:.4f} Ohm")
            print()
            print("Note: Z12 real part (Rm) represents mutual resistance from proximity effect.")
            print()

        return result


def analyze_coil_coupling(coil1, coil2, frequency, n_segments=100, verbose=True):
    """
    Convenience function to analyze coupling between two coils.

    Parameters:
        coil1: First InductionHeatingCoil
        coil2: Second InductionHeatingCoil
        frequency: Operating frequency [Hz]
        n_segments: Segments for Neumann integral
        verbose: Print results

    Returns:
        result: Dict with coupling analysis results
    """
    solver = WPTCoupledSolver(coil1, coil2, frequency)
    return solver.analyze(n_segments=n_segments, verbose=verbose)


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
