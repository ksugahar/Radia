"""
ngsbem_peec_body_coupling.py

Direct coupling of PEEC filaments with FEM eddy current body.

Instead of using mu_eff as intermediary, this module couples PEEC filament
currents DIRECTLY with the FEM eddy current solver.

Architecture:
    PEEC Filaments (wire segments)    Conducting Body (FEM volume mesh)
    ==============================    =================================
    I_k (filament currents)           Hz_scat (eddy current field)
    L (inductance), R (resistance)    sigma, mu_r (material properties)

    Coupling:
    1. Biot-Savart: H_inc(r) = sum_k I_k * h_k(r) at body interior
    2. FEM solve: diffusion eqn with H_inc as source -> Hz_total
    3. Back-reaction: Delta_Z[i,j] = jw*mu * integral h_i * Hz_total_j dV

    Result:
    Z_coupled = Z_PEEC + Delta_Z(omega)

Filament geometry: list of wire segments [(start, end), ...]
Each segment represents a straight current-carrying wire.

Part of Radia project
"""

import numpy as np
import time

MU_0 = 4.0 * np.pi * 1e-7   # H/m



# Biot-Savart functions removed. Use ngsolve.bem MaxwellDL(J*dC) instead.


def create_circular_loop_segments(center, radius, axis, n_segments=36):
    """Create wire segments approximating a circular loop.

    Args:
        center: Loop center [x, y, z]
        radius: Loop radius [m]
        axis: Loop axis direction [ax, ay, az]
        n_segments: Number of straight segments

    Returns:
        segments: List of (p1, p2) tuples
    """
    center = np.asarray(center, dtype=float)
    axis = np.asarray(axis, dtype=float)
    axis = axis / np.linalg.norm(axis)

    if abs(axis[0]) < 0.9:
        t1 = np.cross(axis, [1, 0, 0])
    else:
        t1 = np.cross(axis, [0, 1, 0])
    t1 = t1 / np.linalg.norm(t1)
    t2 = np.cross(axis, t1)

    segments = []
    angles = np.linspace(0, 2 * np.pi, n_segments + 1)
    for k in range(n_segments):
        theta1 = angles[k]
        theta2 = angles[k + 1]
        p1 = center + radius * (np.cos(theta1) * t1 + np.sin(theta1) * t2)
        p2 = center + radius * (np.cos(theta2) * t1 + np.sin(theta2) * t2)
        segments.append((p1.tolist(), p2.tolist()))

    return segments


def create_rectangular_loop_segments(center, width, height, axis='z'):
    """Create wire segments for a rectangular loop.

    Args:
        center: Loop center [x, y, z]
        width: Width in first transverse direction [m]
        height: Height in second transverse direction [m]
        axis: Loop normal axis ('x', 'y', or 'z')

    Returns:
        segments: List of (p1, p2) tuples
    """
    cx, cy, cz = center
    w2 = width / 2
    h2 = height / 2

    if axis == 'z':
        p = [
            [cx - w2, cy - h2, cz],
            [cx + w2, cy - h2, cz],
            [cx + w2, cy + h2, cz],
            [cx - w2, cy + h2, cz],
        ]
    elif axis == 'x':
        p = [
            [cx, cy - w2, cz - h2],
            [cx, cy + w2, cz - h2],
            [cx, cy + w2, cz + h2],
            [cx, cy - w2, cz + h2],
        ]
    elif axis == 'y':
        p = [
            [cx - w2, cy, cz - h2],
            [cx + w2, cy, cz - h2],
            [cx + w2, cy, cz + h2],
            [cx - w2, cy, cz + h2],
        ]
    else:
        raise ValueError(f"axis must be 'x', 'y', or 'z', got '{axis}'")

    segments = []
    for k in range(4):
        segments.append((p[k], p[(k + 1) % 4]))

    return segments


def _project_biot_savart_to_gf(segments, mesh, fes):
    """Project source Hz field onto an H1 GridFunction via BiotSavartCF."""
    from ngsolve import GridFunction
    from radia.biot_savart import biot_savart_wire

    H_cf = biot_savart_wire(segments, current=1.0)
    gf = GridFunction(fes)
    with TaskManager():
        gf.Set(H_cf[2])  # Hz component
        return gf


class CoupledPEECBody:
    """Direct coupling of PEEC filaments with FEM eddy current body.

    Couples wire filament models (PEEC) with a conducting body (FEM)
    by computing the impedance change Delta_Z(omega) due to eddy currents.

    The coupled impedance:
        Z_coupled(omega) = Z_PEEC(omega) + Delta_Z(omega)

    where Delta_Z is the body's contribution, computed by:
    1. Biot-Savart: filament currents -> Hz_inc inside body
    2. FEM solve: diffusion equation with Hz_inc as Dirichlet BC
    3. Reciprocity: Delta_Z[i,j] = jw*mu * integral hz_inc_i * hz_total_j dV
    """

    def __init__(self, body_mesh, body_sigma, body_mu_r=1.0,
                 body_label="conductor", surface_label="surface",
                 fem_order=2):
        """Initialize coupled PEEC-body solver.

        Args:
            body_mesh: NGSolve Mesh (3D volume mesh of conducting body)
            body_sigma: Body conductivity [S/m]
            body_mu_r: Body relative permeability
            body_label: Material label for body volume
            surface_label: Boundary label for body surface
            fem_order: FE order for eddy current solver
        """
        self.body_mesh = body_mesh
        self.body_sigma = body_sigma
        self.body_mu_r = body_mu_r
        self.body_label = body_label
        self.surface_label = surface_label
        self.fem_order = fem_order

        self._filament_groups = []
        self._L_peec = None
        self._R_peec = None
        self._hz_inc_gfs = None  # Cached GridFunctions for Hz_inc

    def add_filament_group(self, segments, L_self=None, R_dc=None):
        """Add a group of wire segments representing one coil/filament.

        Args:
            segments: List of (p1, p2) tuples, each a wire segment
            L_self: Self-inductance of this group [H] (optional)
            R_dc: DC resistance of this group [Ohm] (optional)
        """
        self._filament_groups.append({
            'segments': segments,
            'L_self': L_self,
            'R_dc': R_dc,
        })
        self._hz_inc_gfs = None  # Invalidate cache

    def set_peec_matrices(self, L, R):
        """Set PEEC inductance and resistance matrices directly.

        Args:
            L: Inductance matrix [H], shape (n_groups, n_groups)
            R: Resistance matrix [Ohm], shape (n_groups, n_groups)
        """
        self._L_peec = np.asarray(L, dtype=float)
        self._R_peec = np.asarray(R, dtype=float)

    def prepare_coupling(self):
        """Pre-compute Hz_inc GridFunctions from Biot-Savart.

        Projects the Biot-Savart field from each filament group onto
        H1 GridFunctions for reuse across frequency sweeps.
        """
        from ngsolve import H1, Integrate, CF

        n_groups = len(self._filament_groups)
        if n_groups == 0:
            raise RuntimeError("No filament groups. Call add_filament_group().")

        # Create H1 space for source field projection
        fes = H1(self.body_mesh, order=self.fem_order, complex=True)
        volume = Integrate(CF(1), self.body_mesh).real

        print(f"Preparing coupling for {n_groups} filament group(s)...")
        print(f"  Body volume: {volume:.4e} m^3")
        print(f"  FEM DOFs: {fes.ndof}")

        self._hz_inc_gfs = []
        for g, group in enumerate(self._filament_groups):
            gf = _project_biot_savart_to_gf(
                group['segments'], self.body_mesh, fes)

            # Diagnostic: average Hz_inc
            Hz_avg = Integrate(gf, self.body_mesh) / volume
            print(f"  Group {g}: <Hz_inc> = {Hz_avg:.4e} A/m")

            self._hz_inc_gfs.append(gf)

        self._body_volume = volume
        print("  Done.")

    def compute_delta_Z(self, freq, printrates=False):
        """Compute coupling impedance Delta_Z at a single frequency.

        For each filament group pair (i, j):
        1. Hz_inc_j: Biot-Savart field from group j (pre-computed)
        2. FEM solve: diffusion with Dirichlet BC Hz = Hz_inc_j -> Hz_total_j
        3. Delta_Z[i,j] = jw*mu * integral hz_inc_i * hz_total_j dV

        Args:
            freq: Frequency [Hz]
            printrates: Print diagnostic info

        Returns:
            Delta_Z: Complex coupling impedance matrix (n_groups x n_groups)
        """
        from ngsolve import (H1, BilinearForm, GridFunction, LinearForm,
                              Integrate, CF, dx, grad)

        n_groups = len(self._filament_groups)
        omega = 2.0 * np.pi * freq
        mu = self.body_mu_r * MU_0
        k_sq = 1j * omega * mu * self.body_sigma

        if self._hz_inc_gfs is None:
            self.prepare_coupling()

        Delta_Z = np.zeros((n_groups, n_groups), dtype=complex)

        if omega < 1e-30 or self.body_sigma < 1e-30:
            # DC or non-conducting: no eddy currents
            return Delta_Z

        # Create FEM system (shared across groups)
        fes = H1(self.body_mesh, order=self.fem_order, complex=True,
                  dirichlet=self.surface_label)
        u, v = fes.TnT()
        a = BilinearForm(grad(u) * grad(v) * dx + k_sq * u * v * dx)
        a.Assemble()
        inv = a.mat.Inverse(freedofs=fes.FreeDofs(), inverse="pardiso")

        # For each source group j: solve FEM with Hz_inc_j as Dirichlet BC
        Hz_total_gfs = []
        for j in range(n_groups):
            gfu = GridFunction(fes)

            # Set Dirichlet BC: Hz_total = Hz_inc_j on boundary
            gfu.Set(self._hz_inc_gfs[j],
                    definedon=self.body_mesh.Boundaries(self.surface_label))

            # RHS = 0 (total field formulation)
            f = LinearForm(fes)
            f.Assemble()

            # Solve: A * u_free = -A * u_D
            r = f.vec - a.mat * gfu.vec
            gfu.vec.data += inv * r

            Hz_total_gfs.append(gfu)

            if printrates:
                Hz_avg = Integrate(gfu, self.body_mesh) / self._body_volume
                Hz_inc_avg = Integrate(
                    self._hz_inc_gfs[j], self.body_mesh) / self._body_volume
                mu_eff_j = self.body_mu_r * Hz_avg / Hz_inc_avg
                print(f"  Group {j}: mu_eff = "
                      f"{mu_eff_j.real:.6f}{mu_eff_j.imag:+.6f}j")

        # Compute Delta_Z[i,j] = jw * mu * integral hz_inc_i * hz_total_j dV
        for i in range(n_groups):
            for j in range(n_groups):
                integral = Integrate(
                    self._hz_inc_gfs[i] * Hz_total_gfs[j],
                    self.body_mesh)
                Delta_Z[i, j] = 1j * omega * mu * integral

        return Delta_Z

    def solve_frequency(self, freqs, printrates=False):
        """Compute coupled impedance Z_coupled(f) over a frequency range.

        Z_coupled = Z_PEEC + Delta_Z(omega)

        Args:
            freqs: Array of frequencies [Hz]
            printrates: Print solver info

        Returns:
            Z_coupled: Complex impedance array.
                       If n_groups == 1: shape (len(freqs),)
                       If n_groups > 1:  shape (len(freqs), n_groups, n_groups)
        """
        freqs = np.asarray(freqs, dtype=float)
        n_groups = len(self._filament_groups)
        n_freqs = len(freqs)

        if n_groups == 0:
            raise RuntimeError("No filament groups. Call add_filament_group().")

        if self._hz_inc_gfs is None:
            self.prepare_coupling()

        Z_all = np.zeros((n_freqs, n_groups, n_groups), dtype=complex)

        for k, f in enumerate(freqs):
            omega = 2.0 * np.pi * f
            if printrates:
                print(f"\n--- Frequency {f:.1f} Hz ---")

            # PEEC impedance
            Z_peec = np.zeros((n_groups, n_groups), dtype=complex)
            if self._L_peec is not None and self._R_peec is not None:
                Z_peec = self._R_peec + 1j * omega * self._L_peec
            else:
                for g in range(n_groups):
                    grp = self._filament_groups[g]
                    L = grp.get('L_self', 0.0) or 0.0
                    R = grp.get('R_dc', 0.0) or 0.0
                    Z_peec[g, g] = R + 1j * omega * L

            # Body coupling
            Delta_Z = self.compute_delta_Z(f, printrates=printrates)

            Z_all[k] = Z_peec + Delta_Z

        if n_groups == 1:
            return Z_all[:, 0, 0]

        return Z_all

    def compute_mu_eff(self, freq, group_idx=0):
        """Compute effective permeability from FEM solution (diagnostic).

        mu_eff = mu_r * <Hz_total> / <Hz_inc>

        Args:
            freq: Frequency [Hz]
            group_idx: Which filament group to use as source

        Returns:
            mu_eff: Complex effective permeability
        """
        from ngsolve import (H1, BilinearForm, GridFunction, LinearForm,
                              Integrate, CF, dx, grad)
        from ngsolve import TaskManager

        if self._hz_inc_gfs is None:
            self.prepare_coupling()

        omega = 2.0 * np.pi * freq
        mu = self.body_mu_r * MU_0
        k_sq = 1j * omega * mu * self.body_sigma

        fes = H1(self.body_mesh, order=self.fem_order, complex=True,
                  dirichlet=self.surface_label)
        u, v = fes.TnT()
        a = BilinearForm(grad(u) * grad(v) * dx + k_sq * u * v * dx)
        a.Assemble()
        inv = a.mat.Inverse(freedofs=fes.FreeDofs(), inverse="pardiso")

        gfu = GridFunction(fes)
        gfu.Set(self._hz_inc_gfs[group_idx],
                definedon=self.body_mesh.Boundaries(self.surface_label))
        f = LinearForm(fes)
        f.Assemble()
        r = f.vec - a.mat * gfu.vec
        gfu.vec.data += inv * r

        volume = self._body_volume
        Hz_avg = Integrate(gfu, self.body_mesh) / volume
        Hz_inc_avg = Integrate(
            self._hz_inc_gfs[group_idx], self.body_mesh) / volume

        return self.body_mu_r * Hz_avg / Hz_inc_avg
