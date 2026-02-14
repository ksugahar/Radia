"""
ngbem_coupled.py

Phase 3: Unified PEEC + Core Coupled Solver

Combines:
- Phase 1: ngbem PEEC Loop-Star (conductor coil impedance)
- Phase 2: ngbem FEM/FEM-BEM eddy current (conducting core)
- Radia MMM/MSC (magnetic core, nonlinear support)
- Coupling: Delta_L mutual inductance

Core model selection:
    'fembem'  : FEM-BEM coupling (ngbem Calderon). mu_r=1 ONLY. Unbounded domain.
    'fem'     : FEM with Dirichlet BC. Any mu_r. Bounded domain approximation.
    'radia'   : Radia MMM/MSC. Any mu_r, nonlinear. Unbounded domain.
    None      : No core (PEEC only, or static mu_r without eddy currents).

Block system:
    | Z_LL + Delta_L    M_LS^T |   | I_loop |   | V_port |
    | M_LS               Z_SS  | * | Q_star | = | 0      |

where:
    Z_LL = R_diag + j*omega*L_total
    L_total = L_air + Delta_L * (mu_eff(omega) - 1)
    M_LS = divergence coupling               (charge conservation)
    Z_SS = P / (j*omega)                     (capacitive)

mu_eff(omega) computation:
    FEM/FEMBEM: Solve diffusion in core mesh -> volume-average Hz -> mu_eff
    Radia:      Solve MMM for induced M -> field integration -> mu_eff
    None:       mu_eff = mu_r (static, frequency-independent)

Part of Radia project
"""

import numpy as np
import time

# Physical constants
MU_0 = 4.0 * np.pi * 1e-7
EPS_0 = 8.854187817e-12

# Valid core model choices
VALID_CORE_MODELS = (None, 'fembem', 'fem', 'radia')


class CoupledPEECMMM:
    """Coupled PEEC (conductor) + Core solver.

    Uses Phase 1 (ngbem_peec.py) for conductor Loop-Star matrices
    and a selectable core model for eddy current / magnetic effects.

    Core model options:
        'fembem' : FEM-BEM coupling via ngbem Calderon projector.
                   - Linear materials only (mu_r=1 for scalar Hz formulation)
                   - Unbounded domain (natural exterior BEM)
                   - High-order BEM elements available
                   - Best for: aluminum shields, copper conductors

        'fem'    : FEM-only with Dirichlet BC.
                   - Any linear mu_r
                   - Bounded domain (Dirichlet BC = Hz_inc on boundary)
                   - Simpler, faster, good for validation
                   - Best for: conducting magnetic cores (steel, iron)

        'radia'  : Radia MMM/MSC.
                   - Any mu_r, including nonlinear (MatSatIsoTab)
                   - Unbounded domain (integral equation)
                   - Best for: nonlinear cores, permanent magnets

        None     : No core model. Static mu_r, no eddy currents.
    """

    def __init__(self, peec_solver,
                 core_model=None,
                 core_mesh=None,
                 core_sigma=0.0,
                 mu_r=1.0,
                 mu_r_imag=0.0,
                 radia_core=None,
                 fem_order=2):
        """Initialize coupled solver.

        Args:
            peec_solver: NGBEMPEECSolver instance (assembled)
            core_model: Core solver type ('fembem', 'fem', 'radia', or None)
            core_mesh: NGSolve 3D volume Mesh (for 'fembem' or 'fem')
            core_sigma: Core conductivity [S/m] (for 'fembem' or 'fem')
            mu_r: Relative permeability of core (real part)
            mu_r_imag: Imaginary permeability (loss tangent = mu_r_imag/mu_r)
            radia_core: Radia object handle for magnetic core (for 'radia')
            fem_order: Finite element order for core solver (for 'fembem'/'fem')
        """
        if core_model not in VALID_CORE_MODELS:
            raise ValueError(
                f"core_model must be one of {VALID_CORE_MODELS}, "
                f"got '{core_model}'")

        if core_model == 'fembem' and mu_r != 1.0:
            print(f"  [WARNING] core_model='fembem' uses scalar Hz formulation "
                  f"valid ONLY for mu_r=1.")
            print(f"            Got mu_r={mu_r}. Results may be inaccurate.")
            print(f"            Use core_model='fem' or 'radia' for mu_r != 1.")

        if core_model in ('fembem', 'fem') and core_mesh is None:
            raise ValueError(
                f"core_mesh required for core_model='{core_model}'")

        if core_model == 'radia' and radia_core is None:
            raise ValueError(
                "radia_core required for core_model='radia'")

        self.peec = peec_solver
        self.core_model = core_model
        self.core_mesh = core_mesh
        self.core_sigma = core_sigma
        self.mu_r = mu_r
        self.mu_r_imag = mu_r_imag
        self.radia_core = radia_core
        self.fem_order = fem_order

        # Get PEEC matrices
        mats = peec_solver.get_matrices()
        self.L_air = mats['L']          # n_loop x n_loop
        self.P = mats['P']              # n_star x n_star
        self.M_LS = mats['M_LS']        # n_star x n_loop
        self.R_loop = mats['R_loop']    # n_loop diagonal
        self.n_loop = mats['n_loop']
        self.n_star = mats['n_star']

        # Coupling matrix (computed lazily)
        self.Delta_L = None
        self._coupled = False

        # mu_eff cache (freq -> mu_eff), keyed by (core_model, freq)
        self._mu_eff_cache = {}

        # Results
        self.t_coupling = 0.0

    def compute_coupling_analytically(self, L_core):
        """Set the coupling matrix from an analytically computed value.

        For simple geometries, Delta_L can be computed analytically:
            Delta_L = L_air / (mu_r + 1)  (image method)

        For a single-loop system, Delta_L is a scalar.

        Args:
            L_core: Core coupling matrix or scalar [H]
        """
        L_core = np.atleast_2d(L_core)
        if L_core.shape[0] != self.n_loop:
            if L_core.size == 1:
                self.Delta_L = L_core[0, 0] * np.eye(self.n_loop)
            else:
                raise ValueError(
                    f"L_core shape {L_core.shape} does not match "
                    f"n_loop={self.n_loop}")
        else:
            self.Delta_L = L_core

        self._coupled = True

    def compute_coupling_radia(self):
        """Compute coupling matrix Delta_L using Radia MMM.

        For each loop DOF i, j:
        1. Apply unit current in loop j
        2. Compute H field at core element centers (Biot-Savart)
        3. Solve Radia for induced magnetization M
        4. Compute A field from M at loop i center
        5. Delta_L[i,j] = dot(A, loop_direction_i) * loop_length_i

        This requires the Radia core object to be set up with material.
        """
        try:
            import radia as rad
        except ImportError:
            raise ImportError(
                "Radia module required for compute_coupling_radia()")

        if self.radia_core is None:
            raise RuntimeError("No Radia core object provided")

        t_start = time.perf_counter()

        # Get edge (loop) geometry from PEEC mesh
        mesh = self.peec.mesh
        edge_data = self._extract_edge_geometry(mesh)

        n_loop = self.n_loop
        Delta_L = np.zeros((n_loop, n_loop), dtype=complex)

        # For each source loop j: compute H at core, solve M, compute A
        for j in range(n_loop):
            rad.UtiDelAll()
            # TODO: Use rad.ObjBckg with Biot-Savart field from loop j
            pass  # Placeholder for full Radia coupling

        self.Delta_L = Delta_L
        self._coupled = True
        self.t_coupling = time.perf_counter() - t_start

    def _extract_edge_geometry(self, mesh):
        """Extract edge midpoints and directions from mesh."""
        edge_data = []
        return edge_data

    def _compute_mu_eff_fem(self, omega, mode='fem'):
        """Compute mu_eff from FEM or FEM-BEM eddy current solution.

        Solves the diffusion equation in the core mesh:
            laplacian(Hz) = j*omega*mu*sigma*Hz

        FEM mode: Dirichlet BC Hz = Hz_inc on boundary (bounded domain).
        FEMBEM mode: Calderon projector coupling (unbounded domain).

        Then computes:
            mu_eff = mu_r * <Hz_total>_volume / Hz_inc

        Args:
            omega: Angular frequency [rad/s]
            mode: 'fem' or 'fembem'

        Returns:
            mu_eff: Complex effective relative permeability
        """
        from ngbem_eddy import EddyCurrentFEMBEM
        from ngsolve import Integrate, CF

        freq = omega / (2.0 * np.pi)

        # Check cache
        cache_key = (mode, round(freq, 6))
        if cache_key in self._mu_eff_cache:
            return self._mu_eff_cache[cache_key]

        # Create eddy current solver
        eddy = EddyCurrentFEMBEM(
            self.core_mesh,
            sigma=self.core_sigma,
            mu_r=self.mu_r,
            order=self.fem_order,
            conductor_label="conductor",
            surface_label="surface")

        # Assemble and solve
        B_ext = [0, 0, 1.0]  # 1T z-directed
        if mode == 'fembem':
            eddy.assemble_fembem(freq)
            eddy.solve(B_ext=B_ext, mode='fembem')
        else:
            eddy.assemble_fem(freq)
            eddy.solve(B_ext=B_ext, mode='fem')

        # Get total Hz field
        Hz_total = eddy.get_total_field()
        Hz_inc = B_ext[2] / MU_0

        # Volume average: <Hz> = Integrate(Hz, mesh) / Volume
        volume = Integrate(CF(1), self.core_mesh).real
        Hz_avg = Integrate(Hz_total, self.core_mesh) / volume

        # Effective permeability
        mu_eff = self.mu_r * Hz_avg / Hz_inc

        # Store eddy solver for diagnostics
        self._last_eddy_solver = eddy

        # Cache result
        self._mu_eff_cache[cache_key] = mu_eff

        return mu_eff

    def _compute_mu_eff_radia(self, omega):
        """Compute mu_eff from Radia MMM/MSC solution.

        For conducting magnetic cores, uses Radia to solve for
        magnetization and compute effective permeability.

        Currently returns static mu_r (Radia eddy current coupling
        requires full implementation of compute_coupling_radia).

        Args:
            omega: Angular frequency [rad/s]

        Returns:
            mu_eff: Complex effective relative permeability
        """
        # TODO: Full Radia eddy current computation
        # For now, return static mu_r for Radia core model
        return self.mu_r - 1j * self.mu_r_imag

    def _get_mu_eff(self, omega):
        """Get effective permeability at given frequency.

        Dispatches to the appropriate solver based on core_model.

        Args:
            omega: Angular frequency [rad/s]

        Returns:
            mu_eff: Complex effective relative permeability
        """
        if omega <= 0 or self.core_sigma <= 0:
            return self.mu_r - 1j * self.mu_r_imag

        if self.core_model == 'fembem':
            return self._compute_mu_eff_fem(omega, mode='fembem')
        elif self.core_model == 'fem':
            return self._compute_mu_eff_fem(omega, mode='fem')
        elif self.core_model == 'radia':
            return self._compute_mu_eff_radia(omega)
        else:
            # No core model: static mu_r
            return self.mu_r - 1j * self.mu_r_imag

    def solve_frequency(self, freqs, mode='mqs'):
        """Compute port impedance Z(f) for coupled system.

        Args:
            freqs: Array of frequencies [Hz]
            mode: 'mqs' (Loop-only) or 'full' (Loop-Star)

        Returns:
            Z_port: Complex impedance array
        """
        freqs = np.atleast_1d(freqs)
        Z_port = np.zeros(len(freqs), dtype=complex)

        for i, freq in enumerate(freqs):
            omega = 2.0 * np.pi * freq

            if mode == 'mqs':
                Z_port[i] = self._solve_mqs(omega)
            else:
                Z_port[i] = self._solve_full(omega)

        return Z_port

    def _solve_mqs(self, omega):
        """MQS (Loop-only) impedance with core coupling.

        Z_branch = R_diag + j*omega*(L_air + Delta_L * (mu_eff(omega) - 1))
        """
        R = np.diag(self.R_loop)
        L_total = self.L_air.copy()

        if self._coupled and self.Delta_L is not None:
            mu_eff = self._get_mu_eff(omega)
            L_total = L_total + self.Delta_L * (mu_eff - 1.0)

        Z_branch = R + 1j * omega * L_total

        if self.n_loop == 1:
            return Z_branch[0, 0]
        else:
            try:
                Y_branch = np.linalg.inv(Z_branch)
                e = np.ones(self.n_loop) / self.n_loop
                Z_port = 1.0 / (e @ Y_branch @ e)
                return Z_port
            except np.linalg.LinAlgError:
                return Z_branch[0, 0]

    def _solve_full(self, omega):
        """Full Loop-Star impedance with core coupling.

        Block system:
        | Z_LL    M_LS^T |   | I_L |   | V_port |
        | M_LS    Z_SS   | * | Q_S | = | 0      |

        Z_LL = R + j*omega*(L_air + Delta_L * (mu_eff - 1))
        Z_SS = P / (j*omega)
        """
        R = np.diag(self.R_loop)
        L_total = self.L_air.copy()

        if self._coupled and self.Delta_L is not None:
            mu_eff = self._get_mu_eff(omega)
            L_total = L_total + self.Delta_L * (mu_eff - 1.0)

        Z_LL = R + 1j * omega * L_total

        if omega > 0:
            Z_SS = self.P / (1j * omega)
        else:
            return self._solve_mqs(omega)

        try:
            Z_SS_inv = np.linalg.inv(Z_SS)
            Z_schur = Z_LL - self.M_LS.T @ Z_SS_inv @ self.M_LS
        except np.linalg.LinAlgError:
            Z_schur = Z_LL

        if self.n_loop == 1:
            return Z_schur[0, 0]
        else:
            try:
                Y = np.linalg.inv(Z_schur)
                e = np.ones(self.n_loop) / self.n_loop
                return 1.0 / (e @ Y @ e)
            except np.linalg.LinAlgError:
                return Z_schur[0, 0]

    def get_inductance(self, freq=0):
        """Get effective inductance at given frequency.

        L_eff = Im(Z_port) / omega

        Args:
            freq: Frequency [Hz] (0 = DC inductance)

        Returns:
            L_eff: Effective inductance [H]
        """
        if freq == 0:
            freq = 1.0
        omega = 2.0 * np.pi * freq
        Z = self._solve_mqs(omega)
        return np.imag(Z) / omega

    def get_resistance(self, freq=0):
        """Get effective resistance at given frequency.

        R_eff = Re(Z_port)

        Args:
            freq: Frequency [Hz]

        Returns:
            R_eff: Effective resistance [Ohm]
        """
        if freq == 0:
            freq = 1.0
        omega = 2.0 * np.pi * freq
        Z = self._solve_mqs(omega)
        return np.real(Z)

    def print_summary(self):
        """Print summary of coupled system."""
        print(f"CoupledPEECMMM Summary:")
        print(f"  Loop DOFs: {self.n_loop}")
        print(f"  Star DOFs: {self.n_star}")
        print(f"  Core model: {self.core_model or 'none (static mu_r)'}")
        print(f"  Core mu_r: {self.mu_r:.1f}")
        if self.mu_r_imag > 0:
            print(f"  Core mu_r_imag: {self.mu_r_imag:.1f}")
            print(f"  Loss tangent: {self.mu_r_imag/self.mu_r:.4f}")

        if self.core_model in ('fembem', 'fem'):
            print(f"  Core sigma: {self.core_sigma:.2e} S/m")
            print(f"  FE order: {self.fem_order}")
            if self.core_sigma > 0:
                delta_ref = np.sqrt(2.0 / (2*np.pi*1000 * self.mu_r * MU_0
                                            * self.core_sigma))
                print(f"  Skin depth at 1kHz: {delta_ref*1e3:.2f} mm")
            if self.core_model == 'fembem':
                print(f"  Exterior BC: BEM (Calderon, unbounded)")
            else:
                print(f"  Exterior BC: Dirichlet (bounded)")
        elif self.core_model == 'radia':
            print(f"  Radia core: object #{self.radia_core}")

        print(f"  Coupled: {self._coupled}")

        if self._coupled:
            L_air_diag = np.diag(self.L_air)
            print(f"  L_air (diagonal): min={np.min(L_air_diag):.4e}, "
                  f"max={np.max(L_air_diag):.4e} H")
            if self.Delta_L is not None:
                DL_diag = np.real(np.diag(self.Delta_L))
                print(f"  Delta_L (diagonal): min={np.min(DL_diag):.4e}, "
                      f"max={np.max(DL_diag):.4e} H")
                ratio = np.max(np.abs(DL_diag)) / np.max(np.abs(L_air_diag))
                print(f"  Delta_L / L_air ratio: {ratio:.3f}")


def compute_delta_L(L_air, mu_r):
    """Compute geometric coupling matrix Delta_L.

    Delta_L represents the geometric coupling between coil and core.
    For a flat plate coil above a magnetic half-space:
        Delta_L = L_air / (mu_r + 1)

    This is derived from:
        L_total(DC) = L_air + Delta_L * (mu_r - 1)
                    = L_air * 2*mu_r / (mu_r + 1)  (image method)
    =>  Delta_L = L_air / (mu_r + 1)

    Works for all mu_r including mu_r=1 (non-magnetic conductor).

    Args:
        L_air: Air inductance matrix [H]
        mu_r: Relative permeability

    Returns:
        Delta_L: Coupling matrix [H]
    """
    return L_air / (mu_r + 1.0)
