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
VALID_CORE_MODELS = (None, 'fembem', 'fem', 'radia', 'vector_fembem')


def _biot_savart_finite_filament(p1, p2, obs_point):
    """Compute H-field from a finite current filament carrying I = 1 A.

    Formula: H = (1/4*pi*d) * (cos(alpha1) - cos(alpha2)) * e_perp

    Same implementation as peec_coupled.biot_savart_finite_filament().

    Args:
        p1: Start point [x, y, z] (meters)
        p2: End point [x, y, z] (meters)
        obs_point: Observation point [x, y, z] (meters)

    Returns:
        H-field vector [Hx, Hy, Hz] in A/m (for I = 1 A)
    """
    p1 = np.asarray(p1, dtype=float)
    p2 = np.asarray(p2, dtype=float)
    r = np.asarray(obs_point, dtype=float)

    dl = p2 - p1
    L = np.linalg.norm(dl)
    if L < 1e-20:
        return np.zeros(3)
    e_l = dl / L

    r1 = r - p1
    r2 = r - p2

    cross = np.cross(e_l, r1)
    d = np.linalg.norm(cross)
    if d < 1e-20:
        return np.zeros(3)

    e_perp = cross / d

    r1_mag = np.linalg.norm(r1)
    r2_mag = np.linalg.norm(r2)
    if r1_mag < 1e-20 or r2_mag < 1e-20:
        return np.zeros(3)

    cos_alpha1 = np.dot(e_l, r1) / r1_mag
    cos_alpha2 = np.dot(e_l, r2) / r2_mag

    H = (1.0 / (4.0 * np.pi * d)) * (cos_alpha1 - cos_alpha2) * e_perp
    return H


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

        'vector_fembem' : Vector FEM-BEM (H(curl) + Maxwell SLP).
                   - Any linear mu_r (not restricted to 1)
                   - Unbounded domain (Weggler stabilized BEM)
                   - Full 3D vector eddy currents
                   - Best for: magnetic conducting cores (steel, iron)

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
            print(f"            Use core_model='vector_fembem', 'fem', or "
                  f"'radia' for mu_r != 1.")

        if core_model in ('fembem', 'fem', 'vector_fembem') and core_mesh is None:
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

    def compute_coupling_radia(self, solver_method=0, solver_prec=0.0001,
                                solver_maxiter=1000):
        """Compute coupling matrix Delta_L using Radia MMM.

        For each edge DOF j (unit current along edge j):
        1. Compute H-field at core via Biot-Savart (finite filament)
        2. Set as background field via rad.ObjBckg()
        3. Solve Radia for induced magnetization M
        4. Compute vector potential A from M at each edge center
        5. Delta_L[i,j] = dot(A(center_i), dir_i) * length_i

        This requires the Radia core object to be set up with material.

        Uses the same algorithm as peec_coupled.CoupledPEECSolver but
        adapted for edge-based DOFs (ngbem HDivSurface).

        Args:
            solver_method: Radia solver method (0=LU, 1=BiCGSTAB, 2=HACApK)
            solver_prec: Solver precision
            solver_maxiter: Maximum iterations
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

        if len(edge_data) == 0:
            raise RuntimeError(
                "No edges extracted from mesh. Check mesh validity.")

        n_loop = self.n_loop
        Delta_L = np.zeros((n_loop, n_loop))

        # Ensure radia_core is a list of object handles
        if isinstance(self.radia_core, int):
            core_handles = [self.radia_core]
        elif isinstance(self.radia_core, (list, tuple)):
            core_handles = list(self.radia_core)
        else:
            core_handles = [self.radia_core]

        print(f"  [ngbem_coupled] Computing Delta_L ({n_loop}x{n_loop}) "
              f"via Radia...")

        # Create persistent container for core objects (survives across iterations)
        mag_container = rad.ObjCnt(core_handles)

        for j in range(n_loop):
            p1_j = edge_data[j]['p1']
            p2_j = edge_data[j]['p2']

            # Biot-Savart H-field from unit current along edge j
            def h_field_from_edge_j(point, _p1=p1_j, _p2=p2_j):
                H = _biot_savart_finite_filament(_p1, _p2, point)
                B = MU_0 * H
                return [B[0], B[1], B[2]]

            # Create background field and solve
            bkg = rad.ObjBckg(h_field_from_edge_j)
            container = rad.ObjCnt([mag_container, bkg])
            rad.Solve(container, solver_prec, solver_maxiter, solver_method)

            # Extract A-field at each target edge center
            for i in range(n_loop):
                center_i = edge_data[i]['center'].tolist()
                A_vec = np.array(rad.Fld(container, 'a', center_i))
                dir_i = edge_data[i]['direction']
                len_i = edge_data[i]['length']
                Delta_L[i, j] = np.dot(A_vec, dir_i) * len_i

            # Only delete background and outer container (mag_container survives)
            rad.UtiDel(bkg)
            rad.UtiDel(container)

            if (j + 1) % max(1, n_loop // 10) == 0:
                print(f"    Edge {j+1}/{n_loop} done")

        self.Delta_L = Delta_L
        self._coupled = True
        self.t_coupling = time.perf_counter() - t_start
        print(f"  [ngbem_coupled] Delta_L computed in {self.t_coupling:.2f} s")

    def _extract_edge_geometry(self, mesh):
        """Extract edge midpoints, directions, and endpoints from mesh.

        For order=0 HDivSurface, each DOF corresponds to one mesh edge.

        Args:
            mesh: NGSolve Mesh

        Returns:
            list of dicts: {'center', 'direction', 'length', 'p1', 'p2'}
        """
        from ngbem_interface import extract_edge_geometry

        geom = extract_edge_geometry(mesh)
        n_edges = len(geom['centers'])

        edge_data = []
        for i in range(n_edges):
            edge_data.append({
                'center': geom['centers'][i],
                'direction': geom['directions'][i],
                'length': geom['lengths'][i],
                'p1': geom['p1'][i],
                'p2': geom['p2'][i],
            })
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
        with TaskManager():
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

    def _compute_mu_eff_vector_fembem(self, omega):
        """Compute mu_eff using vector FEM-BEM eddy current solver.

        Uses H(curl) FEM + Maxwell SLP BEM to solve the full 3D vector
        eddy current problem in the core. Handles arbitrary mu_r.

        Args:
            omega: Angular frequency [rad/s]

        Returns:
            mu_eff: Complex effective relative permeability
        """
        from ngbem_eddy import VectorEddyCurrentFEMBEM
        from ngsolve import Integrate, CF
        from ngsolve import TaskManager

        freq = omega / (2.0 * np.pi)

        # Check cache
        cache_key = ('vector_fembem', round(freq, 6))
        if cache_key in self._mu_eff_cache:
            return self._mu_eff_cache[cache_key]

        # Create and assemble vector FEM-BEM solver (lazily)
        if not hasattr(self, '_vector_fembem_solver'):
            self._vector_fembem_solver = VectorEddyCurrentFEMBEM(
                self.core_mesh,
                sigma=self.core_sigma,
                mu_r=self.mu_r,
                order=self.fem_order,
                conductor_label="conductor",
                surface_label="surface")
            self._vector_fembem_solver.assemble()

        solver = self._vector_fembem_solver

        # Solve with uniform z-directed B field
        B_ext = [0, 0, 1.0]
        solver.solve(freq, B_ext=B_ext)

        # Compute stored energy to derive effective permeability
        # For a volume V in field H_inc = B_ext/mu_0:
        # W = 0.5 * mu_eff * mu_0 * |H_inc|^2 * V
        # Also: W = 0.5 * Re(integral 1/mu * |curl A|^2 dV) (from FEM)
        H_inc = B_ext[2] / MU_0
        volume = Integrate(CF(1), self.core_mesh).real

        W = solver.compute_stored_energy()
        # mu_eff = 2 * W / (mu_0 * |H_inc|^2 * V)
        mu_eff = 2.0 * W / (MU_0 * H_inc**2 * volume)

        self._mu_eff_cache[cache_key] = mu_eff
        return mu_eff

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
        elif self.core_model == 'vector_fembem':
            return self._compute_mu_eff_vector_fembem(omega)
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

        For Radia core: Delta_L already includes full material response,
            L_total = L_air + Delta_L
        For FEM/FEMBEM core: Delta_L is geometric coupling,
            L_total = L_air + Delta_L * (mu_eff(omega) - 1)
        """
        R = np.diag(self.R_loop)
        L_total = self.L_air.copy()

        if self._coupled and self.Delta_L is not None:
            if self.core_model == 'radia':
                # Radia Delta_L already includes full material response
                L_total = L_total + self.Delta_L
            else:
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

        Z_LL = R + j*omega*L_total
        Z_SS = P / (j*omega)
        """
        R = np.diag(self.R_loop)
        L_total = self.L_air.copy()

        if self._coupled and self.Delta_L is not None:
            if self.core_model == 'radia':
                L_total = L_total + self.Delta_L
            else:
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
