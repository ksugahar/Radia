"""
ngbem_eddy.py

FEM-BEM Eddy Current Solver using ngsolve/ngbem.

Two solver modes:
1. FEM-BEM coupled: Costabel symmetric coupling (separate H1 + SurfaceL2 spaces)
2. FEM-only: Dirichlet BC from known external field (for validation)

Scalar Hz formulation:
    Interior (conductor): laplacian(Hz) = j*omega*mu*sigma * Hz  (diffusion)
    Exterior (air):       laplacian(Hz) = 0                       (Laplace, via BEM)

Costabel symmetric coupling (separate spaces, BlockMatrix):

    | A_FEM + D    (-0.5*M + K)^T |   | Hz_scat |   | -k^2 * Hz_inc |
    | (-0.5*M + K)     -V         | * | lambda  | = | 0             |

where:
    A_FEM = FEM stiffness + k^2 mass (interior diffusion)
    D = hypersingular operator (exterior BEM)
    V = single layer potential operator (exterior BEM)
    K = double layer potential operator (exterior BEM)
    M = boundary mass matrix (trace coupling)
    k^2 = j*omega*mu*sigma

Physics:
    - External field Hz_inc (uniform or from PEEC filament source)
    - Scattered field Hz_scat satisfies diffusion eqn with source -k^2 * Hz_inc
    - Total field: Hz_total = Hz_scat + Hz_inc
    - Current: J = curl(H) => Jx = dHz/dy, Jy = -dHz/dx
    - Loss: P = (1/(2*sigma)) * integral(|grad Hz_total|^2) dV

Designed for coupling with PEEC filament models:
    - accept_source_field(): set Hz_inc as CoefficientFunction (from Biot-Savart)
    - solve() returns Hz_scat, lambda; total field = Hz_scat + Hz_inc

Part of Radia project
"""

import numpy as np
import time

# Physical constants
MU_0 = 4.0 * np.pi * 1e-7   # H/m
EPS_0 = 8.854187817e-12      # F/m


def create_conductor_mesh(width, height, depth, maxh=0.002,
                          conductor_label="conductor",
                          surface_label="surface"):
    """Create a rectangular conducting block with volume mesh.

    Args:
        width: Block width (x) [m]
        height: Block height (y) [m]
        depth: Block depth (z) [m]
        maxh: Maximum element size [m]
        conductor_label: Material label for conductor volume
        surface_label: Boundary label for conductor surface

    Returns:
        mesh: NGSolve Mesh (3D volume mesh with labeled boundary)
    """
    from netgen.occ import Box, Pnt, OCCGeometry
    from ngsolve import Mesh

    block = Box(Pnt(-width/2, -height/2, -depth/2),
                Pnt(width/2, height/2, depth/2))
    block.solids.name = conductor_label
    block.faces.name = surface_label

    geo = OCCGeometry(block)
    ngmesh = geo.GenerateMesh(maxh=maxh)
    return Mesh(ngmesh)


class EddyCurrentFEMBEM:
    """FEM-BEM solver for eddy currents in conducting bodies.

    Scalar Hz formulation with Costabel symmetric coupling.

    Uses SEPARATE function spaces (H1 + SurfaceL2) and BlockMatrix,
    following the ngsolve BEM tutorial approach.

    Modes:
        'fembem': Full FEM-BEM with Costabel symmetric coupling
        'fem':    FEM-only with Dirichlet BC Hz=Hz_inc on boundary
    """

    def __init__(self, mesh, sigma=5.8e7, mu_r=1.0, order=2,
                 conductor_label="conductor", surface_label="surface"):
        """Initialize FEM-BEM eddy current solver.

        Args:
            mesh: NGSolve Mesh (3D volume mesh of conductor)
            sigma: Electrical conductivity [S/m]
            mu_r: Relative permeability (real part)
            order: Finite element order
            conductor_label: Material label for conductor volume
            surface_label: Boundary label for conductor surface
        """
        self.mesh = mesh
        self.sigma = sigma
        self.mu_r = mu_r
        self.mu = mu_r * MU_0
        self.order = order
        self.conductor_label = conductor_label
        self.surface_label = surface_label

        # State
        self.omega = 0.0
        self.freq = 0.0
        self._assembled = False
        self._solved = False
        self.gfu_h1 = None   # H1 GridFunction (Hz_scat or Hz_total)
        self.gfu_l2 = None   # SurfaceL2 GridFunction (lambda = dHz/dn)
        self._mode = None
        self._Hz_inc = None

        # Custom source field (CoefficientFunction, for PEEC coupling)
        self._Hz_inc_cf = None

        # Timing
        self.t_assemble = 0.0
        self.t_solve = 0.0

    def assemble_fem(self, freq):
        """Assemble FEM-only system (Dirichlet BC) for given frequency.

        Args:
            freq: Frequency [Hz]
        """
        from ngsolve import (H1, BilinearForm, GridFunction, dx, grad)

        t_start = time.perf_counter()

        self.freq = freq
        self.omega = 2.0 * np.pi * freq
        self._mode = 'fem'

        # Skin depth
        if self.omega > 0 and self.sigma > 0:
            self.delta = np.sqrt(2.0 / (self.omega * self.mu * self.sigma))
        else:
            self.delta = float('inf')

        # k^2 = j*omega*mu*sigma
        self._k_sq = 1j * self.omega * self.mu * self.sigma

        # H1 FE space (complex-valued)
        self._fes = H1(self.mesh, order=self.order, complex=True,
                        dirichlet=self.surface_label)
        u, v = self._fes.TnT()

        # Bilinear form: (grad u, grad v) + k^2 (u, v) = 0
        self._a = BilinearForm(
            grad(u) * grad(v) * dx + self._k_sq * u * v * dx)
        self._a.Assemble()

        self._assembled = True
        self.t_assemble = time.perf_counter() - t_start

    def assemble_fembem(self, freq, intorder=12):
        """Assemble FEM-BEM coupled system for given frequency.

        Uses SEPARATE H1 and SurfaceL2 spaces with BlockMatrix,
        following the ngsolve BEM tutorial exactly.

        Args:
            freq: Frequency [Hz]
            intorder: Integration order for BEM singular quadrature
        """
        from ngsolve import (H1, SurfaceL2, BilinearForm, GridFunction,
                              TaskManager, ds, dx, grad)
        from ngsolve.bem import (SingleLayerPotentialOperator,
                                  DoubleLayerPotentialOperator,
                                  HypersingularOperator)

        t_start = time.perf_counter()

        self.freq = freq
        self.omega = 2.0 * np.pi * freq
        self._mode = 'fembem'

        # Skin depth
        if self.omega > 0 and self.sigma > 0:
            self.delta = np.sqrt(2.0 / (self.omega * self.mu * self.sigma))
        else:
            self.delta = float('inf')

        self._k_sq = 1j * self.omega * self.mu * self.sigma

        # --- SEPARATE function spaces (ngsolve BEM tutorial approach) ---
        self._fes_h1 = H1(self.mesh, order=self.order, complex=True)
        self._fes_l2 = SurfaceL2(
            self.mesh, order=self.order - 1,
            complex=True, dual_mapping=True,
            definedon=self.mesh.Boundaries(self.surface_label))

        u, v = self._fes_h1.TnT()
        uL2, vL2 = self._fes_l2.TnT()

        # --- FEM interior bilinear form ---
        self._a_fem = BilinearForm(
            grad(u) * grad(v) * dx + self._k_sq * u * v * dx)
        self._a_fem.Assemble()

        # --- BEM operators (exterior Laplace, separate space API) ---
        print("  Assembling BEM operators...")
        with TaskManager():
            self._V_op = SingleLayerPotentialOperator(
                self._fes_l2, intorder=intorder)
            self._K_op = DoubleLayerPotentialOperator(
                self._fes_h1, self._fes_l2,
                trial_definedon=self.mesh.Boundaries(self.surface_label),
                test_definedon=self.mesh.Boundaries(self.surface_label),
                intorder=intorder)
            self._D_op = HypersingularOperator(
                self._fes_h1,
                definedon=self.mesh.Boundaries(self.surface_label),
                intorder=intorder)
            self._M_bf = BilinearForm(
                self._fes_h1.TrialFunction()
                * self._fes_l2.TestFunction().Trace()
                * ds(self.surface_label)).Assemble()
        print("  Done.")

        # --- Preconditioner blocks ---
        self._pre_h1 = BilinearForm(
            (grad(u) * grad(v) + 1e-6 * u * v) * dx).Assemble()
        self._pre_l2 = BilinearForm(
            uL2 * vL2.Trace() * ds(self.surface_label)).Assemble()

        self._assembled = True
        self.t_assemble = time.perf_counter() - t_start

    def set_source_field(self, Hz_inc_cf):
        """Set custom source field for PEEC coupling.

        Instead of a uniform Hz_inc, use a CoefficientFunction
        representing the incident field from PEEC filaments (Biot-Savart).

        Args:
            Hz_inc_cf: ngsolve CoefficientFunction giving Hz_inc(x,y,z)
        """
        self._Hz_inc_cf = Hz_inc_cf

    def solve(self, B_ext=None, mode='fem', printrates=False):
        """Solve for eddy current distribution under external field.

        Args:
            B_ext: External B field [Bx, By, Bz] in Tesla (uniform).
                   Ignored if set_source_field() was called.
            mode: 'fem' (Dirichlet BC) or 'fembem' (Costabel coupling)
            printrates: Print solver convergence rates

        Returns:
            gfu_h1: GridFunction with Hz solution (scat for fembem, total for fem)
        """
        if not self._assembled:
            raise RuntimeError("Call assemble_fem() or assemble_fembem() first")

        if B_ext is None:
            B_ext = [0, 0, 1.0]
        B_ext = np.asarray(B_ext, dtype=float)

        # Hz_inc = Bz / mu_0 (for uniform field)
        self._Hz_inc = B_ext[2] / MU_0

        t_start = time.perf_counter()

        if self._mode == 'fem':
            self._solve_fem(printrates)
        else:
            self._solve_fembem(printrates)

        self._solved = True
        self.t_solve = time.perf_counter() - t_start
        return self.gfu_h1

    def _solve_fem(self, printrates=False):
        """Solve FEM-only with Dirichlet BC: Hz = Hz_inc on boundary."""
        from ngsolve import (GridFunction, LinearForm, dx, CF)

        u, v = self._fes.TnT()

        self.gfu_h1 = GridFunction(self._fes)

        # Set Dirichlet BC: Hz = Hz_inc on boundary
        if self._Hz_inc_cf is not None:
            self.gfu_h1.Set(self._Hz_inc_cf,
                            definedon=self.mesh.Boundaries(self.surface_label))
        else:
            self.gfu_h1.Set(CF(self._Hz_inc),
                            definedon=self.mesh.Boundaries(self.surface_label))

        # RHS = 0 (total field formulation, no volume source)
        f = LinearForm(self._fes)
        f.Assemble()

        # Modify RHS for Dirichlet BC
        r = f.vec - self._a.mat * self.gfu_h1.vec

        # Solve
        inv = self._a.mat.Inverse(
            freedofs=self._fes.FreeDofs(), inverse="pardiso")
        self.gfu_h1.vec.data += inv * r

    def _solve_fembem(self, printrates=False):
        """Solve FEM-BEM coupled system using BlockMatrix + GMRes.

        Uses separate H1 and SurfaceL2 spaces with Costabel symmetric coupling:

            | A + D          (-0.5*M + K)^T |   | Hz_scat |   | f_H1 |
            | (-0.5*M + K)   -V             | * | lambda  | = | 0    |

        where f_H1 = -k^2 * (Hz_inc, v) is the FEM volume source.
        """
        from ngsolve import (GridFunction, LinearForm, TaskManager,
                              BlockMatrix, BlockVector, dx, CF)
        from ngsolve.solvers import GMRes

        # --- RHS ---
        v = self._fes_h1.TestFunction()
        vL2 = self._fes_l2.TestFunction()

        # FEM RHS: -k^2 * (Hz_inc, v)_volume
        f_H1 = LinearForm(self._fes_h1)
        if self._Hz_inc_cf is not None:
            f_H1 += CF(-self._k_sq) * self._Hz_inc_cf * v * dx
        else:
            f_H1 += CF(-self._k_sq * self._Hz_inc) * v * dx
        f_H1.Assemble()

        # BEM RHS: 0
        f_L2 = LinearForm(self._fes_l2)
        f_L2.Assemble()

        # --- Block system (Costabel symmetric coupling) ---
        # Off-diagonal: (-0.5*M + K)
        offdiag = -0.5 * self._M_bf.mat + self._K_op.mat

        lhs = BlockMatrix([
            [self._a_fem.mat + self._D_op.mat, offdiag.T],
            [offdiag, (-1) * self._V_op.mat]
        ])
        rhs = BlockVector([f_H1.vec, f_L2.vec])

        # --- Block preconditioner ---
        pre = BlockMatrix([
            [self._pre_h1.mat.Inverse(
                freedofs=self._fes_h1.FreeDofs(), inverse="pardiso"),
             None],
            [None,
             self._pre_l2.mat.Inverse(
                 freedofs=self._fes_l2.FreeDofs())]
        ])

        # --- Solve with GMRes ---
        if printrates:
            print(f"  FEMBEM: H1 DOFs={self._fes_h1.ndof}, "
                  f"L2 DOFs={self._fes_l2.ndof}")
        with TaskManager():
            sol = GMRes(A=lhs, b=rhs, pre=pre,
                        tol=1e-8, maxsteps=500, printrates=printrates)

        # --- Extract solution ---
        self.gfu_h1 = GridFunction(self._fes_h1)
        self.gfu_h1.vec[:] = sol[0]  # Hz_scat

        self.gfu_l2 = GridFunction(self._fes_l2)
        self.gfu_l2.vec[:] = sol[1]  # lambda = dHz_scat/dn

    def get_total_field(self):
        """Get total Hz field (scattered + incident).

        For FEM mode: gfu_h1 already contains total field.
        For FEMBEM mode: gfu_h1 contains scattered field, add Hz_inc.

        Returns:
            Hz_total: CoefficientFunction representing total Hz field
        """
        from ngsolve import CF

        if self._mode == 'fem':
            return self.gfu_h1
        else:
            if self._Hz_inc_cf is not None:
                return self.gfu_h1 + self._Hz_inc_cf
            else:
                return self.gfu_h1 + CF(self._Hz_inc)

    def get_scattered_field(self):
        """Get scattered Hz field only.

        For FEM mode: gfu_h1 - Hz_inc.
        For FEMBEM mode: gfu_h1 directly.

        Returns:
            Hz_scat: CoefficientFunction
        """
        from ngsolve import CF

        if self._mode == 'fem':
            if self._Hz_inc_cf is not None:
                return self.gfu_h1 - self._Hz_inc_cf
            else:
                return self.gfu_h1 - CF(self._Hz_inc)
        else:
            return self.gfu_h1

    def get_neumann_data(self):
        """Get Neumann data (dHz/dn on boundary).

        Only available in FEMBEM mode.

        Returns:
            gfu_l2: GridFunction on SurfaceL2 space (lambda = dHz_scat/dn)
        """
        if self._mode != 'fembem':
            raise RuntimeError("Neumann data only available in FEMBEM mode")
        return self.gfu_l2

    def compute_loss(self):
        """Compute eddy current power loss.

        For scalar Hz formulation:
            J = curl H => Jx = dHz/dy, Jy = -dHz/dx
            P = (1/(2*sigma)) * integral(|J|^2) dV
              = (1/(2*sigma)) * integral(|grad_xy Hz|^2) dV

        Returns:
            P_loss: Time-averaged power loss [W]
        """
        from ngsolve import Integrate, Conj, grad, dx

        if not self._solved:
            raise RuntimeError("Call solve() before compute_loss()")

        # grad(Hz_inc) = 0 for uniform field, so grad(total) = grad(scat)
        if self._mode == 'fem':
            u = self.gfu_h1
        else:
            u = self.gfu_h1  # Hz_scat; grad(scat) = grad(total) for uniform inc

        grad_u = grad(u)
        integrand = (1.0 / (2.0 * self.sigma)) * grad_u * Conj(grad_u)
        P_loss = Integrate(integrand, self.mesh).real

        return P_loss

    def compute_stored_energy(self):
        """Compute stored magnetic energy.

        W_m = (mu/2) * integral(|Hz|^2) dV

        Returns:
            W_m: Stored magnetic energy [J]
        """
        from ngsolve import Integrate, Conj

        if not self._solved:
            raise RuntimeError("Call solve() before compute_stored_energy()")

        Hz_total = self.get_total_field()
        W_m = 0.5 * self.mu * Integrate(Hz_total * Conj(Hz_total),
                                         self.mesh).real
        return W_m

    def get_skin_depth(self):
        """Get skin depth at current frequency.

        Returns:
            delta: Skin depth [m]
        """
        return self.delta

    def print_summary(self):
        """Print summary of solver state."""
        print(f"EddyCurrentFEMBEM Summary:")
        print(f"  Mode: {self._mode}")
        print(f"  Frequency: {self.freq:.1f} Hz")
        print(f"  Conductivity: {self.sigma:.2e} S/m")
        print(f"  Permeability: mu_r = {self.mu_r:.1f}")
        print(f"  Skin depth: {self.delta*1e3:.4f} mm")
        print(f"  FE order: {self.order}")
        if self._mode == 'fem':
            print(f"  DOFs: {self._fes.ndof}")
        else:
            print(f"  H1 DOFs (interior):     {self._fes_h1.ndof}")
            print(f"  L2 DOFs (boundary BEM): {self._fes_l2.ndof}")
        print(f"  Assembly time: {self.t_assemble:.3f} s")
        if self._solved:
            print(f"  Solve time: {self.t_solve:.3f} s")
            P = self.compute_loss()
            print(f"  Eddy current loss: {P:.4e} W")
