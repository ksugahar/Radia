"""
ngsbem_eddy.py

Eddy Current Solvers using ngsolve/ngsbem.

Three solver modes:
1. FEM-only: Dirichlet BC from known external field (for validation)
2. FEM-BEM coupled: Costabel symmetric coupling (EddyCurrentFEMBEM)
3. BEM+SIBC: Surface impedance replaces interior FEM (EddyCurrentBEMSIBC)

Scalar Hz formulation:
    Interior (conductor): laplacian(Hz) = j*omega*mu*sigma * Hz  (diffusion)
    Exterior (air):       laplacian(Hz) = 0                       (Laplace, via BEM)

Costabel symmetric coupling (FEM-BEM):

    | A_FEM + D    (-0.5*M + K)^T |   | Hz_scat |   | -k^2 * Hz_inc |
    | (-0.5*M + K)     -V         | * | lambda  | = | 0             |

BEM+SIBC (replaces A_FEM with gamma*M_surf):

    | gamma*M_s + D  (-0.5*M + K)^T |   | Hz_scat |   | -gamma*(Hz_inc,v)_S |
    | (-0.5*M + K)       -V         | * | lambda  | = | 0                   |

    gamma = sqrt(j*omega*mu*sigma) = (1+j)/delta
    SIBC approximates interior DtN map: dHz/dn ~ gamma * Hz

Part of Radia project
"""

import numpy as np
import time
from scipy.linalg import lu_factor, lu_solve

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
                              ds, dx, grad)
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
        from ngsolve import (GridFunction, LinearForm,
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
        # Note: While the block matrix is complex symmetric (where COCG might be
        # mathematically optimal for memory), we use GMRes because COCG often
        # fails to stably converge for these particular BEM-coupled systems without
        # a specialized spectral preconditioner. GMRes guarantees residual minimization.
        if printrates:
            print(f"  FEMBEM: H1 DOFs={self._fes_h1.ndof}, "
                  f"L2 DOFs={self._fes_l2.ndof}")
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
              = (1/(2*sigma)) * integral(|dHz/dx|^2 + |dHz/dy|^2) dV

        Note: We use only the x,y components of grad(Hz), NOT |grad Hz|^2,
        because Jz = 0 always in the scalar Hz formulation. Including
        |dHz/dz|^2 would massively overestimate loss (skin-depth decay
        on top/bottom faces contributes to dHz/dz but NOT to J).

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
        # |J|^2 = |dHz/dx|^2 + |dHz/dy|^2 (NOT |dHz/dz|^2)
        J_sq = (grad_u[0] * Conj(grad_u[0])
                + grad_u[1] * Conj(grad_u[1]))
        integrand = (1.0 / (2.0 * self.sigma)) * J_sq
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


def _extract_dense_complex(mat, ndof):
    """Extract dense complex numpy matrix from NGSolve BaseMatrix.

    Uses COO() extraction for BEM operators (SparseMatrix with 100% fill,
    use_fmm=False), which is ~2500x faster than column-by-column MatVec.
    For FEM sparse matrices, COO also works correctly (sparse -> dense).
    Falls back to MatVec for SumMatrix (use_fmm=True) or other types.
    """
    if hasattr(mat, 'COO'):
        from scipy.sparse import coo_matrix
        rows, cols, vals = mat.COO()
        M = coo_matrix((vals, (rows, cols)),
                        shape=(mat.height, mat.width)).toarray()
        # Restrict to ndof x ndof if matrix is larger (e.g. FEM with extra DOFs)
        M = M[:ndof, :ndof]
        if M.dtype != complex:
            M = M.astype(complex)
        return M

    # Fallback: column-by-column MatVec (slow but universal)
    M = np.zeros((ndof, ndof), dtype=complex)
    for i in range(ndof):
        ei = mat.CreateColVector()
        ei[:] = 0
        ei[i] = 1.0
        col = mat.CreateColVector()
        mat.Mult(ei, col)
        for j in range(ndof):
            M[j, i] = complex(col[j])
    return M


def _extract_square_matrix(operator_mat, n_dof, indices=None):
    """Extract dense matrix from square NGSolve operator.

    Args:
        operator_mat: NGSolve BaseMatrix (square: n_dof x n_dof)
        n_dof: Number of DOFs
        indices: If given, restrict to these DOF indices

    Returns:
        Dense numpy array (complex)
    """
    full = _extract_dense_complex(operator_mat, n_dof)
    if indices is not None:
        return full[np.ix_(indices, indices)]
    return full


def _extract_cross_matrix(operator_mat, n_input, n_output,
                           col_indices=None):
    """Extract dense matrix from non-square NGSolve operator.

    Uses COO() when available for fast extraction.
    Falls back to column-by-column MatVec otherwise.

    Args:
        operator_mat: NGSolve BaseMatrix (n_output x n_input)
        n_input: Trial (column) space DOFs
        n_output: Test (row) space DOFs
        col_indices: If given, restrict columns to these indices

    Returns:
        Dense numpy array (n_output x n_input) or restricted
    """
    if hasattr(operator_mat, 'COO'):
        from scipy.sparse import coo_matrix
        rows, cols, vals = operator_mat.COO()
        M = coo_matrix((vals, (rows, cols)),
                        shape=(operator_mat.height, operator_mat.width)).toarray()
        if M.dtype != complex:
            M = M.astype(complex)
        # Trim to expected dimensions if needed
        M = M[:n_output, :n_input]
        if col_indices is not None:
            return M[:, col_indices]
        return M

    # Fallback: column-by-column MatVec
    rv = operator_mat.CreateRowVector()
    cv = operator_mat.CreateColVector()
    rv_size = len(rv)
    cv_size = len(cv)

    if rv_size == n_input and cv_size == n_output:
        create_input = operator_mat.CreateRowVector
        create_output = operator_mat.CreateColVector
    elif cv_size == n_input and rv_size == n_output:
        create_input = operator_mat.CreateColVector
        create_output = operator_mat.CreateRowVector
    else:
        if rv_size >= n_output:
            create_input = operator_mat.CreateColVector
            create_output = operator_mat.CreateRowVector
        else:
            create_input = operator_mat.CreateRowVector
            create_output = operator_mat.CreateColVector

    M = np.zeros((n_output, n_input), dtype=complex)
    for j in range(n_input):
        x = create_input()
        x[:] = 0
        if j < len(x):
            x[j] = 1.0
        y = create_output()
        operator_mat.Mult(x, y)
        for i in range(min(n_output, len(y))):
            M[i, j] = complex(y[i])

    if col_indices is not None:
        return M[:, col_indices]
    return M


class EddyCurrentBEMSIBC:
    """BEM + SIBC eddy current solver.

    Replaces the interior FEM solve in FEM-BEM with a Surface Impedance
    Boundary Condition (SIBC), eliminating the need for volume PDE solves.

    Formulation (Costabel coupling with SIBC replacing A_FEM):

        | gamma*M_s + D  (-0.5*M + K)^T |   | Hz_scat |   | -gamma*(Hz_inc,v)_S |
        | (-0.5*M + K)       -V         | * | lambda  | = | 0                   |

    where:
        gamma = sqrt(j*omega*mu*sigma) = (1+j)/delta
        M_s = surface mass matrix (H1 trace on boundary)
        D = hypersingular operator (exterior BEM)
        V = single layer potential operator (exterior BEM)
        K = double layer operator (exterior BEM)
        M = H1-SurfaceL2 coupling mass matrix

    SIBC approximation:
        The exact interior DtN map (from solving the diffusion PDE) is
        replaced by the analytical relation: dHz/dn ~ gamma * Hz.
        This is valid when the skin depth delta << object size.

    Advantages over FEM-BEM:
        - BEM operators are frequency-independent (Laplace kernel)
        - Only gamma changes with frequency -> fast frequency sweep
        - No volume PDE assembly/solve per frequency
        - Fewer DOFs when restricted to surface

    Limitations:
        - Requires delta << smallest object dimension
        - Not valid at DC (gamma -> 0, SIBC breaks down)
        - Less accurate than FEM-BEM for thick skin depth

    Usage:
        mesh = create_conductor_mesh(0.01, 0.01, 0.005, maxh=0.003)
        solver = EddyCurrentBEMSIBC(mesh, sigma=5.8e7, order=2)
        solver.assemble(intorder=12)
        solver.solve(freq=1e6, B_ext=[0, 0, 1.0])
        P = solver.compute_loss()

        # Efficient frequency sweep
        freqs = np.logspace(3, 7, 50)
        results = solver.frequency_sweep(freqs, B_ext=[0, 0, 1.0])
    """

    def __init__(self, mesh, sigma=5.8e7, mu_r=1.0, order=2,
                 surface_label="surface"):
        """Initialize BEM+SIBC solver.

        Args:
            mesh: NGSolve Mesh (3D volume mesh of conductor)
            sigma: Electrical conductivity [S/m]
            mu_r: Relative permeability
            order: FE order for H1 and BEM spaces
            surface_label: Boundary label for conductor surface
        """
        self.mesh = mesh
        self.sigma = sigma
        self.mu_r = mu_r
        self.mu = mu_r * MU_0
        self.order = order
        self.surface_label = surface_label

        # State
        self.omega = 0.0
        self.freq = 0.0
        self.delta = float('inf')
        self._assembled = False
        self._solved = False
        self.gfu_h1 = None
        self.gfu_l2 = None
        self._Hz_inc = None
        self._Hz_inc_cf = None

        # Numpy matrices (extracted from BEM operators)
        self._V_np = None       # Single layer (n_l2 x n_l2)
        self._K_np = None       # Double layer (n_l2 x n_surf)
        self._D_np = None       # Hypersingular (n_surf x n_surf)
        self._M_coupling_np = None  # H1-L2 mass (n_l2 x n_surf)
        self._M_surf_np = None  # Surface H1 mass (n_surf x n_surf)
        self._Hz_inc_load = None  # Load vector for unit Hz_inc

        # DOF indices
        self._surf_indices = None  # H1 surface DOF indices
        self._n_surf = 0           # Number of surface H1 DOFs
        self._n_l2 = 0             # Number of SurfaceL2 DOFs

        # Timing
        self.t_assemble = 0.0
        self.t_solve = 0.0

    def assemble(self, intorder=12):
        """Assemble frequency-independent BEM operators.

        Extracts V, K, D, M operators to numpy arrays for efficient
        frequency sweeps. This is the expensive step (done once).

        Args:
            intorder: Integration order for BEM singular quadrature
        """
        from ngsolve import (H1, SurfaceL2, BilinearForm, GridFunction,
                              ds, dx, grad, CF)
        from ngsolve.bem import (SingleLayerPotentialOperator,
                                  DoubleLayerPotentialOperator,
                                  HypersingularOperator)
        try:
            from .ngsbem_peec import extract_dense_matrix
        except ImportError:
            from ngsbem_peec import extract_dense_matrix

        t_start = time.perf_counter()
        label = self.surface_label

        # --- Function spaces ---
        self._fes_h1 = H1(self.mesh, order=self.order, complex=True)
        self._fes_l2 = SurfaceL2(
            self.mesh, order=self.order - 1,
            complex=True, dual_mapping=True,
            definedon=self.mesh.Boundaries(label))

        n_h1 = self._fes_h1.ndof
        self._n_l2 = self._fes_l2.ndof

        # --- Identify surface DOFs in H1 space ---
        surf_dofs_ba = self._fes_h1.GetDofs(self.mesh.Boundaries(label))
        self._surf_indices = [i for i in range(n_h1) if surf_dofs_ba[i]]
        self._n_surf = len(self._surf_indices)

        print(f"  BEM+SIBC: H1 total={n_h1}, surface={self._n_surf}, "
              f"L2={self._n_l2}")

        # --- Assemble BEM operators ---
        print("  Assembling BEM operators...")
        u, v = self._fes_h1.TnT()
        uL2, vL2 = self._fes_l2.TnT()

        V_op = SingleLayerPotentialOperator(
            self._fes_l2, intorder=intorder)
        K_op = DoubleLayerPotentialOperator(
            self._fes_h1, self._fes_l2,
            trial_definedon=self.mesh.Boundaries(label),
            test_definedon=self.mesh.Boundaries(label),
            intorder=intorder)
        D_op = HypersingularOperator(
            self._fes_h1,
            definedon=self.mesh.Boundaries(label),
            intorder=intorder)
        M_bf = BilinearForm(
            self._fes_h1.TrialFunction()
            * self._fes_l2.TestFunction().Trace()
            * ds(label)).Assemble()

        # --- Surface mass matrix (H1 trace) ---
        M_surf_bf = BilinearForm(u * v * ds(label))
        M_surf_bf.Assemble()

        # --- Weighted surface mass matrix for SIBC loss ---
        # SIBC loss: P = Re(Z_s)/2 * integral |n x H|^2 dS
        # For H = (0,0,Hz): |n x H|^2 = |Hz|^2 * (1 - nz^2)
        from ngsolve import specialcf
        n = specialcf.normal(3)
        nz_sq = n[2] * n[2]
        M_loss_bf = BilinearForm(u * v * (1.0 - nz_sq) * ds(label))
        M_loss_bf.Assemble()

        # --- Extract to numpy (restricted to surface DOFs) ---
        print("  Extracting matrices to numpy...")
        si = self._surf_indices

        # V: SurfaceL2 x SurfaceL2 (square, no restriction)
        self._V_np = _extract_dense_complex(V_op.mat, self._n_l2)

        # D: H1 x H1 (square) -> restrict to surface DOFs
        self._D_np = _extract_square_matrix(D_op.mat, n_h1, si)

        # K: H1 -> SurfaceL2 (cross-space, n_l2 x n_h1)
        self._K_np = _extract_cross_matrix(
            K_op.mat, n_h1, self._n_l2, col_indices=si)

        # M_coupling: H1 -> SurfaceL2 (cross-space mass, n_l2 x n_h1)
        self._M_coupling_np = _extract_cross_matrix(
            M_bf.mat, n_h1, self._n_l2, col_indices=si)

        # M_surf: H1 surface mass matrix (square) -> restrict to surface
        self._M_surf_np = _extract_square_matrix(
            M_surf_bf.mat, n_h1, si)

        # M_loss: weighted surface mass matrix for SIBC loss
        self._M_loss_np = _extract_square_matrix(
            M_loss_bf.mat, n_h1, si)

        # --- Precompute load vector for unit Hz_inc ---
        # f[i] = integral(phi_i * 1) dS = sum_j M_surf[i,j]
        # (partition of unity: sum(phi_j) = 1)
        self._Hz_inc_load = np.sum(self._M_surf_np, axis=1)

        # Preconditioner matrices not needed (using numpy direct solve)

        self._assembled = True
        self.t_assemble = time.perf_counter() - t_start
        print(f"  Assembly done in {self.t_assemble:.3f} s")

    def set_source_field(self, Hz_inc_cf):
        """Set custom source field for PEEC coupling.

        Args:
            Hz_inc_cf: ngsolve CoefficientFunction giving Hz_inc(x,y,z)
        """
        self._Hz_inc_cf = Hz_inc_cf

    def solve(self, freq, B_ext=None, printrates=False):
        """Solve at a single frequency using BEM+SIBC.

        Args:
            freq: Frequency [Hz]
            B_ext: External B field [Bx, By, Bz] in Tesla
            printrates: Print solver info

        Returns:
            gfu_h1: GridFunction with Hz_scat on surface
        """
        if not self._assembled:
            raise RuntimeError("Call assemble() first")

        if B_ext is None:
            B_ext = [0, 0, 1.0]
        B_ext = np.asarray(B_ext, dtype=float)

        self.freq = freq
        self.omega = 2.0 * np.pi * freq
        self._Hz_inc = B_ext[2] / MU_0

        t_start = time.perf_counter()

        # --- SIBC parameter ---
        if self.omega > 0 and self.sigma > 0:
            gamma = np.sqrt(1j * self.omega * self.mu * self.sigma)
            self.delta = np.sqrt(2.0 / (self.omega * self.mu * self.sigma))
        else:
            raise ValueError(
                "BEM+SIBC requires freq > 0 and sigma > 0. "
                "Use FEM-BEM for DC problems.")

        ns = self._n_surf
        nl = self._n_l2
        n_total = ns + nl

        # --- Build block system ---
        # Block (1,1): gamma * M_surf + D  (n_surf x n_surf)
        A11 = gamma * self._M_surf_np + self._D_np

        # Block (1,2): (-0.5*M_coupling + K)^T  (n_surf x n_l2)
        offdiag = -0.5 * self._M_coupling_np + self._K_np  # (n_l2 x n_surf)
        A12 = offdiag.T  # (n_surf x n_l2)

        # Block (2,1): (-0.5*M_coupling + K)  (n_l2 x n_surf)
        A21 = offdiag

        # Block (2,2): -V  (n_l2 x n_l2)
        A22 = -self._V_np

        A = np.block([[A11, A12], [A21, A22]])

        # --- RHS ---
        rhs = np.zeros(n_total, dtype=complex)

        # SIBC source: -gamma * (Hz_inc, v)_surface
        if self._Hz_inc_cf is not None:
            # Custom source field: project onto H1 surface DOFs
            from ngsolve import GridFunction, CF
            gfu_tmp = GridFunction(self._fes_h1)
            gfu_tmp.Set(self._Hz_inc_cf)
            Hz_inc_vec = gfu_tmp.vec.FV().NumPy()[self._surf_indices].copy()
            rhs[:ns] = -gamma * self._M_surf_np @ Hz_inc_vec
        else:
            rhs[:ns] = -gamma * self._Hz_inc * self._Hz_inc_load

        # --- Solve ---
        if printrates:
            print(f"  BEM+SIBC: freq={freq:.2e} Hz, "
                  f"gamma={gamma:.4e}, delta={self.delta*1e3:.4f} mm")
            print(f"  System size: {n_total} ({ns} H1_surf + {nl} L2)")
            cond = np.linalg.cond(A)
            print(f"  Condition number: {cond:.2e}")

        sol = np.linalg.solve(A, rhs)

        # --- Extract and store solution ---
        Hz_scat_surf = sol[:ns]
        lambda_vec = sol[ns:]

        # Store Hz_scat as GridFunction (surface DOFs only)
        from ngsolve import GridFunction
        self.gfu_h1 = GridFunction(self._fes_h1)
        self.gfu_h1.vec[:] = 0
        vec_np = self.gfu_h1.vec.FV().NumPy()
        for k, idx in enumerate(self._surf_indices):
            vec_np[idx] = Hz_scat_surf[k]

        self.gfu_l2 = GridFunction(self._fes_l2)
        self.gfu_l2.vec.FV().NumPy()[:] = lambda_vec

        # Store surface solution vector for loss computation
        self._Hz_scat_surf = Hz_scat_surf
        self._gamma = gamma

        self._solved = True
        self.t_solve = time.perf_counter() - t_start
        return self.gfu_h1

    def get_total_field(self):
        """Get total Hz field (scattered + incident).

        Returns:
            Hz_total: CoefficientFunction
        """
        from ngsolve import CF
        if self._Hz_inc_cf is not None:
            return self.gfu_h1 + self._Hz_inc_cf
        else:
            return self.gfu_h1 + CF(self._Hz_inc)

    def get_scattered_field(self):
        """Get scattered Hz field.

        Returns:
            Hz_scat: GridFunction
        """
        return self.gfu_h1

    def get_neumann_data(self):
        """Get Neumann data (dHz/dn on boundary).

        Returns:
            gfu_l2: GridFunction on SurfaceL2
        """
        return self.gfu_l2

    def compute_loss(self):
        """Compute eddy current power loss using SIBC formula.

        SIBC loss for scalar Hz formulation:
            P = (1/2) * Re(Z_s) * integral_S |n x H|^2 dS
        For H = (0,0,Hz): |n x H|^2 = |Hz|^2 * (1 - nz^2)

        This correctly weights the surface integral: side faces (nz=0)
        contribute fully, top/bottom faces (nz=+/-1) contribute zero
        (H is normal to surface there, no tangential component).

        Returns:
            P_loss: Time-averaged power loss [W]
        """
        if not self._solved:
            raise RuntimeError("Call solve() before compute_loss()")

        # Z_s = gamma / sigma
        Z_s = self._gamma / self.sigma

        # Hz_total on surface DOFs
        if self._Hz_inc_cf is not None:
            from ngsolve import GridFunction
            gfu_tmp = GridFunction(self._fes_h1)
            gfu_tmp.Set(self._Hz_inc_cf)
            Hz_inc_surf = gfu_tmp.vec.FV().NumPy()[
                self._surf_indices].copy()
        else:
            Hz_inc_surf = np.full(self._n_surf, self._Hz_inc,
                                  dtype=complex)

        Hz_total_surf = self._Hz_scat_surf + Hz_inc_surf

        # P = (1/2) * Re(Z_s) * <Hz_total, M_loss, Hz_total>
        # M_loss is weighted by (1 - nz^2) for tangential H component
        P = 0.5 * Z_s.real * np.real(
            Hz_total_surf.conj() @ self._M_loss_np @ Hz_total_surf)
        return P

    def compute_stored_energy(self):
        """Compute stored magnetic energy (SIBC approximation).

        W_m ~ (mu * delta / 4) * integral_S |Hz_total|^2 dS

        Returns:
            W_m: Stored magnetic energy [J]
        """
        if not self._solved:
            raise RuntimeError("Call solve() before compute_stored_energy()")

        if self._Hz_inc_cf is not None:
            from ngsolve import GridFunction
            gfu_tmp = GridFunction(self._fes_h1)
            gfu_tmp.Set(self._Hz_inc_cf)
            Hz_inc_surf = gfu_tmp.vec.FV().NumPy()[
                self._surf_indices].copy()
        else:
            Hz_inc_surf = np.full(self._n_surf, self._Hz_inc,
                                  dtype=complex)

        Hz_total_surf = self._Hz_scat_surf + Hz_inc_surf
        quadform = np.real(
            Hz_total_surf.conj() @ self._M_surf_np @ Hz_total_surf)
        W_m = 0.25 * self.mu * self.delta * quadform
        return W_m

    def compute_impedance(self):
        """Get SIBC surface impedance at current frequency.

        Returns:
            Z_s: Complex surface impedance [Ohm] = gamma / sigma
        """
        if not self._solved:
            raise RuntimeError("Call solve() first")
        return self._gamma / self.sigma

    def get_skin_depth(self):
        """Get skin depth at current frequency.

        Returns:
            delta: Skin depth [m]
        """
        return self.delta

    def frequency_sweep(self, freqs, B_ext=None):
        """Efficient multi-frequency sweep.

        BEM operators are assembled once (in assemble()). Only the SIBC
        parameter gamma changes per frequency, so each frequency point
        requires only a matrix solve (no reassembly).

        Args:
            freqs: Array of frequencies [Hz]
            B_ext: External B field [Bx, By, Bz] in Tesla

        Returns:
            dict with:
                'freqs': frequency array [Hz]
                'P_loss': power loss per frequency [W]
                'W_stored': stored energy per frequency [J]
                'Z_s': surface impedance per frequency [Ohm]
                'delta': skin depth per frequency [m]
        """
        if not self._assembled:
            raise RuntimeError("Call assemble() first")

        freqs = np.asarray(freqs, dtype=float)
        n_freq = len(freqs)

        P_loss = np.zeros(n_freq)
        W_stored = np.zeros(n_freq)
        Z_s_arr = np.zeros(n_freq, dtype=complex)
        delta_arr = np.zeros(n_freq)

        for i, f in enumerate(freqs):
            self.solve(f, B_ext)
            P_loss[i] = self.compute_loss()
            W_stored[i] = self.compute_stored_energy()
            Z_s_arr[i] = self.compute_impedance()
            delta_arr[i] = self.delta

        return {
            'freqs': freqs,
            'P_loss': P_loss,
            'W_stored': W_stored,
            'Z_s': Z_s_arr,
            'delta': delta_arr,
        }

    def print_summary(self):
        """Print summary of solver state."""
        print("EddyCurrentBEMSIBC Summary:")
        print(f"  Conductivity: {self.sigma:.2e} S/m")
        print(f"  Permeability: mu_r = {self.mu_r:.1f}")
        print(f"  FE order: {self.order}")
        print(f"  Surface H1 DOFs: {self._n_surf}")
        print(f"  SurfaceL2 DOFs:  {self._n_l2}")
        print(f"  Total system:    {self._n_surf + self._n_l2}")
        print(f"  Assembly time: {self.t_assemble:.3f} s")
        if self._solved:
            print(f"  Frequency: {self.freq:.1f} Hz")
            print(f"  Skin depth: {self.delta*1e3:.4f} mm")
            print(f"  Surface impedance: {self.compute_impedance():.4e} Ohm")
            print(f"  Solve time: {self.t_solve:.3f} s")
            P = self.compute_loss()
            print(f"  Eddy current loss: {P:.4e} W")


# ======================================================================
# Loop-only BEM+SIBC for conducting shields (vector formulation)
# ======================================================================

class LoopBasisBuilder:
    """Build divergence-free basis from surface mesh topology.

    Uses NGSolve's div operator to build the face-edge divergence matrix
    D (n_faces x n_edges), ensuring consistency with the HDivSurface
    basis function convention used by M and V operators.

    The divergence-free subspace is null(D). For a closed surface:
      n_loops = n_edges - rank(D) = V - 1  (for genus-0)

    IMPORTANT: The graph-theoretic cycle space (null of edge-vertex
    incidence) has dimension T-1, which is LARGER than V-1. Those extra
    vectors are NOT divergence-free in the BEM sense and would cause
    the V_2 term of the Maxwell SLP to leak into the loop subspace.
    """

    def __init__(self, mesh, fes):
        """Build divergence-free basis from mesh and HDivSurface FE space.

        Uses NGSolve's own div operator to build D, ensuring the loop
        basis T = null(D) is in the same convention as the BEM operators.

        Args:
            mesh: NGSolve Mesh (with volume elements; boundary = closed surface)
            fes: HDivSurface FE space (order=0)
        """
        from scipy.linalg import null_space
        from ngsolve import BND, SurfaceL2, BilinearForm, div, ds

        # Collect active (boundary) DOFs
        active_dof_set = set()
        for el in mesh.Elements(BND):
            for d_raw in fes.GetDofNrs(el):
                d = d_raw if d_raw >= 0 else ~d_raw
                active_dof_set.add(d)

        self.active_dofs = sorted(active_dof_set)
        self.n_active = len(self.active_dofs)
        self._dof_to_idx = {d: i for i, d in enumerate(self.active_dofs)}

        bnd_elements = list(mesh.Elements(BND))
        n_faces = len(bnd_elements)
        self.n_faces = n_faces

        # Store face geometry for RHS computation
        self._face_verts = []     # list of vertex coords per face
        self._face_centroids = [] # centroids
        self._face_areas = []     # areas

        # DOF info per face (for sign extraction from D later)
        face_dof_indices = []  # list of [(idx, dof_nr), ...] per face

        for face_idx, el in enumerate(bnd_elements):
            verts = [mesh.vertices[v.nr].point for v in el.vertices]
            coords = [np.array([v[0], v[1], v[2]]) for v in verts]

            self._face_verts.append(coords)
            centroid = sum(coords) / len(coords)
            self._face_centroids.append(centroid)

            if len(coords) == 3:
                area = 0.5 * np.linalg.norm(
                    np.cross(coords[1] - coords[0], coords[2] - coords[0]))
            else:
                area = (0.5 * np.linalg.norm(
                    np.cross(coords[1] - coords[0], coords[2] - coords[0]))
                    + 0.5 * np.linalg.norm(
                    np.cross(coords[2] - coords[0], coords[3] - coords[0])))
            self._face_areas.append(area)

            dofs_on_face = []
            for d_raw in fes.GetDofNrs(el):
                d = d_raw if d_raw >= 0 else ~d_raw
                if d in self._dof_to_idx:
                    dofs_on_face.append((self._dof_to_idx[d], d))
            face_dof_indices.append(dofs_on_face)

        # -----------------------------------------------------------------
        # Build D using NGSolve's div operator for consistent convention.
        #
        # The geometric sign computation (cross-product based) can differ
        # from NGSolve's internal RT0 convention for ~50% of edges.
        # Using NGSolve's own div operator ensures D is in the SAME
        # convention as M (mass) and V (SLP), so the loop projection
        # M_LL = T^T M T is consistent.
        # -----------------------------------------------------------------
        fes_l2 = SurfaceL2(mesh, order=0)

        # Build mapping: SurfaceL2 DOF index -> BND element index
        l2_dof_to_face = {}
        for fi, el in enumerate(bnd_elements):
            for d in fes_l2.GetDofNrs(el):
                d_abs = d if d >= 0 else ~d
                l2_dof_to_face[d_abs] = fi

        # Mixed bilinear form: <div(u_hdiv), v_l2> on boundary surface
        u_trial = fes.TrialFunction()
        v_test = fes_l2.TestFunction()

        bf_div = BilinearForm(trialspace=fes, testspace=fes_l2)
        bf_div += div(u_trial.Trace()) * v_test * ds(bonus_intorder=2)
        bf_div.Assemble()

        D_mat = bf_div.mat

        # Use GridFunction vectors (not CreateColVector which may have
        # incompatible vector types for mixed HDivSurface/SurfaceL2 forms)
        from ngsolve import GridFunction as GF
        gf_trial = GF(fes)
        gf_test = GF(fes_l2)

        # Extract D column by column for active DOFs only
        D_raw = np.zeros((fes_l2.ndof, self.n_active))
        for k, dk in enumerate(self.active_dofs):
            gf_trial.vec[:] = 0
            gf_trial.vec[dk] = 1.0
            gf_test.vec.data = D_mat * gf_trial.vec
            for f in range(fes_l2.ndof):
                D_raw[f, k] = gf_test.vec[f]

        # Reorder rows: L2 DOF order -> BND element order
        D = np.zeros((n_faces, self.n_active))
        for l2_dof in range(fes_l2.ndof):
            face_idx = l2_dof_to_face.get(l2_dof, l2_dof)
            D[face_idx, :] = D_raw[l2_dof, :]

        self._D = D

        # Build face_dofs with correct signs from D (NGSolve convention).
        # For RT0 order=0: D[f,k] = +-1 (sign of basis function on face f).
        # This sign is used in _compute_excitation for centroid quadrature:
        #   phi_k(centroid) = sign * (centroid - p_opp) / (2 * area)
        self._face_dofs = []
        for face_idx in range(n_faces):
            face_dof_list = []
            for dof_idx, dof_nr in face_dof_indices[face_idx]:
                d_val = D[face_idx, dof_idx]
                sign = np.sign(d_val) if abs(d_val) > 1e-10 else 1.0
                face_dof_list.append((dof_idx, sign, dof_nr))
            self._face_dofs.append(face_dof_list)

        # Loop basis = null space of D (divergence-free subspace)
        self.T_loop = null_space(D)
        self.n_loops = self.T_loop.shape[1]

        # Verify divergence-free condition
        residual = np.linalg.norm(D @ self.T_loop)
        if residual > 1e-10:
            raise RuntimeError(
                f"Div-free basis verification failed: ||D * T|| = {residual:.2e}")

    def project(self, M_full):
        """Project a full-DOF matrix into loop subspace.

        Args:
            M_full: (n_active x n_active) matrix (real or complex)

        Returns:
            M_LL: (n_loops x n_loops) projected matrix
        """
        if np.iscomplexobj(M_full):
            return self.T_loop.conj().T @ M_full @ self.T_loop
        return self.T_loop.T @ M_full @ self.T_loop

    def project_rhs(self, rhs_full):
        """Project a full-DOF RHS vector into loop subspace.

        Args:
            rhs_full: (n_active,) vector

        Returns:
            rhs_loop: (n_loops,) projected vector
        """
        if np.iscomplexobj(rhs_full):
            return self.T_loop.conj().T @ rhs_full
        return self.T_loop.T @ rhs_full

    def expand(self, J_loop):
        """Expand loop coefficients back to full DOF space.

        Args:
            J_loop: (n_loops,) loop coefficients

        Returns:
            J_full: (n_active,) full DOF coefficients
        """
        return self.T_loop @ J_loop


class ShieldBEMSIBC:
    """Loop-only BEM+SIBC solver for conducting shields.

    Solves eddy currents on a closed conducting surface using:
      (Zs * M_LL + jw * mu_0 * V_LL) * I_loop = -jw * b_loop

    where:
      V_LL: Vectorial SLP projected to div-free loop subspace
            (from Maxwell SLP with kappa=1; V_2 term vanishes)
      M_LL: surface mass matrix projected to loop subspace
      Zs:   surface impedance = (1+j) / (sigma * delta)
      I_loop: loop current coefficients (div-free surface currents)
      b_loop: excitation from incident vector potential A_inc

    The loop basis is the null space of the FACE-EDGE divergence
    matrix (not the edge-vertex incidence). This ensures div(J)=0
    in the BEM/RT0 sense, eliminating the 1/kappa^2 breakdown.

    Key features:
      - V and M assembled once (frequency-independent)
      - Only Zs and omega change per frequency -> fast sweep
      - Couples with PEEC coils via Biot-Savart excitation
    """

    def __init__(self, mesh, sigma, mu_r=1.0):
        """Initialize shield BEM+SIBC solver.

        Args:
            mesh: NGSolve Mesh (must have volume elements for surface extraction)
            sigma: Conductivity [S/m]
            mu_r: Relative permeability (default 1.0 for non-magnetic)
        """
        self.mesh = mesh
        self.sigma = sigma
        self.mu_r = mu_r
        self.mu = mu_r * MU_0

        self._assembled = False
        self._solved = False
        self.t_assemble = 0.0
        self.t_solve = 0.0

    def assemble(self, intorder=4):
        """Assemble BEM operators and build loop basis.

        This is the expensive step (done once). After assembly,
        frequency sweeps only require matrix solves.

        Args:
            intorder: Integration order for BEM quadrature
        """
        from ngsolve import HDivSurface, BilinearForm, InnerProduct, ds
        from ngsolve.bem import MaxwellSingleLayerPotentialOperator

        t0 = time.time()

        # HDivSurface FE space (order=0: one DOF per boundary edge)
        self._fes = HDivSurface(self.mesh, order=0)
        self._n_dof = self._fes.ndof
        u, v = self._fes.TnT()

        # Build loop basis from mesh topology
        self._loop = LoopBasisBuilder(self.mesh, self._fes)
        n_active = self._loop.n_active
        n_loops = self._loop.n_loops
        active_dofs = self._loop.active_dofs

        # Maxwell SLP with kappa=1 (arbitrary; V_2 term vanishes in loop subspace)
        V_op = MaxwellSingleLayerPotentialOperator(
            self._fes, kappa=1.0, intorder=intorder)
        V_ngmat = V_op.mat

        # Surface mass matrix
        M_bf = BilinearForm(self._fes)
        M_bf += InnerProduct(u.Trace(), v.Trace()) * ds(bonus_intorder=2)
        M_bf.Assemble()
        M_ngmat = M_bf.mat

        # Extract to dense numpy (active DOFs only)
        # IMPORTANT: BEM and FEM operators use DIFFERENT vector types.
        # Using BEM vectors (V_ngmat.CreateColVector()) with BilinearForm
        # gives WRONG results. Must use separate vectors for each operator.
        V_full = np.zeros((n_active, n_active), dtype=complex)
        M_full = np.zeros((n_active, n_active))

        # BEM operator vectors (for Maxwell SLP)
        x_v = V_ngmat.CreateColVector()
        y_v = V_ngmat.CreateColVector()

        for i, di in enumerate(active_dofs):
            x_v[:] = 0
            x_v[di] = 1.0
            V_ngmat.Mult(x_v, y_v)
            for j, dj in enumerate(active_dofs):
                V_full[j, i] = y_v[dj]

        # FEM operator vectors (for mass matrix) - MUST be separate
        x_m = M_ngmat.CreateColVector()
        y_m = M_ngmat.CreateColVector()

        for i, di in enumerate(active_dofs):
            x_m[:] = 0
            x_m[di] = 1.0
            M_ngmat.Mult(x_m, y_m)
            for j, dj in enumerate(active_dofs):
                M_full[j, i] = y_m[dj]

        # Project to loop subspace (V_2 term vanishes here!)
        self._V_LL = self._loop.project(V_full)
        self._M_LL = self._loop.project(M_full)

        # Store for excitation computation
        self._V_full = V_full
        self._M_full = M_full

        self._assembled = True
        self.t_assemble = time.time() - t0

    def solve(self, freq, A_inc_func=None):
        """Solve for eddy currents at given frequency.

        System: (Zs * M_LL + jw * mu_0 * V_LL) * I_loop = -jw * b_loop

        Args:
            freq: Frequency [Hz]
            A_inc_func: Callable(points) -> (n_points, 3) vector potential
                from external source (e.g., PEEC coil Biot-Savart).
                If None, uses uniform unit excitation for testing.

        Returns:
            dict with solution info
        """
        if not self._assembled:
            raise RuntimeError("Call assemble() first")

        t0 = time.time()

        self.freq = freq
        omega = 2.0 * np.pi * freq
        delta = np.sqrt(2.0 / (omega * self.mu * self.sigma))
        gamma = (1.0 + 1j) / delta
        Zs = gamma / self.sigma
        self.delta = delta
        self._Zs = Zs
        self._omega = omega

        n_loops = self._loop.n_loops

        # System matrix from EFIE:
        #   E_inc_tan + E_scat_tan = Zs * J   (SIBC boundary condition)
        #   E_scat = -jw*mu_0*V*J              (SLP, loop subspace: no grad Phi)
        #   => (Zs*M + jw*mu_0*V)*J = -jw*<A_inc, phi>
        #
        # PLUS sign on V: Lenz's law requires A_scat opposing A_inc
        A_sys = Zs * self._M_LL + 1j * omega * MU_0 * self._V_LL

        # RHS: -jw * M * A_inc projected to loop space
        # For PEEC coupling: A_inc = vector potential from coil at edge centers
        if A_inc_func is not None:
            rhs_full = self._compute_excitation(A_inc_func, omega)
            rhs_loop = self._loop.project_rhs(rhs_full)
        else:
            # Uniform test excitation
            rhs_loop = np.ones(n_loops, dtype=complex)

        # Solve
        self._J_loop = np.linalg.solve(A_sys, rhs_loop)
        self._J_full = self._loop.expand(self._J_loop)

        self._solved = True
        self.t_solve = time.time() - t0

        return {
            'freq': freq,
            'delta': delta,
            'Zs': Zs,
            'n_loops': n_loops,
            'J_loop_max': np.abs(self._J_loop).max(),
            't_solve': self.t_solve,
        }

    def _compute_excitation(self, A_inc_func, omega):
        """Compute excitation RHS from incident vector potential.

        For EFIE: E_inc = -jw * A_inc on surface
        RHS[i] = -jw * integral_S A_inc(r) . phi_i(r) dS

        Uses centroid quadrature on each triangle. For RT0 order=0,
        the basis function on triangle T opposite vertex p_k is:
            phi_i(r) = sign * (r - p_k) / (2 * area_T)

        So: integral_T A_inc . phi_i dS ~ sign * A_inc(c) . (c - p_k) / 2

        Args:
            A_inc_func: Callable(points_array) -> (n_points, 3) vector potential
            omega: Angular frequency

        Returns:
            rhs_full: (n_active,) excitation vector
        """
        loop = self._loop
        n_active = loop.n_active

        # Evaluate A_inc at all face centroids
        centroids = np.array(loop._face_centroids)
        A_at_centroids = A_inc_func(centroids)  # (n_faces, 3)

        # Build edge-to-vertices map for finding opposite vertex
        edge_vertex_map = {}  # dof_nr -> set of vertex nrs
        for edge in self.mesh.edges:
            if edge.nr in loop._dof_to_idx:
                edge_vertex_map[edge.nr] = {v.nr for v in edge.vertices}

        # Build face vertex number lists
        from ngsolve import BND
        face_vert_nrs = []
        for el in self.mesh.Elements(BND):
            face_vert_nrs.append([v.nr for v in el.vertices])

        # Accumulate RHS using weak form on triangles
        rhs_full = np.zeros(n_active, dtype=complex)

        for face_idx in range(loop.n_faces):
            A_c = A_at_centroids[face_idx]
            centroid = centroids[face_idx]
            verts = loop._face_verts[face_idx]
            vnrs = face_vert_nrs[face_idx]

            for dof_idx, sign, dof_nr in loop._face_dofs[face_idx]:
                # Find opposite vertex: the face vertex NOT on this edge
                edge_vnrs = edge_vertex_map.get(dof_nr, set())
                p_opp = None
                for k, vnr in enumerate(vnrs):
                    if vnr not in edge_vnrs:
                        p_opp = verts[k]
                        break

                if p_opp is None:
                    continue  # degenerate case

                # RT0 basis at centroid: phi(c) = sign * (c - p_opp) / (2*area)
                # integral_T A . phi dS ~ A(c) . phi(c) * area
                #                       = sign * A(c) . (c - p_opp) / 2
                integral = sign * np.dot(A_c, centroid - p_opp) / 2.0
                rhs_full[dof_idx] += -1j * omega * integral

        return rhs_full

    def compute_loss(self):
        """Compute eddy current loss in the shield.

        P = Re(Zs)/2 * integral |J_t|^2 dS
          = Re(Zs)/2 * J^T * M * J  (in full DOF space)

        Returns:
            P_loss: Power loss [W]
        """
        if not self._solved:
            raise RuntimeError("Call solve() first")

        Zs = self._Zs
        J = self._J_full

        # P = Re(Zs)/2 * J^H * M_full * J
        MJ = self._M_full @ J
        P = 0.5 * Zs.real * np.real(J.conj() @ MJ)
        return P

    def compute_stored_energy(self):
        """Compute stored magnetic energy from shield currents.

        W = mu_0/(2*omega) * Im(J^H * V * J)  (in loop subspace)

        Returns:
            W_stored: Stored energy [J]
        """
        if not self._solved:
            raise RuntimeError("Call solve() first")

        J = self._J_loop
        VJ = self._V_LL @ J
        W = MU_0 / (2.0 * self._omega) * np.abs(np.imag(J.conj() @ VJ))
        return W

    def compute_A_scattered(self, obs_points):
        """Compute vector potential from shield eddy currents at external points.

        A_scat(r) = mu_0/(4*pi) * integral_S J(r')/|r-r'| dS'

        Uses NGSolve Integrate with the solved HDivSurface current distribution.
        Same approach as compute_vector_potential() in ngsbem_peec.py.

        Must call solve() first.

        Args:
            obs_points: (N, 3) array of observation points [m]

        Returns:
            A_scat: (N, 3) complex array of vector potential [T*m]
        """
        from ngsolve import GridFunction, x, y, z, sqrt, Integrate

        if not self._solved:
            raise RuntimeError("Call solve() first")

        obs_points = np.atleast_2d(obs_points)
        n_obs = len(obs_points)
        A = np.zeros((n_obs, 3), dtype=complex)

        bnd = self.mesh.Boundaries("surface")
        fes = self._fes

        active_dofs = self._loop.active_dofs

        for part_label, get_part in [('real', np.real), ('imag', np.imag)]:
            gf_J = GridFunction(fes)
            for k in range(len(self._J_full)):
                gf_J.vec[active_dofs[k]] = float(get_part(self._J_full[k]))

            for i in range(n_obs):
                rx, ry, rz = obs_points[i]
                dist = sqrt((x - rx)**2 + (y - ry)**2 + (z - rz)**2 + 1e-30)
                kernel = 1.0 / dist
                for comp in range(3):
                    val = MU_0 / (4 * np.pi) * Integrate(
                        gf_J[comp] * kernel, self.mesh,
                        definedon=bnd, order=5)
                    if part_label == 'real':
                        A[i, comp] += val
                    else:
                        A[i, comp] += 1j * val

        return A

    def _precompute_excitation_geometry(self, topology_dict):
        """Precompute geometry-only part of excitation RHS for all filaments.

        The excitation for filament j at frequency f is:
            rhs_full_j = -jw * rhs_geom_j

        where rhs_geom_j depends only on geometry (Biot-Savart + RT0 quadrature).
        This is computed once and reused for all frequencies.

        Returns:
            rhs_geom: (n_active, n_seg) real/complex array (geometry-only)
        """
        centers = np.asarray(topology_dict['segment_centers'])
        directions = np.asarray(topology_dict['segment_directions'])
        lengths = np.asarray(topology_dict['segment_lengths'])
        n_seg = len(centers)

        loop = self._loop
        n_active = loop.n_active
        centroids = np.array(loop._face_centroids)

        # Build edge-to-vertices map (for finding opposite vertex)
        edge_vertex_map = {}
        for edge in self.mesh.edges:
            if edge.nr in loop._dof_to_idx:
                edge_vertex_map[edge.nr] = {v.nr for v in edge.vertices}

        from ngsolve import BND
        face_vert_nrs = []
        for el in self.mesh.Elements(BND):
            face_vert_nrs.append([v.nr for v in el.vertices])

        rhs_geom = np.zeros((n_active, n_seg), dtype=float)

        for j in range(n_seg):
            # Biot-Savart A at all face centroids for filament j
            A_at_centroids = _biot_savart_A(
                centroids, centers[j], directions[j], lengths[j])

            # RT0 centroid quadrature (same as _compute_excitation but w/o -jw)
            for face_idx in range(loop.n_faces):
                A_c = A_at_centroids[face_idx]
                centroid = centroids[face_idx]
                verts = loop._face_verts[face_idx]
                vnrs = face_vert_nrs[face_idx]

                for dof_idx, sign, dof_nr in loop._face_dofs[face_idx]:
                    edge_vnrs = edge_vertex_map.get(dof_nr, set())
                    p_opp = None
                    for k, vnr in enumerate(vnrs):
                        if vnr not in edge_vnrs:
                            p_opp = verts[k]
                            break
                    if p_opp is None:
                        continue
                    integral = sign * np.dot(A_c, centroid - p_opp) / 2.0
                    rhs_geom[dof_idx, j] += integral

        return rhs_geom

    def _build_A_scattered_operator(self, obs_points):
        """Build coupling matrix from shield DOFs to A_scattered at obs points.

        C[obs_i, comp, dof_k] such that:
            A_scat[obs_i, comp] = sum_k C[obs_i, comp, k] * J_full[k]

        Uses 7-point triangle quadrature (degree 5) for accuracy.

        Args:
            obs_points: (n_obs, 3) array

        Returns:
            C: (n_obs*3, n_active) real coupling matrix
        """
        loop = self._loop
        n_active = loop.n_active
        obs_points = np.atleast_2d(obs_points)
        n_obs = len(obs_points)

        # 7-point quadrature on reference triangle (degree 5)
        # Barycentric coordinates (L1, L2, L3) and weights
        a1, b1 = 0.059715871789770, 0.470142064105115
        a2, b2 = 0.797426985353087, 0.101286507323456
        quad_bary = np.array([
            [1.0/3, 1.0/3, 1.0/3],
            [a1, b1, b1], [b1, a1, b1], [b1, b1, a1],
            [a2, b2, b2], [b2, a2, b2], [b2, b2, a2],
        ])
        quad_weights = np.array([
            0.225,
            0.132394152788506, 0.132394152788506, 0.132394152788506,
            0.125939180544827, 0.125939180544827, 0.125939180544827,
        ])
        n_quad = len(quad_weights)

        # Build edge-to-vertices map
        edge_vertex_map = {}
        for edge in self.mesh.edges:
            if edge.nr in loop._dof_to_idx:
                edge_vertex_map[edge.nr] = {v.nr for v in edge.vertices}

        from ngsolve import BND
        face_vert_nrs = []
        for el in self.mesh.Elements(BND):
            face_vert_nrs.append([v.nr for v in el.vertices])

        # C matrix: (n_obs, 3, n_active)
        C = np.zeros((n_obs, 3, n_active))

        for face_idx in range(loop.n_faces):
            verts = loop._face_verts[face_idx]
            area = loop._face_areas[face_idx]
            vnrs = face_vert_nrs[face_idx]

            if len(verts) != 3:
                continue  # skip non-triangles (shouldn't happen for BEM)

            v0, v1, v2 = np.array(verts[0]), np.array(verts[1]), np.array(verts[2])

            # Physical quadrature points
            r_q = (quad_bary[:, 0:1] * v0[np.newaxis, :]
                   + quad_bary[:, 1:2] * v1[np.newaxis, :]
                   + quad_bary[:, 2:3] * v2[np.newaxis, :])  # (n_quad, 3)

            # Distance from obs points to quad points: (n_obs, n_quad)
            diff = obs_points[:, np.newaxis, :] - r_q[np.newaxis, :, :]  # (n_obs, n_quad, 3)
            dist = np.sqrt(np.sum(diff**2, axis=2) + 1e-30)  # (n_obs, n_quad)

            # kernel = w_q / |r - r'| : (n_obs, n_quad)
            kernel = quad_weights[np.newaxis, :] / dist  # (n_obs, n_quad)

            for dof_idx, sign, dof_nr in loop._face_dofs[face_idx]:
                # Find opposite vertex
                edge_vnrs = edge_vertex_map.get(dof_nr, set())
                p_opp = None
                for k, vnr in enumerate(vnrs):
                    if vnr not in edge_vnrs:
                        p_opp = np.array(verts[k])
                        break
                if p_opp is None:
                    continue

                # phi_k(r_q) = sign * (r_q - p_opp) / (2*area)  : (n_quad, 3)
                phi_q = sign * (r_q - p_opp[np.newaxis, :]) / (2.0 * area)

                # Dunavant quadrature: integral_T f dS = area * sum_q w_q * f(r_q)
                # (weights sum to 1.0, not 0.5)
                # kernel is (n_obs, n_quad), phi_q is (n_quad, 3)
                contrib = area * np.einsum('oq,qc->oc', kernel, phi_q)
                C[:, :, dof_idx] += MU_0 / (4.0 * np.pi) * contrib

        # Reshape to (n_obs*3, n_active) for efficient matrix multiply
        return C.reshape(n_obs * 3, n_active)

    def compute_impedance_matrix(self, freq, topology_dict):
        """Compute shield contribution to PEEC impedance matrix.

        Optimized version using precomputed operators:
          1. Excitation geometry vectors (computed once per topology)
          2. A_scattered coupling matrix (computed once per obs points)
          3. LU factorization (once per frequency)
          4. Batch solve all filament RHS at once

        Args:
            freq: Frequency [Hz]
            topology_dict: Dict with keys:
                segment_centers, segment_directions, segment_lengths

        Returns:
            Delta_Z: (n_seg x n_seg) complex impedance matrix [Ohm]
        """
        if not self._assembled:
            raise RuntimeError("Call assemble() first")

        omega = 2.0 * np.pi * freq

        centers = np.asarray(topology_dict['segment_centers'])
        directions = np.asarray(topology_dict['segment_directions'])
        lengths = np.asarray(topology_dict['segment_lengths'])
        n_seg = len(centers)

        loop = self._loop

        # --- Precompute (cached, done once for all frequencies) ---

        # Excitation geometry: only depends on filament geometry
        if not hasattr(self, '_cached_rhs_geom'):
            t0 = time.time()
            self._cached_rhs_geom = self._precompute_excitation_geometry(
                topology_dict)
            self._t_precompute_rhs = time.time() - t0
            print(f"  [ShieldBEMSIBC] Precomputed excitation geometry "
                  f"({n_seg} filaments): {self._t_precompute_rhs:.3f} s")

        # A_scattered coupling: only depends on obs point geometry
        obs_key = id(topology_dict)
        if not hasattr(self, '_cached_C_ascat') or self._cached_obs_key != obs_key:
            t0 = time.time()
            self._cached_C_ascat = self._build_A_scattered_operator(centers)
            self._cached_obs_key = obs_key
            self._t_precompute_ascat = time.time() - t0
            print(f"  [ShieldBEMSIBC] Precomputed A_scattered operator "
                  f"({n_seg} obs pts x {loop.n_active} DOFs): "
                  f"{self._t_precompute_ascat:.3f} s")

        rhs_geom = self._cached_rhs_geom  # (n_active, n_seg)
        C_ascat = self._cached_C_ascat     # (n_obs*3, n_active)

        # --- Per-frequency computation ---

        # Surface impedance
        delta = np.sqrt(2.0 / (omega * self.mu * self.sigma))
        gamma = (1.0 + 1j) / delta
        Zs = gamma / self.sigma

        # Slab impedance if thickness is set
        if hasattr(self, 'thickness') and self.thickness is not None:
            t = self.thickness
            gt = gamma * t
            if abs(gt) > 30:
                Zs_slab = Zs
            else:
                Zs_slab = Zs * np.cosh(gt) / np.sinh(gt)
            Zs = Zs_slab

        # System matrix (same for all filaments)
        A_sys = Zs * self._M_LL + 1j * omega * MU_0 * self._V_LL

        # Scale RHS by -jw and project to loop subspace
        rhs_full_all = -1j * omega * rhs_geom  # (n_active, n_seg)
        rhs_loop_all = loop.T_loop.conj().T @ rhs_full_all  # (n_loops, n_seg)

        # Batch solve: one LU factorization, n_seg forward/back substitutions
        lu, piv = lu_factor(A_sys)
        J_loop_all = lu_solve((lu, piv), rhs_loop_all)  # (n_loops, n_seg)

        # Expand to full DOF space
        J_full_all = loop.T_loop @ J_loop_all  # (n_active, n_seg)

        # A_scattered via precomputed coupling matrix
        # C_ascat is (n_seg*3, n_active), J_full_all is (n_active, n_seg)
        A_scat_flat = C_ascat @ J_full_all  # (n_seg*3, n_seg)
        A_scat_all = A_scat_flat.reshape(n_seg, 3, n_seg)  # (n_obs, 3, n_src)

        # Delta_Z from Faraday's law
        Delta_Z = np.zeros((n_seg, n_seg), dtype=complex)
        for i in range(n_seg):
            for j in range(n_seg):
                Delta_Z[i, j] = (1j * omega
                                 * np.dot(A_scat_all[i, :, j], directions[i])
                                 * lengths[i])

        return Delta_Z

    def frequency_sweep(self, freqs, A_inc_func=None):
        """Multi-frequency sweep (fast: only Zs, omega change).

        Args:
            freqs: Array of frequencies [Hz]
            A_inc_func: Excitation function (same for all frequencies)

        Returns:
            dict with arrays: freqs, P_loss, W_stored, delta, Zs
        """
        if not self._assembled:
            raise RuntimeError("Call assemble() first")

        freqs = np.asarray(freqs, dtype=float)
        n_freq = len(freqs)

        P_loss = np.zeros(n_freq)
        W_stored = np.zeros(n_freq)
        delta_arr = np.zeros(n_freq)
        Zs_arr = np.zeros(n_freq, dtype=complex)

        t0 = time.time()
        for i, f in enumerate(freqs):
            self.solve(f, A_inc_func)
            P_loss[i] = self.compute_loss()
            W_stored[i] = self.compute_stored_energy()
            delta_arr[i] = self.delta
            Zs_arr[i] = self._Zs
        t_total = time.time() - t0

        return {
            'freqs': freqs,
            'P_loss': P_loss,
            'W_stored': W_stored,
            'delta': delta_arr,
            'Zs': Zs_arr,
            't_total': t_total,
        }

    def print_summary(self):
        """Print solver summary."""
        print("ShieldBEMSIBC Summary:")
        print(f"  Conductivity: {self.sigma:.2e} S/m")
        print(f"  mu_r: {self.mu_r:.1f}")
        print(f"  Total HDivSurface DOFs: {self._n_dof}")
        print(f"  Active boundary DOFs: {self._loop.n_active}")
        print(f"  Loop DOFs (reduced): {self._loop.n_loops}")
        print(f"  Boundary faces: {self._loop.n_faces}")
        print(f"  Assembly time: {self.t_assemble:.3f} s")
        if self._solved:
            print(f"  Frequency: {self.freq:.1f} Hz")
            print(f"  Skin depth: {self.delta*1e3:.4f} mm")
            print(f"  |Zs| = {abs(self._Zs):.4e} Ohm")
            print(f"  Solve time: {self.t_solve:.6f} s")


def _biot_savart_A(points, center, direction, length):
    """Vector potential A from a finite current filament (unit current)."""
    from radia.biot_savart import a_filament
    points = np.atleast_2d(points)
    p1 = np.array(center) - 0.5 * length * np.array(direction)
    p2 = np.array(center) + 0.5 * length * np.array(direction)
    return np.array([a_filament(p1, p2, pt) for pt in points])


# ====================================================================
# Vector Eddy Current FEM-BEM Solver
# ====================================================================

class VectorEddyCurrentFEMBEM:
    """Vector eddy current FEM-BEM solver using H(curl) + Maxwell SLP.

    3D A-formulation:
        Interior (conductor):  curl(1/mu * curl A) + j*omega*sigma*A = 0
        Exterior (air):        Handled by Maxwell SLP (Weggler stabilized)
        Coupling:              Tangential trace from H(curl) to HDivSurface

    Johnson-Nedelec type coupling with Weggler stabilization.
    Uses only Maxwell SLP (available in bundled ngsolve.bem), avoiding
    the need for Maxwell DLP or Hypersingular operators.

    3x3 block system (Johnson-Nedelec coupling):
        | a_FEM    +B^T        0        | | A   |   | f_1 |
        | B        -A_kappa   -Q_kappa^T | | j   | = | f_2 |
        | 0        -Q_kappa   -k^2*V_k  | | rho |   | f_3 |

    where:
        a_FEM = (1/mu * curl u, curl v) + j*omega*sigma*(u,v)
        B     = tangential trace coupling H(curl) -> HDivSurface
        A_k   = HelmholtzSL on HDivSurface (vectorial part)
        V_k   = HelmholtzSL on SurfaceL2 (scalar part)
        Q_k   = HelmholtzSL on div(HDivSurface) x SurfaceL2

    Advantages over scalar Hz formulation:
        - Full 3D vector eddy currents (not just Hz component)
        - Handles arbitrary mu_r (not limited to mu_r=1)
        - Works for thick skin depth (no SIBC approximation)
        - Weggler stabilization prevents 1/kappa^2 breakdown

    Part of Radia project
    """

    def __init__(self, mesh, sigma, mu_r=1.0, order=1,
                 conductor_label="conductor", surface_label="surface"):
        """Initialize vector FEM-BEM eddy current solver.

        Args:
            mesh: NGSolve Mesh (3D volume mesh of conductor)
            sigma: Electrical conductivity [S/m]
            mu_r: Relative permeability (any value, not restricted to 1)
            order: H(curl) FE order (1=Nedelec first kind, lowest order)
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
        self.delta = float('inf')
        self._assembled = False
        self._solved = False
        self.kappa = 0.01

        # Dense numpy matrices (extracted from NGSolve operators)
        self._a_curl_np = None    # (n_hcurl x n_hcurl) curl-curl stiffness
        self._a_mass_np = None    # (n_hcurl x n_hcurl) mass matrix
        self._B_np = None         # (n_hdiv_free x n_hcurl) trace coupling
        self._S_bem_np = None     # (n_hdiv_free+n_l2 x ...) combined BEM
        self._A_k_np = None       # (n_hdiv_free x n_hdiv_free) Maxwell SLP
        self._V_k_np = None       # (n_l2 x n_l2) scalar SLP
        self._Q_k_np = None       # (n_l2 x n_hdiv_free) divergence coupling

        # DOF counts (full space)
        self._n_hcurl = 0
        self._n_hdiv = 0
        self._n_l2 = 0

        # Active (free) DOF indices and counts
        self._active_hdiv = None  # Free HDivSurface DOF indices
        self._n_hdiv_free = 0     # Number of free HDivSurface DOFs

        # Solution state
        self._A_coeffs = None     # H(curl) solution
        self._j_coeffs = None     # HDivSurface current (free DOFs only)
        self._rho_coeffs = None   # SurfaceL2 charge

        # Timing
        self.t_assemble = 0.0
        self.t_solve = 0.0

    def assemble(self, kappa=0.01, intorder=4):
        """Assemble frequency-independent operators.

        This is the expensive step (done once). After assembly, frequency
        sweeps only require updating the mass matrix coefficient and solving.

        Args:
            kappa: Wavenumber for exterior BEM (small positive value
                   for quasi-static Laplace limit). The Weggler stabilization
                   keeps the condition number O(1) for any kappa.
            intorder: Integration order for BEM singular quadrature.
        """
        from ngsolve import (HCurl, HDivSurface, SurfaceL2,
                              BilinearForm, GridFunction,
                              ds, dx, curl, div, BND)
        from ngsolve.bem import HelmholtzSL

        t0 = time.time()
        self.kappa = kappa

        # --- Function spaces ---
        self._fes_hcurl = HCurl(self.mesh, order=self.order, complex=True)
        self._fes_hdiv = HDivSurface(self.mesh, order=0, complex=True)
        self._fes_l2 = SurfaceL2(self.mesh, order=0, complex=True,
                                  dual_mapping=True)

        self._n_hcurl = self._fes_hcurl.ndof
        self._n_hdiv = self._fes_hdiv.ndof
        self._n_l2 = self._fes_l2.ndof

        # --- Determine free (active) DOFs ---
        # HDivSurface may have non-free DOFs (e.g., wire basket constraints)
        fd_hdiv = self._fes_hdiv.FreeDofs()
        self._active_hdiv = [i for i in range(self._n_hdiv) if fd_hdiv[i]]
        self._n_hdiv_free = len(self._active_hdiv)

        u_hc, v_hc = self._fes_hcurl.TnT()

        # --- FEM interior operators (frequency-independent parts) ---
        # Curl-curl stiffness: (1/mu * curl u, curl v)
        a_curl_bf = BilinearForm(self._fes_hcurl)
        a_curl_bf += (1.0 / self.mu) * curl(u_hc) * curl(v_hc) * dx
        a_curl_bf.Assemble()

        # Mass matrix: (u, v) -- will be multiplied by j*omega*sigma at solve
        a_mass_bf = BilinearForm(self._fes_hcurl)
        a_mass_bf += u_hc * v_hc * dx
        a_mass_bf.Assemble()

        # --- Extract FEM matrices to dense numpy ---
        self._a_curl_np = self._extract_matrix(a_curl_bf.mat, self._n_hcurl)
        self._a_mass_np = self._extract_matrix(a_mass_bf.mat, self._n_hcurl)

        # --- Coupling operator B: tangential trace H(curl) -> HDivSurface ---
        # B[i,j] = <gamma_t(phi_j), psi_i>_Gamma
        # where phi_j is H(curl) basis, psi_i are FREE HDivSurface DOFs
        self._B_np = self._build_trace_coupling()

        # --- BEM operators (Weggler stabilized, product space) ---
        # Must use product space fes_hdiv * fes_l2 for correct extraction.
        # Individual spaces give inconsistent vector dimensions.
        print("  [VectorFEMBEM] Assembling BEM operators (Weggler)...")
        fes_bem = self._fes_hdiv * self._fes_l2
        (uHDiv, uL2), (vHDiv, vL2) = fes_bem.TnT()
        n_bem_full = fes_bem.ndof

        # A_kappa: HelmholtzSL on HDivSurface x HDivSurface (vectorial)
        A_k = HelmholtzSL(
            uHDiv.Trace() * ds(bonus_intorder=intorder), kappa
        ) * vHDiv.Trace() * ds(bonus_intorder=intorder)

        # V_kappa: HelmholtzSL on SurfaceL2 x SurfaceL2 (scalar)
        V_k = HelmholtzSL(
            uL2 * ds(bonus_intorder=intorder), kappa
        ) * vL2 * ds(bonus_intorder=intorder)

        # Q_kappa: HelmholtzSL on div(HDivSurface) x SurfaceL2
        Q_k = HelmholtzSL(
            div(uHDiv.Trace()) * ds(bonus_intorder=intorder), kappa
        ) * vL2 * ds(bonus_intorder=intorder)

        # Combined BEM matrix (Weggler stabilized):
        # S_bem = [A_k, Q_k^T; Q_k, kappa^2*V_k]
        S_bem_mat = (A_k.mat + Q_k.mat + Q_k.mat.T
                     + kappa**2 * V_k.mat)

        # --- Extract combined BEM matrix using FREE DOFs only ---
        # HDivSurface may have non-free DOFs (wire basket constraints).
        # Product space DOFs: [HDivSurface(0..n_hdiv-1), SurfaceL2(n_hdiv..)]
        # Active DOFs in product space = active_hdiv + all SurfaceL2
        active_bem = (self._active_hdiv
                      + [self._n_hdiv + i for i in range(self._n_l2)])
        n_bem = len(active_bem)  # n_hdiv_free + n_l2

        self._S_bem_np = np.zeros((n_bem, n_bem), dtype=complex)
        x_b = S_bem_mat.CreateColVector()
        y_b = S_bem_mat.CreateColVector()
        for idx_j, j in enumerate(active_bem):
            x_b[:] = 0
            x_b[j] = 1.0
            S_bem_mat.Mult(x_b, y_b)
            for idx_i, i in enumerate(active_bem):
                self._S_bem_np[idx_i, idx_j] = y_b[i]

        # Extract individual blocks for diagnostics (free DOFs only)
        n2f = self._n_hdiv_free
        n3 = self._n_l2
        self._A_k_np = self._S_bem_np[:n2f, :n2f].copy()
        self._Q_k_np = self._S_bem_np[n2f:, :n2f].copy()
        if abs(kappa) > 1e-15:
            self._V_k_np = self._S_bem_np[n2f:, n2f:].copy() / (kappa**2)
        else:
            self._V_k_np = self._S_bem_np[n2f:, n2f:].copy()

        self._assembled = True
        self.t_assemble = time.time() - t0
        print(f"  [VectorFEMBEM] Assembly done in {self.t_assemble:.2f} s")
        print(f"    H(curl) DOFs: {self._n_hcurl}")
        print(f"    HDivSurface DOFs: {self._n_hdiv_free} free "
              f"(of {self._n_hdiv})")
        print(f"    SurfaceL2 DOFs: {self._n_l2}")
        print(f"    Total system size: "
              f"{self._n_hcurl + self._n_hdiv_free + self._n_l2}")

    def _extract_matrix(self, mat, ndof):
        """Extract dense numpy matrix from NGSolve BaseMatrix.

        Args:
            mat: NGSolve BaseMatrix
            ndof: Number of DOFs

        Returns:
            M: (ndof, ndof) complex numpy array
        """
        M = np.zeros((ndof, ndof), dtype=complex)
        x = mat.CreateColVector()
        y = mat.CreateColVector()
        for i in range(ndof):
            x[:] = 0
            x[i] = 1.0
            mat.Mult(x, y)
            for j in range(ndof):
                M[j, i] = y[j]
        return M

    def _extract_rect_matrix(self, mat, nrows, ncols):
        """Extract rectangular dense matrix from NGSolve BaseMatrix.

        For operators mapping from trial space (ncols) to test space (nrows).

        Args:
            mat: NGSolve BaseMatrix
            nrows: Number of rows (test DOFs)
            ncols: Number of columns (trial DOFs)

        Returns:
            M: (nrows, ncols) complex numpy array
        """
        M = np.zeros((nrows, ncols), dtype=complex)
        x = mat.CreateColVector()
        y = mat.CreateRowVector()
        for i in range(ncols):
            x[:] = 0
            x[i] = 1.0
            mat.Mult(x, y)
            for j in range(nrows):
                M[j, i] = y[j]
        return M

    def _build_trace_coupling(self):
        """Build tangential trace coupling matrix B (free DOFs only).

        B[i,j] = integral_Gamma gamma_t(phi_j) . psi_i dS

        where phi_j are H(curl) basis functions and psi_i are
        FREE HDivSurface basis functions.

        Returns:
            B: (n_hdiv_free, n_hcurl) complex numpy array
        """
        from ngsolve import BilinearForm, ds

        u_hc, _ = self._fes_hcurl.TnT()
        _, v_hd = self._fes_hdiv.TnT()

        # Try mixed BilinearForm approach
        try:
            B_bf = BilinearForm(trialspace=self._fes_hcurl,
                                testspace=self._fes_hdiv)
            B_bf += u_hc.Trace() * v_hd.Trace() * ds(self.surface_label,
                                                       bonus_intorder=2)
            B_bf.Assemble()

            # Extract to dense numpy (free HDivSurface DOFs only)
            B = np.zeros((self._n_hdiv_free, self._n_hcurl), dtype=complex)
            x = B_bf.mat.CreateColVector()
            y = B_bf.mat.CreateRowVector()
            for i in range(self._n_hcurl):
                x[:] = 0
                x[i] = 1.0
                B_bf.mat.Mult(x, y)
                for idx_j, j in enumerate(self._active_hdiv):
                    B[idx_j, i] = y[j]
            return B

        except Exception as e:
            print(f"  [VectorFEMBEM] Mixed BilinearForm failed: {e}")
            print(f"  [VectorFEMBEM] Falling back to column-by-column "
                  f"extraction")
            return self._build_trace_coupling_manual()

    def _build_trace_coupling_manual(self):
        """Build trace coupling B matrix by manual column-by-column extraction.

        For each H(curl) basis function phi_j:
        1. Set GridFunction to unit vector e_j
        2. Evaluate tangential trace on boundary
        3. Project onto HDivSurface basis (free DOFs only)

        Returns:
            B: (n_hdiv_free, n_hcurl) complex numpy array
        """
        from ngsolve import (GridFunction, BilinearForm, LinearForm,
                              InnerProduct, ds, BND)

        B = np.zeros((self._n_hdiv_free, self._n_hcurl), dtype=complex)

        # H(curl) GridFunction for each basis vector
        _, v_hd = self._fes_hdiv.TnT()
        gf_hc = GridFunction(self._fes_hcurl)

        for j in range(self._n_hcurl):
            gf_hc.vec[:] = 0
            gf_hc.vec[j] = 1.0

            # Build RHS: <gamma_t(phi_j), psi_i> for all i
            rhs_lf = LinearForm(self._fes_hdiv)
            rhs_lf += InnerProduct(gf_hc, v_hd.Trace()) * ds(
                self.surface_label, bonus_intorder=2)
            rhs_lf.Assemble()

            # Extract only free HDivSurface DOFs
            for idx_i, i in enumerate(self._active_hdiv):
                B[idx_i, j] = rhs_lf.vec[i]

        return B

    def solve(self, freq, A_inc_func=None, B_ext=None):
        """Solve eddy current problem at given frequency.

        Uses TOTAL FIELD formulation to avoid catastrophic cancellation:
          Row 1: a_FEM * A_total + B^T * j_scat = 0
          Row 2: B * A_total - S * [j_scat; rho_scat] = B * A_inc

        The unknown A is the TOTAL interior field (not scattered).
        The BEM unknowns j, rho are scattered surface quantities.

        Args:
            freq: Frequency [Hz]
            A_inc_func: Callable(points) -> A_inc vectors (n, 3) array.
                        For PEEC coupling: Biot-Savart vector potential
                        from coil filaments.
            B_ext: External B field [Bx, By, Bz] in Tesla (uniform).
                   Used if A_inc_func is None. Computes A_inc = 0.5 * B x r.

        Returns:
            dict with solution fields and diagnostics
        """
        if not self._assembled:
            raise RuntimeError("Call assemble() first")

        t0 = time.time()

        self.freq = freq
        self.omega = 2.0 * np.pi * freq
        self._B_ext_stored = B_ext  # Store for compute_loss_sibc()

        # Skin depth
        if self.omega > 0 and self.sigma > 0:
            self.delta = np.sqrt(2.0 / (self.omega * self.mu * self.sigma))
        else:
            self.delta = float('inf')

        n1 = self._n_hcurl
        n2 = self._n_hdiv_free   # Free HDivSurface DOFs only
        n3 = self._n_l2
        N = n1 + n2 + n3

        # --- Build frequency-dependent FEM block ---
        a_FEM = self._a_curl_np + 1j * self.omega * self.sigma * self._a_mass_np

        # --- Assemble 3x3 block system ---
        A_sys = np.zeros((N, N), dtype=complex)

        # Block (1,1): a_FEM
        A_sys[:n1, :n1] = a_FEM

        # Block (1,2): +B^T  (Johnson-Nedelec coupling)
        A_sys[:n1, n1:n1+n2] = self._B_np.T

        # Block (2,1): +B
        A_sys[n1:n1+n2, :n1] = self._B_np

        # BEM 2x2 block (Weggler stabilized, from product space)
        # S_bem = [A_k, Q_k^T; Q_k, kappa^2*V_k]
        A_sys[n1:, n1:] = -self._S_bem_np

        # --- Build RHS (TOTAL FIELD formulation) ---
        rhs = np.zeros(N, dtype=complex)

        # Compute A_inc coefficients
        A_inc_coeffs = self._project_incident_field(A_inc_func, B_ext)

        if A_inc_coeffs is not None:
            # Total field formulation:
            #   Row 1: a_FEM * A_total + B^T * j = 0  (homogeneous interior)
            #   Row 2: B * A_total - S * [j; rho] = B * A_inc  (BEM drive)
            #
            # Row 1 RHS = 0 (no source inside conductor)
            # Row 2 RHS = +B @ A_inc (incident field drives via BEM coupling)
            rhs[n1:n1+n2] = self._B_np @ A_inc_coeffs

            # Row 3 RHS = 0 (no charge source)

        # --- Symmetric diagonal scaling (FEM-BEM equilibration) ---
        # FEM block O(10^10) vs BEM block O(10^-2) causes cond ~10^14.
        # Scaling: D^{-1} A D^{-1} * (D x) = D^{-1} rhs
        diag_abs = np.abs(np.diag(A_sys))
        diag_abs[diag_abs < 1e-30] = 1.0  # avoid div by zero
        D_scale = np.sqrt(diag_abs)
        D_inv = 1.0 / D_scale

        A_scaled = A_sys * D_inv[np.newaxis, :] * D_inv[:, np.newaxis]
        rhs_scaled = rhs * D_inv

        # --- Solve scaled system ---
        sol_scaled = np.linalg.solve(A_scaled, rhs_scaled)
        sol = sol_scaled * D_inv  # unscale

        # --- Extract solution ---
        # A_coeffs is now A_TOTAL (not scattered)
        self._A_coeffs = sol[:n1]
        self._j_coeffs = sol[n1:n1+n2]
        self._rho_coeffs = sol[n1+n2:]
        self._A_inc_coeffs = A_inc_coeffs

        self._solved = True
        self.t_solve = time.time() - t0

        return {
            't_solve': self.t_solve,
            'A_max': np.max(np.abs(self._A_coeffs)),
            'j_max': np.max(np.abs(self._j_coeffs)),
            'rho_max': np.max(np.abs(self._rho_coeffs)),
            'n_total': N,
            'cond': np.linalg.cond(A_scaled),
        }

    def _project_incident_field(self, A_inc_func, B_ext):
        """Project incident vector potential onto H(curl) FE space.

        Args:
            A_inc_func: Callable(points) -> A_inc (n,3) array, or None.
            B_ext: Uniform B field [Bx, By, Bz] in Tesla, or None.

        Returns:
            A_inc_coeffs: H(curl) coefficient vector, or None if no source.
        """
        from ngsolve import GridFunction, CF, x, y, z

        if A_inc_func is None and B_ext is None:
            return None

        if B_ext is not None:
            B_ext = np.asarray(B_ext, dtype=float)

            # A_inc = 0.5 * B x r (Coulomb gauge, uniform field)
            # A = 0.5 * (B x r) = 0.5 * (By*z - Bz*y, Bz*x - Bx*z, Bx*y - By*x)
            Bx, By, Bz = B_ext
            A_inc_cf = CF((0.5 * (By * z - Bz * y),
                           0.5 * (Bz * x - Bx * z),
                           0.5 * (Bx * y - By * x)))

            gf = GridFunction(self._fes_hcurl)
            gf.Set(A_inc_cf)
            return np.array([gf.vec[i] for i in range(self._n_hcurl)])

        # A_inc_func provided: evaluate and project
        if A_inc_func is not None:
            return self._project_func_to_hcurl(A_inc_func)

        return None

    def _project_func_to_hcurl(self, A_func):
        """Project a vector function onto H(curl) using L2 projection.

        Args:
            A_func: Callable(points) -> (n, 3) array of vector potential values.

        Returns:
            coeffs: H(curl) coefficient vector
        """
        from ngsolve import GridFunction

        # Evaluate A_func at mesh integration points and project
        # Simple approach: evaluate at edge midpoints for order=1 Nedelec
        gf = GridFunction(self._fes_hcurl)

        # Use NGSolve's Set with a CoefficientFunction wrapper
        # Build a piecewise-constant CF from sampled values
        try:
            from ngsolve import CF, x, y, z

            # Create CF that evaluates A_func
            # For numerical integration, we sample at enough points
            class VectorPotentialCF:
                def __init__(self, func):
                    self.func = func

                def __call__(self, mip):
                    pt = np.array([mip.point[0], mip.point[1], mip.point[2]])
                    A = self.func(pt.reshape(1, 3))
                    return tuple(A[0])

            # Direct approach: evaluate on mesh points
            # Sample at vertex locations and create interpolation
            n_vert = self.mesh.nv
            pts = np.zeros((n_vert, 3))
            for i in range(n_vert):
                p = self.mesh.vertices[i].point
                pts[i] = [p[0], p[1], p[2]]

            A_vals = np.atleast_2d(A_func(pts))

            # Create CF from sampled values
            # Use the average value as a constant approximation per element
            A_avg = np.mean(A_vals, axis=0)
            A_cf = CF((A_avg[0], A_avg[1], A_avg[2]))
            gf.Set(A_cf)

        except Exception:
            # Fallback: zero field
            gf.vec[:] = 0

        return np.array([gf.vec[i] for i in range(self._n_hcurl)])

    def compute_loss(self):
        """Compute eddy current loss from volume integral.

        P = 0.5 * Re(integral sigma * |E_total|^2 dV)
          = 0.5 * omega^2 * sigma * Re(A_total^H * M_mass * A_total)

        where E_total = -j*omega*A_total.

        WARNING: This formula requires the mesh to resolve the skin depth.
        When maxh >> delta (skin depth), the FEM produces A_total ~ 0
        (perfect conductor limit) and P ~ 0. Use compute_loss_sibc()
        for a mesh-independent SIBC-based estimate, or use
        ShieldBEMSIBC for thin-skin problems.

        Returns:
            P: Eddy current loss [W]
        """
        if not self._solved:
            raise RuntimeError("Call solve() first")

        # _A_coeffs is already A_total (total field formulation)
        A_total = self._A_coeffs

        # P = 0.5 * sigma * omega^2 * integral |A_total|^2 dV
        P = (0.5 * self.omega**2 * self.sigma
             * np.real(A_total.conj() @ self._a_mass_np @ A_total))
        return P

    def compute_loss_sibc(self):
        """Compute eddy current loss using analytical SIBC on boundary faces.

        For each boundary triangle, computes:
            H_inc_t = H_inc - (H_inc . n) * n   (tangential H)
            dP = 0.5 * Re(Zs) * |H_inc_t|^2 * area

        This is mesh-INDEPENDENT and valid when delta << thickness
        (thin-skin regime where compute_loss() fails on coarse mesh).

        NOTE: This gives the half-space SIBC estimate. It does not account
        for finite-body edge/corner enhancement effects. For accurate
        thin-skin loss, use ShieldBEMSIBC which includes BEM mutual
        coupling effects.

        Returns:
            P: Eddy current loss [W] (analytical SIBC estimate)
        """
        if not self._solved:
            raise RuntimeError("Call solve() first")
        if self.omega <= 0 or self.sigma <= 0:
            return 0.0

        # Surface impedance
        delta = np.sqrt(2.0 / (self.omega * self.mu * self.sigma))
        Zs = (1 + 1j) / (self.sigma * delta)

        # Need B_ext to compute H_inc; extract from stored A_inc
        if not hasattr(self, '_B_ext_stored') or self._B_ext_stored is None:
            return 0.0

        B_ext = np.asarray(self._B_ext_stored, dtype=float)
        H_inc = B_ext / self.mu  # H_inc = B_ext / (mu_r * mu_0)

        from ngsolve import BND

        P_total = 0.0
        for el in self.mesh.Elements(BND):
            # Get triangle vertices
            verts = [self.mesh.vertices[v.nr].point for v in el.vertices]
            v0 = np.array(verts[0])
            v1 = np.array(verts[1])
            v2 = np.array(verts[2])

            # Outward normal and area
            edge1 = v1 - v0
            edge2 = v2 - v0
            cross = np.cross(edge1, edge2)
            area = 0.5 * np.linalg.norm(cross)
            if area < 1e-30:
                continue
            n_hat = cross / (2.0 * area)

            # Tangential H_inc
            H_t = H_inc - np.dot(H_inc, n_hat) * n_hat
            H_t_mag2 = np.dot(H_t, H_t)

            # SIBC loss contribution: 0.5 * Re(Zs) * |H_t|^2 * area
            P_total += 0.5 * Zs.real * H_t_mag2 * area

        return P_total

    def compute_stored_energy(self):
        """Compute magnetic stored energy.

        W = 0.25 * Re(integral (1/mu) * |curl A|^2 dV)
          = 0.25 * Re(A^H * a_curl * A)

        Returns:
            W: Magnetic stored energy [J]
        """
        if not self._solved:
            raise RuntimeError("Call solve() first")

        W = 0.25 * np.real(
            self._A_coeffs.conj() @ self._a_curl_np @ self._A_coeffs)
        return W

    def frequency_sweep(self, freqs, A_inc_func=None, B_ext=None):
        """Run frequency sweep (fast: only FEM mass coefficient changes).

        Args:
            freqs: Array of frequencies [Hz]
            A_inc_func: Incident field function (see solve())
            B_ext: Uniform B field (see solve())

        Returns:
            dict with arrays of P_loss, W_stored, delta per frequency
        """
        if not self._assembled:
            raise RuntimeError("Call assemble() first")

        freqs = np.asarray(freqs, dtype=float)
        n_freq = len(freqs)

        P_loss = np.zeros(n_freq)
        P_loss_sibc = np.zeros(n_freq)
        W_stored = np.zeros(n_freq)
        delta_arr = np.zeros(n_freq)

        t0 = time.time()
        for i, f in enumerate(freqs):
            self.solve(f, A_inc_func, B_ext)
            P_loss[i] = self.compute_loss()
            P_loss_sibc[i] = self.compute_loss_sibc()
            W_stored[i] = self.compute_stored_energy()
            delta_arr[i] = self.delta
        t_total = time.time() - t0

        return {
            'freqs': freqs,
            'P_loss': P_loss,
            'P_loss_sibc': P_loss_sibc,
            'W_stored': W_stored,
            'delta': delta_arr,
            't_total': t_total,
        }
