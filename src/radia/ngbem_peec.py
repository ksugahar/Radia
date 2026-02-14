"""
ngbem_peec.py

PEEC Loop-Star Solver using ngbem (NGSolve BEM)

Computes port impedance Z(f) for conductors using Galerkin BEM with
natural Loop-Star decomposition via HDivSurface / SurfaceL2 function spaces.

Block system (Loop-Star):
    | R + jwL       M_LS^T  |   | I_loop |   | V_port |
    | M_LS      P/(jw)      | * | Q_star | = | 0      |

where:
    L = mu_0 * LaplaceSL(HDivSurface)     - Inductance (edge DOFs)
    P = SingleLayerPotentialOperator(SurfaceL2) / eps_0  - Potential coeff (cell DOFs)
    M_LS = int div(J_edge) * phi_cell dS   - Divergence coupling
    R = surface resistance from sigma, thickness

Modes:
    'mqs':  Loop-only (Z = R + jwL), valid for DC - 1 MHz
    'full': Full Loop-Star with capacitive effects

Usage:
    from ngsolve import Mesh
    from ngbem_peec import NGBEMPEECSolver, create_plate_mesh

    mesh = create_plate_mesh(0.01, 0.01, 0.003, label="conductor")
    solver = NGBEMPEECSolver(mesh, conductor_label="conductor",
                              sigma=5.8e7, thickness=0.035e-3)
    solver.assemble()

    freqs = np.logspace(3, 9, 50)
    Z = solver.solve_frequency(freqs, mode='full')

Part of Radia project
"""

import numpy as np
import time

# Physical constants
MU_0 = 4.0 * np.pi * 1e-7   # H/m
EPS_0 = 8.854187817e-12      # F/m


def extract_dense_matrix(mat, ndof):
    """Extract dense numpy matrix from NGSolve BaseMatrix.

    Args:
        mat: NGSolve BaseMatrix (from IntegralOperator.mat)
        ndof: Number of degrees of freedom

    Returns:
        Dense numpy array (ndof x ndof)
    """
    M = np.zeros((ndof, ndof))
    for i in range(ndof):
        ei = mat.CreateColVector()
        ei[:] = 0
        ei[i] = 1.0
        col = mat.CreateColVector()
        mat.Mult(ei, col)
        for j in range(ndof):
            M[j, i] = col[j]
    return M


def create_plate_mesh(width, height, maxh, label="conductor"):
    """Create rectangular plate surface mesh using Netgen OCC.

    Args:
        width: Plate width [m]
        height: Plate height [m]
        maxh: Maximum element size [m]
        label: Boundary label for the conductor surface

    Returns:
        mesh: NGSolve Mesh (surface in 3D)
    """
    from netgen.occ import WorkPlane, OCCGeometry, Axes, Pnt, Dir
    from ngsolve import Mesh

    wp = WorkPlane(Axes(p=Pnt(0, 0, 0), n=Dir(0, 0, 1), h=Dir(1, 0, 0)))
    rect = wp.Rectangle(width, height).Face()
    rect.name = label

    geo = OCCGeometry(rect)
    ngmesh = geo.GenerateMesh(maxh=maxh)
    return Mesh(ngmesh)


def get_mesh_triangles(mesh):
    """Extract triangle vertices and areas from NGSolve surface mesh.

    Args:
        mesh: NGSolve Mesh

    Returns:
        triangles: List of [[x1,y1,z1], [x2,y2,z2], [x3,y3,z3]]
        areas: numpy array of triangle areas [m^2]
    """
    from ngsolve import BND

    triangles = []
    areas = []
    for el in mesh.Elements(BND):
        verts = []
        for v in el.vertices:
            pt = mesh.vertices[v.nr].point
            verts.append([pt[0], pt[1], pt[2]])
        if len(verts) == 3:
            triangles.append(verts)
            e1 = np.array(verts[1]) - np.array(verts[0])
            e2 = np.array(verts[2]) - np.array(verts[0])
            areas.append(0.5 * np.linalg.norm(np.cross(e1, e2)))

    return triangles, np.array(areas)


class NGBEMPEECSolver:
    """PEEC solver using ngbem Galerkin BEM with Loop-Star decomposition.

    Uses NGSolve BEM (ngbem) for matrix assembly:
    - HDivSurface (edge-based RWG) for inductance L (Loop DOFs)
    - SurfaceL2 (cell-based piecewise constant) for potential P (Star DOFs)
    - Divergence bilinear form for Loop-Star coupling M_LS

    This provides:
    - Galerkin discretization (symmetric matrices by construction)
    - High-order elements (order > 0)
    - Natural Loop-Star decomposition (no explicit basis transformation)
    - FMM acceleration available for large problems
    """

    def __init__(self, mesh, conductor_label="conductor", sigma=5.8e7,
                 thickness=None, order=0, intorder=5):
        """Initialize ngbem PEEC solver.

        Args:
            mesh: NGSolve Mesh (surface mesh of conductor)
            conductor_label: Boundary label for conductor surface
            sigma: Conductivity [S/m] (default: copper 5.8e7)
            thickness: Conductor thickness [m] for resistance computation.
                       If None, resistance is not computed.
            order: Finite element order (0=piecewise constant, 1=linear, ...)
            intorder: Integration order for BEM singular quadrature (5-7 recommended)
        """
        self.mesh = mesh
        self.conductor_label = conductor_label
        self.sigma = sigma
        self.thickness = thickness
        self.order = order
        self.intorder = intorder

        # Matrices (populated by assemble())
        self.L = None       # Inductance [H], (n_loop x n_loop)
        self.P = None       # Potential coefficients [1/F], (n_star x n_star)
        self.M_LS = None    # Loop-Star coupling, (n_star x n_loop)
        self.R_loop = None  # Loop resistance [Ohm], (n_loop,) diagonal

        # DOF counts
        self.n_loop = 0     # Edge DOFs (HDivSurface)
        self.n_star = 0     # Cell DOFs (SurfaceL2)

        # Timing
        self.t_assemble = 0.0

    def assemble(self):
        """Assemble L, P, M_LS, R matrices using ngbem.

        This is the main computation step. Call once, then use
        solve_frequency() for multiple frequency sweeps.
        """
        from ngsolve import HDivSurface, SurfaceL2, TaskManager, ds
        from ngsolve import BilinearForm, BND
        from ngsolve.bem import LaplaceSL, SingleLayerPotentialOperator

        t_start = time.perf_counter()
        label = self.conductor_label

        # --- Function spaces ---
        self._fes_hdiv = HDivSurface(self.mesh, order=self.order)
        self._fes_l2 = SurfaceL2(self.mesh, order=self.order,
                                  dual_mapping=False)

        self.n_loop = self._fes_hdiv.ndof
        self.n_star = self._fes_l2.ndof

        # --- L matrix: mu_0 * vector single layer on HDivSurface ---
        j_trial = self._fes_hdiv.TrialFunction()
        j_test = self._fes_hdiv.TestFunction()

        with TaskManager():
            L_op = LaplaceSL(j_trial * ds(label)) * j_test * ds(label)

        L_dense = extract_dense_matrix(L_op.mat, self.n_loop)
        self.L = MU_0 * L_dense

        # --- P matrix: scalar single layer on SurfaceL2 / eps_0 ---
        with TaskManager():
            V_op = SingleLayerPotentialOperator(self._fes_l2,
                                                intorder=self.intorder)

        V_dense = extract_dense_matrix(V_op.mat, self.n_star)
        self.P = V_dense / EPS_0

        # --- M_LS: divergence coupling (FEM bilinear form, NOT BEM) ---
        self.M_LS = self._build_divergence_coupling()

        # --- R: resistance from conductivity and geometry ---
        if self.thickness is not None and self.sigma > 0:
            self.R_loop = self._compute_edge_resistance()
        else:
            self.R_loop = np.zeros(self.n_loop)

        self.t_assemble = time.perf_counter() - t_start

    def _build_divergence_coupling(self):
        """Build Loop-Star coupling matrix M_LS.

        M_LS[i,j] = int_cell_i div(J_edge_j) dS

        where J_edge_j is the j-th RWG basis function (HDivSurface)
        and cell_i is the i-th triangle (SurfaceL2 DOF).

        This is a pure FEM operation (no BEM kernel).

        Returns:
            M_LS: numpy array (n_star x n_loop)
        """
        from ngsolve import BilinearForm, ds

        # Create product space
        fes_product = self._fes_hdiv * self._fes_l2
        (u_hdiv, u_l2), (v_hdiv, v_l2) = fes_product.TnT()

        # The divergence coupling: div(j) tested against phi
        # In the product space, this couples the first component (HDivSurface)
        # with the second component (SurfaceL2)
        #
        # M_LS[star_i, loop_j] = int div(J_j) * phi_i dS
        #
        # Using the product space bilinear form:
        from ngsolve import div
        bf = BilinearForm(fes_product)
        bf += div(u_hdiv.Trace()) * v_l2 * ds
        bf.Assemble()

        # Extract the off-diagonal block (HDivSurface -> SurfaceL2)
        # The product space has DOFs ordered as [HDivSurface, SurfaceL2]
        # We need the (1,0) block: rows=SurfaceL2, cols=HDivSurface
        n_loop = self.n_loop
        n_star = self.n_star
        n_total = n_loop + n_star

        M_full = extract_dense_matrix(bf.mat, n_total)

        # Extract the (1,0) block: M_LS = M_full[n_loop:, :n_loop]
        M_LS = M_full[n_loop:n_loop + n_star, :n_loop]

        return M_LS

    def _compute_edge_resistance(self):
        """Compute per-edge DC resistance from conductivity and thickness.

        For a surface conductor with thickness t and conductivity sigma,
        the sheet resistance is R_sheet = 1 / (sigma * t) [Ohm/square].

        For an RWG edge basis function with edge length l_e and
        associated area A_e, the resistance contribution is:
            R_edge ~ R_sheet * l_e / w_e
        where w_e ~ A_e / l_e is the effective width.

        Returns:
            R: numpy array (n_loop,) of per-edge resistance [Ohm]
        """
        from ngsolve import BND

        R_sheet = 1.0 / (self.sigma * self.thickness)

        # Get mesh edge lengths
        # For order=0 HDivSurface, each DOF corresponds to one internal edge
        # We approximate edge geometry from the mesh
        triangles, areas = get_mesh_triangles(self.mesh)
        n_tri = len(triangles)

        # For simple resistance estimate: use average element properties
        # More accurate: compute per-edge length from mesh topology
        total_area = np.sum(areas)
        mean_edge_length = np.sqrt(4.0 * total_area / (np.sqrt(3) * n_tri))

        # Each edge has two adjacent triangles with combined area ~ 2*A_mean
        # Effective width ~ 2*A_mean / l_edge
        A_mean = total_area / n_tri
        w_eff = 2.0 * A_mean / mean_edge_length

        R_per_edge = R_sheet * mean_edge_length / w_eff
        return np.full(self.n_loop, R_per_edge)

    def solve_frequency(self, freqs, mode='mqs'):
        """Compute port impedance Z(f) over a frequency range.

        Args:
            freqs: Array of frequencies [Hz]
            mode: 'mqs' for Loop-only (fast), 'full' for Loop-Star

        Returns:
            Z_port: Complex impedance array, shape (len(freqs),)
        """
        if self.L is None:
            raise RuntimeError("Call assemble() before solve_frequency()")

        freqs = np.asarray(freqs, dtype=float)
        Z_port = np.zeros(len(freqs), dtype=complex)

        for i, f in enumerate(freqs):
            omega = 2.0 * np.pi * f

            if mode == 'mqs':
                Z_port[i] = self._solve_mqs(omega)
            else:
                Z_port[i] = self._solve_full_loop_star(omega)

        return Z_port

    def _solve_mqs(self, omega):
        """MQS mode: Loop-only impedance (capacitive effects neglected).

        Z_branch = diag(R) + jw*L
        Z_port = 1 / (e^T * Z_branch^{-1} * e)

        where e is the port excitation vector.
        """
        Z_branch = np.diag(self.R_loop.astype(complex)) + 1j * omega * self.L

        # For simple total impedance (series assumption):
        # Z_port = sum of all Z_branch elements
        # This is valid for a 1-port single-conductor problem.
        #
        # For general multi-port, need proper MNA formulation.
        # Use sum-of-all approximation for now.
        try:
            Z_inv = np.linalg.inv(Z_branch)
            # Port vector: uniform excitation (all edges carry same current)
            e = np.ones(self.n_loop) / self.n_loop
            Z_port = 1.0 / (e @ Z_inv @ e)
        except np.linalg.LinAlgError:
            # Fallback: sum of diagonal
            Z_port = np.sum(self.R_loop) + 1j * omega * np.sum(self.L)

        return Z_port

    def _solve_full_loop_star(self, omega):
        """Full Loop-Star: includes capacitive (Star) DOFs.

        Block system:
            | Z_LL    M_LS^T  |   | I_L |   | V_port |
            | M_LS    Z_SS    | * | Q_S | = | 0      |

        where:
            Z_LL = diag(R) + jw*L   (n_loop x n_loop)
            Z_SS = P / (jw)          (n_star x n_star)
            M_LS                      (n_star x n_loop)

        Port impedance via Schur complement:
            Z_port = e^T * (Z_LL - M_LS^T * Z_SS^{-1} * M_LS)^{-1} * e
        """
        # Loop-Loop block
        Z_LL = np.diag(self.R_loop.astype(complex)) + 1j * omega * self.L

        if omega < 1e-10:
            # DC limit: only resistance
            return np.sum(self.R_loop) + 0j

        # Star-Star block
        Z_SS = self.P / (1j * omega)

        # Schur complement: Z_eff = Z_LL - M_LS^T * Z_SS^{-1} * M_LS
        try:
            Z_SS_inv = np.linalg.inv(Z_SS)
            Schur = Z_LL - self.M_LS.T @ Z_SS_inv @ self.M_LS
        except np.linalg.LinAlgError:
            # If P is singular, fall back to MQS
            Schur = Z_LL

        # Port extraction (uniform excitation approximation)
        try:
            Schur_inv = np.linalg.inv(Schur)
            e = np.ones(self.n_loop) / self.n_loop
            Z_port = 1.0 / (e @ Schur_inv @ e)
        except np.linalg.LinAlgError:
            Z_port = np.sum(self.R_loop) + 1j * omega * np.sum(self.L)

        return Z_port

    def get_matrices(self):
        """Return assembled matrices for external use.

        Returns:
            dict with keys: 'L', 'P', 'M_LS', 'R_loop',
                            'n_loop', 'n_star'
        """
        return {
            'L': self.L,
            'P': self.P,
            'M_LS': self.M_LS,
            'R_loop': self.R_loop,
            'n_loop': self.n_loop,
            'n_star': self.n_star,
        }

    def print_summary(self):
        """Print summary of assembled matrices."""
        if self.L is None:
            print("NGBEMPEECSolver: Not assembled yet. Call assemble() first.")
            return

        print(f"NGBEMPEECSolver Summary:")
        print(f"  Conductor label: {self.conductor_label}")
        print(f"  FE order: {self.order}")
        print(f"  Loop DOFs (edges, HDivSurface): {self.n_loop}")
        print(f"  Star DOFs (cells, SurfaceL2):   {self.n_star}")
        print(f"  Total DOFs:                     {self.n_loop + self.n_star}")
        print(f"  Assembly time: {self.t_assemble:.3f} s")
        print()

        # L matrix properties
        L_eig = np.linalg.eigvalsh(self.L)
        L_sym = np.linalg.norm(self.L - self.L.T) / np.linalg.norm(self.L)
        print(f"  L matrix ({self.n_loop}x{self.n_loop}):")
        print(f"    Diagonal range: [{np.min(np.diag(self.L)):.4e}, "
              f"{np.max(np.diag(self.L)):.4e}] H")
        print(f"    Symmetry ||L-L^T||/||L||: {L_sym:.2e}")
        print(f"    Eigenvalues: {np.sum(L_eig > 0)} positive, "
              f"{np.sum(L_eig <= 0)} non-positive")
        print(f"    Eigenvalue range: [{L_eig[0]:.4e}, {L_eig[-1]:.4e}]")
        print()

        # P matrix properties
        P_eig = np.linalg.eigvalsh(self.P)
        P_sym = np.linalg.norm(self.P - self.P.T) / np.linalg.norm(self.P)
        print(f"  P matrix ({self.n_star}x{self.n_star}):")
        print(f"    Diagonal range: [{np.min(np.diag(self.P)):.4e}, "
              f"{np.max(np.diag(self.P)):.4e}] 1/F")
        print(f"    Symmetry ||P-P^T||/||P||: {P_sym:.2e}")
        print(f"    Eigenvalues: {np.sum(P_eig > 0)} positive, "
              f"{np.sum(P_eig <= 0)} non-positive")
        print(f"    Eigenvalue range: [{P_eig[0]:.4e}, {P_eig[-1]:.4e}]")
        print()

        # M_LS properties
        print(f"  M_LS matrix ({self.n_star}x{self.n_loop}):")
        print(f"    Range: [{np.min(self.M_LS):.4e}, {np.max(self.M_LS):.4e}]")
        print(f"    Nonzeros: {np.count_nonzero(np.abs(self.M_LS) > 1e-15)}"
              f" / {self.M_LS.size}")
        print()

        # Resistance
        print(f"  R_loop: [{np.min(self.R_loop):.4e}, "
              f"{np.max(self.R_loop):.4e}] Ohm")
        if self.thickness:
            R_sheet = 1.0 / (self.sigma * self.thickness)
            print(f"  Sheet resistance: {R_sheet:.4e} Ohm/sq")
