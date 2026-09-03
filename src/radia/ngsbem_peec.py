"""
ngsbem_peec.py

PEEC Loop-Star Solver using ngsbem (NGSolve BEM)

Computes port impedance Z(f) for conductors using Galerkin BEM with
natural Loop-Star decomposition via HDivSurface / SurfaceL2 function spaces.

Block system (Loop-Star):
    | R + jwL       M_LS^T  |   | I_loop |   | V_port |
    | M_LS      P/(jw)      | * | Q_star | = | 0      |

where:
    L = mu_0 * LaplaceSL(HDivSurface)     - Inductance (edge DOFs)
    P = LaplaceSL(SurfaceL2) / eps_0      - Potential coeff (cell DOFs)
    M_LS = int div(J_edge) * phi_cell dS   - Divergence coupling
    R = surface resistance from sigma, thickness

Modes:
    'mqs':  Loop-only (Z = R + jwL), valid for DC - 1 MHz
    'full': Full Loop-Star (reformulated Schur complement, numerically stable)
    'stabilized': Weggler's stabilized BEM (Q_0 from LaplaceSL, most robust)

Usage:
    from ngsolve import Mesh
    from radia.ngsbem_peec import NGBEMPEECSolver, create_plate_mesh

    mesh = create_plate_mesh(0.01, 0.01, 0.003, label="conductor")
    solver = NGBEMPEECSolver(mesh, conductor_label="conductor",
                              sigma=5.8e7, thickness=0.035e-3)
    solver.assemble()

    freqs = np.logspace(1, 6, 50)
    Z = solver.solve_frequency(freqs, mode='full')

Part of Radia project
"""

import numpy as np
import scipy.linalg
import time

# Physical constants
MU_0 = 4.0 * np.pi * 1e-7   # H/m
EPS_0 = 8.854187817e-12      # F/m
C_0 = 1.0 / np.sqrt(MU_0 * EPS_0)  # speed of light [m/s]


def _bem_measure_pair(ds, label, intorder, trial_order, test_order):
    """Return measures matching an absolute BEM integration order.

    NGSolve's variational BEM API adds two and both FE orders to the two
    bonus-order values, so the remaining order belongs on one measure.
    """
    source_bonus = max(
        0, int(intorder) - 2 - int(trial_order) - int(test_order)
    )
    return ds(label, bonus_intorder=source_bonus), ds(label)


# ============================================================
# Dowell formula for skin effect in rectangular conductors
# ============================================================

def dowell_F_R(xi):
    """Dowell AC resistance factor: R_ac/R_dc.

    F_R(xi) = xi * (sinh(2*xi) + sin(2*xi)) / (cosh(2*xi) - cos(2*xi))

    where xi = d / (2*delta), d = conductor thickness, delta = skin depth.

    Valid for rectangular conductors with d << w (thin plate approximation).

    Args:
        xi: Normalized half-thickness d/(2*delta). Scalar or array.

    Returns:
        F_R: AC/DC resistance ratio (>=1). Same shape as xi.
    """
    if np.isscalar(xi):
        if xi < 0.001:
            return 1.0
        xi2 = 2 * xi
        return xi * (np.sinh(xi2) + np.sin(xi2)) / (np.cosh(xi2) - np.cos(xi2))
    else:
        xi = np.asarray(xi, dtype=float)
        result = np.ones_like(xi)
        mask = xi >= 0.001
        xi_m = xi[mask]
        xi2 = 2 * xi_m
        result[mask] = xi_m * (np.sinh(xi2) + np.sin(xi2)) / (np.cosh(xi2) - np.cos(xi2))
        return result


def dowell_F_L(xi):
    """Dowell AC inductance factor: L_ac/L_dc (internal inductance only).

    F_L(xi) = 3*(sinh(2*xi) - sin(2*xi)) / (xi^3 * (cosh(2*xi) - cos(2*xi)))

    Args:
        xi: Normalized half-thickness d/(2*delta). Scalar or array.

    Returns:
        F_L: AC/DC internal inductance ratio (<=1). Same shape as xi.
    """
    if np.isscalar(xi):
        if xi < 0.001:
            return 1.0
        xi2 = 2 * xi
        denom = xi**3 * (np.cosh(xi2) - np.cos(xi2))
        if abs(denom) < 1e-30:
            return 1.0
        return 3.0 * (np.sinh(xi2) - np.sin(xi2)) / denom
    else:
        xi = np.asarray(xi, dtype=float)
        result = np.ones_like(xi)
        mask = xi >= 0.001
        xi_m = xi[mask]
        xi2 = 2 * xi_m
        denom = xi_m**3 * (np.cosh(xi2) - np.cos(xi2))
        valid = np.abs(denom) > 1e-30
        mask2 = mask.copy()
        mask2[mask] &= valid
        result[mask2] = 3.0 * (np.sinh(2 * xi[mask2]) - np.sin(2 * xi[mask2])) / denom[valid]
        return result


def make_dowell_Zs_func(n_loop, R_dc_per_edge, thickness, sigma, mu_r=1.0):
    """Create Dowell Zs_func callable for use with solve_frequency(Zs_func=...).

    The Dowell formula gives the AC resistance factor F_R(xi) where
    xi = d/(2*delta), d = thickness, delta = skin depth = sqrt(2/(w*mu*sigma)).

    Zs[i] = R_dc[i] * (F_R(xi) - 1)  (excess resistance beyond DC)

    Args:
        n_loop: Number of loop DOFs (edges)
        R_dc_per_edge: DC resistance per edge, array (n_loop,) [Ohm]
        thickness: Conductor thickness [m]
        sigma: Conductivity [S/m]
        mu_r: Relative permeability (default 1.0)

    Returns:
        Callable(freq) -> Zs array (n_loop,) complex
    """
    mu = mu_r * MU_0

    def Zs_func(freq):
        if freq <= 0:
            return np.zeros(n_loop, dtype=complex)
        omega = 2.0 * np.pi * freq
        delta = np.sqrt(2.0 / (omega * mu * sigma))
        xi = thickness / (2.0 * delta)
        F_R = dowell_F_R(xi)
        # Excess AC resistance (beyond R_dc which is already in R_loop)
        return R_dc_per_edge * (F_R - 1.0) + 0j

    return Zs_func


def extract_dense_matrix(mat, ndof):
    """Extract dense numpy matrix from NGSolve BaseMatrix.

    Uses COO() extraction when available (SparseMatrix with use_fmm=False),
    which is ~2500x faster than column-by-column MatVec for BEM matrices.
    Falls back to MatVec for SumMatrix (use_fmm=True) or other types.

    Args:
        mat: NGSolve BaseMatrix (from IntegralOperator.mat)
        ndof: Number of degrees of freedom

    Returns:
        Dense numpy array (ndof x ndof)
    """
    try:
        from scipy.sparse import coo_matrix
        rows, cols, vals = mat.COO()
        return coo_matrix((vals, (rows, cols)),
                          shape=(mat.height, mat.width)).toarray()
    except Exception:
        pass

    # Fallback: column-by-column MatVec (slow but universal)
    ei = mat.CreateColVector()
    ei[:] = 0
    ei[0] = 1.0
    col = mat.CreateColVector()
    mat.Mult(ei, col)
    is_complex = isinstance(col[0], complex)

    M = np.zeros((ndof, ndof), dtype=complex if is_complex else float)
    for j in range(ndof):
        M[j, 0] = col[j]
    for i in range(1, ndof):
        ei[:] = 0
        ei[i] = 1.0
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
    """PEEC solver using ngsbem Galerkin BEM with Loop-Star decomposition.

    Uses NGSolve BEM (ngsbem) for matrix assembly:
    - HDivSurface (edge-based RWG) for inductance L (Loop DOFs)
    - SurfaceL2 (cell-based piecewise constant) for potential P (Star DOFs)
    - Divergence bilinear form for Loop-Star coupling M_LS

    This provides:
    - Galerkin discretization (symmetric matrices by construction)
    - High-order elements (order > 0)
    - Natural Loop-Star decomposition (no explicit basis transformation)


    Solver notes:
    - Dense direct: LDL^T (Bunch-Kaufman) for complex symmetric A = A^T
    - Iterative (large-scale): GMRes is recommended over COCG/COCR because
      BEM-discretized systems exhibit erratic COCG convergence without
      specialized spectral preconditioners (SSOR, incomplete Cholesky).
      GMRes guarantees monotonic residual minimization.
    """

    def __init__(self, mesh, conductor_label="conductor", sigma=5.8e7,
                 thickness=None, order=0, intorder=5):
        """Initialize ngsbem PEEC solver.

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
        self.M_LS = None    # Loop-Star coupling (FEM), (n_star x n_loop)
        self.Q_0 = None     # BEM coupling LaplaceSL(div(J), rho), (n_star x n_loop)
        self.V_0 = None     # Raw LaplaceSL on SurfaceL2, (n_star x n_star)
        self.R_loop = None  # Loop resistance [Ohm], (n_loop,) diagonal
        self._P_inv_M_LS = None  # Precomputed P^{-1} @ M_LS for Schur complement

        # DOF counts
        self.n_loop = 0     # Edge DOFs (HDivSurface)
        self.n_star = 0     # Cell DOFs (SurfaceL2)

        # Timing
        self.t_assemble = 0.0

    def assemble(self):
        """Assemble L, P, Q_0, M_LS, R matrices using ngsbem.

        This is the main computation step. Call once, then use
        solve_frequency() for multiple frequency sweeps.

        Assembles:
        - L: mu_0 * LaplaceSL(HDivSurface) — inductance
        - P: LaplaceSL(SurfaceL2) / eps_0 — potential coefficient
        - V_0: LaplaceSL(SurfaceL2) — raw BEM scalar SL
        - Q_0: LaplaceSL(div(HDivSurface), SurfaceL2) — BEM coupling
        - M_LS: int div(J) * phi dS — local FEM coupling
        - R_loop: DC resistance per edge
        """
        from ngsolve import HDivSurface, SurfaceL2, ds
        from ngsolve import BilinearForm, BND
        from ngsolve.bem import LaplaceSL

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

        # .Trace() is required for correct BEM integration.
        # Without it, boundary-edge DOFs get corrupted diagonal
        # entries (e.g. -1e+11) due to incorrect kernel evaluation.
        L_op = LaplaceSL(
            j_trial.Trace() * ds(label)
        ) * j_test.Trace() * ds(label)

        L_dense = extract_dense_matrix(L_op.mat, self.n_loop)
        self.L = MU_0 * L_dense

        # --- V_0 / P: scalar single layer on SurfaceL2 ---
        uL2, vL2 = self._fes_l2.TnT()
        source_ds, test_ds = _bem_measure_pair(
            ds, label, self.intorder, self.order, self.order
        )
        V_op = LaplaceSL(uL2 * source_ds) * vL2 * test_ds

        V_dense = extract_dense_matrix(V_op.mat, self.n_star)
        self.V_0 = V_dense              # Raw LaplaceSL (for stabilized mode)
        self.P = V_dense / EPS_0         # Potential coefficient P = V_0/eps_0

        # --- Q_0: BEM coupling LaplaceSL(div(HDivSurface), SurfaceL2) ---
        self.Q_0 = self._build_bem_coupling()

        # --- M_LS: divergence coupling (FEM bilinear form, NOT BEM) ---
        self.M_LS = self._build_divergence_coupling()

        # --- Precompute P^{-1} @ M_LS for Schur complement (mode='full') ---
        try:
            self._P_inv_M_LS = scipy.linalg.solve(
                self.P, self.M_LS, assume_a='sym')
        except np.linalg.LinAlgError:
            self._P_inv_M_LS = None

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

    def _build_bem_coupling(self):
        """Build BEM coupling matrix Q_0 using LaplaceSL.

        Q_0[i,j] = LaplaceSL(div(J_edge_j), phi_cell_i)
                 = int_S int_S div(J_j(r')) * G_0(r,r') * phi_i(r) dS dS'

        where G_0(r,r') = 1/(4*pi*|r-r'|) is the Laplace Green's function.

        This is a BEM operator (involves Green's function), unlike M_LS
        which is a local FEM bilinear form. Q_0 is needed for the
        stabilized formulation (Weggler).

        Uses the product space (HDivSurface * SurfaceL2) approach
        from Lucy Weggler's demo notebook.

        Returns:
            Q_0: numpy array (n_star x n_loop), real
        """
        from ngsolve import ds, div
        from ngsolve.bem import LaplaceSL

        # Product space for BEM coupling operator
        fes_hdiv_r = type(self._fes_hdiv)(self.mesh, order=self.order)
        fes_l2_r = type(self._fes_l2)(self.mesh, order=self.order,
                                       dual_mapping=False)
        fes_prod = fes_hdiv_r * fes_l2_r
        (uHDiv, uL2), (vHDiv, vL2) = fes_prod.TnT()

        n_total = fes_prod.ndof

        Q_op = LaplaceSL(
            div(uHDiv.Trace()) * ds(bonus_intorder=self.intorder)
        ) * vL2 * ds(bonus_intorder=self.intorder)

        # Extract the full product-space matrix and get the (1,0) block
        Q_full = extract_dense_matrix(Q_op.mat, n_total)

        # Q_op maps HDivSurface trial [0:n_loop] to SurfaceL2 test [n_loop:]
        Q_0 = np.real(Q_full[self.n_loop:, :self.n_loop])

        return Q_0

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

        # Simplified resistance estimate: all edges get the same R value
        # based on average element dimensions. For a more accurate model,
        # compute per-edge length and dual area from mesh topology.
        total_area = np.sum(areas)
        mean_edge_length = np.sqrt(4.0 * total_area / (np.sqrt(3) * n_tri))

        # Each edge has two adjacent triangles with combined area ~ 2*A_mean
        # Effective width ~ 2*A_mean / l_edge
        A_mean = total_area / n_tri
        w_eff = 2.0 * A_mean / mean_edge_length

        R_per_edge = R_sheet * mean_edge_length / w_eff
        return np.full(self.n_loop, R_per_edge)

    def solve_frequency(self, freqs, mode='mqs', Zs_func=None):
        """Compute port impedance Z(f) over a frequency range.

        Args:
            freqs: Array of frequencies [Hz]
            mode: Solver mode:
                'mqs': Loop-only (Z = R + jwL), fastest, valid DC - 1 MHz
                'full': Full Loop-Star with reformulated Schur complement
                'stabilized': Weggler's stabilized BEM (most robust)
            Zs_func: Optional callable(freq) returning Zs array (n_loop,)
                     for frequency-dependent surface impedance (e.g. Dowell).
                     If None, only DC resistance is used.

        Returns:
            Z_port: Complex impedance array, shape (len(freqs),)
        """
        if self.L is None:
            raise RuntimeError("Call assemble() before solve_frequency()")

        freqs = np.asarray(freqs, dtype=float)
        Z_port = np.zeros(len(freqs), dtype=complex)

        for i, f in enumerate(freqs):
            omega = 2.0 * np.pi * f
            Zs = Zs_func(f) if Zs_func is not None else None

            if mode == 'mqs':
                Z_port[i] = self._solve_mqs(omega, Zs=Zs)
            elif mode == 'stabilized':
                Z_port[i] = self._solve_stabilized(omega, Zs=Zs)
            else:  # 'full'
                Z_port[i] = self._solve_full_loop_star(omega, Zs=Zs)

        return Z_port

    def _solve_mqs(self, omega, Zs=None):
        """MQS mode: Loop-only impedance (capacitive effects neglected).

        Z_branch = diag(R + Zs) + jw*L  (complex symmetric: A = A^T)
        Z_port = 1 / (e^T * Z_branch^{-1} * e)

        Uses LDL^T factorization (Bunch-Kaufman) for complex symmetric solve.

        Args:
            omega: Angular frequency [rad/s]
            Zs: Optional surface impedance array (n_loop,) complex.
                Added to DC resistance diagonal.
        """
        R_diag = self.R_loop.astype(complex)
        if Zs is not None:
            R_diag = R_diag + np.asarray(Zs, dtype=complex)
        Z_branch = np.diag(R_diag) + 1j * omega * self.L

        # Port vector: uniform excitation (all edges carry same current)
        e = np.ones(self.n_loop) / self.n_loop

        # Solve Z_branch * x = e via LDL^T (complex symmetric)
        x = scipy.linalg.solve(Z_branch, e, assume_a='sym')
        Z_port = 1.0 / (e @ x)

        return Z_port

    def _solve_full_loop_star(self, omega, Zs=None):
        """Full Loop-Star: includes capacitive (Star) DOFs.

        Reformulated Schur complement (numerically stable):
            Z_eff = Z_LL - jw * M_LS^T * P^{-1} * M_LS

        The key fix: instead of forming Z_SS = P/(jw) (blows up at low freq),
        we precompute P^{-1} * M_LS once and multiply by jw at each frequency.

        This gives the exact same mathematical result but avoids the
        O(1/omega) conditioning issue in the intermediate linear solve.
        """
        # Loop-Loop block
        R_diag = self.R_loop.astype(complex)
        if Zs is not None:
            R_diag = R_diag + np.asarray(Zs, dtype=complex)
        Z_LL = np.diag(R_diag) + 1j * omega * self.L

        if omega < 1e-10:
            # DC limit: only resistance
            e = np.ones(self.n_loop) / self.n_loop
            x = scipy.linalg.solve(np.diag(self.R_loop), e)
            return 1.0 / (e @ x) + 0j

        # Schur complement: Z_eff = Z_LL - jw * M_LS^T @ P^{-1} @ M_LS
        # P^{-1} @ M_LS is precomputed in assemble() as self._P_inv_M_LS
        if self._P_inv_M_LS is not None:
            cap_correction = (1j * omega) * (self.M_LS.T @ self._P_inv_M_LS)
        else:
            # Fallback: solve at each frequency (slower but always works)
            X = scipy.linalg.solve(self.P, self.M_LS, assume_a='sym')
            cap_correction = (1j * omega) * (self.M_LS.T @ X)

        Schur = Z_LL - cap_correction

        # Port extraction
        e = np.ones(self.n_loop) / self.n_loop
        x = scipy.linalg.solve(Schur, e, assume_a='sym')
        Z_port = 1.0 / (e @ x)

        return Z_port

    def _solve_stabilized(self, omega, Zs=None):
        """Stabilized formulation using BEM Q_0 coupling (Weggler).

        Uses Q_0 = LaplaceSL(div(J), rho) instead of local FEM M_LS.
        At MQS (k->0): the system is a saddle-point [L, Q_0; Q_0^T, 0]
        which enforces divergence-free current (no charge accumulation).

        For finite frequency: includes k^2*V_0 capacitive block.

        Reference:
            L. Weggler, "Stabilized DtN formulation for low-frequency BEM"
            https://github.com/Weggler/docu-ngsbem/blob/main/demos/
            Maxwell_DtN_Stabilized.ipynb
        """
        if self.Q_0 is None:
            raise RuntimeError("Q_0 not assembled. Call assemble() first.")

        R_diag = self.R_loop.astype(complex)
        if Zs is not None:
            R_diag = R_diag + np.asarray(Zs, dtype=complex)

        if omega < 1e-10:
            # DC limit: only resistance
            e = np.ones(self.n_loop) / self.n_loop
            x = scipy.linalg.solve(np.diag(self.R_loop), e)
            return 1.0 / (e @ x) + 0j

        # Build full stabilized block system
        n = self.n_loop + self.n_star
        Z_full = np.zeros((n, n), dtype=complex)

        # (1,1): Z_LL = diag(R+Zs) + jw*L
        Z_full[:self.n_loop, :self.n_loop] = (
            np.diag(R_diag) + 1j * omega * self.L)

        # Off-diagonal: mu_0 * Q_0 coupling (BEM, not local FEM)
        # Scale by jw*mu_0 to match EFIE formulation
        scale = 1j * omega * MU_0
        Z_full[:self.n_loop, self.n_loop:] = scale * self.Q_0.T
        Z_full[self.n_loop:, :self.n_loop] = scale * self.Q_0

        # (2,2): jw*mu_0*k^2*V_0 (saddle-point at k=0, stabilized at finite k)
        k = omega / C_0
        Z_full[self.n_loop:, self.n_loop:] = (
            scale * k**2 * self.V_0)

        # RHS: port excitation in loop DOFs
        rhs = np.zeros(n, dtype=complex)
        e = np.ones(self.n_loop) / self.n_loop
        rhs[:self.n_loop] = e

        # Solve full system (well-conditioned by stabilization)
        sol = scipy.linalg.solve(Z_full, rhs)

        # Port impedance from loop DOFs
        I_L = sol[:self.n_loop]
        Z_port = 1.0 / (e @ I_L)

        return Z_port

    def solve_loop_star(self, freq, mode='full'):
        """Solve full Loop-Star system and return solution vectors.

        This method exposes the internal current and charge distributions,
        useful for visualization and post-processing.

        Args:
            freq: Frequency [Hz]
            mode: 'full' (reformulated Schur) or 'stabilized' (Weggler BEM)

        Returns:
            I_loop: Complex array (n_loop,) — edge current coefficients
            Q_star: Complex array (n_star,) — cell charge coefficients
        """
        if self.L is None:
            raise RuntimeError("Call assemble() before solve_loop_star()")

        omega = 2.0 * np.pi * freq
        e = np.ones(self.n_loop) / self.n_loop

        if omega < 1e-10:
            # DC limit: only resistance, no charge
            I_loop = scipy.linalg.solve(np.diag(self.R_loop), e)
            Q_star = np.zeros(self.n_star, dtype=complex)
            return I_loop, Q_star

        if mode == 'stabilized':
            # Solve full stabilized system
            n = self.n_loop + self.n_star
            Z_full = np.zeros((n, n), dtype=complex)
            Z_full[:self.n_loop, :self.n_loop] = (
                np.diag(self.R_loop.astype(complex)) + 1j * omega * self.L)
            scale = 1j * omega * MU_0
            Z_full[:self.n_loop, self.n_loop:] = scale * self.Q_0.T
            Z_full[self.n_loop:, :self.n_loop] = scale * self.Q_0
            k = omega / C_0
            Z_full[self.n_loop:, self.n_loop:] = scale * k**2 * self.V_0

            rhs = np.zeros(n, dtype=complex)
            rhs[:self.n_loop] = e
            sol = scipy.linalg.solve(Z_full, rhs)
            return sol[:self.n_loop], sol[self.n_loop:]

        # mode='full': reformulated Schur complement
        Z_LL = np.diag(self.R_loop.astype(complex)) + 1j * omega * self.L

        # Schur complement: Z_eff = Z_LL - jw * M_LS^T @ P^{-1} @ M_LS
        if self._P_inv_M_LS is not None:
            cap_correction = (1j * omega) * (self.M_LS.T @ self._P_inv_M_LS)
        else:
            X = scipy.linalg.solve(self.P, self.M_LS, assume_a='sym')
            cap_correction = (1j * omega) * (self.M_LS.T @ X)

        Schur = Z_LL - cap_correction

        # Solve for loop currents
        I_loop = scipy.linalg.solve(Schur, e, assume_a='sym')

        # Back-substitute for star charges: Q = -jw * P^{-1} * M_LS * I
        Q_star = -(1j * omega) * scipy.linalg.solve(
            self.P, self.M_LS @ I_loop, assume_a='sym')

        return I_loop, Q_star

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
        print(f"  M_LS matrix ({self.n_star}x{self.n_loop}) [local FEM]:")
        print(f"    Range: [{np.min(self.M_LS):.4e}, {np.max(self.M_LS):.4e}]")
        print(f"    Nonzeros: {np.count_nonzero(np.abs(self.M_LS) > 1e-15)}"
              f" / {self.M_LS.size}")
        print()

        # Q_0 properties
        if self.Q_0 is not None:
            print(f"  Q_0 matrix ({self.Q_0.shape[0]}x{self.Q_0.shape[1]}) "
                  f"[BEM LaplaceSL]:")
            print(f"    Range: [{np.min(self.Q_0):.4e}, "
                  f"{np.max(self.Q_0):.4e}]")
            print(f"    Nonzeros: "
                  f"{np.count_nonzero(np.abs(self.Q_0) > 1e-15)}"
                  f" / {self.Q_0.size}")
            print()

        # Resistance
        print(f"  R_loop: [{np.min(self.R_loop):.4e}, "
              f"{np.max(self.R_loop):.4e}] Ohm")
        if self.thickness:
            R_sheet = 1.0 / (self.sigma * self.thickness)
            print(f"  Sheet resistance: {R_sheet:.4e} Ohm/sq")


# ============================================================
# Visualization functions (netgen.webgui)
# ============================================================

def draw_peec_model(mesh, solver=None, conductor_label="conductor"):
    """Draw PEEC model in netgen.webgui: surface mesh with edges and DOF counts.

    Renders the surface mesh (triangles) with wireframe overlay.
    If a solver is provided, prints the Loop/Star DOF counts.

    Args:
        mesh: NGSolve Mesh (surface mesh)
        solver: NGBEMPEECSolver instance (optional, for DOF info)
        conductor_label: Boundary label

    Returns:
        scene: webgui scene object
    """
    from ngsolve.webgui import Draw

    if solver is not None:
        print(f"PEEC model: {solver.n_loop} Loop DOFs (edges), "
              f"{solver.n_star} Star DOFs (cells)")

    scene = Draw(mesh)
    return scene


def draw_loop_current(solver, freq, mesh):
    """Visualize Loop current magnitude |I| on the mesh.

    Solves the full Loop-Star system at the given frequency and projects
    the edge current magnitudes onto a per-element scalar field.

    Each triangle's current is the mean |I| of its three edge DOFs.
    This shows where the current concentrates on the conductor surface.

    Args:
        solver: NGBEMPEECSolver instance (must be assembled)
        freq: Frequency [Hz]
        mesh: NGSolve Mesh

    Returns:
        gf_current: GridFunction with current magnitude
        scene: webgui scene object
    """
    from ngsolve import SurfaceL2, GridFunction, BND
    from ngsolve.webgui import Draw

    I_loop, Q_star = solver.solve_loop_star(freq)

    # Create scalar field for current magnitude
    fes_l2 = SurfaceL2(mesh, order=0)
    gf_current = GridFunction(fes_l2, name="Loop_current")

    # Map edge DOF magnitudes to element-average current
    for el in mesh.Elements(BND):
        edge_nrs = [e.nr for e in el.edges]
        # Average |I| of the edges belonging to this element
        edge_vals = [abs(I_loop[e]) for e in edge_nrs
                     if e < len(I_loop)]
        if edge_vals:
            gf_current.vec[el.nr] = np.mean(edge_vals)

    scene = Draw(gf_current, mesh, name="|I_loop|")
    return gf_current, scene


def draw_star_charge(solver, freq, mesh):
    """Visualize Star charge density Q on the mesh.

    Solves the full Loop-Star system at the given frequency and displays
    the per-cell charge density. Star DOFs are naturally piecewise constant
    on cells (SurfaceL2), so no projection is needed.

    Positive/negative charge accumulation shows the capacitive behavior
    of the conductor.

    Args:
        solver: NGBEMPEECSolver instance (must be assembled)
        freq: Frequency [Hz]
        mesh: NGSolve Mesh

    Returns:
        gf_charge: GridFunction with charge density
        scene: webgui scene object
    """
    from ngsolve import SurfaceL2, GridFunction
    from ngsolve.webgui import Draw

    I_loop, Q_star = solver.solve_loop_star(freq)

    # Star DOFs map directly to SurfaceL2 cells
    fes_l2 = SurfaceL2(mesh, order=0)
    gf_charge = GridFunction(fes_l2, name="Star_charge")

    for i in range(min(solver.n_star, fes_l2.ndof)):
        gf_charge.vec[i] = Q_star[i].real

    scene = Draw(gf_charge, mesh, name="Q_star")
    return gf_charge, scene


def draw_loop_star(solver, freq, mesh):
    """Visualize both Loop currents and Star charges side by side.

    Convenience wrapper that calls draw_loop_current() and
    draw_star_charge() and returns both results.

    Args:
        solver: NGBEMPEECSolver instance (must be assembled)
        freq: Frequency [Hz]
        mesh: NGSolve Mesh

    Returns:
        (gf_current, gf_charge): Tuple of GridFunctions
    """
    print(f"Loop-Star decomposition at f = {freq:.2e} Hz:")

    gf_current, _ = draw_loop_current(solver, freq, mesh)
    gf_charge, _ = draw_star_charge(solver, freq, mesh)

    # Print statistics
    I_loop, Q_star = solver.solve_loop_star(freq)
    Z_port = solver.solve_frequency(np.array([freq]), mode='full')[0]
    L_nH = np.imag(Z_port) / (2 * np.pi * freq) * 1e9

    print(f"  |I_loop| range: [{np.min(np.abs(I_loop)):.4e}, "
          f"{np.max(np.abs(I_loop)):.4e}]")
    print(f"  |Q_star| range: [{np.min(np.abs(Q_star)):.4e}, "
          f"{np.max(np.abs(Q_star)):.4e}]")
    print(f"  Z_port = {Z_port.real:.4e} + j{Z_port.imag:.4e} Ohm")
    print(f"  L_eff = {L_nH:.2f} nH")

    return gf_current, gf_charge


def compute_vector_potential(solver, I_loop, obs_points):
    """Compute vector potential A at observation points from solved current J.

    This is the key output of the PEEC solver for coupling to external physics:
    a coil's effect on a magnetic core, a shield's back-action on a coil, etc.

        A(r) = mu_0 / (4*pi) * integral_S  J(r') / |r - r'|  dS'

    Uses NGSolve Integrate over boundary elements with the solved
    HDivSurface current distribution.

    Args:
        solver: NGBEMPEECSolver instance (must be assembled)
        I_loop: Complex array (n_loop,) of edge current coefficients
                (from solve_loop_star or solve_frequency)
        obs_points: Array of observation points, shape (N, 3) [m]

    Returns:
        A: Complex array (N, 3) — vector potential [T*m] at each point
    """
    from ngsolve import x, y, z, sqrt, Integrate, GridFunction

    mesh = solver.mesh
    bnd = mesh.Boundaries(solver.conductor_label)

    obs_points = np.asarray(obs_points)
    n_obs = len(obs_points)

    # Process real and imaginary parts separately (Integrate is real-valued)
    A = np.zeros((n_obs, 3), dtype=complex)

    for part_label, get_part in [('real', np.real), ('imag', np.imag)]:
        gf_J = GridFunction(solver._fes_hdiv)
        for i in range(solver.n_loop):
            gf_J.vec[i] = float(get_part(I_loop[i]))

        for i in range(n_obs):
            rx, ry, rz = obs_points[i]
            # Distance with small regularization (obs points should be off-surface)
            dist = sqrt((x - rx)**2 + (y - ry)**2 + (z - rz)**2 + 1e-30)
            kernel = 1.0 / dist

            for comp in range(3):
                val = MU_0 / (4 * np.pi) * Integrate(
                    gf_J[comp] * kernel, mesh, definedon=bnd, order=5)

                if part_label == 'real':
                    A[i, comp] += val
                else:
                    A[i, comp] += 1j * val

    return A
