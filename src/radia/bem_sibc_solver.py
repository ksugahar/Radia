"""
Scalar potential BIE with SIBC for MQS eddy current on conducting surfaces.

Solves for the scalar magnetic potential phi on a conductor surface:
    (1/2*M - DL + gamma * SL * M^{-1} * K) phi = rhs
    gamma = Z_s / (jw * mu_0)

where M = H1 surface mass, K = H1 surface stiffness (Laplace-Beltrami),
DL = Laplace double layer, SL = Laplace single layer (all from ngsolve.bem).

Usage:
    from radia.bem_sibc_solver import ScalarBIESIBCSolver

    solver = ScalarBIESIBCSolver(mesh, order=1)
    result = solver.solve(phi_inc_cf, Z_s=Z_s, omega=omega)

    J_rms = result['H_t_rms']       # Surface current RMS [A/m]
    P = result['P_density']          # Power loss density [W/m^2]
    phi_gf = result['phi']           # Solved potential GridFunction
"""

import math
import time
import numpy as np
from scipy.linalg import solve as scipy_solve

MU_0 = 4e-7 * np.pi


class ScalarBIESIBCSolver:
    """Scalar potential BIE + SIBC solver for conducting surfaces.

    The solver assembles BEM operators once, then can be called repeatedly
    with different Z_s (frequency sweep) or phi_inc (different sources).
    """

    def __init__(self, mesh, order=1, assemble_dense=True,
                  use_intree_bem=False, intree_geom_order=2,
                  intree_singular_n_q=8, intree_regular_quad_degree=11,
                  use_intree_hacapk=False, hacapk_aca_eps=1e-10,
                  hacapk_leaf=64, hacapk_eta=2.0,
                  bnd_label=None, log_fn=None):
        """Initialize solver and assemble BEM operators.

        Args:
            mesh: NGSolve Mesh (surface mesh, dim=2, or volume mesh with BND).
                When ``bnd_label`` is given, ``mesh`` is interpreted as a
                volume (or larger surface) mesh, and only the BND elements
                whose boundary label equals ``bnd_label`` are used.  This
                is the path used for curved scalar SIBC BEM: pass the
                parent ``vol_mesh`` (after ``vol_mesh.Curve(p)``) plus the
                workpiece sideset name, and the curving is preserved via
                ``mesh.GetTrafo(el)`` on the parent.
            order: H1 polynomial order on surface (default 1).
            assemble_dense: if True (default, backward-compat), extract
                ``DL`` and ``SL`` to dense ``ndof x ndof`` numpy arrays
                via N column matvecs.  This is O(N^3) total (~67 min for
                N=2474) but lets ``solve()`` use a single dense scipy
                solve.
                If False, KEEP ``DL_bf`` and ``SL_bf`` as NGSolve
                bilinear forms only and use ``solve_iterative()`` (GMRES
                with LinearOperator wrappers).  ~50x faster for typical
                IH wp meshes -- recommended for new code.
            bnd_label: optional NGSolve boundary label name to filter BND
                element by.  Required when ``mesh`` is the parent volume
                mesh and only a subset of its boundary belongs to the
                workpiece (e.g. ``bnd_label="sibc"`` selects the SIBC
                sideset only).  Only consumed by the in-tree paths
                (``use_intree_bem=True``); the ngsolve.bem path uses
                ``ds`` over all BND.
            log_fn: optional progress callback ``log_fn(tag, msg)`` used
                to surface assembly-phase boundaries to the panel debug
                log + console.  When None (default) the solver is silent
                -- backwards compatible for non-panel callers / tests.
                Caller should pass ``calc_common.progress`` for the IH
                panel pipeline.  Phases logged: SLDL assembly (C++ vs
                Python), mass / stiffness assembly, HACApK compression.
        """
        from ngsolve import (H1, BilinearForm, GridFunction, ds, grad,
                             TaskManager, InnerProduct)

        # Phase log (default no-op so non-panel callers stay silent).
        if log_fn is None:
            def _log_phase(_tag, _msg):
                pass
        else:
            _log_phase = log_fn

        self.mesh = mesh
        self.order = order
        self.use_intree_bem = use_intree_bem
        self._bnd_label = bnd_label
        # Lagrange-P2 in-tree path is taken when use_intree_bem and order >= 2.
        # In that mode self.fes is None (we don't use NGSolve FES at all):
        # SL, DL, M, K are all assembled in Lagrange P2 basis directly
        # (NGSolve uses hierarchical Lobatto for H1 order>=2; mixing the
        # two would introduce a basis-mismatch in the assembled SL/DL).
        self._intree_lagrange_p2 = bool(use_intree_bem and order >= 2)
        self.dof_coords = None   # populated for the Lagrange-P2 path

        if self._intree_lagrange_p2:
            # Skip NGSolve FES creation entirely.
            self.fes = None
            t0 = time.perf_counter()
            _log_phase("BEM",
                f"Lagrange-P2 path: bnd_label={bnd_label!r}, "
                f"geom_order={intree_geom_order}")
            self._init_lagrange_p2_path(
                mesh, order, intree_geom_order,
                intree_singular_n_q, intree_regular_quad_degree,
                use_intree_hacapk, hacapk_aca_eps, hacapk_leaf, hacapk_eta,
                bnd_label=bnd_label, log_fn=log_fn)
            self.t_assembly = time.perf_counter() - t0
            _log_phase("BEM",
                f"Lagrange-P2 path done (ndof={self.ndof}, "
                f"{self.t_assembly:.1f}s)")
            return

        # H1 on surface (P1, or order>=2 with use_intree_bem=False)
        self.fes = H1(mesh, order=order)
        u, v = self.fes.TnT()
        self.ndof = self.fes.ndof

        t0 = time.perf_counter()

        ndof = self.ndof

        if use_intree_bem:
            # In-tree Sauter-Schwab Galerkin BEM (no ngsolve.bem dep).
            # Uses the C++ assembler (radia._radia_pybind._AssembleSLDL_
            # Galerkin) when available -- pure-Python fallback only when
            # the C++ symbol is missing (older wheel without Phase 1.9).
            from radia.bem.sibc_hacapk import extract_surface_curved
            _t_ext = time.perf_counter()
            verts, tris, v_global, tri_p2 = extract_surface_curved(
                mesh, geom_order=intree_geom_order)
            _log_phase("BEM",
                f"extract_surface_curved: {len(tris)} tris, {len(verts)} verts "
                f"(geom_order={intree_geom_order}, "
                f"{time.perf_counter()-_t_ext:.1f}s)")
            try:
                from radia import _radia_pybind as _rpb
                _cpp_assemble = _rpb._AssembleSLDL_Galerkin
            except (ImportError, AttributeError):
                _cpp_assemble = None
            _t_sldl = time.perf_counter()
            if _cpp_assemble is not None:
                _log_phase("BEM",
                    f"SLDL Galerkin assembly: C++ kernel, "
                    f"{len(tris)} tris, quad_deg={intree_regular_quad_degree}, "
                    f"sing_n_q={intree_singular_n_q}")
                v_arr = np.ascontiguousarray(verts, dtype=np.float64)
                t_arr = np.ascontiguousarray(tris, dtype=np.int64)
                p_arr = np.ascontiguousarray(tri_p2, dtype=np.float64)
                SL_loc, DL_loc = _cpp_assemble(
                    v_arr, t_arr, p_arr,
                    intree_regular_quad_degree,
                    intree_singular_n_q,
                    0)   # n_threads=0 -> OpenMP default
            else:
                _log_phase("BEM",
                    f"SLDL Galerkin assembly: pure-Python fallback "
                    f"(C++ _AssembleSLDL_Galerkin not available -- slow!), "
                    f"{len(tris)} tris")
                from radia.bem.sibc_hacapk import (
                    assemble_SL_dense_curved, assemble_DL_dense_curved)
                SL_loc = assemble_SL_dense_curved(
                    verts, tris, tri_p2,
                    regular_quad_degree=intree_regular_quad_degree,
                    singular_n_q=intree_singular_n_q)
                DL_loc = assemble_DL_dense_curved(
                    verts, tris, tri_p2,
                    regular_quad_degree=intree_regular_quad_degree,
                    singular_n_q=intree_singular_n_q)
            _log_phase("BEM",
                f"SLDL Galerkin assembly done "
                f"({time.perf_counter()-_t_sldl:.1f}s)")
            # Lift to full ndof basis (interior vertices contribute zero
            # rows/cols since their hat is 0 on BND).
            self.SL = np.zeros((ndof, ndof))
            self.DL = np.zeros((ndof, ndof))
            self.SL[np.ix_(v_global, v_global)] = SL_loc
            self.DL[np.ix_(v_global, v_global)] = DL_loc
            self._intree_v_global = v_global
            self._DL_bf = None
            self._SL_bf = None

            # Optional: compress SL and DL via HACApK for fast MatVec.
            # When use_intree_hacapk=True the compressed handles are
            # built and stored as self._SL_hacapk / self._DL_hacapk; the
            # caller can then choose between solve() (dense LU) and
            # solve_hacapk() (GMRES via H-matrix MatVec).  At N<5000 on
            # compact wp geometries (~80%+ near-field) compression gives
            # essentially no speedup, but at N>10K or for frequency
            # sweeps with reused operators the compressed MatVec wins.
            self._SL_hacapk = None
            self._DL_hacapk = None
            if use_intree_hacapk:
                _t_hca = time.perf_counter()
                _log_phase("BEM",
                    f"HACApK compress: SL + DL, ndof={ndof}, "
                    f"aca_eps={hacapk_aca_eps}, leaf={int(hacapk_leaf)}, "
                    f"eta={hacapk_eta}")
                from radia import _radia_pybind as _rpb_h
                # Build coordinate array sized for the FULL ndof; for
                # interior vertices that don't appear in BND we use the
                # MESH vertex coordinate (the H-matrix clustering only
                # uses these for spatial bisection, never for kernel).
                coords_full = np.zeros((ndof, 3), dtype=np.float64)
                for i in range(ndof):
                    coords_full[i] = mesh.vertices[i].point
                coords_full = np.ascontiguousarray(coords_full)
                SL_arr = np.ascontiguousarray(self.SL)
                DL_arr = np.ascontiguousarray(self.DL)
                self._SL_hacapk = _rpb_h.HACApKBEMManager(coords_full, SL_arr)
                self._SL_hacapk.BuildHMatrix(
                    aca_eps=hacapk_aca_eps,
                    leaf_size=int(hacapk_leaf),
                    eta=hacapk_eta,
                    max_rank=-1, print_level=0)
                self._DL_hacapk = _rpb_h.HACApKBEMManager(coords_full, DL_arr)
                self._DL_hacapk.BuildHMatrix(
                    aca_eps=hacapk_aca_eps,
                    leaf_size=int(hacapk_leaf),
                    eta=hacapk_eta,
                    max_rank=-1, print_level=0)
                _log_phase("BEM",
                    f"HACApK compress done "
                    f"({time.perf_counter()-_t_hca:.1f}s)")
        else:
            from ngsolve.bem import LaplaceDL, LaplaceSL
            # BEM operators.  Probed 2026-04-29 on N=2477 wp:
            # use_fmm=True/False both run at ~1.67 s/matvec; default
            # H-matrix used.
            # NOTE: Caller MUST be inside `with TaskManager():` per
            # CLAUDE.md "Caller Wraps, Helper Does NOT" (2026-05-27).
            DL_bf = LaplaceDL(u.Trace() * ds) * v.Trace() * ds
            SL_bf = LaplaceSL(u.Trace() * ds) * v.Trace() * ds
            self._DL_bf = DL_bf
            self._SL_bf = SL_bf

            if assemble_dense:
                self.DL = np.zeros((ndof, ndof))
                self.SL = np.zeros((ndof, ndof))
                for j in range(ndof):
                    ej = GridFunction(self.fes)
                    ej.vec[:] = 0
                    ej.vec[j] = 1.0
                    r1 = ej.vec.CreateVector()
                    r1.data = DL_bf.mat * ej.vec
                    self.DL[:, j] = r1.FV().NumPy().copy()
                    r2 = ej.vec.CreateVector()
                    r2.data = SL_bf.mat * ej.vec
                    self.SL[:, j] = r2.FV().NumPy().copy()
            else:
                self.DL = None
                self.SL = None

        # Surface mass M
        _t_mk = time.perf_counter()
        _log_phase("BEM", f"mass + stiffness assembly (ndof={ndof})")
        mass_bf = BilinearForm(self.fes)
        mass_bf += u.Trace() * v.Trace() * ds
        mass_bf.Assemble()
        self.M = np.zeros((ndof, ndof))
        rows, cols, vals = mass_bf.mat.COO()
        for r_, c_, val in zip(rows, cols, vals):
            self.M[int(r_), int(c_)] = val

        # Surface stiffness K (Laplace-Beltrami)
        stiff_bf = BilinearForm(self.fes)
        stiff_bf += InnerProduct(grad(u).Trace(), grad(v).Trace()) * ds
        stiff_bf.Assemble()
        self.K = np.zeros((ndof, ndof))
        rows, cols, vals = stiff_bf.mat.COO()
        for r_, c_, val in zip(rows, cols, vals):
            self.K[int(r_), int(c_)] = val

        _log_phase("BEM",
            f"mass + stiffness done; M_inv "
            f"({time.perf_counter()-_t_mk:.1f}s)")
        self.M_inv = np.linalg.inv(self.M)

        # Gauge vector: <1, v>_S for Lagrange multiplier
        self._c_gauge = self.M @ np.ones(ndof)

        self.t_assembly = time.perf_counter() - t0

    def _init_lagrange_p2_path(self, mesh, order,
                                intree_geom_order,
                                intree_singular_n_q,
                                intree_regular_quad_degree,
                                use_intree_hacapk,
                                hacapk_aca_eps, hacapk_leaf, hacapk_eta,
                                bnd_label=None, log_fn=None):
        """Build SL, DL, M, K, gauge entirely in Lagrange P2 basis.

        Bypasses NGSolve's H1 hierarchical FES so that all matrices are
        consistent in the same basis -- the user passes phi_inc as an
        ndarray of length ``ndof`` whose i-th entry is the scalar
        potential value at the i-th Lagrange P2 node (vertex or edge
        mid-point).  ``self.dof_coords`` exposes those node positions
        so the caller can sample any analytical phi at the right places.

        When ``bnd_label`` is given, only BND elements with that label are
        included.  Combined with passing the parent (already curved)
        volume mesh as ``mesh``, this gives a true curved-Tri6 P2
        scalar BEM on the workpiece surface.
        """
        if order != 2:
            raise NotImplementedError(
                f"Lagrange-P2 in-tree path supports order=2 only "
                f"(got order={order}).  Higher H1 orders need a generalized "
                f"Lagrange basis assembler.")
        from radia.bem.sibc_hacapk import (
            extract_surface_p2_lagrange,
            assemble_SL_dense_curved_p2, assemble_DL_dense_curved_p2,
            assemble_mass_curved_p2, assemble_stiffness_curved_p2,
        )

        verts, tris, v_global, tri_p2, dofs_per_tri, n_dof, dof_coords = \
            extract_surface_p2_lagrange(mesh, bnd_label=bnd_label,
                                        geom_order=intree_geom_order)

        self.ndof = n_dof
        self.dof_coords = dof_coords
        self._intree_v_global = v_global
        self._intree_dofs_per_tri = dofs_per_tri
        self._intree_tri_p2 = tri_p2

        # SL / DL: Python implementation (slow O(n_t^2)).  Replaced with a
        # C++ binding ``_AssembleSLDL_Galerkin_P2`` once the .pyd ships
        # the P2 entry point.
        try:
            from radia import _radia_pybind as _rpb
            _cpp_assemble_p2 = getattr(_rpb, "_AssembleSLDL_Galerkin_P2",
                                       None)
        except ImportError:
            _cpp_assemble_p2 = None
        if _cpp_assemble_p2 is not None:
            v_arr = np.ascontiguousarray(verts, dtype=np.float64)
            t_arr = np.ascontiguousarray(tris, dtype=np.int64)
            p_arr = np.ascontiguousarray(tri_p2, dtype=np.float64)
            d_arr = np.ascontiguousarray(dofs_per_tri, dtype=np.int64)
            self.SL, self.DL = _cpp_assemble_p2(
                v_arr, t_arr, p_arr, d_arr, int(n_dof),
                int(intree_regular_quad_degree),
                int(intree_singular_n_q),
                0)
        else:
            self.SL = assemble_SL_dense_curved_p2(
                verts, tris, tri_p2, dofs_per_tri, n_dof,
                regular_quad_degree=intree_regular_quad_degree,
                singular_n_q=intree_singular_n_q)
            self.DL = assemble_DL_dense_curved_p2(
                verts, tris, tri_p2, dofs_per_tri, n_dof,
                regular_quad_degree=intree_regular_quad_degree,
                singular_n_q=intree_singular_n_q)

        self.M = assemble_mass_curved_p2(
            tri_p2, dofs_per_tri, n_dof,
            regular_quad_degree=intree_regular_quad_degree)
        self.K = assemble_stiffness_curved_p2(
            tri_p2, dofs_per_tri, n_dof,
            regular_quad_degree=intree_regular_quad_degree)

        self.M_inv = np.linalg.inv(self.M)
        self._c_gauge = self.M @ np.ones(n_dof)

        # HACApK compression: feed the manager `dof_coords` (one position
        # per Lagrange DOF -- vertex or edge midpoint).  HACApK is
        # basis-agnostic; it only uses the coords for spatial bisection
        # of the cluster tree.
        self._SL_hacapk = None
        self._DL_hacapk = None
        if use_intree_hacapk:
            from radia import _radia_pybind as _rpb_h
            coords_full = np.ascontiguousarray(dof_coords, dtype=np.float64)
            SL_arr = np.ascontiguousarray(self.SL)
            DL_arr = np.ascontiguousarray(self.DL)
            self._SL_hacapk = _rpb_h.HACApKBEMManager(coords_full, SL_arr)
            self._SL_hacapk.BuildHMatrix(
                aca_eps=hacapk_aca_eps,
                leaf_size=int(hacapk_leaf),
                eta=hacapk_eta,
                max_rank=-1, print_level=0)
            self._DL_hacapk = _rpb_h.HACApKBEMManager(coords_full, DL_arr)
            self._DL_hacapk.BuildHMatrix(
                aca_eps=hacapk_aca_eps,
                leaf_size=int(hacapk_leaf),
                eta=hacapk_eta,
                max_rank=-1, print_level=0)

    def solve(self, phi_inc_cf, Z_s, omega):
        """Solve scalar BIE + SIBC for given incident potential and impedance.

        Args:
            phi_inc_cf: NGSolve CoefficientFunction for incident scalar potential
                        on the surface, OR ndarray of length ndof (nodal values).
            Z_s: Surface impedance.
                 - **complex scalar**: legacy global SIBC (uniform Z_s
                   over the workpiece)
                 - **ndarray of length ndof (complex)**: per-node Z_s
                   (per-panel curvature SIBC; build the array from
                   panel-level Z_s values via vertex averaging or any
                   other H1 projection).
            omega: Angular frequency [rad/s].

        Returns:
            dict with keys:
                phi: GridFunction (complex) - solved surface potential
                phi_vec: ndarray (complex) - coefficient vector
                H_t_rms: float - RMS tangential H [A/m] (= surface current density)
                P_density: float - time-averaged power loss density [W/m^2]
                gamma: complex - Z_s / (jw * mu_0) (or its mean if per-node)
                t_solve: float - solve time [s]
        """
        if not self._intree_lagrange_p2:
            from ngsolve import (LinearForm, GridFunction, Integrate, CF, ds,
                                 grad, BND, InnerProduct)

        t0 = time.perf_counter()
        ndof = self.ndof

        # RHS: <phi_inc, v>_S
        if isinstance(phi_inc_cf, np.ndarray):
            rhs_vec = self.M @ phi_inc_cf
        elif self._intree_lagrange_p2:
            raise TypeError(
                "Lagrange-P2 path requires phi_inc as an ndarray of length "
                f"ndof={ndof}.  Sample your CF at solver.dof_coords (the "
                "physical positions of the Lagrange P2 nodes -- vertices + "
                "edge mid-points) and pass the resulting array.")
        else:
            v_h1 = self.fes.TestFunction()
            lf = LinearForm(self.fes)
            lf += phi_inc_cf * v_h1.Trace() * ds
            lf.Assemble()
            rhs_vec = lf.vec.FV().NumPy().copy()

        # System matrix.
        # Original global-Z_s formulation:
        #   gamma = Z_s / (jw * mu_0)         (scalar)
        #   A = 1/2 M - DL + gamma * (SL M^{-1} K)
        # Per-node Z_s formulation:
        #   Gamma_ii = Z_s_i / (jw * mu_0)    (diagonal)
        #   The Robin term gamma * Delta_s phi has the per-node coefficient
        #   on the test side: <gamma_i grad_s phi, grad_s v> -> the SL M^-1 K
        #   block is left-multiplied by diag(gamma) (per-row scaling).
        # That is, the discrete form (1/2 M - DL) phi + diag(gamma) (SL M^-1 K)
        # phi = M phi_inc.
        if isinstance(Z_s, np.ndarray):
            if Z_s.shape != (ndof,):
                raise ValueError(
                    f"Per-node Z_s must have shape ({ndof},), got {Z_s.shape}")
            if omega <= 0:
                gamma_vec = np.zeros(ndof, dtype=complex)
            else:
                gamma_vec = Z_s.astype(complex) / (1j * omega * MU_0)
            gamma_for_log = complex(np.mean(gamma_vec))
            # diag(gamma_vec) @ M = gamma_vec[:, None] * M (row-scaling)
            robin_block = gamma_vec[:, None] * (self.SL @ self.M_inv @ self.K)
            A_sys = (0.5 * self.M - self.DL + robin_block).astype(complex)
        else:
            gamma = Z_s / (1j * omega * MU_0) if omega > 0 and Z_s != 0 else 0
            gamma_for_log = complex(gamma)
            A_sys = (0.5 * self.M - self.DL
                     + gamma * self.SL @ self.M_inv @ self.K).astype(complex)

        # Solve with gauge (Lagrange multiplier for int phi dS = 0)
        phi_vec = self._solve_with_gauge(A_sys, rhs_vec.astype(complex))
        t_solve = time.perf_counter() - t0

        # Extract H_t_rms = sqrt(<|grad_s phi|^2> / area)
        if self._intree_lagrange_p2:
            # Use our Lagrange-P2 stiffness K to compute the surface-
            # gradient norm of phi:  ||grad_S phi||^2_L2 = phi^T K phi
            # (works for both real and imaginary parts).  The total area
            # is 1^T M 1.
            Hsq_re = float(phi_vec.real @ (self.K @ phi_vec.real))
            Hsq_im = float(phi_vec.imag @ (self.K @ phi_vec.imag))
            area = float(np.ones(ndof) @ (self.M @ np.ones(ndof)))
            H_t_rms = math.sqrt((abs(Hsq_re) + abs(Hsq_im)) / abs(area))
        else:
            gf = GridFunction(self.fes)
            gf.vec.FV().NumPy()[:] = phi_vec.real
            Hsq_re = Integrate(InnerProduct(grad(gf), grad(gf)),
                               self.mesh, BND)
            gf.vec.FV().NumPy()[:] = phi_vec.imag
            Hsq_im = Integrate(InnerProduct(grad(gf), grad(gf)),
                               self.mesh, BND)
            area = Integrate(CF(1), self.mesh, BND)
            H_t_rms = math.sqrt((abs(Hsq_re) + abs(Hsq_im)) / abs(area))

        # Power density: P' = (1/2) Re(Z_s) |J_s|^2 = (1/2) Re(Z_s) H_t_rms^2
        # (time-averaged). For per-node Z_s, the area-averaged Re(Z_s) is
        # the right scalar to report; the caller already integrates the
        # local power per panel via the panel-level R if needed.
        if isinstance(Z_s, np.ndarray):
            Z_s_avg_re = float(np.mean(Z_s.real))
            P_density = 0.5 * Z_s_avg_re * H_t_rms ** 2
        else:
            P_density = 0.5 * Z_s.real * H_t_rms ** 2 if Z_s != 0 else 0

        # GridFunction output (None for the Lagrange-P2 path)
        if self._intree_lagrange_p2:
            gf_phi = None
        else:
            gf_phi = GridFunction(self.fes)
            gf_phi.vec.FV().NumPy()[:] = phi_vec.real  # real part

        return {
            'phi': gf_phi,
            'phi_vec': phi_vec,
            'H_t_rms': float(H_t_rms),
            'P_density': float(P_density),
            'area': float(abs(area)),
            'gamma': gamma_for_log,
            't_solve': round(t_solve, 3),
        }

    def solve_hacapk(self, phi_inc_cf, Z_s, omega, *,
                       tol=1e-8, maxiter=300, restart=50):
        """SIBC solve via scipy GMRES + HACApK H-matrix MatVec.

        Requires use_intree_bem=True and use_intree_hacapk=True at
        construction.  Wraps the in-tree HACApK-compressed SL and DL
        operators in a scipy LinearOperator and solves the SIBC
        scalar BIE without ever assembling A_sys densely.

        For compact workpiece meshes at N ~ 2-5K, this is typically
        comparable in speed to dense LU because the H-matrix has
        near-zero compression (>80% near-field pairs).  For N > 10K
        or frequency sweeps reusing the same SL/DL it can be faster.

        Args / Returns: same as solve().
        """
        if not self._intree_lagrange_p2:
            from ngsolve import (LinearForm, GridFunction, Integrate, CF, ds,
                                 grad, BND, InnerProduct)
        from scipy.sparse.linalg import LinearOperator, gmres
        if self._SL_hacapk is None or self._DL_hacapk is None:
            raise RuntimeError(
                "solve_hacapk: HACApK handles not built.  Construct the "
                "solver with use_intree_bem=True and use_intree_hacapk=True.")

        t0 = time.perf_counter()
        ndof = self.ndof

        # RHS: <phi_inc, v>_S
        if isinstance(phi_inc_cf, np.ndarray):
            rhs_vec = self.M @ phi_inc_cf
        elif self._intree_lagrange_p2:
            raise TypeError(
                "Lagrange-P2 path requires phi_inc as an ndarray of length "
                f"ndof={ndof}.  Sample your CF at solver.dof_coords.")
        else:
            v_h1 = self.fes.TestFunction()
            lf = LinearForm(self.fes)
            lf += phi_inc_cf * v_h1.Trace() * ds
            lf.Assemble()
            rhs_vec = lf.vec.FV().NumPy().copy()

        # Accept scalar OR per-node ndarray Z_s.  The MatVec below
        # uses element-wise multiplication of the SL block by
        # ``gr_vec`` / ``gi_vec`` (length-ndof), so the scalar case
        # is handled by broadcasting a full-length array with the
        # same value at every DOF -- no separate code path.
        if isinstance(Z_s, np.ndarray):
            if Z_s.shape != (ndof,):
                raise ValueError(
                    f"Per-node Z_s must have shape ({ndof},), "
                    f"got {Z_s.shape}")
            if omega <= 0:
                gamma_vec = np.zeros(ndof, dtype=complex)
            else:
                gamma_vec = Z_s.astype(complex) / (1j * omega * MU_0)
            gamma_for_log = complex(np.mean(gamma_vec))
        else:
            gamma_scalar = (Z_s / (1j * omega * MU_0)
                             if (omega > 0 and Z_s != 0) else 0 + 0j)
            gamma_vec = np.full(ndof, gamma_scalar, dtype=complex)
            gamma_for_log = complex(gamma_scalar)
        gr_vec = np.ascontiguousarray(gamma_vec.real)
        gi_vec = np.ascontiguousarray(gamma_vec.imag)

        # Pre-compute M^{-1} @ K (real) once for reuse in MatVec.
        # The Robin term `SL @ M_inv @ K` applied to x is:
        #   (SL @ M_inv @ K) @ x = SL @ (M_inv @ K @ x)
        # so caching M_inv @ K saves one dense matmul per MatVec.
        Minv_K = self.M_inv @ self.K  # ndof x ndof, real

        SL_op = self._SL_hacapk
        DL_op = self._DL_hacapk
        M_full = self.M

        def matvec_complex(x_complex):
            """y = A_sys @ x where A_sys = (1/2)M - DL + diag(gamma)*SL*M^-1*K.

            Per-node Z_s is implemented by element-wise scaling of the
            SL block's output rows by gr_vec/gi_vec; scalar Z_s reduces
            to the same path with gamma_vec broadcast to a constant.
            """
            x_re = np.ascontiguousarray(x_complex.real)
            x_im = np.ascontiguousarray(x_complex.imag)
            # (1/2)*M - DL part (real)
            half_M_re = 0.5 * (M_full @ x_re)
            half_M_im = 0.5 * (M_full @ x_im)
            DL_re = DL_op.MatVec(x_re)
            DL_im = DL_op.MatVec(x_im)
            # SL @ M^-1 @ K @ x part (real intermediate)
            MinvK_re = Minv_K @ x_re
            MinvK_im = Minv_K @ x_im
            SL_re = SL_op.MatVec(np.ascontiguousarray(MinvK_re))
            SL_im = SL_op.MatVec(np.ascontiguousarray(MinvK_im))
            # diag(gamma) * (SL_re + j*SL_im) row-wise:
            # term[i] = (gr_i*SL_re[i] - gi_i*SL_im[i]) + j(gr_i*SL_im[i] + gi_i*SL_re[i])
            term_re = gr_vec * SL_re - gi_vec * SL_im
            term_im = gr_vec * SL_im + gi_vec * SL_re
            y_re = half_M_re - DL_re + term_re
            y_im = half_M_im - DL_im + term_im
            return y_re + 1j * y_im

        # Augment with gauge multiplier
        c_gauge = self._c_gauge

        def matvec_aug(x_aug):
            x = x_aug[:ndof]
            mu = x_aug[ndof]
            y = matvec_complex(x)
            y = y + mu * c_gauge.astype(complex)
            y_aug = np.empty(ndof + 1, dtype=complex)
            y_aug[:ndof] = y
            y_aug[ndof] = c_gauge @ x
            return y_aug

        A_aug = LinearOperator((ndof + 1, ndof + 1), matvec=matvec_aug,
                                dtype=complex)
        rhs_aug = np.zeros(ndof + 1, dtype=complex)
        rhs_aug[:ndof] = rhs_vec
        sol_aug, info = gmres(A_aug, rhs_aug, rtol=tol,
                                atol=0.0, maxiter=maxiter, restart=restart)
        if info != 0:
            print(f"  WARN: gmres did not converge cleanly (info={info})",
                  flush=True)
        phi_vec = sol_aug[:ndof]
        t_solve = time.perf_counter() - t0

        # Compute H_t_rms / P_density (same as solve())
        if self._intree_lagrange_p2:
            Hsq_re = float(phi_vec.real @ (self.K @ phi_vec.real))
            Hsq_im = float(phi_vec.imag @ (self.K @ phi_vec.imag))
            area = float(np.ones(ndof) @ (self.M @ np.ones(ndof)))
            H_t_rms = math.sqrt((abs(Hsq_re) + abs(Hsq_im)) / abs(area))
            gf_phi = None
        else:
            gf = GridFunction(self.fes)
            gf.vec.FV().NumPy()[:] = phi_vec.real
            Hsq_re = Integrate(InnerProduct(grad(gf), grad(gf)),
                               self.mesh, BND)
            gf.vec.FV().NumPy()[:] = phi_vec.imag
            Hsq_im = Integrate(InnerProduct(grad(gf), grad(gf)),
                               self.mesh, BND)
            area = Integrate(CF(1), self.mesh, BND)
            H_t_rms = math.sqrt((abs(Hsq_re) + abs(Hsq_im)) / abs(area))
            gf_phi = GridFunction(self.fes)
            gf_phi.vec.FV().NumPy()[:] = phi_vec.real
        # Per-node Z_s: report the area-averaged Re(Z_s) for the scalar
        # power-density quantity, matching the dense path's convention.
        if isinstance(Z_s, np.ndarray):
            Z_s_avg_re = float(np.mean(Z_s.real))
            P_density = 0.5 * Z_s_avg_re * H_t_rms ** 2
        else:
            P_density = 0.5 * Z_s.real * H_t_rms ** 2 if Z_s != 0 else 0.0
        return {
            'phi': gf_phi,
            'phi_vec': phi_vec,
            'H_t_rms': float(H_t_rms),
            'P_density': float(P_density),
            'area': float(abs(area)),
            'gamma': complex(gamma_for_log),
            't_solve': round(t_solve, 3),
        }

    def _solve_with_gauge(self, A_mat, rhs):
        """Solve A*phi = rhs with gauge constraint int(phi)dS = 0."""
        n = len(rhs)
        A_aug = np.zeros((n + 1, n + 1), dtype=complex)
        A_aug[:n, :n] = A_mat
        A_aug[:n, n] = self._c_gauge
        A_aug[n, :n] = self._c_gauge
        rhs_aug = np.zeros(n + 1, dtype=complex)
        rhs_aug[:n] = rhs
        return scipy_solve(A_aug, rhs_aug)[:n]

    def solve_iterative(self, phi_inc_cf, Z_s, omega, *,
                         tol=1e-8, maxiter=200, restart=50):
        """GMRES-based solve that AVOIDS dense DL/SL extraction.

        Drop-in replacement for ``solve()`` with ~50x speedup on the
        typical IH wp mesh (N=2000-5000 dofs).  The expensive
        column-by-column matvec extraction in ``__init__`` is skipped
        (you pass ``assemble_dense=False``) and the BEM operators
        ``DL_bf`` and ``SL_bf`` are wrapped as scipy ``LinearOperator``
        objects for matrix-free GMRES.

        Math is identical to ``solve()``:
            A_sys = (1/2) M - DL + gamma * SL @ M^{-1} @ K
            A_sys @ phi = M @ phi_inc
            gauge: int(phi) dS = 0  via Lagrange multiplier
        Same convergence properties; the iterative path just doesn't
        store ``A_sys`` as dense numpy.
        """
        from scipy.sparse.linalg import LinearOperator, gmres
        from ngsolve import (LinearForm, GridFunction, BND, ds,
                              Integrate, CF, InnerProduct, grad)

        t0 = time.perf_counter()
        ndof = self.ndof

        # RHS: <phi_inc, v>_S
        if isinstance(phi_inc_cf, np.ndarray):
            rhs_vec = self.M @ phi_inc_cf
        else:
            v_h1 = self.fes.TestFunction()
            lf = LinearForm(self.fes)
            lf += phi_inc_cf * v_h1.Trace() * ds
            lf.Assemble()
            rhs_vec = lf.vec.FV().NumPy().copy()

        # Per-node Z_s NOT supported in iterative path yet (the dense
        # solve has its own row-scaling logic).
        if isinstance(Z_s, np.ndarray):
            raise NotImplementedError(
                "solve_iterative does not support per-node Z_s yet -- "
                "use solve() with assemble_dense=True for that case.")

        gamma = (Z_s / (1j * omega * MU_0)
                 if (omega > 0 and Z_s != 0) else complex(0))

        # Helper: complex matvec on a real bem bilinear form.  NGSolve
        # bem operators are real-valued, so we dispatch real and imag
        # parts separately and recombine.
        fes = self.fes
        DL_bf = self._DL_bf
        SL_bf = self._SL_bf

        def _bem_matvec(bf, x_complex):
            gf_re = GridFunction(fes)
            gf_im = GridFunction(fes)
            gf_re.vec.FV().NumPy()[:] = x_complex.real
            gf_im.vec.FV().NumPy()[:] = x_complex.imag
            out_re = gf_re.vec.CreateVector()
            out_im = gf_im.vec.CreateVector()
            out_re.data = bf.mat * gf_re.vec
            out_im.data = bf.mat * gf_im.vec
            return (out_re.FV().NumPy().astype(complex)
                    + 1j * out_im.FV().NumPy().astype(complex))

        # A_sys @ x  (matrix-free)
        def A_matvec(x):
            x = np.asarray(x, dtype=complex).ravel()
            t1 = 0.5 * (self.M @ x)
            t2 = _bem_matvec(DL_bf, x)
            kx = self.K @ x
            mkx = self.M_inv @ kx
            slmkx = _bem_matvec(SL_bf, mkx)
            return t1 - t2 + gamma * slmkx

        # Augmented system to enforce int(phi) dS = 0:
        #   [ A    c ] [phi]   [b]
        #   [ c^T  0 ] [lam] = [0]
        c_gauge = self._c_gauge.astype(complex)

        def A_aug_matvec(x_aug):
            x = x_aug[:ndof]
            lam = x_aug[ndof]
            top = A_matvec(x) + lam * c_gauge
            bot = c_gauge @ x
            out = np.empty(ndof + 1, dtype=complex)
            out[:ndof] = top
            out[ndof] = bot
            return out

        A_op = LinearOperator((ndof + 1, ndof + 1),
                               matvec=A_aug_matvec, dtype=complex)

        rhs_aug = np.zeros(ndof + 1, dtype=complex)
        rhs_aug[:ndof] = rhs_vec

        sol_aug, info = gmres(A_op, rhs_aug, rtol=tol,
                               maxiter=maxiter, restart=restart)
        phi_vec = sol_aug[:ndof]

        t_solve = time.perf_counter() - t0

        if info > 0:
            print(f"[bem_sibc] WARN: GMRES did not converge in {info} iters")

        # Extract H_t_rms / P_density (same as solve())
        gf = GridFunction(self.fes)
        gf.vec.FV().NumPy()[:] = phi_vec.real
        Hsq_re = Integrate(InnerProduct(grad(gf), grad(gf)),
                            self.mesh, BND)
        gf.vec.FV().NumPy()[:] = phi_vec.imag
        Hsq_im = Integrate(InnerProduct(grad(gf), grad(gf)),
                            self.mesh, BND)
        area = Integrate(CF(1), self.mesh, BND)
        H_t_rms = math.sqrt((abs(Hsq_re) + abs(Hsq_im)) / abs(area))
        P_density = (0.5 * Z_s.real * H_t_rms ** 2
                     if Z_s != 0 else 0.0)

        gf_phi = GridFunction(self.fes)
        gf_phi.vec.FV().NumPy()[:] = phi_vec.real
        return {
            'phi': gf_phi,
            'phi_vec': phi_vec,
            'H_t_rms': float(H_t_rms),
            'P_density': float(P_density),
            'area': float(abs(area)),
            'gamma': complex(gamma),
            't_solve': round(t_solve, 3),
            'gmres_info': int(info),
        }

    # ------------------------------------------------------------------
    # 2-port driven solve (DELETED 2026-05-02): see solve_port_jump
    # ------------------------------------------------------------------
    # An earlier attempt (2026-05-02) implemented `solve_port` /
    # `solve_port_hacapk` with naive Dirichlet phi=V_source on src cap and
    # phi=0 on sink cap.  The implementation was numerically self-
    # consistent (dense LU vs HACApK GMRES agreed to 1e-12) but PHYSICALLY
    # BROKEN: on a closed conductor surface the exterior is simply
    # connected (for ball-topology bodies), so Ampere's law forces
    # enclosed current = 0 around any exterior loop.  The Dirichlet
    # voltage doesn't correspond to an EE port voltage in this framework.
    # The smoke test on a Cu rod gave KCL violation 89% and negative R.
    #
    # Replaced by `solve_port_jump`: closed-loop coil with cut-surface
    # jump constraint psi_src - psi_snk = I (current source driven), L/R
    # via Telegen energy/dissipation post-processing.  See Phase A series
    # in the implementation log for details.
    # ------------------------------------------------------------------

    def _identify_port_dofs(self, source_label, sink_label):
        """Identify body / source / sink surface DOF indices from BND labels.

        For the cut-surface (jump) port formulation, we treat the source
        and sink BND labels as the two "lips" of the topological cut.
        Returns dict with int64 arrays:
          'source_idx', 'sink_idx', 'body_idx', 'surface_idx'.
        Body = surface - source - sink.  Source and sink must be
        DISJOINT at the node level (caps must be geometrically separate).
        """
        from ngsolve import BND
        bnd_labels = list(self.mesh.GetBoundaries())

        src_label_indices = {i for i, n in enumerate(bnd_labels)
                             if n == source_label}
        snk_label_indices = {i for i, n in enumerate(bnd_labels)
                             if n == sink_label}
        if not src_label_indices:
            raise ValueError(
                f"port: source BND label {source_label!r} not found; "
                f"available: {bnd_labels}")
        if not snk_label_indices:
            raise ValueError(
                f"port: sink BND label {sink_label!r} not found; "
                f"available: {bnd_labels}")
        if src_label_indices & snk_label_indices:
            raise ValueError(
                f"port: source and sink labels resolve to the same BND "
                f"index — labels must be distinct.")

        src_set, snk_set, surf_set = set(), set(), set()
        for el in self.mesh.Elements(BND):
            for v in el.vertices:
                surf_set.add(v.nr)
                if el.index in src_label_indices:
                    src_set.add(v.nr)
                elif el.index in snk_label_indices:
                    snk_set.add(v.nr)

        overlap = src_set & snk_set
        if overlap:
            raise ValueError(
                f"port: {len(overlap)} vertices belong to BOTH source "
                f"and sink labels — port caps must be geometrically "
                f"disconnected.")

        body_set = surf_set - src_set - snk_set

        return {
            'source_idx': np.array(sorted(src_set), dtype=np.int64),
            'sink_idx': np.array(sorted(snk_set), dtype=np.int64),
            'body_idx': np.array(sorted(body_set), dtype=np.int64),
            'surface_idx': np.array(sorted(surf_set), dtype=np.int64),
        }


def telegen_extract_coil_LR(solver, phi_vec, current, omega, Z_s):
    """Telegen extraction of self-inductance L and resistance R for a
    closed-loop coil driven by a cut-surface BEM-SIBC solve.

    Inputs:
        solver  : ScalarBIESIBCSolver, with M and K already assembled.
        phi_vec : complex ndarray (ndof,) -- TOTAL psi on the coil
                  surface (= phi_inc + phi_response from solver.solve()).
        current : complex scalar -- prescribed coil current I [A]
                  (the current carried by the topological cut filament
                  used to build phi_inc).
        omega   : angular frequency [rad/s].
        Z_s     : surface impedance (complex scalar) -- the same value
                  used for the BEM-SIBC solve.

    Method:
        Time-averaged magnetic energy stored in the exterior::

            <W_mag> = (mu_0 / 4) * integral_(ext) |H|^2 dV
                    = -(mu_0 / 4) * Re(psi^H M q)        (Green's id.)

        Power dissipated in the body skin layer (Leontovich SIBC)::

            P_diss = (1/2) Re(Z_s) * (psi^H K psi)

        Self-inductance::

            L = 4 <W_mag> / |I|^2 = -mu_0 * Re(psi^H M q) / |I|^2

        Series resistance::

            R = P_diss / |I_rms|^2 = P_diss / (|I|^2 / 2) = 2 P_diss / |I|^2

        (using the convention that ``current`` is the peak phasor
        amplitude.)  Z = R + j omega L.

    The companion ``q`` vector is reconstructed via the workpiece-flow
    Robin closure ``q = gamma * M^-1 K psi`` where
    ``gamma = Z_s / (j omega mu_0)``.

    Returns:
        dict with R_port, L_port, Z_port, W_mag, P_diss, plus the raw
        bilinear-form diagnostics.
    """
    if abs(current) < 1e-30:
        raise ValueError(
            "telegen_extract_coil_LR: prescribed current is zero")
    I2 = abs(current) ** 2

    gamma = (Z_s / (1j * omega * MU_0)
             if (omega > 0 and Z_s != 0) else 0+0j)
    q = gamma * (solver.M_inv @ (solver.K @ phi_vec))

    Mq = solver.M @ q
    psi_H_M_q = np.vdot(phi_vec, Mq)   # complex
    W_mag = -0.25 * MU_0 * psi_H_M_q.real
    L = -MU_0 * psi_H_M_q.real / I2

    Kpsi = solver.K @ phi_vec
    psi_H_K_psi = np.vdot(phi_vec, Kpsi).real
    P_diss = 0.5 * Z_s.real * psi_H_K_psi
    R = 2.0 * P_diss / I2

    Z = R + 1j * omega * L
    return {
        'R_port': float(R),
        'L_port': float(L),
        'Z_port': complex(Z),
        'W_mag': float(W_mag),
        'P_diss': float(P_diss),
        'psi_H_M_q': complex(psi_H_M_q),
        'psi_H_K_psi': float(psi_H_K_psi),
        'gamma': complex(gamma),
    }



def surface_euler_characteristic(mesh):
    """Euler characteristic chi = V - E + F of a closed triangulated surface.

    ``chi == 2`` -> genus 0 (sphere-like).  ``chi == 0`` -> genus 1 (torus).
    ``genus = (2 - chi) / 2`` for a closed orientable surface.

    WHY THIS MATTERS for the scalar BIE + SIBC: the solver's surface
    current is ``J_s = n x (-grad_s phi)`` with a SINGLE-VALUED phi, so
    the net current through ANY cut of the surface is identically zero.
    On a genus-0 workpiece that is exact.  On a genus >= 1 workpiece whose
    handle links the coil flux (e.g. the Takahashi tube: the coil flux
    threads the bore), the physical eddy current contains a NET
    circulating (shorted-transformer-turn) component that this
    representation CANNOT express -- its Lenz screening is lost and the
    solver over-estimates H_t / P_wp.  Measured on Takahashi 7 kHz
    (2026-07-16, workpiece chi = 0): H_t 66.4 kA/m vs FEM 46.1 (x1.44),
    P_wp 38 kW vs 17 kW (x2.2), while the SAME solver matches the
    analytic mu_r-swept sphere (genus 0) benchmark to 0.3%.  Fixing the
    incident potential's branch-cut wall (surface-Poisson reconstruction,
    3% field consistency) moved P_wp by only ~5% -- the missing loop
    degree of freedom is the dominant error.

    Callers should compute this on the extracted workpiece surface and
    fail loud / caveat the output for ``chi != 2``.  The proper fix is a
    cohomology extension (add harmonic-1-form loop DOFs constrained by the
    Faraday EMF condition) -- see the radia.cohomology engine.

    Args:
        mesh: NGSolve surface mesh (or mesh whose BND elements form the
            closed surface).

    Returns:
        int Euler characteristic.
    """
    from ngsolve import BND

    verts = set()
    edges = set()
    n_faces = 0
    for el in mesh.Elements(BND):
        vid = [v.nr for v in el.vertices]
        n_faces += 1
        verts.update(vid)
        a, b, c = vid
        for u, w in ((a, b), (b, c), (c, a)):
            edges.add((u, w) if u < w else (w, u))
    return len(verts) - len(edges) + n_faces


def compute_phi_inc_from_loop(obs_points, loop_center, loop_radius, current,
                              n_quad=30, gap_deg=0):
    """Compute incident scalar magnetic potential from a circular current loop.

    Uses path integration: phi(P) = phi_axis(z) - int_axis^P H.dl
    H field from radia.biot_savart.h_segments (analytical formula).
    """
    from radia.biot_savart import h_segments_batch

    obs = np.asarray(obs_points, dtype=float)
    center = np.asarray(loop_center, dtype=float)
    a = float(loop_radius)
    I = float(current)
    if obs.ndim == 1:
        obs = obs.reshape(1, 3)

    obs_local = obs - center[np.newaxis, :]

    # Build coil wire segments as array (N_seg, 2, 3)
    arc_deg = 360 - gap_deg
    n_seg = max(200, int(arc_deg))
    theta = np.linspace(0, np.radians(arc_deg), n_seg + 1)
    coil_segs = np.zeros((n_seg, 2, 3))
    coil_segs[:, 0, 0] = a * np.cos(theta[:-1])
    coil_segs[:, 0, 1] = a * np.sin(theta[:-1])
    coil_segs[:, 1, 0] = a * np.cos(theta[1:])
    coil_segs[:, 1, 1] = a * np.sin(theta[1:])

    frac = arc_deg / 360.0

    # Gauss-Legendre quadrature for horizontal integration
    t_gl, w_gl = np.polynomial.legendre.leggauss(n_quad)
    t_01 = 0.5 * (t_gl + 1)
    w_01 = 0.5 * w_gl

    n_pts = len(obs_local)
    phi = np.zeros(n_pts)

    xy = obs_local[:, :2]
    z_arr = obs_local[:, 2]
    rho = np.sqrt(xy[:, 0]**2 + xy[:, 1]**2)

    # Analytical phi on z-axis for circular loop
    r_za = np.sqrt(z_arr**2 + a * a)
    phi_axis = (I / 2.0) * (1.0 - z_arr / r_za) * frac

    on_axis = rho < 1e-12 * a
    phi[on_axis] = phi_axis[on_axis]

    off_axis = ~on_axis
    if np.any(off_axis):
        idx_off = np.where(off_axis)[0]
        n_off = len(idx_off)

        # Build all quadrature points at once: (n_off * n_quad, 3)
        dl_vecs = np.zeros((n_off, 3))
        dl_vecs[:, 0] = obs_local[idx_off, 0]
        dl_vecs[:, 1] = obs_local[idx_off, 1]

        all_quad = np.zeros((n_off * n_quad, 3))
        for iq in range(n_quad):
            sl = slice(iq * n_off, (iq + 1) * n_off)
            all_quad[sl, 0] = t_01[iq] * dl_vecs[:, 0]
            all_quad[sl, 1] = t_01[iq] * dl_vecs[:, 1]
            all_quad[sl, 2] = z_arr[idx_off]

        # Vectorized H-field at all quadrature points
        H_all = h_segments_batch(coil_segs, all_quad, current=I)

        # Dot with dl_vec and integrate
        for iq in range(n_quad):
            sl = slice(iq * n_off, (iq + 1) * n_off)
            H_dot_dl = (H_all[sl, 0] * dl_vecs[:, 0]
                        + H_all[sl, 1] * dl_vecs[:, 1])
            phi[idx_off] -= w_01[iq] * H_dot_dl

        phi[idx_off] += phi_axis[idx_off]

    return phi



def compute_phi_inc_from_surface_J(obs_points, src_centroids, src_areas,
                                    src_J_vecs, n_quad=20, chunk_size=512):
    """Compute phi_inc from solved surface current via path integration.

    H field computed as sum of dipole contributions from surface elements.
    Fully vectorised over both the axis sweep (Stage 1) and the
    horizontal path sweep (Stage 2). For ih_sample.vol (170 obs nodes,
    4006 source panels, n_quad=20) this brings the runtime from ~160 s
    (per-node loop) down to ~1 s.

    Args:
        obs_points: (N, 3) observation node coordinates.
        src_centroids: (M, 3) source panel centroids.
        src_areas: (M,) source panel areas.
        src_J_vecs: (M, 3) source panel surface current densities (A/m).
        n_quad: Gauss-Legendre order for the path integrals.
        chunk_size: maximum number of obs nodes processed per chunk
            (memory cap; the inner H field array is
            chunk_size * n_quad * M * 3 floats).

    Returns:
        (N,) phi_inc values.
    """
    obs = np.asarray(obs_points, dtype=float)
    if obs.ndim == 1:
        obs = obs.reshape(1, 3)

    centers = np.asarray(src_centroids, dtype=float)
    areas = np.asarray(src_areas, dtype=float)
    if np.iscomplexobj(src_J_vecs):
        # np.asarray(..., dtype=float) on a complex array silently
        # discards the imaginary part (only a ComplexWarning).  The
        # BEM-A impedance-EFIE coil current IS complex -- callers must
        # bridge Re and Im in two separate calls and combine
        # (phi = phi_re + 1j*phi_im), as calc_inductance does.  Fail
        # fast per CLAUDE.md "No Fallbacks" instead of silently
        # truncating the phasor.
        raise TypeError(
            "compute_phi_inc_from_surface_J is real-only but received a "
            "complex src_J_vecs.  Call it separately on np.real(J) and "
            "np.imag(J) and combine phi_re + 1j*phi_im.")
    J = np.asarray(src_J_vecs, dtype=float)

    INV_4PI = 1.0 / (4.0 * np.pi)

    def _H_at_points(pts):
        """Vectorized H field at *pts* (shape (N, 3)) from all sources.

        Returns (N, 3) H values.  This is the per-chunk hot loop:
        memory ~ N * M * 3 floats for the dx and cross arrays.
        """
        dx = pts[:, None, :] - centers[None, :, :]  # (N, M, 3)
        r2 = np.sum(dx * dx, axis=2)
        r2 = np.maximum(r2, 1e-60)
        r3_inv = areas[None, :] / (r2 * np.sqrt(r2))  # (N, M)
        cross = np.cross(J[None, :, :], dx)  # (N, M, 3)
        return INV_4PI * np.sum(cross * r3_inv[:, :, None], axis=1)

    t_gl, w_gl = np.polynomial.legendre.leggauss(n_quad)
    t_01 = 0.5 * (t_gl + 1)
    w_01 = 0.5 * w_gl

    src_extent = float(np.max(np.abs(centers)))
    z_far = 20.0 * src_extent

    # ------------------------------------------------------------------
    # Stage 1: phi(0, 0, z) for every unique z used by an observation.
    # ------------------------------------------------------------------
    z_unique = np.unique(obs[:, 2])
    n_z = len(z_unique)

    # Quadrature points along the (0,0,z) -> (0,0,z_far) ray:
    #   z(t) = z_far - t * (z_far - z_unique[k])
    # x_quad shape: (n_z, n_quad, 3); we then flatten to feed _H_at_points.
    x_axis_quad = np.zeros((n_z, n_quad, 3))
    x_axis_quad[:, :, 2] = (z_far
                            - t_01[None, :] * (z_far - z_unique[:, None]))

    # Memory-bounded chunking on n_z * n_quad observation rows.
    rows_per_chunk = max(1, chunk_size // max(1, n_quad))
    H_axis = np.empty((n_z, n_quad, 3))
    flat = x_axis_quad.reshape(-1, 3)
    H_flat = np.empty_like(flat)
    n_rows = flat.shape[0]
    for s in range(0, n_rows, chunk_size):
        e = min(s + chunk_size, n_rows)
        H_flat[s:e] = _H_at_points(flat[s:e])
    H_axis = H_flat.reshape(n_z, n_quad, 3)

    # phi_axis[k] = (z_far - z_unique[k]) * int_0^1 H_z(z(t)) dt
    phi_axis_arr = (z_far - z_unique) * np.sum(
        w_01[None, :] * H_axis[:, :, 2], axis=1)
    phi_axis_lookup = dict(zip(z_unique, phi_axis_arr))

    # ------------------------------------------------------------------
    # Stage 2: horizontal path (0,0,z) -> (x,y,z) per observation node,
    # vectorised across all off-axis nodes.
    # ------------------------------------------------------------------
    n_obs = len(obs)
    phi = np.empty(n_obs)

    rho2 = obs[:, 0] ** 2 + obs[:, 1] ** 2
    eps_axis = (1e-12 * src_extent) ** 2
    on_axis = rho2 < eps_axis

    # On-axis nodes: phi(P) = phi_axis(z_P) directly.
    if np.any(on_axis):
        for k in np.where(on_axis)[0]:
            phi[k] = phi_axis_lookup[obs[k, 2]]

    off_idx = np.where(~on_axis)[0]
    if off_idx.size:
        x_off = obs[off_idx, 0]
        y_off = obs[off_idx, 1]
        z_off = obs[off_idx, 2]

        # dl vector per off-axis node (z component is zero).
        dl_x = x_off  # (n_off,)
        dl_y = y_off

        # Phi at the axis end of each path (z = z_off[k]).
        phi_axis_off = np.array([phi_axis_lookup[z] for z in z_off])

        # Process off-axis nodes in chunks to keep peak memory bounded.
        max_nodes_per_chunk = max(1, chunk_size // max(1, n_quad))
        for s in range(0, off_idx.size, max_nodes_per_chunk):
            e = min(s + max_nodes_per_chunk, off_idx.size)
            n_chunk = e - s

            # Build (n_chunk, n_quad, 3) quadrature points along
            # (0,0,z) -> (x,y,z): position(t) = (t*x, t*y, z)
            quad_pts = np.empty((n_chunk, n_quad, 3))
            quad_pts[:, :, 0] = t_01[None, :] * dl_x[s:e, None]
            quad_pts[:, :, 1] = t_01[None, :] * dl_y[s:e, None]
            quad_pts[:, :, 2] = z_off[s:e, None]
            flat_pts = quad_pts.reshape(-1, 3)

            H_flat = _H_at_points(flat_pts)
            H_quad = H_flat.reshape(n_chunk, n_quad, 3)

            # integrand[k, q] = H_x[k, q] * dl_x[k] + H_y[k, q] * dl_y[k]
            integrand = (H_quad[:, :, 0] * dl_x[s:e, None]
                         + H_quad[:, :, 1] * dl_y[s:e, None])
            path_int = np.sum(w_01[None, :] * integrand, axis=1)  # (n_chunk,)
            phi[off_idx[s:e]] = phi_axis_off[s:e] - path_int

    return phi


def _h_segments_complex(segments, obs_points, currents,
                         seg_chunk=4096):
    """Finite-segment Biot-Savart H at obs_points with complex per-segment currents.

    segments:  (N_seg, 2, 3) real endpoints
    obs_points: (N_obs, 3) real
    currents:  (N_seg,) complex
    seg_chunk: max segments processed at once (memory cap, IGNORED in
        C++ path -- the C++ kernel uses per-obs accumulation that is
        constant memory in N_seg)
    returns:   (N_obs, 3) complex

    Calls the C++ kernel `radia._radia_pybind._HFromSegmentsComplex`
    when available (~50-100x faster than the legacy NumPy fallback at
    typical IH coil sizes).  Falls back to the NumPy path when the C++
    symbol is missing (older wheel build).
    """
    obs = np.asarray(obs_points, dtype=float)
    if obs.ndim == 1:
        obs = obs.reshape(1, 3)

    # C++ fast path
    try:
        from radia import _radia_pybind as _rpb_bs
        _cpp_bs = _rpb_bs._HFromSegmentsComplex
    except (ImportError, AttributeError):
        _cpp_bs = None
    if _cpp_bs is not None:
        seg = np.ascontiguousarray(np.asarray(segments, dtype=float))
        I = np.asarray(currents, dtype=complex)
        I_re = np.ascontiguousarray(I.real)
        I_im = np.ascontiguousarray(I.imag)
        obs_c = np.ascontiguousarray(obs)
        H_re, H_im = _cpp_bs(seg, obs_c, I_re, I_im, 0)
        return H_re + 1j * H_im

    # ---- legacy NumPy fallback ----
    INV_4PI = 1.0 / (4.0 * np.pi)
    obs = np.asarray(obs_points, dtype=float)
    if obs.ndim == 1:
        obs = obs.reshape(1, 3)
    n_obs = len(obs)

    seg = np.asarray(segments, dtype=float)
    p1s_all = seg[:, 0, :]
    p2s_all = seg[:, 1, :]
    I_all = np.asarray(currents, dtype=complex)

    H_total = np.zeros((n_obs, 3), dtype=complex)

    # Chunk on segments to cap peak memory at
    # n_obs * seg_chunk * 3 (* 16 B complex).
    for s0 in range(0, len(p1s_all), seg_chunk):
        s1 = min(s0 + seg_chunk, len(p1s_all))
        p1s = p1s_all[s0:s1]  # (M, 3)
        p2s = p2s_all[s0:s1]
        I = I_all[s0:s1]      # (M,)
        M = s1 - s0

        dl = p2s - p1s                             # (M, 3)
        L = np.linalg.norm(dl, axis=1)             # (M,)
        valid_seg = L > 1e-30
        if not np.any(valid_seg):
            continue
        e_l = np.zeros_like(dl)
        e_l[valid_seg] = dl[valid_seg] / L[valid_seg, np.newaxis]

        # r1, r2 shape: (n_obs, M, 3)
        r1 = obs[:, None, :] - p1s[None, :, :]
        r2 = obs[:, None, :] - p2s[None, :, :]

        # cross(e_l, r1) per (obs, seg) -> (n_obs, M, 3)
        cross = np.cross(e_l[None, :, :], r1)
        d = np.linalg.norm(cross, axis=2)          # (n_obs, M)
        r1_mag = np.linalg.norm(r1, axis=2)
        r2_mag = np.linalg.norm(r2, axis=2)

        ok = (d > 1e-30) & (r1_mag > 1e-30) & (r2_mag > 1e-30)
        ok &= valid_seg[None, :]

        # Safe denominators to avoid /0 warnings outside 'ok'
        d_safe = np.where(ok, d, 1.0)
        r1_safe = np.where(ok, r1_mag, 1.0)
        r2_safe = np.where(ok, r2_mag, 1.0)

        cos_a1 = np.sum(r1 * e_l[None, :, :], axis=2) / r1_safe
        cos_a2 = np.sum(r2 * e_l[None, :, :], axis=2) / r2_safe
        geom = np.where(ok, INV_4PI / d_safe * (cos_a1 - cos_a2), 0.0)

        # e_perp = cross / d (zero where not ok)
        e_perp = cross / d_safe[..., None]

        # H contribution per (obs, seg): I[seg] * geom[obs, seg] * e_perp[obs, seg, :]
        weight = I[None, :] * geom                  # (n_obs, M) complex
        H_total += np.sum(weight[..., None] * e_perp, axis=1)

    return H_total


def compute_phi_inc_from_filaments(obs_points, filament_paths, currents,
                                    n_quad=20, chunk_size=512, z_far=None):
    """Complex phi_inc from a bundle of filaments with per-filament currents.

    Path-integration follows the same two-stage pattern as
    ``compute_phi_inc_from_surface_J`` (axis ray from z_far + horizontal
    sweep) but evaluates H via the exact finite-segment Biot-Savart
    formula and supports COMPLEX per-filament currents, so the output
    is a complex phi_inc suitable for ScalarBIESIBCSolver driven by a
    coil with non-trivial AC current distribution (e.g. Tier A skin
    effect).

    The method assumes the z-axis ray (0, 0, z_far) -> (0, 0, z_obs)
    and the horizontal ray (0, 0, z_obs) -> (x_obs, y_obs, z_obs) do
    not pierce any filament wire — true for a ring coil centered on
    the z-axis with the workpiece inside the coil bore.

    Args:
        obs_points: (N, 3) observation points [m].
        filament_paths: list of K filaments; each filament is a list of
            ``(p1, p2)`` endpoint tuples, as returned by
            ``CoilBuilder.to_filaments``.
        currents: length-K array of per-filament currents, real or
            complex [A].
        n_quad: Gauss-Legendre order for the 1-D path integrals.
        chunk_size: memory cap on ``n_obs * n_quad`` rows per H batch.
        z_far: reference height where phi ~ 0. Default: 20x max extent.

    Returns:
        (N,) complex phi values.
    """
    obs = np.asarray(obs_points, dtype=float)
    if obs.ndim == 1:
        obs = obs.reshape(1, 3)
    n_obs = len(obs)

    flat_segs = []
    flat_I = []
    for fil_segs, Ik in zip(filament_paths, currents):
        Ik_c = complex(Ik)
        for (p1, p2) in fil_segs:
            flat_segs.append((p1, p2))
            flat_I.append(Ik_c)

    if not flat_segs:
        return np.zeros(n_obs, dtype=complex)

    segs = np.asarray(flat_segs, dtype=float)
    Iseg = np.asarray(flat_I, dtype=complex)

    extent = float(np.max(np.abs(segs)))
    if z_far is None:
        z_far = 20.0 * max(extent, 1e-3)

    # --- Dipole tail correction: phi(0,0,z_far) ~ m_z / (4 pi z_far^2) ---
    # For the gauge phi(infinity) = 0 to be numerically honored, we seed the
    # axis integration with the dipole value at z_far instead of zero. This
    # cancels the O(1/z_far^2) offset from truncating the axis integral at
    # a finite z_far; residual is quadrupole-order (1/z_far^4).
    r_mid = 0.5 * (segs[:, 0, :] + segs[:, 1, :])
    dl_vec = segs[:, 1, :] - segs[:, 0, :]
    m_z = 0.5 * np.sum(
        Iseg * (r_mid[:, 0] * dl_vec[:, 1] - r_mid[:, 1] * dl_vec[:, 0]))
    phi_tail = m_z / (4.0 * np.pi * z_far ** 2)

    t_gl, w_gl = np.polynomial.legendre.leggauss(n_quad)
    t_01 = 0.5 * (t_gl + 1)
    w_01 = 0.5 * w_gl

    # ---- Stage 1: axis phi at unique z ----
    z_unique = np.unique(obs[:, 2])
    n_z = len(z_unique)
    axis_quad = np.zeros((n_z, n_quad, 3))
    axis_quad[:, :, 2] = z_far - t_01[None, :] * (z_far - z_unique[:, None])
    flat = axis_quad.reshape(-1, 3)

    max_rows = max(1, chunk_size)
    H_axis_flat = np.empty((flat.shape[0], 3), dtype=complex)
    for s in range(0, flat.shape[0], max_rows):
        e = min(s + max_rows, flat.shape[0])
        H_axis_flat[s:e] = _h_segments_complex(segs, flat[s:e], Iseg)
    H_axis = H_axis_flat.reshape(n_z, n_quad, 3)
    phi_axis_arr = phi_tail + (z_far - z_unique) * np.sum(
        w_01[None, :] * H_axis[:, :, 2], axis=1)
    phi_axis_lookup = dict(zip(z_unique, phi_axis_arr))

    # ---- Stage 2: horizontal path (0,0,z) -> (x,y,z) ----
    phi = np.empty(n_obs, dtype=complex)
    rho2 = obs[:, 0] ** 2 + obs[:, 1] ** 2
    eps_axis = (1e-12 * max(extent, 1e-3)) ** 2
    on_axis = rho2 < eps_axis
    if np.any(on_axis):
        for k in np.where(on_axis)[0]:
            phi[k] = phi_axis_lookup[obs[k, 2]]
    off_idx = np.where(~on_axis)[0]
    if off_idx.size:
        x_off = obs[off_idx, 0]
        y_off = obs[off_idx, 1]
        z_off = obs[off_idx, 2]
        phi_axis_off = np.array([phi_axis_lookup[z] for z in z_off],
                                dtype=complex)
        max_nodes_per_chunk = max(1, chunk_size // max(1, n_quad))
        for s in range(0, off_idx.size, max_nodes_per_chunk):
            e = min(s + max_nodes_per_chunk, off_idx.size)
            n_chunk = e - s
            quad_pts = np.empty((n_chunk, n_quad, 3))
            quad_pts[:, :, 0] = t_01[None, :] * x_off[s:e, None]
            quad_pts[:, :, 1] = t_01[None, :] * y_off[s:e, None]
            quad_pts[:, :, 2] = z_off[s:e, None]
            flat_pts = quad_pts.reshape(-1, 3)
            H_flat = _h_segments_complex(segs, flat_pts, Iseg)
            H_quad = H_flat.reshape(n_chunk, n_quad, 3)
            integrand = (H_quad[:, :, 0] * x_off[s:e, None]
                         + H_quad[:, :, 1] * y_off[s:e, None])
            path_int = np.sum(w_01[None, :] * integrand, axis=1)
            phi[off_idx[s:e]] = phi_axis_off[s:e] - path_int

    return phi


def compute_phi_inc_from_filaments_surface_path(
        wp_mesh, filament_paths, currents, *,
        ref_vertex_idx=0):
    """Reconstruct phi_inc on a closed workpiece surface mesh by
    integrating -H_BS·dl along surface edges, where H_BS is the
    Biot-Savart field of the filament bundle.

    Motivation (kubota report, 2026-05-21): the legacy two-stage
    ``compute_phi_inc_from_filaments`` assumes the axis ray
    (0,0,z_far)->(0,0,z_obs) and the horizontal ray
    (0,0,z_obs)->(x_obs,y_obs,z_obs) do not pierce any wire.  For a
    real coil with lead wires (x=129mm on 3turnCoil_work) plus 16
    perimeter filaments forming closed loops, this assumption fails:
    the path integrates through wire crossings and picks up
    spurious branch-cut jumps from the multivalued scalar
    potential.  Symptom: q_surf spatial peak appears on the
    coil-FAR face instead of coil-NEAR.

    This implementation avoids the problem entirely by:
      1. Evaluating H_BS exactly at each surface vertex via
         ``_h_segments_complex`` (the proven-correct primitive).
      2. Building surface-edge adjacency from wp_mesh's boundary
         elements (each triangle/quad contributes its edges).
      3. BFS from ``ref_vertex_idx`` along a spanning tree of the
         surface graph; at each new vertex v reached from parent u,
         set ``phi(v) = phi(u) - 0.5*(H(u)+H(v))·(x_v - x_u)``
         (trapezoidal rule along the edge).

    Why this works topologically: the workpiece outer surface is a
    closed simply-connected surface (sphere topology -- the
    workpiece is a solid cylinder / blob, not a torus).  Any loop
    on the surface bounds a disk on the surface; the disk does NOT
    enclose the exterior coil currents, so ∮H·dl = 0 around any
    such loop and phi is single-valued on the surface.  The BFS
    spanning tree yields one consistent assignment; discretisation
    error from trapezoidal integration is O(h^2) per edge.

    Args:
        wp_mesh: NGSolve Mesh of the workpiece outer surface (the
            extracted 2D mesh that ScalarBIESIBCSolver consumes).
        filament_paths: list of K filaments; each filament is a
            list of ``(p1, p2)`` endpoint tuples.
        currents: length-K array of per-filament currents, real or
            complex [A].
        ref_vertex_idx: vertex used as the gauge (phi=0).  Any
            vertex on the surface is valid; the global gauge is
            absorbed by the BIE's grad-of-phi operator anyway.

    Returns:
        (n_surface_vertices,) complex phi_inc.
    """
    from ngsolve import BND

    n_v = wp_mesh.nv

    # 1. Pre-flatten filament segments + currents (shared with the
    #    legacy two-stage routine).
    flat_segs = []
    flat_I = []
    for fil_segs, Ik in zip(filament_paths, currents):
        Ik_c = complex(Ik)
        for (p1, p2) in fil_segs:
            flat_segs.append((p1, p2))
            flat_I.append(Ik_c)
    if not flat_segs:
        return np.zeros(n_v, dtype=complex)
    segs = np.asarray(flat_segs, dtype=float)
    Iseg = np.asarray(flat_I, dtype=complex)

    # 2. H_BS at every surface vertex (vectorised; single call).
    coords = np.empty((n_v, 3), dtype=float)
    for vi in range(n_v):
        p = wp_mesh.vertices[vi].point
        coords[vi, 0] = p[0]
        coords[vi, 1] = p[1]
        coords[vi, 2] = p[2]
    H_vert = _h_segments_complex(segs, coords, Iseg)  # (n_v, 3) complex

    # 3. Surface adjacency from the mesh's BND elements.  Each
    #    triangle (3 verts) contributes 3 edges; each quad (4)
    #    contributes 4 edges.  Build an EDGE LIST (not just sets)
    #    so the LSQ in step 4 can use each edge as a separate
    #    constraint.  De-duplicate by canonical ordering (u < v).
    edge_set = set()
    for el in wp_mesh.Elements(BND):
        verts = [v.nr for v in el.vertices]
        k = len(verts)
        for i in range(k):
            a, b = verts[i], verts[(i + 1) % k]
            edge_set.add((min(a, b), max(a, b)))
    edges = np.array(sorted(edge_set), dtype=int)
    n_e = len(edges)

    # 4. Least-squares solve for single-valued phi consistent with
    #    the discretised -H·dl integrand along every edge.
    #    System: for each edge (u, v), phi(v) - phi(u) = b_e where
    #    b_e = -0.5*(H(u)+H(v))·(x_v - x_u).
    #    BFS spanning tree gave path-dependent accumulation of the
    #    trapezoidal-rule error and put the q_surf peak at the wp
    #    end caps instead of under the coil.  LSQ over ALL edges
    #    averages the noise out and matches the true |H|^2 peak.
    #
    #    Sparse incidence matrix A (n_e, n_v): A[e, u] = -1,
    #    A[e, v] = +1.  Gauge: drop column ref_vertex_idx (and
    #    enforce phi(ref) = 0).  scipy.sparse.linalg.lsmr solves
    #    A_reduced @ phi_reduced = b efficiently.
    from scipy import sparse
    from scipy.sparse.linalg import lsmr

    rows = np.repeat(np.arange(n_e), 2)
    cols = edges.reshape(-1)
    data = np.tile(np.array([-1.0, 1.0]), n_e)
    A = sparse.csr_matrix((data, (rows, cols)),
                          shape=(n_e, n_v))

    # Trapezoidal -H·dl for each edge (separately for re/im).
    u_idx = edges[:, 0]
    v_idx = edges[:, 1]
    Hmid = 0.5 * (H_vert[u_idx] + H_vert[v_idx])
    dl = coords[v_idx] - coords[u_idx]
    b_full = -(Hmid[:, 0] * dl[:, 0]
               + Hmid[:, 1] * dl[:, 1]
               + Hmid[:, 2] * dl[:, 2])  # complex (n_e,)

    # Drop the gauge column and solve for the remaining DOFs
    keep_cols = np.array([i for i in range(n_v) if i != ref_vertex_idx])
    A_red = A[:, keep_cols]

    # Solve real and imag parts separately (lsmr is real-only).
    phi = np.zeros(n_v, dtype=complex)
    phi_re = lsmr(A_red, b_full.real, atol=1e-10, btol=1e-10)[0]
    phi_im = lsmr(A_red, b_full.imag, atol=1e-10, btol=1e-10)[0]
    phi[keep_cols] = phi_re + 1j * phi_im
    # phi[ref_vertex_idx] stays 0 (gauge)

    return phi


def compute_phi_inc_from_filaments_arbitrary_axis(
        obs_points, filament_paths, currents, *,
        n_quad=20, chunk_size=512,
        far_ref=None, axis_dir=None):
    """Complex phi_inc with arbitrary integration axis (generalises
    compute_phi_inc_from_filaments to non-z-symmetric coils).

    The original ``compute_phi_inc_from_filaments`` assumes the coil is
    a ring centered on the +z axis with the workpiece inside the bore --
    its 2-stage path (axis ray (0,0,z_far)->(0,0,z_obs) + horizontal
    sweep) only stays clear of filament wires under that geometry.

    This variant lets the caller specify:
      * ``axis_dir``  : unit vector for the "axis" leg of the path
                         (default (0,0,1) = z-axis, matching the
                         original behaviour).
      * ``far_ref``    : point at which phi is taken to be ~0
                         (default: ``axis_dir * 20*max_extent``,
                         matching the original ``z_far``).

    Path topology (matches the original 2-stage form, just rotated):
      Stage 1: ``far_ref`` --(along axis_dir)--> ``(obs projected onto
                 the axis line through far_ref)``
      Stage 2: that intermediate --(perpendicular to axis_dir)-->
                 obs

    The dipole tail correction at far_ref is the projection of the
    coil's magnetic dipole moment onto axis_dir (so ``axis_dir = +z``
    recovers the m_z formula).

    For a 3-turn coil whose helix axis is z but is offset to +x (e.g.
    ``3turnCoil_work``), the original function still works because
    z-axis is INSIDE the coil bore.  This generalised form is needed
    when the coil's helix axis is not aligned with z, OR when the
    workpiece sticks out of the coil bore so far that the original
    horizontal sweep would graze a filament wire.

    Args:
        obs_points: (N, 3) observation points [m].
        filament_paths: list of K filaments; each filament is a list of
            ``(p1, p2)`` endpoint tuples.
        currents: length-K array of per-filament currents [A] (real or
            complex).
        n_quad: Gauss-Legendre order for the 1-D path integrals.
        chunk_size: memory cap on ``n_obs * n_quad`` rows per H batch.
        far_ref: (3,) point [m] where phi ~ 0.  Default: 20x max extent
            along axis_dir.
        axis_dir: (3,) unit vector for stage-1 integration direction.
            Default (0, 0, 1).

    Returns:
        (N,) complex phi values.
    """
    obs = np.asarray(obs_points, dtype=float)
    if obs.ndim == 1:
        obs = obs.reshape(1, 3)
    n_obs = len(obs)

    flat_segs = []
    flat_I = []
    for fil_segs, Ik in zip(filament_paths, currents):
        Ik_c = complex(Ik)
        for (p1, p2) in fil_segs:
            flat_segs.append((p1, p2))
            flat_I.append(Ik_c)
    if not flat_segs:
        return np.zeros(n_obs, dtype=complex)
    segs = np.asarray(flat_segs, dtype=float)
    Iseg = np.asarray(flat_I, dtype=complex)

    if axis_dir is None:
        axis_dir = np.array([0.0, 0.0, 1.0])
    else:
        axis_dir = np.asarray(axis_dir, dtype=float)
        n = np.linalg.norm(axis_dir)
        if n < 1e-30:
            raise ValueError("axis_dir must be non-zero")
        axis_dir = axis_dir / n

    extent = float(np.max(np.abs(segs)))
    if far_ref is None:
        far_ref = axis_dir * (20.0 * max(extent, 1e-3))
    far_ref = np.asarray(far_ref, dtype=float)

    # Dipole tail correction at far_ref:
    #   phi_dipole(r) = (m . r_hat) / (4 pi |r|^2)
    # m = (1/2) sum_seg I_seg * (r_mid x dl)
    r_mid = 0.5 * (segs[:, 0, :] + segs[:, 1, :])
    dl_vec = segs[:, 1, :] - segs[:, 0, :]
    cross_md = np.cross(r_mid, dl_vec)  # (N_seg, 3) real
    # m_complex per cartesian component: shape (3,) complex
    m_vec = 0.5 * np.sum(Iseg[:, None] * cross_md, axis=0)
    r_far_mag = np.linalg.norm(far_ref)
    if r_far_mag > 1e-30:
        r_hat = far_ref / r_far_mag
        phi_far = float(np.dot(m_vec.real, r_hat)) / (4.0 * np.pi * r_far_mag ** 2) \
                  + 1j * float(np.dot(m_vec.imag, r_hat)) / (4.0 * np.pi * r_far_mag ** 2)
    else:
        phi_far = 0.0 + 0.0j

    t_gl, w_gl = np.polynomial.legendre.leggauss(n_quad)
    t_01 = 0.5 * (t_gl + 1)  # (n_quad,)
    w_01 = 0.5 * w_gl

    # ---- Stage 1: along axis_dir, from far_ref to (obs projected onto
    # axis line through far_ref).
    # Param: r(t) = far_ref + t * (mid_obs - far_ref)
    # where mid_obs = far_ref + ((obs - far_ref) . axis_dir) * axis_dir
    proj = np.dot(obs - far_ref[None, :], axis_dir)  # (n_obs,) along-axis distance
    mid_obs = far_ref[None, :] + proj[:, None] * axis_dir[None, :]  # (n_obs, 3)
    # Segment vector from far_ref to mid_obs (length = proj along axis_dir)
    s1_vec = mid_obs - far_ref[None, :]  # (n_obs, 3)

    # Quadrature points for stage 1
    quad_pts_s1 = far_ref[None, None, :] + t_01[None, :, None] * s1_vec[:, None, :]
    flat_s1 = quad_pts_s1.reshape(-1, 3)
    H_s1_flat = np.empty((flat_s1.shape[0], 3), dtype=complex)
    max_rows = max(1, chunk_size)
    for s in range(0, flat_s1.shape[0], max_rows):
        e = min(s + max_rows, flat_s1.shape[0])
        H_s1_flat[s:e] = _h_segments_complex(segs, flat_s1[s:e], Iseg)
    H_s1 = H_s1_flat.reshape(n_obs, n_quad, 3)
    integrand_s1 = np.sum(H_s1 * s1_vec[:, None, :], axis=2)  # H . dr/dt
    # phi at mid_obs = phi_far - integral of H.dl from far_ref to mid_obs
    phi_mid = phi_far - np.sum(w_01[None, :] * integrand_s1, axis=1)

    # ---- Stage 2: from mid_obs to obs (perpendicular to axis_dir).
    s2_vec = obs - mid_obs  # (n_obs, 3)
    quad_pts_s2 = mid_obs[:, None, :] + t_01[None, :, None] * s2_vec[:, None, :]
    flat_s2 = quad_pts_s2.reshape(-1, 3)
    H_s2_flat = np.empty((flat_s2.shape[0], 3), dtype=complex)
    for s in range(0, flat_s2.shape[0], max_rows):
        e = min(s + max_rows, flat_s2.shape[0])
        H_s2_flat[s:e] = _h_segments_complex(segs, flat_s2[s:e], Iseg)
    H_s2 = H_s2_flat.reshape(n_obs, n_quad, 3)
    integrand_s2 = np.sum(H_s2 * s2_vec[:, None, :], axis=2)
    phi = phi_mid - np.sum(w_01[None, :] * integrand_s2, axis=1)

    return phi


def extract_H_t_per_dof_grad(phi_vec_complex, wp_mesh):
    """Per-vertex |H_t| via manual triangle-wise P1 gradient of phi.

    H_t = -grad_s(phi) on the workpiece surface.  For each surface
    triangle build the constant gradient of the linear interpolant
    of phi from its 3 vertex values, accumulate area-weighted to
    the 3 vertices, and return ``|H_t|_per_vertex``.

    This is the physically-correct local-gradient extraction.  It
    REPLACES the legacy Galerkin localization ``phi_i * (K @ phi)_i``,
    which is a Laplacian sample (samples the LOCAL Laplacian of phi,
    not the LOCAL gradient norm) and consequently mis-locates the
    spatial peak.  The v4.67.0 release fixed this for the q_surf
    spatial output; this helper applies the same fix to the ESIM
    per-DOF |H_t| extraction (Karl iteration outer loop).

    Parameters
    ----------
    phi_vec_complex : (n_v,) complex
        Magnetic scalar potential on workpiece-surface vertices
        (P1 DOFs in the BIE basis order=1 path).  For order>=2 the
        caller should down-project to vertices before invoking this
        function.
    wp_mesh : ngsolve.Mesh
        Workpiece SURFACE mesh (2-D embedded in 3-D).  We read
        vertices and BND triangle connectivity.

    Returns
    -------
    H_t_per_vertex : (n_v,) float
        |H_t|_i = sqrt(|grad_s phi|^2_real + |grad_s phi|^2_imag)
        at each surface vertex, area-weighted from the incident
        triangles.

    Notes
    -----
    The 3D triangle-P1 gradient formula:
        grad N_0 = (p_2 - p_1) x n_hat / (2 area)
        grad N_1 = (p_0 - p_2) x n_hat / (2 area)
        grad N_2 = (p_1 - p_0) x n_hat / (2 area)
    where ``n_hat`` is the outward unit normal.  These are in-plane
    perpendicular to the opposite edge and have magnitude 1/h_i
    (h_i = altitude from vertex i to the opposite edge).
    """
    import numpy as _np
    from ngsolve import BND as _BND

    n_v = wp_mesh.nv
    tri_v = []
    for el in wp_mesh.Elements(_BND):
        tri_v.append([v.nr for v in el.vertices])
    tri_v = _np.asarray(tri_v, dtype=_np.int64)
    vert_xyz = _np.asarray(
        [list(wp_mesh.vertices[i].point) for i in range(n_v)],
        dtype=float)

    p0 = vert_xyz[tri_v[:, 0]]
    p1 = vert_xyz[tri_v[:, 1]]
    p2 = vert_xyz[tri_v[:, 2]]
    normal = _np.cross(p1 - p0, p2 - p0)
    area2 = _np.linalg.norm(normal, axis=1)
    tri_area = 0.5 * area2
    n_hat = normal / _np.maximum(area2[:, None], 1e-30)

    gN0 = _np.cross(p2 - p1, n_hat) / _np.maximum(area2[:, None], 1e-30)
    gN1 = _np.cross(p0 - p2, n_hat) / _np.maximum(area2[:, None], 1e-30)
    gN2 = _np.cross(p1 - p0, n_hat) / _np.maximum(area2[:, None], 1e-30)

    phi = _np.asarray(phi_vec_complex)
    phi_re = phi.real
    phi_im = phi.imag

    grad_re = (phi_re[tri_v[:, 0:1]] * gN0
               + phi_re[tri_v[:, 1:2]] * gN1
               + phi_re[tri_v[:, 2:3]] * gN2)
    grad_im = (phi_im[tri_v[:, 0:1]] * gN0
               + phi_im[tri_v[:, 1:2]] * gN1
               + phi_im[tri_v[:, 2:3]] * gN2)
    Hsq_tri = (_np.sum(grad_re * grad_re, axis=1)
               + _np.sum(grad_im * grad_im, axis=1))

    vert_Hsq_num = _np.zeros(n_v)
    vert_area_sum = _np.zeros(n_v)
    for j in range(3):
        _np.add.at(vert_Hsq_num, tri_v[:, j], tri_area * Hsq_tri)
        _np.add.at(vert_area_sum, tri_v[:, j], tri_area)
    Hsq_per_vertex = vert_Hsq_num / _np.maximum(vert_area_sum, 1e-30)
    return _np.sqrt(_np.maximum(Hsq_per_vertex, 0.0))


def extract_surface_J_from_phi(mesh, phi_vec_complex, order=1):
    """Extract per-element complex surface current density J_s.

    SIBC scalar BIE: phi on surface is magnetic scalar potential. The
    tangential magnetic field is H_t = -grad_s(phi), and the surface
    current density is J_s = n x H_t = -n x grad_s(phi).

    Args:
        mesh: NGSolve surface mesh (2D or 3D with BND).
        phi_vec_complex: (ndof,) complex coefficient vector for phi.
        order: H1 polynomial order matching the BIE solve (default 1).

    Returns:
        centroids: (M, 3) float, element area-centroids.
        areas: (M,) float, element areas.
        J_s: (M, 3) complex, element-averaged surface current density.
    """
    from ngsolve import (H1, GridFunction, Integrate, CF, BND,
                         grad, specialcf, Cross)

    fes = H1(mesh, order=order)
    gf_re = GridFunction(fes)
    gf_im = GridFunction(fes)
    gf_re.vec.FV().NumPy()[:] = phi_vec_complex.real
    gf_im.vec.FV().NumPy()[:] = phi_vec_complex.imag

    n_cf = specialcf.normal(3)
    Jre_cf = -Cross(n_cf, grad(gf_re).Trace())
    Jim_cf = -Cross(n_cf, grad(gf_im).Trace())

    elem_A = Integrate(CF(1), mesh, VOL_or_BND=BND, element_wise=True)
    Jre = [Integrate(Jre_cf[i], mesh, VOL_or_BND=BND, element_wise=True)
           for i in range(3)]
    Jim = [Integrate(Jim_cf[i], mesh, VOL_or_BND=BND, element_wise=True)
           for i in range(3)]

    centroids, areas, J_s = [], [], []
    for el in mesh.Elements(BND):
        a = abs(elem_A[el.nr])
        if a < 1e-30:
            continue
        verts = [mesh.vertices[v.nr].point for v in el.vertices]
        c = np.mean([(v[0], v[1], v[2]) for v in verts], axis=0)
        Jvec = np.array([(Jre[i][el.nr] + 1j * Jim[i][el.nr]) / a
                         for i in range(3)], dtype=complex)
        centroids.append(c)
        areas.append(a)
        J_s.append(Jvec)
    return np.array(centroids), np.array(areas), np.array(J_s)


def A_from_surface_J(obs_points, src_centroids, src_areas, src_J,
                     chunk_size=4096):
    """Magnetic vector potential at obs points from surface current density.

        A(r) = (mu_0 / 4 pi) * sum_j  J_s[j] * area[j] / |r - c_j|

    Args:
        obs_points: (N, 3) float observation coordinates.
        src_centroids: (M, 3) panel centroids.
        src_areas: (M,) panel areas.
        src_J: (M, 3) complex (or real) surface current density.
        chunk_size: rows of the pair-wise distance matrix per chunk
            (memory cap; each chunk is N_chunk x M).

    Returns:
        (N, 3) complex array (real if src_J is real).
    """
    r = np.asarray(obs_points, dtype=float)
    c = np.asarray(src_centroids, dtype=float)
    a = np.asarray(src_areas, dtype=float)
    J = np.asarray(src_J)
    if r.ndim == 1:
        r = r.reshape(1, 3)
    N = r.shape[0]
    dtype = complex if np.iscomplexobj(J) else float
    MU_0 = 4e-7 * np.pi
    pref = MU_0 / (4.0 * np.pi)
    A = np.zeros((N, 3), dtype=dtype)
    for s in range(0, N, chunk_size):
        e = min(s + chunk_size, N)
        dx = r[s:e, None, 0] - c[None, :, 0]
        dy = r[s:e, None, 1] - c[None, :, 1]
        dz = r[s:e, None, 2] - c[None, :, 2]
        dist = np.sqrt(dx * dx + dy * dy + dz * dz)
        dist[dist < 1e-12] = 1e-12
        weight = pref * (a / dist)                    # (n, M)
        A[s:e] = np.einsum('nm,mc->nc', weight, J)
    return A


def flux_linkage_in_filaments(filament_paths, src_centroids, src_areas, src_J):
    """Phi_k = integral of A_src . dl along each filament polyline.

    Uses midpoint rule on each polyline segment (A evaluated at midpoint,
    dotted with the segment's dl vector).

    Args:
        filament_paths: list of K filament polylines, each a list of
            ``(p1, p2)`` endpoint tuples (3-vectors in meters).
        src_centroids, src_areas, src_J: panel data for A_from_surface_J.

    Returns:
        (K,) complex (or real) array of flux linkages.
    """
    K = len(filament_paths)
    dtype = complex if np.iscomplexobj(src_J) else float
    Phi = np.zeros(K, dtype=dtype)
    for k, path in enumerate(filament_paths):
        mids = np.array([0.5 * (np.asarray(p1, float) + np.asarray(p2, float))
                         for p1, p2 in path])
        dls = np.array([np.asarray(p2, float) - np.asarray(p1, float)
                        for p1, p2 in path])
        A_mid = A_from_surface_J(mids, src_centroids, src_areas, src_J)
        Phi[k] = np.sum(A_mid * dls)
    return Phi


# Legacy alias
def compute_phi_inc_from_coil(obs_points, coil_path, current, n_quad=50):
    """Compute phi_inc using filament approximation (legacy interface)."""
    path = np.asarray(coil_path, dtype=float)
    center = np.mean(path, axis=0)
    center[2] = path[0, 2]
    radius = np.mean(np.linalg.norm(path[:, :2] - center[np.newaxis, :2], axis=1))
    theta_start = np.arctan2(path[0, 1], path[0, 0])
    theta_end = np.arctan2(path[-1, 1], path[-1, 0])
    arc = (theta_end - theta_start) % (2 * np.pi)
    gap_deg = max(0, 360 - np.degrees(arc))
    return compute_phi_inc_from_loop(obs_points, center, radius, current,
                                     n_quad=n_quad, gap_deg=gap_deg)
