"""Reduced vector potential solver (A_r formulation, Radia + NGSolve).

Implements the reduced vector potential method for magnetostatic
field computation combining Radia BEM source fields with NGSolve FEM.

Method:
    B = B_s + curl(A_r)
    H = nu * B
    Weak form: int(nu * curl(A_r) . curl(v)) = int((nu_0 - nu) * B_s . curl(v))

Features:
    - Linear materials (mu_r per domain)
    - Nonlinear materials (B-H curve, Newton with coenergy)
    - Nonlinear Picard iteration with under-relaxation
    - HCurl FE space with nograds=True gauge
    - BDDC+CG iterative solver for large problems

Usage:
    import radia as rad
    from radia.vector_potential_solver import VectorPotentialSolver

    solver = VectorPotentialSolver(mesh, mu_r_dict={'iron': 1000}, order=2)
    solver.set_source_from_radia(coil, resolution=41)
    solver.solve_linear(dirichlet='outer')
    B_cf = solver.get_B()

Comparison with ScalarPotentialSolver (Simkin):
    ScalarPotentialSolver: H = H_s - grad(phi),  source = H_s
    VectorPotentialSolver: B = B_s + curl(A_r),   source = B_s
"""

import numpy as np
from math import sqrt as msqrt

MU_0 = 4 * np.pi * 1e-7


def _build_nu_of_b_interpolator(bh_data):
    """Invert the production monotone PCHIP B(H) law for reduced-A.

    HDiv-MMM and Omega-reduced-Omega consume B(H) directly.  Reduced-A needs
    H(B), so each scalar value is inverted against that same PCHIP curve rather
    than constructing a second interpolation law.  Above the table, the shared
    vacuum-slope continuation has the analytic inverse.
    """
    from scipy.optimize import brentq
    from radia.scalar_potential_solver import _build_bh_interpolator

    bh = np.asarray(bh_data, dtype=float)
    if bh.ndim != 2 or bh.shape[1] < 2 or bh.shape[0] < 2:
        raise ValueError("bh_data must contain at least two [H, B] rows")
    h_values = bh[:, 0]
    b_values = bh[:, 1]
    b_of_h = _build_bh_interpolator(bh)
    h_min = float(h_values[0])
    h_max = float(h_values[-1])
    b_min = float(b_values[0])
    b_max = float(b_values[-1])
    positive = np.flatnonzero((h_values > 0.0) & (b_values > 0.0))
    if positive.size == 0:
        raise ValueError("bh_data must contain a positive H, B sample")
    initial_reluctivity = float(h_values[positive[0]] / b_values[positive[0]])

    def nu_of_b(value):
        b_value = max(float(value), 0.0)
        if b_value <= max(b_min, 1.0e-15):
            return initial_reluctivity
        if b_value <= b_max:
            h_value = brentq(
                lambda candidate: b_of_h(candidate) - b_value,
                h_min,
                h_max,
                xtol=1.0e-10,
                rtol=1.0e-13,
                maxiter=100,
            )
        else:
            h_value = h_max + (b_value - b_max) / MU_0
        return float(h_value / b_value)

    return nu_of_b, b_max


class VectorPotentialSolver:
    """Reduced vector potential magnetostatic solver (Radia + NGSolve).

    Parameters
    ----------
    mesh : ngsolve.Mesh
        FEM mesh with material labels.
    iron_domains : str or list of str
        Material name(s) for iron regions.
    mu_r : float
        Default relative permeability for iron.
    mu_r_dict : dict, optional
        Per-domain mu_r, e.g. {'iron': 1000, 'steel': 500}.
    order : int
        FEM polynomial order (default: 2).
    kelvin_region : str, optional
        Material name of Kelvin exterior domain region (reserved for future use).
    kelvin_radius : float, optional
        Inner radius R of Kelvin exterior domain [m].
    kelvin_center : list/tuple, optional
        Center of Kelvin exterior domain [x, y, z] in meters.
    ams_options : dict, optional
        AMS preconditioner options. Keys:
        - 'chebyshev_degree' (int, default 3): Chebyshev smoother degree.
        - 'num_smooth' (int, default 2): Pre/post smoothing steps.
        - 'h1_solver' (str, default 'auto'): H1 subspace solver.
          'auto' selects direct (sparsecholesky) for ndof_h1 < 100k,
          h1amg otherwise. Can be 'direct' or 'h1amg'.
        - 'eigenratio' (float, default 30.0): Chebyshev smoothing
          interval [lambda_max/eigenratio, lambda_max].
    """

    def __init__(self, mesh, iron_domains='iron', mu_r=1000.0, order=2,
                 mu_r_dict=None, kelvin_region=None, kelvin_radius=None,
                 kelvin_center=None, ams_options=None):
        self.mesh = mesh
        self.order = order

        # Material setup
        if mu_r_dict is not None:
            self._mu_r_dict = dict(mu_r_dict)
            self.iron_domains = list(mu_r_dict.keys())
            self.mu_r = max(mu_r_dict.values())
        else:
            if isinstance(iron_domains, str):
                iron_domains = [iron_domains]
            self.iron_domains = iron_domains
            self.mu_r = float(mu_r)
            self._mu_r_dict = {d: float(mu_r) for d in self.iron_domains}

        # Kelvin transform
        self._kelvin_region = kelvin_region
        self._kelvin_radius = kelvin_radius
        self._kelvin_center = kelvin_center or [0, 0, 0]
        if self._kelvin_region and self._kelvin_radius is None:
            raise ValueError("kelvin_radius is required when kelvin_region is set")

        # AMS preconditioner options
        self._ams_options = {
            'chebyshev_degree': 3,
            'num_smooth': 2,
            'h1_solver': 'auto',
            'eigenratio': 30.0,
        }
        if ams_options:
            self._ams_options.update(ams_options)
        self._ams_cache = {}

        # Source field: B_s (vector, 3 components)
        self._B_source_cf = None
        self._radia_obj = None

        # Solution storage
        self._A_gf = None
        self._B_cf = None
        self._H_cf = None

    # ------------------------------------------------------------------
    # Source field setup
    # ------------------------------------------------------------------

    def set_source_from_radia(self, radia_obj, resolution=41, bbox=None):
        """Set source field B_s from a Radia coil object.

        Computes B_s on a voxel grid via ``rad.Fld(obj, 'b', points)``
        and creates a VoxelCoefficient for efficient FEM integration.

        Parameters
        ----------
        radia_obj : int
            Radia object handle (coil).
        resolution : int
            Voxel grid resolution per dimension.
        bbox : list of [min, max] pairs, optional
            Custom bounding box [[xmin,xmax],[ymin,ymax],[zmin,zmax]].
        """
        import radia as rad
        from ngsolve import VoxelCoefficient, CF

        self._radia_obj = radia_obj

        if bbox is not None:
            pass
        elif self._kelvin_region:
            bbox = self._compute_physical_bbox()
        else:
            pmin, pmax = self.mesh.ngmesh.bounding_box
            pmin = [pmin[i] for i in range(3)]
            pmax = [pmax[i] for i in range(3)]
            margin = 0.05 * max(pmax[i] - pmin[i] for i in range(3))
            bbox = [[pmin[i] - margin, pmax[i] + margin] for i in range(3)]

        nx = ny = nz = resolution
        x = np.linspace(bbox[0][0], bbox[0][1], nx)
        y = np.linspace(bbox[1][0], bbox[1][1], ny)
        z = np.linspace(bbox[2][0], bbox[2][1], nz)

        xx, yy, zz = np.meshgrid(x, y, z, indexing='ij')
        points = np.column_stack([xx.ravel(), yy.ravel(), zz.ravel()])

        B_field = np.asarray(rad.Fld(radia_obj, 'b', points))

        start = (bbox[0][0], bbox[1][0], bbox[2][0])
        end = (bbox[0][1], bbox[1][1], bbox[2][1])

        cfs = []
        for comp in range(3):
            data = B_field[:, comp].reshape(nx, ny, nz)
            data = np.ascontiguousarray(data.transpose(2, 1, 0))
            cfs.append(VoxelCoefficient(start, end, data, linear=True))

        self._B_source_cf = CF(tuple(cfs))

    def set_source_from_callback(self, b_func, resolution=41):
        """Set source field B_s from a Python callback.

        Parameters
        ----------
        b_func : callable
            ``b_func(x, y, z) -> (Bx, By, Bz)`` returning B in Tesla.
            **Contract**: ``b_func`` MUST accept NumPy ndarrays of equal
            shape for x, y, z and return an array-like of shape
            ``(..., 3)`` (the trailing axis is the vector component).
            A scalar-only callback (one that only handles plain floats)
            is detected at runtime and wrapped via ``np.vectorize`` with
            a ``DeprecationWarning`` -- this fallback is ~50-200x slower
            than a properly vectorized callback at resolution=41.
        resolution : int
            Voxel grid resolution per dimension. Default 41 -> ~69k
            points; a non-vectorized callback dominates panel startup.
        """
        import warnings
        from ngsolve import VoxelCoefficient, CF

        pmin, pmax = self.mesh.ngmesh.bounding_box
        pmin = [pmin[i] for i in range(3)]
        pmax = [pmax[i] for i in range(3)]
        margin = 0.05 * max(pmax[i] - pmin[i] for i in range(3))
        bbox = [[pmin[i] - margin, pmax[i] + margin] for i in range(3)]

        nx = ny = nz = resolution
        x = np.linspace(bbox[0][0], bbox[0][1], nx)
        y = np.linspace(bbox[1][0], bbox[1][1], ny)
        z = np.linspace(bbox[2][0], bbox[2][1], nz)

        xx, yy, zz = np.meshgrid(x, y, z, indexing='ij')
        pts_flat = np.column_stack([xx.ravel(), yy.ravel(), zz.ravel()])

        # Probe: detect whether b_func is array-vectorized.  This is a
        # single, documented back-compat shim (see docstring) -- NOT a
        # silent algorithmic fallback (which CLAUDE.md "No Fallbacks"
        # forbids).
        xs = pts_flat[:, 0]
        ys = pts_flat[:, 1]
        zs = pts_flat[:, 2]
        probe = np.array([0.0, 1.0])
        vectorized = True
        try:
            _ = np.asarray(b_func(probe, probe, probe))
        except (TypeError, ValueError):
            vectorized = False

        if vectorized:
            B_field = np.asarray(b_func(xs, ys, zs), dtype=float)
            B_field = B_field.reshape(-1, 3)
        else:
            warnings.warn(
                "set_source_from_callback: b_func is not NumPy-vectorized; "
                "falling back to np.vectorize wrapper "
                "(~50-200x slower at resolution=41). "
                "Make b_func accept ndarrays of x, y, z and return "
                "shape (..., 3) to remove this warning.",
                DeprecationWarning,
                stacklevel=2,
            )
            vec_b = np.vectorize(b_func, signature='(),(),()->(3)')
            B_field = np.asarray(vec_b(xs, ys, zs), dtype=float)

        start = (bbox[0][0], bbox[1][0], bbox[2][0])
        end = (bbox[0][1], bbox[1][1], bbox[2][1])

        cfs = []
        for comp in range(3):
            data = B_field[:, comp].reshape(nx, ny, nz)
            data = np.ascontiguousarray(data.transpose(2, 1, 0))
            cfs.append(VoxelCoefficient(start, end, data, linear=True))

        self._B_source_cf = CF(tuple(cfs))

    def set_source_cf(self, B_source_cf):
        """Set source field B_s directly from an NGSolve CoefficientFunction.

        Parameters
        ----------
        B_source_cf : ngsolve.CoefficientFunction
            Vector CF of dimension 3 giving B_s in Tesla.
        """
        self._B_source_cf = B_source_cf

    # ------------------------------------------------------------------
    # Linear solver
    # ------------------------------------------------------------------

    def solve_linear(self, dirichlet='default', eps=1e-10, solver='auto'):
        """Solve the linear A_r formulation.

        Finds A_r in HCurl satisfying:
            int(nu * curl(A_r) . curl(v)) = int((nu_0 - nu) * B_s . curl(v))

        Then: B = B_s + curl(A_r), H = nu * B.

        Parameters
        ----------
        dirichlet : str
            Boundary label for Dirichlet BC (A x n = 0).
        eps : float
            Gauge regularization strength.
        solver : str
            'auto' (AMS if available, else BDDC for large, direct for small),
            'ams' (Chebyshev AMS + TaskManager), 'bddc', or 'direct'.
        """
        from ngsolve import (GridFunction, BilinearForm, LinearForm,
                             curl, InnerProduct, dx, Preconditioner,
                             TaskManager)
        from ngsolve.krylovspace import CGSolver

        if self._B_source_cf is None:
            raise RuntimeError("Set source field first (set_source_from_radia)")

        nu_cf = self._build_nu_cf()
        nu_0 = 1.0 / MU_0

        fes = self._make_hcurl_space(dirichlet)

        solver = self._select_solver(fes.ndof, solver)

        A = fes.TrialFunction()
        v = fes.TestFunction()

        all_mats = list(dict.fromkeys(self.mesh.GetMaterials()))

        # Bilinear form: nu * curl(A) . curl(v) + eps*nu_0 * A . v.
        # eps is dimensionless; scaling by nu_0 is essential on the Kelvin
        # ball where nu tends to zero at the compactification centre.
        gauge_coeff = eps * nu_0
        a = BilinearForm(fes)
        for mat in all_mats:
            a += nu_cf * InnerProduct(curl(A), curl(v)) * dx(mat)
            a += gauge_coeff * InnerProduct(A, v) * dx(mat)

        pre_bddc = None
        if solver == 'bddc':
            pre_bddc = Preconditioner(a, 'bddc')
        a.Assemble()

        # RHS: (nu_0 - nu) * B_s . curl(v)
        nu_contrast = nu_0 - nu_cf
        f = LinearForm(fes)
        for mat in self.iron_domains:
            f += nu_contrast * InnerProduct(self._B_source_cf, curl(v)) * dx(mat)
        f.Assemble()

        self._A_gf = GridFunction(fes, name='A_reduced')

        # NOTE: Caller MUST be inside `with TaskManager():` per CLAUDE.md
        # "Caller Wraps, Helper Does NOT" (2026-05-27).
        if solver == 'ams':
            pre = self._setup_ams_preconditioner(a.mat, fes, gauge_coeff)
            inv = CGSolver(mat=a.mat, pre=pre, maxiter=2000,
                           tol=1e-10, printrates=False)
            self._A_gf.vec.data = inv * f.vec
        elif solver == 'bddc':
            inv = CGSolver(mat=a.mat, pre=pre_bddc.mat, maxiter=2000,
                           tol=1e-10, printrates=False)
            self._A_gf.vec.data = inv * f.vec
        else:
            self._A_gf.vec.data = (
                a.mat.Inverse(fes.FreeDofs()) * f.vec)

        self._B_cf = self._B_source_cf + curl(self._A_gf)
        self._H_cf = nu_cf * self._B_cf

        return self._A_gf

    # ------------------------------------------------------------------
    # Newton nonlinear solver (coenergy)
    # ------------------------------------------------------------------

    def solve_nonlinear_newton(self, bh_data, tol=1e-4, maxiter=50,
                               dirichlet='default', verbose=True,
                               solver='auto'):
        """Newton iteration with SymbolicEnergy for nonlinear B-H curve.

        Uses magnetic coenergy w*(B) = integral_0^B H(B') dB' and NGSolve
        automatic differentiation for the Jacobian. Energy-based line
        search guarantees monotone convergence.

        Parameters
        ----------
        bh_data : list of [H, B] pairs
            Tabulated B-H curve. H in A/m, B in Tesla.
        tol : float
            Convergence tolerance (energy residual).
        maxiter : int
            Maximum Newton iterations.
        dirichlet : str
            Boundary label for Dirichlet BC.
        verbose : bool
            Print iteration progress.
        solver : str
            'auto', 'ams', 'bddc', or 'direct'.
        """
        from ngsolve import (HCurl, GridFunction, BilinearForm, BSpline,
                             SymbolicEnergy, InnerProduct, curl, dx,
                             sqrt, Preconditioner, TaskManager)
        from ngsolve.krylovspace import CGSolver

        if self._B_source_cf is None:
            raise RuntimeError("Set source field first")
        if self._kelvin_region:
            raise NotImplementedError(
                "Kelvin reduced-A uses solve_nonlinear() (Picard). "
                "The Newton energy path requires a Kelvin-pulled-back source."
            )

        # Build H(B) BSpline and coenergy w*(B) = integral_0^B H(B') dB'
        # Same pattern as ScalarPotentialSolver but with inverted curve.
        bh = np.array(bh_data)
        H_tab, B_tab = bh[:, 0], bh[:, 1]

        # Sort by B (invert the B-H curve)
        sort_idx = np.argsort(B_tab)
        B_sorted = B_tab[sort_idx].tolist()
        H_sorted = H_tab[sort_idx].tolist()

        # Remove duplicate B values
        B_unique, H_unique = [B_sorted[0]], [H_sorted[0]]
        for i in range(1, len(B_sorted)):
            if B_sorted[i] > B_unique[-1]:
                B_unique.append(B_sorted[i])
                H_unique.append(H_sorted[i])

        # Extend to high B with vacuum slope (dH/dB ~ 1/mu_0)
        B_max_ext = max(B_unique[-1] * 3, 5.0)
        H_max_ext = H_unique[-1] + (1.0 / MU_0) * (B_max_ext - B_unique[-1])
        B_tab_ext = B_unique + [B_max_ext]
        H_tab_ext = H_unique + [H_max_ext]

        # BSpline for H(B): len(knots) = len(values) + 1
        h_of_b_spline = BSpline(2, [0] + B_tab_ext, H_tab_ext)
        # Coenergy: w*(B) = integral_0^B H(B') dB'
        w_star_bspline = h_of_b_spline.Integrate()

        # FE space: HCurl with nograds gauge
        if dirichlet == 'default':
            fes = HCurl(self.mesh, order=self.order, dirichlet='.*',
                        nograds=True)
        else:
            fes = HCurl(self.mesh, order=self.order, dirichlet=dirichlet,
                        nograds=True)

        if verbose:
            print(f"  HCurl DOFs: {fes.ndof}")

        A = fes.TrialFunction()

        # B_total = B_s + curl(A_r)
        B_total = self._B_source_cf + curl(A)
        B2 = InnerProduct(B_total, B_total)
        B_mag_safe = sqrt(B2 + 1e-12)

        # Physical materials
        all_mats = list(set(self.mesh.GetMaterials()))
        phys_mats = [m for m in all_mats if m != self._kelvin_region]
        iron_mats = [m for m in phys_mats if m in self._mu_r_dict]
        air_mats = [m for m in phys_mats if m not in self._mu_r_dict]

        # Energy functional (symmetric for Newton)
        a = BilinearForm(fes, symmetric=True)

        # Air: w*(B) = B^2 / (2*mu_0)
        for mat in air_mats:
            a += SymbolicEnergy(
                (1.0 / (2.0 * MU_0)) * InnerProduct(B_total, B_total),
                definedon=self.mesh.Materials(mat))

        # Iron: w*(|B|) from BSpline
        for mat in iron_mats:
            a += SymbolicEnergy(
                w_star_bspline(B_mag_safe),
                definedon=self.mesh.Materials(mat))

        # Source coupling: -nu_0 * B_s . curl(A_r)
        # Without this, the energy stationarity gives curl(H) = 0 (wrong).
        # With it, stationarity gives curl(H - H_s) = 0, i.e. curl(H_r) = 0 (correct).
        # Derivation: E_total = int w*(|curl A|) - J.A;  A = A_s + A_r;
        #   J.A_r = (1/mu_0)*curl(B_s).A_r = (1/mu_0)*B_s.curl(A_r) + boundary
        nu_0 = 1.0 / MU_0
        a += SymbolicEnergy(
            -nu_0 * InnerProduct(self._B_source_cf, curl(A)))

        # Gauge regularization: eps * |A|^2 / 2
        eps = 1e-10
        gauge_coeff = eps * nu_0
        a += SymbolicEnergy(
            gauge_coeff / 2.0 * InnerProduct(A, A))

        # Newton iteration with energy line search
        sol = GridFunction(fes, name='A_newton')
        sol.vec[:] = 0

        au = sol.vec.CreateVector()
        r = sol.vec.CreateVector()
        w = sol.vec.CreateVector()
        sol_new = sol.vec.CreateVector()

        solver = self._select_solver(fes.ndof, solver)

        pre_bddc = None
        if solver == 'bddc':
            pre_bddc = Preconditioner(a, type='bddc',
                                      inverse='sparsecholesky')

        converged = False
        # NOTE: Caller MUST be inside `with TaskManager():` per CLAUDE.md
        # "Caller Wraps, Helper Does NOT" (2026-05-27).
        for it in range(maxiter):
            E0 = a.Energy(sol.vec)

            a.AssembleLinearization(sol.vec)
            a.Apply(sol.vec, au)
            r.data = -au

            if solver == 'ams':
                eps_gauge = 1e-10
                pre = self._setup_ams_preconditioner(
                    a.mat, fes, eps_gauge)
                inv = CGSolver(mat=a.mat, pre=pre,
                               maxiter=2000, tol=1e-10,
                               printrates=False)
                w.data = inv * r
            elif solver == 'bddc':
                inv = CGSolver(mat=a.mat, pre=pre_bddc.mat,
                               maxiter=2000, tol=1e-10,
                               printrates=False)
                w.data = inv * r
            else:
                inv = a.mat.Inverse(fes.FreeDofs(),
                                    inverse='pardiso')
                w.data = inv * r

            err = InnerProduct(w, r)
            if verbose:
                print(f"   Newton {it}: err = {err:.2e}, "
                      f"E = {E0:.6e}")

            if abs(err) < tol:
                converged = True
                if verbose:
                    print(f"   Converged at Newton iteration {it}")
                break

            # Energy line search
            sol_new.data = sol.vec + w
            E = a.Energy(sol_new)
            tau = 1.0
            while E > E0 and tau > 1e-10:
                tau *= 0.5
                sol_new.data = sol.vec + tau * w
                E = a.Energy(sol_new)
                if verbose and tau < 0.5:
                    print(f"     line search: tau = {tau:.2e}, "
                          f"E = {E:.6e}")

            sol.vec.data = sol_new

        if not converged and verbose:
            print(f"   WARNING: Not converged after {maxiter} iterations")

        # Store results
        self._A_gf = sol
        B_cf = self._B_source_cf + curl(sol)
        self._B_cf = B_cf

        # H from inverted B-H: H = nu_eff * B where nu_eff = H(|B|) / |B|
        B_mag_post = sqrt(InnerProduct(B_cf, B_cf) + 1e-30)
        nu_eff = h_of_b_spline(B_mag_post) / B_mag_post
        iron_indicator = self._build_domain_indicator()
        from ngsolve import CF
        nu_total = iron_indicator * nu_eff + (1 - iron_indicator) * (1.0 / MU_0)
        self._H_cf = nu_total * B_cf

        return sol

    # ------------------------------------------------------------------
    # Picard nonlinear solver
    # ------------------------------------------------------------------

    def solve_nonlinear(self, bh_data, tol=1e-4, maxiter=50, relax=0.3,
                        dirichlet='default', verbose=True, solver='auto'):
        """Picard iteration for nonlinear A_r formulation.

        Each iteration: solve linear -> evaluate |B| -> update nu -> repeat.

        Weak form per iteration:
            (nu_k * curl A_{k+1}, curl v) = ((1/mu_0 - nu_k) * B_s, curl v)

        Parameters
        ----------
        bh_data : list of [H, B] pairs
            Tabulated B-H curve. H in A/m, B in Tesla.
        tol : float
            Convergence tolerance (max |dB|/B_sat).
        maxiter : int
            Maximum Picard iterations.
        relax : float
            Under-relaxation for nu update (0.3 = 30% damping).
        dirichlet : str
            Boundary label for Dirichlet BC.
        verbose : bool
            Print iteration progress.
        solver : str
            'auto', 'ams', 'bddc', or 'direct'.
        """
        from ngsolve import (HCurl, L2, GridFunction, BilinearForm, LinearForm,
                             curl, InnerProduct, Norm, dx, VOL, Preconditioner,
                             TaskManager)
        if self._B_source_cf is None:
            raise RuntimeError("Set source field first")

        nu_of_b, B_sat = _build_nu_of_b_interpolator(bh_data)
        nu_iron_init = nu_of_b(0.0)
        nu_air = 1.0 / MU_0

        # FE space
        fes = self._make_hcurl_space(dirichlet)

        if verbose:
            print(f"  HCurl DOFs: {fes.ndof}")

        A_trial = fes.TrialFunction()
        v = fes.TestFunction()

        # Per-element nu (L2, order=0)
        fes_nu = L2(self.mesh, order=0)
        nu_gf = GridFunction(fes_nu)

        # Initialize nu on physical materials. The Kelvin metric remains a
        # continuous CoefficientFunction and is assembled separately.
        for el in self.mesh.Elements(VOL):
            mat = str(el.mat) if hasattr(el, 'mat') else str(
                self.mesh.GetMaterials()[el.nr])
            if mat in self._mu_r_dict:
                nu_gf.vec[el.nr] = nu_iron_init
            elif mat != self._kelvin_region:
                nu_gf.vec[el.nr] = nu_air

        all_mats = list(dict.fromkeys(self.mesh.GetMaterials()))
        phys_mats = [m for m in all_mats if m != self._kelvin_region]
        kelvin_nu = None
        if self._kelvin_region:
            from radia.kelvin_source import kelvin_nu_factor_3d_cf
            kelvin_nu = nu_air * kelvin_nu_factor_3d_cf(
                self._kelvin_center, self._kelvin_radius)

        A_gf = GridFunction(fes, name='A_picard')
        A_gf.vec[:] = 0

        B_old_arr = None
        converged = False
        final_relative_change = None
        maximum_linear_relative_residual = 0.0
        iterations = 0
        eps = 1e-6
        # Keep the dimensionless gauge convention identical to solve_linear().
        # Scaling by the vacuum reluctivity is essential in the Kelvin ball,
        # where the transformed curl-curl coefficient vanishes at its centre.
        gauge_coeff = eps * nu_air
        solver = self._select_solver(fes.ndof, solver)

        if verbose:
            solver_names = {'ams': 'AMS+CG', 'bddc': 'BDDC+CG',
                            'direct': 'PARDISO SPD'}
            print(f"  Picard iteration (solver: {solver_names.get(solver, solver)}):")

        for it in range(maxiter):
            iterations = it + 1
            # Assemble: nu * curl(A) . curl(v) + eps * A . v
            a = BilinearForm(fes, symmetric=True)
            for mat in phys_mats:
                a += nu_gf * InnerProduct(curl(A_trial), curl(v)) * dx(mat)
                a += gauge_coeff * InnerProduct(A_trial, v) * dx(mat)
            if kelvin_nu is not None:
                a += kelvin_nu * InnerProduct(curl(A_trial), curl(v)) * dx(
                    self._kelvin_region)
                a += gauge_coeff * InnerProduct(A_trial, v) * dx(
                    self._kelvin_region)

            # RHS: (1/mu_0 - nu) * B_s . curl(v)
            nu_contrast = 1.0 / MU_0 - nu_gf
            f = LinearForm(fes)
            for mat in self.iron_domains:
                f += nu_contrast * InnerProduct(
                    self._B_source_cf, curl(v)) * dx(mat)

            if solver == 'bddc':
                pre_bddc = Preconditioner(a, 'bddc')

            # NOTE: Caller MUST be inside `with TaskManager():` per CLAUDE.md
            # "Caller Wraps, Helper Does NOT" (2026-05-27).
            a.Assemble()
            f.Assemble()

            if solver == 'ams':
                from ngsolve.krylovspace import CGSolver
                pre = self._setup_ams_preconditioner(
                    a.mat, fes, gauge_coeff)
                inv = CGSolver(a.mat, pre, maxiter=2000, tol=1e-8,
                               printrates=False)
                A_gf.vec.data = inv * f.vec
            elif solver == 'bddc':
                from ngsolve.krylovspace import CGSolver
                inv = CGSolver(a.mat, pre_bddc.mat, maxiter=2000,
                               tol=1e-8, printrates=False)
                A_gf.vec.data = inv * f.vec
            else:
                A_gf.vec.data = a.mat.Inverse(
                    fes.FreeDofs(), inverse='pardisospd') * f.vec

            linear_residual = f.vec.CreateVector()
            linear_residual.data = f.vec - a.mat * A_gf.vec
            rhs_norm = float(Norm(f.vec))
            linear_relative_residual = (
                float(Norm(linear_residual)) / rhs_norm if rhs_norm > 0.0 else 0.0
            )
            maximum_linear_relative_residual = max(
                maximum_linear_relative_residual, linear_relative_residual
            )
            if not np.isfinite(linear_relative_residual):
                raise RuntimeError(
                    "reduced-A linear solve produced a non-finite residual "
                    f"at Picard iteration {it}"
                )
            if linear_relative_residual > 1.0e-4:
                raise RuntimeError(
                    "reduced-A linear solve failed its relative-residual "
                    f"contract at Picard iteration {it}: "
                    f"{linear_relative_residual:.6e}"
                )

            # B_total = B_s + curl(A_r)
            B_total_cf = self._B_source_cf + curl(A_gf)

            # Update nu in iron from |B|
            B_new_list = []
            for el in self.mesh.Elements(VOL):
                mat = str(el.mat) if hasattr(el, 'mat') else str(
                    self.mesh.GetMaterials()[el.nr])
                if mat not in self._mu_r_dict:
                    continue

                cx, cy, cz = self._element_centroid(el)
                try:
                    mip = self.mesh(cx, cy, cz)
                    B_val = B_total_cf(mip)
                    B_mag = float(np.sqrt(
                        B_val[0]**2 + B_val[1]**2 + B_val[2]**2))
                except Exception:
                    B_mag = 0.0

                if B_mag > 1e-10:
                    nu_new = nu_of_b(B_mag)
                else:
                    nu_new = nu_iron_init

                # Under-relaxation
                alpha = relax
                nu_old = nu_gf.vec[el.nr]
                nu_gf.vec[el.nr] = alpha * nu_new + (1.0 - alpha) * nu_old
                B_new_list.append(B_mag)

            B_new_arr = np.array(B_new_list)

            # Convergence check
            if B_old_arr is not None and len(B_old_arr) == len(B_new_arr):
                max_change = np.max(np.abs(B_new_arr - B_old_arr)) / B_sat
                final_relative_change = float(max_change)
                if verbose:
                    mip_0 = self.mesh(0, 0, 0)
                    Bz_now = B_total_cf(mip_0)[2]
                    print(f"   iter {it}: max |dB|/B_sat = {max_change:.2e}, "
                          f"Bz = {Bz_now * 1e3:.1f} mT")
                if max_change < tol:
                    converged = True
                    if verbose:
                        print(f"   Converged at iteration {it}")
                    break
            else:
                if verbose:
                    mip_0 = self.mesh(0, 0, 0)
                    Bz_now = B_total_cf(mip_0)[2]
                    print(f"   iter {it}: initial, "
                          f"Bz = {Bz_now * 1e3:.1f} mT")

            B_old_arr = B_new_arr.copy()

        self._last_nonlinear_stats = {
            'converged': bool(converged),
            'iterations': int(iterations),
            'final_relative_change': final_relative_change,
            'tolerance': float(tol),
            'maximum_iterations': int(maxiter),
            'maximum_linear_relative_residual': maximum_linear_relative_residual,
            'direct_inverse': 'pardisospd' if solver == 'direct' else None,
        }
        if not converged and verbose:
            print(f"   WARNING: Not converged after {maxiter} iterations")

        # Store results
        self._A_gf = A_gf
        B_cf = self._B_source_cf + curl(A_gf)
        self._B_cf = B_cf
        if kelvin_nu is None:
            self._H_cf = nu_gf * B_cf
        else:
            self._H_cf = self.mesh.MaterialCF(
                {self._kelvin_region: kelvin_nu}, default=nu_gf) * B_cf

        return A_gf

    # ------------------------------------------------------------------
    # Hysteresis solver (Hantila polarization method for A_r)
    # ------------------------------------------------------------------

    def solve_hysteresis(self, mat_factory, tol=1e-4, maxiter=50, alpha=500.0,
                         dirichlet='default', verbose=True, solver='auto',
                         relax=0.0):
        """Hantila polarization iteration for A_r with energy hysteresis.

        Uses the Hantila (1975) polarization method adapted for the
        vector potential formulation. The LHS is assembled once (constant
        reluctivity nu_alpha in iron, nu_0 in air). The RHS is updated
        each iteration from the polarization residual R = M - alpha*H.

        H is derived from B without needing Inverse(B->H):
            H = (B/mu_0 - R_prev) / (1 + alpha)

        Parameters
        ----------
        mat_factory : callable
            Factory: mat_factory() -> int (Radia material handle).
        tol : float
            Convergence tolerance (max |dB| / B_sat).
        maxiter : int
            Maximum Picard iterations.
        alpha : float
            Polarization parameter (>= max susceptibility). Default 500.
        dirichlet : str
            Boundary label for Dirichlet BC (A x n = 0).
        verbose : bool
            Print iteration progress.
        solver : str
            'auto', 'ams', 'bddc', or 'direct'.
        relax : float
            Under-relaxation for R update (0.0 = full step, 0.5 = half).
        """
        import radia as rad
        from ngsolve import (HCurl, VectorL2, GridFunction, BilinearForm,
                             LinearForm, curl, InnerProduct, dx, VOL,
                             Preconditioner, TaskManager, CF)

        if self._B_source_cf is None:
            raise RuntimeError("Set source field first")

        nu_0 = 1.0 / MU_0
        nu_alpha = 1.0 / (MU_0 * (1 + alpha))

        # Per-element material handles (persistent across time steps)
        n_iron = 0
        for el in self.mesh.Elements(VOL):
            if str(el.mat) in self._mu_r_dict:
                n_iron += 1

        if not hasattr(self, '_hys_handles') or self._hys_handles is None:
            self._hys_handles = [mat_factory() for _ in range(n_iron)]
            if verbose:
                print(f"   Created {n_iron} per-element material handles")

        handles = self._hys_handles

        # FE space: HCurl with nograds gauge
        if dirichlet == 'default':
            fes = HCurl(self.mesh, order=self.order, dirichlet='.*',
                        nograds=True)
        else:
            fes = HCurl(self.mesh, order=self.order, dirichlet=dirichlet,
                        nograds=True)

        if verbose:
            print(f"   HCurl DOFs: {fes.ndof}")

        A_trial = fes.TrialFunction()
        v = fes.TestFunction()

        # Polarization residual R = M - alpha*H (VectorL2 order=0)
        fes_vec = VectorL2(self.mesh, order=0)
        R_gf = GridFunction(fes_vec, name='R_polar')
        R_gf.vec[:] = 0

        # Material classification
        all_mats = list(set(self.mesh.GetMaterials()))
        phys_mats = [m for m in all_mats if m != self._kelvin_region]
        iron_mats = [m for m in phys_mats if m in self._mu_r_dict]
        air_mats = [m for m in phys_mats if m not in self._mu_r_dict]

        # Iron element indices and centroids
        iron_els = []
        for el in self.mesh.Elements(VOL):
            if str(el.mat) in self._mu_r_dict:
                iron_els.append((el.nr, self._element_centroid(el)))

        # Bilinear form: CONSTANT LHS (assembled once)
        eps = 1e-10
        a = BilinearForm(fes)
        for mat in iron_mats:
            a += nu_alpha * InnerProduct(curl(A_trial), curl(v)) * dx(mat)
        for mat in air_mats:
            a += nu_0 * InnerProduct(curl(A_trial), curl(v)) * dx(mat)
        a += eps * InnerProduct(A_trial, v) * dx  # gauge

        solver_type = self._select_solver(fes.ndof, solver)
        pre_bddc = None
        if solver_type == 'bddc':
            pre_bddc = Preconditioner(a, 'bddc')
        a.Assemble()

        nel = len(list(self.mesh.Elements(VOL)))

        # Save material states before iteration
        saved_states = [rad.MatHysSaveState(h) for h in handles]

        A_gf = GridFunction(fes, name='A_hysteresis')
        A_gf.vec[:] = 0
        B_old = np.zeros(len(iron_els))

        converged = False
        for iteration in range(maxiter):
            # RHS: contrast + polarization residual
            f = LinearForm(fes)
            for mat in iron_mats:
                f += (nu_0 - nu_alpha) * InnerProduct(
                    self._B_source_cf, curl(v)) * dx(mat)
                f += (1.0 / (1 + alpha)) * InnerProduct(
                    R_gf, curl(v)) * dx(mat)
            f.Assemble()

            # Solve (LHS pre-assembled).
            # NOTE: Caller MUST be inside `with TaskManager():` per CLAUDE.md
            # "Caller Wraps, Helper Does NOT" (2026-05-27).
            if solver_type == 'ams':
                from ngsolve.krylovspace import CGSolver
                pre = self._setup_ams_preconditioner(a.mat, fes, eps)
                inv = CGSolver(a.mat, pre, maxiter=2000, tol=1e-10,
                               printrates=False)
                A_gf.vec.data = inv * f.vec
            elif solver_type == 'bddc':
                from ngsolve.krylovspace import CGSolver
                inv = CGSolver(a.mat, pre_bddc.mat, maxiter=2000,
                               tol=1e-10, printrates=False)
                A_gf.vec.data = inv * f.vec
            else:
                A_gf.vec.data = a.mat.Inverse(
                    fes.FreeDofs(), inverse='pardiso') * f.vec

            # B = B_s + curl(A_r) at centroids
            B_total_cf = self._B_source_cf + curl(A_gf)

            # Restore states before Forward(H)
            for h, s in zip(handles, saved_states):
                rad.MatHysRestoreState(h, s)

            # Derive H from B, Forward(H) -> M, update R
            B_new = np.zeros(len(iron_els))
            for idx, (el_nr, centroid) in enumerate(iron_els):
                try:
                    mip = self.mesh(*centroid)
                    B_val = [float(B_total_cf(mip)[i]) for i in range(3)]
                except Exception:
                    B_val = [0.0, 0.0, 0.0]

                # H = (B/mu_0 - R_prev) / (1 + alpha)
                R_prev = [R_gf.vec[d * nel + el_nr] for d in range(3)]
                H_val = [(B_val[i] / MU_0 - R_prev[i]) / (1 + alpha)
                         for i in range(3)]

                # Forward(H) -> M
                M_val = list(rad.MatMvsH(handles[idx], 'm', H_val))

                B_mag = msqrt(B_val[0]**2 + B_val[1]**2 + B_val[2]**2)
                B_new[idx] = B_mag

                # R = M - alpha * H (with optional under-relaxation)
                omega = 1.0 - relax
                for d in range(3):
                    R_target = M_val[d] - alpha * H_val[d]
                    if relax > 0 and iteration > 0:
                        R_gf.vec[d * nel + el_nr] = (
                            (1 - omega) * R_gf.vec[d * nel + el_nr]
                            + omega * R_target)
                    else:
                        R_gf.vec[d * nel + el_nr] = R_target

            # Convergence check
            if iteration > 0:
                max_change = np.max(np.abs(B_new - B_old)) / max(2.0, 1e-10)
                if verbose:
                    print(f"   iter {iteration}: "
                          f"max |dB|/B_sat = {max_change:.2e}")
                if max_change < tol:
                    converged = True
                    if verbose:
                        print(f"   Converged at iteration {iteration}")
                    break
            elif verbose:
                print(f"   iter 0: initial solve")

            B_old[:] = B_new

        if not converged and verbose:
            print(f"   WARNING: Not converged after {maxiter} iterations")

        # Commit converged states for next time step
        if converged:
            for h in handles:
                rad.MatHysCommitState(h)

        # Store results
        self._A_gf = A_gf
        B_cf = self._B_source_cf + curl(A_gf)
        self._B_cf = B_cf

        # H = nu_alpha*B - R/(1+alpha) in iron, nu_0*B in air
        iron_indicator = self._build_domain_indicator()
        H_iron = nu_alpha * B_cf - (1.0 / (1 + alpha)) * R_gf
        H_air = nu_0 * B_cf
        self._H_cf = iron_indicator * H_iron + (1 - iron_indicator) * H_air

        # Save per-element M for to_radia() pipeline
        self._M_per_element = {}
        for idx, (el_nr, centroid) in enumerate(iron_els):
            try:
                mip = self.mesh(*centroid)
                B_val = [float(B_cf(mip)[i]) for i in range(3)]
            except Exception:
                B_val = [0.0, 0.0, 0.0]
            R_prev = [R_gf.vec[d * nel + el_nr] for d in range(3)]
            H_val = [(B_val[i] / MU_0 - R_prev[i]) / (1 + alpha)
                     for i in range(3)]
            M_val = list(rad.MatMvsH(handles[idx], 'm', H_val))
            self._M_per_element[el_nr] = M_val

        return A_gf

    # ------------------------------------------------------------------
    # Result access
    # ------------------------------------------------------------------

    def get_B(self):
        """Get B field CoefficientFunction (Tesla, dim=3).

        B = B_s + curl(A_r).
        """
        if self._B_cf is None:
            raise RuntimeError("Call solve_linear() or solve_nonlinear_newton() first")
        return self._B_cf

    def get_H(self):
        """Get H field CoefficientFunction (A/m, dim=3).

        H = nu * B.
        """
        if self._H_cf is None:
            raise RuntimeError("Call solve_linear() or solve_nonlinear_newton() first")
        return self._H_cf

    def get_A(self):
        """Get reduced vector potential GridFunction (HCurl)."""
        if self._A_gf is None:
            raise RuntimeError("Call solve_linear() or solve_nonlinear_newton() first")
        return self._A_gf

    def project_to_hdiv(self, order=None):
        """Project B field onto HDiv GridFunction (ensures div(B)=0).

        Parameters
        ----------
        order : int, optional
            HDiv polynomial order. Default: same as solver order.

        Returns
        -------
        ngsolve.GridFunction
            B field in HDiv space.
        """
        from ngsolve import HDiv, GridFunction

        if self._B_cf is None:
            raise RuntimeError("Call solve first")

        if order is None:
            order = self.order
        fes_hdiv = HDiv(self.mesh, order=order)
        B_gf = GridFunction(fes_hdiv, name='B')
        B_gf.Set(self._B_cf)
        return B_gf

    def get_M_per_element(self):
        """Get per-element magnetization from hysteresis solve.

        Returns
        -------
        dict
            {element_nr: [Mx, My, Mz]} in A/m.
        """
        if not hasattr(self, '_M_per_element') or self._M_per_element is None:
            raise RuntimeError("Call solve_hysteresis() first")
        return self._M_per_element

    def to_radia(self, coil=None):
        """Convert solved magnetization to Radia objects for analytical field.

        Creates Radia ObjTetrahedron/ObjHexahedron elements with per-element
        magnetization from solve_hysteresis(). The resulting Radia objects
        evaluate B using exact analytical surface charge formulas -- no mesh
        needed in the gap region.

        Parameters
        ----------
        coil : int, optional
            Radia coil object handle to combine with iron contribution.

        Returns
        -------
        int
            Radia container object handle.
        """
        import radia as rad
        from radia.netgen_mesh_import import netgen_mesh_to_radia

        M_dict = self.get_M_per_element()

        def material_func(el_idx):
            if el_idx in M_dict:
                return {'magnetization': M_dict[el_idx]}
            return {'magnetization': [0, 0, 0]}

        iron_domains = list(self._mu_r_dict.keys())
        container = netgen_mesh_to_radia(
            self.mesh, material=material_func,
            material_filter=iron_domains, verbose=False)

        if coil is not None:
            container = rad.ObjCnt([container, coil])

        return container

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _build_nu_cf(self):
        """Build domain-wise nu CoefficientFunction (linear).

        nu = 1/(mu_0 * mu_r) for iron, 1/mu_0 for air.
        """
        nu_0 = 1.0 / MU_0
        nu_dict = {d: nu_0 / mr for d, mr in self._mu_r_dict.items()}
        if self._kelvin_region:
            from radia.kelvin_source import kelvin_nu_factor_3d_cf
            nu_dict[self._kelvin_region] = nu_0 * kelvin_nu_factor_3d_cf(
                self._kelvin_center, self._kelvin_radius)
        return self.mesh.MaterialCF(nu_dict, default=nu_0)

    def _make_hcurl_space(self, dirichlet):
        """Create the NGSolve-native HCurl space for the configured domain."""
        from ngsolve import HCurl, Periodic

        if not self._kelvin_region:
            boundary = '.*' if dirichlet == 'default' else dirichlet
            return HCurl(
                self.mesh, order=self.order, dirichlet=boundary, nograds=True)

        from radia.kelvin_identify_ngsolve import has_kelvin_identification
        if not has_kelvin_identification(self.mesh):
            raise RuntimeError(
                "Kelvin HCurl requires kelvin_int/kelvin_ext point "
                "identifications in the .vol mesh"
            )
        kwargs = {"order": self.order, "nograds": True}
        if dirichlet not in ('default', 'GND'):
            kwargs["dirichlet"] = dirichlet
        return Periodic(HCurl(self.mesh, **kwargs))

    def _build_domain_indicator(self):
        """Build indicator CF: 1.0 in iron, 0.0 elsewhere."""
        iron_dict = {d: 1.0 for d in self.iron_domains}
        return self.mesh.MaterialCF(iron_dict, default=0.0)

    def _compute_physical_bbox(self):
        """Compute bounding box of physical (non-Kelvin) domains."""
        from ngsolve import VOL
        pmin = [1e30, 1e30, 1e30]
        pmax = [-1e30, -1e30, -1e30]
        for el in self.mesh.Elements(VOL):
            mat = (el.mat if hasattr(el, 'mat')
                   else str(self.mesh.GetMaterials()[el.nr]))
            if mat == self._kelvin_region:
                continue
            for v in el.vertices:
                pt = self.mesh.vertices[v.nr].point
                for i in range(3):
                    pmin[i] = min(pmin[i], pt[i])
                    pmax[i] = max(pmax[i], pt[i])
        margin = 0.1 * max(pmax[i] - pmin[i] for i in range(3))
        return [[pmin[i] - margin, pmax[i] + margin] for i in range(3)]

    def _setup_ams_preconditioner(self, a_mat, fes, h1_mass_coeff):
        """Set up AMS preconditioner with Chebyshev smoother.

        Parameters
        ----------
        a_mat : BaseMatrix
            Assembled HCurl system matrix.
        fes : HCurl FESpace
            HCurl finite element space.
        h1_mass_coeff : CoefficientFunction or float
            Coefficient for H1 gradient subspace problem.
            Magnetostatic: eps (gauge regularization).
            Eddy current: eps * nu + omega * sigma.
        """
        import radia.sparsesolv_ngsolve as ssn

        opts = self._ams_options
        if int(self.order) != 1:
            raise ValueError(
                "AMS requires HCurl order=1; use solver='bddc' for order>=2"
            )

        # Cache G_mat and h1_fes (invariant across Newton/Picard iterations)
        if self._ams_cache.get('fes_id') != id(fes):
            G_mat, h1_fes = fes.CreateGradient()
            self._ams_cache = {
                'G_mat': G_mat, 'h1_fes': h1_fes, 'fes_id': id(fes),
            }
        G_mat = self._ams_cache['G_mat']
        h1_fes = self._ams_cache['h1_fes']

        if int(h1_fes.ndof) != int(self.mesh.nv):
            raise RuntimeError(
                "AMS order-1 gradient space must have one H1 DOF per mesh "
                f"vertex (got {h1_fes.ndof} DOFs for {self.mesh.nv} vertices)"
            )
        coordinates = [self.mesh.vertices[index].point
                       for index in range(self.mesh.nv)]
        coord_x = [float(point[0]) for point in coordinates]
        coord_y = [float(point[1]) for point in coordinates]
        coord_z = [float(point[2]) for point in coordinates]

        return ssn.HypreBasedAMSPreconditioner(
            mat=a_mat,
            grad_mat=G_mat,
            freedofs=fes.FreeDofs(),
            coord_x=coord_x,
            coord_y=coord_y,
            coord_z=coord_z,
            num_smooth=opts['num_smooth'],
            cycle_type=1,
            print_level=0,
        )

    def _select_solver(self, fes_ndof, solver):
        """Select solver type based on DOF count and availability."""
        if getattr(self, '_kelvin_region', None):
            if solver in ('ams', 'bddc'):
                raise ValueError(
                    f"solver='{solver}' is not supported for Periodic Kelvin "
                    "HCurl; use solver='direct'"
                )
            return 'direct' if solver == 'auto' else solver
        if solver != 'auto':
            if solver == 'ams' and int(self.order) != 1:
                raise ValueError(
                    "solver='ams' requires HCurl order=1; use solver='bddc' "
                    "for order>=2"
                )
            return solver
        if fes_ndof <= 200_000:
            return 'direct'
        if int(self.order) != 1:
            return 'bddc'
        try:
            import radia.sparsesolv_ngsolve as ssn
            if hasattr(ssn, 'HypreBasedAMSPreconditioner'):
                return 'ams'
            return 'bddc'
        except ImportError:
            return 'bddc'

    def _solve_system(self, a, f, fes, gf, pre, use_iterative):
        """Solve assembled linear system (direct or iterative)."""
        if use_iterative:
            from ngsolve.krylovspace import CGSolver
            inv = CGSolver(mat=a.mat, pre=pre.mat, maxiter=2000,
                           printrates=False, tol=1e-10)
            gf.vec.data = inv * f.vec
        else:
            gf.vec.data = a.mat.Inverse(fes.FreeDofs()) * f.vec

    def _element_centroid(self, el):
        """Compute centroid of a volume element."""
        verts = el.vertices
        coords = np.zeros(3)
        for v in verts:
            pt = self.mesh.vertices[v.nr].point
            coords += np.array([pt[0], pt[1], pt[2]])
        coords /= len(verts)
        return tuple(coords)
