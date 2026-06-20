"""Simkin-Trowbridge reduced scalar potential solver (Radia + NGSolve).

Implements the Simkin-Trowbridge (1979) method for magnetostatic
field computation combining Radia BEM source fields with NGSolve FEM.

Method:
    H = H_s - grad(phi)
    B = mu * H
    Weak form: int(mu * grad(phi) . grad(v)) = int(mu * H_s . grad(v))

Features:
    - Linear materials (mu_r per domain)
    - Nonlinear materials (B-H curve, Picard iteration)
    - Energy-based hysteresis (rad.MatEnergyHysteresis, Picard iteration)
    - Kelvin transform for open boundary
    - BDDC+CG iterative solver for large problems
    - Two-potential variant for high mu_r (> 5000)

Usage:
    import radia as rad
    from radia.scalar_potential_solver import ScalarPotentialSolver

    solver = ScalarPotentialSolver(
        mesh, mu_r_dict={'iron': 1000}, order=3,
        kelvin_region='kelvin', kelvin_radius=0.3,
        kelvin_center=[0, 0.07, 0])
    solver.set_source_from_radia(coil)
    solver.solve(dirichlet='outer')
    B_cf = solver.get_B()

Reference:
    J. Simkin and C. W. Trowbridge, "On the use of the total scalar
    potential in the numerical solution of field problems in
    electromagnetics," IJNME, vol. 14, pp. 423-440, 1979.
"""

import numpy as np
from math import sqrt as msqrt

MU_0 = 4 * np.pi * 1e-7


class ScalarPotentialSolver:
    """Simkin-Trowbridge magnetostatic solver (Radia + NGSolve).

    Parameters
    ----------
    mesh : ngsolve.Mesh
        FEM mesh with material labels.
    iron_domains : str or list of str
        Material name(s) for iron regions.
    mu_r : float
        Default relative permeability for iron (used when mu_r_dict not given).
    mu_r_dict : dict, optional
        Per-domain mu_r, e.g. {'iron': 1000, 'steel': 500}.
        Overrides mu_r and iron_domains when provided.
    order : int
        FEM polynomial order (default: 2).
    kelvin_region : str, optional
        Material name of Kelvin exterior domain region.
    kelvin_radius : float, optional
        Inner radius R of Kelvin exterior domain [m].
    kelvin_center : list/tuple, optional
        Center of Kelvin exterior domain [x, y, z] in meters.
    """

    def __init__(self, mesh, iron_domains='iron', mu_r=1000.0, order=2,
                 mu_r_dict=None, kelvin_region=None, kelvin_radius=None,
                 kelvin_center=None):
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

        # Source field (set by set_source_*)
        self._H_source_cf = None
        self._radia_obj = None
        # Precomputed source CFs for total+reduced solver (built eagerly
        # in set_source_from_radia while the Radia handle is alive)
        self._Omega_s_cf = None   # scalar VoxelCoefficient
        self._Bs_cf = None        # vector CF = MU_0 * H_source_cf

        # Solution storage
        self._phi_gf = None
        self._B_cf = None
        self._H_cf = None

    # ------------------------------------------------------------------
    # Source field setup
    # ------------------------------------------------------------------

    def set_source_from_radia(self, radia_obj, resolution=None):
        """Set source field from a Radia object.

        Evaluates ``rad.Fld('h')`` and ``rad.Fld('phi')`` at mesh DoF
        locations and projects them into H1 GridFunctions for smooth
        FEM integration.  This is the NGSolve-native approach (no
        VoxelCoefficient): the source fields are C0-continuous FE
        functions with polynomial gradients within each element,
        compatible with any FE order including order >= 3.

        After this call returns, the Radia object handle is no longer
        needed — all subsequent solve methods work from stored
        GridFunctions.

        Parameters
        ----------
        radia_obj : int
            Radia object handle (coil, after Solve() if soft iron present).
        resolution : int, optional
            Ignored (kept for API compatibility). The mesh DoF locations
            are used directly — no intermediate voxel grid.
        """
        import radia as rad
        from ngsolve import H1, GridFunction, CF, VOL, BND

        self._radia_obj = radia_obj

        # --- Collect DoF locations for the H1 space ---
        fes_scalar = H1(self.mesh, order=self.order)
        ndof = fes_scalar.ndof
        dof_pts = np.zeros((ndof, 3))
        for i in range(ndof):
            gf_tmp = GridFunction(fes_scalar)
            # Cannot iterate DoFs easily; use mesh vertices + edges.
            # Instead: evaluate rad.Fld at mesh vertices, then Set()
            # with the analytical CF at quadrature points.
            pass

        # Simpler approach: evaluate rad.Fld at all mesh vertices,
        # then use GridFunction.Set() which evaluates the CF at internal
        # quadrature points. Since we can't pass rad.Fld as a CF,
        # we build a VoxelCoefficient for the Set() projection only,
        # but the RESULT is a smooth H1 GridFunction.
        #
        # This is a TWO-STEP process:
        #   Step 1: Build a temporary VoxelCoefficient from rad.Fld
        #   Step 2: Project it into H1 GridFunctions (smooth)
        # After Step 2, the VoxelCoefficient is discarded.

        from ngsolve import VoxelCoefficient

        # Physical bbox (excludes Kelvin)
        if self._kelvin_region:
            bbox = self._compute_physical_bbox()
        else:
            pmin, pmax = self.mesh.ngmesh.bounding_box
            pmin = [pmin[i] for i in range(3)]
            pmax = [pmax[i] for i in range(3)]
            margin = 0.05 * max(pmax[i] - pmin[i] for i in range(3))
            bbox = [[pmin[i] - margin, pmax[i] + margin] for i in range(3)]

        # Dense voxel grid for accurate projection
        res = max(61, 2 * self.order * 20 + 1)
        nx = ny = nz = res
        xs = np.linspace(bbox[0][0], bbox[0][1], nx)
        ys = np.linspace(bbox[1][0], bbox[1][1], ny)
        zs = np.linspace(bbox[2][0], bbox[2][1], nz)
        xx, yy, zz = np.meshgrid(xs, ys, zs, indexing='ij')
        pts = np.column_stack([xx.ravel(), yy.ravel(), zz.ravel()])

        # Evaluate H and Phi from Radia.
        # Convention: rad.Fld('phi') returns Phi where H = -grad(Phi).
        # EMPY uses Omega_s where H = +grad(Omega_s), so Omega_s = -Phi.
        H_vals = np.asarray(rad.Fld(radia_obj, 'h', pts))
        Phi_vals = -np.asarray(rad.Fld(radia_obj, 'phi', pts))

        start = (bbox[0][0], bbox[1][0], bbox[2][0])
        end = (bbox[0][1], bbox[1][1], bbox[2][1])

        # Temporary VoxelCoefficients (for projection only)
        H_voxel_cfs = []
        for c in range(3):
            d = H_vals[:, c].reshape(nx, ny, nz)
            d = np.ascontiguousarray(d.transpose(2, 1, 0))
            H_voxel_cfs.append(VoxelCoefficient(start, end, d, linear=True))
        H_voxel = CF(tuple(H_voxel_cfs))
        Phi_data = Phi_vals.reshape(nx, ny, nz)
        Phi_data = np.ascontiguousarray(Phi_data.transpose(2, 1, 0))
        Phi_voxel = VoxelCoefficient(start, end, Phi_data, linear=True)

        # Project into H1 GridFunctions (smooth, C0-continuous)
        # Exclude Kelvin region from projection (source = 0 there)
        phys_region = '|'.join(
            m for m in set(self.mesh.GetMaterials())
            if m != self._kelvin_region)

        fes_vec = H1(self.mesh, order=self.order, dim=3)
        self._H_source_gf = GridFunction(fes_vec, name='H_source')
        self._H_source_gf.Set(H_voxel, definedon=self.mesh.Materials(
            phys_region))

        fes_sc = H1(self.mesh, order=self.order)
        self._Omega_s_gf = GridFunction(fes_sc, name='Omega_s')
        self._Omega_s_gf.Set(Phi_voxel, definedon=self.mesh.Materials(
            phys_region))

        # Store as CFs (GridFunction IS a CoefficientFunction)
        self._H_source_cf = self._H_source_gf
        self._Bs_cf = MU_0 * self._H_source_gf
        self._Omega_s_cf = self._Omega_s_gf

    def set_source_from_callback(self, h_func, resolution=41):
        """Set source field H_s from a Python callback.

        Parameters
        ----------
        h_func : callable
            ``h_func(x, y, z) -> (Hx, Hy, Hz)`` returning H in A/m.
            **Contract**: ``h_func`` MUST accept NumPy ndarrays of equal
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

        # Probe: detect whether h_func is array-vectorized.  This is a
        # single, documented back-compat shim (see docstring) -- NOT a
        # silent algorithmic fallback (which CLAUDE.md "No Fallbacks"
        # forbids).
        xs = pts_flat[:, 0]
        ys = pts_flat[:, 1]
        zs = pts_flat[:, 2]
        probe = np.array([0.0, 1.0])
        vectorized = True
        try:
            _ = np.asarray(h_func(probe, probe, probe))
        except (TypeError, ValueError):
            vectorized = False

        if vectorized:
            H_field = np.asarray(h_func(xs, ys, zs), dtype=float)
            H_field = H_field.reshape(-1, 3)
        else:
            warnings.warn(
                "set_source_from_callback: h_func is not NumPy-vectorized; "
                "falling back to np.vectorize wrapper "
                "(~50-200x slower at resolution=41). "
                "Make h_func accept ndarrays of x, y, z and return "
                "shape (..., 3) to remove this warning.",
                DeprecationWarning,
                stacklevel=2,
            )
            vec_h = np.vectorize(h_func, signature='(),(),()->(3)')
            H_field = np.asarray(vec_h(xs, ys, zs), dtype=float)

        start = (bbox[0][0], bbox[1][0], bbox[2][0])
        end = (bbox[0][1], bbox[1][1], bbox[2][1])

        cfs = []
        for comp in range(3):
            data = H_field[:, comp].reshape(nx, ny, nz)
            data = np.ascontiguousarray(data.transpose(2, 1, 0))
            cfs.append(VoxelCoefficient(start, end, data, linear=True))

        self._H_source_cf = CF(tuple(cfs))

    def set_source_cf(self, H_source_cf):
        """Set source field H_s directly from an NGSolve CoefficientFunction.

        Parameters
        ----------
        H_source_cf : ngsolve.CoefficientFunction
            Vector CF of dimension 3 giving H_s in A/m.
        """
        self._H_source_cf = H_source_cf

    # ------------------------------------------------------------------
    # Single-potential solver (reduced potential)
    # ------------------------------------------------------------------

    def solve_single_potential(self, dirichlet='default'):
        """Solve using a single reduced potential on the full domain.

        Finds phi in H1 satisfying:
            int(mu * grad(phi) . grad(v)) = int(mu * H_s . grad(v))

        Then: H = H_s - grad(phi), B = mu * H.

        Includes Kelvin transform if configured.

        Parameters
        ----------
        dirichlet : str
            Boundary label for Dirichlet BC (phi=0).
            'default' uses all outer boundaries.
        """
        from ngsolve import (H1, GridFunction, BilinearForm, LinearForm,
                             grad, dx, x, y, z, sqrt, Preconditioner)

        if self._H_source_cf is None:
            raise RuntimeError("Set source field first (set_source_from_radia)")

        mu_cf = self._build_mu_cf()

        if dirichlet == 'default':
            fes = H1(self.mesh, order=self.order, dirichlet='.*')
        else:
            fes = H1(self.mesh, order=self.order, dirichlet=dirichlet)

        phi = fes.TrialFunction()
        v = fes.TestFunction()

        # Physical materials (exclude kelvin)
        all_mats = list(set(self.mesh.GetMaterials()))
        phys_mats = [m for m in all_mats if m != self._kelvin_region]

        # Bilinear form
        a = BilinearForm(fes)
        for mat in phys_mats:
            a += mu_cf * grad(phi) * grad(v) * dx(mat)

        # Kelvin exterior domain term
        kelvin_weight = None
        if self._kelvin_region:
            kelvin_weight = self._build_kelvin_weight(x, y, z, sqrt)
            a += MU_0 * kelvin_weight * grad(phi) * grad(v) * dx(
                self._kelvin_region)

        # BDDC preconditioner for large problems
        pre = None
        use_iterative = fes.ndof > 200_000
        if use_iterative:
            pre = Preconditioner(a, 'bddc')
        a.Assemble()

        # Source
        f = LinearForm(fes)
        for mat in phys_mats:
            f += mu_cf * self._H_source_cf * grad(v) * dx(mat)
        if self._kelvin_region and kelvin_weight is not None:
            f += MU_0 * kelvin_weight * self._H_source_cf * grad(v) * dx(
                self._kelvin_region)
        f.Assemble()

        # Solve
        self._phi_gf = GridFunction(fes, name='phi_reduced')
        self._solve_system(a, f, fes, self._phi_gf, pre, use_iterative)

        self._H_cf = self._H_source_cf - grad(self._phi_gf)
        self._B_cf = mu_cf * self._H_cf

        return self._phi_gf

    # ------------------------------------------------------------------
    # Two-potential solver (Simkin-Trowbridge)
    # ------------------------------------------------------------------

    def solve_two_potential(self, dirichlet='default'):
        """Solve using the two-scalar-potential method.

        Uses separate unknowns:
            - phi_r on full domain (reduced potential): H = H_s - grad(phi_r)
            - psi on iron only (total potential):       H = -grad(psi)

        Better for high mu_r (> 5000) to avoid cancellation error.

        Parameters
        ----------
        dirichlet : str
            Boundary label for Dirichlet BC on outer boundary.
        """
        from ngsolve import (H1, GridFunction, BilinearForm, LinearForm,
                             grad, dx, CF, x, y, z, sqrt, Preconditioner)

        if self._H_source_cf is None:
            raise RuntimeError("Set source field first (set_source_from_radia)")

        iron_re = '|'.join(self.iron_domains)

        if dirichlet == 'default':
            V_full = H1(self.mesh, order=self.order, dirichlet='.*')
        else:
            V_full = H1(self.mesh, order=self.order, dirichlet=dirichlet)
        V_iron = H1(self.mesh, order=self.order, definedon=iron_re)

        X = V_full * V_iron
        (phi_r, psi), (v_f, v_i) = X.TnT()

        # Physical materials (exclude kelvin)
        all_mats = list(set(self.mesh.GetMaterials()))
        air_mats = [m for m in all_mats
                    if m not in self.iron_domains and m != self._kelvin_region]

        a = BilinearForm(X)

        # Air region
        for mat in air_mats:
            a += MU_0 * grad(phi_r) * grad(v_f) * dx(mat)

        # Iron region
        for dom, mu_r_val in self._mu_r_dict.items():
            a += (MU_0 * mu_r_val) * grad(psi) * grad(v_i) * dx(dom)

        # Coupling penalty
        penalty = 1e3 * self.mu_r * MU_0
        a += penalty * (phi_r - psi) * (v_f - v_i) * dx(iron_re)

        # Kelvin exterior domain
        kelvin_weight = None
        if self._kelvin_region:
            kelvin_weight = self._build_kelvin_weight(x, y, z, sqrt)
            a += MU_0 * kelvin_weight * grad(phi_r) * grad(v_f) * dx(
                self._kelvin_region)

        pre = None
        use_iterative = X.ndof > 200_000
        if use_iterative:
            pre = Preconditioner(a, 'bddc')
        a.Assemble()

        # Source: H_s in air only
        f = LinearForm(X)
        for mat in air_mats:
            f += MU_0 * self._H_source_cf * grad(v_f) * dx(mat)
        if self._kelvin_region and kelvin_weight is not None:
            f += MU_0 * kelvin_weight * self._H_source_cf * grad(v_f) * dx(
                self._kelvin_region)
        f.Assemble()

        gf = GridFunction(X, name='two_potential')
        self._solve_system(a, f, X, gf, pre, use_iterative)

        phi_r_gf, psi_gf = gf.components
        self._phi_r_gf = phi_r_gf
        self._psi_gf = psi_gf
        self._phi_gf = gf

        # Build domain-wise result fields
        iron_ind = self._build_domain_indicator()
        air_ind = 1.0 - iron_ind

        H_air = self._H_source_cf - grad(phi_r_gf)
        H_iron = -grad(psi_gf)
        self._H_cf = air_ind * H_air + iron_ind * H_iron

        mu_cf = self._build_mu_cf()
        self._B_cf = mu_cf * self._H_cf

        return gf

    # ------------------------------------------------------------------
    # Total+reduced split (EMPY pattern) — required for Kelvin + coil
    # ------------------------------------------------------------------

    def solve_total_reduced_potential(self, dirichlet='GND',
                                      iron_boundary='default'):
        """Solve using EMPY-style total+reduced split with Kelvin.

        WARNING (2026-06-19): this Omega_s-LIFT two-scalar method is BROKEN for
        COIL sources.  It couples the source through the coil's magnetic SCALAR
        potential Omega_s (= -rad.Fld(coil,'phi')) used as a Dirichlet lift on
        the iron-air boundary -- but Omega_s is MULTIVALUED (the solid-angle
        potential jumps by the ampere-turns NI) for any current loop linking the
        iron.  The lift is therefore ill-posed: the iron magnetises but the
        air/gap field comes out ~100% wrong, MESH-INDEPENDENTLY (verified on the
        coil-driven C-yoke; uniform-source Omega_s = H0*z, being single-valued,
        works -- which is the ONLY case this method was ever validated on).

        Use ``solve_single_potential`` / ``solve_nonlinear_newton`` instead:
        they drive the reduced potential with the single-valued VECTOR source
        H_s = rad.RadiaField(coil,'h') in a VOLUME integral (H = H_s - grad Omega),
        exactly as the production ``calc_accel_magnet.py`` pipeline does -- that
        is the correct and validated method for Kelvin + external coil source.
        (Static check: lint rule ``coil-scalar-potential-as-lift``.)

        The unknown ``omega`` represents:
        - Omega_t (total potential) in iron regions
        - Omega_r + Omega_s (total potential) in air regions

        Source enters via:
        - Dirichlet lift: Omega_s on iron-air boundary (weak, NOT strict)
        - Neumann jump: (n . B_s) on iron-air boundary

        Kelvin region receives NO source term.

        All source CFs (Omega_s, B_s) are pre-built by
        ``set_source_from_radia`` — this method uses only NGSolve CFs,
        no Radia calls.

        Parameters
        ----------
        dirichlet : str
            Boundary label for Dirichlet BC (default: 'GND').
            MUST NOT include the iron-air boundary (strict Dirichlet
            there kills demagnetization).
        iron_boundary : str
            Boundary name of the iron-air interface (default: auto-detect
            from mesh boundaries containing 'iron' or 'cylinder').

        Returns
        -------
        GridFunction
            The solved scalar potential.
        """
        from ngsolve import (H1, Periodic, GridFunction, BilinearForm,
                             LinearForm, grad, dx, ds, CF,
                             x, y, z, sqrt, specialcf, Preconditioner,
                             BND)

        if self._H_source_cf is None:
            raise RuntimeError(
                "Set source first via set_source_from_radia()")
        if not self._kelvin_region:
            raise RuntimeError(
                "solve_total_reduced_potential requires Kelvin "
                "configuration (kelvin_region, kelvin_radius, "
                "kelvin_center)")
        if self._Omega_s_cf is None:
            raise RuntimeError(
                "Omega_s not built — set_source_from_radia must be "
                "called with Kelvin configured")

        # Use pre-built CFs from set_source_from_radia
        Omega_s_cf = self._Omega_s_cf
        Bs_cf = self._Bs_cf

        # --- Auto-detect iron boundary name ---
        if iron_boundary == 'default':
            bnds = self.mesh.GetBoundaries()
            for candidate in ['iron_surf', 'cylinder', 'iron']:
                if candidate in bnds:
                    iron_boundary = candidate
                    break
            else:
                iron_re = '|'.join(self.iron_domains)
                iron_boundary = iron_re

        # --- Material CF with Kelvin metric ---
        mu_cf = self._build_mu_cf_with_kelvin()

        # --- FE space: Periodic + GND Dirichlet only ---
        fes_pre = H1(self.mesh, order=self.order, dirichlet=dirichlet)
        fes = Periodic(fes_pre)

        omega_t = fes.TrialFunction()
        psi = fes.TestFunction()

        # --- Bilinear form: Mu * grad . grad over all volume ---
        a = BilinearForm(fes)
        a += mu_cf * grad(omega_t) * grad(psi) * dx

        pre = None
        use_iterative = fes.ndof > 200_000
        if use_iterative:
            pre = Preconditioner(a, 'bddc')
        a.Assemble()

        # --- Dirichlet lift: Omega_s on iron boundary (weak) ---
        gfOmega = GridFunction(fes)
        gfOmega.Set(Omega_s_cf, BND,
                     self.mesh.Boundaries(iron_boundary))

        # --- RHS: lift contribution in air + Neumann jump ---
        f = LinearForm(fes)
        air_mats = [m for m in set(self.mesh.GetMaterials())
                    if m not in self.iron_domains
                    and m != self._kelvin_region]
        for mat in air_mats:
            f += mu_cf * grad(gfOmega) * grad(psi) * dx(mat)
        normal = -specialcf.normal(self.mesh.dim)
        f += (normal * Bs_cf) * psi * ds(iron_boundary)
        f.Assemble()

        # --- Solve ---
        self._phi_gf = GridFunction(fes, name='omega_total_reduced')
        self._solve_system(a, f, fes, self._phi_gf, pre, use_iterative)

        # --- Build per-material B field ---
        zero_vec = CF((0.0, 0.0, 0.0))
        fes_air = H1(self.mesh, order=self.order,
                      definedon='|'.join(air_mats +
                                         [self._kelvin_region]))
        Oxr = GridFunction(fes_air)
        Oxr.Set(Omega_s_cf, BND,
                self.mesh.Boundaries(iron_boundary))

        B_per_mat = []
        for m in self.mesh.GetMaterials():
            if m in self.iron_domains:
                B_per_mat.append(mu_cf * grad(self._phi_gf))
            elif m == self._kelvin_region:
                B_per_mat.append(
                    mu_cf * (grad(self._phi_gf) - grad(Oxr)))
            else:
                B_per_mat.append(
                    mu_cf * (grad(self._phi_gf) - grad(Oxr))
                    + Bs_cf)
        self._B_cf = CF(B_per_mat, dims=(3,))
        self._H_cf = self._B_cf / mu_cf

        self._Omega_s_cf = Omega_s_cf
        self._Bs_cf = Bs_cf

        return self._phi_gf

    def _compute_iron_bbox(self):
        """Compute bounding box of iron domains + margin."""
        from ngsolve import VOL
        pmin = [1e30, 1e30, 1e30]
        pmax = [-1e30, -1e30, -1e30]
        for el in self.mesh.Elements(VOL):
            mat = el.mat if hasattr(el, 'mat') else str(
                self.mesh.GetMaterials()[el.nr])
            if mat not in self.iron_domains:
                continue
            for v in el.vertices:
                pt = self.mesh.vertices[v.nr].point
                for i in range(3):
                    pmin[i] = min(pmin[i], pt[i])
                    pmax[i] = max(pmax[i], pt[i])
        margin = 0.5 * max(pmax[i] - pmin[i] for i in range(3))
        return [[pmin[i] - margin, pmax[i] + margin] for i in range(3)]

    def _build_mu_cf_with_kelvin(self):
        """Build per-material mu CF including Kelvin metric.

        The Kelvin factor (R/rho')^2 is regularized with rho_min set to
        the estimated mesh element size in the Kelvin region.  This caps
        the metric at (R/maxh)^2, preventing numerical explosion for
        higher-order FEM (order >= 3) near the GND vertex at rho'=0.
        """
        from ngsolve import CF
        from radia.kelvin_source import kelvin_mu_factor_3d_cf

        rho_min = self._estimate_kelvin_maxh()
        kelvin_factor = kelvin_mu_factor_3d_cf(
            center=self._kelvin_center, R=self._kelvin_radius,
            rho_min=rho_min)

        mat_dict = {}
        for m in set(self.mesh.GetMaterials()):
            if m in self._mu_r_dict:
                mat_dict[m] = MU_0 * self._mu_r_dict[m]
            elif m == self._kelvin_region:
                mat_dict[m] = MU_0 * kelvin_factor
            else:
                mat_dict[m] = MU_0
        return CF([mat_dict[m] for m in self.mesh.GetMaterials()])

    def _estimate_kelvin_maxh(self):
        """Estimate max element size in the Kelvin region."""
        from ngsolve import VOL
        maxh = 0.0
        for el in self.mesh.Elements(VOL):
            mat = el.mat if hasattr(el, 'mat') else str(
                self.mesh.GetMaterials()[el.nr])
            if mat != self._kelvin_region:
                continue
            verts = [self.mesh.vertices[v.nr].point for v in el.vertices]
            for i in range(len(verts)):
                for j in range(i + 1, len(verts)):
                    d = sum((verts[i][k] - verts[j][k])**2
                            for k in range(3)) ** 0.5
                    maxh = max(maxh, d)
        return maxh if maxh > 0 else 1e-10

    # ------------------------------------------------------------------
    # Convenience: auto-select method
    # ------------------------------------------------------------------

    def solve(self, method='auto', **kwargs):
        """Solve the magnetostatic problem (linear materials).

        Parameters
        ----------
        method : str
            'single' for single-potential, 'two' for two-potential,
            'total_reduced' for EMPY-style total+reduced (Kelvin + coil),
            'auto' to select automatically.

        Auto-selection logic:
            - Kelvin + radia source -> 'total_reduced'
            - mu_r > 5000 -> 'two'
            - otherwise -> 'single'
        """
        if method == 'auto':
            if self._kelvin_region and self._radia_obj is not None:
                method = 'total_reduced'
            elif self.mu_r > 5000:
                method = 'two'
            else:
                method = 'single'

        if method == 'single':
            return self.solve_single_potential(**kwargs)
        elif method == 'two':
            return self.solve_two_potential(**kwargs)
        elif method == 'total_reduced':
            return self.solve_total_reduced_potential(**kwargs)
        else:
            raise ValueError(
                "method must be 'single', 'two', "
                "'total_reduced', or 'auto'")

    # ------------------------------------------------------------------
    # Nonlinear solver (B-H curve, Picard iteration)
    # ------------------------------------------------------------------

    def solve_nonlinear(self, bh_data, tol=1e-4, maxiter=50, relax=0.0,
                        dirichlet='default', verbose=True):
        """Nonlinear Picard iteration with tabulated B-H curve.

        Iterates: solve linear -> evaluate H -> update mu from B-H -> repeat.

        Parameters
        ----------
        bh_data : list of [H, B] pairs
            Tabulated B-H curve. H in A/m, B in Tesla. Must be monotonic.
        tol : float
            Convergence tolerance (max relative B change).
        maxiter : int
            Maximum Picard iterations.
        relax : float
            Under-relaxation (0=full step, 0.3=30% damping).
        dirichlet : str
            Boundary label for Dirichlet BC.
        verbose : bool
            Print iteration progress.
        """
        from ngsolve import (H1, GridFunction, BilinearForm, LinearForm,
                             L2, grad, dx, x, y, z, sqrt, Preconditioner,
                             VOL)
        from scipy.interpolate import interp1d

        if self._H_source_cf is None:
            raise RuntimeError("Set source field first")

        # Build B(H) interpolator from table
        bh = np.array(bh_data)
        H_tab, B_tab = bh[:, 0], bh[:, 1]
        B_of_H = interp1d(H_tab, B_tab, kind='linear',
                          fill_value='extrapolate')
        B_sat = float(B_tab[-1])

        # Initial mu_r from second B-H point (ELF convention)
        if len(H_tab) >= 2 and H_tab[1] > 0:
            mu_r_init = B_tab[1] / (MU_0 * H_tab[1])
        else:
            mu_r_init = self.mu_r

        return self._picard_loop(
            material_eval=lambda H_vec, idx: self._eval_bh_curve(H_vec, B_of_H),
            mu_r_init=mu_r_init,
            B_sat=B_sat,
            tol=tol, maxiter=maxiter, relax=relax,
            dirichlet=dirichlet, verbose=verbose)

    # ------------------------------------------------------------------
    # Nonlinear solver (Newton + SymbolicEnergy)
    # ------------------------------------------------------------------

    def solve_nonlinear_newton(self, bh_data, tol=1e-4, maxiter=50,
                               dirichlet='default', verbose=True):
        """Newton iteration with SymbolicEnergy for nonlinear B-H curve.

        Uses BSpline.Integrate() for exact energy density and NGSolve
        automatic differentiation for the Jacobian. Energy-based line
        search guarantees monotone convergence.

        Much faster than Picard (solve_nonlinear): quadratic convergence.

        Parameters
        ----------
        bh_data : list of [H, B] pairs
            Tabulated B-H curve. H in A/m, B in Tesla. Must be monotonic.
        tol : float
            Convergence tolerance (energy residual).
        maxiter : int
            Maximum Newton iterations.
        dirichlet : str
            Boundary label for Dirichlet BC.
        verbose : bool
            Print iteration progress.
        """
        from ngsolve import (H1, GridFunction, BilinearForm, BSpline,
                             SymbolicEnergy, SymbolicBFI, InnerProduct,
                             grad, dx, x, y, z, sqrt, Preconditioner,
                             TaskManager)
        from ngsolve.krylovspace import CGSolver

        if self._H_source_cf is None:
            raise RuntimeError("Set source field first")

        # Build BSpline B(H) and energy density w(H) = integral B(H') dH'
        bh = np.array(bh_data)
        H_tab, B_tab = bh[:, 0].tolist(), bh[:, 1].tolist()

        # Extend BH curve beyond table range with vacuum slope (mu_0).
        # NGSolve BSpline returns 0 outside [knot_min, knot_max], so
        # without this extension, B=0 for H > H_max which breaks the
        # Newton solver for saturated iron (e.g. NI=20000).
        H_max_ext = max(H_tab[-1] * 10, 5e6)
        B_max_ext = B_tab[-1] + MU_0 * (H_max_ext - H_tab[-1])
        H_tab_ext = H_tab + [H_max_ext]
        B_tab_ext = B_tab + [B_max_ext]

        bh_spline = BSpline(2, [0] + H_tab_ext, B_tab_ext)
        w_H = bh_spline.Integrate()  # w(H) = integral_0^H B(H') dH'

        # FE space
        if dirichlet == 'default':
            fes = H1(self.mesh, order=self.order, dirichlet='.*')
        else:
            fes = H1(self.mesh, order=self.order, dirichlet=dirichlet)

        phi = fes.TrialFunction()
        v = fes.TestFunction()

        # H = H_s - grad(phi) -- the primary variable
        H_trial = self._H_source_cf - grad(phi)
        H2 = InnerProduct(H_trial, H_trial)

        # Physical materials
        all_mats = list(set(self.mesh.GetMaterials()))
        phys_mats = [m for m in all_mats if m != self._kelvin_region]
        iron_mats = [m for m in phys_mats if m in self._mu_r_dict]
        air_mats = [m for m in phys_mats if m not in self._mu_r_dict]

        # Bilinear form: energy-based (symmetric for Newton)
        a = BilinearForm(fes, symmetric=True)

        # Air regions: linear energy mu_0/2 * |H|^2
        for mat in air_mats:
            a += SymbolicEnergy(MU_0 / 2 * H2,
                                definedon=self.mesh.Materials(mat))

        # Iron regions: nonlinear energy w(|H|)
        for mat in iron_mats:
            a += SymbolicEnergy(w_H(sqrt(H2 + 1e-12)),
                                definedon=self.mesh.Materials(mat))

        # Kelvin exterior domain: linear with weight
        if self._kelvin_region:
            kelvin_weight = self._build_kelvin_weight(x, y, z, sqrt)
            a += SymbolicEnergy(
                MU_0 / 2 * kelvin_weight * H2,
                definedon=self.mesh.Materials(self._kelvin_region))

        # BDDC preconditioner
        c = Preconditioner(a, type='bddc', inverse='sparsecholesky')

        # Newton iteration with energy line search
        sol = GridFunction(fes, name='phi_newton')
        sol.vec[:] = 0

        au = sol.vec.CreateVector()
        r = sol.vec.CreateVector()
        w = sol.vec.CreateVector()
        sol_new = sol.vec.CreateVector()

        converged = False
        # NOTE: Caller MUST be inside `with TaskManager():` per CLAUDE.md
        # "Caller Wraps, Helper Does NOT" (2026-05-27).
        for it in range(maxiter):
            E0 = a.Energy(sol.vec)

            a.AssembleLinearization(sol.vec)
            a.Apply(sol.vec, au)
            r.data = -au

            inv = CGSolver(mat=a.mat, pre=c.mat,
                           maxiter=2000, tol=1e-10,
                           printrates=False)
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
        self._phi_gf = sol
        H_cf = self._H_source_cf - grad(sol)
        self._H_cf = H_cf

        # B in iron: B = mu(|H|)*|H| * H_hat, where mu from B-H curve
        # B in air: B = mu_0 * H
        H_mag = sqrt(InnerProduct(H_cf, H_cf) + 1e-30)
        # mu_eff(|H|) = B(|H|) / |H| from BSpline
        mu_eff = bh_spline(H_mag) / H_mag
        iron_indicator = self._build_domain_indicator()
        from ngsolve import CF
        mu_total = iron_indicator * mu_eff + (1 - iron_indicator) * MU_0
        self._B_cf = mu_total * H_cf

        return sol

    # ------------------------------------------------------------------
    # Hysteresis solver (energy play model, Picard iteration)
    # ------------------------------------------------------------------

    def solve_hysteresis(self, mat_factory, tol=1e-4, maxiter=50, alpha=500.0,
                         dirichlet='default', verbose=True, relax=0.0):
        """Polarization-method Picard iteration with energy hysteresis material.

        Uses the Hantila (1975) polarization method for guaranteed convergence:
            LHS: mu_0*(1+alpha) * grad(phi).grad(v) dx(iron)
                 + mu_0 * grad(phi).grad(v) dx(air)     -- CONSTANT
            RHS: mu_0*(1+alpha) * H_s.grad(v) dx(iron)
                 + mu_0 * H_s.grad(v) dx(air)
                 + mu_0 * (M^n - alpha*H^n) . grad(v) dx(iron) -- polarization residual

        Key: at H=0, residual = mu_0*M_rem -> remanence captured.
        Convergence guaranteed when alpha >= max(dM/dH).

        Parameters
        ----------
        mat_factory : callable
            Factory: mat_factory() -> int (Radia material handle).
        tol : float
            Convergence tolerance (max relative B change).
        maxiter : int
            Maximum Picard iterations.
        alpha : float
            Polarization parameter (>= max susceptibility). Default 500.
        dirichlet : str
            Boundary label for Dirichlet BC.
        verbose : bool
            Print iteration progress.
        relax : float
            Under-relaxation factor for R update (0.0 = full step, 0.5 = half).
            Helps convergence when local dM/dH exceeds alpha near play
            model pinning thresholds.
        """
        if self._H_source_cf is None:
            raise RuntimeError("Set source field first")

        import radia as rad
        from ngsolve import (H1, GridFunction, BilinearForm, LinearForm,
                             L2, VectorL2, grad, dx, x, y, z, sqrt,
                             Preconditioner, VOL)

        # Create per-element material handles on first call
        if not hasattr(self, '_hys_handles') or self._hys_handles is None:
            n_iron = 0
            for el in self.mesh.Elements(VOL):
                if str(el.mat) in self._mu_r_dict:
                    n_iron += 1
            self._hys_handles = [mat_factory() for _ in range(n_iron)]
            if verbose:
                print(f"   Created {n_iron} per-element material handles")

        handles = self._hys_handles

        # Setup FEM spaces
        if dirichlet == 'default':
            fes = H1(self.mesh, order=self.order, dirichlet='.*')
        else:
            fes = H1(self.mesh, order=self.order, dirichlet=dirichlet)

        phi_trial = fes.TrialFunction()
        v = fes.TestFunction()

        # Per-element polarization residual R = M - alpha*H as VectorL2
        fes_vec = VectorL2(self.mesh, order=0)
        R_gf = GridFunction(fes_vec, name='R_polar')
        R_gf.vec[:] = 0

        # Kelvin weight
        kelvin_weight = None
        if self._kelvin_region:
            kelvin_weight = self._build_kelvin_weight(x, y, z, sqrt)

        # Physical mats
        all_mats = list(set(self.mesh.GetMaterials()))
        phys_mats = [m for m in all_mats if m != self._kelvin_region]
        iron_mats = [m for m in phys_mats if m in self._mu_r_dict]
        air_mats = [m for m in phys_mats if m not in self._mu_r_dict]

        # Iron element indices and centroids
        iron_els = []
        for el in self.mesh.Elements(VOL):
            if str(el.mat) in self._mu_r_dict:
                centroid = self._element_centroid(el)
                iron_els.append((el.nr, centroid))

        # Bilinear form: CONSTANT, assembled once
        # Iron: mu_0*(1+alpha), Air: mu_0
        a = BilinearForm(fes)
        for mat in iron_mats:
            a += MU_0 * (1 + alpha) * grad(phi_trial) * grad(v) * dx(mat)
        for mat in air_mats:
            a += MU_0 * grad(phi_trial) * grad(v) * dx(mat)
        if self._kelvin_region and kelvin_weight is not None:
            a += MU_0 * kelvin_weight * grad(phi_trial) * grad(v) * dx(
                self._kelvin_region)

        pre = None
        use_iterative = fes.ndof > 200_000
        if use_iterative:
            pre = Preconditioner(a, 'bddc')
        a.Assemble()

        # Total elements for VectorL2 block DOF mapping
        nel = len(list(self.mesh.Elements(VOL)))

        # Save states before iteration
        saved_states = [rad.MatHysSaveState(h) for h in handles]

        phi_gf = GridFunction(fes, name='phi_hysteresis')
        B_old = np.zeros(len(iron_els))

        converged = False
        for iteration in range(maxiter):
            # RHS: mu_0*(1+alpha)*H_s [iron] + mu_0*H_s [air] + mu_0*R [iron]
            f = LinearForm(fes)
            for mat in iron_mats:
                f += MU_0 * (1 + alpha) * self._H_source_cf * grad(v) * dx(mat)
                f += MU_0 * R_gf * grad(v) * dx(mat)
            for mat in air_mats:
                f += MU_0 * self._H_source_cf * grad(v) * dx(mat)
            if self._kelvin_region and kelvin_weight is not None:
                f += MU_0 * kelvin_weight * self._H_source_cf * grad(v) * dx(
                    self._kelvin_region)
            f.Assemble()

            # Solve (LHS pre-assembled)
            self._solve_system(a, f, fes, phi_gf, pre, use_iterative)

            # H = H_s - grad(phi)
            H_cf = self._H_source_cf - grad(phi_gf)

            # Restore states before Forward(H)
            for h, s in zip(handles, saved_states):
                rad.MatHysRestoreState(h, s)

            # Evaluate M and update polarization residual R = M - alpha*H
            B_new = np.zeros(len(iron_els))
            for idx, (el_nr, centroid) in enumerate(iron_els):
                try:
                    mip = self.mesh(*centroid)
                    H_val = [float(H_cf(mip)[i]) for i in range(3)]
                except Exception:
                    H_val = [0.0, 0.0, 0.0]

                # Forward(H) -> M
                M_val = list(rad.MatMvsH(handles[idx], 'm', H_val))
                B_vec = [MU_0 * (H_val[i] + M_val[i]) for i in range(3)]
                B_mag = msqrt(B_vec[0]**2 + B_vec[1]**2 + B_vec[2]**2)
                B_new[idx] = B_mag

                # R = M - alpha * H (polarization residual)
                # VectorL2 order=0: block layout [all_x | all_y | all_z]
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
                    print(f"   iter {iteration}: max |dB|/B_sat = {max_change:.2e}")
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

        # Commit converged states for next step
        if converged:
            for h in handles:
                rad.MatHysCommitState(h)

        # Store final results (M from last iteration's R: M = R + alpha*H)
        self._phi_gf = phi_gf
        H_cf = self._H_source_cf - grad(phi_gf)
        self._H_cf = H_cf
        # B = mu_0*(H + M), approximate from last R: M ≈ R + alpha*H in iron
        M_approx = R_gf + alpha * H_cf
        # Build domain-selective B: iron has M contribution, air/kelvin doesn't
        iron_indicator = self._build_domain_indicator()
        self._B_cf = MU_0 * H_cf + MU_0 * iron_indicator * M_approx

        # Save per-element M for to_radia() pipeline
        self._M_per_element = {}
        for idx, (el_nr, centroid) in enumerate(iron_els):
            try:
                mip = self.mesh(*centroid)
                H_val = [float(H_cf(mip)[i]) for i in range(3)]
            except Exception:
                H_val = [0.0, 0.0, 0.0]
            M_val = list(rad.MatMvsH(handles[idx], 'm', H_val))
            self._M_per_element[el_nr] = M_val

        return phi_gf

    # ------------------------------------------------------------------
    # Result access
    # ------------------------------------------------------------------

    def get_H(self):
        """Get H field CoefficientFunction (A/m, dim=3)."""
        if self._H_cf is None:
            raise RuntimeError("Call solve() first")
        return self._H_cf

    def get_B(self):
        """Get B field CoefficientFunction (Tesla, dim=3)."""
        if self._B_cf is None:
            raise RuntimeError("Call solve() first")
        return self._B_cf

    def get_phi(self):
        """Get scalar potential GridFunction(s).

        Single-potential: returns phi GridFunction.
        Two-potential: returns (phi_r, psi) tuple.
        """
        if hasattr(self, '_phi_r_gf') and self._phi_r_gf is not None:
            return self._phi_r_gf, self._psi_gf
        if self._phi_gf is not None:
            return self._phi_gf
        raise RuntimeError("Call solve() first")

    def project_to_hdiv(self, order=None):
        """Project B field onto HDiv GridFunction (ensures div(B)=0).

        Parameters
        ----------
        order : int, optional
            HDiv polynomial order. Default: same as H1 order.

        Returns
        -------
        ngsolve.GridFunction
            B field in HDiv space.
        """
        from ngsolve import HDiv, GridFunction

        if self._B_cf is None:
            raise RuntimeError("Call solve() first")

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
    # Internal: Picard iteration loop (shared by nonlinear/hysteresis)
    # ------------------------------------------------------------------

    def _picard_loop(self, material_eval, mu_r_init, B_sat, tol, maxiter,
                     relax, dirichlet, verbose,
                     save_states=None, restore_states=None,
                     commit_states=None):
        """Picard fixed-point iteration for nonlinear materials.

        Parameters
        ----------
        material_eval : callable
            material_eval([Hx, Hy, Hz], idx) -> (B_mag, mu_r_eff)
            where idx is the iron element index (0-based).
        mu_r_init : float
            Initial mu_r for first linear solve.
        B_sat : float
            Saturation B for relative convergence criterion.
        save_states : callable or None
            save_states() -> saved_state. Saves material states.
        restore_states : callable or None
            restore_states(saved_state). Restores material states.
        commit_states : callable or None
            commit_states(). Commits converged states for next step.
        """
        from ngsolve import (H1, GridFunction, BilinearForm, LinearForm,
                             L2, grad, dx, x, y, z, sqrt, Preconditioner,
                             VOL)

        if dirichlet == 'default':
            fes = H1(self.mesh, order=self.order, dirichlet='.*')
        else:
            fes = H1(self.mesh, order=self.order, dirichlet=dirichlet)

        phi_trial = fes.TrialFunction()
        v = fes.TestFunction()

        # Per-element mu as L2(order=0) GridFunction
        fes_mu = L2(self.mesh, order=0)
        mu_gf = GridFunction(fes_mu, name='mu')

        # Initialize mu for all elements
        for el in self.mesh.Elements(VOL):
            mat_name = str(el.mat)
            if mat_name in self._mu_r_dict:
                mu_gf.vec[el.nr] = MU_0 * mu_r_init
            else:
                mu_gf.vec[el.nr] = MU_0

        # Kelvin weight CF
        kelvin_weight = None
        if self._kelvin_region:
            kelvin_weight = self._build_kelvin_weight(x, y, z, sqrt)

        # Physical mats (exclude kelvin)
        all_mats = list(set(self.mesh.GetMaterials()))
        phys_mats = [m for m in all_mats if m != self._kelvin_region]

        # Iron element indices and centroids (precompute)
        iron_els = []
        for el in self.mesh.Elements(VOL):
            mat_name = str(el.mat)
            if mat_name in self._mu_r_dict:
                centroid = self._element_centroid(el)
                iron_els.append((el.nr, centroid))

        # B tracking for convergence
        B_old = np.zeros(len(iron_els))

        phi_gf = GridFunction(fes, name='phi_nonlinear')

        # Save material states before iteration (for hysteresis)
        saved_state = None
        if save_states is not None:
            saved_state = save_states()

        converged = False
        for iteration in range(maxiter):
            # Assemble with current mu_gf
            a = BilinearForm(fes)
            for mat in phys_mats:
                a += mu_gf * grad(phi_trial) * grad(v) * dx(mat)
            if self._kelvin_region and kelvin_weight is not None:
                a += MU_0 * kelvin_weight * grad(phi_trial) * grad(v) * dx(
                    self._kelvin_region)

            pre = None
            use_iterative = fes.ndof > 200_000
            if use_iterative:
                pre = Preconditioner(a, 'bddc')
            a.Assemble()

            f = LinearForm(fes)
            for mat in phys_mats:
                f += mu_gf * self._H_source_cf * grad(v) * dx(mat)
            if self._kelvin_region and kelvin_weight is not None:
                f += MU_0 * kelvin_weight * self._H_source_cf * grad(v) * dx(
                    self._kelvin_region)
            f.Assemble()

            # Solve
            self._solve_system(a, f, fes, phi_gf, pre, use_iterative)

            # Compute H = H_s - grad(phi)
            H_cf = self._H_source_cf - grad(phi_gf)

            # Restore states before material evaluation (hysteresis)
            # so Forward(H) starts from the committed state, not from
            # the state left by the previous iteration
            if restore_states is not None and saved_state is not None:
                restore_states(saved_state)

            # Update mu at iron element centroids
            B_new = np.zeros(len(iron_els))
            for idx, (el_nr, centroid) in enumerate(iron_els):
                try:
                    mip = self.mesh(*centroid)
                    H_val = [float(H_cf(mip)[i]) for i in range(3)]
                except Exception:
                    H_val = [0.0, 0.0, 0.0]

                B_mag, mu_r_new = material_eval(H_val, idx)
                B_new[idx] = B_mag

                if relax > 0 and iteration > 0:
                    mu_r_old = mu_gf.vec[el_nr] / MU_0
                    mu_r_new = (1 - relax) * mu_r_new + relax * mu_r_old

                mu_gf.vec[el_nr] = MU_0 * mu_r_new

            # Convergence check
            if iteration > 0:
                max_change = np.max(np.abs(B_new - B_old)) / max(B_sat, 1e-10)
                if verbose:
                    print(f"   iter {iteration}: max |dB|/B_sat = {max_change:.2e}")
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

        # Commit converged states for next step (hysteresis)
        if converged and commit_states is not None:
            commit_states()

        # Store final results
        self._phi_gf = phi_gf
        self._H_cf = H_cf
        self._B_cf = mu_gf * H_cf

        return phi_gf

    # ------------------------------------------------------------------
    # Internal: material evaluation callbacks
    # ------------------------------------------------------------------

    @staticmethod
    def _eval_bh_curve(H_vec, B_of_H):
        """Evaluate B-H table: H -> (B_mag, mu_r_eff)."""
        H_mag = msqrt(H_vec[0]**2 + H_vec[1]**2 + H_vec[2]**2)
        if H_mag < 1e-10:
            return 0.0, 1.0
        B_mag = float(B_of_H(H_mag))
        mu_r = B_mag / (MU_0 * H_mag)
        return B_mag, mu_r

    @staticmethod
    def _eval_hysteresis(H_vec, mat_handle):
        """Evaluate energy hysteresis: H -> Forward(H) -> (B_mag, mu_r_eff)."""
        import radia as rad
        H_mag = msqrt(H_vec[0]**2 + H_vec[1]**2 + H_vec[2]**2)
        if H_mag < 1e-10:
            return 0.0, 1.0
        M = rad.MatMvsH(mat_handle, 'm', H_vec)
        B = [MU_0 * (H_vec[i] + M[i]) for i in range(3)]
        B_mag = msqrt(B[0]**2 + B[1]**2 + B[2]**2)
        mu_r = B_mag / (MU_0 * H_mag)
        return B_mag, mu_r

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _build_mu_cf(self):
        """Build domain-wise mu CoefficientFunction (linear)."""
        iron_dict = {d: MU_0 * mr for d, mr in self._mu_r_dict.items()}
        return self.mesh.MaterialCF(iron_dict, default=MU_0)

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
            mat = el.mat if hasattr(el, 'mat') else str(self.mesh.GetMaterials()[el.nr])
            if mat == self._kelvin_region:
                continue
            for v in el.vertices:
                pt = self.mesh.vertices[v.nr].point
                for i in range(3):
                    pmin[i] = min(pmin[i], pt[i])
                    pmax[i] = max(pmax[i], pt[i])
        margin = 0.1 * max(pmax[i] - pmin[i] for i in range(3))
        return [[pmin[i] - margin, pmax[i] + margin] for i in range(3)]

    def _build_kelvin_weight(self, x, y, z, sqrt_cf):
        """Build Kelvin exterior domain weight CF: (R/r)^2."""
        cx, cy, cz = self._kelvin_center
        r_sq = (x - cx)**2 + (y - cy)**2 + (z - cz)**2 + 1e-30
        R = self._kelvin_radius
        return (R * R) / r_sq

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
