"""Total scalar potential solver with cohomology cuts (GMSH-FREE).

Implements the total scalar potential formulation for magnetostatics:

    H = -grad(phi) + sum_k I_k h_k

where phi is the scalar potential (H1), h_k are cohomology basis functions (HCurl, curl-free) with unit
circulation around the k-th coil loop, and I_k = N_k * I_coil are the ampere-turns for each coil.

The cohomology basis functions h_k are computed by the pure-Python ``radia.cohomology`` engine (no Gmsh, no
.msh -> .vol transfer): given an NGSolve/Netgen mesh of the multiply-connected domain (air + iron, with each
coil region left as a through-hole), it returns one curl-free, unit-circulation h_k per independent loop.
Each h_k satisfies:

    curl(h_k) = 0                  (irrotational)
    oint_{loop j} h_k . dl = delta_jk   (unit circulation around loop j)

This eliminates the need for:
    - Biot-Savart source field computation
    - VoxelCoefficient interpolation (no Gibbs artifacts at coil surfaces)
    - Explicit cut surface definition (automatic via cohomology)
    - Gmsh (homology solver + mesh transfer) -- the engine is now radia.cohomology

Weak form (I_k known):
    integral(mu * grad(phi) . grad(v)) = sum_k I_k * integral(mu * h_k . grad(v))

This is a standard H1 Poisson problem with a topological source term.

Usage:
    import ngsolve as ng
    from radia.cohomology_cut import CohomologyCutSolver

    mesh = ng.Mesh(...)               # air + iron, each coil region a through-hole (multiply-connected)
    solver = CohomologyCutSolver()
    n = solver.setup_from_mesh(mesh)  # n = number of independent coil loops (= b1)
    solver.solve([100.0], {'iron': 1000.0})  # NI=100, mu_r=1000
    B_cf = solver.get_B()

For MULTIPLE coils, pass ``loop_circles`` -- one (cx, cy, rho[, z]) test circle threading each coil -- so the
basis is tied to the physical coils (oint_{circle j} h_k = delta_jk); then NI_list[j] is the j-th coil.

Reference:
    Pellikka et al., "Homology and Cohomology Computation in Finite Element Modeling," SIAM J. Sci. Comput.,
    2013.  (The radia.cohomology engine realises the same H^1 cut basis via the combinatorial Hodge Laplacian.)
"""

import numpy as np

MU_0 = 4 * np.pi * 1e-7


class CohomologyCutSolver:
    """Total scalar potential solver with automatic cohomology cuts (gmsh-free).

    Workflow:
        1. Build an NGSolve mesh of the domain (air + iron) with each coil region left as a through-hole
           (so the domain is multiply-connected; b1 = number of coils).
        2. Call setup_from_mesh() to compute the cohomology cut basis (via radia.cohomology, no gmsh).
        3. Call solve() with ampere-turns and material properties.
        4. Access results via get_H(), get_B(), get_phi().
    """

    def __init__(self):
        self._mesh = None
        self._h_basis = []
        self._phi_gf = None
        self._H_cf = None
        self._B_cf = None
        self._mu_cf = None
        self._loops = None

    def setup_from_mesh(self, mesh, loop_circles=None):
        """Compute the cohomology cut basis on an NGSolve mesh -- pure Python, no gmsh.

        Parameters
        ----------
        mesh : ngsolve.Mesh
            Mesh of the multiply-connected domain (air + iron, each coil region a through-hole).
        loop_circles : list of (cx, cy, rho[, z]) or None
            One test circle threading each coil.  When given (length must equal b1), the basis is recombined so
            that oint_{circle j} h_k = delta_jk -- i.e. h_j is tied to the j-th physical coil and NI_list[j] is
            its ampere-turns.  When None, the b1 cut functions come in the engine's natural (topological) order
            with unit circulation around an internal homology-generator basis (fine for a single coil).

        Returns
        -------
        int
            Number of cohomology generators (= number of independent coil loops, b1).
        """
        from ngsolve import GridFunction
        from radia.cohomology import cohomology_basis, circulation

        basis, b1, fes, _ctx, loops = cohomology_basis(mesh)
        self._mesh = mesh
        self._loops = loops

        if loop_circles is not None and b1 > 0:
            if len(loop_circles) != b1:
                raise ValueError(
                    f"{len(loop_circles)} loop_circles supplied but b1={b1} cohomology generators were found")
            # C[k][j] = circulation of basis[k] around test circle j; recombine h_j = sum_k a[k,j] basis[k]
            # with a = inv(C^T) so that oint_{circle i} h_j = delta_ij (tie each cut to its physical coil).
            C = np.array([[circulation(basis[k], mesh, *loop_circles[j]) for j in range(b1)]
                          for k in range(b1)])
            a = np.linalg.inv(C.T)
            tied = []
            for j in range(b1):
                gf = GridFunction(fes)
                acc = np.zeros(fes.ndof)
                for k in range(b1):
                    acc += a[k, j] * basis[k].vec.FV().NumPy()
                gf.vec.FV().NumPy()[:] = acc
                tied.append(gf)
            basis = tied

        self._h_basis = basis
        return b1

    def solve(self, NI_list, mu_r_dict=None, order=2,
              dirichlet='outer', kelvin_region=None, kelvin_radius=None,
              kelvin_center=None):
        """Solve the total scalar potential problem.

        Parameters
        ----------
        NI_list : list of float
            Ampere-turns for each coil [NI_1, NI_2, ...].
            Length must match number of cohomology generators.
        mu_r_dict : dict or None
            Material name -> relative permeability, e.g. {'iron': 1000}.
            Regions not listed default to mu_r=1 (air).
        order : int
            H1 polynomial order (default: 2).
        dirichlet : str
            Boundary label for Dirichlet BC (phi=0).
        kelvin_region : str or None
            Material name for the Kelvin transform shell region.
            If set, applies weight (R/r)^2 to the bilinear form there.
        kelvin_radius : float or None
            Kelvin sphere radius R (in meters). Required if kelvin_region is set.
        kelvin_center : tuple or None
            Center of the Kelvin sphere (x0, y0, z0). Default: origin.
        """
        from ngsolve import (H1, GridFunction, BilinearForm, LinearForm,
                             grad, dx, x, y, z, sqrt)

        if self._mesh is None:
            raise RuntimeError("Call setup_from_mesh() first")

        if len(NI_list) != len(self._h_basis):
            raise ValueError(
                f"NI_list has {len(NI_list)} entries but "
                f"{len(self._h_basis)} cohomology generators exist")

        if kelvin_region and not kelvin_radius:
            raise ValueError("kelvin_radius required when kelvin_region is set")

        # Build mu CoefficientFunction
        self._mu_cf = self._build_mu_cf(mu_r_dict)

        # FEM spaces
        fes = H1(self._mesh, order=order, dirichlet=dirichlet)
        phi = fes.TrialFunction()
        v = fes.TestFunction()

        # Choose solver strategy based on problem size
        ndofs = fes.ndof
        use_iterative = ndofs > 200000

        # Bilinear form: a(phi, v) = integral(mu * grad(phi) . grad(v))
        a = BilinearForm(fes)
        if kelvin_region:
            # Physical region: standard bilinear form
            all_mats = list(set(self._mesh.GetMaterials()))
            phys_mats = [m for m in all_mats if m != kelvin_region]
            for mat in phys_mats:
                a += self._mu_cf * grad(phi) * grad(v) * dx(mat)
            # Kelvin exterior domain: weight = (R/r)^2 (3D Kelvin transform)
            cx, cy, cz = kelvin_center if kelvin_center else (0, 0, 0)
            r_cf = sqrt((x - cx)**2 + (y - cy)**2 + (z - cz)**2)
            R = kelvin_radius
            kelvin_weight = (R * R) / (r_cf * r_cf + 1e-30)
            a += MU_0 * kelvin_weight * grad(phi) * grad(v) * dx(kelvin_region)
        else:
            a += self._mu_cf * grad(phi) * grad(v) * dx

        # BDDC preconditioner must be registered before Assemble
        pre = None
        if use_iterative:
            from ngsolve import Preconditioner
            pre = Preconditioner(a, 'bddc')
        a.Assemble()

        # Source: f(v) = sum_k NI_k * integral(mu * h_k . grad(v))
        # Source only in physical region (not in Kelvin exterior domain)
        f = LinearForm(fes)
        if kelvin_region:
            for k, NI_k in enumerate(NI_list):
                if abs(NI_k) > 0:
                    for mat in phys_mats:
                        f += NI_k * self._mu_cf * self._h_basis[k] * grad(v) * dx(mat)
        else:
            for k, NI_k in enumerate(NI_list):
                if abs(NI_k) > 0:
                    f += NI_k * self._mu_cf * self._h_basis[k] * grad(v) * dx
        f.Assemble()

        # Solve
        self._phi_gf = GridFunction(fes, name='phi_total')
        if use_iterative:
            from ngsolve.krylovspace import CGSolver
            print(f"   Using BDDC+CG ({ndofs} DOFs)...", flush=True)
            inv = CGSolver(mat=a.mat, pre=pre.mat, maxiter=2000,
                           printrates=False, tol=1e-10)
            self._phi_gf.vec.data = inv * f.vec
            print(f"   CG converged ({ndofs} DOFs)", flush=True)
        else:
            self._phi_gf.vec.data = a.mat.Inverse(fes.FreeDofs()) * f.vec

        # Build result fields: H = -grad(phi) + sum NI_k * h_k
        from ngsolve import CF
        H_source = CF((0, 0, 0))
        for k, NI_k in enumerate(NI_list):
            if abs(NI_k) > 0:
                H_source = H_source + NI_k * self._h_basis[k]
        self._H_cf = -grad(self._phi_gf) + H_source
        self._B_cf = self._mu_cf * self._H_cf

        return self._phi_gf

    def solve_nonlinear(self, NI_list, bh_data, order=2,
                        dirichlet='outer', iron_domain='iron',
                        tol=1e-4, maxiter=50, relax=1.0):
        """Solve with nonlinear B-H curve (Newton-Raphson).

        Parameters
        ----------
        NI_list : list of float
            Ampere-turns for each coil.
        bh_data : list of [H, B] pairs
            B-H curve data [[H1, B1], [H2, B2], ...].
        order : int
            H1 polynomial order.
        dirichlet : str
            Boundary label for Dirichlet BC.
        iron_domain : str
            Material name for iron region.
        tol : float
            Convergence tolerance.
        maxiter : int
            Maximum iterations.
        relax : float
            Under-relaxation factor (0 < relax <= 1).
        """
        from ngsolve import (H1, GridFunction, BilinearForm, LinearForm,
                             grad, dx, Integrate, InnerProduct, CF)

        if self._mesh is None:
            raise RuntimeError("Call setup_from_mesh() first")

        bh = np.array(bh_data)
        H_data, B_data = bh[:, 0], bh[:, 1]

        fes = H1(self._mesh, order=order, dirichlet=dirichlet)

        self._phi_gf = GridFunction(fes, name='phi_total')
        phi_new = GridFunction(fes)

        mu_r_current = 1000.0  # initial guess

        for iteration in range(maxiter):
            mu_r_dict = {iron_domain: mu_r_current}
            self._mu_cf = self._build_mu_cf(mu_r_dict)

            phi = fes.TrialFunction()
            v = fes.TestFunction()

            a = BilinearForm(fes)
            a += self._mu_cf * grad(phi) * grad(v) * dx
            a.Assemble()

            f = LinearForm(fes)
            for k, NI_k in enumerate(NI_list):
                if abs(NI_k) > 0:
                    f += NI_k * self._mu_cf * self._h_basis[k] * grad(v) * dx
            f.Assemble()

            phi_new.vec.data = a.mat.Inverse(fes.FreeDofs()) * f.vec

            # Under-relaxation
            diff_vec = phi_new.vec.CreateVector()
            diff_vec.data = phi_new.vec - self._phi_gf.vec
            self._phi_gf.vec.data += relax * diff_vec

            # Update H and compute mu_r from B-H curve
            H_source = CF((0, 0, 0))
            for k, NI_k in enumerate(NI_list):
                if abs(NI_k) > 0:
                    H_source = H_source + NI_k * self._h_basis[k]
            H_cf = -grad(self._phi_gf) + H_source
            H_mag_sq = Integrate(InnerProduct(H_cf, H_cf) * dx(iron_domain),
                                 self._mesh)
            vol_iron = Integrate(CF(1) * dx(iron_domain), self._mesh)
            H_avg = np.sqrt(H_mag_sq / max(vol_iron, 1e-30))

            B_interp = np.interp(H_avg, H_data, B_data)
            mu_r_new = B_interp / (MU_0 * H_avg) if H_avg > 1 else 1000.0
            mu_r_new = max(1.0, mu_r_new)

            change = abs(mu_r_new - mu_r_current) / max(mu_r_current, 1.0)
            mu_r_current = mu_r_new

            if change < tol:
                break

        self._mu_cf = self._build_mu_cf({iron_domain: mu_r_current})
        self._H_cf = -grad(self._phi_gf) + H_source
        self._B_cf = self._mu_cf * self._H_cf

        return self._phi_gf

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
        """Get scalar potential GridFunction."""
        if self._phi_gf is None:
            raise RuntimeError("Call solve() first")
        return self._phi_gf

    def get_mesh(self):
        """Get the NGSolve mesh."""
        return self._mesh

    def get_cohomology_basis(self):
        """Get list of HCurl GridFunctions (cohomology basis)."""
        return self._h_basis

    def project_to_hdiv(self, order=None):
        """Project B field onto HDiv GridFunction (ensures div(B)=0).

        Parameters
        ----------
        order : int, optional
            HDiv polynomial order. Default: 2.

        Returns
        -------
        ngsolve.GridFunction
            B field in HDiv space.
        """
        from ngsolve import HDiv, GridFunction

        if self._B_cf is None:
            raise RuntimeError("Call solve() first")
        if order is None:
            order = 2
        fes_hdiv = HDiv(self._mesh, order=order)
        B_gf = GridFunction(fes_hdiv, name='B')
        B_gf.Set(self._B_cf)
        return B_gf

    # ------------------------------------------------------------------
    # Internal methods
    # ------------------------------------------------------------------

    def _build_mu_cf(self, mu_r_dict):
        """Build domain-wise mu CoefficientFunction."""
        if mu_r_dict is None:
            from ngsolve import CF
            return CF(MU_0)

        iron_dict = {name: MU_0 * mu_r for name, mu_r in mu_r_dict.items()}
        return self._mesh.MaterialCF(iron_dict, default=MU_0)
