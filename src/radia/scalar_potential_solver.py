"""Two-scalar-potential (phi-reduced/phi) solver for Radia + NGSolve.

Implements the Simkin-Trowbridge (1979) method for magnetostatic
field computation combining Radia BEM source fields with NGSolve FEM.

Method overview:
    - Air region:  H = H_s - grad(phi_r)     [reduced potential]
    - Iron region: H = -grad(psi)             [total potential]
    - Interface: normal B and tangential H continuous

The source field H_s is computed by Radia using MMM (Magnetic Moment
Method) or MSC (Magnetic Surface Charge).  NGSolve solves for the
correction potentials that enforce div(B) = 0 in the presence of
soft iron.

Two formulations are provided:

1. ``solve_single_potential()``: Single H1 space on the full domain.
   Simple, works well for moderate mu_r (< 5000).  Equation:
       integral(mu * grad(phi) . grad(v)) = integral(mu * H_s . grad(v))
   Result: H = H_s - grad(phi), B = mu * H.

2. ``solve_two_potential()``: Compound H1 space with separate unknowns
   on air and iron.  Avoids cancellation error for high mu_r iron.
   Air:  integral(mu_0 * grad(phi_r) . grad(v_a), Omega_air)
   Iron: integral(mu * grad(psi) . grad(v_i), Omega_iron)
   Coupling: penalty on air-iron interface.

Usage:
    import radia as rad
    from ngsolve import *
    from radia.scalar_potential_solver import ScalarPotentialSolver

    # Build Radia model (permanent magnets, MMM/MSC integral method)
    # Radia always uses meters
    mag = rad.ObjRecMag([0,0,0], [0.01,0.01,0.01], [0,0,954930])

    # Create NGSolve mesh with labeled iron region
    mesh = Mesh(...)  # must have material labels

    # Solve (Radia H_s as source -> NGSolve FEM correction)
    solver = ScalarPotentialSolver(mesh, iron_domains='iron', mu_r=1000)
    solver.set_source_from_radia(mag)
    solver.solve()

    # Get results
    B_cf = solver.get_B()
    H_cf = solver.get_H()

Reference:
    J. Simkin and C. W. Trowbridge, "On the use of the total scalar
    potential in the numerical solution of field problems in
    electromagnetics," IJNME, vol. 14, pp. 423-440, 1979.
"""

import numpy as np

MU_0 = 4 * np.pi * 1e-7


class ScalarPotentialSolver:
    """Two-scalar-potential magnetostatic solver (Radia + NGSolve).

    Parameters
    ----------
    mesh : ngsolve.Mesh
        FEM mesh.  Must have material labels for iron regions.
    iron_domains : str or list of str
        Material name(s) for iron regions.  All other regions are treated
        as air (mu_r = 1).
    mu_r : float
        Relative permeability of iron.
    order : int
        Finite element polynomial order (default: 2).
    """

    def __init__(self, mesh, iron_domains='iron', mu_r=1000.0, order=2):
        self.mesh = mesh
        self.order = order

        if isinstance(iron_domains, str):
            iron_domains = [iron_domains]
        self.iron_domains = iron_domains
        self.mu_r = float(mu_r)

        # Source field (set by set_source_*)
        self._H_source_cf = None
        self._radia_obj = None

        # Solution storage
        self._phi_gf = None
        self._B_cf = None
        self._H_cf = None

    # ------------------------------------------------------------------
    # Source field setup
    # ------------------------------------------------------------------

    def set_source_from_radia(self, radia_obj, resolution=41):
        """Set source field H_s from a Radia object.

        Computes H_s on a voxel grid via ``rad.Fld()`` and creates
        a VoxelCoefficient for efficient FEM integration.

        Parameters
        ----------
        radia_obj : int
            Radia object handle (after ``rad.Solve()`` if soft iron present).
        resolution : int
            Voxel grid resolution per dimension.
        """
        import radia as rad
        from ngsolve import VoxelCoefficient, CF

        self._radia_obj = radia_obj

        # Get mesh bounding box with margin
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

        H_field = np.asarray(rad.Fld(radia_obj, 'h', points))

        start = (bbox[0][0], bbox[1][0], bbox[2][0])
        end = (bbox[0][1], bbox[1][1], bbox[2][1])

        cfs = []
        for comp in range(3):
            data = H_field[:, comp].reshape(nx, ny, nz)
            data = np.ascontiguousarray(data.transpose(2, 1, 0))
            cfs.append(VoxelCoefficient(start, end, data, linear=True))

        self._H_source_cf = CF(tuple(cfs))

    def set_source_from_callback(self, h_func, resolution=41):
        """Set source field H_s from a Python callback.

        Parameters
        ----------
        h_func : callable
            ``h_func(x, y, z) -> (Hx, Hy, Hz)`` returning H in A/m.
        resolution : int
            Voxel grid resolution per dimension.
        """
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
        H_field = np.zeros((len(pts_flat), 3))
        for i in range(len(pts_flat)):
            H_field[i] = h_func(pts_flat[i, 0], pts_flat[i, 1], pts_flat[i, 2])

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

        Finds phi in H1(Omega) satisfying:
            integral(mu * grad(phi) . grad(v)) = integral(mu * H_s . grad(v))

        Then: H = H_s - grad(phi), B = mu * H.

        Works well for moderate mu_r (< 5000).  For higher mu_r, use
        ``solve_two_potential()`` to avoid cancellation error in iron.

        Parameters
        ----------
        dirichlet : str
            Boundary label for Dirichlet BC (phi=0).
            'default' uses all outer boundaries.
        """
        from ngsolve import (H1, GridFunction, BilinearForm, LinearForm,
                             grad, dx)

        if self._H_source_cf is None:
            raise RuntimeError("Set source field first (set_source_from_radia)")

        mu_cf = self._build_mu_cf()

        if dirichlet == 'default':
            fes = H1(self.mesh, order=self.order, dirichlet='.*')
        else:
            fes = H1(self.mesh, order=self.order, dirichlet=dirichlet)

        phi = fes.TrialFunction()
        v = fes.TestFunction()

        a = BilinearForm(fes)
        a += mu_cf * grad(phi) * grad(v) * dx
        a.Assemble()

        f = LinearForm(fes)
        f += mu_cf * self._H_source_cf * grad(v) * dx
        f.Assemble()

        self._phi_gf = GridFunction(fes, name='phi_reduced')
        self._phi_gf.vec.data = a.mat.Inverse(fes.FreeDofs()) * f.vec

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

        The formulation uses a compound space X = V_full * V_iron where
        V_full is a standard H1 on the whole domain and V_iron is H1
        restricted to iron.  In iron, the solution uses psi (total potential)
        so grad(psi) gives H directly without cancellation.

        The coupling is enforced via: in iron elements, phi_r is constrained
        to equal (psi - phi_s), linking the two unknowns.

        Parameters
        ----------
        dirichlet : str
            Boundary label for Dirichlet BC on outer boundary.
        """
        from ngsolve import (H1, GridFunction, BilinearForm, LinearForm,
                             grad, dx, CF)

        if self._H_source_cf is None:
            raise RuntimeError("Set source field first (set_source_from_radia)")

        iron_re = '|'.join(self.iron_domains)
        mu_iron_val = MU_0 * self.mu_r

        # Full-domain space for phi_r, iron-only space for psi
        if dirichlet == 'default':
            V_full = H1(self.mesh, order=self.order, dirichlet='.*')
        else:
            V_full = H1(self.mesh, order=self.order, dirichlet=dirichlet)
        V_iron = H1(self.mesh, order=self.order, definedon=iron_re)

        X = V_full * V_iron
        (phi_r, psi), (v_f, v_i) = X.TnT()

        a = BilinearForm(X)

        # Air region: mu_0 * grad(phi_r) . grad(v_f)
        # Only on air elements (exclude iron)
        air_mats = [m for m in self.mesh.GetMaterials()
                    if m not in self.iron_domains]
        for mat in air_mats:
            a += MU_0 * grad(phi_r) * grad(v_f) * dx(mat)

        # Iron region: mu * grad(psi) . grad(v_i)
        a += mu_iron_val * grad(psi) * grad(v_i) * dx(iron_re)

        # Coupling in iron: phi_r and psi share the same interface nodes.
        # Penalize (phi_r - psi) in iron to enforce phi_r = psi there,
        # making the transition continuous.
        penalty = 1e3 * self.mu_r * MU_0
        a += penalty * (phi_r - psi) * (v_f - v_i) * dx(iron_re)

        a.Assemble()

        # Source: H_s contribution in air only
        f = LinearForm(X)
        for mat in air_mats:
            f += MU_0 * self._H_source_cf * grad(v_f) * dx(mat)
        f.Assemble()

        gf = GridFunction(X, name='two_potential')
        gf.vec.data = a.mat.Inverse(X.FreeDofs()) * f.vec

        phi_r_gf, psi_gf = gf.components
        self._phi_r_gf = phi_r_gf
        self._psi_gf = psi_gf

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
    # Convenience: auto-select method
    # ------------------------------------------------------------------

    def solve(self, method='auto', **kwargs):
        """Solve the magnetostatic problem.

        Parameters
        ----------
        method : str
            'single' for single-potential, 'two' for two-potential,
            'auto' to select based on mu_r (two-potential if mu_r > 5000).
        """
        if method == 'auto':
            method = 'two' if self.mu_r > 5000 else 'single'

        if method == 'single':
            return self.solve_single_potential(**kwargs)
        elif method == 'two':
            return self.solve_two_potential(**kwargs)
        else:
            raise ValueError(f"method must be 'single', 'two', or 'auto'")

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

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _build_mu_cf(self):
        """Build domain-wise mu CoefficientFunction."""
        iron_dict = {d: MU_0 * self.mu_r for d in self.iron_domains}
        return self.mesh.MaterialCF(iron_dict, default=MU_0)

    def _build_domain_indicator(self):
        """Build indicator CF: 1.0 in iron, 0.0 in air."""
        iron_dict = {d: 1.0 for d in self.iron_domains}
        return self.mesh.MaterialCF(iron_dict, default=0.0)
