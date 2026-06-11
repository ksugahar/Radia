"""Clebsch potential representation for magnetostatics.

B = grad(psi) x grad(chi)   automatically satisfies div(B) = 0.

Two solvers:

  AxisymStreamFunctionSolver
      Axisymmetric magnetostatics.  Uses Psi = r * A_phi (the standard
      Stokes flux function).  Weak form::

          int (1/(mu * r)) grad(Psi) . grad(v) dr dz
              = int J_phi * v dr dz

      B_r = -(1/r) dPsi/dz,   B_z = (1/r) dPsi/dr.
      Standard NGSolve H1 + 1/r coefficient.  Simpler and more accurate
      at the axis than the raw A_phi formulation because Psi -> 0 at r=0
      (natural Dirichlet) so the 1/r singularity is suppressed.

  ClebschSolver  (3-D, novel)
      Forward magnetostatic problem: given J, find (psi, chi) in H^1 x H^1
      such that B = grad(psi) x grad(chi) satisfies curl(B/mu) = J.
      Alternating iteration::

          step psi  : fix chi, solve anisotropic-diffusion BVP for psi
          step chi  : fix psi, solve anisotropic-diffusion BVP for chi

      The natural "vector potential" is A = psi grad(chi)  =>  curl A = B.
      No gauge fixing needed for B (div B = 0 by construction).  Small
      gauge-fixing mass term eps_gauge stabilises the null-space of each
      degenerate-elliptic sub-problem.

      Reference: EMPY project (W:\\00_CAE\\NGSolve\\EMPY\\EMPY_Analysis\\
      Clebsch_potential).  3-D NGSolve realisation first reported here;
      WA-O1 CEFC 2026 (Ren et al.) "not yet 3D".

Both solvers are TaskManager-caller-wrapped by the *caller* per the
lab TaskManager policy.
"""
from __future__ import annotations

import math

import numpy as np

__all__ = [
    "AxisymStreamFunctionSolver",
    "ClebschSolver",
]

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _cross_cf(a, b):
    """Cross product of two 3-component NGSolve CoefficientFunctions."""
    from ngsolve import CF
    return CF((a[1]*b[2] - a[2]*b[1],
               a[2]*b[0] - a[0]*b[2],
               a[0]*b[1] - a[1]*b[0]))


# ---------------------------------------------------------------------------
# Axisymmetric stream function  Psi = r * A_phi
# ---------------------------------------------------------------------------

class AxisymStreamFunctionSolver:
    """Axisymmetric magnetostatic solver via stream function Psi = r * A_phi.

    The mesh lives in the r-z half-plane (x = r >= 0, y = z).
    Governing weak form::

        int (1/(mu * r)) grad(Psi) . grad(v) dr dz = int J_phi * v dr dz

    Dirichlet BC ``Psi = 0`` on the axis (r = 0) is imposed by specifying
    the boundary name via *axis_bnd*.  An outer-boundary homogeneous
    Dirichlet can be added via *outer_bnd*.

    Parameters
    ----------
    mesh : ngsolve.Mesh
        2-D mesh in the r-z half-plane (r >= 0).
    J_phi_cf : CoefficientFunction
        Azimuthal current density [A/m^2].  Zero outside conductors.
    mu_cf : CoefficientFunction, optional
        Permeability [H/m].  Defaults to mu_0.
    order : int
        H^1 polynomial order.  Default 2.
    axis_bnd : str
        Boundary name for the axis r = 0.  Defaults to ``"axis"``.
    outer_bnd : str
        Additional Dirichlet boundary (e.g. outer truncation).
        Defaults to ``"outer"``.  Pass ``None`` to skip.
    """

    def __init__(self, mesh, J_phi_cf, mu_cf=None, order=2,
                 axis_bnd="axis", outer_bnd="outer"):
        from ngsolve import H1, GridFunction, CF, x as r_coord
        self.mesh = mesh
        self.J_phi_cf = J_phi_cf
        self.mu_cf = mu_cf if mu_cf is not None else CF(4e-7 * math.pi)
        self._r = r_coord  # x = r in axisymmetric r-z mesh

        bnd_parts = [axis_bnd]
        if outer_bnd:
            bnd_parts.append(outer_bnd)
        dirichlet_str = "|".join(b for b in bnd_parts if b)

        self.fes = H1(mesh, order=order, dirichlet=dirichlet_str)
        self.gf_psi = GridFunction(self.fes, name="Psi")
        self._solved = False

    def assemble_and_solve(self, inverse="sparsecholesky"):
        """Assemble and solve in one call.

        Call this inside a ``with TaskManager():`` block in the calling script.

        Parameters
        ----------
        inverse : str
            NGSolve direct solver.  Default ``"sparsecholesky"``.
        """
        from ngsolve import BilinearForm, LinearForm, grad, InnerProduct, dx

        fes = self.fes
        u, v = fes.TnT()
        mu = self.mu_cf
        r = self._r

        a = BilinearForm(fes, symmetric=True)
        a += (1 / (mu * r)) * InnerProduct(grad(u), grad(v)) * dx

        f = LinearForm(fes)
        f += self.J_phi_cf * v * dx

        a.Assemble()
        f.Assemble()
        self.gf_psi.vec.data = a.mat.Inverse(fes.FreeDofs(),
                                              inverse=inverse) * f.vec
        self._solved = True

    @property
    def B_r_cf(self):
        """B_r = -(1/r) dPsi/dz as a CoefficientFunction."""
        if not self._solved:
            raise RuntimeError("Call assemble_and_solve() first.")
        from ngsolve import grad
        r = self._r
        return -(1 / r) * grad(self.gf_psi)[1]  # [1] = d/dz component

    @property
    def B_z_cf(self):
        """B_z = (1/r) dPsi/dr as a CoefficientFunction."""
        if not self._solved:
            raise RuntimeError("Call assemble_and_solve() first.")
        from ngsolve import grad
        r = self._r
        return (1 / r) * grad(self.gf_psi)[0]   # [0] = d/dr component

    @property
    def B_cf(self):
        """(B_r, B_z) as a 2-component CoefficientFunction."""
        from ngsolve import CF
        return CF((self.B_r_cf, self.B_z_cf))

    @property
    def A_phi_cf(self):
        """A_phi = Psi / r as a CoefficientFunction."""
        if not self._solved:
            raise RuntimeError("Call assemble_and_solve() first.")
        return self.gf_psi / self._r


# ---------------------------------------------------------------------------
# 3-D Clebsch alternating solver
# ---------------------------------------------------------------------------

class ClebschSolver:
    """3-D Clebsch potential solver: B = grad(psi) x grad(chi).

    Forward problem: given *J_cf* (current density), find (*psi*, *chi*) in
    H^1 x H^1 such that ``curl(B / mu) = J`` where ``B = grad(psi) x grad(chi)``.

    The alternating sub-problems are degenerate-elliptic (anisotropic
    diffusion with projection kernel):

    step psi  (fix chi)::

        a_psi(u,v) = (1/mu) [|grad chi|^2 (grad u . grad v)
                               - (grad u . grad chi)(grad v . grad chi)]
                   + eps * u * v
        l_psi(v)   = (J . grad chi) * v

    step chi  (fix psi, IBP of int psi J . grad v dV)::

        a_chi(u,v) = (1/mu) [|grad psi|^2 (grad u . grad v)
                               - (grad u . grad psi)(grad v . grad psi)]
                   + eps * u * v
        l_chi(v)   = -(J . grad psi) * v

    Both systems are symmetric positive semi-definite; the mass term *eps*
    lifts the constant null space.

    Initialise *chi* with the azimuthal angle ``atan2(y, x)`` for z-axis
    coil problems (most common).  For non-axisymmetric sources, supply a
    smooth initial *chi* via :meth:`initialize`.

    Parameters
    ----------
    mesh : ngsolve.Mesh
        3-D volumetric mesh.
    J_cf : CoefficientFunction, shape (3,)
        Current density [A/m^2].  Zero outside conductors.
    mu_cf : CoefficientFunction, optional
        Permeability [H/m].  Defaults to mu_0.
    order : int
        H^1 polynomial order for both psi and chi.  Default 2.
    dirichlet_bnd : str
        Boundary name for outer Dirichlet BC (``psi = 0``, ``chi`` from init).
        Default ``"outer"``.
    eps_gauge : float
        Gauge-fixing mass coefficient.  Stabilises the constant null space
        of each degenerate-elliptic sub-problem.  Default 1e-8.
    """

    def __init__(self, mesh, J_cf, mu_cf=None, order=2,
                 dirichlet_bnd="outer", eps_gauge=1e-8):
        from ngsolve import H1, GridFunction, CF
        self.mesh = mesh
        self.J_cf = J_cf
        self.mu_cf = mu_cf if mu_cf is not None else CF(4e-7 * math.pi)
        self.eps_gauge = float(eps_gauge)
        self.dirichlet_bnd = dirichlet_bnd

        self.fes = H1(mesh, order=order, dirichlet=dirichlet_bnd)
        self.gf_psi = GridFunction(self.fes, name="psi")
        self.gf_chi = GridFunction(self.fes, name="chi")
        self._niter_done = 0

    def initialize(self, chi_init=None, psi_init=None):
        """Set initial potentials.

        Parameters
        ----------
        chi_init : CoefficientFunction, optional
            Initial chi.  Default: ``atan2(y, x)`` — azimuthal angle.
            Good for z-axis coil problems.  The Dirichlet BC on *outer*
            boundary is also set from this expression via ``gf_chi.Set``.
        psi_init : CoefficientFunction, optional
            Initial psi.  Default: 0.

        Call this inside a ``with TaskManager():`` block.
        """
        from ngsolve import atan2, x, y, CF, TaskManager
        if chi_init is None:
            chi_init = atan2(y, x)
        self.gf_chi.Set(chi_init)
        if psi_init is not None:
            self.gf_psi.Set(psi_init)
        else:
            self.gf_psi.vec[:] = 0.0

    def _step_psi(self, inverse="sparsecholesky"):
        """Fix chi, solve for psi."""
        from ngsolve import BilinearForm, LinearForm, grad, InnerProduct, dx
        fes = self.fes
        u, v = fes.TnT()
        mu = self.mu_cf
        gc = grad(self.gf_chi)
        gc_sq = InnerProduct(gc, gc)

        a = BilinearForm(fes, symmetric=True)
        a += (1 / mu) * (gc_sq * InnerProduct(grad(u), grad(v))
                         - InnerProduct(grad(u), gc) * InnerProduct(grad(v), gc)) * dx
        a += self.eps_gauge * u * v * dx

        f = LinearForm(fes)
        f += InnerProduct(self.J_cf, gc) * v * dx

        a.Assemble()
        f.Assemble()
        self.gf_psi.vec.data = a.mat.Inverse(fes.FreeDofs(),
                                              inverse=inverse) * f.vec

    def _step_chi(self, inverse="sparsecholesky"):
        """Fix psi, solve for chi."""
        from ngsolve import BilinearForm, LinearForm, grad, InnerProduct, dx
        fes = self.fes
        u, v = fes.TnT()
        mu = self.mu_cf
        gp = grad(self.gf_psi)
        gp_sq = InnerProduct(gp, gp)

        a = BilinearForm(fes, symmetric=True)
        a += (1 / mu) * (gp_sq * InnerProduct(grad(u), grad(v))
                         - InnerProduct(grad(u), gp) * InnerProduct(grad(v), gp)) * dx
        a += self.eps_gauge * u * v * dx

        f = LinearForm(fes)
        # IBP of int psi J . grad v dV = -int (J . grad psi) v dV + bdry
        f += -InnerProduct(self.J_cf, gp) * v * dx

        a.Assemble()
        f.Assemble()
        self.gf_chi.vec.data = a.mat.Inverse(fes.FreeDofs(),
                                              inverse=inverse) * f.vec

    def step(self, inverse="sparsecholesky"):
        """Perform one psi-then-chi alternating step.

        Call this inside a ``with TaskManager():`` block.

        Returns
        -------
        float
            Relative change in psi: ``||psi_new - psi_old|| / ||psi_new||``.
        """
        psi_old = np.array(self.gf_psi.vec)
        self._step_psi(inverse=inverse)
        self._step_chi(inverse=inverse)
        self._niter_done += 1

        psi_new = np.array(self.gf_psi.vec)
        denom = np.linalg.norm(psi_new)
        return float(np.linalg.norm(psi_new - psi_old) / (denom + 1e-14))

    def solve_axisymmetric(self, inverse="sparsecholesky"):
        """One-shot Clebsch solve for z-axis-symmetric sources.

        For an axisymmetric source J = J_phi(r,z) * e_phi, the azimuthal
        angle chi = atan2(y, x) is the *exact* Clebsch coordinate — no
        iteration over chi is needed.  Only one psi-step is performed::

            chi  = atan2(y, x)   [exact, fixed]
            psi  = solution of the anisotropic BVP  (= r * A_phi)
            B    = grad(psi) x grad(chi)             (exact)

        This IS the Stokes stream-function solve in 3-D guise:
        the psi-step operator reduces to  int (1/(mu*r)) grad(u).grad(v) dV
        which matches the 2-D axisymmetric stream-function equation (times 2*pi).

        Call this inside a ``with TaskManager():`` block.

        Parameters
        ----------
        inverse : str
            NGSolve direct solver.  Default ``"sparsecholesky"``.
        """
        from ngsolve import atan2, x, y, CF

        # Use the projected chi for the psi-solve (proven numerically stable).
        # Near the branch cut (x<0, y~0), the H^1-projected atan2 loses accuracy
        # but the field evaluated FAR from the branch cut is correct.
        self.gf_chi.Set(atan2(y, x))
        self._step_psi(inverse=inverse)

        # Store exact chi gradient for accurate B evaluation everywhere.
        # This avoids the branch-cut errors in grad(gf_chi) without affecting
        # the psi-solve (B = grad(psi) x chi_grad_exact is used for field queries).
        r2 = x**2 + y**2 + 1e-30
        self._chi_grad_exact = CF((-y / r2, x / r2, 0.0))

    @property
    def B_cf_axisym(self):
        """B = grad(psi) x grad(chi_exact) using the analytic chi gradient.

        Preferred over ``B_cf`` after ``solve_axisymmetric()``: avoids
        branch-cut errors in the numerically-projected gf_chi near x < 0.
        """
        from ngsolve import grad
        if not hasattr(self, "_chi_grad_exact"):
            raise RuntimeError(
                "B_cf_axisym is only available after solve_axisymmetric().")
        return _cross_cf(grad(self.gf_psi), self._chi_grad_exact)

    def solve(self, niter=30, tol=1e-8, inverse="sparsecholesky", verbose=True):
        """Alternating Clebsch iteration until convergence or *niter* steps.

        Call this inside a ``with TaskManager():`` block.

        Parameters
        ----------
        niter : int
            Maximum number of alternating iterations.
        tol : float
            Convergence threshold for ``||dpsi||/||psi||``.
        inverse : str
            NGSolve direct solver name.  Default ``"sparsecholesky"``.
        verbose : bool
            Print iteration progress.

        Returns
        -------
        converged : bool
        """
        for it in range(niter):
            delta = self.step(inverse=inverse)
            if verbose:
                print(f"  Clebsch iter {it+1:3d}: |dpsi|/|psi| = {delta:.2e}")
            if delta < tol:
                if verbose:
                    print(f"  Converged after {it+1} iterations.")
                return True
        if verbose:
            print(f"  Warning: did not converge after {niter} iterations "
                  f"(last |dpsi|/|psi| = {delta:.2e}).")
        return False

    @property
    def B_cf(self):
        """B = grad(psi) x grad(chi) as a 3-component CoefficientFunction."""
        from ngsolve import grad
        return _cross_cf(grad(self.gf_psi), grad(self.gf_chi))

    @property
    def A_cf(self):
        """Vector potential A = psi * grad(chi)  (Clebsch gauge).

        Satisfies  curl(A) = B  exactly.
        """
        from ngsolve import grad
        gc = grad(self.gf_chi)
        return self.gf_psi * gc

    def residual_rms(self):
        """Compute RMS of the curl-curl residual: ||curl(B/mu) - J||_L2.

        Requires ``ngsolve.Integrate``.  Used for convergence diagnostics.
        Only meaningful after at least one :meth:`solve` call.
        """
        from ngsolve import Integrate, curl, CF, sqrt, dx
        B = self.B_cf
        mu = self.mu_cf
        H = B / mu
        # curl H - J  (3-component)
        res = curl(H) - self.J_cf
        val = Integrate(res * res, self.mesh)
        return float(np.sqrt(val))
