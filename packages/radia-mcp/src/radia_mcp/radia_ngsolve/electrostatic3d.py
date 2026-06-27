"""3D electrostatics + capacitance on standard NGSolve H1 -- the COMSOL AC/DC
"Computing Capacitance" / Electrostatics interface, in 3D.

The 2D/axisymmetric twin lives in ``scalar_fem2d`` (FEMM csolv analog). This module
is the genuine 3D capability: solve the Laplace/Poisson potential with conductors as
Dirichlet boundaries, then read capacitance from the field ENERGY (the robust route
COMSOL uses for the capacitance matrix):

    -div(eps grad V) = rho ,  E = -grad V ,  W = 1/2 integral eps |grad V|^2 ,  C = 2W/V^2

Validated against the spherical-capacitor closed form C = 4 pi eps0 a b / (b - a).
"""
import math
import numpy as np
from ngsolve import (H1, BilinearForm, LinearForm, GridFunction, grad, dx,
                     Integrate, CoefficientFunction)

EPS0 = 8.8541878128e-12   # vacuum permittivity [F/m]


def solve_electrostatic_3d(mesh, eps_cf, bc_values, order=2, rho=None):
    """Electrostatic potential ``V`` with conductors imposed as Dirichlet boundaries.

    ``eps_cf``    : permittivity CoefficientFunction (``EPS0``, or ``mesh.MaterialCF``
                    for piecewise dielectrics ``eps0*eps_r``).
    ``bc_values`` : dict {boundary-name -> fixed potential [V]} (the conductors).
    ``rho``       : optional space-charge density CoefficientFunction (default 0).

    Returns the H1 GridFunction ``V`` (read ``E = -grad(V)``).
    """
    dirichlet = "|".join(bc_values.keys())
    fes = H1(mesh, order=order, dirichlet=dirichlet)
    gfV = GridFunction(fes)
    gfV.vec[:] = 0.0
    # Set ALL conductor potentials in ONE Set: gf.Set(definedon=...) zeros the whole
    # vector first, so a per-boundary loop would wipe each previous boundary's data.
    bnd_cf = mesh.BoundaryCF(bc_values, default=0.0)
    gfV.Set(bnd_cf, definedon=mesh.Boundaries(dirichlet))
    u, v = fes.TnT()
    a = BilinearForm(fes, symmetric=True)
    a += eps_cf * grad(u) * grad(v) * dx
    f = LinearForm(fes)
    if rho is not None:
        f += rho * v * dx
    a.Assemble(); f.Assemble()
    # Dirichlet lift: solve for the free-DOF correction to the prescribed boundary data
    r = f.vec.CreateVector()
    r.data = f.vec - a.mat * gfV.vec
    gfV.vec.data += a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * r
    return gfV


def capacitance_from_energy(gfV, eps_cf, voltage, mesh=None):
    """Two-conductor capacitance from the electrostatic energy:
    ``W = 1/2 integral eps |grad V|^2`` and ``C = 2 W / V^2`` [F]."""
    if mesh is None:
        mesh = gfV.space.mesh
    W = 0.5 * Integrate(eps_cf * grad(gfV) * grad(gfV) * dx, mesh)
    return 2.0 * W / voltage**2


def spherical_capacitor_C(a, b, eps_r=1.0):
    """Exact capacitance of a spherical capacitor (inner radius a, outer b):
    C = 4 pi eps0 eps_r a b / (b - a) [F]."""
    return 4.0 * math.pi * EPS0 * eps_r * a * b / (b - a)


def capacitance_matrix(mesh, eps_cf, conductors, order=2):
    """Maxwell capacitance matrix C (N x N) for the named conductor boundaries --
    the COMSOL "Computing Capacitance" capacitance matrix, in 3D.

    Energises each conductor to 1 V in turn (all others grounded), solves the
    electrostatic problem, and reads EVERY conductor's charge from the FEM REACTION
    Q_i = <K V, chi_i> (chi_i = the FE function 1 on conductor i) -- more accurate
    and robust than surface-integrating eps dV/dn. C[i, j] = Q_i when V_j = 1 V.

    Physical properties (good sanity checks): symmetric (C_ij = C_ji), diagonal
    C_ii > 0, off-diagonal C_ij < 0 (i != j), and each row sums to the conductor's
    capacitance to ground/infinity (0 for a closed system). ``conductors`` is the
    list of Dirichlet boundary names. Returns an (N, N) numpy array [F].
    """
    names = list(conductors)
    diri = "|".join(names)
    fes = H1(mesh, order=order, dirichlet=diri)
    u, v = fes.TnT()
    a = BilinearForm(fes, symmetric=True)
    a += eps_cf * grad(u) * grad(v) * dx
    a.Assemble()
    inv = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")

    # indicator gridfunction (=1 on each conductor's boundary trace) for the reaction
    chi = []
    for nm in names:
        c = GridFunction(fes); c.vec[:] = 0.0
        c.Set(1.0, definedon=mesh.Boundaries(nm))
        chi.append(c)

    C = np.zeros((len(names), len(names)))
    for j, nm in enumerate(names):
        gfV = GridFunction(fes); gfV.vec[:] = 0.0
        gfV.Set(mesh.BoundaryCF({nm: 1.0}, default=0.0), definedon=mesh.Boundaries(diri))
        res = gfV.vec.CreateVector(); res.data = -a.mat * gfV.vec
        gfV.vec.data += inv * res
        KV = gfV.vec.CreateVector(); KV.data = a.mat * gfV.vec    # reaction (charge) vector
        for i in range(len(names)):
            C[i, j] = KV.InnerProduct(chi[i].vec)
    return C
