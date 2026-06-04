"""2D scalar-potential FEM (FEMM csolv / hsolv analogs) on standard NGSolve H1.

FEMM's electrostatics (``csolv``), heat flow (``hsolv``) and DC current flow are
all the SAME elliptic operator  -div(c grad u) = f  with a different material
coefficient ``c`` and a different post-processing (capacitance / thermal
resistance / conductance). This module ships that shared core plus thin,
physics-named wrappers, each validated against an analytical benchmark:

    electrostatics : -div(eps grad V) = rho ,  E=-grad V , C = 2W/V^2
    heat (steady)  : -div(k   grad T) = q   ,  q_flux=-k grad T
    current flow   : -div(sig grad V) = 0   ,  J=-sig grad V , G=2P/V^2

All are 2D planar (quantities per unit length in the out-of-plane direction).
"""
import math

from ngsolve import (H1, BilinearForm, LinearForm, GridFunction, grad, dx,
                     Integrate, CoefficientFunction, BND, x as _r)

EPS0 = 8.8541878128e-12


def solve_poisson_2d(mesh, coeff, dirichlet_values, source=None, order=2):
    """Core 2D elliptic solve  -div(coeff grad u) = source  with Dirichlet data.

    Parameters
    ----------
    mesh             : 2D mesh with named electrode/wall boundaries.
    coeff            : material coefficient CF (eps, k, or sigma).
    dirichlet_values : {boundary_name: value} fixed-potential boundaries.
    source           : volumetric source CF (rho, heat q) or None.
    order            : H1 order (default 2).

    Returns the H1 GridFunction ``u`` (potential / temperature). The gradient
    ``grad(u)`` gives -E / -q_flux/k / -J/sigma.
    """
    fes = H1(mesh, order=order, dirichlet="|".join(dirichlet_values))
    u, v = fes.TnT()
    a = BilinearForm(fes, symmetric=True)
    a += coeff * grad(u) * grad(v) * dx
    f = LinearForm(fes)
    if source is not None:
        f += source * v * dx
    a.Assemble()
    f.Assemble()
    gfu = GridFunction(fes)
    gfu.Set(mesh.BoundaryCF(dirichlet_values, default=0.0), BND)
    r = f.vec - a.mat * gfu.vec
    gfu.vec.data += a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * r
    return gfu


def solve_poisson_axi(mesh, coeff, dirichlet_values, source=None, order=2):
    """AXISYMMETRIC scalar elliptic solve  -div(coeff grad u) = source on an
    (r, z) half-plane mesh (r = x >= 0). Carries the toroidal Jacobian r:
    int coeff grad(u).grad(v) r dr dz = int source v r dr dz. The symmetry axis
    (r=0) is a natural (Neumann) boundary -- leave it unnamed. Use for FEMM
    axisymmetric electrostatics / heat (sphere, cone, disk electrodes).

    Returns the H1 GridFunction ``u``. (No 1/r singularity here: this is the
    SCALAR potential, unlike the magnetic A_phi which needs axihenrotte.)
    """
    fes = H1(mesh, order=order, dirichlet="|".join(dirichlet_values))
    u, v = fes.TnT()
    a = BilinearForm(fes, symmetric=True)
    a += coeff * grad(u) * grad(v) * _r * dx
    f = LinearForm(fes)
    if source is not None:
        f += source * v * _r * dx
    a.Assemble()
    f.Assemble()
    gfu = GridFunction(fes)
    gfu.Set(mesh.BoundaryCF(dirichlet_values, default=0.0), BND)
    res = f.vec - a.mat * gfu.vec
    gfu.vec.data += a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * res
    return gfu


def capacitance_axi(V, mesh, eps, v_applied):
    """Axisymmetric capacitance [F] (full 3D, not per length): C = 2W/V^2 with
    W = (1/2) int eps |grad V|^2 (2 pi r) dr dz. Validated: concentric-sphere
    capacitor C = 4 pi eps ab/(b-a) to 0.15 % (tests/test_axi_scalar.py)."""
    W = math.pi * Integrate(eps * grad(V) * grad(V) * _r * dx, mesh)
    return 2.0 * W / (v_applied * v_applied)


def solve_electrostatic(mesh, eps, potentials, charge=None, order=2):
    """FEMM ``csolv`` analog: -div(eps grad V) = rho. ``potentials`` =
    {electrode_boundary: volts}. Returns V; field E = -grad(V).

    Validated: coaxial capacitor C/L = 2 pi eps / ln(b/a) to 0.2 %
    (tests/test_electrostatic.py)."""
    return solve_poisson_2d(mesh, eps, potentials, source=charge, order=order)


def capacitance(V, mesh, eps, v_applied):
    """Capacitance per length [F/m] from the field energy: C = 2W/V^2,
    W = 1/2 int eps |grad V|^2 dA. ``v_applied`` = electrode potential difference."""
    W = 0.5 * Integrate(eps * grad(V) * grad(V) * dx, mesh)
    return 2.0 * W / (v_applied * v_applied)


def solve_thermal(mesh, k, temperatures, heat_source=None, order=2):
    """FEMM ``hsolv`` analog (steady state): -div(k grad T) = q.
    ``temperatures`` = {boundary: T}. Returns T; heat flux q = -k grad(T)."""
    return solve_poisson_2d(mesh, k, temperatures, source=heat_source, order=order)


def solve_current_flow(mesh, sigma, potentials, order=2):
    """DC conduction: -div(sigma grad V) = 0. ``potentials`` = {electrode: volts}.
    Returns V; current density J = -sigma grad(V)."""
    return solve_poisson_2d(mesh, sigma, potentials, source=None, order=order)


def conductance(V, mesh, sigma, v_applied):
    """Conductance per length [S/m] from ohmic power: G = 2P/V^2,
    P = 1/2 int sigma |grad V|^2 dA (DC, real)."""
    P = 0.5 * Integrate(sigma * grad(V) * grad(V) * dx, mesh)
    return 2.0 * P / (v_applied * v_applied)
