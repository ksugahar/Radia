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

from ngsolve import (H1, BilinearForm, LinearForm, GridFunction, grad, dx, ds,
                     Integrate, CoefficientFunction, BND, x as _r)

EPS0 = 8.8541878128e-12


def solve_poisson_2d(mesh, coeff, dirichlet_values, source=None, order=2,
                     robin=None):
    """Core 2D elliptic solve  -div(coeff grad u) = source  with Dirichlet and
    optional Robin (mixed / convection) boundary data.

    Parameters
    ----------
    mesh             : 2D mesh with named electrode/wall boundaries.
    coeff            : material coefficient CF (eps, k, or sigma).
    dirichlet_values : {boundary_name: value} fixed-potential boundaries.
    source           : volumetric source CF (rho, heat q) or None.
    order            : H1 order (default 2).
    robin            : {boundary: (h, u_inf)} Robin / convection boundary
                       -coeff du/dn = h (u - u_inf); adds  int_G h u v ds to the
                       stiffness and  int_G h u_inf v ds to the load. Unnamed
                       boundaries stay natural (Neumann, zero flux / insulated).

    Returns the H1 GridFunction ``u`` (potential / temperature). The gradient
    ``grad(u)`` gives -E / -q_flux/k / -J/sigma.
    """
    fes = H1(mesh, order=order, dirichlet="|".join(dirichlet_values))
    u, v = fes.TnT()
    a = BilinearForm(fes, symmetric=True)
    a += coeff * grad(u) * grad(v) * dx
    if robin:
        for bnd, (h, _uinf) in robin.items():
            a += h * u * v * ds(definedon=mesh.Boundaries(bnd))
    f = LinearForm(fes)
    if source is not None:
        f += source * v * dx
    if robin:
        for bnd, (h, uinf) in robin.items():
            f += h * uinf * v * ds(definedon=mesh.Boundaries(bnd))
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


def solve_thermal(mesh, k, temperatures, heat_source=None, order=2,
                  convection=None):
    """FEMM ``hsolv`` analog (steady state): -div(k grad T) = q.
    ``temperatures`` = {boundary: T} fixed-temperature walls.
    ``convection`` = {boundary: (h, T_inf)} convective (Robin) walls,
    -k dT/dn = h (T - T_inf) [h in W/m^2/K]. Unnamed walls are insulated
    (natural Neumann). Returns T; heat flux q = -k grad(T).

    Validated: 1D slab Dirichlet + convection, surface T and flux vs the
    series thermal resistance L/k + 1/h (tests/test_scalar_fem2d_ext.py)."""
    return solve_poisson_2d(mesh, k, temperatures, source=heat_source,
                            order=order, robin=convection)


def solve_current_flow(mesh, sigma, potentials, order=2):
    """DC conduction: -div(sigma grad V) = 0. ``potentials`` = {electrode: volts}.
    Returns V; current density J = -sigma grad(V)."""
    return solve_poisson_2d(mesh, sigma, potentials, source=None, order=order)


def conductance(V, mesh, sigma, v_applied):
    """Conductance per length [S/m] from ohmic power: G = 2P/V^2,
    P = 1/2 int sigma |grad V|^2 dA (DC, real)."""
    P = 0.5 * Integrate(sigma * grad(V) * grad(V) * dx, mesh)
    return 2.0 * P / (v_applied * v_applied)


def solve_current_flow_ac(mesh, sigma, eps, omega, potentials, order=2):
    """FEMM AC current-flow: time-harmonic conduction in a lossy dielectric,
    -div((sigma + j w eps) grad V) = 0 with COMPLEX potential V.

    ``sigma`` (conduction) and ``eps`` (permittivity) are material CFs; ``omega``
    the angular frequency. ``potentials`` = {electrode: volts} (real or complex).
    Returns the complex H1 GridFunction V (E = -grad V; J = (sigma+j w eps) E).
    Terminal admittance via :func:`admittance_ac`.

    Validated: coaxial Y/L = 2 pi (sigma + j w eps) / ln(b/a) = G + j w C
    (tests/test_scalar_fem2d_ext.py)."""
    c = CoefficientFunction(sigma + 1j * omega * eps)
    fes = H1(mesh, order=order, complex=True, dirichlet="|".join(potentials))
    u, v = fes.TnT()
    a = BilinearForm(fes, symmetric=True)
    a += c * grad(u) * grad(v) * dx
    f = LinearForm(fes)
    a.Assemble()
    f.Assemble()
    gfu = GridFunction(fes)
    gfu.Set(mesh.BoundaryCF({b: complex(val) for b, val in potentials.items()},
                            default=0.0), BND)
    r = f.vec - a.mat * gfu.vec
    gfu.vec.data += a.mat.Inverse(fes.FreeDofs(), inverse="umfpack") * r
    return gfu


def admittance_ac(V, mesh, sigma, eps, omega, v_applied):
    """Terminal admittance per length [S/m] of an AC current-flow solve:
    Y = (1/V^2) int (sigma + j w eps) grad V . grad V dA = G + j w C.
    (Non-conjugated product -- this is I/V, not power.)"""
    c = CoefficientFunction(sigma + 1j * omega * eps)
    return Integrate(c * grad(V) * grad(V) * dx, mesh) / (v_applied * v_applied)
