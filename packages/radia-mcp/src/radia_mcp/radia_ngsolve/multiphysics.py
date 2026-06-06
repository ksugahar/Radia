"""Multiphysics couplings for radia-ngsolve -- the COMSOL-class problems
(induction heating EM->thermal, ...) as REUSABLE building blocks rather than
one-off scripts.

These are the executable counterpart of the knowledge in ``multiphysics_usage``.
v1 ships the one-way magneto-thermal (induction-heating) coupling:

    eddy solve (A-Phi harmonic)  ->  joule_loss_density()  ->  solve_heat_steady()

validated against the analytic cylinder-in-axial-AC eddy loss (exact I0/I1
Bessel) and the 1-D radial heat equation (P/L 8.9 %, T_centre 9.6 %; see
examples/comsol_class/induction_heating.py).
"""
from ngsolve import (H1, BilinearForm, LinearForm, GridFunction, InnerProduct,
                     Conj, grad, dx)


def joule_loss_density(gfA, gfPhi, sigma_cf, omega, A0=None):
    """Time-averaged Joule (ohmic) loss density  q = 1/2 sigma |E|^2  [W/m^3]
    from a harmonic eddy-current solution, with  E = -j omega (A_total + grad Phi).

    CRITICAL gotcha (cost real debugging): in a SCATTERED-field eddy solve --
    where the source is a background potential ``A0`` (curl A0 = applied B0) and
    the returned ``gfA`` is the SCATTERED potential (so field probes do
    ``curl(gfA) + B0``) -- the TOTAL electric field needs A0 added back:
    pass ``A0`` (the background CF). Omitting it overestimates the loss by ~10x.
    For a COIL-driven (total-field) solve, ``gfA`` is already total -> A0=None.

    Args:
        gfA, gfPhi : components of the solve_eddy_current_harmonic_APhi result.
        sigma_cf   : conductivity CoefficientFunction (0 outside conductors).
        omega      : angular frequency [rad/s].
        A0         : background vector-potential CF for scattered solves, else None.

    Returns a real CoefficientFunction (0 where sigma=0).
    """
    A_total = gfA if A0 is None else (A0 + gfA)
    E = -1j * omega * (A_total + grad(gfPhi))
    return 0.5 * sigma_cf * InnerProduct(E, Conj(E)).real


def solve_heat_steady(mesh, q, k, conductor, dirichlet, order=2,
                      inverse="sparsecholesky"):
    """Steady heat conduction  -div(k grad T) = q  on the ``conductor`` material
    with  T = 0  on the ``dirichlet`` boundary (temperature RISE above a cooled
    surface). One-way coupled to an EM Joule source ``q`` (induction heating).

    Args:
        mesh       : NGSolve mesh.
        q          : volumetric heat source CF [W/m^3] (e.g. joule_loss_density).
        k          : thermal conductivity [W/(m K)] (scalar or CF).
        conductor  : material name the heat equation is solved on.
        dirichlet  : boundary name held at T=0 (the cooled surface).
        order      : H1 order (default 2).

    Returns the temperature GridFunction (defined on ``conductor``).
    """
    fesT = H1(mesh, order=order, definedon=mesh.Materials(conductor),
              dirichlet=dirichlet)
    T, s = fesT.TnT()
    a = BilinearForm(fesT); a += k * grad(T) * grad(s) * dx
    f = LinearForm(fesT);   f += q * s * dx
    a.Assemble(); f.Assemble()
    gT = GridFunction(fesT)
    gT.vec.data = a.mat.Inverse(fesT.FreeDofs(), inverse=inverse) * f.vec
    return gT


def joule_heat_source(gfV, sigma_cf):
    """DC Joule heat density  q = sigma |grad V|^2 = |J|^2 / sigma  [W/m^3] from a
    conduction solution ``gfV`` (J = -sigma grad V, e.g. from
    ``scalar_fem2d.solve_current_flow``). The DC twin of :func:`joule_loss_density`;
    feed it to :func:`solve_heat_steady` for an ELECTRO-THERMAL (Joule-heating)
    coupling: solve_current_flow -> joule_heat_source -> solve_heat_steady."""
    return sigma_cf * InnerProduct(grad(gfV), grad(gfV))


def joule_bar_temperature_rise(sigma, voltage, length, k, x):
    """Steady temperature RISE of a uniform Joule-heated bar -- length L along x,
    voltage V applied across it, BOTH ends held cold, sides insulated:

        dT(x) = (q/2k) x (L - x),  q = sigma (V/L)^2,  max dT = sigma V^2/(8k) at x=L/2.

    The peak rise is INDEPENDENT of length (q ~ 1/L^2, the L^2 cancels). ``x`` from a
    cold end; scalar or iterable. The exact reference for the electro-thermal coupling
    (examples/comsol_class/joule_heating.py)."""
    q = sigma * (voltage / length) ** 2

    def _dt(xx):
        return q / (2.0 * k) * xx * (length - xx)
    try:
        return [_dt(float(xx)) for xx in x]
    except TypeError:
        return _dt(float(x))
