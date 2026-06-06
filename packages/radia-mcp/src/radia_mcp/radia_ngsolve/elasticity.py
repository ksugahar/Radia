"""2D linear ELASTICITY for radia-ngsolve -- the structural half of the
multiphysics couplings (thermo-mechanical thermal stress, magneto/electro-mechanical
force -> deflection). Small-strain isotropic elasticity, plane stress or plane strain,
with body force, boundary tractions, and an optional thermal pre-strain.

    -div(sigma) = f ,   sigma = 2 mu eps(u) + lam tr(eps(u)) I - sigma_thermal
    eps(u) = sym(grad u) ,   sigma_thermal = 2(mu+lam) alpha dT I  (isotropic).

Plane stress: lam = E nu/(1-nu^2);  plane strain: lam = E nu/((1+nu)(1-2nu));
mu = E/(2(1+nu)) for both.
"""
from ngsolve import (VectorH1, BilinearForm, LinearForm, GridFunction, InnerProduct,
                     CoefficientFunction, grad, dx, ds, Id, Trace)


def _lame(E, nu, plane):
    mu = E / (2.0 * (1.0 + nu))
    lam = E * nu / (1.0 - nu * nu) if plane == "stress" \
        else E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
    return mu, lam


def _eps(w):
    g = grad(w)
    return 0.5 * (g + g.trans)


def solve_linear_elasticity(mesh, E, nu, dirichletx="", dirichlety="", dirichlet="",
                            traction=None, body_force=None, thermal=None,
                            plane="stress", order=2):
    """Solve 2D small-strain isotropic elasticity. Returns the displacement
    VectorH1 GridFunction ``u`` (then strain ``sym(grad u)``, stress via
    :func:`stress_2d`).

    Per-component clamps: ``dirichletx`` fixes u_x = 0 on those boundaries,
    ``dirichlety`` fixes u_y = 0, ``dirichlet`` fixes BOTH (full clamp). Combine as
    needed (e.g. dirichletx="left", dirichlety="bottom" = the statically-determinate
    uniaxial setup). ``traction`` = {boundary: (tx, ty)} applied surface force/area
    [Pa]. ``body_force`` = (fx, fy) CF [N/m^3]. ``thermal`` = (alpha, dT) isotropic
    thermal pre-strain (dT scalar or CF) -> thermal stress 2(mu+lam) alpha dT.
    """
    mu, lam = _lame(E, nu, plane)
    dx_bc = "|".join(s for s in (dirichletx, dirichlet) if s)
    dy_bc = "|".join(s for s in (dirichlety, dirichlet) if s)
    fes = VectorH1(mesh, order=order, dirichletx=dx_bc, dirichlety=dy_bc)
    u, v = fes.TnT()
    Imat = Id(mesh.dim)

    def sig(w):
        e = _eps(w)
        return 2.0 * mu * e + lam * Trace(e) * Imat
    a = BilinearForm(fes)
    a += InnerProduct(sig(u), _eps(v)) * dx
    f = LinearForm(fes)
    if body_force is not None:
        f += InnerProduct(CoefficientFunction(body_force), v) * dx
    if traction is not None:
        for bnd, t in traction.items():
            f += InnerProduct(CoefficientFunction(tuple(t)), v) * ds(definedon=mesh.Boundaries(bnd))
    if thermal is not None:
        alpha, dT = thermal
        f += 2.0 * (mu + lam) * alpha * dT * Trace(_eps(v)) * dx   # thermal-stress load
    a.Assemble()
    f.Assemble()
    gu = GridFunction(fes)
    gu.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec
    return gu


def stress_2d(gu, E, nu, plane="stress", thermal=None):
    """Cauchy stress CF from a displacement solution ``gu``:
    sigma = 2 mu eps(u) + lam tr(eps) I - 2(mu+lam) alpha dT I (the thermal term
    subtracted when ``thermal=(alpha, dT)`` was applied). sigma[0,0]=sigma_xx, etc."""
    mu, lam = _lame(E, nu, plane)
    e = _eps(gu)
    s = 2.0 * mu * e + lam * Trace(e) * Id(2)
    if thermal is not None:
        alpha, dT = thermal
        s = s - 2.0 * (mu + lam) * alpha * dT * Id(2)
    return s


def cantilever_tip_deflection(w, length, E, I_area):
    """Euler-Bernoulli tip deflection of a CANTILEVER (clamped one end, free the
    other) under a UNIFORM transverse load ``w`` [N/m]:  delta = w L^4 / (8 E I),
    with ``I_area`` the second moment of area (rectangular per-unit-depth section
    height h: I = h^3/12). The slender-beam reference for the magneto-mechanical
    chain -- a current-carrying beam in field B0 feels w = I*B0 [N/m] (Lorentz)."""
    return w * length ** 4 / (8.0 * E * I_area)


def constrained_bar_thermal_stress(E, alpha, mean_dT):
    """Axial thermal stress in a bar fully constrained along its axis (both ends
    fixed in x, free laterally / plane stress) at temperature rise field dT(x):
    equilibrium (sigma_xx const) + zero net axial strain give

        sigma_xx = -E alpha <dT>

    where <dT> is the LENGTH-AVERAGE rise (for a parabolic Joule profile, <dT> =
    (2/3) dT_max). Compressive (negative) for heating. The structural end of the
    electro-thermo-mechanical chain (examples/comsol_class/electro_thermo_mech.py)."""
    return -E * alpha * mean_dT
