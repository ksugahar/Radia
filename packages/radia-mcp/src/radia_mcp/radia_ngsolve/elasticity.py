"""2D linear ELASTICITY for radia-ngsolve -- the structural half of the
multiphysics couplings (thermo-mechanical thermal stress, magneto/electro-mechanical
force -> deflection). Small-strain isotropic elasticity, plane stress or plane strain,
with body force, boundary tractions, and an optional thermal pre-strain.

    -div(sigma) = f ,   sigma = 2 mu eps(u) + lam tr(eps(u)) I - sigma_thermal
    eps(u) = sym(grad u) ,   sigma_thermal = 2(mu+lam) alpha dT I  (isotropic).

Plane stress: lam = E nu/(1-nu^2);  plane strain: lam = E nu/((1+nu)(1-2nu));
mu = E/(2(1+nu)) for both.
"""
import math

from ngsolve import (VectorH1, BilinearForm, LinearForm, GridFunction, InnerProduct,
                     CoefficientFunction, grad, dx, ds, Id, Trace, ArnoldiSolver)


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


def elastic_modal_2d(mesh, E, nu, rho, n_modes, dirichletx="", dirichlety="",
                     dirichlet="", plane="stress", order=2, shift=-1.0):
    """Free-vibration (modal) ANGULAR FREQUENCIES omega [rad/s] of a 2D elastic body, by the
    generalised eigenproblem  K phi = omega^2 M phi  with the elasticity stiffness K and the
    consistent mass  M = int rho u.v . The dynamic / modal companion of
    :func:`solve_linear_elasticity` -- the structural-NVH eigenmodes an electromagnetic force
    spectrum can excite. Per-component clamps as in solve_linear_elasticity (constrain enough to
    remove the rigid-body modes, else they appear as omega ~ 0 and are dropped).

    Validated on the fixed-free AXIAL rod (nu=0, transverse motion suppressed via u_y clamps):
    omega_n = (2n-1) pi/(2L) sqrt(E/rho), to ~machine precision; and the fixed-fixed rod
    omega_n = n pi/L sqrt(E/rho). Returns the ascending list of omega (length ``n_modes``)."""
    mu, lam = _lame(E, nu, plane)
    dx_bc = "|".join(s for s in (dirichletx, dirichlet) if s)
    dy_bc = "|".join(s for s in (dirichlety, dirichlet) if s)
    fes = VectorH1(mesh, order=order, dirichletx=dx_bc, dirichlety=dy_bc)
    u, v = fes.TnT()
    Imat = Id(mesh.dim)

    def sig(w):
        e = _eps(w)
        return 2.0 * mu * e + lam * Trace(e) * Imat
    K = BilinearForm(InnerProduct(sig(u), _eps(v)) * dx).Assemble()
    M = BilinearForm(rho * InnerProduct(u, v) * dx).Assemble()
    gf = GridFunction(fes, multidim=n_modes + 6)
    lams = ArnoldiSolver(K.mat, M.mat, fes.FreeDofs(), list(gf.vecs), shift=shift)
    w2 = sorted(l.real for l in lams)
    w2 = [x for x in w2 if x > 1e-6 * max(w2)]
    return [math.sqrt(x) for x in w2[:n_modes]]


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


def thermal_bending_tip_deflection(alpha, dT_thickness, thickness, length):
    """Euler-Bernoulli tip deflection of a CANTILEVER with a LINEAR through-thickness
    temperature gradient ``dT_thickness`` (top-to-bottom difference) over ``thickness`` h.
    It bends STRESS-FREE into the constant curvature kappa = alpha*dT_thickness/h, so

        delta_tip = kappa L^2/2 = alpha * dT_thickness * length^2 / (2 h).

    The thermal-bimorph / bimetallic-strip actuator principle (homogeneous beam, linear
    gradient). Distinct from :func:`constrained_bar_thermal_stress` (AXIAL blocked expansion
    -> stress, no bending). Validated examples/comsol_class/thermal_bending.py to 0.24 %."""
    return alpha * dT_thickness * length ** 2 / (2.0 * thickness)


def bimetal_tip_deflection(alpha1, alpha2, dT, thickness, length):
    """Cantilever tip deflection of a BIMETALLIC STRIP -- two bonded layers of EQUAL
    thickness and modulus but different expansion (alpha1 bottom, alpha2 top) under a
    UNIFORM rise dT (Timoshenko 1925). The strip curves with

        kappa = 3 (alpha2 - alpha1) dT / (2 h)   (h = total thickness),

    so the cantilever tip delta = kappa L^2/2 = 3(alpha2-alpha1) dT L^2/(4 h). The
    thermostat / thermal-switch / thermal-relay principle. Distinct from
    :func:`thermal_bending_tip_deflection` (ONE material, through-thickness GRADIENT):
    here a uniform dT + a material Delta-alpha drives the bending. Pass alpha as a
    region-wise CF to solve_linear_elasticity(thermal=(alpha, dT)). Validated
    examples/comsol_class/bimetallic_strip.py to 1.5 %."""
    kappa = 3.0 * (alpha2 - alpha1) * dT / (2.0 * thickness)
    return kappa * length ** 2 / 2.0


def biaxial_thermal_stress(E, alpha, dT, nu):
    """Biaxial thermal stress of a plate fully constrained in-plane (u=0 on all edges)
    under a UNIFORM rise dT, PLANE STRESS:

        sigma_x = sigma_y = -E alpha dT / (1 - nu)     [Pa]  (compressive).

    The 2D constraint factor 1/(1-nu) -- vs the 1D uniaxial blocked bar
    -E alpha dT (:func:`constrained_bar_thermal_stress`). Exact (u=0 is the solution for a
    uniform eigenstrain with all edges clamped). (PLANE STRAIN fully constrained is
    -E alpha dT/(1-2 nu); NOTE radia's plane='strain' thermal load is currently (1+nu) too
    small -- flagged for fix -- so this helper/validation covers plane STRESS.)"""
    return -E * alpha * dT / (1.0 - nu)


def rule_of_mixtures_cte(moduli, areas, alphas):
    """Effective (stiffness-weighted) coefficient of thermal expansion of a bonded
    composite whose layers co-strain axially:

        alpha_eff = sum(E_i A_i alpha_i) / sum(E_i A_i).

    ``areas`` are the cross-sectional shares A_i (per unit depth: layer widths). The
    common free-expansion strain of the laminate is alpha_eff*dT. Stiffer / larger
    layers pull alpha_eff toward their own alpha (Voigt / iso-strain mixing rule)."""
    den = sum(E * A for E, A in zip(moduli, areas))
    return sum(E * A * al for E, A, al in zip(moduli, areas, alphas)) / den


def laminate_residual_stress(E_i, alpha_i, alpha_eff, dT):
    """INTERNAL (residual) axial stress in layer i of a FREE bonded composite under a
    uniform temperature rise ``dT`` -- no external constraint, the stress comes purely
    from the CTE mismatch forcing a common strain alpha_eff*dT (:func:`rule_of_mixtures_cte`):

        sigma_i = E_i (alpha_eff - alpha_i) dT.

    Layers with alpha_i > alpha_eff (want to expand more) end in COMPRESSION, the rest in
    tension; the cross-section is self-equilibrated, sum(A_i sigma_i) = 0. The composite-
    laminate / die-attach / coating residual-stress and bimetal-internal-stress workhorse;
    the FREE-bar counterpart of the externally blocked :func:`constrained_bar_thermal_stress`.
    Pass region-wise E and alpha CFs to solve_linear_elasticity(thermal=(alpha, dT)) and read
    sigma_xx per layer. Validated examples/comsol_class/laminate_thermal_stress.py to 0.00 %."""
    return E_i * (alpha_eff - alpha_i) * dT
