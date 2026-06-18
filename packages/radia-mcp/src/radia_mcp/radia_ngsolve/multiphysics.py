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
import math

from ngsolve import (H1, BilinearForm, LinearForm, GridFunction, InnerProduct,
                     Conj, CoefficientFunction, grad, dx, ds, BND)


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


def solve_heat_transient(mesh, T0, k, rho_cp, dirichlet, t_end, n_steps,
                         source=None, conductor=None, order=2, theta=0.5,
                         inverse="sparsecholesky"):
    """TRANSIENT heat conduction  rho_cp dT/dt = div(k grad T) + source  by the theta time
    scheme -- the time-domain complement of :func:`solve_heat_steady`. Homogeneous Dirichlet
    ``T = 0`` on the ``dirichlet`` boundary (the initial field ``T0`` must vanish there).

    The theta-rule discretisation (theta = 0.5 = Crank-Nicolson, 2nd order in dt; theta = 1 =
    implicit Euler, unconditionally stable, 1st order) reads, with mass M = int rho_cp u v and
    stiffness A = int k grad u . grad v,

        (M + theta dt A) T^{n+1} = (M - (1-theta) dt A) T^n + dt F .

    Args:
        mesh      : NGSolve mesh.
        T0        : initial-temperature CoefficientFunction (set at t=0).
        k         : thermal conductivity [W/(m K)] (scalar or CF).
        rho_cp    : volumetric heat capacity rho*c_p [J/(m^3 K)] (scalar or CF); the
                    diffusivity is alpha = k / rho_cp.
        dirichlet : boundary held at T=0.
        t_end, n_steps : final time and number of equal steps (dt = t_end/n_steps).
        source    : volumetric heat source CF [W/m^3] (steady; e.g. a Joule density), or None.
        conductor : restrict the solve to this material (default: whole mesh).
        order, theta, inverse : H1 order, time-scheme weight, sparse factoriser.

    Returns the temperature GridFunction at ``t_end``. Verified on the single-mode decay of the
    unit square (T0 = sin(pi x) sin(pi y) -> exp(-alpha 2 pi^2 t), Crank-Nicolson 2nd-order in dt).
    """
    kw = {"definedon": mesh.Materials(conductor)} if conductor else {}
    fesT = H1(mesh, order=order, dirichlet=dirichlet, **kw)
    T, s = fesT.TnT()
    A = BilinearForm(k * grad(T) * grad(s) * dx); A.Assemble()
    M = BilinearForm(rho_cp * T * s * dx); M.Assemble()
    gT = GridFunction(fesT); gT.Set(T0)
    dt = t_end / n_steps
    mstar = M.mat.CreateMatrix()
    mstar.AsVector().data = M.mat.AsVector() + theta * dt * A.mat.AsVector()
    inv = mstar.Inverse(fesT.FreeDofs(), inverse=inverse)
    fvec = LinearForm(source * s * dx).Assemble().vec if source is not None else None
    rhs = gT.vec.CreateVector()
    for _ in range(n_steps):
        rhs.data = M.mat * gT.vec - (1.0 - theta) * dt * (A.mat * gT.vec)
        if fvec is not None:
            rhs.data += dt * fvec
        gT.vec.data = inv * rhs
    return gT


def thermal_penetration_depth(k, rho_cp, omega):
    """Thermal-wave penetration depth  delta = sqrt(2 alpha/omega),  alpha = k/rho_cp [m] -- the
    1/e decay length of a temperature oscillation at angular frequency ``omega`` (the thermal
    analog of the EM skin depth). Used by :func:`solve_heat_harmonic` / lock-in thermography."""
    return math.sqrt(2.0 * (k / rho_cp) / omega)


def solve_heat_harmonic(mesh, k, rho_cp, omega, dirichlet_values, source=None, order=2,
                        inverse="umfpack"):
    """FREQUENCY-DOMAIN (harmonic) heat conduction -- the COMPLEX temperature phasor T(x) for a
    boundary oscillating at angular frequency ``omega``:

        k grad^2 T = i omega rho_cp T  (+ source) .

    The frequency-domain sibling of :func:`solve_heat_transient`; it shares the complex-diffusion
    structure of the AC eddy skin solve (``eddy2d``). It launches damped THERMAL WAVES of
    penetration depth delta = sqrt(2 alpha/omega) (alpha = k/rho_cp,
    :func:`thermal_penetration_depth`): a semi-infinite solid driven T = T0 on the surface gives
    T(x) = T0 exp(-(1+i) x/delta), i.e. |T| = T0 exp(-x/delta) with phase lag x/delta -- the basis
    of lock-in thermography / photothermal depth profiling.

    ``dirichlet_values`` = {boundary: complex T}. Returns the complex temperature GridFunction."""
    fes = H1(mesh, order=order, complex=True, dirichlet="|".join(dirichlet_values))
    T, w = fes.TnT()
    a = BilinearForm(k * grad(T) * grad(w) * dx + 1j * omega * rho_cp * T * w * dx)
    f = LinearForm(fes)
    if source is not None:
        f += source * w * dx
    a.Assemble(); f.Assemble()
    gfu = GridFunction(fes)
    gfu.Set(mesh.BoundaryCF(dirichlet_values, default=0.0), BND)
    r = f.vec - a.mat * gfu.vec
    gfu.vec.data += a.mat.Inverse(fes.FreeDofs(), inverse=inverse) * r
    return gfu


def joule_heated_cylinder_peak_dT(q, a, k):
    """Peak (centre) temperature RISE of a round conductor (radius ``a``) carrying a UNIFORM
    DC Joule source ``q`` [W/m^3] (q = J^2/sigma, uniform current density), with the surface
    held at the coolant temperature. Steady radial conduction -div(k grad T)=q gives the
    PARABOLIC profile dT(r) = q (a^2 - r^2)/(4 k), so the centre rise is

        dT_peak = q a^2 / (4 k).

    The cylindrical sibling of the slab electro-thermal (#19, dT_peak = sigma V^2/(8 k)); the
    classic current-carrying-wire / ampacity temperature rise. Validated LIVE 3-way (NGSolve
    :func:`solve_heat_steady` on a disk == this closed form == the reference commercial FE
    code, all 0.00 %) -- examples/comsol_class/cylinder_joule_heating.py,
    tests/test_cylinder_joule_heating.py."""
    return q * a * a / (4.0 * k)


def solve_heat_steady_robin(mesh, q, k, conductor, robin, h, T_inf=0.0, order=2,
                            dirichlet="", dirichlet_value=0.0, inverse="sparsecholesky"):
    """Steady heat conduction with a CONVECTIVE (Robin / Newton-cooling) boundary on ``robin``:

        -k dT/dn = h (T - T_inf)        (film coefficient ``h`` [W/(m^2 K)]),

    i.e. the cooled surface is NOT held at a fixed temperature -- it FLOATS to balance the
    convected heat, so the surface sits a film rise q*(half-thickness)/h above the coolant.
    The realistic cooling BC (vs the idealised T=0 of :func:`solve_heat_steady`). Weak form
    adds ``h T s`` to the stiffness and ``h T_inf s`` to the load on ``robin``; the remaining
    boundaries are natural (insulated) unless named in ``dirichlet`` (T=0). Returns the
    temperature (rise above ``T_inf``, which defaults to 0).
    """
    fesT = H1(mesh, order=order, definedon=mesh.Materials(conductor), dirichlet=dirichlet)
    T, s = fesT.TnT()
    a = BilinearForm(fesT)
    a += k * grad(T) * grad(s) * dx
    a += h * T * s * ds(definedon=mesh.Boundaries(robin))
    f = LinearForm(fesT)
    f += q * s * dx
    if T_inf != 0.0:
        f += h * T_inf * s * ds(definedon=mesh.Boundaries(robin))
    a.Assemble(); f.Assemble()
    gT = GridFunction(fesT)
    if dirichlet and dirichlet_value != 0.0:               # inhomogeneous Dirichlet (e.g. fin base)
        gT.Set(CoefficientFunction(dirichlet_value), definedon=mesh.Boundaries(dirichlet))
        r = f.vec.CreateVector()
        r.data = f.vec - a.mat * gT.vec
        gT.vec.data += a.mat.Inverse(fesT.FreeDofs(), inverse=inverse) * r
    else:
        gT.vec.data = a.mat.Inverse(fesT.FreeDofs(), inverse=inverse) * f.vec
    return gT


def convective_slab_peak_dT(q, thickness, k, h):
    """Centre temperature RISE of a slab (thickness ``L``) with a UNIFORM volumetric source
    ``q`` [W/m^3] cooled by CONVECTION (film coefficient ``h``) on BOTH faces:

        dT_centre = q L^2/(8 k)  +  q L/(2 h),

    the CONDUCTION parabola (#19/#41, q L^2/8k) PLUS the convective FILM rise q L/(2 h) (the
    surface floats this far above the coolant). The surface rise alone is q L/(2 h); as
    h -> inf the film vanishes and the fixed-T result q L^2/8k is recovered. The realistic
    Joule-heated / cooled-device temperature rise. Validated LIVE 3-way (NGSolve
    :func:`solve_heat_steady_robin` == this closed form == the reference commercial FE code)
    -- examples/comsol_class/convective_electrothermal.py, tests/test_convective_electrothermal.py."""
    return q * thickness * thickness / (8.0 * k) + q * thickness / (2.0 * h)


def thermal_resistance_series(thicknesses, conductivities):
    """Series thermal resistance per unit area of a MULTILAYER stack (1-D conduction):

        R_total = sum_i L_i / k_i   [m^2 K / W],

    the thermal-circuit analog of series electrical resistances. The temperature DROP across
    layer i carrying a heat flux Q [W/m^2] is Q*L_i/k_i, so the junction-to-base rise of a
    die stack is Q * R_total (plus q*L1^2/(2 k1) for self-heating WITHIN a uniformly
    heat-generating top layer). The electronics-cooling / multilayer-wall workhorse:
    high-k layers (Cu spreader) drop little, low-k layers (solder, TIM) dominate. Validated
    LIVE 3-way (NGSolve solve_heat_steady region-wise k == this network == the reference FE
    code) -- examples/comsol_class/die_stack_thermal.py, tests/test_die_stack_thermal.py."""
    return sum(L / k for L, k in zip(thicknesses, conductivities))


def fin_parameter(h, k, thickness):
    """Fin parameter m = sqrt(h P /(k A)) for a thin straight fin convecting from BOTH faces
    (per unit depth: perimeter P = 2*depth, cross-section A = thickness*depth), so

        m = sqrt(2 h /(k * thickness))   [1/m].

    Sets the conduction-vs-convection length scale 1/m; valid when the Biot number
    h*thickness/k << 1 (the fin is thin enough to be ~isothermal across its thickness)."""
    import math as _m
    return _m.sqrt(2.0 * h / (k * thickness))


def fin_tip_temperature_rise(dT_base, m, length):
    """Tip temperature RISE of a straight fin with an ADIABATIC (insulated) tip:

        dT(x) = dT_base * cosh(m (L - x))/cosh(m L)   ->   dT_tip = dT_base / cosh(m L).

    ``dT_base`` = base rise above coolant, ``m`` = :func:`fin_parameter`. The fin temperature
    falls from base to tip as the surface convects heat away."""
    import math as _m
    return dT_base / _m.cosh(m * length)


def fin_efficiency(m, length):
    """Straight-fin EFFICIENCY (adiabatic tip):

        eta_fin = tanh(m L)/(m L)

    = actual heat dissipated / the heat if the WHOLE fin were at the base temperature. eta -> 1
    for a short/fat/conductive fin (mL -> 0) and -> 1/(mL) -> 0 for a long/thin/insulating one
    -- the classic heat-sink design trade-off. Validated LIVE 3-way (NGSolve
    :func:`solve_heat_steady_robin` with a base Dirichlet + convective faces == this closed form
    == the reference FE code) -- examples/comsol_class/cooling_fin.py, tests/test_cooling_fin.py."""
    import math as _m
    return _m.tanh(m * length) / (m * length)


def solve_heat_steady_nonlinear(mesh, q_of_T, k, conductor, dirichlet, order=2,
                                relax=1.0, max_iter=60, tol=1e-7):
    """Steady heat with a TEMPERATURE-DEPENDENT volumetric source -- Picard fixed point on
    ``-div(k grad T) = q(T)``.  ``q_of_T`` is callable(T_cf) -> q_cf (the source given the
    current temperature), e.g. a TWO-WAY electro-thermal sigma(T) Joule density
    ``q = J^2/sigma(T) = (J^2/sigma0)(1+alpha T)`` (electrical resistivity rises with T, so
    Joule heating feeds back into the temperature), or a radiative/reaction volumetric term.

    The 2-way twin of :func:`solve_heat_steady`.  Returns the converged temperature
    GridFunction (T=0 on ``dirichlet``).  Validated EXACT (0.00 %) against the closed form
    for sigma(T)=sigma0/(1+alpha T): T_max=(1/alpha)(sec(sqrt(b) L/2)-1), b=alpha J^2/(sigma0 k)
    (examples/comsol_class/electrothermal_sigmaT.py)."""
    from ngsolve import CoefficientFunction, Integrate, dx as _dx
    Tprev = CoefficientFunction(0.0)
    gT = None
    last = None
    dom = _dx(definedon=mesh.Materials(conductor))
    for it in range(max_iter):
        gT = solve_heat_steady(mesh, q_of_T(Tprev), k, conductor, dirichlet, order=order)
        cur = Integrate(gT * gT * dom, mesh)
        if last is not None and abs(cur - last) <= tol * max(1.0, abs(cur)) and it >= 2:
            break
        last = cur
        Tprev = gT if relax >= 1.0 else (1.0 - relax) * Tprev + relax * gT
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


LORENZ = 2.44e-8          # Wiedemann-Franz Lorenz number [V^2/K^2] (k/(sigma T))


def phi_theta_local_temperature(V, U, T0=293.15, lorenz=LORENZ):
    """The phi-theta (voltage-temperature) relation -- the temperature at the point of local
    potential ``V`` in a current-carrying conductor whose two terminals (at potentials U and 0)
    are both held at ``T0``:

        T(V)^2 = T0^2 + V (U - V)/L                  (L = Wiedemann-Franz Lorenz number).

    A REMARKABLE result: under the Wiedemann-Franz law (k = L sigma T) the temperature is a
    function of the POTENTIAL ALONE -- INDEPENDENT of the conductor's geometry, conductivity, or
    current. Derived by the Kirchhoff transform psi = L sigma T^2/2 (linearises the coupled
    electro-thermal problem: -lap(psi) = sigma|grad V|^2, with psi quadratic in V)."""
    import math as _m
    return _m.sqrt(T0 * T0 + V * (U - V) / lorenz)


def phi_theta_max_temperature(U, T0=293.15, lorenz=LORENZ):
    """Maximum (hot-spot) temperature of a current-carrying conductor / electrical CONSTRICTION
    with both terminals at ``T0`` and total voltage ``U`` between them:

        T_max^2 - T0^2 = U^2/(4 L),   at the V = U/2 surface.

    The constriction "supertemperature" -- geometry-INDEPENDENT (Kohlrausch / Diesselhorst /
    Holm contact theory). The Wiedemann-Franz-exact generalisation of the constant-property bar
    peak sigma V^2/8k (:func:`joule_bar_temperature_rise`): at small rise
    T_max^2-T0^2 ~ 2 T0 (T_max-T0) -> T_max-T0 ~ U^2/(8 L T0). Use :func:`melting_voltage` for the
    inverse (the contact voltage that softens/melts the metal)."""
    return phi_theta_local_temperature(0.5 * U, U, T0, lorenz)


def melting_voltage(T_melt, T0=293.15, lorenz=LORENZ):
    """The geometry-INDEPENDENT contact VOLTAGE at which an electrical constriction reaches the
    temperature ``T_melt`` (e.g. the softening or melting point of the metal):

        U = sqrt(4 L (T_melt^2 - T0^2)).

    The inverse of :func:`phi_theta_max_temperature`. The classic contact result that a metal
    softens/melts at a CHARACTERISTIC voltage set only by the material (Lorenz number + melting
    point), NOT by the contact size or force -- e.g. copper softens ~0.12 V, melts ~0.43 V; the
    basis of contact-voltage limits in connectors, relays, and motor terminals."""
    import math as _m
    return _m.sqrt(4.0 * lorenz * (T_melt * T_melt - T0 * T0))
