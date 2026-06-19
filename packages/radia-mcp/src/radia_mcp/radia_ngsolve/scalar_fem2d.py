"""2D scalar-potential FEM (FEMM csolv / hsolv analogs) on standard NGSolve H1.

FEMM's electrostatics (``csolv``), heat flow (``hsolv``) and DC current flow are
all the SAME elliptic operator  -div(c grad u) = f  with a different material
coefficient ``c`` and a different post-processing (capacitance / thermal
resistance / conductance). This module ships that shared core plus thin,
physics-named wrappers, each validated against an analytical benchmark:

    electrostatics  : -div(eps grad V)    = rho ,  E=-grad V    , C = 2W/V^2
    heat (steady)   : -div(k   grad T)    = q   ,  q_flux=-k grad T , G_th=2P/dT^2
    current flow    : -div(sig grad V)    = 0   ,  J=-sig grad V , G = 2P/V^2
    magnetic scalar : -div(mu  grad phi_m)= 0   ,  H=-grad phi_m, B=mu H , P = 2W/F^2

The last is the current-free magnetic scalar potential (COMSOL 'Magnetic Fields,
No Currents' / mfnc) -- the magnetic-circuit / reluctance-network primitive, with
mu = mu0*mu_r in place of eps / k / sigma. On the coaxial annulus all four give the
SAME radial-Laplace lumped value  2 pi c / ln(b/a)  (c = eps, k, sigma, mu).

All are 2D planar (quantities per unit length in the out-of-plane direction).
"""
import math

from ngsolve import (H1, BilinearForm, LinearForm, GridFunction, grad, dx, ds,
                     Integrate, CoefficientFunction, BND, x as _r)

EPS0 = 8.8541878128e-12
SIGMA_SB = 5.670374419e-8     # Stefan-Boltzmann constant [W/m^2/K^4]


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


def resistance_axi(V, mesh, sigma, v_applied):
    """Axisymmetric DC RESISTANCE [Ohm] (full 3D, not per length) between the electrodes that fix
    the solved potential ``V`` (from :func:`solve_poisson_axi` with coeff = sigma): by the power
    method  R = v_applied^2 / P  with  P = int sigma |grad V|^2 (2 pi r) dr dz  the total ohmic
    dissipation. The current-flow / resistance dual of :func:`capacitance_axi` (which is the
    electrostatic 2W/V^2). Validated on the uniform cylinder R = h/(sigma pi b^2) (exact, 1-D
    current) and the disk-contact SPREADING resistance R = 1/(4 sigma a) (Holm, ~1 % on a finite
    domain -- the contact-rim current singularity J ~ 1/sqrt(a^2-r^2) limits the rate)."""
    P = 2.0 * math.pi * Integrate(sigma * grad(V) * grad(V) * _r * dx, mesh)
    return v_applied * v_applied / P


def solve_electrostatic(mesh, eps, potentials, charge=None, order=2):
    """FEMM ``csolv`` analog: -div(eps grad V) = rho. ``potentials`` =
    {electrode_boundary: volts}. Returns V; field E = -grad(V).

    Validated: coaxial capacitor C/L = 2 pi eps / ln(b/a) to 0.2 %
    (tests/test_electrostatic.py).  The same coaxial C/L is independently cross-confirmed by an
    external FE electrostatics solver (~0.04 %) and by a curved high-order-hex 3D FE solve
    (~0 % at geometry order 3) -- a multi-method capacitance consistency."""
    return solve_poisson_2d(mesh, eps, potentials, source=charge, order=order)


def capacitance(V, mesh, eps, v_applied):
    """Capacitance per length [F/m] from the field energy: C = 2W/V^2,
    W = 1/2 int eps |grad V|^2 dA. ``v_applied`` = electrode potential difference."""
    W = 0.5 * Integrate(eps * grad(V) * grad(V) * dx, mesh)
    return 2.0 * W / (v_applied * v_applied)


def solve_thermal(mesh, k, temperatures, heat_source=None, order=2,
                  convection=None, radiation=None, relax=1.0, max_iter=60,
                  tol=1e-8):
    """FEMM ``hsolv`` analog (steady state): -div(k grad T) = q.
    ``temperatures`` = {boundary: T} fixed-temperature walls.
    ``convection`` = {boundary: (h, T_inf)} convective (Robin) walls,
    -k dT/dn = h (T - T_inf) [h in W/m^2/K].
    ``radiation`` = {boundary: (eps, T_inf)} radiative walls,
    -k dT/dn = eps*sigma_SB*(T^4 - T_inf^4) [T in KELVIN]. This is NONLINEAR;
    it is solved by Picard, linearizing each sweep as an effective Robin film
    coefficient  h_rad(T) = eps*sigma_SB*(T^2 + T_inf^2)*(T + T_inf)  evaluated
    at the previous temperature (combined with any convection on the same wall).
    Unnamed walls are insulated (natural Neumann). Returns T; q = -k grad(T).

    Validated: 1D slab Dirichlet + convection vs L/k + 1/h (+0.003%); Dirichlet
    + radiation vs the implicit conduction/radiation balance solved by brentq
    (tests/test_scalar_fem2d_ext.py)."""
    if radiation is None:
        return solve_poisson_2d(mesh, k, temperatures, source=heat_source,
                                order=order, robin=convection)

    base = dict(convection) if convection else {}
    gfu = solve_poisson_2d(mesh, k, temperatures, source=heat_source,
                           order=order, robin=base or None)
    prev = None
    for _ in range(max_iter):
        robin = dict(base)
        for bnd, (eps, Tinf) in radiation.items():
            h_rad = eps * SIGMA_SB * (gfu * gfu + Tinf * Tinf) * (gfu + Tinf)
            if bnd in robin:                       # convection + radiation wall
                h_c, T_c = robin[bnd]
                h_tot = h_c + h_rad
                robin[bnd] = (h_tot, (h_c * T_c + h_rad * Tinf) / h_tot)
            else:
                robin[bnd] = (h_rad, Tinf)
        gnew = solve_poisson_2d(mesh, k, temperatures, source=heat_source,
                                order=order, robin=robin)
        if relax != 1.0:
            gnew.vec.data = (1.0 - relax) * gfu.vec + relax * gnew.vec
        cur = Integrate(gnew * gnew * dx, mesh)
        gfu = gnew
        if prev is not None and abs(cur - prev) < tol * max(abs(cur), 1e-30):
            break
        prev = cur
    return gfu


def solve_current_flow(mesh, sigma, potentials, order=2):
    """DC conduction: -div(sigma grad V) = 0. ``potentials`` = {electrode: volts}.
    Returns V; current density J = -sigma grad(V)."""
    return solve_poisson_2d(mesh, sigma, potentials, source=None, order=order)


def conductance(V, mesh, sigma, v_applied):
    """Conductance per length [S/m] from ohmic power: G = 2P/V^2,
    P = 1/2 int sigma |grad V|^2 dA (DC, real)."""
    P = 0.5 * Integrate(sigma * grad(V) * grad(V) * dx, mesh)
    return 2.0 * P / (v_applied * v_applied)


def thermal_conductance(T, mesh, k, delta_T):
    """Thermal conductance per length [W/m/K] from the conduction dissipation: G_th = 2P/dT^2,
    P = 1/2 int k |grad T|^2 dA -- the heat member of the Laplace-operator triad alongside electric
    :func:`conductance` and :func:`capacitance` (same -div(c grad u)=0, c = sigma / eps / k).
    Validated on the coaxial annulus: G_th/L = 2 pi k / ln(b/a) (tests/test_scalar_fem2d.py)."""
    P = 0.5 * Integrate(k * grad(T) * grad(T) * dx, mesh)
    return 2.0 * P / (delta_T * delta_T)


def solve_magnetic_scalar(mesh, mu, scalar_potentials, order=2):
    """Current-free magnetic scalar potential (COMSOL 'Magnetic Fields, No Currents' / mfnc analog):
    -div(mu grad phi_m) = 0, with H = -grad(phi_m) and B = mu H (mu = mu0*mu_r). ``scalar_potentials``
    = {boundary: magnetomotive potential [A]} fixed-MMF boundaries. The fourth member of the FEMM-style
    scalar-Laplace family (electrostatic / current-flow / thermal): the reduced scalar potential for
    CURRENT-FREE magnetics (gaps, PM exteriors) where a full vector potential A is unnecessary -- the
    magnetic-circuit / reluctance-network primitive. Returns phi_m; the permeance via :func:`permeance`."""
    return solve_poisson_2d(mesh, mu, scalar_potentials, source=None, order=order)


def permeance(phi, mesh, mu, mmf):
    """Magnetic permeance per length [H/m] from the field co-energy: P' = 2W/F^2,
    W = 1/2 int mu |grad phi_m|^2 dA -- the magnetic member of the scalar-Laplace family alongside
    electric :func:`conductance`, :func:`capacitance` and :func:`thermal_conductance` (same
    -div(c grad u)=0, c = sigma / eps / k / mu). For a finite axial length L the magnetic-circuit
    reluctance is R_m = 1/(P' L).  Validated on the coaxial annulus: P'/L = 2 pi mu / ln(b/a)
    (tests/test_scalar_fem2d.py) -- the radial-gap reluctance used in magnetic-equivalent-circuit models."""
    W = 0.5 * Integrate(mu * grad(phi) * grad(phi) * dx, mesh)
    return 2.0 * W / (mmf * mmf)


def solve_magnetostatic_az(mesh, nu, currents, dirichlet_values, order=2):
    """In-plane (2D planar) magnetostatic vector potential A_z (FEMM 'magnetics' /
    COMSOL 'mf' analog): -div(nu grad A_z) = J_z, with reluctivity nu = 1/(mu0 mu_r),
    out-of-plane current density J_z [A/m^2] and flux density B = (dA_z/dy, -dA_z/dx)
    so |B| = |grad A_z|. This is the CURRENT-carrying (vector-potential) member of the
    same elliptic operator -div(c grad u) = f as the scalar potentials in this module --
    the magnetostatic primitive for slot conductors, coils, busbars and the per-length
    inductance below (whereas :func:`solve_magnetic_scalar` is the current-free mfnc member).

    ``currents``         = {region_name: J_z [A/m^2]} drive density per material region
                           (0 in regions not listed).
    ``dirichlet_values`` = {boundary: A_z}: a far boundary A_z=0 (flux-parallel, FEMM's
                           default outer boundary) or a perfectly-conducting flux-return shell.

    Returns the H1 GridFunction A_z; the inductance via :func:`inductance_2d`. Validated:
    solid-conductor coax L' = mu0/(2pi)(ln(b/a)+1/4) and the magnetic-fill variant
    L' = mu0/(2pi)(mu_r ln(b/a)+1/4) (tests/test_inductance_2d.py)."""
    src = mesh.MaterialCF(dict(currents), default=0.0)
    return solve_poisson_2d(mesh, nu, dirichlet_values, source=src, order=order)


def solve_magnetostatic_az_magnet(mesh, nu, magnetization, dirichlet_values,
                                  currents=None, order=2):
    """2D planar magnetostatic A_z with a PERMANENT-MAGNET (magnetization) source -- the
    PM-magnetics extension of :func:`solve_magnetostatic_az` (which carries free currents only).

    A uniform magnetization ``M = (Mx, My)`` [A/m] in a region is the curl-of-M source of
    ``-div(nu grad A_z) = (curl M)_z``; integrated by parts it enters the weak form as

        int_region ( Mx dv/dy - My dv/dx ) dx      (the equivalent surface current M x n),

    which captures the magnet's bound surface current automatically for a uniform M. With
    ``nu = 1/mu0`` (a mu_r=1 magnet) the solution is the free-space magnet field; pass ``currents``
    too for combined coil+PM problems. ``magnetization`` = {region: (Mx, My)};
    ``currents`` = {region: J_z [A/m^2]} (optional). Flux density B = (dA_z/dy, -dA_z/dx).
    Validated against the closed form ``radia.analytical_formulas.rect_magnet_2d`` (2D bar magnet)."""
    Mx = mesh.MaterialCF({r: m[0] for r, m in magnetization.items()}, default=0.0)
    My = mesh.MaterialCF({r: m[1] for r, m in magnetization.items()}, default=0.0)
    fes = H1(mesh, order=order, dirichlet="|".join(dirichlet_values))
    u, v = fes.TnT()
    a = BilinearForm(nu * grad(u) * grad(v) * dx, symmetric=True).Assemble()
    gv = grad(v)
    lf = (Mx * gv[1] - My * gv[0]) * dx                       # magnetization (curl M) source
    if currents:
        Jz = mesh.MaterialCF(dict(currents), default=0.0)
        lf = lf + Jz * v * dx
    f = LinearForm(lf).Assemble()
    gfu = GridFunction(fes)
    gfu.Set(mesh.BoundaryCF(dirichlet_values, default=0.0), BND)
    res = f.vec - a.mat * gfu.vec
    gfu.vec.data += a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * res
    return gfu


def inductance_2d(A, mesh, nu, total_current):
    """Inductance per length [H/m] from the magnetic field energy: L' = 2W/I^2,
    W = 1/2 int nu |grad A_z|^2 dA = 1/2 int |B|^2/mu dA -- the magnetostatic
    (vector-potential) member of the energy-method family alongside electric
    :func:`capacitance` (C=2W/V^2), :func:`conductance` (G=2P/V^2) and magnetic-scalar
    :func:`permeance` (P'=2W/F^2). ``nu`` = 1/(mu0 mu_r) reluctivity field: a *frozen* nu
    field (the converged reluctivity of a saturated nonlinear solve, held fixed) makes this
    the frozen-permeability incremental inductance used for machine Ld/Lq; a constant mu_r is
    the linear limit. ``total_current`` = net conductor current I [A]. Validated against the
    solid-coax closed form mu0/(2pi)(ln(b/a)+1/4) to 0.16% (tests/test_inductance_2d.py)."""
    W = 0.5 * Integrate(nu * grad(A) * grad(A) * dx, mesh)
    return 2.0 * W / (total_current * total_current)


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
