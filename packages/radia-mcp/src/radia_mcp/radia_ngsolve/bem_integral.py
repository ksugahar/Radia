# -*- coding: utf-8 -*-
r"""ngsolve.bem INTEGRAL-EQUATION (BEM) solver helpers -- the OPEN boundary-element / integral solver,
the radia-ngsolve counterpart of a commercial high-frequency integral-equation (MoM) solver.

The boundary element method solves an exterior PDE (Laplace / Helmholtz / Maxwell) with unknowns on the
BODY SURFACE only (free space implicit to infinity) -- the same surface-charge / surface-current idea as
radia's MMM and the ELF BEM, lifted to a full integral-equation formulation.  ngsolve.bem ships the
Laplace, Helmholtz (scalar wave, wavenumber kappa) and Maxwell (full-wave EFIE/MFIE) single/double-layer
potential operators, curved + high-order, FMM/H-matrix accelerated.

This module starts the radia-ngsolve BEM capability with the canonical ELECTROSTATIC verification (the
kappa=0 static limit of the high-frequency integral solver): the CAPACITANCE of a conductor via the
Laplace single-layer indirect method.  The high-frequency members (Helmholtz Mie scattering, Maxwell RCS)
build on the same operator pattern.

Reference (analytic): an isolated conducting sphere of radius R has capacitance C = 4 pi eps0 R; in the
normalised units used here (eps0 = 1) the unit sphere gives C = 4 pi.
"""
import cmath
import math
import numpy as np
import ngsolve as ng
import ngsolve.bem as bem
from ngsolve import GridFunction, LinearForm, BilinearForm, ds, Integrate, specialcf, Trace
from ngsolve.krylovspace import CGSolver
from netgen.occ import Sphere, Pnt, OCCGeometry
import netgen.occ as occ


def single_layer_capacitance(mesh, order=3, intorder=12, curve_order=4, tol=1e-9, maxiter=600):
    """Capacitance of the conductor(s) bounded by ``mesh`` (a closed surface) by the Laplace SINGLE-LAYER
    indirect BEM.  Sets the conductor potential to 1, solves the first-kind integral equation V sigma = 1
    for the surface charge density sigma (CG with an L2-mass preconditioner -- V is SPD), and returns the
    total charge Q = integral(sigma) = C (since the potential is unity).  ``eps0 = 1`` (normalised)."""
    fes = ng.SurfaceL2(mesh, order=order, dual_mapping=True)
    u, v = fes.TnT()
    V = bem.SingleLayerPotentialOperator(fes, intorder=intorder)          # Laplace single layer (SPD)
    rhs = LinearForm(1.0 * v.Trace() * ds(bonus_intorder=intorder - 4)).Assemble()
    mass = BilinearForm(u.Trace() * v.Trace() * ds(bonus_intorder=intorder - 4)).Assemble().mat
    pre = mass.Inverse(freedofs=fes.FreeDofs())
    sigma = GridFunction(fes)
    inv = CGSolver(V.mat, pre, printrates=False, maxiter=maxiter, tol=tol)
    sigma.vec.data = inv * rhs.vec
    Q = Integrate(sigma, mesh, definedon=mesh.Boundaries(".*"))
    return {"ndof": fes.ndof, "capacitance": float(Q)}


def sphere_capacitance(R=1.0, maxh=0.4, order=3, **kw):
    """Capacitance of an isolated conducting sphere of radius R via the Laplace single-layer BEM, vs the
    analytic C = 4 pi eps0 R (eps0=1).  Returns ``{ndof, capacitance, analytic, rel_err}``."""
    mesh = ng.Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), R)).GenerateMesh(maxh=maxh)).Curve(kw.pop("curve_order", 4))
    res = single_layer_capacitance(mesh, order=order, **kw)
    res["analytic"] = 4.0 * np.pi * R
    res["rel_err"] = abs(res["capacitance"] - res["analytic"]) / res["analytic"]
    return res


def spheroid_capacitance_analytic(a_polar, b_equatorial):
    r"""The EXACT capacitance of an isolated conducting SPHEROID -- polar semi-axis ``a_polar`` (along z),
    equatorial semi-axis ``b_equatorial`` -- in the eps0=1 normalisation.  With the focal distance
    c_f = sqrt(|a^2 - b^2|):

        PROLATE (a > b):  C = 4 pi c_f / artanh(c_f / a),
        OBLATE  (a < b):  C = 4 pi c_f / arcsin(c_f / b),
        SPHERE  (a = b):  C = 4 pi a   (the smooth limit of both).

    The oblate disk limit (a -> 0) gives the conducting-disk capacitance C = 8 b; the prolate needle limit
    gives the slender-rod log.  Returns the capacitance (float)."""
    a = float(a_polar); b = float(b_equatorial)
    if a <= 0.0 or b <= 0.0:
        raise ValueError("semi-axes must be > 0 (got a=%r, b=%r)" % (a_polar, b_equatorial))
    if abs(a - b) < 1e-12 * max(a, b):
        return 4.0 * np.pi * a
    if a > b:
        cf = math.sqrt(a * a - b * b)
        return 4.0 * np.pi * cf / math.atanh(cf / a)
    cf = math.sqrt(b * b - a * a)
    return 4.0 * np.pi * cf / math.asin(cf / b)


def spheroid_capacitance(a_polar, b_equatorial, maxh=0.18, order=3, **kw):
    """Capacitance of an isolated conducting SPHEROID (polar semi-axis ``a_polar``, equatorial
    ``b_equatorial``) via the Laplace single-layer BEM, vs the exact closed form
    (:func:`spheroid_capacitance_analytic`).  Extends :func:`sphere_capacitance` to a non-spherical
    conductor -- the prolate/oblate generalisation of C = 4 pi eps0 R.  Returns
    ``{ndof, capacitance, analytic, rel_err, kind}``."""
    a = float(a_polar); b = float(b_equatorial)
    ax = occ.Axes(Pnt(0, 0, 0), occ.Z, occ.X)
    ell = occ.Ellipsoid(ax, b, b, a)                                # r1=r2=b (equatorial), r3=a (polar, z)
    mesh = ng.Mesh(OCCGeometry(ell).GenerateMesh(maxh=maxh)).Curve(kw.pop("curve_order", 4))
    res = single_layer_capacitance(mesh, order=order, **kw)
    res["analytic"] = spheroid_capacitance_analytic(a, b)
    res["rel_err"] = abs(res["capacitance"] - res["analytic"]) / res["analytic"]
    res["kind"] = "sphere" if abs(a - b) < 1e-12 * max(a, b) else ("prolate" if a > b else "oblate")
    return res


def conductor_polarizability(mesh, axis="z", order=3, intorder=12, maxiter=600):
    """Static polarizability (induced-dipole response to a UNIFORM applied field) of the
    conductor(s) bounded by ``mesh``, via the Laplace single-layer indirect BEM. A neutral
    conductor in a unit applied field along ``axis`` sits at potential 0, so its induced
    surface charge must produce a single-layer potential equal to that coordinate's trace
    on the surface:  V sigma = r_axis|_Gamma. The induced dipole component is then

        alpha = integral_Gamma sigma * r_axis dS         (eps0 = 1; alpha = p / E0).

    This is the applied-field / multipole-RESPONSE complement to
    :func:`single_layer_capacitance` (which is the fixed-potential self-capacitance, datum
    1 instead of a coordinate). ``qnet`` returns the net induced charge (~0 = neutral)."""
    coord = {"x": ng.x, "y": ng.y, "z": ng.z}[axis]
    fes = ng.SurfaceL2(mesh, order=order, dual_mapping=True)
    u, v = fes.TnT()
    V = bem.SingleLayerPotentialOperator(fes, intorder=intorder)
    rhs = LinearForm(coord * v.Trace() * ds(bonus_intorder=intorder - 4)).Assemble()
    mass = BilinearForm(u.Trace() * v.Trace() * ds(bonus_intorder=intorder - 4)).Assemble().mat
    pre = mass.Inverse(freedofs=fes.FreeDofs())
    sigma = GridFunction(fes)
    inv = CGSolver(V.mat, pre, printrates=False, maxiter=maxiter)
    sigma.vec.data = inv * rhs.vec
    bnd = mesh.Boundaries(".*")
    alpha = Integrate(sigma * coord, mesh, definedon=bnd)
    qnet = Integrate(sigma, mesh, definedon=bnd)
    return {"ndof": fes.ndof, "alpha": float(alpha), "qnet": float(qnet)}


def sphere_polarizability(R=1.0, maxh=0.3, order=3, curve_order=4, **kw):
    """Polarizability of an isolated conducting sphere via the Laplace single-layer BEM, vs
    the exact isotropic value  alpha = 4 pi eps0 R^3  (eps0 = 1; the induced surface charge is
    sigma = 3 eps0 E0 cos theta).  The applied-field-response twin of :func:`sphere_capacitance`.
    Returns ``{ndof, alpha, qnet, analytic, rel_err}``."""
    mesh = ng.Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), R)).GenerateMesh(maxh=maxh)).Curve(curve_order)
    res = conductor_polarizability(mesh, axis="z", order=order, **kw)
    res["analytic"] = 4.0 * np.pi * R**3
    res["rel_err"] = abs(res["alpha"] - res["analytic"]) / res["analytic"]
    return res


def helmholtz_exterior_point_source(kappa, a=1.0, maxh=0.4, order=3, rext=3.0, intorder=12):
    """HIGH-FREQUENCY Helmholtz single-layer BEM (wavenumber ``kappa``), verified by the MANUFACTURED
    SOLUTION of an interior point source -- the high-frequency (kappa>0) generalisation of the static
    (kappa=0) capacitance, the open counterpart of a commercial high-frequency integral (MoM) solver.

    A point source at the origin radiates the Helmholtz field u(x) = exp(i kappa |x|)/(4 pi |x|); on the
    sphere |x|=a its trace is the CONSTANT datum g = exp(i kappa a)/(4 pi a).  Solving the exterior
    Dirichlet problem with the Helmholtz single layer (V sigma = g; GMRES, V is complex symmetric) and
    evaluating the single-layer potential u_BEM = INT_Gamma sigma(y) exp(i kappa|x-y|)/(4 pi|x-y|) dGamma
    at an exterior point x=(0,0,rext) MUST reproduce the point-source field there (exterior-Dirichlet
    uniqueness).  Returns ``{ndof, kappa, ka, u_bem, u_exact, rel_err}`` (rel_err ~ machine precision).
    """
    from ngsolve import x, y, z, exp, sqrt
    from ngsolve.krylovspace import GMResSolver
    import cmath
    mesh = ng.Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), a)).GenerateMesh(maxh=maxh)).Curve(4)
    fes = ng.SurfaceL2(mesh, order=order, complex=True, dual_mapping=True)
    u, v = fes.TnT()
    V = bem.HelmholtzSingleLayerPotentialOperator(fes, fes, kappa, intorder=intorder)
    g = cmath.exp(1j * kappa * a) / (4.0 * np.pi * a)                      # constant interior-source trace
    rhs = LinearForm(g * v.Trace() * ds(bonus_intorder=intorder - 4)).Assemble()
    mass = BilinearForm(u.Trace() * v.Trace() * ds(bonus_intorder=intorder - 4)).Assemble().mat
    pre = mass.Inverse(freedofs=fes.FreeDofs())
    sigma = GridFunction(fes)
    inv = GMResSolver(V.mat, pre, printrates=False, maxiter=1000, tol=1e-10)
    sigma.vec.data = inv * rhs.vec
    R = sqrt(x * x + y * y + (z - rext) ** 2)                             # |x_ext - y|, x_ext=(0,0,rext)
    u_bem = Integrate(sigma * exp(1j * kappa * R) / (4.0 * np.pi * R) * ds(bonus_intorder=intorder - 4), mesh)
    u_exact = cmath.exp(1j * kappa * rext) / (4.0 * np.pi * rext)
    return {"ndof": fes.ndof, "kappa": kappa, "ka": kappa * a, "u_bem": complex(u_bem),
            "u_exact": u_exact, "rel_err": abs(u_bem - u_exact) / abs(u_exact)}


def mie_soundsoft_scattered(kappa, a, px, py, pz, nmax=40):
    r"""Exact MIE-SERIES scattered field of a SOUND-SOFT (Dirichlet) sphere of radius ``a`` at the point
    ``(px,py,pz)``, for the incident plane wave ``u_inc = exp(i kappa z)``.  Sound-soft means
    ``u_inc + u_scat = 0`` on the sphere, so the partial-wave coefficients are ``-j_n(ka)/h_n^{(1)}(ka)``:

        u_scat(r,theta) = -sum_n (2n+1) i^n [j_n(ka)/h_n(ka)] h_n(k r) P_n(cos theta),  h_n = j_n + i y_n.

    This is the canonical analytic gate for an exterior Helmholtz scattering solver (the scalar / acoustic
    Mie series).  ``j_n``/``y_n`` are spherical Bessel functions; the sum converges by ``n ~ ka + 10``."""
    from scipy.special import spherical_jn, spherical_yn, eval_legendre
    r = (px * px + py * py + pz * pz) ** 0.5
    costh = pz / r
    s = 0.0 + 0.0j
    for n in range(nmax + 1):
        hn_a = spherical_jn(n, kappa * a) + 1j * spherical_yn(n, kappa * a)
        hn_r = spherical_jn(n, kappa * r) + 1j * spherical_yn(n, kappa * r)
        s += (2 * n + 1) * (1j ** n) * (spherical_jn(n, kappa * a) / hn_a) * hn_r * eval_legendre(n, costh)
    return -s


def helmholtz_soundsoft_sphere_scattering(kappa, a=1.0, maxh=0.35, order=3, rext=2.0, intorder=12):
    r"""HIGH-FREQUENCY Helmholtz exterior SCATTERING by a SOUND-SOFT (Dirichlet) sphere via the ngsolve.bem
    single-layer BEM, verified against the exact MIE SERIES -- a TRUE scattering problem (incident plane
    wave on a real obstacle), the high-frequency scattering counterpart of a commercial integral (MoM)
    solver, beyond the #manufactured-solution check (helmholtz_exterior_point_source).

    Incident plane wave ``u_inc = exp(i kappa z)``; the sound-soft condition ``u_inc + u_scat = 0`` on the
    sphere gives the Dirichlet datum ``u_scat = -u_inc`` on Gamma.  The indirect single-layer ansatz
    ``u_scat(x) = INT_Gamma G_kappa(x,y) sigma(y) dGamma`` then solves the first-kind equation
    ``V sigma = -u_inc|_Gamma`` (GMRES; V complex symmetric).  Evaluating the single-layer potential at the
    exterior on-axis point ``(0,0,rext)`` and equatorial point ``(rext,0,0)`` must reproduce the Mie field.

    NB: the sound-soft single layer is singular at the interior Dirichlet eigenvalues ``j_n(ka)=0`` (lowest
    ``ka=pi``) -- this helper is for ``ka < pi``; a combined-field (Brakhage-Werner) formulation lifts that
    cap for higher ``ka``.  Returns ``{ndof, kappa, ka, axial:{u_bem,u_mie,rel_err}, equator:{...}}``.
    """
    from ngsolve.krylovspace import GMResSolver
    mesh = ng.Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), a)).GenerateMesh(maxh=maxh)).Curve(4)
    fes = ng.SurfaceL2(mesh, order=order, complex=True, dual_mapping=True)
    u, v = fes.TnT()
    V = bem.HelmholtzSingleLayerPotentialOperator(fes, fes, kappa, intorder=intorder)
    uinc = ng.exp(1j * kappa * ng.z)                                      # incident plane wave (+z)
    rhs = LinearForm(-uinc * v.Trace() * ds(bonus_intorder=intorder - 4)).Assemble()
    mass = BilinearForm(u.Trace() * v.Trace() * ds(bonus_intorder=intorder - 4)).Assemble().mat
    pre = mass.Inverse(freedofs=fes.FreeDofs())
    sigma = GridFunction(fes)
    inv = GMResSolver(V.mat, pre, printrates=False, maxiter=2000, tol=1e-10)
    sigma.vec.data = inv * rhs.vec
    res = {"ndof": fes.ndof, "kappa": kappa, "ka": kappa * a}
    for label, (px, py, pz) in [("axial", (0.0, 0.0, rext)), ("equator", (rext, 0.0, 0.0))]:
        R = ng.sqrt((ng.x - px) ** 2 + (ng.y - py) ** 2 + (ng.z - pz) ** 2)
        u_bem = complex(Integrate(sigma * ng.exp(1j * kappa * R) / (4.0 * np.pi * R)
                                  * ds(bonus_intorder=intorder - 4), mesh))
        u_mie = mie_soundsoft_scattered(kappa, a, px, py, pz)
        res[label] = {"u_bem": u_bem, "u_mie": u_mie, "rel_err": abs(u_bem - u_mie) / abs(u_mie)}
    return res


def mie_soundhard_scattered(kappa, a, px, py, pz, nmax=40):
    r"""Exact MIE-SERIES scattered field of a SOUND-HARD (rigid / Neumann) sphere of radius ``a`` at the
    point ``(px,py,pz)``, for the incident plane wave ``u_inc = exp(i kappa z)`` -- the NEUMANN DUAL of the
    sound-soft (Dirichlet) :func:`mie_soundsoft_scattered`.  Sound-hard means the TOTAL normal derivative
    vanishes on the sphere, ``d/dr (u_inc + u_scat) = 0`` at r=a, so the partial-wave coefficients carry the
    DERIVATIVE Riccati-Bessel ratio ``-j_n'(ka)/h_n^{(1)'}(ka)`` (vs ``-j_n(ka)/h_n(ka)`` for sound-soft):

        u_scat(r,theta) = sum_n (2n+1) i^n [-j_n'(ka)/h_n'(ka)] h_n(k r) P_n(cos theta),  h_n = j_n + i y_n.

    This is the canonical analytic gate for an exterior sound-hard (rigid) Helmholtz scattering solver (the
    rigid-sphere acoustic Mie series; the j_n'/h_n' primes are derivatives w.r.t. the argument)."""
    from scipy.special import spherical_jn, spherical_yn, eval_legendre
    r = (px * px + py * py + pz * pz) ** 0.5
    costh = pz / r
    s = 0.0 + 0.0j
    for n in range(nmax + 1):
        jnp = spherical_jn(n, kappa * a, derivative=True)
        ynp = spherical_yn(n, kappa * a, derivative=True)
        hnp_a = jnp + 1j * ynp
        hn_r = spherical_jn(n, kappa * r) + 1j * spherical_yn(n, kappa * r)
        s += (2 * n + 1) * (1j ** n) * (-jnp / hnp_a) * hn_r * eval_legendre(n, costh)
    return s


def helmholtz_soundhard_sphere_scattering(kappa, a=1.0, maxh=0.35, order=3, rext=2.0, intorder=12):
    r"""HIGH-FREQUENCY Helmholtz exterior SCATTERING by a SOUND-HARD (rigid / Neumann) sphere via the
    ngsolve.bem DIRECT (Calderon) BEM, verified against the exact rigid-sphere MIE series -- the NEUMANN DUAL
    of the sound-soft single-layer scattering (:func:`helmholtz_soundsoft_sphere_scattering`), the rigid
    acoustic obstacle a commercial high-frequency integral (MoM) solver targets.

    Incident plane wave ``u_inc = exp(i kappa z)``; the sound-hard (rigid) condition ``d/dn(u_inc+u_scat)=0``
    on Gamma gives the KNOWN Neumann datum ``sigma = du_scat/dn = -du_inc/dn = -i kappa n_z exp(i kappa z)``.
    The scattered field's Cauchy data on Gamma satisfy the EXTERIOR Calderon identity (the same NGSolve sign
    as :func:`laplace_exterior_dtn`): ``V sigma = (-1/2 M + K) u`` with K the double-layer, M the mass, u the
    unknown Dirichlet trace.  So the DIRECT second-kind solve is ``(K - 1/2 M) u = V sigma`` (GMRES; K, V the
    Helmholtz double/single-layer H-matrices).  The scattered field is reconstructed by the exterior
    representation ``u_scat(x) = INT_Gamma [dG/dn_y(x,y) u(y) - G(x,y) sigma(y)] dGamma`` (double-layer minus
    single-layer potential) and compared to the rigid Mie series at the exterior axial / equatorial points --
    to ~1e-7.

    NB: the direct (K - 1/2 M) operator is rank-deficient at the interior DIRICHLET eigenvalues (lowest
    ka=pi); this helper is for ka < pi (a Burton-Miller combined formulation lifts that cap).  Returns
    ``{ndof, kappa, ka, axial:{u_bem,u_mie,rel_err}, equator:{...}}``.
    """
    from ngsolve.krylovspace import GMResSolver
    from ngsolve import specialcf
    mesh = ng.Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), a)).GenerateMesh(maxh=maxh)).Curve(4)
    fes = ng.SurfaceL2(mesh, order=order, complex=True, dual_mapping=True)
    u, v = fes.TnT()
    V = bem.HelmholtzSingleLayerPotentialOperator(fes, fes, kappa, intorder=intorder)
    K = bem.HelmholtzDoubleLayerPotentialOperator(fes, fes, kappa, intorder=intorder)
    mass = BilinearForm(u.Trace() * v.Trace() * ds(bonus_intorder=intorder - 4)).Assemble().mat
    nrm = specialcf.normal(3)
    g = -1j * kappa * nrm[2] * ng.exp(1j * kappa * ng.z)                  # Neumann datum sigma = -du_inc/dn
    pre = mass.Inverse(freedofs=fes.FreeDofs())
    bg = LinearForm(g * v.Trace() * ds(bonus_intorder=intorder - 4)).Assemble()
    sigma = GridFunction(fes)
    sigma.vec.data = pre * bg.vec                                          # project the Neumann data to coeffs
    rhs = (V.mat * sigma.vec).Evaluate()                                  # V sigma  (single-layer H-matrix)
    uD = GridFunction(fes)                                                 # Dirichlet trace u_scat|_Gamma
    inv = GMResSolver(K.mat + (-0.5) * mass, pre, printrates=False, maxiter=3000, tol=1e-10)
    uD.vec.data = inv * rhs                                                # (K - 1/2 M) u = V sigma
    res = {"ndof": fes.ndof, "kappa": kappa, "ka": kappa * a}
    for label, (px, py, pz) in [("axial", (0.0, 0.0, rext)), ("equator", (rext, 0.0, 0.0))]:
        R = ng.sqrt((ng.x - px) ** 2 + (ng.y - py) ** 2 + (ng.z - pz) ** 2)
        E = ng.exp(1j * kappa * R)
        dGdn = (1j * kappa * R - 1.0) * E / (4.0 * np.pi * R ** 3) * \
            ((ng.x - px) * nrm[0] + (ng.y - py) * nrm[1] + (ng.z - pz) * nrm[2])
        SL = sigma * E / (4.0 * np.pi * R)                                # single-layer potential of sigma
        u_bem = complex(Integrate((uD * dGdn - SL) * ds(bonus_intorder=intorder - 4), mesh))  # DL(u) - SL(sigma)
        u_mie = mie_soundhard_scattered(kappa, a, px, py, pz)
        res[label] = {"u_bem": u_bem, "u_mie": u_mie, "rel_err": abs(u_bem - u_mie) / abs(u_mie)}
    return res


def sphere_radiation_ratio(kappa, a, n):
    r"""The exact surface-field / normal-derivative-datum ratio ``u/sigma = h_n(ka)/(k h_n^{(1)'}(ka))`` of a
    sphere driven in the n-th spherical RADIATION mode (n=0 monopole, n=1 dipole, ...) -- the analytic gate
    for the velocity-driven radiation impedance (:func:`helmholtz_sphere_radiation_impedance`), the RADIATION
    sibling of the scattering Mie ratio.  The monopole reduces to ``a/(i ka - 1)``; the specific radiation
    impedance ``z = j omega rho (u/sigma)`` has the famous resistance ``Re(z) = rho c (ka)^2/(1+(ka)^2)``
    (monopole).  ``h_n = j_n + i y_n``; primes are derivatives w.r.t. the argument.  Returns complex."""
    from scipy.special import spherical_jn, spherical_yn
    x = kappa * a
    hn = spherical_jn(n, x) + 1j * spherical_yn(n, x)
    hnp = spherical_jn(n, x, derivative=True) + 1j * spherical_yn(n, x, derivative=True)
    return hn / (kappa * hnp)


def helmholtz_sphere_radiation_impedance(kappa, a=1.0, maxh=0.3, order=3, intorder=12):
    r"""The acoustic RADIATION impedance of a sphere driven by a prescribed surface NORMAL-VELOCITY datum, via
    the ngsolve.bem velocity-driven (Neumann) direct Calderon BEM -- the RADIATION complement of the exterior
    SCATTERING family (sound-soft / sound-hard Mie); the open counterpart of an antenna / transducer
    radiation solve.

    A radiating sphere with normal-derivative datum ``sigma = du/dn`` on Gamma (the source) has its surface
    field ``u|_Gamma`` related to sigma by the SAME exterior Calderon identity as the sound-hard scattering
    (:func:`helmholtz_soundhard_sphere_scattering`): ``(K - 1/2 M) u = V sigma`` (GMRES; K, V the Helmholtz
    double / single-layer H-matrices).  Driving the MONOPOLE (``sigma = 1``, uniform) and DIPOLE
    (``sigma = cos theta = z/a``) modes and extracting the matching surface coefficient gives the radiation
    ratio ``u/sigma``, compared to the exact ``h_n(ka)/(k h_n'(ka))`` (:func:`sphere_radiation_ratio`) -- to
    ~1e-7.  The monopole real part is the radiation resistance ``R/(rho c) = (ka)^2/(1+(ka)^2)``.  Returns
    ``{ndof, kappa, ka, monopole:{ratio_bem,ratio_analytic,rel_err}, dipole:{...}}``."""
    from ngsolve.krylovspace import GMResSolver
    mesh = ng.Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), a)).GenerateMesh(maxh=maxh)).Curve(4)
    fes = ng.SurfaceL2(mesh, order=order, complex=True, dual_mapping=True)
    u, v = fes.TnT()
    V = bem.HelmholtzSingleLayerPotentialOperator(fes, fes, kappa, intorder=intorder)
    K = bem.HelmholtzDoubleLayerPotentialOperator(fes, fes, kappa, intorder=intorder)
    mass = BilinearForm(u.Trace() * v.Trace() * ds(bonus_intorder=intorder - 4)).Assemble().mat
    pre = mass.Inverse(freedofs=fes.FreeDofs())
    area = Integrate(ng.CoefficientFunction(1.0) * ds, mesh).real
    res = {"ndof": fes.ndof, "kappa": kappa, "ka": kappa * a}
    for n, sig_cf, label in [(0, ng.CoefficientFunction(1.0 + 0j), "monopole"),
                             (1, ng.CoefficientFunction(ng.z / a + 0j), "dipole")]:
        bg = LinearForm(sig_cf * v.Trace() * ds(bonus_intorder=intorder - 4)).Assemble()
        sigvec = GridFunction(fes)
        sigvec.vec.data = pre * bg.vec                                # project the Neumann datum to coeffs
        rhs = (V.mat * sigvec.vec).Evaluate()
        uD = GridFunction(fes)
        inv = GMResSolver(K.mat + (-0.5) * mass, pre, printrates=False, maxiter=3000, tol=1e-10)
        uD.vec.data = inv * rhs                                       # (K - 1/2 M) u = V sigma
        if n == 0:
            ratio_bem = complex(Integrate(uD * ds, mesh)) / area      # <u> / <sigma=1>
        else:
            num = complex(Integrate(uD * (ng.z / a) * ds, mesh))      # <u cos> / <cos^2>
            den = complex(Integrate((ng.z / a) ** 2 * ds, mesh))
            ratio_bem = num / den
        an = sphere_radiation_ratio(kappa, a, n)
        res[label] = {"ratio_bem": ratio_bem, "ratio_analytic": an, "rel_err": abs(ratio_bem - an) / abs(an)}
    return res


def helmholtz_cfie_soundsoft_sphere(kappa, a=1.0, maxh=0.35, order=3, rext=2.0, intorder=12):
    r"""COMBINED-FIELD (Brakhage-Werner) sound-soft sphere scattering via ngsolve.bem -- the IRREGULAR-
    FREQUENCY-ROBUST formulation that fixes the single-layer's ka<pi limitation
    (helmholtz_soundsoft_sphere_scattering).  The indirect ansatz couples the double- and single-layer
    potentials, u_scat = (D - i kappa S) psi, which is UNIQUELY solvable for ALL real kappa (no spurious
    interior resonances).  The boundary integral equation is (1/2 I + K - i kappa V) psi = -u_inc|_Gamma,
    assembled here as ``(1/2 M + HelmholtzCombinedFieldOperator) psi = rhs`` (the operator excludes the 1/2
    jump, added via the mass matrix M).  The scattered field is evaluated as

        u_scat(x) = INT_Gamma psi(y) [ dG(x,y)/dn_y - i kappa G(x,y) ] dGamma(y),

    with dG/dn_y = (i kappa R - 1) exp(i kappa R)/(4 pi R^3) (y-x).n(y).  Verified against the scalar MIE
    series at the exterior axial and equatorial points -- to ~machine precision (~1e-6) for ALL ka INCLUDING
    ka=pi (where the single-layer's operator is rank-deficient; see helmholtz_irregular_frequency_condition).
    Returns ``{ndof, kappa, ka, axial:{u_bem,u_mie,rel_err}, equator:{...}}``.
    """
    from ngsolve.krylovspace import GMResSolver
    from ngsolve import specialcf
    mesh = ng.Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), a)).GenerateMesh(maxh=maxh)).Curve(4)
    fes = ng.SurfaceL2(mesh, order=order, complex=True, dual_mapping=True)
    u, v = fes.TnT()
    cfo = bem.HelmholtzCombinedFieldOperator(fes, fes, kappa=kappa, intorder=intorder)
    mass = BilinearForm(u.Trace() * v.Trace() * ds(bonus_intorder=intorder - 4)).Assemble().mat
    rhs = LinearForm(-ng.exp(1j * kappa * ng.z) * v.Trace() * ds(bonus_intorder=intorder - 4)).Assemble()
    psi = GridFunction(fes)
    pre = mass.Inverse(freedofs=fes.FreeDofs())
    inv = GMResSolver(cfo.mat + 0.5 * mass, pre, printrates=False, maxiter=3000, tol=1e-10)
    psi.vec.data = inv * rhs.vec
    nrm = specialcf.normal(3)
    res = {"ndof": fes.ndof, "kappa": kappa, "ka": kappa * a}
    for label, (px, py, pz) in [("axial", (0.0, 0.0, rext)), ("equator", (rext, 0.0, 0.0))]:
        R = ng.sqrt((ng.x - px) ** 2 + (ng.y - py) ** 2 + (ng.z - pz) ** 2)
        E = ng.exp(1j * kappa * R)
        dGdn = (1j * kappa * R - 1.0) * E / (4.0 * np.pi * R ** 3) * \
            ((ng.x - px) * nrm[0] + (ng.y - py) * nrm[1] + (ng.z - pz) * nrm[2])
        G = E / (4.0 * np.pi * R)
        u_bem = complex(Integrate(psi * (dGdn - 1j * kappa * G) * ds(bonus_intorder=intorder - 4), mesh))
        u_mie = mie_soundsoft_scattered(kappa, a, px, py, pz)
        res[label] = {"u_bem": u_bem, "u_mie": u_mie, "rel_err": abs(u_bem - u_mie) / abs(u_mie)}
    return res


def helmholtz_irregular_frequency_condition(kappa, a=1.0, maxh=0.9, order=0, intorder=10):
    r"""Condition numbers of the single-layer V vs the combined-field (1/2 M + CFO) operator at wavenumber
    ``kappa`` -- the rigorous "why the combined field" demonstration.  At the interior Dirichlet eigenvalues
    (lowest ka=pi, where j_0(ka)=0) the single-layer V is RANK-DEFICIENT (cond -> inf in the continuum), so
    cond(V) SPIKES; the Brakhage-Werner combined operator is uniquely solvable for all kappa, so its
    condition number stays bounded and flat.  Uses a coarse, low-order surface space so the H-matrix
    operators can be densified (applied to each unit vector) and fed to ``numpy.linalg.cond``.
    Returns ``{ndof, kappa, ka, cond_single_layer, cond_combined_field}``."""
    mesh = ng.Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), a)).GenerateMesh(maxh=maxh)).Curve(3)
    fes = ng.SurfaceL2(mesh, order=order, complex=True, dual_mapping=True)
    u, v = fes.TnT()
    nd = fes.ndof
    V = bem.HelmholtzSingleLayerPotentialOperator(fes, fes, kappa, intorder=intorder)
    cfo = bem.HelmholtzCombinedFieldOperator(fes, fes, kappa=kappa, intorder=intorder)
    mass = BilinearForm(u.Trace() * v.Trace() * ds(bonus_intorder=intorder - 2)).Assemble().mat

    def _dense(op):
        x = op.CreateColVector(); y = op.CreateColVector()
        M = np.zeros((nd, nd), dtype=complex)
        for j in range(nd):
            x[:] = 0.0; x[j] = 1.0
            y.data = op * x
            M[:, j] = y.FV().NumPy()[:]
        return M
    Vd = _dense(V.mat)
    Cd = _dense(cfo.mat) + 0.5 * _dense(mass)
    return {"ndof": nd, "kappa": kappa, "ka": kappa * a,
            "cond_single_layer": float(np.linalg.cond(Vd)),
            "cond_combined_field": float(np.linalg.cond(Cd))}


def mie_soundsoft_farfield(kappa, a, theta, nmax=40):
    r"""Exact scalar MIE far-field scattering amplitude f(theta) of a sound-soft sphere (incident
    exp(i kappa z)): the scattered field is ``u_scat ~ f(theta) exp(i k r)/r`` as r->inf, with

        f(theta) = (i/kappa) sum_n (2n+1) [j_n(ka)/h_n(ka)] P_n(cos theta).

    The analytic gate for the RCS / far-field of an exterior Helmholtz solver."""
    from scipy.special import spherical_jn, spherical_yn, eval_legendre
    import cmath
    ct = math.cos(theta); s = 0.0 + 0.0j
    for n in range(nmax + 1):
        hn = spherical_jn(n, kappa * a) + 1j * spherical_yn(n, kappa * a)
        s += (2 * n + 1) * (spherical_jn(n, kappa * a) / hn) * eval_legendre(n, ct)
    return 1j * s / kappa


def mie_soundsoft_cross_section(kappa, a, nmax=40):
    r"""Total scattering cross-section of a sound-soft sphere, ``sigma_sca = (4 pi/k^2) sum (2n+1)
    |j_n(ka)/h_n(ka)|^2`` -- the analytic gate for the optical theorem (for a lossless scatterer
    sigma_t = sigma_sca)."""
    from scipy.special import spherical_jn, spherical_yn
    s = 0.0
    for n in range(nmax + 1):
        hn = spherical_jn(n, kappa * a) + 1j * spherical_yn(n, kappa * a)
        s += (2 * n + 1) * abs(spherical_jn(n, kappa * a) / hn) ** 2
    return 4.0 * math.pi / kappa ** 2 * s


def helmholtz_soundsoft_far_field(kappa, a=1.0, maxh=0.3, order=3, intorder=12,
                                  angles_deg=(0.0, 45.0, 90.0, 135.0, 180.0)):
    r"""FAR-FIELD scattering amplitude f(theta) + the OPTICAL THEOREM for the sound-soft sphere, on the
    combined-field (Brakhage-Werner) density -- the RCS apparatus of the open ngsolve.bem integral solver
    (the open counterpart of a commercial MoM solver's RCS).  The EXACT far-field of the combined potential
    ``u_scat = INT psi[dG/dn - i k G]`` is

        f(xhat) = (-i k/4 pi) INT_Gamma psi(y) [(xhat . n_y) + 1] exp(-i k xhat.y) dGamma(y)

    (no 1/r truncation -- a surface integral against the far-field kernel).  Verified vs the scalar Mie
    far-field ``mie_soundsoft_farfield`` to ~machine precision, AND the OPTICAL THEOREM
    ``sigma_t = (4 pi/k) Im f(0) == sigma_sca`` (Mie cross-section) -- a forward-scattering self-consistency
    gate independent of the angular comparison.  ``angles_deg`` are the scattering angles theta.
    Returns ``{ndof, kappa, ka, angles_deg, f_bem, f_mie, rel_err, sigma_optical, sigma_mie, optical_rel_err}``.
    """
    from ngsolve.krylovspace import GMResSolver
    from ngsolve import specialcf
    mesh = ng.Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), a)).GenerateMesh(maxh=maxh)).Curve(4)
    fes = ng.SurfaceL2(mesh, order=order, complex=True, dual_mapping=True)
    u, v = fes.TnT()
    cfo = bem.HelmholtzCombinedFieldOperator(fes, fes, kappa=kappa, intorder=intorder)
    mass = BilinearForm(u.Trace() * v.Trace() * ds(bonus_intorder=intorder - 4)).Assemble().mat
    rhs = LinearForm(-ng.exp(1j * kappa * ng.z) * v.Trace() * ds(bonus_intorder=intorder - 4)).Assemble()
    psi = GridFunction(fes)
    pre = mass.Inverse(freedofs=fes.FreeDofs())
    inv = GMResSolver(cfo.mat + 0.5 * mass, pre, printrates=False, maxiter=3000, tol=1e-10)
    psi.vec.data = inv * rhs.vec
    nrm = specialcf.normal(3)

    def farfield(thx, thz):
        xn = thx * nrm[0] + thz * nrm[2]
        phase = ng.exp(-1j * kappa * (thx * ng.x + thz * ng.z))
        return complex((-1j * kappa / (4.0 * np.pi))
                       * Integrate(psi * (xn + 1.0) * phase * ds(bonus_intorder=intorder - 4), mesh))
    f_bem, f_mie, rel = [], [], []
    for thd in angles_deg:
        th = math.radians(thd)
        fb = farfield(math.sin(th), math.cos(th)); fm = mie_soundsoft_farfield(kappa, a, th)
        f_bem.append(fb); f_mie.append(fm); rel.append(abs(fb - fm) / abs(fm))
    f0 = farfield(0.0, 1.0)
    sig_opt = 4.0 * math.pi / kappa * f0.imag
    sig_mie = mie_soundsoft_cross_section(kappa, a)
    return {"ndof": fes.ndof, "kappa": kappa, "ka": kappa * a, "angles_deg": list(angles_deg),
            "f_bem": f_bem, "f_mie": f_mie, "rel_err": rel, "sigma_optical": sig_opt,
            "sigma_mie": sig_mie, "optical_rel_err": abs(sig_opt - sig_mie) / sig_mie}


def mie_soundhard_farfield(kappa, a, theta, nmax=40):
    r"""Exact scalar MIE far-field scattering amplitude f(theta) of a SOUND-HARD (rigid / Neumann) sphere
    (incident ``exp(i kappa z)``) -- the NEUMANN DUAL of :func:`mie_soundsoft_farfield`: the ``j_n/h_n`` ratio
    is replaced by the DERIVATIVE ratio ``j_n'(ka)/h_n'(ka)``,

        f(theta) = (i/kappa) sum_n (2n+1) [j_n'(ka)/h_n'(ka)] P_n(cos theta).

    The analytic gate for the rigid-body RCS / far-field of an exterior Helmholtz solver."""
    from scipy.special import spherical_jn, spherical_yn, eval_legendre
    ct = math.cos(theta); s = 0.0 + 0.0j
    for n in range(nmax + 1):
        jnp = spherical_jn(n, kappa * a, derivative=True)
        ynp = spherical_yn(n, kappa * a, derivative=True)
        hnp = jnp + 1j * ynp
        s += (2 * n + 1) * (jnp / hnp) * eval_legendre(n, ct)
    return 1j * s / kappa


def mie_soundhard_cross_section(kappa, a, nmax=40):
    r"""Total scattering cross-section of a SOUND-HARD (rigid) sphere, ``sigma_sca = (4 pi/k^2) sum (2n+1)
    |j_n'(ka)/h_n'(ka)|^2`` -- the rigid-body analytic gate for the optical theorem (lossless: sigma_t =
    sigma_sca)."""
    from scipy.special import spherical_jn, spherical_yn
    s = 0.0
    for n in range(nmax + 1):
        jnp = spherical_jn(n, kappa * a, derivative=True)
        ynp = spherical_yn(n, kappa * a, derivative=True)
        hnp = jnp + 1j * ynp
        s += (2 * n + 1) * abs(jnp / hnp) ** 2
    return 4.0 * math.pi / kappa ** 2 * s


def helmholtz_soundhard_far_field(kappa, a=1.0, maxh=0.3, order=3, intorder=12,
                                  angles_deg=(0.0, 45.0, 90.0, 135.0, 180.0)):
    r"""FAR-FIELD scattering amplitude f(theta) + the OPTICAL THEOREM for the SOUND-HARD (rigid / Neumann)
    sphere via the ngsolve.bem DIRECT (Calderon) solve -- the RIGID-body RCS apparatus of the open
    ngsolve.bem integral solver (the open counterpart of a commercial high-frequency MoM solver's RCS), the
    NEUMANN DUAL of :func:`helmholtz_soundsoft_far_field`.

    The direct second-kind solve ``(K - 1/2 M) u = V sigma`` (sigma = -du_inc/dn the known Neumann datum, u
    the unknown Dirichlet trace, as in :func:`helmholtz_soundhard_sphere_scattering`) gives the scattered
    Cauchy data; the FAR-FIELD of the exterior representation ``u_scat = INT[u dG/dn_y - sigma G]`` is

        f(xhat) = (1/4 pi) INT_Gamma [ -i k (xhat . n_y) u(y) - sigma(y) ] exp(-i k xhat.y) dGamma(y)

    (no 1/r truncation -- a surface integral against the far-field kernel).  Verified vs the rigid Mie
    far-field :func:`mie_soundhard_farfield` AND the OPTICAL THEOREM ``sigma_t = (4 pi/k) Im f(0) ==
    sigma_sca`` (rigid Mie cross-section).  ka < pi (the direct operator's interior-Dirichlet eigenvalue cap).
    Returns ``{ndof, kappa, ka, angles_deg, f_bem, f_mie, rel_err, sigma_optical, sigma_mie,
    optical_rel_err}``.
    """
    from ngsolve.krylovspace import GMResSolver
    from ngsolve import specialcf
    mesh = ng.Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), a)).GenerateMesh(maxh=maxh)).Curve(4)
    fes = ng.SurfaceL2(mesh, order=order, complex=True, dual_mapping=True)
    u, v = fes.TnT()
    V = bem.HelmholtzSingleLayerPotentialOperator(fes, fes, kappa, intorder=intorder)
    K = bem.HelmholtzDoubleLayerPotentialOperator(fes, fes, kappa, intorder=intorder)
    mass = BilinearForm(u.Trace() * v.Trace() * ds(bonus_intorder=intorder - 4)).Assemble().mat
    nrm = specialcf.normal(3)
    g = -1j * kappa * nrm[2] * ng.exp(1j * kappa * ng.z)                  # Neumann datum sigma = -du_inc/dn
    pre = mass.Inverse(freedofs=fes.FreeDofs())
    bg = LinearForm(g * v.Trace() * ds(bonus_intorder=intorder - 4)).Assemble()
    sigma = GridFunction(fes)
    sigma.vec.data = pre * bg.vec                                          # project the Neumann data to coeffs
    rhs = (V.mat * sigma.vec).Evaluate()                                  # V sigma
    uD = GridFunction(fes)                                                 # Dirichlet trace u_scat|_Gamma
    inv = GMResSolver(K.mat + (-0.5) * mass, pre, printrates=False, maxiter=3000, tol=1e-10)
    uD.vec.data = inv * rhs                                                # (K - 1/2 M) u = V sigma

    def farfield(thx, thz):
        xn = thx * nrm[0] + thz * nrm[2]
        phase = ng.exp(-1j * kappa * (thx * ng.x + thz * ng.z))
        return complex((1.0 / (4.0 * np.pi))
                       * Integrate(((-1j * kappa * xn) * uD - sigma) * phase * ds(bonus_intorder=intorder - 4), mesh))
    f_bem, f_mie, rel = [], [], []
    for thd in angles_deg:
        th = math.radians(thd)
        fb = farfield(math.sin(th), math.cos(th)); fm = mie_soundhard_farfield(kappa, a, th)
        f_bem.append(fb); f_mie.append(fm); rel.append(abs(fb - fm) / abs(fm))
    f0 = farfield(0.0, 1.0)
    sig_opt = 4.0 * math.pi / kappa * f0.imag
    sig_mie = mie_soundhard_cross_section(kappa, a)
    return {"ndof": fes.ndof, "kappa": kappa, "ka": kappa * a, "angles_deg": list(angles_deg),
            "f_bem": f_bem, "f_mie": f_mie, "rel_err": rel, "sigma_optical": sig_opt,
            "sigma_mie": sig_mie, "optical_rel_err": abs(sig_opt - sig_mie) / sig_mie}


def mie_impedance_farfield(kappa, a, beta, theta, nmax=40):
    r"""Exact scalar MIE far-field scattering amplitude f(theta) of an IMPEDANCE (Robin) sphere of radius
    ``a`` with the surface condition ``du/dn + i kappa beta u = 0`` (normalized admittance ``beta``) -- the
    INTERPOLATION between the sound-soft (Dirichlet, beta -> inf) and sound-hard (Neumann, beta = 0)
    spheres.  Applying the Robin BC to each partial wave gives the scattering ratio
    ``s_n = (j_n'(ka) + i beta j_n(ka)) / (h_n'(ka) + i beta h_n(ka))`` (primes = derivatives w.r.t. the
    argument), so

        f(theta) = (i/kappa) sum_n (2n+1) s_n P_n(cos theta).

    beta = 0 recovers :func:`mie_soundhard_farfield` (j_n'/h_n'); beta -> inf recovers
    :func:`mie_soundsoft_farfield` (j_n/h_n).  The analytic gate for an impedance / coated-obstacle
    integral-equation solve."""
    from scipy.special import spherical_jn, spherical_yn, eval_legendre
    ct = math.cos(theta); s = 0.0 + 0.0j
    for n in range(nmax + 1):
        jn = spherical_jn(n, kappa * a); yn = spherical_yn(n, kappa * a)
        jnp = spherical_jn(n, kappa * a, derivative=True); ynp = spherical_yn(n, kappa * a, derivative=True)
        hn = jn + 1j * yn; hnp = jnp + 1j * ynp
        s += (2 * n + 1) * ((jnp + 1j * beta * jn) / (hnp + 1j * beta * hn)) * eval_legendre(n, ct)
    return 1j * s / kappa


def helmholtz_impedance_far_field(kappa, a=1.0, beta=0.5, maxh=0.3, order=3, intorder=12,
                                  angles_deg=(0.0, 45.0, 90.0, 135.0, 180.0)):
    r"""FAR-FIELD scattering amplitude f(theta) of an IMPEDANCE (Robin) sphere ``du/dn + i kappa beta u = 0``
    via the ngsolve.bem DIRECT (Calderon) solve -- the impedance / coated-obstacle member of the open
    ngsolve.bem scattering family, INTERPOLATING the sound-soft (:func:`helmholtz_soundsoft_far_field`) and
    sound-hard (:func:`helmholtz_soundhard_far_field`) limits.

    The Robin BC on the scattered field is ``sigma + i kappa beta u = g`` with ``sigma = du_scat/dn``, ``u =
    u_scat`` and the known datum ``g = -(du_inc/dn + i kappa beta u_inc) = -i kappa (n_z + beta) exp(i kappa
    z)``.  Substituting ``sigma = g - i kappa beta u`` into the EXTERIOR Calderon identity ``V sigma =
    (-1/2 M + K) u`` (the same NGSolve sign as :func:`helmholtz_soundhard_sphere_scattering`) gives the
    single second-kind solve

        (K - 1/2 M + i kappa beta V) u = V g,

    reusing the SAME V (single-layer), K (double-layer), M (mass) operators as the sound-hard solve (just the
    added ``i kappa beta V`` term and the Robin datum).  The scattered field is reconstructed by ``u_scat =
    INT[u dG/dn_y - sigma G]`` with ``sigma = g - i kappa beta u`` and its far-field is

        f(xhat) = (1/4 pi) INT_Gamma [ -i k (xhat . n_y) u - sigma ] exp(-i k xhat.y) dGamma.

    Verified vs the impedance Mie far-field :func:`mie_impedance_farfield`.  For a finite (radiating) beta the
    impedance is dissipative, which regularises the operator (no interior-resonance cap).  Returns
    ``{ndof, kappa, ka, beta, angles_deg, f_bem, f_mie, rel_err}``.
    """
    from ngsolve.krylovspace import GMResSolver
    from ngsolve import specialcf
    mesh = ng.Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), a)).GenerateMesh(maxh=maxh)).Curve(4)
    fes = ng.SurfaceL2(mesh, order=order, complex=True, dual_mapping=True)
    u, v = fes.TnT()
    V = bem.HelmholtzSingleLayerPotentialOperator(fes, fes, kappa, intorder=intorder)
    K = bem.HelmholtzDoubleLayerPotentialOperator(fes, fes, kappa, intorder=intorder)
    mass = BilinearForm(u.Trace() * v.Trace() * ds(bonus_intorder=intorder - 4)).Assemble().mat
    nrm = specialcf.normal(3)
    g = -1j * kappa * (nrm[2] + beta) * ng.exp(1j * kappa * ng.z)         # Robin datum -(dn u_inc + i k beta u_inc)
    pre = mass.Inverse(freedofs=fes.FreeDofs())
    bg = LinearForm(g * v.Trace() * ds(bonus_intorder=intorder - 4)).Assemble()
    gproj = GridFunction(fes)
    gproj.vec.data = pre * bg.vec                                          # L2 projection of the Robin datum g
    rhs = (V.mat * gproj.vec).Evaluate()                                  # V g
    uD = GridFunction(fes)
    op = K.mat + (-0.5) * mass + (1j * kappa * beta) * V.mat               # K - 1/2 M + i k beta V
    inv = GMResSolver(op, pre, printrates=False, maxiter=3000, tol=1e-10)
    uD.vec.data = inv * rhs
    sigma = gproj - (1j * kappa * beta) * uD                               # sigma = g - i k beta u  (CF combo)

    def farfield(thx, thz):
        xn = thx * nrm[0] + thz * nrm[2]
        phase = ng.exp(-1j * kappa * (thx * ng.x + thz * ng.z))
        return complex((1.0 / (4.0 * np.pi))
                       * Integrate(((-1j * kappa * xn) * uD - sigma) * phase * ds(bonus_intorder=intorder - 4), mesh))
    f_bem, f_mie, rel = [], [], []
    for thd in angles_deg:
        th = math.radians(thd)
        fb = farfield(math.sin(th), math.cos(th)); fm = mie_impedance_farfield(kappa, a, beta, th)
        f_bem.append(fb); f_mie.append(fm); rel.append(abs(fb - fm) / abs(fm))
    return {"ndof": fes.ndof, "kappa": kappa, "ka": kappa * a, "beta": beta,
            "angles_deg": list(angles_deg), "f_bem": f_bem, "f_mie": f_mie, "rel_err": rel}


# ---------------------------------------------------------------------------
# ELECTROMAGNETIC (vector Maxwell) Mie -- analytic RCS gate for a future
# MaxwellSingleLayer EFIE/CFIE PEC-sphere solve.  The vector counterpart of the
# scalar/acoustic Mie helpers above; this is the canonical analytic gate a
# commercial HF integral-equation (MoM) RCS solver is checked against.
# ---------------------------------------------------------------------------

def mie_pec_coefficients(ka, nmax=None):
    r"""Electromagnetic Mie scattering coefficients of a PERFECTLY CONDUCTING (PEC) sphere, argument
    ``ka`` (size parameter, k = 2 pi / lambda, a = radius).  Using the Riccati-Bessel functions
    ``psi_n(x) = x j_n(x)`` and ``xi_n(x) = x h_n^{(1)}(x)``, the PEC limit of the Bohren-Huffman
    coefficients is

        a_n = psi_n'(ka) / xi_n'(ka)     (TM / electric multipole),
        b_n = psi_n(ka)  / xi_n(ka)      (TE / magnetic multipole).

    Returns ``(n, a_n, b_n)`` numpy arrays, n = 1 .. nmax.  ``nmax`` defaults to the Wiscombe
    rule ``ka + 4 ka^{1/3} + 10`` (well-converged for ka up to ~100)."""
    from scipy.special import spherical_jn, spherical_yn
    x = float(ka)
    if nmax is None:
        nmax = int(x + 4.0 * x ** (1.0 / 3.0) + 10.0)
    n = np.arange(1, nmax + 1)
    jn = spherical_jn(n, x); jnp = spherical_jn(n, x, derivative=True)
    yn = spherical_yn(n, x); ynp = spherical_yn(n, x, derivative=True)
    hn = jn + 1j * yn; hnp = jnp + 1j * ynp
    psi = x * jn; psip = jn + x * jnp                 # psi_n = x j_n, psi_n' = j_n + x j_n'
    xi = x * hn; xip = hn + x * hnp                   # xi_n  = x h_n, xi_n'  = h_n + x h_n'
    a_n = psip / xip                                  # TM (electric)
    b_n = psi / xi                                    # TE (magnetic)
    return n, a_n, b_n


def mie_pec_sphere_rcs(ka, nmax=None):
    r"""Closed-form ELECTROMAGNETIC radar cross-sections of a PEC sphere, normalised by the geometric
    cross-section ``pi a^2``.  From the Mie coefficients (:func:`mie_pec_coefficients`):

        total scattering   sigma_sca = (2 pi / k^2) sum_{n>=1} (2n+1) (|a_n|^2 + |b_n|^2),
        backscatter (RCS)  sigma_b   = (pi / k^2) | sum_{n>=1} (2n+1) (-1)^n (a_n - b_n) |^2.

    Two analytic limits gate the result (verified in the test suite):
      * RAYLEIGH (ka -> 0):     sigma_b / (pi a^2) -> 9 (ka)^4   (electric + magnetic dipole).
      * GEOMETRIC (ka -> inf):  sigma_b / (pi a^2) -> 1  (specular)  and
                                sigma_sca / (pi a^2) -> 2  (extinction paradox).

    Returns ``{ka, nmax, sigma_sca_over_pa2, sigma_back_over_pa2}`` (both dimensionless, i.e. divided by
    pi a^2; multiply by pi a^2 for absolute cross-sections)."""
    x = float(ka)
    n, a_n, b_n = mie_pec_coefficients(x, nmax=nmax)
    sca = (2.0 * np.pi / x ** 2) * np.sum((2 * n + 1) * (np.abs(a_n) ** 2 + np.abs(b_n) ** 2))
    back = (np.pi / x ** 2) * np.abs(np.sum((2 * n + 1) * (-1.0) ** n * (a_n - b_n))) ** 2
    return {"ka": x, "nmax": int(n[-1]),
            "sigma_sca_over_pa2": float(sca / np.pi),
            "sigma_back_over_pa2": float(back / np.pi)}


def _mie_angular_functions(mu, nmax):
    """Bohren-Huffman angular functions pi_n(mu), tau_n(mu), n=1..nmax (mu = cos theta)."""
    pi = np.zeros(nmax + 1); tau = np.zeros(nmax + 1)
    pi[1] = 1.0; tau[1] = mu
    for n in range(2, nmax + 1):
        pi[n] = ((2 * n - 1) / (n - 1)) * mu * pi[n - 1] - (n / (n - 1)) * pi[n - 2]
        tau[n] = n * mu * pi[n] - (n + 1) * pi[n - 1]
    return pi[1:], tau[1:]


def mie_pec_amplitude(ka, theta, nmax=None):
    r"""Complex scattering AMPLITUDE functions of a PEC sphere -- S1(theta) (perpendicular / H-plane)
    and S2(theta) (parallel / E-plane) -- from the Mie coefficients (:func:`mie_pec_coefficients`) and
    the Bohren-Huffman angular functions pi_n, tau_n:

        S1 = sum_{n>=1} (2n+1)/(n(n+1)) [a_n pi_n + b_n tau_n],
        S2 = sum_{n>=1} (2n+1)/(n(n+1)) [a_n tau_n + b_n pi_n].

    This is the full BISTATIC angular scattering pattern (forward theta=0 -> backscatter theta=pi), the
    vector-EM counterpart of the SCALAR :func:`mie_soundsoft_farfield`; the existing EM Mie helpers
    gave only the total and backscatter cross-sections. Returns ``(S1, S2)`` (complex)."""
    n, a_n, b_n = mie_pec_coefficients(ka, nmax=nmax)
    pin, taun = _mie_angular_functions(math.cos(theta), len(n))
    fac = (2 * n + 1) / (n * (n + 1))
    S1 = complex(np.sum(fac * (a_n * pin + b_n * taun)))
    S2 = complex(np.sum(fac * (a_n * taun + b_n * pin)))
    return S1, S2


def mie_pec_bistatic_rcs(ka, theta, nmax=None):
    r"""Bistatic radar cross-section sigma(theta)/(pi a^2) of a PEC sphere for the two principal planes,
    sigma = (4 pi/k^2)|S|^2  ->  /(pi a^2) = (4/x^2)|S|^2  with x = ka:  the E-plane (incident E in the
    scattering plane) uses S2, the H-plane uses S1.  At theta=pi both reduce to the monostatic backscatter
    (:func:`mie_pec_sphere_rcs`); the forward (theta=0) amplitude obeys the OPTICAL THEOREM
    ``Q_ext = (4/x^2) Re[S(0)] = Q_sca`` (lossless PEC).  Returns
    ``{theta, e_plane_over_pa2, h_plane_over_pa2}``."""
    x = float(ka)
    S1, S2 = mie_pec_amplitude(x, theta, nmax=nmax)
    return {"theta": float(theta),
            "e_plane_over_pa2": (4.0 / x ** 2) * abs(S2) ** 2,
            "h_plane_over_pa2": (4.0 / x ** 2) * abs(S1) ** 2}


def maxwell_efie_pec_sphere_rcs(ka, a=1.0, maxh=0.3, order=2, intorder=12):
    r"""VECTOR-MAXWELL EFIE solve: the monostatic (backscatter) radar cross-section of a PERFECTLY
    CONDUCTING (PEC) sphere via the ngsolve.bem ``MaxwellSingleLayerPotentialOperator`` (the electric-field
    integral equation), verified against the analytic EM Mie RCS (:func:`mie_pec_sphere_rcs`).  This is the
    ACTUAL full-wave vector solve the scalar/acoustic Mie helpers above only gated analytically -- the open
    counterpart of a commercial high-frequency integral (MoM) RCS solver.

    Formulation.  The PEC condition gamma_t E = 0 gives the EFIE for the surface current J in the
    div-conforming (RWG) surface space ``HDivSurface``:  V[J] = -i gamma_t E_inc, where V is the NGBEM
    Maxwell single-layer operator (which carries the kappa scaling internally -- the leading ``-i`` and the
    far-field ``kappa^2`` below were identified from a 1/ka^2 ratio scan, then verified ka-independent).  The
    incident wave is x-polarised, +z propagating, ``E_inc = xhat exp(i kappa z)``.  Solve V J = rhs (GMRES,
    V complex symmetric, L2-mass preconditioner).  The backscatter far-field amplitude (xhat = -zhat) is the
    transverse part of ``A = INT_Gamma J(y) exp(i kappa z) dGamma``, and (Z0 = 1)

        sigma_b = (kappa^2 / 4 pi) (|A_x|^2 + |A_y|^2),   sigma_b/(pi a^2)  vs  mie_pec_sphere_rcs.

    Verified to ~1e-4 across ka = 0.3 .. 3 (Rayleigh -> resonance -> optical onset).  Returns
    ``{ndof, ka, sigma_back_over_pa2, mie, rel_err}``."""
    from ngsolve.krylovspace import GMResSolver
    kappa = ka / a
    mesh = ng.Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), a)).GenerateMesh(maxh=maxh)).Curve(4)
    fes = ng.HDivSurface(mesh, order=order, complex=True)
    u, v = fes.TnT()
    V = bem.MaxwellSingleLayerPotentialOperator(fes, kappa, intorder=intorder)
    Einc = ng.CF((ng.exp(1j * kappa * ng.z), 0, 0))                       # x-pol, +z plane wave
    rhs = LinearForm((-1j) * (Einc * v.Trace()) * ds(bonus_intorder=intorder - 4)).Assemble()
    mass = BilinearForm(u.Trace() * v.Trace() * ds(bonus_intorder=intorder - 4)).Assemble().mat
    pre = mass.Inverse(freedofs=fes.FreeDofs())
    J = GridFunction(fes)
    inv = GMResSolver(V.mat, pre, printrates=False, maxiter=4000, tol=1e-9)
    J.vec.data = inv * rhs.vec
    ph = ng.exp(1j * kappa * ng.z)                                        # exp(-i k xhat.y), xhat = -zhat
    Ax = complex(Integrate(J[0] * ph * ds(bonus_intorder=intorder - 4), mesh))
    Ay = complex(Integrate(J[1] * ph * ds(bonus_intorder=intorder - 4), mesh))
    sigma_b = (kappa ** 2 / (4.0 * np.pi)) * (abs(Ax) ** 2 + abs(Ay) ** 2)
    bem_val = float(sigma_b / (np.pi * a * a))
    mie = mie_pec_sphere_rcs(ka)["sigma_back_over_pa2"]
    return {"ndof": fes.ndof, "ka": float(ka), "sigma_back_over_pa2": bem_val,
            "mie": float(mie), "rel_err": abs(bem_val - mie) / mie}


# ---------------------------------------------------------------------------
# GIBC -- Generalized Impedance Boundary Condition (curvature-corrected SIBC)
# ---------------------------------------------------------------------------
#
# The Leontovich SIBC assumes a FLAT conductor surface.  For curved surfaces the
# Dirichlet-to-Neumann map of the conductor interior has curvature corrections
# (Mitzner 1968; Senior & Volakis 1995):
#
#   order=0  Leontovich:   Z_s = sqrt(jω μ_c / σ_c)
#   order=1  Mitzner:      Z_s^(1) = Z_s (1 + H / k_c)
#   order=2  Senior:       Z_s^(2) = Z_s (1 + H/k_c + (H² − K/2)/k_c²)
#
# where H = mean curvature = (1/R₁ + 1/R₂)/2,  K = Gaussian curvature = 1/(R₁R₂),
# k_c = sqrt(jω μ_c σ_c) (complex conductor wave number, = (1+j)/δ for good conductor),
# and H, K are obtained from specialcf.Weingarten(3) which requires a curved mesh
# (geom_order ≥ 2).
#
# For a sphere of radius R: H = 1/R, K = 1/R², so
#   Z_s^(1) = Z_s (1 + 1/(R k_c)),  Z_s^(2) = Z_s (1 + 1/(R k_c) + 1/(2 R² k_c²))
#
# These are CoefficientFunctions that vary over the surface -- the same Z_s CF
# is used directly in any FEM Robin BC: a += (1/Z_s_cf) * inner(n×H_t, v) * ds.
#
# DtN Schur interpretation: SIBC = leading-order (0th) approximation to the exact
# conductor DtN (Schur complement of the conductor subdomain); GIBC = higher-order
# asymptotic approximation, converging to the exact DtN as |k_c R| → ∞.


def gibc_impedance_cf(omega, mu_c, sigma_c, order=1):
    """GIBC surface impedance CoefficientFunction, valid on curved boundary meshes.

    Parameters
    ----------
    omega : float   angular frequency [rad/s]
    mu_c  : float   conductor permeability [H/m]  (= mu_r * mu_0)
    sigma_c : float conductor conductivity [S/m]
    order : int     0 = Leontovich SIBC (scalar), 1 = Mitzner curvature correction,
                    2 = Senior-Volakis second-order correction

    Returns a boundary CoefficientFunction for Z_s [Ω] (complex).  For order > 0
    the mesh must be curved (geom_order ≥ 2) so Weingarten is non-zero.

    DtN Schur hierarchy:
      order=0 (SIBC)  : local Leontovich approximation to conductor DtN
      order=1 (GIBC1) : + curvature correction H/k_c
      order=2 (GIBC2) : + (H² − K/2)/k_c² (Gaussian curvature)"""
    Z_s0 = cmath.sqrt(1j * omega * mu_c / sigma_c)       # Leontovich (complex scalar)
    if order == 0:
        return ng.CoefficientFunction(Z_s0)
    k_c = cmath.sqrt(1j * omega * mu_c * sigma_c)        # = j omega mu_c / Z_s0
    W = specialcf.Weingarten(3)
    two_H = Trace(W)                                      # = 1/R_1 + 1/R_2
    H = two_H * 0.5                                       # mean curvature
    # order=1: Z_s0 * (1 + H / k_c)
    cf = ng.CoefficientFunction(Z_s0) * (1.0 + H * (1.0 / k_c))
    if order == 1:
        return cf
    # order=2: add Gaussian-curvature correction K = ((2H)^2 - Trace(W^2)) / 2
    trW2 = Trace(W * W)
    K = (two_H * two_H - trW2) * 0.5
    cf2 = cf + ng.CoefficientFunction(Z_s0) * (H * H - K * 0.5) * (1.0 / k_c ** 2)
    return cf2


def sphere_gibc_benchmark(R=1.0, omega=1e6, mu_r=1.0, sigma=5.8e7, maxh=0.3, order=3):
    """GIBC benchmark on a sphere of radius R: compare SIBC (order 0), GIBC1 (order 1),
    and GIBC2 (order 2) impedances against the exact surface impedance from the
    spherical harmonic DtN.

    For a solid conducting sphere, the exact surface impedance for the dipole mode
    (n=1 spherical harmonic) is derived from the interior Mie series:

        Z_s_exact = j ω μ_c / k_c  ×  [j₀(k_c R) / j₁(k_c R)]  ... (TE/TM, l=1)

    but for the purpose of the GIBC *asymptotic* validation, we compare the relative
    correction term H/k_c against its mean value on the sphere (should be 1/R/k_c).

    Returns dict with:
      R, omega, k_c, skin_depth, kg = |k_c|*R (normalised radius),
      Z_s_sibc  (scalar, Leontovich),
      Z_s_gibc1 (spatial average over sphere, order 1),
      Z_s_gibc2 (spatial average over sphere, order 2),
      gibc1_correction = mean|Z_s_gibc1/Z_s_sibc - 1| (should ≈ 1/(|k_c|R)),
      gibc2_correction = mean|Z_s_gibc2/Z_s_sibc - 1| (should ≈ 1/(|k_c|R) + 1/(2|k_c|²R²))
    """
    mu_c = mu_r * 4e-7 * math.pi
    mesh = ng.Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), R)).GenerateMesh(maxh=maxh)).Curve(3)
    Z_s0_scalar = cmath.sqrt(1j * omega * mu_c / sigma)
    k_c = cmath.sqrt(1j * omega * mu_c * sigma)
    delta = math.sqrt(2.0 / (omega * mu_c * sigma))      # skin depth
    area = float(Integrate(ng.CoefficientFunction(1.0), mesh, ng.BND).real)
    # order=1 GIBC: Z_s0 * (1 + H/k_c), H=1/R on sphere
    cf1 = gibc_impedance_cf(omega, mu_c, sigma, order=1)
    avg1 = complex(Integrate(cf1, mesh, ng.BND)) / area
    # order=2 GIBC: additionally (H^2 - K/2)/k_c^2 = 1/(2R^2 k_c^2) on sphere
    cf2 = gibc_impedance_cf(omega, mu_c, sigma, order=2)
    avg2 = complex(Integrate(cf2, mesh, ng.BND)) / area
    corr1 = abs(avg1 / Z_s0_scalar - 1.0)
    corr2 = abs(avg2 / Z_s0_scalar - 1.0)
    return {"R": R, "omega": omega, "k_c": k_c, "skin_depth": delta,
            "kg": abs(k_c) * R,
            "Z_s_sibc": Z_s0_scalar, "Z_s_gibc1": avg1, "Z_s_gibc2": avg2,
            "gibc1_correction": corr1, "gibc2_correction": corr2,
            "expected_gibc1_correction": 1.0 / (abs(k_c) * R),
            "expected_gibc2_correction": 1.0 / (abs(k_c) * R) + 1.0 / (2.0 * abs(k_c) ** 2 * R ** 2)}


# ---------------------------------------------------------------------------
# BEM exterior DtN (Laplace) -- Schur complement of the exterior subdomain
# ---------------------------------------------------------------------------
#
# The EXTERIOR Laplace DtN Λ_ext maps the Dirichlet data u|_Γ to its outward
# Neumann data ∂u/∂n|_Γ (n outward from conductor = outward from interior).
#
# BIE (first-kind, exterior, direct representation) -- NGSolve sign convention:
#
#   NGSolve's DoubleLayerPotentialOperator uses n_y = outward from the enclosed
#   volume (standard Sauter-Schwab convention).  The Calderon identity gives
#
#     Interior DtN:  V σ_int = (+½ M + K) u   →  σ = V⁻¹ (+½ M + K) u,  eigenvalue +n/R
#     Exterior DtN:  V σ_ext = (-½ M + K) u   →  σ = V⁻¹ (-½ M + K) u,  eigenvalue -(n+1)/R
#
# where σ = ∂u/∂n with n = outward from the enclosed domain (= +r̂ for sphere).
# The sign difference is the exterior-side jump of the double-layer (−½) vs the
# interior-side jump (+½); K uses the OUTWARD normal n_y throughout.
#
# For the spherical harmonic Y_n^m: Λ_ext Y_n = −(n+1)/R Y_n, i.e.
# the exterior DtN eigenvalue is −(n+1)/R.  For the dipole (n=1): Λ f = −2/R f.
#
# DtN Schur connection:
#   - FEM interior: K_FEM u_I + K_IB u_B = f  (interior + boundary DOFs)
#   - Schur of FEM: S_FEM = K_BB − K_BI K_II⁻¹ K_IB  (FEM DtN on boundary)
#   - BEM exterior: Λ_ext = V⁻¹ (−½M + K)              (BEM DtN on boundary)
#   - Coupled system: (S_FEM − Λ_ext) u_B = rhs_condensed
#     (the weak boundary term is −∮(∂u/∂n)v = −⟨Λ_ext u,v⟩; Λ_ext is
#      negative-definite, so −Λ_ext ADDS a positive-definite term.  In the FEM
#      assembly the weak DtN is −P^T M_bnd Λ_ext P -- the SurfaceL2 mass M_bnd
#      turns the DtN DOF-map into the weak form; see laplace_fem_bem_schur.)
#   This gives an EXACT open boundary without Kelvin transformation or mesh truncation.
#
# The function below validates the BEM DtN on a sphere against the known eigenvalue,
# and assembles the dense DtN matrix for use in FEM-BEM coupling.


def _dense_mat(op, ndof):
    """Densify a linear operator op: R^ndof → R^ndof by applying it to each unit vector."""
    x = op.CreateColVector()
    y = op.CreateColVector()
    M = np.zeros((ndof, ndof), dtype=complex)
    for j in range(ndof):
        x[:] = 0.0; x[j] = 1.0
        y.data = op * x
        M[:, j] = y.FV().NumPy()[:]
    return M


def laplace_exterior_dtn(fes_bnd, intorder=12):
    """Assemble the exterior Laplace DtN matrix Λ = V⁻¹ (-½ M + K) for the closed
    surface mesh encoded by ``fes_bnd`` (a SurfaceL2 space).

    The matrix maps the DOF vector of the boundary potential u|_Γ to the DOF vector
    of the outward normal derivative σ = ∂u/∂n|_Γ (outward from the enclosed volume).

    The exterior first-kind BIE is V σ = (-½ M + K) u (exterior-side jump −½ of the
    double-layer, same K as the interior BIE V σ = (+½ M + K) u but with negative mass).
    This follows the Sauter-Schwab Calderon identity with n_y = outward from the enclosed
    domain (the NGSolve DoubleLayerPotentialOperator convention).

    Returns the dense DtN matrix as a complex numpy array (ndof × ndof).
    For a sphere of radius R the DtN eigenvalue for dipole mode Y₁^m is −2/R
    (verified by :func:`sphere_exterior_dtn_eigenvalue`).

    Sign in FEM-BEM Schur coupling: the weak boundary term is −∮(∂u/∂n)v dS =
    −⟨Λ_ext u, v⟩, so the contribution to the FEM stiffness is −P^T Λ_ext P
    (Λ_ext is negative-definite, so this ADDS a positive-definite term; see
    :func:`fem_bem_coupling.laplace_fem_bem_schur`).

    ``fes_bnd`` may be a SurfaceL2 restricted to one boundary of a VOLUME mesh
    (``definedon=mesh.Boundaries("outer")``).  Such a space allocates DOFs on the
    other boundaries too and marks them non-free, which would make the BEM single
    layer V singular (zero rows); the assembly below therefore solves on the FREE
    sub-block and scatters zeros back.  For a standalone closed surface every DOF
    is free and this is a no-op."""
    u, v = fes_bnd.TnT()
    V = bem.SingleLayerPotentialOperator(fes_bnd, intorder=intorder)
    K = bem.DoubleLayerPotentialOperator(fes_bnd, fes_bnd, intorder=intorder)
    M = BilinearForm(u.Trace() * v.Trace() * ds(bonus_intorder=intorder - 4)).Assemble().mat
    ndof = fes_bnd.ndof
    V_d = _dense_mat(V.mat, ndof)
    K_d = _dense_mat(K.mat, ndof)
    M_d = _dense_mat(M, ndof)
    rhs_mat = -0.5 * M_d + K_d              # Λ_ext = V⁻¹ (-½ M + K)  [exterior BIE]

    free = np.array([fes_bnd.FreeDofs()[i] for i in range(ndof)], dtype=bool)
    if free.all():
        return np.linalg.solve(V_d, rhs_mat)
    idx = np.where(free)[0]
    dtn_free = np.linalg.solve(V_d[np.ix_(idx, idx)], rhs_mat[np.ix_(idx, idx)])
    dtn = np.zeros((ndof, ndof), dtype=dtn_free.dtype)
    dtn[np.ix_(idx, idx)] = dtn_free
    return dtn


def sphere_exterior_dtn_eigenvalue(R=1.0, maxh=0.4, order=1, intorder=10):
    """Validate the exterior Laplace DtN on a sphere of radius R against the exact
    dipole eigenvalue −(n+1)/R = −2/R (spherical harmonic mode n=1).

    Prescribes f = z/R = cos θ (the l=1 dipole) as a boundary potential, applies the
    assembled DtN matrix, and checks that σ_BEM ≈ (−2/R) × f pointwise.  The test
    uses a low-order space (order=1) so the DtN densification is fast enough to run
    as a gate check.

    Returns ``{ndof, R, dtn_eigenvalue_expected, dtn_eigenvalue_bem, rel_err}``."""
    mesh = ng.Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), R)).GenerateMesh(maxh=maxh)).Curve(3)
    fes = ng.SurfaceL2(mesh, order=order, dual_mapping=True)
    u, v = fes.TnT()
    dtn = laplace_exterior_dtn(fes, intorder=intorder)
    # L2 project f = z/R onto fes via mass-matrix solve (mass = <u,v> on the surface mesh)
    mass_bf = BilinearForm(u.Trace() * v.Trace() * ds(bonus_intorder=intorder - 4)).Assemble()
    mass = mass_bf.mat
    rhs_f = LinearForm((ng.z / R) * v.Trace() * ds(bonus_intorder=intorder - 4)).Assemble()
    gf_f = GridFunction(fes)
    gf_f.vec.data = mass.Inverse(freedofs=fes.FreeDofs()) * rhs_f.vec
    f_vec = gf_f.vec.FV().NumPy()[:]
    sigma_vec = dtn @ f_vec                              # σ = Λ f
    # Exact: σ = −2/R × f for the n=1 harmonic
    expected_vec = (-2.0 / R) * f_vec
    # Mass-weighted norms and Rayleigh-quotient DtN eigenvalue
    M_d = _dense_mat(mass, fes.ndof)
    err_vec = sigma_vec - expected_vec
    err_norm = float(np.sqrt(abs(np.real(err_vec.conj() @ M_d @ err_vec))))
    f_norm   = float(np.sqrt(abs(np.real(expected_vec.conj() @ M_d @ expected_vec))))
    lam_bem = float(np.real(sigma_vec.conj() @ M_d @ f_vec) /
                    np.real(f_vec.conj() @ M_d @ f_vec))
    return {"ndof": fes.ndof, "R": R,
            "dtn_eigenvalue_expected": -2.0 / R,
            "dtn_eigenvalue_bem": lam_bem,
            "rel_err": abs(lam_bem - (-2.0 / R)) / (2.0 / R),
            "l2_rel_err": err_norm / f_norm if f_norm > 0 else float("nan")}


# ---------------------------------------------------------------------------
# DtN spectrum vs mesh size -- the spectral explanation of coarse-mesh accuracy
# ---------------------------------------------------------------------------
#
# Every open-boundary closure (Kelvin transform, BEM coupling, PML, infinite
# elements, asymptotic Robin BC) is, on the truncation surface Γ, an
# approximation of ONE continuous operator: the exterior Laplace DtN Λ_ext.
# On a sphere of radius R it is diagonalised by the spherical harmonics with the
# MESH-INDEPENDENT eigenvalue ladder λ_n = −(n+1)/R (degree n, multiplicity
# 2n+1).  The discrete DtN matrix Λ_h reproduces the LOW-degree eigenvalues
# accurately and almost independently of h even on coarse meshes (the
# eigenfunctions are smooth low harmonics); the per-mode error grows with n.
#
# Because a compact source's exterior field is dominated by the low multipoles
# (the n-th decays like r^{−(n+1)}), the coarse mesh already resolves everything
# that matters -- which is exactly the coarse-mesh accuracy that Kameari showed
# empirically for the Kelvin transformation, here read off the matrix spectrum.
#
# Knowledge: dtn_coarse_mesh.  Companion: laplace_exterior_dtn (the operator).


def exterior_dtn_spectrum(R=1.0, maxh=0.6, order=1, intorder=10, nmax=4):
    """Eigenvalue spectrum of the dense exterior Laplace DtN matrix Λ_h on a
    sphere of radius R, matched degree-by-degree to the analytic ladder
    λ_n = −(n+1)/R (spherical-harmonic degree n, multiplicity 2n+1).

    The eigenvalues of Λ_ext are a property of the CONTINUOUS operator (no mesh);
    this measures how well the matrix at mesh size ``maxh`` reproduces them.
    The leading negative eigenvalues (sorted toward zero first) are bucketed by
    the analytic multiplicity 2n+1 and each bucket mean compared to −(n+1)/R.

    Parameters
    ----------
    R        : float  sphere radius.
    maxh     : float  surface mesh size.
    order    : int    SurfaceL2 order (1 keeps the densification tractable).
    intorder : int    BEM quadrature order.
    nmax     : int    highest spherical-harmonic degree to match.

    Returns dict ``{ndof, R, order, maxh, eigenvalues, modes}`` where ``modes``
    is a list of ``{n, multiplicity, lambda_exact, lambda_mean, rel_err,
    spread}`` (spread = max−min within the degenerate bucket)."""
    mesh = ng.Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), R)).GenerateMesh(maxh=maxh)).Curve(3)
    fes = ng.SurfaceL2(mesh, order=order, dual_mapping=True)
    dtn = laplace_exterior_dtn(fes, intorder=intorder)

    w = np.linalg.eigvals(dtn)
    # Physical DtN eigenvalues are real and negative; drop numerical noise and
    # any spurious non-negative eigenvalues, then sort toward zero (n=0 first).
    real_mask = np.abs(np.imag(w)) <= 1e-6 * (1.0 + np.abs(w))
    w = np.real(w[real_mask])
    w = np.sort(w[w < 0.0])[::-1]                 # descending: −1/R (n=0) leads

    modes = []
    idx = 0
    for n in range(nmax + 1):
        mult = 2 * n + 1
        if idx + mult > len(w):
            break
        grp = w[idx:idx + mult]
        lam_exact = -(n + 1) / R
        lam_mean = float(np.mean(grp))
        modes.append({
            "n": n,
            "multiplicity": mult,
            "lambda_exact": lam_exact,
            "lambda_mean": lam_mean,
            "rel_err": abs(lam_mean - lam_exact) / abs(lam_exact),
            "spread": float(np.max(grp) - np.min(grp)),
        })
        idx += mult
    return {"ndof": fes.ndof, "R": R, "order": order, "maxh": maxh,
            "eigenvalues": [float(x) for x in w[:(nmax + 1) ** 2]],
            "modes": modes}


def dtn_spectrum_vs_mesh(R=1.0, maxh_list=(0.5, 0.4, 0.3),
                         order=1, intorder=10, nmax=4, band_tol=0.005,
                         low_nmax=2):
    """Sweep :func:`exterior_dtn_spectrum` across mesh sizes and summarise the
    spectral facts that explain the coarse-mesh accuracy of open-boundary
    methods (Kameari's Kelvin-transform observation, read off the matrix).

    NB on a sphere the surface mesh FLOORS at a minimum triangulation: for R=1,
    every ``maxh >= 0.5`` yields the same coarsest mesh.  Pass distinct meshes
    by using ``maxh <= 0.5`` (e.g. 0.5, 0.4, 0.3 give ndof 336, 564, 924 at
    order=1).  The DENSIFICATION is O(ndof^2); maxh=0.3 already takes minutes.

    Returns dict:
      per_mesh                : list of exterior_dtn_spectrum() results
      maxh_list / ndof_list   : the meshes used and their boundary DOF counts
      error_table             : error_table[n] = [rel_err of degree n at each
                                mesh] (None where the mesh cannot resolve n)
      coarse_low_mode_max_err : max rel_err over degrees n<=low_nmax on the
                                COARSEST (first) mesh -- the headline "already
                                accurate on a coarse mesh" number
      low_mode_max_err        : same, but max over ALL meshes (the low-mode floor)
      degree_monotonic        : per mesh, True if rel_err strictly increases with
                                degree n -- the spectral signature (accurate at
                                the bottom, degraded toward the resolution limit)
      degree_growth           : per mesh, rel_err(highest matched n)/rel_err(n=0)
                                -- how strongly error is ordered by degree
      accurate_band           : per mesh, highest degree n with rel_err < band_tol
                                -- grows as maxh shrinks (refinement WIDENS the
                                band of well-resolved modes)

    The spectral story (what the numbers show, vs a naive "it's mesh-independent"):
      * the LOW modes are already accurate to a small floor on the COARSEST mesh
        (coarse_low_mode_max_err is small);
      * at EVERY mesh the error is ordered by degree (degree_monotonic / growth);
      * refinement lowers every mode steeply (~O(h^4) order-1 Galerkin eigenvalue
        superconvergence, measured rate p~3.9) and widens accurate_band, but the
        low modes -- which dominate a compact source's field -- were already
        below engineering tolerance before refining.
    That is Kameari's coarse-mesh accuracy expressed as a property of Lambda_h.
    """
    per_mesh = [exterior_dtn_spectrum(R=R, maxh=h, order=order,
                                      intorder=intorder, nmax=nmax)
                for h in maxh_list]

    error_table = {}
    for n in range(nmax + 1):
        error_table[n] = [next((mm["rel_err"] for mm in res["modes"]
                                if mm["n"] == n), None)
                          for res in per_mesh]

    def _low_errs(results):
        return [mm["rel_err"] for res in results for mm in res["modes"]
                if mm["n"] <= low_nmax]

    coarse_low = _low_errs([per_mesh[0]]) if per_mesh else []
    all_low = _low_errs(per_mesh)
    coarse_low_mode_max_err = max(coarse_low) if coarse_low else float("nan")
    low_mode_max_err = max(all_low) if all_low else float("nan")

    degree_monotonic, degree_growth, accurate_band = [], [], []
    for res in per_mesh:
        errs = [m["rel_err"] for m in res["modes"]]
        degree_monotonic.append(all(errs[i] < errs[i + 1]
                                    for i in range(len(errs) - 1)))
        degree_growth.append(errs[-1] / errs[0] if errs and errs[0] > 0 else float("nan"))
        band = -1
        for m in res["modes"]:
            if m["rel_err"] < band_tol:
                band = m["n"]
            else:
                break
        accurate_band.append(band)

    return {"per_mesh": per_mesh,
            "maxh_list": list(maxh_list),
            "ndof_list": [res["ndof"] for res in per_mesh],
            "error_table": error_table,
            "coarse_low_mode_max_err": coarse_low_mode_max_err,
            "low_mode_max_err": low_mode_max_err,
            "degree_monotonic": degree_monotonic,
            "degree_growth": degree_growth,
            "accurate_band": accurate_band,
            "band_tol": band_tol,
            "low_nmax": low_nmax}
