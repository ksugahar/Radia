"""Time-harmonic A-V eddy currents in a terminal-driven 3-D conductor.

This is the independent route the fin PEEC gate 3 asks for.  Both existing
routes -- the BEM-A reference and the surface PEEC -- impose a Leontovich
surface impedance, so neither can check the other's central assumption.  Here
the conductor's interior is solved, the skin profile comes out of the solution
rather than going into it, and nothing is assumed about the ratio of skin
depth to feature size.

Formulation (A-V, "modified" gauge).  With ``exp(j omega t)`` and
``E = -j omega A - grad V``:

    find A in HCurl(Omega), V in H1(Omega_c) such that
        int_Omega  nu (curl A).(curl A')
      + int_Omega_c  j omega sigma (A + grad V).A'                = 0
        int_Omega_c  j omega sigma (A + grad V).(grad V')         = 0

``HCurl`` is built with ``nograds=True``: the curl-free part of the electric
field lives in ``grad V`` inside the conductor, and outside it there is no
current to represent, so removing the gradient shapes both fixes the gauge and
makes the air block non-singular.  The two equations are the same bilinear
form tested against the two fields, so the assembled system is symmetric.

Driving.  ``V`` is Dirichlet: zero on the sink cap, one on the source cap.
That is a unit potential drive, not a unit current, so every returned quantity
is rescaled to the terminal current the solution produces.  Driving by
potential keeps the system linear and square; driving by current would need a
constraint row and gives the same answer after the same rescale.

Scope.  Linear isotropic media, one conductor, one terminal pair, no imposed
external source field.  The outer boundary is a flux-parallel truncation, so
its adequacy is a mesh question the caller has to answer by enlarging it.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

MU0 = 4e-7 * np.pi


@dataclass(frozen=True)
class EddyAPhiResult:
    """One solved operating point, rescaled to unit terminal current."""

    frequency_hz: float
    sigma_S_per_m: float
    terminal_current_A: complex
    terminal_impedance_ohm: complex
    resistance_ohm: float
    inductance_H: float
    ohmic_loss_W: float
    ndof: int
    ndof_A: int
    ndof_V: int
    skin_depth_m: float
    gf_A: object
    gf_V: object
    mesh: object
    conductor: str
    scale: complex

    def current_density(self):
        """``J = -j omega sigma (A + grad Vt)`` at unit terminal current."""
        import ngsolve as ng

        omega = 2.0 * np.pi * self.frequency_hz
        return (-1j * omega * self.sigma_S_per_m
                * (self.gf_A + ng.grad(self.gf_V)) * self.scale)

    def flux_density(self):
        """``B = curl A`` as a CoefficientFunction."""
        import ngsolve as ng

        return ng.curl(self.gf_A) * self.scale


def solve_eddy_aphi(mesh, *, conductor, source, sink, frequency_hz, sigma,
                    order=2, mu_r=1.0, outer="outer", condense=True,
                    solver="umfpack", air_sigma_ratio=1e-8,
                    print_dofs=False):
    """Solve one terminal-driven eddy-current operating point.

    Args:
        mesh: NGSolve mesh whose volume regions include ``conductor``.
        conductor: material name of the conducting region.
        source, sink: boundary names of the two terminal caps.
        frequency_hz: drive frequency; must be positive.
        sigma: conductivity of ``conductor`` in S/m.
        order: HCurl / H1 polynomial order.  Two is the useful minimum; the
            skin profile is exponential and low order needs many more
            elements to reach the same error.
        mu_r: relative permeability of every region.
        outer: boundary name carrying the flux-parallel truncation.
        condense: use static condensation.
        solver: direct sparse solver passed to ``.Inverse``.

    Returns:
        EddyAPhiResult rescaled so the terminal current is exactly 1 A.
    """
    import ngsolve as ng

    if frequency_hz <= 0:
        raise ValueError("frequency_hz must be positive; "
                         "this is a time-harmonic formulation")
    if sigma <= 0:
        raise ValueError("sigma must be positive")
    omega = 2.0 * np.pi * float(frequency_hz)
    nu = 1.0 / (MU0 * float(mu_r))
    delta = np.sqrt(2.0 / (omega * MU0 * float(mu_r) * float(sigma)))

    region = mesh.Materials(conductor)
    if region.Mask().NumSet() == 0:
        raise ValueError(f"no elements carry the material name {conductor!r}")
    for name in (source, sink, outer):
        if mesh.Boundaries(name).Mask().NumSet() == 0:
            raise ValueError(f"no boundary carries the name {name!r}")

    # Unknowns are A and the SCALED potential Vt = V / (j omega).  The scaling
    # is what makes the two coupling blocks equal, so the assembled system is
    # symmetric; every result therefore carries a j omega that has to be put
    # back -- the applied terminal voltage is j omega * (Vt_source - Vt_sink),
    # and J = -j omega sigma (A + grad Vt).
    #
    # A must be gauged.  With full HCurl the pair (A, Vt) is redundant --
    # A -> A + grad phi, Vt -> Vt - phi leaves E untouched and costs nothing in
    # the curl-curl term -- and the solve runs away along exactly that mode.
    # `nograds=True` removes the gradient shapes, which removes the redundancy
    # itself rather than penalising it.
    sigma_air = float(sigma) * float(air_sigma_ratio)
    sigma_cf = mesh.MaterialCF({conductor: float(sigma)}, default=sigma_air)

    # A is constrained on the WHOLE outer surface, terminal caps included.
    # The caps carry a name of their own because V is prescribed there, but
    # they are still part of the truncation, and leaving them out gives A the
    # natural curl-curl condition -- tangential H = 0.  On the end plane of a
    # conductor carrying current along its axis H is azimuthal and tangential
    # there, so that condition is false, and it cost 19% in the DC limit where
    # the answer is simply L / (sigma A).
    a_dirichlet = "|".join((outer, source, sink))
    fes_a = ng.HCurl(mesh, order=order, complex=True, nograds=True,
                     dirichlet=a_dirichlet)
    fes_v = ng.H1(mesh, order=order, complex=True, definedon=region,
                  dirichlet=f"{source}|{sink}")
    fes = ng.FESpace([fes_a, fes_v])
    (a_trial, v_trial), (a_test, v_test) = fes.TnT()

    form = ng.BilinearForm(fes, symmetric=True, condense=condense)
    form += nu * ng.curl(a_trial) * ng.curl(a_test) * ng.dx
    form += 1j * omega * sigma_cf * a_trial * a_test * ng.dx
    drive = 1j * omega * float(sigma)
    form += drive * ng.grad(v_trial) * a_test * ng.dx(definedon=region)
    form += drive * (a_trial + ng.grad(v_trial)) * ng.grad(v_test) * ng.dx(
        definedon=region)

    gfu = ng.GridFunction(fes)
    gf_a, gf_v = gfu.components
    # Unit potential drive: the caps are the only place V is prescribed.
    gf_v.Set(ng.CF(1.0 + 0j), definedon=mesh.Boundaries(source))
    # The same shape serves as the test function for the terminal current
    # below -- one on the source cap, zero on the sink cap.
    indicator = ng.GridFunction(fes_v)
    indicator.Set(ng.CF(1.0 + 0j), definedon=mesh.Boundaries(source))

    rhs = ng.LinearForm(fes)
    form.Assemble()
    rhs.Assemble()
    residual = rhs.vec.CreateVector()
    residual.data = rhs.vec - form.mat * gfu.vec
    if condense:
        residual.data += form.harmonic_extension_trans * residual
    inverse = form.mat.Inverse(freedofs=fes.FreeDofs(condense),
                               inverse=solver)
    gfu.vec.data += inverse * residual
    if condense:
        gfu.vec.data += form.harmonic_extension * gfu.vec
        gfu.vec.data += form.inner_solve * residual

    current = terminal_current(mesh, gf_a, gf_v, indicator, sigma=sigma,
                               omega=omega, region=region)
    if current == 0:
        raise RuntimeError("the solution carries no terminal current; check "
                           "that the caps bound the conductor")
    # Vt = 1 was prescribed, so the applied voltage is j omega volts.
    applied = 1j * omega
    impedance = applied / current
    scale = 1.0 / current
    loss = 0.5 * impedance.real  # P = |I|^2 R / 2 at I = 1 A peak

    if print_dofs:
        print(f"eddy A-V: ndof={fes.ndof} (A {fes_a.ndof}, V {fes_v.ndof}) "
              f"delta={delta * 1e3:.4f} mm", flush=True)

    return EddyAPhiResult(
        frequency_hz=float(frequency_hz), sigma_S_per_m=float(sigma),
        terminal_current_A=complex(current),
        terminal_impedance_ohm=complex(impedance),
        resistance_ohm=float(impedance.real),
        inductance_H=float(impedance.imag / omega),
        ohmic_loss_W=float(loss), ndof=int(fes.ndof), ndof_A=int(fes_a.ndof),
        ndof_V=int(fes_v.ndof), skin_depth_m=float(delta), gf_A=gf_a,
        gf_V=gf_v, mesh=mesh, conductor=conductor, scale=complex(scale))


def terminal_current(mesh, gf_a, gf_v, indicator, *, sigma, omega, region):
    """Terminal current by flux recovery, as a volume integral.

    Integrating ``J . n`` over the cap looks more direct and does not work:
    ``grad`` of an H1 space restricted to the conductor evaluates to zero on a
    boundary facet, because the facet does not resolve its neighbouring volume
    element, and the cap integral comes back exactly zero.

    Since ``div J = 0`` in the conductor, testing against any ``w`` that is one
    on the source cap and zero on the sink cap recovers the same number from
    the interior:

        int_Omega_c J . grad w  =  oint w J.n  =  -I

    and it is exact for the discrete solution rather than a facet
    interpolation of it.
    """
    import ngsolve as ng

    j_field = -1j * omega * sigma * (gf_a + ng.grad(gf_v))
    return -complex(ng.Integrate(j_field * ng.grad(indicator),
                                 mesh, definedon=region))


def ohmic_loss(result):
    """Volume loss ``int |J|^2 / (2 sigma)`` at the rescaled unit current."""
    import ngsolve as ng

    j_field = result.current_density()
    return float(ng.Integrate(
        (j_field * ng.Conj(j_field)).real / (2.0 * result.sigma_S_per_m),
        result.mesh, definedon=result.mesh.Materials(result.conductor)).real)
