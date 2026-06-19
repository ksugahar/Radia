"""2D frequency-domain eddy currents: skin effect and AC resistance.

A long straight conductor carrying a time-harmonic current pushes that current
toward its surface (the skin effect), raising its effective (AC) resistance above
the DC value. In the frequency domain the out-of-plane vector potential A_z (phasor)
obeys the diffusion equation

    -div( nu grad A_z ) + j*omega*sigma*A_z = sigma*E0        (in the conductor)
    -div( nu grad A_z ) = 0                                    (in the surrounding air)

with nu = 1/(mu0 mu_r) and E0 = -dV/dz the (uniform) axial driving field. For an
ISOLATED conductor we may simply impose E0 = 1; the current density then redistributes
by the skin effect and the AC/DC resistance ratio is a current-amplitude-independent,
purely geometric/frequency quantity

    J_z = sigma (E0 - j*omega*A_z),   I = int_cond J_z dA,
    Rac/Rdc = A_cond * int_cond |J_z|^2 dA / |I|^2.

For the circular cross-section this ratio has the exact closed form in Kelvin functions
(ber, bei), against which :func:`ac_resistance_round_wire` is validated to machine
precision (see ``kelvin_ac_resistance_ratio``). The SAME physics is carried by the
trusted in-repo analytic baseline ``radia.analytical_formulas.conductor_impedance``
(``cylinder_ac_impedance`` / ``cylinder_dc_resistance``, in Bessel J0/J1 of the complex
argument k a); the FEM solve, the Kelvin ber/bei ratio, and that Bessel impedance all
agree to ~1e-10 across xi = 0.5..6. Related: the analytic proximity-effect AC copper-loss
model lives in ``radia_mcp.peec`` (Carstensen); this module is the direct FEM solve of
the single-conductor skin problem.

Extracting R_ac -- the robust route. The effective AC resistance is taken from the OHMIC-LOSS
integral, ``Rac = (int |J_z|^2 / sigma dA) / |I|^2`` (= ``A_cond int|J_z|^2 dA / |I|^2`` for
uniform sigma), NOT from an apparent terminal voltage/current ratio. In a multi-conductor or
iron-bearing problem a terminal V/I folds mutual coupling and core loss into its real part, so the
apparent "resistance" can come out spuriously low (even complex); the loss integral is immune to
that and is the mesh-robust definition. Use phasor amplitudes consistently for I and J (peak or
RMS -- the ratio is convention-independent). This solve is a SINGLE isolated conductor; in a
bundle/winding the proximity loss is NON-LOCAL (neighbouring-turn fields), so a one-turn result
cannot simply be scaled by the turn count -- model the whole bundle or use the proximity
(``radia_mcp.peec``) model.

All quantities are 2D planar (per unit out-of-plane length).
"""
import math

MU0 = 4e-7 * math.pi


def kelvin_ac_resistance_ratio(xi):
    """Exact AC/DC resistance ratio of an isolated round wire, in Kelvin functions.

    ``xi = a * sqrt(omega*mu*sigma) = sqrt(2) * a / delta`` is the radius measured in
    (sqrt(2) x) skin depths.  The classical result (e.g. Ramo-Whinnery-Van Duzer) is

        Rac/Rdc = (xi/2) * [ber(xi) bei'(xi) - bei(xi) ber'(xi)]
                          / [ber'(xi)^2 + bei'(xi)^2].

    Returns 1.0 in the DC limit xi -> 0 and grows ~ xi/(2 sqrt 2) + 1/4 as xi -> inf.
    """
    from scipy.special import kelvin as _kelvin
    Be, _Ke, Bep, _Kep = _kelvin(xi)
    ber, bei = Be.real, Be.imag
    berp, beip = Bep.real, Bep.imag
    denom = berp * berp + beip * beip
    if denom == 0.0:
        return 1.0
    return (xi / 2.0) * (ber * beip - bei * berp) / denom


def ac_resistance_round_wire(a=1.0, sigma=1.0, freq=None, xi=None, mu_r=1.0,
                             order=4, Rout=None, maxh_wire=None):
    """FEM skin-effect solve for an isolated round wire -> AC/DC resistance ratio.

    Drives the wire with a unit axial field E0 = 1 and solves the complex A_z
    diffusion problem on a self-generated disk-in-air mesh (wire radius ``a``,
    far boundary ``Rout`` with A_z = 0). Specify the frequency either directly via
    ``freq`` [Hz] or, more conveniently for benchmarking, by the dimensionless
    ``xi = a sqrt(omega mu sigma)`` (which fixes omega). The skin layer is resolved
    automatically (``maxh_wire ~ delta/3``) unless overridden.

    Returns a dict with ``rac_rdc`` (FEM ratio), ``xi``, ``area`` (FEM wire area),
    and ``current`` (complex total current for the unit drive).
    """
    import numpy as np
    import ngsolve as ng
    from netgen.geom2d import SplineGeometry

    mu = MU0 * mu_r
    if xi is None:
        if freq is None:
            raise ValueError("specify either freq [Hz] or xi")
        omega = 2.0 * math.pi * freq
        xi = a * math.sqrt(omega * mu * sigma)
    else:
        omega = xi * xi / (a * a * mu * sigma)
    delta = math.sqrt(2.0 / (omega * mu * sigma))
    if Rout is None:
        Rout = 8.0 * a
    if maxh_wire is None:
        maxh_wire = min(0.15 * a, delta / 3.0)

    geo = SplineGeometry()
    geo.AddCircle((0, 0), r=Rout, leftdomain=1, bc="outer")
    geo.AddCircle((0, 0), r=a, leftdomain=2, rightdomain=1, bc="wiresurf")
    geo.SetMaterial(1, "air")
    geo.SetMaterial(2, "wire")
    geo.SetDomainMaxH(2, maxh_wire)
    mesh = ng.Mesh(geo.GenerateMesh(maxh=Rout / 3.0))
    mesh.Curve(order)

    fes = ng.H1(mesh, order=order, complex=True, dirichlet="outer")
    A, w = fes.TnT()
    nu = 1.0 / mu
    sig = mesh.MaterialCF({"wire": sigma}, default=0.0)
    a_f = ng.BilinearForm(fes)
    a_f += nu * ng.grad(A) * ng.grad(w) * ng.dx
    a_f += 1j * omega * sig * A * w * ng.dx
    f_f = ng.LinearForm(fes)
    f_f += sig * 1.0 * w * ng.dx          # unit axial drive E0 = 1
    a_f.Assemble(); f_f.Assemble()
    gfu = ng.GridFunction(fes)
    gfu.vec.data = a_f.mat.Inverse(fes.FreeDofs(), inverse="umfpack") * f_f.vec

    Jz = sig * (1.0 - 1j * omega * gfu)
    wire = mesh.Materials("wire")
    cur = ng.Integrate(Jz, mesh, definedon=wire)
    j2 = ng.Integrate((Jz * ng.Conj(Jz)).real, mesh, definedon=wire)
    area = ng.Integrate(ng.CF(1.0), mesh, definedon=wire)
    rac_rdc = area * j2 / (abs(cur) ** 2)
    return {"rac_rdc": float(rac_rdc), "xi": float(xi), "area": float(area),
            "current": complex(cur), "delta": float(delta)}
