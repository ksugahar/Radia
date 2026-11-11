"""Time-harmonic eddy currents in a voltage-driven axisymmetric ring.

The beak fin's delivery metrics are mid-span quantities of a section that does
not change along the conductor.  Revolving that section about an axis gives a
ring whose meridian IS the beak profile, so the skin physics under test is
present in a two-dimensional problem -- and a two-dimensional mesh can resolve
a 0.17 mm skin over an 8 mm section without the element count that defeats the
three-dimensional solve.

This builds on `radia.axifem`, whose stiffness and sigma-mass operators are
already validated against the Kameari sphere, the Hiruma disk and Stoll.  What
is added here is the drive and the terminal quantities: those operators solve
the induced problem, and a fin carries an impressed current.

Formulation.  The axihenrotte unknown is ``A_phi`` itself, and both
operators carry the revolved volume measure ``2 pi r dr dz``:

    K_ij = int nu curl N_i . curl N_j 2 pi r dr dz,
    M_ij = int sigma N_i N_j 2 pi r dr dz.

A loop voltage ``V`` around the ring is the source field ``V / (2 pi r)``, so

    J_phi = sigma (V / (2 pi r) - j omega A_phi)
    (K + j omega M) A = V * int sigma N_i dr dz

-- the ``2 pi r`` of the measure cancels the ``1 / (2 pi r)`` of the source.
The drive is NOT ``M`` applied to a constant: that would impose a uniform
``E_phi``, which no loop voltage produces, and would leave the field and the
reconstructed current describing different problems.

The ring's impedance is ``V / I``.  Dividing by the mean circumference gives a
per-unit-length figure comparable with a straight conductor, which is exact in
the limit of a large radius and is how the straight-fin comparison is made.

Exterior.  The mesh ends on a semicircle ``outer`` of radius ``R`` centred on
the axis at ``z = z_c``; the axis is the boundary ``axis``.  Outside it the
field is source-free, ``A_phi = sum_n a_n rho^-(n+1) P_n^1(cos theta)``, and
eliminating the exterior gives the exact Dirichlet-to-Neumann term

    a_ext(u, v) = sum_n c_n g_n(u) g_n(v),
    g_n(u) = int_outer r^2 / R dP_n/dx(x) u ds / (N_n R^2),  x = (z - z_c)/R,
    c_n = 2 pi R n N_n / mu,  N_n = 2 n (n + 1) / (2 n + 1),

a rank-``modes`` update applied by the Woodbury identity on one factorisation.
A pinned box in its place overstates the inductance's sensitivity to the
truncation (+3.5 % at three box radii); the exterior term makes the result
independent of ``R`` up to the mode cut and the mesh.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

MU0 = 4e-7 * np.pi


@dataclass(frozen=True)
class AxisymRingResult:
    """One solved operating point, rescaled to unit terminal current."""

    frequency_hz: float
    sigma_S_per_m: float
    terminal_current_A: complex
    terminal_impedance_ohm: complex
    resistance_ohm: float
    inductance_H: float
    ohmic_loss_W: float
    mean_radius_m: float
    ndof: int
    skin_depth_m: float
    gf_a: object
    mesh: object
    conductor: str
    scale: complex

    @property
    def impedance_per_metre(self):
        """Ring impedance divided by the mean circumference."""
        return self.terminal_impedance_ohm / (2.0 * np.pi * self.mean_radius_m)

    def current_density(self):
        """``J_phi`` as a CoefficientFunction, at unit terminal current."""
        import ngsolve as ng

        omega = 2.0 * np.pi * self.frequency_hz
        # V = 1 was applied; `scale` carries the whole solution to unit
        # terminal current, the drive included.
        return (ng.CF(self.scale) * self.sigma_S_per_m
                * (1.0 / (2.0 * np.pi * ng.x) - 1j * omega * self.gf_a))

    def loss_density(self):
        """``|J|^2 / (2 sigma)``: time-average watts per cubic metre."""
        import ngsolve as ng

        j = self.current_density()
        return (j * ng.Conj(j)).real / (2.0 * self.sigma_S_per_m)

    def total_loss(self):
        """Dissipated power in the whole ring, in watts.

        Integrated in the revolved volume measure -- see
        :mod:`radia.axisym_measure` for why that has to be said out loud.
        """
        from radia.axisym_measure import axi_volume_integral

        return float(axi_volume_integral(
            self.loss_density(), self.mesh,
            definedon=self.mesh.Materials(self.conductor)).real)

    def loss_beyond(self, r_m):
        """Dissipated power in the part of the ring beyond radius ``r_m``."""
        from radia.axisym_measure import axi_volume_beyond

        return float(axi_volume_beyond(
            self.loss_density(), self.mesh, r_m,
            definedon=self.mesh.Materials(self.conductor)).real)


_MAX_EXTERIOR_MODES = 40
_EXTERIOR_TAIL_LIMIT = 1e-6


def _legendre_derivatives(x, modes):
    """``dP_n/dx`` for n = 1..modes as Horner chains in the CF ``x``.

    A three-term recurrence would share subexpressions that NGSolve's
    CoefficientFunction tree re-evaluates, which grows exponentially in n.
    The monomial coefficients grow like 2^n, which is why the mode count is
    capped and a distant boundary is preferred to a high mode count.
    """
    import ngsolve as ng
    from numpy.polynomial import legendre

    out = []
    for n in range(1, modes + 1):
        coeffs = legendre.leg2poly(legendre.legder([0] * n + [1]))
        cf = ng.CF(float(coeffs[-1]))
        for a in coeffs[-2::-1]:
            cf = cf * x + float(a)
        out.append(cf)
    return out


_ARC_TOLERANCE = 1e-5


def _check_exterior_arc(mesh, radius, center_z):
    """Fail unless the mapped ``outer`` trace is the stated semicircle.

    Vertices alone are not enough: a chordal (``Curve(1)``) mesh has the same
    vertices and a different boundary.  The check integrates over the mapped
    trace, so it sees the curved geometry the solve integrates on: the trace
    must cover half the circumference, stay on the circle in the mean-square
    sense, and end on the axis at both poles.
    """
    import ngsolve as ng

    names = set(mesh.GetBoundaries())
    if "outer" not in names:
        raise ValueError("the mesh has no boundary named 'outer'")
    if "axis" not in names:
        raise ValueError("the mesh has no boundary named 'axis'")
    outer = mesh.Boundaries("outer")
    order = 2 * max(int(mesh.GetCurveOrder()), 1) + 4
    length = float(ng.Integrate(ng.CF(1.0), mesh, ng.BND, definedon=outer,
                                order=order))
    rho = ng.sqrt(ng.x * ng.x + (ng.y - center_z) * (ng.y - center_z))
    miss = float(ng.Integrate((rho - radius) ** 2, mesh, ng.BND,
                              definedon=outer, order=order))
    rms = (max(miss, 0.0) / max(length, 1e-300)) ** 0.5 / radius
    cover = length / (np.pi * radius) - 1.0
    poles = []
    for el in mesh.Elements(ng.BND):
        if el.mat == "axis":
            poles.extend(mesh[v].point for v in el.vertices)
    poles = np.array(poles, dtype=float)
    reach = (float(np.max(np.abs(poles[:, 1] - center_z))) / radius - 1.0
             if poles.size else -1.0)
    if (rms > _ARC_TOLERANCE or abs(cover) > _ARC_TOLERANCE
            or abs(reach) > _ARC_TOLERANCE
            or float(np.max(np.abs(poles[:, 0]))) > _ARC_TOLERANCE * radius):
        raise ValueError(
            f"'outer' is not the mapped semicircle of radius {radius} about "
            f"z = {center_z}: rms miss {rms:.2e}, coverage {cover:+.2e}, "
            f"axis reach {reach:+.2e} (tolerance {_ARC_TOLERANCE:.0e}); "
            "a chordal or partial boundary is not the exterior the DtN term "
            "assumes -- curve the mesh")


def solve_axisym_ring(mesh, *, conductor, frequency_hz, sigma,
                      exterior_radius, exterior_center_z=0.0,
                      exterior_modes=30, order=2, mu_r=1.0,
                      solver="umfpack"):
    """Solve one voltage-driven axisymmetric eddy operating point.

    Args:
        mesh: two-dimensional NGSolve mesh in ``(r, z)``; ``r`` is the first
            coordinate and the axis is at ``r = 0``.
        conductor: material name of the conducting region.
        frequency_hz: drive frequency; must be positive.
        sigma: conductivity of ``conductor`` in S/m.
        order: axihenrotte element order.
        exterior_radius: radius of the semicircular boundary ``outer``.
        exterior_center_z: axial position of its centre.
        exterior_modes: number of exterior multipoles kept.
        mu_r: relative permeability everywhere, exterior included.

    Returns:
        AxisymRingResult rescaled so the terminal current is exactly 1 A.
    """
    import ngsolve as ng

    from radia.axifem import (AxiHenrotteSigmaMassBFI,
                              AxiHenrotteStiffnessBFI, H1Henrotte)

    if frequency_hz <= 0:
        raise ValueError("frequency_hz must be positive; "
                         "this is a time-harmonic formulation")
    if sigma <= 0:
        raise ValueError("sigma must be positive")
    if mesh.dim != 2:
        raise ValueError("the axisymmetric solve needs a 2-D (r, z) mesh")
    if not 1 <= exterior_modes <= _MAX_EXTERIOR_MODES:
        raise ValueError(f"exterior_modes must be in 1..{_MAX_EXTERIOR_MODES}; "
                         "the monomial form loses digits beyond that, so "
                         "move 'outer' outward instead")
    _check_exterior_arc(mesh, float(exterior_radius), float(exterior_center_z))
    region = mesh.Materials(conductor)
    if region.Mask().NumSet() == 0:
        raise ValueError(f"no elements carry the material name {conductor!r}")

    omega = 2.0 * np.pi * float(frequency_hz)
    delta = np.sqrt(2.0 / (omega * MU0 * float(mu_r) * float(sigma)))
    mu_cf = ng.CF(MU0 * float(mu_r))
    sigma_cf = mesh.MaterialCF({conductor: float(sigma)}, default=0.0)

    # The axis is pinned by symmetry; 'outer' is natural and closed below.
    fes = H1Henrotte(mesh, order=order, complex=True, dirichlet="axis")
    stiffness = ng.BilinearForm(fes, symmetric=True)
    stiffness += AxiHenrotteStiffnessBFI(mu_cf)
    mass = ng.BilinearForm(fes, symmetric=True)
    mass += AxiHenrotteSigmaMassBFI(sigma_cf)

    gf_a = ng.GridFunction(fes)
    stiffness.Assemble()
    mass.Assemble()
    test = fes.TestFunction()
    drive = ng.LinearForm(fes)
    drive += sigma_cf * test * ng.dx      # V = 1; 2 pi r cancels 1/(2 pi r)
    drive.Assemble()
    system = stiffness.mat.CreateMatrix()
    system.AsVector().data = (stiffness.mat.AsVector()
                              + 1j * omega * mass.mat.AsVector())
    inverse = system.Inverse(freedofs=fes.FreeDofs(), inverse=solver)
    gf_a.vec.data = inverse * drive.vec

    # Exterior: Woodbury on sum_n c_n g_n g_n^T.
    big_r = float(exterior_radius)
    derivs = _legendre_derivatives((ng.y - float(exterior_center_z)) / big_r,
                                   int(exterior_modes))
    g_vecs, y_vecs, c_inv = [], [], []
    for n, dp in enumerate(derivs, start=1):
        norm_n = 2.0 * n * (n + 1) / (2 * n + 1)
        form = ng.LinearForm(fes)
        form += (ng.x * ng.x / big_r * dp / (norm_n * big_r * big_r)
                 * test.Trace() * ng.ds("outer"))
        form.Assemble()
        y = form.vec.CreateVector()
        y.data = inverse * form.vec
        g_vecs.append(form.vec)
        y_vecs.append(y)
        c_inv.append(MU0 * float(mu_r) / (2.0 * np.pi * big_r * n * norm_n))
    k = len(g_vecs)
    capacitance = np.diag(np.asarray(c_inv, dtype=complex))
    for i in range(k):
        for j in range(k):
            capacitance[i, j] += ng.InnerProduct(y_vecs[j], g_vecs[i],
                                                 conjugate=False)
    projections = np.array([ng.InnerProduct(gf_a.vec, g, conjugate=False)
                            for g in g_vecs])
    weights = np.linalg.solve(capacitance, projections)
    for w, y in zip(weights, y_vecs):
        gf_a.vec.data -= complex(w) * y
    # The multipole series converges as (r_source / R)^n; a cut that leaves
    # energy in the last mode is a truncation, not an exterior.
    final = np.array([ng.InnerProduct(gf_a.vec, g, conjugate=False)
                      for g in g_vecs])
    mode_energy = np.abs(final) ** 2 / np.asarray(c_inv)
    tail = float(mode_energy[-1] / mode_energy.sum())
    if tail > _EXTERIOR_TAIL_LIMIT:
        raise RuntimeError(
            f"the exterior series is not converged: mode {k} holds "
            f"{tail:.2e} of the exterior energy (limit "
            f"{_EXTERIOR_TAIL_LIMIT:.0e}); move 'outer' away from the "
            "conductor or keep more modes")

    from radia.axisym_measure import axi_section_integral

    radius = ng.x
    j_phi = sigma_cf * (1.0 / (2.0 * np.pi * radius) - 1j * omega * gf_a)
    # Terminal current and section geometry are meridian-measure quantities:
    # the current flows THROUGH this cross-section.  Dissipation is not, and
    # goes through the volume measure instead.
    current = complex(axi_section_integral(j_phi, mesh, definedon=region))
    if current == 0:
        raise RuntimeError("the solution carries no ring current")
    area = float(axi_section_integral(ng.CF(1.0), mesh,
                                      definedon=region).real)
    mean_radius = float(axi_section_integral(
        radius, mesh, definedon=region).real) / area

    scale = 1.0 / current
    impedance = 1.0 / current          # V = 1 was applied
    return AxisymRingResult(
        frequency_hz=float(frequency_hz), sigma_S_per_m=float(sigma),
        terminal_current_A=complex(current),
        terminal_impedance_ohm=complex(impedance),
        resistance_ohm=float(impedance.real),
        inductance_H=float(impedance.imag / omega),
        ohmic_loss_W=float(0.5 * impedance.real),
        mean_radius_m=mean_radius, ndof=int(fes.ndof),
        skin_depth_m=float(delta), gf_a=gf_a, mesh=mesh,
        conductor=conductor, scale=complex(scale))
