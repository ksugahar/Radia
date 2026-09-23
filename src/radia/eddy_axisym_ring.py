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

Formulation.  With ``psi = 2 pi r A_phi`` and a loop voltage ``V`` applied
around the ring, the source electric field is ``V / (2 pi r)`` and

    J_phi = sigma (V / (2 pi r) - j omega A_phi)
          = sigma / (2 pi r) * (V - j omega psi)

so the discrete system is

    (K + j omega M_sigma) psi = V * M_sigma * 1

where the right-hand side is the sigma-mass operator applied to the
coefficient vector of the constant function ``psi = 1``.  Building it that way
rather than by summing rows matters: the axihenrotte basis is polynomial in
``s = r^2``, so its shape functions are not a partition of unity and a row sum
would not be the interpolant of a constant.

The ring's impedance is ``V / I``.  Dividing by the mean circumference gives a
per-unit-length figure comparable with a straight conductor, which is exact in
the limit of a large radius and is how the straight-fin comparison is made.
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
    gf_psi: object
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
        applied = ng.CF(1.0 + 0j) - 1j * omega * self.gf_psi
        return (ng.CF(self.scale) * self.sigma_S_per_m
                / (2.0 * np.pi * ng.x) * applied)

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


def solve_axisym_ring(mesh, *, conductor, frequency_hz, sigma, order=2,
                      mu_r=1.0, dirichlet="axis|outer", solver="umfpack"):
    """Solve one voltage-driven axisymmetric eddy operating point.

    Args:
        mesh: two-dimensional NGSolve mesh in ``(r, z)``; ``r`` is the first
            coordinate and the axis is at ``r = 0``.
        conductor: material name of the conducting region.
        frequency_hz: drive frequency; must be positive.
        sigma: conductivity of ``conductor`` in S/m.
        order: axihenrotte element order.
        mu_r: relative permeability everywhere.
        dirichlet: boundary names where ``psi`` is pinned.  The axis belongs
            here: ``A_phi`` vanishes there by symmetry.

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
    region = mesh.Materials(conductor)
    if region.Mask().NumSet() == 0:
        raise ValueError(f"no elements carry the material name {conductor!r}")

    omega = 2.0 * np.pi * float(frequency_hz)
    delta = np.sqrt(2.0 / (omega * MU0 * float(mu_r) * float(sigma)))
    mu_cf = ng.CF(MU0 * float(mu_r))
    sigma_cf = mesh.MaterialCF({conductor: float(sigma)}, default=0.0)

    fes = H1Henrotte(mesh, order=order, complex=True, dirichlet=dirichlet)
    stiffness = ng.BilinearForm(fes, symmetric=True)
    stiffness += AxiHenrotteStiffnessBFI(mu_cf)
    mass = ng.BilinearForm(fes, symmetric=True)
    mass += AxiHenrotteSigmaMassBFI(sigma_cf)

    gf_one = ng.GridFunction(fes)
    gf_psi = ng.GridFunction(fes)
    stiffness.Assemble()
    mass.Assemble()
    # The interpolant of psi = 1, not a row sum: the basis is polynomial
    # in s = r^2 and is not a partition of unity.
    gf_one.Set(ng.CF(1.0 + 0j))
    rhs = gf_psi.vec.CreateVector()
    rhs.data = mass.mat * gf_one.vec
    system = stiffness.mat.CreateMatrix()
    system.AsVector().data = (stiffness.mat.AsVector()
                              + 1j * omega * mass.mat.AsVector())
    gf_psi.vec.data = system.Inverse(freedofs=fes.FreeDofs(),
                                     inverse=solver) * rhs

    from radia.axisym_measure import axi_section_integral

    radius = ng.x
    j_phi = sigma_cf / (2.0 * np.pi * radius) * (1.0 - 1j * omega * gf_psi)
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
        skin_depth_m=float(delta), gf_psi=gf_psi, mesh=mesh,
        conductor=conductor, scale=complex(scale))
