# -*- coding: utf-8 -*-
r"""Air-gap harmonic element (AGE) -- analytic core + gate.

The mesh-free annular coupling for a rotating-machine air gap.  In the current-free
gap (ri < r < ro) the 2D vector potential is harmonic, so each Fourier mode is

    A_n(r, theta) = Re[ (alpha r^n + beta r^-n) e^{i n theta} ]   (Laplace, annulus)

with complex coefficients (alpha, beta) fixed by the rotor-ring phasor A_n(ri)=A_in
and the stator-ring phasor A_n(ro)=A_out.  Coupling the two FE regions through this
ANALYTIC transfer (instead of meshing the gap) gives the AGE's two gems:

  * the gap is never discretised -- rotor rotation is a phase shift of A_in;
  * the transmitted torque is computed from the gap harmonics in CLOSED FORM and is
    RADIUS-INDEPENDENT (mesh-independent), unlike a Maxwell-stress contour.

This module is the analytic CORE: the harmonic transfer, the radial-trace
(Dirichlet-to-Neumann / Steklov) action that borders the NGSolve system, and the
closed-form harmonic torque -- all validated against direct Maxwell-stress quadrature
in tests/test_airgap_element.py.  The reusable NGSolve assembly that couples a ROTOR and a
STATOR FE region across the un-meshed gap -- real OR eddy-current (jw*mu*sigma), with rotor
rotation as a pure phase and the mesh-free closed-form torque -- is :mod:`airgap_machine`,
which builds on this core (validated to machine precision in tests/test_airgap_eddy_machine.py).

Open method (Abdel-Razek/Konrad 1982; Davat 1985 air-gap macro-element); analytic,
publishable.
"""
from __future__ import annotations

import cmath
import math

MU0 = 4e-7 * math.pi


def annular_harmonic_coeffs(n, ri, ro, A_in, A_out):
    """Complex (alpha, beta) of A_n(r) = alpha r^n + beta r^-n in the current-free
    annulus ri<r<ro, from the ring phasors A_n(ri)=A_in (rotor), A_n(ro)=A_out (stator).
    Solves [[ri^n, ri^-n],[ro^n, ro^-n]] [alpha; beta] = [A_in; A_out]."""
    rin, rio = ri ** n, ri ** -n
    ron, roo = ro ** n, ro ** -n
    det = rin * roo - rio * ron
    alpha = (A_in * roo - A_out * rio) / det
    beta = (A_out * rin - A_in * ron) / det
    return alpha, beta


def annular_field(n, ri, ro, A_in, A_out, r, theta):
    """The harmonic-n potential A_n(r, theta) [real] at a gap point (analytic)."""
    alpha, beta = annular_harmonic_coeffs(n, ri, ro, A_in, A_out)
    return ((alpha * r ** n + beta * r ** -n) * cmath.exp(1j * n * theta)).real


def annular_radial_trace(n, ri, ro, A_in, A_out):
    """The Dirichlet-to-Neumann (Steklov) ACTION of the gap for harmonic n: given the
    ring potentials (A_in, A_out), return the radial derivatives (dA/dr|ri, dA/dr|ro)
    -- the Neumann data the gap imposes back on the two FE regions (the AGE coupling
    that will border the NGSolve stiffness).  dA_n/dr = n alpha r^{n-1} - n beta r^{-n-1}."""
    alpha, beta = annular_harmonic_coeffs(n, ri, ro, A_in, A_out)
    d_ri = n * alpha * ri ** (n - 1) - n * beta * ri ** (-n - 1)
    d_ro = n * alpha * ro ** (n - 1) - n * beta * ro ** (-n - 1)
    return d_ri, d_ro


def annular_dtn_matrix(n, ri, ro):
    """The 2x2 Dirichlet-to-Neumann (Steklov) matrix M of the gap for harmonic n,

        [dA/dr|ri ; dA/dr|ro] = M @ [A_in ; A_out] ,

    returned as ((M11, M12), (M21, M22)).  The stator-ring row (M21, M22) is what
    borders the NGSolve stiffness in the AGE assembly: M22 is the Robin self-term on
    the gap-facing boundary (coefficient of A there), M21 couples the rotor-ring
    phasor A_in into that boundary's Neumann data.  Linear, so read off from the two
    unit responses of :func:`annular_radial_trace`."""
    d_ri_in, d_ro_in = annular_radial_trace(n, ri, ro, 1.0, 0.0)   # A_in=1, A_out=0
    d_ri_out, d_ro_out = annular_radial_trace(n, ri, ro, 0.0, 1.0)  # A_in=0, A_out=1
    return ((d_ri_in, d_ri_out), (d_ro_in, d_ro_out))


def airgap_dtn_modes(ri, ro, harmonics):
    """Stator-ring (M21, M22) of the gap DtN for each harmonic n in ``harmonics``.

    The gap couples each Fourier mode INDEPENDENTLY (it is current-free / linear), so
    the multi-harmonic AGE boundary block is the sum over n of the rank-1 spectral
    projector ``(M22(n)/norm) c_n c_n^T`` (c_n = the boundary integral of cos/sin n*theta),
    with rotor-ring forcing ``-M21(n) A_in,n c_n``.  This is what lets a SINGLE FE solve
    carry the full air-gap spectrum (slot / pole harmonics) across the un-meshed gap.
    Returns ``{n: (M21, M22)}``."""
    return {n: annular_dtn_matrix(n, ri, ro)[1] for n in harmonics}


def airgap_harmonic_torque(n, ri, ro, A_in, A_out, axial_length=1.0, mu0=MU0):
    """Mesh-independent torque [N*m] transmitted across the gap by harmonic n:

        T_n = -(2 pi n^2 L / mu0) * Im(alpha * conj(beta)).

    This is RADIUS-INDEPENDENT (the torque is conserved across the current-free gap --
    the AGE gem), depends only on the rotor/stator phasor PHASE (zero when in phase,
    maximal at 90 deg), and needs no Maxwell-stress contour or air-gap mesh."""
    alpha, beta = annular_harmonic_coeffs(n, ri, ro, A_in, A_out)
    return -(2.0 * math.pi * n * n * axial_length / mu0) * (alpha * beta.conjugate()).imag


def maxwell_torque_numeric(n, ri, ro, A_in, A_out, r, axial_length=1.0, mu0=MU0, npts=720):
    """Independent check: direct weighted Maxwell-stress torque at radius r by
    theta-quadrature, T(r) = (r^2 L / mu0) * integral_0^2pi B_r B_theta dtheta, with
    B_r = (1/r) dA/dtheta, B_theta = -dA/dr.  The harmonic formula must equal this at
    EVERY r in the gap (the mesh-independence the test locks)."""
    alpha, beta = annular_harmonic_coeffs(n, ri, ro, A_in, A_out)
    Ahat = alpha * r ** n + beta * r ** -n
    dAhat = n * alpha * r ** (n - 1) - n * beta * r ** (-n - 1)
    Br_hat = 1j * n / r * Ahat            # phasor of B_r = (1/r) dA/dtheta
    Bth_hat = -dAhat                      # phasor of B_theta = -dA/dr
    dth = 2.0 * math.pi / npts
    acc = 0.0
    for k in range(npts):
        e = cmath.exp(1j * n * k * dth)
        acc += (Br_hat * e).real * (Bth_hat * e).real * dth
    return r * r * axial_length / mu0 * acc
