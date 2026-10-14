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
which builds on this core (validated to machine precision in validation_test/radia_mcp/test_airgap_eddy_machine.py).

Open method (Abdel-Razek/Konrad 1982; Davat 1985 air-gap macro-element); analytic,
publishable.
"""
from __future__ import annotations

import cmath
import math

MU0 = 4e-7 * math.pi


def annular_rotation_phase(harmonic, angle_rad):
    """Return the rotor-trace phase for harmonic ``n`` after rotation ``angle_rad``.

    With ``A_n(theta) ~ exp(i*n*theta)``, replacing ``theta`` by
    ``theta-angle_rad`` multiplies the rotor trace by ``exp(-i*n*angle_rad)``.
    The AGE mesh and DtN operator are unchanged.
    """

    if isinstance(harmonic, bool) or int(harmonic) != harmonic or int(harmonic) <= 0:
        raise ValueError("AGE harmonic must be a positive integer")
    angle = float(angle_rad)
    if not math.isfinite(angle):
        raise ValueError("AGE rotation angle must be finite")
    return cmath.exp(-1j * int(harmonic) * angle)


def planar_translation_phase(wavenumber, displacement):
    """Return the moving-trace phase for a planar AGE translation.

    With ``A_k(x) ~ exp(i*k*x)``, replacing ``x`` by ``x-displacement``
    multiplies the moving trace by ``exp(-i*k*displacement)`` without remeshing.
    """

    wave = float(wavenumber)
    shift = float(displacement)
    if not math.isfinite(wave) or wave <= 0.0:
        raise ValueError("AGE wavenumber must be positive and finite")
    if not math.isfinite(shift):
        raise ValueError("AGE displacement must be finite")
    return cmath.exp(-1j * wave * shift)


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


# ---------------------------------------------------------------------------
# Planar gap (linear motor / axial-flux machine -- the planar analogue of the annular AGE)
# ---------------------------------------------------------------------------
#
# In the current-free planar gap 0 < z < g each transverse Fourier mode with
# wavenumber k satisfies Laplace → d²A/dz² = k²A, so
#
#     A(z)  =  a sinh(k z) + b sinh(k (g-z)),    S = sinh(k g)
#     A(0) = b S  →  b = A_bottom / S
#     A(g) = a S  →  a = A_top    / S
#
# DtN (Steklov) matrix M,  [dA/dz|_0 ; dA/dz|_g] = M @ [A_bottom ; A_top]:
#
#     M = (k/S) [[ -cosh(k g),  1          ],
#                [ -1,          cosh(k g)  ]]
#
# Thrust force (spatial average over one period in x, on the bottom plane z=0):
#
#     F_x = (k² x_period y_depth / (2 μ₀ S)) Im(A_bottom conj(A_top))
#
# POSITION-INDEPENDENT: verified by maxwell_thrust_numeric at every z ∈ [0, g]
# (the planar gem, exactly analogous to the radius-independence of the rotary torque).
#
# Correspondence with the annular AGE:
#   r^n  ↔  sinh(k z)         (basis function in the gap)
#   angular index n  ↔  wavenumber k = n π / W (or 2π n / W)
#   rotary torque T  ↔  linear thrust F_x
#
# Open method (Abdel-Razek/Konrad 1982 generalised); analytic, publishable.


def planar_harmonic_coeffs(k, g, A_bottom, A_top):
    """Complex (a, b) with A(z) = a sinh(kz) + b sinh(k(g-z)), normalised by S = sinh(kg).

    From A(0) = b S = A_bottom and A(g) = a S = A_top:
        a = A_top    / sinh(k g)
        b = A_bottom / sinh(k g)"""
    S = cmath.sinh(k * g)
    return A_top / S, A_bottom / S


def planar_field(k, g, A_bottom, A_top, z):
    """The planar harmonic potential A(z) [real] at gap position z ∈ [0, g]."""
    a, b = planar_harmonic_coeffs(k, g, A_bottom, A_top)
    return (a * cmath.sinh(k * z) + b * cmath.sinh(k * (g - z))).real


def planar_radial_trace(k, g, A_bottom, A_top):
    """DtN action: (dA/dz|_{z=0}, dA/dz|_{z=g}) from (A_bottom, A_top).

    dA/dz|_0 = k (a - b cosh(kg)) = k/S (A_top - A_bottom cosh(kg))
    dA/dz|_g = k (a cosh(kg) - b) = k/S (A_top cosh(kg) - A_bottom)"""
    a, b = planar_harmonic_coeffs(k, g, A_bottom, A_top)
    ckg = cmath.cosh(k * g)
    d_bottom = k * (a - b * ckg)
    d_top    = k * (a * ckg - b)
    return d_bottom, d_top


def planar_dtn_matrix(k, g):
    """2×2 Dirichlet-to-Neumann (Steklov) matrix M for the planar gap,

        [dA/dz|_0 ; dA/dz|_g]  =  M @ [A_bottom ; A_top]

        M = (k/S) [[ -cosh(kg),  1         ],
                   [ -1,          cosh(kg) ]]

    Returned as ((M11, M12), (M21, M22)).  Direct planar analogue of
    :func:`annular_dtn_matrix`; assembles into the NGSolve stiffness exactly
    as the annular DtN does (same Woodbury low-rank structure).

    The bottom-plane row (M[0]) borders the lower FE region; the top-plane
    row (M[1]) borders the upper FE region."""
    d_b_in,  d_t_in  = planar_radial_trace(k, g, 1.0, 0.0)  # A_bottom=1, A_top=0
    d_b_out, d_t_out = planar_radial_trace(k, g, 0.0, 1.0)  # A_bottom=0, A_top=1
    return ((d_b_in, d_b_out), (d_t_in, d_t_out))


def planar_dtn_modes(g, wavenumbers):
    """Bottom-plane (M[0][0], M[0][1]) of the gap DtN for each wavenumber k in
    ``wavenumbers``.  Returns ``{k: (M01, M00)}`` (coupling, self-term) -- direct
    analogue of :func:`airgap_dtn_modes`."""
    return {k: planar_dtn_matrix(k, g)[0] for k in wavenumbers}


def planar_harmonic_thrust(k, g, A_bottom, A_top, x_period, y_depth=1.0, mu0=MU0):
    """Position-independent thrust force [N] on the bottom plane (z=0) in the +x direction,
    over the area x_period × y_depth:

        F_x = (k² x_period y_depth / (2 μ₀ sinh(kg)))  Im(A_bottom conj(A_top)).

    POSITION-INDEPENDENT: the same force results at every z in the gap (verified by
    :func:`maxwell_thrust_numeric`).  Sign: F_x > 0 when Im(A_bottom conj(A_top)) > 0,
    i.e. when A_top leads A_bottom in phase by an angle δ ∈ (π, 2π) (or lags by δ ∈ (0,π)).
    This is the exact planar analogue of :func:`airgap_harmonic_torque`."""
    S = math.sinh(k * g)
    return (k * k * x_period * y_depth / (2.0 * mu0 * S)) * (A_bottom * A_top.conjugate()).imag


def maxwell_thrust_numeric(k, g, A_bottom, A_top, z, x_period, y_depth=1.0, mu0=MU0):
    """Independent check: Maxwell-stress thrust at height z ∈ [0, g] via closed-form
    spatial average.

    <T_xz>_x = (k / 2μ₀) Im[A(z) conj(dA/dz)(z)]

    and the total force F_x = <T_xz>_x × x_period × y_depth.  This must equal
    :func:`planar_harmonic_thrust` at EVERY z (the planar gem / position-independence).
    Uses the exact phasor formula (no quadrature), so the comparison is analytic."""
    a, b = planar_harmonic_coeffs(k, g, A_bottom, A_top)
    Az  = a * cmath.sinh(k * z) + b * cmath.sinh(k * (g - z))
    dAz = k * (a * cmath.cosh(k * z) - b * cmath.cosh(k * (g - z)))
    # <Bx Bz>_x = (k/2) Im[Az dAz*], force on bottom = <T_xz> = (1/mu0)<Bx Bz>
    return (k * x_period * y_depth / (2.0 * mu0)) * (Az * dAz.conjugate()).imag
