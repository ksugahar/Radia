# -*- coding: utf-8 -*-
r"""Air-gap harmonic element (AGE) -- the reusable rotating-machine FE solver, eddy-capable.

Builds on the analytic core :mod:`airgap_element` (harmonic DtN + closed-form torque).  This
module packages the NGSolve assembly that couples a ROTOR FE region and a STATOR FE region
across the UN-MESHED air gap, for a real OR complex (eddy-current) machine, with rotor rotation
as a pure phase and a mesh-free closed-form torque:

    coupling = airgap_coupling(fes, ri, ro, "rotor_ring", "stator_ring", harmonics)
    factor   = airgap_factorize(a.mat, coupling, fes.FreeDofs())   # K fixed across rotation
    gfu      = airgap_solve(factor, fes, dirichlet_cf=excite(theta_r))   # one cheap solve / angle
    T        = airgap_torque(coupling, gfu, axial_length=L)        # radius-independent

The gap (ri<r<ro) is current-free, so even in the time-harmonic (eddy) regime it stays Laplace
and the SAME real 2x2 DtN couples a CONDUCTING (jw*mu*sigma) FE region across the un-meshed gap.
Each Fourier harmonic n couples independently as the low-rank block ``U G U^T`` (U = the cos/sin
ring-mode boundary vectors), added to the FE stiffness by the Woodbury identity -- so one
factorization serves a whole rotation sweep.  cos AND sin modes (rank-4 per harmonic) carry
arbitrary rotor phase; rotation is a phase of the rotor excitation, no remesh.

This is the eddy x mesh-free-gap rotating-machine core that a magnetostatic air-gap element
cannot give and that mainstream codes only reach with a meshed sliding band.  Validated to
machine precision against a fully-meshed complex reference (validation_test/radia_mcp/test_airgap_eddy_machine.py).

Open method (Abdel-Razek/Konrad 1982; Davat 1985 air-gap macro-element).

IMPORTANT (complex/eddy): NGSolve ``InnerProduct(a, b)`` is HERMITIAN (= a . conj(b)).  The
Woodbury projections here are bilinear (a . b, no conjugation), so every InnerProduct uses
``conjugate=False`` -- harmless when real, but conjugating it would silently strip the eddy
phase from the coupled solution.
"""
from __future__ import annotations

import math

import numpy as np
from ngsolve import GridFunction, LinearForm, ds, x, y, atan2, cos, sin, BND, InnerProduct

from .airgap_element import annular_dtn_matrix, airgap_harmonic_torque


def _inv(mat, freedofs):
    for solver in ("umfpack", "pardiso", "sparsecholesky"):
        try:
            return mat.Inverse(freedofs, inverse=solver)
        except Exception:
            continue
    return mat.Inverse(freedofs)


def airgap_coupling(fes, ri, ro, rotor_ring, stator_ring, harmonics):
    """Assemble the AGE gap coupling for the un-meshed annulus ri<r<ro bordered by the
    ``rotor_ring`` (r=ri) and ``stator_ring`` (r=ro) boundaries of ``fes``.

    Returns a dict with the per-harmonic ring-mode vectors ``U`` (4 per harmonic, in the order
    rotor-cos, stator-cos, rotor-sin, stator-sin), the block-diagonal inverse coupling ``Ginv``,
    and the ring norms.  Each harmonic n couples independently via its 2x2 DtN
    (:func:`airgap_element.annular_dtn_matrix`) as the rank-4 block ``U G U^T``; the cos & sin
    modes let the rotor harmonic sit at any phase (rotation).  Works for real or complex ``fes``."""
    mesh = fes.mesh
    v = fes.TestFunction()
    norm_ri, norm_ro = math.pi * ri, math.pi * ro
    U, Gblocks = [], []
    for n in harmonics:
        cn, sn = cos(n * atan2(y, x)), sin(n * atan2(y, x))
        for cf, bc in ((cn, rotor_ring), (cn, stator_ring), (sn, rotor_ring), (sn, stator_ring)):
            L = LinearForm(fes)
            L += cf * v * ds(definedon=mesh.Boundaries(bc))
            L.Assemble()
            vec = L.vec.CreateVector(); vec.data = L.vec
            U.append(vec)
        (M11, M12), (M21, M22) = annular_dtn_matrix(n, ri, ro)
        G = np.array([[-M11 / norm_ri, -M12 / norm_ro], [M21 / norm_ri, M22 / norm_ro]])
        Gfull = np.zeros((4, 4)); Gfull[0:2, 0:2] = G; Gfull[2:4, 2:4] = G   # cos block + sin block
        Gblocks.append(Gfull)
    nh = len(list(harmonics))
    Gbig = np.zeros((4 * nh, 4 * nh))
    for k, Gf in enumerate(Gblocks):
        Gbig[4 * k:4 * k + 4, 4 * k:4 * k + 4] = Gf
    return {"harmonics": list(harmonics), "ri": ri, "ro": ro, "U": U,
            "Ginv": np.linalg.inv(Gbig), "norm_ri": norm_ri, "norm_ro": norm_ro}


def airgap_factorize(Kmat, coupling, freedofs):
    """Factor the AGE-coupled operator ``K + U G U^T`` ONCE (K is fixed across rotor rotation).
    Returns the Woodbury pieces (Kinv, w = Kinv U, smallinv = (Ginv + U^T Kinv U)^-1) reused by
    every per-angle :func:`airgap_solve`.  ``conjugate=False`` keeps the projections bilinear so
    the complex (eddy) case is exact."""
    U = coupling["U"]
    Kinv = _inv(Kmat, freedofs)
    w = [u.CreateVector() for u in U]
    for i, u in enumerate(U):
        w[i].data = Kinv * u
    UtKU = np.array([[InnerProduct(U[i], w[j], conjugate=False) for j in range(len(U))]
                     for i in range(len(U))])
    smallinv = np.linalg.inv(coupling["Ginv"] + UtKU)
    return {"Kmat": Kmat, "Kinv": Kinv, "w": w, "smallinv": smallinv, "U": U}


def airgap_solve(factor, fes, dirichlet_cf=None, source_lf=None):
    """Solve ``(K + U G U^T) u = f`` for one rotor position (Dirichlet lift + optional volume
    source ``source_lf``).  Reuses the cached factorization, so a rotation sweep is one
    :func:`airgap_factorize` + many cheap solves.  Returns the ``GridFunction``."""
    U, w, smallinv = factor["U"], factor["w"], factor["smallinv"]
    gfu = GridFunction(fes)
    if dirichlet_cf is not None:
        gfu.Set(dirichlet_cf, BND)
    rhs = gfu.vec.CreateVector()
    rhs.data = -(factor["Kmat"] * gfu.vec)          # Dirichlet lift
    if source_lf is not None:
        rhs.data += source_lf.vec
    Kr = rhs.CreateVector(); Kr.data = factor["Kinv"] * rhs
    beta = smallinv @ np.array([InnerProduct(U[i], Kr, conjugate=False) for i in range(len(U))])
    delta = Kr.CreateVector(); delta.data = Kr
    for i in range(len(U)):
        coefficient = complex(beta[i])
        if not delta.is_complex:
            scale = max(1.0, abs(coefficient.real))
            if abs(coefficient.imag) > 1.0e-12 * scale:
                raise ValueError("real AGE space received a complex Woodbury coefficient")
            coefficient = float(coefficient.real)
        delta.data -= coefficient * w[i]
    gfu.vec.data += delta
    return gfu


def airgap_ring_phasors(coupling, gfu):
    """Recover the complex gap-ring harmonic phasors ``{n: (a_RI, a_RO)}`` from an FE solution
    (cos - i*sin amplitudes on the rotor/stator rings).  These are the un-meshed gap's interface
    data -- what the closed-form torque and any post-processing read."""
    U, nri, nro = coupling["U"], coupling["norm_ri"], coupling["norm_ro"]
    out = {}
    for k, n in enumerate(coupling["harmonics"]):
        rc, sc, rs, ss = U[4 * k], U[4 * k + 1], U[4 * k + 2], U[4 * k + 3]
        a_RI = (InnerProduct(rc, gfu.vec, conjugate=False) - 1j * InnerProduct(rs, gfu.vec, conjugate=False)) / nri
        a_RO = (InnerProduct(sc, gfu.vec, conjugate=False) - 1j * InnerProduct(ss, gfu.vec, conjugate=False)) / nro
        out[n] = (a_RI, a_RO)
    return out


def airgap_torque(coupling, gfu, axial_length=1.0):
    """Mesh-independent electromagnetic torque [N*m] from the recovered gap-ring phasors --
    the sum over harmonics of the closed-form AGE torque
    (:func:`airgap_element.airgap_harmonic_torque`).  Radius-independent (no Maxwell-stress
    contour, no air-gap mesh).  For complex (eddy) ``gfu`` the phasor cross-term gives the
    time-average rotating-field / drag torque."""
    T = 0.0
    for n, (a_RI, a_RO) in airgap_ring_phasors(coupling, gfu).items():
        T += airgap_harmonic_torque(n, coupling["ri"], coupling["ro"], a_RI, a_RO, axial_length)
    return T
