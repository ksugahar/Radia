# -*- coding: utf-8 -*-
r"""FEM-BEM Schur coupling -- exact open boundary for interior FEM problems.

The EXTERIOR Laplace DtN (Dirichlet-to-Neumann map) assembled from BEM
operators replaces the unbounded exterior domain with a boundary operator on
the coupling surface, giving an EXACT open boundary without mesh truncation,
PML, or Kelvin transformation.

Coupled system (FEM interior + BEM exterior DtN on outer boundary Γ):

    (K_FEM + P^T Λ_ext P) u_free = f_rhs

where:
  K_FEM  = FEM Laplacian stiffness (H1 interior + boundary rows/cols)
  Λ_ext  = V⁻¹ (-½ M_bnd + K_bnd)  [dense ndof_bnd × ndof_bnd]
  P      = SurfaceL2 DOFs ← H1 volume DOFs  (coupling / projection matrix)
  u_free = free DOF vector after Dirichlet condensation

The coupling matrix P is assembled as the L2-mass-weighted restriction of H1
basis traces to the SurfaceL2 space:

    M_bnd P = R,  R[i,j] = ∫_Γ φ_j^{H1} ψ_i^{SL2} dS

where φ_j^{H1} is the j-th H1 basis function, ψ_i^{SL2} the i-th SurfaceL2
test function, and M_bnd the SurfaceL2 mass matrix.

See :mod:`bem_integral` for the BEM DtN assembly.
See :mod:`airgap_machine` for the AGE gap element (a different Schur coupling).

Reference: Steinbach (2007) Numerical Approx. Methods for Elliptic BVPs.
"""
from __future__ import annotations

import math
import numpy as np
import ngsolve as ng
from ngsolve import (
    GridFunction, LinearForm, BilinearForm,
    ds, dx, BND, x as ng_x, y as ng_y, z as ng_z, sqrt as ng_sqrt,
)
from netgen.occ import Sphere, Pnt, OCCGeometry

from .bem_integral import laplace_exterior_dtn


def _dense_ng(mat, ndof):
    """Densify an NGSolve operator into a numpy array (ndof × ndof)."""
    e = mat.CreateColVector()
    y = mat.CreateColVector()
    D = np.zeros((ndof, ndof), dtype=complex)
    for j in range(ndof):
        e[:] = 0.0; e[j] = 1.0
        y.data = mat * e
        D[:, j] = y.FV().NumPy()[:]
    return D


def _coupling_matrix(fes_vol, fes_bnd, intorder, bnd_outer):
    """Build the coupling matrix P (ndof_bnd × ndof_vol): M_bnd P = R.

    R[i, j] = ∫_Γ φ_j^{H1}|_Γ  ψ_i^{SL2} dS  (H1 trace projected to SurfaceL2).

    Assembled as a rectangular NGSolve BilinearForm, then solved with the
    SurfaceL2 mass matrix.

    Returns ``(P, M_bnd)``.  M_bnd (the dense SurfaceL2 mass) is returned because
    the WEAK DtN boundary term needs it: ∮_Γ(∂u/∂n)v dS = P^T (M_bnd Λ_ext) P.

    NB the mixed form uses the restricted SurfaceL2 as the TRIAL space (volume H1
    as TEST).  Using the ``definedon``-restricted SurfaceL2 as the TEST space of a
    mixed form on a VOLUME mesh corrupts the heap on teardown -> nondeterministic
    NGSolve crash (exit 127/139); as the TRIAL space it is stable (isolated by the
    bem-crash-isolate experiment: swapped/full-boundary stable, definedon-as-test
    crashes regardless of dual_mapping or integration region)."""
    ndof_bnd = fes_bnd.ndof
    ndof_vol = fes_vol.ndof

    # Mixed bilinear form: SurfaceL2 (restricted) as TRIAL, volume H1 as TEST.
    # mat is (ndof_vol [test] × ndof_bnd [trial]); mat * e_i (e_i a bnd unit
    # vector) = column i = R[i, :] (the i-th SurfaceL2 basis fn against all H1).
    u_bnd = fes_bnd.TrialFunction()
    v_vol = fes_vol.TestFunction()
    r_bf = BilinearForm(trialspace=fes_bnd, testspace=fes_vol)
    r_bf += u_bnd.Trace() * v_vol * ds(bnd_outer, bonus_intorder=intorder - 4)
    r_bf.Assemble()

    e_bnd = r_bf.mat.CreateRowVector()    # size ndof_bnd (trial)
    y_vol = r_bf.mat.CreateColVector()    # size ndof_vol (test)
    R = np.zeros((ndof_bnd, ndof_vol))
    for i in range(ndof_bnd):
        e_bnd[:] = 0.0; e_bnd[i] = 1.0
        y_vol.data = r_bf.mat * e_bnd
        R[i, :] = y_vol.FV().NumPy()[:]

    # Mass matrix on SurfaceL2 boundary (real-valued; drop the complex dtype that
    # _dense_ng allocates so the P solve stays real).
    u_bnd, v_bnd2 = fes_bnd.TnT()
    mass_bf = BilinearForm(u_bnd.Trace() * v_bnd2.Trace() * ds(bonus_intorder=intorder - 4))
    mass_bf.Assemble()
    M_d = np.real(_dense_ng(mass_bf.mat, ndof_bnd))

    # When fes_bnd is a SurfaceL2 restricted to one boundary of a VOLUME mesh, the
    # DOFs on the other boundaries are non-free (zero mass rows) -> M_d singular.
    # Solve on the FREE sub-block; non-free P rows stay 0 (those basis traces
    # vanish on Γ anyway).  Standalone closed surface: all free -> plain solve.
    # Returns (P, M_d): the surface mass M_d is needed to assemble the WEAK DtN
    # boundary term  ∮(∂u/∂n)v dS = P^T (M_bnd Λ) P  (Λ = laplace_exterior_dtn).
    free = np.array([fes_bnd.FreeDofs()[i] for i in range(ndof_bnd)], dtype=bool)
    if free.all():
        return np.linalg.solve(M_d, R), M_d
    idx = np.where(free)[0]
    P = np.zeros((ndof_bnd, ndof_vol), dtype=R.dtype)
    P[idx, :] = np.linalg.solve(M_d[np.ix_(idx, idx)], R[idx, :])
    return P, M_d


def laplace_fem_bem_schur(
    mesh,
    bnd_outer,
    bnd_dirichlet="",
    dirichlet_cf=None,
    source_cf=None,
    order=2,
    intorder=8,
    bem_intorder=None,
):
    """Solve -Δu = source inside ``mesh`` with exact BEM DtN on ``bnd_outer``.

    Assembles the fully coupled (dense) system:

        (K_FEM − P^T M_bnd Λ_ext P) u_free = f_rhs

    where Λ_ext is the exterior Laplace DtN from :func:`bem_integral.laplace_exterior_dtn`,
    M_bnd is the SurfaceL2 boundary mass, and P maps H1 volume DOFs to SurfaceL2
    boundary DOFs via L2 projection.  The boundary term is the WEAK DtN
    −∮_Γ(∂u/∂n)v dS = −⟨Λ_ext u, v⟩ = −P^T M_bnd Λ_ext P (the surface mass M_bnd
    turns the DtN DOF map Λ_ext into the weak form -- omitting it is a magnitude
    bug).  The MINUS sign with Λ_ext negative-definite (eigenvalue −(n+1)/R) gives
    a positive-definite term, matching the analytic Robin BC +(2/R)∮uv of
    :func:`sphere_shell_analytic_dtn` for the dipole.

    Validated end-to-end on the spherical-shell dipole at order=2: L2 rel_err
    5.5e-2 (maxh=0.5) / 5.4e-2 (maxh=0.4), matching the analytic-DtN Robin path
    (:func:`sphere_shell_analytic_dtn`) to 3 sig figs -- the FEM interior error
    dominates at this resolution.  Deterministic over repeated runs.

    .. note::
        The coupling mixed form (:func:`_coupling_matrix`) uses the ``definedon``-
        restricted SurfaceL2 as the TRIAL space (volume H1 as TEST).  Using it as
        the TEST space on a VOLUME mesh corrupts the heap on teardown -> a
        nondeterministic NGSolve crash (exit 127/139, any dual_mapping).  This was
        isolated (swapped/full-boundary stable; definedon-as-test crashes
        regardless of integration region) and fixed by the swap; do NOT revert it.

    For large problems use a sparse iterative solver with the DtN as a shell operator;
    this dense reference implementation is intended for verification up to ~1000 DOFs.

    Parameters
    ----------
    mesh          : ng.Mesh  Volume mesh; ``bnd_outer`` must be a closed surface.
    bnd_outer     : str      Outer coupling boundary name.
    bnd_dirichlet : str      Dirichlet boundary name(s) (may be empty).
    dirichlet_cf  : CF       Dirichlet datum (set on BND before solve).
    source_cf     : CF       RHS body source (optional).
    order         : int      H1 polynomial order (default 2).
    intorder      : int      FEM quadrature order.
    bem_intorder  : int      BEM quadrature order (default intorder+2).

    Returns
    -------
    gfu : GridFunction   Solution on H1(mesh, order)."""
    if bem_intorder is None:
        bem_intorder = intorder + 2

    fes_vol = ng.H1(mesh, order=order, dirichlet=bnd_dirichlet)
    fes_bnd = ng.SurfaceL2(mesh, order=max(order - 1, 0), dual_mapping=True,
                           definedon=mesh.Boundaries(bnd_outer))

    u, v = fes_vol.TnT()
    a = BilinearForm(ng.grad(u) * ng.grad(v) * dx(bonus_intorder=intorder))
    a.Assemble()
    K = _dense_ng(a.mat, fes_vol.ndof)

    f = LinearForm(fes_vol)
    if source_cf is not None:
        f += source_cf * v * dx(bonus_intorder=intorder)
    f.Assemble()

    gfu = GridFunction(fes_vol)
    if dirichlet_cf is not None and bnd_dirichlet:
        gfu.Set(dirichlet_cf, BND)

    # BEM DtN and coupling matrix (P maps vol→bnd DOFs; M_bnd = SurfaceL2 mass)
    dtn = laplace_exterior_dtn(fes_bnd, intorder=bem_intorder)
    P, M_bnd = _coupling_matrix(fes_vol, fes_bnd, intorder, bnd_outer)

    # Coupled stiffness: weak boundary term −∮(∂u/∂n)v dS = −P^T (M_bnd Λ_ext) P
    # (the surface mass M_bnd turns the DtN DOF map Λ_ext into the weak form;
    # for the dipole −P^T M_bnd Λ P = +(2/R)∮uv, the analytic Robin term).
    K_coupled = K - P.T @ (M_bnd @ dtn) @ P
    # the physical problem is real; the BEM densification leaves ~machine-zero
    # imaginary parts in Λ_ext -- drop them so the H1 solution stays real.
    K_coupled = np.real(K_coupled)

    # RHS: lift Dirichlet, solve on free DOFs
    rhs = np.real(f.vec.FV().NumPy()[:]) - K_coupled @ np.real(gfu.vec.FV().NumPy()[:])
    free = np.array([fes_vol.FreeDofs()[i] for i in range(fes_vol.ndof)], dtype=bool)

    sol = np.linalg.solve(K_coupled[np.ix_(free, free)], rhs[free])
    full = np.real(gfu.vec.FV().NumPy()[:])
    full[free] = sol
    gfu.vec.FV().NumPy()[:] = full
    return gfu


def sphere_shell_analytic_dtn(R_inner=0.5, R_outer=1.0, maxh_vol=0.3,
                              order=2, intorder=8):
    """Verification: Laplace on spherical shell with ANALYTIC DtN Robin BC on outer sphere.

    Inner BC: u = z/R_inner  (dipole datum, Y₁¹ mode, = cos θ at r=R_inner)
    Outer BC: ∂u/∂n = (−2/R_outer) u  (analytic DtN eigenvalue for dipole)

    Exact solution:  u(r) = R_inner² z / r³   (exterior dipole, decays ~ 1/r²)
    which satisfies ∇²u = 0 and:
      - u(R_inner,θ) = R_inner² cos θ / R_inner² = cos θ = z/R_inner ✓
      - ∂u/∂r|_{R_outer} = -2 R_inner² z / R_outer⁴ = (-2/R_outer) × (R_inner²/R_outer²) cos θ ✓
        (Λ_ext eigenvalue = -2/R_outer for dipole)

    Exact solution: u(r,θ) = R_inner² cos θ / r²  = R_inner² z / r³

    This test validates the FEM+DtN Robin BC framework.  The BEM exterior DtN is
    separately validated by ``sphere_exterior_dtn_eigenvalue`` (eigenvalue −2/R for
    the dipole mode).  Together they confirm:
        FEM + BEM-DtN = FEM + analytic-DtN (Robin) = exact open BC solution.

    Returns dict with L2 relative error and DOF count."""
    from netgen.occ import Sphere as NGSphere, Pnt as NGPnt, OCCGeometry as NGOCC
    outer_s = NGSphere(NGPnt(0, 0, 0), R_outer); outer_s.bc("outer")
    inner_s = NGSphere(NGPnt(0, 0, 0), R_inner); inner_s.bc("inner")
    mesh = ng.Mesh(NGOCC(outer_s - inner_s).GenerateMesh(maxh=maxh_vol)).Curve(3)

    fes = ng.H1(mesh, order=order, dirichlet="inner")
    u, v = fes.TnT()
    a = BilinearForm(fes)
    a += ng.grad(u) * ng.grad(v) * dx(bonus_intorder=intorder)
    a += (2.0 / R_outer) * u * v * ds("outer", bonus_intorder=intorder)
    a.Assemble()

    gfu = ng.GridFunction(fes)
    gfu.Set(ng_z / R_inner, ng.BND)
    rhs = gfu.vec.CreateVector()
    rhs.data = -(a.mat * gfu.vec)
    gfu.vec.data += a.mat.Inverse(freedofs=fes.FreeDofs()) * rhs

    r2 = ng_x**2 + ng_y**2 + ng_z**2
    r3 = r2 * ng_sqrt(r2)
    u_exact = R_inner**2 * ng_z / r3

    err_cf = gfu - u_exact
    l2_err = math.sqrt(abs(float(ng.Integrate(err_cf * err_cf * dx, mesh))))
    l2_ref = math.sqrt(abs(float(ng.Integrate(u_exact * u_exact * dx, mesh))))
    return {
        "ndof": fes.ndof,
        "R_inner": R_inner, "R_outer": R_outer,
        "l2_err": l2_err,
        "l2_ref": l2_ref,
        "rel_err": l2_err / l2_ref if l2_ref > 0 else float("nan"),
    }


# ---------------------------------------------------------------------------
# Kelvin-closure effective DtN eigenvalue -- the VOLUME-FEM companion to the
# BEM exterior_dtn_spectrum (bem_integral.py)
# ---------------------------------------------------------------------------
#
# The Kelvin transformation realises the exterior open boundary by a VOLUME FEM
# solve on the inverted ball.  This function MEASURES the effective exterior DtN
# eigenvalue that the Kelvin closure presents on the truncation sphere Γ for a
# pure spherical-harmonic mode n, and compares it to the analytic −(n+1)/R.
#
# Mechanism (why Kelvin is coarse-mesh accurate, MEASURED not argued): the 3D
# Kelvin inversion r ↦ R²/ρ with the conformal weight (R/ρ) maps the slow-decay
# exterior mode r^−(n+1) Y_n to the SOLID HARMONIC ρ^n Y_n on the ball -- a
# degree-n polynomial.  Order-p Lagrange FEM represents a degree-n polynomial
# EXACTLY iff p ≥ n, so:
#   - dipole  (n=1) at order ≥ 1: exact up to GEOMETRY (curved-sphere) error,
#     which converges fast under refinement -- accurate on a coarse mesh;
#   - quadrupole (n=2): order=1 cannot represent the quadratic image (large
#     error); order ≥ 2 is exact up to geometry error;
#   - mode n needs order ≥ n.
# A compact source's far field is dipole-dominated, and the dipole inverts to a
# LINEAR field that even order-1 coarse FEM resolves -- which IS the coarse-mesh
# accuracy Kameari demonstrated, here read off the realised DtN eigenvalue.
#
# This is the COMPLEMENT of the BEM spectrum: BEM Λ_h is accurate for all smooth
# low surface-harmonics (error grows smoothly with n); the Kelvin VOLUME closure
# is accurate exactly up to its polynomial order (sharp order threshold n ≤ p).
#
# Knowledge: dtn_coarse_mesh.  BEM side: bem_integral.exterior_dtn_spectrum.


def _solid_harmonic(n):
    """Zonal solid harmonic h_n = ρ^n P_n(cosθ) as a Cartesian-polynomial CF, for
    ARBITRARY degree n ≥ 0, built by the Legendre recurrence

        (k+1) h_{k+1} = (2k+1) z h_k − k ρ² h_{k−1},   h_0 = 1, h_1 = z, ρ²=x²+y²+z².

    Its trace on |x|=R is the degree-n boundary datum (a spherical harmonic), and it
    is itself the EXACT Kelvin image u* of the exterior mode r^{−(n+1)}Y_n (a degree-n
    harmonic polynomial, Δ h_n = 0).  The (Legendre) normalization is immaterial: every
    consumer here reports a SCALE-INVARIANT quantity (the effective-DtN energy quotient,
    or a relative error in which the datum and the exact solution share the same h_n).
    Lifting the former n ≤ 3 cap lets the n=1..9 spectrum be computed from the library."""
    if n < 0:
        raise ValueError("solid harmonic degree must be n >= 0, got %r" % n)
    r2 = ng_x * ng_x + ng_y * ng_y + ng_z * ng_z
    h = [ng.CoefficientFunction(1.0), ng_z]
    for k in range(1, n):
        h.append(((2 * k + 1) * ng_z * h[k] - k * r2 * h[k - 1]) / (k + 1))
    return h[n]


def _solid_harmonic_2d(n):
    """2D solid harmonic h_n = r^n cos(nθ) as a Cartesian-polynomial CF (Δ h_n = 0),
    for ARBITRARY n ≥ 0, by the Chebyshev-type recurrence

        h_{k+1} = 2 x h_k − (x²+y²) h_{k−1},   h_0 = 1, h_1 = x.

    The 2D analogue of :func:`_solid_harmonic`: the exact Kelvin image of the exterior
    mode r^{−n} e^{inθ} under the 2D (circle) inversion.  Used by
    :func:`kelvin_dtn_eigenvalue` with ``dim=2``.  (n=0..3 reproduce 1, x, x²−y²,
    x³−3xy² exactly.)"""
    if n < 0:
        raise ValueError("2D solid harmonic degree must be n >= 0, got %r" % n)
    r2 = ng_x * ng_x + ng_y * ng_y
    h = [ng.CoefficientFunction(1.0), ng_x]
    for k in range(1, n):
        h.append(2 * ng_x * h[k] - r2 * h[k - 1])
    return h[n]


def kelvin_dtn_eigenvalue(R=1.0, degree=1, maxh=0.4, order=1, intorder=10, dim=3, curve=None):
    """Measure the effective exterior DtN eigenvalue the KELVIN closure realises on
    the truncation surface for harmonic mode ``degree`` (n), vs the analytic value.

    The Kelvin transformation IS a spherical inversion (r ↦ R²/r), so its truncation
    surface is a sphere (``dim=3``) or circle (``dim=2``).  The exterior mode is
    mapped to the solid harmonic ρ^n Y_n -- a degree-n harmonic polynomial -- and we
    solve Laplace on the inverted Kelvin ball/disk with that datum, then form the
    variationally consistent effective DtN eigenvalue from the energy quotient:

        λ_eff = offset − ∫_Ω |∇u*|²  /  ∮_Γ u*² ,     offset = −(dim−2)/R
        λ_exact = −(n + dim − 2)/R   ( 3D: −(n+1)/R ;  2D: −n/R )

    The ``offset`` is the conformal weight (R/r)^{d−2} of the Kelvin map (present in
    3D, ABSENT in 2D).  At FEM ``order ≥ n`` the degree-n polynomial image is captured
    in the REFERENCE space, so the POLYNOMIAL-approximation error vanishes (the dipole
    inverts to a LINEAR field -> order-1 coarse accurate); ``order < n`` is
    polynomial-approximation limited.  NB on the curved 3D sphere the REALISED accuracy
    then FLOORS at the isoparametric-geometry + conformal-weight error (~5-6 digits at
    order = n, e.g. n=2 -> 2.5e-3, n=3 -> 8.4e-4, both -> ~1e-5 as order rises) -- this
    is Kameari's "5-6 digits", NOT machine-exact; in 2D (no conformal weight) the floor
    is deeper (~1e-7…1e-9).  ``dim=2`` is the cross-section case for static apparatus /
    rotating machines.

    ``curve`` sets the isoparametric (geometry) order of the curved truncation sphere;
    ``None`` keeps the historical default ``min(order+1, 3)`` (3D) / ``min(order+1, 4)``
    (2D).  Because the residual floor IS this curved-geometry error, raising ``curve``
    lowers it ~30×/order (e.g. 3D order=9: Curve 3 ~1e-5 → Curve 5 ~1e-7); ``degree`` is
    no longer capped (the solid-harmonic image is built by recurrence for any n).

    Returns dict ``{ndof, degree, order, maxh, dim, curve, lam, lam_exact, rel_err}``."""
    if dim not in (2, 3):
        raise ValueError("dim must be 2 or 3, got %r" % dim)
    if dim == 3:
        cv = curve if curve is not None else min(order + 1, 3)
        mesh = ng.Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), R)).GenerateMesh(maxh=maxh)).Curve(cv)
        datum = _solid_harmonic(degree)
    else:
        from netgen.geom2d import SplineGeometry
        cv = curve if curve is not None else min(order + 1, 4)
        geo = SplineGeometry(); geo.AddCircle((0, 0), R, bc="kelvin_circle")
        mesh = ng.Mesh(geo.GenerateMesh(maxh=maxh)); mesh.Curve(cv)
        datum = _solid_harmonic_2d(degree)

    fes = ng.H1(mesh, order=order, dirichlet=".*")
    u, v = fes.TnT()
    a = BilinearForm(ng.grad(u) * ng.grad(v) * dx)
    a.Assemble()

    gfu = GridFunction(fes)
    gfu.Set(datum, BND)                                   # u* = ρ^n Y_n on Γ
    rhs = gfu.vec.CreateVector()
    rhs.data = -(a.mat * gfu.vec)
    gfu.vec.data += a.mat.Inverse(freedofs=fes.FreeDofs()) * rhs

    # λ_eff = offset − ∫|∇u*|² / ∮u*²   (weak, variationally consistent)
    offset = -(dim - 2) / R                               # 3D: −1/R ; 2D: 0
    energy = float(ng.Integrate(ng.grad(gfu) * ng.grad(gfu) * dx(bonus_intorder=intorder), mesh))
    bmass = float(ng.Integrate(gfu * gfu * ds(bonus_intorder=intorder), mesh))
    lam = offset - energy / bmass if bmass > 0 else float("nan")
    lam_exact = -(degree + dim - 2) / R
    return {
        "ndof": fes.ndof, "degree": degree, "order": order, "maxh": maxh, "dim": dim,
        "curve": cv, "lam": lam, "lam_exact": lam_exact,
        "rel_err": abs(lam - lam_exact) / abs(lam_exact),
    }


# ---------------------------------------------------------------------------
# Isolated open-BC error: Kelvin (and BEM) vs the EXACT open boundary, with the
# interior FEM discretisation error CANCELLED
# ---------------------------------------------------------------------------
#
# The total error of an open-BC field solve is  interior-FEM error + open-BC
# (operator) error.  A field-refinement study (Kameari) sees the SUM.  To
# evaluate the KELVIN method's OWN error one must cancel the interior-FEM part.
#
# On a SHARED interior mesh, swap ONLY the boundary operator on Γ:
#   exact open BC : Robin ∂u/∂n = λ_exact u,  λ_exact = −(n+1)/R   → u_exactDtN
#   Kelvin closure: Robin ∂u/∂n = λ_Kelvin u (λ_Kelvin = kelvin_dtn_eigenvalue,
#                   the Kelvin closure's measured effective DtN)    → u_Kelvin
#   BEM DtN       : the validated laplace_fem_bem_schur             → u_BEM
# For a single-mode (multipole) shell problem the Kelvin closure acts on Γ as
# multiplication by its effective DtN eigenvalue, so the Robin-with-λ_Kelvin solve
# IS the Kelvin field (for that mode) -- no separate Kelvin mesh/periodic-BC.
#
# u_exactDtN uses the EXACT open BC, so its error vs the analytic field is PURE
# interior FEM error.  Then  ||u_Kelvin − u_exactDtN||  is the isolated Kelvin
# open-BC error (interior error cancels: same mesh, only Γ-operator differs).


def _shell_robin_solve(mesh, R_inner, degree, lam, order, intorder):
    """Solve −Δu=0 on the shell with inner Dirichlet = the mode-n solid harmonic
    and a Robin BC ∂u/∂n = lam·u on "outer" (the weak term −∮(∂u/∂n)v = −lam∮uv)."""
    inner_datum = _solid_harmonic(degree) / R_inner ** degree
    fes = ng.H1(mesh, order=order, dirichlet="inner")
    u, v = fes.TnT()
    a = BilinearForm(fes)
    a += ng.grad(u) * ng.grad(v) * dx(bonus_intorder=intorder)
    a += (-lam) * u * v * ds("outer", bonus_intorder=intorder)
    a.Assemble()
    gfu = GridFunction(fes)
    gfu.Set(inner_datum, BND)
    rhs = gfu.vec.CreateVector()
    rhs.data = -(a.mat * gfu.vec)
    gfu.vec.data += a.mat.Inverse(freedofs=fes.FreeDofs()) * rhs
    return gfu


def kelvin_vs_exact_open_bc_error(R_inner=0.5, R_outer=1.0, maxh=0.4, order=2,
                                  degree=1, kelvin_maxh=0.6, kelvin_order=1,
                                  intorder=8, include_bem=True):
    """Evaluate the KELVIN method's open-boundary error, ISOLATED from the
    interior FEM discretisation error, on the spherical-shell multipole problem.

    Solves the SAME shell interior mesh three ways (exact-DtN Robin, Kelvin-DtN
    Robin, BEM-DtN Schur) and compares.  Because the exact-DtN solve carries ZERO
    open-BC error, ``||u_method − u_exactDtN||`` is the method's pure open-BC error
    (the shared interior FEM error cancels).

    Parameters
    ----------
    R_inner, R_outer : float  shell radii (Γ = outer sphere, the truncation surface).
    maxh, order      : interior shell discretisation (SHARED by all methods).
    degree           : spherical-harmonic degree n of the multipole (1=dipole, …).
    kelvin_maxh, kelvin_order : the KELVIN closure's own (independent) discretisation
                       -- coarse by default, to show coarse-mesh accuracy.
    include_bem      : also run laplace_fem_bem_schur (slower dense BEM).

    Returns dict:
      ndof, degree, lam_exact, lam_kelvin, kelvin_dtn_rel_err,
      interior_fem_error  -- ||u_exactDtN − u_analytic|| (shared FEM error; open BC exact),
      kelvin_total_error  -- ||u_Kelvin − u_analytic||,
      kelvin_openbc_error -- ||u_Kelvin − u_exactDtN||   (ISOLATED Kelvin open-BC error),
      bem_total_error / bem_openbc_error  (if include_bem).

    The headline: ``kelvin_openbc_error`` is far below ``interior_fem_error`` even
    on a coarse Kelvin mesh -- the Kelvin closure is NOT the accuracy bottleneck,
    which is Kameari's coarse-mesh accuracy stated as a separated error budget."""
    outer = Sphere(Pnt(0, 0, 0), R_outer); outer.bc("outer")
    inner = Sphere(Pnt(0, 0, 0), R_inner); inner.bc("inner")
    mesh = ng.Mesh(OCCGeometry(outer - inner).GenerateMesh(maxh=maxh)).Curve(3)

    r2 = ng_x ** 2 + ng_y ** 2 + ng_z ** 2
    r = ng_sqrt(r2)
    u_exact = R_inner ** (degree + 1) * _solid_harmonic(degree) / r ** (2 * degree + 1)

    def l2(cf):
        return math.sqrt(abs(float(ng.Integrate(cf * cf * dx, mesh))))

    lam_exact = -(degree + 1) / R_outer
    kel = kelvin_dtn_eigenvalue(R=R_outer, degree=degree, maxh=kelvin_maxh,
                                order=kelvin_order, intorder=10)
    lam_kelvin = kel["lam"]

    u_exactDtN = _shell_robin_solve(mesh, R_inner, degree, lam_exact, order, intorder)
    u_kelvin = _shell_robin_solve(mesh, R_inner, degree, lam_kelvin, order, intorder)

    ref = l2(u_exact)
    out = {
        "ndof": u_exactDtN.space.ndof, "degree": degree,
        "R_inner": R_inner, "R_outer": R_outer, "maxh": maxh, "order": order,
        "kelvin_maxh": kelvin_maxh, "kelvin_order": kelvin_order,
        "kelvin_ndof": kel["ndof"],
        "lam_exact": lam_exact, "lam_kelvin": lam_kelvin,
        "kelvin_dtn_rel_err": kel["rel_err"],
        "interior_fem_error": l2(u_exactDtN - u_exact) / ref,
        "kelvin_total_error": l2(u_kelvin - u_exact) / ref,
        "kelvin_openbc_error": l2(u_kelvin - u_exactDtN) / ref,
    }
    if include_bem:
        inner_datum = _solid_harmonic(degree) / R_inner ** degree
        u_bem = laplace_fem_bem_schur(mesh, bnd_outer="outer", bnd_dirichlet="inner",
                                      dirichlet_cf=inner_datum, order=order, intorder=intorder)
        out["bem_total_error"] = l2(u_bem - u_exact) / ref
        out["bem_openbc_error"] = l2(u_bem - u_exactDtN) / ref
    return out


def kelvin_openbc_error_vs_exterior_mesh(kelvin_maxh_list=(0.7, 0.5, 0.35, 0.25),
                                         R_inner=0.5, R_outer=1.0, maxh=0.4, order=2,
                                         degree=1, kelvin_order=1, intorder=8):
    """Verify open-boundary accuracy as a function of the EXTERIOR-region mesh size.

    In the Kelvin transformation the unbounded exterior is the Kelvin ball, so the
    "exterior mesh size" is ``kelvin_maxh``.  This sweeps it with the INTERIOR shell
    mesh held FIXED, so the (fixed) interior FEM error cancels and the reported
    open-BC error is purely the exterior-discretisation contribution.

    This is Kameari's exterior-refinement study, ISOLATED: refine the exterior and
    watch the open-BC error converge -- while the interior FEM error stays put.

    Returns dict:
      interior_fem_error : the FIXED interior FEM error (open BC exact); independent
                           of the exterior mesh (the cancellation check).
      per_mesh           : list of one dict per exterior mesh, each with
                           {kelvin_maxh, kelvin_ndof, kelvin_dtn_rel_err,
                            kelvin_openbc_error, ratio_fem_over_openbc}
      converges          : True if kelvin_openbc_error is non-increasing as the
                           exterior mesh refines (distinct meshes only).
      always_below_fem   : True if kelvin_openbc_error < interior_fem_error at EVERY
                           exterior mesh (the exterior is never the bottleneck).

    Signature fact: the open-BC error CONVERGES with exterior refinement but stays
    far below the fixed interior FEM error at every level -- so a coarse exterior
    (Kelvin) mesh already suffices.  That is the quantitative, error-isolated form
    of Kameari's coarse-mesh accuracy, read off the EXTERIOR mesh size."""
    per_mesh = []
    interior_fem_error = None
    for kmaxh in kelvin_maxh_list:
        r = kelvin_vs_exact_open_bc_error(
            R_inner=R_inner, R_outer=R_outer, maxh=maxh, order=order, degree=degree,
            kelvin_maxh=kmaxh, kelvin_order=kelvin_order, intorder=intorder,
            include_bem=False)
        interior_fem_error = r["interior_fem_error"]   # fixed (interior mesh fixed)
        ob = r["kelvin_openbc_error"]
        per_mesh.append({
            "kelvin_maxh": kmaxh,
            "kelvin_ndof": r["kelvin_ndof"],
            "kelvin_dtn_rel_err": r["kelvin_dtn_rel_err"],
            "kelvin_openbc_error": ob,
            "ratio_fem_over_openbc": (interior_fem_error / ob) if ob > 0 else float("inf"),
        })

    # convergence over DISTINCT exterior meshes (the Kelvin ball floors for coarse maxh)
    distinct = []
    seen = set()
    for m in per_mesh:
        if m["kelvin_ndof"] not in seen:
            seen.add(m["kelvin_ndof"]); distinct.append(m)
    converges = all(distinct[i + 1]["kelvin_openbc_error"] <= distinct[i]["kelvin_openbc_error"]
                    for i in range(len(distinct) - 1))
    always_below_fem = all(m["kelvin_openbc_error"] < interior_fem_error for m in per_mesh)

    return {"interior_fem_error": interior_fem_error, "per_mesh": per_mesh,
            "converges": converges, "always_below_fem": always_below_fem,
            "degree": degree, "maxh": maxh, "order": order}


def kelvin_twosphere_shell_dipole(R_inner=0.5, R_outer=1.0, offset=3.0,
                                  maxh=0.25, order=2, intorder=6):
    """The lab's REAL Kelvin transformation (two offset spheres + periodic BC +
    material modulation), solved for the shell dipole, vs the exact u = R_in²z/r³.

    Confirms that the genuine Kelvin coupling -- NOT the single-ball effective-DtN
    equivalent used elsewhere -- reproduces the analytic exterior dipole.  Setup
    (Nagamine convention, see knowledge: kelvin_transformation):
      * interior shell R_in<r<R_out, inner Dirichlet u = z/R_in (the dipole datum);
      * a SEPARATE Kelvin sphere of the SAME radius R_out, OFFSET in space; the
        shell's outer face and the Kelvin sphere are PERIODICALLY identified
        (translation) -- this is the true Kelvin inversion, not a concentric shell;
      * material modulation μ' = (R_out/r')² in the Kelvin ball (r' = distance from
        its centre = the image of infinity), GND (u=0) at that centre;
      * solve ∫ μ ∇u·∇v over shell ∪ Kelvin-ball with the Periodic() H1 space.

    Returns ``{ndof, maxh, order, rel_err}`` -- the shell-region L2 error vs the
    exact dipole (total error; the Kelvin open-BC part is tiny, so this is mostly
    the interior FEM error).  ~1.9 % at maxh=0.18 order=2."""
    import netgen.occ as occ
    from netgen.occ import Sphere as _Sph, Pnt as _Pnt, Vec as _Vec, IdentificationType
    from ngsolve import Periodic, grad as _grad, x as _x, y as _y, z as _z, sqrt as _sqrt

    outer = _Sph(_Pnt(0, 0, 0), R_outer)
    inner = _Sph(_Pnt(0, 0, 0), R_inner)
    for f in outer.faces:
        f.name = "kelvin_int"
    for f in inner.faces:
        f.name = "inner"
    shell = (outer - inner); shell.mat("shell")
    kball = _Sph(_Pnt(offset, 0, 0), R_outer)
    for f in kball.faces:
        f.name = "kelvin_ext"
    kball.mat("kelvin")
    gnd = occ.Vertex(_Pnt(offset, 0, 0)); gnd.name = "GND"

    fi = [f for f in shell.faces if f.name == "kelvin_int"][0]
    fe = [f for f in kball.faces if f.name == "kelvin_ext"][0]
    fi.Identify(fe, "kelvin", IdentificationType.PERIODIC,
                occ.gp_Trsf.Translation(_Vec(offset, 0, 0)))

    shape = occ.Glue([shell, kball, gnd])
    mesh = ng.Mesh(OCCGeometry(shape).GenerateMesh(maxh=maxh)).Curve(3)

    fes = Periodic(ng.H1(mesh, order=order, dirichlet="inner|GND"))
    rp2 = (_x - offset) ** 2 + _y ** 2 + _z ** 2 + 1e-20
    mu = mesh.MaterialCF({"shell": 1.0, "kelvin": R_outer ** 2 / rp2}, default=1.0)

    u, v = fes.TnT()
    a = BilinearForm(mu * _grad(u) * _grad(v) * dx(bonus_intorder=intorder))
    a.Assemble()
    gfu = GridFunction(fes)
    gfu.Set(_z / R_inner, BND, definedon=mesh.Boundaries("inner"))
    rhs = gfu.vec.CreateVector()
    rhs.data = -(a.mat * gfu.vec)
    gfu.vec.data += a.mat.Inverse(freedofs=fes.FreeDofs(), inverse="sparsecholesky") * rhs

    r2 = _x * _x + _y * _y + _z * _z
    u_exact = R_inner ** 2 * _z / (r2 * _sqrt(r2))
    err = math.sqrt(abs(float(ng.Integrate((gfu - u_exact) ** 2 * dx("shell"), mesh))))
    ref = math.sqrt(abs(float(ng.Integrate(u_exact ** 2 * dx("shell"), mesh))))
    return {"ndof": fes.ndof, "maxh": maxh, "order": order,
            "rel_err": err / ref if ref > 0 else float("nan")}
