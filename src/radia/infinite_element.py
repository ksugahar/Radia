# -*- coding: utf-8 -*-
"""radia.infinite_element -- NGSolve-native static / low-frequency INFINITE-ELEMENT open boundary.

The (static) infinite element (IE) closes an exterior Laplace problem on a truncation surface by
expanding the decaying exterior field as (surface FE) x (radial decay basis) and adding the radial
levels as EXTRA surface DOFs that NGSolve assembles + solves monolithically with the interior FE
system.  On a SPHERE this is identical to the Kelvin transformation (same exterior polynomial space;
see `docs/open_boundary/INFINITE_ELEMENT_SOTA.md` and
`validation_test/open_boundary/test_infinite_element.py`).

DESIGN (NGSolve-native -- "complement NGSolve, do not reimplement"):
  * the radial decay operators are tiny (P x P) and live in numpy here (the well-conditioned vertex +
    integrated-Legendre basis -- NEVER naive monomials, per the act7_28 conditioning lesson);
  * the P-1 radial "bubble" levels are boundary-only NGSolve H1 spaces, compounded with the interior
    H1; the exterior energy ``sum_kl [ R1_kl int_Gamma u_k u_l ds + R0_kl int_Gamma grad_S u_k . grad_S
    u_l ds ]`` is added as boundary integrals, so NGSolve assembles the augmented SPARSE system and
    solves it monolithically (no scipy, no dense Schur condensation);
  * the exterior field is recovered everywhere as the radial expansion ``phi(r,s) = sum_k U_k(s)
    N_k(a/r)`` from the solved level GridFunctions (:func:`exterior_field`).

The exterior energy is the TENSOR PRODUCT of a surface part and a radial part: for a single mode n the
block reduces to ``R1 + n(n+1) R0`` (act7_25).  Note the ENERGY operator ``S`` (the condensed surface
stiffness added to the FE system) has per-mode eigenvalue ``eig(S, Mtil) = (n+1)*a`` (the exterior
energy / radial Steklov), while the DtN eigenvalue (du/dr over u at r=a) is ``-(n+1)/a``; these
coincide only at ``a=1`` (the default).  The FE coupling uses the energy operator (correct for all a).

Public API
----------
``radial_operators(P, a)``                  -> (R1, R0, g)  the P x P radial operators + trace vector.
``dtn_surface_matrix(MS, KS, P, a)``         -> S            condensed DtN matrix (numpy, for analysis).
``ie_compound_space(mesh, P, order)``        -> X            compound [interior H1, P-1 surface H1].
``add_exterior_ie(a_bf, X, P, a)``           -> None         add the IE boundary terms to a BilinearForm.
``exterior_field(gf, P, a, points)``         -> ndarray      evaluate the exterior radial expansion.
"""
import numpy as np

__all__ = ["radial_operators", "dtn_surface_matrix", "ie_compound_space",
           "add_exterior_ie", "exterior_field"]


# ===========================================================================
# radial decay operators (numpy) -- orthogonal nodal basis (vertex + integrated-Legendre bubbles)
# ===========================================================================
def _legval(j, xi):
    c = np.zeros(j + 1); c[j] = 1.0
    return np.polynomial.legendre.legval(xi, c)


def _radial_eval(P, t):
    """nodal basis on t = a/r in (0,1]: N_1=t (vertex/trace), N_k=int-Legendre bubble (0 at t=0,1)."""
    t = np.asarray(t, float)
    N = np.zeros((P, t.size)); Np = np.zeros((P, t.size))
    N[0] = t; Np[0] = np.ones_like(t)
    xi = 2.0 * t - 1.0
    for k in range(2, P + 1):
        N[k - 1] = (_legval(k, xi) - _legval(k - 2, xi)) / (2.0 * k - 1.0)
        Np[k - 1] = _legval(k - 1, xi) * 2.0
    return N, Np


def radial_operators(P, a=1.0, nq=160):
    """Radial decay operators (sphere radius ``a``): ``R1`` (P x P) radial stiffness
    ``int rho_k' rho_l' r^2 dr``, ``R0`` (P x P) radial mass ``int rho_k rho_l dr``, and the trace
    vector ``g`` (= e_1; only the vertex function is nonzero at the surface).
    """
    x, w = np.polynomial.legendre.leggauss(nq)
    t = 0.5 * (x + 1.0); w = 0.5 * w
    N, Np = _radial_eval(P, t)
    R1 = a * (Np * w) @ Np.T
    R0 = a * (N / t ** 2 * w) @ N.T
    g = np.zeros(P); g[0] = 1.0
    return R1, R0, g


def dtn_surface_matrix(MS, KS, P, a=1.0, nq=160):
    """Condensed DtN surface stiffness ``S`` (N x N) from the UNIT-SPHERE surface mass ``MS`` and
    Laplace-Beltrami ``KS`` (N x N, symmetric).  Builds the P-level tensor blocks ``R1_kl*MS +
    R0_kl*KS`` and statically condenses the radial bubble levels onto the trace (numpy).  Useful for
    the discrete Steklov spectrum: ``eig(S, MS)`` gives the per-mode exterior ENERGY ``(n+1)*a`` (the
    condensed radial stiffness), which equals the DtN magnitude ``(n+1)/a`` only at ``a=1``.  The
    production solve keeps the levels explicit (:func:`add_exterior_ie`) and does NOT condense.
    """
    MS = np.ascontiguousarray(MS, float); KS = np.ascontiguousarray(KS, float)
    R1, R0, _ = radial_operators(P, a, nq)
    block = lambda k, l: R1[k, l] * MS + R0[k, l] * KS
    A11 = block(0, 0)
    if P == 1:
        return 0.5 * (A11 + A11.T)
    b = list(range(1, P))
    A1b = np.hstack([block(0, l) for l in b])
    Abb = np.vstack([np.hstack([block(k, l) for l in b]) for k in b])
    S = A11 - A1b @ np.linalg.solve(Abb, A1b.T)
    return 0.5 * (S + S.T)


# ===========================================================================
# NGSolve-native coupling: compound space + boundary IE terms + exterior recovery
# ===========================================================================
def ie_compound_space(mesh, P, order=2, definedon=None):
    """Compound FESpace ``[interior H1, S_1, ..., S_{P-1}]`` for the IE: the interior H1 (its boundary
    trace is radial level 0) plus ``P-1`` boundary-only H1 "bubble" levels on the truncation surface.
    ``definedon`` restricts the surface levels (default: all boundaries).
    """
    import ngsolve as ng
    if definedon is None:
        definedon = mesh.Boundaries(".*")
    Vint = ng.H1(mesh, order=order)
    Ss = [ng.H1(mesh, order=order, definedon=definedon) for _ in range(P - 1)]
    return ng.FESpace([Vint] + Ss)


def add_exterior_ie(a_bf, X, P, a=1.0, nq=160, definedon=None):
    """Add the static IE exterior energy to a BilinearForm ``a_bf`` on the compound space ``X``
    (from :func:`ie_compound_space`).  Appends the boundary integrals
    ``sum_kl [ (R1_kl/a^2) int_Gamma u_k u_l ds + R0_kl int_Gamma grad_S u_k . grad_S u_l ds ]``
    (level 0 = the interior trace, levels 1..P-1 = the surface bubbles), so the augmented SPARSE
    system closes the exterior.  The caller adds the interior physics + RHS to the same ``X``.

    ``definedon`` restricts the IE to a specific truncation boundary (e.g. an OUTER sphere enclosing a
    body, so an internal material interface is NOT treated as the open boundary).  It MUST match the
    ``definedon`` passed to :func:`ie_compound_space` (the bubble levels live there).  Default = all
    boundaries.
    """
    import ngsolve as ng
    R1, R0, _ = radial_operators(P, a, nq)
    trial = X.TrialFunction(); test = X.TestFunction()
    n = ng.specialcf.normal(X.mesh.dim)
    inv_a2 = 1.0 / (a * a)
    ds = ng.ds if definedon is None else ng.ds(definedon=definedon)

    def val(fns, k):
        return fns[0].Trace() if k == 0 else fns[k]

    def sgrad(fns, k):
        g = ng.grad(fns[k]).Trace()
        return g - (g * n) * n        # tangential (surface) gradient; boundary levels already tangential

    for k in range(P):
        for l in range(P):
            if R1[k, l] != 0.0:
                a_bf += (inv_a2 * R1[k, l]) * val(trial, k) * val(test, l) * ds
            if R0[k, l] != 0.0:
                a_bf += R0[k, l] * sgrad(trial, k) * sgrad(test, l) * ds


def exterior_field(gf, P, a, points):
    """Evaluate the IE exterior field ``phi(x) = sum_k U_k(s) N_k(a/|x|)`` at exterior ``points``
    (|x| >= a), where ``U_k(s)`` is level-k of the solved compound GridFunction ``gf`` at the surface
    point ``s = a x/|x|``.  Returns a 1-D ndarray (one value per point).  This makes the exterior
    field available EVERYWHERE from the NGSolve solution (the IE is not interior-only).
    """
    import ngsolve as ng
    mesh = gf.space.mesh
    pts = np.atleast_2d(np.asarray(points, float))
    out = np.empty(len(pts))
    for i, x in enumerate(pts):
        r = float(np.linalg.norm(x))
        s = x / r * a                          # project onto the truncation sphere (radius a)
        mp = mesh(float(s[0]), float(s[1]), float(s[2]), ng.BND)
        Nv, _ = _radial_eval(P, np.array([a / r]))
        out[i] = sum(float(gf.components[k](mp)) * Nv[k, 0] for k in range(P))
    return out
