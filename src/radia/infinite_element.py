# -*- coding: utf-8 -*-
"""radia.infinite_element -- static / low-frequency INFINITE-ELEMENT open-boundary closure.

The (static) infinite element (IE) closes an exterior Laplace problem on a truncation surface by
expanding the decaying exterior field as (surface FE) x (radial decay basis) and statically
condensing the radial DOFs onto the surface trace, giving a DtN surface stiffness ``S_Gamma`` that
the caller adds to the interior FE system.  On a SPHERE this is identical to the Kelvin
transformation (same exterior polynomial space; see
``examples/kelvin_transformation/DtN_spectrum/act7_28_ie_vs_kelvin_fair_dtn.py``).

The numerical kernel is C++ (``src/core/rad_infinite_element.cpp``, the port of the Python prototypes
``act7_32``/``act7_33``): the well-conditioned radial decay operators (R1 radial stiffness, R0 radial
mass, in the vertex + integrated-Legendre basis -- NEVER naive monomials, per the act7_28 conditioning
lesson) and the static condensation.  The generic surface FE matrices (mass + Laplace-Beltrami) come
from NGSolve (this module does not reimplement them -- "complement NGSolve").

Public API
----------
``radial_operators(P, a=1.0)``           -> (R1, R0, g)   the P x P radial operators + trace vector.
``dtn_surface_operator(Mtil, Ktil, P, a)`` -> S            condense given unit-sphere surface matrices.
``surface_dtn_from_mesh(mesh, P, a, order)`` -> (S, bnd)   the full closure on an NGSolve mesh's
                                                            spherical boundary (lazy NGSolve import).
"""
import numpy as np

__all__ = ["radial_operators", "dtn_surface_operator", "surface_dtn_from_mesh"]


def radial_operators(P, a=1.0, nq=160):
    """Orthogonal nodal radial decay operators of the static IE (sphere radius ``a``).

    Returns ``(R1, R0, g)``: ``R1`` (P x P) radial stiffness ``int rho_k' rho_l' r^2 dr``, ``R0``
    (P x P) radial mass ``int rho_k rho_l dr`` in the well-conditioned vertex + integrated-Legendre
    basis, and the trace vector ``g`` (= e_1).
    """
    from radia import _radia_pybind as _rp
    d = _rp._ie_radial_operators(int(P), float(a), int(nq))
    R1 = np.asarray(d["R1"], float).reshape(P, P)
    R0 = np.asarray(d["R0"], float).reshape(P, P)
    g = np.asarray(d["g"], float)
    return R1, R0, g


def dtn_surface_operator(Mtil, Ktil, P, a=1.0, nq=160):
    """Condensed DtN surface stiffness ``S_Gamma`` (N x N) of the static IE.

    ``Mtil`` / ``Ktil`` are the UNIT-SPHERE surface mass and Laplace-Beltrami matrices (N x N,
    symmetric).  For a truncation sphere of radius ``a``, ``Mtil = M_physical / a^2`` and
    ``Ktil = K_physical`` (see :func:`surface_dtn_from_mesh`, which does the scaling).  Builds the
    P-level tensor blocks ``R1_kl*Mtil + R0_kl*Ktil`` and condenses the radial bubble levels onto
    the trace.  ``eig(S, Mtil) -> the analytic Steklov ladder (n+1)/a``.
    """
    from radia import _radia_pybind as _rp
    Mtil = np.ascontiguousarray(Mtil, float)
    Ktil = np.ascontiguousarray(Ktil, float)
    N = Mtil.shape[0]
    if Mtil.shape != (N, N) or Ktil.shape != (N, N):
        raise ValueError(f"Mtil/Ktil must be square N x N (got {Mtil.shape}, {Ktil.shape})")
    d = _rp._ie_dtn_operator(Mtil.ravel().tolist(), Ktil.ravel().tolist(),
                             int(N), int(P), float(a), int(nq))
    if d["info"] != 0:
        raise RuntimeError(f"IE condensation LAPACK solve failed (info={d['info']})")
    return np.asarray(d["S"], float).reshape(N, N)


def surface_dtn_from_mesh(mesh, P, a=1.0, order=None, nq=160):
    """Build the IE DtN closure on the SPHERICAL boundary of an NGSolve volume ``mesh`` (radius ``a``).

    Assembles the boundary surface mass ``M^S`` and Laplace-Beltrami ``K^S`` (NGSolve, via the
    tangential projection of the volume H1 gradient trace), restricts to the boundary DOFs, scales to
    the unit sphere (``Mtil = M^S/a^2``, ``Ktil = K^S``), and condenses the radial DOFs in C++.

    Returns ``(S, bnd)``: the dense DtN surface stiffness ``S`` (n_bnd x n_bnd, the exterior energy to
    add to the interior system) and the boundary DOF indices ``bnd`` (into the volume H1 space).
    ``order`` defaults to the mesh's curve order / 2 if available, else 2.  The caller wraps the
    NGSolve assembly in ``with TaskManager():`` per the lab TaskManager policy.
    """
    import ngsolve as ng
    import scipy.sparse as sp

    if order is None:
        order = 2
    fes = ng.H1(mesh, order=order)
    u, v = fes.TnT()
    n = ng.specialcf.normal(3)
    bm = ng.BilinearForm(fes, symmetric=True, check_unused=False)
    bm += u * v * ng.ds
    bm.Assemble()
    gu, gv = ng.grad(u).Trace(), ng.grad(v).Trace()
    gut = gu - (gu * n) * n
    gvt = gv - (gv * n) * n
    bk = ng.BilinearForm(fes, symmetric=True, check_unused=False)
    bk += (gut * gvt) * ng.ds
    bk.Assemble()

    def _csr(m):
        r, c, val = m.COO()
        return sp.csr_matrix((np.asarray(val), (np.asarray(r), np.asarray(c))),
                             shape=(m.height, m.height))

    bnd_ba = fes.GetDofs(mesh.Boundaries(".*"))
    bnd = np.array([i for i in range(fes.ndof) if bnd_ba[i]], dtype=int)
    MS = _csr(bm.mat)[np.ix_(bnd, bnd)].toarray()
    KS = _csr(bk.mat)[np.ix_(bnd, bnd)].toarray()
    MS = 0.5 * (MS + MS.T)
    KS = 0.5 * (KS + KS.T)
    Mtil = MS / (a * a)        # unit-sphere mass (physical/a^2); Ktil = K^S (a-factors cancel)
    S = dtn_surface_operator(Mtil, KS, P, a=a, nq=nq)
    return S, bnd
