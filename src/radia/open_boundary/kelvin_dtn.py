# -*- coding: utf-8 -*-
"""Kelvin-built exterior DtN open boundary -> CLN (the material-aware / non-separable
companion of the closed-form `dtn_cln`).

WHAT THIS IS.  Where `dtn_cln` realises the EXACT closed-form DtN of a separable
(spherical) homogeneous exterior, this module **BUILDS** the exterior DtN by a
Kelvin-transformed FEM and Schur-condensation -- so it carries an arbitrary
(NON-separable) truncation shape AND a MATERIAL exterior (iron shield / layered /
inhomogeneous), which a closed-form symbol cannot.  The physical DtN ladder is the
generalised Steklov eigenproblem `(S, M_Gamma)`; for a non-separable shape it is a
convergent BAND approximation in `q = sqrt(s)` (a few-stage CLN), verified by the
point-group SPLITTING of the sphere's degeneracies (square: C4v; cube: O_h).

HONEST PROVENANCE (the crux -- 3 layers, do not overclaim; full record in the MCP
topic kelvin_transformation(topic="material_exterior")):
  * Kelvin open boundary with an AIR exterior + arbitrary INTERIOR material/anisotropy
    is CLASSICAL -- Freeman & Lowther (IEEE T-Magn 1988/89); FEMM ships it.
  * Transforming sigma / eps / mu under a coordinate (conformal) map is CLASSICAL --
    transformation optics / Ward-Pendry 1996 (incl. DC-sigma cloaks); a PML is the
    same CTM family.  So "transform sigma under the Kelvin map" is NOT new in the
    abstract.
  * Sugahara's OWN validated FUSIONS are the contributions: (i) Kelvin INVERSION as an
    EXACT OPEN BOUNDARY with the sigma-CONFORMAL transform so a CONDUCTOR crosses the
    truncation, eddy-current-testing-validated (IEEE Magnetics 2022 -- the formulation
    basis for the (a/r)^4 sigma / (a/r)^2 mu weights this module's nu, sigma carry);
    (ii) the Kelvin material-aware DtN as an INVERSE-DESIGN kernel (SF-with-iron, conf
    ~0.83, self-check pending).  This module is the verified reusable operator, not a
    paper claim.

VERIFIED (`validation_test/open_boundary/test_kelvin_dtn.py`, ported from the
Kelvin DtN research demos and maintained as executable validation):
  * the radial Kelvin (R/rho')^2-weighted ball reproduces the closed-form
    `dtn_cln.eddy_dtn` per multipole (the "Kelvin BUILDS the exact DtN" check);
  * the generalised Steklov ladder of a NON-separable cube is O_h-split (the l=2
    quintet -> E_g doublet + T_2g triplet = 2+3, the l=1 dipole stays a degenerate
    triplet) and mesh-convergent -- an analytic-free correctness proof;
  * a few-stage CLN in sqrt(s) reduces the built DtN over the DC->evanescent band.

SCOPE.  Same island as `dtn_cln`: COMPACT / quasi-spherical MQS.  Kelvin is
sphere-locked (Liouville), so an elongated truncation wastes the spherical shell
(a box CFS-PML hugs better); genuine wave radiation is outside radia's MQS scope.
See docs/open_boundary/OPEN_BOUNDARY_MAP.md.
"""
import numpy as np

__all__ = [
    "kelvin_fem_radial_dtn",
    "kelvin_dtn_matrix",
    "steklov_spectrum",
    "band_cln_fit",
]

# 3-point Gauss-Legendre on [-1, 1]
_GP = np.array([-np.sqrt(3.0 / 5.0), 0.0, np.sqrt(3.0 / 5.0)])
_GW = np.array([5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0])


# ---------------------------------------------------------------------------
# radial Kelvin-FEM build of the exact eddy DtN per multipole (pure numpy)
# ---------------------------------------------------------------------------
def _assemble_inner(nodes, n):
    """Inner radial K (stiffness) + M (mass), r^2-measure: energy
    int ( r^2 u'^2 + n(n+1) u^2 + s r^2 u^2 ) dr over [R0, Rmid]."""
    N = nodes.size
    K = np.zeros((N, N))
    M = np.zeros((N, N))
    cent = n * (n + 1)
    for e in range(N - 1):
        a, b = nodes[e], nodes[e + 1]
        d = b - a
        rg = 0.5 * (a + b) + 0.5 * d * _GP
        jac = 0.5 * d
        N0 = (b - rg) / d
        N1 = (rg - a) / d
        dN = (-1.0 / d, 1.0 / d)
        Ns = (N0, N1)
        for p in range(2):
            for qq in range(2):
                K[e + p, e + qq] += np.sum(_GW * jac * (rg ** 2 * dN[p] * dN[qq] + cent * Ns[p] * Ns[qq]))
                M[e + p, e + qq] += np.sum(_GW * jac * rg ** 2 * Ns[p] * Ns[qq])
    return K, M


def _assemble_kelvin(nodes, n, R):
    """Kelvin-ball static stiffness for the compactified exterior [Rmid, inf) -> the
    ball r' in [0, Rmid].  3-D Kelvin maps the exterior energy to the ball with the
    conformal MATERIAL weight mu' = (R/r')^2 (radia Omega/H1 convention):
        energy = int R^2 [ v'^2 + n(n+1)/r'^2 v^2 ] dr'.
    (The 1/r'^2 term is integrable with the GND v(0)=0 at the Kelvin centre.)"""
    N = nodes.size
    K = np.zeros((N, N))
    cent = n * (n + 1)
    for e in range(N - 1):
        a, b = nodes[e], nodes[e + 1]
        d = b - a
        rg = 0.5 * (a + b) + 0.5 * d * _GP
        jac = 0.5 * d
        N0 = (b - rg) / d
        N1 = (rg - a) / d
        dN = (-1.0 / d, 1.0 / d)
        Ns = (N0, N1)
        for p in range(2):
            for qq in range(2):
                K[e + p, e + qq] += np.sum(_GW * jac * R ** 2
                                           * (dN[p] * dN[qq] + cent * Ns[p] * Ns[qq] / rg ** 2))
    return K


def kelvin_fem_radial_dtn(n, s, R0=1.0, Rmid=3.0, h_in=0.01, h_kel=0.02):
    """Eddy DtN eigenvalue at R0 for multipole n, BUILT by a radial Kelvin-FEM:
    inner (K + s M) on [R0, Rmid] + a (R/rho')^2-weighted Kelvin-ball static tail
    on r' in [0, Rmid] (the compactified exterior), GND at the Kelvin centre.
    Reproduces the closed-form `dtn_cln.eddy_dtn` -- the "Kelvin BUILDS the exact
    DtN" check (separable / homogeneous case)."""
    inner = np.linspace(R0, Rmid, int(round((Rmid - R0) / h_in)) + 1)
    Ki, Mi = _assemble_inner(inner, n)
    Ni = inner.size
    kel = np.linspace(0.0, Rmid, int(round(Rmid / h_kel)) + 1)
    Kk = _assemble_kelvin(kel, n, Rmid)
    Nk = kel.size
    Ng = Ni + Nk - 1
    A = np.zeros((Ng, Ng), dtype=complex)
    A[:Ni, :Ni] += Ki + complex(s) * Mi
    kmap = np.empty(Nk, dtype=int)
    kmap[:Nk - 1] = np.arange(Ni, Ni + Nk - 1)
    kmap[Nk - 1] = Ni - 1                                  # shared Rmid interface node
    for i in range(Nk):
        for j in range(Nk):
            A[kmap[i], kmap[j]] += Kk[i, j]
    gnd = Ni                                               # Kelvin centre v(0)=0 (GND)
    fixed = [0, gnd]
    free = [k for k in range(Ng) if k not in fixed]
    u = np.zeros(Ng, dtype=complex)
    u[0] = 1.0                                             # u(R0)=1 (Dirichlet trace)
    rhs = -A[np.ix_(free, fixed)] @ u[fixed]
    u[free] = np.linalg.solve(A[np.ix_(free, free)], rhs)
    reaction = A[0, :] @ u                                 # weak-form flux at R0
    return -complex(reaction)


# ---------------------------------------------------------------------------
# arbitrary-shape, MATERIAL-aware Kelvin/FEM DtN matrix via Schur (needs NGSolve)
# ---------------------------------------------------------------------------
def kelvin_dtn_matrix(mesh, order, s, gamma="gamma", dirichlet="outer",
                      nu=1.0, sigma=1.0):
    """Build the eddy DtN MATRIX on the truncation surface `gamma` by Schur-condensing
    the Kelvin/FEM exterior A(s) = K + s M onto Gamma, with K = int nu grad.grad and
    M = int sigma u v (so a non-unit / spatially-varying `nu`, `sigma` carries the
    Kelvin (R/rho')^2 weight and/or a MATERIAL exterior).  Returns (S, Mg, g_idx):
    S = Schur complement (the eddy DtN flux functional), Mg = the Gamma boundary mass,
    g_idx = the Gamma DOF indices.  The physical DtN ladder is the GENERALISED
    eigenproblem (S, Mg) -- see steklov_spectrum (NOT the raw eigenvalues of S).

    Requires NGSolve (imported lazily) + scipy.  Caller wraps in TaskManager."""
    from ngsolve import H1, BilinearForm, grad, dx, ds
    import scipy.sparse as sp
    import scipy.sparse.linalg as spla

    fes = H1(mesh, order=order, dirichlet=dirichlet)
    u, v = fes.TnT()
    K = BilinearForm(nu * grad(u) * grad(v) * dx); K.Assemble()
    M = BilinearForm(sigma * u * v * dx); M.Assemble()
    Mgf = BilinearForm(u * v * ds(gamma), check_unused=False); Mgf.Assemble()
    gam = fes.GetDofs(mesh.Boundaries(gamma))
    free = fes.FreeDofs()
    g_idx = [i for i in range(fes.ndof) if gam[i] and free[i]]
    i_idx = [i for i in range(fes.ndof) if free[i] and not gam[i]]

    def _csr(mat):
        i, j, val = mat.COO()
        return sp.csr_matrix((np.array(val), (np.array(i), np.array(j))), shape=(fes.ndof, fes.ndof))

    A = _csr(K.mat) + complex(s) * _csr(M.mat)
    Agg = A[np.ix_(g_idx, g_idx)].toarray()
    Agi = A[np.ix_(g_idx, i_idx)].toarray()
    Aii = A[np.ix_(i_idx, i_idx)].tocsc()
    Aig = A[np.ix_(i_idx, g_idx)]
    lu = spla.splu(Aii)
    S = Agg - Agi @ lu.solve(Aig.toarray())
    Mg = _csr(Mgf.mat)[np.ix_(g_idx, g_idx)].toarray().real
    return S, Mg, g_idx


def steklov_spectrum(S, Mg):
    """The physical DtN ladder of a built DtN matrix: the GENERALISED eigenproblem
    (sym(S), Mg) with the Gamma boundary mass.  Returns (eigenvalues ascending,
    eigenvectors).  For a non-separable shape the spectrum splits the sphere's
    l-fold degeneracies by the truncation's point group (square C4v, cube O_h)."""
    from scipy.linalg import eigh
    S = np.asarray(S)
    return eigh(0.5 * (S + S.T).real, np.asarray(Mg).real)


def band_cln_fit(s_band, dtn_values, stages):
    """Fit a built (non-separable) DtN over a band to a rational function in
    q = sqrt(s) of the given number of `stages` (the convergent band-CLN).  Returns
    (fit_values, nrmse).  Exact closed-form ladders (separable) are in dtn_cln; this
    is the NON-separable convergent approximation."""
    s_band = np.asarray(s_band, dtype=complex)
    G = np.asarray(dtn_values, dtype=complex)
    q = np.sqrt(s_band)
    Vd = np.vstack([q ** k for k in range(stages)]).T
    Amat = np.hstack([Vd, -(G[:, None]) * Vd[:, 1:]])
    coef, *_ = np.linalg.lstsq(np.vstack([Amat.real, Amat.imag]),
                               np.concatenate([(G * Vd[:, 0]).real, (G * Vd[:, 0]).imag]),
                               rcond=None)
    fit = (Vd @ coef[:stages]) / (Vd @ np.concatenate([[1.0], coef[stages:]]))
    nrmse = float(np.sqrt(np.mean(np.abs(fit - G) ** 2)) / np.sqrt(np.mean(np.abs(G) ** 2)))
    return fit, nrmse
