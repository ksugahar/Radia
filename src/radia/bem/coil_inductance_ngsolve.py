"""Coil self-inductance via ngsolve.bem LaplaceSL on HDivSurface RT₀.

Production replacement for the intree (Python) BEM-A path retired
2026-05-03 after benchmarking showed ngsolve.bem was 50-60x faster
than our pure-Python double-loop assembler at every measured N (302
to 1014 triangles).  ngsolve.bem stores the dense Galerkin matrix
internally and exposes it via COO for O(1) extraction; intree's
``radia.bem.efie_rwg.solve_inductance_source_sink_intree`` was deleted.

Solves the constrained impedance-EFIE on the coil surface (the SOLE
AC formulation since 2026-07-02):

  [jw*mu0*SL + Z_s*M   D^T] [J]   [0]
  [D                    0 ] [p] = [g]      (complex, omega > 0)

where SL = LaplaceSL (single-layer BEM operator on surface currents),
M = HDivSurface mass matrix carrying the complex Leontovich surface
impedance Z_s = (1+j)/(sigma*delta), D = divergence matrix
(HDivSurface -> SurfaceL2), g = source/sink current injection
(+1/A_src at source, -1/A_snk at sink), J = the FINITE-IMPEDANCE
surface current, p = Lagrange multiplier enforcing current
conservation.  At omega == 0 the system reduces to the real vacuum
saddle [SL, D^T; D, 0] (DC limit, R = 0).

External inductance: ``L = mu_0 * (J^H @ SL @ J)``.
AC SIBC resistance:  ``R = Re(Z_s) * (J^H @ M @ J)``.

The historical PEC formulation (real perfect-conductor saddle, R
evaluated post-hoc from the PEC current) was REMOVED 2026-07-02: it
over-estimated R ~3x on tightly-wound coils because the PEC J
concentrates singularly at near-contact gaps/edges where the
Leontovich integral breaks down (kubota 3-turn: 15.14 mOhm vs the
physical 4.63).  See docs/peec/VOLUME_PEEC_DESIGN.md.

Per-triangle J sampling: ``compute_centroids_areas_J(...)`` returns
the centroid + area + averaged J vector per BND triangle, ready to
feed into ``radia.bem_sibc_solver.compute_phi_inc_from_surface_J`` for
the workpiece weak-coupling bridge.

Requires: NGSolve >= 6.2.2603 with ngsolve.bem, scipy.
"""
from __future__ import annotations

import time

import numpy as np
from scipy.linalg import solve as scipy_solve
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import minres as scipy_minres
from scipy.sparse.linalg import gmres as scipy_gmres


MU_0 = 4e-7 * np.pi


def _solve_saddle(K, rhs, method, tol=1e-8, maxiter=2000):
    """Solve the saddle-point system K x = rhs.

    method:
      "lu"     dense LAPACK direct solve (scipy.linalg.solve, partial-
               pivot LU) with overwrite_a=True so LAPACK reuses K's
               storage for the LU factors (saves N^2 memory vs the
               default copy).  Best for ndof <~ 25000.  Memory O(N^2),
               work O(N^3).
      "minres" Krylov solver for REAL symmetric indefinite systems
               (scipy.sparse.linalg.minres) -- valid ONLY for the
               omega=0 real vacuum-L saddle.  REJECTED (ValueError)
               for the complex AC impedance-EFIE saddle: minres
               assumes a Hermitian operator, but that system is
               complex-symmetric (K^T = K, not Hermitian), so minres
               silently returns a wrong solution.
      "gmres"  General complex-capable Krylov
               (scipy.sparse.linalg.gmres).  The large-N choice for
               the AC impedance-EFIE saddle.

    NOTE: with method="lu", K is overwritten with its LU factors in
    place.  Callers must not reuse K after this call.
    """
    if method == "minres" and np.iscomplexobj(K):
        # scipy MINRES assumes a real-symmetric / Hermitian operator.
        # The impedance-EFIE saddle is complex-SYMMETRIC (K^T = K, NOT
        # Hermitian), so minres silently returns a WRONG solution
        # (measured: for K=(1+1j)*I it returns (K^H)^-1 b, not K^-1 b,
        # emitting only ComplexWarnings).  Fail fast per CLAUDE.md
        # "No Fallbacks"; use lu or gmres.
        raise ValueError(
            "MINRES cannot solve the complex impedance-EFIE saddle: "
            "scipy minres assumes a Hermitian operator, but this system "
            "is complex-symmetric (K^T = K), so it silently returns a "
            "wrong solution.  Use solver='lu' (dense direct) or "
            "solver='gmres'.")
    if method == "lu":
        x = scipy_solve(K, rhs, overwrite_a=True)
        # K has been overwritten -- cannot compute K@x for residual cheaply
        # without keeping a copy.  Skip the residual check for "lu".
        info = {"method": "lu", "iterations": 1, "residual": 0.0}
        return x, info
    if method == "minres":
        x, code = scipy_minres(K, rhs, rtol=tol, maxiter=maxiter)
        if code != 0:
            raise RuntimeError(
                f"MINRES did not converge (info={code}, "
                f"maxiter={maxiter}, tol={tol})")
        info = {"method": "minres", "iterations": maxiter,
                "residual": float(np.linalg.norm(K @ x - rhs))}
        return x, info
    if method == "gmres":
        x, code = scipy_gmres(K, rhs, rtol=tol, maxiter=maxiter,
                              restart=50)
        if code != 0:
            raise RuntimeError(
                f"GMRES did not converge (info={code}, "
                f"maxiter={maxiter}, tol={tol})")
        info = {"method": "gmres", "iterations": maxiter,
                "residual": float(np.linalg.norm(K @ x - rhs))}
        return x, info
    raise ValueError(
        f"Unknown saddle-point solver: {method!r}.  "
        f"Choices: lu, minres, gmres.")


def _to_dense(mat):
    """Extract dense NumPy array from NGSolve BaseMatrix via COO.

    NGSolve BEM operators (LaplaceSL etc.) store the full dense matrix
    internally as a SparseMatrix with 100% fill.  The built-in
    ToDense() is O(N) slower than necessary (~144 s vs 0.06 s at
    N=5085) because it performs N column-by-column MatVecs instead
    of a direct memory copy.

    This function extracts the COO triplets and converts via scipy,
    which is ~2500x faster.  Note: requires the operator to have
    been constructed with ``use_fmm=False`` so the underlying
    storage is dense (FMM = lazy operator, no COO).
    """
    rows, cols, vals = mat.COO()
    return coo_matrix((vals, (rows, cols)),
                      shape=(mat.height, mat.width)).toarray()


_LOOP_COCR_TOL = 1e-8
_LOOP_COCR_MAXITER = 2000


def _cocr(matvec, b, tol=_LOOP_COCR_TOL, maxiter=_LOOP_COCR_MAXITER):
    """COCR (Sogabe-Zhang 2007) for complex-SYMMETRIC A = A^T (NOT Hermitian).

    Uses UNCONJUGATED inner products (numpy ``@`` on 1-D arrays does not
    conjugate), the correct short-recurrence Krylov method for the complex-
    symmetric loop operator ``Pi A11 Pi``.  Contrast: scipy GMRES restart
    cripples this system (~1e5 matvecs), and COCR on the raw indefinite saddle
    breaks down structurally -- the div-free reduction is what makes COCR fit.
    Valid for real-symmetric A too.  Returns ``(x, iters, status)`` with status
    in {converged, breakdown, maxiter}.  1 matvec per iteration.
    """
    nb = float(np.linalg.norm(b))
    x = np.zeros_like(b)
    if nb == 0.0:
        return x, 0, "converged"
    r = b.copy()
    Ar = matvec(r)
    p = r.copy()
    Ap = Ar.copy()
    rAr = r @ Ar
    for k in range(maxiter):
        ApAp = Ap @ Ap
        if ApAp == 0:
            return x, k, "breakdown"
        alpha = rAr / ApAp
        x = x + alpha * p
        r = r - alpha * Ap
        if float(np.linalg.norm(r)) / nb < tol:
            return x, k + 1, "converged"
        Ar = matvec(r)
        rAr_new = r @ Ar
        if rAr == 0:
            return x, k + 1, "breakdown"
        beta = rAr_new / rAr
        rAr = rAr_new
        p = r + beta * p
        Ap = Ar + beta * Ap
    return x, maxiter, "maxiter"


def _edge_midpoint_coords(mesh, n_J):
    """Edge-midpoint coordinates for the HACApK cluster tree, aligned with the
    RT0 HDivSurface DOF numbering (DOF i == edge.nr i, verified for order-0
    HDivSurface).  Only valid for fes_order == 0; raises otherwise."""
    coords = np.zeros((n_J, 3))
    n_edges = 0
    for e in mesh.edges:
        vs = [mesh.vertices[v.nr].point for v in mesh[e].vertices]
        coords[e.nr] = 0.5 * (np.asarray(vs[0]) + np.asarray(vs[1]))
        n_edges += 1
    if n_edges != n_J:
        raise ValueError(
            f"HACApK loop matvec needs fes_order==0 (RT0): got {n_edges} "
            f"edges but n_J={n_J}.  Use solver='cocr' (dense matvec) for "
            "higher order.")
    return coords


class _LoopReducedSaddle:
    """Build-once / solve-many divergence-free reduction of the BEM saddle.

    Reduces ``[[A11, D_red^T],[D_red, 0]] [J;p] = [rhs_J; g_red]`` to the
    divergence-free (loop / stream-function) subspace via the sparse
    orthogonal projector onto ker(D_red)::

        Pi = I - D_red^T (D_red D_red^T)^-1 D_red            (sparse splu)
        J  = J_p + J_loop,  J_p = D_red^T (D_red D_red^T)^-1 g_red,
        (Pi A11 Pi) J_loop = Pi (rhs_J - A11 J_p)            (J_loop in ker D_red)

    with ``A11 = jw*mu0*SL + Z_s*M`` (omega>0) or ``SL`` (omega==0).  ``Pi A11
    Pi`` is complex-symmetric and A11-like (SPD-dominated single layer), so
    COCR converges in ~20-30 MESH-INDEPENDENT iterations, replacing the dense
    saddle LU (O(N^3) work, O(N^2) memory) with ~24 matvecs.  The reduction is
    EXACT: on the gapped-torus fixture it reproduces the dense-LU R/L to <0.01%.

    The projector splu, ``J_p`` and (optionally HACApK-compressed) A11 are all
    built ONCE in ``__init__``; ``solve(rhs_J)`` then reuses them for any
    J-block right-hand side.  This is what the coupled Picard loop
    (CoupledBEMSolver) needs -- the coil H-matrix is built once and re-solved
    each iteration against the workpiece back-reaction ``rhs_J = -f_back`` --
    and what the standalone impedance-EFIE (``_loop_cocr_solve``, rhs_J=0) uses.

    Why the loop reduction (not the raw saddle): COCR breaks down on the
    indefinite saddle -- rhs = [0; g] lives entirely in the constraint block so
    the initial ``r^T A r = 0`` -- and diverges on the Schur complement; the
    div-free subspace is where the short-recurrence method actually fits.
    ``D_red`` (one redundant row dropped) is REQUIRED: the full D is
    rank-deficient and ``(D D^T)^-1`` then amplifies g's null component.

    ``matvec_backend``: ``"dense"`` (``SL @ v``, reuses the already-assembled
    dense SL) or ``"hacapk"`` (compress SL to an O(N log N) H-matrix via
    ``HACApKBEMManager``, the ``bem_sibc_solver`` pattern; fes_order==0 only,
    needs ``coords``).  HACApK MatVec accuracy is ~3e-7 (>> the FMM backend's
    ~1e-3), so it preserves the accuracy-sensitive R.
    """

    def __init__(self, SL, M, D_red, g_red, omega, Z_s, matvec_backend,
                 coords=None, hacapk_aca_eps=1e-8, hacapk_leaf=64,
                 hacapk_eta=2.0, log_fn=None):
        import scipy.sparse as sp
        import scipy.sparse.linalg as spla
        _log = log_fn if log_fn is not None else (lambda _t, _m: None)
        self.ac = omega > 0.0
        jwm = 1j * omega * MU_0

        if matvec_backend == "hacapk":
            if coords is None:
                raise ValueError(
                    "hacapk loop matvec requires edge-midpoint coords")
            from radia import _radia_pybind as _rpb
            t0 = time.perf_counter()
            SL_h = _rpb.HACApKBEMManager(np.ascontiguousarray(coords),
                                         np.ascontiguousarray(SL))
            SL_h.BuildHMatrix(aca_eps=hacapk_aca_eps,
                              leaf_size=int(hacapk_leaf),
                              eta=hacapk_eta, max_rank=-1, print_level=0)
            st = SL_h.GetStats()
            _log("BEMA",
                f"loop-COCR HACApK H-matrix built "
                f"(compression={st['compression']:.3f}, O(N log N) matvec, "
                f"{time.perf_counter()-t0:.1f}s)")

            def _sl(v):
                return (SL_h.MatVec(np.ascontiguousarray(v.real))
                        + 1j * SL_h.MatVec(np.ascontiguousarray(v.imag)))
        elif matvec_backend == "dense":
            def _sl(v):
                return SL @ v
        else:
            raise ValueError(
                f"unknown loop matvec backend {matvec_backend!r} "
                f"(choices: dense, hacapk)")

        if self.ac:
            def _a11(v):
                return jwm * _sl(v) + Z_s * (M @ v)
        else:
            def _a11(v):
                return _sl(v)
        self._a11 = _a11

        Dr = sp.csr_matrix(D_red)
        Drt = Dr.T.tocsr()
        lu_DDt = spla.splu((Dr @ Drt).tocsc())

        def _proj(v):
            Dv = Dr @ v
            if np.iscomplexobj(v):
                y = lu_DDt.solve(Dv.real) + 1j * lu_DDt.solve(Dv.imag)
            else:
                y = lu_DDt.solve(Dv)
            return v - (Drt @ y)
        self._proj = _proj

        self.J_p = Drt @ lu_DDt.solve(g_red)
        if self.ac:
            self.J_p = self.J_p.astype(complex)
        self._c0 = _a11(self.J_p)                 # A11 J_p, reused every solve

        def _mop(v):
            return _proj(_a11(_proj(v)))
        self._mop = _mop
        self._backend = matvec_backend

    def a11(self, v):
        """Apply the (H-matrix or dense) A11 operator -- used for the
        magnetic energy ``J^H A11 J`` in the coupled solver."""
        return self._a11(v)

    def solve(self, rhs_J=None, include_particular=True,
              tol=_LOOP_COCR_TOL, maxiter=_LOOP_COCR_MAXITER):
        """Solve the saddle for a J-block right-hand side ``rhs_J`` (default 0).

        ``include_particular`` selects the constraint-block RHS: True (default)
        uses ``g_red`` (the terminal drive -> adds the particular solution
        ``J_p``); False uses ``0`` (a zero-net-current solve, e.g. the imaginary
        back-reaction current whose terminal current is zero).

        Returns ``(J, iters, status)``; raises RuntimeError on non-convergence
        (No-Fallbacks: a stall means a degenerate mesh or mislabelled ports).
        """
        if include_particular:
            b_red = -self._proj(self._c0)
        else:
            b_red = np.zeros_like(self.J_p)
        if rhs_J is not None:
            b_red = b_red + self._proj(rhs_J)
        J_loop, iters, status = _cocr(self._mop, b_red, tol=tol,
                                      maxiter=maxiter)
        if status != "converged":
            raise RuntimeError(
                f"loop-COCR did not converge (status={status!r}, "
                f"iters={iters}, tol={tol}).  The loop operator is "
                f"SPD-dominated and converges in ~20-30 iters on a well-formed "
                f"coil surface; a stall indicates a degenerate mesh or "
                f"mislabelled source/sink.")
        J = (self.J_p + J_loop) if include_particular else J_loop
        return J, iters, status


def _loop_cocr_solve(SL, M, D_red, g_red, omega, Z_s, matvec_backend,
                     coords=None, tol=_LOOP_COCR_TOL,
                     maxiter=_LOOP_COCR_MAXITER, hacapk_aca_eps=1e-8,
                     hacapk_leaf=64, hacapk_eta=2.0, log_fn=None):
    """Loop-reduced COCR solve of the impedance-EFIE saddle (rhs = [0; g_red]).

    Thin wrapper over :class:`_LoopReducedSaddle` for the standalone coil
    solver.  See that class for the full derivation.  Returns ``(J, info)``.
    """
    saddle = _LoopReducedSaddle(
        SL, M, D_red, g_red, omega, Z_s, matvec_backend, coords=coords,
        hacapk_aca_eps=hacapk_aca_eps, hacapk_leaf=hacapk_leaf,
        hacapk_eta=hacapk_eta, log_fn=log_fn)
    J, iters, _status = saddle.solve(tol=tol, maxiter=maxiter)
    if not saddle.ac:
        J = J.real.copy()
    info = {"method": f"cocr[{matvec_backend}]",
            "iterations": int(iters), "residual": 0.0,
            "matvec": matvec_backend}
    return J, info


def compute_inductance_source_sink(
        mesh, source_label="source", sink_label="sink",
        fes_order=0, solver="lu", omega=0.0, Z_s_complex=None,
        log_fn=None):
    """Coil impedance via the ngsolve.bem impedance-EFIE saddle (SOLE formulation).

    For ``omega > 0`` this solves the **impedance-EFIE**: the Leontovich
    surface impedance sits INSIDE the saddle system, so the recovered J
    is the *finite-impedance* surface current::

        [ jω μ0 SL + Z_s M   D^T ] [J]   [0]
        [ D                   0  ] [p] = [g]      (complex)

    with R = Re(Z_s)·(Jᴴ M J) (Leontovich SIBC dissipation of the
    finite-impedance J) and L = μ0·(Jᴴ SL J) (EXTERNAL inductance; the
    internal surface reactance Im(Z_s)·(Jᴴ M J)/ω is deliberately NOT
    folded into L, keeping L a pure geometry quantity).

    The historical "PEC post-hoc" formulation (real perfect-conductor
    saddle + R = Re(Z_s)·JᵀMJ evaluated afterwards) was REMOVED
    2026-07-02: the perfect-conductor J concentrates singularly at
    near-touching surfaces and edges, where it varies below the skin
    depth and the Leontovich integral breaks down -- on the kubota
    3-turn pancake it over-estimated R 3× (15.14 mΩ vs the physical
    4.63 mΩ confirmed by volume PEEC / perimeter PEEC / analytic
    proximity).  On smooth geometry (isolated straight wire) the
    impedance-EFIE reproduces the closed-form Bessel R (0.16-0.5 %),
    i.e. nothing was lost by the removal.  See
    docs/peec/VOLUME_PEEC_DESIGN.md "2026-07-02 outcome".

    For ``omega == 0`` (DC) the AC impedance degenerates and the solve
    reduces to the real vacuum-inductance saddle ``[SL, D^T; D, 0]``
    with R = 0 -- the physical DC limit, selected by the caller passing
    ``omega=0`` (frequency=0), not a fallback.

    Args:
        mesh: NGSolve Mesh -- must be a PURE SURFACE mesh (no volume
            elements).  A volume mesh's internal tets add saddle null
            modes the D[:-1,:] deflation cannot remove and the LU
            reports a singular matrix; extract the boundary first, as
            the panel does via
            ``surface_mesh_extract._extract_surface_mesh_filtered``
            (see validation_test/bem/test_coil_bem_a_volume_vol.py).
        source_label, sink_label: BND labels for current injection /
            extraction faces.  Set in the OCC face.name = "source"/"sink"
            BEFORE meshing, or rely on the panel's smallest-2-PLANE
            auto-detection.
        fes_order: HDivSurface polynomial order.  0 = RT₀ (= RWG, the
            production setting); 1+ for higher-order is supported by
            NGSolve but not validated against radia goldens here.
        solver: saddle-point solver.
            - "lu" (dense direct, default): O(N^3) LU of the full complex
              saddle; best for small N (<~5000).
            - "cocr": SCALABLE path -- reduce the saddle to the divergence-
              free (loop / stream-function) subspace and solve the complex-
              symmetric ``Pi A11 Pi`` with COCR in ~24 MESH-INDEPENDENT
              iterations (exact vs LU to <0.01 %), using the dense ``SL @ v``
              matvec.  The recommended solver for medium/large coils.
            - "hacapk_cocr": the same COCR, but the SL matvec is an
              O(N log N) ``HACApKBEMManager``-compressed H-matrix (the
              ``bem_sibc_solver`` pattern; MatVec accuracy ~3e-7, fes_order==0
              only).  Identical R/L to "cocr"; cuts the per-iteration matvec
              cost on larger meshes.
            - "gmres": scipy GMRES on the dense saddle (unpreconditioned --
              stalls on large saddles, kept for comparison).
            - "minres" is accepted only for the ``omega == 0`` real solve:
              scipy's MINRES silently discards the imaginary part of a complex
              system, so it is REJECTED for the AC impedance-EFIE (fail fast).
        omega: angular frequency [rad/s].  0 selects the DC vacuum-L
            solve (R = 0).
        Z_s_complex: complex Leontovich surface impedance
            ``Z_s = (1+j)/(σ δ)`` [Ohm/sq]; REQUIRED when ``omega > 0``.

    Returns:
        dict with keys
        ``L``         : float [H], external self-inductance at unit terminal current
        ``R``         : float [Ω], SIBC dissipation of the impedance-EFIE J
                        (0.0 for the omega=0 DC solve)
        ``J``         : (n_J,) HDivSurface coefficients (complex for omega>0)
        ``gf_J``      : GridFunction(HDivSurface) holding Re(J)
        ``gf_J_im``   : GridFunction(HDivSurface) holding Im(J)
                        (all-zero for the omega=0 DC solve) -- sample BOTH
                        to reconstruct the complex per-triangle current
        ``SL``, ``D`` : dense Galerkin matrices (post-processing)
        ``A_source``, ``A_sink`` : port areas [m²]
        ``residual``  : max|D J - g|, machine precision for LU solve
        ``n_J``, ``n_f`` : DOF counts
        ``t_assembly``, ``t_solve``, ``t_total`` : timings [s]
    """
    from ngsolve import (HDivSurface, SurfaceL2, TaskManager, ds, BND,
                         BilinearForm, LinearForm, div, GridFunction)
    from ngsolve.bem import LaplaceSL

    # Optional progress log (default no-op).
    _log = log_fn if log_fn is not None else (lambda _t, _m: None)

    # --- Parameter validation FIRST (before the expensive dense
    #     LaplaceSL assembly): fail fast on a bad omega/Zs contract
    #     instead of burning minutes of assembly first. ---
    if solver not in {"lu", "gmres", "minres", "cocr", "hacapk_cocr"}:
        raise ValueError(
            f"Unknown saddle-point solver: {solver!r}.  "
            "Choices: lu, minres, gmres, cocr, hacapk_cocr.")
    if omega > 0.0:
        if Z_s_complex is None or Z_s_complex == 0:
            raise ValueError(
                "omega > 0 requires the complex Leontovich surface "
                "impedance Z_s_complex = (1+1j)/(sigma*delta) "
                f"(got omega={omega}, Z_s_complex={Z_s_complex!r}).  "
                "For a pure vacuum-L solve pass omega=0.")
    elif omega != 0.0:
        # Negative or NaN omega must not silently degrade to the DC
        # vacuum solve (a sign typo in --frequency would return R=0
        # with Z_s_complex ignored).  Fail fast.
        raise ValueError(
            f"omega must be > 0 (AC impedance-EFIE) or exactly 0 "
            f"(DC vacuum-L solve); got omega={omega!r}.")

    t_start = time.perf_counter()

    fes_J = HDivSurface(mesh, order=fes_order)
    fes_L2 = SurfaceL2(mesh, order=max(0, fes_order - 1))
    n_J = fes_J.ndof
    n_f = fes_L2.ndof
    if solver == "hacapk_cocr" and fes_order != 0:
        raise ValueError(
            "solver='hacapk_cocr' currently requires fes_order==0 (RT0); "
            "use solver='cocr' (dense matvec) for higher-order HDivSurface.")
    _log("BEMA",
        f"FES built: n_J={n_J} (HDivSurface RT{fes_order}), "
        f"n_f={n_f} (SurfaceL2 P{max(0,fes_order-1)})")

    # --- Divergence matrix D: n_f x n_J ---
    _t_D = time.perf_counter()
    u_J = fes_J.TrialFunction()
    q = fes_L2.TestFunction()
    bf_D = BilinearForm(trialspace=fes_J, testspace=fes_L2)
    bf_D += div(u_J.Trace()) * q * ds
    bf_D.Assemble()
    D = _to_dense(bf_D.mat)
    _log("BEMA",
        f"div matrix D assembled "
        f"({n_f}x{n_J}, {time.perf_counter()-_t_D:.1f}s)")

    # --- LaplaceSL matrix: n_J x n_J ---
    t0 = time.perf_counter()
    _log("BEMA",
        f"LaplaceSL assembly start (dense, n_J={n_J}, "
        f"memory ~{(n_J*n_J*8)/1e9:.1f} GB)")
    jt, jv = fes_J.TnT()
    # NOTE: Caller MUST be inside `with TaskManager():` per CLAUDE.md
    # "Caller Wraps, Helper Does NOT" (2026-05-27).
    V_op = LaplaceSL(jt.Trace() * ds, use_fmm=False) * jv.Trace() * ds
    SL = _to_dense(V_op.mat)
    t_assembly = time.perf_counter() - t0
    _log("BEMA", f"LaplaceSL assembled ({t_assembly:.1f}s)")

    # --- Source/sink RHS ---
    f_src = LinearForm(fes_L2)
    f_src += q * ds(source_label)
    f_src.Assemble()
    g_src = f_src.vec.FV().NumPy().copy()
    A_src = np.sum(g_src)

    f_snk = LinearForm(fes_L2)
    f_snk += q * ds(sink_label)
    f_snk.Assemble()
    g_snk = f_snk.vec.FV().NumPy().copy()
    A_snk = np.sum(g_snk)

    if A_src < 1e-30 or A_snk < 1e-30:
        return {"error": (f"Source/sink faces not found "
                          f"(A_src={A_src}, A_snk={A_snk}).  "
                          f"Check BND labels {source_label!r} and "
                          f"{sink_label!r}.  Available: "
                          f"{list(set(mesh.GetBoundaries()))}")}

    g = g_src / A_src - g_snk / A_snk

    # --- Saddle point system (drop one constraint to deflate the
    #     1-D null space of the constant function on closed surfaces) ---
    D_red = D[:-1, :]
    g_red = g[:-1]
    n_constraint = n_f - 1

    t0 = time.perf_counter()
    _saddle_n = n_J + n_constraint
    _zero_c = np.zeros((n_constraint, n_constraint))

    if omega > 0.0:
        # --- Impedance-EFIE (the sole AC formulation): Z_s sits INSIDE
        #     the system so J is the finite-impedance (not perfect-
        #     conductor) current.  The (1,1) block is jω μ0 SL + Z_s M
        #     (complex); the resistive Z_s M term penalises the singular
        #     edge/near-contact concentration that made the removed PEC
        #     post-hoc integral over-estimate R 3x on tightly-wound
        #     coils.  (omega/Z_s_complex contract already validated at
        #     the top of the function, before assembly.) ---
        # Mass matrix on HDivSurface (∫_S jt·jv dS) for the Z_s M block.
        bf_M = BilinearForm(fes_J)
        bf_M += jt.Trace() * jv.Trace() * ds
        bf_M.Assemble()
        M = _to_dense(bf_M.mat)
        if solver in ("cocr", "hacapk_cocr"):
            # Scalable path: reduce the saddle to the div-free (loop)
            # subspace and solve the complex-symmetric Pi A11 Pi with COCR
            # (~24 mesh-independent iters) instead of the O(N^3) dense LU.
            mv = "hacapk" if solver == "hacapk_cocr" else "dense"
            coords = (_edge_midpoint_coords(mesh, n_J)
                      if mv == "hacapk" else None)
            _log("BEMA",
                f"impedance-EFIE {solver} solve (div-free reduction, "
                f"matvec={mv})")
            J, solve_info = _loop_cocr_solve(
                SL, M, D_red, g_red, omega, Z_s_complex, mv,
                coords=coords, log_fn=_log)
        else:
            _log("BEMA",
                f"impedance-EFIE saddle assembly (complex, "
                f"{_saddle_n}x{_saddle_n}, ~{(_saddle_n*_saddle_n*16)/1e9:.1f} GB)")
            A11 = 1j * omega * MU_0 * SL + Z_s_complex * M
            Kc = np.block([
                [A11,                          D_red.T.astype(complex)],
                [D_red.astype(complex),        _zero_c.astype(complex)]
            ])
            del A11   # np.block copied it; free ~n_J^2 x 16 B before the LU
            rhs = np.zeros(n_J + n_constraint, dtype=complex)
            rhs[n_J:] = g_red
            _log("BEMA", f"impedance-EFIE solve start (method={solver})")
            x, solve_info = _solve_saddle(Kc, rhs, method=solver)
            J = x[:n_J]                       # complex current
        t_lu = time.perf_counter() - t0
        JHMJ = float(np.real(np.conj(J) @ M @ J))
        JHSLJ = float(np.real(np.conj(J) @ SL @ J))
        # R = Re(Zs)*(J^H M J): the Leontovich SIBC dissipation of the
        # finite-impedance J -- physical, unlike the removed PEC
        # post-hoc integral.
        R_coil = float(Z_s_complex.real * JHMJ)
        # L = mu0*(J^H SL J): EXTERNAL inductance (pure geometry).  The
        # internal surface reactance Im(Zs)*(J^H M J)/omega is
        # deliberately NOT folded into L: L stays the geometry quantity
        # every downstream consumer and golden expects, and only R
        # carries the conductor physics.
        L = float(MU_0 * JHSLJ)
        residual = float(np.max(np.abs(D @ J - g)))
        _log("BEMA",
            f"impedance-EFIE done "
            f"({solve_info.get('method', solver)}, {t_lu:.1f}s): "
            f"R={R_coil*1e3:.4f} mOhm, L={L*1e9:.2f} nH")
        gf_J = GridFunction(fes_J)
        gf_J.vec.FV().NumPy()[:] = np.ascontiguousarray(J.real)
        gf_J_im = GridFunction(fes_J)
        gf_J_im.vec.FV().NumPy()[:] = np.ascontiguousarray(J.imag)
    else:
        # --- omega == 0: DC / vacuum-inductance solve (real saddle,
        #     R = 0).  The physical DC limit selected by frequency=0;
        #     the AC surface impedance does not exist here.  (Negative
        #     or NaN omega already rejected at the top of the
        #     function.) ---
        if solver in ("cocr", "hacapk_cocr"):
            mv = "hacapk" if solver == "hacapk_cocr" else "dense"
            coords = (_edge_midpoint_coords(mesh, n_J)
                      if mv == "hacapk" else None)
            _log("BEMA",
                f"DC vacuum-L {solver} solve (div-free reduction, "
                f"matvec={mv})")
            J, solve_info = _loop_cocr_solve(
                SL, None, D_red, g_red, 0.0, None, mv,
                coords=coords, log_fn=_log)
        else:
            _log("BEMA",
                f"DC vacuum-L saddle assembly: K is {_saddle_n}x{_saddle_n} "
                f"(~{(_saddle_n*_saddle_n*8)/1e9:.1f} GB), solver={solver}")
            K = np.block([
                [SL,              D_red.T],
                [D_red, _zero_c]
            ])
            rhs = np.zeros(n_J + n_constraint)
            rhs[n_J:] = g_red
            _log("BEMA",
                f"DC saddle solve start (method={solver}, ndof={_saddle_n})")
            x, solve_info = _solve_saddle(K, rhs, method=solver)
            J = x[:n_J]
        t_lu = time.perf_counter() - t0
        _log("BEMA",
            f"DC solve done ({solve_info.get('method', solver)}, "
            f"iters={solve_info.get('iterations','-')}, "
            f"residual={solve_info.get('residual',0.0):.2e}, {t_lu:.1f}s)")
        L = MU_0 * J @ SL @ J
        residual = np.max(np.abs(D @ J - g))
        R_coil = 0.0
        gf_J = GridFunction(fes_J)
        gf_J.vec.FV().NumPy()[:] = J
        gf_J_im = GridFunction(fes_J)
        gf_J_im.vec[:] = 0.0

    t_total = time.perf_counter() - t_start

    return {
        'L': float(L),
        'R': R_coil,
        'n_J': n_J,
        'n_f': n_f,
        'A_source': float(A_src),
        'A_sink': float(A_snk),
        't_assembly': round(t_assembly, 2),
        't_solve': round(t_lu, 2),
        't_total': round(t_total, 2),
        'residual': float(residual),
        'solver': solve_info,
        'J': J,
        'SL': SL,
        'D': D,
        'gf_J': gf_J,
        'gf_J_im': gf_J_im,
    }


def compute_centroids_areas_J(mesh, gf_J):
    """Sample the HDivSurface coil current at each BND triangle's centroid.

    Returns (centroids, areas, J_per_tri) -- (n_t, *)-shaped float64
    arrays.  ``gf_J`` is a REAL GridFunction; since the impedance-EFIE
    current is complex, callers reconstruct the complex per-triangle
    current by sampling BOTH GridFunctions returned by
    ``compute_inductance_source_sink``::

        cen, ar, J_re = compute_centroids_areas_J(mesh, res["gf_J"])
        _,   _,  J_im = compute_centroids_areas_J(mesh, res["gf_J_im"])
        J_per_tri = J_re + 1j * J_im

    (``radia.bem_sibc_solver.compute_phi_inc_from_surface_J`` is
    real-only, so downstream keeps bridging Re and Im separately.)

    Per-element averaging via NGSolve ``Integrate(..., element_wise=True)``
    over each component of the HDivSurface field.  This avoids the
    HDivSurface-specific reference-coord evaluation API and Just Works
    for any fes_order (RT₀ or higher).
    """
    from ngsolve import Integrate, BND, CF
    n_bnd = mesh.GetNE(BND)
    # Per-element area (denominator for the average).
    elem_A = Integrate(CF(1), mesh, VOL_or_BND=BND, element_wise=True)
    # Per-element integral of each J component.
    elem_J = [
        Integrate(gf_J[i], mesh, VOL_or_BND=BND, element_wise=True)
        for i in range(3)
    ]
    centroids = np.zeros((n_bnd, 3), dtype=np.float64)
    areas = np.zeros(n_bnd, dtype=np.float64)
    J_per_tri = np.zeros((n_bnd, 3), dtype=np.float64)
    for i, el in enumerate(mesh.Elements(BND)):
        a = abs(elem_A[el.nr])
        if a < 1e-30:
            continue
        pts = [mesh.vertices[v.nr].point for v in el.vertices]
        centroids[i] = np.mean(pts, axis=0)
        areas[i] = a
        J_per_tri[i] = [elem_J[k][el.nr] / a for k in range(3)]
    return centroids, areas, J_per_tri
