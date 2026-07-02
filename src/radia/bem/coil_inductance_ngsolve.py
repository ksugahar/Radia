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
        solver: saddle-point solver — "lu" (dense direct, default) or
            "gmres".  "minres" is accepted only for the ``omega == 0``
            real solve: scipy's MINRES silently discards the imaginary
            part of a complex system (solves the real part only), so it
            is REJECTED for the AC impedance-EFIE (fail fast).
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
            f"impedance-EFIE done ({solver}, {t_lu:.1f}s): "
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
            f"DC saddle solve done ({solver}, "
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
