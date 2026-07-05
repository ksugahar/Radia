"""(ACA+)+TSVD accelerated least-norm solver -- generic and kernel-agnostic.

Solves the underdetermined least-norm system  A phi = B  (M field points x
N basis sources, M < N) via a TSVD-regularized pseudo-inverse.  A ~= C D^T is
factored with ACA+ (HACApK's cHACApK_acaplus), then only the small factors are
TSVD'd -> ~ (M/k)^2 faster than naive O(N M^2) TSVD.

KERNEL-AGNOSTIC.  The matrix entry A(i,j) is supplied by the caller as a
callable ``entry(i, j) -> float`` (0-based).  The SAME machinery therefore
serves any Radia source family using Radia's *already-implemented* field
computation -- there is no field kernel baked into this module:

    - coils                : Biot-Savart H/A from filaments
    - permanent magnets,   : MMM / MSC field from magnetization
      soft iron

For convenience, ``radia_field_kernel`` builds such a callback directly from
Radia object handles via ``radia.Fld`` (works for coils AND magnetic materials).

C++ core: src/core/rad_stream_function.cpp (ACA+ delegated to HACApK; only the
TSVD recompression, manuscript Method 2/3 of IEEJ SA-25-020, lives there).

Example
-------
    import numpy as np, radia as rad
    from radia.stream_function import aca_tsvd, pseudo_inverse_solve, radia_field_kernel

    obs = np.random.rand(M, 3)          # observation points
    sources = [rad.ObjRecMag(...), ...] # N Radia objects (coils or magnets)

    entry = radia_field_kernel(obs, sources, component=2)   # A(i,j) = Bz(obs_i, src_j)
    res = aca_tsvd(M, N, entry, modes=40)
    phi = pseudo_inverse_solve(res, B_target, k_mode=30)
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from radia._radia_pybind import _stream_aca_tsvd as _cpp_aca_tsvd

__all__ = [
    "StreamTSVD",
    "RegularizedTSVD",
    "aca_tsvd",
    "pseudo_inverse_solve",
    "pseudo_inverse_solve_regularized",
    "solve",
    "radia_field_kernel",
]


@dataclass
class StreamTSVD:
    """Recompressed truncated SVD of A:  A ~= U diag(S) V^T (truncated to modes).

    Attributes
    ----------
    U : ndarray, shape (M, modes), row-major (C-contiguous)
    S : ndarray, shape (modes,)
    V : ndarray, shape (N, modes), row-major (C-contiguous)
    k_aca : int
        ACA+ rank found before TSVD truncation.
    method : str
        Which method produced this SVD: ``"aca_qr_tsvd"`` (ACA + QR + TSVD, the
        fast default) or ``"dense"`` (direct dense TSVD, the exact reference).
        (The legacy manuscript Method 2/3 were removed 2026-07-05, JIAM-2026-36.)
    """

    U: np.ndarray
    S: np.ndarray
    V: np.ndarray
    k_aca: int
    method: str

    @property
    def M(self) -> int:
        return int(self.U.shape[0])

    @property
    def N(self) -> int:
        return int(self.V.shape[0])

    @property
    def modes(self) -> int:
        return int(self.S.shape[0])


def aca_tsvd(M, N, entry, modes=None, kmax=None,
             aca_eps=1.0e-4, method="aca_qr_tsvd") -> StreamTSVD:
    """Truncated SVD of an M x N matrix A supplied by the callback ``entry(i,j)``.

    Two methods -- and ONLY two (peer review JIAM-2026-36):

    - ``method="aca_qr_tsvd"`` (default; short aliases ``"qr"`` / ``"aca"``):
      ACA+ factors ``A ~= C D^T`` (only O(k_aca*(M+N)) entry evaluations), then
      the standard "SVD of a low-rank product" recompression -- QR each
      tall-skinny factor + ONE small ``k_aca x k_aca`` TSVD.  Fast; the
      production path.  (Name = the full pipeline: ACA + QR + TSVD.)
    - ``method="dense"`` (alias ``"tsvd"``): the plain/direct TSVD -- materialise
      the full A via the callback and take its dense SVD (``numpy.linalg.svd``),
      then truncate.  EXACT (no ACA approximation), the trusted reference /
      validation baseline, but O(M*N) entry calls + an O(N*M^2) SVD -- use for
      SMALL problems only.

    Parameters
    ----------
    M, N : int
        Rows (field / observation points) and columns (basis sources); usually N>M.
    entry : callable
        ``entry(i, j) -> float`` returning A(i,j) for 0-based i in [0,M), j in [0,N).
    modes : int, optional
        Singular triplets to return (clamped to the rank).  Defaults to ``kmax``.
    kmax : int, optional
        Maximum ACA+ rank (``aca_qr_tsvd`` only).  Defaults to ``min(M, N)``.
    aca_eps : float, optional
        ACA+ stopping tolerance (``aca_qr_tsvd`` only).  Default 1e-4.
    method : {"aca_qr_tsvd", "dense"}, optional
        "aca_qr_tsvd" (default; aliases "qr"/"aca") = ACA + QR + TSVD;
        "dense"/"tsvd" = direct dense SVD (exact reference).  Legacy integers
        2/3 and None map to "aca_qr_tsvd" (the manuscript Method 2/3 were removed
        2026-07-05; see ``memory/aca_tsvd_qr_recompression.md``).

    Returns
    -------
    StreamTSVD
    """
    M = int(M)
    N = int(N)
    if M <= 0 or N <= 0:
        raise ValueError(f"M and N must be positive, got M={M}, N={N}")
    if not callable(entry):
        raise TypeError("entry must be callable: entry(i, j) -> float")

    # Only two methods are supported (No-Fallback: an unknown value RAISES).
    # Legacy int 2/3 and None -> "aca_qr_tsvd" (the manuscript Method 2/3 were removed).
    _m = "aca_qr_tsvd" if method in (None, 2, 3) else str(method).strip().lower()
    if _m in ("aca_qr_tsvd", "aca+qr+tsvd", "acaqrtsvd", "qr", "aca"):
        _m = "aca_qr_tsvd"
    elif _m in ("dense", "tsvd", "direct", "svd", "full"):
        _m = "dense"
    else:
        raise ValueError(
            f"method must be 'aca_qr_tsvd' (default; ACA+QR+TSVD) or 'dense' "
            f"(direct TSVD); got {method!r}")

    if kmax is None:
        kmax = min(M, N)
    kmax = int(min(kmax, M, N))
    if modes is None:
        modes = kmax
    modes = int(min(modes, kmax))

    if _m == "dense":
        # 通常のTSVD: materialise A via the callback, take its dense SVD, truncate.
        A = np.fromiter((float(entry(i, j)) for i in range(M) for j in range(N)),
                        dtype=float, count=M * N).reshape(M, N)
        Uf, Sf, Vtf = np.linalg.svd(A, full_matrices=False)
        m = int(min(modes, Sf.shape[0]))
        return StreamTSVD(U=np.ascontiguousarray(Uf[:, :m]),
                          S=np.ascontiguousarray(Sf[:m]),
                          V=np.ascontiguousarray(Vtf[:m, :].T),
                          k_aca=int(min(M, N)), method="dense")

    # method="aca_qr_tsvd": ACA+ then the QR-of-a-low-rank-product recompression (C++).
    def _entry(i, j):                       # C++ side always gets plain floats
        return float(entry(i, j))
    U, S, V, k_aca = _cpp_aca_tsvd(
        M, N, _entry, int(modes), int(kmax), float(aca_eps))
    return StreamTSVD(U=U, S=S, V=V, k_aca=int(k_aca), method="aca_qr_tsvd")


def pseudo_inverse_solve(result: StreamTSVD, B, k_mode=None) -> np.ndarray:
    """Least-norm pseudo-inverse solve:  phi = V diag(1/S) U^T B.

    Uses the first ``k_mode`` singular triplets (TSVD regularization).

    Parameters
    ----------
    result : StreamTSVD
    B : array_like, shape (M,)
        Target field values.
    k_mode : int, optional
        Number of modes to use (<= result.modes).  Default = result.modes.

    Returns
    -------
    ndarray, shape (N,)
        Basis coefficients phi.
    """
    B = np.asarray(B, dtype=float).ravel()
    if B.shape[0] != result.M:
        raise ValueError(f"B length {B.shape[0]} != M {result.M}")
    k = result.modes if k_mode is None else int(min(k_mode, result.modes))
    k = max(1, k)
    U = result.U[:, :k]   # (M, k)
    S = result.S[:k]      # (k,)
    V = result.V[:, :k]   # (N, k)
    with np.errstate(divide="ignore", invalid="ignore"):
        c = (U.T @ B) / S
    c = np.where(S > 0.0, c, 0.0)
    return V @ c


# --------------------------------------------------------------------------
# Regularisation-folded ACA+TSVD pseudo-inverse
# --------------------------------------------------------------------------
@dataclass
class RegularizedTSVD:
    """Folds an SPD stiffness ``S`` into the ACA+TSVD pseudo-inverse.

    Solves the constrained min-norm problem

        min  phi^T S phi    s.t.    A phi = B

    with the closed-form Lagrangian solution

        phi = S^-1 A^T (A S^-1 A^T)^-1 B.

    Substituting the truncated SVD ``A = U Sigma V^T`` (orthonormal V
    columns from ACA+TSVD) gives the FACTORED form

        phi = (S^-1 V) . W^-1 . Sigma^-1 . U^T B,    W = V^T S^-1 V (k x k).

    The factorisation ``(S^-1 V, W^-1)`` is precomputed once; subsequent
    ``solve(B)`` calls cost only ``O(k * (M + N))`` floating-point ops
    (no further sparse / dense solves).  This makes the Path-A
    compensated iteration's inner solve cheap regardless of the FE
    DOF count.

    ``solve(B, alpha)`` with ``alpha > 0`` is the **generalised Tikhonov**
    (soft data fit) ``min ||A phi - B||^2 + alpha . phi^T S phi`` via the
    SAME cached factorisation -- only an ``alpha I`` is added inside the
    ``k x k`` core ``(alpha I + Sigma^2 W)``.  So an entire L-curve /
    Pareto alpha-sweep re-uses one factorisation (the ``+ alpha I`` is the
    whole cost of moving along the front); ``alpha = 0`` recovers the
    exact-fit form above.

    Special case ``S = I``: ``W = V^T V = I_k`` since ACA+TSVD V columns
    are orthonormal, so ``phi = V Sigma^-1 U^T B`` -- reduces exactly to
    the standard L2 pseudo-inverse (see ``pseudo_inverse_solve``).

    Parameters
    ----------
    base : StreamTSVD
        ACA+TSVD factorisation of A (M x N, M <= k).
    Sinv_V : ndarray, shape (N, k)
        Precomputed ``S^-1 V`` (one k-RHS solve with S).
    W_inv : ndarray, shape (k, k)
        Inverse of ``W = V^T Sinv_V``.

    Notes
    -----
    For an N x N stiffness ``S`` restricted to free DOFs, the caller
    must pass ``S_free`` and the matching ``V`` (= base.V on free
    DOFs).  Padding back to full DOF (e.g. Dirichlet zero on boundary)
    is the caller's responsibility.

    See ``docs/stream_function/theory.ipynb`` and the docs-local
    ``demo_planar_uniform_fem_psi.py`` helper for the canonical Path-A usage
    pattern.
    """

    base: StreamTSVD
    Sinv_V: np.ndarray         # (N, k)
    W_inv: np.ndarray          # (k, k)
    W: np.ndarray = None       # (k, k) = V^T S^-1 V; needed for alpha > 0

    @classmethod
    def from_stiffness(cls, base: "StreamTSVD", S) -> "RegularizedTSVD":
        """Precompute ``Sinv_V`` and ``W_inv`` from an SPD stiffness ``S``.

        Parameters
        ----------
        base : StreamTSVD
            ACA+TSVD factorisation of the system matrix A.
        S : ndarray or scipy.sparse matrix, shape (N, N)
            SPD regularisation matrix (e.g. H1 stiffness on free DOFs).
            Dense numpy arrays go through ``np.linalg.solve``; scipy sparse
            matrices are factored ONCE with ``scipy.sparse.linalg.splu`` and
            all columns of ``V`` are back-solved together.

        Raises
        ------
        ValueError
            If ``S`` is singular / not SPD (No-Fallback: the dense
            ``LinAlgError`` and the sparse ``RuntimeError`` are normalised into
            one clear error), or if the inverted core ``W = V^T S^-1 V`` is
            catastrophically ill-conditioned (``cond > 1e12``, where the
            exact-fit solve would silently lose most of its digits).  Add a
            ridge / mass term to ``S`` or use a better-conditioned seminorm.

        Returns
        -------
        RegularizedTSVD
            Cached factorisation; pass ``B`` to ``.solve(B)``.
        """
        V = base.V                                 # (N, k)
        try:
            if hasattr(S, "tocsc"):
                # scipy sparse path -- factor S ONCE (sparse LU) then back-solve
                # all k columns of V together.  (Was k separate spsolve calls,
                # each re-factorising S; this is O(1) factorisation + k solves --
                # the win that lets N=15k-DOF designs build in ~1 s, not ~13 s.)
                from scipy.sparse.linalg import splu
                Sinv_V = splu(S.tocsc()).solve(np.asarray(V, dtype=float))
            else:
                S_dense = np.asarray(S, dtype=float)
                if S_dense.shape != (V.shape[0], V.shape[0]):
                    raise ValueError(
                        f"S shape {S_dense.shape} must be ({V.shape[0]}, "
                        f"{V.shape[0]}) to match base.V")
                Sinv_V = np.linalg.solve(S_dense, V)
        except (np.linalg.LinAlgError, RuntimeError) as e:
            # No-Fallback: the dense path raises numpy.linalg.LinAlgError and the
            # sparse splu raises RuntimeError on a singular S -- normalise BOTH
            # into ONE clear, actionable error rather than leaking two types.
            raise ValueError(
                f"regularisation matrix S is singular / not SPD "
                f"({type(e).__name__}: {e}). Add a ridge / mass term to S (the "
                f"_seminorm helper adds +1e-10 mass for exactly this), or use "
                f"--confine abe/on to remove the constant-psi null space.") from e
        W = V.T @ Sinv_V                            # (k, k) -- the inverted core
        # No-Fallback: a catastrophically ill-conditioned core makes the exact-
        # fit solve SILENTLY lose most of its digits.  Production paths sit at
        # cond(W) ~ 1e2 (l2/abe) .. 1e7 (h1/off); only a genuinely broken S
        # exceeds 1e12 (< ~4 good digits) -- raise rather than return garbage.
        cond_W = float(np.linalg.cond(W))
        if not np.isfinite(cond_W) or cond_W > 1e12:
            raise ValueError(
                f"regularised core W = V^T S^-1 V is too ill-conditioned "
                f"(cond ~ {cond_W:.1e}); the exact-fit solve would lose most "
                f"digits. Add a ridge / mass term to S, use --confine abe, or a "
                f"better-conditioned --regularize.")
        W_inv = np.linalg.inv(W)               # back-compat; W is sound past here
        return cls(base=base, Sinv_V=Sinv_V, W_inv=W_inv, W=W)

    def solve(self, B, k_mode=None, alpha=0.0) -> np.ndarray:
        """Apply the cached regularised pseudo-inverse to ``B``.

        ``alpha = 0`` (default) -- exact data fit, minimum ``S`` seminorm::

            phi = Sinv_V . W_inv . diag(1/Sigma) . U^T B

        ``alpha > 0`` -- generalised **Tikhonov** (soft data fit), trading
        misfit against the seminorm,
        ``min ||A phi - B||^2 + alpha . phi^T S phi``::

            phi(alpha) = Sinv_V . (alpha I + Sigma^2 W)^-1 . Sigma . U^T B
                       == (A^T A + alpha S)^-1 A^T B

        This is the SAME cached factorisation with a single ``alpha I``
        added inside the ``k x k`` core, so an L-curve / Pareto alpha-sweep
        re-solves only the small core (no re-factorisation).  Special case
        ``S = I``: the core reduces to the classic Tikhonov filter factors
        ``sigma / (sigma^2 + alpha)``.  See
        ``docs/stream_function/regularization.md`` and its docs-local
        ``demo_pareto_tikhonov_aca.py`` helper.

        Parameters
        ----------
        B : array_like, shape (M,)
            Target field values.
        k_mode : int, optional
            Truncate to first ``k_mode`` SVD modes (re-inverts the k_mode
            x k_mode top-left block of W since W_inv depends on k).
            Default = base.modes (use full cached factorisation).
        alpha : float, optional
            Tikhonov weight.  ``0`` (default) = exact-fit min-seminorm;
            ``> 0`` = soft-fit Tikhonov in the ``S`` metric.  Composes
            with ``k_mode`` (hard spectral truncation + smooth damping).

        Returns
        -------
        ndarray, shape (N,)
            ``alpha = 0``: ``phi`` satisfying ``A phi = B`` with minimum
            ``phi^T S phi``.  ``alpha > 0``: the Tikhonov solution
            ``(A^T A + alpha S)^-1 A^T B``.
        """
        B = np.asarray(B, dtype=float).ravel()
        if B.shape[0] != self.base.M:
            raise ValueError(f"B length {B.shape[0]} != M {self.base.M}")
        k_full = self.base.modes
        k = k_full if k_mode is None else int(min(k_mode, k_full))
        k = max(1, k)

        U = self.base.U[:, :k]                      # (M, k)
        Sigma = self.base.S[:k]                     # (k,)

        if k == k_full:
            Sinv_V = self.Sinv_V                    # (N, k_full)
            W_inv = self.W_inv                      # (k_full, k_full)
            W = self.W if self.W is not None else np.linalg.inv(W_inv)
        else:
            Sinv_V = self.Sinv_V[:, :k]             # (N, k)
            V_k = self.base.V[:, :k]                # (N, k)
            W = V_k.T @ Sinv_V                      # (k, k)

        UtB = U.T @ B                               # (k,)

        if alpha <= 0.0:
            # exact-fit min-seminorm: phi = Sinv_V . W^-1 . Sigma^-1 . U^T B.
            # Solve W y = c (factored) instead of caching/applying inv(W): a
            # factored solve is strictly more stable on an ill-conditioned W
            # (No-Fallback: do not silently lose digits to an explicit inverse).
            with np.errstate(divide="ignore", invalid="ignore"):
                c = UtB / Sigma                     # (k,)
            c = np.where(Sigma > 0.0, c, 0.0)
            return Sinv_V @ np.linalg.solve(W, c)

        # generalised Tikhonov: (alpha I + Sigma^2 W) y = Sigma U^T B
        #   phi = Sinv_V . y == (A^T A + alpha S)^-1 A^T B
        core = alpha * np.eye(k) + (Sigma ** 2)[:, None] * W
        y = np.linalg.solve(core, Sigma * UtB)
        return Sinv_V @ y


def pseudo_inverse_solve_regularized(result: StreamTSVD, B, S, k_mode=None):
    """One-shot regularised pseudo-inverse solve.

    Builds ``RegularizedTSVD.from_stiffness(result, S)`` and applies
    it once.  Use ``RegularizedTSVD`` directly when the same ``S`` is
    reused across multiple ``B`` (= Path-A compensated iteration).

    Parameters
    ----------
    result : StreamTSVD
        ACA+TSVD factorisation of A.
    B : array_like, shape (M,)
        Target field values.
    S : ndarray or scipy.sparse matrix, shape (N, N)
        SPD regularisation matrix (e.g. H1 stiffness on free DOFs).
    k_mode : int, optional
        Truncate to first ``k_mode`` SVD modes (default = result.modes).

    Returns
    -------
    ndarray, shape (N,)
        ``phi`` satisfying ``A phi = B`` with minimum ``phi^T S phi``.
    """
    reg = RegularizedTSVD.from_stiffness(result, S)
    return reg.solve(B, k_mode=k_mode)


def solve(M, N, entry, B, modes=None, k_mode=None,
          kmax=None, aca_eps=1.0e-4, method="aca_qr_tsvd"):
    """Convenience: (ACA+)+TSVD decompose then pseudo-inverse solve.
    ``method`` selects "aca_qr_tsvd" (ACA+QR+TSVD, default) or "dense" (direct
    TSVD); see :func:`aca_tsvd`.

    Returns
    -------
    (phi, result) : (ndarray shape (N,), StreamTSVD)
    """
    result = aca_tsvd(M, N, entry, modes=modes, kmax=kmax,
                      aca_eps=aca_eps, method=method)
    phi = pseudo_inverse_solve(result, B, k_mode=k_mode)
    return phi, result


def radia_field_kernel(obs_points, sources, component=2, field="b"):
    """Build a matrix-entry callback from Radia's existing field computation.

    ``A(i, j)`` = (``component`` of the ``field``) at observation point i
    produced by Radia source object ``sources[j]``, evaluated via ``radia.Fld``.
    This is the generic reuse of Radia's already-implemented kernels and works
    for ANY source family -- coils (Biot-Savart), permanent magnets / soft iron
    (MMM / MSC) -- without embedding a new field formula.

    Parameters
    ----------
    obs_points : array_like, shape (M, 3)
        Observation points.
    sources : sequence of int
        Length-N sequence of Radia object handles.
    component : int, optional
        Field component index 0/1/2 for x/y/z (default 2 = z).
    field : str, optional
        Radia field id passed to ``radia.Fld`` (default "b"; "h", "a", ...).

    Returns
    -------
    callable
        ``entry(i, j) -> float`` for use with ``aca_tsvd`` / ``solve``
        (M = len(obs_points), N = len(sources)).
    """
    import radia as rad

    obs = np.ascontiguousarray(obs_points, dtype=float).reshape(-1, 3)
    src = [int(s) for s in sources]
    comp = int(component)

    def entry(i, j):
        B = rad.Fld(src[j], field, obs[i].tolist())
        return float(np.asarray(B, dtype=float).ravel()[comp])

    return entry
