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
    "aca_tsvd",
    "pseudo_inverse_solve",
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
    method : int
        Recompression method used (2 or 3).
    """

    U: np.ndarray
    S: np.ndarray
    V: np.ndarray
    k_aca: int
    method: int

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
             aca_eps=1.0e-4, method=3) -> StreamTSVD:
    """(ACA+)+TSVD recompressed truncated SVD of an M x N matrix A.

    Parameters
    ----------
    M : int
        Number of rows (field / observation points).
    N : int
        Number of columns (basis sources).  Typically N > M (underdetermined).
    entry : callable
        ``entry(i, j) -> float`` returning A(i,j) for 0-based i in [0,M),
        j in [0,N).  Supplied by the caller from Radia's field computation
        (see ``radia_field_kernel``).  Called on demand by ACA+
        (O(k_aca * (M + N)) evaluations, not the full M*N).
    modes : int, optional
        Number of singular triplets to return.  Clamped to the ACA+ rank
        k_aca.  Defaults to ``kmax``.
    kmax : int, optional
        Maximum ACA+ rank.  Defaults to ``min(M, N)``.
    aca_eps : float, optional
        ACA+ stopping tolerance (absolute pivot/row/col threshold).  Default 1e-4.
    method : int, optional
        3 (default) = improved 2-SVD recompression (f90 method_aca_tsvd_2);
        2 = full re-SVD of both factors (f90 method_aca_tsvd_1).

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
    if method not in (2, 3):
        raise ValueError(f"method must be 2 or 3, got {method}")

    if kmax is None:
        kmax = min(M, N)
    kmax = int(min(kmax, M, N))
    if modes is None:
        modes = kmax
    modes = int(min(modes, kmax))

    # Wrap so the C++ side always receives plain Python floats.
    def _entry(i, j):
        return float(entry(i, j))

    U, S, V, k_aca = _cpp_aca_tsvd(
        M, N, _entry, int(modes), int(kmax), float(aca_eps), int(method))
    return StreamTSVD(U=U, S=S, V=V, k_aca=int(k_aca), method=int(method))


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


def solve(M, N, entry, B, modes=None, k_mode=None,
          kmax=None, aca_eps=1.0e-4, method=3):
    """Convenience: (ACA+)+TSVD decompose then pseudo-inverse solve.

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
