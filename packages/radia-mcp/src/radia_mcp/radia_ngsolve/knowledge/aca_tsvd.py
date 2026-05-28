"""(ACA+)+TSVD knowledge -- accelerated kernel-agnostic least-norm solver.

Read this when:
* Solving an underdetermined field-synthesis / inverse-source problem
  ``A phi = B`` (M field points x N basis sources, M < N) in the Radia stack.
* Designing a coil (stream function method) or reconstructing a magnetization.
* Wondering how to accelerate a TSVD pseudo-inverse whose matrix entries are
  Radia field evaluations (Biot-Savart / MMM / MSC).
* Pairing the linear solve with an outer CMA-ES (Optuna) design loop
  (the IEEJ SA-25-020 "ACA stream function + CMA-ES" workflow).

Production module: ``radia.stream_function`` (src/radia/stream_function.py),
C++ core src/core/rad_stream_function.cpp.  ACA+ is delegated to the in-repo
HACApK C library (cHACApK_acaplus).  Docs: docs/stream_function.md.  Examples:
examples/stream_function/.  Tests: tests/test_stream_function.py.

The MCP server exposes this via ``aca_tsvd(topic=...)``.  Topics: overview,
method, api, kernel_agnostic, performance, cmaes, validation, all.
"""


ACA_TSVD_OVERVIEW = r"""
# (ACA+)+TSVD least-norm solver

Accelerated, kernel-agnostic least-norm solver -- the numerical core of the
stream function method of coil design, generalised to ANY Radia source family.

## Problem

    A phi = B            A in R^{M x N},  usually M < N  (underdetermined)
    A(i,j) = (a field component) at observation i produced by basis source j
    B(i)   = desired field at observation i
    phi(j) = unknown strength of basis source j

Least-norm solution = truncated-SVD (TSVD) pseudo-inverse:
    A ~= U diag(S) V^T  ;  phi = V diag(1/S) U^T B   (truncated to k modes)

Truncation at k regularises the ill-posed inverse (small singular values
amplify noise).  Sweeping k traces the L-curve (residual vs solution norm).

## Why (ACA+)+TSVD

Forming the dense A and a full SVD is O(N M^2).  For smooth kernels A is
numerically LOW RANK, so:
  1. ACA+ factors A ~= C D^T with rank k_aca << min(M,N), evaluating only
     ~k_aca*(M+N) entries of A (never the full M*N).
  2. TSVD recompresses the small factors -> SVD of the rank-k_aca approximation.

Net ~ (M/k_aca)^2 cheaper than dense.  See topic "performance" for measured
numbers (~10x at N=2048, growing with N).

## Single source of truth for ACA+

ACA+ is NOT re-implemented: it is HACApK's cHACApK_acaplus (src/ext/HACApK),
fed an arbitrary matrix-entry callback via the HACApK_set_entry_func override
(default behaviour, the MMM/MSC system matrix, is unchanged when the override
is null).  Only the TSVD recompression (manuscript Method 2/3, SA-25-020) lives
in rad_stream_function.cpp.
"""


ACA_TSVD_METHOD = r"""
# Method: ACA+ then TSVD recompression

## Step 1 -- ACA+ (HACApK cHACApK_acaplus)

Adaptive Cross Approximation with pivoting builds A ~= C D^T,
C in R^{M x k_aca}, D in R^{N x k_aca}, evaluating O(k_aca(M+N)) entries.
Stops when the relative block norm < aca_eps.  Parameters used by
rad_stream_function.cpp to reproduce the validated reference exactly:
  param[61]=1  (absolute ACA_EPS, apxnorm = first block norm)
  param[64]=1  (minimum-rank guard)
  eps      = 1e-12 (relative convergence tolerance)
  pACA_EPS = aca_eps (user absolute pivot threshold)
Output zaa(M,kmax)=C and zab(N,kmax)=D are column-major.

## Step 2 -- TSVD recompression (two equivalent methods, IEEJ SA-25-020)

method=3 (DEFAULT, f90 method_aca_tsvd_2):
  SVD(C) -> Uc,Sc,VTc ; E = diag(Sc) VTc D^T ; SVD(E) -> UE,SE,VTE
  U = Uc UE ;  S = SE ;  V = VTE^T            (only TWO SVDs)

method=2 (f90 method_aca_tsvd_1):
  SVD(C), SVD(D), Middle = Sc (VTc VTd^T) Sd, SVD(Middle), combine.

Both give A ~= U diag(S) V^T truncated to `modes` (<= k_aca).  LAPACKE_dgesdd
(SVD) + cblas_dgemm (products), column-major internally; U/V converted to
row-major for NumPy.
"""


ACA_TSVD_API = r"""
# Python API (radia.stream_function)

    from radia.stream_function import (
        aca_tsvd, pseudo_inverse_solve, solve, radia_field_kernel, StreamTSVD,
    )

## aca_tsvd(M, N, entry, modes=None, kmax=None, aca_eps=1e-4, method=3) -> StreamTSVD
  entry(i, j) -> float : A(i,j), 0-based i in [0,M), j in [0,N).  Called on
  demand by ACA+ (O(k_aca(M+N)) calls, not M*N).
  Returns StreamTSVD: U (M,modes) row-major, S (modes,), V (N,modes) row-major,
  k_aca, method.  modes clamped to k_aca; kmax default min(M,N).

## pseudo_inverse_solve(result, B, k_mode=None) -> phi
  phi = V diag(1/S) U^T B using the first k_mode modes (default result.modes).
  For truncated-Tikhonov damping do it yourself from result.U/S/V:
    f = S/(S^2 + lam^2);  phi = V[:, :k] @ (f[:k] * (U[:, :k].T @ B))

## solve(M, N, entry, B, modes=None, k_mode=None, ...) -> (phi, result)
  Convenience: aca_tsvd then pseudo_inverse_solve.

## radia_field_kernel(obs_points, sources, component=2, field="b") -> entry
  Builds entry(i,j) = (component of field) at obs_points[i] from Radia object
  sources[j] via radia.Fld.  Reuses Radia's existing field for ANY source
  (coils, permanent magnets, soft iron).  component 0/1/2 = x/y/z.

## Example

    entry = radia_field_kernel(obs, loops, component=2)   # A(i,j)=Bz_i(loop_j)
    res   = aca_tsvd(len(obs), len(loops), entry, modes=20)
    phi   = pseudo_inverse_solve(res, B_target, k_mode=10) # loop currents
"""


ACA_TSVD_KERNEL_AGNOSTIC = r"""
# Kernel-agnostic design

The solver embeds NO field kernel.  A(i,j) is a caller callback, so the same
machinery serves every Radia source family using Radia's already-implemented
field computation:

  coils (thin wires)            -> Biot-Savart (ObjFlmCur, ObjArcCur)
  permanent magnets / soft iron -> MMM / MSC surface-charge field
                                   (ObjRecMag, ObjHexahedron, ...)

Use radia_field_kernel(obs, sources, ...) to build the callback from Radia
object handles -- there is no coil/magnet-specific code in the solver.

Why this matters (history, 2026-05-29):
- v1 hand-ported ACA+ from coil_solver.f90 -> rejected (2-fold maintenance of
  ACA+; HACApK already has the C implementation).  Fix: delegate to
  cHACApK_acaplus via HACApK_set_entry_func.
- v2 still embedded a coil-specific mirrored-rectangular-loop Biot-Savart kernel
  -> rejected (duplicates Radia's Biot-Savart; not reusable for magnets).
  Fix: kernel becomes a generic callback; coil Biot-Savart lives only in the
  test that cross-checks the f90 reference.

CLAUDE.md alignment: "Use HACApK only" (no custom H-matrix algorithms) and the
no-duplication principle.
"""


ACA_TSVD_PERFORMANCE = r"""
# Performance (measured)

examples/stream_function/bench_aca_vs_dense.py: (ACA+)+TSVD vs naive dense TSVD
(build full A via the SAME per-call kernel, then numpy.linalg.svd).  Smooth
1/(1+alpha r^2) kernel, M = N/4, LAB 2026-05-29.  The kernel is numerically low
rank so k_aca stays ~30 while N grows; the eval-count reduction
M*N -> ~k_aca(M+N) widens with N and, since both methods call the same kernel,
shows up almost one-for-one in wall-clock time:

  N     M    k_aca  kernel evals (naive->ACA)  eval cut  time (naive->ACA)  speedup
  256   64   31     16 384 -> 9 892            1.7x      37 ms -> 23 ms     1.6x
  512   128  33     65 536 -> 21 781           3.0x      141 ms -> 48 ms    2.9x
  1024  256  34     262 144 -> 45 588          5.8x      554 ms -> 96 ms    5.8x
  2048  512  36     1 048 576 -> 98 469        10.6x     2217 ms -> 204 ms  10.8x

Singular values match the dense SVD to ~2e-9.

Takeaways:
- Speedup GROWS with N (~ N/(5 k_aca) for M=N/4): dense is O(N M^2); ACA is
  O(k_aca(M+N)) evals + tiny factor SVDs.
- The win is largest when the kernel is EXPENSIVE (Biot-Savart / MMM-MSC /
  radia.Fld) -- every avoided A(i,j) is an avoided field evaluation.
- The matrix entry is a Python callback (per-call overhead).  Both methods pay
  it equally, so the RATIO is eval-count-driven; a C++ kernel lowers both
  absolute times by the same factor without changing the speedup.
- The flip side of low rank: ACA compresses min(M,N) to ~k_aca (e.g. 25 -> 2
  for a planar magnet array seen from afar).  Great for compression, but it
  means the array cannot support a high-DOF inverse from far sensors.
"""


ACA_TSVD_CMAES = r"""
# Optimisation layer: linear (TSVD) vs nonlinear (CMA-ES)

The (ACA+)+TSVD solve is the LINEAR design layer: source amplitudes phi with
fixed directions/positions.  When design variables enter the field NONLINEARLY
(magnetization directions/angles, magnet positions, coil-region geometry), it
is no longer a linear least-norm solve -> use a black-box optimiser.  This is
the "ACA stream function + CMA-ES" split (IEEJ SA-25-020): fast linear inner
solve = (ACA+)+TSVD; nonlinear outer search = CMA-ES.

Use Optuna's CmaEsSampler (do NOT re-implement CMA-ES):

    import optuna
    study = optuna.create_study(direction="minimize",
                                sampler=optuna.samplers.CmaEsSampler(seed=42))
    study.optimize(objective, n_trials=200)   # objective builds Radia magnets,
                                              # evaluates the field, returns scalar

CMA-ES is for continuous, mid-dimension (10-300) BBO; cast int/categorical with
care.  See the optuna_* MCP tools (optuna_algorithm topic="samplers",
optuna_recipes_advanced) for sampler choice, multi-objective (NSGA-II),
pruning, and lab BBO recipes.

Example: examples/stream_function/demo_cmaes_magnet_design.py optimises 16
magnetization angles for a uniform transverse field (16-D CMA-ES, ~3x objective
reduction).  Practical note (magnetics): the planar-array field operator is
either well-conditioned (close sensors -> regularisation idle) or numerically
rank ~1 (far sensors -> reconstruction hopeless), so the natural CMA-ES use is
nonlinear design (angles/geometry), while (ACA+)+TSVD owns the linear amplitude
solve.
"""


ACA_TSVD_VALIDATION = r"""
# Validation

tests/test_stream_function.py (10 tests):
- Reconstruction vs the TRUE dense A: machine precision when k_aca=min(M,N);
  tracks aca_eps otherwise.  Near-field full-rank + far-field low-rank cases.
- vs Fortran reference coil_solver.f90 (method_aca_tsvd_1/2 -- a faithful port
  of the same HACApK ACA+): identical k_aca and
  ||S_f90 - S_radia|| / ||S_f90|| ~ 1e-15 for BOTH methods.  Runs in a fresh
  subprocess to dodge conftest's DLL-search pollution (the f2py module bundles
  its own Intel/MKL DLLs).  LAB-only (skipif W: drive missing).
- Methods 2 and 3 agree to ~1e-9.
- Least-norm solve recovers B in range(A); validates B-length.
- Generic Radia-field path: factors a permanent-magnet array MMM/MSC coupling
  to < 1e-5 (test_radia_field_kernel_magnets).

f2py reference (LAB): W:\04_..\046_伊藤海人\2026_01_06_f2py_matlab比較\f2py.
"""


def get_aca_tsvd_knowledge(topic: str = "overview") -> str:
    """Return (ACA+)+TSVD knowledge for the requested topic."""
    t = (topic or "overview").strip().lower()
    table = {
        "overview": ACA_TSVD_OVERVIEW,
        "method": ACA_TSVD_METHOD,
        "api": ACA_TSVD_API,
        "kernel_agnostic": ACA_TSVD_KERNEL_AGNOSTIC,
        "performance": ACA_TSVD_PERFORMANCE,
        "cmaes": ACA_TSVD_CMAES,
        "validation": ACA_TSVD_VALIDATION,
    }
    aliases = {
        "kernel": "kernel_agnostic", "generic": "kernel_agnostic",
        "speed": "performance", "speedup": "performance", "benchmark": "performance",
        "optuna": "cmaes", "cma-es": "cmaes", "cma_es": "cmaes",
        "stream_function": "overview", "stream": "overview", "tsvd": "method",
        "aca": "method", "validate": "validation", "f90": "validation",
    }
    t = aliases.get(t, t)
    if t == "all":
        return "\n\n".join([ACA_TSVD_OVERVIEW, ACA_TSVD_METHOD, ACA_TSVD_API,
                            ACA_TSVD_KERNEL_AGNOSTIC, ACA_TSVD_PERFORMANCE,
                            ACA_TSVD_CMAES, ACA_TSVD_VALIDATION])
    if t in table:
        return table[t]
    return (f"Unknown topic '{topic}'.  Available: overview, method, api, "
            "kernel_agnostic, performance, cmaes, validation, all.\n\n"
            + ACA_TSVD_OVERVIEW)
