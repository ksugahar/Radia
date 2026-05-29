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
method, api, kernel_agnostic, performance, cmaes, validation, literature
(stream function method theory + Turner/Peeren/Abe-DUCAS lineage), workflow
(SF -> CAD/STEP -> PEEC -> field demos), single_stroke (one-continuous-wire
connection: Kuijpers/Lomonova prior art + our crossover-compensation research),
all.
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


ACA_TSVD_LITERATURE = r"""
# Stream function method (SFM) -- theory and literature lineage

Our solver is the numerical core of the STREAM FUNCTION METHOD (a.k.a. CURRENT
POTENTIAL method) of coil / surface-current design.  The lineage below is the
citation backbone for any paper built on radia.stream_function.

## The method in one paragraph

A divergence-free surface current K on a meshed source surface S is written from
a single scalar STREAM FUNCTION (current potential) psi:
    K = grad(psi) x n            (n = surface normal)
Equal-spaced iso-contours of psi ARE the wires: between levels psi=(n-1)*dI and
psi=n*dI flows current dI, so placing one conductor at psi=(n-0.5)*dI per band
reproduces the continuous K.  The field is linear in psi:  B = A psi  (A from
Biot-Savart on each FE / node basis).  Designing a coil = solving the
underdetermined inverse for psi given a target B over a region (FOV / DSV).
Regularise by truncated SVD (TSVD).  radia.stream_function = an accelerated,
kernel-agnostic re-implementation of exactly this CP+TSVD solve.

## Lineage (chronological, with roles)

- TARGET FIELD METHOD -- R. Turner, J. Phys. D 19 (1986) L147 (+ review Magn.
  Reson. Imaging 11 (1993) 903).  Cylinder, Fourier-Bessel inversion of B_target
  -> surface current j(phi,z); continuity div j = 0; convergence needs apodising
  (Gaussian) the target.  Discretisation = "integrated-current contour lines at
  equal-current intervals" -- the contour==wire principle we use.  m=0 target ->
  axisymmetric (Gz, azimuthal rings); m=+-1 -> transverse Golay double-saddle
  (Gx).  This is the origin of cylindrical gradient/shim coil design.

- CURRENT POTENTIAL IN FUSION -- A. Kameari, J. Comput. Phys. 42 (1981) 124
  (thin-conductor current potential; CP==SF identity).  P. Merkel, Nucl. Fusion
  27 (1987) 867 (NESCOIL: stellarator modular coils from a CP on a winding
  surface).  A.H. Boozer, Phys. Plasmas 7 (2000) 629.  Same math as MRI SFM but
  on a toroidal winding surface -- the "arbitrary surface" generalisation.

- STREAM FUNCTION MODELLED DIRECTLY -- G.N. Peeren, J. Comput. Phys. 191 (2003)
  305 (TU/e / Philips).  Argues to optimise psi DIRECTLY (not current density
  then reconstruct psi); casts shape/field-synthesis as quadratic programming
  with linear constraints.  Foundational TU/e line that leads to Lomonova /
  Kuijpers single-stroke work (see topic "single_stroke").

- SFM FOR GRADIENT COILS -- A.L. Lemdiasov & R. Ludwig, Concepts Magn. Reson. B
  26 (2005) 67.  Poole / Crozier / Lopez et al. (equivalent magnetization
  current; the "fingerprint" transverse patterns), e.g. IEEE TM 45 (2009) 767;
  minimax current density, J. Phys. D 43 (2010) 095001.  Z. Liu, J. Hennig,
  J.G. Korvink, IEEE TM 48 (2012) 1179 (discretised SF with HIGH-ORDER
  SMOOTHNESS -- directly relevant to smooth single-stroke wiring).

- CP + TSVD = OUR EXACT METHOD (DUCAS) -- M. Abe et al.:
  * Phys. Plasmas 10 (2003) 1022 -- DUCAS ("Design tool Using Current potential
    And SVD") for stellarator modular coils.
  * IEEE TM 49 (2013) 5645 -- DUCAS for MRI gradient coils (weighted nodes,
    initial CP).
  * IEEE TM 50 (2014) 5100911 -- active-shield gradient coil (two current
    surfaces solved iteratively).
  DUCAS formulation (== radia.stream_function):
    B = A T        T = node CP (== psi) vector,  A = Biot-Savart response matrix
    j = grad(T) x n                              (current between nodes i,j = Ti-Tj)
    TSVD pseudo-inverse truncated to MD eigenmodes (MD == our `modes`/k):
        (W_B A R W_I^-1)* = sum_i v_i u_i^T / lambda_i  (i=1..MD)
    node weights W_I suppress peaked currents; eigenmode strength
    D_i = u_i^T W_B B_TG; pick MD so residual peak-to-peak < eps.
  Mapping to us: Abe MD <-> our truncation `k_mode`; Abe A response matrix <->
  our radia_field_kernel; Abe node weights <-> optional W in the solve.  We ADD
  ACA+ acceleration (Abe forms A densely) and make A a generic callback so the
  same code drives coils AND magnets/iron.  The lab f90 coil_solver.f90
  (method_aca_tsvd_1/2) we validate against is a DUCAS-lineage code.

- NONLINEAR OUTER LOOP -- D. Tomasi, Magn. Reson. Med. 45 (2001) 505 (SF +
  simulated annealing to optimise short self-shielded gradient coils).  The
  modern analogue is our CMA-ES outer loop (topic "cmaes"): SA/CMA-ES handles the
  nonlinear DOF, TSVD owns the linear amplitude solve.

## Source folder (LAB, owner-password PDFs; decrypt with pikepdf)
W:\03_文献..\..\流れ関数法 (Turner target-field subfolder; Truncated Singular
Value Decomposition subfolder = the Abe / NESCOIL / TSVD line).
"""


ACA_TSVD_WORKFLOW = r"""
# End-to-end workflow: SF design -> CAD(STEP) -> PEEC -> field

The stream function solve is only the first stage.  examples/stream_function/
carries it all the way to a manufacturable, field-verified conductor:

## 1D axisymmetric -- demo_coil_design_gz.py (cylindrical Gz gradient)
Target Bz = Gz*z is axisymmetric -> surface current is purely azimuthal -> psi
reduces to psi(z) and the basis is FULL azimuthal RINGS at z_1..z_N (the correct
reduced model; also keeps every field eval on Radia's reliable full-ring path).
On-axis + OFF-AXIS volume target (rho in {0, 0.4a, 0.7a} at a single azimuth --
safe because full rings are axisymmetric) -> (ACA+)+TSVD ring currents I(z) ->
psi(z) = cumulative integral -> equal-current contour -> wire rings (generalised
Maxwell pair).  Verifies on-axis dBz/dz + DSV volume nonuniformity.

## 2D transverse -- demo_coil_design_gx.py (cylindrical Gx gradient)
Target Bz = Gx*x is NOT axisymmetric -> genuine 2D surface psi(phi,z) fingerprint.
Per-node quad loops on the cylinder; A(i,j)=Bz via an EXACT numpy Biot-Savart
segment kernel (NOT rad.ObjFlmCur -- it has a ~10x-wrong-Bz bug for small TILTED
loops at phi=90 deg; full rings are unaffected).  marching-squares contour ->
~68 saddle loops; reconstructed Bz matches Gx*x to ~0.8% RMS.  k_aca reaches
min(M,N) here (the transverse target is full-rank, unlike the low-rank Gz case).
Single-stroke connection of the nested fingerprint loops = future work
(topic "single_stroke").

## Full chain -- demo_sf_to_peec_gz.py (--with-peec)
SF design -> SINGLE-STROKE smooth helix (coaxial equal-current rings joined into
one continuous wire; cosine-blended axial-ramp crossovers; handedness flips at
current-sign changes) -> CAD STEP (build123d Spline centerline + Frenet swept
solid; wire radius auto < half the min turn spacing so the tube cannot
self-intersect) -> PEEC (peec_matrices.PEECBuilder -> PEECCircuitSolver ->
port impedance; L = Im(Z)/(2 pi f), R = Re(Z)) -> exact Biot-Savart field ->
verify dBz/dz vs the design Gz*z.  ~16-turn run: ~15 m single wire,
dBz/dz ~ 0.99, ~2.6% nonlinearity, L ~ 38 uH, R ~ 16 mOhm @ 1 kHz.  Confirms the
SF design survives the single-stroke manufacturing constraint.

## Helper APIs used
radia.biot_savart.h_segments_batch(segments, obs, current) -- exact segment
Biot-Savart (B = MU0 * H), the trustworthy kernel.  radia.coil_from_cad (STEP
centerline).  peec_matrices.PEECBuilder / peec_topology.PEECCircuitSolver.
CoilBuilder is for PLANAR racetrack/saddle coils; a solenoidal helix uses the
smooth-helix + Spline path, not CoilBuilder arcs.

## Stage-2 panel CLI
src/radia/panels/calc_stream_coil.py wraps the Gz design as a headless Layer-4
script (argparse in, JSON out, no Cubit/PySide6); locked by
tests/panels/test_stream_coil_golden.py (fitted_dBdz in [0.9,1.1],
gradient_nonlinearity < 0.05).  Stage-3 PySide6 panel = future work.
"""


ACA_TSVD_SINGLE_STROKE = r"""
# Single-stroke (one continuous wire) coil design

## The problem
SFM iso-contours are SEPARATE closed loops (Gz: rings; Gx: nested fingerprint
loops).  A real wound coil must be ONE continuous conductor driven by one current
source.  Connecting the loops adds CROSSOVER / connection segments that carry the
full current and produce a PARASITIC field that was not part of the design ->
degrades the target-region field.  How to connect with minimal damage is the
single-stroke problem.

## Prior art -- Kuijpers, Jansen, Lomonova (TU/e EPE), Compumag 2023 (Kyoto)
Poster [525] "Comparison of Discretization Methods for Continuous Stream-Function
Distributions" (TU/e Electromechanics & Power Electronics; the Peeren SFM line).
Approach: fit a B-SPLINE SURFACE to the stream-function distribution (SFD) ->
psi value AND gradient at any (x,y) on the source surface.  Three connection
algorithms producing a CONNECTED (single-stroke) coil:
  - Method 1: cut-and-couple iso-contours over segments of LINEAR psi-decrease,
    using a 4th-ORDER POLYNOMIAL BLENDING FUNCTION to interpolate the joins.
  - Method 2: descend the fitted B-spline SFD surface point-to-point with a
    different constant gradient per step; user chooses start/end points.
  - Method 3: user-defined start/end + blending function over the NON-linear SFD
    (combines 1 and 2; most flexible).
Reported (planar coil, target at z=6 mm): Bx/By/Bz rms errors ~5-13%; methods 1
and 3 beat method 2 by >= 0.8 / 4.9 / 1.9 % in x/y/z.
KEY FINDING (the opening for us): "the area where the connected curve DEVIATES
from the iso-contour lines is coupled to the area in the target region with the
ERROR in the flux density."  They OBSERVE the connection -> field-error coupling
and SELECT the least-bad connection; they do NOT compensate it in the solve.
Related smoothness prior art: Liu/Hennig/Korvink, IEEE TM 48 (2012) 1179
(high-order-smooth discretised SF).

## Our research direction (paper idea, 2026-05-29)
CORE CLAIM: make the single-stroke connection PART of the field design, not a
post-processing step -- enabled by the ACA+-accelerated, kernel-agnostic solve
(cheap re-solve; works with magnetic materials / arbitrary surfaces that
free-space methods cannot treat).  Why ACA+ is essential: for material/expensive
kernels it cuts entry evaluations O(MN)->O(k(M+N)), making an outer
connection-optimisation loop tractable.

  A (backbone) -- CROSSOVER-COMPENSATED ITERATED LEAST-NORM.  Fix a connection
    topology tau -> its crossover field B_c(tau) is known; solve
    A psi = B_target - B_c(tau) and iterate (fixed point).  The parasitic
    connection field is designed away in the pattern.  Each iteration = one fast
    (ACA+)+TSVD solve.  Goal: single-stroke field == idealised disconnected-loop
    field.
  D (theory) -- MULTIVALUED POTENTIAL / BRANCH CUT.  Single-stroke == psi with a
    dI jump across a branch cut: psi~ = psi + (dI/2pi)*theta_tau.  Connection
    topology = choice of branch cut; the spiral-ramp join (psi - alpha*u) is the
    special case.  "Minimal multivalued correction making the contour a single
    connected curve with minimal field perturbation."
  B/C (search) -- field-cost combinatorial routing (min-cost Eulerian / TSP over
    crossovers, cost = parasitic DSV field / length / inductance) + BILEVEL
    (inner = ACA+TSVD amplitude, outer = CMA-ES/combinatorial connection).  This
    is the natural sequel to SA-25-020 "ACA stream function + CMA-ES".
  E (differentiator) -- KERNEL-AGNOSTIC generality: single-stroke design WITH
    magnetic materials (iron yoke / active shield) and on NON-cylindrical /
    conformal surfaces, using radia_field_kernel (MMM/MSC).  Free-space
    Biot-Savart SFM tools (Turner/Peeren/Kuijpers) cannot do this.

Positioning (honest): SFM, spiral/blending connection, and the
connection->error observation are PRIOR ART (Kuijpers et al.).  Novelty = folding
the connection field INTO the solve x ACA+ acceleration x materials/arbitrary
surfaces.  Confirm against literature before claiming.

## Status
Gz single-stroke = DONE (smooth helix, demo_sf_to_peec_gz.py).  Gx fingerprint
single-stroke = OPEN (future work).  Proposed demonstrators: Gx (compensated vs
naive vs ideal), shielded (iron / active shield), biplanar/saddle (surface
generality).
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
        "literature": ACA_TSVD_LITERATURE,
        "workflow": ACA_TSVD_WORKFLOW,
        "single_stroke": ACA_TSVD_SINGLE_STROKE,
    }
    aliases = {
        "kernel": "kernel_agnostic", "generic": "kernel_agnostic",
        "speed": "performance", "speedup": "performance", "benchmark": "performance",
        "optuna": "cmaes", "cma-es": "cmaes", "cma_es": "cmaes",
        "stream_function": "literature", "stream": "literature",
        "sfm": "literature", "lineage": "literature", "lit": "literature",
        "turner": "literature", "abe": "literature", "ducas": "literature",
        "peeren": "literature", "current_potential": "literature",
        "tsvd": "method", "aca": "method",
        "validate": "validation", "f90": "validation",
        "cad": "workflow", "peec": "workflow", "gz": "workflow", "gx": "workflow",
        "demo": "workflow", "pipeline": "workflow",
        "single-stroke": "single_stroke", "one_stroke": "single_stroke",
        "onestroke": "single_stroke", "winding": "single_stroke",
        "crossover": "single_stroke", "connect": "single_stroke",
        "kuijpers": "single_stroke", "fingerprint": "single_stroke",
        "spiral": "single_stroke", "manufacturable": "single_stroke",
    }
    t = aliases.get(t, t)
    if t == "all":
        return "\n\n".join([ACA_TSVD_OVERVIEW, ACA_TSVD_METHOD, ACA_TSVD_API,
                            ACA_TSVD_KERNEL_AGNOSTIC, ACA_TSVD_PERFORMANCE,
                            ACA_TSVD_CMAES, ACA_TSVD_VALIDATION,
                            ACA_TSVD_LITERATURE, ACA_TSVD_WORKFLOW,
                            ACA_TSVD_SINGLE_STROKE])
    if t in table:
        return table[t]
    return (f"Unknown topic '{topic}'.  Available: overview, method, api, "
            "kernel_agnostic, performance, cmaes, validation, literature, "
            "workflow, single_stroke, all.\n\n"
            + ACA_TSVD_OVERVIEW)
