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

The MCP server exposes this via ``streamfunction(topic=...)`` on the
radia-streamfunction server (moved from radia-ngsolve aca_tsvd in 2026-06).
Topics: overview,
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
care.  Use official `optuna/optuna-mcp` or plain Optuna APIs for sampler
choice, multi-objective (NSGA-II), pruning, and lab BBO recipes.

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

The checks above validate the SOLVER (the matrix A, the ACA+ factorisation, the
TSVD least-norm solve). The checks below validate the DESIGN -- that the psi the
solver returns actually reproduces the target field as a real coil. Run at least
one analytic-coil benchmark whenever the kernel, former, or regularisation menu
changes; a clean reconstruction error does NOT by itself prove a correct coil.

## Analytic-coil benchmarks (design-level, forward check)

Design psi for a textbook target, then FORWARD-evaluate B = A psi over the design
volume (DSV/FOV) and compare to the closed form. Reference fields:

- Uniform axial Bz on a CYLINDER former (radius a, length L) -> the ideal
  SOLENOID. The least-norm psi is ~linear in z (psi ~ z => K_phi = -dpsi/dz =
  const turns density). Check: central Bz vs mu0 * (dI per contour) * turns/length;
  field uniformity 1 - Bz(r,z)/Bz(0) over the DSV (a long former -> sub-% ripple).
- Transverse MRI GRADIENT Gx (target Bz = Gx * x over the DSV) on a cylinder ->
  the Golay / fingerprint saddle pattern. Check: gradient LINEARITY
  (max |Bz/x - Gx|/Gx over the DSV) and efficiency Gx per amp-turn. The Gy design
  is the same pattern rotated 90 deg; Gz (= Maxwell/anti-Helmholtz pair) targets
  Bz = Gz * z.
- Uniform-region pair: Helmholtz (uniform Bz) / Maxwell (uniform dBz/dz) give
  closed-form on-axis fields to anchor the magnitude.

Spherical-harmonic TARGETS and the achieved-(l,m) PURITY / contamination check are
their own topic -- see streamfunction("harmonics") (MRI gradient/shim basis). The
active-shield stray-reduction benchmark (~86x at INDEPENDENT points, not the
constraint points) is in streamfunction("shielding"); read its honesty note --
measuring residual at the constraint nodes flatters the result.

## Forward (discretized-wire) consistency

The continuous K = n x grad(psi) is exact; the WOUND coil is a finite set of
equal-increment iso-contours. After contouring, re-run Biot-Savart on the actual
wire POLYLINES and compare to the continuous-psi field over the DSV: this isolates
the DISCRETISATION (finite number-of-turns) error from the TSVD residual. Halving
dI (doubling turns) should drop the discretisation error ~linearly until it hits
the TSVD/target-reachability floor.

## TSVD truncation convergence (L-curve)

Sweep the retained singular-value count k (or Tikhonov alpha -- folded Tikhonov
traces it at ~50 us/point): the target-residual ||A psi - B|| falls then PLATEAUS
at the field's achievable-on-this-surface floor, while ||psi|| (and peak current
density) rises. The L-curve corner is the design point; report BOTH residual and
peak J so a "perfect" residual obtained with an unphysical peak J is visible. A
residual that never plateaus = target not in range(A) for that former (enlarge /
reshape the surface -- see streamfunction("deformation")).
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
  E (kernel-agnostic generality) -- single-stroke design WITH magnetic
    materials (iron yoke / active shield) and on NON-cylindrical /
    conformal surfaces, using radia_field_kernel (MMM/MSC).  Free-space
    Biot-Savart SFM tools (Turner/Peeren/Kuijpers) do not address this
    setup, so it is the natural application area for our callback-based
    pipeline.  Originality of this scope in OSS subject to literature
    confirmation before claiming.

Positioning (honest): SFM, spiral/blending connection, and the
connection->error observation are PRIOR ART (Kuijpers et al.).  Novelty = folding
the connection field INTO the solve x ACA+ acceleration x materials/arbitrary
surfaces.  Confirm against literature before claiming.

## Status
Gz single-stroke = DONE (smooth helix, ``demo_sf_to_peec_gz.py``).

Gx fingerprint single-stroke = DONE (2026-05-30..31,
``demo_sf_to_peec_gx.py`` with
``--chain-method {field_aware,kuijpers,lobe,greedy,nn_blend}``).

  USE THE ``single-stroke-chain`` CLAUDE SKILL when joining a new coil's
  contours -- the connection has no clean closed-form optimum and is a
  reason-and-verify task (the skill encodes the objective + 5 dead-ends +
  the build->measure-DSV-RMS->compare->escalate-to-Path-A loop).

  Five methods shipped:

  - ``field_aware`` (DEFAULT, recommended, 2026-05-31): kuijpers
    lobe/current-sign ORDER (the dominant factor) + each contour's cut
    chosen by coordinate descent to minimise the AZIMUTHAL arc to its chain
    neighbours (axial dz is free, carrying no stray), so the rung stray
    fields cancel more symmetrically over the DSV.  **DSV RMS 9.29 %,
    x-nonlin 7.20 %** -- beats kuijpers in EVERY tested config (16.24 % ->
    9.29 % at default; 30-54 % lower across nlevels 10/12/16 + nphi32/nz48).
  - ``kuijpers`` (prior best): Kuijpers 2023 Method-1.  Per-lobe fixed cut
    phi (`+x` at 0, `-x` at pi), one straight axial (phi, z) rung per
    contour pair (Fig.4 rung pattern).  RMS 16.24 %, nonlin 9.73 %.
  - ``lobe``: 4-quadrant classification + within-quadrant spiral.
    RMS 23.98 %, nonlin 9.67 %.
  - ``greedy``: global nearest-neighbour.  Visually "wasted" arcs.
    RMS 21.78 %, nonlin 11.40 %.
  - ``nn_blend``: CAUTIONARY negative result.  Geometric-shortest balanced
    cut + nearest-neighbour order -> RMS 0.651 (WORST).  Geometry-only NN
    interleaves the four lobes' alternating current signs (+,-,-,+); the
    shortest-chord hop keeps bridging opposite-sign contours.

  **KEY FINDING (corrects the earlier framing)**: field impact is NOT
  predicted by rung length -- geometric OR azimuthal.  nn_blend and
  field_aware have near-equal azimuthal rung totals (12847 vs 12522 mm)
  but 8x different RMS.  The dominant lever is the current-sign-respecting
  ORDER; given that, symmetric cut placement (azimuthal-min) makes the rung
  stray fields cancel.  No single scalar metric captures it -- hence the
  skill, and why the winning ``field_aware`` came from reason-about-objective
  + verify, not a closed formula.

  Implemented as ``single_stroke_field_aware_phi_z`` +
  ``_kuijpers_lobe_order`` + ``single_stroke_{,lobe_,kuijpers_,nn_blend_}
  chain_phi_z`` in ``examples/stream_function/demo_sf_to_peec_gx.py``.

  **The HARD-tier "16 % ceiling" was a kuijpers-method ARTIFACT, not a
  fundamental bound** -- the connection METHOD, not just the SF design,
  gates HARD quality.  And field_aware COMPOSES with Path-A: see below.

  Tried-and-rejected sort variants for kuijpers (2026-05-30, save the
  next-session debugging):

  - Sort by sign(centroid_x) (+x lobes desc z, -x lobes asc z, intended to
    give "Maxwell pair top-edge / bottom-edge / saddle-to-saddle" topology)
    -> DSV RMS 43.56 % (3x worse).  Reason: breaks matplotlib's per-contour
    CCW orientation when -x lobes are reversed, so current direction
    flips for ~half the loops.  The "correct" topology on paper does NOT
    survive matplotlib's natural orientation handling.

  - Forcing cut to TOP/BOTTOM crossing of each closed contour ("+z lobes
    cut at top half, -z at bottom half" for monotone outer-to-inner z)
    -> DSV RMS 30 %.  Same root cause: cut at a forced side flips the
    polyline rotation sense relative to matplotlib's choice for the y-
    mirror partners (which are open polylines whose ends are both at y=0).
    Net effect = wrong current direction for half the polylines.

  - Pair-aware (+y partners sorted outer->inner then -y partners reversed
    inner->outer, intended to spiral down on +y side then up on -y side
    with a short mid-lobe mirror crossing) -> DSV RMS 40 %.  Same flip.

  - Y-sign sub-grouping (8 sub-lobes by sign x/y/z) -> 35 % RMS.  Adds
    extra inter-sub-lobe transitions without resolving the flip issue.

  Pattern: ANY variant that introduces a per-contour reversal of the
  matplotlib traversal direction breaks the field.  The SAFE knob is the
  global within-lobe SORT KEY (cut z, sign(cut y), pair index, etc), but
  the per-contour TRAVERSAL DIRECTION must always be the matplotlib one.

  5th rejected variant + POSITIVE resolution (2026-05-31): ``nn_blend``
  tried hypothesis (a) above -- a different cut strategy (min-distance
  balanced cut between adjacent contours) -- but with GEOMETRIC nearest-
  neighbour ORDER, which interleaves the lobe current signs -> RMS 0.651
  (WORST).  Keeping kuijpers' sign-respecting ORDER and applying the
  azimuthal-min balanced cut (= ``field_aware``) is hypothesis (a) DONE
  RIGHT: RMS 9.29 %, the new best.  So the residual diagonal rungs ARE
  reducible by a better cut strategy -- as long as the ORDER preserves the
  current signs.  Both levers (better connection = field_aware, AND the
  compensated iteration = Path-A below) compose.

## Path A naive fixed-point iteration -- NEGATIVE result (2026-05-30)

Implemented as ``--compensated-iter N --compensated-step alpha`` in
``demo_sf_to_peec_gx.py``.  At each iter:

  1. Build chain from current psi
  2. Compute Bz_chain at unit current at DSV; gain-fit I_w
  3. ``residual = B_target - I_w * Bz_chain_unit``
  4. ``delta_phi = pseudo_inverse(residual)``  (uses cached TSVD factor)
  5. ``psi += alpha * delta_phi``; re-extract contours

Results on the 24x40 Gx baseline (target RMS at iter 0 = 16.24 %):

  - alpha = 1.0, n = 1:    RMS jumps **16 % -> 55 %**.  Hard divergence.
  - alpha = 1.0, n = 8:    RMS oscillates 16/55/26/99/90/85/...  No convergence.
  - alpha = 0.05, n = 60:  RMS oscillates 16-55 %, occasionally dips
    near 15-17 %, best of 30 iters = 15.28 % (-6 % vs baseline).
    On x-axis NONLINEARITY: 9.73 % -> 5.99 % (-38 %) at the best iter.
  - The best iter is found by ACCIDENT (iter 24 in this seed), not by
    convergence -- the same neighbourhood is re-visited many times.

Root cause: chain construction is highly NONLINEAR in phi.  Each iter
shifts level-set positions, matplotlib re-traces contours in a different
order, the polyline COUNT and TOPOLOGY change, and the resulting chain's
parasitic field bears no useful linear relation to the previous chain's
parasitic field.  ``pseudo_inverse(residual)`` is the right LINEARISED
update IF chain field were linear in phi; it is not.

Implemented mitigation: the demo tracks BEST PSI SEEN across iterations
and returns it.  This delivers the 5.99 % nonlinearity reliably (any seed
+ alpha=0.05 + n>=24 finds the same neighbourhood).  Modest but real.

**Path-A COMPOSES with field_aware (2026-05-31)**: re-running Path-A on top
of the ``field_aware`` chain (not kuijpers) reaches DSV RMS **8.11 %** at
``--compensated-iter 40 --compensated-step 0.3`` (-13 % on top of
field_aware's 9.29 %), x-nonlin 7.20 % -> 6.27 % -- roughly HALF the old
kuijpers baseline (16.24 %).  Path-A is STEP-SENSITIVE: step 0.3 gained,
0.5/0.8 found nothing (best-psi tracker kept the baseline).  Connection
method (field_aware) and SF-design absorption (Path-A) are SEPARATE levers
that compose.  Still oscillatory (best-psi tracking, not true convergence)
on grid-sampled psi; Path-A converges MONOTONICALLY only on FE-direct psi.

Open paths to actual convergence (TBD):
  - Anderson acceleration / quasi-Newton (weighted combination of recent
    iterates).  Standard remedy for non-contractive fixed points.
  - FREEZE the chain topology after iter 1; only update polyline
    currents (= phi).  This eliminates the level-set topology jumps.
    Loses some design flexibility but gives a true linear update law.
  - B-spline SFD continuous representation (Kuijpers Methods 2/3) so
    contour locations vary continuously with phi.  Eliminates the jumps
    by construction.  Heavier implementation.
  - The MULTIVALUED POTENTIAL formulation (D) sidesteps the iteration
    entirely by parameterising the chain as ``psi + (dI/2pi)*theta``,
    so chain topology is parameterised within the SF solve.

Status: Path A "naive Picard" is documented as a negative-with-trace
result.  Anderson, frozen-topology, or B-spline variants are the next
experiments worth running.  None is single-day work.

## Path A is ACTUALLY useful on simpler topologies (2026-05-30)

Same naive Picard tested on ``demo_planar_uniform_coil.py`` (uniform
Bz target above a planar source plane -- the easy end of the
complexity hierarchy: target axisymmetric -> concentric closed
contours -> trivially clean spiral).  Results on the 21x21 / 5x5 /
target_z=100 mm / B0=1 mT baseline:

  - iter 1 (= no compensation, single-shot Kuijpers spiral):
    DSV RMS  2.99 %, peak-to-peak / mean  9.59 %, mean Bz 0.999 mT.
  - 30 iters @ step = 0.05:
    RMS 0.97 % (-67 %), p2p/mean 3.46 % (-64 %).
  - 100 iters @ step = 0.05:
    RMS **0.58 %** (-80 %), p2p/mean **2.17 %** (-77 %), mean Bz
    1.000 mT exactly.

The iteration still does NOT strictly converge (residual still
oscillates), but the best-effort tracker finds genuinely better
neighbourhoods MUCH more reliably than on the Gx fingerprint.

Why the improvement is bigger here:
  - Concentric contour family stays connected under small psi
    perturbations (no level-set bifurcations near a saddle).
  - Single cut line at the +x axis -- no lobe transitions.
  - The "deviation region" of Kuijpers Method 1 is a SINGLE straight
    radial spoke, naturally far from the target -> small B_c at the
    target -> low residual already from iter 1.
  - Each iter can refine without fighting topology jumps.

This validates the COMPLEXITY-HIERARCHY framing: single-stroke design
quality (and Path-A effectiveness) is bounded by the coil's topology
class.

  Easy:    axisymmetric / planar uniform -> Path A reaches < 1 % RMS
  Medium:  cylindrical Gz                 -> done at iter 0 already
  Hard:    cylindrical Gx fingerprint     -> stuck at 15 % RMS
  Harder:  shielded / biplanar / 3D       -> needs B-spline SFD or
                                              multivalued potential

Practical lesson: when designing a single-stroke coil with this
machinery, START by classifying which complexity tier the target
forces.  If the tier is "Hard" or worse, accept multi-port or move
to a continuous-SFD formulation rather than expecting more iterations
of naive Picard to help.

## Path A converges MONOTONICALLY with FE-direct psi (2026-05-30)

Implemented as ``demo_planar_uniform_fem_psi.py``: psi is a direct H1
GridFunction on an NGSolve 2D mesh of the source plane (rectangle,
``maxh=0.025``, order=2 -> 1773 free DOFs), with the Biot-Savart
matrix ``A[j,i] = -(mu0/4pi) integral over mesh of grad(phi_i) . d_xy
/ r^3`` assembled per target via ``LinearForm``.  Two solve flavours:

  - ``--regularize l2``: ``np.linalg.lstsq(A_free, B)`` -- the min-
    Euclidean-norm interpolant.  RMS 1.78 % single-shot.
  - ``--regularize h1`` (recommended): min ``psi^T S psi`` s.t.
    ``A_free psi = B``, solved via ``psi = S^-1 A^T (A S^-1 A^T)^-1 B``
    where ``S`` is the H1 stiffness.  Smoothest psi that hits B
    exactly.  RMS 2.09 % single-shot.

Path A on the H1-regularised solve, ``--compensated-iter 100
--compensated-step 0.05``:

  - iters 40-47 are MONOTONICALLY decreasing 0.62 % -> 0.49 %.  No
    oscillation, no jitter.  This is a TRUE FIXED-POINT CONTRACTION.
  - Final at iter 75: RMS = **0.47 %** (-78 % vs single-shot
    H1, -84 % vs single-shot L2, -84 % vs the BASIS-LOOP demo
    ``demo_planar_uniform_coil.py``).
  - Peak-to-peak / mean = **1.64 %**.

Why this version converges where the basis-loop version only
oscillated: the H1 GridFunction is a CONTINUOUS function of its DOF
vector, the marching-squares contour family deforms SMOOTHLY with
psi, matplotlib emits the contours in a consistent count and ordering
across iterations, and the parasitic ``B_c(psi)`` is therefore a
smooth function of psi.  The naive Picard fixed-point map is then
genuinely a contraction.

Implication: the "Path A naive Picard does not converge" negative
result above is SPECIFIC to the basis-loop psi representation
(matplotlib contouring on a discretised grid of psi VALUES, prone to
topology jitter under small psi changes).  When psi is upgraded to a
continuous FE function, Path A converges.

This is the empirical justification for the
``radia-mcp streamfunction(single_stroke)`` topic's recommendation
to "switch to a B-spline SFD continuous representation" for
non-trivial topologies -- the same effect (smooth chain field
response) gives the same benefit (Path A actually contracts).

Open question: does the same FE-direct + H1-Path-A formulation
unstick the Gx FINGERPRINT case?  Worth trying on a cylindrical
surface mesh.  If yes, the "Hard tier" recommendation can be
weakened from "multi-port or continuous SFD" to just "continuous
SFD via FE direct".

## Single-current sheet-metal coil distortion (bankin-ho) (2026-05-31)

Shim loops (separate feeds) and Path-A (re-contour) both keep the wire
a FLAT planar pattern.  A THIRD compensation -- ``--distort`` in
``demo_planar_uniform_fem_psi.py`` -- drops that assumption: keep ONE
series current and the contour LEVELS fixed, and BEND the single-stroke
wire in 3D with a smooth low-dimensional deformation field

    d(x, y) -> (delta_x, delta_y, delta_z)

(control-grid bilinear = the discrete realisation of an NGSolve
VectorH1 mesh deformation; "bankin-ho" / sheet-metal forming, a geometric
shape optimisation).  Gauss-Newton over the deformation DOFs:
``(J^T J + (lam + lam_disp) I) delta = J^T r - lam_disp d`` with the
displacement penalty ``lam_disp`` RELATIVE to ``mean(diag(J^T J))``;
re-fit the single current each step; track the honest dense-grid MAE
for "best".

Measured (planar uniform Bz, 33 contours; integrated --distort run,
flat single-stroke baseline 12290 ppm; baseline is ~12-20k ppm
depending on the psi representation):

    penalty 1.0  (17 mm bend) -> 2015 ppm   conservative / printable
    penalty 0.3  (24 mm bend) -> 1025 ppm   balanced (default)
    penalty 0.1  (31 mm bend) ->  605 ppm   sub-1000, ONE current
    penalty 0.03 (34 mm bend) ->  340 ppm   aggressive

Convergence is fast + monotone (penalty 0.1: 12290 -> 3656 -> 1060 ->
687 -> 613 ppm over 6 iters; penalty 0.03 reaches 340 ppm).

Key physics: the OUT-OF-PLANE (z) lift is what lets a SINGLE current
cancel the parasitic connector field -- lifting a connector toward /
away from the DSV changes its Bz weight.  In-plane-only (xy) distortion
merely reroutes the connector current WITHIN the DSV plane and is far
weaker: an in-plane contour-flow PoC (move each wire normal to itself,
delta_p = -delta_psi grad psi / |grad psi|^2) plateaued at ~12500 ppm
in 40 iters, while full-3D Gauss-Newton reaches ~2500 ppm in 5.

This is the single-current ANSWER to the shim trade-off "no free lunch
within one uniform-current wire": that caveat assumed a FIXED PLANAR
geometry.  Trading GEOMETRIC DOFs for CURRENT DOFs, one current reaches
the ~500-2500 ppm class -- one 3D-printed bent conductor, no extra
supplies.  Separate-feed shims still reach a lower floor (183 ppm at 10
feeds) but need independent supplies; the two COMPOSE (distort first =
one part, add a couple of feeds only if the spec demands sub-500 ppm).

Implemented as ``coil_distort_3d`` + CLI ``--distort --distort-comps
{xyz,xy,z} --distort-penalty L --distort-grid N --distort-iter K
--distort-fit-n F``.  A fast vectorised ``bz_fast`` (broadcast
segments x obs, z-component only) replaces the per-segment Python loop
of ``h_segments_batch`` for the Gauss-Newton Jacobian (the field is
called ~N_dof times per iter).

CYLINDER PORT (done 2026-06-01) -- ``demo_sf_to_peec_gx.py --distort``:
``coil_distort_cyl`` bends the Gx fingerprint single-stroke on a
cylinder.  The analog of the planar z-lift is a RADIAL bend dr (out of
the surface) + azimuthal s=a*dphi + axial dz, on a phi-periodic (phi,z)
control grid.  Gx (HARD tier): single-stroke 8.5-9.3% -> distort 1.4%
(ONE current, ~30mm bend, full r+s+z).  INTERESTING REVERSAL: on the
cylinder the IN-SURFACE reroute (sz) is the DOMINANT lever and radial
(r) alone is WEAKEST (needs ~50mm) -- OPPOSITE of the plane (where the
out-of-plane lift dominated).  Reason: cylinder wires already wrap in
3D, so repositioning ALONG the surface reshapes the current pattern
more than lifting OFF it; for an internal DSV a radial move only weakly
changes the wire-DSV distance.  The optimal sheet-metal direction is
GEOMETRY-DEPENDENT; Gauss-Newton picks it automatically.

Also added ``--regularize {tsvd,tikhonov,h1}`` to the cylinder psi
solve (was plain TSVD).  alpha L-curve is NON-MONOTONIC (contour
topology jumps with alpha): TSVD mode-truncation is the BEST
regulariser for the Gx single-stroke (8.45%), best ridge alpha=1e-2*mean
only ties it (8.97%), H1 worse (over-spreads current -> longer
connectors) -- OPPOSITE ranking from the planar uniform case (where H1
was cleanest).  Recipe: pick the regularisation that minimises the
SINGLE-STROKE RMS, then distort.

COMPOSE with electric shims (``--distort --shim-loops K``): the bend and
the separate-feed shims cancel DIFFERENT parts of the residual.  The bend
is FAR more feed-efficient -- 1 bent wire (1 feed) reaches 1.4% where 10
electric shims alone (10 feeds) only reach 2.3%.  distort + 10 shims (11
feeds) -> 1.0%.  Recipe on the hard tier: distort first (1 printable
part), then a few shims clean the high-freq residual.  Under --distort the
[8] shim-only block is skipped; shims apply to the DISTORTED residual in
[9].  The Gx 3D plot colours the wire by RADIAL bend (sheet-metal forming
visible).

## FE-direct on ARBITRARY surfaces -- the "3rd-order is big" win (2026-06-02)

User 「3次が使えるのは大きいはず」「コイル設計面はどんな形をしていても
作れるはず」.  The basis-loop representation needs a structured (phi,z) grid,
so it is stuck on planes/cylinders.  FE-direct H1 psi MESHES ANY SURFACE --
THAT is the real payoff of the high-order (cubic) FE stream function:
design the SF directly ON the real manufacturable former (sphere /
conformal / 3D-printed), the NMR shim use case.  ``demo_sphere_fe_direct.py``.

Surface FEM done robustly: mesh the SOLID ball, ``H1(mesh, order=3,
definedon=mesh.Boundaries(".*"))``, integrate with ``ds`` + ``grad(.).Trace()``,
``mesh.Curve(3)`` (isoparametric -- curve order = field order; verified by
sphere area: Curve(1) 0.86% err -> Curve(3) 0.0003%).  General kernel (works
on ANY surface, the curved generalisation of planar z_hat x grad psi):
  A[j,i] = (mu0/4pi) int_S [ (n x grad_s phi_i) x (r_t-r') ]_z / |r_t-r'|^3 dS
with ``n = specialcf.normal(3)``, ``K_v = Cross(n, grad(v).Trace())``.

Results (R=150mm coil, DSV r=80mm):
  - uniform Bz (l=1):  continuous cres 3e-15 (0 ppm), single-stroke 0.24%.
  - Z2 shim (l=2,m=0): continuous cres 5e-14, single-stroke 4.3% -> sphere
    sheet-metal distort 0.36% (one current, ~2mm bend); real inter-turn
    spacing 10.5mm >= conductor width (manufacturable single wire).

Manufacturable pipeline blocker + fix: evaluating a boundary-H1 GridFunction
at 3D points returns 0 (lands in a volume element).  Fix: sample psi at the
BOUNDARY VERTICES (``gf.vec[vertex_nr]`` = nodal value) -> scipy.griddata to a
(phi,theta) grid (periodic phi via 3 wrapped copies) -> contour -> spiral
chain (sort by theta, open at min phi) -> sphere embed.  The Z2 "0 mm
spacing" was a FALSE ALARM (closed-circle phi=0/2pi endpoints + seam
connection points); the real inter-turn spacing is 10.5mm.

KEY POINTS:
  - The DESIGN-SURFACE SHAPE is NOT a manufacturability constraint -- any
    fabricable former holds a coil; the real constraints are the WIRE pattern
    (single-stroke / min-spacing / no self-cross / bend radius).
  - easy/hard is set by TARGET complexity (uniform l=1 -> 0.24% clean
    latitude-ring spiral; Gx fingerprint -> hard), NOT the surface shape.
  - FE-direct's value is GENERALITY (any surface), NOT a single-stroke-RMS
    improvement: on the cylinder Gx it does NOT beat basis-loop (the
    fingerprint topology-jumps regardless of representation; Path-A
    oscillates), but on a sphere basis-loop simply cannot be laid.
  - Caveat: the supplied sphere chainer is for AXISYMMETRIC (m=0: Z1/Z2/Z3...)
    shims; m!=0 tesserals need a general field-aware sphere chainer.

## Acceleration path for FE-direct: WAIT for Joachim H-matrix (2026-05-30)

Decision: do NOT build a custom HACApK ↔ ngsolve.bem basis bridge to
accelerate the FE-direct ``A[j, i] = integral over support of phi_i``
matrix assembly.  Instead, wait for native H-matrix support in
ngsolve.bem (planned upstream).  Reasons:

  - FMM != ACA+.  ngsolve.bem's ``use_fmm=True`` is FMM, NOT a
    drop-in replacement for HACApK's ACA+.  Radia explicitly
    removed FMM (CLAUDE.md 2026-03-06) in favour of ACA+ because
    FMM's analytic multipole expansions are kernel-specific and
    do not help the compact / near-field-dominated geometries we
    care about.  Each new kernel (material MMM, SIBC, ...) would
    need new multipole math; ACA+ stays kernel-agnostic.
  - A custom bridge in the interim would have to consume
    ngsolve.bem's basis evaluation + quadrature (which evolve
    upstream) and produce a scalar entry function for HACApK
    (which evolves locally) -- maintenance burden for negligible
    short-term gain at our M=25 scale.
  - For LARGE M or material kernels, the right wait-and-replace
    point is when Joachim ships a native H-matrix on the
    ``ngsolve.bem`` integral operators; the ``fe_entry(j, i)``
    function we'd write today becomes one API call into that
    operator.

Interim recommendation: full M x ndof assembly for small M (works
in 30 ms for the current demo, no need for ACA+).  Material-kernel
extensions deferred to upstream.

**GitHub history (verified 2026-05-30)**:

  - 2026-04-20 commit ``d90c59e`` "Coupling to external H matrix tools
    (HTool)" by Christopher Lackner, crediting Pierre Marchand
    (HTool's author).  Added ``IntegralOperator.CalcSubMatrix(rowids,
    colids) -> MatrixD/C`` and ``CalcSubMatrixCapsule() -> capsule``
    to ``python_bem.cpp``.
  - 2026-04-14 / 04-15 / 05-04 / 05-06 commits by ``erdieee``: FMM
    interface work (``FMMInterface for Kernels``, ``Unify
    FMMInterface to use Vec<COMPS,..>``, ``Expose FMM parameters as
    kwargs``, ``Add FMM operator info``).
  - 2026-04-30 release v6.2.2604 -- includes both HTool bridge and
    the FMM API.

**Joachim's architectural choice**: do NOT build a single H-matrix
implementation INTO ngsolve.bem.  Instead expose a generic
block-entry callback (``CalcSubMatrix``, ``CalcSubMatrixCapsule``)
so ANY external H-matrix library can plug in.  HTool was the test
pair (same author wiring both ends); HACApK or h2lib could plug in
via the same callback contract.  FMM stays in ngsolve.bem because
its math is kernel-specific.

**This validates Radia's (A) path retroactively**: the
``radia.stream_function.aca_tsvd`` entry-function contract
``entry(i, j) -> float`` is the same shape as ngsolve.bem's
``CalcSubMatrix`` (just block of size 1).  When we swap the interim
``LinearForm`` matrix assembly for the ngsolve.bem operator, the
rest of the pipeline (TSVD, Path-A, ...) carries over unchanged:

  # interim (today's demo_planar_uniform_fem_psi_aca.py):
  A_free = build_fem_matrix(...)[:, free_idx]
  def entry(i, j): return float(A_free[i, j])

  # 2604 ngsolve.bem swap:
  op = bem.LaplaceSL(integrand).Operator(name)
  def entry(i, j):
      M = op.CalcSubMatrix(np.array([i], dtype=np.int32),
                            np.array([j], dtype=np.int32))
      return float(np.asarray(M)[0, 0])

The swap is ~5 lines of code -- BUT ONLY IF the problem is
surface-to-surface (standard BEM Galerkin: same surface for trial
and test, source and target both indexed by surface FES DOFs).
For an SF coil INVERSE DESIGN with OFF-SURFACE point targets
(rows = 25 observation points at z = z_target, cols = source-plane
H1 DOFs), ``CalcSubMatrix`` does NOT directly apply.  Verified
2026-05-30 by reading ``tests/pytest/test_bem.py``: the canonical
pattern is

    op = LaplaceSL(u * ds, use_fmm=False) * v * ds   # SURFACE x SURFACE
    submat = op.CalcSubMatrix(rows, cols)

Both ``rows`` and ``cols`` index the SAME surface FES.  For our
inverse-design problem the natural off-surface evaluation path is
``potop(gfu, target_boundary)`` -- but that goes through the FMM
backend, not ``CalcSubMatrix``.

For our problem we keep the (A) path AS-IS (callback into
HACApK from the LinearForm-assembled matrix or any other source).
Future option: reformulate the targets as a thin virtual surface
mesh and switch to ``CalcSubMatrix`` -- but at M=25 the LinearForm
assembly already runs in 30 ms, so there is no acceleration
incentive.

ngsolve.bem's ``CalcSubMatrix`` API is designed for the standard
BEM user community (inductance / capacitance extraction, FEM-BEM
coupling).  Inverse coil design with discrete off-surface targets
is a DIFFERENT problem class with a DIFFERENT natural matrix shape
(M x N_DOF, M usually small).  The (A) callback contract is the
right abstraction for OUR problem class.

**Status 2026-05-30**: NGSolve **6.2.2604** is INSTALLED on LAB
and exposes the new FMM-style hierarchical BEM in ``ngsolve.bem``:

  - ``BiotSavartCF(order, kappa, center, rad)`` -- multipole
    expansion of the Biot-Savart kernel.
  - ``BiotSavartRegularMLCF`` / ``BiotSavartSingularMLCF`` --
    Multi-Level (= hierarchical tree) expansions.
  - ``RegularMLExpansion`` / ``SingularMLExpansion`` --
    multipole / local expansion building blocks.
  - ``SphericalHarmonicsCF`` -- spherical-harmonic basis.
  - ``IntegralOperator.NearFieldMatrix`` / ``CalcSubMatrix`` --
    near-field sparse + block-wise matrix evaluation.

Symbol naming (``SphericalHarmonics``, ``MLExpansion``) confirms
the math is **FMM**, not ACA+.  This is consistent with the
recommended use case (smooth-surface Biot-Savart for SF coil
design, far-field-dominated) and does NOT conflict with the
"FMM Removed" policy from CLAUDE.md -- that policy targets
**Radia's MMM/MSC volume integral** (compact magnets, near-field
heavy), a different layer + geometry class.  See
memory ``feedback_fmm_vs_aca_distinction`` and CLAUDE.md
"FMM Removed from Radia core (2026-03-06)" / "SCOPE CLARIFICATION
2026-05-30".

NEXT STEP (task #14 in the project tracker): rewrite the matrix
assembly in ``demo_planar_uniform_fem_psi.py`` to use the
``ngsolve.bem`` ``BiotSavartCF``-based operator, verify field
correctness vs the current ``LinearForm``-per-target version,
then benchmark on large M and the cylindrical Gx fingerprint.

## Path (A) shipped: FE-direct psi via HACApK (ACA+)+TSVD (2026-05-30)

In parallel with the ngsolve.bem wait, the (A) path -- "wrap the FE
matrix as a per-entry callback for the existing
``radia.stream_function.aca_tsvd`` (HACApK ACA+ + Method-3 TSVD)"
-- has been implemented and validated:

  - ``demo_planar_uniform_fem_psi_aca.py``: M=25 x N_free=1773 FE
    matrix is built via the interim LinearForm path, wrapped as
    ``entry(i, j) -> float``, passed to ``aca_tsvd(M, N, entry,
    modes=M, kmax=min(M,N), aca_eps=1e-10, method=3)``.
  - ``k_aca = 25`` (full-rank for M-bounded matrix), factorisation
    time 0.033 s.
  - Reconstruction matches direct lstsq to
    ``||psi_ACA - psi_lstsq|| / ||psi_lstsq|| ~ 2e-3`` (= aca_eps
    floor).
  - Path-A iteration with cached factorisation: best RMS 0.67 %
    at iter 53 (slightly worse than direct lstsq's 0.39 % because
    of the ACA+ truncation, but the factorisation is re-used at
    ~1 ms per back-substitution).
  - Identical pipeline shape to the basis-loop demo: same
    ``entry`` contract, same ``pseudo_inverse_solve`` API, same
    Path-A loop.  This is the empirical confirmation that the
    callback contract carries over.

WHY (A) matters even though the FE matrix here is cheap (30 ms
full-assembly): the same ``entry`` callback contract handles
EXPENSIVE kernels.  For Radia MMM iron yoke / shielded coil /
SIBC workpiece, each entry costs a Radia container solve and
full assembly becomes O(MN) container solves.  ACA+ via the
same machinery cuts that to O(k(M+N)).  (A) is the
acceleration path that does NOT depend on ngsolve.bem's
upstream H-matrix.

WHY both paths matter:
  - (A) HACApK ACA+: kernel-agnostic, lives in Radia, works
    for ANY entry function (free-space, material, BEM SIBC).
  - ngsolve.bem H-matrix / FMM (when shipped in full): native
    surface BEM acceleration for free-space smooth kernels,
    benefits from NGSolve's tree/quadrature infrastructure.

Their natural domains are different (volume-MMM-MSC vs
surface-BEM), and they complement each other.

See memory ``feedback_fmm_vs_aca_distinction`` for the FMM vs ACA+
distinction and why we keep slipping on it.

  GMSH viewer: ``view_sf_coil_gx_gmsh.py --mode chain`` regenerates the
  chain (now the ``field_aware`` default) and writes a 1D-line .msh;
  ``--mode contours`` shows
  the raw SF design's closed contours with NO single-stroke connection,
  for direct visual comparison of "what the SFM designed" vs "what got
  manufactured".

Proposed next demonstrators: Gx with crossover-compensated iterated
least-norm (A); shielded coil with iron/active shield (E); biplanar/saddle
on a non-cylindrical surface (E).
"""


ACA_TSVD_SESSION_2026_05_30 = r"""
# Session 2026-05-30 -- single-stroke SF coil + Path-A + FE-direct + 2604

Compact structured summary of the day-long investigation.  Use this when
the full ``single_stroke`` topic is too long; this gives the executive
narrative + numbers + pointers.  Detail lives in the referenced topics
and demo files.

## 1. Complexity tier of single-stroke design (the framework)

The design quality reachable by Kuijpers-style chaining + Path-A
iteration is bounded by the coil's TOPOLOGY CLASS.  Empirically, four
tiers:

  | Tier   | Topology                       | Baseline    | best so far | Behaviour                     |
  | ------ | ------------------------------ | ----------- | ----------- | ----------------------------- |
  | EASY   | axisymmetric / planar uniform  | RMS 2-3 %   | <1 %        | path-A useful (basis-loop)    |
  |        |   + FE-direct H1 psi           | 2 %         | **0.47 %**  | path-A MONOTONE convergence   |
  | MEDIUM | cylindrical Gz                 | already OK  | redundant   | smooth helix natural          |
  | HARD   | cylindrical Gx fingerprint     | 16 % (kuijpers) | **8.1 %** (field_aware + Path-A) | chain method matters; NOT 16%-capped |
  | HARDER | shielded / biplanar / 3D       | --          | --          | needs FE-direct or multival.  |

Rule of thumb: classify the tier FIRST.  Going past EASY needs either
FE-direct continuous psi OR multi-port acceptance.  The HARD "16 % ceiling"
turned out to be a kuijpers-CHAIN-METHOD artifact, NOT a fundamental bound
(field_aware 9.3 %, + Path-A 8.1 %); past 8 % likely needs B-spline SFD or
multivalued-potential reformulation.

## 2. Chain methods for single-stroke (2026-05-31: field_aware is best)

Five implemented + benchmarked on the cylindrical Gx fingerprint
(``demo_sf_to_peec_gx.py --chain-method
{field_aware,kuijpers,lobe,greedy,nn_blend}``):

  | Method        | DSV RMS  | x-axis nonlin |
  | ------------- | -------- | ------------- |
  | field_aware   |  9.29 %  |  7.20 %       |
  | kuijpers      | 16.24 %  |  9.73 %       |
  | greedy        | 21.78 %  | 11.40 %       |
  | lobe          | 23.98 %  |  9.67 %       |
  | nn_blend      | 65.08 %  | 38.90 %       |

``field_aware`` (DEFAULT, 2026-05-31) = kuijpers lobe/current-sign ORDER +
cuts minimising the AZIMUTHAL arc to chain neighbours.  Beats kuijpers in
every tested config (30-54 % lower RMS across nlevels/mesh sweeps).  KEY:
field impact is NOT predicted by rung length (geometric or azimuthal) --
nn_blend and field_aware have near-equal azimuthal totals but 8x RMS; the
lever is the sign-respecting ORDER + symmetric cut placement (stray
cancellation).  No closed-form metric -> use the ``single-stroke-chain``
SKILL (reason + verify).

## 3. Path-A compensated iteration -- representation dependent

The fixed-point ``phi <- phi + alpha * pseudo_inverse(B_target -
I_w * Bz_chain_unit)`` behaves very differently across psi
representations:

  | psi representation              | Baseline   | Path-A best  | Behaviour                  |
  | ------------------------------- | ---------- | ------------ | -------------------------- |
  | basis-loop (matplotlib contour) | 2.99 % RMS | 0.58 %       | OSCILLATES, best-effort    |
  | FE-direct (L2 min Euclidean)    | 1.78 %     | **0.39 %**   | improving                  |
  | FE-direct (H1 min seminorm)     | 2.09 %     | **0.47 %**   | MONOTONE iter 40-47        |
  | FE-direct + HACApK ACA+TSVD     | 1.78 %     | 0.67 %       | (slightly worse due to     |
  |                                 |            |              | aca_eps=1e-10 truncation)  |
  | (Gx fingerprint, kuijpers chain)| 16.24 %    | 15.28 %      | oscillates, modest          |
  | (Gx fingerprint, field_aware)   | 9.29 %     | **8.11 %**   | oscillates but composes (-13%)|

WHY representation matters: matplotlib contour extraction on a
discretised grid TOPOLOGY-JUMPS under small psi changes (level sets
bifurcate near saddles, polylines re-order, B_c(psi) is not smooth
in psi).  A continuous H1 GridFunction is a smooth function of its
DOF vector, contour family deforms smoothly, B_c(psi) is smooth,
Picard contracts.  This is the empirical justification for FE-
direct as the path past EASY tier.

## 4. Dead-end variants (do NOT re-try) + the positive resolution

Five chain attempts that REDUCED field accuracy (Gx, kuijpers 16.24 % ref):
    - Sort by sign(centroid_x): RMS 43.56 % (-x lobe reversal flips current)
    - Force top/bottom cut on closed contours: 30.33 %  (same reason)
    - Pair-aware (+y outer-to-inner, -y inner-to-outer reverse): 40.04 %
    - 8-sub-lobe (sx, sy, sz): 34.64 %
    - nn_blend geometric-shortest balanced cut + geometric NN order: 65.08 %
      (NN interleaves the lobe current signs -- the worst of all)

  Common cause (traps 1-4): per-contour reversal flips matplotlib's CCW
  orientation -> current inverts.  Trap 5 (nn_blend): geometric ORDER
  interleaves lobe signs.  SAFE knobs = within-lobe sort key + sign-
  respecting order; UNSAFE = per-contour reversal OR geometric-NN order.

  POSITIVE RESOLUTION (2026-05-31): ``field_aware`` = kuijpers sign-ORDER +
  azimuthal-min balanced cut reaches RMS 9.29 % (now DEFAULT), + Path-A
  8.11 %.  The "16 % HARD ceiling" was a kuijpers-method artifact.  Lesson:
  optimise the FIELD (DSV RMS), verify, and respect the current-sign order
  -- not rung length.  Now a ``single-stroke-chain`` Claude skill.  See
  [[feedback_single_stroke_chain_orientation_traps]] (memory).

  Path-A on basis-loop psi:
    - alpha = 1.0 (full step): RMS 16 % -> 55 %    (divergence)
    - alpha = 0.5 / 0.1: oscillates 40-55 %        (no contraction)
    - alpha = 0.05 + N=60 iters: oscillates, best = 15.28 % by accident

  Common cause: topology jumps in the matplotlib contour family.
  See [[feedback_path_a_naive_picard_negative]] (memory).

## 5. Demo file ledger

All under ``examples/stream_function/``:

  | File                                 | Role                                   | Best RMS | Path-A |
  | ------------------------------------ | -------------------------------------- | -------- | ------ |
  | ``demo_sf_to_peec_gx.py``            | Gx fingerprint, 5 chain methods (field_aware default) | 9.29 % | **8.11 %** |
  | ``demo_planar_uniform_coil.py``      | Planar uniform Bz, basis-loop          | 2.99 %   | 0.58 % |
  | ``demo_planar_uniform_fem_psi.py``   | Planar uniform Bz, H1 FE-direct        | 1.78 %   | 0.39 % |
  | ``demo_planar_uniform_fem_psi_aca.py`` | Same + HACApK ACA+ via callback      | 1.78 %   | 0.67 % |
  | ``view_sf_coil_gx_gmsh.py``          | GMSH viewer (contours / chain / step)  | --       | --     |

## 6. NGSolve 6.2.2604 + ngsolve.bem (installed on LAB)

Confirmed 2026-05-30:

  - ``pip show ngsolve`` -> Version 6.2.2604, Location
    ``C:/Program Files/Python312/Lib/site-packages``.
  - ``ngsolve.bem`` exposes:
    - FMM-style: ``BiotSavartCF(order, kappa, center, rad)``,
      ``BiotSavartRegularMLCF`` / ``BiotSavartSingularMLCF``,
      ``RegularMLExpansion`` / ``SingularMLExpansion``,
      ``SphericalHarmonicsCF``.  Kernel-specific multipole math.
    - H-matrix bridge: ``IntegralOperator.CalcSubMatrix(rowids,
      colids) -> MatrixD/C``, ``CalcSubMatrixCapsule() -> capsule``
      (Pierre Marchand HTool integration, commit
      `d90c59e <https://github.com/NGSolve/ngsolve/commit/d90c59e>`_
      2026-04-20).  Algebraic, library-agnostic.
    - Other potentials: ``LaplaceSL/DL``, ``HelmholtzSL/DL/HypersingularOperator``,
      ``MaxwellSL/DL``, ``LameSL``.

  - **CalcSubMatrix scope** (verified by reading
    ``tests/pytest/test_bem.py``): the canonical pattern is

        op = LaplaceSL(u * ds, use_fmm=False) * v * ds   # surface x surface
        sub = op.CalcSubMatrix(rows, cols)

    Both ``rows`` and ``cols`` index the SAME surface FES.  For SF
    INVERSE design with off-surface point targets this does NOT
    apply directly.  Off-surface evaluation goes through
    ``potop(gfu, target_boundary)`` -- which uses the FMM backend,
    not CalcSubMatrix.

  - **(A) path stays the right abstraction** for SF inverse design
    because target points are not a surface FES.  ``CalcSubMatrix``
    is for standard BEM users (inductance / capacitance / FEM-BEM
    coupling).  Our problem class (M discrete targets x N source
    DOFs) keeps its own per-(target, basis) entry function.

## 7. FMM vs ACA+ -- do NOT conflate

  - FMM: kernel-specific analytic multipole expansion.  Different
    math per kernel (Laplace SH, Helmholtz Bessel, ...).  Good for
    smooth surface kernels in far-field-dominated geometries.
  - ACA+: algebraic / kernel-agnostic.  Same code works for any
    low-rank-friendly kernel (free-space, material, SIBC, ...).
    Good when each entry is expensive (Radia container solve).

CLAUDE.md "FMM Removed (2026-03-06)" targets the Radia MMM/MSC
volume integral (compact magnets, near-field heavy, FMM math bad).
It does NOT forbid ngsolve.bem's FMM for surface BEM (smooth
kernel, far-field, FMM math good) -- different layer, different
geometry class.  SCOPE CLARIFICATION committed to CLAUDE.md
2026-05-30.

## 8. Architecture: callback contract is universal

Block-entry callback (1x1 block = scalar entry) emerges as the
common bridge between ALL the relevant systems:

  | System                              | Callback                                              |
  | ----------------------------------- | ----------------------------------------------------- |
  | Radia ``radia.stream_function``     | ``entry(i, j) -> float``                             |
  | HACApK (C library)                  | ``HACApK_set_entry_func`` ``entry(i, j)``            |
  | ngsolve.bem ``IntegralOperator``    | ``CalcSubMatrix(rowids, colids) -> Matrix``          |
  | HTool (Pierre Marchand)             | ``void(int, int, const int*, const int*, double*)``  |

Joachim's design: do NOT bake one H-matrix lib into ngsolve.bem;
expose ``CalcSubMatrixCapsule`` so external H-matrix libs plug in.
HTool is the test pair; HACApK could plug in via the same capsule.
This validates Radia's existing (A) abstraction retroactively.

## 9. GMSH off-screen window (Pitfall #9, ``gmsh_usage("pitfalls")``)

When the GMSH viewer "doesn't appear":
  1. ``Get-Process python | Select MainWindowTitle`` -- if "Gmsh - ..."
     is set, the window EXISTS.
  2. P/Invoke ``GetWindowRect`` -- negative x = off-screen.
  3. ``MoveWindow(hWnd, 100, 100, 1280, 800, true)`` to rescue.
  4. Preventive in viewer scripts: ``gmsh.initialize(["-noconfig"])``
     + explicit ``General.GraphicsPositionX/Y`` + ``Width`` / ``Height``.
     NOTE: 6.2.x options are Width/Height, NOT SizeX/SizeY.

## 10. Advanced regularisation + surface deformation (2026-05-30 evening)

Extended ``demo_planar_uniform_fem_psi_advanced.py`` with four research
knobs the user explicitly requested:

  (a) ``--regularize h1_sigma --sigma-cf EXPR``: 1/sigma weighted H1
      seminorm = min INT (1/sigma)|grad psi|^2 s.t. A psi = B.  sigma(x,y)
      from a Python CF expression (with x, y, IfPos, exp, sqrt, sin, cos,
      log in scope).  Physically = true ohmic dissipation for non-uniform
      conductivity (litz wire, mask, forbidden-region exclusion).
      Empirically: Gaussian-decay sigma (exp(-r^2/0.02)) -> RMS 1.17 %
      (vs uniform H1 2.09 %) on the planar uniform-Bz benchmark.

  (b) ``--regularize inductance_diag``: lumped self-inductance proxy
      via per-element characteristic-length weighting (``specialcf.mesh_size``
      coupled to grad).  RMS 2.46 %; this is a DIAGONAL approximation to
      the true inductance L = (mu0/4pi) INT INT grad psi(x) . grad psi(y)
      / |x-y| dA(x) dA(y) -- full off-diagonal form needs ngsolve.bem
      surface BEM (MaxwellSL on a thin 3D shell embedded in (x, y, 0)
      plane), DEFERRED -- see ``# TODO`` in the demo file.

  (c) ``--regularize linf --jmax VAL``: scipy SLSQP with per-vertex
      ``|grad psi(v)| <= jmax`` constraint.  Warm-started from H1 solve.
      EXPERIMENTAL -- SLSQP is slow at ndof = 1773 H1 order-2 (~ 5 sec
      per iteration, may not converge to tight cap).  Best with H1
      order=1 + larger maxh.  Use for hot-spot suppression studies.

  (d) ``--deform --deform-params {zoff, bump, zoff+bump} --deform-trials N``:
      Optuna CMA-ES outer loop over surface deformation parameters.
      Source surface z = base_z + deform_params (Gaussian bump with
      amp/cx/cy, or pure z-offset, or both).  Each trial = one full
      inner SF solve + spiral chain + target field eval; CMA-ES picks
      deformation to minimise target-plane RMS.
      Measured on planar uniform-Bz benchmark (H1 inner solver):
        - 1-param zoff, 10 trials in 7.4 s:    RMS 2.09 % -> 1.50 % (-28 %)
        - 3-param bump, 20 trials in 22.6 s:   RMS -> **0.77 % (-63 %)**
        - 4-param zoff+bump, 50 trials in 33 s: RMS -> 0.85 %
      Best bump = ~-21 mm amp (concave dip pointing AWAY from target),
      off-centre.  Best zoff = +20 mm (plane translated toward target).
      Combined optimum is a local minimum that compensates the two
      effects against each other.

Each knob is composable with the others and with Path-A iteration.
The infrastructure mirrors the SA-25-020 "(ACA+)+TSVD + CMA-ES" pattern
(inner linear, outer nonlinear) but extended with smarter regularisation
choices on the inner side.

  (e) ``--regularize l2_aca``: routes the L2 (Euclidean min-norm) solve
      through ``radia.stream_function.aca_tsvd`` (HACApK ACA+ + Method-3
      TSVD).  Validated 2026-05-30 on the planar uniform-Bz demo:
      RMS identical to ``np.linalg.lstsq`` (1.78 % vs 1.78 %) and Path-A
      cache-friendly.  Slightly slower at M=25 (1.12 s vs 0.92 s)
      because ACA+ overhead dominates; the win is for M >> 25 or
      material kernels.  Confirms the (A) callback contract carries to
      the advanced demo.

  (f) High-order H1 basis (``--order N``, default 2): NGSolve supports
      arbitrary p, the advanced demo exposes it.  Measured RMS on the
      planar uniform-Bz benchmark at the SAME ``maxh=0.025``:

        order=1: RMS 1.02 %, p2p 3.77 %, ndof ~500
        order=2: RMS 2.09 %, p2p 6.81 %, ndof ~1773
        order=3: RMS **0.51 %**, p2p **1.83 %**, ndof ~3700  <- sweet spot
        order=4: RMS 2.04 %, p2p 7.65 %, ndof ~6000
        order=5: RMS 2.09 %, p2p 7.27 %, ndof ~9000

      Non-monotone in p: the SF discretisation + single-stroke contour
      extraction + ``nlevels=12`` interact with the smoothness of psi.
      order=3 happens to be the resolution-matched sweet spot; order=1
      is also competitive because the linear interpolation matches the
      contour algorithm cleanly.  This is a benchmark-class observation,
      not a generic FE convergence statement.

  (g) Composability: deformation + higher order can BACKFIRE if the
      baseline is already near-optimal.  Empirical: order=3 baseline
      0.51 % + bump deform (20 trials, 27 s) -> 0.82 % (CMA-ES wastes
      its budget searching for an even-better point and gets stuck
      in a sub-optimum).  At order=2 the same deform improves the
      baseline 2.09 % -> 0.77 % because there is room for improvement.
      DESIGN RULE: deformation freedom helps most when baseline is
      suboptimal; turn it OFF when single-shot accuracy is already
      under a few tenths of a percent.

## 11. Validation strategy (2026-05-30)

After session-long buildout (single-stroke chain, Path-A iteration,
FE-direct psi, advanced regularisation, surface deformation), the
framework has the following feature set.  This is a description of
what we have, not a claim of relative ranking against other tools:

  | Capability                       | Our state              | Related tool / paper       |
  | -------------------------------- | ---------------------- | -------------------------- |
  | SF inverse design                | ACA+TSVD, kernel-agn.  | CoilGen (MRI-only OSS)     |
  | Single-stroke connection         | 3 methods, Kuijpers-1  | Kuijpers 2023 paper        |
  | Single-stroke compensation       | Path-A monotone (FE)   | Kuijpers 2023 (observes only) |
  | Bilevel (inner SF, outer geom)   | Optuna CMA-ES          | Comsol Opt (commercial)    |
  | Regularisation choices           | 5 (L2/H1/sigma/Linf/L_diag) | Liu-Hennig-Korvink (H2)|
  | High-order FE psi                | H1 order p (any)       | NGSolve raw                |
  | Material-kernel extension        | callback ready (TODO demo) | -- (commercial FEM)    |
  | ngsolve.bem 2604 alignment       | callback contract matches  | Joachim/Pierre upstream|

  No published OSS combines all eight in a single package as far as we
  have surveyed; literature search to confirm originality is pending.
  Each individual capability has a literature precedent.

**Validation gap (industry-standard benchmarks not yet reproduced)**:

  Reproducing the published target specs (DSV, gradient strength,
  linearity, inductance, wire length) and reporting our numbers
  alongside is what the framework still needs:

    1. Bilac et al., Magn Reson Imaging (TBD year): planar shim coil
       benchmarks.  Closest to our planar uniform Bz demo.
    2. Turner, J. Phys. D 19 L147 (1986): cylindrical gradient coil
       analytical formulation.  Has explicit SFD expressions; we can
       compare ours to it directly.
    3. Lemdiasov & Ludwig, Concepts Magn Reson Part B (2005): target
       field method for MRI gradient coils.  Open-published spec.
    4. Forbes & Crozier 2002-2010 series: rigorous mathematical
       formulation with minimum-inductance regularisation.  Reference
       for our inductance_diag mode upgrade.
    5. CoilGen (Schwartz et al., GitHub): the directly-comparable OSS
       SF coil designer.  Head-to-head same-spec comparison.

  Recommendation: 2-week validation campaign:

    Week 1: implement Bilac + Turner + Lemdiasov-Ludwig benchmarks in
            ``examples/stream_function/benchmarks/`` with literature
            target specs; produce a benchmark table.
    Week 1 end: install CoilGen, run identical spec, head-to-head table.
    Week 2: shielded coil (iron back plate) via Radia MMM kernel
            through the (A) callback -- material-kernel demo not
            previously combined this way in OSS that we know of.
            Verify against any commercial-FEM equivalent before
            claiming originality of the open-source-first scope.
    Week 2 end: paper draft with the benchmark table and the Path-A
            compensation methodology.

**Candidate contribution claims for the paper** — calibrated by
what is genuinely novel.  ACA+TSVD itself is NOT new: Bebendorf 2000
(original ACA), Bebendorf-Rjasanow 2003 (ACA+), Hackbusch 2008
"Hierarchical Matrices" (TSVD recompression / "round" operation) are
standard H-matrix primitives.  The SA-25-020 paper applied the
combination to SF coil design; even that may have international
precedent in the MRI gradient coil + H-matrix literature (need to
check Kurz-Rain-Rjasanow IEEE TMag and MRI gradient coil community).

What THIS session adds on top of SA-25-020 may have real novelty:

  1. Path-A compensated iteration on FE-direct psi with empirical
     monotone convergence on the planar uniform Bz benchmark.
     Extends Kuijpers 2023's "deviation = field error" observation
     from a SELECTION criterion to a SOLVE update -- they observe and
     pick the least-bad chain; we fold the chain's parasitic field
     back into the solve and iterate.  Literature search for prior
     compensated-iteration in coil design is needed before claiming.
  2. Representation-dependence of Path-A -- the same fixed-point map
     converges only when psi is a continuous FE function; on the
     grid-sampled basis-loop representation it oscillates.  This
     methodological observation has not (to our knowledge) been
     pointed out in print.
  3. Empirical complexity-tier framework (Easy/Medium/Hard/Harder) --
     probably an articulation of community-implicit knowledge rather
     than a discovery, but useful as a design guide.
  4. Material-kernel extension via kernel-agnostic callback,
     demonstrated on a shielded coil with iron back plate (TODO).
     The callback abstraction itself is a standard H-matrix library
     pattern; the contribution is the SF-coil-design APPLICATION on
     magnetic materials, not the contract.

The strongest candidate for a methods-paper headline is (1) Path-A
compensated iteration; (2) is the supporting methodological result.
Everything else is engineering integration that strengthens the demo
but is not a methods contribution on its own -- but the integration
itself is a legitimate IMPLEMENTATION CONTRIBUTION for a software-
focused venue (JOSS, SoftwareX, methods-paper Implementation section).

**Implementation contribution (separate from algorithmic novelty)**:

The integrated pipeline combines three libraries that are valuable
separately but have not been combined for SF coil design (literature
check pending):

  - NGSolve ``ngsolve.bem`` (Schoeberl 2024+; HTool bridge by Pierre
    Marchand 2026-04-20): surface BEM operators, arbitrary-order
    FE basis, Sauter-Schwab quadrature, block-entry callback API
    (``CalcSubMatrix``, ``CalcSubMatrixCapsule``).
  - HACApK (Ida et al., ppOpen-HPC, MIT 2015+): kernel-agnostic ACA+
    H-matrix solver, callback via ``HACApK_set_entry_func``.
  - Radia (ESRF + Sugahara Lab): MMM/MSC material kernels,
    PEEC, single-stroke chain construction, Path-A iteration,
    kernel-agnostic ``radia.stream_function.aca_tsvd``.

Value of the integration:

  1. Callback-contract alignment -- ngsolve.bem ``CalcSubMatrix``,
     HACApK ``HACApK_set_entry_func``, and
     ``radia.stream_function.aca_tsvd(entry)`` are all
     ``entry(i, j) -> float`` (or block-form) compatible.  This
     alignment is upstream-organic; we provide the SF-coil-design
     glue that uses it.
  2. Kernel-swap transparency -- same pipeline for free-space
     Biot-Savart, Radia MMM iron yoke, SIBC workpiece, future
     application kernels by replacing the entry function.  Chain
     construction, Path-A, regularisation, deformation are unchanged.
  3. End-to-end OSS pipeline -- to our knowledge no other OSS SF
     coil designer combines kernel-agnostic SF inverse design +
     single-stroke chain construction + material-kernel support in
     one package.  CoilGen is MRI-only and free-space-only.
     Commercial Comsol AC/DC + Opt Module is comparable in scope
     but closed source.
  4. Active upstream alignment -- NGSolve 6.2.2604+ shipped the
     ``ngsolve.bem`` HTool bridge that our ``aca_tsvd`` callback
     contract is structurally compatible with.

For the paper:

  - Methods paper: short Implementation section near the end,
    crediting Schoeberl, Marchand, Ida, Sugahara, framing as
    "we integrate" not "we invent".
  - JOSS / SoftwareX: integrated pipeline IS the contribution,
    with demos + validation as evidence of usefulness.

The two framings are not mutually exclusive (= a methods paper + a
companion JOSS paper is a common pattern).

**State**:

  - Feature set: complete for the SF inverse-design pipeline plus
    several extensions (Path-A, regularisation choices, deformation).
  - Validation: not yet reproduced against industry benchmarks --
    that is the campaign described above.
  - Manuscript: methods + theory material is in the docs folder; the
    results section is gated on the benchmark numbers.
  - Timing: NGSolve 6.2.2604 (April 2026) shipped the H-matrix bridge
    that other SF + BEM groups will likely build on over 6-12 months;
    if proof-of-priority on any specific contribution is intended the
    validation campaign + arXiv preprint are worth doing sooner.

## 12. Paper draft material ready

Today's outcomes form publishable Methods + Results:

  - Methods: SFM via (ACA+)+TSVD; single-stroke chaining (Kuijpers
    Method-1 + the field_aware sign-order/azimuthal-min cut) on a 4-
    lobe fingerprint; Path-A iteration; FE-direct H1 psi with min-
    seminorm regularisation.
  - Results: complexity tier table (numbers above); Gx chain
    16.24 % (kuijpers) -> 9.29 % (field_aware) -> 8.11 % (+ Path-A);
    planar 0.47 % monotone-converged FE-direct; dead-end catalogue
    + the field_aware positive resolution as supplementary material.
  - Open extensions: cylindrical Gx FE-direct (Hard tier test);
    material kernel via (A) callback; multivalued potential D-path.

References to cite:
  - Kuijpers, Jansen, Lomonova, Compumag 2023 [525] (TU/e EPE).
  - Liu, Hennig, Korvink, IEEE TM 48 (2012) 1179 (H2-smooth SFD).
  - Sugahara Lab, SA-25-020 ((ACA+)+TSVD lineage).
  - HACApK (ppOpen-HPC, MIT, ``src/ext/HACApK/``).
"""


ACA_TSVD_REGULARIZED = r"""
# Regularisation folded into the ACA+TSVD factorisation (2026-05-31)

The constrained min-norm problem

    min   psi^T S psi      s.t.    A psi = B

(S any SPD regularisation matrix on the free DOFs) admits a clean
closed form on top of the ACA+TSVD ``A ~= U Sigma V^T`` (orthonormal V
columns from method 3 recompression):

    psi = (S^-1 V) . W^-1 . Sigma^-1 . U^T B,    W = V^T S^-1 V (k x k).

Shipped as ``radia.stream_function.RegularizedTSVD`` (radia 4.83.0+):

    @dataclass
    class RegularizedTSVD:
        base: StreamTSVD            # ACA+TSVD of A (kernel-agnostic)
        Sinv_V: ndarray             # (N, k) precomputed S^-1 V
        W_inv: ndarray              # (k, k) precomputed (V^T S^-1 V)^-1

        @classmethod
        def from_stiffness(cls, base, S) -> RegularizedTSVD: ...

        def solve(self, B, k_mode=None) -> ndarray:
            return self.Sinv_V @ (self.W_inv @ ((self.base.U.T @ B)
                                                 / self.base.S))

Plus a one-shot helper ``pseudo_inverse_solve_regularized(result, B,
S)`` for non-iterated use.

## Derivation

    A^T          = V Sigma U^T
    S^-1 A^T     = (S^-1 V) Sigma U^T
    A S^-1 A^T   = U Sigma (V^T S^-1 V) Sigma U^T = U Sigma W Sigma U^T
    (A S^-1 A^T) pinv = U Sigma^-1 W^-1 Sigma^-1 U^T   (rank-k Moore-Penrose)
    psi          = S^-1 A^T (A S^-1 A^T) pinv B
                 = (S^-1 V) Sigma U^T . U Sigma^-1 W^-1 Sigma^-1 U^T B
                 = (S^-1 V) . W^-1 . Sigma^-1 . U^T B    (Uses U^T U = I_k)

## Sanity (verified 2026-05-31)

  - S = I: W = V^T V = I_k (because V has orthonormal cols), so
    formula reduces to psi = V Sigma^-1 U^T B == ``pseudo_inverse_solve``.
    Numerical match: 9.96e-16 (machine precision).
  - Random SPD S: matches the direct Lagrangian solution
    psi = S^-1 A^T (A S^-1 A^T)^-1 B to 9.3e-15.
  - Constraint: A psi = B to 5.9e-15.
  - Error scales with aca_eps: at eps=1e-13 the cached psi matches
    direct to 4.9e-9; at eps=1e-10 the cached form has ACA+ noise of
    relative size ~1.5e-3 in the (highly-underdetermined) nullspace
    but seminorm matches to 2.2e-6 (the optimisation target IS
    near-machine-precision unchanged).

## Tikhonov vs TSVD vs generalised-Tikhonov -- ONE factorisation, three filters

The folded form above is the CONSTRAINED (equality, lambda -> 0) limit
min psi^T S psi s.t. A psi = B.  The whole Tikhonov family lives on the
SAME ACA+TSVD ``A ~= U Sigma V^T`` by inserting a spectral FILTER phi_n
on the modes -- the three classic regularisers differ ONLY in phi_n:

    TSVD (truncate at k)        phi_n = 1 (n<=k), else 0          hard cut
    Tikhonov (ridge lambda)     phi_n = sigma_n^2/(sigma_n^2+lambda^2)  smooth
    constrained min-seminorm    phi_n = 1 for every kept mode     (lambda -> 0)

STANDARD (identity) Tikhonov is just S = I:

    psi(lambda) = V diag(sigma_n/(sigma_n^2+lambda^2)) U^T B
                = sum_n [sigma_n^2/(sigma_n^2+lambda^2)] (u_n^T B / sigma_n) v_n

i.e. ``pseudo_inverse_solve`` damped by f = Sigma/(Sigma^2 + lambda^2)
(the api topic gives this one-liner).  GENERALISED (S != I) Tikhonov
min ||A psi - B||^2 + lambda^2 ||L psi||^2 with S = L^T L is classically
the GSVD of the pair (A, L); the folded ``RegularizedTSVD`` is an
ACA-FRIENDLY IMPLICIT GSVD -- it never forms the (dense) GSVD, only
``S^-1 V`` (a SPARSE solve when S is the L2 / H1 stiffness; dense only
for the inductance S) and the k x k ``W = V^T S^-1 V``.  So "ridge alpha"
and "min-energy seminorm S" are the SAME machine: lambda is a smooth
spectral filter, S re-shapes WHICH norm is minimised.  Tikhonov with
lambda ~ sigma_k is the SOFT version of TSVD truncated at rank k.

## Three distinct dials -- do not conflate them (the aca_eps pitfall)

ACA+TSVD-with-Tikhonov has THREE separate truncation / regularisation
dials; mixing them up regularises SILENTLY:

  1. aca_eps -> k_aca   : the ACA+ rank = the APPROXIMATION of A
                          (compression accuracy).  NOT a regulariser;
                          keep it TIGHT.
  2. k_mode (<= k_aca)  : the TSVD regularisation rank (hard filter).
  3. lambda             : the Tikhonov ridge (smooth filter).

PITFALL (the aca_eps sanity above is the evidence): a LOOSE aca_eps lets
the ACA+ truncation ITSELF act as an uncontrolled regulariser -- it drops
small singular directions you may have wanted to keep + regularise
explicitly.  Rule: set aca_eps TIGHTER than the regularisation you intend
(eps ~ 1e-10..1e-13), THEN steer the design with k_mode / lambda / S --
the auditable dials.  At eps=1e-10 the cached psi already carries ~1.5e-3
nullspace noise; the optimisation TARGET (seminorm + constraint) stays
near machine precision, but the specific least-norm REPRESENTATIVE drifts.

## Why the DUAL form (M << N)

SF design is massively UNDER-determined (M target rows << N surface DOFs).
The primal normal equations (A^T A + lambda S) are N x N (huge); the
folded form solves the DUAL ``A S^-1 A^T`` (M x M, small) -- the
constrained-min-seminorm Lagrangian psi = S^-1 A^T (A S^-1 A^T)^-1 B.
Re-solving at a new lambda / B reuses ``S^-1 V`` and ``W^-1`` and redoes
only the k-dimensional filter (O(k^2 + k M)) plus one (N, k) expansion
matvec -- THAT is the ~50 us "k x k core" re-solve the Pareto sweep rides
on.  The regularised inverse lives in the k-dimensional ACA range, never
the N-dimensional DOF space.

## Choosing lambda / k_mode (L-curve, Morozov, and the SF-specific meaning)

  - L-curve: plot (||A psi - B|| , ||psi|| or peak J) over a lambda / k
    sweep; the CORNER (max curvature) is the design point.  ``lcurve()``
    in ``topology_optimization.linear_inverse`` returns the point list.
  - Morozov discrepancy: stop where ||A psi - B|| == the noise level.  SF
    has NO measurement noise -- the floor is the FIELD UNREACHABILITY on
    this winding surface (the residual PLATEAU), so Morozov here == "stop
    at the achievable-on-this-surface floor".
  - GCV (parameter-free): not wired in; use the L-curve corner.
  HONEST: the SF alpha-L-curve is NON-MONOTONIC -- the iso-contour
  TOPOLOGY jumps with alpha (a saddle merges / splits a current region),
  so report BOTH residual AND peak J and pick the corner, never residual
  alone.  Measured: for the Gx single-stroke plain TSVD mode-truncation
  beat the best ridge alpha (8.45 % vs 8.97 %) -- the best regulariser is
  target-dependent, so the menu (tsvd / tikhonov / h1 / inductance) earns
  its keep.

## Conditioning -- the constant-psi null space

Both A and a GRADIENT seminorm S = L^T L annihilate an additive CONSTANT
in psi: grad(const) = 0, and a uniform current potential carries no
current (K = n x grad psi = 0).  So S is SPD only on the quotient and
``W = V^T S^-1 V`` can blow up.  Two fixes, both already in the code:
  (a) a tiny ridge ``S += (1e-9 * mean_diag) I`` (``_inductance_seminorm``
      ridge=1e-9) lifts the harmless null space;
  (b) ``--confine abe / on`` GROUNDS it via the reduction matrix R (one
      free constant per physical boundary component).
Without either, W^-1 amplifies the constant mode -- a DC offset in psi
that the equal-increment contouring happens to ignore, but which wrecks
any ||psi||-based L-curve selection.  This is the regularisation-side
reason the Verify-First FES checks (slaved DOFs, constant-mode grounding)
matter before trusting an L-curve.

## Why this matters

  1. ACA+ amortisation -- the heavy O(k * (M + N)) ACA+ runs ONCE,
     regardless of how many B vectors / how many S choices / how many
     Path-A iterations follow.
  2. Path-A inner solve becomes O(k * (M + N)) -- a single (k x N)
     matvec, one k x k inverse-apply, one (M, k) matvec.  Measured
     0.04 ms / iter at N=1773, k=25; vs ~10 ms / iter for the direct
     H1 path-solve.  100-iter Path-A wall time drops 50s -> ~1s.
  3. Regularisation sweep is cheap -- once ACA+ is paid, scanning L2 /
     H1 / sigma-weighted H1 / inductance-diagonal / Lawson-IRLS L_inf
     each costs only one ``Sinv_V`` build + one k x k inverse.
  4. Material-kernel lift -- when each A(i, j) is a Radia ``rad.Solve()``
     + ``rad.Fld()`` (= MMM iron yoke, shielded coil), the ACA+
     amortisation is the dominant cost; the regularisation-folded
     form means all design choices share that cost.

## Path-A cached integration (production pattern)

```python
from radia.stream_function import aca_tsvd, RegularizedTSVD

# Once per design problem
res = aca_tsvd(M, N, entry, modes=M, aca_eps=1e-10)
reg = RegularizedTSVD.from_stiffness(res, S_free)
psi_initial = reg.solve(B_target)

# In the Path-A outer loop:
for it in range(n_iter):
    # ... compute residual from current chain ...
    delta_psi = reg.solve(residual)             # 0.04 ms each
    psi = psi + step * delta_psi
    # ... recompute chain from psi ...
```

The base ``res`` and the cache ``reg`` survive the loop; the per-iter
work shrinks to a few matvecs.

## IRLS for L_inf on top of the cache

``solve_linf_irls`` in
``examples/stream_function/demo_planar_uniform_fem_psi_advanced.py``
rebuilds S^(k) each IRLS iter (vertex weights that depend on
|grad psi^(k-1)|) but ALWAYS folds into the SAME ACA+TSVD base:

    res = aca_tsvd(...)                # ONCE
    for k in range(n_iter_irls):
        w_v = bounded_active_set_weight(|grad psi^(k-1)|_v)
        S_k = assemble_weighted_h1_stiffness(w_v)
        reg_k = RegularizedTSVD.from_stiffness(res, S_k)
        psi_k = reg_k.solve(B)

Bounded-ratio (``weight_ratio_max=20``) + damped (``damping=0.5``) +
active-set boost (``hot_quantile=0.85, hot_boost=3``) keeps the system
well-conditioned where pure Lawson IRLS oscillates.  Optional final
scalar rescale enforces ``max_v |grad psi| <= jmax`` exactly.

On the planar uniform-Bz benchmark the H1 baseline already has
peak/mean gradient ~= 1.47 (near-L_inf-uniform), so IRLS gives 10-30%
peak reduction.  On harder geometries with hot spots (Gx fingerprint
chains, surface-deformed coils) the framework has more bite without
re-factorising A.

## Demo

``examples/stream_function/demo_regularized_aca.py`` exercises all 5
modes routed through the cached form (l2 / h1 / h1_sigma /
inductance_diag / linf_irls) and reports per-mode fold time, solve
time, continuous residual, single-stroke chain RMS, and best-fit
single current I_w on the planar uniform-Bz benchmark.  All four
closed-form modes hit ||A psi - B|| / ||B|| < 1e-14 (machine-precision
constraint satisfaction); IRLS does not because its jmax rescale is
binding (RMS 1.12 % is the field-fit cost of the L_inf cap).

## Constrained reg-aware Optuna loop -- "min reg s.t. RMS <= eps"

The cached ``RegularizedTSVD`` lifts naturally into the standard
MRI-shim / gradient-coil bilevel formulation: minimise the
regularisation norm (= inductance / dissipation / surface current L2,
depending on which S) subject to a field-uniformity tolerance.

Shipped in ``demo_planar_uniform_fem_psi_advanced.py`` with two new
CLI flags:

    --minimize-reg          # switch the Optuna objective from RMS to reg_norm
    --eps-rms 0.025         # RMS tolerance (default 2 %)
    --reg-penalty 100       # quadratic-penalty weight (auto-scaled to reg_norm)

Implementation: scalar single-objective CMA-ES (same sampler as the
RMS-min loop) with penalty form

    cost = reg_norm + reg_scale * penalty_weight * max(0, RMS/eps - 1)^2

where ``reg_scale`` is set from the FIRST trial's reg_norm so units
match.  ``best_feasible`` (the lowest reg_norm among trials satisfying
the constraint) is tracked separately on ``study.user_attrs``.

### Measured Pareto trade-off (planar uniform Bz, 8 trials, eps=2.5 %)

| Mode                          | RMS    | psi^T S psi (H1)  | Note                          |
|-------------------------------|--------|--------------------|-------------------------------|
| Flat baseline (no deform)     | 2.09 % | 2.59e+06           | reference                     |
| min(RMS) (CMA-ES)             | 0.77 % | 2.74e+06           | accuracy-opt; energy UP from baseline |
| min(reg) s.t. RMS <= 2.5 %    | 1.82 % | 2.15e+06           | meets spec, 28 % LOWER energy vs min(RMS) |

The accuracy-optimal point uses MORE energy than the flat baseline
because squeezing RMS from 2 % to 0.8 % requires localised current
concentrations.  The reg-min constrained point trades some
(in-spec) accuracy for ~25-30 % lower surface current L2 -- the
classical MRI-shim Pareto answer.

### Cache reuse

Because ``A`` changes per trial when the surface deforms, the cached
``RegularizedTSVD`` is rebuilt per trial.  What IS amortised:
``_compute_reg_norm`` builds ``S`` once per trial (~50 ms) then a
single matvec.

## Optimising the regularisation SHAPE (fixed surface, ACA+ reused)

The dual of the deformation loop: keep the surface FIXED, optimise
``sigma(x, y)`` (where current flows) in the 1/sigma-weighted H1
seminorm.  Now ``A`` is CONSTANT, so the ACA+TSVD base is computed
ONCE and reused across every Optuna trial -- only the cheap fold
``S^-1 V + W = V^T S^-1 V`` is rebuilt per trial.  This is the
genuine demonstration of the folded form's amortisation (the
deformation loop cannot show it because A changes there).

Shipped in ``examples/stream_function/demo_reg_hyperparam_aca.py``:
CMA-ES over a Gaussian conductivity feature

    sigma(x, y) = 1 + amp * exp(-((x-cx)^2 + (y-cy)^2)/width)

inner solve ``min INT (1/sigma)|grad psi|^2 s.t. A psi = B`` folded
through the cached base.

### Result (planar uniform Bz, 30 trials, 15 s)

| Quantity                | Uniform sigma (plain H1) | Optimised sigma |
|-------------------------|--------------------------|-----------------|
| single-stroke chain RMS | 2.09 %                   | **0.73 %**      |
| wire length             | 18.98 m                  | 18.91 m         |

Optimiser finds amp ~= -0.81 (an ~80 % conductivity DIP at centre)
which pushes current toward the periphery -- the Helmholtz-like
distribution that produces a more uniform central Bz.  2.9x RMS drop
with NO geometry change, NO new target, purely by reshaping the
regularisation penalty.  Comparable to the surface deformation loop
(0.77 %) via an orthogonal, physically cheaper DOF.

### Timing breakdown (proves the amortisation)

    ACA+TSVD base (U, Sigma, V):  ONCE, ~33 ms
    per trial:  build S(sigma) -> S^-1 V -> W -> solve  (~43 ms fold)
                ACA+ base NOT recomputed

At FE scale the ~33 ms base is small, but with a MATERIAL kernel
(rad.Solve()+rad.Fld() per entry, ~100x slower) the ACA+ base would
dominate and reusing it across the 30-trial sweep is the entire
performance argument.

## Multi-objective NSGA-II Pareto front -- ``--pareto``

The penalty form gives ONE point.  To trace the full Pareto curve
of ``(RMS, psi^T S psi)`` in a single Optuna run, switch to
multi-objective NSGA-II:

    --pareto                # requires --deform; mutually exclusive with --minimize-reg

Implementation:

  - ``optuna.samplers.NSGAIISampler(seed=42, population_size=auto)``
  - ``optuna.create_study(directions=["minimize", "minimize"])``
  - ``objective(trial) -> (rms, reg_norm)`` tuple
  - ``study.best_trials`` = non-dominated set; demo deduplicates
    (NSGA-II re-evaluating same params across generations is common)
    and sorts by RMS

Outputs:

  - ``demo_pareto_results.json`` -- all trials + Pareto front
  - ``demo_pareto_plot.png`` -- scatter + Pareto curve

### Measured Pareto front (planar uniform Bz, 50 trials, H1, bump)

| RMS    | psi^T S psi  | bump (amp, cx, cy)        | Note         |
|--------|--------------|----------------------------|--------------|
| 0.58 % | 2.80e+06     | (-0.032, -0.059, +0.007)   | accuracy end |
| 0.61 % | 2.46e+06     | (+0.021, -0.144, +0.141)   | **knee**     |
| 1.04 % | 2.29e+06     | (+0.045, +0.140, -0.121)   |              |
| 1.39 % | 2.27e+06     | (+0.045, -0.059, +0.146)   |              |
| 1.67 % | 2.27e+06     | (+0.045, -0.062, +0.146)   |              |
| 2.13 % | 2.20e+06     | (+0.044, +0.099, +0.081)   | energy end   |

Wall time: **33 s** (50 inner SF solves).  Sharp KNEE at ~0.6 % RMS --
12 % reg drop for only 0.03 % accuracy loss.  Past the knee the curve
flattens.  Plain ``min(RMS)`` CMA-ES collapses to (0.77 %, 2.74e+06)
on the same problem; NSGA-II's accuracy end is BETTER on field but
slightly worse on energy -- expected behaviour, scalar minimisation
finds one specific point while NSGA-II explores the front.

### Mode selection table

| Situation                                              | Use            |
|--------------------------------------------------------|----------------|
| Spec known up-front ("RMS must be <= 2 %")             | --minimize-reg --eps-rms 0.02 |
| Choosing the spec (operating point unknown)            | --pareto       |
| Single number for a DoE                                | --minimize-reg |
| Pareto-front figure for a paper                        | --pareto       |
| Budget < 15 trials                                     | --minimize-reg |
| Budget >= 30 trials                                    | either; --pareto gives more info/trial |

## References / cross-links

  - ``docs/stream_function/regularization.md`` -- user-facing docs
  - ``docs/stream_function/deformation.md`` -- constrained reg-aware
    loop section + Pareto trade-off table
  - ``radia.stream_function.RegularizedTSVD`` -- implementation
  - ``streamfunction(performance)`` -- ACA+ amortisation numbers
  - ``streamfunction(session_2026_05_30)`` -- Section 10 (regularisation
    sweep + deformation) is the prior, non-folded form
  - Hansen 1998 "Rank-Deficient and Discrete Ill-Posed Problems"
    SIAM -- standard reference for TSVD + Tikhonov; the specific
    Path-A x cached folded form is the SF-coil-design-specific
    composition implemented here.


================================================================
PARETO FRONT (field homogeneity vs peak current density) + the
SHEET-METAL (板金) lever   [2026-06, radia 4.89.1+]
================================================================

1. FOLDED TIKHONOV (the "+ alpha I" core).  RegularizedTSVD.solve(B,
   alpha=0.0) was extended: for alpha > 0 the SAME cached ACA+TSVD
   factorisation gains a single alpha*I inside the k x k core,

       psi(alpha) = S^-1 V . (alpha I + Sigma^2 W)^-1 . Sigma . U^T B
                  == (A^T A + alpha S)^-1 A^T B,   W = V^T S^-1 V.

   This UNIFIES, on one factorisation: TSVD (hard truncation via
   k_mode), Tikhonov (smooth filter sigma/(sigma^2+alpha); S=I gives the
   classic filter factors to machine precision), and the generalised
   seminorm S.  alpha=0 recovers the exact-fit min-seminorm form.  An
   entire L-curve / Pareto alpha-sweep re-solves only the small core
   (~50 us/point), no re-factorisation.  Verified vs dense
   (A^T A + alpha S)^-1 A^T B (machine precision for S=I; ACA tolerance
   for general S).  4 golden tests in tests/test_stream_function.py.

2. PARETO FRONT.  Objectives f1 = ||A psi - B||/||B|| (homogeneity),
   f2 = max_v |grad psi| (peak surface current density = ||K||_inf).
   Tikhonov is the weighted-sum scalarisation of the convex bi-objective
   (misfit, seminorm), so the alpha-sweep traces the WHOLE convex Pareto
   front.  demo_pareto_tikhonov_aca.py: --front {energy, peak, both}.

3. FOUR STACKABLE LEVERS to push the front (planar gradient hot-spot,
   order-3 psi):
     1. Tikhonov alpha  -- moves ALONG the front.
     2. L-inf seminorm (bounded-ratio IRLS) -- redistributes current
        within a FIXED A: -18% peak (median) at matched homogeneity.
     3. geometry -- changes A: planar former size -34% (monotone,
        diminishing returns); cylinder LENGTH has an OPTIMUM
        (~50 cm, -37%; longer stops helping once it covers the DSV).
     4. sheet-metal (板金) deformation -- FORM the surface (see 4).
   NSGA-II over (geometry, alpha) traces the 3-objective
   (homogeneity, peak, former-size) surface: low peak REQUIRES a larger
   former.  demo_pareto_geometry_nsga.py, demo_pareto_cylinder.py.

4. SHEET-METAL (板金) -- the lever DIRECTION is geometry-dependent, and
   the honest claim requires separating STANDOFF from BENDING:
     * STANDOFF = forming the sheet CLOSER (smaller average gap): large
       but achievable by repositioning a FLAT former -> NOT genuine.
     * BENDING  = reshaping at FIXED average standoff -> genuine forming.
   - PLANAR: out-of-surface bending z=f(x,y) is the lever.  Genuine
     (zero-mean, fixed avg standoff, |z|<=5cm) pushes the whole front
     -5..-18% (best -17% near exact homogeneity, per-alpha CMA-ES,
     warm-started).  Standoff-allowed reaches -53% but is
     standoff-dominated.  Distance-dominated field => bending bounded.
     demo_pareto_deform.py (--zero-mean default / --allow-standoff).
   - CYLINDER: the lever FLIPS.  Out-of-surface (radial) forming is WEAK
     (~-3%, optimiser barely uses the |dr| budget); the IN-SURFACE axial
     bending is DOMINANT.  Length-preserving axial reparametrisation
     Z(z)=z+sum_k b_k sin(k pi (z+L/2)/L) redistributes loop spacing at
     FIXED radius => 100% genuine forming, NO standoff component (no
     zero-mean trick needed).  Whole Gx-fingerprint front -10..-25%
     (best -25%).  Done RIGHT: local-spacing peak (non-uniform axial
     derivative) + spacing-WEIGHTED graph Laplacian seminorm (z-edge
     a*dphi/dZ, phi-edge dZ/(a*dphi)) + monotone dZ/dz>0 penalty.
     demo_pareto_cylinder_deform.py.
   - 2D in-surface (axial Z(z) + azimuthal Phi(phi)): the AZIMUTHAL
     component matters ONLY when the target has azimuthal structure
     (m>=2).  Per-target exact-homogeneity test (axial-only vs 2D
     axial+azimuthal, azimuthal benefit = 2D vs axial):
       Gx        (m=1)  +0.0%  (smooth cos phi, no azimuthal hot spot)
       Z2        (m=0)  -1.8%  (axisymmetric; axial dominates at -20%)
       S2=xy     (m=2)  -1.2%
       C2=x^2-y^2(m=2)  -5.4%  (axial-only is only -5.0% here, so 2D
                                ~DOUBLES the in-surface benefit to ~-10%)
     So the in-surface lever DIRECTION follows the target symmetry:
     low-m / axisymmetric -> AXIAL; high-m azimuthal targets (C2 ellipse,
     S2) -> AZIMUTHAL is a comparable lever.  Axial-only
     (demo_pareto_cylinder_deform.py) is the right default for Gx/Z2; 2D
     azimuthal pays off for m>=2 shims.  (Earlier "azimuthal marginal"
     was Gx-specific.)
   Same plane-vs-cylinder reversal as the single-stroke distortion lever
   (out-of-surface lift on the plane, in-surface reroute on the
   cylinder) -- see topic "single_stroke".

DEMOS (examples/stream_function/):
  demo_pareto_tikhonov_aca.py     -- folded-Tikhonov (homogeneity, peak) front
  demo_pareto_geometry_nsga.py    -- geometry lever + NSGA-II joint front
  demo_pareto_cylinder.py         -- cylinder Gx, length lever (optimum)
  demo_pareto_deform.py           -- planar sheet-metal (zero-mean genuine)
  demo_pareto_cylinder_deform.py  -- cylinder in-surface axial bending (genuine)
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
        "session_2026_05_30": ACA_TSVD_SESSION_2026_05_30,
        "regularized": ACA_TSVD_REGULARIZED,
    }
    aliases = {
        # Pareto front + folded Tikhonov + sheet-metal (板金) lever
        "pareto": "regularized", "pareto_front": "regularized",
        "pareto-front": "regularized", "front": "regularized",
        "peak_current_density": "regularized", "peak_current": "regularized",
        "peak_j": "regularized", "peak": "regularized", "homogeneity": "regularized",
        "tikhonov": "regularized", "tikhonov_alpha": "regularized",
        "alpha_sweep": "regularized", "l_curve": "regularized",
        "lcurve": "regularized", "l-curve": "regularized",
        "filter_factors": "regularized", "geometry_lever": "regularized",
        "former_size": "regularized", "cylinder_length": "regularized",
        "nsga": "regularized", "nsga2": "regularized", "nsga-ii": "regularized",
        "joint_front": "regularized", "linf": "regularized", "l_inf": "regularized",
        "minimax": "regularized", "irls": "regularized",
        "pareto_bankin": "regularized", "bankin_pareto": "regularized",
        "sheet_metal_pareto": "regularized", "surface_forming": "regularized",
        "in_surface": "regularized", "in-surface": "regularized",
        "axial_bending": "regularized", "radial_forming": "regularized",
        "lever_flip": "regularized", "deform_pareto": "regularized",
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
        "field_aware": "single_stroke", "field-aware": "single_stroke",
        "nn_blend": "single_stroke", "chain": "single_stroke",
        "rung": "single_stroke", "azimuthal": "single_stroke",
        "distort": "single_stroke", "distortion": "single_stroke",
        "sheet_metal": "single_stroke", "sheet-metal": "single_stroke",
        "bankin": "single_stroke", "wire_distortion": "single_stroke",
        "coil_distortion": "single_stroke", "deform_wire": "single_stroke",
        "mesh_deformation": "single_stroke", "3d_print": "single_stroke",
        "sphere": "single_stroke", "arbitrary_surface": "single_stroke",
        "arbitrary-surface": "single_stroke", "curved_former": "single_stroke",
        "former": "single_stroke", "fe_direct": "single_stroke",
        "fe-direct": "single_stroke", "surface_fem": "single_stroke",
        "nmr_shim": "single_stroke", "shim_coil": "single_stroke",
        "conformal": "single_stroke", "z2": "single_stroke",
        # 2026-05-30 session summary
        "session": "session_2026_05_30",
        "summary": "session_2026_05_30",
        "today": "session_2026_05_30",
        "2026_05_30": "session_2026_05_30",
        "2026-05-30": "session_2026_05_30",
        # Specific session sub-topics route to the session summary
        "complexity_tier": "session_2026_05_30",
        "tier": "session_2026_05_30",
        "chain_methods": "session_2026_05_30",
        "path_a": "session_2026_05_30",
        "path-a": "session_2026_05_30",
        "compensated": "session_2026_05_30",
        "fe_direct": "session_2026_05_30",
        "fe-direct": "session_2026_05_30",
        "fem_psi": "session_2026_05_30",
        "dead_ends": "session_2026_05_30",
        "ngsbem": "session_2026_05_30",
        "ngsolve_bem": "session_2026_05_30",
        "ngsolve.bem": "session_2026_05_30",
        "2604": "session_2026_05_30",
        "htool": "session_2026_05_30",
        "calcsubmatrix": "session_2026_05_30",
        "fmm_vs_aca": "session_2026_05_30",
        "callback_contract": "session_2026_05_30",
        "planar_uniform": "session_2026_05_30",
        # 2026-05-31: regularisation-folded ACA+TSVD closed form
        "regularised": "regularized",
        "regularization": "regularized",
        "regularisation": "regularized",
        "reg_tsvd": "regularized",
        "regularizedtsvd": "regularized",
        "regularized_tsvd": "regularized",
        "folded": "regularized",
        "closed_form": "regularized",
        "closed-form": "regularized",
        "lawson": "regularized",
        "irls": "regularized",
        "linf_irls": "regularized",
        "h1_cache": "regularized",
        "sinv_v": "regularized",
        "w_inv": "regularized",
        "demo_regularized_aca": "regularized",
        "path_a_cache": "regularized",
        # Constrained reg-aware Optuna loop (2026-05-31)
        "minimize_reg": "regularized",
        "minimise_reg": "regularized",
        "min_reg": "regularized",
        "eps_rms": "regularized",
        "pareto": "regularized",
        "constrained": "regularized",
        "mri_shim": "regularized",
        "shim": "regularized",
        "gradient_coil": "regularized",
        "reg_aware": "regularized",
        "reg-aware": "regularized",
        # NSGA-II multi-objective Pareto (2026-05-31)
        "nsga": "regularized",
        "nsga2": "regularized",
        "nsga-ii": "regularized",
        "nsgaii": "regularized",
        "multiobjective": "regularized",
        "multi_objective": "regularized",
        "multi-objective": "regularized",
        "pareto_front": "regularized",
        "pareto-front": "regularized",
        "trade_off": "regularized",
        "tradeoff": "regularized",
        "trade-off": "regularized",
        # sigma-shape reg-hyperparameter sweep (fixed surface, ACA+ reused)
        "sigma_shape": "regularized",
        "sigma-shape": "regularized",
        "reg_hyperparam": "regularized",
        "reg_hyperparameter": "regularized",
        "hyperparameter": "regularized",
        "hyperparam": "regularized",
        "sigma_optuna": "regularized",
        "fixed_surface": "regularized",
        "demo_reg_hyperparam_aca": "regularized",
        "aca_reuse": "regularized",
        "aca_amortisation": "regularized",
        "aca_amortization": "regularized",
    }
    t = aliases.get(t, t)
    if t == "all":
        return "\n\n".join([ACA_TSVD_OVERVIEW, ACA_TSVD_METHOD, ACA_TSVD_API,
                            ACA_TSVD_KERNEL_AGNOSTIC, ACA_TSVD_PERFORMANCE,
                            ACA_TSVD_CMAES, ACA_TSVD_VALIDATION,
                            ACA_TSVD_LITERATURE, ACA_TSVD_WORKFLOW,
                            ACA_TSVD_SINGLE_STROKE,
                            ACA_TSVD_SESSION_2026_05_30,
                            ACA_TSVD_REGULARIZED])
    if t in table:
        return table[t]
    return (f"Unknown topic '{topic}'.  Available: overview, method, api, "
            "kernel_agnostic, performance, cmaes, validation, literature, "
            "workflow, single_stroke, session_2026_05_30, regularized, all.\n\n"
            + ACA_TSVD_OVERVIEW)
