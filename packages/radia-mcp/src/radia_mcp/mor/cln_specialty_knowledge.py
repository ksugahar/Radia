"""CLN lab-specialty methods: termination, Hiruma, BEM+FEM, Nagamine error.

Covers public-safe curated corpus, 12_比留間法, 13_BEM+FEM,
15_誤差論@長嶺氏}. These are the lab-internal techniques not yet
fully published in the canonical Kameari-Ebrahimi-Sugahara-Shindo-
Matsuo 2018 IEEE TMAG paper or its standard extensions.

The four specialty topics, in order of increasing scope:

  08 / termination      - what to put at the end of a finite ladder
  12 / Hiruma method    - non-symmetric Lanczos (Hokkaido / Igarashi lab)
                          giving CLN directly from (G + s C) form
  13 / BEM + FEM        - coupling Kameari-CLN on the conductor with a
                          BEM exterior, multi-mode multi-port ladder
  15 / Nagamine error   - rigorous a-posteriori bound on truncation
                          and on finite-element discretisation of L_n, R_n

Companion module: radia_mcp.mor.cln_knowledge (the canonical Kameari-
Ebrahimi-Sugahara-Shindo-Matsuo recursion, multi-expansion-point
extension, applications).  This file is the *internal* layer.
"""

# ---------------------------------------------------------------------------
# OVERVIEW
# ---------------------------------------------------------------------------

OVERVIEW = """
# CLN lab-specialty methods (Group D)

The canonical CLN paper (Kameari-Ebrahimi-Sugahara-Shindo-Matsuo, IEEE
TMAG 54(3):7201804, 2018) defines a recursion that produces a ladder
of (R, L) stages.  Four research threads inside Sugahara Lab refine or
extend that recursion in directions that have NOT (yet) become part
of the canonical literature:

| Folder | Topic | Lead | Status |
|--------|-------|------|--------|
| 08_終端問題 | Ladder TERMINATION impedance | Sugahara, 2016-2018 | internal note |
| 12_比留間法 | Hiruma's non-symmetric Lanczos -> CLN | Hiruma (Hokkaido / Igarashi) | published abroad, internal MATLAB ports |
| 13_BEM+FEM | Coupling open-boundary BEM with conductor-CLN | Sugahara, 2018 | multi-mode multi-port; ongoing |
| 15_誤差論@長嶺氏 | Verified error bound on R_n, L_n | Nagamine (collaborator) | 4 internal PDFs, EMC 2022 paper |

Why these are "specialty" and not "core":

  * 08 (termination) - The canonical paper computes Z(s) by direct
    continued-fraction evaluation.  At a finite ladder length N, the
    response loses accuracy at high frequency.  Replacing the open or
    short termination with a frequency-dependent Z_term(s) recovers
    accuracy without enlarging N.  Sugahara's 2018-03-17 note shows
    the construction and the residual error.

  * 12 (Hiruma) - The canonical CLN recursion needs the differential
    operators (curl curl, etc.) of the original PDE.  Hiruma derived
    a non-symmetric Lanczos algorithm that ingests the ALGEBRAIC
    system (G + s C) x = b directly and produces the same Cauer
    ladder constants.  This makes CLN reachable from any black-box
    FEM matrix exporter.  It also unifies CLN with PVL / SyPVL /
    PRIMA (the PEEC MOR family) on the algebraic side.

  * 13 (BEM+FEM) - The canonical CLN is interior-domain only.
    Coupling with an exterior BEM (eddy current in a finite
    conductor, free space outside) requires changing the conductor
    boundary condition from Dirichlet to Robin and re-deriving the
    Kameari recursion.  Result: a MULTI-PORT, MULTI-MODE ladder
    where each H_m mode contributes one port.

  * 15 (Nagamine) - The canonical paper proves CONVERGENCE of the
    ladder as N -> infinity but does not bound the FE discretisation
    error in any individual R_n, L_n.  Nagamine (visiting Matsuo Lab
    at Kyoto) closes this gap with a verified-arithmetic computation
    of the exact R_n, L_n for the square cross-section, and a
    cut-off-frequency / characteristic-length diagnostic for mesh
    adequacy.

All four threads share the same starting point: the basic Kameari
recursion in cln_knowledge.py.  They diverge in WHERE they cut into
that recursion (termination at the tail, algebraic ingestion at the
front, BEM coupling at the side, error bound underneath).
"""

# ---------------------------------------------------------------------------
# 08: TERMINATION
# ---------------------------------------------------------------------------

TERMINATION_THEORY = """
# Ladder TERMINATION (08_終端問題)

## The problem

The Kameari recursion at stage N has computed R_0, L_1, R_2, ..., L_{2N-1}
(or up to R_{2N}, depending on whether the ladder ends on R or L).  The
truncated impedance is then evaluated by

   Z_N(s) = R_0 + 1 / (s L_1 + 1 / (R_2 + 1 / (s L_3 + ...))).

For a 1-D magnetic sheet in the H-mode (Shindo-Kameari analytic case),
the exact impedance is

   Z(s) = (1/sigma d) * tanh(sqrt(s sigma mu) d / 2)
                       / (sqrt(s sigma mu) d / 2)

and the truncated Z_N(s) ALTERNATES between under- and over-shoot of
the exact value as N grows (the "tail-flapping" or shippo-furi pattern,
in Sugahara's 2016-12-23 note).

## The fix: terminate with a frequency-dependent impedance

Instead of stopping the ladder open (Z_end = infinity) or shorted
(Z_end = 0), insert a small impedance Z_term(s) at the tail.  If
Z_term(s) is chosen to approximate the rest of the ladder (the
infinite tail from stage N+1 onward), then Z_N(s) becomes accurate
over a much wider frequency range.

Sugahara 2018-03-17 internal note proposes:

   Z_term^E(s) = R_{2N} + s L_{2N+1}    (E-mode, R-then-L terminator)
   Z_term^H(s) = s L_{2N+1} + R_{2N+2}  (H-mode, L-then-R terminator)

where R_{2N+2}, L_{2N+3} are taken from the ANALYTIC asymptotic form
of R_n, L_n (here, R_n = sigma d / (2n+1), L_n = mu d / (2(2n+1))
for the magnetic sheet) -- they are NOT computed by FEM.

## Empirical findings (2016-12-23 note)

For the round-trip 2-conductor model submitted to Compumag 2017:

  * Without termination, Z_N requires N >= 30 for 1 % accuracy over
    1 kHz to 1 MHz.

  * With ideal R-termination, Z_6 has the same accuracy as Z_{N>15}
    in the mid-frequency band; at high frequency the high-stage Z_N
    still wins because the truncated termination has its own error.

  * Conclusion: termination IMPROVES the lower-mid band but is NOT a
    substitute for adding more stages at high frequency.  In Sugahara's
    own words (2016-12-23):

      "Adding stages is more effective than optimising the terminator."

The note recommends terminator-optimisation as a POLISH after fixing
N, not as a way to reduce N.

## Open question (2018 note, section 3)

Is there a closed-form Z_term(s) that makes Z_N(s) exact for s in some
finite band [s_min, s_max]?  The infinite-ladder limit is known:
Z_term^infinity(s) = (1/sigma d) tanh(...)/(...) -- but this requires
the analytic solution, which CLN is trying to avoid.  Looking for a
RATIONAL approximation of Z_term^infinity matched to N is an open
direction.
"""

TERMINATION_CODE = """
# Termination: representative code (08_終端問題)

## Driver: main.m

The 2018-03-13 driver compares the ladder with and without termination
optimisation.  Skeleton:

```matlab
% main.m   (public-safe curated corpus)
global y z N
f   = logspace(3, 7, 200);
sig = 5.8e7;  mu = 4*pi*1e-7;  d = 0.001;
z   = sqrt(-1j*2*pi*f*sig*mu) * d;
Z0  = tan(z) ./ z;             % analytic reference

for N = 1 : 5
  % case 1 -- no termination
  err = cf([0, 0]);            % cf returns y in global, evaluating
                               % the truncated continued fraction
  Z1  = y;
  % case 2 -- optimised termination
  for nf = 1 : length(f)
    z = sqrt(-1j*2*pi*f(nf)*sig*mu) * d;
    p = fminsearch('cf', [0, 0]);   % search for best (Re, Im) of Z_term
    Z2(nf) = y;
  end
end
```

## Helper: cf.m -- the truncated continued fraction WITH a terminator

```matlab
function err = cf(p)
  global y z N
  y = z.^2 / ((2*N + 1) + complex(p(1), p(2)));   % terminator
  for n = N : -1 : 1
    nn = 2*n - 1;
    y  = z.^2 ./ (nn - y);                        % back-up the ladder
  end
  y   = y ./ z.^2;
  err = sum(abs(y - tan(z) ./ z) .^ 2);           % L2 vs analytic
end
```

The terminator slot is `(2*N + 1) + complex(p(1), p(2))`.  When `p=0`
this is the open-termination ladder of length N; when `p` is optimised
by fminsearch the ladder matches the analytic tan(z)/z over a wider
band of `z = sqrt(-j omega sigma mu) d`.

## Termination form: termination.m

```matlab
% termination.m (analytic study; no FEM, symbolic ladder)
for N = 2
  switch nCASE
    case 1   % polynomial terminator with cubic / quartic corrections
      y = z.^2 / (2*N + 1) + 0.00002 * z.^3 + 0.002 * z.^4 ...
                            - 0.0003 * z.^5;
    case 2   % bare 1 / (2N + 1) terminator
      y = z.^2 / (2*N + 1);
  end
  for n = N : -1 : 1
    nn = 2*n - 1;  y = z.^2 ./ (nn - y);
  end
end
```

Case 1 (polynomial-corrected terminator) gives ~3 orders of magnitude
lower error than case 2 over the band z in (0, pi/2), at the cost of
4 extra coefficients per terminator.

## Hiruma-style infinite-ladder test: hiruma_infinity.m

The same folder contains a driver that builds the FULL ladder length
2N+2 from a synthetic symmetric (G, C) pair and compares the ladder
expansion against the direct solve.  This is the algebraic prototype
that is generalised in 12_比留間法 below.
"""

# ---------------------------------------------------------------------------
# 12: HIRUMA METHOD
# ---------------------------------------------------------------------------

HIRUMA_METHOD = """
# Hiruma method (12_比留間法)

## Who is Hiruma

Shingo Hiruma (比留間 真悟) was a Ph.D. student in Hajime Igarashi's
laboratory at Hokkaido University, Graduate School of Information
Science and Technology.  He is now an Assistant Professor at Kyoto
University, Matsuo Lab (same Matsuo as the canonical CLN paper),
which is what brings his work directly into the Sugahara Lab CLN
catalogue.

His committee-paper series in the Japanese committee meeting
(電気学会 静止器・回転機合同研究会) is the origin of the Hiruma
method as Sugahara Lab uses the term:

  * #9-04   "対称 Lanczos アルゴリズムについて"   (symmetric Lanczos)
  * #11-08  "(...) の Lanczos アルゴリズムについて その2"  (extension)
  * clnb_0309 (Hiruma's CLN follow-up note)

## The core idea

Canonical Kameari CLN starts from the PDE (curl(nu curl A) + s sigma A
= J_0) and builds basis fields H_n, E_n by recursive static FEM
solves.  Each step needs the differential operators of the PDE.

Hiruma starts from the discretised algebraic system

   (G + s C) x = B u                                          (1)
   y = c^T x                                                  (2)

with G the FE conductance matrix and C the FE capacitance/mass matrix
(or, in the magnetic case, K + s N as Hiruma writes it: K the stiffness,
N the conductivity mass).  The transfer function is

   H(s) = c^T (G + s C)^{-1} B u
        = l^T (I - s A)^{-1} r                                (3)

where  A = -G^{-1} C,  l = c,  r = G^{-1} B u.

Hiruma's NON-SYMMETRIC LANCZOS expands H(s) in a continued fraction
directly from (3) without ever touching the PDE.

## The recursion (Hiruma, scalar input/output)

Let w_1 = G^{-1} l = G^{-1} r (initial vector).  Set lambda_1 = w_1^T G w_1.
Then for n = 1, 2, ...:

   w_{n+1}^new = w_{n-1} - (A w_n) / lambda_n     if n odd  (or even - mirror)
   lambda_{n+1} = w_{n+1}^T M_{n+1} w_{n+1}       (M alternates G and C)

After N steps the recursion gives lambda_1, lambda_2, ..., lambda_{2N+1}
and these are EXACTLY the ladder constants R_0 = 1/lambda_1,
L_1 = lambda_2, R_2 = 1/lambda_3, L_3 = lambda_4, ...

The non-symmetric form (l != r) handles input-different-from-output
transfer functions; the same Krylov subspace K(A, u_0) is generated
twice as K(A, u_0), K(A^T, v_0) but Hiruma SHORTENS this back to one
Krylov subspace K(A, u_0) when G is SPD by using G-orthogonality:

   (u_i, v_j)_G = u_i^T G v_j = 0   for  0 <= i < j <= N-1.

This is the key insight: the symmetric inner product (.,.)_G replaces
the two-sided non-symmetric inner product of vanilla non-symmetric
Lanczos.  The result is equivalent to SyPVL (Symmetric Pade via
Lanczos) but generates the CAUER (continued fraction) form instead of
the FOSTER (partial fraction) form of PVL.

## What Hiruma's method gains over canonical CLN

  1. Black-box FEM: any FEM solver that exports (G, C) works.  No
     need to expose H_n, E_n basis fields or curl operators.

  2. Unifies CLN with the PEEC MOR family (PVL, SyPVL, PRIMA, Arnoldi).
     Hiruma showed (paper clnb_0309) that the Kameari basis fields are
     a Lanczos sequence in disguise.  This is the bridge between the
     FEM-CLN literature and the PEEC-MOR literature.

  3. Multi-port extension is immediate: if b is replaced by a matrix
     [b1, b2, ...], the recursion runs block-wise and gives a multi-
     port Cauer network.  See hiruma_matrix.m.

  4. Eddy-current loss in a coil layer (sell-mesh / Yee-cell sub-
     problem) can be CLN-expanded directly without re-deriving the
     PDE in the cell.

## What is LOST relative to canonical CLN

  * The basis fields H_n, E_n no longer have a physical
    interpretation -- they are abstract Krylov vectors.  Visualising
    the eddy-current modal pattern is harder.

  * The algorithm is conditional on a starting vector w_1 = G^{-1} l.
    Different starting vectors give different ladder constants;
    this is the "expansion point" question and is examined separately
    in 12/2019_05_29_展開点 (see HIRUMA_EXPANSION_POINT below).

  * Lanczos breakdown (lambda_n = 0 for some n) is possible and must
    be checked at each step.

## Hiruma + multi-expansion point

The 2019 expansion-point notes (Sugahara 2019-06-14, 2019-08-03)
restart Hiruma's Krylov sequence at s0 != 0: define K_{s0} = G + s0 C,
w_1 = K_{s0}^{-1} l, lambda_n alternating between K_{s0} and C
inner products.  Rewriting (G + s C) x = b as (K_{s0} + (s-s0) C) x = b
makes the truncated ladder exact at s = s0 and accurate over a wider
band than the s0 = 0 ladder -- the algebraic counterpart of the
canonical Kuriyama 2019 multi-expansion-point CLN.

Sugahara's experimental finding (2019-06-14):  imaginary-axis shifts
(s0 = j omega_target) widen the accurate frequency band; real-axis
shifts improve the matched-frequency accuracy.  Surprisingly, the
algorithm RUNS even when K_{s0} is not SPD -- an open theoretical
question as of 2019.

## Hiruma + sub-cell unit problem ("sell-ho", #11-08 section 5)

Hiruma applies the algebraic Cauer expansion to a UNIT CELL of a wire
layer (a single conductor cross-section in periodic BC).  Given a
uniform external B_dot, the cell's (K + s N) system is Hiruma-
expanded, the truncated ladder gives an effective mu_eff(s) for the
layer.  This is the "homogenisation" application that motivates the
algebraic approach: each unit cell shape (round, rectangular, litz)
gets its CLN without re-deriving the Kameari PDE recursion.  It is
the algebraic counterpart of Shindo-Kameari-Matsuo 2017's
analytical-sheet Cauer circuit, ready for any cell geometry.
"""

HIRUMA_CODE = """
# Hiruma method: representative code (12_比留間法)

## Scalar 1-port driver: hiruma_scalar.m

```matlab
% public-safe curated corpus
N = 6;
G = symmetric_random(N);    % SPD stiffness
C = symmetric_random(N);    % SPD mass
b1 = rand(N, 1);
b2 = rand(N, 1);
A  = -G \\ C;
s  = linspace(1, 10, 100);

% Reference direct solve
for k = 1 : length(s)
  Y = [b1, b2].' * ((G + s(k)*C) \\ [b1, b2]);
end

% Hiruma recursion for 2 inputs (b1, b2): non-symmetric two-sided
for nb = 1 : 2
  if nb == 1, w1 = G \\ b1; else, w1 = G \\ b2; end
  lambda(1) = w1.' * G * w1;
  w(2).v    = - w1 / lambda(1);
  for n = 1 : N
    lambda(2*n)   = w(2*n  ).v.' * C * w(2*n  ).v;
    w(2*n+1).v    = w(2*n-1).v - (A * w(2*n).v) / lambda(2*n);
    lambda(2*n+1) = w(2*n+1).v.' * G * w(2*n+1).v;
    w(2*n+2).v    = w(2*n).v - w(2*n+1).v / lambda(2*n+1);
  end
end
```

Output: ladder constants (1/lambda_1, lambda_2, 1/lambda_3, ...) =
(R_0, L_1, R_2, ...) of the Cauer network.  Comparing the truncated
continued-fraction Y(s) against the direct solve verifies that
Hiruma's recursion reproduces the reference within ladder-truncation
error.

## Multi-port driver: hiruma_matrix.m

```matlab
% Same setup, but b is now a matrix [b1, b2, b3, b4]
b = rand(N, 4);
w(1).v       = G \\ b;            % matrix of Krylov starting vectors
yy(1).G      = w(1).v.' * G * w(1).v;       % 4x4 inner product
w(2).v       = -(yy(1).G \\ w(1).v.').';
for n = 1 : N
  yy(2*n).G   = w(2*n  ).v.' * C * w(2*n  ).v;
  w(2*n+1).v  = w(2*n-1).v - (yy(2*n).G \\ (A*w(2*n).v).').';
  yy(2*n+1).G = w(2*n+1).v.' * G * w(2*n+1).v;
  w(2*n+2).v  = w(2*n).v   - (yy(2*n+1).G \\ w(2*n+1).v.').';
end
```

The yy(n).G are now block-Lanczos inner-product matrices.  The
resulting Cauer ladder is multi-port (4 in, 4 out).  Truncated Y(s)
matches the direct solve over the chosen s-range -- this is the
"multi-port Hiruma method" used in 13_BEM+FEM below for multi-mode
expansion (H_1 mode, H_2 mode, ... each as a separate port).

## Expansion-point variant: hiruma3.m, hiruma4.m, hiruma5.m

```matlab
% public-safe curated corpus
for nCASE = [1, 2, 3]
  switch nCASE
    case 1, s0 = 0;        % canonical CLN expansion at s = 0
    case 2, s0 = 15;       % shift along the real axis
    case 3, s0 = 15j;      % shift along the imaginary axis
  end
  K  = R + s0 * L;         % shifted stiffness
  w(1).v       = K \\ V;
  lambda(1)    = w(1).v.' * K * w(1).v;
  w(2).v       = -w(1).v / lambda(1);
  for n = 1 : N
    lambda(2*n)   = w(2*n  ).v.' * L * w(2*n  ).v;
    w(2*n+1).v   = w(2*n-1).v + (K \\ (L*w(2*n).v)) / lambda(2*n);
    lambda(2*n+1) = w(2*n+1).v.' * K * w(2*n+1).v;
    w(2*n+2).v   = w(2*n).v - w(2*n+1).v / lambda(2*n+1);
  end
  % evaluate Y(s) over a wider s-band using (s - s0) shifted variables
end
```

s0 = 0 is the baseline.  s0 = 15j gives accuracy improvement on the
imaginary axis (physical frequencies omega = -s/j); s0 = 15 gives
accuracy improvement on the real axis (relaxation eigenvalues).

## Symbolic prototype: hiruma_check_eq126.m, hiruma_check0.m

The hiruma_check series of files verifies, on small (N=2) symbolic
examples, that the Hiruma algorithm reproduces the partial-fraction
expansion of a transfer function H(s) = (9 + 9/5 s) / (1 + 7/10 s -
2/25 s^2) when written as a continued fraction with lambda constants
1/9, 1/18, 162/25, 2/81.  This is the smallest non-trivial test
case for Hiruma.

```matlab
H = 1 ./ (1/9 + 1 ./ (1./(s/18) + 1 ./ (1/(162/25) - 1 ./ (1./(2/81*s)))));
```

(corresponds to lambda_1=9, lambda_2=18, lambda_3=162/25, lambda_4=2/81).
"""

# ---------------------------------------------------------------------------
# 13: BEM + FEM
# ---------------------------------------------------------------------------

BEM_FEM_COUPLED = """
# BEM + FEM coupled CLN (13_BEM+FEM)

## The motivation (Sugahara 2018-07-18 note)

Canonical Kameari CLN puts the eddy-current conductor inside a closed
FEM domain with Dirichlet A=0 on its faces.  For a round-conductor
2-D problem the H_m mode current is then

   J_m(r) = -p_m * sigma * j omega * J_m(k r) * cos(m theta)    (1)

with J_m the m-th Bessel function and k = sqrt(-j omega mu sigma);
Kameari recursion with Dirichlet recovers (1).

For an OPEN-BOUNDARY problem (conductor in infinite air), only the
conductor is meshed; the air is handled by BEM (surface Green's
function).  The conductor boundary is now ROBIN (tangential H and
normal B continuity), giving

   J_m(r) = J_m(k r) * (m/k - mu*(J_{m-1} - J_{m+1})/2)^{-1}
            * (n . (H_inc + Sum_n p_n H_n^scat))   (2)

Same r-dependence J_m(k r) but a different (s-dependent, complex)
prefactor.  Sugahara 2018-07-18: "the integration-by-parts BOUNDARY
TERM picks up a different contribution under Robin than under
Dirichlet, so the lambda_n acquire s-dependence -- no longer the
clean integer sequence of equation (1)."

## The structure of the coupled system

Let the conductor occupy domain Omega_c with surface Gamma; the air
is the unbounded exterior.

  FEM part   (interior, in Omega_c):
     A_z + Robin BC on Gamma:    n . curl A = f(A|_Gamma, H_inc)
     -> for each mode H_m, one Kameari recursion -> one port

  BEM part   (exterior, on Gamma):
     Surface integral equation for the secondary current/magnetisation
     -> a matrix P_BEM(s) coupling all surface DoFs

  Coupling:
     Continuity of A_z and (1/mu) dA_z/dn across Gamma
     -> a multi-mode multi-port Cauer ladder with mode index m and
        port index k (one port per conductor / per mode).

The BEM part is set up using a standard Bio-Savart 2-D kernel (see
Bio2D.m in 13_BEM+FEM/); the FEM part can be done with FEMM (free
2-D FEM, see GenFEMM.m drivers) or FreeFEM (see FreeFEM_orig.edp in
TSVD/円筒導体/).

The result is a MULTI-PORT Hiruma-style Cauer ladder: each H_m mode
contributes one port, the BEM coupling between conductors contributes
the off-diagonal blocks of the multi-port Y matrix.

## Sub-experiment: TSVD basis reduction (TSVD/ folder)

Direct multi-port CLN with M conductors and N_m modes is M*N_m ports.
For M = 6, N_m = 36 (cosine modes on a circular conductor) this is
216 ports -- intractable.  TSVD (Truncated SVD, gen_svd_mode.m)
builds the BEM "interaction matrix" Az between all boundary surface
points and a single source point, does SVD, keeps the top K = 5-15
left singular vectors as a reduced port basis.  The 2018-06-25
internal note 連続関数のSVD derives the continuum analog (KL
expansion of the smooth BEM kernel has exponentially decaying
singular values), which justifies the truncation.

## Sub-experiments: hoverboard CLN + FEMM_*_coil benchmarks

  * ホバーボード_CLN/{Scalar,Matrix}/  -- PM hoverboard over a
    conducting plate, scalar (1-port lumped) vs matrix (multi-port
    H_m-mode) BEM+FEM CLN.  Matrix form gives 5-10x B-field accuracy
    improvement at the cost of larger Hiruma recursion.

  * FEMM_{one,two,six}_coil, FEMM_{dipole,monopole,quadrupole,
    sextupole} -- benchmark recipes (GenFEMM.m / GenFigure1.m /
    GenFigure2.m) covering successively richer coil symmetries vs
    Bessel-mode MoM reference (MoM_A.m).  All verify that the Hiruma
    BEM+FEM CLN converges to the analytic MoM as N (stages) and K
    (TSVD modes) grow.

## State of the work (2018-07-18, latest note)

  * H_1 mode: ladder converges to MoM reference at N=10, K=5.
  * H_2 mode: ladder converges at N=15, K=8.
  * Higher modes: BEM-side discretisation dominates; refine the BEM
    surface before refining N or K.

Open extensions (not yet done): 3-D version, coupling with the PEEC
PRIMA solver in radia.lanczos_reduction, validation against
ngsolve.bem on a 3-D test case.
"""

BEM_FEM_CODE = """
# BEM + FEM CLN: representative code (13_BEM+FEM)

## BEM kernel: Bio2D.m

The Biot-Savart 2-D kernel for a line current at the origin:

```matlab
% public-safe curated corpus
function [hx, hy, az] = Bio2D(x0, y0, X, Y)
  % field at (x0,y0) from a polygonal current loop with vertices (X,Y)
  for k = 1 : length(X)-1
    dx = X(k+1) - X(k);  dy = Y(k+1) - Y(k);
    L  = sqrt(dx^2 + dy^2);
    % ... line-source contribution to (hx, hy, az) ...
  end
end
```

Used by gen_svd_mode.m to build the surface-to-surface Az interaction
matrix that gets SVD-truncated.

## SVD-mode generation: gen_svd_mode.m

```matlab
% rectangle surface discretisation, 4 sides of N points each
Nx = 1001; Ny = 501; dl = 0.01/1000;
edge_x = (-(Nx-1)/2 : (Nx-1)/2) * dl;
edge_y = (-(Ny-1)/2 : (Ny-1)/2) * dl;
x = [edge_x, lx*ones(1,Ny), -edge_x, -lx*ones(1,Ny)];
y = [-ly*ones(1,Nx), edge_y, ly*ones(1,Nx), -edge_y];

% Source: a unit circle of N=360 polygonal points
N_src = 360;
T = linspace(0, 2*pi, N_src+1);  T = T(1:N_src);
X = a*cos(T);  Y = a*sin(T);

% Interaction matrix: A_{m,k} = az from source point k at obs point m
for m = 1 : length(x)
  [hx, hy, at] = Bio2D(x(m), y(m), X, Y);
  Az(m, :) = at;
end
[u, s, v] = svd(Az);       % TSVD of the BEM kernel
```

The top K columns of u (left singular vectors) become the reduced
basis for the BEM coupling matrix; this basis is fed back into the
Hiruma ladder recursion as the multi-port input.

## Bessel-mode multi-coil reference: MoM_A.m (FEMM_two_coils/)

For two parallel cylindrical conductors, the analytic eddy-current
distribution is a sum of Bessel J_m modes inside each conductor and
1/r^m harmonics outside.  Coupling between the two is captured by
addition-theorem-style integrals.  The driver:

```matlab
% FEMM_two_coils/MoM_A.m  (Method of Moments reference)
k = sqrt(-1j * omega * mu0 * mu * sig);   % skin-depth wavenumber
% Build cross-coupling matrix A of size 2N x 2N (N modes per conductor)
for n = 1 : N
  % qn = exterior mode coefficient transformer at the conductor surface
  qn = a^n * (n*besselj(n, k*a) - 1/mu * a*k/2 *  ...
              (besselj(n-1, k*a) - besselj(n+1, k*a))) /  ...
             (n*besselj(n, k*a) + 1/mu * a*k/2 *  ...
              (besselj(n-1, k*a) - besselj(n+1, k*a)));
  % field at conductor-2 surface from a unit mode n at conductor-1
  for nt = 1 : length(theta)
    r = sqrt((a*cos(theta(nt))+xs(2) - xs(1))^2 +  ...
             (a*sin(theta(nt))+ys(2) - ys(1))^2);
    az(nt) = -qn * cos(n*theta(nt)) / r^n;
  end
  A(N+m, n) = -trapz(theta, az.*cos(m*theta)) / pi;
end
Af = A \\ b;       % solve for mode coefficients
```

The script produces the analytic J_z(x, y) and A_z(x, y) maps that
serve as goldenset for the Hiruma BEM+FEM CLN run.

## Cylindrical-conductor reference with Kelvin transformation

TSVD/円筒導体_Kelvin/check_basis_analytic.m drives a symbolic_Robin
helper `H1(r, theta, a, sig, mu)` that returns the analytic CLN basis
A_n(r), J_n(r) for the H_1 mode under a Robin boundary condition.
This is the gold standard against which any Hiruma-style algebraic
CLN driver is validated.  The companion TSVD/円筒導体/ folder runs
the same problem with a FreeFEM .edp interior solve (each TSVD mode
becomes a Dirichlet boundary on Gamma; the interior A is then fed
back into the Kameari recursion through surface integration -- the
2-D Kelvin-transformation route for unbounded eddy-current problems).
"""

# ---------------------------------------------------------------------------
# 15: NAGAMINE ERROR THEORY
# ---------------------------------------------------------------------------

NAGAMINE_ERROR_THEORY = """
# Nagamine error theory (15_誤差論@長嶺氏)

## Who is Nagamine

Hideaki Nagamine (長嶺 英朗) is a researcher / Ph.D. student in
Tetsuji Matsuo's Lab at Kyoto University -- the same Matsuo as the
canonical Kameari-Ebrahimi-Sugahara-Shindo-Matsuo 2018 paper.  He
specialises in VERIFIED NUMERICAL COMPUTATION (interval arithmetic,
affine arithmetic) and applies it to the CLN method to give
rigorous bounds on the circuit constants R_n, L_n.

His committee paper is:

  * Nagamine et al., "Convergence analysis of the Cauer ladder
    network method using eigenfunction expansion", 2022 IEEE 20th
    Biennial Conference on Electromagnetic Field Computation
    (CEFC 2022).  This is the PUBLISHED form of the lab notes
    referenced below.

  * His doctoral / committee output appears as Nag+22a, Nag+22b in
    the bibliographies of his lab notes.

## The 4 PDFs in 15_誤差論@長嶺氏

### PDF 1: 2023-03-25-clnsqkv.pdf (Nagamine, 2023-03-25)

"Numerical solution of CLN for square cross-section" (正方形領域の
CLN 法の数値解), 26 pages.  For a 1 mm x 1 mm infinite square
conductor (mu = 100 mu_0, sigma = 1e4 S/m), computes the EXACT
R_n, L_n by interval arithmetic (IA) and affine arithmetic (AA) via
the Stieltjes transform + quotient-difference algorithm (商差法) --
analytic, not FEM.  Runtime 9 min 13 s for L_1 through R_18.  Affine
arithmetic (Comba-Stolfi correlation tracking) tightens the intervals
up to ~4x vs plain IA (R_18 ratio 1.5, L_17 ratio 3.2).  Also
compares Lanczos vs Arnoldi: Lanczos needs SPD G, Arnoldi handles
non-symmetric G at higher cost; both succeed for this problem.

### PDF 2: 2023_01_24_有限要素解析とCLN の回路定数の精度.pdf

"Precision of CLN circuit constants from finite-element analysis
(shared with Sugahara)", 9 pages.  Compares FEM-computed R_n, L_n
(FreeFEM++, rect.cln.edp, mesh quality MQ = 100..1500) against the
verified-arithmetic exact values from PDF 1.  Figure 1 plots the
relative error vs n for each MQ; the error is NOT monotonically
decreasing in n -- stages 8-9 onward diverge from exact for MQ < 1500.
Quote (page 5): "the error of circuit constants tends to grow as
the number of stages increases".  The convergent-in-N CLN ladder is
NOT convergent for any finite mesh.  The FreeFEM driver is included
verbatim in the appendix:

```cpp
// rect.cln.edp (excerpt, page 4)
Vh dA, w;
solve Prob(dA, w, eps = SolverEps)
  = int2d(Th)(Nu * rotz(w)' * rotz(dA))
  - int2d(Th)(Sigma * w * Eb[n] / G[n])
  + on(BoxWallLabelDirichlet, dA = 0);
Ab[n+1] = Ab[n] + dA;
L[n+1]  = int2d(Th)(Nu * rotz(Ab[n+1])' * rotz(Ab[n+1]));
Eb[n+1] = Eb[n] - Ab[n+1] / L[n+1];
G[n+1]  = int2d(Th)(Sigma * Eb[n+1] * Eb[n+1]);
```

Conclusion: each stage's static FEM solve depends on the previous
stage's basis functions, so FEM error in A_n is amplified through the
recursion.

### PDF 3: 2023_02_02_CLN 法の回路定数の有限要素離散化誤差.pdf

"Proposal on the FE discretisation error of CLN circuit constants
(shared with Sugahara)", 30 pages.  Defines:

  * characteristic frequency f_n = max_i (R_{2i} / (2 pi L_{2i-1}))
    over the n cut-off frequencies of the Foster equivalent of the
    n-stage Cauer ladder (skin-depth analog for multi-frequency CLN)
  * characteristic length delta_n = 1 / sqrt(pi f_n sigma mu)

For the 1-D magnetic-sheet E-mode problem (Shindo-Kameari-Matsuo
2017a):  R_{2n} = sigma d / (4n+1), L_{2n+1} = mu d / (4*(4n+3)).
Numerically, f_n ~ n^{3.8-3.9}, matching the SKM17b O(n^4) effective-
frequency-range scaling.

MAIN EMPIRICAL FINDING (page 4): relative error of R_n, L_n vs the
ratio delta_n / Delta_x (mesh size):

   delta_n / Delta_x = 100  ->  relative error 0.001 %
   delta_n / Delta_x =  10  ->  relative error 0.1 %
   delta_n / Delta_x =   1  ->  relative error several %

i.e. delta_n shrinks 10x -> relative error grows 100x.  Same scaling
holds for the cylindrical conductor (R_{2n} = (2n+1)/(sigma pi a^2),
L_{2n+1} = mu/(8(n+1)pi), Shindo-Matsuo 2019).

PRACTICAL RULE: trust only stages with Delta_x <= delta_n / 10.
This is the FIRST quantitative mesh-adequacy rule for CLN.

Anomaly: for Delta_x < 1 um the scaling breaks down -- relative
error grows without delta_n shrinking.  Nagamine speculates floating-
point round-off in the recursive FE solve.  Open question.

### PDF 4: foster回路表現導出(Rm)_021726_watai.pdf (Watai, 2026)

Author appears to be Watai (different person, not Nagamine).  5 pages
but extracted text is EMPTY -- likely an image-only scan.  Title
suggests "Derivation of Foster circuit representation R_m" -- this is
the Foster (partial fraction) dual of the Cauer (continued fraction)
ladder representation, used in Nagamine's characteristic-frequency
definition above.

## Why this matters for the lab

Before Nagamine's work, the CLN method was used with the IMPLICIT
assumption that all R_n, L_n are correct as long as the mesh is
"reasonable".  Nagamine showed:

  1. Even at MQ = 1500 (FreeFEM mesh size 1/1500), stages 8+ have
     >10 % error in R, L for the square cross-section.

  2. There is a quantitative DESIGN RULE: pick the maximum stage n
     such that delta_n >= 10 * Delta_x.  Beyond that the FEM
     discretisation dominates.

  3. The characteristic length delta_n is a CHEAP proxy: it is
     computed from the (already-available) Cauer ladder constants,
     no extra FEM run is needed.

  4. The Foster equivalent (partial-fraction form) is more numerically
     stable than the Cauer (continued-fraction form) at high stage
     numbers, because the cut-off frequencies f_{ni} converge to
     fixed values as n grows.

This closes a 5-year gap between Kameari 2018 (existence and
convergence of the ladder) and the practical question "how many
stages should I trust for my mesh?".
"""

# ---------------------------------------------------------------------------
# LAB LESSONS
# ---------------------------------------------------------------------------

LAB_LESSONS = """
# Lab lessons across the four specialty methods

1. TERMINATION (08) is a polish, not a substitute for stages.
   Sugahara 2016-12-23: optimising the terminator gives the same
   accuracy as adding 2-3 stages.  Add stages first; polish later.

2. HIRUMA (12) is the right entry point for a black-box FEM.
   Given (K, N) matrices from any solver (FreeFEM, ngsolve, FEMM,
   PEEC exporter), Hiruma's algebraic Lanczos produces the Cauer
   ladder without re-deriving Kameari's PDE recursion per cell shape.

3. HIRUMA EXPANSION-POINT shifts: imaginary-axis first.
   Sugahara 2019: real-axis s0 improves accuracy at the matched
   relaxation frequency but degrades it on the imaginary (physical
   omega) axis.  Imaginary-axis s0 = j omega_target widens accuracy
   over a band of physical frequencies -- the default choice.

4. Non-SPD shifted system "just works" for Hiruma -- empirically.
   K_{s0} = G + s0 C may not be SPD for s0 != 0, yet Hiruma's
   Lanczos still converges.  Open theoretical question; treat as a
   useful surprise, not a guarantee.

5. BEM+FEM (13) needs Robin BC, not Dirichlet.
   Canonical Kameari assumes closed Dirichlet A=0.  BEM coupling
   requires Robin (tangential H + normal B continuity); lambda_n
   change and acquire s-dependence.  Do NOT reuse a Dirichlet driver.

6. TSVD truncation makes multi-port BEM+FEM tractable.
   M conductors x N_m modes -> SVD-truncated to top K = 5-15
   effective ports.  Beyond K the BEM discretisation error dominates.

7. NAGAMINE (15): trust only delta_n >= 10 * Delta_x.
   delta_n is computed in O(N) from already-available R_n, L_n
   (Foster cut-off).  Do NOT trust the first n stages just because
   they look converged.

8. CLN ladder is NOT monotone-convergent in n for a finite mesh.
   Both R and L can grow in error vs n for fixed Delta_x (Nagamine
   2023-01-24).  Do NOT use "the error stopped decreasing" as a
   stopping criterion; use delta_n >= 10 * Delta_x.

9. Foster (partial fraction) is more stable than Cauer (continued
   fraction) at high n.  Foster cut-off frequencies converge; the
   continued-fraction tail does not.  Use Foster for diagnostics
   (NOT for SPICE export -- Cauer is the canonical export form).

10. Nagamine's delta_n definition is one CHOICE.  Other proxies
    (highest pole of Z, last-stage R_n/L_n) oscillate with n;
    Foster cut-off converges.  Cite Nagamine's definition explicitly.
"""

# ---------------------------------------------------------------------------
# KEY FILES
# ---------------------------------------------------------------------------

KEY_FILES = """
# Key files cited (paths relative to public-safe curated corpus)

## 08_終端問題  (Termination)

  08_終端問題/2016_12_23_Cauer梯子回路の終端の検討.docx
      Original 2016 note, Compumag 2017 supporting material.  Key
      finding: "adding stages > optimising terminator".
  08_終端問題/2018_03_13_CLN_termination/main.m
      Main driver: 1-5 stages, with/without fminsearch terminator,
      vs analytic tan(z)/z over 1 kHz - 10 MHz.
  08_終端問題/2018_03_13_CLN_termination/cf.m
      Continued-fraction evaluator with terminator slot (8 lines).
  08_終端問題/2018_03_13_CLN_termination/termination.m
      Symbolic study of polynomial-correction terminators.
  08_終端問題/2018_03_13_CLN_termination/2018_03_17_CLN回路の終端
      インピーダンス.docx
      Full 2018-03-17 note: Z_term construction for E- and H-mode.

## 12_比留間法  (Hiruma)

  12_比留間法/2018_03_30_比留間法/hiruma_scalar.m
      Scalar 1-port Hiruma driver, 2-input verification.
  12_比留間法/2018_03_30_比留間法/hiruma_matrix.m
      Block multi-input Hiruma driver (4 ports).
  12_比留間法/2018_03_30_比留間法/hiruma_check_eq126.m
      Symbolic N=2 verification of Hiruma equation (126).
  12_比留間法/2018_03_30_比留間法/先行文献/{9-04, 11-08}.pdf
      Hiruma's two committee papers (symmetric / non-symmetric
      Lanczos derivations of the Cauer ladder).
  12_比留間法/2018_03_30_比留間法/先行文献/clnb_0309.pdf
      Hiruma's follow-up: "Kameari fields are Lanczos vectors".
  12_比留間法/2019_05_29_展開点に関する検討/hiruma{3,4,5}.m
      Expansion-point variants: s0 = 0 (canonical), s0 = 15 (real
      shift), s0 = 15j (imaginary shift).
  12_比留間法/2019_05_29_展開点に関する検討/2019_08_03_展開点に関する
      検討.pdf
      Sugahara 2019-08 note: multi-expansion-point Hiruma for the
      cylindrical E0 mode.

## 13_BEM+FEM  (BEM + FEM coupled CLN)

  13_BEM+FEM/2018_07_18_現在の検討状況のまとめ.docx
      July 2018 state summary: Dirichlet -> Robin BC issue, the
      s-dependence of lambda_n in the coupled case.
  13_BEM+FEM/Bio2D.m, Gauss2D.m
      2-D BEM kernel and quadrature (used by all 13_ drivers).
  13_BEM+FEM/TSVD/gen_svd_mode.m
      SVD-mode basis generation + FreeFEM .edp emission per mode.
  13_BEM+FEM/TSVD/円筒導体_Kelvin/check_basis_analytic.m
      H_1 mode CLN basis vs analytic Bessel (Kelvin transform).
  13_BEM+FEM/FEMM_two_coils/MoM_A.m
      MoM analytic reference for 2-coil + 1-conductor benchmark.
  13_BEM+FEM/ホバーボード_CLN/{Scalar,Matrix}/plt_B_map.m
      Two-port vs multi-port Hiruma BEM+FEM for hoverboard.
  13_BEM+FEM/参考文献/Efficient circuit representation of
      eddy-current fields.pdf
      Shindo-Kameari-Matsuo COMPEL 2017 (canonical sheet / wire
      reference).
  13_BEM+FEM/参考文献/固定部と可動部のCLN 結合.pdf
      Stator-rotor CLN coupling (motivation for multi-port block
      extension).

## 15_誤差論@長嶺氏  (Nagamine error theory)

  15_誤差論@長嶺氏/2023-03-25-clnsqkv.pdf  (26 pages)
      Verified arithmetic L_1..R_18 for square cross-section.
  15_誤差論@長嶺氏/2023_02_07_長嶺氏_CLN誤差/2023_01_24_有限要素解析
      とCLN の回路定数の精度.pdf  (9 pages)
      FreeFEM-vs-exact: non-monotone error + rect.cln.edp source.
  15_誤差論@長嶺氏/2023_02_07_長嶺氏_CLN誤差/2023_02_02_CLN 法の回路
      定数の有限要素離散化誤差.pdf  (30 pages)
      Characteristic frequency / length, delta_n >= 10 Delta_x rule.
  15_誤差論@長嶺氏/foster回路表現導出(Rm)_021726_watai.pdf  (5 pages)
      Watai's Foster-form derivation (image-only scan).
"""

# ---------------------------------------------------------------------------
# CITATIONS
# ---------------------------------------------------------------------------

CITATIONS = """
# Citations (BibTeX-ish)

## Hiruma method primary sources

@misc{Hiruma2017IEEJa,
  author = {Shingo Hiruma and Hajime Igarashi},
  title  = {On a Lanczos algorithm related to Kameari's method (#9-04)},
  note   = {Committee paper, IEEJ Joint Technical Meeting on Static
            Apparatus and Rotating Machinery (静止器・回転機合同研究会).
            Hokkaido University, Igarashi Lab.  Archived as 9-04.pdf},
  year   = {TBD}
}

@misc{Hiruma2017IEEJb,
  author = {Shingo Hiruma and Hajime Igarashi},
  title  = {On a Lanczos algorithm related to Kameari's method, part 2
            (#11-08)},
  note   = {Extends with eddy-current sub-cell (sell-ho) application.
            Archived as 11-08.pdf},
  year   = {TBD}
}

@misc{Hiruma2018CLNb,
  author = {Shingo Hiruma and Hajime Igarashi},
  title  = {CLN follow-up note (clnb_0309)},
  note   = {2018-03-09.  Identifies Kameari basis fields with Lanczos
            vectors.  Archived as clnb_0309.pdf},
  year   = {2018}
}

@article{Hiruma2022CLNerror,
  author  = {Shingo Hiruma and others},
  title   = {Error estimator for Cauer ladder network representation},
  journal = {IEEE Transactions on Magnetics},
  volume  = {58}, number = {9}, year = {2022}, pages = {1--4},
  note    = {Art no. 7500904},
  doi     = {10.1109/TMAG.2022.3165943}
}

## Nagamine error theory primary sources

@misc{Nagamine2023a,
  author = {Hideaki Nagamine},
  title  = {Precision of CLN circuit constants from finite-element
            analysis (shared with Sugahara)},
  note   = {2023-01-24 lab note, 9 pages.  Kyoto, Matsuo Lab.
            Archived in public-safe curated corpus
            CLN誤差/},
  year   = {2023}
}

@misc{Nagamine2023b,
  author = {Hideaki Nagamine},
  title  = {Proposal on the FE discretisation error of CLN circuit
            constants (shared with Sugahara)},
  note   = {2023-02-02 lab note, 30 pages.  Defines f_n, delta_n;
            derives the delta_n >= 10 Delta_x rule},
  year   = {2023}
}

@misc{Nagamine2023c,
  author = {Hideaki Nagamine},
  title  = {Numerical solution of CLN for square cross-section},
  note   = {2023-03-25 lab note, 26 pages.  Verified arithmetic
            L_1..R_18 with interval + affine arithmetic.  Archived as
            2023-03-25-clnsqkv.pdf},
  year   = {2023}
}

@inproceedings{Nagamine2022CEFC,
  author    = {Hideaki Nagamine and others},
  title     = {Convergence analysis of the Cauer ladder network method
               using eigenfunction expansion},
  booktitle = {Proc. 2022 IEEE 20th Biennial CEFC},
  year      = {2022},
  note      = {Published form; Nag+22a, Nag+22b in his lab notes}
}

@misc{Watai2026,
  author = {Watai},
  title  = {Derivation of Foster circuit representation R_m},
  note   = {Internal lab note (file timestamp 2026 inferred from
            021726).  5 pages, image-only scan.  Year TBD},
  year   = {TBD}
}

## BEM+FEM coupling references

@article{Shindo2017COMPEL,
  author  = {Yuji Shindo and Akihisa Kameari and Tetsuji Matsuo},
  title   = {Efficient circuit representation of eddy-current fields},
  journal = {COMPEL},
  year    = {2017},
  doi     = {10.1108/COMPEL-02-2017-0084},
  note    = {Cauer ladder for magnetic sheet (E and H mode) from
             continued fraction of tanh(kd/2) AND from Legendre-
             polynomial expansion of H_y.  Required reading for
             13_BEM+FEM; also Nagamine's test-problem reference}
}

@misc{Sugahara2018BEMFEM,
  author = {Ken Sugahara},
  title  = {Current state of BEM-FEM CLN coupling (internal note)},
  note   = {2018-07-18.  Identifies Dirichlet -> Robin BC and the
            s-dependence of Kameari lambda_n.  Archived as
            2018_07_18_現在の検討状況のまとめ.docx},
  year   = {2018}
}

@misc{LabBEMFEMCoupling,
  author = {Sugahara Lab (internal)},
  title  = {固定部と可動部のCLN 結合 (Stator-rotor CLN coupling)},
  note   = {Internal reference PDF, motivation for the multi-port
            Hiruma block extension},
  year   = {TBD}
}

## Termination internal notes

@misc{Sugahara2016Termination,
  author = {Ken Sugahara},
  title  = {Cauer 梯子回路の終端の検討},
  note   = {2016-12-23.  Compumag 2017 supporting material.
            Conclusion: adding stages > optimising terminator},
  year   = {2016}
}

@misc{Sugahara2018Termination,
  author = {Ken Sugahara},
  title  = {CLN 回路の終端インピーダンス},
  note   = {2018-03-17.  Defines Z_term^E(s), Z_term^H(s) from
            analytic R_n, L_n of the magnetic sheet},
  year   = {2018}
}

## Base CLN paper (parent module: cln_knowledge.py)

@article{Kameari2018CLN,
  author  = {Akihisa Kameari and S. Ebrahimi and Ken Sugahara and
             Yuji Shindo and Tetsuji Matsuo},
  title   = {Cauer ladder network representation of eddy-current fields
             for model order reduction using finite-element method},
  journal = {IEEE Transactions on Magnetics},
  volume  = {54}, number = {3}, year = {2018}, pages = {1--4},
  note    = {Art no. 7201804},
  doi     = {10.1109/TMAG.2017.2743224}
}

@article{Shindo2019Cylinder,
  author = {Yuji Shindo and Tetsuji Matsuo},
  title  = {Cauer ladder representation of the cylindrical conductor
            eddy-current field},
  note   = {Cited as SM19 in Nagamine 2023-02-02.  Source of
            R_{2n} = (2n+1)/(sigma pi a^2), L_{2n+1} = mu/(8(n+1)pi).
            Year TBD},
  year   = {2019}
}
"""

# ---------------------------------------------------------------------------
# Dispatch
# ---------------------------------------------------------------------------

_TOPICS = {
    "overview": OVERVIEW,
    "termination_theory": TERMINATION_THEORY,
    "termination_code": TERMINATION_CODE,
    "hiruma_method": HIRUMA_METHOD,
    "hiruma_code": HIRUMA_CODE,
    "bem_fem_coupled": BEM_FEM_COUPLED,
    "bem_fem_code": BEM_FEM_CODE,
    "nagamine_error_theory": NAGAMINE_ERROR_THEORY,
    "lab_lessons": LAB_LESSONS,
    "key_files": KEY_FILES,
    "citations": CITATIONS,
    "all": "",
}
_TOPICS["all"] = "\n\n".join(v for k, v in _TOPICS.items() if k != "all")


def get_cln_specialty_documentation(topic: str = "all") -> str:
    """Dispatch by topic.

    Topics:
      "all"                    - all sections concatenated
      "overview"               - big picture across the 4 specialty methods
      "termination_theory"     - 08, why terminate, what to terminate with
      "termination_code"       - 08, representative MATLAB drivers
      "hiruma_method"          - 12, Hiruma's algebraic non-symmetric Lanczos
      "hiruma_code"            - 12, MATLAB drivers (scalar, matrix,
                                 expansion-point variants)
      "bem_fem_coupled"        - 13, BEM exterior + FEM interior coupled CLN
      "bem_fem_code"           - 13, MATLAB / FreeFEM drivers
      "nagamine_error_theory"  - 15, verified arithmetic, characteristic
                                 frequency / length, delta_n >= 10 Delta_x rule
      "lab_lessons"            - 10 cross-cutting lessons
      "key_files"              - file inventory with paths
      "citations"              - BibTeX-ish references
    """
    topic = topic.lower().strip()
    if topic in _TOPICS:
        return _TOPICS[topic]
    return ("Unknown topic '" + topic + "'. Available: "
            + ", ".join(sorted(_TOPICS.keys())))
