#!/usr/bin/env wolframscript
(* ::Package:: *)

(*
  rect2D_schur_F.wls

  3-axis ROM evaluation for the 2D infinite-length rectangular prism
  scalar diffusion problem.  Parallels cuboid3D_schur_F.wls (3D)
  and slab_cpe_schur.wls (1D).  This is the journal section V.0
  verification, cross-checking against Nagamine et al., "Verified
  Numerical Computations of the Cauer Network Representation of a
  Square Prism Conductor," J. Magn. Magn. Mater. (2026, JJIAM
  submission). The Nagamine reference provides interval-arithmetic
  bounds on Cauer rungs {L_n, R_n} for the square case (a = b).

  Reference problem (per-unit-length, infinite z extent):
     -Laplacian u + s mu sigma u = -s mu sigma   on  Omega = [0,a]x[0,b]
     u = 0  on  dOmega
     u_ext = 1

  Foster spectrum (odd m,n only; even modes uncouple from uniform src):
     lambda_{mn} = (m pi/a)^2 + (n pi/b)^2
     tau_{mn}    = mu sigma / lambda_{mn}
     g_{mn}^2    = 64 sigma a b / ((m n)^2 pi^4)     (per-unit-length)

  Foster admittance (per-unit-length, S/m):
     Y(s) = Sum_{odd m,n} g_{mn}^2 / (1 + s tau_{mn})
     Y(0) = sigma a b                         (DC, Nagamine eq.20)
     Y(s) -> K_SIBC / sqrt(s) at large s
     K_SIBC = 2(a+b) sqrt(sigma/mu)           (perimeter formula)

  ROMs compared on three axes (A wall / B Taylor / C asymptote):
     CLN_N3 / CLN_N10 at s_0 = 0  (canonical Cauer Pade)
     Schur-F (CLN_N3 + 1 sqrt(s) block, basis 4) -- the test

  Run: wolframscript -file rect2D_schur_F.wls
*)

ClearAll[a, b, sigma, mu, nu, Mmax];

(* Geometry: 2D square prism cross-section, a = b = 5 mm Cu.
   Matches Nagamine et al. 2026 JJIAM "square case" with a = 5 mm. *)
a       = 5*^-3;
b       = 5*^-3;
sigma   = 5.96*^7;       (* Cu (matches digest cylinder; Nagamine uses sigma = 5.8e7) *)
mu      = 4 Pi 1*^-7;
nu      = 1/mu;

(* SIBC constant via cross-section perimeter *)
Pperim   = 2 (a + b);
KSIBC2D  = Pperim Sqrt[sigma/mu];

Print["=== rect2D_schur_F.wls ==="];
Print["a = ", N[a 1000], " mm, b = ", N[b 1000], " mm"];
Print["P = 2(a+b) = ", N[Pperim 1000, 6], " mm"];
Print["sigma = ", sigma, " S/m, mu = ", N[mu, 6], " H/m"];
Print["K_SIBC^{2D-rect} = P sqrt(sigma/mu) = ", N[KSIBC2D, 8]];

(* === Foster spectrum truncation ======================================== *)
(* Odd modes m,n = 1,3,5,...,Mmax. Total modes = ((Mmax+1)/2)^2 *)
Mmax = 99;   (* 50 x 50 = 2500 modes; DC converges to ~5 digits and
                f=1e9 Hz asymptote captures the bulk of the SIBC tail.
                For full 12-digit accuracy use Mmax=999 (250000 modes;
                ~30 s runtime), tracked in Nagamine cross-check. *)
oddRange = Range[1, Mmax, 2];

(* Foster eigenvalues, time constants, couplings (per-unit-length).
   g^2 ABSORBS sigma to match slab/cuboid convention:
   Y(s) = Sum g^2/(1 + s*tau) *)
lamMN[m_, n_] := N[(m Pi/a)^2 + (n Pi/b)^2];
tauMN[m_, n_] := mu sigma / lamMN[m, n];
g2MN[m_, n_]  := N[64 sigma a b / ((m n)^2 Pi^4)];

(* Enumerate modes (flat index) *)
modeList = Flatten[Table[{m, n}, {m, oddRange}, {n, oddRange}], 1];
Nmodes   = Length[modeList];
lamVec   = N[Table[lamMN[mn[[1]], mn[[2]]], {mn, modeList}]];
tauVec   = mu sigma / lamVec;
g2Vec    = N[Table[g2MN[mn[[1]], mn[[2]]], {mn, modeList}]];

Print["Nmodes = ", Nmodes, " (m,n odd up to ", Mmax, ")"];

(* Verify DC: Y(0) = Sum g^2 should equal sigma a b *)
yDC = Total[g2Vec];
Print["Y(0) Foster sum   = ", N[yDC, 8]];
Print["Y(0) analytic     = sigma a b = ", N[sigma a b, 8]];
Print["DC convergence    = ", N[yDC/(sigma a b), 8], " (target 1)"];

(* === Modal operators in Foster basis ==================================== *)
LambdaM  = N[DiagonalMatrix[lamVec]];
sigmaM   = N[sigma IdentityMatrix[Nmodes]];
Imat     = N[IdentityMatrix[Nmodes]];

Kmodal[s_] := nu LambdaM + sigma s Imat;
K0matrix   = Kmodal[0];

(* Source vector: uniform u_ext = 1 -> RHS proportional to projection
   onto each mode.  Match slab/cuboid convention: b_k ~ sqrt(g_k^2/tau_k) *)
bvec = N[Sqrt[g2Vec / tauVec]];

(* Exact Foster admittance (truncated double sum, the "ground truth") *)
YExactFoster[omega_] := Module[{s = I omega},
  Total[g2Vec / (1 + s tauVec)]];

(* === K_0 inner product utilities ======================================= *)
K0Inner[u_, v_]     := u . K0matrix . v;
K0MGS[vIn_, basis_] := Module[{v = vIn, ai, cc},
  Do[ai = basis[[i]];
     cc = K0Inner[ai, v];
     v  = v - cc ai, {i, 1, Length[basis]}];
  v];
K0Normalise[v_] := Module[{nrm = Sqrt[Re[K0Inner[v, v]]]},
  If[nrm < 1*^-300, v, v/nrm]];

(* === CLN Krylov basis at s = s0 (canonical Cauer Pade at s=0) ========= *)
buildCLNbasis[s0_, M_] := Module[{Q = {}, vk, KsInv},
  KsInv = Inverse[Kmodal[s0]];
  vk    = KsInv . (sigmaM . bvec);
  vk    = K0Normalise[vk];
  AppendTo[Q, vk];
  Do[vk = KsInv . (sigmaM . Last[Q]);
     vk = K0MGS[vk, Q];
     vk = K0Normalise[vk];
     AppendTo[Q, vk], {M - 1}];
  Q];

(* === Reduced admittance (no Schur) ==================================== *)
YReducedFromBasis[Qlist_, omega_] := Module[{Q, Kr, Qb},
  Q  = Qlist;
  Kr = Q . Kmodal[I omega] . Transpose[Q];
  Qb = Q . bvec;
  sigma (Qb . LinearSolve[Kr, Qb])];

(* === Schur-F 1x1 termination ========================================== *)
(* Y_R(s) = Y_CLN(s) + K_SIBC sqrt(s) / (s + d) *)
YReducedSchurF[Qlist_, omega_, d_] := Module[{yCLN, s, yExtra},
  yCLN   = YReducedFromBasis[Qlist, omega];
  s      = I omega;
  yExtra = KSIBC2D Sqrt[s] / (s + d);
  yCLN + yExtra];

(* === Sanity: DC and asymptote of exact Foster ========================= *)
Print[""];
Print["=== Sanity check on exact Foster Y(s) ==="];
yHF = YExactFoster[2 Pi 1*^9];
KempHF = Abs[yHF] Sqrt[2 Pi 1*^9];
Print["|Y(2pi 1e9)| sqrt(omega) num = ", N[KempHF, 8]];
Print["K_SIBC^{2D-rect} formula     = ", N[KSIBC2D, 8]];
Print["asymptote ratio (num/formula) = ", N[KempHF / KSIBC2D, 6]];
Print["  (NOTE: 121 modes is not enough at 1e9 Hz; this checks the trend)"];

(* === Extract first few Cauer rungs (L_n, R_n) for Nagamine cross-check  *)
(* Lanczos tridiagonalization of (K0^{-1} sigma) on b = bvec yields the
   tridiagonal (alpha_i, beta_i), and the Cauer rungs come from the
   Stieltjes recursion:
       L_{2k-1} = 1 / (alpha_k - L_{2k-3}^{-1} (1/R_{2k-2}))   [odd]
       R_{2k}   = beta_k^2 L_{2k-1} L_{2k+1}                     [even]
   For verification against Nagamine, we use the simpler "moment-based"
   extraction: Taylor coefficients c_k of Y(s) at s=0 -> Hankel-Pade
   QD algorithm -> Cauer rungs.
*)
Print[""];
Print["=== Cauer rung extraction (for Nagamine cross-check) ==="];

(* Compute Taylor moments c_k = (-1)^k * Sum_m g_m^2 * tau_m^k *)
moments[K_] := Table[
  N[Sum[g2Vec[[i]] (-tauVec[[i]])^k, {i, 1, Nmodes}], 30],
  {k, 0, K-1}];

cMoments = moments[20];
Print["First 5 Taylor coeffs c_k of Y(s) at s=0:"];
Do[Print["  c_", k, " = ", N[cMoments[[k+1]], 8]], {k, 0, 4}];

(* Cauer rungs via QD algorithm on the Hankel matrix of moments.
   For verified comparison with Nagamine: he extracts L_1, R_2, L_3, R_4,
   ... using interval arithmetic on the same Foster moments.  Here we
   use machine-precision Lanczos as the cross-check baseline.
   FULL QD implementation in qd_extract.m4 (TODO: port to Mathematica).
   For now, output the moments themselves -- Nagamine's Table 2 lists
   moment bounds we can directly compare. *)
Export["rect2D_schur_F_moments.csv",
  Prepend[Table[{k, cMoments[[k+1]]}, {k, 0, Length[cMoments]-1}],
         {"k", "c_k"}], "CSV"];
Print["Wrote rect2D_schur_F_moments.csv (", Length[cMoments],
      " moments for Nagamine cross-check)"];

(* === Build ROMs ======================================================== *)
Print[""];
Print["=== Building ROMs ==="];
QCLN3   = buildCLNbasis[0, 3];
QCLN10  = buildCLNbasis[0, 10];
Print["  CLN_N3   basis size = ", Length[QCLN3]];
Print["  CLN_N10  basis size = ", Length[QCLN10]];

(* === Frequency sweep + CSV export ===================================== *)
fHzList = N[10^Subdivide[0, 9, 90]];   (* 91 points, 1 Hz to 1e9 Hz *)

dSchur = N[2 Pi 1*^5];   (* digest crossover; tune if wall-band fit is poor *)
Print[""];
Print["=== Frequency sweep (1 Hz - 1e9 Hz, 91 points) ==="];
Print["Schur-F crossover d/(2pi) = ", N[dSchur/(2 Pi), 6], " Hz"];

results = Table[
  Module[{omega = 2 Pi f, yEx, yCLN3, yCLN10, ySchur, kref},
    omega  = 2 Pi f;
    yEx    = YExactFoster[omega];
    yCLN3  = YReducedFromBasis[QCLN3,  omega];
    yCLN10 = YReducedFromBasis[QCLN10, omega];
    ySchur = YReducedSchurF[QCLN3, omega, dSchur];
    kref   = KSIBC2D / Sqrt[omega];
    {f, Abs[yEx], Abs[yCLN3], Abs[yCLN10], Abs[ySchur], kref}],
  {f, fHzList}];

Export["rect2D_schur_F.csv",
  Prepend[results,
   {"f_Hz", "Y_exact", "Y_CLN_N3", "Y_CLN_N10", "Y_Schur_F_b4", "K_div_sqrts"}],
  "CSV"];
Print["Wrote rect2D_schur_F.csv (", Length[results], " rows)"];

(* === 3-axis summary ==================================================== *)
Print[""];
Print["=== 3-axis summary (square case a=b=5mm Cu) ==="];

(* A-axis: wall band peak rel err *)
wallMask = Select[Range[Length[fHzList]],
  1*^5 <= fHzList[[#]] <= 1*^7 &];
walls = {};
Do[Module[{label, idx, name = labelpair[[2]], col = labelpair[[1]]},
  idx = col;
  Module[{errs = Table[
       Abs[results[[k, idx]] - results[[k, 2]]] / Abs[results[[k, 2]]],
       {k, wallMask}]},
    AppendTo[walls,
      {name, N[100 Max[errs], 4]}]]],
  {labelpair, {{3, "CLN_N3"}, {4, "CLN_N10"}, {5, "Schur-F b4"}}}];
Print[StringPadRight["ROM", 16] <> "wall peak err"];
Do[Print[StringPadRight[w[[1]], 16] <> ToString[w[[2]]] <> " %"], {w, walls}];

(* C-axis: asymptote ratio at 1e9 Hz *)
idxC = First[Position[fHzList, _?(# >= 1*^9 &), 1, 1]][[1]];
Print[""];
Print["C-axis asymptote ratio r(1e9):"];
Print["  CLN_N3       = ", N[results[[idxC, 3]] / results[[idxC, 6]], 4]];
Print["  CLN_N10      = ", N[results[[idxC, 4]] / results[[idxC, 6]], 4]];
Print["  Schur-F b4   = ", N[results[[idxC, 5]] / results[[idxC, 6]], 4],
      "    (target 1)"];

Print[""];
Print["=== rect2D_schur_F.wls DONE ==="];
