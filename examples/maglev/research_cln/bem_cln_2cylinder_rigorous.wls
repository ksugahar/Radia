#!/usr/bin/env wolframscript
(* ::Package:: *)

(*
  bem_cln_2cylinder_rigorous.wls

  Phase 2.5 rigorous derivation of the inter-element coupling for the
  BEM-CLN multi-cylinder assembly.  The phenomenological saturating
  g(s, D) of Phase 1 is replaced by the standard 2D dipole-dipole
  Green's function 1/D^2, with the per-cylinder building block
  re-expressed as a polarizability alpha(s) rather than the
  IGTE-style admittance Y_cyl(s).  The polarizability vanishes at
  DC and saturates at high frequency by construction, so the assembly
  is rigorous without any phenomenological factor.

  Single-cylinder polarizability derivation (Bessel exact)
  ---------------------------------------------------------
  For an infinite Cu cylinder of radius a in uniform transverse
  applied field B_0 y_hat, the matching of A_z at r = a between
    Helmholtz inside: A_z = D I_1(gamma r) cos(theta)
    Laplace outside: A_z = -B_0 r cos(theta) + c cos(theta)/r
  yields
    c(s) = a^2 [1 - 2 I_1(gamma a) / (gamma a I_0(gamma a))] B_0,
    gamma = sqrt(s mu sigma).
  We DEFINE the polarizability alpha(s) by  c(s) = alpha(s) B_0:
    alpha(s) = a^2 [1 - 2 I_1(gamma a) / (gamma a I_0(gamma a))].
  Limits:
    alpha(0) = a^2 [1 - 1] = 0         (no eddy current at DC)
    alpha(infinity) = a^2 [1 - 0] = a^2 (PEC diamagnetic shielding limit)
  Crucially, alpha(s) ranges in [0, a^2] for all real-positive s, so
  any coupling matrix g_ij ~ 1/D_ij^2 acting on alpha is bounded:
  alpha/D^2 <= a^2/D^2 < 1 provided D > a.  No phenomenological
  saturation factor is needed.

  Relation to the IGTE-style Y_cyl(s)
  ------------------------------------
    Y_cyl(s) = pi a^2 sigma * 2 I_1(gamma a) / (gamma a I_0(gamma a))
             = pi a^2 sigma * (1 - alpha(s)/a^2)
  Equivalently:
    alpha(s) = a^2 [1 - Y_cyl(s)/(pi a^2 sigma)]
            = a^2 - Y_cyl(s) / (pi sigma).
  This conversion lets a Schur-F-augmented Y_cyl^Schur(s) build a
  Schur-F-augmented polarizability alpha^Schur(s) automatically.

  Coupling: 2D dipole-dipole field at distance D
  -----------------------------------------------
  A 2D dipole c*cos(theta)/r at the origin produces a perturbation
  of the transverse B_y component at position (D, 0):
    Delta B_y = c / D^2.
  So the local field at cylinder n is
    B_y^local_n = B_0 + sum_{m != n} c_m / D_{nm}^2.
  With c_n = alpha(s) B_y^local_n, this gives
    c_n - alpha(s) sum_{m != n} c_m / D_{nm}^2 = alpha(s) B_0,
  i.e.
    (I - alpha(s) * G_{nm}) . vec{c} = alpha(s) B_0 vec{1},
    G_{nm} = 1/D_{nm}^2 (off-diagonal), G_{nn} = 0.
  The total induced polarizability of the N-cylinder system:
    alpha_total(s, D) = sum_n c_n / B_0.

  Run: wolframscript -file bem_cln_2cylinder_rigorous.wls
*)

ClearAll[a, sigma, mu, dParam, basisN];

a       = 5*^-3;
sigma   = 5.96*^7;
mu      = 4 Pi 1*^-7;
dParam  = 2 Pi 1*^5;
basisN  = 4;

Print["=== bem_cln_2cylinder_rigorous.wls ==="];
Print["Rigorous polarizability coupling (no phenomenological factor)"];
Print["a = ", a 1000, " mm, sigma = ", sigma, " S/m"];
Print[""];

KsibcCyl = 2 Pi a Sqrt[sigma/mu];

(* === Per-cylinder admittance + polarizability ========================== *)
YcylBessel[s_] := Module[{gamma},
  gamma = Sqrt[s mu sigma];
  Pi a^2 sigma * 2 BesselI[1, gamma a] / (gamma a BesselI[0, gamma a])];

alphaBessel[s_] := a^2 - YcylBessel[s] / (Pi sigma);

alphaModes = N[BesselJZero[0, Range[1, 60]], 30];
gVec2      = 4 / alphaModes^2;
tauVec     = mu sigma a^2 / alphaModes^2;

YcylFosterN[s_, n_] := Pi a^2 sigma Total[gVec2[[1 ;; n]] / (1 + s tauVec[[1 ;; n]])];

YcylSchur[s_, n_] := YcylFosterN[s, n] + KsibcCyl Sqrt[s] / (s + dParam);

alphaSchur[s_, n_] := a^2 - YcylSchur[s, n] / (Pi sigma);

(* === Sanity: alpha limits =============================================== *)
Print["alpha(s) sanity checks:"];
Print["  alpha(s -> 0+)            = ", N[Re[alphaBessel[I 2 Pi 0.01]], 6],
      " (target 0)"];
Print["  alpha(s -> infinity)      = ", N[Re[alphaBessel[I 2 Pi 1*^14]], 6],
      " (target a^2 = ", N[a^2, 6], ")"];
Print["  alpha(s) Schur at f=0     = ", N[Re[alphaSchur[I 2 Pi 0.01, basisN]], 6]];
Print["  alpha(s) Schur at f=1e12  = ", N[Re[alphaSchur[I 2 Pi 1*^12, basisN]], 6]];
Print[""];

(* === N-cylinder rigorous BEM-CLN solver ================================= *)
couplingMatrix[Nc_, D_] := Table[
  If[i == j, 0, 1/(Abs[i - j] D)^2],
  {i, 1, Nc}, {j, 1, Nc}];

solveBemCLN[s_, Nc_, D_, alphaFun_] := Module[
  {alf, G, sys, rhs, cvec, ones},
  alf = alphaFun[s];
  G   = couplingMatrix[Nc, D];
  ones = ConstantArray[1, Nc];
  sys  = IdentityMatrix[Nc] - alf G;     (* (I - alpha G) *)
  rhs  = alf ones;
  cvec = LinearSolve[sys, rhs];
  Total[cvec]                              (* sum of induced c_n / B_0 *)
];

alphaTotalBessel[s_, Nc_, D_] := solveBemCLN[s, Nc, D, alphaBessel];
alphaTotalSchur[s_, Nc_, D_]  := solveBemCLN[s, Nc, D, alphaSchur[#, basisN] &];

(* === Verification table: N = 2 to 10 at D = 15 mm ====================== *)
Print["Verification: linear chain, D = 15 mm, f = 10^6 Hz"];
Print[""];
Print[StringPadRight["N", 4] <>
      StringPadRight["DOF Schur", 12] <>
      StringPadRight["alpha_total (Bessel)", 24] <>
      StringPadRight["alpha_total (Schur)", 24] <>
      "Schur/Bessel"];
Print[StringRepeat["-", 80]];

With[{Dprobe = 15*^-3, fProbe = 1*^6, omegaProbe = 2 Pi 1*^6},
  With[{sProbe = I omegaProbe},
    Do[
      With[{aBes = Re[alphaTotalBessel[sProbe, Nc, Dprobe]],
            aSch = Re[alphaTotalSchur[sProbe, Nc, Dprobe]]},
        Print[StringPadRight[ToString[Nc], 4],
              StringPadRight[ToString[Nc (basisN + 1)], 12],
              StringPadRight[ToString[N[aBes, 6]], 24],
              StringPadRight[ToString[N[aSch, 6]], 24],
              N[aSch / aBes, 6]]
      ],
      {Nc, {2, 3, 5, 10}}]
  ]];

(* === DC and PEC limits ================================================== *)
Print[""];
Print["DC limit (s -> 0):  alpha_total(0) should be 0 (no eddy current)"];
With[{Nc = 5, Dprobe = 15*^-3, s0 = I 2 Pi 0.001},
  Print["  alpha_total(s=0, N=5, D=15mm) = ", N[Re[alphaTotalBessel[s0, Nc, Dprobe]], 8]];
];

Print[""];
Print["High-frequency PEC limit:  alpha -> N a^2 / (1 - sum coupling)"];
With[{Nc = 5, Dprobe = 15*^-3, sInf = I 2 Pi 1*^14},
  With[{alphaPECsingle = a^2,
        Gmat = couplingMatrix[Nc, Dprobe]},
    With[{eigvalsG = Sort[Eigenvalues[Gmat], Greater]},
      Print["  largest eigenvalue of G = ", N[eigvalsG[[1]], 6]];
      Print["  a^2 * (largest eigval)   = ", N[alphaPECsingle eigvalsG[[1]], 6],
            " (must be < 1 for stability)"];
      Print["  alpha_total(s -> infty)  = ",
            N[Re[alphaTotalBessel[sInf, Nc, Dprobe]], 8]]
    ]]];

Print[""];
Print["Far-field (D -> infty): alpha_total -> N alpha_single"];
With[{Nc = 5, Dfar = 2*^5, sProbe = I 2 Pi 1*^6},
  With[{alf = Re[alphaBessel[sProbe]],
        aTot = Re[alphaTotalBessel[sProbe, Nc, Dfar]]},
    Print["  N alpha_single at f=1e6 = ", N[Nc alf, 8]];
    Print["  alpha_total at D=200m   = ", N[aTot, 8]];
    Print["  ratio                    = ", N[aTot/(Nc alf), 8]]
  ]];

(* === SIBC tail check ==================================================== *)
Print[""];
Print["SIBC convergence (N=5, D=15mm): alpha approaches PEC at high s"];
Print["  Expected at high s: alpha_single -> a^2; alpha_total -> ..."];
Print[""];
Print[StringPadRight["f [Hz]", 12] <>
      StringPadRight["alpha_single", 20] <>
      "alpha_total / (N a^2)"];
Print[StringRepeat["-", 65]];

With[{Nc5 = 5, Dprobe = 15*^-3},
  Do[
    With[{omega = 2 Pi 10.0^logf},
      With[{s = I omega},
        With[{alf  = Re[alphaBessel[s]],
              aTot = Re[alphaTotalBessel[s, Nc5, Dprobe]]},
          Print[StringPadRight[ToString["f=1e" <> ToString[logf]], 12],
                StringPadRight[ToString[N[alf, 6]], 20],
                N[aTot / (Nc5 a^2), 6]]
        ]]],
    {logf, {2, 4, 6, 8, 10, 12, 14}}]
];

(* === Convert back to admittance scale for comparison with Phase 1 ======= *)
Print[""];
Print["Conversion: total admittance Y_total = pi sigma (N a^2 - alpha_total)"];
With[{Nc = 2, Dprobe = 15*^-3},
  Do[
    With[{omega = 2 Pi 10.0^logf},
      With[{s = I omega},
        With[{aTotBes = Re[alphaTotalBessel[s, Nc, Dprobe]],
              aTotSch = Re[alphaTotalSchur[s, Nc, Dprobe]]},
          With[{yTotBes = Pi sigma (Nc a^2 - aTotBes),
                yTotSch = Pi sigma (Nc a^2 - aTotSch)},
            Print["f=1e", logf, " Hz: |Y_total| Bessel = ", N[yTotBes, 5],
                  ", Schur = ", N[yTotSch, 5]]
          ]]]],
    {logf, {2, 4, 6, 8, 10, 12}}]
];

(* === Write CSV ========================================================== *)
Print[""];
Print["Writing rigorous BEM-CLN CSV..."];
freqList = Table[10^k, {k, 2, 14, 0.5}];
NList = {2, 3, 5, 10};

csvData = {{"N", "f_Hz", "alpha_total_Bessel", "alpha_total_Schur",
            "Y_total_Bessel", "Y_total_Schur"}};
Do[
  Do[
    With[{omega = 2 Pi freq, s = I 2 Pi freq},
      With[{aTotBes = Re[alphaTotalBessel[s, Nc, 15*^-3]],
            aTotSch = Re[alphaTotalSchur[s, Nc, 15*^-3]]},
        AppendTo[csvData,
          {Nc, N[freq, 6], N[aTotBes, 8], N[aTotSch, 8],
           N[Pi sigma (Nc a^2 - aTotBes), 8],
           N[Pi sigma (Nc a^2 - aTotSch), 8]}]
      ]]
    , {freq, freqList}]
  , {Nc, NList}];

Export["bem_cln_rigorous.csv", csvData, "CSV"];
Print["Wrote bem_cln_rigorous.csv (", (Length[csvData] - 1), " rows)"];

Print[""];
Print["=== bem_cln_2cylinder_rigorous.wls DONE ==="];
