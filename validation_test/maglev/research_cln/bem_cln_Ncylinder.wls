#!/usr/bin/env wolframscript
(* ::Package:: *)

(*
  bem_cln_Ncylinder.wls

  Phase 2: N-cylinder sparse-block BEM-CLN assembly.  Generalises the
  Phase 1 two-cylinder demo (bem_cln_2cylinder.wls) to an arbitrary
  number N of identical infinite Cu cylinders arranged in a linear
  chain at equal spacing D.  The applied field is uniform transverse
  B_0 in y.

  System (N x N) at each frequency s:
    (I + diag(Y_cyl(s)) . G(s)) . m = Y_cyl(s) . B_0 . 1
  where
    Y_cyl(s)     scalar per-cylinder admittance (Bessel or Schur-F);
    G_ij(s) = (mu_0 / (2 pi |r_i - r_j|^2)) * s / (s + s_sat) for i != j;
    G_ii = 0;
    m_n          induced dipole moment of cylinder n.

  Total system admittance:
    Y_total(s, N, D) = (1 / B_0) sum_n m_n(s).

  DOF count (Schur-F per cylinder, basis N_C):
    total = N * (N_C + 1).      Memory scales O(N^2) for coupling matrix.

  Verifies:
   (a) DC limit: Y_total(0, N, D) -> N * Y_cyl(0)
   (b) Far-field limit (D -> infinity): Y_total -> N * Y_cyl
   (c) SIBC tail: |Y_total(j omega)| sqrt(omega) -> N * K_SIBC^(cyl)
       as omega -> infinity (via per-cylinder Schur-F)
   (d) Bessel-exact vs Schur-F-basis-4 building blocks: agreement in wall band

  Run: wolframscript -file bem_cln_Ncylinder.wls
*)

ClearAll[a, sigma, mu, dParam, basisN, omegaSat, sSat];

a        = 5*^-3;
sigma    = 5.96*^7;
mu       = 4 Pi 1*^-7;
dParam   = 2 Pi 1*^5;
basisN   = 4;
omegaSat = 2 Pi 1*^5;
sSat     = I omegaSat;

Print["=== bem_cln_Ncylinder.wls ==="];
Print["N-cylinder linear-chain BEM-CLN assembly"];
Print["a = ", a 1000, " mm, sigma = ", sigma, " S/m, Schur-F basis ", basisN];
Print[""];

KsibcCyl = 2 Pi a Sqrt[sigma/mu];
Print["K_SIBC^(cyl) = ", N[KsibcCyl, 6]];

(* === Per-cylinder building blocks ====================================== *)
YcylBessel[s_] := Module[{gamma},
  gamma = Sqrt[s mu sigma];
  Pi a^2 sigma * 2 BesselI[1, gamma a] / (gamma a BesselI[0, gamma a])];

alphaBessel = N[BesselJZero[0, Range[1, 60]], 30];
gVec2       = 4 / alphaBessel^2;
tauVec      = mu sigma a^2 / alphaBessel^2;

YcylFosterN[s_, n_] := Pi a^2 sigma Total[gVec2[[1 ;; n]] / (1 + s tauVec[[1 ;; n]])];

YcylSchur[s_, n_] := YcylFosterN[s, n] + KsibcCyl Sqrt[s] / (s + dParam);

(* === Coupling matrix G(s) for linear chain of N cylinders ============== *)
(* Cylinder n at position (n - 1) * D along x.  Distance D_ij = |i - j| D. *)
gPair[s_, dist_] := (mu / (2 Pi dist^2)) (s / (s + sSat));

couplingMatrix[s_, Nc_, D_] := Table[
  If[i == j, 0, gPair[s, Abs[i - j] D]],
  {i, 1, Nc}, {j, 1, Nc}];

(* === Solve N x N BEM-CLN system at frequency s ========================= *)
solveBemCLN[s_, Nc_, D_, ycylFun_] := Module[
  {yc, G, sys, rhs, mvec, ones},
  yc = ycylFun[s];
  G = couplingMatrix[s, Nc, D];
  ones = ConstantArray[1, Nc];
  sys = IdentityMatrix[Nc] + yc G;
  rhs = yc ones;
  mvec = LinearSolve[sys, rhs];
  Total[mvec]
];

YtotalBemBessel[s_, Nc_, D_] := solveBemCLN[s, Nc, D, YcylBessel];
YtotalBemSchur[s_, Nc_, D_]  := solveBemCLN[s, Nc, D, YcylSchur[#, basisN] &];

(* === Verification table: N = 2, 3, 5, 10 at D = 15 mm ================== *)
Print[""];
Print["Verification: linear chain, D = 15 mm, at f = 10^6 Hz wall band"];
Print[""];
Print[StringPadRight["N", 4] <>
      StringPadRight["DOF (Schur)", 14] <>
      StringPadRight["|Y_tot Bessel|", 18] <>
      StringPadRight["|Y_tot Schur|", 18] <>
      StringPadRight["Schur/Bessel", 14] <>
      "N*Y_cyl_uncoupled"];
Print[StringRepeat["-", 90]];

With[{Dprobe = 15*^-3, fProbe = 1*^6, omegaProbe = 2 Pi 1*^6},
  With[{sProbe = I omegaProbe},
    With[{ycylUncoupled = Abs[YcylBessel[sProbe]]},
      Do[
        With[{ytBes  = Abs[YtotalBemBessel[sProbe, Nc, Dprobe]],
              ytSch  = Abs[YtotalBemSchur[sProbe, Nc, Dprobe]]},
          Print[StringPadRight[ToString[Nc], 4],
                StringPadRight[ToString[Nc (basisN + 1)], 14],
                StringPadRight[ToString[N[ytBes, 5]], 18],
                StringPadRight[ToString[N[ytSch, 5]], 18],
                StringPadRight[ToString[N[ytSch / ytBes, 6]], 14],
                N[Nc ycylUncoupled, 5]]
        ],
        {Nc, {2, 3, 5, 10}}]
    ]]];

(* === SIBC tail check for N = 5 chain =================================== *)
Print[""];
Print["SIBC tail (N = 5, D = 15 mm): r(omega) = |Y_total| sqrt(omega) / (N K_SIBC^cyl)"];
Print["  Expected r -> 1 at high omega via per-cylinder Schur-F"];
Print[""];
Print[StringPadRight["f [Hz]", 10] <>
      StringPadRight["r (Bessel)", 14] <>
      "r (Schur-F)"];
Print[StringRepeat["-", 50]];

With[{Nc5 = 5, Dprobe = 15*^-3, Ktot5 = 5 KsibcCyl},
  Do[
    With[{omega = 2 Pi (10.0^logf)},
      With[{s = I omega},
        With[{yBes  = Abs[YtotalBemBessel[s, Nc5, Dprobe]],
              ySchur = Abs[YtotalBemSchur[s, Nc5, Dprobe]]},
          Print["f = 1e", logf, " Hz   ",
                "r_Bessel = ", N[yBes Sqrt[omega] / Ktot5, 5], "   ",
                "r_Schur = ",  N[ySchur Sqrt[omega] / Ktot5, 5]]
        ]]],
    {logf, {4, 5, 6, 7, 8, 10, 12}}]
];

(* === DC and parallel limits for sanity check =========================== *)
Print[""];
Print["Sanity: DC and parallel limits"];
With[{omega0 = 2 Pi 0.01, omegaInfD = 2 Pi 1*^6, Nc = 5, Dprobe = 15*^-3, Dfar = 2*^5},
  With[{s0 = I omega0, sInfD = I omegaInfD},
    With[{ycyl0 = Abs[Re[YcylBessel[s0]]],
          ytotDC = Abs[Re[YtotalBemBessel[s0, Nc, Dprobe]]],
          ytotFar = Abs[YtotalBemBessel[sInfD, Nc, Dfar]],
          ycylWall = Abs[YcylBessel[sInfD]]},
      Print["  DC: Y_total(0) / (N Y_cyl(0))  = ", N[ytotDC / (Nc ycyl0), 8],
            "  (target 1.000000)"];
      Print["  Far-field D=200m: Y_tot / (N Y_cyl) at f=1e6 = ",
            N[ytotFar / (Nc ycylWall), 8],
            "  (target 1.000000)"]
    ]]];

(* === Sweep CSV output =================================================== *)
Print[""];
Print["Writing N-scaling sweep CSV..."];
Dout = 15*^-3;
freqList = Table[10^k, {k, 4, 12, 1}];
NList = {2, 3, 5, 10};

csvData = {{"N", "f_Hz", "Y_tot_Bessel_mag", "Y_tot_Schur_mag",
            "DOF_total_Schur", "N_KsibcCyl_div_sqrt_omega"}};
Do[
  Do[
    With[{omega = 2 Pi freq, s = I 2 Pi freq},
      With[{yBes = Abs[YtotalBemBessel[s, Nc, Dout]],
            ySchur = Abs[YtotalBemSchur[s, Nc, Dout]],
            asympt = Nc KsibcCyl / Sqrt[omega]},
        AppendTo[csvData,
          {Nc, N[freq, 6], N[yBes, 8], N[ySchur, 8],
           Nc (basisN + 1), N[asympt, 8]}]
      ]]
    , {freq, freqList}]
  , {Nc, NList}];

Export["bem_cln_Ncylinder.csv", csvData, "CSV"];
Print["Wrote bem_cln_Ncylinder.csv (", (Length[csvData] - 1), " rows)"];

Print[""];
Print["=== bem_cln_Ncylinder.wls DONE ==="];
