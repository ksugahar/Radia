#!/usr/bin/env wolframscript
(* ::Package:: *)

(*
  bem_cln_2cuboid_3D.wls

  Phase 3 Task B foundation: rigorous 3D 2-cuboid BEM-CLN assembly.
  Extends the Phase 2.5 polarizability framework from 2D infinite
  cylinders to 3D finite cuboids, using the cuboid Foster eigenmode
  decomposition (cuboid3D_schur_F.wls building block) and the 3D
  dipole-dipole Green's function.

  Per-cuboid building block (scalar diffusion model, V.B journal):
    Y_cuboid(s) = sigma * Sum_{m,n,p odd} g^2_mnp / (1 + s tau_mnp)
    g^2_mnp     = 512 a b c / ((mnp)^2 pi^6)
    lambda_mnp  = (m pi/a)^2 + (n pi/b)^2 + (p pi/c)^2
    tau_mnp     = mu sigma / lambda_mnp

  Schur-F per-cuboid:
    Y_cuboid_S(s) = Y_cuboid_N(s) + K_SIBC^(cuboid) sqrt(s) / (s + d_param)
    K_SIBC^(cuboid) = 2 (ab + bc + ca) sqrt(sigma/mu)

  Polarizability for 3D cuboid (volume susceptibility):
    alpha_3D(s) = V * (1 - Y_cuboid(s) / (sigma V))
                = V - Y_cuboid(s) / sigma
  Limits:
    alpha(0) = 0           (DC, no eddy currents)
    alpha(inf) = V         (PEC, full volume excluded)

  3D BEM-CLN coupling (3D dipole-dipole Green's function):
    A point dipole m_y at origin produces dB_y at distance D along x:
      dB_y = -(mu_0 / (4 pi D^3)) m_y
    For two cuboids n=1, 2 with c_n = alpha(s) B^local_n:
      B^local_n = B_0 - sum_{m != n} (mu_0 / (4 pi D_nm^3)) c_m
    System:
      (I + alpha(s) G_3D) c = alpha(s) B_0 1_vec
      G_ij = mu_0 / (4 pi D_ij^3)  for i != j, 0 on diagonal
    Note the +sign of alpha G (opposite to 2D's -sign), reflecting
    the destructive 3D dipole-dipole interaction (3 cos^2 theta - 1
    angular factor).

  Total polarizability: alpha_total = sum c_n / B_0.

  Run: wolframscript -file bem_cln_2cuboid_3D.wls
*)

ClearAll[a, b, c, sigma, mu, dParam, basisN];

a       = 5*^-3;
b       = 2*^-3;
c       = 1*^-3;
sigma   = 5.80*^7;
mu      = 4 Pi 1*^-7;
dParam  = 2 Pi 1*^5;
basisN  = 3;

VCub    = a b c;
SCub    = 2 (a b + b c + c a);
KsibcCub = SCub Sqrt[sigma/mu];

Print["=== bem_cln_2cuboid_3D.wls ==="];
Print["3D 2-cuboid BEM-CLN: rigorous polarizability framework"];
Print["Cuboid: ", a 1000, "x", b 1000, "x", c 1000, " mm Cu"];
Print["  V = ", N[VCub, 6], " m^3"];
Print["  S = ", N[SCub, 6], " m^2"];
Print["  K_SIBC^(cuboid) = S sqrt(sigma/mu) = ", N[KsibcCub, 6]];
Print[""];

(* === Per-cuboid Foster eigenmode building block ======================== *)
Mmax = 13;
modeList = Flatten[Table[{m, n, p},
  {m, 1, Mmax, 2}, {n, 1, Mmax, 2}, {p, 1, Mmax, 2}], 2];
Nmodes = Length[modeList];

g2MNP[m_, n_, p_] := 512 a b c / ((m n p)^2 Pi^6);
lambdaMNP[m_, n_, p_] := (m Pi/a)^2 + (n Pi/b)^2 + (p Pi/c)^2;
tauMNP[m_, n_, p_]    := mu sigma / lambdaMNP[m, n, p];

g2Vec  = N[Table[g2MNP[mn[[1]], mn[[2]], mn[[3]]], {mn, modeList}], 30];
tauVec = N[Table[tauMNP[mn[[1]], mn[[2]], mn[[3]]], {mn, modeList}], 30];

(* Foster Y_cuboid(s) *)
YcubFoster[s_] := sigma Total[g2Vec / (1 + s tauVec)];

(* Schur-F Y_cuboid *)
YcubSchur[s_] := YcubFoster[s] + KsibcCub Sqrt[s] / (s + dParam);

(* Polarizability *)
alphaCubFoster[s_] := VCub - YcubFoster[s] / sigma;
alphaCubSchur[s_]  := VCub - YcubSchur[s] / sigma;

(* Sanity: DC and PEC limits *)
Print["Per-cuboid alpha(s) sanity:"];
Print["  Y_cuboid(0) Foster   = ", N[YcubFoster[0], 8],
      "  (target sigma V = ", N[sigma VCub, 8], ")"];
Print["  alpha(0)             = ", N[Re[alphaCubFoster[I 2 Pi 0.01]], 8]];
Print["  alpha(inf) Foster    = ", N[Re[alphaCubFoster[I 2 Pi 1*^14]], 8],
      "  (target V = ", N[VCub, 8], ")"];
Print["  alpha(inf) Schur     = ", N[Re[alphaCubSchur[I 2 Pi 1*^14]], 8]];
Print[""];

(* === 3D BEM-CLN coupling =============================================== *)
(* G_ij = mu_0 / (4 pi D_ij^3) for i != j *)
gCoupling3D[D_] := mu / (4 Pi D^3);

coupling3DMatrix[Nc_, D_] := Table[
  If[i == j, 0, gCoupling3D[Abs[i - j] D]],
  {i, 1, Nc}, {j, 1, Nc}];

solveBemCLN3D[s_, Nc_, D_, alphaFun_] := Module[
  {alf, G, sys, rhs, cvec, ones},
  alf  = alphaFun[s];
  G    = coupling3DMatrix[Nc, D];
  ones = ConstantArray[1, Nc];
  sys  = IdentityMatrix[Nc] + alf G;     (* I + alpha G (3D sign) *)
  rhs  = alf ones;
  cvec = LinearSolve[sys, rhs];
  Total[cvec]
];

alphaTotal3DFoster[s_, Nc_, D_] := solveBemCLN3D[s, Nc, D, alphaCubFoster];
alphaTotal3DSchur[s_, Nc_, D_]  := solveBemCLN3D[s, Nc, D, alphaCubSchur];

(* === Verification: 2-cuboid pair at various D ========================== *)
Print["3D 2-cuboid BEM-CLN: alpha_total(s, D) at f = 10^6 Hz"];
Print[""];
Print[StringPadRight["D [mm]", 10] <>
      StringPadRight["alpha (per cuboid)", 22] <>
      StringPadRight["alpha_total (Foster)", 22] <>
      "alpha_total/(2 alpha)"];
Print[StringRepeat["-", 78]];

With[{omegaProbe = 2 Pi 1*^6},
  With[{sProbe = I omegaProbe},
    With[{alfSingle = Re[alphaCubFoster[sProbe]]},
      Do[
        With[{Dval = Dmm 1*^-3},
          With[{aTot = Re[alphaTotal3DFoster[sProbe, 2, Dval]]},
            Print[StringPadRight[ToString[N[Dmm, 4]], 10],
                  StringPadRight[ToString[N[alfSingle, 6]], 22],
                  StringPadRight[ToString[N[aTot, 6]], 22],
                  N[aTot/(2 alfSingle), 8]]
          ]],
        {Dmm, {5, 10, 15, 30, 100, 1000}}]
    ]]];

Print[""];

(* === DC and PEC sanity check =========================================== *)
Print["Sanity checks at D = 15 mm:"];
With[{Dprobe = 15*^-3, Nc = 2},
  With[{s0 = I 2 Pi 0.01, sInf = I 2 Pi 1*^14},
    With[{alphaDC   = Re[alphaTotal3DFoster[s0, Nc, Dprobe]],
          alphaPEC  = Re[alphaTotal3DFoster[sInf, Nc, Dprobe]],
          alphaSingleInf = Re[alphaCubFoster[sInf]]},
      Print["  DC limit alpha_total(s=0)            = ",
            N[alphaDC, 8], "  (target 0)"];
      Print["  PEC limit alpha_total(s=infinity)    = ",
            N[alphaPEC, 8]];
      Print["  N alpha_single (PEC, uncoupled)      = ",
            N[Nc alphaSingleInf, 8], "  (target ~V)"];
      Print["  alpha_total / N alpha_single (PEC)   = ",
            N[alphaPEC/(Nc alphaSingleInf), 8],
            "  (= 1 / (1 + V mu_0/(4 pi D^3) coupling))"]
    ]]];

(* === Compare Foster vs Schur-F per-cuboid in 3D BEM-CLN ================ *)
Print[""];
Print["Foster vs Schur-F per-cuboid building blocks in 2-cuboid BEM-CLN:"];
Print[""];
Print[StringPadRight["f [Hz]", 10] <>
      StringPadRight["alpha_tot Foster", 22] <>
      StringPadRight["alpha_tot Schur", 22] <>
      "ratio"];
Print[StringRepeat["-", 70]];

With[{Dprobe = 15*^-3},
  Do[
    With[{omega = 2 Pi 10.0^logf},
      With[{s = I omega},
        With[{aFos = Re[alphaTotal3DFoster[s, 2, Dprobe]],
              aSch = Re[alphaTotal3DSchur[s, 2, Dprobe]]},
          Print[StringPadRight[ToString["1e" <> ToString[logf]], 10],
                StringPadRight[ToString[N[aFos, 6]], 22],
                StringPadRight[ToString[N[aSch, 6]], 22],
                N[aSch/aFos, 6]]
        ]]],
    {logf, {3, 4, 5, 6, 7, 8, 10, 12}}]
];

(* === Far-field test ===================================================== *)
Print[""];
With[{Dfar = 100, sProbe = I 2 Pi 1*^6, Nc = 2},
  With[{alfSingle = Re[alphaCubFoster[sProbe]],
        aTotFar = Re[alphaTotal3DFoster[sProbe, Nc, Dfar]]},
    Print["Far-field D = 100 m: alpha_total/(N alpha_single) = ",
          N[aTotFar/(Nc alfSingle), 12], "  (target 1.0)"]
  ]];

(* === CSV output ========================================================= *)
Print[""];
Print["Writing CSV..."];
csvData = {{"f_Hz", "D_mm", "alpha_per_cuboid", "alpha_total_Foster", "alpha_total_Schur"}};
Dout = 15*^-3;
Do[
  With[{omega = 2 Pi (10.0^logf), s = I 2 Pi (10.0^logf)},
    With[{alf = Re[alphaCubFoster[s]],
          aFos = Re[alphaTotal3DFoster[s, 2, Dout]],
          aSch = Re[alphaTotal3DSchur[s, 2, Dout]]},
      AppendTo[csvData,
        {10.0^logf, 15.0, N[alf, 8], N[aFos, 8], N[aSch, 8]}]
    ]],
  {logf, {2, 3, 4, 5, 6, 7, 8, 10, 12}}];

Export["bem_cln_2cuboid_3D.csv", csvData, "CSV"];
Print["Wrote bem_cln_2cuboid_3D.csv"];

Print[""];
Print["=== bem_cln_2cuboid_3D.wls DONE ==="];
