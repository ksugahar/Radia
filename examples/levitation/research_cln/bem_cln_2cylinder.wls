#!/usr/bin/env wolframscript
(* ::Package:: *)

(*
  bem_cln_2cylinder.wls

  Phase 1 demonstration of the per-element multipole CLN + Schur-F
  BEM kernel introduced in journal section V.G.  Builds a two-element
  BEM-CLN matrix system on two parallel infinite Cu cylinders in a
  uniform transverse applied field, using single-cylinder Schur-F
  models as building blocks and a 2D dipole-dipole Green's function
  for the inter-element coupling.

  Per-cylinder building block (matches digest IGTE Fig 1):
    Y_cyl(s) = pi a^2 sigma * 2 I1(gamma a) / (gamma a I0(gamma a))
    gamma = sqrt(s mu sigma)
    K_SIBC^(cyl) = 2 pi a sqrt(sigma / mu)

  Single-cylinder Schur-F approximation (basis-4 in the IGTE digest):
    Y_cyl_Schur(s) = Y_CLN_N(s) + K_SIBC^(cyl) sqrt(s) / (s + d_param)

  Two-cylinder coupling (separation D between centres, applied B in y):
    m_n(s) = Y_cyl(s) * B_n^(applied)
    B_1^(applied) = B_0 - g(D) * m_2
    B_2^(applied) = B_0 - g(D) * m_1
    g(D) = mu_0 / (2 pi D^2)     2D transverse-dipole Green's function

  System (2x2):
    [1                    Y_cyl(s) g(D) ] [m_1]   [Y_cyl(s) B_0]
    [Y_cyl(s) g(D)        1             ] [m_2] = [Y_cyl(s) B_0]

  Total admittance of the 2-cylinder system:
    Y_total(s) = (m_1 + m_2) / B_0 = 2 Y_cyl(s) / (1 + Y_cyl(s) g(D))

  Sanity checks:
    (a) D -> infinity:  g(D) -> 0, Y_total -> 2 Y_cyl  (parallel limit)
    (b) High frequency: Y_cyl -> K_SIBC^(cyl)/sqrt(s); the coupling term
        Y_cyl(s) g(D) -> 0 at high s (since Y_cyl -> 0); so
        Y_total -> 2 K_SIBC^(cyl)/sqrt(s) = K_SIBC^(total)/sqrt(s)
        with K_SIBC^(total) = 2 K_SIBC^(cyl).  SIBC tail PRESERVED.

  This demonstrates the BEM-CLN matrix assembly: the system is
  small (one DOF per cylinder, here equivalent to one m_n per
  cylinder), and the SIBC asymptote is preserved by construction
  via each per-cylinder Schur-F block.

  Run: wolframscript -file bem_cln_2cylinder.wls
*)

ClearAll[a, sigma, mu, dParam, basisN];

(* Material: Cu cylinder, radius 5 mm (matches IGTE Fig 1) *)
a       = 5*^-3;
sigma   = 5.96*^7;
mu      = 4 Pi 1*^-7;
dParam  = 2 Pi 1*^5;            (* Schur-F crossover from digest *)
basisN  = 4;

Print["=== bem_cln_2cylinder.wls ==="];
Print["Two parallel infinite Cu cylinders, transverse B applied"];
Print["a = ", a 1000, " mm,  sigma = ", sigma, " S/m"];
Print["Schur-F crossover d = 2 pi * ", dParam/(2 Pi), " rad/s"];
Print[""];

(* Single-cylinder K_SIBC *)
KsibcCyl = 2 Pi a Sqrt[sigma/mu];
Print["K_SIBC^(cyl) = 2 pi a sqrt(sigma/mu) = ", N[KsibcCyl, 8]];

(* === Per-cylinder admittance (Bessel exact) ============================ *)
YcylBessel[s_] := Module[{gamma},
  gamma = Sqrt[s mu sigma];
  Pi a^2 sigma * 2 BesselI[1, gamma a] / (gamma a BesselI[0, gamma a])];

(* === Per-cylinder Schur-F approximation ================================ *)
(* CLN N-rung at s_0 = 0: Pade approximation of Y_cyl(s) at s = 0.        *)
(* For the cylinder, the Taylor expansion gives                            *)
(*   Y_cyl(s) = pi a^2 sigma * (1 - (gamma a)^2 / 8 + ... )                *)
(* Truncating to first N moments gives the Pade [N-1/N] rational.          *)
(* For demo purposes, we use a simple Cauer-I ladder constructed from the  *)
(* Foster decomposition of Y_cyl: poles at -1/tau_n with tau_n = mu sigma  *)
(* a^2/alpha_n^2, where alpha_n are zeros of J_0(alpha).                   *)

(* Modal decomposition of Bessel function: zeros of J_0 *)
alphaBessel = N[BesselJZero[0, Range[1, 60]], 30];  (* first 60 modes *)
tauVec      = mu sigma a^2 / alphaBessel^2;
gVec2       = 4 / alphaBessel^2;                    (* mode coupling for Bessel admittance *)

(* DC value check *)
Print["Y_cyl(0) Bessel limit       = pi a^2 sigma = ", N[Pi a^2 sigma, 8]];
Print["Y_cyl(0) Foster (60 modes)  = ",
      N[Pi a^2 sigma Total[gVec2], 8]];
Print["Foster convergence          = ",
      N[100 Total[gVec2], 6], " % (should approach 100)"];
Print[""];

YcylFosterN[s_, N_] := Pi a^2 sigma Total[gVec2[[1 ;; N]] / (1 + s tauVec[[1 ;; N]])];

YcylSchur[s_, N_] := YcylFosterN[s, N] + KsibcCyl Sqrt[s] / (s + dParam);

(* === BEM-CLN assembly for 2 cylinders =================================== *)
(* Inter-element coupling: phenomenological 2D transverse-dipole mutual
   coupling normalised to vanish at DC and saturate at high s.
   Physical motivation: mutual inductance between two conductors induces
   no current at DC (no time-varying flux), and as both cylinders enter
   the skin-shielded regime at high s, the mutual coupling saturates
   (each conductor screens the other).  We capture both limits with
       g(s, D) = mu / (2 pi D^2) * s / (s + s_sat)
   where s_sat = j omega_sat is a saturation frequency (taken at the
   wall-band corner, ~10^5 Hz for Cu cylinder a = 5 mm). *)

omegaSat = 2 Pi 1*^5;
sSat     = I omegaSat;

gCoupling[s_, D_] := (mu / (2 Pi D^2)) * (s / (s + sSat));

(* Y_total(s, D) using Schur-F per-cylinder building blocks *)
YtotalSchur[s_, D_, N_] := Module[{Yc, g},
  Yc = YcylSchur[s, N];
  g = gCoupling[s, D];
  2 Yc / (1 + Yc g)];

(* Y_total(s, D) using Bessel exact per-cylinder building blocks *)
YtotalBessel[s_, D_] := Module[{Yc, g},
  Yc = YcylBessel[s];
  g = gCoupling[s, D];
  2 Yc / (1 + Yc g)];

(* === Sanity check 1: D -> infinity, parallel limit ======================= *)
Print["Sanity check 1: D -> infinity, parallel limit Y_total = 2 Y_cyl"];
Print[StringPadRight["D [mm]", 12] <>
      StringPadRight["Y_total/2Y_cyl at f=1e6 Hz", 32] <>
      "deviation from 1"];
Print[StringRepeat["-", 70]];
Do[
  Module[{D = Dmm 1*^-3, f = 1*^6, omega, s, ycyl, ytot, ratio},
    omega = 2 Pi f; s = I omega;
    ycyl = YcylBessel[s];
    ytot = YtotalBessel[s, D];
    ratio = ytot / (2 ycyl);
    Print[StringPadRight[ToString[N[Dmm, 4]], 12],
          StringPadRight[ToString[N[Abs[ratio], 6]], 32],
          N[Abs[ratio] - 1, 4]]
  ],
  {Dmm, {15, 50, 200, 2000, 200000}}];

Print[""];
Print["Sanity check 2: SIBC tail of 2-cylinder system at moderate D"];
Print["  expected: Y_total -> 2 K_SIBC^(cyl)/sqrt(s) = K_SIBC^(total)/sqrt(s)"];
Print["  with K_SIBC^(total) = 2 K_SIBC^(cyl) = ", N[2 KsibcCyl, 8]];
Print[""];
Print["f [Hz]  | |Y_total Bessel|  | ratio (Bes)  | |Y_total Schur-F|  | ratio (Schur)"];
Print[StringRepeat["-", 80]];
With[{Dprobe = 15*^-3, KsibcTotal = 2 KsibcCyl},
  Do[
    With[{omega = 2 Pi freq},
      With[{s = I omega},
        With[{yBes  = YtotalBessel[s, Dprobe],
              ySchur = YtotalSchur[s, Dprobe, basisN]},
          Print["f = ", N[freq, 3],
                "  |Y_B| = ", N[Abs[yBes], 5],
                "  rB = ",   N[Abs[yBes] Sqrt[omega] / KsibcTotal, 4],
                "  |Y_S| = ", N[Abs[ySchur], 5],
                "  rS = ",   N[Abs[ySchur] Sqrt[omega] / KsibcTotal, 4]]
        ]]],
    {freq, {1*^4, 1*^5, 1*^6, 1*^7, 1*^8, 1*^10, 1*^12}}]
];

Print[""];
Print["Sanity check 3: DC limit at moderate D"];
Print["  expected: Y_total(0) = 2 Y_cyl(0) (coupling vanishes at s=0)"];
Module[{Dprobe = 15*^-3, omega = 2 Pi * 0.01, s, yBes, ytot},
  s = I omega;
  yBes = YcylBessel[s];
  ytot = YtotalBessel[s, Dprobe];
  Print["  Y_cyl(s -> 0)        = ", N[Re[yBes], 6]];
  Print["  2 Y_cyl(s -> 0)      = ", N[Re[2 yBes], 6]];
  Print["  Y_total(s -> 0)      = ", N[Re[ytot], 6]];
  Print["  ratio                = ", N[Re[ytot]/Re[2 yBes], 6]]
];

(* === Frequency sweep CSV output ========================================= *)
Print[""];
Print["Writing frequency sweep CSV..."];
Dout = 15*^-3;
freqList = Table[10^k, {k, 2, 12, 0.25}];
resultsTable = {{"f_Hz", "Y_cyl_Bessel_mag", "Y_cyl_Schur_mag",
                 "Y_total_Bessel_mag", "Y_total_Schur_mag",
                 "K_SIBC_total_div_sqrt_omega"}};
KsibcTotal = 2 KsibcCyl;
Do[
  Module[{omega = 2 Pi f, s},
    s = I omega;
    AppendTo[resultsTable,
      {N[f, 6],
       N[Abs[YcylBessel[s]], 6],
       N[Abs[YcylSchur[s, basisN]], 6],
       N[Abs[YtotalBessel[s, Dout]], 6],
       N[Abs[YtotalSchur[s, Dout, basisN]], 6],
       N[KsibcTotal / Sqrt[omega], 6]}]
  ],
  {f, freqList}];

Export["bem_cln_2cylinder.csv", resultsTable, "CSV"];
Print["Wrote bem_cln_2cylinder.csv (", Length[freqList], " frequency points)"];

Print[""];
Print["=== bem_cln_2cylinder.wls DONE ==="];
