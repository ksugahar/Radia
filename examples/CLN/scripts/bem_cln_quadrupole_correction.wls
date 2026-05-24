#!/usr/bin/env wolframscript
(* ::Package:: *)

(*
  bem_cln_quadrupole_correction.wls

  Phase 4 final-2: quadrupole-order coupling correction to the
  Phase 2.5 rigorous BEM-CLN framework for 2D long cylinder pairs.

  Phase 2.5 includes only dipole-dipole coupling (1/D^2).  For
  finite ratio a/D the field from one cylinder is NOT uniform
  across the other cylinder's cross-section: there is a gradient
  ~1/D^3 that induces a *quadrupole* moment in the receiving
  cylinder.  The correction is O((a/D)^2) of the dipole result.

  For 5 mm cylinder pair at D = 15 mm: (a/D)^2 = (1/3)^2 = 0.111,
  so ~11 % correction expected.

  The 4 x 4 coupled system for 2 cylinders:
    (c_dip_1, c_quad_1, c_dip_2, c_quad_2) coupled by 2x2 sub-blocks
    [G_dd, G_dq; G_qd, G_qq]

  Coupling tensor entries (2D dipole/quadrupole pair, applied B in y,
  cylinders separated by D along x):
    G_dd(D)  = 1/D^2        (dipole field at distance D)
    G_dq(D)  = grad of dipole field / D^3 ~ 2/D^3 in magnitude
    G_qd(D)  = quadrupole field at distance D ~ 2/D^3
    G_qq(D)  = grad of quadrupole field ~ 6/D^4

  Per-cylinder polarizabilities:
    alpha_dipole(s) = a^2 [1 - 2 I_1(gamma a)/(gamma a I_0(gamma a))]
                                                       (Phase 2.5)
    alpha_quad(s)   ~ a^4 [1 - I_2(gamma a)/(... I_1(gamma a))]
                                                       (next Bessel order)

  For a simplified derivation we use scaling estimates:
    alpha_quad(s) / alpha_dip(s) = a^2 * scaling_factor
  where scaling_factor approaches 1 at PEC (perfect shielding of all
  multipoles) and 0 at DC.  A useful approximation is
    alpha_quad(s) ~= a^2 * alpha_dip(s)
  which captures the correct DC and PEC limits and the leading
  s-dependence.

  Run: wolframscript -file bem_cln_quadrupole_correction.wls
*)

ClearAll[a, sigma, mu];

a     = 5*^-3;
sigma = N[5.96*^7, 30];
mu    = N[4 Pi 1*^-7, 30];

Print["=== bem_cln_quadrupole_correction.wls ==="];
Print["Quadrupole correction to Phase 2.5 dipole-only coupling"];
Print["Cylinder a = ", a 1000, " mm, Cu"];
Print[""];

(* === Bessel exact polarizability (per Phase 2.5) ======================== *)
YcylBessel[s_?NumericQ] := Module[{gamma},
  gamma = Sqrt[s mu sigma];
  N[Pi a^2 sigma * 2 BesselI[1, gamma a] / (gamma a BesselI[0, gamma a]), 20]];

alphaDip[s_?NumericQ] := N[a^2 - YcylBessel[s] / (Pi sigma), 20];

(* === Quadrupole polarizability (limits-correct approximation) =========== *)
(* The exact Bessel-derived n=2 polarizability is complicated.  For an
   estimate that captures the correct DC and PEC limits and the leading
   s-dependence, we use the simplified scaling
        alpha_quad(s) = (a^2 / 2) * alpha_dip(s)
   which has:
     alpha_quad(0)        = 0                 (no eddy current at DC)
     alpha_quad(infinity) = a^4 / 2           (n=2 PEC shielding)
   The exact n=2 Bessel form has the same limits but different shape at
   intermediate frequencies; for the (a/D)^2 magnitude estimate this
   approximation is sufficient.  Refinement to exact n=2 Bessel via
   I_2(gamma a) matching is left as a small future improvement.        *)
alphaQuad[s_?NumericQ] := N[(a^2 / 2) alphaDip[s], 20];

(* Sanity: limits *)
Print["alpha_dip limits:"];
Print["  alpha_dip(0+)        = ", N[Re[alphaDip[I 2 Pi 0.01]], 6]];
Print["  alpha_dip(infinity)  = ", N[Re[alphaDip[I 2 Pi 1*^14]], 6],
      "  (target a^2 = ", N[a^2, 6], ")"];
Print[""];
Print["alpha_quad limits:"];
Print["  alpha_quad(0+)       = ", N[Re[alphaQuad[I 2 Pi 0.01]], 6]];
Print["  alpha_quad(infinity) = ", N[Re[alphaQuad[I 2 Pi 1*^14]], 6],
      "  (target a^4 = ", N[a^4, 6], ")"];
Print[""];

(* === 2D coupling tensor (dipole + quadrupole) =========================== *)
(* For cylinders separated by D along x, applied B in y *)
gDD[D_] := 1/D^2;           (* dipole-dipole *)
gDQ[D_] := 2/D^3;           (* dipole field gradient -> drives quadrupole *)
gQD[D_] := 2/D^3;           (* by reciprocity *)
gQQ[D_] := 6/D^4;           (* quadrupole-quadrupole *)

(* === 2-cylinder 4 x 4 system =========================================== *)
solve2cylQuad[s_?NumericQ, D_?NumericQ, B0_?NumericQ] := Module[
  {ad, aq, gdd, gdq, gqd, gqq, sysMat, rhs, sol},
  ad = alphaDip[s];
  aq = alphaQuad[s];
  gdd = gDD[D]; gdq = gDQ[D]; gqd = gQD[D]; gqq = gQQ[D];
  (* System: (c_d1, c_q1, c_d2, c_q2) coupled
     c_d1 - ad gdd c_d2 - ad gdq c_q2 = ad B0
     c_q1 - aq gqd c_d2 - aq gqq c_q2 = 0       (no quadrupole drive from uniform B)
     -ad gdd c_d1 - ad gdq c_q1 + c_d2 = ad B0
     -aq gqd c_d1 - aq gqq c_q1 + c_q2 = 0     *)
  sysMat = {
    {1,         0,         -ad gdd,   -ad gdq},
    {0,         1,         -aq gqd,   -aq gqq},
    {-ad gdd,  -ad gdq,    1,          0     },
    {-aq gqd,  -aq gqq,    0,          1     }
  };
  rhs = {ad B0, 0, ad B0, 0};
  sol = LinearSolve[sysMat, rhs];
  Association["c_d1" -> sol[[1]], "c_q1" -> sol[[2]],
              "c_d2" -> sol[[3]], "c_q2" -> sol[[4]]]];

(* === Phase 2.5 reference (dipole only) ================================ *)
solve2cylDipOnly[s_?NumericQ, D_?NumericQ, B0_?NumericQ] := Module[
  {ad, gdd},
  ad = alphaDip[s];
  gdd = gDD[D];
  ad B0 / (1 - ad gdd)];

(* === Test: 2 cylinders, a = 5 mm, D = 15 mm =========================== *)
Dpair = 15*^-3;
aOverD = a / Dpair;
Print["a/D = ", N[aOverD, 4]];
Print["(a/D)^2 = ", N[aOverD^2, 4], " (expected quadrupole correction magnitude)"];
Print[""];

Print["Comparison: dipole-only (Phase 2.5) vs dipole+quadrupole at f = 10^6 Hz"];
Print[""];
With[{sProbe = I 2 Pi 1*^6, B0 = 1.0},
  Module[{cDipOnly, sol4x4, cTotal4x4, correctionPct},
    cDipOnly = solve2cylDipOnly[sProbe, Dpair, B0];
    sol4x4 = solve2cylQuad[sProbe, Dpair, B0];
    cTotal4x4 = sol4x4["c_d1"];     (* dipole moment of cyl 1 *)
    correctionPct = 100 (Abs[cTotal4x4] - Abs[cDipOnly]) / Abs[cDipOnly];
    Print["  c_dip Phase 2.5 only           = ", N[cDipOnly, 5]];
    Print["  c_dip with quadrupole coupling = ", N[cTotal4x4, 5]];
    Print["  |c_quad| (induced quadrupole)  = ", N[Abs[sol4x4["c_q1"]], 5]];
    Print["  Correction (%)                 = ", N[correctionPct, 4]];
    Print["  Expected ~ (a/D)^2 = ", N[100 aOverD^2, 4], " %"]
  ]];

(* === Frequency sweep ==================================================== *)
Print[""];
Print["Quadrupole correction across wall band:"];
Print[""];
Print[StringPadRight["f [Hz]", 12] <>
      StringPadRight["|c_dip Phase2.5|", 22] <>
      StringPadRight["|c_dip with quad|", 22] <>
      "correction %"];
Print[StringRepeat["-", 70]];
Do[
  With[{omega = 2 Pi (10.0^logf), s = I 2 Pi (10.0^logf), B0 = 1.0},
    With[{cDip = solve2cylDipOnly[s, Dpair, B0],
          sol = solve2cylQuad[s, Dpair, B0]},
      Print["f=1e", logf, "       ",
            StringPadRight[ToString[N[Abs[cDip], 5]], 22],
            StringPadRight[ToString[N[Abs[sol["c_d1"]], 5]], 22],
            N[100 (Abs[sol["c_d1"]] - Abs[cDip])/Abs[cDip], 4]]
    ]],
  {logf, {3, 4, 5, 6, 7, 8}}];

Print[""];
Print["=== bem_cln_quadrupole_correction.wls DONE ==="];
