(* Stage 5 (gamma): Cauer-I extraction from Stage 5c Galerkin moments.

   Available moments (from Stage 5c PEC + dipole):
     alpha_1 = -mu_0 sigma (a^2+b^2)/48   (closed form)
     alpha_2 = +4.633960934739281e-10
     alpha_3 = -7.312818820589975e-15
     alpha_4 = +1.217153366239987e-19

   With 4 moments we can extract up to 2 R-L stages (4 unknowns).
   Extraction via continued-fraction expansion of Z(s) = -1/chi(s).
*)
$HistoryLength = 0;
prec = 60;

Print["=== Stage 5 (gamma): Cauer-I extraction ==="];
Print[""];

sigmaN = SetPrecision[58000000, prec];
mu0N = SetPrecision[4 Pi/10^7, prec];
aN = SetPrecision[5/1000, prec];
bN = SetPrecision[2/1000, prec];
cN = SetPrecision[1/1000, prec];
PiN = SetPrecision[Pi, prec];

(* Load Stage 5c moments *)
alpha1 = -mu0N sigmaN (aN^2 + bN^2)/48;
alpha2 = SetPrecision[4.633960934739281399894378775*^-10, prec];
alpha3 = SetPrecision[-7.3128188205899750960242786617291870572383071*^-15, prec];
alpha4 = SetPrecision[1.217153366239987350795349309870076568962*^-19, prec];

(* TODO: extend alpha_5..alpha_N from same Galerkin matrices *)
(* For now, use 4 moments. Full implementation re-runs Stage 5c augmented system
   with extended series expansion (no extra NIntegrate needed). *)

Print["chi(s) Taylor coefficients (Stage 5c PEC + dipole):"];
Print["  alpha_1 = ", N[alpha1, 16], "  (closed form)"];
Print["  alpha_2 = ", N[alpha2, 16]];
Print["  alpha_3 = ", N[alpha3, 16]];
Print["  alpha_4 = ", N[alpha4, 16]];
Print[""];

(* Note: alpha_5+ would extend the Pade table to more pairs.
   With 8 moments -> 4 pairs, 16 moments -> 8 pairs, etc.
   Hankel conditioning kappa ~ 4^N/N limits practical N to ~50 with 60-digit input.
   Basis truncation (8 interior modes + 1 dipole) limits accurate moments to ~10. *)


(* Z(s) = -1/chi(s) Laurent expansion at s=0 *)
(* chi(s) = alpha_1 s + alpha_2 s^2 + alpha_3 s^3 + alpha_4 s^4 + O(s^5)
            = alpha_1 s [1 + (alpha_2/alpha_1) s + (alpha_3/alpha_1) s^2 + ...]
   Z(s) = -1/chi = -1/(alpha_1 s) [1 + ... ]^-1
        = -1/(alpha_1 s) + (alpha_2/alpha_1^2)
          + (alpha_3/alpha_1^2 - alpha_2^2/alpha_1^3) s
          + ...
*)
chiSer = alpha1 sv + alpha2 sv^2 + alpha3 sv^3 + alpha4 sv^4;
ZLaurent = Series[-1/chiSer, {sv, 0, 3}];
Print["Z(s) = -1/chi(s) Laurent series (first 4 terms):"];
Print["  Z(s) = ", Normal[ZLaurent]];
Print[""];

(* Cauer-I form: Z(s) = 1/(s C_0) + R_0 + 1/(s C_1 + 1/(R_1 + 1/(s C_2 + ...)))
   But for eddy current: actual ladder has R-L pairs, not R-C.

   For CLN equivalent magnetic-circuit form:
     chi(s) corresponds to inductive Foster ladder
     Z_M(s) = 1/(s + 1/(R_0/L_0 + 1/(s + 1/(R_1/L_1 + ...))))
     ladder: R_0, L_0, R_1, L_1, ...

   Standard Cauer-I extraction: alternating L (in series) and R (in parallel).
   Assume Z(s) = -1/chi(s) starts with 1/s pole (1/(s C_0) form).
*)

cauerExtract[laurent_, svar_, maxSteps_] := Module[
  {Zcur, vals, k, leading, residual, leadCoeff},
  Zcur = laurent;
  vals = {};
  For[k = 0, k < maxSteps, k++,
    (* Get leading s^p term *)
    leadCoeff = SeriesCoefficient[Zcur, {svar, 0, -1}];
    If[Abs[N[leadCoeff, 30]] > 10^-30,
      (* leading 1/s term -> "L_k" element *)
      AppendTo[vals, {2 k, "1/(sC_" <> ToString[k] <> ")", 1/leadCoeff}];
      Zcur = Zcur - leadCoeff/svar,
      (* no 1/s term -- check constant *)
      leadCoeff = SeriesCoefficient[Zcur, {svar, 0, 0}];
      If[Abs[N[leadCoeff, 30]] > 10^-30,
        AppendTo[vals, {2 k, "R_" <> ToString[k], leadCoeff}];
        Zcur = Zcur - leadCoeff,
        Break[]]];
    (* Now invert remainder for next stage *)
    Zcur = Series[1/Zcur, {svar, 0, maxSteps - k}];
    leadCoeff = SeriesCoefficient[Zcur, {svar, 0, -1}];
    If[Abs[N[leadCoeff, 30]] > 10^-30,
      AppendTo[vals, {2 k + 1, "1/(sC_" <> ToString[k] <> "')", 1/leadCoeff}];
      Zcur = Zcur - leadCoeff/svar,
      leadCoeff = SeriesCoefficient[Zcur, {svar, 0, 0}];
      If[Abs[N[leadCoeff, 30]] > 10^-30,
        AppendTo[vals, {2 k + 1, "R_" <> ToString[k] <> "'", leadCoeff}];
        Zcur = Zcur - leadCoeff,
        Break[]]];
    If[Length[Zcur[[3]]] < 2, Break[]];
    Zcur = Series[1/Zcur, {svar, 0, maxSteps - k - 1}];
    ];
  vals];

stages = cauerExtract[ZLaurent, sv, 4];
Print["Cauer-I ladder elements:"];
Print["Index | Element              | Value (16 digits)"];
Print[StringJoin[Table["-", {70}]]];
Do[Module[{idx, type, val},
  {idx, type, val} = stages[[i]];
  Print[StringPadLeft[ToString[idx], 5], " | ",
    StringPadRight[type, 22], " | ",
    ToString[CForm[N[val, 16]]]]],
  {i, Length[stages]}];
Print[""];

(* Direct extraction for first 2 stages of standard R-L Cauer:
   Cauer-I: Z(s) = sL_0 + 1/(G_1 + 1/(sL_1 + ...))   with G = 1/R
   But since Z(s) starts with 1/s (capacitive-like leading), the "L_0" interpretation
   doesn't fit. Instead, use the Foster admittance form:

   Y(s) := -chi(s) (effective magnetic admittance)
   Y(s) = -alpha_1 s - alpha_2 s^2 - alpha_3 s^3 - ...
        = |alpha_1| s + (-alpha_2) s^2 + ... (alpha_1 < 0)
   For Foster Y(s) = sum a_n/(1+s tau_n) with a_n, tau_n > 0:
     Y(s) Taylor: a_n - a_n tau_n s + a_n tau_n^2 s^2 - ... summed
                = M_0 - M_1 s + M_2 s^2 - M_3 s^3 + ...
                where M_k = sum a_n tau_n^k

   But our Y starts at s (not constant) because chi(0) = 0. So
   Y(s)/s = |alpha_1| - alpha_2 s + ... corresponds to "Y/s = M_0 - M_1 s + ..." with M_0 = |alpha_1|.

   With M_0, M_1, M_2, M_3 known, can solve for 2 Foster pairs (a_1, tau_1, a_2, tau_2)
   via Pade. *)

Print["=== Foster pairs from chi(s)/s moments ==="];
(* Foster Y(s)/s = sum a_n/(1 + s tau_n) Taylor:
   = sum a_n - s sum a_n tau_n + s^2 sum a_n tau_n^2 - s^3 sum a_n tau_n^3 + ...
   With Y = -chi, Y/s = -alpha_1 - alpha_2 s - alpha_3 s^2 - alpha_4 s^3
   Match: sum a_n = -alpha_1 (= M_0)
          -sum a_n tau_n = -alpha_2 -> sum a_n tau_n = alpha_2 (= M_1)
          +sum a_n tau_n^2 = -alpha_3 -> sum a_n tau_n^2 = -alpha_3 (= M_2)
          -sum a_n tau_n^3 = -alpha_4 -> sum a_n tau_n^3 = alpha_4 (= M_3) *)
M0 = -alpha1; M1 = alpha2; M2 = -alpha3; M3 = alpha4;
Print["M_0 = ", N[M0, 16]];
Print["M_1 = ", N[M1, 16]];
Print["M_2 = ", N[M2, 16]];
Print["M_3 = ", N[M3, 16]];
Print[""];

(* Pade [1,1] approximant: (M_0 - M_0 tau s + ...) ~ a/(1 + s tau) form
   With moments [M_0, M_1, M_2, M_3]:
   Y/s = M_0 - M_1 s + M_2 s^2 - M_3 s^3 + ...
       ~ sum a_n/(1+s tau_n)

   For 2 Foster pairs: 4 unknowns, 4 moments -- solve.
   System:
     a_1 + a_2 = M_0
     a_1 tau_1 + a_2 tau_2 = M_1
     a_1 tau_1^2 + a_2 tau_2^2 = M_2
     a_1 tau_1^3 + a_2 tau_2^3 = M_3

   tau_1, tau_2 are roots of det([[M_1 M_0],[M_2 M_1]]-tau[[M_0 ?],...]) characteristic poly.

   Simpler approach: build Hankel matrix and use Prony's method.
   The poles tau_n are eigenvalues of T = inv(H1) H2 where H1 = [[M_0,M_1],[M_1,M_2]],
   H2 = [[M_1,M_2],[M_2,M_3]]. Then a_n via Vandermonde.
*)

H1 = {{M0, M1}, {M1, M2}};
H2 = {{M1, M2}, {M2, M3}};
tauMat = Inverse[H1].H2;
tauVals = Eigenvalues[tauMat];
Print["Foster time constants (Pade [2,2]):"];
Print["  tau_1, tau_2 = ", N[tauVals, 16]];
Print[""];

(* Extract residues a_n via Vandermonde *)
Vmat = {{1, 1}, {tauVals[[1]], tauVals[[2]]}};
aVals = LinearSolve[Vmat, {M0, M1}];
Print["Foster residues:"];
Print["  a_1, a_2 = ", N[aVals, 16]];
Print[""];

Print["Foster representation: chi(s) = -s sum a_n/(1+s tau_n) = "];
Print["  -s ", N[aVals[[1]], 8], "/(1+s ", N[tauVals[[1]], 8], ")"];
Print["  -s ", N[aVals[[2]], 8], "/(1+s ", N[tauVals[[2]], 8], ")"];
Print[""];

(* Convert Foster -> Cauer-I via standard recursion *)
(* For chi(s)/s = sum a_n/(1+s tau_n), the corresponding "admittance" Y(s) = chi(s)/s
   has 2-pole structure. Cauer-I: Y(s) = 1/(L_0 s + 1/(R_0 + 1/(L_1 s + 1/(R_1)))).
   But we need to rewrite in standard form.

   Direct interpretation:
   chi(s) = -|alpha_1| s + alpha_2 s^2 - ...
   "current" J ~ s chi(s) ... unclear without physical model.

   Skip detailed conversion -- just report Foster pairs as the meaningful output. *)

Print["=== Time constants compared to sphere reference ==="];
sphTau1 = mu0N sigmaN (5(aN^2 + bN^2)/24)/PiN^2;  (* tau_1 of sphere R_eff *)
sphTau2 = mu0N sigmaN (5(aN^2 + bN^2)/24)/(2 PiN)^2;
Print["Sphere R_eff first tau (tau_1 = mu0 sigma R_eff^2/pi^2): ",
  N[sphTau1 * 10^6, 8], " us"];
Print["Sphere R_eff second tau (tau_2 = mu0 sigma R_eff^2/(2pi)^2): ",
  N[sphTau2 * 10^6, 8], " us"];
Print["Galerkin extracted tau_1: ", N[tauVals[[1]] * 10^6, 8], " us"];
Print["Galerkin extracted tau_2: ", N[tauVals[[2]] * 10^6, 8], " us"];

Print[""];
Print["=== Stage 5 (gamma) summary ==="];
Print["With 4 Taylor moments (Stage 5c PEC+dipole), Pade [2,2] gives:"];
Print["  2 Foster pairs (tau_1, a_1) and (tau_2, a_2)"];
Print["  These approximate the true Foster series truncated to 2 dominant modes"];
Print["More accurate moments (Stage 5d, 5e) -> more Foster pairs / Cauer stages"];
