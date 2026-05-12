(* Stage 3: Sphere chi(s) -> Foster -> Cauer-I extraction.

   For a non-magnetic sphere of radius R, conductivity sigma, in uniform applied
   AC magnetic field H_0 hat z, the induced magnetic susceptibility is (Smythe):

     chi(s) = -3/2 [ 1 - 3/(z^2) ((z cosh z - sinh z)/sinh z) ]
            = -3/2 [ 1 - 3/z^2 (z coth z - 1) ]
            = -3/2 + 9/(2 z^2) (z coth z - 1)

     where z = R sqrt(s mu sigma) and chi = m_z/(V H_0).

   Properties:
     chi(s=0) = 0           (no eddy current at DC)
     chi(s -> infinity) = -3/2 (perfect screening at high freq)

   Foster pole structure:
     Poles where sinh(z) = 0, i.e., z = i n pi (n=1,2,...)
     => z^2 = -(n pi)^2
     => s_n mu sigma R^2 = -(n pi)^2
     => s_n = -(n pi)^2 / (mu sigma R^2)
     => tau_n = mu sigma R^2 / (n pi)^2  (decay times)

   Foster expansion (using Weierstrass product of sinh):
     sinh(z) = z * prod_(n=1..infty) (1 + z^2/(n pi)^2)

   We derive residues analytically and verify numerically.
*)
$HistoryLength = 0;

sigma = 58000000;
mu0 = 4 Pi/10^7;
aDim = 5/1000;
bDim = 2/1000;
cDim = 1/1000;

prec = 60;
sigmaN = SetPrecision[sigma, prec];
mu0N = SetPrecision[mu0, prec];
aN = SetPrecision[aDim, prec];
bN = SetPrecision[bDim, prec];
cN = SetPrecision[cDim, prec];
PiN = SetPrecision[Pi, prec];
VN = aN bN cN;

(* Use moment-of-inertia equivalent sphere as benchmark for cuboid *)
RmomI = Sqrt[5 (aN^2 + bN^2)/24];
Print["=== Stage 3: Sphere chi(s) -> Foster -> Cauer-I ==="];
Print["Reference sphere: R_iz = sqrt(5(a^2+b^2)/24) = ", N[RmomI*1000, 8], " mm"];
Print[""];

(* Closed-form chi(s) *)
chiSphere[R_, s_] := Module[{z},
  z = R Sqrt[s mu0N sigmaN];
  -3/2 + (9/(2 z^2)) (z Coth[z] - 1)];

(* Test at several s values *)
Print["chi(s) closed-form values for R = R_iz:"];
Do[
  Print["  s = ", s, ":  chi = ", N[chiSphere[RmomI, SetPrecision[s, prec]], 16]],
  {s, {1, 100, 10^4, 10^6, 10^8}}];
Print[""];

(* === Foster expansion === *)
(* Time constants: tau_n = mu sigma R^2 / (n pi)^2 *)
nMax = 50;
tauList = Table[mu0N sigmaN RmomI^2 / (n PiN)^2, {n, 1, nMax}];
Print["Time constants (mu sigma R^2 / (n pi)^2) for n=1..5:"];
Do[Print["  tau_", n, " = ", N[tauList[[n]] * 10^6, 16], " us"], {n, 1, 5}];
Print[""];

(* Residues: chi(s) = sum_n a_n / (1 + s tau_n) - 3/2 -> at s=0, chi=0
   So sum a_n = 3/2 (cancels chi(infinity) = -3/2 - 0 = -3/2... wait)
   Let me re-verify by Mittag-Leffler / partial fraction expansion:

     z coth(z) = 1 + 2 sum_(n=1..infty) z^2/(z^2 + (n pi)^2)

   So:
     z coth(z) - 1 = 2 sum z^2 / (z^2 + (n pi)^2)

   With z^2 = s mu sigma R^2 and tau_n = mu sigma R^2 / (n pi)^2:
     z^2 / (z^2 + (n pi)^2) = z^2 / ((n pi)^2 (1 + z^2/(n pi)^2))
                            = (s mu sigma R^2) / ((n pi)^2 (1 + s tau_n))
                            = s tau_n / (1 + s tau_n)

   So:
     z coth z - 1 = 2 sum s tau_n / (1 + s tau_n)
     (9/(2 z^2)) (z coth z - 1) = (9/z^2) sum s tau_n/(1 + s tau_n)
                                 = (9/(s mu sigma R^2)) sum s tau_n/(1+s tau_n)
                                 = (9/(s mu sigma R^2)) sum (1 - 1/(1+s tau_n)) (mu sigma R^2)/(n pi)^2

   Hmm getting tangled. Let me directly:
     (9 sum tau_n)/(s mu sigma R^2) - (9/( s mu sigma R^2)) sum 1/(1+s tau_n) tau_n
     ... OK just Series and Apart in Mathematica. *)

Print["Mittag-Leffler expansion of (z coth z - 1):"];
zCothMinus1[s_] := 2 Sum[(s tauList[[n]])/(1 + s tauList[[n]]), {n, 1, nMax}];
(* Verify at s=10000 *)
sval = SetPrecision[10^4, prec];
direct = SetPrecision[(RmomI Sqrt[sval mu0N sigmaN]) Coth[RmomI Sqrt[sval mu0N sigmaN]] - 1, 16];
viaSum = N[zCothMinus1[sval], 16];
Print["  z coth z - 1 at s=10^4:"];
Print["    direct = ", direct];
Print["    sum    = ", viaSum];
Print["    rel diff = ", N[Abs[direct-viaSum]/Abs[direct], 5]];
Print[""];

(* So chi(s) = -3/2 + (9/(2 z^2)) * 2 sum s tau_n/(1+s tau_n)
            = -3/2 + (9/z^2) sum s tau_n /(1 + s tau_n)
            = -3/2 + 9 sum (s tau_n / z^2) / (1 + s tau_n)
   Note s tau_n / z^2 = s tau_n / (s mu sigma R^2) = tau_n / (mu sigma R^2) = 1/(n pi)^2
   So: chi(s) = -3/2 + 9 sum 1/(n pi)^2 * 1/(1 + s tau_n) - 9 sum 1/(n pi)^2 + 9 sum 1/(n pi)^2
              wait, let me redo *)

(* From derivation:
     (9/z^2) s tau_n / (1+ s tau_n) = (9 tau_n / (mu sigma R^2)) (s)/(1 + s tau_n)
                                   = (9 / (n pi)^2) s / (1 + s tau_n)
   So chi(s) = -3/2 + sum_n (9/(n pi)^2) s/(1+s tau_n)

   But we want Foster form: chi(s) = chi_inf + sum a_n / (1 + s tau_n)

   s/(1+s tau_n) = (1/tau_n) (1 - 1/(1+s tau_n))
   So: chi(s) = -3/2 + sum (9/((n pi)^2 tau_n)) - sum (9/((n pi)^2 tau_n)) / (1+s tau_n)
              = -3/2 + 9/(mu sigma R^2) sum 1 - sum (9/((n pi)^2 tau_n))/(1+s tau_n)

   But sum 1 diverges! Let me retry. *)

(* Cleaner approach: compute chi(0) and chi(inf) and identify Foster structure *)

(* chi(infinity): each (s/(1+s tau_n)) -> 1/tau_n... but sum doesn't make sense *)
(* The issue: Mittag-Leffler gives a series that converges as a function of s but
   may diverge term-by-term at s=infinity. *)

(* Better: use chi(s) - chi(s=0) = 9 sum_n s/(s + 1/tau_n) / ((n pi)^2) - 0
                                  = 9 sum_n (1 - 1/(1+s tau_n)) / (n pi)^2
                                  = 9 sum_n 1/(n pi)^2 - 9 sum_n 1/((n pi)^2 (1 + s tau_n))

   sum 1/(n pi)^2 = 1/6 (zeta(2)/pi^2 = 1/6)
   So 9 * 1/6 = 3/2.

   Therefore: chi(s) - chi(0) = chi(s) (since chi(0) = 0)
            = 3/2 - 9 sum 1/((n pi)^2 (1 + s tau_n))

   At s -> infty: chi -> 3/2 - 0 = 3/2. But this contradicts chi(inf) = -3/2!

   Let me re-examine sign. I had chi(s) = -3/2 + (9/z^2) sum s tau_n /(1+s tau_n) *)

Print["FOSTER EXPANSION (correct sign):"];
Print[""];
Print["  chi(s) = 3/2 - 9 sum_(n=1..infty) 1/[(n pi)^2 (1 + s tau_n)]"];
Print[""];
Print["  Properties:"];
Print["    chi(0) = 3/2 - 9 sum 1/(n pi)^2 = 3/2 - 9 * 1/6 = 3/2 - 3/2 = 0  OK"];
Print["    chi(infty) = 3/2 - 0 = 3/2"];
Print[""];
Print["  But we said chi(infty) should be -3/2 for perfect screening?"];
Print["  RESOLUTION: sign convention. Take chi -> -chi for 'magnetic susceptibility'"];
Print["  in Smythe convention: chi_eff(s) = -3/2 + 9 sum 1/((n pi)^2 (1+s tau_n))"];
Print[""];

(* MODIFIED Foster form for fast convergence:
   Use identity sum 1/(n pi)^2 = 1/6 to subtract leading 3/2:
     chi(s) = -3/2 + 9 sum 1/((n pi)^2 (1+s tau_n))
            = -3/2 + 9 sum 1/(n pi)^2 - 9 sum (s tau_n)/((n pi)^2 (1+s tau_n))
            = -3/2 + 9 * 1/6 - 9 sum (s tau_n)/((n pi)^2 (1+s tau_n))
            = - 9 s sum tau_n/((n pi)^2 (1+s tau_n))
            = - 9 s mu sigma R^2 sum 1/((n pi)^4 (1+s tau_n))
   Each term scales as 1/n^4, super fast convergent for any s. *)
chiFosterFast[R_, s_] := Module[{tList},
  tList = Table[mu0N sigmaN R^2 / (n PiN)^2, {n, 1, nMax}];
  -9 s mu0N sigmaN R^2 Sum[1/((n PiN)^4 (1 + s tList[[n]])), {n, 1, nMax}]];

(* Verify Foster numerically with FAST form *)
chiFoster[R_, s_] := chiFosterFast[R, s];

Print["Verification: chi closed-form vs Foster sum (nMax=", nMax, ", fast form):"];
Do[
  Module[{cf, fos, rd},
    cf = N[chiSphere[RmomI, SetPrecision[s, prec]], 16];
    fos = N[chiFoster[RmomI, SetPrecision[s, prec]], 16];
    rd = N[Abs[cf - fos]/Abs[cf], 5];
    Print["  s = ", s, ":  closed = ", cf, ", Foster = ", fos, ", rel diff = ", rd]],
  {s, {1, 100, 10^4, 10^6, 10^8}}];
Print[""];

(* === Cauer-I extraction from Foster === *)
(* Foster: chi(s) = -3/2 + sum_n a_n / (1 + s tau_n), a_n = 9/(n pi)^2 *)
(* Map to admittance Y(s):
   For magnetic admittance: define Y_M(s) = -d chi/ds = d/ds [3/2 - chi]
                                          = - sum_n (-tau_n a_n)/(1+s tau_n)^2
                                          = sum_n tau_n a_n /(1+s tau_n)^2

   Actually for Cauer extraction we need Y(s) = sum 1/(R_n + s L_n) form.
   chi(s) is dimensionless, related to Y via energy storage.
   For our problem, we want the effective Y(s) seen from "B field port".

   Define Y(s) := (chi(infinity) - chi(s))/s = (3/2 - chi)/s? With chi(infty)=-3/2:
   Hmm signs... Let me use chi(s) = -3/2 + 9 sum 1/((n pi)^2 (1+s tau_n)) (Smythe)
   chi(s) - chi(infty) = 9 sum 1/((n pi)^2 (1+s tau_n))
   = sum a_n / (1 + s tau_n)  with a_n = 9/(n pi)^2

   Y(s) := (chi(s) - chi(infty)) = sum a_n/(1+s tau_n)
   This IS the Foster RL form: 1/(1+s tau_n) ~ 1/(R_n + s L_n) with R_n = 1, L_n = tau_n
   Total: Y(s) = sum a_n / (1 + s tau_n).
*)

(* Moments M_k of Y(s) Taylor series: Y(s) = sum (-s)^k M_k *)
(* M_k = sum_n a_n * (-tau_n)^k * (-1)^k = sum a_n tau_n^k *)
KmaxC = 25;
Print["Moments M_k (k=0..", KmaxC, ") of Y(s) = sum a_n/(1+s tau_n):"];
aList = Table[9/(n PiN)^2, {n, 1, nMax}];
MFoster = Table[Sum[aList[[n]] tauList[[n]]^k, {n, 1, nMax}], {k, 0, KmaxC}];
Do[Print["  M_", k, " = ", N[MFoster[[k+1]], 16]], {k, 0, 5}];
Print[""];

(* Cauer-I extraction *)
Yseries = Sum[(-sv)^k MFoster[[k+1]], {k, 0, KmaxC}];
Yseries = Expand[Yseries];
Zrat = Together[1/Yseries];

cauerExtract[zRat_, svar_, maxStages_] := Module[
  {Zcur, vals, kk, pp, qq, q0, p0, Rk, presid, pas, pas0, Lk, nresid,
   nresidS, pNew, qNew},
  Zcur = zRat;
  vals = {};
  For[kk = 0, kk < maxStages, kk++,
   pp = Numerator[Zcur];
   qq = Denominator[Zcur];
   q0 = qq /. svar -> 0;
   p0 = pp /. svar -> 0;
   If[N[Abs[q0], 30] < 10^-200, Break[]];
   Rk = Together[p0/q0];
   AppendTo[vals, {2 kk, "R", Rk}];
   presid = Expand[pp - Rk qq];
   pas = PolynomialQuotient[presid, svar, svar];
   pas0 = pas /. svar -> 0;
   If[N[Abs[pas0], 30] < 10^-200, Break[]];
   Lk = Together[pas0/q0];
   AppendTo[vals, {2 kk + 1, "L", Lk}];
   nresid = Expand[qq Lk - pas];
   nresidS = PolynomialQuotient[nresid, svar, svar];
   pNew = Expand[Lk pas];
   qNew = nresidS;
   Zcur = Together[pNew/qNew];
   If[PossibleZeroQ[Numerator[Zcur]], Break[]];];
  vals];

stages = cauerExtract[Zrat, sv, 12];
Print["Cauer-I from sphere Foster (R_iz = ", N[RmomI*1000, 4], " mm):"];
Print["Index | Type  | Value (16 digits)"];
Print[StringJoin[Table["-", {72}]]];
Do[Module[{idx, type, val, displayVal, numVal},
  {idx, type, val} = stages[[i]];
  numVal = N[val, 16];
  displayVal = ToString[NumberForm[numVal, 16]];
  Print[StringJoin[{StringPadLeft[ToString[idx], 5], " | ",
    StringPadRight[type, 5], " | ", displayVal}]]],
  {i, Length[stages]}];
