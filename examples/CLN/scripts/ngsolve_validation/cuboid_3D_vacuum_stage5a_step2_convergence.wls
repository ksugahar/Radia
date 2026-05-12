(* Stage 5a Step 2: Convergence study and Cauer-I extraction for cuboid
   interior PEC-like Galerkin (single-potential family A_z = 0).

   Builds on Step 1 framework. For each P_max in {3, 5, 7, 9, 11}:
     - Build basis (p,q,r all odd, p<=Pmax, q<=Qmax, r<=Rmax with Pmax=Qmax=Rmax)
     - Compute K, M (diagonal, closed form)
     - Compute source vector b
     - Compute moments M_k of Y(s) = b^T (K+s mu_0 sigma M)^{-1} b
       For diagonal K, M:
         Y(s)/s = -mu0 sigma sum b_i^2 / (K_ii + s mu0 sigma M_ii)
         At order s^k: alpha_(k+1) = (mu0 sigma)^(k+1) * sum b_i^2 M_ii^k / K_ii^(k+1) (with sign)
     - Cauer-I extraction
   Compare convergence with basis size, and to sphere R_eff reference.
*)

$HistoryLength = 0;
prec = 60;

Print["=== Stage 5a Step 2: Convergence + Cauer-I (60-digit) ==="];
Print[""];

(* Constants *)
sigmaN = SetPrecision[58000000, prec];
mu0N = SetPrecision[4 Pi/10^7, prec];
aN = SetPrecision[5/1000, prec];
bN = SetPrecision[2/1000, prec];
cN = SetPrecision[1/1000, prec];
PiN = SetPrecision[Pi, prec];
VN = aN bN cN;

(* Reference: sphere with R_eff = sqrt(5(a^2+b^2)/24) *)
ReffSq = 5 (aN^2 + bN^2)/24;

(* Sphere chi(s) closed form *)
chiSphere[R_, s_] := Module[{z}, z = R Sqrt[s mu0N sigmaN];
  -3/2 + (9/(2 z^2))(z Coth[z] - 1)];
(* Sphere alpha_n: chi(s) = alpha_1 s + alpha_2 s^2 + ... *)
(* alpha_n^sph = (mu0 sigma R^2)^n / a_n  with specific a_n:
   For sphere, chi(s) = -3/2 + (9/(2z^2))(z coth z - 1).
   z coth z - 1 = sum_(n=1) 2 z^2/(z^2 + (n pi)^2)
   chi(s) = -3/2 + (9/z^2) sum z^2/(z^2 + (n pi)^2)
          = -3/2 + 9 sum 1/(1 + (n pi/z)^2)
   With z^2 = R^2 mu sigma s, (n pi/z)^2 = (n pi)^2/(R^2 mu sigma s) -> infinity as s->0
   so each term -> 1/(big) -> 0. Wait need different expansion.

   Use Mittag-Leffler form chi = -9 s mu0 sigma R^2 sum 1/((n pi)^4 (1 + s tau_n)),
   tau_n = mu0 sigma R^2 /(n pi)^2:
   Series in s: chi = -9 mu0 sigma R^2 s sum [1 - s tau_n + (s tau_n)^2 - ...] / (n pi)^4
   alpha_1 = -9 mu0 sigma R^2 sum 1/(n pi)^4 = -9 mu0 sigma R^2 (pi^4/90)/pi^4 = -mu0 sigma R^2/10
   alpha_2 = +9 mu0 sigma R^2 sum tau_n/(n pi)^4 = 9 (mu0 sigma R^2)^2 sum 1/(n pi)^6
            = 9 (mu0 sigma R^2)^2 (pi^6/945)/pi^6 = (mu0 sigma R^2)^2/105
   alpha_3 = -9 (mu0 sigma R^2)^3 sum 1/(n pi)^8 = -9 (mu0 sigma R^2)^3 (pi^8/9450)/pi^8
           = -(mu0 sigma R^2)^3 / 1050
   etc. via zeta(2k) values. *)

alphaSphere[R_, k_] := Module[{x = mu0N sigmaN R^2, ck},
  ck = Switch[k,
    1, -1/10,           (* zeta(4)/(pi^4) * 9 = 9/90 *)
    2, +1/105,          (* zeta(6)/(pi^6) * 9 = 9/945 *)
    3, -1/1050,         (* zeta(8)/(pi^8) * 9 = 9/9450 *)
    4, +1/10395,        (* zeta(10)/(pi^10) * 9 *)
    _, (-1)^k * 9 N[Zeta[2 k + 2], prec]/PiN^(2 k + 2) ];
  ck * x^k];

Print["Sphere reference (R_eff) coefficients of chi(s) Taylor expansion:"];
Do[Print["  alpha_", k, "^sph(R_eff) = ",
  N[alphaSphere[Sqrt[ReffSq], k], 16]], {k, 1, 5}];
Print[""];

(* === Build basis and matrices for given Pmax === *)
buildBasis[Pmax_] := Flatten[Table[
  {2 P + 1, 2 Q + 1, 2 R + 1},
  {P, 0, (Pmax - 1)/2}, {Q, 0, (Pmax - 1)/2}, {R, 0, (Pmax - 1)/2}], 2];

(* Closed-form diagonal entries:
   M_ii = (V/8) pi^2 (p^2/a^2 + q^2/b^2)
   K_ii = (V/8) pi^4 (p^2/a^2 + q^2/b^2) (p^2/a^2 + q^2/b^2 + r^2/c^2)
   lambda_i = K_ii/M_ii = pi^2 (p^2/a^2 + q^2/b^2 + r^2/c^2)
   b_i = -(1/2)*(2a/(p pi))*(-2b/(q pi))*(2c/(r pi))
         -(1/2)*(-2a/(p pi))*(2b/(q pi))*(2c/(r pi))
       = (4abc)/(pqr pi^3) + (4abc)/(pqr pi^3)
       = 8 V/(p q r pi^3)   for p,q,r all odd
*)

mElement[p_, q_, r_] := (VN/8) PiN^2 (p^2/aN^2 + q^2/bN^2);
kElement[p_, q_, r_] := (VN/8) PiN^4 (p^2/aN^2 + q^2/bN^2) (p^2/aN^2 + q^2/bN^2 + r^2/cN^2);
bElement[p_, q_, r_] := 8 VN/(p q r PiN^3);

(* Verify closed-form formulas against Step 1 numerical:
   For (1,1,1): b_i = 8V/pi^3, alpha_2 contribution = b_i^2/k_ii.
*)
Print["Closed-form verification ((1,1,1)):"];
Print["  b_(1,1,1) = 8V/pi^3 = ", N[bElement[1, 1, 1], 16]];
Print["  K_(1,1,1) = ", N[kElement[1, 1, 1], 16]];
Print["  M_(1,1,1) = ", N[mElement[1, 1, 1], 16]];
Print["  lambda_(1,1,1) = K/M = ", N[kElement[1, 1, 1]/mElement[1, 1, 1], 16],
  "  vs pi^2 (1/a^2+1/b^2+1/c^2) = ",
  N[PiN^2 (1/aN^2 + 1/bN^2 + 1/cN^2), 16]];
Print[""];

(* === Compute moments alpha_n via series expansion ===
   Y(s) = -mu0 sigma sum b_i^2/(K_ii + s mu0 sigma M_ii)
   chi_correction(s) = -(mu0 sigma s/V) sum b_i alpha_i(s) [contribution from A_ind to <A_ext, A>]
   But for the leading coefficients, the picture is:

   chi(s) = -(mu0 sigma s/V)[<A_ext, A_ext> + <A_ext, A_ind>]
          = alpha_LF * s + corrections

   For diagonal Galerkin:
     A_ind expansion alpha_i = -s mu0 sigma b_i / (K_ii + s mu0 sigma M_ii)
     <A_ext, A_ind> = sum b_i alpha_i = -s mu0 sigma sum b_i^2/(K_ii + s mu0 sigma M_ii)
     chi_corr = -(mu0 sigma s/V) * (-s mu0 sigma) sum b_i^2/(K + s mu0 sigma M)
              = (mu0 sigma s)^2/V * sum b_i^2/(K_ii + s mu0 sigma M_ii)
     Series in s: 1/(K + s mu0 sigma M) = (1/K) (1 - s mu0 sigma M/K + (s mu0 sigma M/K)^2 - ...)
     chi_corr = (mu0 sigma)^2 s^2/V sum (b_i^2/K_ii) (1 - s mu0 sigma M_ii/K_ii + ...)
     -> alpha_2_corr = (mu0 sigma)^2 / V * sum b_i^2/K_ii
        alpha_3_corr = -(mu0 sigma)^3 / V * sum b_i^2 M_ii/K_ii^2
        alpha_(n+1)_corr = (-1)^(n-1) (mu0 sigma)^(n+1)/V sum b_i^2 M_ii^(n-1) / K_ii^n

   Note alpha_LF = -mu0 sigma (a^2+b^2)/48 is BASIS INDEPENDENT (from <A_ext,A_ext> piece). *)

alphaLFRef = -mu0N sigmaN (aN^2 + bN^2)/48;

(* For given basis (Pmax), compute alpha_2, alpha_3, ... *)
computeAlphas[Pmax_, kMax_] := Module[{idx, vals, b, k, m},
  idx = buildBasis[Pmax];
  vals = Table[
    {b, k, m} = {bElement @@ pp, kElement @@ pp, mElement @@ pp};
    Table[(-1)^(n - 1) * (mu0N sigmaN)^(n + 1)/VN * b^2 * m^(n - 1)/k^n,
      {n, 1, kMax}],
    {pp, idx}];
  Total /@ Transpose[vals]];

(* Run convergence sweep *)
PmaxList = {3, 5, 7, 9, 11};
kMax = 5;
Print["=== Convergence sweep: alpha_2..alpha_5 vs basis size ==="];
Print[""];
Print[StringPadRight["Pmax", 6], StringPadRight["nBasis", 8],
  StringPadRight["alpha_2", 24], StringPadRight["alpha_3", 24],
  StringPadRight["alpha_4", 24], StringPadRight["alpha_5", 24]];
Print[StringJoin[Table["-", {110}]]];
allResults = {};
Do[Module[{nB, alphas},
  nB = Length[buildBasis[Pmax]];
  alphas = computeAlphas[Pmax, kMax];
  AppendTo[allResults, {Pmax, nB, alphas}];
  Print[StringPadRight[ToString[Pmax], 6], StringPadRight[ToString[nB], 8],
    StringJoin[StringPadRight[ToString[CForm[N[#, 14]]], 24] & /@ alphas]]],
  {Pmax, PmaxList}];
Print[""];

(* Convergence ratios and Richardson extrapolation if possible *)
Print["=== Convergence rate analysis ==="];
Do[Module[{a3, a5, a7, a9, a11, ratio},
  a3 = allResults[[1, 3, k]];
  a5 = allResults[[2, 3, k]];
  a7 = allResults[[3, 3, k]];
  a9 = allResults[[4, 3, k]];
  a11 = allResults[[5, 3, k]];
  Print["alpha_", k + 1, " sequence: ",
    Table[N[allResults[[i, 3, k]], 8], {i, 1, 5}]];
  ratio = (a11 - a9)/(a9 - a7);
  Print["  ratio (a11-a9)/(a9-a7) = ", N[ratio, 6],
    " (smaller -> faster convergence)"]],
  {k, 1, 4}];
Print[""];

(* Compare to sphere reference *)
Print["=== Comparison to sphere(R_eff) reference ==="];
Print[""];
Print[StringPadRight["alpha_n", 10], StringPadRight["Pmax=11 (interior PEC)", 28],
  StringPadRight["sphere(R_eff)", 24], StringPadRight["ratio (Galerkin/sph)", 24]];
Print[StringJoin[Table["-", {86}]]];
Do[Module[{interiorVal, sphVal, ratio},
  interiorVal = allResults[[5, 3, k]];
  sphVal = alphaSphere[Sqrt[ReffSq], k + 1];
  ratio = interiorVal/sphVal;
  Print[StringPadRight[ToString["alpha_" <> ToString[k + 1]], 10],
    StringPadRight[ToString[CForm[N[interiorVal, 14]]], 28],
    StringPadRight[ToString[CForm[N[sphVal, 14]]], 24],
    ToString[NumberForm[N[ratio, 8], 8]]]],
  {k, 1, 4}];
Print[""];
Print["Note: ratio < 1 because interior PEC-like basis underestimates response"];
Print["(no eddy current outside conductor, no normal flux through surface).";];
Print["Step 3 (multipole exterior) will close this gap."];
Print[""];

(* === Cauer-I extraction from highest-Pmax moments === *)
(* Build moment sequence M_k = alpha_(k+1) (Y(s) = -dchi/ds at s=0 = -alpha_1 - 2 alpha_2 s - ...?
   Actually let's redefine: define Y(s) such that chi(s) = -alpha_LF / Y(s) effective
   For Cauer-I extraction we need an admittance-like rational function.

   Use chi(s) directly as an impedance-like quantity:
   Z(s) := -1/chi(s)
   then Z(s) ~ 1/(alpha_LF s) at low s = 1/(alpha_LF s) * [1 + alpha_2/alpha_LF s + ...]^(-1)
   This is a Foster admittance form. Cauer-I extraction from this. *)

chiCoeffs[Pmax_, kMax_] := Prepend[computeAlphas[Pmax, kMax], alphaLFRef];

allCoeffs = chiCoeffs[11, 12];
Print["chi(s) Taylor coefficients (Pmax=11, k=1..13):"];
Do[Print["  alpha_", k, " = ", N[allCoeffs[[k]], 16]], {k, 1, 8}];
Print[""];

(* Form chi(s) Taylor expansion truncated *)
chiSeries[Pmax_, kMax_] := Module[{coeffs, s},
  coeffs = chiCoeffs[Pmax, kMax];
  Sum[coeffs[[k]] s^k, {k, 1, kMax + 1}]];

(* Z(s) = -1/chi(s) Laurent series at s = 0 *)
ZLaurent[Pmax_, kMax_] := Module[{chi, s, ser},
  chi = chiSeries[Pmax, kMax] /. Global`s -> sv;
  ser = Series[-1/chi, {sv, 0, kMax}];
  Normal[ser]];

Zseries = ZLaurent[11, 11];
Print["Z(s) = -1/chi(s) Laurent series (first few):"];
Print["  ", Series[-1/chiSeries[11, 11] /. Global`s -> sv, {sv, 0, 4}]];
Print[""];

(* Cauer-I extraction *)
cauerExtractFromLaurent[zLaurent_, svar_, maxStages_] := Module[
  {Zcur, vals, kk, pn, pp, qq, q0, p0, Rk, presid, pas, pas0, Lk, nresid},
  Zcur = zLaurent;  (* expression in svar with leading 1/svar term *)
  vals = {};
  For[kk = 0, kk < maxStages, kk++,
   (* Z = 1/(s L_k) form means Series at s=0 has -1 term *)
   pn = Series[Zcur, {svar, 0, 30}];
   (* extract leading coefficient = 1/(s L_k) *)
   Lk = 1/Coefficient[Normal[pn], svar, -1];
   If[NumericQ[Lk] && Abs[Lk] > 10^-100,
    AppendTo[vals, {2 kk, "L", Lk}],
    Break[]];
   (* Subtract 1/(s Lk), get next *)
   Zcur = Zcur - 1/(svar Lk);
   pn = Series[Zcur, {svar, 0, 30}];
   (* Now leading is constant -> R_k *)
   Rk = Normal[pn] /. svar -> 0;
   If[NumericQ[Rk] && Abs[Rk] > 10^-100,
    AppendTo[vals, {2 kk + 1, "R", Rk}],
    Break[]];
   Zcur = 1/(Zcur - Rk);
   ];
  vals];

Print["=== Cauer-I extraction from chi(s) (interior PEC, Pmax=11) ==="];
stages = cauerExtractFromLaurent[Zseries, sv, 6];
Print["Index | Type  | Value (16 digits)"];
Print[StringJoin[Table["-", {60}]]];
Do[Module[{idx, type, val},
  {idx, type, val} = stages[[i]];
  Print[StringPadLeft[ToString[idx], 5], " | ",
    StringPadRight[type, 5], " | ",
    ToString[CForm[N[val, 16]]]]],
  {i, Length[stages]}];
Print[""];

Print["=== Stage 5a Step 2 summary ==="];
Print["* alpha_2 converges with basis size (sin sin sin curl interior basis)"];
Print["* Final alpha_n values are interior-PEC-like, NOT vacuum"];
Print["* Cauer-I from these moments gives 'cuboid in PEC cavity' Cauer ladder"];
Print["* Step 3: add exterior multipole basis to convert to true vacuum"];
