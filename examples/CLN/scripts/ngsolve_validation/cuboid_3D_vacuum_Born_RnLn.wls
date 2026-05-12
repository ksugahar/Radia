(* Final: extract (R_n, L_n) Cauer ladder from Born-derived alpha_n.
   Combines closed-form alpha_1, alpha_2 with numerical alpha_3+ from various sources. *)
$HistoryLength = 0;
prec = 60;

Print["=== Born CLN extraction: R_n, L_n ==="];
Print[""];

sigCu = SetPrecision[58000000, prec];
mu0 = SetPrecision[4 Pi/10^7, prec];
aN = SetPrecision[5/1000, prec];
bN = SetPrecision[2/1000, prec];
cN = SetPrecision[1/1000, prec];
piN = SetPrecision[Pi, prec];
vol = aN bN cN;

(* Born-derived alpha_n *)
(* alpha_1: closed form *)
alpha1Born = -mu0 sigCu (aN^2 + bN^2)/48;

(* alpha_2: closed form via Born series, computed earlier *)
alpha2Born = SetPrecision[7.4783240535434074691212813798857319216568135873611*^-10, prec];

(* alpha_3: from Born T_xx + T_yy 3D integrals (MC 3-digit precision) *)
alpha3BornMC = SetPrecision[-1.435266483988741*^-14, prec];

(* Stage 5c Galerkin alpha_n (for cross-check) *)
alpha2Galerkin = SetPrecision[4.633960934739281399894378775*^-10, prec];
alpha3Galerkin = SetPrecision[-7.3128188205899750960242786617*^-15, prec];

Print["Born series moments:"];
Print["  alpha_1 (closed)       = ", N[alpha1Born, 16]];
Print["  alpha_2 (closed)       = ", N[alpha2Born, 16]];
Print["  alpha_3 (MC ~3 digits) = ", N[alpha3BornMC, 8], "  <-- update with actual"];
Print[""];

Print["Stage 5c Galerkin (closed PEC + 1 dipole, comparison):"];
Print["  alpha_2_Galerkin       = ", N[alpha2Galerkin, 16],
  "  ratio Born/Gal = ", N[alpha2Born/alpha2Galerkin, 6]];
Print["  alpha_3_Galerkin       = ", N[alpha3Galerkin, 16]];
Print[""];

(* === Foster moments and Cauer-I extraction === *)
(* Y(s) = -chi(s), Y/s = sum a_n/(1+s tau_n)
   Y/s Taylor: M_0 - M_1 s + M_2 s^2 - M_3 s^3 + ...
   Identification:
     M_0 = -alpha_1 = +mu_0 sigma (a^2+b^2)/48
     M_1 = +alpha_2 = +(mu_0 sigma)^2 ...
     M_2 = -alpha_3 = ... (positive)
     M_k = (-1)^(k+1) alpha_(k+1)
   All M_k > 0 for valid Foster. *)

M0 = -alpha1Born;
M1 = alpha2Born;
M2 = -alpha3BornMC;

Print["Foster moments M_k = (-1)^(k+1) alpha_(k+1):"];
Print["  M_0 = ", N[M0, 16], If[N[M0, 30] > 0, " (+ ✓)", " (-) ✗"]];
Print["  M_1 = ", N[M1, 16], If[N[M1, 30] > 0, " (+ ✓)", " (-) ✗"]];
Print["  M_2 = ", N[M2, 16], If[N[M2, 30] > 0, " (+ ✓)", " (-) ✗"]];
Print[""];

(* Pade [1,1] gives 1 Foster pair from 3 moments: M_0, M_1, M_2 *)
(* a/(1+s tau) Taylor: a - a tau s + a tau^2 s^2 ...
   M_0 = a, M_1 = a tau, M_2 = a tau^2
   tau = M_1/M_0
   a = M_0
   Verify: M_2 should = M_1^2/M_0 if single pole *)

tauPade1 = M1/M0;
aPade1 = M0;

Print["Pade [1,1] (1 Foster pair):"];
Print["  tau_1 = M_1/M_0 = ", N[tauPade1 * 10^6, 12], " us"];
Print["  a_1 = M_0      = ", N[aPade1, 16]];
Print["  Self-check: M_1^2/M_0 = ", N[M1^2/M0, 12], " vs M_2 = ", N[M2, 12]];
Print[""];

(* Pade [2,2] from 4 moments would be ideal. We have M_0, M_1, M_2 (3 moments).
   So at most 1 Foster pair properly determined. *)

(* Single R-L stage from single Foster pair:
   chi_1pair(s) = -a_1/(1+s tau_1)
   For circuit equivalent:
     R_0 = 1/a_1 * V (some unit normalization), L_0 = R_0 * tau_1 *)

(* Actually for CLN ladder R_n, L_n the conversion from Foster requires
   continued fraction extraction. With 1 pair, just get 1 stage:
     R_0, L_0 with L_0/R_0 = tau_1, magnitude from a_1. *)

Print["Single-stage equivalent CLN:"];
Print["  tau_1 = ", N[tauPade1 * 10^6, 8], " us"];
Print["  R/L ratio = 1/tau = ", N[1/tauPade1, 8], " s^-1"];
Print[""];

(* For CLN time constant: tau_n = mu_0 sigma <kappa_n> = some L/R ratio *)
(* Actual R, L values depend on circuit normalization (Tanimoto convention). *)

Print["=== Comparison to Stage 5c Galerkin Cauer extraction ==="];
Print["From earlier Pade[2,2] of Galerkin moments (4 moments):"];
Print["  tau_1 = 16.74 us, tau_2 = 1.63 us"];
Print[""];
Print["Pade[3,3] from Galerkin gave converged tau_1 = 16.96 us"];
Print[""];

(* Assess Born vs Galerkin *)
Print["=== Assessment ==="];
Print["Born tau_1 (single pair): ", N[tauPade1 * 10^6, 6], " us"];
Print["Galerkin tau_1 (closed PEC + 1 dipole, Pade[3,3]): 16.96 us"];
Print["Ratio Born/Gal: ", N[tauPade1 * 10^6 / 16.96, 6]];
Print[""];

Print["NOTE: with only 3 moments (alpha_1, 2, 3), Pade[1,1] gives 1 pair."];
Print["For multi-pair extraction, need alpha_4, alpha_5, ... from Born series."];
Print["These require T_xx, T_yy at higher Born orders -- 12D, 15D, ... integrals."];
Print[""];
Print["Practical: combine Born alpha_1, alpha_2 (closed form, exact) +"];
Print["          Galerkin alpha_3, alpha_4 (cheap, basis-truncated) for full ladder."];
