(* High-precision Cauer extraction from spectral BEM Foster moments.

   Uses moments alpha_n from WP=50 P=6 spectral BEM run (hardcoded here).
   Extracts CLN-I (s=0) continued fraction at WP=50 and shows:
     - Cauer p_k coefficients to many digits
     - Cauer truncation -> recovered tau_k vs Foster top
     - Convergence as K increases

   Date: 2026-05-01
*)

$HistoryLength = 0;
prec = 50;

Print["================================================================"];
Print["High-precision Cauer extraction (WP=50)"];
Print["from spectral BEM P=6 moments"];
Print["================================================================"];
Print[""];

(* Hardcode WP=50 P=6 alpha values (from spectral_v3_wp50.log) *)
alpha1 = SetPrecision[-0.0000440346570278169352257847180889677, prec];  (* exact via Parseval *)
alpha2 = SetPrecision[7.4783240521932061573535486933*^-10, prec];
alpha3 = SetPrecision[-1.38987444698506133728520045089107882524730075*^-14, prec];
alpha4 = SetPrecision[2.671196698538943965550262235483672977442*^-19, prec];
alpha5 = SetPrecision[-5.2088290905449930546402512476415822*^-22, prec];
(* ... up to alpha_n we have. P=6 had 12 moments computed, but we don't have them all in log.
   For now use 5; extend when P=8 completes. *)

alphasList = {alpha1, alpha2, alpha3, alpha4, alpha5};
NMoments = Length[alphasList];

Print["Moments (alpha_n, n=1..", NMoments, "):"];
Do[
  Print["  alpha_", n, " = ", N[alphasList[[n]], 16]],
  {n, NMoments}];
Print[""];

(* === Convert alphas to Taylor coefficients of chi(s) ===
   chi(s) = alpha_0 + alpha_1 s + alpha_2 s^2 + ...
   But alpha_0 = sum a_k = M_0 not directly given by alpha_1, alpha_2, ...
   Trick: Pade[N-1, N] from alpha_1..alpha_{2N-1} reconstructs alpha_0 implicitly.

   For Cauer continued fraction starting from M_0:
   We need alpha_0. From Foster: alpha_0 = sum a_k.
   For our problem, we need the explicit M_0 from BEM.
   For now, derive M_0 from Pade[1,1]: M_0_pade = alpha_1^2/alpha_2 (under-estimates).
   Better: take it from BEM directly. Hardcoding from earlier output:
*)
M0 = SetPrecision[3.05897258746351459800805934293582037956, prec];  (* from Phase 7 WP=30 *)

Print["alpha_0 = M_0 (from BEM) = ", N[M0, 16]];
Print[""];

(* Taylor coefficients c_n of chi(s) = sum c_n s^n where c_n = alpha_n
   (alpha_0 = M_0; alpha_1, alpha_2, ... already in alphasList) *)
cTaylor = Prepend[alphasList, M0];

(* === CLN-I (s=0) continued fraction extraction ===
   chi(s) = 1 / (p_1 + s/(p_2 + s/(p_3 + ...)))
   Iteration: invert Taylor, peel leading term, divide by s. *)
invertTaylor[c_List] := Module[{n = Length[c], e},
  e = ConstantArray[0, n];
  e[[1]] = 1/c[[1]];
  Do[
    e[[k+1]] = -Sum[c[[j+1]] e[[k-j+1]], {j, 1, k}]/c[[1]],
    {k, 1, n - 1}];
  e];

cauerExtract[cTaylorIn_List, maxStages_] := Module[
  {c = cTaylorIn, p = {}, e},
  Do[
    If[Length[c] < 2 || Abs[c[[1]]] < 10^(-prec + 5), Break[]];
    e = invertTaylor[c];
    AppendTo[p, e[[1]]];
    c = Drop[e, 1],
    {k, maxStages}];
  p];

pCauer = cauerExtract[cTaylor, NMoments];
Print["Cauer p_k (", Length[pCauer], " stages):"];
Do[
  Print["  p_", k, " = ", N[pCauer[[k]], 18]],
  {k, Length[pCauer]}];
Print[""];

(* === Cauer truncation -> recovered tau_k === *)
(* Build polynomial: chi^(K)(s) = N(s)/D(s) from K Cauer stages.
   chi^(K) = 1/(p_1 + s/(p_2 + s/(... + s/p_K)))
   Build inside out. *)
cauerToRational[pList_List, K_Integer] := Module[
  {N, D, k, sN, padded, max, Dnew},
  N = {1};
  D = {pList[[K]]};
  Do[
    sN = Prepend[N, 0];
    Dnew = pList[[k]] D;
    max = Max[Length[Dnew], Length[sN]];
    padded = ConstantArray[0, max];
    padded[[1 ;; Length[Dnew]]] += Dnew;
    padded[[1 ;; Length[sN]]] += sN;
    {N, D} = {D, padded},
    {k, K - 1, 1, -1}];
  {N, D}];

extractRealNegPoles[Dcoeffs_] := Module[{s, poly, roots, real},
  poly = Sum[Dcoeffs[[k+1]] s^k, {k, 0, Length[Dcoeffs] - 1}];
  roots = NSolve[poly == 0, s, WorkingPrecision -> prec];
  real = s /. roots;
  Cases[real, x_ /; NumericQ[x] && Re[x] < 0 && Abs[Im[x]] < 10^(-prec/3) (Abs[Re[x]] + 1)]];

Print["Cauer truncation -> Pade poles (recovered tau_k):"];
Print[StringPadRight["K", 4], StringPadRight["recovered tau [us]", 80]];
Do[
  If[K > Length[pCauer], Break[]];
  Module[{Nc, Dc, poles, taus},
    {Nc, Dc} = cauerToRational[pCauer, K];
    poles = extractRealNegPoles[Dc];
    taus = Sort[-1/Re[#] & /@ poles, Greater];
    Print[StringPadRight[ToString[K], 4],
      StringJoin[
        Riffle[ToString[N[# 10^6, 12]] & /@ taus, ", "]]]],
  {K, 2, Length[pCauer]}];
Print[""];

(* === Compare with closed Pade[1,1] ===*)
tau1Pade11 = alpha2/(-alpha1);
Print["Reference: tau_1 (Pade[1,1]) = ", N[tau1Pade11 10^6, 16], " us"];
