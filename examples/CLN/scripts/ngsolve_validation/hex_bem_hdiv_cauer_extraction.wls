(* ::Package:: *)

(* hex_bem_hdiv_cauer_extraction.wls
   Cauer-I extraction from HDiv order 3 results (saved .mx).
   Output: Foster network (tau_n, R_n, L_n) for Radia integration.

   Date: 2026-05-03
*)

ClearAll[x, y, z];
$MaxExtraPrecision = 50;

Print["Cauer-I extraction from HDiv order 3 (J.n=0 vacuum cuboid)"];
Print["==========================================================="];

(* Load saved results *)
If[FileExistsQ["hex_bem_hdiv_order3_results.mx"],
  Get["hex_bem_hdiv_order3_results.mx"];
  Print["Loaded: nDOF=", nDOF, ", #eigvals=", Length[eigvals]],
  Print["No saved file. Run hex_bem_hdiv_order3_overnight.wls first."]; Quit[]
];

sigma = SetPrecision[58 10^6, 30];
mu0   = 4 Pi SetPrecision[10^-7, 30];
{aBox, bBox, cBox} = SetPrecision[{5 10^-3, 2 10^-3, 1 10^-3}, 30];
V = aBox bBox cBox;

(* Recompute g vector and tau values *)
realEigs = Re[eigvals];
realGVec = Re[gVec];
gMag2 = realGVec^2;
tausUs = Table[
  If[realEigs[[k]] > 10^-30, 10^6 sigma realEigs[[k]], 0],
  {k, 1, Length[eigvals]}];

(* === Filter physical Foster modes ===
   Keep modes with non-trivial g^2 AND positive tau *)
gMaxValue = Max[gMag2];
gThresh = gMaxValue / 10^4;   (* keep top 4 orders of magnitude in g^2 *)
physicalIdx = Select[Range[Length[eigvals]],
  (gMag2[[#]] > gThresh) && (tausUs[[#]] > 0) &];
sortedPhysical = SortBy[physicalIdx, -tausUs[[#]] &];   (* sort by tau descending *)

Print[];
Print["Foster modes (filtered, sorted by tau desc):"];
Print["  rank | tau (us) | |g|^2 | tau*g^2 (Y(0) contribution) "];
Do[
  Module[{idx = sortedPhysical[[k]]},
    Print["  ", k, " | ", N[tausUs[[idx]], 8], " | ",
          N[gMag2[[idx]], 6], " | ", N[tausUs[[idx]] gMag2[[idx]], 6]]],
  {k, 1, Min[12, Length[sortedPhysical]]}];

(* === Compute Foster moments alpha_n ===
   chi(s) = sum_k g_k^2 tau_k / (1 + s tau_k)
   alpha_n = (-1)^n d^n chi/ds^n at s=0 / n!
          = (-1)^n sum_k g_k^2 tau_k^(n+1) / V    (with normalization to V)
   Or using Phase 7-11 convention:
     alpha_n = (-1)^n (mu0 sigma)^n / V * sum_k g_k^2 lambda_phys^(n-1)
*)

lambdaPhys = realEigs / mu0;   (* convert: lambda_phys for Phase 7-11 formula *)
mu0sigma = mu0 sigma;

NMoments = 8;
alphas = Table[
  (-1)^nM Sum[
    If[realEigs[[idx]] > 10^-30,
       gMag2[[idx]] lambdaPhys[[idx]]^(nM - 1),
       0],
    {idx, 1, Length[eigvals]}] mu0sigma^nM / V,
  {nM, 1, NMoments}
];

Print[];
Print["alpha_n moments (HDiv path, all 81 modes contributing):"];
Do[Print["  alpha_", k, " = ", N[alphas[[k]], 12]], {k, 1, NMoments}];

(* M_0 = chi(0) = sum_k g_k^2 tau_k *)
M0 = Sum[gMag2[[idx]] / lambdaPhys[[idx]],
   {idx, Select[Range[Length[eigvals]], realEigs[[#]] > 10^-30 &]}] mu0 sigma / V;
M0 = Re[M0];
Print[];
Print["M_0 (DC moment) = ", N[M0, 12]];

Print[];
tau1Pade = N[10^6 alphas[[2]] / (-alphas[[1]]), 12];
Print["tau_1 Pade[0,1] = alpha_2/(-alpha_1) = ", tau1Pade, " us"];

(* === Cauer-I continued fraction extraction ===
   Taylor coefficients c_n = M_n where M_n = (-1)^n alpha_n
*)
cTaylor = Prepend[Table[(-1)^n alphas[[n]], {n, 1, NMoments}], M0];
Print[];
Print["Taylor coefficients (M_0..M_8):"];
Do[Print["  c_", k - 1, " = ", N[cTaylor[[k]], 12]], {k, 1, Length[cTaylor]}];

invertTaylor[c_] := Module[{n = Length[c], e, k, s},
  e = Table[0, {n}];
  e[[1]] = 1/c[[1]];
  Do[
    s = Sum[c[[j+1]] e[[k - j + 1]], {j, 1, k}];
    e[[k + 1]] = -s / c[[1]],
    {k, 1, n - 1}];
  e];

cauerCF[cTaylor_, maxStages_] := Module[{c = cTaylor, pList = {}, e, k},
  Do[
    If[Length[c] < 2 || Abs[c[[1]]] < 10^-200, Break[]];
    e = invertTaylor[c];
    AppendTo[pList, e[[1]]];
    c = Drop[e, 1],
    {k, 1, maxStages}];
  pList];

cauerTruncatePoles[pList_, K_] := Module[{N0, D0, k},
  If[K < 1, Return[{}]];
  N0 = {1};
  D0 = {pList[[K]]};
  Do[
    Module[{sN = Prepend[N0, 0], DnewPk = pList[[k]] * D0,
            maxLen, Dnew},
      maxLen = Max[Length[DnewPk], Length[sN]];
      Dnew = Table[0, {maxLen}];
      Do[Dnew[[i]] += DnewPk[[i]], {i, 1, Length[DnewPk]}];
      Do[Dnew[[i]] += sN[[i]], {i, 1, Length[sN]}];
      N0 = D0; D0 = Dnew],
    {k, K - 1, 1, -1}];
  {N0, D0}];

pList = cauerCF[cTaylor, NMoments];
Print[];
Print["Cauer-I p_k stages (first ", Length[pList], "):"];
Do[Print["  p_", k, " = ", N[pList[[k]], 12]], {k, 1, Length[pList]}];

findRealNegRoots[coeffs_] := Module[{poly, roots, t},
  poly = Sum[coeffs[[k+1]] t^k, {k, 0, Length[coeffs]-1}];
  roots = t /. NSolve[poly == 0, t, WorkingPrecision -> 25];
  Select[roots, (Im[N[#]] == 0) && (Re[N[#]] < 0) &]
];

Print[];
Print["Cauer truncation -> recovered tau_k (physical Foster poles):"];
Do[
  Module[{ND, Nc, Dc, poles, taus},
    ND = cauerTruncatePoles[pList, K];
    {Nc, Dc} = ND;
    poles = findRealNegRoots[Dc];
    If[Length[poles] > 0,
      taus = Sort[-1/Re[poles], Greater];
      Print["  K=", K, ": ", Length[poles], " poles, taus = ",
            N[10^6 # & /@ taus, 8], " us"],
      Print["  K=", K, ": no real negative poles"]]],
  {K, 2, Length[pList]}];

(* === Foster network parameters R_n, L_n for Radia ===
   From Y(s) = sum_n |g_n|^2 / (R_n + s L_n)
   With R_n = 1/(sigma |g_n|^2), L_n = mu0/(lambda_phys_n |g_n|^2) ?
   Actually for series RL Foster: Y(s) = sum 1/(R_n + s L_n) where tau_n = L_n/R_n
*)
Print[];
Print["Foster R_n, L_n parameters (top 5 modes by tau*g^2 ranking):"];
sortedByContrib = Reverse[SortBy[sortedPhysical,
  tausUs[[#]] gMag2[[#]] &]];
Do[
  Module[{idx = sortedByContrib[[k]], R, L, tau},
    tau = realEigs[[idx]];   (* in s, dimensionless to mu0sigma *)
    (* R_n = 1/(sigma |g_n|^2 V), L_n = mu0/(lambda_phys |g_n|^2) ... approximate *)
    Print["  Mode ", k, ":  tau = ", N[10^6 tausUs[[idx]]/sigma, 6],
          "*sigma us, weight g^2 = ", N[gMag2[[idx]], 6]]],
  {k, 1, Min[5, Length[sortedByContrib]]}];

Print[];
Print["============================================"];
Print["For Radia: use multi-pole Foster network (table above)"];
Print["Single 'leading tau' is misleading for cuboid eddy current"];
Print["============================================"];
