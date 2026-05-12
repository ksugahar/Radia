(* ::Package:: *)

(* hex_bem_hdiv_redo.wls
   Cuboid eddy current redo with CORRECT BC: J.n=0 on all 6 faces (vacuum).
   Run HDiv interior div-free orders 1, 2, 3, 4 for convergence confirmation.

   This corrects the methodology: the previous Phase 7-11 result (19.77 us)
   was for parity-restricted subspace, not isolated-in-vacuum BC.
   Correct Radia-use answer: see HDiv result below.

   Date: 2026-05-03
*)

ClearAll[x, y, z, xp, yp, zp];
$MaxExtraPrecision = 50;
$HistoryLength = 5;

Print[StringRepeat["=", 70]];
Print["Cuboid 5x2x1 mm Cu eddy current REDO with isolated-in-vacuum BC"];
Print["BC: J.n = 0 on all 6 faces (correct physics for Radia integration)"];
Print["Started: ", DateString[]];
Print[StringRepeat["=", 70]];

sigma = SetPrecision[58 10^6, 30];
mu0   = 4 Pi SetPrecision[10^-7, 30];
{aBox, bBox, cBox} = SetPrecision[{5 10^-3, 2 10^-3, 1 10^-3}, 30];

DuCrossDv[u_, v_, vars_] := Cross[D[u, {vars}], D[v, {vars}]];
wDuCrossDv[u_, v_, w_, vars_] := w Cross[D[u, {vars}], D[v, {vars}]];

IntegratedLegendre[n_, x_] := Module[{p3, p2, p1, j, result},
  result = {};
  If[n < 2, Return[{}]];
  p3 = 0; p2 = -1; p1 = x;
  Do[p3 = p2; p2 = p1;
     p1 = (2 j - 3)/j x p2 - (j - 3)/j p3;
     AppendTo[result, p1], {j, 2, n}];
  result];

LegendrePoly[n_, x_] := Module[{P, k},
  P = Table[0, {n + 1}];
  P[[1]] = 1; If[n >= 1, P[[2]] = x];
  Do[P[[k + 1]] = (2 - 1/k) x P[[k]] + (1/k - 1) P[[k - 1]], {k, 2, n}];
  P];

ShapesHDivHexDivFree[order_, x_, y_, z_] := Module[
  {p, xi, eta, zeta, LXi, LEta, LZeta, PXi, PEta, PZeta, shapes, vars},
  vars = {x, y, z};
  p = order;
  xi = 2 x - 1; eta = 2 y - 1; zeta = 2 z - 1;
  LXi = IntegratedLegendre[p + 1, xi];
  LEta = IntegratedLegendre[p + 1, eta];
  LZeta = IntegratedLegendre[p + 1, zeta];
  PXi = LegendrePoly[p, xi];
  PEta = LegendrePoly[p, eta];
  PZeta = LegendrePoly[p, zeta];
  shapes = {};
  Do[AppendTo[shapes, wDuCrossDv[LXi[[i+1]] LZeta[[k+1]], y, -4 PEta[[j+2]], vars]],
    {i, 0, p-1}, {j, 0, p-1}, {k, 0, p-1}];
  Do[AppendTo[shapes, wDuCrossDv[LEta[[j+1]] LZeta[[k+1]], x, 4 PXi[[i+2]], vars]],
    {i, 0, p-1}, {j, 0, p-1}, {k, 0, p-1}];
  Do[AppendTo[shapes, DuCrossDv[LEta[[j+1]] LZeta[[k+1]], 2 xi, vars]],
    {j, 0, p-1}, {k, 0, p-1}];
  Do[AppendTo[shapes, DuCrossDv[LXi[[i+1]] LZeta[[k+1]], 2 eta, vars]],
    {i, 0, p-1}, {k, 0, p-1}];
  Do[AppendTo[shapes, DuCrossDv[LXi[[i+1]] LEta[[j+1]], 2 zeta, vars]],
    {i, 0, p-1}, {j, 0, p-1}];
  shapes];

fastIntegrate3D[poly_, vars_] := Module[{p = Expand[poly], rules},
  If[p === 0, Return[0]];
  rules = CoefficientRules[p, vars];
  Total[Function[r, r[[2]] / Times @@ Map[(# + 1) &, r[[1]]]] /@ rules]];


KentryRobust[i_, j_, shapesRef_] := Module[{vi, vj, dotProduct, kernel, result},
  vi = shapesRef[[i]];
  vj = shapesRef[[j]] /. {x -> xp, y -> yp, z -> zp};
  dotProduct = vi[[1]] vj[[1]] aBox^2 + vi[[2]] vj[[2]] bBox^2 + vi[[3]] vj[[3]] cBox^2;
  kernel = 1 / Sqrt[aBox^2 (x - xp)^2 + bBox^2 (y - yp)^2 + cBox^2 (z - zp)^2];
  result = Quiet[NIntegrate[dotProduct kernel,
    {x, 0, 1}, {y, 0, 1}, {z, 0, 1}, {xp, 0, 1}, {yp, 0, 1}, {zp, 0, 1},
    Method -> {"GlobalAdaptive", "SymbolicProcessing" -> 0},
    WorkingPrecision -> 15,
    PrecisionGoal -> 3, AccuracyGoal -> 5, MaxRecursion -> 6],
    {NIntegrate::slwcon, NIntegrate::ncvb, NIntegrate::eincr, NIntegrate::precw}];
  (mu0 / (4 Pi)) result
];

runOrder[ord_] := Module[
  {shapesRef, n, MrefXX, MrefYY, MrefZZ, Mphys, K, t0, totalTime, counter, nUnique,
   eigvals, eigvecs, eigvecsNorm, bVec, gVec, gMag2, realEigs, tausUs,
   sortedByG, gMaxValue, gThresh, physicalIdx, sortedPhysical, Kmat, leadingTau},
  Print[];
  Print[StringRepeat["-", 60]];
  Print["HDiv div-free order ", ord, " (J.n=0 on all faces)"];
  Print[StringRepeat["-", 60]];

  shapesRef = ShapesHDivHexDivFree[ord, x, y, z];
  n = Length[shapesRef];
  Print["DOFs: ", n];

  t0 = AbsoluteTime[];
  MrefXX = Table[fastIntegrate3D[Expand[shapesRef[[i, 1]] shapesRef[[j, 1]]], {x, y, z}], {i, n}, {j, n}];
  MrefYY = Table[fastIntegrate3D[Expand[shapesRef[[i, 2]] shapesRef[[j, 2]]], {x, y, z}], {i, n}, {j, n}];
  MrefZZ = Table[fastIntegrate3D[Expand[shapesRef[[i, 3]] shapesRef[[j, 3]]], {x, y, z}], {i, n}, {j, n}];
  Mphys = N[(aBox/(bBox cBox)) MrefXX + (bBox/(aBox cBox)) MrefYY + (cBox/(aBox bBox)) MrefZZ, 30];
  Print["M matrix in ", Round[AbsoluteTime[] - t0, 0.1], " s"];

  nUnique = n (n + 1) / 2;
  Print["K matrix: ", nUnique, " unique entries ..."];
  K = Table[0, {n}, {n}];
  totalTime = 0;
  counter = 0;
  Do[
    counter++;
    t0 = AbsoluteTime[];
    K[[i, j]] = KentryRobust[i, j, shapesRef];
    K[[j, i]] = K[[i, j]];
    totalTime += AbsoluteTime[] - t0;
    If[Mod[counter, 500] == 0,
      Print["  ", counter, "/", nUnique, "; ", Round[totalTime/60, 0.1], " min"]],
    {i, 1, n}, {j, i, n}
  ];
  Print["K matrix in ", Round[totalTime/60, 0.1], " min"];

  Kmat = N[K, 30];
  {eigvals, eigvecs} = Eigensystem[{Kmat, Mphys}];

  eigvecsNorm = Table[
    Module[{v = eigvecs[[k]], vMv},
      vMv = Re[v . Mphys . v];
      If[vMv > 10^-50, v / Sqrt[vMv], 0 v]],
    {k, 1, Length[eigvals]}];

  bVec = Table[
    Module[{phi = shapesRef[[i]]},
      (-aBox/2) fastIntegrate3D[Expand[phi[[1]] y aBox bBox cBox / (bBox cBox)], {x, y, z}] +
      (bBox/2)  fastIntegrate3D[Expand[phi[[2]] x aBox bBox cBox / (aBox cBox)], {x, y, z}]],
    {i, 1, n}];
  bVec = N[bVec, 30];
  gVec = Table[Re[eigvecsNorm[[k]] . bVec], {k, 1, Length[eigvals]}];
  gMag2 = gVec^2;

  realEigs = Re[eigvals];
  tausUs = Table[
    If[realEigs[[k]] > 10^-30, 10^6 sigma realEigs[[k]], 0],
    {k, 1, Length[eigvals]}];

  (* Find leading by |g|^2 (the physical excitation amplitude) *)
  sortedByG = Reverse[Ordering[gMag2]];
  Print[];
  Print["Top 5 modes by |g|^2 (physical excitation strength):"];
  Do[
    Module[{idx = sortedByG[[k]]},
      Print["  rank ", k, ":  tau = ", N[tausUs[[idx]], 12], " us, |g|^2 = ",
            N[gMag2[[idx]], 6]]],
    {k, 1, Min[5, Length[sortedByG]]}];

  leadingTau = tausUs[[sortedByG[[1]]]];
  Print[];
  Print["LEADING tau (max |g|^2): ", N[leadingTau, 12], " us"];

  {ord, n, leadingTau, gMag2[[sortedByG[[1]]]]}
];


(* === Run orders 1, 2, 3, 4 === *)
results = {};
Do[
  results = Append[results, runOrder[ord]],
  {ord, 1, 4}
];

Print[];
Print[StringRepeat["=", 70]];
Print["CONVERGENCE STUDY: HDiv div-free (J.n=0) cuboid 5x2x1 mm Cu"];
Print[StringRepeat["=", 70]];
Print["  order | DOFs |  Leading tau (us)  | |g|^2 (excitation)"];
Do[
  Print["    ", results[[k, 1]], "   |  ", results[[k, 2]], "  | ",
        N[results[[k, 3]], 12], " | ", N[results[[k, 4]], 4]],
  {k, 1, Length[results]}];

Print[];
Print["This is the CORRECT physical setup for Radia integration (vacuum BC)"];
Print["Phase 7-11 reference 19.77 us was for parity-restricted subspace ;"];
Print["it does NOT correspond to isolated-conductor-in-vacuum eddy current."];
Print[];
Print["Finished: ", DateString[]];
