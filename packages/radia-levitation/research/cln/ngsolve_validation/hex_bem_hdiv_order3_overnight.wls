(* ::Package:: *)

(* hex_bem_hdiv_order3_overnight.wls
   Overnight HDiv order 3 div-free BEM run on cuboid 5x2x1 mm Cu.
   Robust setup with periodic checkpointing.

   Order 3: 2*27 + 3*9 = 81 div-free interior DOFs (= Phase 7-11 P=3 equivalent).
   K matrix: 81x81 = 3321 unique entries.

   At WP=18, ~5-15 sec/entry expected -> 5-12 hours total.

   Strategy:
   - Each entry: 6D NIntegrate WP=18 with timeout fallback to WP=15
   - Save K matrix to .mx every 200 entries (so we can resume if needed)
   - At end: full Foster + Cauer analysis

   Date: 2026-05-03 (overnight run)
*)

ClearAll[x, y, z, xp, yp, zp];
$MaxExtraPrecision = 60;
$HistoryLength = 5;   (* save memory *)

Print[StringRepeat["=", 70]];
Print["HDiv order 3 div-free BEM overnight run (WP=18, ~hours)"];
Print["Started at: ", DateString[]];
Print[StringRepeat["=", 70]];

sigma = SetPrecision[58 10^6, 30];
mu0   = 4 Pi SetPrecision[10^-7, 30];
{aBox, bBox, cBox} = SetPrecision[{5 10^-3, 2 10^-3, 1 10^-3}, 30];

(* === Helpers === *)
DuCrossDv[u_, v_, vars_] := Cross[D[u, {vars}], D[v, {vars}]];
wDuCrossDv[u_, v_, w_, vars_] := w Cross[D[u, {vars}], D[v, {vars}]];

IntegratedLegendre[n_, x_] := Module[{p3, p2, p1, j, result},
  result = {};
  If[n < 2, Return[{}]];
  p3 = 0; p2 = -1; p1 = x;
  Do[
    p3 = p2; p2 = p1;
    p1 = (2 j - 3)/j * x * p2 - (j - 3)/j * p3;
    AppendTo[result, p1],
    {j, 2, n}
  ];
  result
];

LegendrePoly[n_, x_] := Module[{P, k},
  P = Table[0, {n + 1}];
  P[[1]] = 1;
  If[n >= 1, P[[2]] = x];
  Do[
    P[[k + 1]] = (2 - 1/k) x P[[k]] + (1/k - 1) P[[k - 1]],
    {k, 2, n}
  ];
  P
];

ShapesHDivHexDivFree[order_, x_, y_, z_] := Module[
  {p, xi, eta, zeta, LXi, LEta, LZeta, PXi, PEta, PZeta, shapes, vars},
  vars = {x, y, z};
  p = order;
  xi = 2 x - 1; eta = 2 y - 1; zeta = 2 z - 1;
  LXi   = IntegratedLegendre[p + 1, xi];
  LEta  = IntegratedLegendre[p + 1, eta];
  LZeta = IntegratedLegendre[p + 1, zeta];
  PXi   = LegendrePoly[p, xi];
  PEta  = LegendrePoly[p, eta];
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
  shapes
];

fastIntegrate3D[poly_, vars_] := Module[{p = Expand[poly], rules},
  If[p === 0, Return[0]];
  rules = CoefficientRules[p, vars];
  Total[Function[r, r[[2]] / Times @@ Map[(# + 1) &, r[[1]]]] /@ rules]
];

(* === Setup order 3 === *)
ord = 3;
shapesRef = ShapesHDivHexDivFree[ord, x, y, z];
nDOF = Length[shapesRef];
Print["HDiv order ", ord, " div-free DOFs: ", nDOF];

(* Mass matrix (analytical, fast) *)
t0 = AbsoluteTime[];
Print[];
Print["Computing reference component mass matrices ..."];
MrefXX = Table[fastIntegrate3D[Expand[shapesRef[[i, 1]] shapesRef[[j, 1]]], {x, y, z}],
               {i, nDOF}, {j, nDOF}];
MrefYY = Table[fastIntegrate3D[Expand[shapesRef[[i, 2]] shapesRef[[j, 2]]], {x, y, z}],
               {i, nDOF}, {j, nDOF}];
MrefZZ = Table[fastIntegrate3D[Expand[shapesRef[[i, 3]] shapesRef[[j, 3]]], {x, y, z}],
               {i, nDOF}, {j, nDOF}];
Mphys = N[(aBox/(bBox cBox)) MrefXX + (bBox/(aBox cBox)) MrefYY + (cBox/(aBox bBox)) MrefZZ, 30];
Print["  Mass done in ", Round[AbsoluteTime[] - t0, 0.1], " s"];

(* === K matrix entry with robust NIntegrate === *)
KentryRobust[i_, j_] := Module[{vi, vj, dotProduct, kernel, result},
  vi = shapesRef[[i]];
  vj = shapesRef[[j]] /. {x -> xp, y -> yp, z -> zp};
  dotProduct = vi[[1]] vj[[1]] aBox^2 + vi[[2]] vj[[2]] bBox^2 + vi[[3]] vj[[3]] cBox^2;
  kernel = 1 / Sqrt[aBox^2 (x - xp)^2 + bBox^2 (y - yp)^2 + cBox^2 (z - zp)^2];
  (* WP=15 with relaxed goals - matches HCurl 12-sec/1485-entries baseline *)
  result = Quiet[NIntegrate[dotProduct kernel,
    {x, 0, 1}, {y, 0, 1}, {z, 0, 1}, {xp, 0, 1}, {yp, 0, 1}, {zp, 0, 1},
    Method -> {"GlobalAdaptive", "SymbolicProcessing" -> 0},
    WorkingPrecision -> 15,
    PrecisionGoal -> 3, AccuracyGoal -> 5, MaxRecursion -> 6],
    {NIntegrate::slwcon, NIntegrate::ncvb, NIntegrate::eincr, NIntegrate::precw}
  ];
  (mu0 / (4 Pi)) result
];

(* === Build K matrix with periodic save === *)
nUnique = nDOF (nDOF + 1) / 2;
Print[];
Print["Building K matrix: ", nUnique, " unique entries (overnight expected)"];
Print["Started at: ", DateString[]];

K = Table[0, {nDOF}, {nDOF}];
totalTime = 0;
counter = 0;
saveFile = "K_hdiv_order3_partial.mx";

Do[
  counter++;
  t0 = AbsoluteTime[];
  K[[i, j]] = KentryRobust[i, j];
  K[[j, i]] = K[[i, j]];
  totalTime += AbsoluteTime[] - t0;
  If[Mod[counter, 100] == 0,
    Print["  ", counter, "/", nUnique, "; ",
          Round[totalTime/60, 0.1], " min; ",
          DateString[]];
    (* Save partial results *)
    DumpSave[saveFile, {K, counter, totalTime}];
  ],
  {i, 1, nDOF}, {j, i, nDOF}
];

Print[];
Print["K matrix done in ", Round[totalTime/60, 0.1], " min (", DateString[], ")"];
DumpSave["K_hdiv_order3_final.mx", K];

(* === Eigensystem === *)
Kmat = N[K, 30];
Print[];
Print["Solving Eigensystem ..."];
{eigvals, eigvecs} = Eigensystem[{Kmat, Mphys}];
Print["Done. ", Length[eigvals], " modes."];

eigvecsNorm = Table[
  Module[{v = eigvecs[[k]], vMv},
    vMv = Re[v . Mphys . v];
    If[vMv > 10^-50, v / Sqrt[vMv], 0 v]],
  {k, 1, Length[eigvals]}];

(* Project onto A_ext *)
bVec = Table[
  Module[{phi = shapesRef[[i]]},
    (-aBox/2) fastIntegrate3D[Expand[phi[[1]] y aBox bBox cBox / (bBox cBox)], {x, y, z}] +
    (bBox/2)  fastIntegrate3D[Expand[phi[[2]] x aBox bBox cBox / (aBox cBox)], {x, y, z}]
  ],
  {i, 1, nDOF}
];
bVec = N[bVec, 30];
gVec = Table[Re[eigvecsNorm[[k]] . bVec], {k, 1, Length[eigvals]}];
gMag2 = gVec^2;

realEigs = Re[eigvals];
tausUs = Table[
  If[realEigs[[k]] > 10^-30, 10^6 sigma realEigs[[k]], 0],
  {k, 1, Length[eigvals]}];

(* Filter: significant g_k AND positive tau *)
gMaxValue = Max[gMag2];
gThresh = gMaxValue / 10^6;
physicalIdx = Select[Range[Length[eigvals]],
  (gMag2[[#]] > gThresh) && (tausUs[[#]] > 0) &];
sortedPhysical = SortBy[physicalIdx, -tausUs[[#]] &];

Print[];
Print["===== Top 15 PHYSICAL modes (g^2 > max/10^6 and tau > 0) ====="];
Do[
  Module[{idx = sortedPhysical[[k]]},
    Print["  rank ", k, ":  tau = ", N[tausUs[[idx]], 12], " us,  |g|^2 = ",
          N[gMag2[[idx]], 6]]],
  {k, 1, Min[15, Length[sortedPhysical]]}
];

Print[];
Print[StringRepeat["=", 70]];
Print["LEADING physical mode: tau = ",
      N[tausUs[[sortedPhysical[[1]]]], 12], " us"];
Print["Phase 7-11 reference (P=4 K=8): 19.71 us"];
Print["Phase 7-11 reference (P=6):     19.77 us"];
Print[StringRepeat["=", 70]];
Print["Finished at: ", DateString[]];

(* Save final results *)
DumpSave["hex_bem_hdiv_order3_results.mx",
  {nDOF, K, Mphys, eigvals, eigvecs, gVec, tausUs, sortedPhysical}];
Print["Saved to hex_bem_hdiv_order3_results.mx"];
