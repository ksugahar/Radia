(* ::Package:: *)

(* hex_bem_order2.wls
   HCurl order 2 BEM on cuboid 5x2x1 mm Cu -> target tau_1 = 19.77 us.
   Adapted from hex_bem_full_hcurl.wls with ord=2.

   Order 2: 144 DOFs (vs 54 at order 1).
   K matrix: 144x144 = 10440 unique entries.

   Date: 2026-05-03
*)

Get["w:\\30_CauerLadderNetwork\\2026_04_01_\:9577\:65b9\:5f62CLN\\ngsolve_validation\\schoberl_zaglmayr_basis.wls"];

ClearAll[x, y, z, xp, yp, zp];
$MaxExtraPrecision = 50;

Print[StringRepeat["=", 70]];
Print["HCurl order 2 BEM on cuboid 5x2x1 mm Cu  (-> tau_1 target = 19.77 us)"];
Print[StringRepeat["=", 70]];

sigma = SetPrecision[58 10^6, 30];
mu0   = 4 Pi SetPrecision[10^-7, 30];
{aBox, bBox, cBox} = SetPrecision[{5 10^-3, 2 10^-3, 1 10^-3}, 30];

ord = 2;
shapesRef = SchoberlZaglmayrBasis`ShapesHEX[ord, Range[0, 7], x, y, z];
nDOF = Length[shapesRef];
Print["HCurl order ", ord, ": ", nDOF, " DOFs"];

fastIntegrate3D[poly_, vars_] := Module[{p = Expand[poly], rules},
  If[p === 0, Return[0]];
  rules = CoefficientRules[p, vars];
  Total[Function[r, r[[2]] / Times @@ Map[(# + 1) &, r[[1]]]] /@ rules]
];

componentMassMatrix[shapes_, comp_] := Module[{n = Length[shapes], M, exprs},
  exprs = Table[Expand[shapes[[i, 1, comp]]], {i, 1, n}];
  M = Table[
    If[j < i, 0, fastIntegrate3D[exprs[[i]] exprs[[j]], {x, y, z}]],
    {i, 1, n}, {j, 1, n}
  ];
  M + Transpose[M] - DiagonalMatrix[Diagonal[M]]
];

t0 = AbsoluteTime[];
MrefXX = componentMassMatrix[shapesRef, 1];
MrefYY = componentMassMatrix[shapesRef, 2];
MrefZZ = componentMassMatrix[shapesRef, 3];
Mphys = N[(bBox cBox / aBox) MrefXX +
          (aBox cBox / bBox) MrefYY +
          (aBox bBox / cBox) MrefZZ, 30];
Print["Mass matrices in ", Round[AbsoluteTime[] - t0, 0.1], " s"];

Kentry[i_, j_, WP_:15] := Module[{vi, vj, integrand, factor, scaled},
  vi = shapesRef[[i, 1]];
  vj = shapesRef[[j, 1]] /. {x -> xp, y -> yp, z -> zp};
  integrand = (vi[[1]] vj[[1]] / aBox^2 +
               vi[[2]] vj[[2]] / bBox^2 +
               vi[[3]] vj[[3]] / cBox^2);
  factor = aBox^2 bBox^2 cBox^2;
  scaled = factor / Sqrt[aBox^2 (x - xp)^2 + bBox^2 (y - yp)^2 + cBox^2 (z - zp)^2];
  (mu0 / (4 Pi)) NIntegrate[
    integrand scaled,
    {x, 0, 1}, {y, 0, 1}, {z, 0, 1},
    {xp, 0, 1}, {yp, 0, 1}, {zp, 0, 1},
    Method -> "GlobalAdaptive", WorkingPrecision -> WP,
    PrecisionGoal -> 4, AccuracyGoal -> 6, MaxRecursion -> 8
  ]
];

Print[];
Print["Building K matrix at WP=15 ..."];
nUnique = nDOF (nDOF + 1) / 2;
Print["  Unique entries: ", nUnique];
K = Table[0, {nDOF}, {nDOF}];
totalTime = 0;
counter = 0;
Do[
  counter++;
  t0 = AbsoluteTime[];
  K[[i, j]] = Kentry[i, j, 15];
  K[[j, i]] = K[[i, j]];
  totalTime += AbsoluteTime[] - t0;
  If[Mod[counter, 1000] == 0,
    Print["  ", counter, "/", nUnique, "; ", Round[totalTime/60, 0.1], " min"]
  ],
  {i, 1, nDOF}, {j, i, nDOF}
];
Print["K matrix in ", Round[totalTime/60, 0.1], " min"];
Kmat = N[K, 30];

Print[];
Print["Solving Eigensystem[{K, M}] ..."];
t0 = AbsoluteTime[];
{eigvals, eigvecs} = Eigensystem[{Kmat, Mphys}];
Print["Done in ", Round[AbsoluteTime[] - t0, 0.1], " s"];

eigvecsNorm = Table[
  Module[{v = eigvecs[[k]], vMv},
    vMv = Re[v . Mphys . v];
    If[vMv > 10^-50, v / Sqrt[vMv], 0 v]],
  {k, 1, Length[eigvals]}];

bVec = Table[
  -bBox cBox/2 * fastIntegrate3D[Expand[shapesRef[[i, 1, 1]] y], {x, y, z}] +
   aBox cBox/2 * fastIntegrate3D[Expand[shapesRef[[i, 1, 2]] x], {x, y, z}],
  {i, 1, nDOF}
];
gVec = Table[Re[eigvecsNorm[[k]] . bVec], {k, 1, Length[eigvals]}];
gMag2 = gVec^2;

realEigs = Re[eigvals];
tausUs = Table[
  If[realEigs[[k]] > 10^-30, 10^6 sigma realEigs[[k]], 0],
  {k, 1, Length[eigvals]}];

sortedByG = Reverse[Ordering[gMag2]];
Print[];
Print["Top 10 modes by |g_k|^2:"];
Do[
  Module[{idx = sortedByG[[k]]},
    Print["  rank ", k, ":  tau = ", N[tausUs[[idx]], 8], " us,  |g|^2 = ",
          N[gMag2[[idx]], 6]]],
  {k, 1, Min[10, Length[sortedByG]]}
];

Print[];
Print["============================================"];
Print["LEADING physical mode (HCurl order 2): tau = ",
      N[tausUs[[sortedByG[[1]]]], 10], " us"];
Print["Phase 7-11 reference (P=4 K=8 Cauer):    19.71 us"];
Print["Phase 7-11 reference (P=6):              19.77 us"];
Print["HCurl order 1 baseline (54 DOFs):        19.43 us"];
Print["============================================"];
