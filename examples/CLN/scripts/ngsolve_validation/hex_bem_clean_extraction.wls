(* ::Package:: *)

(* hex_bem_clean_extraction.wls
   Clean HCurl BEM analysis: filter spurious modes via projection magnitude,
   identify leading physical mode, attempt Cauer-I.

   Date: 2026-05-03
*)

Get["w:\\30_CauerLadderNetwork\\2026_04_01_\:9577\:65b9\:5f62CLN\\ngsolve_validation\\schoberl_zaglmayr_basis.wls"];

ClearAll[x, y, z, xp, yp, zp];
$MaxExtraPrecision = 50;

Print[StringRepeat["=", 70]];
Print["HCurl order 1 BEM: clean leading-mode extraction"];
Print[StringRepeat["=", 70]];

sigma = SetPrecision[58 10^6, 30];
mu0   = 4 Pi SetPrecision[10^-7, 30];
{aBox, bBox, cBox} = SetPrecision[{5 10^-3, 2 10^-3, 1 10^-3}, 30];
V = aBox bBox cBox;

ord = 1;
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
Print["M matrix in ", Round[AbsoluteTime[] - t0, 0.1], " s"];

Kentry[i_, j_, WP_] := Module[{vi, vj, integrand, factor, scaled},
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
Print["Building K matrix (WP=15) ..."];
nUnique = nDOF (nDOF + 1) / 2;
K = Table[0, {nDOF}, {nDOF}];
totalTime = 0;
counter = 0;
Do[
  counter++;
  t0 = AbsoluteTime[];
  K[[i, j]] = Kentry[i, j, 15];
  K[[j, i]] = K[[i, j]];
  totalTime += AbsoluteTime[] - t0,
  {i, 1, nDOF}, {j, i, nDOF}
];
Print["K matrix done in ", Round[totalTime/60, 0.1], " min"];
Kmat = N[K, 30];

Print[];
Print["Solving Eigensystem[{K, M}] ..."];
{eigvals, eigvecs} = Eigensystem[{Kmat, Mphys}];
Print["Done. Total ", Length[eigvals], " modes."];

(* M-orthonormalize eigenvectors *)
eigvecsNorm = Table[
  Module[{v = eigvecs[[k]], vMv},
    vMv = Re[v . Mphys . v];
    If[vMv > 10^-50, v / Sqrt[vMv], 0 v]
  ],
  {k, 1, Length[eigvals]}
];

(* Project A_ext = (-y/2, x/2, 0) onto each mode *)
bVec = Table[
  -bBox cBox/2 * fastIntegrate3D[Expand[shapesRef[[i, 1, 1]] y], {x, y, z}] +
   aBox cBox/2 * fastIntegrate3D[Expand[shapesRef[[i, 1, 2]] x], {x, y, z}],
  {i, 1, nDOF}
];

gVec = Table[Re[eigvecsNorm[[k]] . bVec], {k, 1, Length[eigvals]}];

(* Convert lambda to tau *)
realEigs = Re[eigvals];
tauUs = Table[
  If[realEigs[[k]] > 10^-30, 10^6 sigma * realEigs[[k]], 0],
  {k, 1, Length[eigvals]}
];

(* Sort by |g_k|^2 (excitation strength) *)
gMag2 = gVec^2;
sortedByG = Reverse[Ordering[gMag2]];

Print[];
Print["Top 20 modes ranked by |g_k|^2 (excitation strength):"];
Print["  rank   tau (us)        |g_k|^2          tau*g^2 (Foster weight)"];
Do[
  Module[{idx = sortedByG[[k]], tau = tauUs[[sortedByG[[k]]]], g2 = gMag2[[sortedByG[[k]]]]},
    Print["  ", k, "    ", N[tau, 8], "   ", N[g2, 6], "   ", N[tau g2, 6]]
  ],
  {k, 1, Min[20, Length[sortedByG]]}
];

(* Identify leading PHYSICAL mode: largest tau among those with non-trivial g *)
threshold = N[Max[gMag2] / 100, 30];   (* keep modes with |g|^2 > 1% of max *)
physicalModes = Select[
  Range[Length[eigvals]],
  (gMag2[[#]] > threshold) && (tauUs[[#]] > 0) &
];
sortedPhysical = SortBy[physicalModes, -tauUs[[#]] &];
Print[];
Print["Physical modes (|g|^2 > 1% of max):"];
Print["  Number of physical modes: ", Length[sortedPhysical]];
Print["  Top 5 by tau:"];
Do[
  Module[{idx = sortedPhysical[[k]]},
    Print["    tau_", k, " = ", N[tauUs[[idx]], 10], " us, |g|^2 = ", N[gMag2[[idx]], 6]]
  ],
  {k, 1, Min[5, Length[sortedPhysical]]}
];

Print[];
Print["============================================"];
Print["Leading physical tau (HCurl order 1 BEM, this run): ",
      If[Length[sortedPhysical] > 0,
         N[tauUs[[sortedPhysical[[1]]]], 8], "no physical mode"], " us"];
Print["Phase 7-11 reference (P=4 K=8 Cauer): 19.71 us"];
Print["Phase 7-11 reference (P=6): 19.77 us"];
Print["============================================"];

(* Foster decomposition for plotting Y(s) *)
Print[];
Print["Foster decomposition Y(s) = sigma Σ |g_k|^2 / (1 + s tau_k):"];
Print["  (only physical modes shown, scaled by max(g^2))"];
gNorm = gMag2 / Max[gMag2];
Do[
  Module[{idx = sortedPhysical[[k]]},
    Print["  k=", k, ":  tau = ", N[tauUs[[idx]], 8], " us,  weight = ", N[gNorm[[idx]], 6]]
  ],
  {k, 1, Min[10, Length[sortedPhysical]]}
];
