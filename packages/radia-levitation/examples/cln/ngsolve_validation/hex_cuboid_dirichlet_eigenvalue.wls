(* ::Package:: *)

(* hex_cuboid_dirichlet_eigenvalue.wls
   Single hex element with Dirichlet H(curl) BC: tangential A = 0 on all 6 faces.
   This removes ALL edge and face DOFs, leaving only cell-interior DOFs.

   Test convergence: order 1, 2, 3, 4 — leading tau_1 should converge to ~19.77 us.

   Cell DOF count formula (HEX H(curl) order p):
     gradient cell:   p^3
     rotations cell:  2 p^3
     missing cell:    3 p^2
     total:           3 p^3 + 3 p^2

   Date: 2026-05-02
*)

Get["w:\\30_CauerLadderNetwork\\2026_04_01_\:9577\:65b9\:5f62CLN\\ngsolve_validation\\schoberl_zaglmayr_basis.wls"];

ClearAll[x, y, z];
$MaxExtraPrecision = 100;

(* === Fast monomial integration === *)
fastIntegrate3D[poly_] := Module[{p = Expand[poly], rules},
  If[p === 0, Return[0]];
  rules = CoefficientRules[p, {x, y, z}];
  Total[
    Function[r, r[[2]] / ((r[[1, 1]] + 1) (r[[1, 2]] + 1) (r[[1, 3]] + 1))]
      /@ rules
  ]
];

componentMatrixSym[shapes_, slot_, comp_] := Module[{n = Length[shapes], M, exprs},
  exprs = Table[Expand[shapes[[i, slot, comp]]], {i, 1, n}];
  M = Table[
    If[j < i, 0, fastIntegrate3D[exprs[[i]] exprs[[j]]]],
    {i, 1, n}, {j, 1, n}
  ];
  M + Transpose[M] - DiagonalMatrix[Diagonal[M]]
];

(* === Identify cell-only (interior) DOFs ===
   In ShapesHEX ordering for vnums = {0..7}:
     indices 1..12  = lowest-order edges (12)
     12 high-order edge gradients per p (12 p)
     6 faces × ... face shapes
     cell shapes at the END
   Cell DOF count at order p: 3 p^3 + 3 p^2
   Total DOF count at order p: complicated; safer to compute cell DOFs from formula.
*)

cellDofCount[p_] := 3 p^3 + 3 p^2;

(* === Solve for given order === *)
solveOrder[ord_] := Module[{shapes, n, MrefX, MrefY, MrefZ, SrefX, SrefY, SrefZ,
                            sigma, mu0, aBox, bBox, cBox, abc,
                            Mphys, Sphys, nCell, cellIdx, Mc, Sc, eigvals,
                            sortedEigs, t0},
  Print[StringRepeat["=", 60]];
  Print["HEX Dirichlet eigenvalue, order = ", ord];
  Print[StringRepeat["=", 60]];

  shapes = SchoberlZaglmayrBasis`ShapesHEX[ord, Range[0, 7], x, y, z];
  n = Length[shapes];
  Print["  Total DOFs: ", n];

  nCell = cellDofCount[ord];
  cellIdx = Range[n - nCell + 1, n];   (* cell DOFs are at the end *)
  Print["  Cell DOFs (Dirichlet kept): ", nCell, " (indices ", n - nCell + 1, "..", n, ")"];

  t0 = AbsoluteTime[];
  MrefX = componentMatrixSym[shapes, 1, 1];
  MrefY = componentMatrixSym[shapes, 1, 2];
  MrefZ = componentMatrixSym[shapes, 1, 3];
  SrefX = componentMatrixSym[shapes, 2, 1];
  SrefY = componentMatrixSym[shapes, 2, 2];
  SrefZ = componentMatrixSym[shapes, 2, 3];
  Print["  All 6 component matrices built in ", Round[AbsoluteTime[] - t0, 0.1], " s"];

  sigma = SetPrecision[58 10^6, 30];
  mu0   = 4 Pi SetPrecision[10^-7, 30];
  {aBox, bBox, cBox} = SetPrecision[{5 10^-3, 2 10^-3, 1 10^-3}, 30];
  abc = aBox bBox cBox;

  Mphys = (bBox cBox / aBox) MrefX + (aBox cBox / bBox) MrefY + (aBox bBox / cBox) MrefZ;
  Sphys = (aBox / (bBox cBox)) SrefX + (bBox / (aBox cBox)) SrefY + (cBox / (aBox bBox)) SrefZ;

  (* Restrict to cell DOFs *)
  Mc = N[Mphys[[cellIdx, cellIdx]], 30];
  Sc = N[Sphys[[cellIdx, cellIdx]], 30];

  t0 = AbsoluteTime[];
  eigvals = Eigenvalues[{Sc, Mc}];
  Print["  Eigensystem solved in ", Round[AbsoluteTime[] - t0, 0.1], " s"];

  (* Filter: keep real positive *)
  sortedEigs = Sort[Select[Re[eigvals], # > 0 &]];

  Print["  Number of positive real eigenvalues: ", Length[sortedEigs]];
  Print["  First few lambda (ascending):"];
  Do[Print["    lambda_", i, " = ", N[sortedEigs[[i]], 12]],
     {i, 1, Min[5, Length[sortedEigs]]}];
  Print["  First few tau = sigma mu0 / lambda [us]:"];
  Do[Print["    tau_", i, " = ", N[10^6 sigma mu0 / sortedEigs[[i]], 12], " us"],
     {i, 1, Min[5, Length[sortedEigs]]}];
  Print[];

  Return[{ord, nCell, sortedEigs}];
];

(* === Run convergence test === *)
results = {};
Do[
  AppendTo[results, solveOrder[ord]],
  {ord, 1, 4}
];

(* === Convergence summary === *)
Print[];
Print[StringRepeat["=", 60]];
Print["Convergence of leading tau_1 toward Phase 7-11 reference"];
Print[StringRepeat["=", 60]];
Print["Reference (Phase 7-11 spectral, parity-restricted): tau_1 = 19.77 us"];
Print[];
sigmaUsed = SetPrecision[58 10^6, 30];
mu0Used = 4 Pi SetPrecision[10^-7, 30];
Do[
  Module[{ord = results[[i, 1]], nCell = results[[i, 2]],
          eigs = results[[i, 3]], tau1},
    tau1 = If[Length[eigs] >= 1, N[10^6 sigmaUsed mu0Used / eigs[[1]], 8], "no positive eig"];
    Print["  order ", ord, " (", nCell, " cell DOFs): tau_1 = ", tau1, " us"]
  ],
  {i, 1, Length[results]}
];
