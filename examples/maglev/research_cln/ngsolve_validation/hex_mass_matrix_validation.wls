(* ::Package:: *)

(* hex_mass_matrix_validation.wls
   Compute HEX H(curl) mass matrix M_ij = int_[0,1]^3 phi_i . phi_j dV
   symbolically using Mathematica, for validation against NGSolve numerical.

   Output: Table of diagonal entries, sparsity pattern, and a few off-diagonals.

   Date: 2026-05-02
*)

Get["w:\\30_CauerLadderNetwork\\2026_04_01_\:9577\:65b9\:5f62CLN\\ngsolve_validation\\schoberl_zaglmayr_basis.wls"];

ClearAll[x, y, z];

(* === HEX Mass matrix at given order === *)
HexMassMatrix[order_, opts___] := Module[{shapes, n, M, t0},
  Print[StringRepeat["=", 64]];
  Print["HEX H(curl) Mass matrix: order = ", order];
  Print[StringRepeat["=", 64]];

  shapes = SchoberlZaglmayrBasis`ShapesHEX[order, Range[0, 7], x, y, z];
  n = Length[shapes];
  Print["  Number of basis functions: ", n];

  t0 = AbsoluteTime[];
  M = Table[
    If[j < i,
      0,   (* fill below diagonal later via symmetry *)
      Integrate[shapes[[i, 1]] . shapes[[j, 1]], {x, 0, 1}, {y, 0, 1}, {z, 0, 1}]
    ],
    {i, 1, n}, {j, 1, n}
  ];
  Print["  Computed upper triangle in ", Round[(AbsoluteTime[] - t0), 0.1], " seconds"];

  M = M + Transpose[M] - DiagonalMatrix[Diagonal[M]];
  Return[M];
];

(* === HEX Stiffness (curl-curl) matrix === *)
HexStiffnessMatrix[order_] := Module[{shapes, n, S, t0},
  Print[StringRepeat["-", 64]];
  Print["HEX H(curl) curl-curl Stiffness matrix: order = ", order];
  Print[StringRepeat["-", 64]];

  shapes = SchoberlZaglmayrBasis`ShapesHEX[order, Range[0, 7], x, y, z];
  n = Length[shapes];
  Print["  Number of basis functions: ", n];

  t0 = AbsoluteTime[];
  S = Table[
    If[j < i, 0,
      Integrate[shapes[[i, 2]] . shapes[[j, 2]], {x, 0, 1}, {y, 0, 1}, {z, 0, 1}]
    ],
    {i, 1, n}, {j, 1, n}
  ];
  Print["  Computed upper triangle in ", Round[(AbsoluteTime[] - t0), 0.1], " seconds"];
  S = S + Transpose[S] - DiagonalMatrix[Diagonal[S]];
  Return[S];
];


(* === Run order 0 first (fastest) === *)
ord = 0;
allowOrder0 = False;
M0 = Quiet[Check[
  HexMassMatrix[ord],
  Module[{shapes, n, M},
    Print["  Order 0 not allowed by ShapesHEX constraint; trying order 1..."];
    allowOrder0 = False;
    Null
  ]]];

If[M0 === Null,
  Print[];
  ord = 1;
  M1 = HexMassMatrix[1];
  S1 = HexStiffnessMatrix[1];
  Print[];
  Print["===== Mass matrix structure (order 1, 54x54) ====="];
  n = Length[M1];
  Print["  Diagonal entries:"];
  Do[Print["    M[", i, "][", i, "] = ", M1[[i, i]]], {i, 1, Min[12, n]}];
  Print["  Number of nonzero entries (upper tri): ",
        Count[Flatten[Table[M1[[i, j]], {i, 1, n}, {j, i, n}]], _?(# =!= 0 &)]];
  Print["  Number of total entries (upper tri): ", n (n + 1)/2];

  Print[];
  Print["===== Stiffness matrix structure (order 1) ====="];
  Print["  Diagonal entries:"];
  Do[Print["    S[", i, "][", i, "] = ", S1[[i, i]]], {i, 1, Min[12, n]}];
  Print["  Number of nonzero entries (upper tri): ",
        Count[Flatten[Table[S1[[i, j]], {i, 1, n}, {j, i, n}]], _?(# =!= 0 &)]];
];
