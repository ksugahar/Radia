#!/usr/bin/env wolframscript
(*
  extract_K_reference.wls — pull selected K matrix entries from the
  overnight .mx file and dump them to JSON for the Phase F-3 quadrature
  validation.

  Output: tests/K_reference_values.json

  The K matrix in the .mx file already includes the (mu0/(4 Pi)) prefactor
  (see hex_bem_hdiv_order3_overnight.wls line 125: KentryRobust returns
  (mu0/(4 Pi)) * NIntegrate[...]). For Phase F-3 comparison we strip the
  prefactor by dividing by mu0/(4 Pi) = 1.0e-7.
*)

ClearAll[x, y, z];

Print["Loading hex_bem_hdiv_order3_results.mx ..."];
Get["W:\\30_CauerLadderNetwork\\2026_04_01_\:9577\:65b9\:5f62CLN\\ngsolve_validation\\hex_bem_hdiv_order3_results.mx"];
Print["  K dimensions: ", Dimensions[K]];

mu0Over4Pi = 4 Pi * 10^-7 / (4 Pi);   (* = 10^-7 *)

(* Pick a spread of indices: leading + middle + tail; physical (max g²) mode is rank 4 of sortedPhysical *)
testIdx = {{1, 1}, {1, 2}, {2, 2}, {1, 8}, {8, 8}, {7, 7}, {41, 41}};

values = Table[
  Module[{ii = pair[[1]], jj = pair[[2]], v},
    v = N[Re[K[[ii, jj]]] / mu0Over4Pi, 16];
    <|
      "i" -> ii - 1,    (* convert to 0-based for Python *)
      "j" -> jj - 1,
      "K_ij_bare" -> v,
      "K_ij_with_prefactor" -> N[Re[K[[ii, jj]]], 16]
    |>
  ],
  {pair, testIdx}
];

result = <|
  "method" -> "Mathematica NIntegrate WP=15 (overnight, sym positive entries)",
  "geometry" -> <|"a_mm" -> 5, "b_mm" -> 2, "c_mm" -> 1|>,
  "order" -> 3,
  "n_dofs" -> 81,
  "mu0_over_4pi" -> N[mu0Over4Pi],
  "entries" -> values
|>;

scriptDir = DirectoryName[$InputFileName];
outPath = FileNameJoin[{scriptDir, "..", "tests", "K_reference_values.json"}];
Export[outPath, result, "JSON"];
Print["Saved: ", outPath];
Do[
  Print["  K[", e["i"], ",", e["j"], "] bare = ", e["K_ij_bare"]],
  {e, values}
];
