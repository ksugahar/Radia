#!/usr/bin/env wolframscript
(*
  hex_bem_hdiv_order3_export_json.wls — read the overnight .mx and export
  eigenvalues + Foster amplitudes to JSON for downstream Python Cauer
  extraction at arbitrary moment count + precision.

  Mirrors `bem_disk_axisym_cauer.json` schema so that the same Python
  Cauer-extraction script structure can run on cuboid data.
*)

ClearAll[x, y, z];
$MaxExtraPrecision = 60;

Print["Loading hex_bem_hdiv_order3_results.mx ..."];
Get["W:\\30_CauerLadderNetwork\\2026_04_01_\:9577\:65b9\:5f62CLN\\ngsolve_validation\\hex_bem_hdiv_order3_results.mx"];
Print["  nDOF = ", nDOF];
Print["  eigvals length = ", Length[eigvals]];
Print["  gVec   length = ", Length[gVec]];
Print["  sortedPhysical (top 5) = ", sortedPhysical[[1;;5]]];

(* Constants used in this benchmark *)
sigma = 58 10^6;
mu0 = 4 Pi 10^-7;
{aBox, bBox, cBox} = {5 10^-3, 2 10^-3, 1 10^-3};

(* Use sortedPhysical: g_k^2 > g_max/10^6 AND tau_k > 0, sorted by tau desc.
   This filters out the symmetry-forbidden modes (g_k = 0) that would
   otherwise contribute spurious alpha_n moments at arbitrary precision. *)
keepIdx = sortedPhysical;

eigvalsSorted = N[Re[eigvals[[keepIdx]]], 30];
gVecSorted    = N[Re[gVec[[keepIdx]]],    30];
tausUsSorted  = 10^6 * sigma * eigvalsSorted;

Print["  physical modes kept (sortedPhysical) = ", Length[keepIdx]];
Print["  tau_1..tau_6 us    = ", N[tausUsSorted[[1;;Min[6, Length[tausUsSorted]]]], 8]];

result = <|
  "method" -> "BEM HDiv order 3 div-free interior basis (81 DOFs), J.n=0 isolated cuboid",
  "params" -> <|"a_mm" -> aBox 1000, "b_mm" -> bBox 1000, "c_mm" -> cBox 1000,
                "sigma" -> N[sigma], "mu0" -> N[mu0],
                "nDOF" -> nDOF, "n_modes" -> Length[keepIdx]|>,
  "bem_tau_us"          -> tausUsSorted,
  "bem_eigvals_lambda"  -> eigvalsSorted,
  "foster_amplitudes_g" -> gVecSorted,
  "leading_tau_us"      -> tausUsSorted[[1]]
|>;

Export["hex_bem_hdiv_order3_export.json", result, "JSON"];
Print["Saved: hex_bem_hdiv_order3_export.json (", Length[keepIdx], " modes)"];
