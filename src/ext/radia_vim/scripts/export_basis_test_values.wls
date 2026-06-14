#!/usr/bin/env wolframscript
(*
  export_basis_test_values.wls — dump reference values for the HDiv div-free
  hex basis at a few test points, for use as ground truth in
  tests/test_basis_vs_mathematica.py.

  The reference implementation is the canonical Mathematica
  ShapesHDivHexDivFree (from hex_bem_hdiv_order3_overnight.wls), reproduced
  here verbatim so the test does not depend on .mx files.

  Output: tests/basis_test_values.json
*)

ClearAll[x, y, z];

(* === Reference impl (verbatim from hex_bem_hdiv_order3_overnight.wls) === *)

DuCrossDv[u_, v_, vars_] := Cross[D[u, {vars}], D[v, {vars}]];
wDuCrossDv[u_, v_, w_, vars_] := w Cross[D[u, {vars}], D[v, {vars}]];

IntegratedLegendre[n_, xv_] := Module[{p3, p2, p1, j, result},
  result = {};
  If[n < 2, Return[{}]];
  p3 = 0; p2 = -1; p1 = xv;
  Do[
    p3 = p2; p2 = p1;
    p1 = (2 j - 3)/j * xv * p2 - (j - 3)/j * p3;
    AppendTo[result, p1],
    {j, 2, n}
  ];
  result
];

LegendrePoly[n_, xv_] := Module[{P, k},
  P = Table[0, {n + 1}];
  P[[1]] = 1;
  If[n >= 1, P[[2]] = xv];
  Do[
    P[[k + 1]] = (2 - 1/k) xv P[[k]] + (1/k - 1) P[[k - 1]],
    {k, 2, n}
  ];
  P
];

ShapesHDivHexDivFree[order_, xv_, yv_, zv_] := Module[
  {p, xi, eta, zeta, LXi, LEta, LZeta, PXi, PEta, PZeta, shapes, vars},
  vars = {xv, yv, zv};
  p = order;
  xi = 2 xv - 1; eta = 2 yv - 1; zeta = 2 zv - 1;
  LXi   = IntegratedLegendre[p + 1, xi];
  LEta  = IntegratedLegendre[p + 1, eta];
  LZeta = IntegratedLegendre[p + 1, zeta];
  PXi   = LegendrePoly[p, xi];
  PEta  = LegendrePoly[p, eta];
  PZeta = LegendrePoly[p, zeta];
  shapes = {};
  Do[AppendTo[shapes, wDuCrossDv[LXi[[i+1]] LZeta[[k+1]], yv, -4 PEta[[j+2]], vars]],
    {i, 0, p-1}, {j, 0, p-1}, {k, 0, p-1}];
  Do[AppendTo[shapes, wDuCrossDv[LEta[[j+1]] LZeta[[k+1]], xv, 4 PXi[[i+2]], vars]],
    {i, 0, p-1}, {j, 0, p-1}, {k, 0, p-1}];
  Do[AppendTo[shapes, DuCrossDv[LEta[[j+1]] LZeta[[k+1]], 2 xi, vars]],
    {j, 0, p-1}, {k, 0, p-1}];
  Do[AppendTo[shapes, DuCrossDv[LXi[[i+1]] LZeta[[k+1]], 2 eta, vars]],
    {i, 0, p-1}, {k, 0, p-1}];
  Do[AppendTo[shapes, DuCrossDv[LXi[[i+1]] LEta[[j+1]], 2 zeta, vars]],
    {i, 0, p-1}, {j, 0, p-1}];
  shapes
];

(* === Test points === *)
testPoints = {
  {0.5,    0.5,    0.5},    (* center, mostly zero by symmetry *)
  {0.3,    0.7,    0.2},    (* asymmetric *)
  {0.1,    0.9,    0.8},
  {0.999,  0.001,  0.5},    (* near boundary *)
  {0.25,   0.5,    0.75}
};

(* === Run for orders 1, 2, 3 === *)
allResults = Table[
  Module[{shapes = ShapesHDivHexDivFree[ord, x, y, z], values, ndof},
    ndof = Length[shapes];
    values = Table[
      Module[{px = pt[[1]], py = pt[[2]], pz = pt[[3]], rule},
        rule = {x -> px, y -> py, z -> pz};
        N[shapes /. rule, 16]
      ],
      {pt, testPoints}
    ];
    <|
      "order" -> ord,
      "n_dofs" -> ndof,
      "values" -> values
    |>
  ],
  {ord, 1, 3}
];

result = <|
  "test_points" -> testPoints,
  "by_order" -> allResults
|>;

scriptDir = DirectoryName[$InputFileName];
outPath = FileNameJoin[{scriptDir, "..", "tests", "basis_test_values.json"}];
Export[outPath, result, "JSON"];

Print["Saved: ", outPath];
Print["  orders 1..3, ", Length[testPoints], " test points each"];
Print["  order 3 ndof = ", allResults[[3]]["n_dofs"]];
