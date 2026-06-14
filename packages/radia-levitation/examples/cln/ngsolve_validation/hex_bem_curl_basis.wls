(* ::Package:: *)

(* hex_bem_curl_basis.wls
   BEM eigenvalue with H(curl) basis + Newton kernel:
     K[i,j] = (mu0/4pi) int_cell int_cell phi_i(r) . phi_j(r') / |r-r'| dV dV'
     M[i,j] = int_cell phi_i . phi_j dV
     Generalized eigenproblem: K v = lambda M v -> tau = lambda * sigma
   This corresponds to Phase 7-11's BEM (vacuum-coupled) eddy current operator.

   Strategy: start with HCurl order 1 LOWEST 12 EDGE basis only (sub-block).
   Each K entry = 6D NIntegrate at moderate WP. Estimate scaling, then expand.

   Date: 2026-05-02
*)

Get["w:\\30_CauerLadderNetwork\\2026_04_01_\:9577\:65b9\:5f62CLN\\ngsolve_validation\\schoberl_zaglmayr_basis.wls"];

ClearAll[x, y, z, xp, yp, zp];
$MaxExtraPrecision = 50;

(* === Physical cuboid === *)
sigma = SetPrecision[58 10^6, 30];
mu0   = 4 Pi SetPrecision[10^-7, 30];
{aBox, bBox, cBox} = SetPrecision[{5 10^-3, 2 10^-3, 1 10^-3}, 30];

(* === Reference HCurl basis on [0,1]^3 === *)
ord = 1;
shapes = SchoberlZaglmayrBasis`ShapesHEX[ord, Range[0, 7], x, y, z];
nFull = Length[shapes];
Print["HCurl order 1: ", nFull, " total basis"];

(* Take only LOWEST-ORDER edge functions (indices 1..12, before high-order edge added).
   In my current ordering, EVERY edge has 1 lowest + p high-order, in pairs.
   For order 1: indices 1, 3, 5, 7, ..., 23 are lowest; indices 2, 4, ..., 24 are high.
   But to be safe, explicitly construct the lowest-only basis at order 0 by
   manually constructing the lowest 12 from C++ formula:
*)
lowestEdges = Module[{lami, sigmaSym, edges, e, xi, lamE},
  lami = {(1-x)(1-y)(1-z), x(1-y)(1-z), x*y*(1-z), (1-x)*y*(1-z),
          (1-x)(1-y)*z,    x(1-y)*z,    x*y*z,    (1-x)*y*z};
  sigmaSym = {(1-x)+(1-y)+(1-z), x+(1-y)+(1-z), x+y+(1-z), (1-x)+y+(1-z),
              (1-x)+(1-y)+z,     x+(1-y)+z,     x+y+z,     (1-x)+y+z};
  Table[
    e = SchoberlZaglmayrBasis`HexEdges[[i+1]];
    xi = sigmaSym[[e[[2]]+1]] - sigmaSym[[e[[1]]+1]];
    lamE = lami[[e[[1]]+1]] + lami[[e[[2]]+1]];
    SchoberlZaglmayrBasis`uDvOp[0.5 lamE, xi, {x, y, z}],
    {i, 0, 11}
  ]
];

Print["Constructed ", Length[lowestEdges], " lowest-order edge basis functions"];

(* Apply Piola covariant scaling: phi_phys_alpha = phi_ref_alpha / d_alpha
   where d_alpha = aBox, bBox, cBox for alpha = x, y, z. *)
piolaShape[shapeRef_, a_, b_, c_] := Module[{val = shapeRef[[1]]},
  {{val[[1]]/a, val[[2]]/b, val[[3]]/c}, shapeRef[[2]]}   (* curl unchanged for now *)
];

physShapes = piolaShape[#, aBox, bBox, cBox] & /@ lowestEdges;

(* === Mass matrix M[i,j] = int phi_i . phi_j dV (in physical coords) === *)
fastIntegrate3D[poly_, {a_, b_, c_}] := Module[{p = Expand[poly], rules},
  If[p === 0, Return[0]];
  rules = CoefficientRules[p, {x, y, z}];
  (* Integrate over [0, a/aBox] x [0, b/bBox] x [0, c/cBox] = [0, 1]^3 in ref *)
  (* Map to physical coords: x_phys = a x_ref, dV_phys = a*b*c dV_ref *)
  (* For polynomial in (x_ref, y_ref, z_ref), int over [0,1]^3 then multiply by a*b*c *)
  a b c Total[Function[r, r[[2]] / ((r[[1,1]]+1)(r[[1,2]]+1)(r[[1,3]]+1))] /@ rules]
];

n = Length[physShapes];
Print["Computing mass matrix M (", n, "x", n, ") in physical coords ..."];
t0 = AbsoluteTime[];
M = Table[
  fastIntegrate3D[
    physShapes[[i, 1]] . physShapes[[j, 1]], {aBox, bBox, cBox}
  ],
  {i, 1, n}, {j, 1, n}
];
Print["  M done in ", Round[AbsoluteTime[] - t0, 0.1], " s"];

Print["M diagonal (first 6):"];
Do[Print["  M[", i, "][", i, "] = ", N[M[[i, i]], 12]], {i, 1, 6}];

(* === BEM K matrix entry via 6D NIntegrate ===
   Single entry test first to estimate timing.
*)
Kentry[shapeI_, shapeJ_] := Module[{vi, vj, integrand},
  vi = shapeI[[1]];
  vj = shapeJ[[1]] /. {x -> xp, y -> yp, z -> zp};
  integrand = (vi . vj) / Sqrt[(x-xp)^2 + (y-yp)^2 + (z-zp)^2];
  NIntegrate[integrand,
    {x, 0, aBox}, {y, 0, bBox}, {z, 0, cBox},
    {xp, 0, aBox}, {yp, 0, bBox}, {zp, 0, cBox},
    Method -> {"GlobalAdaptive", "MaxErrorIncreases" -> 100000,
               "SymbolicProcessing" -> 0},
    WorkingPrecision -> 25, PrecisionGoal -> 6, AccuracyGoal -> 8,
    MaxRecursion -> 20
  ] (mu0 / (4 Pi))
];

Print[];
Print["===== Single entry test: K[1,1] (lowest-order edge 0-1, self-coupling) ====="];
t0 = AbsoluteTime[];
K11 = Kentry[physShapes[[1]], physShapes[[1]]];
Print["  K[1,1] = ", N[K11, 8], "  (", Round[AbsoluteTime[] - t0, 0.1], " s)"];

(* If this completes in reasonable time, do a few more entries *)
Print[];
Print["===== Sample 4 diagonal entries at higher WP ====="];
Do[
  t0 = AbsoluteTime[];
  Module[{Kii = Kentry[physShapes[[i]], physShapes[[i]]]},
    Print["  K[", i, "][", i, "] = ", N[Kii, 10], "  (", Round[AbsoluteTime[] - t0, 0.1], " s)"]
  ],
  {i, {1, 3, 5, 9}}
];
