(* Derive cuboid Newton potential by symbolic triple integration.
   Use the result as ground truth, validate corner-sum formula against it. *)
$HistoryLength = 0;

ClearAll[xpp, ypp, zpp, xv, yv, zv, cubA, cubB, cubC, uu, vv, ww];

Print["=== Cuboid Newton potential: symbolic derivation ==="];
Print[""];

(* Test simple case: a=b=c=1, field point at origin *)
Print["Test 1: cube [-1,1]^3 at field point (0,0,0)"];
Print["  N = integral 1/sqrt(x'^2+y'^2+z'^2) over cube"];
sym1 = Integrate[1/Sqrt[xpp^2 + ypp^2 + zpp^2], {xpp, -1, 1}];
Print["  Inner xpp at (yp, zp): ", sym1];
sym2 = Integrate[sym1, {ypp, -1, 1}];
Print["  After ypp: ", Short[sym2, 3]];
sym3 = Integrate[sym2, {zpp, -1, 1}];
Print["  Final cube Newton at origin: ", sym3];
Print["    Numeric: ", N[sym3, 16]];
Print[""];

(* Compare to direct NIntegrate *)
ndi = NIntegrate[1/Sqrt[xpp^2 + ypp^2 + zpp^2],
  {xpp, -1, 1}, {ypp, -1, 1}, {zpp, -1, 1},
  WorkingPrecision -> 30, AccuracyGoal -> 22, MaxRecursion -> 12];
Print["  NIntegrate: ", N[ndi, 16]];
Print[""];

(* Now derive corner formula by indefinite triple integration *)
Print["=== Indefinite triple integral (anti-derivative) ==="];
Print["  G(u,v,w) such that d^3 G / du dv dw = 1/sqrt(u^2+v^2+w^2)"];
Print[""];

(* Step 1: integrate 1/r wrt u *)
g1 = Integrate[1/Sqrt[uu^2 + vv^2 + ww^2], uu,
  Assumptions -> {uu \[Element] Reals, vv \[Element] Reals, ww \[Element] Reals}];
Print["  After u-integration: ", g1];
Print[""];

(* Step 2: integrate wrt v *)
g2 = Integrate[g1, vv,
  Assumptions -> {uu \[Element] Reals, vv \[Element] Reals, ww \[Element] Reals}];
Print["  After v-integration: ", g2];
Print[""];

(* Step 3: integrate wrt w *)
g3 = Integrate[g2, ww,
  Assumptions -> {uu \[Element] Reals, vv \[Element] Reals, ww \[Element] Reals}];
Print["  After w-integration (G(u,v,w)): ", g3];
Print["  LeafCount = ", LeafCount[g3]];
Print[""];

(* Now use this g3 to compute N at field point (xv, yv, zv) for cuboid [-cubA, cubA]^3.
   N = G(cubA - xv, cubB - yv, cubC - zv) - G(-cubA - xv, ...) etc.
   = sum over corners of G evaluated at (corner - field point), with sign factor. *)

cornerSum[xv0_, yv0_, zv0_, ca_, cb_, cc_] :=
  Sum[Sign * (g3 /. {uu -> uVal, vv -> vVal, ww -> wVal}),
    {Sign, uVal, vVal, wVal},
    {{1, ca - xv0}, {-1, -ca - xv0}},
    {{1, cb - yv0}, {-1, -cb - yv0}},
    {{1, cc - zv0}, {-1, -cc - zv0}}];

Print["=== Test corner sum at (0,0,0) for cube [-1,1]^3 ==="];
test1 = cornerSum[0, 0, 0, 1, 1, 1];
Print["  Corner sum result: ", test1];
Print["  Numeric: ", N[test1, 16]];
Print["  vs NIntegrate: ", N[ndi, 16]];
Print["  rel diff: ", N[Abs[(test1 - ndi)/ndi], 6]];
