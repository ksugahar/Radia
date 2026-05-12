(* Phase A v6: ComplexExpand to clean up complex logs.

   v5 failed because Re[g3 evaluated at numeric values] gives wrong result
   due to independent branch-cut issues in complex log terms.
   Solution: ComplexExpand g3 to combine complex logs into atan/log of real args. *)
$HistoryLength = 0;

ClearAll[uu, vv, ww];

Print["=== Phase A v6: ComplexExpand antiderivative ==="];
Print[""];

(* Derive g3 *)
g1 = Integrate[1/Sqrt[uu^2 + vv^2 + ww^2], uu,
  Assumptions -> {uu \[Element] Reals, vv \[Element] Reals, ww \[Element] Reals}];
g2 = Integrate[g1, vv,
  Assumptions -> {uu \[Element] Reals, vv \[Element] Reals, ww \[Element] Reals}];
g3 = Integrate[g2, ww,
  Assumptions -> {uu \[Element] Reals, vv \[Element] Reals, ww \[Element] Reals}];
Print["g3 raw LeafCount = ", LeafCount[g3]];

(* ComplexExpand to handle complex logs *)
g3Real = ComplexExpand[g3, TargetFunctions -> {Re, Im, Log, ArcTan}];
g3Real = Simplify[g3Real,
  Assumptions -> {uu \[Element] Reals, vv \[Element] Reals, ww \[Element] Reals,
                  uu^2 + vv^2 + ww^2 > 0}];
Print["g3 after ComplexExpand+Simplify LeafCount = ", LeafCount[g3Real]];
Print["  ", Short[g3Real, 6]];
Print[""];

(* Check if g3Real is real for real arguments *)
testU = 0.5; testV = 0.7; testW = 0.3;
testVal = g3Real /. {uu -> testU, vv -> testV, ww -> testW};
Print["g3Real at (0.5, 0.7, 0.3) = ", N[testVal, 16]];
Print["  Im part: ", N[Im[testVal], 16]];
Print[""];

(* Define gFunc based on g3Real (no Re[] needed if real) *)
gFunc[uval_?NumericQ, vval_?NumericQ, wval_?NumericQ] := Block[{result},
  result = g3Real /. {uu -> SetPrecision[uval, 50],
                       vv -> SetPrecision[vval, 50],
                       ww -> SetPrecision[wval, 50]};
  N[result, 30]];

newtonPot[xv_?NumericQ, yv_?NumericQ, zv_?NumericQ,
          aHalf_?NumericQ, bHalf_?NumericQ, cHalf_?NumericQ] := Block[
  {result = 0, sx, sy, sz},
  Do[Do[Do[
    result += sx sy sz gFunc[sx aHalf - xv, sy bHalf - yv, sz cHalf - zv],
    {sz, {1, -1}}], {sy, {1, -1}}], {sx, {1, -1}}];
  Re[result]];

(* Test on unit cube *)
Print["Test 1: cube [-1,1]^3 at origin"];
testN1 = newtonPot[0., 0., 0., 1., 1., 1.];
ndi1 = NIntegrate[1/Sqrt[xpp^2 + ypp^2 + zpp^2],
  {xpp, -1, 1}, {ypp, -1, 1}, {zpp, -1, 1},
  WorkingPrecision -> 30, AccuracyGoal -> 22, MaxRecursion -> 12];
Print["  corner sum = ", N[testN1, 16]];
Print["  NIntegrate = ", N[ndi1, 16]];
Print["  rel diff   = ", N[Abs[(testN1 - ndi1)/ndi1], 6]];
Print[""];

(* Test on actual cuboid *)
Print["Test 2: cuboid 5x2x1 mm"];
prec = 50;
aN = SetPrecision[5/1000, prec];
bN = SetPrecision[2/1000, prec];
cN = SetPrecision[1/1000, prec];
testPts = {{aN/4, bN/4, cN/4}, {0, 0, 0}, {aN/3, bN/5, 0}};
Do[Module[{r = pt, ncl, ndir, x0, y0, z0},
  {x0, y0, z0} = r;
  ncl = newtonPot[x0, y0, z0, aN/2, bN/2, cN/2];
  ndir = NIntegrate[1/Sqrt[(x0 - xpp)^2 + (y0 - ypp)^2 + (z0 - zpp)^2],
    {xpp, -aN/2, aN/2}, {ypp, -bN/2, bN/2}, {zpp, -cN/2, cN/2},
    WorkingPrecision -> 30, AccuracyGoal -> 22, MaxRecursion -> 12];
  Print["  r = ", N[r*1000, 6], " mm"];
  Print["    closed form = ", N[ncl, 16]];
  Print["    NIntegrate  = ", N[ndir, 16]];
  Print["    rel diff    = ", N[Abs[(ncl - ndir)/ndir], 6]]],
  {pt, testPts}];
