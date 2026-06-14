(* Phase A v5: use Mathematica-derived g3 antiderivative. *)
$HistoryLength = 0;

ClearAll[uu, vv, ww];

Print["=== Phase A v5: Mathematica-derived antiderivative ==="];
Print[""];

(* Derive g3 = triple antiderivative of 1/sqrt(u^2+v^2+w^2) *)
g1 = Integrate[1/Sqrt[uu^2 + vv^2 + ww^2], uu,
  Assumptions -> {uu \[Element] Reals, vv \[Element] Reals, ww \[Element] Reals}];
g2 = Integrate[g1, vv,
  Assumptions -> {uu \[Element] Reals, vv \[Element] Reals, ww \[Element] Reals}];
g3 = Integrate[g2, ww,
  Assumptions -> {uu \[Element] Reals, vv \[Element] Reals, ww \[Element] Reals}];
Print["g3 LeafCount = ", LeafCount[g3]];
Print[""];

(* g3 evaluation function: substitute numeric values and take real part *)
gFunc[uval_?NumericQ, vval_?NumericQ, wval_?NumericQ] := Re[
  N[g3 /. {uu -> SetPrecision[uval, 50], vv -> SetPrecision[vval, 50],
            ww -> SetPrecision[wval, 50]}, 30]];

(* Newton potential by 8-corner sum. cuboid centered at origin: [-A, A]x[-B, B]x[-C, C] *)
newtonPot[xv_?NumericQ, yv_?NumericQ, zv_?NumericQ,
          aHalf_?NumericQ, bHalf_?NumericQ, cHalf_?NumericQ] := Block[
  {result = 0, sx, sy, sz},
  Do[Do[Do[
    result += sx sy sz gFunc[sx aHalf - xv, sy bHalf - yv, sz cHalf - zv],
    {sz, {1, -1}}], {sy, {1, -1}}], {sx, {1, -1}}];
  result];

(* Test on unit cube *)
Print["Test 1: cube [-1,1]^3 at origin"];
testN1 = newtonPot[0., 0., 0., 1., 1., 1.];
ndi1 = NIntegrate[1/Sqrt[xpp^2 + ypp^2 + zpp^2],
  {xpp, -1, 1}, {ypp, -1, 1}, {zpp, -1, 1},
  WorkingPrecision -> 30, AccuracyGoal -> 22, MaxRecursion -> 12];
Print["  corner sum  = ", N[testN1, 16]];
Print["  NIntegrate  = ", N[ndi1, 16]];
Print["  rel diff    = ", N[Abs[(testN1 - ndi1)/ndi1], 6]];
Print[""];

(* Test 2: our actual cuboid 5x2x1 mm *)
Print["Test 2: cuboid 5x2x1 mm at various field points"];
prec = 50;
aN = SetPrecision[5/1000, prec];
bN = SetPrecision[2/1000, prec];
cN = SetPrecision[1/1000, prec];
testPts = {{aN/4, bN/4, cN/4}, {aN/3, bN/5, 0}, {0, 0, 0},
           {aN/2 - cN/10, 0, 0}, {0, bN/3, 0}};
Do[Module[{r = pt, ncl, ndir, t1, t2, x0, y0, z0},
  {x0, y0, z0} = r;
  t1 = AbsoluteTime[];
  ncl = newtonPot[x0, y0, z0, aN/2, bN/2, cN/2];
  t1 = AbsoluteTime[] - t1;
  t2 = AbsoluteTime[];
  ndir = NIntegrate[1/Sqrt[(x0 - xpp)^2 + (y0 - ypp)^2 + (z0 - zpp)^2],
    {xpp, -aN/2, aN/2}, {ypp, -bN/2, bN/2}, {zpp, -cN/2, cN/2},
    WorkingPrecision -> 30, AccuracyGoal -> 22, MaxRecursion -> 12];
  t2 = AbsoluteTime[] - t2;
  Print["  r = ", N[r*1000, 6], " mm"];
  Print["    closed form  = ", N[ncl, 16], "  (", N[t1, 4], "s)"];
  Print["    direct 3D NI = ", N[ndir, 16], "  (", N[t2, 4], "s)"];
  Print["    rel diff     = ", N[Abs[(ncl - ndir)/ndir], 6]]],
  {pt, testPts}];
Print[""];

(* === If verified, build B_x via IBP === *)
Print["=== Test B_x via IBP: B_x = xv*N(r) + Fp(r) - Fm(r) ==="];

facePlus[xv_?NumericQ, yv_?NumericQ, zv_?NumericQ] := NIntegrate[
  Sqrt[(xv - aN/2)^2 + (yv - ypp)^2 + (zv - zpp)^2],
  {ypp, -bN/2, bN/2}, {zpp, -cN/2, cN/2},
  WorkingPrecision -> 30, AccuracyGoal -> 22, MaxRecursion -> 10];

faceMinus[xv_?NumericQ, yv_?NumericQ, zv_?NumericQ] := NIntegrate[
  Sqrt[(xv + aN/2)^2 + (yv - ypp)^2 + (zv - zpp)^2],
  {ypp, -bN/2, bN/2}, {zpp, -cN/2, cN/2},
  WorkingPrecision -> 30, AccuracyGoal -> 22, MaxRecursion -> 10];

bxIBP[xv_?NumericQ, yv_?NumericQ, zv_?NumericQ] :=
  xv * newtonPot[xv, yv, zv, aN/2, bN/2, cN/2]
    + facePlus[xv, yv, zv] - faceMinus[xv, yv, zv];

Do[Module[{r = pt, fcl, fdir, t1, t2, x0, y0, z0},
  {x0, y0, z0} = r;
  t1 = AbsoluteTime[];
  fcl = bxIBP[x0, y0, z0];
  t1 = AbsoluteTime[] - t1;
  t2 = AbsoluteTime[];
  fdir = NIntegrate[
    xpp/Sqrt[(x0 - xpp)^2 + (y0 - ypp)^2 + (z0 - zpp)^2],
    {xpp, -aN/2, aN/2}, {ypp, -bN/2, bN/2}, {zpp, -cN/2, cN/2},
    WorkingPrecision -> 30, AccuracyGoal -> 22, MaxRecursion -> 12];
  t2 = AbsoluteTime[] - t2;
  Print["  r = ", N[r*1000, 6], " mm"];
  Print["    via IBP+face = ", N[fcl, 16], "  (", N[t1, 4], "s)"];
  Print["    direct 3D NI = ", N[fdir, 16], "  (", N[t2, 4], "s)"];
  Print["    rel diff     = ", N[Abs[(fcl - fdir)/fdir], 6]]],
  {pt, testPts}];
