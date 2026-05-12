(* Phase A v2: B_x(r) closed form via IBP decomposition.
   Fix: avoid underscores in variable names (Mathematica pattern conflict).

   B_x(r) = x_R N(r) + Fp(r) - Fm(r)
   where:
     N(r) = integral_Omega 1/|r-rpp| dVpp                (cuboid Newton potential)
     Fp(r) = integral_face(xpp=a/2)  |r-(a/2,ypp,zpp)|  dypp dzpp
     Fm(r) = integral_face(xpp=-a/2) |r-(-a/2,ypp,zpp)| dypp dzpp
*)
$HistoryLength = 0;
prec = 60;

Print["=== Phase A v2: B_x closed form (camelCase fix) ==="];
Print[""];

ClearAll[xR, yR, zR, xpp, ypp, zpp, asym, bsym, csym, Rsq];

(* === A1: Newton potential N(r) === *)
Print["A1: Newton potential newton(r)"];
Print["  A1a: xpp integral"];
newtonXpp = TimeConstrained[
  Integrate[1/Sqrt[(xR - xpp)^2 + Rsq], {xpp, -asym/2, asym/2},
    Assumptions -> {asym > 0, Rsq > 0, xR \[Element] Reals}],
  60, $Failed];
If[newtonXpp === $Failed, Print["    TIMEOUT"]; Quit[]];
Print["    LeafCount = ", LeafCount[newtonXpp]];

newtonXppSubst = newtonXpp /. Rsq -> (yR - ypp)^2 + (zR - zpp)^2;
Print["  A1b: ypp integral (timeout 600s)"];
newtonYpp = TimeConstrained[
  Integrate[newtonXppSubst, {ypp, -bsym/2, bsym/2},
    Assumptions -> {asym > 0, bsym > 0, csym > 0,
      xR \[Element] Reals, yR \[Element] Reals, zR \[Element] Reals,
      zpp \[Element] Reals}],
  600, $Failed];
If[newtonYpp === $Failed, Print["    TIMEOUT"]; Quit[]];
Print["    LeafCount = ", LeafCount[newtonYpp]];

Print["  A1c: zpp integral (timeout 1200s)"];
newtonFull = TimeConstrained[
  Integrate[newtonYpp, {zpp, -csym/2, csym/2},
    Assumptions -> {asym > 0, bsym > 0, csym > 0,
      xR \[Element] Reals, yR \[Element] Reals, zR \[Element] Reals}],
  1200, $Failed];
If[newtonFull === $Failed, Print["    TIMEOUT"]; Quit[]];
Print["    LeafCount = ", LeafCount[newtonFull]];
Print["  newton(r) closed form FOUND"];
Print[""];

(* === A2: Face integrals === *)
Print["A2: Face integrals fPlus, fMinus"];
Print["  A2a: fPlus inner ypp"];
fPlusYpp = TimeConstrained[
  Integrate[Sqrt[(xR - asym/2)^2 + (yR - ypp)^2 + (zR - zpp)^2],
    {ypp, -bsym/2, bsym/2},
    Assumptions -> {asym > 0, bsym > 0, csym > 0,
      xR \[Element] Reals, yR \[Element] Reals, zR \[Element] Reals,
      zpp \[Element] Reals}],
  300, $Failed];
If[fPlusYpp === $Failed, Print["    TIMEOUT"]; Quit[]];
Print["    LeafCount = ", LeafCount[fPlusYpp]];

Print["  A2b: fPlus zpp"];
fPlus = TimeConstrained[
  Integrate[fPlusYpp, {zpp, -csym/2, csym/2},
    Assumptions -> {asym > 0, bsym > 0, csym > 0,
      xR \[Element] Reals, yR \[Element] Reals, zR \[Element] Reals}],
  600, $Failed];
If[fPlus === $Failed, Print["    TIMEOUT"]; Quit[]];
Print["    fPlus FOUND, LeafCount = ", LeafCount[fPlus]];

(* fMinus: replace asym/2 -> -asym/2 *)
fMinus = fPlus /. {(xR - asym/2) -> (xR + asym/2)};
Print["    fMinus LeafCount = ", LeafCount[fMinus]];
Print[""];

(* === A3: Combine === *)
Print["A3: bxClosed = xR newton + fPlus - fMinus"];
bxClosed = xR newtonFull + fPlus - fMinus;
Print["  bxClosed (raw) LeafCount = ", LeafCount[bxClosed]];

bxClosedSimp = TimeConstrained[
  Simplify[bxClosed,
    Assumptions -> {asym > 0, bsym > 0, csym > 0,
      xR \[Element] Reals, yR \[Element] Reals, zR \[Element] Reals}],
  300, bxClosed];
Print["  bxClosed (simplified) LeafCount = ", LeafCount[bxClosedSimp]];
Print[""];

(* === A4: Numerical verification === *)
Print["A4: numerical verification"];
sigmaN = SetPrecision[58000000, prec];
mu0N = SetPrecision[4 Pi/10^7, prec];
aN = SetPrecision[5/1000, prec];
bN = SetPrecision[2/1000, prec];
cN = SetPrecision[1/1000, prec];
PiN = SetPrecision[Pi, prec];

bxFunc[xRv_?NumberQ, yRv_?NumberQ, zRv_?NumberQ] := Module[{result},
  result = bxClosed /. {asym -> aN, bsym -> bN, csym -> cN,
    xR -> SetPrecision[xRv, prec],
    yR -> SetPrecision[yRv, prec],
    zR -> SetPrecision[zRv, prec]};
  Re[N[result, 30]]];

bxDirect[xRv_?NumberQ, yRv_?NumberQ, zRv_?NumberQ] := NIntegrate[
  xpp/Sqrt[(xRv - xpp)^2 + (yRv - ypp)^2 + (zRv - zpp)^2],
  {xpp, -aN/2, aN/2}, {ypp, -bN/2, bN/2}, {zpp, -cN/2, cN/2},
  WorkingPrecision -> 25, AccuracyGoal -> 18, PrecisionGoal -> 16,
  MaxRecursion -> 12];

testPts = {{aN/4, bN/4, cN/4}, {aN/3, bN/5, 0}, {aN/2 - aN/100, 0, 0}};
Do[Module[{r, fcl, fdir, t1, t2},
  r = pt;
  t1 = AbsoluteTime[];
  fcl = bxFunc @@ r;
  t1 = AbsoluteTime[] - t1;
  t2 = AbsoluteTime[];
  fdir = bxDirect @@ r;
  t2 = AbsoluteTime[] - t2;
  Print["  r = ", N[r, 6]];
  Print["    closed form  = ", N[fcl, 14], "  (", N[t1, 4], "s)"];
  Print["    direct 3D NI = ", N[fdir, 14], "  (", N[t2, 4], "s)"];
  Print["    rel diff = ", N[Abs[(fcl - fdir)/fdir], 6]]],
  {pt, testPts}];
Print[""];

Print["=== Saving ==="];
DumpSave["bx_closed_form.mx", {bxClosed, bxClosedSimp, newtonFull, fPlus, fMinus}];
Print["  Saved to bx_closed_form.mx"];
Print[""];

Print["=== Phase A v2 summary ==="];
Print["* bxClosed LeafCount = ", LeafCount[bxClosed]];
Print["* bxClosedSimp LeafCount = ", LeafCount[bxClosedSimp]];
Print["* Verification ratios above"];
