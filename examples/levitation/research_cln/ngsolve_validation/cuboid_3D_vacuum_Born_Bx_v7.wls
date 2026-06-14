(* Phase A v7: practical hybrid -- 1D closed (xpp) + 2D NIntegrate (ypp, zpp).

   Avoid 3D closed-form branch-cut issues. Use:
     newtonInnerXpp(xR, R) = closed form of integral 1/sqrt((xR-xpp)^2 + R^2) dxpp
     newtonPot(r) = 2D NIntegrate of newtonInnerXpp over (ypp, zpp)
     bxInnerXpp(xR, R) = closed form of integral xpp/sqrt((xR-xpp)^2 + R^2) dxpp
     B_x(r) = 2D NIntegrate of bxInnerXpp over (ypp, zpp)
*)
$HistoryLength = 0;

ClearAll[uu, ypp, zpp, xR, yR, zR, rsq, aH, bH, cH];

Print["=== Phase A v7: 1D closed + 2D NIntegrate ==="];
Print[""];

(* Closed form: integral 1/sqrt((xR-xpp)^2 + rsq) dxpp from -aH to aH *)
(* = arcsinh((aH-xR)/sqrt(rsq)) + arcsinh((aH+xR)/sqrt(rsq))   for xR in [-aH, aH] *)
newtonInnerXpp[xRv_, rsqv_, aHv_] :=
  ArcSinh[(aHv - xRv)/Sqrt[rsqv]] + ArcSinh[(aHv + xRv)/Sqrt[rsqv]];

(* Closed form: integral xpp/sqrt((xR-xpp)^2 + rsq) dxpp from -aH to aH *)
(* Substitute u = xpp - xR: integral (u + xR)/sqrt(u^2 + rsq) du
   = sqrt(u^2 + rsq) + xR * arcsinh(u/sqrt(rsq)), evaluated *)
bxInnerXpp[xRv_, rsqv_, aHv_] := Block[{r1, r2},
  r1 = Sqrt[(aHv - xRv)^2 + rsqv];
  r2 = Sqrt[(aHv + xRv)^2 + rsqv];
  (r1 - r2) + xRv * (ArcSinh[(aHv - xRv)/Sqrt[rsqv]] + ArcSinh[(aHv + xRv)/Sqrt[rsqv]])];

(* Test inner closed forms numerically *)
Print["Test inner xpp closed forms:"];
prec = 50;
aN = SetPrecision[5/1000, prec];
bN = SetPrecision[2/1000, prec];
cN = SetPrecision[1/1000, prec];

testInner = Block[{xR0 = aN/4, yR0 = bN/4, zR0 = cN/4, ypp0 = bN/8, zpp0 = cN/8, rsq0,
                   ncl, nNI, bcl, bNI},
  rsq0 = (yR0 - ypp0)^2 + (zR0 - zpp0)^2;
  ncl = newtonInnerXpp[xR0, rsq0, aN/2];
  nNI = NIntegrate[1/Sqrt[(xR0 - xpp)^2 + rsq0], {xpp, -aN/2, aN/2},
    WorkingPrecision -> 30, AccuracyGoal -> 25, MaxRecursion -> 12];
  bcl = bxInnerXpp[xR0, rsq0, aN/2];
  bNI = NIntegrate[xpp/Sqrt[(xR0 - xpp)^2 + rsq0], {xpp, -aN/2, aN/2},
    WorkingPrecision -> 30, AccuracyGoal -> 25, MaxRecursion -> 12];
  Print["  newtonInnerXpp closed = ", N[ncl, 16]];
  Print["  newtonInnerXpp NI     = ", N[nNI, 16]];
  Print["  rel diff = ", N[Abs[(ncl - nNI)/nNI], 6]];
  Print["  bxInnerXpp closed = ", N[bcl, 16]];
  Print["  bxInnerXpp NI     = ", N[bNI, 16]];
  Print["  rel diff = ", N[Abs[(bcl - bNI)/bNI], 6]];
  ];
Print[""];

(* Define newton potential and B_x via 2D NIntegrate *)
newtonPot2D[xRv_?NumericQ, yRv_?NumericQ, zRv_?NumericQ] :=
  NIntegrate[newtonInnerXpp[xRv, (yRv - ypp)^2 + (zRv - zpp)^2, aN/2],
    {ypp, -bN/2, bN/2}, {zpp, -cN/2, cN/2},
    WorkingPrecision -> 30, AccuracyGoal -> 22, PrecisionGoal -> 20,
    MaxRecursion -> 12];

bxPot2D[xRv_?NumericQ, yRv_?NumericQ, zRv_?NumericQ] :=
  NIntegrate[bxInnerXpp[xRv, (yRv - ypp)^2 + (zRv - zpp)^2, aN/2],
    {ypp, -bN/2, bN/2}, {zpp, -cN/2, cN/2},
    WorkingPrecision -> 30, AccuracyGoal -> 22, PrecisionGoal -> 20,
    MaxRecursion -> 12];

(* Test against direct 3D NIntegrate *)
Print["Test newtonPot2D and bxPot2D:"];
testPts = {{aN/4, bN/4, cN/4}, {0, 0, 0}, {aN/3, bN/5, 0},
           {aN/2 - aN/100, 0, 0}};
Do[Module[{r = pt, ncl, nDir, bcl, bDir, t1, t2, t3, t4, x0, y0, z0},
  {x0, y0, z0} = r;
  t1 = AbsoluteTime[];
  ncl = newtonPot2D[x0, y0, z0];
  t1 = AbsoluteTime[] - t1;
  t2 = AbsoluteTime[];
  nDir = NIntegrate[1/Sqrt[(x0 - xpp)^2 + (y0 - ypp)^2 + (z0 - zpp)^2],
    {xpp, -aN/2, aN/2}, {ypp, -bN/2, bN/2}, {zpp, -cN/2, cN/2},
    WorkingPrecision -> 30, AccuracyGoal -> 22, MaxRecursion -> 12];
  t2 = AbsoluteTime[] - t2;
  t3 = AbsoluteTime[];
  bcl = bxPot2D[x0, y0, z0];
  t3 = AbsoluteTime[] - t3;
  t4 = AbsoluteTime[];
  bDir = NIntegrate[xpp/Sqrt[(x0 - xpp)^2 + (y0 - ypp)^2 + (z0 - zpp)^2],
    {xpp, -aN/2, aN/2}, {ypp, -bN/2, bN/2}, {zpp, -cN/2, cN/2},
    WorkingPrecision -> 30, AccuracyGoal -> 22, MaxRecursion -> 12];
  t4 = AbsoluteTime[] - t4;
  Print["  r = ", N[r*1000, 6], " mm"];
  Print["    Newton:  closed+2D = ", N[ncl, 12], " (", N[t1, 4], "s)"];
  Print["             3D NI    = ", N[nDir, 12], " (", N[t2, 4], "s)"];
  Print["             rel diff = ", N[Abs[(ncl - nDir)/nDir], 6]];
  Print["    B_x:     closed+2D = ", N[bcl, 12], " (", N[t3, 4], "s)"];
  Print["             3D NI    = ", N[bDir, 12], " (", N[t4, 4], "s)"];
  If[Abs[bDir] > 10^-30,
    Print["             rel diff = ", N[Abs[(bcl - bDir)/bDir], 6]],
    Print["             diff = ", N[bcl - bDir, 12], " (bDir near zero)"]]],
  {pt, testPts}];
Print[""];

Print["=== Summary ==="];
Print["* If above rel diffs ~ 1e-10, we have a working B_x(r)."];
Print["* Then T_xx = integral B_x^2 dV (3D NIntegrate) gives alpha_3."];
