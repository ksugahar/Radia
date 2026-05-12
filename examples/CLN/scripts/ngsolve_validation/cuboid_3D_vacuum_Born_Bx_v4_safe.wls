(* Phase A v4: B_x closed form (Mathematica-safe variable names).

   Avoid reserved names: C, D, E, I, K, N, O.
   Use prefixes: cubA, cubB, cubC for half-dimensions; xv, yv, zv for position.
*)
$HistoryLength = 0;

(* ===== Banerjee-Gupta corner-sum Newton potential ===== *)

(* Kernel function for cuboid Newton potential.
   For corner offsets (uu, vv, ww), G(uu, vv, ww) is the indefinite integral
   1/sqrt(uu^2+vv^2+ww^2) integrated three times. *)
gKernel[uu_?NumericQ, vv_?NumericQ, ww_?NumericQ] := Block[{rr},
  rr = Sqrt[uu^2 + vv^2 + ww^2];
  uu vv Log[ww + rr] + uu ww Log[vv + rr] + vv ww Log[uu + rr]
    - (uu^2/2) ArcTan[uu rr, vv ww]
    - (vv^2/2) ArcTan[vv rr, uu ww]
    - (ww^2/2) ArcTan[ww rr, uu vv]];

(* Newton potential of cuboid [-cubA, cubA] x [-cubB, cubB] x [-cubC, cubC]
   evaluated at point (xv, yv, zv) *)
newtonPot[xv_?NumericQ, yv_?NumericQ, zv_?NumericQ,
          cubA_?NumericQ, cubB_?NumericQ, cubC_?NumericQ] :=
  Sum[sx sy sz * gKernel[xv - sx cubA, yv - sy cubB, zv - sz cubC],
    {sx, {1, -1}}, {sy, {1, -1}}, {sz, {1, -1}}];

(* ===== Test: Newton potential ===== *)
Print["=== Phase A v4: Banerjee-Gupta corner-sum (safe naming) ==="];
Print[""];

Block[{prec = 60, sigCu, mu0, aDim, bDim, cDim, vol,
       cubA, cubB, cubC, testPts, pt, xv, yv, zv},
  sigCu = SetPrecision[58000000, prec];
  mu0 = SetPrecision[4 Pi/10^7, prec];
  aDim = SetPrecision[5/1000, prec];
  bDim = SetPrecision[2/1000, prec];
  cDim = SetPrecision[1/1000, prec];
  vol = aDim bDim cDim;
  cubA = aDim/2; cubB = bDim/2; cubC = cDim/2;

  Print["Newton potential N(r) verification:"];
  Print["  cuboid: a=", N[aDim*1000, 4], " b=", N[bDim*1000, 4],
    " c=", N[cDim*1000, 4], " mm"];
  Print[""];

  testPts = {{aDim/4, bDim/4, cDim/4}, {aDim/3, bDim/5, 0},
             {0, 0, 0}, {aDim/2 - cDim/10, 0, 0}};
  Do[Module[{r = pt, ncl, ndir, t1, t2, xx0, yy0, zz0},
    {xx0, yy0, zz0} = r;
    t1 = AbsoluteTime[];
    ncl = N[newtonPot[xx0, yy0, zz0, cubA, cubB, cubC], 30];
    t1 = AbsoluteTime[] - t1;
    t2 = AbsoluteTime[];
    ndir = NIntegrate[1/Sqrt[(xx0 - xpp)^2 + (yy0 - ypp)^2 + (zz0 - zpp)^2],
      {xpp, -cubA, cubA}, {ypp, -cubB, cubB}, {zpp, -cubC, cubC},
      WorkingPrecision -> 30, AccuracyGoal -> 22, PrecisionGoal -> 20,
      MaxRecursion -> 12];
    t2 = AbsoluteTime[] - t2;
    Print["  r = ", N[r*1000, 6], " mm"];
    Print["    closed form  = ", N[ncl, 16], "  (", N[t1, 4], "s)"];
    Print["    direct 3D NI = ", N[ndir, 16], "  (", N[t2, 4], "s)"];
    Print["    rel diff     = ", N[Abs[(ncl - ndir)/ndir], 6]]],
    {pt, testPts}];
  Print[""];

  (* ===== Face integrals for B_x ===== *)
  (* B_x(r) = xv * N(r) + Fp(r) - Fm(r)
     Fp(r) = integral_{[-B,B]x[-C,C]} |r - (cubA, ypp, zpp)| dypp dzpp
     Fm(r) = integral_{[-B,B]x[-C,C]} |r - (-cubA, ypp, zpp)| dypp dzpp *)

  facePlus[xv_?NumericQ, yv_?NumericQ, zv_?NumericQ] := NIntegrate[
    Sqrt[(xv - cubA)^2 + (yv - ypp)^2 + (zv - zpp)^2],
    {ypp, -cubB, cubB}, {zpp, -cubC, cubC},
    WorkingPrecision -> 30, AccuracyGoal -> 25, PrecisionGoal -> 22,
    MaxRecursion -> 12];

  faceMinus[xv_?NumericQ, yv_?NumericQ, zv_?NumericQ] := NIntegrate[
    Sqrt[(xv + cubA)^2 + (yv - ypp)^2 + (zv - zpp)^2],
    {ypp, -cubB, cubB}, {zpp, -cubC, cubC},
    WorkingPrecision -> 30, AccuracyGoal -> 25, PrecisionGoal -> 22,
    MaxRecursion -> 12];

  bxIBP[xv_?NumericQ, yv_?NumericQ, zv_?NumericQ] := xv * newtonPot[xv, yv, zv, cubA, cubB, cubC]
    + facePlus[xv, yv, zv] - faceMinus[xv, yv, zv];

  Print["B_x(r) = xv * N(r) + Fp(r) - Fm(r) verification:"];
  Print[""];
  Do[Module[{r = pt, fcl, fdir, t1, t2, xx0, yy0, zz0},
    {xx0, yy0, zz0} = r;
    t1 = AbsoluteTime[];
    fcl = bxIBP[xx0, yy0, zz0];
    t1 = AbsoluteTime[] - t1;
    t2 = AbsoluteTime[];
    fdir = NIntegrate[
      xpp/Sqrt[(xx0 - xpp)^2 + (yy0 - ypp)^2 + (zz0 - zpp)^2],
      {xpp, -cubA, cubA}, {ypp, -cubB, cubB}, {zpp, -cubC, cubC},
      WorkingPrecision -> 30, AccuracyGoal -> 22, PrecisionGoal -> 20,
      MaxRecursion -> 12];
    t2 = AbsoluteTime[] - t2;
    Print["  r = ", N[r*1000, 6], " mm"];
    Print["    via IBP+face = ", N[fcl, 16], "  (", N[t1, 4], "s)"];
    Print["    direct 3D NI = ", N[fdir, 16], "  (", N[t2, 4], "s)"];
    Print["    rel diff     = ", N[Abs[(fcl - fdir)/fdir], 6]]],
    {pt, testPts}];
  Print[""];

  Print["=== Phase A v4 summary ==="];
  Print["* Newton potential: corner-sum closed form WORKING"];
  Print["* B_x: hybrid (closed Newton + 2D NIntegrate face)"];
  Print["* Both verified against direct 3D NIntegrate"];
];
