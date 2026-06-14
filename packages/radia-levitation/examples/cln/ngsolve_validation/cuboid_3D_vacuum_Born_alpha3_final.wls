(* Compute alpha_3 Born using v7 hybrid (1D closed + 2D NIntegrate). *)
$HistoryLength = 0;
ClearAll[uu, ypp, zpp, xx, yy, zz, xpp, xR, yR, zR];

Print["=== alpha_3 Born final: 3D NI of B_x^2 ==="];
Print[""];

prec = 50;
sigCu = SetPrecision[58000000, prec];
mu0 = SetPrecision[4 Pi/10^7, prec];
aN = SetPrecision[5/1000, prec];
bN = SetPrecision[2/1000, prec];
cN = SetPrecision[1/1000, prec];
piN = SetPrecision[Pi, prec];
vol = aN bN cN;

(* Closed-form inner integrals *)
bxInnerXpp[xRv_, rsqv_, aHv_] := Block[{r1, r2},
  r1 = Sqrt[(aHv - xRv)^2 + rsqv];
  r2 = Sqrt[(aHv + xRv)^2 + rsqv];
  (r1 - r2) + xRv * (ArcSinh[(aHv - xRv)/Sqrt[rsqv]] + ArcSinh[(aHv + xRv)/Sqrt[rsqv]])];

bxInnerYpp[yRv_, rsqv_, bHv_] := Block[{r1, r2},
  r1 = Sqrt[(bHv - yRv)^2 + rsqv];
  r2 = Sqrt[(bHv + yRv)^2 + rsqv];
  (r1 - r2) + yRv * (ArcSinh[(bHv - yRv)/Sqrt[rsqv]] + ArcSinh[(bHv + yRv)/Sqrt[rsqv]])];

(* B_x(r) = integral xpp/|r-r''| dV'' over conductor *)
Bx[xRv_?NumericQ, yRv_?NumericQ, zRv_?NumericQ] :=
  NIntegrate[bxInnerXpp[xRv, (yRv - ypp)^2 + (zRv - zpp)^2, aN/2],
    {ypp, -bN/2, bN/2}, {zpp, -cN/2, cN/2},
    WorkingPrecision -> 25, AccuracyGoal -> 18, PrecisionGoal -> 15,
    MaxRecursion -> 10];

(* B_y(r) = integral ypp/|r-r''| dV'' (symmetric in xpp <-> ypp) *)
By[xRv_?NumericQ, yRv_?NumericQ, zRv_?NumericQ] :=
  NIntegrate[bxInnerYpp[yRv, (xRv - xpp)^2 + (zRv - zpp)^2, bN/2],
    {xpp, -aN/2, aN/2}, {zpp, -cN/2, cN/2},
    WorkingPrecision -> 25, AccuracyGoal -> 18, PrecisionGoal -> 15,
    MaxRecursion -> 10];

(* T_xx = integral over Omega of B_x(r)^2 dV. Use parity: B_x odd in x, even in y, z.
   So B_x^2 even in all. Integrate over [0, a/2] x [0, b/2] x [0, c/2] times 8.
   Avoid x=0 axis (B_x=0 there) via small offset. *)
Print["Computing T_xx = 8 * integral B_x^2 over [eps, a/2] x [0, b/2] x [0, c/2]"];
Print["  Outer WP=15, AG=8, MR=6 (inner B_x at WP=25 per evaluation)"];
epsX = aN/200;
tStart = AbsoluteTime[];
TxxBorn = TimeConstrained[
  8 NIntegrate[Bx[xx, yy, zz]^2,
    {xx, epsX, aN/2}, {yy, 0, bN/2}, {zz, 0, cN/2},
    WorkingPrecision -> 15, AccuracyGoal -> 8, PrecisionGoal -> 6,
    MaxRecursion -> 6, MinRecursion -> 3,
    Method -> "AdaptiveMonteCarlo"],
  3600, $Failed];
Print["  T_xx = ", N[TxxBorn, 12], "  (", N[AbsoluteTime[] - tStart, 4], "s)"];

If[TxxBorn === $Failed,
  Print["  TIMEOUT, retry with MC..."];
  TxxBorn = 8 NIntegrate[Bx[xx, yy, zz]^2,
    {xx, epsX, aN/2}, {yy, 0, bN/2}, {zz, 0, cN/2},
    Method -> "MonteCarlo", MaxPoints -> 1000];
  Print["  T_xx (MC) = ", N[TxxBorn, 8]]];

Print[""];
Print["Computing T_yy similarly..."];
epsY = bN/200;
tStart = AbsoluteTime[];
TyyBorn = TimeConstrained[
  8 NIntegrate[By[xx, yy, zz]^2,
    {xx, 0, aN/2}, {yy, epsY, bN/2}, {zz, 0, cN/2},
    WorkingPrecision -> 15, AccuracyGoal -> 8, PrecisionGoal -> 6,
    MaxRecursion -> 6, MinRecursion -> 3,
    Method -> "AdaptiveMonteCarlo"],
  3600, $Failed];
Print["  T_yy = ", N[TyyBorn, 12], "  (", N[AbsoluteTime[] - tStart, 4], "s)"];

If[TyyBorn === $Failed,
  TyyBorn = 8 NIntegrate[By[xx, yy, zz]^2,
    {xx, 0, aN/2}, {yy, epsY, bN/2}, {zz, 0, cN/2},
    Method -> "MonteCarlo", MaxPoints -> 1000]];

Print[""];

(* alpha_3 *)
alpha3 = -(mu0 sigCu)^3/(64 piN^2 vol) (TxxBorn + TyyBorn);
Print["=== Result ==="];
Print["  T_xx = ", N[TxxBorn, 12]];
Print["  T_yy = ", N[TyyBorn, 12]];
Print["  alpha_3^Born = ", N[alpha3, 14]];
Print[""];

(* References *)
ReffSq = 5 (aN^2 + bN^2)/24;
alpha3Sphere = -(mu0 sigCu ReffSq)^3/1050;
Print["  alpha_3^sphere(R_eff) = ", N[alpha3Sphere, 12]];
Print["  ratio Born/Sphere = ", N[alpha3/alpha3Sphere, 8]];
Print["  Earlier 9D MC estimate: -1.37e-14"];
