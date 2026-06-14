(* Stage 5 Born alpha_3 v3: use Step 1a closed-form xpp integral.

   Avoid the slow 3D inner NIntegrate by hardcoding inner1 (xpp) closed form:
     inner1(xR, R) = sqrt((a/2-xR)^2 + R^2) - sqrt((a/2+xR)^2 + R^2)
                  + xR (arcsinh((a/2-xR)/R) + arcsinh((a/2+xR)/R))
     where R^2 = (yR-ypp)^2 + (zR-zpp)^2.

   Then B_x(r) = 2D NIntegrate (ypp, zpp) -- much faster.
   T_xx = NIntegrate[B_x(r)^2, r in Omega] -- 3D outer.
*)
$HistoryLength = 0;
prec = 60;

Print["=== Stage 5 Born alpha_3 v3: closed inner1 + 2D NIntegrate ==="];
Print[""];

(* Constants *)
sigmaN = SetPrecision[58000000, prec];
mu0N = SetPrecision[4 Pi/10^7, prec];
aN = SetPrecision[5/1000, prec];
bN = SetPrecision[2/1000, prec];
cN = SetPrecision[1/1000, prec];
PiN = SetPrecision[Pi, prec];
VN = aN bN cN;

(* Closed form for inner1: ∫ x''/sqrt((x-x'')^2 + R^2) dx'' over [-a/2, a/2] *)
inner1Closed[xR_, R2_, a_] :=
  Sqrt[(a/2 - xR)^2 + R2] - Sqrt[(a/2 + xR)^2 + R2] +
  xR (ArcSinh[(a/2 - xR)/Sqrt[R2]] + ArcSinh[(a/2 + xR)/Sqrt[R2]]);

(* Verify against direct numerical *)
Print["Verification of inner1 closed form:"];
testR = {(aN/4) //N, (bN/4) //N, (cN/4) //N};
{xR, yR, zR} = testR;
xppDirect = NIntegrate[xpp/Sqrt[(xR - xpp)^2 + (yR - ypp)^2 + (zR - zpp)^2],
  {xpp, -aN/2, aN/2}, {ypp, -bN/2, bN/2}, {zpp, -cN/2, cN/2},
  WorkingPrecision -> 25, AccuracyGoal -> 18];
xppViaClosed = NIntegrate[
  inner1Closed[xR, (yR - ypp)^2 + (zR - zpp)^2, aN],
  {ypp, -bN/2, bN/2}, {zpp, -cN/2, cN/2},
  WorkingPrecision -> 25, AccuracyGoal -> 18];
Print["  B_x(", N[testR, 6], "):"];
Print["    direct 3D NIntegrate = ", N[xppDirect, 12]];
Print["    via inner1 closed + 2D = ", N[xppViaClosed, 12]];
Print["    relative diff = ", N[Abs[(xppDirect - xppViaClosed)/xppDirect], 6]];
Print[""];

(* Now use this for fast B_x evaluation *)
Bx[xR_?NumericQ, yR_?NumericQ, zR_?NumericQ] :=
  NIntegrate[inner1Closed[xR, (yR - ypp)^2 + (zR - zpp)^2, aN],
    {ypp, -bN/2, bN/2}, {zpp, -cN/2, cN/2},
    WorkingPrecision -> 30, AccuracyGoal -> 22, PrecisionGoal -> 20,
    MaxRecursion -> 12];

By[xR_?NumericQ, yR_?NumericQ, zR_?NumericQ] :=
  NIntegrate[inner1Closed[yR, (xR - xpp)^2 + (zR - zpp)^2, bN],
    {xpp, -aN/2, aN/2}, {zpp, -cN/2, cN/2},
    WorkingPrecision -> 30, AccuracyGoal -> 22, PrecisionGoal -> 20,
    MaxRecursion -> 12];

(* Sample timings *)
Print["Sample B_x timings (2D NIntegrate per call):"];
Do[Module[{r = pt, t, val},
  t = AbsoluteTime[];
  val = Bx @@ r;
  Print["  B_x(", N[r, 6], ") = ", N[val, 14], "  (", N[AbsoluteTime[] - t, 4], " s)"]],
  {pt, {{aN/4, bN/4, cN/4}, {aN/3, 0, cN/4}, {0, bN/3, 0}, {aN/2 - 0.0001, 0, 0}}}];
Print[""];

(* === T_xx = integral B_x^2 over Omega ===
   Use parity: B_x odd in x, even in y, z. So B_x^2 even in all.
   Integrate over [0, a/2] x [0, b/2] x [0, c/2] times 8.
   Avoid x=0 axis (B_x=0 there) via lower bound > 0. *)

Print["Computing T_xx = 8 * integral_{[eps, a/2]x[0, b/2]x[0, c/2]} B_x^2 dV"];
Print["  outer WP=20, inner Bx WP=30"];
Print[""];

epsX = SetPrecision[10^(-7), prec];  (* avoid axis where B_x = 0 *)
tStart = AbsoluteTime[];
TxxBornV3 = TimeConstrained[
  8 NIntegrate[Bx[xx, yy, zz]^2,
    {xx, epsX, aN/2}, {yy, 0, bN/2}, {zz, 0, cN/2},
    WorkingPrecision -> 20, AccuracyGoal -> 12, PrecisionGoal -> 10,
    MaxRecursion -> 8, MinRecursion -> 4,
    Method -> "AdaptiveMonteCarlo"],
  1800, $Failed];
Print["  T_xx = ", N[TxxBornV3, 12], "  (", N[AbsoluteTime[] - tStart, 4], " s)"];

If[TxxBornV3 === $Failed,
  Print["  T_xx TIMEOUT, try Monte Carlo with fewer points"];
  TxxBornV3 = NIntegrate[Bx[xx, yy, zz]^2,
    {xx, epsX, aN/2}, {yy, 0, bN/2}, {zz, 0, cN/2},
    Method -> "MonteCarlo", MaxPoints -> 2000];
  TxxBornV3 = 8 TxxBornV3;
  Print["  T_xx (MC) = ", N[TxxBornV3, 8]]];

Print[""];
Print["Computing T_yy = 8 * integral B_y^2 dV (similar)"];
tStart = AbsoluteTime[];
TyyBornV3 = TimeConstrained[
  8 NIntegrate[By[xx, yy, zz]^2,
    {xx, 0, aN/2}, {yy, epsX, bN/2}, {zz, 0, cN/2},
    WorkingPrecision -> 20, AccuracyGoal -> 12, PrecisionGoal -> 10,
    MaxRecursion -> 8, MinRecursion -> 4,
    Method -> "AdaptiveMonteCarlo"],
  1800, $Failed];
Print["  T_yy = ", N[TyyBornV3, 12], "  (", N[AbsoluteTime[] - tStart, 4], " s)"];

If[TyyBornV3 === $Failed,
  Print["  T_yy TIMEOUT, MC fallback"];
  TyyBornV3 = 8 NIntegrate[By[xx, yy, zz]^2,
    {xx, 0, aN/2}, {yy, epsX, bN/2}, {zz, 0, cN/2},
    Method -> "MonteCarlo", MaxPoints -> 2000];
  Print["  T_yy (MC) = ", N[TyyBornV3, 8]]];

Print[""];
alpha3Born = -(mu0N sigmaN)^3/(64 PiN^2 VN) (TxxBornV3 + TyyBornV3);
Print["=== alpha_3 Born ==="];
Print["  T_xx = ", N[TxxBornV3, 12]];
Print["  T_yy = ", N[TyyBornV3, 12]];
Print["  alpha_3^Born = ", N[alpha3Born, 14]];
Print[""];

(* Reference *)
ReffSq = 5 (aN^2 + bN^2)/24;
alpha3Sphere = -(mu0N sigmaN ReffSq)^3/1050;
Print["  alpha_3^sphere(R_eff) = ", N[alpha3Sphere, 14]];
Print["  ratio Born/Sphere = ", N[alpha3Born/alpha3Sphere, 8]];
Print["  (alpha_2 ratio was 0.405; alpha_2^2 = 0.164)"];
Print[""];

(* Save Foster-relevant info *)
Print["=== Updated Foster moments (Born series) ==="];
alphaLF = -mu0N sigmaN (aN^2 + bN^2)/48;
alpha2Born = SetPrecision[7.4783240535434074691212813798857319216568135873611*^-10, prec];
Print["  alpha_1^Born (closed) = ", N[alphaLF, 16]];
Print["  alpha_2^Born (closed) = ", N[alpha2Born, 16]];
Print["  alpha_3^Born (numeric)= ", N[alpha3Born, 16]];
