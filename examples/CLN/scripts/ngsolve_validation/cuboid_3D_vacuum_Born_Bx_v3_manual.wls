(* Phase A v3: B_x via manual classical closed forms.

   Don't rely on Mathematica auto-Integrate; instead use known formulas:
   - Newton potential of cuboid (Banerjee-Gupta, Nagy): corner sum of G(u,v,w)
   - B_x(r) = x_R * N(r) + face integrals, face integrals also classical.

   Verify against direct 3D NIntegrate.
*)
$HistoryLength = 0;
prec = 60;

Print["=== Phase A v3: manual closed forms ==="];
Print[""];

(* Constants *)
sigmaN = SetPrecision[58000000, prec];
mu0N = SetPrecision[4 Pi/10^7, prec];
aN = SetPrecision[5/1000, prec];
bN = SetPrecision[2/1000, prec];
cN = SetPrecision[1/1000, prec];
PiN = SetPrecision[Pi, prec];
VN = aN bN cN;

(* === Cuboid Newton potential (Banerjee-Gupta 1977 / Nagy 1966 form) === *)
(* For cuboid [-A, A] x [-B, B] x [-C, C], N(r) = integral 1/|r-r'| dV'

   N(x, y, z) = sum_{sx,sy,sz in {+,-}} sx sy sz G(x - sx A, y - sy B, z - sz C)
   where G(u, v, w) is the "kernel" function.

   Standard form (from gravity literature):
   G(u, v, w) = u v log(w + r) + u w log(v + r) + v w log(u + r)
              - (u^2/2) atan(v w / (u r))
              - (v^2/2) atan(u w / (v r))
              - (w^2/2) atan(u v / (w r))
   r = sqrt(u^2 + v^2 + w^2).
*)

gKernel[u_, v_, w_] := Module[{r = Sqrt[u^2 + v^2 + w^2]},
  u v Log[w + r] + u w Log[v + r] + v w Log[u + r]
    - (u^2/2) ArcTan[v w/(u r)]
    - (v^2/2) ArcTan[u w/(v r)]
    - (w^2/2) ArcTan[u v/(w r)]];

(* Newton potential: corner sum *)
newtonClosed[x_, y_, z_, A_, B_, C_] := Sum[
  sx sy sz gKernel[x - sx A, y - sy B, z - sz C],
  {sx, {1, -1}}, {sy, {1, -1}}, {sz, {1, -1}}];

(* Verify against direct numerical *)
Print["Newton potential verification:"];
testPts = {{aN/4, bN/4, cN/4}, {aN/3, bN/5, 0}, {0, 0, 0}};
A = aN/2; B = bN/2; C = cN/2;
Do[Module[{r, ncl, ndir, t1, t2},
  r = pt;
  t1 = AbsoluteTime[];
  ncl = N[newtonClosed[r[[1]], r[[2]], r[[3]], A, B, C], 30];
  t1 = AbsoluteTime[] - t1;
  t2 = AbsoluteTime[];
  ndir = NIntegrate[1/Sqrt[(r[[1]] - xpp)^2 + (r[[2]] - ypp)^2 + (r[[3]] - zpp)^2],
    {xpp, -A, A}, {ypp, -B, B}, {zpp, -C, C},
    WorkingPrecision -> 25, AccuracyGoal -> 18, PrecisionGoal -> 16,
    MaxRecursion -> 12];
  t2 = AbsoluteTime[] - t2;
  Print["  r = ", N[r, 6]];
  Print["    closed form  = ", N[ncl, 14], "  (", N[t1, 4], "s)"];
  Print["    direct 3D NI = ", N[ndir, 14], "  (", N[t2, 4], "s)"];
  Print["    rel diff = ", N[Abs[(ncl - ndir)/ndir], 6]]],
  {pt, testPts}];
Print[""];

(* === Face integrals: F_plus, F_minus === *)
(* F_plus(r) = integral_{[-B,B]x[-C,C]} |r - (A, ypp, zpp)| dypp dzpp
   This is a 2D integral with classical closed form (similar Banerjee-style).

   Define h_kernel(p, q, w) = closed form of int int sqrt(p^2 + (y-ypp)^2 + (z-zpp)^2) dypp dzpp
                                                                         from corners.

   For x_R fixed, p = x_R - A (or +A).

   2D integral has form:
   integral_{u_1}^{u_2} du integral_{v_1}^{v_2} dv sqrt(p^2 + u^2 + v^2)
   = corner sum of H(u, v, p) where
   H(u, v, p) = (uv/3)(u² + v² + 4p²) ... wait no, that's wrong.

   Actually integral integral sqrt(p^2 + u^2 + v^2) du dv has known form:
   = (uv r/3) - (p²/3) [u atan(uv/(p r)) + v atan(vu/(p r))]   ... approximate

   Let me use Mathematica to derive: *)
ClearAll[uv, vv, pv];
hInner = Integrate[Sqrt[uv^2 + vv^2 + pv^2], uv,
  Assumptions -> {uv \[Element] Reals, vv \[Element] Reals, pv > 0}];
Print["Inner integral (one variable):"];
Print["  ", Simplify[hInner]];

hFull = TimeConstrained[
  Integrate[hInner, vv,
    Assumptions -> {uv \[Element] Reals, vv \[Element] Reals, pv > 0}],
  60, $Failed];
If[hFull === $Failed,
  Print["  hFull TIMEOUT — try numerical instead"],
  Print["  hFull (one antiderivative) LeafCount = ", LeafCount[hFull]]];
Print[""];

(* Face integral via NIntegrate (works fine in 2D) *)
fPlusNum[x_?NumberQ, y_?NumberQ, z_?NumberQ] := NIntegrate[
  Sqrt[(x - aN/2)^2 + (y - ypp)^2 + (z - zpp)^2],
  {ypp, -bN/2, bN/2}, {zpp, -cN/2, cN/2},
  WorkingPrecision -> 30, AccuracyGoal -> 25, PrecisionGoal -> 22,
  MaxRecursion -> 12];

fMinusNum[x_?NumberQ, y_?NumberQ, z_?NumberQ] := NIntegrate[
  Sqrt[(x + aN/2)^2 + (y - ypp)^2 + (z - zpp)^2],
  {ypp, -bN/2, bN/2}, {zpp, -cN/2, cN/2},
  WorkingPrecision -> 30, AccuracyGoal -> 25, PrecisionGoal -> 22,
  MaxRecursion -> 12];

(* B_x via IBP: B_x(r) = x_R * N(r) + F_plus(r) - F_minus(r) *)
bxFunc[x_?NumberQ, y_?NumberQ, z_?NumberQ] := Module[{nVal, fpVal, fmVal},
  nVal = newtonClosed[x, y, z, aN/2, bN/2, cN/2];
  fpVal = fPlusNum[x, y, z];
  fmVal = fMinusNum[x, y, z];
  x * nVal + fpVal - fmVal];

(* Verify B_x *)
Print["B_x verification (closed-form Newton + 2D NIntegrate face):"];
Do[Module[{r, fcl, fdir, t1, t2},
  r = pt;
  t1 = AbsoluteTime[];
  fcl = N[bxFunc[r[[1]], r[[2]], r[[3]]], 30];
  t1 = AbsoluteTime[] - t1;
  t2 = AbsoluteTime[];
  fdir = NIntegrate[xpp/Sqrt[(r[[1]] - xpp)^2 + (r[[2]] - ypp)^2 + (r[[3]] - zpp)^2],
    {xpp, -aN/2, aN/2}, {ypp, -bN/2, bN/2}, {zpp, -cN/2, cN/2},
    WorkingPrecision -> 25, AccuracyGoal -> 18, PrecisionGoal -> 16,
    MaxRecursion -> 12];
  t2 = AbsoluteTime[] - t2;
  Print["  r = ", N[r, 6]];
  Print["    via IBP+face = ", N[fcl, 14], "  (", N[t1, 4], "s)"];
  Print["    direct 3D NI = ", N[fdir, 14], "  (", N[t2, 4], "s)"];
  Print["    rel diff = ", N[Abs[(fcl - fdir)/fdir], 6]]],
  {pt, testPts}];
Print[""];

Print["=== Phase A v3 summary ==="];
Print["* Newton potential closed form: corner-sum classical formula ✓"];
Print["* Face integrals: 2D NIntegrate (much faster than 3D)"];
Print["* B_x = xR*N + Fp - Fm: hybrid closed-form + 2D NI"];
Print[""];
Print["Next: try to make face integrals also closed form for full B_x closed."];
