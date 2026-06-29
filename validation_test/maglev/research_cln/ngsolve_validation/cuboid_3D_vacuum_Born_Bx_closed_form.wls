(* Phase A: derive B_x(r) closed form via integration-by-parts decomposition.

   B_x(r) = integral_Omega x''/|r-r''| dV''
          = x_R * N(r) + F_plus(r) - F_minus(r)
   where:
     N(r)       = integral_Omega 1/|r-r''| dV''                          [classical]
     F_plus(r)  = integral_face(x''=a/2)  |r - (a/2, y'', z'')|  dy'' dz''
     F_minus(r) = integral_face(x''=-a/2) |r - (-a/2, y'', z'')| dy'' dz''

   Each component should have closed form. Try Mathematica Integrate[].
*)
$HistoryLength = 0;
prec = 60;

Print["=== Phase A: B_x(r) via IBP decomposition ==="];
Print[""];

ClearAll[xR, yR, zR, xpp, ypp, zpp, asym, bsym, csym];

(* === Step A1: Cuboid Newton potential N(r) === *)
Print["A1: Newton potential N(r) = integral_Omega 1/|r-rpp| dVpp"];
Print[""];

(* Inner xpp integral: ∫ 1/sqrt((xR-xpp)^2 + R^2) dxpp = arcsinh((xpp-xR)/R) evaluated *)
Print["  A1a: xpp integral"];
N_xpp = TimeConstrained[
  Integrate[1/Sqrt[(xR - xpp)^2 + Rsq], {xpp, -asym/2, asym/2},
    Assumptions -> {asym > 0, Rsq > 0, xR \[Element] Reals}],
  60, $Failed];
If[N_xpp === $Failed, Print["  TIMEOUT"]; Quit[]];
Print["    LeafCount = ", LeafCount[N_xpp]];

Print["  A1b: ypp integral"];
N_xppSubst = N_xpp /. Rsq -> (yR - ypp)^2 + (zR - zpp)^2;
N_ypp = TimeConstrained[
  Integrate[N_xppSubst, {ypp, -bsym/2, bsym/2},
    Assumptions -> {asym > 0, bsym > 0, csym > 0,
      xR \[Element] Reals, yR \[Element] Reals, zR \[Element] Reals,
      zpp \[Element] Reals}],
  300, $Failed];
If[N_ypp === $Failed,
  Print["    ypp TIMEOUT"]; Quit[]];
Print["    LeafCount = ", LeafCount[N_ypp]];

Print["  A1c: zpp integral"];
N_zpp = TimeConstrained[
  Integrate[N_ypp, {zpp, -csym/2, csym/2},
    Assumptions -> {asym > 0, bsym > 0, csym > 0,
      xR \[Element] Reals, yR \[Element] Reals, zR \[Element] Reals}],
  600, $Failed];
If[N_zpp === $Failed,
  Print["    zpp TIMEOUT"]; Quit[]];
Print["    LeafCount = ", LeafCount[N_zpp]];
Print["  N(r) closed form FOUND"];
Print[""];

(* === Step A2: Face integral F_plus and F_minus === *)
Print["A2: Face integrals F_plus, F_minus"];
Print[""];

Print["  A2a: F_plus inner ypp integral of |r - (a/2, ypp, zpp)|"];
F_plus_ypp = TimeConstrained[
  Integrate[Sqrt[(xR - asym/2)^2 + (yR - ypp)^2 + (zR - zpp)^2],
    {ypp, -bsym/2, bsym/2},
    Assumptions -> {asym > 0, bsym > 0, csym > 0,
      xR \[Element] Reals, yR \[Element] Reals, zR \[Element] Reals,
      zpp \[Element] Reals}],
  120, $Failed];
If[F_plus_ypp === $Failed, Print["    F_plus ypp TIMEOUT"]; Quit[]];
Print["    LeafCount = ", LeafCount[F_plus_ypp]];

Print["  A2b: F_plus zpp integral"];
F_plus = TimeConstrained[
  Integrate[F_plus_ypp, {zpp, -csym/2, csym/2},
    Assumptions -> {asym > 0, bsym > 0, csym > 0,
      xR \[Element] Reals, yR \[Element] Reals, zR \[Element] Reals}],
  300, $Failed];
If[F_plus === $Failed, Print["    F_plus zpp TIMEOUT"]; Quit[]];
Print["    F_plus closed form FOUND"];
Print["    LeafCount = ", LeafCount[F_plus]];

(* F_minus by symmetry: replace asym/2 -> -asym/2, i.e. (xR - asym/2) -> (xR + asym/2) *)
F_minus = F_plus /. (xR - asym/2) -> (xR + asym/2) // Simplify;
Print["    F_minus (via x sign flip) LeafCount = ", LeafCount[F_minus]];
Print[""];

(* === Step A3: Combine === *)
Print["A3: B_x(r) = x_R * N(r) + F_plus(r) - F_minus(r)"];
B_x_closed = xR * N_zpp + F_plus - F_minus;
B_x_closed = Simplify[B_x_closed,
  Assumptions -> {asym > 0, bsym > 0, csym > 0,
    xR \[Element] Reals, yR \[Element] Reals, zR \[Element] Reals}];
Print["  B_x(r) closed form LeafCount = ", LeafCount[B_x_closed]];
Print[""];

(* === Step A4: Numerical verification === *)
Print["A4: numerical verification"];
sigmaN = SetPrecision[58000000, prec];
mu0N = SetPrecision[4 Pi/10^7, prec];
aN = SetPrecision[5/1000, prec];
bN = SetPrecision[2/1000, prec];
cN = SetPrecision[1/1000, prec];
PiN = SetPrecision[Pi, prec];

(* Define numerical B_x function *)
B_x_func[xR_?NumberQ, yR_?NumberQ, zR_?NumberQ] := Re[
  N[B_x_closed /. {asym -> aN, bsym -> bN, csym -> cN,
    xR -> SetPrecision[xR, prec],
    yR -> SetPrecision[yR, prec],
    zR -> SetPrecision[zR, prec]}, 30]];

(* Direct numerical for comparison *)
B_x_direct[xR_?NumberQ, yR_?NumberQ, zR_?NumberQ] := NIntegrate[
  xpp/Sqrt[(xR - xpp)^2 + (yR - ypp)^2 + (zR - zpp)^2],
  {xpp, -aN/2, aN/2}, {ypp, -bN/2, bN/2}, {zpp, -cN/2, cN/2},
  WorkingPrecision -> 30, AccuracyGoal -> 22, PrecisionGoal -> 20,
  MaxRecursion -> 12];

testPts = {{aN/4, bN/4, cN/4}, {aN/3, bN/5, 0}, {aN/2 - aN/100, 0, 0}};
Do[Module[{r = pt, fcl, fdir, t1, t2},
  t1 = AbsoluteTime[];
  fcl = B_x_func @@ r;
  t1 = AbsoluteTime[] - t1;
  t2 = AbsoluteTime[];
  fdir = B_x_direct @@ r;
  t2 = AbsoluteTime[] - t2;
  Print["  r = ", N[r, 6]];
  Print["    closed form B_x = ", N[fcl, 14], "  (", N[t1, 4], "s)"];
  Print["    direct 3D NI    = ", N[fdir, 14], "  (", N[t2, 4], "s)"];
  Print["    rel diff = ", N[Abs[(fcl - fdir)/fdir], 6]]],
  {pt, testPts}];
Print[""];

(* === Step A5: Save closed form for use in alpha_3 computation === *)
Print["A5: saving closed form to file"];
Export["B_x_closed_form.m", B_x_closed];
Print["  Saved to B_x_closed_form.m"];
Print[""];

Print["=== Phase A summary ==="];
Print["* B_x(r) full closed form OBTAINED"];
Print["* LeafCount = ", LeafCount[B_x_closed]];
Print["* Numerical verification: see ratios above"];
Print["* Next: T_xx = integral B_x(r)^2 dV via 3D NIntegrate (single layer)"];
