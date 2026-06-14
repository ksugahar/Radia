(* Stage 5 Born alpha_3 v2: T_xx = ||B_x||^2 reduction.

   Key insight: 9D integral reduces to 3D self-energy:
     T_xx = integral_Omega [B_x(r)]^2 dV
     B_x(r) = integral_Omega x''/|r-r''| dV''  [Aharoni-type 3D]

   Strategy:
     (1) Get B_x(r) closed form via Mathematica Integrate (harder but cleaner)
         OR partial closed form with NIntegrate fallback
     (2) Compute T_xx = NIntegrate over Omega of B_x(r)^2
     (3) alpha_3 = -(mu_0 sigma)^3/(64 pi^2 V) (T_xx + T_yy)
*)
$HistoryLength = 0;
prec = 60;

Print["=== Stage 5 Born alpha_3 v2: T_xx = ||B_x||^2 ==="];
Print[""];

(* Constants *)
sigmaN = SetPrecision[58000000, prec];
mu0N = SetPrecision[4 Pi/10^7, prec];
aN = SetPrecision[5/1000, prec];
bN = SetPrecision[2/1000, prec];
cN = SetPrecision[1/1000, prec];
PiN = SetPrecision[Pi, prec];
VN = aN bN cN;

(* === Step 1: B_x(r) closed form via iterated Integrate === *)
Print["Step 1: deriving B_x(r) closed form"];
Print[""];

ClearAll[xx, yy, zz, xpp, ypp, zpp, asym, bsym, csym, R, Rsq];

(* Inner xpp closed form *)
inner1 = TimeConstrained[
  Integrate[xpp/Sqrt[(xx - xpp)^2 + Rsq], {xpp, -asym/2, asym/2},
    Assumptions -> {asym > 0, Rsq > 0, xx \[Element] Reals,
                    -asym/2 < xx, xx < asym/2}],
  90, $Failed];
If[inner1 === $Failed, Print["xpp closed: TIMEOUT"]; Quit[]];
Print["Step 1a (xpp closed form): FOUND"];
inner1Subst = inner1 /. Rsq -> (yy - ypp)^2 + (zz - zpp)^2;
Print["  inner1 LeafCount = ", LeafCount[inner1Subst]];
Print[""];

(* Try ypp with longer timeout and more careful assumptions *)
(* Re-write to make Mathematica's life easier *)
inner1Simplified = Simplify[inner1Subst,
  Assumptions -> {asym > 0, bsym > 0, csym > 0, ypp \[Element] Reals,
                  zpp \[Element] Reals, yy \[Element] Reals, zz \[Element] Reals}];
Print["After Simplify, LeafCount = ", LeafCount[inner1Simplified]];

Print["Step 1b: trying ypp Integrate (300s timeout)..."];
inner2 = TimeConstrained[
  Integrate[inner1Simplified, {ypp, -bsym/2, bsym/2},
    Assumptions -> {asym > 0, bsym > 0, csym > 0,
                    yy \[Element] Reals, zz \[Element] Reals, zpp \[Element] Reals,
                    -bsym/2 < yy, yy < bsym/2,
                    -csym/2 < zz, zz < csym/2,
                    -csym/2 < zpp, zpp < csym/2}],
  300, $Failed];

If[inner2 === $Failed,
  Print["  ypp closed form: TIMEOUT (300s)"];
  Print["  Falling back to numerical B_x(r)"];
  bxClosed = False,
  Print["  ypp closed form: FOUND"];
  Print["  inner2 LeafCount = ", LeafCount[inner2]];
  bxClosed = "ypp_done"];
Print[""];

If[bxClosed === "ypp_done",
  Print["Step 1c: trying zpp Integrate (300s timeout)..."];
  inner3 = TimeConstrained[
    Integrate[inner2, {zpp, -csym/2, csym/2},
      Assumptions -> {asym > 0, bsym > 0, csym > 0,
                      xx \[Element] Reals, yy \[Element] Reals, zz \[Element] Reals,
                      -asym/2 < xx, xx < asym/2,
                      -bsym/2 < yy, yy < bsym/2,
                      -csym/2 < zz, zz < csym/2}],
    300, $Failed];
  If[inner3 === $Failed,
    Print["  zpp: TIMEOUT (300s) -- B_x has 1+1D closed form, zpp numerical"];
    bxClosed = "ypp_done",
    Print["  zpp: FULL closed form FOUND"];
    Print["  inner3 LeafCount = ", LeafCount[inner3]];
    bxClosed = "full"]];
Print[""];

(* === Step 2: Compute T_xx = integral B_x^2 dV === *)
Print["Step 2: T_xx computation"];
Print[""];

(* Numerical B_x via NIntegrate with available closed form layers *)
Bx[xR_?NumberQ, yR_?NumberQ, zR_?NumberQ] := Module[{result},
  result = NIntegrate[
    xpp/Sqrt[(xR - xpp)^2 + (yR - ypp)^2 + (zR - zpp)^2],
    {xpp, -aN/2, aN/2}, {ypp, -bN/2, bN/2}, {zpp, -cN/2, cN/2},
    WorkingPrecision -> 30, AccuracyGoal -> 22, PrecisionGoal -> 20,
    MaxRecursion -> 12];
  result];

(* Sample evaluation *)
Print["Sample B_x(r):"];
sampleR = {{aN/4, 0, 0}, {0, bN/4, 0}, {0, 0, cN/4}, {aN/4, bN/4, cN/4}};
Do[Module[{r, val}, r = sampleR[[i]]; val = Bx @@ r;
  Print["  B_x(", N[r, 6], ") = ", N[val, 12]]],
  {i, Length[sampleR]}];
Print[""];

(* Use parity: B_x odd in x, even in y, z. So integrate B_x^2 over [0,a/2] x [0,b/2] x [0,c/2]
   and multiply by 8. *)
Print["Computing T_xx = 8 * integral_{[0,a/2]x[0,b/2]x[0,c/2]} B_x(r)^2 dV"];
Print["  WP=30 outer, WP=30 inner (each B_x evaluation = 3D NIntegrate)"];
Print[""];

tStart = AbsoluteTime[];
TxxNum = TimeConstrained[
  8 NIntegrate[Bx[xx, yy, zz]^2,
    {xx, 0, aN/2}, {yy, 0, bN/2}, {zz, 0, cN/2},
    WorkingPrecision -> 25, AccuracyGoal -> 18, PrecisionGoal -> 15,
    MaxRecursion -> 8, MinRecursion -> 4],
  3600, $Failed];

If[TxxNum === $Failed,
  Print["T_xx computation TIMEOUT (3600s)"]; Quit[],
  Print["T_xx (3D outer + nested NIntegrate B_x) = ", N[TxxNum, 16],
    "  (", N[AbsoluteTime[] - tStart, 4], " s)"]];
Print[""];

(* Now T_yy via swap a <-> b *)
By[xR_?NumberQ, yR_?NumberQ, zR_?NumberQ] := Module[{result},
  result = NIntegrate[
    ypp/Sqrt[(xR - xpp)^2 + (yR - ypp)^2 + (zR - zpp)^2],
    {xpp, -aN/2, aN/2}, {ypp, -bN/2, bN/2}, {zpp, -cN/2, cN/2},
    WorkingPrecision -> 30, AccuracyGoal -> 22, PrecisionGoal -> 20,
    MaxRecursion -> 12];
  result];

Print["Computing T_yy = 8 * integral B_y^2 dV"];
tStart = AbsoluteTime[];
TyyNum = TimeConstrained[
  8 NIntegrate[By[xx, yy, zz]^2,
    {xx, 0, aN/2}, {yy, 0, bN/2}, {zz, 0, cN/2},
    WorkingPrecision -> 25, AccuracyGoal -> 18, PrecisionGoal -> 15,
    MaxRecursion -> 8, MinRecursion -> 4],
  3600, $Failed];

If[TyyNum === $Failed,
  Print["T_yy TIMEOUT"]; Quit[],
  Print["T_yy = ", N[TyyNum, 16], "  (", N[AbsoluteTime[] - tStart, 4], " s)"]];
Print[""];

(* === Step 3: alpha_3 from formula === *)
alpha3Born = -(mu0N sigmaN)^3/(64 PiN^2 VN) (TxxNum + TyyNum);
Print["=== alpha_3 Born ==="];
Print["  T_xx                 = ", N[TxxNum, 12]];
Print["  T_yy                 = ", N[TyyNum, 12]];
Print["  T_xx/T_yy ratio      = ", N[TxxNum/TyyNum, 8]];
Print["  alpha_3^Born         = ", N[alpha3Born, 16]];
Print[""];

(* Compare to expectations *)
ReffSq = 5 (aN^2 + bN^2)/24;
alpha3Sphere = -(mu0N sigmaN ReffSq)^3/1050;
Print["  alpha_3^sphere(R_eff) = ", N[alpha3Sphere, 16]];
Print["  ratio Born/Sphere     = ", N[alpha3Born/alpha3Sphere, 8]];
Print["  (alpha_2 ratio was 0.405; alpha_2^2 = ", N[0.405^2, 6],
    ")"];
Print[""];

(* Compare to MC estimate *)
Print["  Earlier 9D MC estimate (T_yy approx via T_xx/(a/b)^4): -1.37e-14"];
Print["  This v2 result: ", N[alpha3Born, 8], "  (more accurate)"];
