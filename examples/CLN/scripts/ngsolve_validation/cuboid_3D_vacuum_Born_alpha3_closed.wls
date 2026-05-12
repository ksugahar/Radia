(* Stage 5 (Born alpha_3 closed form attempt).

   Target: alpha_3 = -(mu_0 sigma)^3/(64 pi^2 V) (T_xx + T_yy)
   T_xx = integral over Omega^3 of xx''/(|r-r'| |r'-r''|) dV dV' dV''

   3-step analytical reduction:
     (1) B_x(r) := integral_Omega x''/|r-r''| dV''        [3D Aharoni-type]
     (2) G_2(r) := integral_Omega B_x(r')/|r-r'| dV'      [convolve with G]
     (3) T_xx = integral_Omega x * G_2(r) dV              [final 3D]

   Each step may be closed form via Mathematica Integrate[].
*)

$HistoryLength = 0;
prec = 60;

Print["=== Stage 5 Born alpha_3 closed-form attempt ==="];
Print[""];

(* Constants *)
sigmaN = SetPrecision[58000000, prec];
mu0N = SetPrecision[4 Pi/10^7, prec];
aN = SetPrecision[5/1000, prec];
bN = SetPrecision[2/1000, prec];
cN = SetPrecision[1/1000, prec];
PiN = SetPrecision[Pi, prec];
VN = aN bN cN;

(* === Step 1: Compute B_x(r) closed form === *)
ClearAll[xx, yy, zz, xpp, ypp, zpp, asym, bsym, csym];

Print["Step 1: B_x(r) = integral xpp/|r-rpp| dV''"];
Print[""];

(* Inner over xpp: xpp/sqrt((xx-xpp)^2 + (yy-ypp)^2 + (zz-zpp)^2) integrated over [-a/2, a/2] *)
(* Substitute u = xpp - xx, then xpp = u + xx: *)
inner1 = TimeConstrained[
  Integrate[xpp/Sqrt[(xx - xpp)^2 + Rsq2], {xpp, -asym/2, asym/2},
    Assumptions -> {asym > 0, Rsq2 > 0, xx \[Element] Reals}],
  90, $Failed];
If[inner1 === $Failed,
  Print["Step 1 inner xpp: TIMEOUT"]; Quit[]];
Print["Step 1a (inner xpp): closed form FOUND"];
Print[""];

(* Substitute Rsq2 -> (yy-ypp)^2 + (zz-zpp)^2 *)
inner1Subst = inner1 /. Rsq2 -> (yy - ypp)^2 + (zz - zpp)^2;

(* Outer over ypp *)
inner2 = TimeConstrained[
  Integrate[inner1Subst, {ypp, -bsym/2, bsym/2},
    Assumptions -> {asym > 0, bsym > 0, csym > 0,
      yy \[Element] Reals, zz \[Element] Reals, zpp \[Element] Reals,
      Abs[yy] < bsym, Abs[zz] < csym}],
  120, $Failed];
If[inner2 === $Failed,
  Print["Step 1b ypp: TIMEOUT, fall back to NIntegrate"];
  $stepStatus = "ypp_failed",
  Print["Step 1b (ypp): closed form FOUND"]; $stepStatus = "ok"];

If[$stepStatus === "ok",
  inner3 = TimeConstrained[
    Integrate[inner2, {zpp, -csym/2, csym/2},
      Assumptions -> {asym > 0, bsym > 0, csym > 0,
        xx \[Element] Reals, yy \[Element] Reals, zz \[Element] Reals,
        Abs[xx] < asym, Abs[yy] < bsym, Abs[zz] < csym}],
    180, $Failed];
  If[inner3 === $Failed,
    Print["Step 1c zpp: TIMEOUT"]; $stepStatus = "zpp_failed",
    Print["Step 1c (zpp): closed form FOUND -- B_x(r) FULLY CLOSED"]; $stepStatus = "step1_ok"]];
Print[""];

If[$stepStatus =!= "step1_ok",
  Print["B_x(r) closed form not fully obtained. Falling back to NIntegrate."];
  Print[""];
  (* Numerical evaluation only -- skip Step 2,3 closed form attempts *)
  Bxnum[xR_?NumberQ, yR_?NumberQ, zR_?NumberQ] := NIntegrate[
    xpp/Sqrt[(xR - xpp)^2 + (yR - ypp)^2 + (zR - zpp)^2],
    {xpp, -aN/2, aN/2}, {ypp, -bN/2, bN/2}, {zpp, -cN/2, cN/2},
    WorkingPrecision -> 30, AccuracyGoal -> 22, PrecisionGoal -> 20,
    MaxRecursion -> 12];

  Print["Computing alpha_3 via 9D-reduced numerical:"];
  Print["  T_xx = integral_Omega^2 [x B_x(r) integral_Omega 1/|r-r'| dV'] -- still 6D"];
  Print[""];

  (* Direct estimation: 9D Monte Carlo via NIntegrate Method->MonteCarlo *)
  Print["Attempting 9D NIntegrate Monte Carlo..."];
  tMC = AbsoluteTime[];
  TxxMC = TimeConstrained[
    NIntegrate[
      xx * xpp/(Sqrt[(xx - xp)^2 + (yy - yp)^2 + (zz - zp)^2] *
                Sqrt[(xp - xpp)^2 + (yp - ypp)^2 + (zp - zpp)^2]),
      {xx, -aN/2, aN/2}, {yy, -bN/2, bN/2}, {zz, -cN/2, cN/2},
      {xp, -aN/2, aN/2}, {yp, -bN/2, bN/2}, {zp, -cN/2, cN/2},
      {xpp, -aN/2, aN/2}, {ypp, -bN/2, bN/2}, {zpp, -cN/2, cN/2},
      Method -> "MonteCarlo", MaxPoints -> 10^6,
      WorkingPrecision -> 20],
    600, $Failed];
  If[TxxMC === $Failed,
    Print["  9D MC: TIMEOUT"],
    Print["  T_xx (9D MC, 10^6 pts) = ", N[TxxMC, 8],
      "  (", N[AbsoluteTime[] - tMC, 4], " s)"]];
  Print[""];

  (* Compute alpha_3 from MC if obtained *)
  If[TxxMC =!= $Failed,
    (* Symmetric: T_yy from a<->b swap; for our box use approximate symmetry *)
    Tratio = N[(aN/bN)^4, 6];  (* rough D_xx/D_yy ratio guess *)
    TyyMC = TxxMC/Tratio;  (* very rough *)
    alpha3Born = -(mu0N sigmaN)^3/(64 PiN^2 VN) * (TxxMC + TyyMC);
    Print["alpha_3 Born (rough MC) = ", N[alpha3Born, 8]]];

  Print["=== Step 1 not fully closed form. Need more careful Mathematica derivation. ==="];
  Quit[]];

(* If Step 1 succeeded, store B_x(r) as closed form expression *)
BxClosed = inner3;
Print["B_x(r) closed form expression size: ", LeafCount[BxClosed]];
Print[""];

(* Step 2: G_2(r) := integral_Omega B_x(r')/|r-r'| dV' *)
Print["Step 2: G_2(r) = integral B_x(rp)/|r-rp| dVp"];
Print[""];

BxRpClosed = BxClosed /. {xx -> xrp, yy -> yrp, zz -> zrp};
integrandStep2 = BxRpClosed/Sqrt[(xx - xrp)^2 + (yy - yrp)^2 + (zz - zrp)^2];

Print["Trying Integrate[] over xrp..."];
g2_xrp = TimeConstrained[
  Integrate[integrandStep2, {xrp, -asym/2, asym/2},
    Assumptions -> {asym > 0, bsym > 0, csym > 0,
      xx \[Element] Reals, yy \[Element] Reals, zz \[Element] Reals,
      yrp \[Element] Reals, zrp \[Element] Reals,
      Abs[xx] < asym/2, Abs[yy] < bsym/2, Abs[zz] < csym/2}],
  300, $Failed];

If[g2_xrp === $Failed,
  Print["  Step 2a xrp: TIMEOUT"];
  Print["  Falling back to NIntegrate for Step 2..."],
  Print["  Step 2a (xrp): FOUND"]];

(* Stop here for now; further steps require Step 2 success *)
Print[""];
Print["=== Stage 5 alpha_3 closed-form: progress checkpoint ==="];
Print["* Step 1 (B_x(r) closed form): ",
  If[$stepStatus === "step1_ok", "DONE", "PARTIAL"]];
Print["* Step 2 (G_2(r) = B_x conv G): ",
  If[g2_xrp === $Failed, "BLOCKED at xrp", "in progress"]];
Print[""];
Print["Next: continue Step 2 reduction or switch to numerical 6D"];
