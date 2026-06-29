(* Mathematica K matrix exporter for Python.

   Computes K_mon[(p1,q1,r1),(p2,q2,r2)] for spectral basis at WP=50,
   stores as JSON for Python to load.

   This is the AUTHORITATIVE K matrix calculation, since Mathematica's
   NIntegrate GlobalAdaptive handles the corner log-singularity correctly
   (mpmath.quad limited to ~15 digits).

   Output: bem_foster_K_export_P{P_max}.json with K elements at 30 digit precision.

   Date: 2026-05-02
*)

$HistoryLength = 0;
prec = 50;
PmaxArg = 4;  (* default; set via wolframscript -file ... -PmaxArg N *)
If[Length[$ScriptCommandLine] >= 2,
  PmaxArg = ToExpression[$ScriptCommandLine[[2]]]];

Print["================================================================"];
Print["Mathematica K matrix exporter (WP=", prec, ", P=", PmaxArg, ")"];
Print["================================================================"];

aN = SetPrecision[5/1000, prec];
bN = SetPrecision[2/1000, prec];
cN = SetPrecision[1/1000, prec];

(* Inner u-integration *)
ClearAll[gAxisPos];
gAxisPos[p1_Integer, p2_Integer, h_, xi_] := gAxisPos[p1, p2, h, xi] =
  Expand[Integrate[(u + xi/2)^p1 (u - xi/2)^p2,
    {u, xi/2 - h/2, h/2 - xi/2},
    Assumptions -> {xi > 0, xi < h, h > 0}]];

(* 1D-closed inner xi integration *)
ClearAll[IxiClosed];
IxiClosed[p1_Integer, p2_Integer, h_, R_] := IxiClosed[p1, p2, h, R] = Module[
  {xi, gPoly, result},
  gPoly = gAxisPos[p1, p2, h, xi];
  result = Integrate[gPoly/Sqrt[xi^2 + R^2], {xi, 0, h},
    Assumptions -> {h > 0, R > 0}];
  Simplify[result, Assumptions -> {h > 0, R > 0}]];

(* K element via 2D NIntegrate *)
ClearAll[KMon];
KMon[{p1_, q1_, r1_}, {p2_, q2_, r2_}] := KMon[{p1, q1, r1}, {p2, q2, r2}] = Module[
  {eta, zeta, R, ixi, gyExpr, gzExpr, integrand, val},
  R = Sqrt[eta^2 + zeta^2];
  ixi = IxiClosed[p1, p2, aN, R];
  gyExpr = gAxisPos[q1, q2, bN, eta];
  gzExpr = gAxisPos[r1, r2, cN, zeta];
  integrand = ixi * gyExpr * gzExpr;
  val = 8 NIntegrate[integrand,
    {eta, 0, bN}, {zeta, 0, cN},
    WorkingPrecision -> prec, AccuracyGoal -> prec - 5,
    PrecisionGoal -> prec - 8, MaxRecursion -> 16,
    Method -> "GlobalAdaptive"];
  val];

(* Build basis (e,o,e) and (o,e,e) parities *)
buildBasis[Pmax_, parity_] := Module[{px, py, pz, basis},
  {px, py, pz} = parity;
  basis = {};
  Do[
    If[Mod[p, 2] == px && Mod[q, 2] == py && Mod[r, 2] == pz && p + q + r <= Pmax,
      AppendTo[basis, {p, q, r}]],
    {p, 0, Pmax}, {q, 0, Pmax - p}, {r, 0, Pmax - p - q}];
  basis];

basisX = buildBasis[PmaxArg, {0, 1, 0}];
basisY = buildBasis[PmaxArg, {1, 0, 0}];
nX = Length[basisX]; nY = Length[basisY];
Print["x-block (e,o,e): ", nX, " functions"];
Print["y-block (o,e,e): ", nY, " functions"];

(* Compute K matrices *)
Print[""];
Print["Building K_mon for x-block..."];
t0 = AbsoluteTime[];
KX = ConstantArray[0, {nX, nX}];
Do[
  Module[{val = KMon[basisX[[i]], basisX[[j]]]},
    KX[[i, j]] = val;
    KX[[j, i]] = val],
  {i, nX}, {j, i, nX}];
Print["  done: ", N[AbsoluteTime[] - t0, 4], " s"];

Print["Building K_mon for y-block..."];
t0 = AbsoluteTime[];
KY = ConstantArray[0, {nY, nY}];
Do[
  Module[{val = KMon[basisY[[i]], basisY[[j]]]},
    KY[[i, j]] = val;
    KY[[j, i]] = val],
  {i, nY}, {j, i, nY}];
Print["  done: ", N[AbsoluteTime[] - t0, 4], " s"];

(* Export to JSON with high-precision string format *)
KEntryToJSON[i_, j_, val_] :=
  <|"i" -> i - 1, "j" -> j - 1, "value" -> ToString[CForm[N[val, 30]]]|>;

KMatrixToJSON[Kmat_, n_] :=
  Flatten[Table[KEntryToJSON[i, j, Kmat[[i, j]]], {i, n}, {j, n}]];

filename = "bem_foster_K_export_P" <> ToString[PmaxArg] <> ".json";

Export[filename,
  <|"P_max" -> PmaxArg,
    "WP" -> prec,
    "a" -> ToString[N[aN, 30], InputForm],
    "b" -> ToString[N[bN, 30], InputForm],
    "c" -> ToString[N[cN, 30], InputForm],
    "x_block" -> <|
      "basis" -> basisX,
      "K_elements" -> KMatrixToJSON[KX, nX]
    |>,
    "y_block" -> <|
      "basis" -> basisY,
      "K_elements" -> KMatrixToJSON[KY, nY]
    |>
  |>];

Print[""];
Print["Saved K matrices to: ", filename];
Print["  x-block: ", nX, "x", nX, " = ", nX*nX, " entries"];
Print["  y-block: ", nY, "x", nY, " = ", nY*nY, " entries"];
