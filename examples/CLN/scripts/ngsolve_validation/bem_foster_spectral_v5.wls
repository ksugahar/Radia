(* Spectral BEM v5: Legendre basis + 1D-closed + 2D-NIntegrate.

   Combines best of all approaches:
   - Legendre basis: mass diagonal, condition ~ 25 (vs monomial 10^21)
   - 1D xi closed-form: arcsinh + sqrt + powers, exact integration
   - 2D NIntegrate at WP=50: smooth integrand, ~50 sec/element, 17-18 digit precision

   Performance baseline (Task C verification):
     D_yy single element: 17 digit agreement with v1 closed at 50 sec
     alpha_2 single-element: 1.74e-18 relative error vs new closed form

   This script: full K matrix at P=6 -> alpha_n to 8 orders -> Cauer extraction.

   Date: 2026-05-02
*)

$HistoryLength = 0;
prec = 50;

Print["================================================================"];
Print["Spectral BEM v5: Legendre + 1D-closed + 2D-NIntegrate"];
Print["WorkingPrecision = ", prec];
Print["================================================================"];
Print[""];

sigmaN = SetPrecision[58 10^6, prec];
mu0N = SetPrecision[4 Pi/10^7, prec];
PiN = SetPrecision[Pi, prec];
aN = SetPrecision[5/1000, prec];
bN = SetPrecision[2/1000, prec];
cN = SetPrecision[1/1000, prec];
VN = aN bN cN;
mu0sigmaN = mu0N sigmaN;

alpha1Ref = -mu0sigmaN (aN^2 + bN^2)/48;
alpha2RefNew = SetPrecision[7.4783240535434074691212813799*^-10, 28];
tau1Ref = alpha2RefNew/(-alpha1Ref);

Print["References (NEW closed form):"];
Print["  alpha_1 = ", N[alpha1Ref, 18]];
Print["  alpha_2 = ", N[alpha2RefNew, 18]];
Print["  tau_1   = ", N[tau1Ref 10^6, 18], " us"];
Print[""];

(* === Inner u-integration symbolic (cached) === *)
ClearAll[gAxisPos];
gAxisPos[p1_Integer, p2_Integer, h_, xi_] := gAxisPos[p1, p2, h, xi] =
  Expand[Integrate[(u + xi/2)^p1 (u - xi/2)^p2,
    {u, xi/2 - h/2, h/2 - xi/2},
    Assumptions -> {xi > 0, xi < h, h > 0}]];

(* === 1D-closed xi integration (cached) === *)
ClearAll[IxiClosed];
IxiClosed[p1_Integer, p2_Integer, h_, R_] := IxiClosed[p1, p2, h, R] = Module[
  {xi, gPoly, result},
  gPoly = gAxisPos[p1, p2, h, xi];
  result = Integrate[gPoly/Sqrt[xi^2 + R^2], {xi, 0, h},
    Assumptions -> {h > 0, R > 0}];
  Simplify[result, Assumptions -> {h > 0, R > 0}]];

(* === K element (monomial) via 2D NIntegrate === *)
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

(* === Legendre coefficients (cached) === *)
ClearAll[legendreT];
legendreT[Pmax_, h_] := legendreT[Pmax, h] = Module[{T, n, k, c},
  T = ConstantArray[0, {Pmax + 1, Pmax + 1}];
  Do[
    c = CoefficientList[LegendreP[n, x], x];
    Do[T[[n + 1, k + 1]] = c[[k + 1]] (2/h)^k, {k, 0, n}],
    {n, 0, Pmax}];
  T];

(* === Build basis (e,o,e) === *)
buildBasis[Pmax_, parity_] := Module[{px, py, pz, basis},
  {px, py, pz} = parity;
  basis = {};
  Do[
    If[Mod[p, 2] == px && Mod[q, 2] == py && Mod[r, 2] == pz && p + q + r <= Pmax,
      AppendTo[basis, {p, q, r}]],
    {p, 0, Pmax}, {q, 0, Pmax - p}, {r, 0, Pmax - p - q}];
  basis];

(* === Run for given Pmax (one parity block at a time, then both) === *)
runBlock[Pmax_, parity_, label_] := Module[
  {basis, n, t0, KMonMat, Tx, Ty, Tz, TblockL, KLeg, massDiag,
   KLegStandard, eigvals, eigvecs},
  basis = buildBasis[Pmax, parity];
  n = Length[basis];
  Print["  Block ", label, " at P=", Pmax, ": ", n, " functions"];

  Print["    Building K_mon matrix (", n (n + 1)/2, " unique)..."];
  t0 = AbsoluteTime[];
  KMonMat = ConstantArray[0, {n, n}];
  Do[
    Module[{val = KMon[basis[[i]], basis[[j]]]},
      KMonMat[[i, j]] = val;
      KMonMat[[j, i]] = val],
    {i, n}, {j, i, n}];
  Print["      K_mon time: ", N[AbsoluteTime[] - t0, 4], " s"];

  Tx = legendreT[Pmax, aN];
  Ty = legendreT[Pmax, bN];
  Tz = legendreT[Pmax, cN];
  TblockL = Table[
    Module[{P, Q, R, p, q, r},
      {P, Q, R} = basis[[i]];
      {p, q, r} = basis[[j]];
      Tx[[P + 1, p + 1]] Ty[[Q + 1, q + 1]] Tz[[R + 1, r + 1]]],
    {i, n}, {j, n}];

  KLeg = TblockL.KMonMat.Transpose[TblockL];

  massDiag = Table[
    Module[{P, Q, R}, {P, Q, R} = basis[[i]];
      aN/(2 P + 1) bN/(2 Q + 1) cN/(2 R + 1)],
    {i, n}];

  KLegStandard = Table[KLeg[[i, j]]/Sqrt[massDiag[[i]] massDiag[[j]]],
    {i, n}, {j, n}];

  Print["    Solving standard eigenvalue..."];
  t0 = AbsoluteTime[];
  {eigvals, eigvecs} = Eigensystem[N[KLegStandard, prec]];
  Print["      eigvalue solve: ", N[AbsoluteTime[] - t0, 4], " s"];

  {basis, eigvals, eigvecs, massDiag}];

(* Project A_ext on M-orthonormal Legendre basis *)
projectAext[basis_, massDiag_, axis_] := Module[{vec, idx, coef},
  vec = ConstantArray[0, Length[basis]];
  Switch[axis,
    "x",
      idx = Position[basis, {0, 1, 0}, 1];
      If[Length[idx] > 0,
        coef = -bN/4;
        vec[[idx[[1, 1]]]] = coef Sqrt[massDiag[[idx[[1, 1]]]]]],
    "y",
      idx = Position[basis, {1, 0, 0}, 1];
      If[Length[idx] > 0,
        coef = aN/4;
        vec[[idx[[1, 1]]]] = coef Sqrt[massDiag[[idx[[1, 1]]]]]]];
  vec];

runForPmax[Pmax_] := Module[
  {bX, eX, vX, mX, bY, eY, vY, mY, aExtX, aExtY, gXk, gYk,
   lambdaX, lambdaY, alphas, NMoments, tau1BEM},
  Print["=== P_max = ", Pmax, " ==="];

  {bX, eX, vX, mX} = runBlock[Pmax, {0, 1, 0}, "x-coupling (e,o,e)"];
  {bY, eY, vY, mY} = runBlock[Pmax, {1, 0, 0}, "y-coupling (o,e,e)"];

  aExtX = projectAext[bX, mX, "x"];
  aExtY = projectAext[bY, mY, "y"];

  gXk = vX.aExtX;
  gYk = vY.aExtY;

  lambdaX = eX/(4 PiN);
  lambdaY = eY/(4 PiN);

  NMoments = 8;
  alphas = Table[
    (-1)^n (mu0sigmaN)^n/VN *
    (Total[gXk^2 lambdaX^(n - 1)] + Total[gYk^2 lambdaY^(n - 1)]),
    {n, 1, NMoments}];

  Print[""];
  Print["  alpha_n vs reference (NEW closed):"];
  Print["    alpha_1 = ", N[alphas[[1]], 18],
        "  err = ", N[Abs[(alphas[[1]] - alpha1Ref)/alpha1Ref], 6]];
  Print["    alpha_2 = ", N[alphas[[2]], 18],
        "  err = ", N[Abs[(alphas[[2]] - alpha2RefNew)/alpha2RefNew], 6]];
  Do[Print["    alpha_", n, " = ", N[alphas[[n]], 18]], {n, 3, NMoments}];

  tau1BEM = alphas[[2]]/(-alphas[[1]]);
  Print["    tau_1   = ", N[tau1BEM 10^6, 18], " us"];
  Print["    tau_1 err = ", N[Abs[(tau1BEM - tau1Ref)/tau1Ref], 6]];
  Print[""];
  alphas];

(* === Run P=4 first (smaller, less memory risk) === *)
result4 = runForPmax[4];
Print["----------------------------------------------------------------"];
