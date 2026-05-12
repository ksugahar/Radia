(* Task C: Spectral BEM v4 = Legendre basis + 1D-closed + 2D-NIntegrate.

   Combines:
   - v1 closed alpha_2 strategy (1D xi integrated symbolically) -> 30 digit precision in 3 min
   - Legendre basis (mass diagonal, no Cholesky, condition number ~ 25 vs monomial 1e21)
   - All K elements via this strategy -> general spectral BEM at 30 digit precision

   Key insight: for monomial weights, the inner xi integral
     int_0^{h_x} xi^k * h_poly(xi) / sqrt(xi^2 + R^2) dxi
   has closed form (arcsinh + sqrt + powers), via Mathematica Integrate[] with
   Assumptions {h > 0, R > 0}.

   We compute this closed form symbolically once per (p1, p2), then 2D NIntegrate
   over (eta, zeta).

   Date: 2026-05-02
*)

$HistoryLength = 0;
prec = 50;

Print["================================================================"];
Print["Task C: Spectral BEM v4 (Legendre + 1D-closed + 2D-NIntegrate)"];
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

(* === Inner u-integration (xi > 0, monomial basis), Expanded === *)
ClearAll[gAxisPos];
gAxisPos[p1_Integer, p2_Integer, h_, xi_] := gAxisPos[p1, p2, h, xi] =
  Expand[Integrate[(u + xi/2)^p1 (u - xi/2)^p2,
    {u, xi/2 - h/2, h/2 - xi/2},
    Assumptions -> {xi > 0, xi < h, h > 0}]];

(* === 1D-closed inner xi integration ===
   For each (p1, p2), compute symbolically:
     IxiClosed[p1, p2, h, R] = int_0^h gAxisPos[p1, p2, h, xi] / sqrt(xi^2 + R^2) dxi
   This is a polynomial-times-1/sqrt integral, has closed form via Mathematica Integrate[].
*)
ClearAll[IxiClosed];
IxiClosed[p1_Integer, p2_Integer, h_, R_] := IxiClosed[p1, p2, h, R] = Module[
  {xi, gPoly, result},
  gPoly = gAxisPos[p1, p2, h, xi];
  result = Integrate[gPoly/Sqrt[xi^2 + R^2], {xi, 0, h},
    Assumptions -> {h > 0, R > 0}];
  Simplify[result, Assumptions -> {h > 0, R > 0}]];

(* Sanity: print first few *)
Print["Sanity checks (closed-form inner xi):"];
Print["  IxiClosed[0,0,h,R] = ", IxiClosed[0, 0, h, R]];
Print["  IxiClosed[1,1,h,R] = ", IxiClosed[1, 1, h, R]];
Print["  IxiClosed[2,2,h,R] = ", IxiClosed[2, 2, h, R] // Short[#, 8] &];
Print[""];

(* === K element via 2D NIntegrate ===
   K = 8 * int_0^b int_0^c IxiClosed[p1, p2, aN, sqrt(eta^2 + zeta^2)] *
                            gAxisPos[q1, q2, bN, eta] * gAxisPos[r1, r2, cN, zeta]
                            d(eta) d(zeta)
*)
ClearAll[KMon];
KMon[{p1_, q1_, r1_}, {p2_, q2_, r2_}] := KMon[{p1, q1, r1}, {p2, q2, r2}] = Module[
  {eta, zeta, R, ixi, gyExpr, gzExpr, integrand},
  R = Sqrt[eta^2 + zeta^2];
  ixi = IxiClosed[p1, p2, aN, R];
  gyExpr = gAxisPos[q1, q2, bN, eta];
  gzExpr = gAxisPos[r1, r2, cN, zeta];
  integrand = ixi * gyExpr * gzExpr;
  8 NIntegrate[integrand,
    {eta, 0, bN}, {zeta, 0, cN},
    WorkingPrecision -> prec, AccuracyGoal -> prec - 5,
    PrecisionGoal -> prec - 8, MaxRecursion -> 16,
    Method -> "GlobalAdaptive"]];

(* Quick test: K[(0,1,0),(0,1,0)] = D_yy (single matrix element) *)
Print["Test: K[(0,1,0),(0,1,0)] = D_yy (single element)"];
t0 = AbsoluteTime[];
DyyTest = KMon[{0, 1, 0}, {0, 1, 0}];
Print["  D_yy (Task C 1D-closed + 2D-NIntegrate) = ", N[DyyTest, 30]];
Print["  Time: ", N[AbsoluteTime[] - t0, 4], " s"];

(* The reference D_yy from v1 closed alpha_2 wls: 5.249182700352793100999037232e-21 *)
DyyRef = SetPrecision[5.24918270035279310099903723219621782347721846299*^-21, 30];
Print["  Reference D_yy (v1 closed at WP=60) = ", N[DyyRef, 30]];
Print["  Relative diff = ", N[Abs[(DyyTest - DyyRef)/DyyRef], 6]];
Print[""];

(* Quick test: K[(1,0,0),(1,0,0)] = D_xx *)
Print["Test: K[(1,0,0),(1,0,0)] = D_xx"];
t0 = AbsoluteTime[];
DxxTest = KMon[{1, 0, 0}, {1, 0, 0}];
Print["  D_xx = ", N[DxxTest, 30]];
Print["  Time: ", N[AbsoluteTime[] - t0, 4], " s"];
DxxRef = SetPrecision[6.55125437660237267362579845589291716611473968610625*^-20, 30];
Print["  Reference D_xx (v1 closed) = ", N[DxxRef, 30]];
Print["  Relative diff = ", N[Abs[(DxxTest - DxxRef)/DxxRef], 6]];
Print[""];

(* Compute alpha_2 from these *)
alpha2Test = (mu0sigmaN)^2/(16 PiN VN) (DxxTest + DyyTest);
Print["alpha_2 (Task C, single-element) = ", N[alpha2Test, 30]];
Print["  vs alpha_2 closed (NEW) = ", N[alpha2RefNew, 30]];
Print["  Relative diff = ", N[Abs[(alpha2Test - alpha2RefNew)/alpha2RefNew], 6]];
