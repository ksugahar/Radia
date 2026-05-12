(* Stage 5 (Born): Closed-form alpha_2 via free-space Green's function.

   Theory:
     chi(s) at low frequency, Born series in s mu_0 sigma:
       A_ind^(1)(r) = -(s sigma mu_0/4 pi) integral A_ext(r')/|r-r'| dV'
       chi_corr at order s^2: alpha_2 = (mu_0 sigma)^2/(16 pi V) (D_xx + D_yy)
       D_xx = double-integral over cuboid^2 of x x'/|r-r'|

   Reduction to 3D:
     Substitute xi = x - x', eta = y - y', zeta = z - z'.
     For polynomial xx' = X^2 - xi^2/4 with X = (x+x')/2:
       integral over X (length a - |xi|) gives
       g_x(xi) = (a-|xi|)/12 * [(a-|xi|)^2 - 3 xi^2]
     For y, z (no polynomial): g_y(eta) = b-|eta|, g_z(zeta) = c-|zeta|.

   D_xx = integral over [-a,a] x [-b,b] x [-c,c] of g_x(xi) g_y(eta) g_z(zeta)/r dV
        = 8 integral over [0,a] x [0,b] x [0,c] (using symmetry).

   This is an Aharoni-type cuboid integral. Try Mathematica symbolic Integrate[]
   first; fall back to high-precision NIntegrate.
*)
$HistoryLength = 0;
prec = 60;

Print["=== Born series alpha_2 (closed-form attempt) ==="];
Print[""];

(* Constants *)
sigmaN = SetPrecision[58000000, prec];
mu0N = SetPrecision[4 Pi/10^7, prec];
aN = SetPrecision[5/1000, prec];
bN = SetPrecision[2/1000, prec];
cN = SetPrecision[1/1000, prec];
PiN = SetPrecision[Pi, prec];
VN = aN bN cN;

(* === Derived integrand reduction === *)
(* g_x(xi) for xi in [0, a] (use symmetry to halve integration range) *)
(* g_x(xi) = (a - xi)/12 * ((a - xi)^2 - 3 xi^2) for xi >= 0 *)

Print["Step 1: try symbolic Integrate over inner xi (with x-poly)"];
Print[""];

(* Try fully symbolic Integrate *)
clearVars = ClearAll[xi, eta, zeta, asym, bsym, csym];
gxSym[xi_, asym_] := (asym - Abs[xi])/12 * ((asym - Abs[xi])^2 - 3 xi^2);
gySym[eta_, bsym_] := bsym - Abs[eta];
gzSym[zeta_, csym_] := csym - Abs[zeta];

(* Attempt 1D symbolic on the inner xi integral with fixed eta, zeta *)
inner1Try = TimeConstrained[
  Integrate[gxSym[xi, asym]/Sqrt[xi^2 + R2sq], {xi, -asym, asym},
    Assumptions -> {asym > 0, R2sq > 0}],
  60, $Failed];
Print["Symbolic xi integral closed form? ", If[inner1Try === $Failed, "TIMEOUT/FAIL", "FOUND"]];
If[inner1Try =!= $Failed,
  Print["  Result = ", inner1Try]];
Print[""];

(* If symbolic xi works, then build the 2D outer integral *)
(* Else go straight to numerical with very high precision *)

Print["Step 2: high-precision NIntegrate over the 3D reduced integral"];
Print[""];
Print["D_xx = integral_{-a}^{a} dxi integral_{-b}^{b} deta integral_{-c}^{c} dzeta ..."];
Print["     = 8 * integral_{0,a,b,c} (after symmetry)"];
Print[""];

(* Symmetric (factor 8) reduced integral *)
gxNum[xi_?NumberQ] := (aN - xi)/12 * ((aN - xi)^2 - 3 xi^2);

tStart = AbsoluteTime[];
DxxInt = 8 NIntegrate[
  gxNum[xi] (bN - eta)(cN - zeta)/Sqrt[xi^2 + eta^2 + zeta^2],
  {xi, 0, aN}, {eta, 0, bN}, {zeta, 0, cN},
  WorkingPrecision -> 50, AccuracyGoal -> 40, PrecisionGoal -> 35,
  MaxRecursion -> 14];
Print["D_xx = ", N[DxxInt, 30], "  (", N[AbsoluteTime[] - tStart, 4], " s)"];

(* Symmetry under (a, b) -> (b, a) gives D_yy *)
gyNum[eta_?NumberQ] := (bN - eta)/12 * ((bN - eta)^2 - 3 eta^2);
tStart = AbsoluteTime[];
DyyInt = 8 NIntegrate[
  (aN - xi) gyNum[eta] (cN - zeta)/Sqrt[xi^2 + eta^2 + zeta^2],
  {xi, 0, aN}, {eta, 0, bN}, {zeta, 0, cN},
  WorkingPrecision -> 50, AccuracyGoal -> 40, PrecisionGoal -> 35,
  MaxRecursion -> 14];
Print["D_yy = ", N[DyyInt, 30], "  (", N[AbsoluteTime[] - tStart, 4], " s)"];
Print[""];

(* alpha_2 from Born formula *)
alpha2Born = (mu0N sigmaN)^2/(16 PiN VN) * (DxxInt + DyyInt);
Print["alpha_2 (Born) = (mu_0 sigma)^2/(16 pi V) (D_xx + D_yy) = ",
  N[alpha2Born, 16]];
Print[""];

(* Reference values *)
ReffSq = 5 (aN^2 + bN^2)/24;
alpha2Sphere = (mu0N sigmaN ReffSq)^2/105;
alpha2Stage5c = SetPrecision[4.633960934739281399894378775*^-10, prec];

Print["=== Comparison ==="];
Print["alpha_2 Born (free-space)            = ", N[alpha2Born, 16]];
Print["alpha_2 Galerkin (Stage 5c, PEC+dip) = ", N[alpha2Stage5c, 16]];
Print["alpha_2 Sphere (R_eff)               = ", N[alpha2Sphere, 16]];
Print[""];
Print["Ratio Born/Galerkin = ", N[alpha2Born/alpha2Stage5c, 8]];
Print["Ratio Born/Sphere   = ", N[alpha2Born/alpha2Sphere, 8]];
Print[""];

If[Abs[N[alpha2Born/alpha2Sphere, 8]] > 0.5 && Abs[N[alpha2Born/alpha2Sphere, 8]] < 2.0,
  Print["==> Born is comparable to sphere reference — promising!"],
  Print["==> Born differs significantly from sphere; check formulation"]];
