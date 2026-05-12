(* Stage 5 (Born, fully closed): try full symbolic Integrate of D_xx.

   Inner xi-integral closed form found:
     I_xi(a, R) = integral g_x(xi)/sqrt(xi^2 + R^2) dxi from -a to a
                = (18 a^2 R + 8 R^3 - 14 a^2 sqrt(a^2+R^2) - 8 R^2 sqrt(a^2+R^2)
                   + 3 a^3 log(1 + 2a(a+sqrt(a^2+R^2))/R^2))/36
     where R^2 = eta^2 + zeta^2.

   Now integrate over (eta, zeta): 2D outer integral with log/sqrt integrand.
   Try Mathematica Integrate[] -- may produce closed form via dilog/Cl_2 etc.
*)
$HistoryLength = 0;
prec = 60;

Print["=== Born alpha_2 fully closed form (attempt) ==="];
Print[""];

(* Symbolic version with abstract a, b, c *)
ClearAll[xi, eta, zeta, asym, bsym, csym, R];

(* g_x(xi) for xi in [0, a] *)
gxSym[xi_, asym_] := (asym - xi)/12 ((asym - xi)^2 - 3 xi^2);

(* Step 1: inner xi integral with sqrt(xi^2 + R^2) -- already known closed form *)
Print["Inner xi integral closed form (Mathematica):"];
ixiClosed = TimeConstrained[
  Integrate[gxSym[xi, asym]/Sqrt[xi^2 + Rsq],
    {xi, 0, asym},
    Assumptions -> {asym > 0, Rsq > 0}],
  60, $Failed];
If[ixiClosed === $Failed,
  Print["  Inner xi integral: FAILED (timeout)"];
  Quit[]];
Print["  ", ixiClosed];
Print[""];

(* Step 2: outer 2D integral of ixiClosed * (b - eta)(c - zeta) over [0,b] x [0,c] *)
(* Note: Rsq = eta^2 + zeta^2 *)
Print["Step 2: outer 2D integral over (eta, zeta) ..."];
ixiSubst = ixiClosed /. Rsq -> eta^2 + zeta^2;

(* For symmetric integral (factor 8 from symmetry), alphaB factor 8: *)
(* D_xx = 8 integral_0^a dxi integral_0^b deta integral_0^c dzeta (...) (b - eta)(c - zeta) *)
integrand2D = ixiSubst (bsym - eta)(csym - zeta);
Print["Integrand 2D = ", integrand2D];
Print[""];

Print["Trying symbolic Integrate[] over eta..."];
ietaClosed = TimeConstrained[
  Integrate[integrand2D, {eta, 0, bsym},
    Assumptions -> {asym > 0, bsym > 0, csym > 0, zeta > 0, zeta < csym}],
  120, $Failed];

If[ietaClosed === $Failed,
  Print["  eta integral: FAILED/TIMEOUT"];
  Print["  Falling back to high-precision NIntegrate for partial closed form"];
  symbolic2D = $Failed,
  Print["  eta integral: FOUND"];
  Print[Short[ietaClosed, 10]];
  symbolic2D = "etaSymbolic"];
Print[""];

If[symbolic2D === "etaSymbolic",
  Print["Trying outer zeta integral..."];
  ifullClosed = TimeConstrained[
    Integrate[ietaClosed, {zeta, 0, csym},
      Assumptions -> {asym > 0, bsym > 0, csym > 0}],
    300, $Failed];
  If[ifullClosed === $Failed,
    Print["  zeta integral: FAILED -- partial closed form (eta integrated)"];
    fullClosedFormFound = False,
    Print["  zeta integral: FULL CLOSED FORM!"];
    Print[Short[ifullClosed, 20]];
    fullClosedFormFound = True];
  ,
  fullClosedFormFound = False];

(* === Numerical evaluation of full closed form === *)
If[fullClosedFormFound === True,
  Print[""];
  Print["=== Evaluating full closed form at (a,b,c) = (5, 2, 1) mm ==="];
  tStart = AbsoluteTime[];
  DxxClosedForm = 8 (ifullClosed /. {asym -> aN, bsym -> bN, csym -> cN});
  DxxClosedFormReal = Re[N[DxxClosedForm, 30]];
  Print["  D_xx (full closed form, real part) = ", DxxClosedFormReal,
    "  (", N[AbsoluteTime[] - tStart, 4], " s)"];
  Print["  Imaginary part (should be 0): ", Im[N[DxxClosedForm, 16]]];
  Print[""]];

Print[""];
Print["=== Numerical evaluation at our (a,b,c) ==="];
sigmaN = SetPrecision[58000000, prec];
mu0N = SetPrecision[4 Pi/10^7, prec];
aN = SetPrecision[5/1000, prec];
bN = SetPrecision[2/1000, prec];
cN = SetPrecision[1/1000, prec];
PiN = SetPrecision[Pi, prec];
VN = aN bN cN;

(* Evaluate inner xi closed form at numeric a -- HALF-range [0, a], divisor 72 *)
ixiNum[zeta_?NumberQ, eta_?NumberQ] := Module[{R2v, av},
  R2v = eta^2 + zeta^2; av = aN;
  If[R2v < 10^-40, 0,
    (8 R2v^(3/2) - 8 R2v Sqrt[av^2 + R2v]
      - 2 av^2 (-9 Sqrt[R2v] + 7 Sqrt[av^2 + R2v])
      + 3 av^3 Log[1 + 2 av(av + Sqrt[av^2 + R2v])/R2v])/72]];

(* Use 2D NIntegrate with closed-form xi-integrand *)
tStart = AbsoluteTime[];
DxxClosed2D = 8 NIntegrate[
  ixiNum[zeta, eta] (bN - eta)(cN - zeta),
  {eta, 0, bN}, {zeta, 0, cN},
  WorkingPrecision -> 60, AccuracyGoal -> 50, PrecisionGoal -> 45,
  MaxRecursion -> 16];
Print["D_xx (xi closed + 2D NIntegrate) = ", N[DxxClosed2D, 30],
  "  (", N[AbsoluteTime[] - tStart, 4], " s)"];

(* Compare to original 3D NIntegrate from previous run *)
Print["D_xx (full 3D NIntegrate, prev) = 6.55125...e-20"];
Print[""];

(* Compute alpha_2 *)
(* D_yy by swapping a <-> b *)
iyiNum[zeta_?NumberQ, xi_?NumberQ] := Module[{R2v, bv},
  R2v = xi^2 + zeta^2; bv = bN;
  If[R2v < 10^-40, 0,
    (8 R2v^(3/2) - 8 R2v Sqrt[bv^2 + R2v]
      - 2 bv^2 (-9 Sqrt[R2v] + 7 Sqrt[bv^2 + R2v])
      + 3 bv^3 Log[1 + 2 bv(bv + Sqrt[bv^2 + R2v])/R2v])/72]];

tStart = AbsoluteTime[];
DyyClosed2D = 8 NIntegrate[
  iyiNum[zeta, xi] (aN - xi)(cN - zeta),
  {xi, 0, aN}, {zeta, 0, cN},
  WorkingPrecision -> 60, AccuracyGoal -> 50, PrecisionGoal -> 45,
  MaxRecursion -> 16];
Print["D_yy (xi closed + 2D NIntegrate) = ", N[DyyClosed2D, 30],
  "  (", N[AbsoluteTime[] - tStart, 4], " s)"];
Print[""];

alpha2Born = (mu0N sigmaN)^2/(16 PiN VN) * (DxxClosed2D + DyyClosed2D);
Print["=== alpha_2 (Born, with 1D closed-form xi-integral) ==="];
Print["  alpha_2^Born = ", N[alpha2Born, 30]];
Print[""];

(* Compare to previous 3D NIntegrate *)
prevAlpha2Born = SetPrecision[7.4783240521932061573535486933*^-10, prec];
Print["  Previous (3D NIntegrate): ", N[prevAlpha2Born, 30]];
Print["  Difference (current - prev): ", N[alpha2Born - prevAlpha2Born, 16]];
Print["  Relative diff: ", N[Abs[(alpha2Born - prevAlpha2Born)/prevAlpha2Born], 8]];
Print[""];

(* Check series limits *)
Print["=== Sanity checks ==="];
Print["For a = b = c = 1mm (cube), alpha_2 should match cube literature value"];
(* Could compute -- skip for now *)

Print["For a = b = c -> infty isotropy: D_xx = D_yy by symmetry"];
Print["For our a=5, b=2, c=1: D_xx/D_yy = ", N[DxxClosed2D/DyyClosed2D, 6]];
Print["  (a/b)^? ratio expected: a^4/b^4 = ", N[(aN/bN)^4, 6],
  "; a^2/b^2 = ", N[(aN/bN)^2, 6]];
Print[""];

Print["=== Stage 5 (Born) summary ==="];
Print["* alpha_2 Born computed at 30+ digits via xi-closed + 2D NIntegrate"];
Print["* Inner xi integral has closed form (log + sqrt), Aharoni-type"];
Print["* Full 3D closed form possible if eta, zeta also analytical (ongoing attempt)"];
