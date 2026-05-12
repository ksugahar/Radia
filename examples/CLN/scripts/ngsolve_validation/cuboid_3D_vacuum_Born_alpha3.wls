(* Stage 5 (Born alpha_3): 9D->6D reduction via closed-form A_ind^(1).

   alpha_3 = -(mu_0 sigma)^3/(64 pi^2 V) * (T_xx + T_yy)
   T_xx = integral over Omega^3 of xx''/(|r-r'| |r'-r''|) dV dV' dV''

   Using A_ind^(1)(r) = -(s mu_0 sigma/4 pi) F(r) where:
     F_x(r) = integral_Omega (-y'')/|r-r''| dV''     (free-space vector potential
     F_y(r) = integral_Omega x''/|r-r''| dV''         from impressed A_ext source)
   These F_x, F_y have CLOSED FORM (Aharoni-type), already derivable.

   Then:
     T_xx = integral_Omega^2 x x'' / (|r-r'||r'-r''|) ...
          = -integral_Omega A_ext_y(r') × A_ind_y^(1)(r')... (with right substitution)

   Cleaner expression: alpha_3 directly via
     alpha_3 = -(mu_0 sigma)^3/(V) * (1/(4 pi)^2) *
                integral_Omega^2 A_ext(r) . A_ext(r')... wait no

   Let me re-derive. Define:
     B_corr(r) := integral_Omega A_ext(r')/|r-r'| dV'    [vector function of r]
   This is closed-form integrable.

   Then A_ind^(1)(r) = -(s mu_0 sigma/4 pi) B_corr(r).
   <A_ext, A_ind^(1)> = -(s mu_0 sigma/4 pi) integral A_ext . B_corr dV
                      = -(s mu_0 sigma/4 pi) (D_xx + D_yy)/4
   ... so 16 pi V denominator for alpha_2.

   For alpha_3:
     A_ind^(2)(r) = -(s mu_0 sigma/4 pi) integral_Omega A_ind^(1)(r')/|r-r'| dV'
                  = (s mu_0 sigma/4 pi)^2 integral_Omega B_corr(r')/|r-r'| dV'
                  =: (s mu_0 sigma/4 pi)^2 C_corr(r)

   where C_corr(r) := integral_Omega B_corr(r')/|r-r'| dV'.

   <A_ext, A_ind^(2)> = (s mu_0 sigma/4 pi)^2 integral A_ext . C_corr dV
   alpha_3 = -mu_0 sigma/V * coefficient at order s^2 of <A_ext, A_ind^(2)>/s^2
           = -(mu_0 sigma/V) (mu_0 sigma/4 pi)^2 integral A_ext . C_corr dV
           = -(mu_0 sigma)^3/(16 pi^2 V) integral_Omega A_ext . C_corr dV

   Hmm coefficient differs from earlier. Let me trust this derivation:
     alpha_3 = -(mu_0 sigma)^3/(16 pi^2 V) * E_3
     E_3 := integral_Omega A_ext(r) . C_corr(r) dV
          = integral_Omega^2 A_ext(r).B_corr(r')/|r-r'| dV dV'
          = integral_Omega^3 A_ext(r) . A_ext(r'')/(|r-r'| |r'-r''|) dV dV' dV''
          = (1/4)(T_xx + T_yy)

   So alpha_3 = -(mu_0 sigma)^3/(64 pi^2 V) (T_xx + T_yy)  -- consistent.
*)

$HistoryLength = 0;
prec = 60;

Print["=== Stage 5 Born alpha_3 (6D reduction via closed-form A_ind^(1)) ==="];
Print[""];

(* Constants *)
sigmaN = SetPrecision[58000000, prec];
mu0N = SetPrecision[4 Pi/10^7, prec];
aN = SetPrecision[5/1000, prec];
bN = SetPrecision[2/1000, prec];
cN = SetPrecision[1/1000, prec];
PiN = SetPrecision[Pi, prec];
VN = aN bN cN;

(* Step 1: Define closed-form B_corr_y(r) = integral_Omega x''/|r-r''| dV''.
   This is "vector potential at r due to uniform x-current density in cuboid". *)
ClearAll[xx, yy, zz, xpp, ypp, zpp, asym, bsym, csym];

(* For our cuboid centered at origin: integration over [-a/2, a/2] etc. *)
(* B_corr_y(r) = integral_{-a/2}^{a/2} dx'' integral... x''/sqrt((x-x'')^2 + ...) dx'' dy'' dz''

   Symbolic Mathematica try: *)
Print["Step 1: try symbolic closed form for inner integral"];
Print["  J_x(r; a, b, c) := integral_{-a/2}^{a/2}_x x''/sqrt((x-x'')^2 + (y-y'')^2 + (z-z'')^2) dx'' dy'' dz''"];

(* Inner over x'' with polynomial x'' *)
innerXpp = TimeConstrained[
  Integrate[xpp/Sqrt[(xx - xpp)^2 + RxOpp2], {xpp, -asym/2, asym/2},
    Assumptions -> {asym > 0, RxOpp2 > 0, xx \[Element] Reals}],
  60, $Failed];
If[innerXpp === $Failed,
  Print["  inner xpp integral: TIMEOUT"]; Quit[],
  Print["  inner xpp integral: FOUND"]];
Print[""];

(* Now integrate over y'' and z'' (no polynomial factors there) *)
(* This is similar to D_xx structure but with different polynomial *)
(* Skip full closed form derivation -- use NIntegrate at high precision *)

(* === Step 2: numerical 9D-reduced integral === *)
(* T_xx = integral_Omega^3 xx''/(|r-r'| |r'-r''|) dV dV' dV''
   Reduce by symmetry: the integrand is even in (y, z, y', z', y'', z'') under
   appropriate sign flips.

   For numerical evaluation, write as nested 3D integrals:
     T_xx = integral_Omega dr * x * integral_Omega dr' /|r-r'| *
            integral_Omega dr'' x'' / |r'-r''|
          = integral_Omega^2 dr dr' (x/|r-r'|) * J_x(r')
     where J_x(r') = integral_Omega x''/|r'-r''| dV''

   J_x is our F_y component (closed form attainable via Aharoni-type).

   For a quick numerical approximation, use 6D NIntegrate with cleverly nested order.
*)

Print["Step 2: numerical 9D-reduced computation"];
Print["  Method A: 9D Monte-Carlo (slow, low precision)"];
Print["  Method B: 6D nested adaptive (better but still hard)"];
Print["  Method C: pre-tabulate J_x(r') closed form, then 6D outer adaptive"];
Print[""];

(* Pre-compute J_x as a function via NIntegrate *)
(* J_x(xR_, yR_, zR_) = integral_{cuboid} x''/sqrt((xR-x'')^2 + (yR-y'')^2 + (zR-z'')^2) dV'' *)

(* Use closed form for inner (x'') integral if possible, then 2D NIntegrate *)
(* Inner integral closed form is well-known: *)
Jpp[xR_, yR_, zR_, asym_, bsym_, csym_] := Module[{xa = asym/2, xb = -asym/2, ya = bsym/2, yb = -bsym/2, za = csym/2, zb = -csym/2},
  (* integral xpp / sqrt((xR-xpp)^2 + (yR-ypp)^2 + (zR-zpp)^2) dxpp dypp dzpp *)
  (* Closed form too complex; use NIntegrate *)
  NIntegrate[xpp/Sqrt[(xR - xpp)^2 + (yR - ypp)^2 + (zR - zpp)^2],
    {xpp, xb, xa}, {ypp, yb, ya}, {zpp, zb, za},
    WorkingPrecision -> 30, AccuracyGoal -> 20, PrecisionGoal -> 18,
    MaxRecursion -> 12]];

(* Sample evaluation *)
Print["Sample J_x(r=origin) = ", N[Jpp[0, 0, 0, aN, bN, cN], 12]];
Print[""];

(* Outer 6D: alpha_3 from <A_ext, C_corr> = integral A_ext(r).C_corr(r) dV *)
(* C_corr(r) = integral B_corr(r')/|r-r'| dV', where B_corr is itself a 3D integral *)
(* This is essentially 9D nested. Computing each Jpp at each (r, r') pair is too slow. *)

Print["Direct 9D NIntegrate would take too long. Skipping."];
Print["Need structured approach: closed-form B_corr -> 6D NIntegrate with tabulation."];
Print[""];

(* PRACTICAL: just compute alpha_3 directly via 9D Monte-Carlo for verification *)
(* Skip for now -- too slow without proper implementation *)

Print["=== Skipping numerical alpha_3 here. ==="];
Print[""];

(* Fall back to using Pade [3,3] from the moments alpha_1, alpha_2_Born, alpha_3, alpha_4 *)
(* ... but we don't have alpha_3, alpha_4 from Born series yet. *)

Print["=== Alternative: direct 6D NIntegrate of T_xx with separated structure ==="];
Print[""];

(* T_xx = integral_Omega^3 xx''/(|r-r'| |r'-r''|) dV dV' dV''
   Substitute (xi, eta, zeta) = r-r' -- then r' integration is over Omega.
   r = r' + (xi, eta, zeta) must be in Omega -- complex bound shape.

   Skip the analytical reduction. Just declare scope:
   alpha_3 closed form: requires careful Mathematica work, not 1-shot.

   For now, document the structure and commit to future work.
*)

Print["=== Structural observations ==="];
Print[""];
Print["alpha_3 Born formula: -(mu_0 sigma)^3/(64 pi^2 V) (T_xx + T_yy)"];
Print["  T_xx = integral Omega^3 xx''/(|r-r'||r'-r''|) dV dV' dV''"];
Print[""];
Print["This 9D integral requires either:"];
Print["  (a) Closed form via iterated Aharoni reduction (potentially possible)"];
Print["  (b) 6D NIntegrate with closed-form inner (3D) Aharoni potential"];
Print["  (c) Direct 9D adaptive (too slow / impossible)"];
Print[""];
Print["Recommended: implement (b) -- compute B_corr(r) symbolically"];
Print["  (using same techniques as D_xx: inner xi-integral closed form +"];
Print["   2D outer NIntegrate, all at WP=60)"];
Print[""];
Print["Estimated implementation: 1-2 days of careful Mathematica work."];
Print[""];

Print["=== Stage 5 alpha_3 summary ==="];
Print["* alpha_2 Born: closed form COMPLETE (7.4783...e-10, 30 digits) ✓"];
Print["* alpha_3 Born: 9D->6D structure identified; full derivation deferred"];
Print["* For paper: have alpha_1 closed, alpha_2 closed, alpha_3+ via Pade extrapolation"];
