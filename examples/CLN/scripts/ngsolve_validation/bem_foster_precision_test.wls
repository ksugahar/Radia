(* Verify Mathematica K element precision at WP=50 vs WP=70.
   If results agree to 50+ digits, WP=50 NIntegrate truly delivers 18+ digit
   accuracy. If they differ at digit 18, WP=50 is overstating precision.
*)

$HistoryLength = 0;

testKelement[prec_] := Module[
  {sigmaN, mu0N, PiN, aN, bN, cN, gx, gy, gz, xi, eta, zeta,
   ixi, integrand, val, t0},
  sigmaN = SetPrecision[58 10^6, prec];
  mu0N = SetPrecision[4 Pi/10^7, prec];
  PiN = SetPrecision[Pi, prec];
  aN = SetPrecision[5/1000, prec];
  bN = SetPrecision[2/1000, prec];
  cN = SetPrecision[1/1000, prec];

  (* gAxis_(1,1)(h, xi) for x-direction (used for D_xx-like K element) *)
  (* gAxis_(0,0)(h, xi) = h - xi for y, z *)

  (* K[(0,1,0), (0,1,0)] = D_yy: 8 int gAxis_00(a, xi) gAxis_11(b, eta) gAxis_00(c, zeta) / r *)

  (* Inner xi closed form for gAxis_00 (= h - xi) on [0, h]:
     I = int_0^h (h - xi)/sqrt(xi^2+R^2) dxi
       = h*arcsinh(h/R) - sqrt(h^2+R^2) + R *)

  Print["== WP = ", prec, " =="];
  t0 = AbsoluteTime[];
  ixi[R_] := aN ArcSinh[aN/R] - Sqrt[aN^2 + R^2] + R;
  gy[eta_] := (bN^3 - 3 bN^2 eta + 2 eta^3)/12;
  gz[zeta_] := cN - zeta;

  val = 8 NIntegrate[
    ixi[Sqrt[eta^2 + zeta^2]] gy[eta] gz[zeta],
    {eta, 0, bN}, {zeta, 0, cN},
    WorkingPrecision -> prec, AccuracyGoal -> prec - 5,
    PrecisionGoal -> prec - 8, MaxRecursion -> 18,
    Method -> "GlobalAdaptive"];
  Print["  K[(0,1,0),(0,1,0)] = D_yy = ", N[val, 30]];
  Print["  Time: ", N[AbsoluteTime[] - t0, 4], " s"];
  val];

DyyWP50 = testKelement[50];
DyyWP70 = testKelement[70];

Print[""];
Print["Comparison:"];
Print["  WP=50: ", N[DyyWP50, 35]];
Print["  WP=70: ", N[DyyWP70, 35]];
Print[""];
Print["  Relative diff: ", N[Abs[(DyyWP50 - DyyWP70)/DyyWP70], 6]];
Print[""];
Print["  -> If diff < 1e-30, WP=50 truly has 30+ digit precision"];
Print["  -> If diff ~ 1e-16, WP=50 was only ~16 digit (matches Python)"];
Print["  -> Anything in between: WP=50 has that many true digits"];
