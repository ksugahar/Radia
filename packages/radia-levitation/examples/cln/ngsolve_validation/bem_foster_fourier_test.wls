(* Task A: Fourier-based K element — verification.

   For vacuum 3D problem, Foster eigenmode integrals are:
     K_{(p_1q_1r_1),(p_2q_2r_2)} = int int x^{p_1}y^{q_1}z^{r_1} x'^{p_2}y'^{q_2}z'^{r_2} / |r-r'| dV dV'

   After variable change u = r - r', polynomial integration over R = (r+r')/2,
   reduces to 3D integral over u with weight g_x(xi) g_y(eta) g_z(zeta) divided by |u|.

   Fourier path:
     1/|u| = (1/(2 pi^2)) int_R^3 e^{i k.u} / k^2 d^3 k
   So K = 8 int_{[0,h_x]x[0,h_y]x[0,h_z]} g_x g_y g_z / |u| du
        = (8/(2 pi^2)) int_R^3 G_x(k_x) G_y(k_y) G_z(k_z) / k^2 d^3 k
   where G_alpha(k_alpha) = int_0^{h_alpha} g_alpha(xi) cos(k_alpha xi) d xi
   (taking real part since g_alpha and 1/|u| are real).

   Test 1: K[(0,0,0),(0,0,0)] = 8 * int_0^a int_0^b int_0^c 1/r du dv dw
   This equals 8 * Aharoni0 (octant self-integral). We have closed form at WP=50.

   Date: 2026-05-02
*)

$HistoryLength = 0;
prec = 50;

Print["================================================================"];
Print["Task A: Fourier-based K element verification"];
Print["WorkingPrecision = ", prec];
Print["================================================================"];
Print[""];

aN = SetPrecision[5/1000, prec];
bN = SetPrecision[2/1000, prec];
cN = SetPrecision[1/1000, prec];
PiN = SetPrecision[Pi, prec];

(* === Direct: existing Aharoni0 closed form (v3 verified) === *)
t0 = AbsoluteTime[];
Aharoni0Direct = NIntegrate[1/Sqrt[xi^2 + eta^2 + zeta^2],
  {xi, 0, aN}, {eta, 0, bN}, {zeta, 0, cN},
  WorkingPrecision -> prec, AccuracyGoal -> prec - 5,
  PrecisionGoal -> prec - 8, MaxRecursion -> 16];
Print["Direct Aharoni0 (3D NIntegrate at WP=", prec, "):"];
Print["  ", N[Aharoni0Direct, 30], "  (", N[AbsoluteTime[] - t0, 4], " s)"];
Print[""];

(* === Fourier path === *)
(* For p=q=0, the weight g_alpha = 1 on [0, h_alpha].
   Fourier: G(k) = int_0^h cos(k xi) dxi = sin(kh)/k. *)
ClearAll[Gconst];
Gconst[h_, k_] := Sin[k h]/k;

(* K[(0,0,0),(0,0,0)] / 8 = (1/(2 pi^2)) int_R^3 G_x G_y G_z / k^2 d^3 k
   In spherical coords: d^3 k = k^2 sin(theta) dk dtheta dphi
   k_x = k sin(theta) cos(phi), k_y = k sin(theta) sin(phi), k_z = k cos(theta).

   K_self_octant = (1/(2 pi^2)) int_R^3 G_x(k_x) G_y(k_y) G_z(k_z) / k^2 d^3 k
                 = (1/(2 pi^2)) int_0^inf dk int_0^pi sin(theta) dtheta int_0^{2pi} dphi G_x G_y G_z
                   (k^2 cancels with k^2 sin(theta))

   Wait: d^3 k = k^2 sin(theta) dk dtheta dphi, divided by k^2 = sin(theta) dk dtheta dphi.
*)
Print["Fourier path verification:"];
Print["  G_x(k) = Sin[k a] / k for unit-weight g_x"];
Print["  G_y(k) = Sin[k b] / k"];
Print["  G_z(k) = Sin[k c] / k"];
Print[""];

(* Try integrating over (theta, phi, k) numerically first *)
Print["Step 1: 3D k-integral via NIntegrate (verify physics)"];
t0 = AbsoluteTime[];
AharoniFourier = (1/(2 PiN^2)) NIntegrate[
  Gconst[aN, k Sin[theta] Cos[phi]] *
  Gconst[bN, k Sin[theta] Sin[phi]] *
  Gconst[cN, k Cos[theta]] *
  Sin[theta],
  {k, 0, Infinity}, {theta, 0, PiN}, {phi, 0, 2 PiN},
  WorkingPrecision -> 30, AccuracyGoal -> 20,
  PrecisionGoal -> 15, MaxRecursion -> 12,
  Method -> "GlobalAdaptive"];
Print["  K_self_octant (Fourier path NIntegrate at WP=30) = ",
  N[AharoniFourier, 16],
  "  (", N[AbsoluteTime[] - t0, 4], " s)"];
Print[""];
Print["  Compare to direct Aharoni0 = ", N[Aharoni0Direct, 16]];
Print["  Ratio = ", N[AharoniFourier/Aharoni0Direct, 16]];
Print[""];

(* === Try closed form via Mathematica Integrate === *)
Print["Step 2: Try Integrate[] for full closed form (5 min timeout)"];
t0 = AbsoluteTime[];
AharoniClosedTry = TimeConstrained[
  Integrate[
    (1/(2 PiN^2)) *
    (Sin[k Sin[theta] Cos[phi] aN]/(k Sin[theta] Cos[phi])) *
    (Sin[k Sin[theta] Sin[phi] bN]/(k Sin[theta] Sin[phi])) *
    (Sin[k Cos[theta] cN]/(k Cos[theta])) *
    Sin[theta],
    {k, 0, Infinity}, {theta, 0, PiN}, {phi, 0, 2 PiN}],
  300, $Failed];
Print["  Time: ", N[AbsoluteTime[] - t0, 4], " s"];
If[AharoniClosedTry === $Failed,
  Print["  Integrate timed out. Closed form not directly achievable."],
  Print["  Closed form: ", AharoniClosedTry];
  Print["  Numeric: ", N[AharoniClosedTry, 30]]];
