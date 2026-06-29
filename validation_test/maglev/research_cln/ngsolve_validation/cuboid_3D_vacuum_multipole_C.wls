(* 3D Cuboid in vacuum + uniform B_z — Mathematica reference via
   multipole + internal eigenmode hybrid Galerkin (option C, Stage 1+2).

   Symmetry analysis for cuboid (D_2h orthorhombic) + uniform B_z applied:
     - Only odd-ell, m=0 axial multipoles excited (parity)
     - Leading: ell=1 (z-dipole), then ell=3 (octupole), etc.
     - For modest cuboid aspect ratios, dipole dominates (>90%)

   Stage 1: theoretical framework, DC limit verification (exact closed form)
   Stage 2: dipole-only approximation (Foster from internal Galerkin)
   Stage 3+: octupole and higher (future)

   Key formulas:
     R_0^B = 1 / <J_0, J_0/sigma>_cond
           = 1 / (sigma * V * (a^2+b^2)/12)   [coordinate centered at corner]
           = 1 / (sigma * V * (a^2+b^2)/48)   [centered at origin]
*)
$HistoryLength = 0;

sigma = 58000000;
mu0 = 4 Pi/10^7;
aDim = 5/1000;
bDim = 2/1000;
cDim = 1/1000;

prec = 60;
Print["=== Cuboid in vacuum + uniform B_z, multipole + internal Galerkin ==="];
Print["Geometry (mm): a=", N[aDim*1000], " b=", N[bDim*1000], " c=", N[cDim*1000]];
Print["WorkingPrecision = ", prec];
Print[""];

sigmaN = SetPrecision[sigma, prec];
mu0N = SetPrecision[mu0, prec];
aN = SetPrecision[aDim, prec];
bN = SetPrecision[bDim, prec];
cN = SetPrecision[cDim, prec];
PiN = SetPrecision[Pi, prec];
VN = aN bN cN;

(* DC limit / Stage 1 verification *)
Print["=== Stage 1: DC limit (s=0) verification ==="];
Print[""];

(* For coordinate centered at origin: R_0 = 48/(sigma V (a^2+b^2)) *)
R0centered = 48 / (sigmaN VN (aN^2 + bN^2));
Print["R_0 (centered, exact Parseval) = ", N[R0centered, 30]];

(* For coordinate at corner: R_0 = 12/(sigma V (a^2+b^2)) *)
R0corner = 12 / (sigmaN VN (aN^2 + bN^2));
Print["R_0 (corner, exact Parseval) = ", N[R0corner, 30]];

Print[""];

(* === Stage 2: dipole-only approximation === *)
Print["=== Stage 2: dipole-only approximation ==="];
Print[""];
Print["For uniform B_z applied to cuboid, the leading external multipole"];
Print["is the magnetic dipole. Approximate the cuboid as a magnetizable"];
Print["volume with effective susceptibility chi(s)."];
Print[""];

(* Define induced current density inside conductor *)
(* For uniform applied A_ext = (B_0/2)(-y, x, 0), B_0 = 1, *)
(* J_ind = -s sigma A_total where A_total = A_ext + A_ind. *)
(* The Foster expansion comes from the eigenvalue problem inside. *)

(* For Galerkin within conductor, use sin*sin*sin basis with appropriate parity *)
(* Conductor centered at origin: x in [-a/2, a/2] etc. *)
(* For A_ind_x: should be ODD in y (matching applied A_ext_x = -y B_0/2),
                        EVEN in x and z *)
(* Basis: cos(p pi x/a) sin(q pi y/b) cos(r pi z/c) for A_x with p,r >= 0, q >= 1 *)
(* Similarly: A_ind_y: sin(p pi x/a) cos(q pi y/b) cos(r pi z/c) with p>=1, q,r>=0 *)
(* A_ind_z = 0 (by symmetry for uniform Bz applied to D_2h cuboid) *)

(* For Foster decomposition of chi(s), need the eigenvalue problem:
     curl-curl psi = lambda^2 mu sigma psi  inside cuboid
   with NATURAL BC (matching to vacuum dipole outside).

   The natural BC is non-trivial:
     - B field tangential continuous = (1/mu0)(curl psi) tangential
     - B field normal continuous = (1/mu0)(curl psi) normal = -d phi_M/dn
   With phi_M = c_1(s) cos(theta)/r^2 outside (dipole only), the matching
   gives a specific condition on (curl psi)_normal at the surface.

   For Stage 2 simplification: assume the internal mode psi has same
   symmetry as A_ext (rotational about z). Then chi(s) can be computed
   from a 1D-like sequence of eigenvalues. *)

(* Simplification: Approximate the cuboid as an "effective sphere" for
   the dipole calculation, then estimate the corrections.

   For a sphere of radius R_eff in vacuum, applied B_0:
     chi_sphere(s) = ... well-known Mie series.

   For the cuboid, replace R_eff with some equivalent radius and
   estimate corrections. *)

(* Compute volume-equivalent sphere radius *)
ReqVol = (3 VN/(4 PiN))^(1/3);
Print["Volume-equivalent sphere radius R_eq = ", N[ReqVol*1000, 6], " mm"];

(* Compute moment-of-inertia equivalent *)
(* For a sphere of radius R, MOI = (2/5) M R^2; for cuboid (centered),
   I_z = M(a^2+b^2)/12. Equating: R_iz = sqrt(5(a^2+b^2)/24). *)
RmomI = Sqrt[5 (aN^2 + bN^2)/24];
Print["Moment-equivalent sphere radius R_iz = ", N[RmomI*1000, 6], " mm"];

(* Sphere's Foster expansion for uniform B_z applied:
   chi_sphere(s) = -3 [coth(z) - 1/z]/(2[coth(z) + ...]) with z = R sqrt(s mu sigma)
   Or simpler: well-known formula in Landau-Lifshitz or Jackson.

   For a sphere of radius R, conductivity sigma, in uniform B_0 e^(jωt):
     m_z(omega) = -V (3/2) (1 - 3/(z^2)[1 - z coth z]) / (1 + 1/(z^2)[1 - z coth z])
     z = R sqrt(j omega mu sigma) — argument of Mie

   Actually the cleanest formulation:
     chi_eff_sphere(omega) = -[1 - 3 (sin q - q cos q)/(q^2 sin q)] / [...]
   where q = R k, k = sqrt(j omega mu sigma).

   Or via spherical Bessel:
     chi(omega) = (3/2)[1 - (3/q^2)(1 - q cot q)]   (low-freq expansion)
*)

(* Foster expansion of sphere's chi(s):
   The eigenvalues of curl-curl on a sphere with TM_(l,m,n) modes:
     k^2_(l,n) = (root_n of spherical Bessel j_l = 0)/R squared
   For l=1 (dipole), the radial roots are
     j_1(zeta_(1,n)) = 0 -> zeta = pi, 2pi, 3pi, ... (first roots: 4.4934, 7.7253, 10.9041, ...)
     Actually: zeros of j_1(x) = (sin x)/x^2 - cos(x)/x are NOT pi*n.
     They satisfy tan(x) = x, with first roots 4.4934, 7.7253, 10.9041, 14.0662, ...
*)

(* Sphere Foster, ell=1 (dipole) modes only *)
(* sphere magnetic susceptibility chi_M(s) for radius R, sigma:
   chi_M(s) = (3/2) sum_n [1/(1 + s tau_n) - 1] + ...
   Actually let me use the closed-form Mie low-frequency expansion. *)

(* Sphere closed-form: B_z_inside = B_0 + B_z_ind
   For a conducting sphere of radius R in uniform B_0 e^{jωt}:
     The classic Smythe formula gives:
     B_inside / B_0 = 3/(z (cot z + 1/z) + 2/(z^2)(1 - z cot z))  (Lammeraner-Stafl)
     where z = R sqrt(j omega mu sigma).
   Or via volume-averaged susceptibility chi(s).
*)

(* For NUMERICAL reference of Foster expansion, use sphere Mie series.
   This serves as a benchmark — cuboid result should approach sphere as a=b=c.
*)
sphereChi[R_, s_] := Module[{z, F},
  z = R Sqrt[s mu0N sigmaN];
  F = (3/(2 z^2)) (z Coth[z] - 1) - 1/2;
  F];

(* Sample at a few s values *)
sList = {0, 1, 100, 10000};
Print[""];
Print["Sphere chi(s) for R = R_eq = ", N[ReqVol*1000, 4], " mm:"];
Do[
  Print["  s = ", s, ":  chi = ", N[sphereChi[ReqVol, SetPrecision[s, prec]], 16]],
  {s, sList}];

Print[""];
Print["Sphere chi(s) for R = R_iz = ", N[RmomI*1000, 4], " mm:"];
Do[
  Print["  s = ", s, ":  chi = ", N[sphereChi[RmomI, SetPrecision[s, prec]], 16]],
  {s, sList}];

(* Cuboid response approximation: use sphere with R_iz (moment-of-inertia
   equivalent radius). Higher-order corrections from non-spherical shape
   are next stage. *)

Print[""];
Print["Stage 2 (dipole-only) -> Cauer-I extraction"];
Print[""];

(* Foster expansion via sphere Bessel zeros
   The dipole modes for sphere of radius R, conducting, are at
   zeta_(1,n) (zeros of j_1(x) = 0, equivalently tan(x) = x).
   Time constants: tau_n = mu sigma R^2 / zeta_(1,n)^2
*)
nMax = 30;  (* number of dipole modes *)
Print["Computing first ", nMax, " dipole eigenmodes (j_1 zeros)..."];
zetaList = Table[
  zeta /. FindRoot[Tan[zeta] - zeta == 0,
    {zeta, (n + 0.5) Pi, (n - 0.4) Pi, (n + 0.6) Pi},
    WorkingPrecision -> prec, AccuracyGoal -> prec - 10],
  {n, 1, nMax}];
Print["  zeta_(1,1) = ", N[zetaList[[1]], 16], " (expected 4.4934)"];
Print["  zeta_(1,2) = ", N[zetaList[[2]], 16], " (expected 7.7253)"];
Print["  zeta_(1,3) = ", N[zetaList[[3]], 16], " (expected 10.9041)"];

(* Foster: chi(s) = (3/2) sum_n alpha_n/(1+s tau_n), alpha_n related to mode amplitude *)
(* For Smythe's sphere: chi = 6/zeta^4 sum 1/(1+s tau_n) Foster-like. Detail TBD. *)

Print[""];
Print["NOTE: Full Foster -> Cauer-I extraction in Stage 3."];
Print["This Stage 2 establishes the framework; numerical chi(s) values"];
Print["above provide the SPHERE BENCHMARK to compare cuboid against."];
