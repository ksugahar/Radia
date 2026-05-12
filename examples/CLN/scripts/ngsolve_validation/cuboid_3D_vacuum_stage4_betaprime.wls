(* Stage 4 (beta'): Semi-analytical hybrid for cuboid in vacuum + uniform B_z.

   Strategy:
     1. Closed-form low-frequency leading coefficient for cuboid:
          chi_cub(s) = alpha_LF^cub * s + O(s^2)
          alpha_LF^cub = -mu_0 sigma (a^2+b^2)/48     (centered convention)
     2. Closed-form low-freq leading for sphere of radius R:
          chi_sph(s) = alpha_LF^sph * s + O(s^2)
          alpha_LF^sph = -mu_0 sigma R^2 / 10
     3. Match alpha_LF: R_eff defined by mu_0 sigma R_eff^2/10 = mu_0 sigma (a^2+b^2)/48
        => R_eff^2 = 5(a^2+b^2)/24
        => R_eff = sqrt(5(a^2+b^2)/24)
        This is exactly the moment-of-inertia-equivalent radius (R_iz) used in Stage 3.
        => Stage 3 R_iz choice was the *unique* low-frequency matching.
     4. The remaining shape correction is a higher-order (s^2 and beyond) effect.
        Compare:
          eta_2(s) := chi_cub(s) / chi_sph(s; R_eff) as a function of s.
        - eta_2(s -> 0) = 1 by construction.
        - eta_2(s -> infty) -> ratio of "perfect screening volumes"
          (limit values give shape-screening correction).
     5. Produce a reference table chi_sph(s; R_eff) at decade-spaced s values.
        NGSolve cuboid_521_vacuum_caseB.py will sweep the same s values to
        provide chi_cub(s); ratio = sphere reference to NGSolve cuboid result.
*)
$HistoryLength = 0;

sigma = 58000000;
mu0 = 4 Pi/10^7;
aDim = 5/1000;
bDim = 2/1000;
cDim = 1/1000;

prec = 60;
sigmaN = SetPrecision[sigma, prec];
mu0N = SetPrecision[mu0, prec];
aN = SetPrecision[aDim, prec];
bN = SetPrecision[bDim, prec];
cN = SetPrecision[cDim, prec];
PiN = SetPrecision[Pi, prec];
VN = aN bN cN;

Print["=== Stage 4 (beta'): semi-analytical hybrid ==="];
Print["Geometry (mm): a=", N[aDim*1000], " b=", N[bDim*1000], " c=", N[cDim*1000]];
Print["Volume V = ", N[VN, 8], " m^3"];
Print[""];

(* === Step 1: cuboid low-frequency closed form === *)
Print["Step 1: Cuboid low-frequency leading coefficient"];
Print["  Derivation (centered cuboid, applied B_0 hat z):"];
Print["    A_ext = (B_0/2)(-y, x, 0)"];
Print["    J_ind ~ -s sigma A_ext  (DC limit of full E = -dA/dt)"];
Print["    m_z = (1/2) integral (r x J)_z dV"];
Print["        = -(s sigma B_0/4) integral (x^2+y^2) dV"];
Print["    integral_centered (x^2+y^2) dV = V (a^2+b^2)/12"];
Print["    => m_z = -s sigma B_0 V (a^2+b^2)/48"];
Print["    => chi_cub(s) = m_z mu_0 / (V B_0) = -s mu_0 sigma (a^2+b^2)/48 + O(s^2)"];
Print[""];
alphaLFcub = -mu0N sigmaN (aN^2 + bN^2)/48;
Print["  alpha_LF^cub = -mu_0 sigma (a^2+b^2)/48 = ", N[alphaLFcub, 16]];
Print[""];

(* === Step 2: sphere low-frequency closed form === *)
Print["Step 2: Sphere low-frequency leading coefficient"];
Print["  From chi_sph(s) = -3/2 + (9/(2 z^2))(z coth z - 1), z = R sqrt(s mu sigma):"];
Print["    Series expansion at s=0:"];
Print["      z coth z - 1 = z^2/3 - z^4/45 + ..."];
Print["      (9/(2 z^2))(z coth z - 1) = 3/2 - z^2/10 + z^4/(2*45) - ..."];
Print["      chi_sph(s) = -z^2/10 + ... = -mu_0 sigma R^2 / 10 * s + O(s^2)"];
Print["    => alpha_LF^sph = -mu_0 sigma R^2 / 10"];
Print[""];

(* Verify by symbolic series *)
chiSphereExact[R_, s_] := Module[{z}, z = R Sqrt[s mu0N sigmaN];
  -3/2 + (9/(2 z^2)) (z Coth[z] - 1)];
sSym = sv;
seriesCheck = Series[-3/2 + (9/(2 z^2)) (z Coth[z] - 1), {z, 0, 6}];
Print["  Series verification: (9/(2 z^2))(z coth z - 1) - 3/2 ="];
Print["    ", Normal[seriesCheck]];
Print["  z^2 = R^2 mu sigma s, so coefficient of s is R^2 mu sigma / 10 with sign -1"];
Print[""];

(* === Step 3: alpha matching => R_eff === *)
Print["Step 3: alpha matching alpha_LF^sph(R_eff) = alpha_LF^cub"];
Print["  -mu_0 sigma R_eff^2/10 = -mu_0 sigma (a^2+b^2)/48"];
Print["  R_eff^2 = 10 (a^2+b^2)/48 = 5(a^2+b^2)/24"];
Print["  R_eff = sqrt(5(a^2+b^2)/24)"];
Print[""];
Reff = Sqrt[5 (aN^2 + bN^2)/24];
Print["  R_eff = ", N[Reff*1000, 12], " mm"];
Print["  (equals R_iz used in Stage 3 - confirms moment-of-inertia is the"];
Print["   unique alpha-matching radius)"];
Print[""];

(* Verify alpha matching numerically *)
alphaLFsph = -mu0N sigmaN Reff^2 / 10;
Print["  Verification: alpha_LF^sph(R_eff) = ", N[alphaLFsph, 20]];
Print["                alpha_LF^cub        = ", N[alphaLFcub, 20]];
Print["                rel diff            = ",
  N[Abs[alphaLFsph - alphaLFcub]/Abs[alphaLFcub], 5]];
Print[""];

(* === Step 4: higher-order coefficient for shape correction analysis === *)
Print["Step 4: Higher-order coefficients for shape correction analysis"];
Print[""];
Print["  Sphere full series at R_eff:"];
seriesEff = Series[chiSphereExact[Reff, sv], {sv, 0, 4}];
Print["    chi_sph(s; R_eff) = ", Normal[seriesEff]];
Print[""];

(* Cuboid second-order coefficient (closed form via internal Galerkin would be needed;
   for now, document that this requires Stage 5 internal eigenmode work) *)
Print["  Cuboid 2nd-order coefficient alpha_2^cub:"];
Print["    Requires solving curl-curl A_1 / mu_0 = -s sigma A_ext on full domain"];
Print["    and computing m_z[A_1].  This is a one-shot static Poisson-like solve."];
Print["    Closed form NOT generally available for cuboid."];
Print["    => Stage 5 task: NGSolve solves it numerically; compare to sphere prediction."];
Print[""];

(* === Step 5: reference table for NGSolve comparison === *)
Print["Step 5: Reference table chi_sph(s; R_eff)"];
Print[""];

(* Use Foster fast form for stability across decade-spaced s *)
nMax = 80;
chiSphereFoster[R_, s_] := Module[{tList},
  tList = Table[mu0N sigmaN R^2 / (n PiN)^2, {n, 1, nMax}];
  -9 s mu0N sigmaN R^2 Sum[1/((n PiN)^4 (1 + s tList[[n]])), {n, 1, nMax}]];

(* Decade-spaced s values *)
sValues = {0, 1, 10, 100, 1000, 10^4, 10^5, 10^6, 10^7, 10^8};

Print["  s [rad/s]          | chi_sph(s; R_eff)             | |chi|"];
Print[StringJoin[Table["-", {78}]]];
Do[Module[{sv2, chiVal, chiAbs, sStr},
   sv2 = If[s === 0, SetPrecision[0, prec], SetPrecision[s, prec]];
   chiVal = If[s === 0, 0, N[chiSphereFoster[Reff, sv2], 16]];
   chiAbs = N[Abs[chiVal], 16];
   sStr = If[s === 0, "0", "10^" <> ToString[Round[Log10[s]]]];
   Print[StringJoin[{StringPadLeft[sStr, 18], " | ",
     StringPadRight[ToString[NumberForm[chiVal, 16]], 28], " | ",
     ToString[NumberForm[chiAbs, 12]]}]]],
  {s, sValues}];
Print[""];

(* DC asymptote check *)
Print["  Low-frequency check: chi/s = ", N[alphaLFsph, 16], " at s=0+"];
Do[Module[{sv2, chiVal, ratio},
  sv2 = SetPrecision[s, prec];
  chiVal = N[chiSphereFoster[Reff, sv2], 20];
  ratio = N[chiVal/sv2, 16];
  Print["  s = ", s, ":  chi(s)/s = ", ratio,
    "  (rel err vs alpha_LF: ", N[Abs[ratio - alphaLFsph]/Abs[alphaLFsph], 5], ")"]],
  {s, {1, 10, 100, 1000}}];
Print[""];

(* High-frequency screening asymptote *)
Print["  High-frequency screening:"];
Do[Module[{sv2, chiVal},
  sv2 = SetPrecision[s, prec];
  chiVal = N[chiSphereFoster[Reff, sv2], 16];
  Print["  s = 10^", Round[Log10[s]], ":  chi = ", chiVal,
    "  (-> -3/2 at infinity)"]],
  {s, {10^6, 10^7, 10^8, 10^9, 10^10}}];
Print[""];

(* Export table as JSON-friendly format *)
Print["=== Exporting reference table ==="];
sExport = {0, 1, 10, 100, 1000, 10^4, 10^5, 10^6, 10^7, 10^8};
chiExport = Table[
  If[s === 0, 0, N[chiSphereFoster[Reff, SetPrecision[s, prec]], 16]],
  {s, sExport}];
ReffExport = N[Reff, 16];
alphaLFExport = N[alphaLFcub, 16];

(* Export to JSON-like text file *)
fout = OpenWrite["cuboid_3D_stage4_betaprime_sphere_reference.dat"];
WriteString[fout, "# Stage 4 (beta') sphere reference for cuboid comparison\n"];
WriteString[fout, "# Geometry: a=", N[aDim*1000], "mm  b=", N[bDim*1000],
  "mm  c=", N[cDim*1000], "mm\n"];
WriteString[fout, "# R_eff = sqrt(5(a^2+b^2)/24) = ",
  ToString[N[Reff*1000, 16]], " mm\n"];
WriteString[fout, "# alpha_LF^cub = -mu_0 sigma (a^2+b^2)/48 = ",
  ToString[N[alphaLFcub, 16]], "\n"];
WriteString[fout, "# Columns: s [rad/s]   chi_sph(s; R_eff)\n"];
Do[Module[{sv3, cv},
  sv3 = sExport[[i]];
  cv = chiExport[[i]];
  WriteString[fout, ToString[sv3], "  ", ToString[N[cv, 16]], "\n"]],
  {i, Length[sExport]}];
Close[fout];
Print["Reference table written to cuboid_3D_stage4_betaprime_sphere_reference.dat"];
Print[""];

Print["=== Stage 4 (beta') summary ==="];
Print["  alpha_LF matching gives R_eff = R_iz = sqrt(5(a^2+b^2)/24) = ",
  N[Reff*1000, 8], " mm"];
Print["  Sphere chi(s; R_eff) is exact reference at low s."];
Print["  Higher-order shape correction eta_2(s) = chi_cub(s)/chi_sph(s; R_eff)"];
Print["  -> 1 at s=0 (by construction)"];
Print["  -> shape-dependent constant at s=infty (perfect screening volume ratio)"];
Print["  Stage 5: NGSolve sweep at the same s values to extract chi_cub(s)."];
