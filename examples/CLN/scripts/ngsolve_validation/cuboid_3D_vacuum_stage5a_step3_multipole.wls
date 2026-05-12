(* Stage 5a Step 3: enriched Galerkin with exterior dipole multipole.

   Internal basis (Step 1+2): A^I_i = curl(z_hat psi_i),
     psi_i = sin(p pi (x+a/2)/a) sin(q pi (y+b/2)/b) sin(r pi (z+c/2)/c),
     p,q,r all odd.  Vanishes at conductor surface.

   External basis (NEW): single mode
     Phi^E(r) = (-y, x, 0)/(r^2 + R_0^2)^(3/2)
   This:
     * is regular at origin (no singularity)
     * decays as 1/r^3 at infinity (dipole tail)
     * has same parity as A_ext (D_2h + uniform B_z)
     * curl Phi^E approaches dipole B-field at large r

   Closed form: K^EE = integral|curl Phi^E|^2 dV = 15 pi^2/(32 R_0^3)
   NIntegrate needed for M^IE, K^IE, M^EE, b^E (all over conductor only).

   Augmented system (n_int + 1 unknowns):
     [K^II  K^IE]  alpha^I    M^II  M^IE   alpha^I       b^I
     [K^EI  K^EE] [alpha^E] + s mu0 sigma [M^EI  M^EE] [alpha^E] = -s mu0 sigma [b^E]

   Solve for chi(s) Taylor coefficients. *)

$HistoryLength = 0;
prec = 60;

Print["=== Stage 5a Step 3: enriched Galerkin (interior + dipole) ==="];
Print[""];

(* Constants *)
sigmaN = SetPrecision[58000000, prec];
mu0N = SetPrecision[4 Pi/10^7, prec];
aN = SetPrecision[5/1000, prec];
bN = SetPrecision[2/1000, prec];
cN = SetPrecision[1/1000, prec];
PiN = SetPrecision[Pi, prec];
VN = aN bN cN;
ReffSq = 5 (aN^2 + bN^2)/24;

(* Sphere reference alpha_n *)
alphaSphere[R_, k_] := Module[{x = mu0N sigmaN R^2, ck},
  ck = Switch[k, 1, -1/10, 2, +1/105, 3, -1/1050, 4, +1/10395, 5, -1/135135,
    _, (-1)^k * 9 * SetPrecision[Zeta[2 k + 2], prec]/PiN^(2 k + 2)];
  ck * x^k];

Print["Sphere reference (R_eff^2 = 5(a^2+b^2)/24):"];
Do[Print["  alpha_", k, "^sph = ", N[alphaSphere[Sqrt[ReffSq], k], 16]],
  {k, 1, 5}];
Print[""];

(* === Interior basis (closed form, from Step 1+2) === *)
Pmax = 3;
basisIndices = Flatten[Table[{2 P + 1, 2 Q + 1, 2 R + 1},
  {P, 0, (Pmax - 1)/2}, {Q, 0, (Pmax - 1)/2}, {R, 0, (Pmax - 1)/2}], 2];
nInt = Length[basisIndices];
Print["Interior basis: Pmax=", Pmax, ", nInt=", nInt];
Print["Indices: ", basisIndices];
Print[""];

(* Closed-form interior matrices (diagonal) *)
mIIElement[p_, q_, r_] := (VN/8) PiN^2 (p^2/aN^2 + q^2/bN^2);
kIIElement[p_, q_, r_] := (VN/8) PiN^4 (p^2/aN^2 + q^2/bN^2)*
                          (p^2/aN^2 + q^2/bN^2 + r^2/cN^2);
bIElement[p_, q_, r_] := 8 VN/(p q r PiN^3);

(* Build interior matrices *)
MII = DiagonalMatrix[Table[mIIElement @@ idx, {idx, basisIndices}]];
KII = DiagonalMatrix[Table[kIIElement @@ idx, {idx, basisIndices}]];
bIVec = Table[bIElement @@ idx, {idx, basisIndices}];

(* === External basis: regularized dipole === *)
(* Choose R_0: comparable to half the smallest cuboid dimension *)
R0 = SetPrecision[1/1000, prec];  (* 1 mm = b/2 *)
Print["External basis: Phi^E(r) = (-y, x, 0)/(r^2 + R_0^2)^(3/2)"];
Print["  R_0 = ", N[R0*1000, 6], " mm"];
Print[""];

(* Closed form K^EE *)
KEE = 15 PiN^2 / (32 R0^3);
Print["K^EE (closed form) = 15 pi^2/(32 R_0^3) = ", N[KEE, 16]];
Print[""];

(* === NIntegrate parameters === *)
wp = 50;  (* WorkingPrecision *)
ag = 40;  (* AccuracyGoal *)
pg = 35;  (* PrecisionGoal *)
maxRecursion = 12;

(* Helper: Phi^E components and curl Phi^E components, for use in NIntegrate *)
PhiEx[xx_, yy_, zz_] := -yy/(xx^2 + yy^2 + zz^2 + R0^2)^(3/2);
PhiEy[xx_, yy_, zz_] := xx/(xx^2 + yy^2 + zz^2 + R0^2)^(3/2);
(* PhiEz = 0 *)

(* curl Phi^E *)
curlPhiEx[xx_, yy_, zz_] := 3 xx zz/(xx^2 + yy^2 + zz^2 + R0^2)^(5/2);
curlPhiEy[xx_, yy_, zz_] := 3 yy zz/(xx^2 + yy^2 + zz^2 + R0^2)^(5/2);
curlPhiEz[xx_, yy_, zz_] := (2 (zz^2 + R0^2) - (xx^2 + yy^2))/
                             (xx^2 + yy^2 + zz^2 + R0^2)^(5/2);

(* Interior basis psi_i and curl A^I_i = (d_xz psi, d_yz psi, -(d_xx+d_yy) psi) *)
psiSym[p_, q_, r_, xx_, yy_, zz_] :=
  Sin[p PiN (xx + aN/2)/aN] * Sin[q PiN (yy + bN/2)/bN] * Sin[r PiN (zz + cN/2)/cN];

curlAIxSym[p_, q_, r_, xx_, yy_, zz_] :=
  (p PiN/aN)(r PiN/cN) Cos[p PiN (xx + aN/2)/aN] *
    Sin[q PiN (yy + bN/2)/bN] * Cos[r PiN (zz + cN/2)/cN];

curlAIySym[p_, q_, r_, xx_, yy_, zz_] :=
  (q PiN/bN)(r PiN/cN) Sin[p PiN (xx + aN/2)/aN] *
    Cos[q PiN (yy + bN/2)/bN] * Cos[r PiN (zz + cN/2)/cN];

curlAIzSym[p_, q_, r_, xx_, yy_, zz_] :=
  -(p^2 PiN^2/aN^2 + q^2 PiN^2/bN^2) Sin[p PiN (xx + aN/2)/aN] *
    Sin[q PiN (yy + bN/2)/bN] * Sin[r PiN (zz + cN/2)/cN];

AIxSym[p_, q_, r_, xx_, yy_, zz_] :=
  (q PiN/bN) Sin[p PiN (xx + aN/2)/aN] * Cos[q PiN (yy + bN/2)/bN] *
    Sin[r PiN (zz + cN/2)/cN];

AIySym[p_, q_, r_, xx_, yy_, zz_] :=
  -(p PiN/aN) Cos[p PiN (xx + aN/2)/aN] * Sin[q PiN (yy + bN/2)/bN] *
    Sin[r PiN (zz + cN/2)/cN];

(* === Compute NIntegrate matrices === *)

(* M^EE = integral_cond (Phi^E . Phi^E) dV = integral (x^2+y^2)/(r^2+R_0^2)^3 dV *)
Print["Computing M^EE (1 NIntegrate)..."];
tStart = AbsoluteTime[];
MEE = NIntegrate[
  (xx^2 + yy^2)/(xx^2 + yy^2 + zz^2 + R0^2)^3,
  {xx, -aN/2, aN/2}, {yy, -bN/2, bN/2}, {zz, -cN/2, cN/2},
  WorkingPrecision -> wp, AccuracyGoal -> ag, PrecisionGoal -> pg,
  MaxRecursion -> maxRecursion];
Print["  M^EE = ", N[MEE, 16], "  (", N[AbsoluteTime[] - tStart, 4], " s)"];

(* b^E = integral_cond (Phi^E . A_ext) dV = (1/2) integral (x^2+y^2)/(r^2+R_0^2)^(3/2) dV *)
Print["Computing b^E (1 NIntegrate)..."];
tStart = AbsoluteTime[];
bE = (1/2) NIntegrate[
  (xx^2 + yy^2)/(xx^2 + yy^2 + zz^2 + R0^2)^(3/2),
  {xx, -aN/2, aN/2}, {yy, -bN/2, bN/2}, {zz, -cN/2, cN/2},
  WorkingPrecision -> wp, AccuracyGoal -> ag, PrecisionGoal -> pg,
  MaxRecursion -> maxRecursion];
Print["  b^E = ", N[bE, 16], "  (", N[AbsoluteTime[] - tStart, 4], " s)"];

(* M^IE_i = integral_cond Phi^I_i . Phi^E dV
        = integral (A^I_x * (-y)/(r^2+R_0^2)^(3/2) + A^I_y * x/(r^2+R_0^2)^(3/2)) dV *)
Print["Computing M^IE (", nInt, " NIntegrates)..."];
tStart = AbsoluteTime[];
MIEVec = Table[Module[{p, q, r, integrand, val},
  {p, q, r} = basisIndices[[i]];
  integrand = (AIxSym[p, q, r, xx, yy, zz] (-yy) +
               AIySym[p, q, r, xx, yy, zz] xx)/
              (xx^2 + yy^2 + zz^2 + R0^2)^(3/2);
  val = NIntegrate[integrand,
    {xx, -aN/2, aN/2}, {yy, -bN/2, bN/2}, {zz, -cN/2, cN/2},
    WorkingPrecision -> wp, AccuracyGoal -> ag, PrecisionGoal -> pg,
    MaxRecursion -> maxRecursion];
  Print["  M^IE[", basisIndices[[i]], "] = ", N[val, 12]];
  val],
  {i, nInt}];
Print["  M^IE total time: ", N[AbsoluteTime[] - tStart, 4], " s"];

(* K^IE_i = integral_cond curl A^I_i . curl Phi^E dV  (interior support only) *)
Print["Computing K^IE (", nInt, " NIntegrates)..."];
tStart = AbsoluteTime[];
KIEVec = Table[Module[{p, q, r, integrand, val},
  {p, q, r} = basisIndices[[i]];
  integrand = (curlAIxSym[p, q, r, xx, yy, zz] curlPhiEx[xx, yy, zz] +
               curlAIySym[p, q, r, xx, yy, zz] curlPhiEy[xx, yy, zz] +
               curlAIzSym[p, q, r, xx, yy, zz] curlPhiEz[xx, yy, zz]);
  val = NIntegrate[integrand,
    {xx, -aN/2, aN/2}, {yy, -bN/2, bN/2}, {zz, -cN/2, cN/2},
    WorkingPrecision -> wp, AccuracyGoal -> ag, PrecisionGoal -> pg,
    MaxRecursion -> maxRecursion];
  Print["  K^IE[", basisIndices[[i]], "] = ", N[val, 12]];
  val],
  {i, nInt}];
Print["  K^IE total time: ", N[AbsoluteTime[] - tStart, 4], " s"];
Print[""];

(* === Build augmented system === *)
(* Block structure:
     M = [[ M^II   M^IE ]
          [ M^EI   M^EE ]]   (n_int+1) x (n_int+1)
   Same for K.
   Source vector b = [b^I; b^E].
*)
nTot = nInt + 1;
Print["Building augmented system (n_total = ", nTot, ")"];

MAug = Table[0, {nTot}, {nTot}];
KAug = Table[0, {nTot}, {nTot}];
bAug = Table[0, {nTot}];

Do[MAug[[i, i]] = MII[[i, i]]; KAug[[i, i]] = KII[[i, i]]; bAug[[i]] = bIVec[[i]],
  {i, nInt}];
Do[MAug[[i, nTot]] = MIEVec[[i]]; MAug[[nTot, i]] = MIEVec[[i]];
   KAug[[i, nTot]] = KIEVec[[i]]; KAug[[nTot, i]] = KIEVec[[i]],
  {i, nInt}];
MAug[[nTot, nTot]] = MEE;
KAug[[nTot, nTot]] = KEE;
bAug[[nTot]] = bE;
Print[""];

(* === Compute Taylor coefficients alpha_n of chi(s) ===
   Galerkin equation: (K + s mu_0 sigma M) alpha = -s mu_0 sigma b
   chi(s) = -(s mu0 sigma/V)[<A_ext, A_ext> + <A_ext, A_ind>]
          = -(s mu0 sigma/V)[V(a^2+b^2)/48 + b . alpha]
          = alpha_LF * s + (mu0 sigma)^2 s^2/V * b . (K^{-1}) b - ... etc.

   Let nu = s mu0 sigma. alpha(nu) satisfies (K + nu M) alpha = -nu b.
   Series alpha = -nu K^{-1} b + nu^2 K^{-1} M K^{-1} b - nu^3 K^{-1} M K^{-1} M K^{-1} b + ...
   <b, alpha> = -nu b^T K^{-1} b + nu^2 b^T K^{-1} M K^{-1} b - ...
   chi_corr = -(nu/V) <b, alpha> = (nu^2/V) b^T K^{-1} b - (nu^3/V) b^T K^{-1} M K^{-1} b + ...
   alpha_2 = (mu0 sigma)^2/V * b^T K^{-1} b
   alpha_3 = -(mu0 sigma)^3/V * b^T (K^{-1} M K^{-1}) b
   etc.
*)

KAugInv = Inverse[KAug];

Print["=== Taylor coefficients of chi(s) ==="];
alphaLFRef = -mu0N sigmaN (aN^2 + bN^2)/48;
Print["alpha_1 = -mu0 sigma (a^2+b^2)/48 = ", N[alphaLFRef, 16],
  "  (closed form, basis-independent)"];

(* alpha_n for n >= 2 *)
KinvMpow[n_] := Module[{cur = KAugInv}, Do[cur = cur . MAug . KAugInv, {n - 1}]; cur];
alphaN[n_] := Module[{matn},
  matn = If[n == 1, KAugInv, Nest[# . MAug . KAugInv &, KAugInv, n - 1]];
  (-1)^(n - 1) (mu0N sigmaN)^n / VN * (bAug . matn . bAug)];

(* Note: alpha_2 = +(mu0 sigma)^2/V b K^{-1} b *)
Do[Print["alpha_", k, " = ", N[alphaN[k], 16],
  "    (sphere ref: ", N[alphaSphere[Sqrt[ReffSq], k], 12],
  ", ratio: ", N[alphaN[k]/alphaSphere[Sqrt[ReffSq], k], 8], ")"],
  {k, 2, 5}];
Print[""];

(* === Compare to interior PEC-only baseline === *)
Print["=== Comparison ==="];
alpha2InteriorOnly = (mu0N sigmaN)^2 / VN * Sum[bIVec[[i]]^2/KII[[i, i]], {i, nInt}];
Print["alpha_2 interior PEC only (Step 1+2): ", N[alpha2InteriorOnly, 16]];
Print["alpha_2 enriched (interior + dipole): ", N[alphaN[2], 16]];
Print["alpha_2 sphere R_eff:                 ",
  N[alphaSphere[Sqrt[ReffSq], 2], 16]];
Print[""];
Print["Improvement ratio (enriched/interior_only): ",
  N[alphaN[2]/alpha2InteriorOnly, 8]];
Print["Gap to sphere ref (enriched/sphere):       ",
  N[alphaN[2]/alphaSphere[Sqrt[ReffSq], 2], 8]];
Print[""];

Print["=== Stage 5a Step 3 summary ==="];
Print["* Single dipole exterior basis enriches interior PEC by factor ",
  N[alphaN[2]/alpha2InteriorOnly, 4]];
Print["* Remaining gap to sphere reference: ratio = ",
  N[alphaN[2]/alphaSphere[Sqrt[ReffSq], 2], 4]];
Print["* Higher multipoles (octupole etc) and R_0 sweep would further close gap"];
Print[""];

(* Save numerical results *)
results = <|
  "Pmax" -> Pmax,
  "nInt" -> nInt,
  "R0_mm" -> N[R0*1000, 6],
  "alphaLF" -> N[alphaLFRef, 30],
  "alpha2_interiorOnly" -> N[alpha2InteriorOnly, 30],
  "alpha2_enriched" -> N[alphaN[2], 30],
  "alpha2_sphereRef" -> N[alphaSphere[Sqrt[ReffSq], 2], 30],
  "alpha3_enriched" -> N[alphaN[3], 30],
  "alpha4_enriched" -> N[alphaN[4], 30]
|>;
Print["Results: ", results];
