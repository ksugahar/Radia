(* Stage 5c: PEC TE_z basis (correct H(curl)-conforming) + exterior dipole.

   This corrects Stage 5a Step 3 by replacing the bogus sin-basis with the
   proper PEC TE_z basis from Stage 5b.

   Interior basis (PEC TE_z):
     A^I_{p,q,r} = curl(z_hat phi^TE_{p,q,r})
     phi^TE_{p,q,r} = cos(p pi (x+a/2)/a) cos(q pi (y+b/2)/b) sin(r pi (z+c/2)/c)
   Parity-restricted: p even >= 0, q even >= 0, r odd >= 1, AND only edge
   modes (p=0 or q=0 but not both) couple to uniform B_z source.

   Exterior basis (1 dipole mode):
     Phi^E(r) = (-y, x, 0)/(r^2 + R_0^2)^(3/2)
   Same as Step 3 but now interior is correctly H(curl)-conforming so K^IE
   and M^IE are well-defined as conductor integrals (no surface delta).

   Augmented Galerkin: n_edge + 1 unknowns.
*)
$HistoryLength = 0;
prec = 60;

Print["=== Stage 5c: PEC TE_z + dipole multipole ==="];
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
alphaLFRef = -mu0N sigmaN (aN^2 + bN^2)/48;

(* Sphere reference *)
alphaSphere[R_, k_] := Module[{x = mu0N sigmaN R^2, ck},
  ck = Switch[k, 1, -1/10, 2, +1/105, 3, -1/1050, 4, +1/10395, _,
    (-1)^k * 9 * SetPrecision[Zeta[2 k + 2], prec]/PiN^(2 k + 2)];
  ck * x^k];

Print["alpha_LF       = ", N[alphaLFRef, 16]];
Print["alpha_2^sph    = ", N[alphaSphere[Sqrt[ReffSq], 2], 16]];
Print[""];

(* === Edge mode indices === *)
(* (0, q, r) family: q in {2, 4, ..., 2 Q_max}, r in {1, 3, ..., 2 R_max + 1} *)
QMaxPrime = 3;  (* q' goes up to 3 -> q = 2,4,6 *)
RMax = 2;  (* r goes up to 2*RMax+1 = 5 *)
edge0qr = Flatten[Table[{0, 2 qp, 2 rp + 1}, {qp, 1, QMaxPrime},
    {rp, 0, RMax}], 1];
edgep0r = Flatten[Table[{2 pp, 0, 2 rp + 1}, {pp, 1, QMaxPrime},
    {rp, 0, RMax}], 1];
edgeIndices = Join[edge0qr, edgep0r];
nInt = Length[edgeIndices];
Print["Edge mode indices (Q_max'=", QMaxPrime, ", R_max=", RMax, "):"];
Print["  (0, q, r) family: ", edge0qr];
Print["  (p, 0, r) family: ", edgep0r];
Print["  total nInt = ", nInt];
Print[""];

(* === Interior matrices (closed form) === *)
(* For (0, q, r): M_ii = (V/4)(q pi/b)^2, K_ii = lambda M_ii where
   lambda = pi^2((q/b)^2 + (r/c)^2). b_i = -V/(r pi). *)
edgeMassDiag[idx_] := Module[{p, q, r, lam},
  {p, q, r} = idx;
  lam = If[p == 0, (q PiN/bN)^2, (p PiN/aN)^2];
  (VN/4) lam];

edgeStiffDiag[idx_] := Module[{p, q, r, lam, lamtot},
  {p, q, r} = idx;
  lam = If[p == 0, (q PiN/bN)^2, (p PiN/aN)^2];
  lamtot = PiN^2 ((p/aN)^2 + (q/bN)^2 + (r/cN)^2);
  (VN/4) lam * lamtot];

edgeSrc[idx_] := Module[{p, q, r}, {p, q, r} = idx; -VN/(r PiN)];

MIIDiag = Table[edgeMassDiag[idx], {idx, edgeIndices}];
KIIDiag = Table[edgeStiffDiag[idx], {idx, edgeIndices}];
bIVec = Table[edgeSrc[idx], {idx, edgeIndices}];

Print["Interior diagonal matrices (closed form, 60-digit):"];
Print["  M^II first 3: ", N[MIIDiag[[1 ;; 3]], 12]];
Print["  K^II first 3: ", N[KIIDiag[[1 ;; 3]], 12]];
Print["  b^I first 3:  ", N[bIVec[[1 ;; 3]], 12]];
Print[""];

(* Interior-only alpha_2 verification *)
alpha2IntOnly = (mu0N sigmaN)^2/VN * Sum[bIVec[[i]]^2/KIIDiag[[i]], {i, nInt}];
Print["alpha_2 interior PEC only (truncated to ", nInt, " edge modes): ",
  N[alpha2IntOnly, 16]];
Print["  (Stage 5b infinite series gave: 1.674e-10 -- this truncated value should be smaller)"];
Print[""];

(* === Exterior basis === *)
R0 = SetPrecision[1/1000, prec];  (* 1 mm *)
KEE = 15 PiN^2 / (32 R0^3);
Print["Exterior basis Phi^E = (-y, x, 0)/(r^2+R_0^2)^(3/2), R_0 = ", N[R0*1000, 6], " mm"];
Print["  K^EE (closed form) = 15 pi^2/(32 R_0^3) = ", N[KEE, 16]];
Print[""];

(* === NIntegrate setup === *)
wp = 50; ag = 35; pg = 30;
mr = 12;

(* Interior basis explicit forms (curl A^I and A^I) *)
AIxTE[p_, q_, r_, xx_, yy_, zz_] := -(q PiN/bN) Cos[p PiN (xx + aN/2)/aN] *
    Sin[q PiN (yy + bN/2)/bN] Sin[r PiN (zz + cN/2)/cN];
AIyTE[p_, q_, r_, xx_, yy_, zz_] := (p PiN/aN) Sin[p PiN (xx + aN/2)/aN] *
    Cos[q PiN (yy + bN/2)/bN] Sin[r PiN (zz + cN/2)/cN];
(* curl A^I *)
curlAIxTE[p_, q_, r_, xx_, yy_, zz_] := -(p PiN/aN)(r PiN/cN) Sin[p PiN (xx + aN/2)/aN] *
    Cos[q PiN (yy + bN/2)/bN] Cos[r PiN (zz + cN/2)/cN];
curlAIyTE[p_, q_, r_, xx_, yy_, zz_] := -(q PiN/bN)(r PiN/cN) Cos[p PiN (xx + aN/2)/aN] *
    Sin[q PiN (yy + bN/2)/bN] Cos[r PiN (zz + cN/2)/cN];
curlAIzTE[p_, q_, r_, xx_, yy_, zz_] := -((p PiN/aN)^2 + (q PiN/bN)^2) *
    Cos[p PiN (xx + aN/2)/aN] Cos[q PiN (yy + bN/2)/bN] Sin[r PiN (zz + cN/2)/cN];

(* curl Phi^E *)
curlPhiEx[xx_, yy_, zz_] := 3 xx zz/(xx^2 + yy^2 + zz^2 + R0^2)^(5/2);
curlPhiEy[xx_, yy_, zz_] := 3 yy zz/(xx^2 + yy^2 + zz^2 + R0^2)^(5/2);
curlPhiEz[xx_, yy_, zz_] := (2 (zz^2 + R0^2) - (xx^2 + yy^2))/
    (xx^2 + yy^2 + zz^2 + R0^2)^(5/2);

(* Phi^E *)
PhiEx[xx_, yy_, zz_] := -yy/(xx^2 + yy^2 + zz^2 + R0^2)^(3/2);
PhiEy[xx_, yy_, zz_] := xx/(xx^2 + yy^2 + zz^2 + R0^2)^(3/2);

(* === Compute NIntegrate matrices === *)
Print["Computing M^EE (1 NIntegrate, expect ~30s)..."];
tStart = AbsoluteTime[];
MEE = NIntegrate[
  (xx^2 + yy^2)/(xx^2 + yy^2 + zz^2 + R0^2)^3,
  {xx, -aN/2, aN/2}, {yy, -bN/2, bN/2}, {zz, -cN/2, cN/2},
  WorkingPrecision -> wp, AccuracyGoal -> ag, PrecisionGoal -> pg,
  MaxRecursion -> mr];
Print["  M^EE = ", N[MEE, 16], "  (", N[AbsoluteTime[] - tStart, 4], " s)"];

Print["Computing b^E (1 NIntegrate)..."];
tStart = AbsoluteTime[];
bE = (1/2) NIntegrate[
  (xx^2 + yy^2)/(xx^2 + yy^2 + zz^2 + R0^2)^(3/2),
  {xx, -aN/2, aN/2}, {yy, -bN/2, bN/2}, {zz, -cN/2, cN/2},
  WorkingPrecision -> wp, AccuracyGoal -> ag, PrecisionGoal -> pg,
  MaxRecursion -> mr];
Print["  b^E = ", N[bE, 16], "  (", N[AbsoluteTime[] - tStart, 4], " s)"];

Print["Computing M^IE (", nInt, " NIntegrates)..."];
tStart = AbsoluteTime[];
MIEVec = Table[Module[{p, q, r, integrand, val},
  {p, q, r} = edgeIndices[[i]];
  integrand = (AIxTE[p, q, r, xx, yy, zz] PhiEx[xx, yy, zz] +
               AIyTE[p, q, r, xx, yy, zz] PhiEy[xx, yy, zz]);
  val = NIntegrate[integrand,
    {xx, -aN/2, aN/2}, {yy, -bN/2, bN/2}, {zz, -cN/2, cN/2},
    WorkingPrecision -> wp, AccuracyGoal -> ag, PrecisionGoal -> pg,
    MaxRecursion -> mr];
  Print["  M^IE[", edgeIndices[[i]], "] = ", N[val, 12]];
  val],
  {i, nInt}];
Print["  M^IE total time: ", N[AbsoluteTime[] - tStart, 4], " s"];

Print["Computing K^IE (", nInt, " NIntegrates)..."];
tStart = AbsoluteTime[];
KIEVec = Table[Module[{p, q, r, integrand, val},
  {p, q, r} = edgeIndices[[i]];
  integrand = (curlAIxTE[p, q, r, xx, yy, zz] curlPhiEx[xx, yy, zz] +
               curlAIyTE[p, q, r, xx, yy, zz] curlPhiEy[xx, yy, zz] +
               curlAIzTE[p, q, r, xx, yy, zz] curlPhiEz[xx, yy, zz]);
  val = NIntegrate[integrand,
    {xx, -aN/2, aN/2}, {yy, -bN/2, bN/2}, {zz, -cN/2, cN/2},
    WorkingPrecision -> wp, AccuracyGoal -> ag, PrecisionGoal -> pg,
    MaxRecursion -> mr];
  Print["  K^IE[", edgeIndices[[i]], "] = ", N[val, 12]];
  val],
  {i, nInt}];
Print["  K^IE total time: ", N[AbsoluteTime[] - tStart, 4], " s"];
Print[""];

(* === Build augmented system === *)
nTot = nInt + 1;
MAug = Table[0, {nTot}, {nTot}];
KAug = Table[0, {nTot}, {nTot}];
bAug = Join[bIVec, {bE}];
Do[MAug[[i, i]] = MIIDiag[[i]]; KAug[[i, i]] = KIIDiag[[i]], {i, nInt}];
Do[MAug[[i, nTot]] = MIEVec[[i]]; MAug[[nTot, i]] = MIEVec[[i]];
   KAug[[i, nTot]] = KIEVec[[i]]; KAug[[nTot, i]] = KIEVec[[i]],
  {i, nInt}];
MAug[[nTot, nTot]] = MEE;
KAug[[nTot, nTot]] = KEE;

KAugInv = Inverse[KAug];

(* Sanity check: K should be symmetric positive definite *)
Print["K matrix Schur complement (>0 for PD):"];
schur = KEE - Sum[KIEVec[[i]]^2/KIIDiag[[i]], {i, nInt}];
Print["  K^EE - K^EI K^II^(-1) K^IE = ", N[schur, 12]];
If[N[schur, 12] > 0, Print["  -> K augmented is positive definite ✓"],
  Print["  -> K augmented is NOT positive definite (problem)"]];
Print[""];

(* alpha_n via series *)
alphaN[n_] := (-1)^n (mu0N sigmaN)^n/VN *
  bAug.Nest[# . MAug . KAugInv &, KAugInv, n - 2].bAug;

Print["=== Augmented Galerkin alpha_n ==="];
Do[Module[{val},
  val = alphaN[k];
  Print["alpha_", k, " = ", N[val, 16],
    "  (sphere ref: ", N[alphaSphere[Sqrt[ReffSq], k], 12],
    ", ratio: ", N[val/alphaSphere[Sqrt[ReffSq], k], 8], ")"]],
  {k, 2, 4}];
Print[""];

Print["=== Improvement summary ==="];
Print["alpha_2 interior only (truncated):  ", N[alpha2IntOnly, 16]];
Print["alpha_2 enriched (interior+dipole): ", N[alphaN[2], 16]];
Print["alpha_2 sphere R_eff:               ",
  N[alphaSphere[Sqrt[ReffSq], 2], 16]];
Print[""];
Print["Improvement ratio (enriched/interior_only): ",
  N[alphaN[2]/alpha2IntOnly, 8]];
Print["Gap to sphere ref (enriched/sphere):       ",
  N[alphaN[2]/alphaSphere[Sqrt[ReffSq], 2], 8]];
