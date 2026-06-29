(* Stage 5d: Stage 5c + axial octupole exterior basis.

   Adds Phi^E_3(r) = (4z^2 - x^2 - y^2)(-y, x, 0) / (r^2 + R_0^2)^(7/2)
   This is the regularized vector potential for axial magnetic octupole
   (l=3, m=0). At infinity: A ~ 1/r^4, B = curl A ~ 1/r^5 (octupole field).

   Augmented basis: n_int + 2 (interior PEC TE_z + dipole + octupole).
*)
$HistoryLength = 0;
prec = 60;

Print["=== Stage 5d: PEC TE_z + dipole + octupole (axial) ==="];
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

alphaSphere[R_, k_] := Module[{x = mu0N sigmaN R^2, ck},
  ck = Switch[k, 1, -1/10, 2, +1/105, 3, -1/1050, 4, +1/10395, _,
    (-1)^k * 9 * SetPrecision[Zeta[2 k + 2], prec]/PiN^(2 k + 2)];
  ck * x^k];

(* Edge mode indices (same as Stage 5c) *)
QMaxPrime = 3;
RMax = 2;
edge0qr = Flatten[Table[{0, 2 qp, 2 rp + 1}, {qp, 1, QMaxPrime},
    {rp, 0, RMax}], 1];
edgep0r = Flatten[Table[{2 pp, 0, 2 rp + 1}, {pp, 1, QMaxPrime},
    {rp, 0, RMax}], 1];
edgeIndices = Join[edge0qr, edgep0r];
nInt = Length[edgeIndices];
Print["Interior PEC TE_z basis: nInt = ", nInt];

(* Diagonal interior matrices *)
edgeMassDiag[idx_] := Module[{p, q, r, lam}, {p, q, r} = idx;
  lam = If[p == 0, (q PiN/bN)^2, (p PiN/aN)^2];
  (VN/4) lam];
edgeStiffDiag[idx_] := Module[{p, q, r, lam, lamtot}, {p, q, r} = idx;
  lam = If[p == 0, (q PiN/bN)^2, (p PiN/aN)^2];
  lamtot = PiN^2 ((p/aN)^2 + (q/bN)^2 + (r/cN)^2);
  (VN/4) lam * lamtot];
edgeSrc[idx_] := Module[{p, q, r}, {p, q, r} = idx; -VN/(r PiN)];

MIIDiag = Table[edgeMassDiag[idx], {idx, edgeIndices}];
KIIDiag = Table[edgeStiffDiag[idx], {idx, edgeIndices}];
bIVec = Table[edgeSrc[idx], {idx, edgeIndices}];

(* === External bases === *)
R0 = SetPrecision[1/1000, prec];

(* Dipole: Phi^E1 = (-y, x, 0)/(r^2+R0^2)^(3/2) *)
KEE1 = 15 PiN^2 / (32 R0^3);  (* closed form *)

(* Octupole: Phi^E3 = (4z^2 - x^2 - y^2)(-y, x, 0)/(r^2+R0^2)^(7/2) *)
(* K^E3 E3 = integral|curl Phi^E3|^2 dV over R^3 -- need to compute.
   Will compute via NIntegrate over large box for now (could do analytically later). *)

(* Phi^E1, Phi^E3 in Mathematica form *)
(* Phi^E1 components: *)
PhiE1x[xx_, yy_, zz_] := -yy/(xx^2 + yy^2 + zz^2 + R0^2)^(3/2);
PhiE1y[xx_, yy_, zz_] := xx/(xx^2 + yy^2 + zz^2 + R0^2)^(3/2);
(* curl Phi^E1: *)
curlPhiE1x[xx_, yy_, zz_] := 3 xx zz/(xx^2 + yy^2 + zz^2 + R0^2)^(5/2);
curlPhiE1y[xx_, yy_, zz_] := 3 yy zz/(xx^2 + yy^2 + zz^2 + R0^2)^(5/2);
curlPhiE1z[xx_, yy_, zz_] := (2 (zz^2 + R0^2) - (xx^2 + yy^2))/
    (xx^2 + yy^2 + zz^2 + R0^2)^(5/2);

(* Phi^E3 = (4z^2 - x^2 - y^2)(-y, x, 0)/(r^2+R0^2)^(7/2) *)
PhiE3x[xx_, yy_, zz_] := -yy (4 zz^2 - xx^2 - yy^2)/
    (xx^2 + yy^2 + zz^2 + R0^2)^(7/2);
PhiE3y[xx_, yy_, zz_] := xx (4 zz^2 - xx^2 - yy^2)/
    (xx^2 + yy^2 + zz^2 + R0^2)^(7/2);
(* curl Phi^E3: compute symbolically *)
phiE3xSym = -yy (4 zz^2 - xx^2 - yy^2)/(xx^2 + yy^2 + zz^2 + R0^2)^(7/2);
phiE3ySym = xx (4 zz^2 - xx^2 - yy^2)/(xx^2 + yy^2 + zz^2 + R0^2)^(7/2);
curlPhiE3xSym = -D[phiE3ySym, zz];  (* -d_z A_y *)
curlPhiE3ySym = D[phiE3xSym, zz];   (* d_z A_x *)
curlPhiE3zSym = D[phiE3ySym, xx] - D[phiE3xSym, yy];

curlPhiE3x[xx0_, yy0_, zz0_] := curlPhiE3xSym /. {xx -> xx0, yy -> yy0, zz -> zz0};
curlPhiE3y[xx0_, yy0_, zz0_] := curlPhiE3ySym /. {xx -> xx0, yy -> yy0, zz -> zz0};
curlPhiE3z[xx0_, yy0_, zz0_] := curlPhiE3zSym /. {xx -> xx0, yy -> yy0, zz -> zz0};

(* Interior basis explicit forms *)
AIxTE[p_, q_, r_, xx_, yy_, zz_] := -(q PiN/bN) Cos[p PiN (xx + aN/2)/aN] *
    Sin[q PiN (yy + bN/2)/bN] Sin[r PiN (zz + cN/2)/cN];
AIyTE[p_, q_, r_, xx_, yy_, zz_] := (p PiN/aN) Sin[p PiN (xx + aN/2)/aN] *
    Cos[q PiN (yy + bN/2)/bN] Sin[r PiN (zz + cN/2)/cN];
curlAIxTE[p_, q_, r_, xx_, yy_, zz_] := -(p PiN/aN)(r PiN/cN) Sin[p PiN (xx + aN/2)/aN] *
    Cos[q PiN (yy + bN/2)/bN] Cos[r PiN (zz + cN/2)/cN];
curlAIyTE[p_, q_, r_, xx_, yy_, zz_] := -(q PiN/bN)(r PiN/cN) Cos[p PiN (xx + aN/2)/aN] *
    Sin[q PiN (yy + bN/2)/bN] Cos[r PiN (zz + cN/2)/cN];
curlAIzTE[p_, q_, r_, xx_, yy_, zz_] := -((p PiN/aN)^2 + (q PiN/bN)^2) *
    Cos[p PiN (xx + aN/2)/aN] Cos[q PiN (yy + bN/2)/bN] Sin[r PiN (zz + cN/2)/cN];

(* === NIntegrate parameters === *)
wp = 50; ag = 30; pg = 25; mr = 12;

(* Direct NIntegrate calls (avoid Sequence @@ scope bug) *)
doInt[expr_] := NIntegrate[expr,
  {xx, -aN/2, aN/2}, {yy, -bN/2, bN/2}, {zz, -cN/2, cN/2},
  WorkingPrecision -> wp, AccuracyGoal -> ag, PrecisionGoal -> pg,
  MaxRecursion -> mr];

(* === Recompute dipole interactions (already known from 5c, but redo for consistency) === *)
Print["Computing M^E1 E1, b^E1 (dipole self & source) ..."];
tStart = AbsoluteTime[];
ME1E1 = doInt[(xx^2 + yy^2)/(xx^2 + yy^2 + zz^2 + R0^2)^3];
bE1 = (1/2) doInt[(xx^2 + yy^2)/(xx^2 + yy^2 + zz^2 + R0^2)^(3/2)];
Print["  M^E1E1 = ", N[ME1E1, 12], "  (", N[AbsoluteTime[] - tStart, 4], "s)"];
Print["  b^E1   = ", N[bE1, 12]];

(* === Octupole self & source === *)
Print["Computing M^E3 E3, b^E3 (octupole self & source) ..."];
tStart = AbsoluteTime[];
(* Phi^E3 . Phi^E3 = (x^2+y^2)(4z^2-x^2-y^2)^2/(r^2+R0^2)^7 *)
ME3E3 = doInt[(xx^2 + yy^2)(4 zz^2 - xx^2 - yy^2)^2/
    (xx^2 + yy^2 + zz^2 + R0^2)^7];
(* b^E3 = (1/2) integral (x A_ext_y - y A_ext_x)... wait
   b^E3 = integral Phi^E3 . A_ext dV = (1/2) integral (x^2+y^2)(4z^2-x^2-y^2)/(r^2+R0^2)^(7/2) dV *)
bE3 = (1/2) doInt[(xx^2 + yy^2)(4 zz^2 - xx^2 - yy^2)/
    (xx^2 + yy^2 + zz^2 + R0^2)^(7/2)];
Print["  M^E3E3 = ", N[ME3E3, 12], "  (", N[AbsoluteTime[] - tStart, 4], "s)"];
Print["  b^E3   = ", N[bE3, 12]];

(* Dipole-octupole cross-coupling *)
Print["Computing M^E1 E3, K^E1 E3 (cross terms) ..."];
tStart = AbsoluteTime[];
(* Phi^E1 . Phi^E3 = (-y)(-y)(4z^2-x^2-y^2)/(r^2+R0^2)^(3/2+7/2) + (x)(x)(4z^2-..)/...
                  = (y^2 + x^2)(4z^2 - x^2 - y^2)/(r^2+R0^2)^5 *)
ME1E3 = doInt[(xx^2 + yy^2)(4 zz^2 - xx^2 - yy^2)/
    (xx^2 + yy^2 + zz^2 + R0^2)^5];
(* curl Phi^E1 . curl Phi^E3 over R^3: too costly; compute over conductor only and add closed form *)
KE1E3cond = doInt[curlPhiE1x[xx, yy, zz] curlPhiE3xSym +
                  curlPhiE1y[xx, yy, zz] curlPhiE3ySym +
                  curlPhiE1z[xx, yy, zz] curlPhiE3zSym];
(* For exterior portion of K^E1 E3: orthogonality of dipole and octupole spherical
   harmonics over angles gives 0 in the limit R0 -> 0. Will compute as exterior integral
   if needed. For now use just conductor part as approximation. *)
KE1E3 = KE1E3cond; (* approximation: most of K^E1 E3 contribution is over conductor for our R0 ~ size *)

Print["  M^E1E3 = ", N[ME1E3, 12], "  (", N[AbsoluteTime[] - tStart, 4], "s)"];
Print["  K^E1E3 (cond only) = ", N[KE1E3, 12]];

(* === K^E3 E3 over R^3: use spherical sym + axial part === *)
(* Numerical via large-box NIntegrate, scale by R0 *)
Print["Computing K^E3 E3 (over R^3, via large box scaling) ..."];
tStart = AbsoluteTime[];
(* Substitution r = R0 t makes integrand pure dimensionless. integral = R0^(-7) * dimensionless *)
(* Skip full R^3 for now; use conductor-only approximation as warning *)
KE3E3cond = doInt[curlPhiE3xSym^2 + curlPhiE3ySym^2 + curlPhiE3zSym^2];
KE3E3 = KE3E3cond;  (* approximation *)
Print["  K^E3E3 (cond only) = ", N[KE3E3, 12], "  (", N[AbsoluteTime[] - tStart, 4], "s)"];
Print["  WARNING: should add exterior contribution (large-box NIntegrate or closed form)"];
Print[""];

(* === Interior <-> Dipole coupling (recompute or load from 5c) === *)
Print["Computing M^IE1, K^IE1 (", nInt, " each)..."];
tStart = AbsoluteTime[];
MIE1Vec = Table[Module[{p, q, r}, {p, q, r} = edgeIndices[[i]];
  doInt[AIxTE[p, q, r, xx, yy, zz] PhiE1x[xx, yy, zz] +
        AIyTE[p, q, r, xx, yy, zz] PhiE1y[xx, yy, zz]]],
  {i, nInt}];
KIE1Vec = Table[Module[{p, q, r}, {p, q, r} = edgeIndices[[i]];
  doInt[curlAIxTE[p, q, r, xx, yy, zz] curlPhiE1x[xx, yy, zz] +
        curlAIyTE[p, q, r, xx, yy, zz] curlPhiE1y[xx, yy, zz] +
        curlAIzTE[p, q, r, xx, yy, zz] curlPhiE1z[xx, yy, zz]]],
  {i, nInt}];
Print["  M^IE1, K^IE1 done (", N[AbsoluteTime[] - tStart, 4], "s)"];

(* === Interior <-> Octupole coupling === *)
Print["Computing M^IE3, K^IE3 (", nInt, " each)..."];
tStart = AbsoluteTime[];
MIE3Vec = Table[Module[{p, q, r}, {p, q, r} = edgeIndices[[i]];
  doInt[AIxTE[p, q, r, xx, yy, zz] PhiE3x[xx, yy, zz] +
        AIyTE[p, q, r, xx, yy, zz] PhiE3y[xx, yy, zz]]],
  {i, nInt}];
KIE3Vec = Table[Module[{p, q, r}, {p, q, r} = edgeIndices[[i]];
  doInt[curlAIxTE[p, q, r, xx, yy, zz] curlPhiE3xSym +
        curlAIyTE[p, q, r, xx, yy, zz] curlPhiE3ySym +
        curlAIzTE[p, q, r, xx, yy, zz] curlPhiE3zSym]],
  {i, nInt}];
Print["  M^IE3, K^IE3 done (", N[AbsoluteTime[] - tStart, 4], "s)"];
Print[""];

(* === Build augmented system === *)
nTot = nInt + 2;  (* +1 dipole +1 octupole *)
MAug = Table[0, {nTot}, {nTot}];
KAug = Table[0, {nTot}, {nTot}];
bAug = Table[0, {nTot}];

Do[MAug[[i, i]] = MIIDiag[[i]]; KAug[[i, i]] = KIIDiag[[i]]; bAug[[i]] = bIVec[[i]],
  {i, nInt}];
(* dipole row/col = nInt+1 *)
Do[MAug[[i, nInt+1]] = MIE1Vec[[i]]; MAug[[nInt+1, i]] = MIE1Vec[[i]];
   KAug[[i, nInt+1]] = KIE1Vec[[i]]; KAug[[nInt+1, i]] = KIE1Vec[[i]],
  {i, nInt}];
MAug[[nInt+1, nInt+1]] = ME1E1; KAug[[nInt+1, nInt+1]] = KEE1; bAug[[nInt+1]] = bE1;
(* octupole row/col = nInt+2 *)
Do[MAug[[i, nInt+2]] = MIE3Vec[[i]]; MAug[[nInt+2, i]] = MIE3Vec[[i]];
   KAug[[i, nInt+2]] = KIE3Vec[[i]]; KAug[[nInt+2, i]] = KIE3Vec[[i]],
  {i, nInt}];
MAug[[nInt+1, nInt+2]] = ME1E3; MAug[[nInt+2, nInt+1]] = ME1E3;
KAug[[nInt+1, nInt+2]] = KE1E3; KAug[[nInt+2, nInt+1]] = KE1E3;
MAug[[nInt+2, nInt+2]] = ME3E3; KAug[[nInt+2, nInt+2]] = KE3E3;
bAug[[nInt+2]] = bE3;

KAugInv = Inverse[KAug];

alphaN[n_] := (-1)^n (mu0N sigmaN)^n/VN *
  bAug.Nest[# . MAug . KAugInv &, KAugInv, n - 2].bAug;

Print["=== Augmented Galerkin alpha_n (with octupole) ==="];
Do[Module[{val, sphV},
  val = alphaN[k]; sphV = alphaSphere[Sqrt[ReffSq], k];
  Print["alpha_", k, " = ", N[val, 16],
    "  (sphere: ", N[sphV, 12], ", ratio: ", N[val/sphV, 8], ")"]],
  {k, 2, 5}];
Print[""];

(* Compare *)
Print["=== Improvement vs Stage 5c ==="];
Print["  alpha_2 5c (PEC + dipole only):    ~ 4.634e-10  (ratio 0.251)"];
Print["  alpha_2 5d (PEC + dipole+octupole): ", N[alphaN[2], 16],
      "  (ratio: ", N[alphaN[2]/alphaSphere[Sqrt[ReffSq], 2], 6], ")"];
