(* Stage 5 (B): Extract more moments from Stage 5c augmented matrices.

   Reconstruct K_aug, M_aug, b_aug from Stage 5c output (NIntegrate values
   at 12-digit precision, no extra integration needed).
   Then compute alpha_n for n = 2..15 and Pade-extract multiple Foster pairs.

   Limitations:
     - 12-digit input precision -> ~10 reliable Pade pairs at most
     - 9-mode basis (8 interior + 1 dipole) -> physical info saturates ~3-4 pairs
*)
$HistoryLength = 0;
prec = 60;

Print["=== Stage 5 (B): Extended moments from Stage 5c data ==="];
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

(* Edge indices same as Stage 5c *)
QMaxPrime = 3; RMax = 2;
edge0qr = Flatten[Table[{0, 2 qp, 2 rp + 1}, {qp, 1, QMaxPrime},
    {rp, 0, RMax}], 1];
edgep0r = Flatten[Table[{2 pp, 0, 2 rp + 1}, {pp, 1, QMaxPrime},
    {rp, 0, RMax}], 1];
edgeIndices = Join[edge0qr, edgep0r];
nInt = Length[edgeIndices];

(* Closed-form interior matrix entries *)
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

(* Hardcoded Stage 5c NIntegrate values (12-digit precision) *)
R0 = SetPrecision[1/1000, prec];
KEE = 15 PiN^2/(32 R0^3);  (* closed form *)
MEE = SetPrecision[687.93020035713306087260104809528524595982, prec];
bE = SetPrecision[1.5537698392973272334753285573472*^-6, prec];

(* M^IE values per edge index (order matches edgeIndices) *)
MIEData = SetPrecision[{
  -1.64128315536591046026960833728223305143,    (* {0,2,1} *)
  -0.47159780549998160817930691474725126699,    (* {0,2,3} *)
  -0.28240304832860431622792094953520901138,    (* {0,2,5} *)
  -1.06418360685479744539883271963967350656,    (* {0,4,1} *)
  -0.32731682906247031267939278713996206915,    (* {0,4,3} *)
  -0.19529056952797727052792829226217485677,    (* {0,4,5} *)
  -1.09517102291662498071043222480821315278,    (* {0,6,1} *)
  -0.32982827278791484993278259713796577571,    (* {0,6,3} *)
  -0.19719905946178936300919896235421951032,    (* {0,6,5} *)
  -1.16143879312752107610234014819494701172,    (* {2,0,1} *)
  -0.34639104797515909020391778485781075838,    (* {2,0,3} *)
  -0.20717989740343753609722445028735892963,    (* {2,0,5} *)
   0.35851032349081182083067554826992061849,    (* {4,0,1} *)
   0.07413031476517436088665772482743094875,    (* {4,0,3} *)
   0.04460854306397915971593006742013012527,    (* {4,0,5} *)
  -0.58612432191866345889248446878654704834,    (* {6,0,1} *)
  -0.16087879360116839208529116839737833662,    (* {6,0,3} *)
  -0.09689788551996370862094499432418230440     (* {6,0,5} *)
  }, prec];

KIEData = SetPrecision[{
   3.3655480957527208726113301*^6,
   1.319872829207292384460137*^6,
   760350.7471501474983941067,
  -1.4264720420822998825925008*^6,
  -133980.8712952001968432643,
  -80832.9981976664330932619,
   324348.5249799661441201139,
   256302.5926924511939446871,
   149840.9695090602836232325,
   799162.7092004847737178193,
   435129.6976967181811553122,
   250249.1895948848190981751,
  -3.3190815872779968601463487*^6,
  -1.2090838683928092213340739*^6,
  -697648.3861909722151371335,
   3.360539108048216242806438*^6,
   1.043735846299343694576735*^6,
   600722.3358986344421484987
  }, prec];

(* Build augmented system *)
nTot = nInt + 1;
MAug = Table[0, {nTot}, {nTot}];
KAug = Table[0, {nTot}, {nTot}];
Do[MAug[[i, i]] = MIIDiag[[i]]; KAug[[i, i]] = KIIDiag[[i]], {i, nInt}];
Do[MAug[[i, nTot]] = MIEData[[i]]; MAug[[nTot, i]] = MIEData[[i]];
   KAug[[i, nTot]] = KIEData[[i]]; KAug[[nTot, i]] = KIEData[[i]],
  {i, nInt}];
MAug[[nTot, nTot]] = MEE;
KAug[[nTot, nTot]] = KEE;
bAug = Append[bIVec, bE];

KAugInv = Inverse[KAug];

alphaLFRef = -mu0N sigmaN (aN^2 + bN^2)/48;
alphaSphere[R_, k_] := Module[{x = mu0N sigmaN R^2, ck},
  ck = Switch[k, 1, -1/10, 2, +1/105, 3, -1/1050, 4, +1/10395, 5, -1/135135,
    6, +1/1351350, 7, -1/14148810, _,
    (-1)^k * 9 * SetPrecision[Zeta[2 k + 2], prec]/PiN^(2 k + 2)];
  ck * x^k];

(* Compute alpha_2 ... alpha_NMax *)
NMaxAlpha = 15;
alphaN[n_] := (-1)^n (mu0N sigmaN)^n/VN *
  bAug.Nest[# . MAug . KAugInv &, KAugInv, n - 2].bAug;

Print["=== alpha_n from Stage 5c augmented matrices (n=2..", NMaxAlpha, ") ==="];
Print[""];
alphaList = {alphaLFRef};  (* alpha_1 *)
Do[AppendTo[alphaList, alphaN[n]], {n, 2, NMaxAlpha}];

Print[StringPadRight["n", 4], StringPadRight["alpha_n (Galerkin)", 30],
  StringPadRight["alpha_n^sph (R_eff)", 28], "ratio"];
Print[StringJoin[Table["-", {90}]]];
Do[Module[{aGal, aSph, ratio},
  aGal = alphaList[[k]];
  aSph = If[k == 1, alphaLFRef, alphaSphere[Sqrt[ReffSq], k]];
  ratio = aGal/aSph;
  Print[StringPadRight[ToString[k], 4],
    StringPadRight[ToString[CForm[N[aGal, 14]]], 30],
    StringPadRight[ToString[CForm[N[aSph, 14]]], 28],
    ToString[NumberForm[N[ratio, 8], 8]]]],
  {k, 1, NMaxAlpha}];
Print[""];

(* === Convert to Foster pairs via Pade =====================================
   With moments M_k = (-1)^k * coefficient of s^(k+1) in chi(s)... actually:
   Y(s)/s = sum_n a_n/(1+s tau_n)
   Y/s Taylor: M_0 - M_1 s + M_2 s^2 - M_3 s^3 + ...   (M_k = sum a_n tau_n^k)

   chi(s) = -Y(s) = -alpha_1 s - alpha_2 s^2 - ...
   chi/s = -alpha_1 - alpha_2 s - alpha_3 s^2 - ...
   Y/s = -chi/s = alpha_1 + alpha_2 s + alpha_3 s^2 + ... (with our alpha_1 < 0!)

   Wait -- our alpha_1 < 0 but Y/s should be positive at s=0. Let me re-derive:
   Y = -chi. chi(s) -> 0 as s -> 0 from below (alpha_1 < 0 means slope is negative,
   chi decreases). So Y = -chi increases from 0. dY/ds at 0 = -alpha_1 > 0. Y > 0 for small s.
   Y/s at s=0+ = -alpha_1 > 0. ✓

   So Y(s)/s Taylor:
   Y/s = (-alpha_1) + (-alpha_2)*s + (-alpha_3)*s^2 + (-alpha_4)*s^3 + ...
   Match (-1)^k M_k:
     k=0: M_0 = -alpha_1
     k=1: -M_1 = -alpha_2  =>  M_1 = alpha_2
     k=2: +M_2 = -alpha_3  =>  M_2 = -alpha_3
     k=3: -M_3 = -alpha_4  =>  M_3 = alpha_4
     k=4: +M_4 = -alpha_5  =>  M_4 = -alpha_5
   General: M_k = (-1)^k (-alpha_(k+1)) = (-1)^(k+1) alpha_(k+1)

   For our Galerkin signs: alpha_1<0, alpha_2>0, alpha_3<0, alpha_4>0, alpha_5<0,...
   So M_k = (-1)^(k+1) alpha_(k+1) ALL POSITIVE. ✓ Stieltjes condition satisfied. *)

momList = Table[(-1)^(k + 1) alphaList[[k + 1]], {k, 0, NMaxAlpha - 1}];
Print["Foster moments M_k (k=0..", NMaxAlpha - 1, "):"];
Do[Print["  M_", k - 1, " = ", N[momList[[k]], 16],
  If[N[momList[[k]], 30] > 0, " ✓positive", " ✗NEGATIVE!"]],
  {k, 1, Length[momList]}];
Print[""];

(* Pade [n,n] from 2n moments via Hankel diagonalization *)
extractFoster[mom_, n_] := Module[{H1, H2, tauVals, aVals, V},
  H1 = Table[mom[[i + j - 1]], {i, n}, {j, n}];
  H2 = Table[mom[[i + j]], {i, n}, {j, n}];
  tauVals = Eigenvalues[Inverse[H1].H2];
  V = Table[tauVals[[j]]^(i - 1), {i, n}, {j, n}];
  aVals = LinearSolve[V, mom[[1 ;; n]]];
  {tauVals, aVals}];

Print["=== Pade Foster extraction (multiple n) ==="];
Print[""];
Do[Module[{tauL, aL, sortedPairs},
  {tauL, aL} = extractFoster[momList, n];
  sortedPairs = SortBy[Transpose[{tauL, aL}], -#[[1]] &];
  Print["Pade [", n, ",", n, "] (", n, " Foster pairs from ", 2 n, " moments):"];
  Do[Print["  tau_", j, " = ", N[sortedPairs[[j, 1]] * 10^6, 12], " us, ",
    "a_", j, " = ", N[sortedPairs[[j, 2]], 12]],
    {j, 1, Length[sortedPairs]}];
  Print[""]],
  {n, 2, Min[Floor[NMaxAlpha/2], 7]}];

(* Compare to sphere R_eff Foster *)
Print["=== Sphere R_eff Foster reference ==="];
Print["  tau_n^sph = mu0 sigma R_eff^2 / (n pi)^2"];
Reff = Sqrt[ReffSq];
sphTaus = Table[mu0N sigmaN Reff^2/(n PiN)^2, {n, 1, 7}];
sphAs = Table[9/(n PiN)^2, {n, 1, 7}];  (* approximate, leading order *)
Do[Print["  tau_", n, " = ", N[sphTaus[[n]] * 10^6, 12], " us, a_", n, " = ",
  N[sphAs[[n]] * mu0N sigmaN Reff^2, 12], " (approx)"],
  {n, 1, 7}];
Print[""];

(* Cauer-I from Foster pairs (analytical recursion) *)
Print["=== Cauer-I ladder from Pade [4,4] (8 Foster pairs from 8 moments) ==="];
Print["[abandoned: 8 moments gives 4 pairs not 8 -- use Pade [4,4]]"];

(* Use Pade [4,4] -- 4 Foster pairs from 8 moments *)
{tauP4, aP4} = extractFoster[momList, 4];
pairsP4 = SortBy[Transpose[{tauP4, aP4}], -#[[1]] &];
Print["4 Foster pairs (Pade [4,4]):"];
Do[Print["  tau_", j, " = ", N[pairsP4[[j, 1]] * 10^6, 12], " us, ",
  "a_", j, " = ", N[pairsP4[[j, 2]], 12]],
  {j, 1, 4}];
Print[""];

(* Convert Foster pairs to Cauer-I via continued-fraction expansion *)
(* Y(s)/s = sum a_n/(1 + s tau_n) -- as rational function in s *)
Print["=== Cauer-I conversion ==="];
yOverS = Sum[pairsP4[[n, 2]]/(1 + sv pairsP4[[n, 1]]), {n, 1, 4}];
yOverS = Together[yOverS];
yFn = sv * yOverS;  (* Y(s) = s * (Y/s) *)
zFn = -1/yFn;  (* Z = -1/Y = -1/-chi = 1/chi -- wait Y = -chi so 1/Y = -1/chi *)
(* Z(s) for ladder extraction *)
Print["chi(s) (4-pole rational from Pade):"];
Print["  ", -yFn // Together];

(* Cauer-I extraction: Z(s) = sL_0 + 1/(G_1 + 1/(sL_1 + 1/(G_2 + ...))) *)
(* For our chi(s), 1/chi has 1/s leading term (since chi ~ alpha_1 s) *)
(* Cauer-I starts with C-element (1/sC_0): *)

cauerI[zFn_, svar_, maxStages_] := Module[
  {Zcur, vals, kk, num, den, num0, denCoeff, c, r},
  Zcur = Together[zFn];
  vals = {};
  For[kk = 0, kk < maxStages, kk++,
    num = Numerator[Zcur]; den = Denominator[Zcur];
    (* Leading 1/s in Z(s)? -- check if den has factor of s *)
    den0 = den /. svar -> 0;
    num0 = num /. svar -> 0;
    If[N[Abs[den0], 30] < 10^-30,
      (* den has zero at s=0, leading 1/s present *)
      denOverS = PolynomialQuotient[den, svar, svar];
      denOverS0 = denOverS /. svar -> 0;
      c = denOverS0/num0;  (* 1/(sC_k), so coefficient = 1/c *)
      AppendTo[vals, {kk, "C_" <> ToString[kk], c}];
      Zcur = Zcur - 1/(svar c);
      Zcur = Together[Zcur],
      (* leading constant -> R_k *)
      r = num0/den0;
      AppendTo[vals, {kk, "R_" <> ToString[kk], r}];
      Zcur = Zcur - r;
      Zcur = Together[Zcur]];
    If[NumericQ[Zcur] && N[Abs[Zcur], 30] < 10^-30, Break[]];
    Zcur = Together[1/Zcur];
    ];
  vals];

(* Apply *)
ZRat = Together[-1/yFn];
stages = cauerI[ZRat, sv, 8];
Print["Cauer-I ladder elements:"];
Print[StringPadLeft["k", 3], "  ", StringPadRight["type", 8], "  value"];
Print[StringJoin[Table["-", {60}]]];
Do[Module[{kk, type, val},
  {kk, type, val} = stages[[i]];
  Print[StringPadLeft[ToString[kk], 3], "  ",
    StringPadRight[type, 8], "  ", N[val, 12]]],
  {i, Length[stages]}];

Print[""];
Print["=== Stage 5 (B) summary ==="];
Print["* From 9-mode Galerkin (Stage 5c), extracted ", NMaxAlpha,
  " Taylor moments alpha_n"];
Print["* Pade [4,4] gives 4 Foster pairs (R-L stages)"];
Print["* Pade [7,7] possible from 14 moments -- check stability"];
Print["* Cauer-I ladder: ", Length[stages], " elements extracted"];
