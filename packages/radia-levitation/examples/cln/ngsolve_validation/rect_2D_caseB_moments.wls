(* 2D rect bar 5x2 mm — Cauer-I for Case B (uniform H_0 applied).

   DUALITY WITH CASE A:
   ====================
   Both Case A (uniform E_z impressed) and Case B (uniform H_0 applied)
   reduce to the same scalar Dirichlet diffusion equation:

     -nabla^2 phi + s mu sigma phi = 0   inside conductor
     phi |_boundary = const               (= -E_z/s for A, = H_0 for B)

   The "skin function" g(s) = (1/ab) integral phi/const dxdy is shared:

     g(s) = (64/pi^4) sum_{odd mx,my} 1/(mx^2 my^2) * 1/(1 + s tau)

   Port admittances differ:
     Y_A(s) = sigma a b * g(s)             (current per unit length per E)
     Y_B(s) = (1 - g(s)) / s                (magnetic admittance)

   The MOMENT RELATION (Y_X(s) = sum (-s)^k M_k^X) is a simple SHIFT:

     M_k^B = M_{k+1}^A / M_0^A

   So we can re-use Case A's S_k computation directly. The Nagamine
   eq.52-style closed form for S_1 (used in Case A) gives Case B's M_0^B
   exactly, S_2 brute force gives M_1^B, etc.

   References:
     Lammeraner-Stafl chapter on plate eddy currents
     Kitao et al. RM-25-045 (IPMSM modified-conductivity, eq.1)
     Nagamine et al. 2026 for the verified-moment / Cauer extraction.
*)
$HistoryLength = 0;

sigma = 58000000;
mu0 = 4 Pi/10^7;
aDim = 5/1000;
bDim = 2/1000;

prec = 100;
Print["=== 2D rect 5x2 mm — Cauer-I Case B (uniform H_0 applied) ==="];
Print["Strategy: derive moments by SHIFT from Case A (M_k^B = M_{k+1}^A / M_0^A)"];
Print[""];

sigmaN = SetPrecision[sigma, prec];
mu0N = SetPrecision[mu0, prec];
aN = SetPrecision[aDim, prec];
bN = SetPrecision[bDim, prec];
PiN = SetPrecision[Pi, prec];

beta2N[mx_, my_] := 64 aN bN / (mx^2 my^2 PiN^4);
lambda2N[mx_, my_] := PiN^2 (mx^2/aN^2 + my^2/bN^2);

(* M_0^A = sigma a b (Parseval) — we will need this for the shift normalization *)
M0A = sigmaN aN bN;
Print["M_0^A (Case A) = sigma * a * b = ", N[M0A, 30]];

(* Compute S_k^A for k = 1..Kmax+1 (one extra needed for M_k^B = M_{k+1}^A / M_0^A) *)
KmaxB = 25;
KmaxAneeded = KmaxB + 1;

(* S_1: Nagamine eq.52 closed form, generalized to a != b *)
zeta5 = SetPrecision[Zeta[5], prec];
expRatio = aN/bN;
M1bar = 31;
S1Tail = Sum[1/(my^5 (1 + Exp[PiN my expRatio])), {my, 1, M1bar, 2}];
S1 = aN bN^3/12 - (16 bN^4)/PiN^5 *
       ((31/32) zeta5 - 2 S1Tail);
Print["S_1^A (Nagamine closed) = ", N[S1, 30]];

(* S_k^A for k >= 2: brute force at Mmax=999, prec=100 *)
Mmax = 999;
oddR = Range[1, Mmax, 2];

Print["S_k^A (k=2..", KmaxAneeded, ") via direct sum, Mmax=", Mmax, " ..."];
modes2D = Flatten[Table[{mx, my}, {mx, oddR}, {my, oddR}], 1];
betaList = beta2N[#[[1]], #[[2]]] & /@ modes2D;
lambdaList = lambda2N[#[[1]], #[[2]]] & /@ modes2D;

SkAll = ConstantArray[0, KmaxAneeded];
SkAll[[1]] = S1;
Do[
  SkAll[[k]] = Total[betaList / (lambdaList^k)];
  , {k, 2, KmaxAneeded}];
Print["  S_2 = ", N[SkAll[[2]], 30]];
Print["  S_5 = ", N[SkAll[[5]], 30]];
Print["  S_10 = ", N[SkAll[[10]], 30]];
Print["  S_26 = ", N[SkAll[[26]], 30]];

(* Case A moments: M_k^A = sigma^(k+1) mu^k S_k^A for k>=1; M_0^A = sigma a b *)
MkA = Prepend[Table[sigmaN^(k+1) mu0N^k SkAll[[k]], {k, 1, KmaxAneeded}], M0A];
Print[""];
Print["Case A moments (verification, M_k^A for k=0..3):"];
Do[Print["  M_", k, "^A = ", N[MkA[[k+1]], 30]], {k, 0, 3}];
Print[""];

(* Case B moments via shift: M_k^B = M_{k+1}^A / M_0^A *)
MkB = Table[MkA[[k+2]] / MkA[[1]], {k, 0, KmaxB-1}];
(* MkB has KmaxB entries: M_0^B, M_1^B, ..., M_{KmaxB-1}^B *)
Print["Case B moments via shift M_k^B = M_{k+1}^A / M_0^A:"];
Do[Print["  M_", k, "^B = ", N[MkB[[k+1]], 30]], {k, 0, 5}];
Print[""];

(* R_0^B and L_0^B prediction from first two moments *)
R0B = 1/MkB[[1]];
L0Bpred = R0B^2 * MkB[[2]];
Print["R_0^B = 1/M_0^B = ", N[R0B, 30]];
Print["  units of R_0^B: depend on port normalization."];
Print["  Y_B(0) = (1-g(0))/s |_{s=0} = -g'(0) = sum tau_n * 64/(pi^4 mx^2 my^2)"];
Print["L_0^B = R_0^B^2 * M_1^B = ", N[L0Bpred, 30]];
Print[""];

(* Build polynomial Y_B(s) = sum (-s)^k M_k^B and extract Cauer-I *)
Yseries = Sum[(-sv)^k MkB[[k+1]], {k, 0, KmaxB-1}];
Yseries = Expand[Yseries];
Zrat = Together[1/Yseries];
zNum = Numerator[Zrat]; zDen = Denominator[Zrat];
Print["Z_B(s) = 1/Y_B(s) deg num = ", Exponent[zNum, sv], ", deg den = ", Exponent[zDen, sv]];
Print[""];

cauerExtract[zRat_, svar_, maxStages_] := Module[
  {Zcur, vals, kk, pp, qq, q0, p0, Rk, presid, pas, pas0, Lk, nresid,
   nresidS, pNew, qNew},
  Zcur = zRat;
  vals = {};
  For[kk = 0, kk < maxStages, kk++,
   pp = Numerator[Zcur];
   qq = Denominator[Zcur];
   q0 = qq /. svar -> 0;
   p0 = pp /. svar -> 0;
   Print["[k=", kk, "] deg(p)=", Exponent[pp, svar], " deg(q)=", Exponent[qq, svar]];
   If[N[Abs[q0], 30] < 10^-200, Print["  q(0)~0, stop"]; Break[]];
   Rk = Together[p0/q0];
   AppendTo[vals, {2 kk, "R", Rk}];
   presid = Expand[pp - Rk qq];
   pas = PolynomialQuotient[presid, svar, svar];
   pas0 = pas /. svar -> 0;
   If[N[Abs[pas0], 30] < 10^-200, Print["  pas(0)~0, stop"]; Break[]];
   Lk = Together[pas0/q0];
   AppendTo[vals, {2 kk + 1, "L", Lk}];
   nresid = Expand[qq Lk - pas];
   nresidS = PolynomialQuotient[nresid, svar, svar];
   pNew = Expand[Lk pas];
   qNew = nresidS;
   Zcur = Together[pNew/qNew];
   If[PossibleZeroQ[Numerator[Zcur]], Break[]];
   ];
  vals];

stages = cauerExtract[Zrat, sv, 12];
Print[""];
Print["Cauer-I Case B (", Length[stages], " entries):"];
Print["Note: R_n^B, L_n^B are in 'magnetic-admittance port' units (different from Case A)"];
Print["Index | Type  | Value (16 digits)"];
Print[StringJoin[Table["-", {72}]]];
Do[Module[{idx, type, val, displayVal, numVal},
  {idx, type, val} = stages[[i]];
  If[type == "R",
    numVal = N[val, 16];
    displayVal = ToString[NumberForm[numVal, 16]]];
  If[type == "L",
    numVal = N[val, 16];
    displayVal = ToString[NumberForm[numVal, 16]]];
  Print[StringJoin[{StringPadLeft[ToString[idx], 5], " | ",
    StringPadRight[type, 5], " | ", displayVal}]]],
  {i, Length[stages]}];
