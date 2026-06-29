(* 2D rect bar 5x2 mm — Cauer-I via verified high-precision moments.

   Strategy (matching 3D cuboid_521_cauerI_moments.wls pattern):
     M_0 = sigma * a * b           (Parseval, EXACT closed form)
     S_1 : Nagamine eq.52-style 1D reduction (Poisson + tanh identity)
            generalized for rectangular cross-section a != b
            -> exponentially convergent residual sum,
               ~50 digit accuracy in <30 terms
     S_k (k>=2) : direct double sum at Mmax=999 with WorkingPrecision=100
                  tail ~ 1/M^(2k+1) is fast: k=2 gives 1e-15 relative
     Cauer-I extraction at WorkingPrecision=100

   Eigenfunctions: phi_(mx,my) = (2/sqrt(ab)) sin(mx pi x/a) sin(my pi y/b)
   Eigenvalues:    lambda^2 = pi^2 (mx^2/a^2 + my^2/b^2)
   beta^2:         64 a b / (mx^2 my^2 pi^4)  (odd mx, my)

   Per-unit-length analysis (no z-dimension), so M_0 = sigma*a*b (no c factor).

   Reference: Nagamine et al. eq. 50-55, generalized to a != b.
   Identity used (Poisson summation):
     sum_(mx odd) 1/(mx^2 b^2 + Y^2) = (Pi/(4 b Y)) tanh(Pi Y/(2b))
*)
$HistoryLength = 0;

sigma = 58000000;
mu0 = 4 Pi/10^7;
aDim = 5/1000;
bDim = 2/1000;

prec = 100;
Print["=== 2D rect 5x2 mm — Cauer-I via verified moments (Nagamine, prec=", prec, ") ==="];
Print[""];

sigmaN = SetPrecision[sigma, prec];
mu0N = SetPrecision[mu0, prec];
aN = SetPrecision[aDim, prec];
bN = SetPrecision[bDim, prec];
PiN = SetPrecision[Pi, prec];

beta2N[mx_, my_] := 64 aN bN / (mx^2 my^2 PiN^4);
lambda2N[mx_, my_] := PiN^2 (mx^2/aN^2 + my^2/bN^2);

(* M_0 closed form: σab *)
M0 = sigmaN aN bN;
Print["M_0 = σ a b (EXACT) = ", N[M0, 30]];
Print["R_0 = 1/M_0 = ", N[1/M0, 30], " Ω/m"];
Print[""];

(* === S_1 via Nagamine eq.52 generalized to 2D rectangle a!=b ===
   S_1_2D = a b^3 / 12
          - (16 b^4)/pi^5 * [ (31/32) zeta(5)
                              - 2 sum_(my odd) 1/(my^5 (1+exp(pi my a/b))) ]

   Note: 2D has no c factor (compare 3D cuboid where every term is multiplied by c).
*)
zeta5 = SetPrecision[Zeta[5], prec];
expRatio = aN/bN;          (* = 5/2 = 2.5 *)
M1bar = 31;
S1Tail = Sum[1/(my^5 (1 + Exp[PiN my expRatio])), {my, 1, M1bar, 2}];
S1 = aN bN^3/12 - (16 bN^4)/PiN^5 *
       ((31/32) zeta5 - 2 S1Tail);

Print["S_1 (Nagamine-style closed form, 2D):"];
Print["  zeta(5) = ", N[zeta5, 30]];
Print["  exp residual sum (my=1,3,..,", M1bar, ") = ", N[S1Tail, 30]];
Print["  S_1 = ", N[S1, 30]];
S1TailBound = N[1/(2 (M1bar + 2)^4 Exp[PiN (M1bar + 2) expRatio]), 5];
Print["  truncation bound for residual sum: < ", S1TailBound];
Print[""];

(* === S_k (k>=2) via direct double sum at Mmax=999, prec=100 === *)
Mmax = 999;
oddR = Range[1, Mmax, 2];
Kmax = 25;

Print["S_k (k>=2) by direct sum, Mmax=", Mmax, " (per direction)..."];
modes2D = Flatten[Table[{mx, my}, {mx, oddR}, {my, oddR}], 1];
betaList = beta2N[#[[1]], #[[2]]] & /@ modes2D;
lambdaList = lambda2N[#[[1]], #[[2]]] & /@ modes2D;
Print["  ", Length[modes2D], " modes evaluated."];

SkValues = ConstantArray[0, Kmax];
SkValues[[1]] = S1;
Do[
  SkValues[[k]] = Total[betaList / (lambdaList^k)];
  , {k, 2, Kmax}];

Print["  S_1 (Nagamine closed)  = ", N[SkValues[[1]], 30]];
Print["  S_2 (direct sum)       = ", N[SkValues[[2]], 30]];
Print["  S_5                    = ", N[SkValues[[5]], 30]];
Print["  S_10                   = ", N[SkValues[[10]], 30]];

(* Cross-check S_1 via brute force *)
S1bf = Total[betaList / lambdaList];
Print[""];
Print["Cross-check S_1:"];
Print["  Nagamine closed = ", N[SkValues[[1]], 30]];
Print["  Brute force     = ", N[S1bf, 30]];
Print["  Relative diff   = ", N[Abs[(SkValues[[1]] - S1bf)/SkValues[[1]]], 5]];
Print[""];

(* Moments *)
MkValues = Prepend[Table[sigmaN^(k+1) mu0N^k SkValues[[k]], {k, 1, Kmax}], M0];
Print["Moments M_k (k=0..5):"];
Do[Print["  M_", k, " = ", N[MkValues[[k+1]], 30]], {k, 0, 5}];
Print[""];

R0pred = 1/M0;
L0pred = R0pred^2 * MkValues[[2]];
Print["L_0 = R_0^2 * M_1 = ", N[L0pred * 10^9, 30], " nH/m"];
Print[""];

(* Cauer-I extraction *)
Yseries = Sum[(-sv)^k MkValues[[k+1]], {k, 0, Kmax}];
Yseries = Expand[Yseries];
Zrat = Together[1/Yseries];
zNum = Numerator[Zrat]; zDen = Denominator[Zrat];
Print["Z(s) deg num = ", Exponent[zNum, sv], ", deg den = ", Exponent[zDen, sv]];
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
Print["Cauer-I result (", Length[stages], " entries):"];
Print["Index | Type  | Value (16 digits)"];
Print[StringJoin[Table["-", {72}]]];
Do[Module[{idx, type, val, displayVal, numVal},
  {idx, type, val} = stages[[i]];
  If[type == "R",
    numVal = N[val, 16];
    displayVal = ToString[NumberForm[numVal, 16]] <> " Ohm/m"];
  If[type == "L",
    numVal = N[val * 10^9, 16];
    displayVal = ToString[NumberForm[numVal, 16]] <> " nH/m"];
  Print[StringJoin[{StringPadLeft[ToString[idx], 5], " | ",
    StringPadRight[type, 5], " | ", displayVal}]]],
  {i, Length[stages]}];
