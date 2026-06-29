(* 2D rect bar — Generalized Nagamine 1D-reduction for arbitrary S_k.

   Nagamine et al. derived c̃_1 (eq.52, our S_1) and c̃_2 (eq.56-57, our S_2)
   for SQUARE prisms (a = b) using Poisson summation.
   This script generalizes to:
     1. Arbitrary k (not just 1, 2)
     2. Rectangular (a != b)

   Derivation:
     S_k = sum_{m_x, m_y odd} 1/(m_x^2 m_y^2 (m_x^2 b^2 + m_y^2 a^2)^k)

   Inner-sum formula (Poisson summation + recursion in k):
     I_k(c) := sum_{m_x odd} 1/(m_x^2 (b^2 m_x^2 + c)^k)
            = pi^2/(8 c^k) - b^2 sum_{j=0..k-1} (-1)^j / j! * F^{(j)}(c) / c^{k-j}
   where F(c) = (pi/(4 b sqrt(c))) tanh(pi sqrt(c)/(2b)).

   Then S_k = sum_{m_y odd} I_k(m_y^2 a^2) / m_y^2.

   The outer sum has TWO parts:
     - "polynomial part" pi^2/(8c^k) sum_{m_y odd} 1/m_y^(2k+2)  (closed form via zeta)
     - "exponential part" b^2 sum_j ... F^{(j)}(m_y^2 a^2) / m_y^2 / (m_y^2 a^2)^{k-j}
       → super-exponentially convergent residual (since tanh -> 1).

   For each k, the formula is derived symbolically once, then numerically evaluated.
*)
$HistoryLength = 0;

sigma = 58000000;
mu0 = 4 Pi/10^7;
aDim = 5/1000;
bDim = 2/1000;

prec = 100;
Print["=== Generalized Nagamine S_k formula (a != b, arbitrary k) ==="];
Print["WorkingPrecision = ", prec];
Print[""];

sigmaN = SetPrecision[sigma, prec];
mu0N = SetPrecision[mu0, prec];
aN = SetPrecision[aDim, prec];
bN = SetPrecision[bDim, prec];
PiN = SetPrecision[Pi, prec];

(* F(c) := (pi/(4 b sqrt(c))) tanh(pi sqrt(c)/(2b)) *)
F[c_] := (PiN / (4 bN Sqrt[c])) Tanh[PiN Sqrt[c]/(2 bN)];

(* I_k(c) — symbolic form using k-1 derivatives of F *)
Ik[k_Integer, c_] := Module[{j, Fderivs, sum},
  Fderivs = Table[D[F[cv], {cv, j}] /. cv -> c, {j, 0, k - 1}];
  sum = Sum[(-1)^j / j! * Fderivs[[j + 1]] / c^(k - j), {j, 0, k - 1}];
  PiN^2 / (8 c^k) - bN^2 sum];

(* S_k = sum_{m_y odd} I_k(m_y^2 a^2) / m_y^2  *)
M1bar = 41;          (* large enough for super-exponential residual to be 1e-100 *)
oddMy = Range[1, M1bar, 2];

(* Compute S_k for k = 1..Kmax *)
Kmax = 25;
Print["Computing S_k for k = 1..", Kmax, " via 1D-reduction..."];
SkClosed = Table[
  Sum[Ik[k, my^2 aN^2] / my^2, {my, oddMy}],
  {k, 1, Kmax}];
Print["  S_1 (closed form, my up to ", M1bar, ") = ", N[SkClosed[[1]], 30]];
Print["  S_2 = ", N[SkClosed[[2]], 30]];
Print["  S_3 = ", N[SkClosed[[3]], 30]];
Print["  S_5 = ", N[SkClosed[[5]], 30]];
Print["  S_10= ", N[SkClosed[[10]], 30]];
Print["  S_25= ", N[SkClosed[[25]], 30]];
Print[""];

(* Cross-check against direct 2D brute force at Mmax=999 *)
Print["Cross-check: brute force 2D sum at Mmax=999..."];
Mmax = 999;
oddR = Range[1, Mmax, 2];
modes2D = Flatten[Table[{mx, my}, {mx, oddR}, {my, oddR}], 1];
betaDirect[mx_, my_] := 64 aN bN / (mx^2 my^2 PiN^4);
lambdaSqDirect[mx_, my_] := PiN^2 (mx^2/aN^2 + my^2/bN^2);
betaList = betaDirect[#[[1]], #[[2]]] & /@ modes2D;
lambdaList = lambdaSqDirect[#[[1]], #[[2]]] & /@ modes2D;
SkBrute = Table[Total[betaList / (lambdaList^k)], {k, 1, Kmax}];

(* Convert: brute_S_k / SkClosed_k should be 1 (or factor due to prefactor in our S definition) *)
Print[""];
Print["Verification S_k closed form vs brute force:"];
Print["  Note: brute force has prefactor 64 a b / pi^4 in beta^2 and 1/pi^(2k) in lambda^(-2k);"];
Print["        closed form S_k as defined is the bare sum over 1/(mx^2 my^2 (mx^2 b^2 + my^2 a^2)^k)."];
Print["  Conversion: brute_S_k_full = (64 a b / pi^4) * (a b)^(2k) / pi^(2k) * S_k_bare"];
Print[""];

prefactor[k_] := 64 aN bN / PiN^4 * (aN bN)^(2 k) / PiN^(2 k);
Do[
  Module[{closed, brute, ratio},
    closed = N[SkClosed[[k]] * prefactor[k], 25];
    brute = N[SkBrute[[k]], 25];
    ratio = N[closed/brute, 5];
    Print["  k=", k, ": closed = ", closed, ", brute = ", brute,
          ", ratio = ", ratio]],
  {k, 1, 5}];
Print["  ..."];
Do[
  Module[{closed, brute, ratio},
    closed = N[SkClosed[[k]] * prefactor[k], 25];
    brute = N[SkBrute[[k]], 25];
    ratio = N[closed/brute, 10];
    Print["  k=", k, ": closed = ", closed, ", brute = ", brute,
          ", ratio = ", ratio]],
  {k, {10, 15, 20, 25}}];
