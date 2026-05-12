(* Task B: Pure Legendre spectral BEM for vacuum 3D problem.

   Use shifted Legendre polynomials L_p(2x/a) directly as basis (no monomial detour).
   Key benefits:
   - Mass matrix DIAGONAL: M_{(P,Q,R),(P',Q',R')} = (a/(2P+1))(b/(2Q+1))(c/(2R+1)) delta
   - No Cholesky needed, just element-wise division by sqrt(mass)
   - Better numerical conditioning at high P (vs monomial which has cond ~ 10^21 at P=8)

   Method:
   - Cell-cell K_legendre[(P,Q,R),(P',Q',R')] = int int L_P(2x/a) L_Q(2y/b) L_R(2z/c)
       L_{P'}(2x'/a) L_{Q'}(2y'/b) L_{R'}(2z'/c) / |r-r'| dV dV'
   - Compute via Legendre's coefficient list -> sum of monomial K elements
     K_leg[I,J] = sum_{i,j} T_I[i] T_J[j] K_mon[i,j] where T_I[i] = coef of x^p y^q z^r in L_P L_Q L_R

   Run at P=8, P=10 to demonstrate scaling behavior.

   Date: 2026-05-02
*)

$HistoryLength = 0;
prec = 50;

Print["================================================================"];
Print["Task B: Pure Legendre spectral BEM (vacuum 3D)"];
Print["WorkingPrecision = ", prec];
Print["================================================================"];
Print[""];

sigmaN = SetPrecision[58 10^6, prec];
mu0N = SetPrecision[4 Pi/10^7, prec];
PiN = SetPrecision[Pi, prec];
aN = SetPrecision[5/1000, prec];
bN = SetPrecision[2/1000, prec];
cN = SetPrecision[1/1000, prec];
VN = aN bN cN;
mu0sigmaN = mu0N sigmaN;

alpha1Ref = -mu0sigmaN (aN^2 + bN^2)/48;
alpha2RefNew = SetPrecision[7.4783240535434074691212813799*^-10, 28];
tau1Ref = alpha2RefNew/(-alpha1Ref);

Print["Reference values (NEW closed form, 1D-closed + 2D-NIntegrate WP=60):"];
Print["  alpha_2 = ", alpha2RefNew];
Print["  tau_1   = ", N[tau1Ref 10^6, 16], " us"];
Print[""];

(* === Inner u-integration (monomial basis) === *)
ClearAll[gAxisPos];
gAxisPos[p1_Integer, p2_Integer, h_, xi_] := gAxisPos[p1, p2, h, xi] =
  Expand[Integrate[(u + xi/2)^p1 (u - xi/2)^p2,
    {u, xi/2 - h/2, h/2 - xi/2},
    Assumptions -> {xi > 0, xi < h, h > 0}]];

ClearAll[KMon];
KMon[{p1_, q1_, r1_}, {p2_, q2_, r2_}] := KMon[{p1, q1, r1}, {p2, q2, r2}] =
  Module[{xi, eta, zeta, gxE, gyE, gzE},
    gxE = gAxisPos[p1, p2, aN, xi];
    gyE = gAxisPos[q1, q2, bN, eta];
    gzE = gAxisPos[r1, r2, cN, zeta];
    8 NIntegrate[gxE gyE gzE/Sqrt[xi^2 + eta^2 + zeta^2],
      {xi, 0, aN}, {eta, 0, bN}, {zeta, 0, cN},
      WorkingPrecision -> prec, AccuracyGoal -> prec - 5,
      PrecisionGoal -> prec - 8, MaxRecursion -> 16,
      Method -> "GlobalAdaptive"]];

(* === Legendre coefficient extraction ===
   T[P, p] = coefficient of x^p in L_P(2x/h) for x in [-h/2, h/2]. *)
ClearAll[legendreT];
legendreT[Pmax_, h_] := legendreT[Pmax, h] = Module[{T, n, k, c},
  T = ConstantArray[0, {Pmax + 1, Pmax + 1}];
  Do[
    c = CoefficientList[LegendreP[n, x], x];
    Do[T[[n + 1, k + 1]] = c[[k + 1]] (2/h)^k, {k, 0, n}],
    {n, 0, Pmax}];
  T];

(* === Build basis (e,o,e) === *)
buildBasis[Pmax_, parity_] := Module[{px, py, pz, basis},
  {px, py, pz} = parity;
  basis = {};
  Do[
    If[Mod[p, 2] == px && Mod[q, 2] == py && Mod[r, 2] == pz && p + q + r <= Pmax,
      AppendTo[basis, {p, q, r}]],
    {p, 0, Pmax}, {q, 0, Pmax - p}, {r, 0, Pmax - p - q}];
  basis];

(* === Run for given Pmax === *)
runLegendre[Pmax_] := Module[
  {basis, n, t0, KMonMat, Tx, Ty, Tz, TblockL, KLeg, massDiag,
   KLegStandard, eigvals, eigvecs, gXk, lambdaX, alphas, NMoments,
   tau1BEM, condM},
  basis = buildBasis[Pmax, {0, 1, 0}];
  n = Length[basis];
  Print["=== P_max = ", Pmax, " ==="];
  Print["  Basis (e,o,e): ", n, " functions"];

  Print["  Building K_mon (", n (n + 1)/2, " elements)..."];
  t0 = AbsoluteTime[];
  KMonMat = Table[KMon[basis[[i]], basis[[j]]], {i, n}, {j, n}];
  Print["    K_mon: ", N[AbsoluteTime[] - t0, 4], " s"];

  Tx = legendreT[Pmax, aN];
  Ty = legendreT[Pmax, bN];
  Tz = legendreT[Pmax, cN];

  TblockL = Table[
    Module[{P, Q, R, p, q, r},
      {P, Q, R} = basis[[i]];
      {p, q, r} = basis[[j]];
      Tx[[P + 1, p + 1]] Ty[[Q + 1, q + 1]] Tz[[R + 1, r + 1]]],
    {i, n}, {j, n}];

  KLeg = TblockL.KMonMat.Transpose[TblockL];

  (* Diagonal mass: M[I, I] = a/(2P+1) * b/(2Q+1) * c/(2R+1) *)
  massDiag = Table[
    Module[{P, Q, R}, {P, Q, R} = basis[[i]];
      aN/(2 P + 1) bN/(2 Q + 1) cN/(2 R + 1)],
    {i, n}];

  condM = Max[massDiag]/Min[massDiag];
  Print["  Mass condition (Legendre, diagonal): ", N[condM, 6]];

  (* Standard eigenvalue problem: K_norm[I,J] = K_leg[I,J] / sqrt(massDiag[I] * massDiag[J]) *)
  KLegStandard = Table[KLeg[[i, j]]/Sqrt[massDiag[[i]] massDiag[[j]]], {i, n}, {j, n}];

  Print["  Solving standard eigenvalue (no Cholesky)..."];
  t0 = AbsoluteTime[];
  {eigvals, eigvecs} = Eigensystem[N[KLegStandard, prec]];
  Print["    eigvalue solve: ", N[AbsoluteTime[] - t0, 4], " s"];

  (* Foster amplitudes: project A_ext^x = -B_0 y/2 *)
  (* In Legendre: y = (b/2) L_1(2y/b), so A_ext^x has coefficient -b/4 at (0, 1, 0). *)
  (* In M-orthonormal Legendre: phi_I = L_I/sqrt(mass), so coef = -b/4 * sqrt(mass[(0,1,0)]). *)
  Module[{aExt, idx010, b010},
    aExt = ConstantArray[0, n];
    idx010 = Position[basis, {0, 1, 0}, 1][[1, 1]];
    aExt[[idx010]] = -bN/4 Sqrt[massDiag[[idx010]]];
    gXk = eigvecs.aExt];

  lambdaX = eigvals/(4 PiN);
  NMoments = 8;
  alphas = Table[
    (* x-block contribution only; would double for full but err shows trend *)
    (-1)^n (mu0sigmaN)^n/VN Total[gXk^2 lambdaX^(n - 1)],
    {n, 1, NMoments}];

  Print["  alpha_n (x-block, Legendre):"];
  Do[Print["    alpha_", n, " = ", N[alphas[[n]], 16]], {n, NMoments}];

  alphas];

(* === Run P=6 first to compare with monomial result === *)
result6 = runLegendre[6];

Print[""];
Print["----------------------------------------------------------------"];
Print["Comparison: monomial spectral_v1 P=6 WP=50 (x-block) gave:"];
Print["  Note: full alpha_2 = sum (gx^2 + gy^2)*lambda needs both blocks"];
Print["  Pure Legendre at P=6 should give same x-block contribution"];
