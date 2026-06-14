(* Spectral BEM with monomial basis on cuboid.

   Basis: psi_{pqr}(r) = x^p y^q z^r, total order p+q+r <= P_max.
   Galerkin discretization of Newton operator G with kernel 1/(4 pi |r-r'|).

   Strategy:
   - Center cuboid at origin: Omega = [-a/2, a/2] x [-b/2, b/2] x [-c/2, c/2]
   - Cuboid D_2h symmetry => K, M block-diagonal in axis-parities.
   - For uniform B_z excitation:
       A_ext^x = -B_0 y/2  (parity: x even, y odd, z even)
       A_ext^y = +B_0 x/2  (parity: x odd, y even, z even)
     Two independent blocks; sum Foster contributions.
   - Mass matrix: closed form (monomial integrals).
   - Newton matrix: 6D -> 3D via center-of-mass change of variables.
       Inner u-integration: closed form polynomial in xi (separation along axis).
       Outer 3D integral over (xi, eta, zeta) of g_x g_y g_z / sqrt(xi^2 + eta^2 + zeta^2):
       evaluated with NIntegrate at WP=30.
   - Solve generalized eigvalue K phi = lambda M phi per block.
   - Compute Foster sum alpha_n, compare with closed form.

   Expected: spectral convergence in P_max:
     P=4 -> 4 digits in alpha_2
     P=6 -> 6 digits
     P=8 -> 8 digits
     P=10 -> 10+ digits (matches Nagamine 2D in spirit)

   Date: 2026-05-01
*)

$HistoryLength = 0;
prec = 50;

Print["================================================================"];
Print["BEM-Foster Phase 8: spectral basis (monomial) on cuboid"];
Print["WorkingPrecision = ", prec];
Print["================================================================"];
Print[""];

(* === Constants === *)
sigmaN = SetPrecision[58 10^6, prec];
mu0N = SetPrecision[4 Pi/10^7, prec];
PiN = SetPrecision[Pi, prec];
aN = SetPrecision[5/1000, prec];
bN = SetPrecision[2/1000, prec];
cN = SetPrecision[1/1000, prec];
VN = aN bN cN;
mu0sigmaN = mu0N sigmaN;

alpha1Ref = -mu0sigmaN (aN^2 + bN^2)/48;
alpha2Ref = SetPrecision[7.4783240535300000146605449742*^-10, prec];
tau1Ref = alpha2Ref/(-alpha1Ref);

Print["alpha_1 (closed) = ", N[alpha1Ref, 16]];
Print["alpha_2 (closed) = ", N[alpha2Ref, 16]];
Print["tau_1   (closed) = ", N[tau1Ref 10^6, 16], " us"];
Print[""];

(* === Inner u-integration (closed form polynomial), assuming xi >= 0 === *)
(* For two intervals [-h/2, h/2] each, weight (u + xi/2)^p1 (u - xi/2)^p2 over
   |x_1| < h/2, |x_2| < h/2, x_1 - x_2 = xi (xi >= 0).
   Result: g(xi; p1, p2, h), polynomial in (h, xi), no Abs[].
   The K matrix integrand is even in xi for p1+p2 even, odd for p1+p2 odd.
   We integrate xi in [0, h] and use octant symmetry. *)

ClearAll[gAxisPos];
gAxisPos[p1_Integer, p2_Integer, h_, xi_] := gAxisPos[p1, p2, h, xi] = Module[
  {u, integrand, result},
  integrand = (u + xi/2)^p1 (u - xi/2)^p2;
  result = Integrate[integrand, {u, xi/2 - h/2, h/2 - xi/2},
                     Assumptions -> {xi > 0, xi < h, h > 0}];
  Expand[result]];

(* Sanity (xi > 0 form, no Abs[]): *)
Print["Sanity: gAxisPos[0,0,h,xi] = ", gAxisPos[0, 0, h, xi]];
Print["Sanity: gAxisPos[1,1,h,xi] = ", Simplify[gAxisPos[1, 1, h, xi]]];
Print["Sanity: gAxisPos[2,0,h,xi] = ", Simplify[gAxisPos[2, 0, h, xi]]];
Print[""];

(* Octant integration parity:
   gAxisPos[p1, p2, h, xi] is for xi > 0.
   For xi < 0: integrand swaps p1 <-> p2, so g(xi<0) = gAxisPos[p2, p1, h, |xi|].
   For the K matrix integrand:
     g_x(xi) g_y(eta) g_z(zeta) / sqrt(...)
   xi -> -xi: g_x(p1,p2,xi) -> g_x(p2,p1,xi).
   So integrand at (-xi) with same p's is gx(p2, p1, xi). Can sum:
   For symmetric case (or by relabeling), integrate xi >= 0 and add p1<->p2 contribution.
   If Kelement is computed with p1<->p2 swap symmetric (we already have), can just
   integrate full range including xi=0 by extending gAxisPos. *)

(* === Newton matrix element K[(p1,q1,r1), (p2,q2,r2)] === *)
(* K = int int x^p1 y^q1 z^r1 x'^p2 y'^q2 z'^r2 / |r-r'| dV dV'
     = int dxi int deta int dzeta
        gAxis[p1,p2,a,xi] gAxis[q1,q2,b,eta] gAxis[r1,r2,c,zeta]
        / sqrt(xi^2 + eta^2 + zeta^2)
*)

ClearAll[Kelement];
Kelement[{p1_, q1_, r1_}, {p2_, q2_, r2_}] := Kelement[{p1, q1, r1}, {p2, q2, r2}] =
  Module[{gxExpr, gyExpr, gzExpr, xi, eta, zeta},
    (* Within parity block (p1+p2 even, q1+q2 even, r1+r2 even):
       gAxisPos is symmetric, integrand is even in each axis.
       Octant integration: domain (0, h_alpha) per axis, multiply by 8. *)
    gxExpr = gAxisPos[p1, p2, aN, xi];
    gyExpr = gAxisPos[q1, q2, bN, eta];
    gzExpr = gAxisPos[r1, r2, cN, zeta];
    8 NIntegrate[
      gxExpr gyExpr gzExpr/Sqrt[xi^2 + eta^2 + zeta^2],
      {xi, 0, aN}, {eta, 0, bN}, {zeta, 0, cN},
      WorkingPrecision -> prec, AccuracyGoal -> prec - 5,
      PrecisionGoal -> prec - 8, MaxRecursion -> 16,
      Method -> "GlobalAdaptive"]];

(* === Mass matrix element M[(p1,q1,r1), (p2,q2,r2)] === *)
(* M = int x^{p1+p2} y^{q1+q2} z^{r1+r2} dV
     = product over axes of int_{-h/2}^{h/2} x^n dx
       (= h^{n+1}/((n+1) 2^n) if n even, 0 if odd)
*)
massAxis[n_, h_] := If[OddQ[n], 0, h^(n + 1)/((n + 1) 2^n)];

ClearAll[Melement];
Melement[{p1_, q1_, r1_}, {p2_, q2_, r2_}] := Melement[{p1, q1, r1}, {p2, q2, r2}] =
  massAxis[p1 + p2, aN] massAxis[q1 + q2, bN] massAxis[r1 + r2, cN];

(* === Build basis for given parity and total order === *)
(* parity = {px, py, pz} where each is 0 (even) or 1 (odd) *)
(* basis index (p, q, r) with p mod 2 == px, q mod 2 == py, r mod 2 == pz, p+q+r <= Pmax *)
buildBasis[Pmax_, parity_] := Module[{px, py, pz, basis},
  {px, py, pz} = parity;
  basis = {};
  Do[
    If[Mod[p, 2] == px && Mod[q, 2] == py && Mod[r, 2] == pz && p + q + r <= Pmax,
      AppendTo[basis, {p, q, r}]],
    {p, 0, Pmax}, {q, 0, Pmax - p}, {r, 0, Pmax - p - q}];
  basis];

(* === Solve Foster eigvalue per parity block === *)
solveSpectralBlock[basis_, label_] := Module[
  {n, K, M, eigsys, eigvals, eigvecs, t0},
  n = Length[basis];
  Print["  Block ", label, ": ", n, " basis functions"];

  Print["    Building K (", n (n + 1)/2, " unique elements)..."];
  t0 = AbsoluteTime[];
  K = ConstantArray[0, {n, n}];
  Do[
    Module[{val = Kelement[basis[[i]], basis[[j]]]},
      K[[i, j]] = val;
      K[[j, i]] = val],
    {i, n}, {j, i, n}];
  Print["    K assembly: ", N[AbsoluteTime[] - t0, 4], " s"];

  M = ConstantArray[0, {n, n}];
  Do[
    M[[i, j]] = Melement[basis[[i]], basis[[j]]];
    M[[j, i]] = M[[i, j]],
    {i, n}, {j, i, n}];

  Print["    Solving generalized eigenvalue problem..."];
  t0 = AbsoluteTime[];
  (* For Hermitian K, M with M positive definite:
     transform to standard form via Cholesky.
     M = L L^T => K phi = lambda M phi  <=>
       (L^{-1} K L^{-T}) (L^T phi) = lambda (L^T phi)
  *)
  Module[{L, Linv, Khat, eigvalsRaw, eigvecsRaw},
    L = CholeskyDecomposition[M];  (* upper-triangular L^T such that M = L^T L *)
    Linv = Inverse[L];
    Khat = Transpose[Linv].K.Linv;
    {eigvalsRaw, eigvecsRaw} = Eigensystem[N[Khat, prec]];
    eigvals = eigvalsRaw;
    (* phi = L^{-1} (eigenvector of Khat) *)
    eigvecs = (Linv.#) & /@ eigvecsRaw];
  Print["    done: ", N[AbsoluteTime[] - t0, 4], " s"];

  (* Eigvecs: rows are eigenvectors in original (M-orthonormal) basis
     Verify: phi.M.phi should be ~ 1 *)
  {eigvals, eigvecs}];

(* === Foster amplitudes from spectral basis === *)
(* A_ext^x = -B_0 y/2; in basis: coefficient of psi_{0,1,0} is -B_0/2.
   So expansion vector aExt^x: 1 in slot for (0,1,0), 0 elsewhere. *)
projectAext[basis_, axis_] := Module[{vec, idx, p},
  vec = ConstantArray[0, Length[basis]];
  Switch[axis,
    "x", (* A^x = -y/2 -> coefficient -1/2 on psi_{0,1,0} *)
    idx = Position[basis, {0, 1, 0}, 1];
    If[Length[idx] > 0, vec[[idx[[1, 1]]]] = -1/2],
    "y", (* A^y = +x/2 -> coefficient +1/2 on psi_{1,0,0} *)
    idx = Position[basis, {1, 0, 0}, 1];
    If[Length[idx] > 0, vec[[idx[[1, 1]]]] = 1/2]];
  vec];

(* === Main: solve for both blocks, compute alpha_n === *)
runForPmax[Pmax_] := Module[
  {basisX, basisY, eigvalsX, eigvecsX, eigvalsY, eigvecsY,
   aExtX, aExtY, gXk, gYk, lambdaXk, lambdaYk, taus, gks2,
   alphas},
  Print["=== Spectral basis, P_max = ", Pmax, " ==="];

  basisX = buildBasis[Pmax, {0, 1, 0}];
  basisY = buildBasis[Pmax, {1, 0, 0}];

  {eigvalsX, eigvecsX} = solveSpectralBlock[basisX, "x-coupling (e,o,e)"];
  {eigvalsY, eigvecsY} = solveSpectralBlock[basisY, "y-coupling (o,e,e)"];

  aExtX = projectAext[basisX, "x"];
  aExtY = projectAext[basisY, "y"];

  (* Project A_ext on M-orthonormal eigenvectors:
     phi_k coefficients in monomial basis: eigvecs[[k]]
     g_k = <phi_k, A_ext>_M = phi_k^T M aExt = (in M-orthonormal basis) phi_k^T M aExt
     Since eigvecs are M-orthonormal: phi_k^T M phi_l = delta_kl
     Then g_k = phi_k^T M aExt *)
  Module[{Mx, My},
    Mx = Table[
      Melement[basisX[[i]], basisX[[j]]], {i, Length[basisX]}, {j, Length[basisX]}];
    My = Table[
      Melement[basisY[[i]], basisY[[j]]], {i, Length[basisY]}, {j, Length[basisY]}];
    gXk = (eigvecsX.Mx.aExtX);
    gYk = (eigvecsY.My.aExtY)];

  (* True Newton eigenvalues lambda_k = (eigval) / (4 pi)
     Note: K matrix here = int int psi_i psi_j / |r-r'| dVdV' (no 4pi factor in our K).
     The Newton operator G has kernel 1/(4pi|r-r'|), so
       K_matrix = 4 pi <G operator>
     Eigvalues of K phi = lambda M phi -> lambdaK_matrix (for K matrix problem).
     Foster eigvalue lambda = lambdaK_matrix / (4 pi).
     tau = mu0 sigma lambda. *)
  lambdaXk = eigvalsX/(4 PiN);
  lambdaYk = eigvalsY/(4 PiN);

  Print[""];
  Print["  Foster modes (top 5 each block):"];
  Module[{tauXshow, tauYshow, gXshow, gYshow, ord},
    tauXshow = mu0sigmaN lambdaXk;
    tauYshow = mu0sigmaN lambdaYk;
    ord = Reverse[Ordering[N[tauXshow, 6]]];
    Print["    x-block, top 5 (tau [us], g_k):"];
    Do[
      Print["      ", N[tauXshow[[ord[[k]]]] 10^6, 8], "  ",
            ScientificForm[N[gXk[[ord[[k]]]], 6]]],
      {k, Min[5, Length[lambdaXk]]}];
    ord = Reverse[Ordering[N[tauYshow, 6]]];
    Print["    y-block, top 5 (tau [us], g_k):"];
    Do[
      Print["      ", N[tauYshow[[ord[[k]]]] 10^6, 8], "  ",
            ScientificForm[N[gYk[[ord[[k]]]], 6]]],
      {k, Min[5, Length[lambdaYk]]}]];

  (* Foster sum: combine both blocks.
     alpha_n = (-1)^n (mu0sigma)^n / V * sum over both blocks of g_k^2 lambda^{n-1}
  *)
  alphas = Table[
    (-1)^n (mu0sigmaN)^n/VN *
    (Total[gXk^2 lambdaXk^(n - 1)] + Total[gYk^2 lambdaYk^(n - 1)]),
    {n, 1, 5}];

  Print[""];
  Print["  alpha_n (spectral) vs closed form:"];
  Print["    alpha_1 = ", N[alphas[[1]], 16], "  closed = ", N[alpha1Ref, 16],
        "  err = ", N[Abs[(alphas[[1]] - alpha1Ref)/alpha1Ref], 6]];
  Print["    alpha_2 = ", N[alphas[[2]], 16], "  closed = ", N[alpha2Ref, 16],
        "  err = ", N[Abs[(alphas[[2]] - alpha2Ref)/alpha2Ref], 6]];
  Print["    alpha_3 = ", N[alphas[[3]], 16]];
  Print["    alpha_4 = ", N[alphas[[4]], 16]];
  Print["    alpha_5 = ", N[alphas[[5]], 16]];

  (* Pade tau_1 *)
  Module[{tau1BEM},
    tau1BEM = alphas[[2]]/(-alphas[[1]]);
    Print["    tau_1   = ", N[tau1BEM 10^6, 14], " us  vs closed = ",
          N[tau1Ref 10^6, 14], " us"];
    Print["    err_tau_1 = ", N[Abs[(tau1BEM - tau1Ref)/tau1Ref], 6]]];

  alphas];

(* === Run for increasing P_max === *)
Do[
  runForPmax[Pmax];
  Print["----------------------------------------------------------------"],
  {Pmax, {6, 8}}];  (* WP=50 P=8 to verify alpha_3, alpha_4 convergence *)
