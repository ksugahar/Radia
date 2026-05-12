(* Demo: extract Foster eigenmodes as analytic polynomials from spectral BEM.

   Run a small P=4 spectral BEM (fast, ~30s at WP=30), then:
   - Print the dominant Foster mode as polynomial in (x, y, z)
   - Compute curl-curl (= -laplacian for Coulomb gauge) as polynomial
   - Demonstrate analytic Lorentz force integration
   - Save mode coefficients

   Date: 2026-05-01
*)

$HistoryLength = 0;
prec = 30;

Print["================================================================"];
Print["Foster eigenmode demo: analytic polynomial form"];
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

(* Reuse octant gAxisPos, Kelement from spectral_v1 *)
ClearAll[gAxisPos];
gAxisPos[p1_Integer, p2_Integer, h_, xi_] := gAxisPos[p1, p2, h, xi] =
  Expand[Integrate[(u + xi/2)^p1 (u - xi/2)^p2,
    {u, xi/2 - h/2, h/2 - xi/2},
    Assumptions -> {xi > 0, xi < h, h > 0}]];

ClearAll[Kelement];
Kelement[{p1_, q1_, r1_}, {p2_, q2_, r2_}] := Kelement[{p1, q1, r1}, {p2, q2, r2}] =
  Module[{xi, eta, zeta, gxE, gyE, gzE},
    gxE = gAxisPos[p1, p2, aN, xi];
    gyE = gAxisPos[q1, q2, bN, eta];
    gzE = gAxisPos[r1, r2, cN, zeta];
    8 NIntegrate[gxE gyE gzE/Sqrt[xi^2 + eta^2 + zeta^2],
      {xi, 0, aN}, {eta, 0, bN}, {zeta, 0, cN},
      WorkingPrecision -> prec, AccuracyGoal -> prec - 5,
      PrecisionGoal -> prec - 8, MaxRecursion -> 16,
      Method -> "GlobalAdaptive"]];

massAxis[n_, h_] := If[OddQ[n], 0, h^(n + 1)/((n + 1) 2^n)];
ClearAll[Melement];
Melement[{p1_, q1_, r1_}, {p2_, q2_, r2_}] := Melement[{p1, q1, r1}, {p2, q2, r2}] =
  massAxis[p1 + p2, aN] massAxis[q1 + q2, bN] massAxis[r1 + r2, cN];

(* === Build basis (e,o,e) at P=4 — only x-coupling block for demo === *)
Pmax = 4;
basis = {};
Do[
  If[Mod[p, 2] == 0 && Mod[q, 2] == 1 && Mod[r, 2] == 0 && p + q + r <= Pmax,
    AppendTo[basis, {p, q, r}]],
  {p, 0, Pmax}, {q, 0, Pmax - p}, {r, 0, Pmax - p - q}];

n = Length[basis];
Print["Basis (e,o,e) at P=", Pmax, ": ", n, " functions"];
Print["  ", basis];
Print[""];

(* Build K, M *)
Print["Building K, M..."];
t0 = AbsoluteTime[];
K = Table[Kelement[basis[[i]], basis[[j]]], {i, n}, {j, n}];
M = Table[Melement[basis[[i]], basis[[j]]], {i, n}, {j, n}];
Print["  done: ", N[AbsoluteTime[] - t0, 3], " s"];

(* Solve generalized eigenvalue *)
L = CholeskyDecomposition[M];
Linv = Inverse[L];
Khat = Transpose[Linv].K.Linv;
{eigvals, eigvecsRaw} = Eigensystem[N[Khat, prec]];
eigvecs = (Linv.#) & /@ eigvecsRaw;  (* M-orthonormal basis-coefficient vectors *)

(* Foster time constants *)
lambdaK = eigvals/(4 PiN);
taus = mu0sigmaN lambdaK;

Print[""];
Print["Foster time constants (sorted descending):"];
ord = Reverse[Ordering[N[taus, 6]]];
Do[
  Print["  tau_", k, " = ", N[taus[[ord[[k]]]] 10^6, 8], " us"],
  {k, n}];
Print[""];

(* === Construct polynomial form of leading Foster mode === *)
modeIndex = ord[[1]];  (* largest tau *)
phiCoeffs = eigvecs[[modeIndex]];

Print["Leading Foster mode (tau = ", N[taus[[modeIndex]] 10^6, 8], " us):"];
Print["  As polynomial in (x, y, z):"];

(* Build symbolic expression *)
ClearAll[x, y, z];
phiPoly = Sum[
  phiCoeffs[[i]] x^basis[[i, 1]] y^basis[[i, 2]] z^basis[[i, 3]],
  {i, n}];
Print["  phi_1(x, y, z) = ", N[phiPoly, 6]];
Print[""];

(* Simpler view: dominant terms *)
absCoeffs = Abs[phiCoeffs];
maxC = Max[absCoeffs];
Print["  Dominant terms (|coeff| > 0.01 max):"];
Do[
  If[absCoeffs[[i]] > 0.01 maxC,
    Print["    ", N[phiCoeffs[[i]], 6], "  *  ",
          x^basis[[i, 1]] y^basis[[i, 2]] z^basis[[i, 3]] // TraditionalForm]],
  {i, n}];
Print[""];

(* === Compute Laplacian (related to current density for Coulomb gauge) === *)
laplPhi = D[phiPoly, {x, 2}] + D[phiPoly, {y, 2}] + D[phiPoly, {z, 2}];
Print["Laplacian of leading mode (related to induced current source):"];
Print["  Laplacian phi_1 = ", N[laplPhi, 6]];
Print[""];

(* === Compute mode "gradient norm" integrated over cuboid: ===
   E_grad = int |grad phi|^2 dV = closed form polynomial integral *)
gradPhi = {D[phiPoly, x], D[phiPoly, y], D[phiPoly, z]};
gradNorm2 = Total[# ^2 & /@ gradPhi];

(* Integrate analytically over cuboid *)
EGrad = Integrate[gradNorm2,
  {x, -aN/2, aN/2}, {y, -bN/2, bN/2}, {z, -cN/2, cN/2}];
Print["Analytic integration of |grad phi_1|^2 over cuboid:"];
Print["  E_grad = ", N[EGrad, 12]];
Print[""];

(* === Save coefficients to JSON === *)
Export["bem_foster_eigenmode_demo.json",
  <|"basis" -> basis,
    "phiCoeffs" -> N[phiCoeffs, 16],
    "tau1" -> N[taus[[modeIndex]] 10^6, 16],
    "P" -> Pmax|>];
Print["Mode coefficients saved to bem_foster_eigenmode_demo.json"];

(* === Bonus: Lorentz-like integral (demonstration) === *)
(* For a "current" J_x = -laplacian phi (Coulomb gauge proxy), and B_ext = B_0 z_hat,
   Lorentz force has y-component proportional to int J_x dV (very rough).
   Compute this as analytic integral. *)
(* Actually for J = -sigma dot A, full computation requires multi-mode coupling.
   For demo, just show int phi over cuboid (volume average): *)
volAvgPhi = Integrate[phiPoly,
  {x, -aN/2, aN/2}, {y, -bN/2, bN/2}, {z, -cN/2, cN/2}]/VN;
Print[""];
Print["Volume average of phi_1 over cuboid (demo of analytic integration):"];
Print["  <phi_1>_Omega = ", N[volAvgPhi, 12]];
