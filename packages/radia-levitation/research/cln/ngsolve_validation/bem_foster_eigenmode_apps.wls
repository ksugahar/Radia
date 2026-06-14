(* P3 + P4: Eigenmode applications.

   Use spectral basis (monomial, P=4, WP=30) to obtain Foster eigenmodes,
   then demonstrate:
     P3: Evaluate eigenmode along corner-to-corner diagonal,
         check for corner singularity (decay rate of monomial coefficients)
     P4: Lorentz force / magnetic moment from analytic integration

   Date: 2026-05-01
*)

$HistoryLength = 0;
prec = 30;

Print["================================================================"];
Print["P3 + P4: eigenmode polynomial applications"];
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

(* === Reused: gAxisPos, Kelement, mass element === *)
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

(* === Build basis (e,o,e) — x-coupled === *)
Pmax = 6;
basis = {};
Do[
  If[Mod[p, 2] == 0 && Mod[q, 2] == 1 && Mod[r, 2] == 0 && p + q + r <= Pmax,
    AppendTo[basis, {p, q, r}]],
  {p, 0, Pmax}, {q, 0, Pmax - p}, {r, 0, Pmax - p - q}];

n = Length[basis];
Print["Basis (e,o,e) at P=", Pmax, ": ", n, " functions"];
Print["  ", basis];
Print[""];

Print["Building K, M matrices..."];
t0 = AbsoluteTime[];
K = Table[Kelement[basis[[i]], basis[[j]]], {i, n}, {j, n}];
M = Table[Melement[basis[[i]], basis[[j]]], {i, n}, {j, n}];
Print["  done: ", N[AbsoluteTime[] - t0, 3], " s"];

L = CholeskyDecomposition[M];
Linv = Inverse[L];
Khat = Transpose[Linv].K.Linv;
{eigvals, eigvecsRaw} = Eigensystem[N[Khat, prec]];
eigvecs = (Linv.#) & /@ eigvecsRaw;

lambdaK = eigvals/(4 PiN);
taus = mu0sigmaN lambdaK;

ord = Reverse[Ordering[N[taus, 6]]];
Print[""];
Print["Foster time constants (sorted descending):"];
Do[
  Print["  tau_", k, " = ", N[taus[[ord[[k]]]] 10^6, 8], " us"],
  {k, n}];
Print[""];

(* === Helper: build polynomial in (x, y, z) from eigvec coefficients === *)
ClearAll[xs, ys, zs];
modePoly[idx_] := Sum[
  eigvecs[[idx, i]] xs^basis[[i, 1]] ys^basis[[i, 2]] zs^basis[[i, 3]],
  {i, n}];

(* === P3: Corner singularity check === *)
Print["================================================================"];
Print["P3: Corner singularity check"];
Print["================================================================"];
Print[""];

domIdx = ord[[1]];
phi1 = modePoly[domIdx];
Print["Dominant mode tau_1 = ", N[taus[[domIdx]] 10^6, 10], " us"];
Print["Polynomial form (showing |coef| > 1% max):"];
absC = Abs[eigvecs[[domIdx]]];
maxC = Max[absC];
Do[
  If[absC[[i]] > 0.01 maxC,
    Print["    ", N[eigvecs[[domIdx, i]], 6], " * x^", basis[[i,1]],
          " y^", basis[[i,2]], " z^", basis[[i,3]]]],
  {i, n}];
Print[""];

Print["3a. Coefficient magnitudes by total monomial order:"];
Print["   (rapid decay = analytic eigenfunction; slow decay = singular)"];
Print[""];
maxOrderInBasis = Pmax;
Do[
  Module[{maxAtOrder = 0, idx},
    Do[
      If[Total[basis[[i]]] == ordTarget && absC[[i]] > maxAtOrder,
        maxAtOrder = absC[[i]]],
      {i, n}];
    Print["   order ", ordTarget, ": max |coef| = ", N[maxAtOrder, 6]]],
  {ordTarget, 1, Pmax, 2}];
Print[""];

Print["3b. Eigenmode value along corner-to-corner diagonal (x = t a/2, y = t b/2, z = t c/2):"];
Print["   t  -> +1 hits corner (a/2, b/2, c/2)"];
Print[""];
Do[
  Module[{tNum = SetPrecision[t, prec], val},
    val = phi1 /. {xs -> tNum aN/2, ys -> tNum bN/2, zs -> tNum cN/2};
    Print["   t = ", N[t, 4],
          "  (x,y,z)=(", N[t aN 1000/2, 4], ",", N[t bN 1000/2, 4], ",", N[t cN 1000/2, 4], ") mm:",
          "  phi_1 = ", N[val, 8]]],
  {t, {-1, -0.99, -0.9, -0.5, 0, 0.5, 0.9, 0.99, 1}}];
Print[""];

Print["3c. Eigenmode value along x-axis (y=0, z=0): should be 0 by odd-y symmetry"];
Do[
  Module[{tNum = SetPrecision[t, prec], val},
    val = phi1 /. {xs -> tNum aN/2, ys -> 0, zs -> 0};
    Print["   x = ", N[t aN 1000/2, 4], " mm: phi_1 = ", N[val, 4]]],
  {t, {0, 0.5, 1.0}}];
Print[""];

Print["3d. Eigenmode value along y-axis (x=0, z=0): linear in y dominates"];
Do[
  Module[{tNum = SetPrecision[t, prec], val},
    val = phi1 /. {xs -> 0, ys -> tNum bN/2, zs -> 0};
    Print["   y = ", N[t bN 1000/2, 4], " mm: phi_1 = ", N[val, 8]]],
  {t, {-1, -0.5, 0, 0.5, 1}}];
Print[""];

(* === P4: Lorentz force and magnetic moment via analytic integration === *)
Print["================================================================"];
Print["P4: Lorentz force / magnetic moment from analytic integration"];
Print["================================================================"];
Print[""];

Print["Setup:"];
Print["  External excitation: B_ext = B_0 z_hat (uniform, transient)"];
Print["  Vector potential: A_ext = (-B_0 y/2, +B_0 x/2, 0)"];
Print["  At frequency omega, induced J via dominant mode: J_x = -i omega sigma phi_1(r)"];
Print[""];

Print["P4a. Volume integral of phi_1 (= 0 by parity):"];
volPhi1 = Integrate[phi1,
  {xs, -aN/2, aN/2}, {ys, -bN/2, bN/2}, {zs, -cN/2, cN/2}];
Print["  int phi_1 dV = ", N[volPhi1, 8]];
Print[""];

Print["P4b. Magnetic moment integral m_z ~ -int y phi_1 dV (non-zero, key physical quantity):"];
yPhi1 = Integrate[ys phi1,
  {xs, -aN/2, aN/2}, {ys, -bN/2, bN/2}, {zs, -cN/2, cN/2}];
Print["  int y phi_1 dV = ", N[yPhi1, 14]];
Print[""];

Print["P4c. Energy integral |grad phi_1|^2 (closed-form polynomial integral):"];
gradPhi1 = {D[phi1, xs], D[phi1, ys], D[phi1, zs]};
EgradInt = Integrate[Total[# ^2 & /@ gradPhi1],
  {xs, -aN/2, aN/2}, {ys, -bN/2, bN/2}, {zs, -cN/2, cN/2}];
Print["  int |grad phi_1|^2 dV = ", N[EgradInt, 14]];
Print[""];

Print["P4d. Laplacian (induced current source in Coulomb gauge):"];
laplPhi1 = D[phi1, {xs, 2}] + D[phi1, {ys, 2}] + D[phi1, {zs, 2}];
Print["  Laplacian phi_1 = ", N[Simplify[laplPhi1], 6]];
Print[""];

Print["P4e. Surface flux through y=b/2 face (eddy current loop signature):"];
fluxY = Integrate[D[phi1, ys] /. ys -> bN/2,
  {xs, -aN/2, aN/2}, {zs, -cN/2, cN/2}];
Print["  int (d_y phi_1)|_{y=b/2} dx dz = ", N[fluxY, 12]];
Print[""];

(* === Summary === *)
Print["================================================================"];
Print["Summary"];
Print["================================================================"];
Print[""];
Print["* All quantities above computed by SYMBOLIC INTEGRATION OF POLYNOMIALS."];
Print["* No numerical integration in eigenmode application stage."];
Print["* Foster eigenmodes are explicit polynomials in (x, y, z)."];
Print["* Lorentz force, magnetic moment, surface flux, energy:"];
Print["  all closed form once eigenmode coefficients are known."];
Print[""];
Print["This is the key advantage of spectral basis vs piecewise const BEM."];
Print["For Radia integration / accelerator magnet applications:"];
Print["  - Eigenmode databases pre-computed (one-time cost)"];
Print["  - Run-time queries (Lorentz force at given frequency, etc.)"];
Print["    require only polynomial evaluation/integration"];
