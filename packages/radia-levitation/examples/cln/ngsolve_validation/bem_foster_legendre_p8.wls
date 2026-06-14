(* Issue 5 fix: re-test Legendre vs monomial at higher P to verify benefit.

   At P=8 the monomial basis has condition number issues from scale variation
   (x^p magnitudes range over many decades for cuboid 5x2x1 mm).
   Legendre polynomials should mitigate this and improve numerical stability.

   We compare:
   - Monomial (v1 spectral approach)
   - Legendre (v3 approach)
   at P=8 WP=30 (faster than WP=50 to make P=8 feasible).

   Date: 2026-05-01
*)

$HistoryLength = 0;
prec = 30;

Print["================================================================"];
Print["Issue 5: Legendre vs monomial at P=8 WP=30"];
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
alpha2Ref = SetPrecision[7.4783240521932061573535486933*^-10, 28];
tau1Ref = alpha2Ref/(-alpha1Ref);

(* === Inner u-integration === *)
ClearAll[gAxisPos];
gAxisPos[p1_Integer, p2_Integer, h_, xi_] := gAxisPos[p1, p2, h, xi] =
  Expand[Integrate[(u + xi/2)^p1 (u - xi/2)^p2,
    {u, xi/2 - h/2, h/2 - xi/2},
    Assumptions -> {xi > 0, xi < h, h > 0}]];

ClearAll[KelementMon];
KelementMon[{p1_, q1_, r1_}, {p2_, q2_, r2_}] := KelementMon[{p1, q1, r1}, {p2, q2, r2}] =
  Module[{xi, eta, zeta, gxE, gyE, gzE},
    gxE = gAxisPos[p1, p2, aN, xi];
    gyE = gAxisPos[q1, q2, bN, eta];
    gzE = gAxisPos[r1, r2, cN, zeta];
    8 NIntegrate[gxE gyE gzE/Sqrt[xi^2 + eta^2 + zeta^2],
      {xi, 0, aN}, {eta, 0, bN}, {zeta, 0, cN},
      WorkingPrecision -> prec, AccuracyGoal -> prec - 5,
      PrecisionGoal -> prec - 8, MaxRecursion -> 14,
      Method -> "GlobalAdaptive"]];

massAxis[n_, h_] := If[OddQ[n], 0, h^(n + 1)/((n + 1) 2^n)];
ClearAll[Melement];
Melement[{p1_, q1_, r1_}, {p2_, q2_, r2_}] := Melement[{p1, q1, r1}, {p2, q2, r2}] =
  massAxis[p1 + p2, aN] massAxis[q1 + q2, bN] massAxis[r1 + r2, cN];

(* === Build basis === *)
Pmax = 8;
basisX = {};
Do[
  If[Mod[p, 2] == 0 && Mod[q, 2] == 1 && Mod[r, 2] == 0 && p + q + r <= Pmax,
    AppendTo[basisX, {p, q, r}]],
  {p, 0, Pmax}, {q, 0, Pmax - p}, {r, 0, Pmax - p - q}];
nX = Length[basisX];
Print["Basis (e,o,e) at P=", Pmax, ": ", nX, " functions"];
Print[""];

(* === Build K (monomial) === *)
Print["Building K_mon (one-time)..."];
t0 = AbsoluteTime[];
KMon = Table[KelementMon[basisX[[i]], basisX[[j]]], {i, nX}, {j, nX}];
Print["  K_mon assembly: ", N[AbsoluteTime[] - t0, 4], " s"];

MMon = Table[Melement[basisX[[i]], basisX[[j]]], {i, nX}, {j, nX}];

(* === Method 1: Monomial + Cholesky === *)
Print[""];
Print["=== Method 1: Monomial basis + Cholesky ==="];
t0 = AbsoluteTime[];
condMon = N[Sqrt[Max[Eigenvalues[N[MMon, prec]]]/Min[Eigenvalues[N[MMon, prec]]]], 6];
Print["  M_mon condition number: ", condMon];
LMon = CholeskyDecomposition[MMon];
LMonInv = Inverse[LMon];
KMonHat = Transpose[LMonInv].KMon.LMonInv;
{eigvalsMon, eigvecsMonRaw} = Eigensystem[N[KMonHat, prec]];
eigvecsMon = (LMonInv.#) & /@ eigvecsMonRaw;
Print["  eigvalue solve: ", N[AbsoluteTime[] - t0, 4], " s"];

(* === Method 2: Legendre (mass diagonal, no Cholesky) === *)
Print[""];
Print["=== Method 2: Legendre basis (mass diagonal) ==="];
t0 = AbsoluteTime[];

(* Legendre transformation: T[P,p] = coefficient of x^p in P_n(2x/h) *)
ClearAll[legendreT];
legendreT[Pmax_, h_] := legendreT[Pmax, h] = Module[{T, n, k, c},
  T = ConstantArray[0, {Pmax + 1, Pmax + 1}];
  Do[
    c = CoefficientList[LegendreP[n, x], x];
    Do[T[[n + 1, k + 1]] = c[[k + 1]] (2/h)^k, {k, 0, n}],
    {n, 0, Pmax}];
  T];

Tx = legendreT[Pmax, aN]; Ty = legendreT[Pmax, bN]; Tz = legendreT[Pmax, cN];

(* Build basis transformation: TblockI[I, J] = T_x[P_I, p_J] T_y[Q_I, q_J] T_z[R_I, r_J]
   indexing same basis tuple, treated as Legendre orders for I, monomial orders for J. *)
TblockL = Table[
  Module[{P, Q, R, p, q, r},
    {P, Q, R} = basisX[[i]];
    {p, q, r} = basisX[[j]];
    Tx[[P + 1, p + 1]] Ty[[Q + 1, q + 1]] Tz[[R + 1, r + 1]]],
  {i, nX}, {j, nX}];

(* K in Legendre = TblockL . KMon . Transpose[TblockL] *)
KLeg = TblockL.KMon.Transpose[TblockL];

(* Mass diagonal in Legendre *)
massDiag = Table[
  Module[{P, Q, R}, {P, Q, R} = basisX[[i]];
    aN/(2 P + 1) bN/(2 Q + 1) cN/(2 R + 1)],
  {i, nX}];

(* M-orthonormal: K_norm[I, J] = K_Leg[I, J] / sqrt(mass[I] mass[J]) *)
KLegNorm = Table[KLeg[[i, j]]/Sqrt[massDiag[[i]] massDiag[[j]]], {i, nX}, {j, nX}];

condLeg = N[Max[massDiag]/Min[massDiag], 6];
Print["  M_leg condition number (= max_M/min_M): ", condLeg];

{eigvalsLeg, eigvecsLegRaw} = Eigensystem[N[KLegNorm, prec]];
(* Transform back to monomial basis *)
eigvecsLeg = Table[
  Transpose[TblockL].(eigvecsLegRaw[[k]]/Sqrt[massDiag]),
  {k, nX}];
Print["  eigvalue solve: ", N[AbsoluteTime[] - t0, 4], " s"];

(* === Compare results === *)
Print[""];
Print["================================================================"];
Print["Comparison: Monomial vs Legendre at P=8"];
Print["================================================================"];

(* Project A_ext on each *)
aExtMon = ConstantArray[0, nX];
idx010 = Position[basisX, {0, 1, 0}, 1][[1, 1]];
aExtMon[[idx010]] = -1/2;

(* Foster amplitudes *)
gMon = (eigvecsMon.MMon.aExtMon);
gLeg = (eigvecsLeg.MMon.aExtMon);

lambdaMon = eigvalsMon/(4 PiN);
lambdaLeg = eigvalsLeg/(4 PiN);

(* Compute alpha_n via both *)
NMoments = 6;
alphasMon = Table[
  (-1)^n (mu0sigmaN)^n/VN Total[gMon^2 lambdaMon^(n - 1)],
  {n, 1, NMoments}];
alphasLeg = Table[
  (-1)^n (mu0sigmaN)^n/VN Total[gLeg^2 lambdaLeg^(n - 1)],
  {n, 1, NMoments}];

Print["alpha_n (x-block only):"];
Do[
  Print["  alpha_", n, " mon = ", N[alphasMon[[n]], 16],
        "  leg = ", N[alphasLeg[[n]], 16],
        "  diff = ", N[alphasMon[[n]] - alphasLeg[[n]], 6]],
  {n, NMoments}];
Print[""];

(* Eigvalue precision: how close are largest few eigenvalues? *)
ordMon = Reverse[Ordering[N[Re[eigvalsMon], 6]]];
ordLeg = Reverse[Ordering[N[Re[eigvalsLeg], 6]]];
Print["Largest 5 eigenvalues:"];
Do[
  Print["  k=", k, " mon = ", N[eigvalsMon[[ordMon[[k]]]], 16],
        "  leg = ", N[eigvalsLeg[[ordLeg[[k]]]], 16]],
  {k, 5}];
