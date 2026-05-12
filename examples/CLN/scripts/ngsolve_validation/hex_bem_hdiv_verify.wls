(* ::Package:: *)

(* hex_bem_hdiv_verify.wls
   Sanity checks on HDiv order 3 BEM result:
   1. K matrix symmetric, M positive definite
   2. Eigenvalue residuals ||K v - lambda M v||
   3. M-orthonormality of eigenvectors
   4. All eigenvalues positive (physical)?
   5. Y(0) consistency check
*)

ClearAll[x, y, z];
$MaxExtraPrecision = 50;

If[!FileExistsQ["hex_bem_hdiv_order3_results.mx"],
  Print["No saved file. Run hex_bem_hdiv_order3_overnight.wls first."]; Quit[]];

Get["hex_bem_hdiv_order3_results.mx"];

sigma = SetPrecision[58 10^6, 30];
mu0 = 4 Pi SetPrecision[10^-7, 30];
{aBox, bBox, cBox} = SetPrecision[{5 10^-3, 2 10^-3, 1 10^-3}, 30];
V = aBox bBox cBox;

Print["==========================================="];
Print["HDiv order 3 BEM result verification"];
Print["==========================================="];
Print["nDOF = ", nDOF];
Print["#eigvals = ", Length[eigvals]];

(* === 1. K, M symmetric === *)
Kasym = Max[Abs[K - Transpose[K]] / (Max[Abs[K]] + 10^-100)];
Masym = Max[Abs[Mphys - Transpose[Mphys]] / Max[Abs[Mphys]]];
Print[];
Print["1. Matrix symmetry:"];
Print["   K asymmetry (relative): ", N[Kasym, 6]];
Print["   M asymmetry (relative): ", N[Masym, 6]];

(* === 2. M positive definite (all eigenvalues of M positive) === *)
Print[];
Print["2. M matrix eigenvalues (should all be positive):"];
mEigs = Sort[Re[Eigenvalues[Mphys]]];
Print["   Smallest M eigenvalue: ", N[mEigs[[1]], 6]];
Print["   Largest M eigenvalue:  ", N[mEigs[[-1]], 6]];
Print["   Number negative: ", Count[mEigs, _?(# < -10^-30 &)]];

(* === 3. K eigenvalues from generalized eigen problem === *)
realEigs = Re[eigvals];
Print[];
Print["3. K v = lambda M v eigenvalues:"];
Print["   Number total: ", Length[eigvals]];
Print["   Number positive (lambda > 1e-30): ", Count[realEigs, _?(# > 10^-30 &)]];
Print["   Number negative: ", Count[realEigs, _?(# < -10^-30 &)]];
Print["   Number near zero: ", Count[realEigs, _?(Abs[#] <= 10^-30 &)]];
posEigs = Sort[Select[realEigs, # > 10^-30 &], Greater];
Print["   Top 10 positive: ", N[posEigs[[1 ;; Min[10, Length[posEigs]]]], 8]];
tausTop = 10^6 sigma * posEigs[[1 ;; Min[10, Length[posEigs]]]];
Print["   Top 10 tau (us): ", N[tausTop, 8]];

(* === 4. Eigenvalue residuals: ||K v - lambda M v||/||K v|| === *)
Print[];
Print["4. Eigenvalue residual check (top 10 by tau, sample 5):"];
sortedIdx = Reverse[Ordering[realEigs]];
Do[
  Module[{idx = sortedIdx[[k]], v, lam, residual, normKv, normMv},
    v = eigvecs[[idx]];
    lam = eigvals[[idx]];
    residual = K . v - lam Mphys . v;
    normKv = Sqrt[Re[v . K . v]];
    Print["   k=", k, ": tau=", N[10^6 sigma Re[lam], 8], " us,  residual norm = ",
          N[Sqrt[Re[residual . residual]] / (Abs[normKv] + 10^-100), 4]]],
  {k, 1, Min[5, Length[sortedIdx]]}];

(* === 5. M-orthogonality of eigenvectors === *)
Print[];
Print["5. M-orthogonality check (top 5 eigenvecs):"];
Do[
  Module[{idx1 = sortedIdx[[i]], idx2 = sortedIdx[[j]], v1, v2, dot},
    v1 = eigvecs[[idx1]] / Sqrt[Re[eigvecs[[idx1]] . Mphys . eigvecs[[idx1]]]];
    v2 = eigvecs[[idx2]] / Sqrt[Re[eigvecs[[idx2]] . Mphys . eigvecs[[idx2]]]];
    dot = Re[v1 . Mphys . v2];
    Print["   <v_", i, ", v_", j, ">_M = ", N[dot, 6]]],
  {i, 1, 3}, {j, 1, 3}];

(* === 6. Y(0) DC value === *)
Print[];
Print["6. Y(0) DC physical check:"];
gMag2 = Re[gVec]^2;
Y0 = Sum[
  If[realEigs[[idx]] > 10^-30,
     gMag2[[idx]] * realEigs[[idx]] sigma,
     0],
  {idx, 1, Length[eigvals]}];
Print["   Y(0) = sigma * sum(g^2 lambda_mine) = ", N[Y0, 8]];

(* Phase 7-11 reference Y(0) approximation *)
(* From cauer_depth_test.py: M_0 = 3.058 (DC moment in their convention) *)
Print["   Phase 7-11 M_0 (in their normalization): 3.058"];
Print["   Direct ratio (my Y(0) / 3.058) = ", N[Y0/3.058, 6]];
Print["   (likely needs normalization adjustment due to basis difference)"];

(* === 7. Mode classification by g^2 distribution === *)
Print[];
Print["7. Mode g^2 distribution analysis:"];
gSorted = Reverse[Sort[gMag2]];
Print["   Largest 10 g^2 values:"];
Do[Print["    g_", k, "^2 = ", N[gSorted[[k]], 6]], {k, 1, Min[10, Length[gSorted]]}];

(* Cumulative coverage *)
totalG2 = Total[gMag2];
cumulative = Accumulate[gSorted] / totalG2;
Print["   Cumulative |g|^2 fraction:"];
Do[Print["     top ", k, " modes: ",
         N[cumulative[[k]], 4]], {k, {1, 3, 5, 10, 20, 50, Length[cumulative]}}];

(* === 8. Eigenvalue spectrum density === *)
Print[];
Print["8. Eigenvalue spectrum density (positive eigvals):"];
posEigsAll = Sort[Select[realEigs, # > 10^-30 &]];
tausAll = 10^6 sigma * posEigsAll;
Print["   Range of tau (us): [", N[Min[tausAll], 6], ", ",
      N[Max[tausAll], 6], "]"];
Print["   Median tau (us): ", N[Median[tausAll], 6]];

Print[];
Print["==========================================="];
Print["Verification complete"];
Print["==========================================="];
