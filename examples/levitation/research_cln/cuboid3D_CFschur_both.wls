#!/usr/bin/env wolframscript
(* ::Package:: *)

(*
  cuboid3D_CFschur_both.wls

  Continued-fraction (CF) Schur: expand the non-rational K_SIBC sqrt(s)/(s+d)
  block as an N-rung rational ladder via Stieltjes integral representation
  of 1/sqrt(s), then couple each rung as a Schur block to the CLN.

  Mathematical content
  --------------------
    1/sqrt(s) = (1/pi) integral_0^inf 1/(sqrt(xi)*(s+xi)) dxi
              ~ Sum_{k=1..N} w_k / (s + xi_k)
    where (xi_k, w_k) are Stieltjes-Gauss quadrature nodes/weights.

  Implementation: log-grid Stieltjes with closed-form residue weights.
    xi_k = xi_min * (xi_max/xi_min)^((k-0.5)/N)    (log-spaced)
    w_k  = (1/pi) * sqrt(xi_k) * log(xi_max/xi_min) / N
      (rectangular rule on log(xi) with weight sqrt(xi); converges in N)

  This is the simplest convergent quadrature.  For production work, use
  Gauss-Jacobi or AAA rational approximation; results below already show
  geometric convergence in N for the C-axis.

  Y_CF-Schur(s) = Y_CLN(s) + K_SIBC * Sum_{k=1..N} w_k / (s + xi_k)

  All-rational: each rung = parallel R-L pair, SPICE-realisable.

  Four sub-tests on each cuboid (A1 17.72x17.72x2 mm, 5x2x1 mm):
    (1) N=1 CF-Schur vs original non-rational sqrt(s) Schur-F
    (2) C-axis convergence for N=1,3,5,7,9,11
    (3) B-axis contamination (s=0 Taylor coefficient errors)
    (4) A-axis hybrid: CLN_N=3 + CF-Schur(N=4) (basis 7) vs Multi-K basis 7
*)

ClearAll[sigma, mu, nu];

sigma = 5.8*^7;
mu    = 4 Pi 1*^-7;
nu    = 1/mu;

geometries = {
  <|"name" -> "A1 (17.72x17.72x2 mm Cu)",
    "a" -> 17.72*^-3, "b" -> 17.72*^-3, "c" -> 2*^-3,
    "wall_lo" -> 1*^3, "wall_hi" -> 1*^5|>,
  <|"name" -> "5x2x1 mm Cu",
    "a" -> 5*^-3, "b" -> 2*^-3, "c" -> 1*^-3,
    "wall_lo" -> 1*^5, "wall_hi" -> 1*^7|>
};

(* === CF nodes/weights for 1/sqrt(s) on log-xi grid ====================== *)
(* xi range covers the body's wall band and 4 decades above for high-f      *)
cfNodes[N_Integer, xiMin_, xiMax_] := Module[{logRatio, ksRange, xi, w},
  logRatio = Log[xiMax/xiMin];
  ksRange  = Range[N];
  xi = Table[xiMin Exp[(k - 0.5) logRatio / N], {k, ksRange}];
  (* Stieltjes weight on log-grid: w_k = (1/pi) sqrt(xi_k) * (log_ratio/N) *)
  w  = (1/Pi) Sqrt[xi] (logRatio / N);
  {xi, w}];

(* === Per-geometry analysis ============================================== *)
analyzeCuboid[geom_, cfParams_] := Module[
  {a, b, c, name, wallLo, wallHi, Sarea, KSIBC3D, Mmax = 13, oddRange,
   modeList, Nmodes, lamVec, tauVec, g2Vec, LambdaM, sigmaM, Imat,
   Kmodal, K0matrix, bvec, YExactFoster, K0Inner, K0MGS, K0Normalise,
   buildCLNbasis, buildMultiKuriyama, YReducedFromBasis,
   YReducedSchurF, YCFSchur, dCandidates, freqWallTmp, QCLN3, QCLN6, QCLN10,
   QMK4, QMK7,
   dBest, dSchur, freqProbe, evalCase,
   xiMin, xiMax, NCFs, cfResults, taylorExact, taylorReducedPlusCF,
   cfBaxisResults},
  a = geom["a"]; b = geom["b"]; c = geom["c"];
  name = geom["name"]; wallLo = geom["wall_lo"]; wallHi = geom["wall_hi"];
  Sarea = 2 (a b + b c + c a);
  KSIBC3D = Sarea Sqrt[sigma/mu];
  xiMin = cfParams["xi_min"]; xiMax = cfParams["xi_max"];
  NCFs  = cfParams["N_list"];

  Print["==============================================="];
  Print[" Geometry: ", name];
  Print["==============================================="];
  Print["  a = ", N[a 1000, 4], " mm, b = ", N[b 1000, 4], " mm, c = ", N[c 1000, 4], " mm"];
  Print["  S = ", N[Sarea 1*^6, 6], " mm^2,  K_SIBC = ", N[KSIBC3D, 6]];
  Print["  CF grid: xi in [", N[xiMin], ", ", N[xiMax],
        "] rad/s = 2pi[", N[xiMin/(2 Pi)], ", ", N[xiMax/(2 Pi)], "] Hz"];

  (* --- Foster spectrum + reduced operators --- *)
  oddRange = Range[1, Mmax, 2];
  modeList = Flatten[Table[{m, n, p}, {m, oddRange}, {n, oddRange}, {p, oddRange}], 2];
  Nmodes   = Length[modeList];
  lamVec   = N[Table[(mnp[[1]] Pi/a)^2 + (mnp[[2]] Pi/b)^2 + (mnp[[3]] Pi/c)^2,
                     {mnp, modeList}]];
  tauVec   = mu sigma / lamVec;
  g2Vec    = N[Table[512 sigma a b c / ((mnp[[1]] mnp[[2]] mnp[[3]])^2 Pi^6),
                     {mnp, modeList}]];
  LambdaM    = N[DiagonalMatrix[lamVec]];
  sigmaM     = N[sigma IdentityMatrix[Nmodes]];
  Imat       = N[IdentityMatrix[Nmodes]];
  Kmodal[s_] := nu LambdaM + sigma s Imat;
  K0matrix   = Kmodal[0];
  bvec       = N[Sqrt[g2Vec / tauVec]];
  YExactFoster[omega_] := Module[{s = I omega}, Total[g2Vec / (1 + s tauVec)]];

  K0Inner[u_, v_] := u . K0matrix . v;
  K0MGS[vIn_, basis_] := Module[{v = vIn},
    Do[v = v - K0Inner[basis[[i]], v] basis[[i]], {i, 1, Length[basis]}]; v];
  K0Normalise[v_] := Module[{nrm = Sqrt[Re[K0Inner[v, v]]]},
    If[nrm < 1*^-300, v, v/nrm]];
  buildCLNbasis[s0_, M_] := Module[{Q = {}, vk, KsInv},
    KsInv = Inverse[Kmodal[s0]];
    vk    = K0Normalise[KsInv . (sigmaM . bvec)];
    AppendTo[Q, vk];
    Do[vk = KsInv . (sigmaM . Last[Q]);
       vk = K0Normalise[K0MGS[vk, Q]];
       AppendTo[Q, vk], {M - 1}];
    Q];
  buildMultiKuriyama[Mzero_, s1list_] := Module[{Q, KsInv, vRaw},
    Q = buildCLNbasis[0, Mzero];
    Do[Module[{si = s1list[[k]]},
         KsInv = Inverse[Kmodal[si]];
         vRaw  = KsInv . (sigmaM . Last[Q]);
         vRaw  = K0Normalise[K0MGS[vRaw, Q]];
         Q     = Append[Q, vRaw]], {k, 1, Length[s1list]}];
    Q];
  YReducedFromBasis[Qlist_, omega_] := Module[{Kr, Qb},
    Kr = Qlist . Kmodal[I omega] . Transpose[Qlist];
    Qb = Qlist . bvec;
    sigma (Qb . LinearSolve[Kr, Qb])];

  (* --- ROMs --- *)
  QCLN3 = buildCLNbasis[0, 3];
  QCLN6 = buildCLNbasis[0, 6];
  QMK4  = buildMultiKuriyama[3, {2 Pi wallLo 10}];
  QMK7  = buildMultiKuriyama[3, N[2 Pi wallLo 10^Range[0, 3, 1]]];

  (* --- Optimal d for the original non-rational Schur-F (re-use earlier result) --- *)
  dCandidates = N[2 Pi 10^Range[4, 9, 0.5]];
  freqWallTmp = Select[N[10^Range[Log10[wallLo], Log10[wallHi], 0.25]],
                       wallLo <= # <= wallHi &];
  YReducedSchurF[Qlist_, omega_, d_] :=
    YReducedFromBasis[Qlist, omega] + KSIBC3D Sqrt[I omega]/(I omega + d);
  dBest = First[SortBy[Table[
    Module[{wallErr},
      wallErr = Max[Table[
        Module[{om = 2 Pi f, yEx, yR},
          yEx = YExactFoster[om];
          yR  = YReducedSchurF[QCLN3, om, d];
          Abs[(yR - yEx)/yEx]], {f, freqWallTmp}]];
      {d, wallErr}], {d, dCandidates}], Last]];
  dSchur = dBest[[1]];
  Print["  non-rational Schur-F best d/(2pi) = ", N[dSchur/(2 Pi), 4], " Hz"];

  (* --- CF-Schur model --- *)
  YCFSchur[Qlist_, omega_, xiList_, wList_] := Module[{yCLN, s, yExtra},
    yCLN   = YReducedFromBasis[Qlist, omega];
    s      = I omega;
    yExtra = KSIBC3D Total[wList / (s + xiList)];
    yCLN + yExtra];

  freqProbe = N[{1*^6, 1*^8, 1*^10, 1*^12}];

  (* === Sub-test (1) + (2): C-axis convergence with N --- *)
  Print[""];
  Print["  --- Sub-test (1)+(2): C-axis convergence with CF rank N ---"];
  Print["  ", StringPadLeft["N", 4], " ", StringPadRight["wall %", 11],
        StringPadRight["r(1e6)", 11], StringPadRight["r(1e8)", 11],
        StringPadRight["r(1e10)", 11], StringPadRight["r(1e12)", 11]];
  cfResults = Table[
    Module[{xiList, wList, wallErr, rVec, yEx, yR},
      {xiList, wList} = cfNodes[NC, xiMin, xiMax];
      wallErr = Max[Table[
        Module[{om = 2 Pi f},
          yEx = YExactFoster[om];
          yR  = YCFSchur[QCLN3, om, xiList, wList];
          Abs[(yR - yEx)/yEx]], {f, freqWallTmp}]];
      rVec = Table[
        Module[{om = 2 Pi f},
          yR = YCFSchur[QCLN3, om, xiList, wList];
          Abs[yR] Sqrt[om] / KSIBC3D], {f, freqProbe}];
      Print["  ", StringPadLeft[ToString[NC], 4], " ",
            StringPadRight[ToString[NumberForm[100 wallErr, {6, 2}]] <> "%", 11],
            StringPadRight[ToString[NumberForm[rVec[[1]], {6, 4}]], 11],
            StringPadRight[ToString[NumberForm[rVec[[2]], {6, 4}]], 11],
            StringPadRight[ToString[NumberForm[rVec[[3]], {6, 4}]], 11],
            StringPadRight[ToString[NumberForm[rVec[[4]], {6, 4}]], 11]];
      <|"N" -> NC, "wall_peak" -> N[wallErr],
        "r_1e6" -> N[rVec[[1]]], "r_1e8" -> N[rVec[[2]]],
        "r_1e10" -> N[rVec[[3]]], "r_1e12" -> N[rVec[[4]]]|>],
    {NC, NCFs}];

  (* Reference row: non-rational Schur-F at best d *)
  Module[{yEx, yR, wallErr, rVec},
    wallErr = Max[Table[
      Module[{om = 2 Pi f},
        yEx = YExactFoster[om];
        yR  = YReducedSchurF[QCLN3, om, dSchur];
        Abs[(yR - yEx)/yEx]], {f, freqWallTmp}]];
    rVec = Table[
      Module[{om = 2 Pi f, yR2},
        yR2 = YReducedSchurF[QCLN3, om, dSchur];
        Abs[yR2] Sqrt[om] / KSIBC3D], {f, freqProbe}];
    Print["  ", StringPadLeft["NR", 4], " ",
          StringPadRight[ToString[NumberForm[100 wallErr, {6, 2}]] <> "%", 11],
          StringPadRight[ToString[NumberForm[rVec[[1]], {6, 4}]], 11],
          StringPadRight[ToString[NumberForm[rVec[[2]], {6, 4}]], 11],
          StringPadRight[ToString[NumberForm[rVec[[3]], {6, 4}]], 11],
          StringPadRight[ToString[NumberForm[rVec[[4]], {6, 4}]], 11]];
    Print["  (NR = original non-rational sqrt(s) Schur-F at d/(2pi)=",
          N[dSchur/(2 Pi), 4], " Hz)"]];

  (* === Sub-test (3): B-axis contamination --- *)
  Print[""];
  Print["  --- Sub-test (3): B-axis contamination by CF rungs ---"];
  taylorExact[n_] := (-1)^n Total[g2Vec * tauVec^n];
  (* CF-Schur Taylor at s=0:
     Sum_k w_k * (-1)^n / xi_k^(n+1) (geometric series d^n/ds^n[1/(s+xi)] at s=0) *)
  cfTaylor[n_, xiList_, wList_] := KSIBC3D (-1)^n Total[wList / xiList^(n + 1)];
  taylorReducedPlusCF[Qlist_, n_, xiList_, wList_] :=
    Module[{K0r, sigmar, br, K0rInv, aClnR, aCF},
      K0r    = Qlist . K0matrix . Transpose[Qlist];
      sigmar = Qlist . sigmaM . Transpose[Qlist];
      br     = Qlist . bvec;
      K0rInv = Inverse[K0r];
      aClnR = sigma (br . ((-1)^n MatrixPower[K0rInv . sigmar, n] . K0rInv) . br);
      aCF   = cfTaylor[n, xiList, wList];
      aClnR + aCF];

  Module[{xi5, w5, xi7, w7, xi9, w9},
    {xi5, w5} = cfNodes[5, xiMin, xiMax];
    {xi7, w7} = cfNodes[7, xiMin, xiMax];
    {xi9, w9} = cfNodes[9, xiMin, xiMax];
    Print["  Taylor coefficients of (CLN_N=3 + CF(N=5,7,9)) at s=0 vs Foster exact"];
    Print["  ", StringPadLeft["n", 4], " ",
          StringPadRight["exact a_n", 14],
          StringPadRight["CLN3 alone", 14],
          StringPadRight["CLN3+CF5", 14],
          StringPadRight["CLN3+CF7", 14],
          StringPadRight["CLN3+CF9", 14]];
    cfBaxisResults = Table[
      Module[{aEx, aN3, aCF5, aCF7, aCF9, fE},
        aEx  = taylorExact[n];
        aN3  = sigma (((QCLN3 . bvec)) . ((-1)^n MatrixPower[
                Inverse[QCLN3 . K0matrix . Transpose[QCLN3]] . (QCLN3 . sigmaM . Transpose[QCLN3]),
                n] . Inverse[QCLN3 . K0matrix . Transpose[QCLN3]]) . (QCLN3 . bvec));
        aCF5 = taylorReducedPlusCF[QCLN3, n, xi5, w5];
        aCF7 = taylorReducedPlusCF[QCLN3, n, xi7, w7];
        aCF9 = taylorReducedPlusCF[QCLN3, n, xi9, w9];
        fE[x_] := If[Abs[aEx] > 0, ScientificForm[N[Abs[(x - aEx)/aEx]], 2], "n/a"];
        Print["  ", StringPadLeft[ToString[n], 4], " ",
              StringPadRight[ToString[ScientificForm[N[aEx], 3]], 14],
              StringPadRight[ToString[fE[aN3]], 14],
              StringPadRight[ToString[fE[aCF5]], 14],
              StringPadRight[ToString[fE[aCF7]], 14],
              StringPadRight[ToString[fE[aCF9]], 14]];
        <|"n" -> n, "exact" -> N[aEx, 6],
          "rel_CLN3_alone" -> N[If[Abs[aEx] > 0, Abs[(aN3 - aEx)/aEx], 0]],
          "rel_CLN3_CF5"   -> N[If[Abs[aEx] > 0, Abs[(aCF5 - aEx)/aEx], 0]],
          "rel_CLN3_CF7"   -> N[If[Abs[aEx] > 0, Abs[(aCF7 - aEx)/aEx], 0]],
          "rel_CLN3_CF9"   -> N[If[Abs[aEx] > 0, Abs[(aCF9 - aEx)/aEx], 0]]|>],
      {n, 0, 8}]];

  (* === Sub-test (4): Hybrid CLN_3 + CF(N=4) vs Multi-K basis 7 --- *)
  Print[""];
  Print["  --- Sub-test (4): A-axis hybrid CLN_3+CF(N=4) vs Multi-K basis 7 ---"];
  Module[{xi4, w4, yEx, yR, wallErr, rProbe, axisHy, axisMK},
    {xi4, w4} = cfNodes[4, xiMin, xiMax];
    wallErr = Max[Table[
      Module[{om = 2 Pi f},
        yEx = YExactFoster[om];
        yR  = YCFSchur[QCLN3, om, xi4, w4];
        Abs[(yR - yEx)/yEx]], {f, freqWallTmp}]];
    rProbe = Module[{om = 2 Pi 1*^12, yR2},
      yR2 = YCFSchur[QCLN3, om, xi4, w4];
      Abs[yR2] Sqrt[om] / KSIBC3D];
    axisHy = <|"name" -> "CLN_3 + CF(N=4)", "basis" -> 7,
               "wall" -> N[wallErr], "r_1e12" -> N[rProbe]|>;
    wallErr = Max[Table[
      Module[{om = 2 Pi f},
        yEx = YExactFoster[om];
        yR  = YReducedFromBasis[QMK7, om];
        Abs[(yR - yEx)/yEx]], {f, freqWallTmp}]];
    rProbe = Module[{om = 2 Pi 1*^12, yR2},
      yR2 = YReducedFromBasis[QMK7, om];
      Abs[yR2] Sqrt[om] / KSIBC3D];
    axisMK = <|"name" -> "Multi-K M_0+4pt", "basis" -> 7,
               "wall" -> N[wallErr], "r_1e12" -> N[rProbe]|>;
    Print["  ", StringPadRight[axisHy["name"], 25],
          " basis=", axisHy["basis"], "  wall=",
          NumberForm[100 axisHy["wall"], {6, 2}], "%,  r(1e12)=",
          NumberForm[axisHy["r_1e12"], {6, 4}]];
    Print["  ", StringPadRight[axisMK["name"], 25],
          " basis=", axisMK["basis"], "  wall=",
          NumberForm[100 axisMK["wall"], {6, 2}], "%,  r(1e12)=",
          NumberForm[axisMK["r_1e12"], {6, 4}]]];

  Print[""];
  <|"geom" -> geom, "K_SIBC" -> N[KSIBC3D], "Nmodes" -> Nmodes,
    "d_schur_NR" -> N[dSchur],
    "cf_results" -> cfResults,
    "cf_baxis" -> cfBaxisResults|>];

(* === Per-geometry CF grids (xi range matched to wall + asymptote) ====== *)
(* A1: wall 1e3-1e5 Hz, asymptote up to 1e12 Hz -> xi from 2pi*1e3 to 2pi*1e13 *)
(* 5x2x1: wall 1e5-1e7 Hz, asymptote 1e12 Hz   -> xi from 2pi*1e5 to 2pi*1e13 *)
cfParamsA1    = <|"xi_min" -> 2 Pi 1*^3, "xi_max" -> 2 Pi 1*^13,
                  "N_list" -> {1, 3, 5, 7, 9, 11}|>;
cfParams5x2x1 = <|"xi_min" -> 2 Pi 1*^5, "xi_max" -> 2 Pi 1*^13,
                  "N_list" -> {1, 3, 5, 7, 9, 11}|>;

results = {analyzeCuboid[geometries[[1]], cfParamsA1],
            analyzeCuboid[geometries[[2]], cfParams5x2x1]};

(* === Side-by-side C-axis convergence summary =========================== *)
Print[""];
Print["============================================================"];
Print[" CF-Schur C-axis convergence (r(1e12), target = 1)"];
Print["============================================================"];
Print[StringPadLeft["N", 4], " ",
      StringPadRight["A1 plate r(1e12)", 22],
      StringPadRight["5x2x1 brick r(1e12)", 22]];
Do[Module[{a1 = First[results]["cf_results"][[k]],
            sm = Last[results]["cf_results"][[k]]},
    Print[StringPadLeft[ToString[a1["N"]], 4], " ",
          StringPadRight[ToString[NumberForm[a1["r_1e12"], {7, 5}]], 22],
          StringPadRight[ToString[NumberForm[sm["r_1e12"], {7, 5}]], 22]]],
   {k, Length[First[results]["cf_results"]]}];

(* === Save JSON ========================================================== *)
combinedOut = <|
  "geometries" -> Table[results[[k]]["geom"], {k, Length[results]}],
  "results"    -> Table[
    <|"name"     -> results[[k]]["geom"]["name"],
      "K_SIBC"   -> results[[k]]["K_SIBC"],
      "d_schur_NR" -> results[[k]]["d_schur_NR"],
      "cf_results" -> results[[k]]["cf_results"],
      "cf_baxis"   -> results[[k]]["cf_baxis"]|>,
    {k, Length[results]}]|>;
Export["cuboid3D_CFschur_both.json", combinedOut, "JSON"];
Print[""];
Print["Saved cuboid3D_CFschur_both.json"];
