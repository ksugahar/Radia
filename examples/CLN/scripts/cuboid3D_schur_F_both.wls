#!/usr/bin/env wolframscript
(* ::Package:: *)

(*
  cuboid3D_schur_F_both.wls

  Run the 3-axis Schur-F verification on two cuboid geometries:
    (A) A1 cuboid: 17.72 x 17.72 x 2 mm Cu  -- plate-like, matches
                   Nagamine A1 reference problem (phase1 VIM scan)
    (B) 5 x 2 x 1 mm Cu                       -- rectangular brick

  Both use the scalar 3D Dirichlet diffusion model (cf.
  kuriyama_cuboid3D_Ytail.wls).  The point: K_SIBC = S sqrt(sigma/mu)
  with S = total surface area is universal across aspect ratios as
  different as 9:9:1 (A1) and 5:2:1.

  For each geometry, computes:
    - Foster spectrum truncated to (m,n,p odd, max 13) -> 343 modes
    - CLN N=3, N=6, N=10 at s_0 = 0
    - Multi-Kuriyama M_0=3 + k log-spaced points
    - Schur-F (CLN_N=3 + 1 sqrt(s) block) at optimised d
  Reports 3-axis table (wall band peak err, asymptote ratio).
*)

ClearAll[sigma, mu, nu];

sigma = 5.8*^7;
mu    = 4 Pi 1*^-7;
nu    = 1/mu;

geometries = {
  <|"name" -> "A1 (17.72x17.72x2 mm Cu)",
    "a" -> 17.72*^-3, "b" -> 17.72*^-3, "c" -> 2*^-3,
    "wall_lo" -> 1*^3, "wall_hi" -> 1*^5|>,  (* A1 wall band is lower (larger geometry) *)
  <|"name" -> "5x2x1 mm Cu",
    "a" -> 5*^-3, "b" -> 2*^-3, "c" -> 1*^-3,
    "wall_lo" -> 1*^5, "wall_hi" -> 1*^7|>
};

(* === Per-geometry analysis ============================================ *)
analyzeCuboid[geom_] := Module[
  {a, b, c, name, wallLo, wallHi, Sarea, KSIBC3D, Mmax = 13, oddRange,
   modeList, Nmodes, lamVec, tauVec, g2Vec, LambdaM, sigmaM, Imat,
   Kmodal, K0matrix, bvec, YExactFoster, K0Inner, K0MGS, K0Normalise,
   buildCLNbasis, buildMultiKuriyama, YReducedFromBasis,
   YReducedSchurF, dCandidates, freqWallTmp, QCLN3, QCLN6, QCLN10,
   QMK4, QMK7, QMK10, dSweepResults, dBest, dSchur,
   freqProbe, evalCase, axisAC, taylorExact, taylorReduced, bTable},
  a = geom["a"]; b = geom["b"]; c = geom["c"];
  name = geom["name"]; wallLo = geom["wall_lo"]; wallHi = geom["wall_hi"];
  Sarea = 2 (a b + b c + c a);
  KSIBC3D = Sarea Sqrt[sigma/mu];

  Print["==============================================="];
  Print[" Geometry: ", name];
  Print["==============================================="];
  Print["  a = ", N[a 1000, 4], " mm, b = ", N[b 1000, 4], " mm, c = ", N[c 1000, 4], " mm"];
  Print["  S = 2(ab+bc+ca) = ", N[Sarea 1*^6, 6], " mm^2"];
  Print["  K_SIBC = S sqrt(sigma/mu) = ", N[KSIBC3D, 6]];
  Print["  wall band: ", N[wallLo], " .. ", N[wallHi], " Hz"];

  oddRange = Range[1, Mmax, 2];
  modeList = Flatten[Table[{m, n, p}, {m, oddRange}, {n, oddRange}, {p, oddRange}], 2];
  Nmodes   = Length[modeList];
  lamVec   = N[Table[(mnp[[1]] Pi/a)^2 + (mnp[[2]] Pi/b)^2 + (mnp[[3]] Pi/c)^2,
                     {mnp, modeList}]];
  tauVec   = mu sigma / lamVec;
  g2Vec    = N[Table[512 sigma a b c / ((mnp[[1]] mnp[[2]] mnp[[3]])^2 Pi^6),
                     {mnp, modeList}]];
  Print["  Nmodes = ", Nmodes];

  LambdaM    = N[DiagonalMatrix[lamVec]];
  sigmaM     = N[sigma IdentityMatrix[Nmodes]];
  Imat       = N[IdentityMatrix[Nmodes]];
  Kmodal[s_] := nu LambdaM + sigma s Imat;
  K0matrix   = Kmodal[0];
  bvec       = N[Sqrt[g2Vec / tauVec]];

  YExactFoster[omega_] := Module[{s = I omega}, Total[g2Vec / (1 + s tauVec)]];

  K0Inner[u_, v_] := u . K0matrix . v;
  K0MGS[vIn_, basis_] := Module[{v = vIn, cc},
    Do[cc = K0Inner[basis[[i]], v];
       v  = v - cc basis[[i]], {i, 1, Length[basis]}];
    v];
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

  YReducedSchurF[Qlist_, omega_, d_] :=
    YReducedFromBasis[Qlist, omega] + KSIBC3D Sqrt[I omega]/(I omega + d);

  QCLN3   = buildCLNbasis[0, 3];
  QCLN6   = buildCLNbasis[0, 6];
  QCLN10  = buildCLNbasis[0, 10];
  (* Multi-K expansion points span 2 decades above wallLo *)
  QMK4    = buildMultiKuriyama[3, {2 Pi wallLo 10}];
  QMK7    = buildMultiKuriyama[3, N[2 Pi wallLo 10^Range[0, 3, 1]]];
  QMK10   = buildMultiKuriyama[3, N[2 Pi wallLo 10^Range[-1, 5, 1]]];

  (* d sweep: spans 10^4 to 10^9 Hz *)
  dCandidates = N[2 Pi 10^Range[4, 9, 0.5]];
  freqWallTmp = Select[N[10^Range[Log10[wallLo], Log10[wallHi], 0.25]],
                       wallLo <= # <= wallHi &];
  dSweepResults = Table[
    Module[{wallErr},
      wallErr = Max[Table[
        Module[{om = 2 Pi f, yEx, yR},
          yEx = YExactFoster[om];
          yR  = YReducedFromBasis[QCLN3, om] + KSIBC3D Sqrt[I om]/(I om + d);
          Abs[(yR - yEx)/yEx]], {f, freqWallTmp}]];
      <|"d" -> d, "wall_peak" -> N[wallErr]|>],
    {d, dCandidates}];
  dBest  = First[SortBy[dSweepResults, #["wall_peak"] &]];
  dSchur = dBest["d"];
  Print["  Best d/(2pi) = ", N[dSchur/(2 Pi), 4], " Hz, Schur-F wall = ",
        NumberForm[100 dBest["wall_peak"], {6, 2}], "%"];

  freqProbe = N[{1*^8, 1*^10, 1*^12}];
  evalCase[caseName_, modelFn_] := Module[{wallErr, rProbe},
    wallErr = Max[Table[
      Module[{om = 2 Pi f, yEx, yR},
        yEx = YExactFoster[om];
        yR  = modelFn[om];
        Abs[(yR - yEx)/yEx]], {f, freqWallTmp}]];
    rProbe = Table[
      Module[{om = 2 Pi f, yR},
        yR = modelFn[om];
        Abs[yR] Sqrt[om] / KSIBC3D], {f, freqProbe}];
    <|"name" -> caseName, "wall_peak" -> N[wallErr],
      "r_1e8" -> N[rProbe[[1]]], "r_1e10" -> N[rProbe[[2]]],
      "r_1e12" -> N[rProbe[[3]]]|>];

  axisAC = {
    evalCase["CLN_N3",          YReducedFromBasis[QCLN3, #] &],
    evalCase["CLN_N6",          YReducedFromBasis[QCLN6, #] &],
    evalCase["CLN_N10",         YReducedFromBasis[QCLN10, #] &],
    evalCase["Multi-K M3+1pt",  YReducedFromBasis[QMK4, #] &],
    evalCase["Multi-K M3+4pt",  YReducedFromBasis[QMK7, #] &],
    evalCase["Multi-K M3+7pt",  YReducedFromBasis[QMK10, #] &],
    evalCase["Schur-F (CLN3+1)", YReducedSchurF[QCLN3, #, dSchur] &]};

  Print[""];
  Print["  --- A/C axes ---"];
  Print["  ", StringPadRight["case", 22], StringPadRight["basis", 7],
        StringPadRight["wall %", 11], StringPadRight["r(1e8)", 11],
        StringPadRight["r(1e10)", 11], StringPadRight["r(1e12)", 11]];
  basisSizes = {3, 6, 10, 4, 7, 10, 4};
  Do[Module[{r = axisAC[[k]]},
      Print["  ", StringPadRight[r["name"], 22],
            StringPadRight[ToString[basisSizes[[k]]], 7],
            StringPadRight[ToString[NumberForm[100 r["wall_peak"], {6, 2}]] <> "%", 11],
            StringPadRight[ToString[NumberForm[r["r_1e8"], {6, 4}]], 11],
            StringPadRight[ToString[NumberForm[r["r_1e10"], {6, 4}]], 11],
            StringPadRight[ToString[NumberForm[r["r_1e12"], {6, 4}]], 11]]],
     {k, Length[axisAC]}];

  (* B axis (Taylor at s=0) *)
  taylorExact[n_]      := (-1)^n Total[g2Vec * tauVec^n];
  taylorReduced[Qlist_, n_] := Module[{K0r, sigmar, br, K0rInv},
    K0r    = Qlist . K0matrix . Transpose[Qlist];
    sigmar = Qlist . sigmaM . Transpose[Qlist];
    br     = Qlist . bvec;
    K0rInv = Inverse[K0r];
    sigma (br . ((-1)^n MatrixPower[K0rInv . sigmar, n] . K0rInv) . br)];
  Print[""];
  Print["  --- B axis (Taylor coefficient rel-err vs exact at s=0) ---"];
  Print["  ", StringPadLeft["n", 4], " ",
        StringPadRight["CLN3", 12], StringPadRight["CLN6", 12],
        StringPadRight["MK M3+4pt", 12], StringPadRight["MK M3+7pt", 12]];
  bTable = Table[
    Module[{aEx, aN3, aN6, aMK7, aMK10, fE},
      aEx   = taylorExact[n];
      aN3   = taylorReduced[QCLN3, n];
      aN6   = taylorReduced[QCLN6, n];
      aMK7  = taylorReduced[QMK7, n];
      aMK10 = taylorReduced[QMK10, n];
      fE[x_] := If[Abs[aEx] > 0, ScientificForm[N[Abs[(x - aEx)/aEx]], 2], "n/a"];
      Print["  ", StringPadLeft[ToString[n], 4], " ",
            StringPadRight[ToString[fE[aN3]], 12],
            StringPadRight[ToString[fE[aN6]], 12],
            StringPadRight[ToString[fE[aMK7]], 12],
            StringPadRight[ToString[fE[aMK10]], 12]];
      <|"n" -> n, "exact" -> N[aEx, 6],
        "rel_CLN_N3" -> N[If[Abs[aEx] > 0, Abs[(aN3 - aEx)/aEx], 0]],
        "rel_CLN_N6" -> N[If[Abs[aEx] > 0, Abs[(aN6 - aEx)/aEx], 0]],
        "rel_MultiK_M3p4" -> N[If[Abs[aEx] > 0, Abs[(aMK7 - aEx)/aEx], 0]],
        "rel_MultiK_M3p7" -> N[If[Abs[aEx] > 0, Abs[(aMK10 - aEx)/aEx], 0]]|>],
    {n, 0, 12}];

  Print[""];
  <|"geom" -> geom, "K_SIBC" -> N[KSIBC3D], "Nmodes" -> Nmodes,
    "d_best" -> N[dSchur], "axisAC" -> axisAC, "bTable" -> bTable|>];

(* === Run both geometries ============================================= *)
results = Map[analyzeCuboid, geometries];

(* === Side-by-side summary ============================================ *)
Print[""];
Print["============================================================"];
Print[" Universality of K_SIBC = S sqrt(sigma/mu) and Schur-F"];
Print["============================================================"];
Print[""];
Print[StringPadRight["", 22],
      StringPadRight["A1 (17.72x17.72x2 mm)", 25],
      "5x2x1 mm"];
Do[Module[{aA1 = First[results]["axisAC"][[k]],
            aSm = Last[results]["axisAC"][[k]]},
    Print[StringPadRight[aA1["name"], 22],
          StringPadRight[
            "wall=" <> ToString[NumberForm[100 aA1["wall_peak"], {6, 2}]] <>
              "%, r(1e12)=" <> ToString[NumberForm[aA1["r_1e12"], {6, 4}]],
            25],
          "wall=" <> ToString[NumberForm[100 aSm["wall_peak"], {6, 2}]] <>
            "%, r(1e12)=" <> ToString[NumberForm[aSm["r_1e12"], {6, 4}]]]],
   {k, 7}];

Print[""];
Print["K_SIBC^{A1}    = ", N[First[results]["K_SIBC"], 6]];
Print["K_SIBC^{5x2x1} = ", N[Last[results]["K_SIBC"], 6]];
Print["ratio          = ", N[First[results]["K_SIBC"]/Last[results]["K_SIBC"], 6],
      " (matches surface area ratio ",
      N[(2(17.72^2 + 17.72*2 + 17.72*2))/(2*(5*2 + 2*1 + 5*1)), 6], ")"];

(* === Save combined JSON ============================================== *)
combinedOut = <|
  "geometries" -> Table[results[[k]]["geom"], {k, Length[results]}],
  "results"    -> Table[
    <|"name"     -> results[[k]]["geom"]["name"],
      "K_SIBC"   -> results[[k]]["K_SIBC"],
      "d_best"   -> results[[k]]["d_best"],
      "axisAC"   -> results[[k]]["axisAC"],
      "bTable"   -> results[[k]]["bTable"]|>,
    {k, Length[results]}]|>;
Export["cuboid3D_schur_F_both.json", combinedOut, "JSON"];
Print[""];
Print["Saved cuboid3D_schur_F_both.json"];
