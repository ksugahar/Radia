#!/usr/bin/env wolframscript
(* ::Package:: *)

(*
  cuboid3D_schur_F.wls

  3-axis ROM evaluation for the 3D cuboid scalar diffusion problem,
  parallel to slab_cpe_schur.wls.

  Reference problem (scalar Dirichlet, uniform source):
     -Laplacian u + s mu sigma u = -s mu sigma   on  Omega = [0,a]x[0,b]x[0,c]
     u = 0  on  dOmega
     u_ext = 1

  Foster spectrum (odd m,n,p only; even modes uncouple):
     lambda_{mnp} = (m pi/a)^2 + (n pi/b)^2 + (p pi/c)^2
     tau_{mnp}    = mu sigma / lambda_{mnp}
     g_{mnp}^2    = 512 a b c / ((m n p)^2 pi^6)

  Foster admittance:
     Y(s) = sigma Sum_{odd m,n,p} g_{mnp}^2 / (1 + s tau_{mnp})

  Asymptote (verified analytically in kuriyama_cuboid3D_Ytail.wls):
     Y -> K_SIBC^{3D} / sqrt(s)
     K_SIBC^{3D} = 2(ab+bc+ca) sqrt(sigma/mu)     (surface area * eta_s)

  ROMs compared on three axes:
     A: wall-band peak rel-err (1e5..1e7 Hz)
     B: Taylor coefficients of Y_R(s) at s=0 vs exact
     C: asymptote ratio r(f) = |Y_R| / (K_SIBC^{3D}/sqrt(omega))

  Models:
     CLN_N3 / CLN_N6 / CLN_N10 at s_0 = 0  (canonical Cauer Pade)
     Multi-Kuriyama M0=3 + k log-spaced expansion points
     Schur-F (CLN_N3 + 1 sqrt(s) block, basis 4) -- the test

  Run: wolframscript -file cuboid3D_schur_F.wls
*)

ClearAll[a, b, c, sigma, mu, nu, Mmax];

(* Geometry: cuboid 5 x 2 x 1 mm Cu (matches kuriyama_cuboid3D_Ytail.wls) *)
a       = 5*^-3;
b       = 2*^-3;
c       = 1*^-3;
sigma   = 5.8*^7;
mu      = 4 Pi 1*^-7;
nu      = 1/mu;

(* Surface SIBC constant (scalar trace) *)
KSIBC3D = 2 (a b + b c + c a) Sqrt[sigma/mu];
Sarea   = 2 (a b + b c + c a);

Print["=== cuboid3D_schur_F.wls ==="];
Print["a = ", N[a 1000], " mm, b = ", N[b 1000], " mm, c = ", N[c 1000], " mm"];
Print["Surface S = 2(ab+bc+ca) = ", N[Sarea 1*^6, 6], " mm^2"];
Print["sigma = ", sigma, " S/m, mu = ", N[mu, 6], " H/m"];
Print["K_SIBC^{3D} = S sqrt(sigma/mu) = ", N[KSIBC3D, 8]];

(* === Foster spectrum truncation ======================================== *)
(* odd modes m,n,p = 1,3,5,...,Mmax. Total modes = ((Mmax+1)/2)^3 *)
Mmax = 13;   (* 7 x 7 x 7 = 343 modes, plenty *)
oddRange = Range[1, Mmax, 2];

(* Foster eigenvalues, time constants, couplings *)
(* g^2 ABSORBS sigma to match slab convention: Y(s) = Sum g^2/(1+s*tau)     *)
lamMNP[m_, n_, p_] := N[(m Pi/a)^2 + (n Pi/b)^2 + (p Pi/c)^2];
tauMNP[m_, n_, p_] := mu sigma / lamMNP[m, n, p];
g2MNP[m_, n_, p_]  := N[512 sigma a b c / ((m n p)^2 Pi^6)];

(* Enumerate modes (flat index) *)
modeList = Flatten[Table[{m, n, p}, {m, oddRange}, {n, oddRange}, {p, oddRange}], 2];
Nmodes   = Length[modeList];
lamVec   = N[Table[lamMNP[mnp[[1]], mnp[[2]], mnp[[3]]], {mnp, modeList}]];
tauVec   = mu sigma / lamVec;
g2Vec    = N[Table[g2MNP[mnp[[1]], mnp[[2]], mnp[[3]]], {mnp, modeList}]];

Print["Nmodes = ", Nmodes, " (m,n,p odd up to ", Mmax, ")"];

(* === Modal operators in Foster basis ==================================== *)
LambdaM  = N[DiagonalMatrix[lamVec]];
sigmaM   = N[sigma IdentityMatrix[Nmodes]];
Imat     = N[IdentityMatrix[Nmodes]];

Kmodal[s_] := nu LambdaM + sigma s Imat;
K0matrix   = Kmodal[0];

(* Source vector: uniform u_ext = 1 -> RHS proportional to projection onto each mode *)
(* For scalar Dirichlet with uniform source, b_k ~ sqrt(g_k^2/tau_k) (analogous to slab) *)
bvec = N[Sqrt[g2Vec / tauVec]];

(* Exact Foster admittance (truncated triple sum, our "ground truth") *)
(* sigma absorbed into g2; Y(s) = Sum g^2/(1+s*tau)                         *)
YExactFoster[omega_] := Module[{s = I omega},
  Total[g2Vec / (1 + s tauVec)]];

(* === K_0 inner product utilities ======================================= *)
K0Inner[u_, v_]     := u . K0matrix . v;
K0MGS[vIn_, basis_] := Module[{v = vIn, ai, cc},
  Do[ai = basis[[i]];
     cc = K0Inner[ai, v];
     v  = v - cc ai, {i, 1, Length[basis]}];
  v];
K0Normalise[v_] := Module[{nrm = Sqrt[Re[K0Inner[v, v]]]},
  If[nrm < 1*^-300, v, v/nrm]];

(* === CLN Krylov basis at s = s0 ======================================== *)
buildCLNbasis[s0_, M_] := Module[{Q = {}, vk, KsInv},
  KsInv = Inverse[Kmodal[s0]];
  vk    = KsInv . (sigmaM . bvec);
  vk    = K0Normalise[vk];
  AppendTo[Q, vk];
  Do[vk = KsInv . (sigmaM . Last[Q]);
     vk = K0MGS[vk, Q];
     vk = K0Normalise[vk];
     AppendTo[Q, vk], {M - 1}];
  Q];

(* === Multi-Kuriyama: s0=0 plus k log-spaced finite points =============== *)
buildMultiKuriyama[Mzero_, s1list_] := Module[{Q, KsInv, vRaw},
  Q = buildCLNbasis[0, Mzero];
  Do[Module[{si = s1list[[k]]},
       KsInv = Inverse[Kmodal[si]];
       vRaw  = KsInv . (sigmaM . Last[Q]);
       vRaw  = K0MGS[vRaw, Q];
       vRaw  = K0Normalise[vRaw];
       Q     = Append[Q, vRaw]], {k, 1, Length[s1list]}];
  Q];

(* === Reduced admittance (no Schur) ==================================== *)
YReducedFromBasis[Qlist_, omega_] := Module[{Q, Kr, Qb},
  Q  = Qlist;
  Kr = Q . Kmodal[I omega] . Transpose[Q];
  Qb = Q . bvec;
  sigma (Qb . LinearSolve[Kr, Qb])];

(* === Schur-F 1x1 termination ========================================== *)
(* Y_R(s) = Y_CLN(s) + K_SIBC sqrt(s) / (s + d)                          *)
YReducedSchurF[Qlist_, omega_, d_] := Module[{yCLN, s, yExtra},
  yCLN   = YReducedFromBasis[Qlist, omega];
  s      = I omega;
  yExtra = KSIBC3D Sqrt[s] / (s + d);
  yCLN + yExtra];

(* === Sanity: DC and asymptote of exact ================================ *)
Print[""];
Print["=== Sanity check on exact Foster Y(s) ==="];
yDC = YExactFoster[0.0];
Print["Y(0) (Foster sum, truncated)  = ", N[yDC, 8]];
yHF = YExactFoster[2 Pi 1*^12];
KempHF = Abs[yHF] Sqrt[2 Pi 1*^12];
Print["|Y(2pi 1e12)| sqrt(omega) num = ", N[KempHF, 8]];
Print["K_SIBC^{3D} formula           = ", N[KSIBC3D, 8]];
Print["asymptote ratio (num/formula) = ", N[KempHF / KSIBC3D, 6]];

(* === Build ROMs ======================================================== *)
Print[""];
Print["=== Building ROMs ==="];
QCLN3   = buildCLNbasis[0, 3];
QCLN6   = buildCLNbasis[0, 6];
QCLN10  = buildCLNbasis[0, 10];
QMK4    = buildMultiKuriyama[3, {2 Pi 1*^6}];
QMK7    = buildMultiKuriyama[3, N[2 Pi 10^Range[5, 8, 1]]];      (* 4 pts: 1e5..1e8 *)
QMK10   = buildMultiKuriyama[3, N[2 Pi 10^Range[4, 10, 1]]];     (* 7 pts: 1e4..1e10 *)
Print["  CLN_N3  basis size = ", Length[QCLN3]];
Print["  CLN_N6  basis size = ", Length[QCLN6]];
Print["  CLN_N10 basis size = ", Length[QCLN10]];
Print["  Multi-K M3+1pt basis size = ", Length[QMK4]];
Print["  Multi-K M3+4pt basis size = ", Length[QMK7]];
Print["  Multi-K M3+7pt basis size = ", Length[QMK10]];

(* Sweep d to find Schur-F sweet spot for THIS geometry *)
dCandidates = N[2 Pi 10^Range[4, 9, 0.5]];
Print[""];
Print["=== d sweep (wall band fit for Schur-F basis 4) ==="];
Print[StringPadRight["d/(2pi) [Hz]", 16] <> StringPadRight["wall peak %", 14] <>
      StringPadRight["r(1e12)", 12]];
freqWallTmp = Select[N[10^Range[5, 7, 0.25]], 1*^5 <= # <= 1*^7 &];
dSweepResults = Table[
  Module[{wallErr, rAsymp},
    wallErr = Max[Table[
      Module[{om = 2 Pi f, yEx, yR},
        yEx = YExactFoster[om];
        yR  = YReducedFromBasis[QCLN3, om] + KSIBC3D Sqrt[I om]/(I om + d);
        Abs[(yR - yEx)/yEx]], {f, freqWallTmp}]];
    rAsymp = Module[{om = 2 Pi 1*^12, yR},
      yR = YReducedFromBasis[QCLN3, om] + KSIBC3D Sqrt[I om]/(I om + d);
      Abs[yR] Sqrt[om] / KSIBC3D];
    Print[StringPadRight[ToString[NumberForm[d/(2 Pi), {6, 0}]], 16],
          StringPadRight[ToString[NumberForm[100 wallErr, {6, 2}]] <> "%", 14],
          StringPadRight[ToString[NumberForm[N[rAsymp], {6, 4}]], 12]];
    <|"d_over_2pi" -> N[d/(2 Pi)], "wall_peak" -> N[wallErr], "r_1e12" -> N[rAsymp]|>],
  {d, dCandidates}];

dBest = First[SortBy[dSweepResults, #["wall_peak"] &]];
Print["Best d/(2pi) = ", dBest["d_over_2pi"], " Hz, wall = ",
      NumberForm[100 dBest["wall_peak"], {6, 2}], "%"];
dSchur = 2 Pi dBest["d_over_2pi"];

(* === A-axis: wall-band peak rel-err =================================== *)
(* Frequency grid covering wall (1e5..1e7) + asymptote probes 1e8,1e10,1e12 *)
freqList = N[10^Range[1, 12, 0.25]];
freqWall = Select[freqList, 1*^5 <= # <= 1*^7 &];
freqProbe = N[{1*^8, 1*^10, 1*^12}];

evalCase[name_, modelFn_] := Module[{wallErr, rProbe},
  wallErr = Max[Table[
    Module[{om = 2 Pi f, yEx, yR},
      yEx = YExactFoster[om];
      yR  = modelFn[om];
      Abs[(yR - yEx)/yEx]], {f, freqWall}]];
  rProbe = Table[
    Module[{om = 2 Pi f, yR},
      yR = modelFn[om];
      Abs[yR] Sqrt[om] / KSIBC3D], {f, freqProbe}];
  <|"name" -> name, "wall_peak" -> N[wallErr],
    "r_1e8" -> N[rProbe[[1]]],
    "r_1e10" -> N[rProbe[[2]]],
    "r_1e12" -> N[rProbe[[3]]]|>];

axisAC = {
  evalCase["CLN_N3",         YReducedFromBasis[QCLN3, #] &],
  evalCase["CLN_N6",         YReducedFromBasis[QCLN6, #] &],
  evalCase["CLN_N10",        YReducedFromBasis[QCLN10, #] &],
  evalCase["Multi-K M3+1pt", YReducedFromBasis[QMK4, #] &],
  evalCase["Multi-K M3+4pt", YReducedFromBasis[QMK7, #] &],
  evalCase["Multi-K M3+7pt", YReducedFromBasis[QMK10, #] &],
  evalCase["Schur-F (CLN3+1)", YReducedSchurF[QCLN3, #, dSchur] &]};

Print[""];
Print["=== A-axis (wall 1e5..1e7) + C-axis (asymptote ratio r(f)) ==="];
Print[StringPadRight["case", 22] <> StringPadRight["basis", 7] <>
      StringPadRight["wall %", 11] <> StringPadRight["r(1e8)", 11] <>
      StringPadRight["r(1e10)", 11] <> StringPadRight["r(1e12)", 11]];
basisSizes = {3, 6, 10, 4, 7, 10, 4};
Do[Module[{r = axisAC[[k]]},
    Print[StringPadRight[r["name"], 22],
          StringPadRight[ToString[basisSizes[[k]]], 7],
          StringPadRight[ToString[NumberForm[100 r["wall_peak"], {6, 2}]] <> "%", 11],
          StringPadRight[ToString[NumberForm[r["r_1e8"], {6, 4}]], 11],
          StringPadRight[ToString[NumberForm[r["r_1e10"], {6, 4}]], 11],
          StringPadRight[ToString[NumberForm[r["r_1e12"], {6, 4}]], 11]]],
   {k, Length[axisAC]}];

(* === B-axis: Taylor coefficients at s=0 vs exact ======================= *)
(* a_n = (-1)^n sum_k g_k^2 tau_k^n (Foster expansion at s=0)             *)
(* For reduced basis: a_n^R = b_r^T (K_0,r^-1 sigma_r)^n K_0,r^-1 b_r *)

taylorExact[n_] := (-1)^n Total[g2Vec * tauVec^n];

taylorReduced[Qlist_, n_] := Module[{K0r, sigmar, br, K0rInv, M, v},
  K0r   = Qlist . K0matrix . Transpose[Qlist];
  sigmar = Qlist . sigmaM . Transpose[Qlist];
  br    = Qlist . bvec;
  K0rInv = Inverse[K0r];
  M     = (-1)^n MatrixPower[K0rInv . sigmar, n] . K0rInv;
  sigma (br . M . br)];

Print[""];
Print["=== B-axis: Taylor coefficients of Y_R(s) at s=0 ==="];
Print["Order | exact a_n         | CLN_N3 rel-err | CLN_N6 rel-err | Multi-K M3+4pts rel-err | Multi-K M3+7pts rel-err | Schur-F rel-err"];
schurTaylor[n_] := taylorReduced[QCLN3, n] + If[n == 0, 0, 0];
(* Schur-F adds K_SIBC sqrt(s)/(s+d): at s=0 this is 0, derivatives non-analytic *)
(* For B-axis comparison, Schur-F Taylor = CLN_N3 Taylor (sqrt(s) gives no integer Taylor at s=0) *)

bTable = Table[
  Module[{aEx, aN3, aN6, aMK7, aMK10, aSchur, relN3, relN6, relMK7, relMK10, relSchur},
    aEx     = taylorExact[n];
    aN3     = taylorReduced[QCLN3, n];
    aN6     = taylorReduced[QCLN6, n];
    aMK7    = taylorReduced[QMK7, n];
    aMK10   = taylorReduced[QMK10, n];
    aSchur  = aN3;  (* Schur-F sqrt(s) block contributes 0 at all integer Taylor orders at s=0 *)
    relN3    = If[Abs[aEx] > 0, Abs[(aN3   - aEx)/aEx], 0];
    relN6    = If[Abs[aEx] > 0, Abs[(aN6   - aEx)/aEx], 0];
    relMK7   = If[Abs[aEx] > 0, Abs[(aMK7  - aEx)/aEx], 0];
    relMK10  = If[Abs[aEx] > 0, Abs[(aMK10 - aEx)/aEx], 0];
    relSchur = If[Abs[aEx] > 0, Abs[(aSchur - aEx)/aEx], 0];
    <|"n" -> n, "exact" -> N[aEx, 6],
      "rel_CLN_N3" -> N[relN3], "rel_CLN_N6" -> N[relN6],
      "rel_MultiK_M3p4" -> N[relMK7], "rel_MultiK_M3p7" -> N[relMK10],
      "rel_Schur_F" -> N[relSchur]|>],
  {n, 0, 12}];

Do[Module[{r = bTable[[k+1]]},
    Print[StringPadLeft[ToString[r["n"]], 5], " | ",
          StringPadRight[ToString[ScientificForm[r["exact"], 4]], 18], " | ",
          StringPadRight[ToString[ScientificForm[r["rel_CLN_N3"], 3]], 14], " | ",
          StringPadRight[ToString[ScientificForm[r["rel_CLN_N6"], 3]], 14], " | ",
          StringPadRight[ToString[ScientificForm[r["rel_MultiK_M3p4"], 3]], 22], " | ",
          StringPadRight[ToString[ScientificForm[r["rel_MultiK_M3p7"], 3]], 22], " | ",
          StringPadRight[ToString[ScientificForm[r["rel_Schur_F"], 3]], 12]]],
   {k, 0, 12}];

(* === Full frequency sweep for figure data ============================= *)
Print[""];
Print["=== Full frequency sweep (for figure) ==="];
sweepResults = Association[];
modelMap = <|
  "CLN_N3"           -> (YReducedFromBasis[QCLN3, #] &),
  "CLN_N10"          -> (YReducedFromBasis[QCLN10, #] &),
  "Multi_K_M3p7"     -> (YReducedFromBasis[QMK10, #] &),
  "Schur_F_CLN3p1"   -> (YReducedSchurF[QCLN3, #, dSchur] &)|>;

Do[Module[{name = First[entry], modelFn = Last[entry], data},
    data = Table[
      Module[{om = 2 Pi f, yEx, yR, relErr},
        yEx    = YExactFoster[om];
        yR     = modelFn[om];
        relErr = N[Abs[(yR - yEx)/yEx]];
        <|"f_Hz" -> N[f], "Y_abs" -> N[Abs[yR]],
          "Y_re" -> N[Re[yR]], "Y_im" -> N[Im[yR]],
          "rel_err" -> relErr,
          "r_asymp" -> N[Abs[yR] Sqrt[om] / KSIBC3D]|>], {f, freqList}];
    sweepResults[name] = data], {entry, Normal[modelMap]}];

exactSweep = Table[
  Module[{om = 2 Pi f, yEx},
    yEx = YExactFoster[om];
    <|"f_Hz" -> N[f], "Y_abs" -> N[Abs[yEx]],
      "Y_re" -> N[Re[yEx]], "Y_im" -> N[Im[yEx]],
      "r_asymp" -> N[Abs[yEx] Sqrt[om] / KSIBC3D]|>], {f, freqList}];

(* === Save JSON ========================================================= *)
out = <|
  "params" -> <|
    "a_m" -> N[a], "b_m" -> N[b], "c_m" -> N[c],
    "surface_area_m2" -> N[Sarea],
    "sigma_S_per_m" -> N[sigma], "mu_H_per_m" -> N[mu],
    "K_SIBC_3D" -> N[KSIBC3D],
    "Nmodes" -> Nmodes, "Mmax" -> Mmax,
    "d_Schur_rad_s" -> N[dSchur]|>,
  "A_C_axes" -> axisAC,
  "B_axis_Taylor" -> bTable,
  "sweep" -> sweepResults,
  "exact_sweep" -> exactSweep|>;

Export["cuboid3D_schur_F.json", out, "JSON"];

(* Long-form CSV for plotting *)
Module[{fp},
  fp = OpenWrite["cuboid3D_schur_F_sweep.csv"];
  WriteString[fp, "case,f_Hz,Y_abs,Y_re,Y_im,rel_err,r_asymp\n"];
  Do[Module[{name = nm},
    Do[Module[{r = sweepResults[name][[i]]},
      WriteString[fp,
        name, ",",
        ToString[CForm[N[r["f_Hz"]]]], ",",
        ToString[CForm[N[r["Y_abs"]]]], ",",
        ToString[CForm[N[r["Y_re"]]]], ",",
        ToString[CForm[N[r["Y_im"]]]], ",",
        ToString[CForm[N[r["rel_err"]]]], ",",
        ToString[CForm[N[r["r_asymp"]]]], "\n"]
       ], {i, 1, Length[sweepResults[name]]}]],
    {nm, Keys[sweepResults]}];
  Close[fp];
  fp = OpenWrite["cuboid3D_schur_F_exact.csv"];
  WriteString[fp, "f_Hz,Y_abs,Y_re,Y_im,r_asymp\n"];
  Do[Module[{r = exactSweep[[i]]},
    WriteString[fp,
      ToString[CForm[N[r["f_Hz"]]]], ",",
      ToString[CForm[N[r["Y_abs"]]]], ",",
      ToString[CForm[N[r["Y_re"]]]], ",",
      ToString[CForm[N[r["Y_im"]]]], ",",
      ToString[CForm[N[r["r_asymp"]]]], "\n"]
     ], {i, 1, Length[exactSweep]}];
  Close[fp]];

Print[""];
Print["Wrote cuboid3D_schur_F.json, cuboid3D_schur_F_sweep.csv, cuboid3D_schur_F_exact.csv"];
Print[""];
Print["==================================================="];
Print[" Summary"];
Print["==================================================="];
Print[" 3D cuboid scalar diffusion problem, ", N[a 1000], " x ", N[b 1000], " x ", N[c 1000], " mm Cu"];
Print[" K_SIBC^{3D} = 2(ab+bc+ca) sqrt(sigma/mu) = ", N[KSIBC3D, 8]];
Print[" 3-axis ROM evaluation: see A_C_axes (wall band + r(f)) and B_axis_Taylor in JSON"];
Print[" Expected: Multi-K wins A; CLN_N10 == Multi-K M3+7pt on B; Schur-F uniquely keeps r(1e12)~1 on C"];
Print["==================================================="];
