#!/usr/bin/env wolframscript
(* ::Package:: *)

(*
  cuboid3D_AAA_CFschur.wls

  AAA-flavoured rational approximation of the Schur block
    z_block(s) = K_SIBC * sqrt(s) / (s + d)
  using a DC-preserving pole-residue ansatz:
    g_N(s) = K_SIBC * sum_{k=1}^N r_k * (1/p_k - 1/(s + p_k))
           = K_SIBC * sum_{k=1}^N r_k * s / (p_k * (s + p_k))

  At s = 0: g_N(0) = 0  (DC PRESERVED by construction)
  At s -> oo: g_N -> K_SIBC * (sum r_k/p_k) -> finite (BANDLIMITED, not 1/sqrt(s))

  This solves the B-axis contamination of the naive Stieltjes CF
  (which had non-zero DC) while accepting that no finite-N rational
  function can reach 1/sqrt(s) asymptotically -- the BANDLIMIT.

  The poles p_k are chosen log-spaced over [p_min, p_max] and the
  residues r_k solve a complex least-squares problem matching
  z_target(jw) at K log-spaced frequency samples on [w_min, w_max].

  Comparison against:
    - Naive Stieltjes CF (B-axis contaminated, A/C axes failed)
    - Non-rational sqrt(s) Schur-F (gold standard, exact at all s)

  Run: wolframscript -file cuboid3D_AAA_CFschur.wls
*)

ClearAll[sigma, mu, nu];
sigma = 5.8*^7;
mu    = 4 Pi 1*^-7;
nu    = 1/mu;

geometries = {
  <|"name" -> "A1 (17.72x17.72x2 mm Cu)",
    "a" -> 17.72*^-3, "b" -> 17.72*^-3, "c" -> 2*^-3,
    "wall_lo" -> 1*^3, "wall_hi" -> 1*^5,
    "d_schur" -> 2 Pi 1*^6,
    "p_min" -> 2 Pi 1*^4, "p_max" -> 2 Pi 1*^12|>,
  <|"name" -> "5x2x1 mm Cu",
    "a" -> 5*^-3, "b" -> 2*^-3, "c" -> 1*^-3,
    "wall_lo" -> 1*^5, "wall_hi" -> 1*^7,
    "d_schur" -> 2 Pi 1*^9,
    "p_min" -> 2 Pi 1*^6, "p_max" -> 2 Pi 1*^14|>
};

(* === DC-preserving rational ansatz with LS residues =================== *)
(* g(s) = K_SIBC * sum_{k=1..N} r_k * s / (p_k * (s + p_k))                *)
(* Target: z_target(s) = K_SIBC * sqrt(s) / (s + d)                       *)
(* LS over N_sample log-spaced frequencies in [w_lo, w_hi]                *)

fitDCPreservingPoles[Nrank_, pMin_, pMax_, dSchur_, KSIBC_, wLo_, wHi_, nSample_:60] :=
  Module[{poles, samples, omegas, sValues, Amat, bvecLS, rVec, basisFn},
    poles  = N[pMin Exp[Log[pMax/pMin] Subdivide[0, 1, Nrank - 1]]];
    omegas = N[Exp[Log[wLo] + (Log[wHi] - Log[wLo]) Subdivide[0, 1, nSample - 1]]];
    sValues = I omegas;
    basisFn[s_, p_] := s / (p (s + p));
    Amat = Table[
      Table[basisFn[sValues[[i]], poles[[k]]], {k, 1, Nrank}],
      {i, 1, nSample}];
    bvecLS = Table[KSIBC Sqrt[sValues[[i]]] / (sValues[[i]] + dSchur), {i, 1, nSample}];
    (* LS: solve A . r = b in least-squares sense *)
    rVec = LeastSquares[Amat, bvecLS];
    rVec = Re[rVec];  (* enforce real residues (passive) *)
    {poles, rVec}];  (* K_SIBC is already absorbed into bvecLS target *)

(* AAA-fit evaluator: g_AAA(s) = sum_k r_k * s / (p_k * (s + p_k))      *)
gAAA[s_, poles_, rEff_] :=
  Total[rEff (s / (poles (s + poles)))];

(* === Per-geometry analysis ============================================= *)
analyzeGeom[geom_] := Module[
  {a, b, c, name, wallLo, wallHi, dSchur, pMin, pMax, Sarea, KSIBC,
   Mmax = 13, oddRange, modeList, Nmodes, lamVec, tauVec, g2Vec,
   LambdaM, sigmaM, Imat, Kmodal, K0matrix, bvec, YExactFoster,
   K0Inner, K0MGS, K0Normalise, buildCLNbasis,
   YReducedFromBasis, YSchurNR, YAAA, QCLN3, freqWall, freqProbe,
   results, evalAAA, evalNR, NCFs},
  a = geom["a"]; b = geom["b"]; c = geom["c"];
  name = geom["name"]; wallLo = geom["wall_lo"]; wallHi = geom["wall_hi"];
  dSchur = geom["d_schur"]; pMin = geom["p_min"]; pMax = geom["p_max"];
  Sarea = 2 (a b + b c + c a);
  KSIBC = Sarea Sqrt[sigma/mu];

  Print["==============================================="];
  Print[" Geometry: ", name];
  Print["==============================================="];
  Print["  K_SIBC = ", N[KSIBC, 6]];
  Print["  Schur-F d/(2pi) = ", N[dSchur/(2 Pi), 4], " Hz"];
  Print["  AAA pole grid [p_min, p_max]/(2pi) = [",
        N[pMin/(2 Pi), 3], ", ", N[pMax/(2 Pi), 3], "] Hz"];

  (* --- Foster spectrum (343 modes) --- *)
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

  YReducedFromBasis[Qlist_, omega_] := Module[{Kr, Qb},
    Kr = Qlist . Kmodal[I omega] . Transpose[Qlist];
    Qb = Qlist . bvec;
    sigma (Qb . LinearSolve[Kr, Qb])];

  YSchurNR[Qlist_, omega_] :=
    YReducedFromBasis[Qlist, omega] + KSIBC Sqrt[I omega]/(I omega + dSchur);

  QCLN3 = buildCLNbasis[0, 3];

  freqWall = Select[N[10^Range[Log10[wallLo], Log10[wallHi], 0.25]],
                    wallLo <= # <= wallHi &];
  freqProbe = N[{1*^6, 1*^8, 1*^10, 1*^12}];

  NCFs = {4, 8, 16};

  Print[""];
  Print["  --- AAA CF-Schur (DC-preserving rational) ---"];
  Print["  ", StringPadLeft["N", 4], " ",
        StringPadRight["wall %", 11], StringPadRight["DC offset", 14],
        StringPadRight["r(1e6)", 11], StringPadRight["r(1e8)", 11],
        StringPadRight["r(1e10)", 11], StringPadRight["r(1e12)", 11]];
  results = Table[
    Module[{poles, rEff, YAAAfn, wallErr, dcOffset, rVec, yEx, yR},
      {poles, rEff} = fitDCPreservingPoles[NC, pMin, pMax, dSchur, KSIBC,
                        2 Pi wallLo, 2 Pi 1*^11, 60];
      YAAAfn[om_] := YReducedFromBasis[QCLN3, om] + gAAA[I om, poles, rEff];
      wallErr = Max[Table[
        Module[{om = 2 Pi f},
          yEx = YExactFoster[om];
          yR  = YAAAfn[om];
          Abs[(yR - yEx)/yEx]], {f, freqWall}]];
      dcOffset = gAAA[0., poles, rEff];   (* should be 0 by construction *)
      rVec = Table[
        Module[{om = 2 Pi f},
          yR = YAAAfn[om];
          Abs[yR] Sqrt[om] / KSIBC], {f, freqProbe}];
      Print["  ", StringPadLeft[ToString[NC], 4], " ",
            StringPadRight[ToString[NumberForm[100 wallErr, {6, 2}]] <> "%", 11],
            StringPadRight[ToString[ScientificForm[N[dcOffset], 3]], 14],
            StringPadRight[ToString[NumberForm[rVec[[1]], {6, 4}]], 11],
            StringPadRight[ToString[NumberForm[rVec[[2]], {6, 4}]], 11],
            StringPadRight[ToString[NumberForm[rVec[[3]], {6, 4}]], 11],
            StringPadRight[ToString[NumberForm[rVec[[4]], {6, 4}]], 11]];
      <|"N" -> NC, "poles" -> N[poles], "rEff" -> N[rEff],
        "wall_peak" -> N[wallErr],
        "dc_offset" -> N[dcOffset],
        "r_1e6"  -> N[rVec[[1]]], "r_1e8"  -> N[rVec[[2]]],
        "r_1e10" -> N[rVec[[3]]], "r_1e12" -> N[rVec[[4]]]|>],
    {NC, NCFs}];

  (* Non-rational reference *)
  Module[{wallErr, rVec, yR},
    wallErr = Max[Table[
      Module[{om = 2 Pi f, yEx},
        yEx = YExactFoster[om];
        yR  = YSchurNR[QCLN3, om];
        Abs[(yR - yEx)/yEx]], {f, freqWall}]];
    rVec = Table[
      Module[{om = 2 Pi f, yR2},
        yR2 = YSchurNR[QCLN3, om];
        Abs[yR2] Sqrt[om] / KSIBC], {f, freqProbe}];
    Print["  ", StringPadLeft["NR", 4], " ",
          StringPadRight[ToString[NumberForm[100 wallErr, {6, 2}]] <> "%", 11],
          StringPadRight["(0 exact)", 14],
          StringPadRight[ToString[NumberForm[rVec[[1]], {6, 4}]], 11],
          StringPadRight[ToString[NumberForm[rVec[[2]], {6, 4}]], 11],
          StringPadRight[ToString[NumberForm[rVec[[3]], {6, 4}]], 11],
          StringPadRight[ToString[NumberForm[rVec[[4]], {6, 4}]], 11]];
    Print["  (NR = non-rational sqrt(s)/(s+d) Schur-F at d/(2pi)=",
          N[dSchur/(2 Pi), 4], " Hz)"]];

  Print[""];
  <|"geom" -> geom, "K_SIBC" -> N[KSIBC],
    "AAA_results" -> results|>];

(* === Run both geometries =============================================== *)
allResults = Map[analyzeGeom, geometries];

(* === Save JSON (strip symbolic Pi etc.) ================================ *)
toReal[x_] := If[NumericQ[x], N[x],
                 If[ListQ[x], Map[toReal, x],
                    If[AssociationQ[x], Association[Map[# -> toReal[x[#]] &, Keys[x]]],
                       x]]];
out = <|
  "geometries" -> toReal[Table[allResults[[k]]["geom"], {k, Length[allResults]}]],
  "AAA_results" -> Table[
    <|"name"     -> allResults[[k]]["geom"]["name"],
      "K_SIBC"   -> N[allResults[[k]]["K_SIBC"]],
      "results"  -> toReal[allResults[[k]]["AAA_results"]]|>,
    {k, Length[allResults]}]|>;
Export["cuboid3D_AAA_CFschur.json", out, "JSON"];
Print[""];
Print["Saved cuboid3D_AAA_CFschur.json"];

(* === Headline summary ================================================= *)
Print[""];
Print["============================================================"];
Print[" Headline: AAA-CF-Schur vs naive Stieltjes vs non-rational"];
Print[" (compare with naive CF-Schur results in cuboid3D_CFschur_both.wls)"];
Print[" Key improvements over naive CF:"];
Print["   - DC offset = 0 by ansatz (B-axis preserved)"];
Print["   - Wall-band accuracy controlled by LS fit"];
Print["   - C-axis: still rational, still bandlimited (cannot reach 1/sqrt)"];
Print["   - For SPICE realisation: passive RLC ladder with real poles"];
Print["============================================================"];
