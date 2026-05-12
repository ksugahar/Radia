(* BEM-Foster Phase 3: full Foster eigendecomposition.

   Strategy:
   - Reuse cell-cell Newton matrix K from Phase 2.
   - Solve generalized eigenvalue problem K phi = mu V_cell phi.
     => Foster eigenvalues lambda_k = mu_k / (4 pi),
        time constants tau_k = mu0 sigma lambda_k.
   - M-orthonormalize eigenvectors: <phi_j, phi_k>_M = delta_jk.
   - Project A_ext onto modes:
       g_k^x = <phi_k, A_ext^x>_M = sum_i phi_{k,i} (-B_0 y_i/2) V_cell
       g_k^y = <phi_k, A_ext^y>_M = sum_i phi_{k,i} (+B_0 x_i/2) V_cell
   - Foster sum (B_0 = 1):
       alpha_n = (-1)^n * (mu0 sigma)^n / V *
                 sum_k (g_k^x^2 + g_k^y^2) * lambda_k^(n-1)
   - Compare alpha_1 (closed), alpha_2 (closed), alpha_3 (Born MC reference).
   - Inspect dominant tau_k and compare with closed-form tau_1 = 16.98 us.

   Date: 2026-05-01
*)

$HistoryLength = 0;
prec = 25;

Print["================================================================"];
Print["BEM-Foster Phase 3: full Foster eigendecomposition"];
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

(* === References === *)
alpha1Ref = -mu0sigmaN (aN^2 + bN^2)/48;       (* closed form *)
alpha2Ref = SetPrecision[7.47832405353*^-10, 16]; (* closed form *)
alpha3RefMC = -1.435*^-14;                       (* 3-digit MC, progress.tex *)
tau1Ref = SetPrecision[16.98281*^-6, 12];        (* alpha_2 / (-alpha_1) *)

Print["Material : Cu, sigma = ", N[sigmaN, 6], " S/m"];
Print["Geometry : 5 x 2 x 1 mm  (V = ", N[VN, 6], " m^3)"];
Print[""];
Print["References:"];
Print["  alpha_1 (closed)    = ", N[alpha1Ref, 14]];
Print["  alpha_2 (closed)    = ", alpha2Ref];
Print["  alpha_3 (3-dig MC)  = ", alpha3RefMC];
Print["  tau_1   (closed)    = ", N[tau1Ref, 12], " s = 16.98 us"];
Print[""];

(* === Cell-cell Newton kernel === *)
cellNewtonOff[dx_, dy_, dz_, hx_, hy_, hz_] :=
  NIntegrate[
    (hx - Abs[u])(hy - Abs[v])(hz - Abs[w])/
      Sqrt[(u + dx)^2 + (v + dy)^2 + (w + dz)^2],
    {u, -hx, hx}, {v, -hy, hy}, {w, -hz, hz},
    WorkingPrecision -> 20, AccuracyGoal -> 14, PrecisionGoal -> 12,
    MaxRecursion -> 14, Method -> "GlobalAdaptive"];

cellNewtonSelf[hx_, hy_, hz_] :=
  8 NIntegrate[
    (hx - u)(hy - v)(hz - w)/Sqrt[u^2 + v^2 + w^2],
    {u, 0, hx}, {v, 0, hy}, {w, 0, hz},
    WorkingPrecision -> 20, AccuracyGoal -> 14, PrecisionGoal -> 12,
    MaxRecursion -> 16, Method -> "GlobalAdaptive"];

cellNewton[dx_, dy_, dz_, hx_, hy_, hz_] :=
  If[dx == 0 && dy == 0 && dz == 0,
    cellNewtonSelf[hx, hy, hz],
    cellNewtonOff[dx, dy, dz, hx, hy, hz]];

(* === Build K matrix and run Foster decomposition === *)
fosterSolve[Nx_, Ny_, Nz_] := Module[
  {hx, hy, hz, Vcell, ncells, centers, K, t0,
   cache, key, val, di,
   eigvals, eigvecs, phis, lambdaK, taus,
   xs, ys, gxlist, gylist, alphas},
  hx = aN/Nx; hy = bN/Ny; hz = cN/Nz;
  Vcell = hx hy hz;
  ncells = Nx Ny Nz;

  centers = Flatten[
    Table[
      {(i - 1/2)hx - aN/2, (j - 1/2)hy - bN/2, (k - 1/2)hz - cN/2},
      {i, Nx}, {j, Ny}, {k, Nz}],
    2];

  Print["  Grid (", Nx, ", ", Ny, ", ", Nz, "): ",
        ncells, " cells, V_cell = ", N[Vcell, 6]];

  (* Build K with translation-invariance cache *)
  Print["  Building K..."];
  t0 = AbsoluteTime[];
  cache = <||>;
  K = ConstantArray[0, {ncells, ncells}];
  Do[
    di = centers[[i]] - centers[[j]];
    key = Abs[Round[di/{hx, hy, hz}]];
    val = Lookup[cache, Key[key], $Failed];
    If[val === $Failed,
      val = cellNewton[di[[1]], di[[2]], di[[3]], hx, hy, hz];
      cache[key] = val];
    K[[i, j]] = val;
    K[[j, i]] = val,
    {i, ncells}, {j, i, ncells}];
  Print["    K assembly: ", N[AbsoluteTime[] - t0, 4], " s, ",
        Length[cache], " unique values"];

  (* Eigendecomposition: K c = mu V_cell c *)
  Print["  Eigendecomposition (", ncells, " modes)..."];
  t0 = AbsoluteTime[];
  {eigvals, eigvecs} = Eigensystem[N[K/Vcell, prec]];
  Print["    done: ", N[AbsoluteTime[] - t0, 4], " s"];

  (* M-orthonormalize: <phi_k, phi_k>_M = sum_i phi_{k,i}^2 V_cell = 1
     => phi_k = (raw eigvec) / Sqrt[V_cell * Norm[raw]^2] *)
  phis = (#/Sqrt[Vcell #.#]) & /@ eigvecs;

  (* True Newton eigenvalues lambda_k = mu_k / (4 pi)
     (since K_ij = 4 pi * Newton kernel matrix element) *)
  lambdaK = eigvals/(4 PiN);
  taus = mu0sigmaN lambdaK;

  (* Project A_ext (B_0 = 1) *)
  xs = centers[[All, 1]];
  ys = centers[[All, 2]];
  gxlist = (phis.(-ys/2)) Vcell;
  gylist = (phis.(+xs/2)) Vcell;

  (* Foster sum (B_0 = 1, so V*B_0^2 = V) *)
  alphas = Table[
    (-1)^n mu0sigmaN^n / VN Total[(gxlist^2 + gylist^2) lambdaK^(n - 1)],
    {n, 1, 5}];

  <|"taus" -> taus, "alphas" -> alphas,
    "gx" -> gxlist, "gy" -> gylist,
    "lambdaK" -> lambdaK, "K" -> K, "centers" -> centers|>];

(* === Run on (8, 4, 2) grid === *)
Print["=== Solving (8, 4, 2) grid ==="];
res = fosterSolve[8, 4, 2];
Print[""];

(* === Foster spectrum (top 10 modes by tau) === *)
Print["=== Foster spectrum: top 10 modes by tau_k ==="];
Print["  k   tau_k [us]            g_k^x           g_k^y           a_k"];

ord = Reverse[Ordering[res["taus"]]];
Do[
  Module[{ki = ord[[k]], tk, gxk, gyk, ak},
    tk = res["taus"][[ki]];
    gxk = res["gx"][[ki]];
    gyk = res["gy"][[ki]];
    ak = (gxk^2 + gyk^2)/(VN res["lambdaK"][[ki]]);
    Print["  ", PaddedForm[k, 2], "  ",
          ScientificForm[N[tk 10^6, 8], 7], "  ",
          ScientificForm[N[gxk, 6], 5], "  ",
          ScientificForm[N[gyk, 6], 5], "  ",
          ScientificForm[N[ak, 6], 5]]],
  {k, Min[10, Length[res["taus"]]]}];
Print[""];

(* === alpha_n via Foster sum === *)
Print["=== alpha_n via Foster sum (B_0 = 1) ==="];
refs = {alpha1Ref, alpha2Ref, alpha3RefMC, "?", "?"};
Do[
  Module[{a, ref, err},
    a = res["alphas"][[n]];
    ref = refs[[n]];
    Print["  alpha_", n, " = ", N[a, 12]];
    If[NumberQ[ref],
      err = Abs[(a - ref)/ref];
      Print["    ref     = ", N[ref, 12], "   rel.err = ", N[err, 4]]]],
  {n, 5}];
Print[""];

(* === Comparison summary === *)
Print["================================================================"];
Print["Summary: BEM-Foster on (8,4,2) grid"];
Print["================================================================"];
Print["* tau_1 (BEM dominant)  = ", N[res["taus"][[ord[[1]]]] 10^6, 8], " us"];
Print["* tau_1 (closed form)   = 16.98281 us"];
Print["* alpha_1 BEM/closed    = ",
      N[res["alphas"][[1]]/alpha1Ref, 8]];
Print["* alpha_2 BEM/closed    = ",
      N[res["alphas"][[2]]/alpha2Ref, 8]];
Print["* alpha_3 BEM/MC-ref    = ",
      N[res["alphas"][[3]]/alpha3RefMC, 6]];
Print[""];

(* Tau spectrum range *)
taus = res["taus"];
Print["* tau_k range: ", N[Max[taus] 10^6, 6], " us  ...  ",
      N[Min[taus] 10^9, 6], " ns"];
Print["* Number of modes: ", Length[taus]];
Print[""];

(* Save Foster pairs *)
Export["bem_foster_phase3_results.json",
  <|"grid" -> {8, 4, 2},
    "taus" -> N[res["taus"], 16],
    "gx" -> N[res["gx"], 16],
    "gy" -> N[res["gy"], 16],
    "alphas" -> N[res["alphas"], 16],
    "alpha1Closed" -> N[alpha1Ref, 16],
    "alpha2Closed" -> N[alpha2Ref, 16],
    "alpha3MC" -> alpha3RefMC|>];
Print["* Results saved to bem_foster_phase3_results.json"];
