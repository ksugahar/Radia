(* BEM-Foster Phase 2: alpha_2 validation on Cartesian hex grid.

   Strategy:
   - Discretize cuboid Omega = [-a/2,a/2] x [-b/2,b/2] x [-c/2,c/2]
     into Nx x Ny x Nz uniform hex cells with piecewise-constant basis.
   - Cell-cell Newton matrix K_{ij} = int_{Omega_i} int_{Omega_j} 1/|r-r'| dVdV'
     reduced to 3D via change of variables u = r - r':
       K_ij(d) = int_{|u_a|<h_a} prod_a(h_a - |u_a|) / |u + d| du,  d = r_i - r_j
   - alpha_2^BEM = (mu0 sigma)^2 / (16 pi V) * sum_{ij} (x_i x_j + y_i y_j) K_{ij}
     (midpoint quadrature; expected at least O(h) convergence, possibly O(h^2)
     by symmetric-error cancellation).
   - Compare to closed form alpha_2 = 7.4783e-10
     (cuboid_3D_vacuum_Born_alpha2_closed.wls, 30+ digit Aharoni 6D).

   Date: 2026-05-01
*)

$HistoryLength = 0;
prec = 25;

Print["================================================================"];
Print["BEM-Foster Phase 2: alpha_2 validation (Cartesian hex grid)"];
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

Print["Material : Cu, sigma = ", N[sigmaN, 6], " S/m"];
Print["Geometry : ", N[aN 1000, 4], " x ", N[bN 1000, 4],
      " x ", N[cN 1000, 4], " mm  (V = ", N[VN, 6], " m^3)"];
Print[""];

alpha2Closed = SetPrecision[7.47832405353*^-10, 16];
Print["alpha_2 (closed form, Born/Aharoni 6D): ", alpha2Closed];
Print[""];

(* === Cell-cell Newton kernel === *)
(* Off-diagonal (d != 0): smooth integrand. *)
cellNewtonOff[dx_, dy_, dz_, hx_, hy_, hz_] :=
  NIntegrate[
    (hx - Abs[u])(hy - Abs[v])(hz - Abs[w])/
      Sqrt[(u + dx)^2 + (v + dy)^2 + (w + dz)^2],
    {u, -hx, hx}, {v, -hy, hy}, {w, -hz, hz},
    WorkingPrecision -> 20, AccuracyGoal -> 14, PrecisionGoal -> 12,
    MaxRecursion -> 14, Method -> "GlobalAdaptive"];

(* Self-cell (d = 0): exploit octant symmetry, 1/r corner singularity. *)
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

(* === Sanity check: unit cube self-energy === *)
Print["=== Sanity check ==="];
t0 = AbsoluteTime[];
KselfUnit = cellNewtonSelf[1, 1, 1];
Print["  K_self(unit cube) = ", N[KselfUnit, 14],
      "  (", N[AbsoluteTime[] - t0, 3], " s)"];
Print["  Expected: ~1.88231264338  (literature, e.g. Hummer 1996)"];
Print[""];

(* === Build BEM matrix and compute alpha_2 === *)
buildAndCompute[Nx_, Ny_, Nz_] := Module[
  {hx, hy, hz, ncells, centers, K, di, t0, xs, ys, Dxx, Dyy,
   alpha2BEM, err, cache, key, val, dt},
  hx = aN/Nx; hy = bN/Ny; hz = cN/Nz;
  ncells = Nx Ny Nz;

  centers = Flatten[
    Table[
      {(i - 1/2)hx - aN/2, (j - 1/2)hy - bN/2, (k - 1/2)hz - cN/2},
      {i, Nx}, {j, Ny}, {k, Nz}],
    2];

  Print["  Grid (", Nx, ", ", Ny, ", ", Nz, "): ",
        ncells, " cells, ", ncells (ncells + 1)/2, " pairs"];
  t0 = AbsoluteTime[];

  (* Translation-invariance cache: K_ij depends on (|delta_i|, |delta_j|, |delta_k|)
     since the kernel integrand is even in each component of d. *)
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
  dt = AbsoluteTime[] - t0;
  Print["    K assembly: ", N[dt, 4], " s, ",
        Length[cache], " unique kernel values"];

  xs = centers[[All, 1]];
  ys = centers[[All, 2]];
  Dxx = xs.K.xs;
  Dyy = ys.K.ys;
  alpha2BEM = (mu0N sigmaN)^2/(16 PiN VN) (Dxx + Dyy);
  err = Abs[(alpha2BEM - alpha2Closed)/alpha2Closed];
  {alpha2BEM, err, dt}];

(* === Convergence study === *)
Print["=== Convergence study ==="];
Print[""];
Print["Grid                alpha_2^BEM            rel. error"];
Print["----------------------------------------------------------"];

results = Reap[
  Do[
    Module[{a2, err, dt},
      {a2, err, dt} = buildAndCompute @@ config;
      Print["  ", config, "    ", N[a2, 12], "   ", N[err, 4],
            "   (", N[dt, 4], " s)"];
      Sow[<|"config" -> config, "alpha2" -> N[a2, 16],
            "relErr" -> N[err, 8], "time" -> N[dt, 6]|>]],
    {config, {{2, 2, 2}, {4, 2, 2}, {4, 4, 2}, {6, 4, 2}, {8, 4, 2}}}]
][[2, 1]];

Print[""];
Print["================================================================"];
Print["Summary"];
Print["================================================================"];
Print["* alpha_2 (closed form): ", alpha2Closed];
Print["* alpha_2 (BEM, finest): ", results[[-1, "alpha2"]]];
Print["* relative error       : ", results[[-1, "relErr"]]];
Print[""];

If[Length[results] >= 2,
  Module[{r1, r2, h1, h2, p},
    r1 = results[[-2]]; r2 = results[[-1]];
    h1 = (VN/Times @@ r1["config"])^(1/3);
    h2 = (VN/Times @@ r2["config"])^(1/3);
    p = Log[r1["relErr"]/r2["relErr"]]/Log[h1/h2];
    Print["* Estimated convergence rate (err ~ h^p): p = ", N[p, 4]]]];

Export["bem_foster_alpha2_results.json",
  <|"alpha2Closed" -> N[alpha2Closed, 16],
    "results" -> results|>];
Print["* Results saved to bem_foster_alpha2_results.json"];
