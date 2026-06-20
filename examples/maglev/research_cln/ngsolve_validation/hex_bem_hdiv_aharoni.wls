(* ::Package:: *)

(* hex_bem_hdiv_aharoni.wls
   HDiv div-free BEM with Aharoni 1D-closed inner reduction:
   K[i,j] entry = sum over component pairs (xx, yy, zz) of:
     - Polynomial expansion of phi_i_alpha * phi_j_alpha
     - Aharoni-extended closed form for one direction (say x for xx component)
     - 4D NIntegrate for remaining (y, y', z, z')

   Total speedup: ~10-100x vs full 6D NIntegrate.

   Date: 2026-05-03
*)

ClearAll[x, y, z, xp, yp, zp, u, up, len, rho];
$MaxExtraPrecision = 60;

Print[StringRepeat["=", 70]];
Print["HDiv div-free BEM with Aharoni 1D-closed inner (fast)"];
Print["Started: ", DateString[]];
Print[StringRepeat["=", 70]];

sigma = SetPrecision[58 10^6, 30];
mu0   = 4 Pi SetPrecision[10^-7, 30];
{aBox, bBox, cBox} = SetPrecision[{5 10^-3, 2 10^-3, 1 10^-3}, 30];

(* === Helpers === *)
DuCrossDv[u_, v_, vars_] := Cross[D[u, {vars}], D[v, {vars}]];
wDuCrossDv[u_, v_, w_, vars_] := w Cross[D[u, {vars}], D[v, {vars}]];

IntegratedLegendre[n_, x_] := Module[{p3, p2, p1, j, result},
  result = {};
  If[n < 2, Return[{}]];
  p3 = 0; p2 = -1; p1 = x;
  Do[p3 = p2; p2 = p1;
     p1 = (2 j - 3)/j x p2 - (j - 3)/j p3;
     AppendTo[result, p1], {j, 2, n}];
  result];

LegendrePoly[n_, x_] := Module[{P, k},
  P = Table[0, {n + 1}];
  P[[1]] = 1; If[n >= 1, P[[2]] = x];
  Do[P[[k + 1]] = (2 - 1/k) x P[[k]] + (1/k - 1) P[[k - 1]], {k, 2, n}];
  P];

ShapesHDivHexDivFree[order_, x_, y_, z_] := Module[
  {p, xi, eta, zeta, LXi, LEta, LZeta, PXi, PEta, PZeta, shapes, vars},
  vars = {x, y, z};
  p = order;
  xi = 2 x - 1; eta = 2 y - 1; zeta = 2 z - 1;
  LXi = IntegratedLegendre[p + 1, xi];
  LEta = IntegratedLegendre[p + 1, eta];
  LZeta = IntegratedLegendre[p + 1, zeta];
  PXi = LegendrePoly[p, xi];
  PEta = LegendrePoly[p, eta];
  PZeta = LegendrePoly[p, zeta];
  shapes = {};
  Do[AppendTo[shapes, wDuCrossDv[LXi[[i+1]] LZeta[[k+1]], y, -4 PEta[[j+2]], vars]],
    {i, 0, p-1}, {j, 0, p-1}, {k, 0, p-1}];
  Do[AppendTo[shapes, wDuCrossDv[LEta[[j+1]] LZeta[[k+1]], x, 4 PXi[[i+2]], vars]],
    {i, 0, p-1}, {j, 0, p-1}, {k, 0, p-1}];
  Do[AppendTo[shapes, DuCrossDv[LEta[[j+1]] LZeta[[k+1]], 2 xi, vars]],
    {j, 0, p-1}, {k, 0, p-1}];
  Do[AppendTo[shapes, DuCrossDv[LXi[[i+1]] LZeta[[k+1]], 2 eta, vars]],
    {i, 0, p-1}, {k, 0, p-1}];
  Do[AppendTo[shapes, DuCrossDv[LXi[[i+1]] LEta[[j+1]], 2 zeta, vars]],
    {i, 0, p-1}, {j, 0, p-1}];
  shapes];

fastIntegrate3D[poly_, vars_] := Module[{p = Expand[poly], rules},
  If[p === 0, Return[0]];
  rules = CoefficientRules[p, vars];
  Total[Function[r, r[[2]] / Times @@ Map[(# + 1) &, r[[1]]]] /@ rules]];

fastIntegrate2D[poly_, vars_] := Module[{p = Expand[poly], rules},
  If[p === 0, Return[0]];
  rules = CoefficientRules[p, vars];
  Total[Function[r, r[[2]] / Times @@ Map[(# + 1) &, r[[1]]]] /@ rules]];

(* === Aharoni-extended F_kl(rho; len) ===
   F_kl = int_0^len int_0^len x^k x'^l / sqrt((x-x')^2 + rho^2) dx dx'
   Compute symbolically once for k, l = 0..maxK, store in lookup.
*)
maxK = 7;   (* HDiv order p has shapes with polynomial degree up to ~p+2; for order 3, up to 5 *)

Print[];
Print["Precomputing Aharoni-extended F_kl(rho; len) for k, l = 0..", maxK, " ..."];
t0 = AbsoluteTime[];

FklTable = Association[];
Do[
  Module[{Fkl, t1 = AbsoluteTime[]},
    Fkl = Integrate[u^k up^ll / Sqrt[(u - up)^2 + rho^2],
      {u, 0, len}, {up, 0, len},
      Assumptions -> {len > 0, rho > 0}];
    FklTable[{k, ll}] = Fkl;
    If[Mod[k * (maxK + 1) + ll, 5] == 0,
      Print["  F[", k, ",", ll, "] in ", Round[AbsoluteTime[] - t1, 0.1], " s"]]
  ],
  {k, 0, maxK}, {ll, 0, maxK}
];
Print["F_kl precompute done in ", Round[AbsoluteTime[] - t0, 0.1], " s"];

(* === K matrix entry via Aharoni reduction ===
   For component pair (alpha, beta), consider integrand:
     phi_i_alpha(x, y, z) * phi_j_beta(x', y', z') / kernel
   where kernel = sqrt(a^2 (x-x')^2 + b^2 (y-y')^2 + c^2 (z-z')^2).

   We reduce the ALPHA direction (e.g. x for x-x component) analytically.
   Polynomial expansion in (x, x') -> sum of x^k x'^l terms.
   Each term × F_kl(rho_yz/a; 1) × inner-direction scaling
*)

(* Apply Aharoni reduction in x-direction:
   integrand(x, x', y, y', z, z') with kernel = sqrt(a^2 dx^2 + ρ²)
   = sum_kl c_kl(y,y',z,z') x^k x'^l / sqrt(a^2 (x-x')^2 + ρ²)
   where ρ² = b² (y-y')² + c² (z-z')²

   Substitute u = x*a, x → u/a, dx → du/a, limits 0 to a:
   ... actually let me think. Reference: x ∈ [0, 1], so x' ∈ [0, 1].
   kernel sqrt(a²(x-x')² + ρ²) where x, x' in [0,1].

   ∫_0^1 ∫_0^1 x^k x'^l / sqrt(a²(x-x')² + ρ²) dx dx'
   Substitute u = a x, u' = a x': u ∈ [0, a], dx = du/a
   = (1/a²) ∫_0^a ∫_0^a (u/a)^k (u'/a)^l / sqrt((u-u')² + ρ²) du du'
   = (1/(a^(k+l+2))) F_kl(ρ; a)
*)

(* Build single K[i,j] entry *)
KentryAharoni[shi_, shj_] := Module[
  {comp, totals = 0, alpha, dirVars, otherVars, polyAlpha,
   coeffsAlpha, kpair, cKL, lenAlpha, rhoForm,
   integrand4D, val},

  (* For each component direction, integrate that direction analytically *)

  (* === x-x component: phi_i_x(x,y,z) phi_j_x(x',y',z') / kernel === *)
  Module[{phiIX = shi[[1]], phiJX = shj[[1]] /. {x -> xp, y -> yp, z -> zp}},
    If[phiIX =!= 0 && phiJX =!= 0,
      polyAlpha = Expand[phiIX phiJX];
      (* Expand polynomial in (x, xp); collect coefficients in (y, z, yp, zp) *)
      coeffsAlpha = CoefficientRules[polyAlpha, {x, xp}];
      (* Each entry: {{k, l}, c(y, z, yp, zp)} *)
      Do[
        kpair = entry[[1]];
        cKL = entry[[2]];
        (* Aharoni F_kl with len = aBox, rho^2 = bBox^2 (y-yp)^2 + cBox^2 (z-zp)^2 *)
        (* In ref [0,1] coords: F_kl_ref = (1/aBox^(k+l+2)) F_kl(rho_phys; aBox) *)
        (* But our K formula has aBox^2 prefactor for x-x component:
           K_xx = (mu0/4pi) (abc)^2 ∫∫ a^2 phi_x phi_x' / kernel dV_ref dV_ref'
                = (mu0/4pi) a^2(abc)^2 ∫∫∫∫ [c_kl(y,z,y',z') × ∫∫ x^k x'^l / kernel dx dx'] dy dy' dz dz' *)
        Module[{Fklval, integrandF, prefactor},
          Fklval = FklTable[{kpair[[1]], kpair[[2]]}];
          (* Substitute len -> aBox/aBox = 1? No.
             Our reference x range is [0,1], physical kernel uses a (x-xp).
             Let me redo: kernel = sqrt(aBox^2 (x-xp)^2 + bBox^2 (y-yp)^2 + cBox^2 (z-zp)^2)
                        with x, xp in [0, 1].
             ∫∫ x^k xp^l / sqrt(aBox^2 (x-xp)^2 + r^2) dx dxp over [0,1]^2
             Let u = aBox x: u in [0, aBox], du = aBox dx, x = u/aBox
             = (1/aBox)^2 ∫∫ (u/aBox)^k (up/aBox)^l / sqrt((u-up)^2 + r^2) du dup
             = (1/aBox^(k+l+2)) F_kl(r; aBox)
             where F_kl is what we precomputed.
             Substitute len -> aBox in F_kl, rho -> sqrt(b^2(y-yp)^2 + c^2(z-zp)^2)
          *)
          prefactor = 1 / aBox^(kpair[[1]] + kpair[[2]] + 2);
          integrandF = prefactor * (Fklval /. {len -> aBox,
                                                rho -> Sqrt[bBox^2 (y - yp)^2 + cBox^2 (z - zp)^2]});
          integrand4D = cKL * integrandF;
          val = NIntegrate[integrand4D,
            {y, 0, 1}, {z, 0, 1}, {yp, 0, 1}, {zp, 0, 1},
            Method -> "GlobalAdaptive", WorkingPrecision -> 18,
            PrecisionGoal -> 5, AccuracyGoal -> 7, MaxRecursion -> 10];
          totals += aBox^2 val
        ],
        {entry, coeffsAlpha}
      ];
    ];
  ];

  (* === y-y component: phi_i_y phi_j_y / kernel ===
     Reduce y-direction analytically. *)
  Module[{phiIY = shi[[2]], phiJY = shj[[2]] /. {x -> xp, y -> yp, z -> zp}},
    If[phiIY =!= 0 && phiJY =!= 0,
      polyAlpha = Expand[phiIY phiJY];
      coeffsAlpha = CoefficientRules[polyAlpha, {y, yp}];
      Do[
        kpair = entry[[1]];
        cKL = entry[[2]];
        Module[{Fklval, integrandF, prefactor},
          Fklval = FklTable[{kpair[[1]], kpair[[2]]}];
          prefactor = 1 / bBox^(kpair[[1]] + kpair[[2]] + 2);
          integrandF = prefactor * (Fklval /. {len -> bBox,
                                                rho -> Sqrt[aBox^2 (x - xp)^2 + cBox^2 (z - zp)^2]});
          integrand4D = cKL * integrandF;
          val = NIntegrate[integrand4D,
            {x, 0, 1}, {z, 0, 1}, {xp, 0, 1}, {zp, 0, 1},
            Method -> "GlobalAdaptive", WorkingPrecision -> 18,
            PrecisionGoal -> 5, AccuracyGoal -> 7, MaxRecursion -> 10];
          totals += bBox^2 val
        ],
        {entry, coeffsAlpha}
      ];
    ];
  ];

  (* === z-z component === *)
  Module[{phiIZ = shi[[3]], phiJZ = shj[[3]] /. {x -> xp, y -> yp, z -> zp}},
    If[phiIZ =!= 0 && phiJZ =!= 0,
      polyAlpha = Expand[phiIZ phiJZ];
      coeffsAlpha = CoefficientRules[polyAlpha, {z, zp}];
      Do[
        kpair = entry[[1]];
        cKL = entry[[2]];
        Module[{Fklval, integrandF, prefactor},
          Fklval = FklTable[{kpair[[1]], kpair[[2]]}];
          prefactor = 1 / cBox^(kpair[[1]] + kpair[[2]] + 2);
          integrandF = prefactor * (Fklval /. {len -> cBox,
                                                rho -> Sqrt[aBox^2 (x - xp)^2 + bBox^2 (y - yp)^2]});
          integrand4D = cKL * integrandF;
          val = NIntegrate[integrand4D,
            {x, 0, 1}, {y, 0, 1}, {xp, 0, 1}, {yp, 0, 1},
            Method -> "GlobalAdaptive", WorkingPrecision -> 18,
            PrecisionGoal -> 5, AccuracyGoal -> 7, MaxRecursion -> 10];
          totals += cBox^2 val
        ],
        {entry, coeffsAlpha}
      ];
    ];
  ];

  (* prefactor: (mu0/4pi) * (abc)^2 (from Piola dV_phys^2) * (1/(abc))^2 (Piola components)
     = (mu0/4pi). The aBox^2, bBox^2, cBox^2 factors are absorbed in dotProduct.
  *)
  (mu0 / (4 Pi)) totals
];

(* === Run order test === *)
runHDivAharoni[order_] := Module[
  {shapesRef, n, MrefXX, MrefYY, MrefZZ, Mphys, K, t0, totalTime, counter, nUnique,
   eigvals, eigvecs, eigvecsNorm, bVec, gVec, gMag2, realEigs, tausUs,
   sortedByG, gMaxValue, gThresh, physicalIdx, sortedPhysical, Kmat},
  Print[];
  Print[StringRepeat["-", 60]];
  Print["HDiv div-free BEM (Aharoni-reduced), order ", order];
  Print[StringRepeat["-", 60]];

  shapesRef = ShapesHDivHexDivFree[order, x, y, z];
  n = Length[shapesRef];
  Print["DOFs: ", n];

  t0 = AbsoluteTime[];
  MrefXX = Table[fastIntegrate3D[Expand[shapesRef[[i, 1]] shapesRef[[j, 1]]], {x, y, z}], {i, n}, {j, n}];
  MrefYY = Table[fastIntegrate3D[Expand[shapesRef[[i, 2]] shapesRef[[j, 2]]], {x, y, z}], {i, n}, {j, n}];
  MrefZZ = Table[fastIntegrate3D[Expand[shapesRef[[i, 3]] shapesRef[[j, 3]]], {x, y, z}], {i, n}, {j, n}];
  Mphys = N[(aBox/(bBox cBox)) MrefXX + (bBox/(aBox cBox)) MrefYY + (cBox/(aBox bBox)) MrefZZ, 30];
  Print["M matrix in ", Round[AbsoluteTime[] - t0, 0.1], " s"];

  nUnique = n (n + 1) / 2;
  Print["Building K matrix (", nUnique, " unique entries) ..."];
  K = Table[0, {n}, {n}];
  totalTime = 0;
  counter = 0;
  Do[
    counter++;
    t0 = AbsoluteTime[];
    K[[i, j]] = KentryAharoni[shapesRef[[i]], shapesRef[[j]]];
    K[[j, i]] = K[[i, j]];
    totalTime += AbsoluteTime[] - t0;
    If[Mod[counter, 50] == 0,
      Print["  ", counter, "/", nUnique, "; ", Round[totalTime/60, 0.1], " min"]],
    {i, 1, n}, {j, i, n}
  ];
  Print["K matrix in ", Round[totalTime/60, 0.1], " min"];

  Kmat = N[K, 30];
  Print["Solving Eigensystem ..."];
  {eigvals, eigvecs} = Eigensystem[{Kmat, Mphys}];

  eigvecsNorm = Table[
    Module[{v = eigvecs[[k]], vMv},
      vMv = Re[v . Mphys . v];
      If[vMv > 10^-50, v / Sqrt[vMv], 0 v]],
    {k, 1, Length[eigvals]}];

  bVec = Table[
    Module[{phi = shapesRef[[i]]},
      (-aBox/2) fastIntegrate3D[Expand[phi[[1]] y aBox bBox cBox / (bBox cBox)], {x, y, z}] +
      (bBox/2)  fastIntegrate3D[Expand[phi[[2]] x aBox bBox cBox / (aBox cBox)], {x, y, z}]],
    {i, 1, n}];
  bVec = N[bVec, 30];
  gVec = Table[Re[eigvecsNorm[[k]] . bVec], {k, 1, Length[eigvals]}];
  gMag2 = gVec^2;

  realEigs = Re[eigvals];
  tausUs = Table[
    If[realEigs[[k]] > 10^-30, 10^6 sigma realEigs[[k]], 0],
    {k, 1, Length[eigvals]}];

  gMaxValue = Max[gMag2];
  gThresh = gMaxValue / 10^6;
  physicalIdx = Select[Range[Length[eigvals]],
    (gMag2[[#]] > gThresh) && (tausUs[[#]] > 0) &];
  sortedPhysical = SortBy[physicalIdx, -tausUs[[#]] &];

  Print[];
  Print["Top physical modes:"];
  Do[
    Module[{idx = sortedPhysical[[k]]},
      Print["  rank ", k, ":  tau = ", N[tausUs[[idx]], 12], " us, |g|^2 = ",
            N[gMag2[[idx]], 6]]],
    {k, 1, Min[10, Length[sortedPhysical]]}];
  Print[];
  Print["LEADING physical: tau = ",
        If[Length[sortedPhysical] > 0,
           N[tausUs[[sortedPhysical[[1]]]], 12], "none"], " us"];

  {n, If[Length[sortedPhysical] > 0, tausUs[[sortedPhysical[[1]]]], 0]}
];

(* === Run orders 1, 2, 3 === *)
results = {};
Do[
  Module[{r = runHDivAharoni[ord]},
    AppendTo[results, {ord, r[[1]], r[[2]]}]],
  {ord, 1, 3}
];

(* === Summary === *)
Print[];
Print[StringRepeat["=", 70]];
Print["HDiv div-free BEM convergence (Aharoni-reduced):"];
Print[StringRepeat["=", 70]];
Print["  order | DOFs | Leading tau (us)"];
Do[
  Print["    ", results[[k, 1]], "   |  ", results[[k, 2]], "  | ",
        N[results[[k, 3]], 12]],
  {k, 1, Length[results]}];
Print[];
Print["Phase 7-11 reference (P=4 K=8): 19.71 us"];
Print["Phase 7-11 reference (P=6):     19.77 us"];
Print[];
Print["Finished: ", DateString[]];
