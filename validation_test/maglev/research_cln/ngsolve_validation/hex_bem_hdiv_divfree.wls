(* ::Package:: *)

(* hex_bem_hdiv_divfree.wls
   BEM with HDiv Schoberl-Zaglmayr Type 1 INTERIOR div-free basis (loop subspace).
   For cuboid with all 6 faces on boundary, only interior DOFs matter.
   These are constructively divergence-free (volumetric loop basis).

   Reference: hdivhofe_impl.hpp lines 919-940 (HEX, ho_div_free=True, Type 1).

   For order p, interior div-free DOF count:
     2*p^3 (wDu_Cross_Dv shapes) + 3*p^2 (Du_Cross_Dv shapes) = 2p^3 + 3p^2
   Order 1: 2+3 = 5 DOFs
   Order 2: 16+12 = 28 DOFs
   Order 3: 54+27 = 81 DOFs

   Date: 2026-05-03
*)

ClearAll[x, y, z, xp, yp, zp];
$MaxExtraPrecision = 50;

Print[StringRepeat["=", 70]];
Print["BEM with HDiv interior div-free basis (volumetric loop)"];
Print[StringRepeat["=", 70]];

sigma = SetPrecision[58 10^6, 30];
mu0   = 4 Pi SetPrecision[10^-7, 30];
{aBox, bBox, cBox} = SetPrecision[{5 10^-3, 2 10^-3, 1 10^-3}, 30];

(* === HDiv helpers (return 3D vector) ===
   Du_Cross_Dv(u, v) = grad u x grad v
   wDu_Cross_Dv(u, v, w) = w * (grad u x grad v)
   Both are div-free (curl-form constructions).
*)
DuCrossDv[u_, v_, vars_] := Cross[D[u, {vars}], D[v, {vars}]];
wDuCrossDv[u_, v_, w_, vars_] := w Cross[D[u, {vars}], D[v, {vars}]];

(* === IntegratedLegendre (NGSolve naming) ===
   Returns L_2(x), L_3(x), ..., L_{n+1}(x). Initial values P_0 = -1, P_1 = x.
   These vanish at x = +/-1 for k >= 2.
*)
IntegratedLegendre[n_, x_] := Module[{p3, p2, p1, j, result},
  result = {};
  If[n < 2, Return[{}]];
  p3 = 0; p2 = -1; p1 = x;
  Do[
    p3 = p2; p2 = p1;
    p1 = (2 j - 3)/j * x * p2 - (j - 3)/j * p3;
    AppendTo[result, p1],
    {j, 2, n}
  ];
  result
];

(* === LegendrePolynomial (returns P_0, P_1, ..., P_n: n+1 entries) === *)
LegendrePoly[n_, x_] := Module[{P, k},
  P = Table[0, {n + 1}];
  P[[1]] = 1;
  If[n >= 1, P[[2]] = x];
  Do[
    P[[k + 1]] = (2 - 1/k) x P[[k]] + (1/k - 1) P[[k - 1]],
    {k, 2, n}
  ];
  P
];


(* === HDiv HEX interior div-free shapes (Type 1) ===
   For order_inner = (p0, p1, p2), in reference [0,1]^3.
   Mapping: xi = 2x-1, eta = 2y-1, zeta = 2z-1.

   From hdivhofe_impl.hpp lines 919-940:
     L_xi[k+2] etc are IntegratedLegendre values starting at index k+2 = 2..p+1.
     pol_xi[i+1] etc are LegendrePolynomial values starting at index i+1 = 1..p.

   Type 1 div-free interior shapes (using ho_div_free=True):
   1. wDu_Cross_Dv(L_xi[i+2]*L_zeta[k+2], y, -4 P_eta[j+1])  for i, j, k = 0..p-1
   2. wDu_Cross_Dv(L_eta[j+2]*L_zeta[k+2], x, +4 P_xi[i+1])
   3. Du_Cross_Dv(L_eta[j+2]*L_zeta[k+2], 2*xi)
   4. Du_Cross_Dv(L_xi[i+2]*L_zeta[k+2], 2*eta)
   5. Du_Cross_Dv(L_xi[i+2]*L_eta[j+2], 2*zeta)
*)

ShapesHDivHexDivFree[order_, x_, y_, z_] := Module[
  {p, xi, eta, zeta, LXi, LEta, LZeta, PXi, PEta, PZeta, shapes, vars},
  vars = {x, y, z};
  p = order;
  xi = 2 x - 1;
  eta = 2 y - 1;
  zeta = 2 z - 1;

  (* Need IntegratedLegendre indices 2..p+1, so call with n = p+1, returns p entries *)
  LXi   = IntegratedLegendre[p + 1, xi];   (* p entries: L_2, L_3, ..., L_{p+1} *)
  LEta  = IntegratedLegendre[p + 1, eta];
  LZeta = IntegratedLegendre[p + 1, zeta];

  (* LegendrePolynomial indices 1..p, so call with n = p, returns p+1 entries (0..p) *)
  PXi   = LegendrePoly[p, xi];   (* p+1 entries; we use index [i+1] = 1..p *)
  PEta  = LegendrePoly[p, eta];
  PZeta = LegendrePoly[p, zeta];

  shapes = {};

  (* Type 1a: wDu_Cross_Dv(L_xi[i+2]*L_zeta[k+2], y, -4*P_eta[j+1]) for i,j,k = 0..p-1 *)
  Do[
    AppendTo[shapes,
      wDuCrossDv[LXi[[i+1]] LZeta[[k+1]], y, -4 PEta[[j+2]], vars]],
    {i, 0, p-1}, {j, 0, p-1}, {k, 0, p-1}
  ];

  (* Type 1b: wDu_Cross_Dv(L_eta[j+2]*L_zeta[k+2], x, +4*P_xi[i+1]) *)
  Do[
    AppendTo[shapes,
      wDuCrossDv[LEta[[j+1]] LZeta[[k+1]], x, 4 PXi[[i+2]], vars]],
    {i, 0, p-1}, {j, 0, p-1}, {k, 0, p-1}
  ];

  (* Type 1c: Du_Cross_Dv(L_eta[j+2]*L_zeta[k+2], 2*xi) for j, k = 0..p-1 *)
  Do[
    AppendTo[shapes,
      DuCrossDv[LEta[[j+1]] LZeta[[k+1]], 2 xi, vars]],
    {j, 0, p-1}, {k, 0, p-1}
  ];

  (* Type 1d: Du_Cross_Dv(L_xi[i+2]*L_zeta[k+2], 2*eta) *)
  Do[
    AppendTo[shapes,
      DuCrossDv[LXi[[i+1]] LZeta[[k+1]], 2 eta, vars]],
    {i, 0, p-1}, {k, 0, p-1}
  ];

  (* Type 1e: Du_Cross_Dv(L_xi[i+2]*L_eta[j+2], 2*zeta) *)
  Do[
    AppendTo[shapes,
      DuCrossDv[LXi[[i+1]] LEta[[j+1]], 2 zeta, vars]],
    {i, 0, p-1}, {j, 0, p-1}
  ];

  shapes
];

(* === Test: order 1 should give 5 DOFs, order 2 should give 28 === *)
Print["Order 1 div-free interior DOFs: ", Length[ShapesHDivHexDivFree[1, x, y, z]]];
Print["Order 2 div-free interior DOFs: ", Length[ShapesHDivHexDivFree[2, x, y, z]]];

(* === Verify div = 0 for one sample shape === *)
Print[];
Print["Div-free verification (one sample shape from order 2):"];
sample = ShapesHDivHexDivFree[2, x, y, z][[1]];
divCheck = Simplify[D[sample[[1]], x] + D[sample[[2]], y] + D[sample[[3]], z]];
Print["  div(sample) = ", divCheck];

(* === Build BEM K and M matrices === *)
fastIntegrate3D[poly_, vars_] := Module[{p = Expand[poly], rules},
  If[p === 0, Return[0]];
  rules = CoefficientRules[p, vars];
  Total[Function[r, r[[2]] / Times @@ Map[(# + 1) &, r[[1]]]] /@ rules]
];


runHDivBEM[order_] := Module[
  {shapesRef, n, MrefXX, MrefYY, MrefZZ, Mphys, Kentry, K, Kmat,
   eigvals, eigvecs, eigvecsNorm, bVec, gVec, gMag2, realEigs, tausUs,
   sortedByG, t0, totalTime, counter, nUnique},
  Print[];
  Print[StringRepeat["-", 60]];
  Print["HDiv div-free BEM, order ", order];
  Print[StringRepeat["-", 60]];

  shapesRef = ShapesHDivHexDivFree[order, x, y, z];
  n = Length[shapesRef];
  Print["DOFs (interior div-free): ", n];

  (* Mass matrix - using HDiv contravariant Piola: phi_phys = (1/detJ) J phi_ref
     For diagonal J = diag(a, b, c): phi_phys_x = phi_ref_x / (bc), etc.
     M_phys[i,j] = (a/(bc)) Mref_xx + (b/(ac)) Mref_yy + (c/(ab)) Mref_zz *)
  t0 = AbsoluteTime[];
  MrefXX = Table[fastIntegrate3D[Expand[shapesRef[[i, 1]] shapesRef[[j, 1]]], {x, y, z}], {i, n}, {j, n}];
  MrefYY = Table[fastIntegrate3D[Expand[shapesRef[[i, 2]] shapesRef[[j, 2]]], {x, y, z}], {i, n}, {j, n}];
  MrefZZ = Table[fastIntegrate3D[Expand[shapesRef[[i, 3]] shapesRef[[j, 3]]], {x, y, z}], {i, n}, {j, n}];
  Mphys = N[(aBox/(bBox cBox)) MrefXX + (bBox/(aBox cBox)) MrefYY + (cBox/(aBox bBox)) MrefZZ, 30];
  Print["Mass matrices in ", Round[AbsoluteTime[] - t0, 0.1], " s"];

  (* K matrix: BEM Newton kernel
     K[i,j] = (mu0/(4pi)) (abc)^2 int_[0,1]^6 (vi.vj_with_piola) / |r_phys - r'_phys| dV dV'
     |r_phys - r'_phys| = sqrt(a^2 dx^2 + b^2 dy^2 + c^2 dz^2)
     v_i_phys . v_j_phys = (a^2/(abc)^2) v_ix v_jx + (b^2/(abc)^2) v_iy v_jy + (c^2/(abc)^2) v_iz v_jz
                         = (1/(bc)^2) v_ix v_jx + ...
  *)
  Kentry[i_, j_] := Module[{vi, vj, dotProduct, kernel, factor},
    vi = shapesRef[[i]];
    vj = shapesRef[[j]] /. {x -> xp, y -> yp, z -> zp};
    dotProduct = vi[[1]] vj[[1]] aBox^2 + vi[[2]] vj[[2]] bBox^2 + vi[[3]] vj[[3]] cBox^2;
    kernel = 1 / Sqrt[aBox^2 (x - xp)^2 + bBox^2 (y - yp)^2 + cBox^2 (z - zp)^2];
    factor = aBox^2 bBox^2 cBox^2 / (aBox bBox cBox)^2;   (* simplifies to 1, kept for clarity *)
    (mu0 / (4 Pi)) NIntegrate[
      dotProduct kernel,
      {x, 0, 1}, {y, 0, 1}, {z, 0, 1}, {xp, 0, 1}, {yp, 0, 1}, {zp, 0, 1},
      Method -> "GlobalAdaptive", WorkingPrecision -> 15,
      PrecisionGoal -> 4, AccuracyGoal -> 6, MaxRecursion -> 8
    ]
  ];

  Print[];
  nUnique = n (n + 1) / 2;
  Print["Building K matrix (", nUnique, " unique entries) ..."];
  K = Table[0, {n}, {n}];
  totalTime = 0;
  counter = 0;
  Do[
    counter++;
    t0 = AbsoluteTime[];
    K[[i, j]] = Kentry[i, j];
    K[[j, i]] = K[[i, j]];
    totalTime += AbsoluteTime[] - t0;
    If[Mod[counter, 100] == 0,
      Print["  ", counter, "/", nUnique, "; ", Round[totalTime/60, 0.1], " min"]],
    {i, 1, n}, {j, i, n}
  ];
  Print["K matrix in ", Round[totalTime/60, 0.1], " min"];
  Kmat = N[K, 30];

  Print[];
  Print["Solving Eigensystem ..."];
  {eigvals, eigvecs} = Eigensystem[{Kmat, Mphys}];

  eigvecsNorm = Table[
    Module[{v = eigvecs[[k]], vMv},
      vMv = Re[v . Mphys . v];
      If[vMv > 10^-50, v / Sqrt[vMv], 0 v]],
    {k, 1, Length[eigvals]}];

  (* Project onto A_ext = (-y/2, x/2, 0) *)
  bVec = Table[
    Module[{phi = shapesRef[[i]]},
      (* phi_phys . A_ext_phys = (a/(bc)) phi_x_ref (-b y/2) + (b/(ac)) phi_y_ref (a x/2) *)
      (-aBox/2) fastIntegrate3D[Expand[phi[[1]] y aBox bBox cBox / (bBox cBox)], {x, y, z}] +
      (bBox/2)  fastIntegrate3D[Expand[phi[[2]] x aBox bBox cBox / (aBox cBox)], {x, y, z}]
    ],
    {i, 1, n}
  ];
  bVec = N[bVec, 30];
  gVec = Table[Re[eigvecsNorm[[k]] . bVec], {k, 1, Length[eigvals]}];
  gMag2 = gVec^2;

  realEigs = Re[eigvals];
  tausUs = Table[
    If[realEigs[[k]] > 10^-30, 10^6 sigma realEigs[[k]], 0],
    {k, 1, Length[eigvals]}];

  sortedByG = Reverse[Ordering[gMag2]];
  Print[];
  Print["Top modes by |g_k|^2:"];
  Do[
    Module[{idx = sortedByG[[k]]},
      Print["  rank ", k, ":  tau = ", N[tausUs[[idx]], 10], " us,  |g|^2 = ",
            N[gMag2[[idx]], 6]]],
    {k, 1, Min[8, Length[sortedByG]]}
  ];

  Print[];
  Print["LEADING physical mode: tau = ", N[tausUs[[sortedByG[[1]]]], 10], " us"];

  {n, tausUs[[sortedByG[[1]]]]}
];

(* === Run for orders 1, 2, 3 === *)
results = {};
Do[
  Module[{r = runHDivBEM[ord]},
    AppendTo[results, {ord, r[[1]], r[[2]]}]],
  {ord, 1, 3}
];

(* === Summary === *)
Print[];
Print[StringRepeat["=", 70]];
Print["HDiv div-free BEM convergence summary:"];
Print[StringRepeat["=", 70]];
Print[];
Print["  order | DOFs | Leading tau (us)"];
Do[
  Print["    ", results[[k, 1]], "   |  ", results[[k, 2]], "  | ",
        N[results[[k, 3]], 10]],
  {k, 1, Length[results]}
];
Print[];
Print["Phase 7-11 reference (P=4 K=8): 19.71 us"];
Print["Phase 7-11 reference (P=6):     19.77 us"];
Print["HCurl order 1 (54 DOFs):        19.43 us  (~2% off, no div-free)"];
