(* ::Package:: *)

(* hex_bem_aharoni.wls
   BEM with HCurl lowest-order edge basis using Aharoni 1D-closed inner trick.

   For x-direction edges (4 of them), basis is phi = (p(y,z)/a, 0, 0) with
   p a bilinear in (y,z). The kernel 1/|r-r'| with x integrated out:
     F_x(dy, dz; a) = int_0^a int_0^a 1/sqrt((x-x')^2 + dy^2 + dz^2) dx dx'
                   = 2 [a sinh^-1(a/rho) + rho - sqrt(a^2 + rho^2)]
   where rho = sqrt(dy^2 + dz^2). Logarithmic singularity at rho=0.

   K_xx[i,j] = (mu0 / 4pi) (1/a^2) int_y int_z int_yp int_zp
                p_i(y,z) p_j(yp,zp) F_x(y-yp, z-zp; a) dy dz dyp dzp

   M_xx[i,j] = (1/a) int p_i(y,z) p_j(y,z) dy dz   (analytical)

   Similar for y- and z-direction edges. Total 12x12 K, M matrices,
   block-diagonal (3 blocks of 4x4).

   Solve generalized eigenproblem -> tau_n = mu0 sigma lambda_n / (4 pi).

   Date: 2026-05-02
*)

ClearAll[x, y, z, xp, yp, zp, dy, dz, du, dv];
$MaxExtraPrecision = 50;

(* Physical cuboid 5 x 2 x 1 mm Cu *)
sigma = SetPrecision[58 10^6, 30];
mu0   = 4 Pi SetPrecision[10^-7, 30];
{aBox, bBox, cBox} = SetPrecision[{5 10^-3, 2 10^-3, 1 10^-3}, 30];

(* === Aharoni 1D-closed kernel === *)
(* F[du, dv, len] = int_0^len int_0^len 1/sqrt((x-x')^2 + du^2 + dv^2) dx dx'
                 = 2[len ArcSinh(len/rho) + rho - sqrt(len^2 + rho^2)]
   where rho = sqrt(du^2 + dv^2). *)
Fkernel[du_, dv_, len_] := Module[{rho2 = du^2 + dv^2, rho},
  rho = Sqrt[rho2];
  2 (len ArcSinh[len/rho] + rho - Sqrt[len^2 + rho^2])
];

(* === HCurl lowest-order edge polynomial weights ===
   x-direction edges (in physical coords, after Piola scaling 1/a):
     edge 0 (V0-V1): p = (1-y/b)(1-z/c)
     edge 1 (V2-V3): p = -(y/b)(1-z/c)    [sign from edge orientation]
     edge 4 (V4-V5): p = (1-y/b)(z/c)
     edge 5 (V6-V7): p = -(y/b)(z/c)
*)

edgePolysX = {
  (1 - y/bBox) (1 - z/cBox),
  -(y/bBox) (1 - z/cBox),
  (1 - y/bBox) (z/cBox),
  -(y/bBox) (z/cBox)
};

edgePolysY = {
  -(1 - x/aBox) (1 - z/cBox),    (* edge 2 (V3-V0): -y direction *)
  (x/aBox) (1 - z/cBox),         (* edge 3 (V1-V2): +y direction *)
  -(1 - x/aBox) (z/cBox),        (* edge 6 (V7-V4) *)
  (x/aBox) (z/cBox)              (* edge 7 (V5-V6) *)
};

edgePolysZ = {
  (1 - x/aBox) (1 - y/bBox),     (* edge 8 (V0-V4) *)
  (x/aBox) (1 - y/bBox),         (* edge 9 (V1-V5) *)
  (x/aBox) (y/bBox),             (* edge 10 (V2-V6) *)
  (1 - x/aBox) (y/bBox)          (* edge 11 (V3-V7) *)
};

(* === Mass matrix blocks (analytical) ===
   M_xx[i,j] = (1/a) int_0^b int_0^c p_i(y,z) p_j(y,z) dy dz  *)
massBlock[polys_, varA_, lenA_, varB_, lenB_, lenPerp_] := Table[
  Integrate[polys[[i]] polys[[j]], {varA, 0, lenA}, {varB, 0, lenB}] / lenPerp,
  {i, 1, Length[polys]}, {j, 1, Length[polys]}
];

Mxx = massBlock[edgePolysX, y, bBox, z, cBox, aBox];
Myy = massBlock[edgePolysY, x, aBox, z, cBox, bBox];
Mzz = massBlock[edgePolysZ, x, aBox, y, bBox, cBox];

Print["Mass matrix blocks computed."];
Print["Mxx diagonal: ", Diagonal[N[Mxx, 12]]];
Print["Myy diagonal: ", Diagonal[N[Myy, 12]]];
Print["Mzz diagonal: ", Diagonal[N[Mzz, 12]]];

(* === K matrix blocks via 4D NIntegrate (Aharoni-reduced) === *)
kBlock[polys_, varA_, lenA_, varB_, lenB_, lenPerp_, WP_:20] := Module[
  {n = Length[polys], K, t0, total = 0, tot = 0},
  K = Table[0, {n}, {n}];
  Do[
    t0 = AbsoluteTime[];
    K[[i, j]] = (mu0 / (4 Pi)) (1/lenPerp^2) NIntegrate[
      polys[[i]] *
        (polys[[j]] /. {varA -> du + varA, varB -> dv + varB}) *
        Fkernel[-du, -dv, lenPerp] *
        DiracDelta[0],   (* placeholder - rewrite below *)
      Method -> "GlobalAdaptive", WorkingPrecision -> WP
    ];
    K[[j, i]] = K[[i, j]];
    Print["  K[", i, "][", j, "] in ", Round[AbsoluteTime[] - t0, 0.1], " s"],
    {i, 1, n}, {j, i, n}
  ];
  K
];

(* Better: direct 4D NIntegrate with substitution variables *)
kBlockDirect[polys_, varA_, lenA_, varAp_, varB_, lenB_, varBp_, lenPerp_, WP_:20] := Module[
  {n = Length[polys], K, t0},
  K = Table[0, {n}, {n}];
  Do[
    t0 = AbsoluteTime[];
    K[[i, j]] = (mu0 / (4 Pi)) (1/lenPerp^2) NIntegrate[
      polys[[i]] *
        (polys[[j]] /. {varA -> varAp, varB -> varBp}) *
        Fkernel[varA - varAp, varB - varBp, lenPerp],
      {varA, 0, lenA}, {varB, 0, lenB},
      {varAp, 0, lenA}, {varBp, 0, lenB},
      Method -> "GlobalAdaptive", WorkingPrecision -> WP, MaxRecursion -> 20,
      AccuracyGoal -> 8
    ];
    K[[j, i]] = K[[i, j]];
    Print["  K[", i, "][", j, "] = ", N[K[[i, j]], 8], "  (", Round[AbsoluteTime[] - t0, 0.1], " s)"],
    {i, 1, n}, {j, i, n}
  ];
  K
];

Print[];
Print["===== Computing K_xx (4x4) ====="];
t0 = AbsoluteTime[];
Kxx = kBlockDirect[edgePolysX, y, bBox, yp, z, cBox, zp, aBox, 20];
Print["K_xx done in ", Round[AbsoluteTime[] - t0, 0.1], " s"];

Print[];
Print["===== Computing K_yy (4x4) ====="];
t0 = AbsoluteTime[];
Kyy = kBlockDirect[edgePolysY, x, aBox, xp, z, cBox, zp, bBox, 20];
Print["K_yy done in ", Round[AbsoluteTime[] - t0, 0.1], " s"];

Print[];
Print["===== Computing K_zz (4x4) ====="];
t0 = AbsoluteTime[];
Kzz = kBlockDirect[edgePolysZ, x, aBox, xp, y, bBox, yp, cBox, 20];
Print["K_zz done in ", Round[AbsoluteTime[] - t0, 0.1], " s"];

(* === Solve generalized eigenproblems for each block === *)
solveBlock[K_, M_, label_] := Module[{eigvals},
  Print[];
  Print[label, ":"];
  eigvals = Eigenvalues[{K, M}];
  Print["  lambda eigenvalues: ", N[Sort[Re[eigvals]], 10]];
  Print["  tau (us) = sigma * lambda  (mu0/4pi already in K):"];
  Do[Print["    tau = ", N[10^6 sigma * eigvals[[k]], 8], " us"], {k, 1, Length[eigvals]}];
  eigvals
];

eigsX = solveBlock[Kxx, Mxx, "X-block (4 modes)"];
eigsY = solveBlock[Kyy, Myy, "Y-block (4 modes)"];
eigsZ = solveBlock[Kzz, Mzz, "Z-block (4 modes)"];

(* === Combined results === *)
Print[];
Print["============================================================"];
Print["All 12 lowest-order edge BEM eigenvalues -> tau_n (us):"];
Print["============================================================"];
allEigs = Sort[Re[Join[eigsX, eigsY, eigsZ]]];
coef = mu0 sigma / (4 Pi);
Do[Print["  tau_", k, " = ", N[10^6 coef * allEigs[[k]], 10], " us"], {k, 1, Length[allEigs]}];

Print[];
Print["Largest tau (leading mode): ", N[10^6 coef * Max[Re[Join[eigsX, eigsY, eigsZ]]], 8], " us"];
Print["Phase 7-11 reference: tau_1 = 19.77 us"];
