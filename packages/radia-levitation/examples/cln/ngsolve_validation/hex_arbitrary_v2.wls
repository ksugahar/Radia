(* ::Package:: *)

(* hex_arbitrary_v2.wls
   Arbitrary hex BEM via reference-coord component matrices + J-induced metric.
   Reuses fast approach from hex_bem_full_hcurl.wls (cuboid working).

   For arbitrary hex with constant Jacobian J (parallelepiped, including cuboid):
     - Metric G = J^T J  (3x3 constant matrix)
     - |r_phys - r'_phys|^2 = (xi - xi')^T G (xi - xi')
     - phi_phys = J^{-T} phi_ref  (covariant)
     - phi_phys . phi_phys' = phi_ref^T (J^{-1} J^{-T}) phi_ref'
     - Mass: M[i,j] = det J * sum_{ab} (J^{-1} J^{-T})_{ab} * Mref_ab[i,j]
     - K: similarly but with kernel scaling

   For TRUE trilinear hex (J varies with xi,eta,zeta), this framework needs
   numerical quadrature; out of scope for this PoC.

   Date: 2026-05-03
*)

Get["w:\\30_CauerLadderNetwork\\2026_04_01_\:9577\:65b9\:5f62CLN\\ngsolve_validation\\schoberl_zaglmayr_basis.wls"];

ClearAll[x, y, z, xp, yp, zp];
$MaxExtraPrecision = 50;

Print[StringRepeat["=", 70]];
Print["Arbitrary hex BEM via constant-J framework (cuboid + parallelepiped)"];
Print[StringRepeat["=", 70]];

sigma = SetPrecision[58 10^6, 30];
mu0   = 4 Pi SetPrecision[10^-7, 30];

(* === HCurl order 1 reference basis === *)
ord = 1;
shapesRef = SchoberlZaglmayrBasis`ShapesHEX[ord, Range[0, 7], x, y, z];
nDOF = Length[shapesRef];
Print["HCurl order ", ord, ": ", nDOF, " DOFs"];

(* === Reference component mass matrices === *)
fastIntegrate3D[poly_, vars_] := Module[{p = Expand[poly], rules},
  If[p === 0, Return[0]];
  rules = CoefficientRules[p, vars];
  Total[Function[r, r[[2]] / Times @@ Map[(# + 1) &, r[[1]]]] /@ rules]
];

componentMassRef[shapes_, comp1_, comp2_] := Module[{n = Length[shapes], M},
  Table[
    fastIntegrate3D[Expand[shapes[[i, 1, comp1]] shapes[[j, 1, comp2]]], {x, y, z}],
    {i, 1, n}, {j, 1, n}
  ]
];

t0 = AbsoluteTime[];
MrefXX = componentMassRef[shapesRef, 1, 1];
MrefYY = componentMassRef[shapesRef, 2, 2];
MrefZZ = componentMassRef[shapesRef, 3, 3];
MrefXY = componentMassRef[shapesRef, 1, 2];
MrefXZ = componentMassRef[shapesRef, 1, 3];
MrefYZ = componentMassRef[shapesRef, 2, 3];
Print["6 component reference mass matrices in ", Round[AbsoluteTime[] - t0, 0.1], " s"];

(* Symmetrize cross terms *)
MrefXYsym = (MrefXY + Transpose[MrefXY])/2;
MrefXZsym = (MrefXZ + Transpose[MrefXZ])/2;
MrefYZsym = (MrefYZ + Transpose[MrefYZ])/2;

(* === Build M_phys for given Jacobian J === *)
buildMphys[J_] := Module[
  {detJ = Det[J], JinvJinvT, M},
  JinvJinvT = Inverse[J] . Transpose[Inverse[J]];   (* J^{-1} J^{-T} *)
  M = detJ * (
    JinvJinvT[[1, 1]] MrefXX +
    JinvJinvT[[2, 2]] MrefYY +
    JinvJinvT[[3, 3]] MrefZZ +
    2 JinvJinvT[[1, 2]] MrefXYsym +
    2 JinvJinvT[[1, 3]] MrefXZsym +
    2 JinvJinvT[[2, 3]] MrefYZsym
  );
  N[M, 30]
];

(* === Build K_phys via 6D NIntegrate with constant-J metric === *)
buildKphys[J_, WP_:15] := Module[
  {detJ, JinvJinvT, G, K, n = nDOF, t0, total = 0, counter = 0,
   vi, vj, integrand, kernel, factor, dotProduct},
  detJ = Det[J];
  JinvJinvT = Inverse[J] . Transpose[Inverse[J]];
  G = Transpose[J] . J;   (* metric: |r_phys - r'_phys|^2 = (xi-xi')^T G (xi-xi') *)

  K = Table[0, {n}, {n}];
  Do[
    counter++;
    t0 = AbsoluteTime[];
    vi = shapesRef[[i, 1]];
    vj = shapesRef[[j, 1]] /. {x -> xp, y -> yp, z -> zp};
    (* dot product with J^{-1} J^{-T} metric *)
    dotProduct = vi . JinvJinvT . vj;
    factor = detJ^2;
    kernel = 1 / Sqrt[
      G[[1, 1]] (x - xp)^2 + G[[2, 2]] (y - yp)^2 + G[[3, 3]] (z - zp)^2 +
      2 G[[1, 2]] (x - xp) (y - yp) + 2 G[[1, 3]] (x - xp) (z - zp) +
      2 G[[2, 3]] (y - yp) (z - zp)
    ];
    integrand = dotProduct factor kernel;
    K[[i, j]] = (mu0 / (4 Pi)) NIntegrate[integrand,
      {x, 0, 1}, {y, 0, 1}, {z, 0, 1},
      {xp, 0, 1}, {yp, 0, 1}, {zp, 0, 1},
      Method -> "GlobalAdaptive", WorkingPrecision -> WP,
      PrecisionGoal -> 4, AccuracyGoal -> 6, MaxRecursion -> 8];
    K[[j, i]] = K[[i, j]];
    total += AbsoluteTime[] - t0;
    If[Mod[counter, 200] == 0,
      Print["  ", counter, "/", n(n+1)/2, "; ", Round[total/60, 0.1], " min"]],
    {i, 1, n}, {j, i, n}
  ];
  N[K, 30]
];

(* === Project A_ext = (-y_phys/2, x_phys/2, 0) onto basis ===
   r_phys = J . xi_ref. So x_phys = J11*x + J12*y + J13*z, y_phys = J21*x + J22*y + J23*z.
   A_ext_phys = (-y_phys/2, x_phys/2, 0).
   A_ext_ref via covariant: A_ref = J^T A_phys.
   <phi_i, A_ext>_phys = int phi_phys . A_phys det J dV_ref
                       = int (J^{-T} phi_ref) . A_phys det J dV_ref
                       = det J * int phi_ref . (J^{-1} A_phys) dV_ref *)
projectAext[J_] := Module[
  {detJ = Det[J], Jinv = Inverse[J], rPhys, AextPhys, AextRef, b, integrand},
  rPhys = J . {x, y, z};
  AextPhys = {-rPhys[[2]]/2, rPhys[[1]]/2, 0};
  AextRef = Jinv . AextPhys;
  b = Table[
    integrand = shapesRef[[i, 1]] . AextRef * detJ;
    fastIntegrate3D[Expand[integrand], {x, y, z}],
    {i, 1, nDOF}
  ];
  N[b, 30]
];


(* === Run for given Jacobian === *)
runJ[J_, label_] := Module[
  {Mphys, Kphys, eigvals, eigvecs, eigvecsNorm, bVec, gVec, gMag2, realEigs, tausUs,
   sortedByG, t0},
  Print[];
  Print[StringRepeat["-", 60]];
  Print[label];
  Print["J = ", N[J, 6]];
  Print["det J = ", N[Det[J], 12], " m^3"];
  Print[StringRepeat["-", 60]];

  t0 = AbsoluteTime[];
  Mphys = buildMphys[J];
  Print["M_phys assembled in ", Round[AbsoluteTime[] - t0, 0.1], " s"];

  t0 = AbsoluteTime[];
  Kphys = buildKphys[J, 15];
  Print["K_phys assembled in ", Round[(AbsoluteTime[] - t0)/60, 0.1], " min"];

  t0 = AbsoluteTime[];
  {eigvals, eigvecs} = Eigensystem[{Kphys, Mphys}];
  Print["Eigensystem in ", Round[AbsoluteTime[] - t0, 0.1], " s"];

  eigvecsNorm = Table[
    Module[{v = eigvecs[[k]], vMv},
      vMv = Re[v . Mphys . v];
      If[vMv > 10^-50, v / Sqrt[vMv], 0 v]],
    {k, 1, Length[eigvals]}];

  bVec = projectAext[J];
  gVec = Table[Re[eigvecsNorm[[k]] . bVec], {k, 1, Length[eigvals]}];
  gMag2 = gVec^2;

  realEigs = Re[eigvals];
  tausUs = Table[
    If[realEigs[[k]] > 10^-30, 10^6 sigma realEigs[[k]], 0],
    {k, 1, Length[eigvals]}];

  sortedByG = Reverse[Ordering[gMag2]];
  Print[];
  Print["Top 5 modes by |g_k|^2:"];
  Do[
    Module[{idx = sortedByG[[k]]},
      Print["  rank ", k, ":  tau = ", N[tausUs[[idx]], 8], " us,  |g|^2 = ",
            N[gMag2[[idx]], 6]]],
    {k, 1, Min[5, Length[sortedByG]]}
  ];
  Print[];
  Print["LEADING physical mode: tau = ", N[tausUs[[sortedByG[[1]]]], 8], " us"];
  tausUs[[sortedByG[[1]]]]
];


(* === Test cases === *)

aBox = SetPrecision[5 10^-3, 30];
bBox = SetPrecision[2 10^-3, 30];
cBox = SetPrecision[1 10^-3, 30];

(* Cuboid: J = diag(a, b, c) *)
JCuboid = DiagonalMatrix[{aBox, bBox, cBox}];
tauCuboid = runJ[JCuboid, "Test A: Cuboid 5x2x1 mm (axis-aligned, baseline)"];

(* Parallelepiped: top shifted +0.5mm in x. J has off-diagonal (1,3) = skew/cBox * cBox = skew? *)
(* Actually for parallelepiped with verts as I defined: r_phys = (a*xi + skew*zeta, b*eta, c*zeta).
   J = d r_phys / d (xi, eta, zeta) = [[a, 0, skew], [0, b, 0], [0, 0, c]]. *)
skew = SetPrecision[0.5 10^-3, 30];
JParallel = {{aBox, 0, skew}, {0, bBox, 0}, {0, 0, cBox}};
tauParallel = runJ[JParallel, "Test B: Parallelepiped (top shifted +0.5mm in x)"];

(* Larger skew test *)
skew2 = SetPrecision[1.0 10^-3, 30];
JParallel2 = {{aBox, 0, skew2}, {0, bBox, 0}, {0, 0, cBox}};
tauParallel2 = runJ[JParallel2, "Test C: Parallelepiped (top shifted +1.0mm = 20% deformation)"];


Print[];
Print[StringRepeat["=", 70]];
Print["Comparison summary (LEADING tau in us):"];
Print[StringRepeat["=", 70]];
Print["  Cuboid baseline:           ", N[tauCuboid, 8], " us  (Phase 7-11 ref: 19.77)"];
Print["  Parallelepiped 0.5mm skew: ", N[tauParallel, 8], " us  (NEW)"];
Print["  Parallelepiped 1.0mm skew: ", N[tauParallel2, 8], " us  (NEW)"];
Print[];
Print["HCurl × BEM × isoparametric mapping → arbitrary hex eddy current."];
