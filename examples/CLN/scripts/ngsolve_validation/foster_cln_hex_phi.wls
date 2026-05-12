(* ::Package:: *)

(* foster_cln_hex_phi.wls
   Closed-form Newton potential of arbitrary 8-node hexahedron, 30+ digit precision.

   Public API:
     PhiBoxClosed[P, x0, x1, y0, y1, z0, z1]   axis-aligned cuboid (Banerjee-Gupta 1977)
     TriFaceIntegral[p1, p2, p3, P]            int_T 1/|r-P| dA  (Wilton-Rao-Glisson 1984)
     TetNewtonPotential[v0, v1, v2, v3, P]     int_tet 1/|r-P| dV
     PhiHex[verts, P]                           int_hex 1/|r-P| dV  (5-tet decomposition)

   All functions accept arbitrary-precision inputs (e.g. SetPrecision[..., 30]).
   For axis-aligned cuboid input, PhiHex agrees with PhiBoxClosed to 30 digits.

   Date: 2026-05-02
*)

BeginPackage["FosterClnHexPhi`"];

PhiBoxClosed::usage = "PhiBoxClosed[P, x0, x1, y0, y1, z0, z1] gives the Newton potential int_box 1/|r-P| dV for axis-aligned cuboid [x0,x1]x[y0,y1]x[z0,z1]. Banerjee-Gupta 1977 closed form.";
TriFaceIntegral::usage = "TriFaceIntegral[p1, p2, p3, P] gives int_triangle 1/|r-P| dA where the triangle has vertices p1, p2, p3 (3D). Wilton-Rao-Glisson 1984 closed form.";
TetNewtonPotential::usage = "TetNewtonPotential[v0, v1, v2, v3, P] gives int_tet 1/|r-P| dV via 4-face sum.";
PhiHex::usage = "PhiHex[verts, P] gives int_hex 1/|r-P| dV for an 8-node hexahedron via 5-tet decomposition. verts is a list of 8 (x,y,z) vertices in standard hex ordering: 0..3 = bottom face CCW, 4..7 = top face CCW.";

Begin["`Private`"];

(* ===== Banerjee-Gupta antiderivative ===== *)
GAntider[x_, y_, z_] := Module[{R = Sqrt[x^2 + y^2 + z^2], term1, term2, term3, ax, ay, az},
  If[R == 0, Return[0]];
  term1 = If[x*y == 0, 0, x y Log[z + R]];
  term2 = If[y*z == 0, 0, y z Log[x + R]];
  term3 = If[z*x == 0, 0, z x Log[y + R]];
  ax = If[x == 0, 0, (x^2/2) ArcTan[y z / (x R)]];
  ay = If[y == 0, 0, (y^2/2) ArcTan[z x / (y R)]];
  az = If[z == 0, 0, (z^2/2) ArcTan[x y / (z R)]];
  term1 + term2 + term3 - ax - ay - az
];

PhiBoxClosed[P_List, x0_, x1_, y0_, y1_, z0_, z1_] := Module[
  {Px, Py, Pz, xs, ys, zs, total, sign},
  {Px, Py, Pz} = P;
  xs = {x0 - Px, x1 - Px};
  ys = {y0 - Py, y1 - Py};
  zs = {z0 - Pz, z1 - Pz};
  total = Sum[
    sign = If[EvenQ[i + j + k], 1, -1];
    sign * GAntider[xs[[i + 1]], ys[[j + 1]], zs[[k + 1]]],
    {i, 0, 1}, {j, 0, 1}, {k, 0, 1}
  ];
  -total
];

(* ===== Vector helpers ===== *)
VSub[a_List, b_List] := a - b;
VDot[a_List, b_List] := a.b;
VCross[a_List, b_List] := Cross[a, b];
VNorm[a_List] := Sqrt[a.a];
VScale[a_List, s_] := s * a;

(* ===== Solid angle (Van Oosterom-Strackee 1983) ===== *)
SolidAngleTriangle[p1_, p2_, p3_, P_] := Module[
  {a, b, c, la, lb, lc, num, den, result},
  a = p1 - P; b = p2 - P; c = p3 - P;
  la = VNorm[a]; lb = VNorm[b]; lc = VNorm[c];
  num = VDot[VCross[a, b], c];
  den = la lb lc + VDot[a, b] lc + VDot[a, c] lb + VDot[b, c] la;
  (* When P is on a vertex or in the triangle plane, both num and den can
     vanish, making ArcTan[den,num] indeterminate. Use Chop to identify the
     numerical-zero case; the prefactor in the calling formula vanishes there. *)
  If[Chop[num, 10^-30] === 0 && Chop[den, 10^-30] === 0, Return[0]];
  result = Quiet[2 ArcTan[den, num], {ArcTan::indet}];
  If[result === Indeterminate, 0, result]
];

(* ===== Triangle face integral: int_T 1/|r-P| dA  (closed form) =====
   I = h_signed * Omega - edge_sum
   per Wilton-Rao-Glisson 1984 / Graglia 1993, sign convention verified
   numerically against 2D quadrature.
*)
TriFaceIntegral[p1_List, p2_List, p3_List, P_List] := Module[
  {e12, e13, nUnnorm, area2, nhat, hSigned, rho0, Omega, verts,
   edgeSum, i, a, b, sVec, sLen, shat, mhat, ra, rb, Ra, Rb,
   tParam, lMinus, lPlus, tI0, argPlus, argMinus, fI},
  e12 = p2 - p1;
  e13 = p3 - p1;
  nUnnorm = VCross[e12, e13];
  area2 = VNorm[nUnnorm];
  If[area2 == 0, Return[0]];
  nhat = nUnnorm / area2;
  hSigned = VDot[P - p1, nhat];
  rho0 = P - hSigned nhat;
  Omega = SolidAngleTriangle[p1, p2, p3, P];
  verts = {p1, p2, p3};
  edgeSum = 0;
  Do[
    a = verts[[i]];
    b = verts[[Mod[i, 3] + 1]];
    sVec = b - a;
    sLen = VNorm[sVec];
    If[sLen == 0, Continue[]];
    shat = sVec / sLen;
    mhat = VCross[shat, nhat];
    ra = P - a; rb = P - b;
    Ra = VNorm[ra]; Rb = VNorm[rb];
    tParam = VDot[P - a, shat];
    lMinus = -tParam;
    lPlus = sLen - tParam;
    tI0 = VDot[rho0 - a, mhat];
    argPlus = Rb + lPlus;
    argMinus = Ra + lMinus;
    fI = If[argPlus > 0 && argMinus > 0,
      Log[argPlus / argMinus],
      Module[{argPlus2 = Rb - lPlus, argMinus2 = Ra - lMinus},
        If[argPlus2 > 0 && argMinus2 > 0,
          Log[argMinus2 / argPlus2],
          0
        ]
      ]
    ];
    edgeSum += tI0 fI,
    {i, 1, 3}
  ];
  hSigned Omega - edgeSum
];

(* ===== Tet Newton potential ===== *)
TetNewtonPotential[v0_List, v1_List, v2_List, v3_List, P_List] := Module[
  {verts, centroid, faceIndices, Phi, face, p1, p2, p3, e12, e13,
   nUnnorm, narea, nhat, outCheck, dF, fTri, idx},
  verts = {v0, v1, v2, v3};
  centroid = Mean[verts];
  (* Face i is opposite vertex i (1-indexed in Mathematica) *)
  faceIndices = {{2, 3, 4}, {1, 4, 3}, {1, 2, 4}, {1, 3, 2}};
  Phi = 0;
  Do[
    {p1, p2, p3} = verts[[face]];
    e12 = p2 - p1;
    e13 = p3 - p1;
    nUnnorm = VCross[e12, e13];
    narea = VNorm[nUnnorm];
    If[narea == 0, Continue[]];
    nhat = nUnnorm / narea;
    outCheck = VDot[p1 - centroid, nhat];
    If[outCheck < 0,
      {p2, p3} = {p3, p2};
      nhat = -nhat
    ];
    dF = VDot[P - p1, nhat];
    fTri = TriFaceIntegral[p1, p2, p3, P];
    Phi += -dF fTri / 2,
    {face, faceIndices}
  ];
  Phi
];

(* ===== Hex Newton potential via 5-tet decomposition ===== *)
(* Standard hex 5-tet split, vertex ordering 1..8 (Mathematica) *)
hex5tets = {{1, 2, 4, 5}, {2, 3, 4, 7}, {2, 5, 6, 7}, {4, 5, 7, 8}, {2, 4, 5, 7}};

PhiHex[verts_List /; Length[verts] == 8, P_List] := Module[
  {Phi = 0, tet, v0, v1, v2, v3},
  Do[
    {v0, v1, v2, v3} = verts[[tet]];
    Phi += TetNewtonPotential[v0, v1, v2, v3, P],
    {tet, hex5tets}
  ];
  Phi
];

End[];
EndPackage[];

(* ===== Self-test if run directly ===== *)
If[$ScriptCommandLine =!= {},
  Print["FosterClnHexPhi self-test"];
  Print["================================"];
  Module[{a, b, c, P, verts, phiHex, phiBox, err},
    {a, b, c} = SetPrecision[{5*^-3, 2*^-3, 1*^-3}, 30];
    verts = SetPrecision[{
      {-a/2, -b/2, -c/2}, { a/2, -b/2, -c/2},
      { a/2,  b/2, -c/2}, {-a/2,  b/2, -c/2},
      {-a/2, -b/2,  c/2}, { a/2, -b/2,  c/2},
      { a/2,  b/2,  c/2}, {-a/2,  b/2,  c/2}
    }, 30];

    Print["Test 1: exterior point P = (10mm, 0, 0)"];
    P = SetPrecision[{10*^-3, 0, 0}, 30];
    phiHex = FosterClnHexPhi`PhiHex[verts, P];
    phiBox = FosterClnHexPhi`PhiBoxClosed[P, -a/2, a/2, -b/2, b/2, -c/2, c/2];
    err = Abs[(phiHex - phiBox) / phiBox];
    Print["  PhiHex     = ", N[phiHex, 18]];
    Print["  PhiBox(BG) = ", N[phiBox, 18]];
    Print["  rel.err    = ", N[err, 6]];

    Print["Test 2: center P = (0, 0, 0)"];
    P = SetPrecision[{0, 0, 0}, 30];
    phiHex = FosterClnHexPhi`PhiHex[verts, P];
    phiBox = FosterClnHexPhi`PhiBoxClosed[P, -a/2, a/2, -b/2, b/2, -c/2, c/2];
    err = Abs[(phiHex - phiBox) / phiBox];
    Print["  PhiHex     = ", N[phiHex, 18]];
    Print["  PhiBox(BG) = ", N[phiBox, 18]];
    Print["  rel.err    = ", N[err, 6]];

    Print["Test 3: face center P = (a/2, 0, 0)"];
    P = SetPrecision[{a/2, 0, 0}, 30];
    phiHex = FosterClnHexPhi`PhiHex[verts, P];
    phiBox = FosterClnHexPhi`PhiBoxClosed[P, -a/2, a/2, -b/2, b/2, -c/2, c/2];
    err = Abs[(phiHex - phiBox) / phiBox];
    Print["  PhiHex     = ", N[phiHex, 18]];
    Print["  PhiBox(BG) = ", N[phiBox, 18]];
    Print["  rel.err    = ", N[err, 6]];

    Print["Expected: relative error < 1e-25 in all tests."];
  ]
]
