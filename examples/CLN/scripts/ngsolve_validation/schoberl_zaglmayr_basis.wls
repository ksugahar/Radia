(* ::Package:: *)

(* schoberl_zaglmayr_basis.wls
   High-order H(curl) Nedelec basis (Schöberl-Zaglmayr 2005), ported from NGSolve
   src/fem/hcurlhofe_impl.hpp and recursive_pol*.hpp at arbitrary precision.

   Status (2026-05-02 initial commit):
     [DONE] Helper operators (Du, uDv, uDv_minus_vDu, wuDv_minus_wvDu)
     [DONE] LegendrePolynomial (3-term recursion)
     [DONE] IntLegNoBubble    (3-term recursion)
     [DONE] DubinerBasis 2D (triangle interior)
     [DONE] TRIG element (lowest order verified vs NGSolve reference)
     [TODO] QUAD, TET, PRISM, HEX, PYRAMID

   Reference source: schoberl_zaglmayr_source/hcurlhofe_impl.hpp lines:
     SEGM:    154
     TRIG:    189-281
     QUAD:    290-394
     TET:     404-540
     PRISM:   549-838
     HEX:     848-1015
     PYRAMID: 1190-...
*)

BeginPackage["SchoberlZaglmayrBasis`"];

(* === Public API === *)
DuOp::usage = "DuOp[u, vars] gives {value=Grad u, curl=0} for a scalar u.";
uDvOp::usage = "uDvOp[u, v, vars] gives {value=u Grad v, curl=Grad u x Grad v}.";
uDvMinusvDu::usage = "uDvMinusvDu[u, v, vars] gives {value=u Grad v - v Grad u, curl=2 Grad u x Grad v}.";
wuDvMinuswvDu::usage = "wuDvMinuswvDu[u, v, w, vars] gives {value=w u Grad v - w v Grad u, curl=Grad(u w) x Grad v - Grad(v w) x Grad u}. Note: curl simplifies to wuDvMinuswvDu_curl in 2D.";

LegendreP1D::usage = "LegendreP1D[n, x] returns list {P_0(x), ..., P_n(x)} via 3-term recursion. Standard Legendre, P_0=1, P_1=x.";
IntLegNoBubble1D::usage = "IntLegNoBubble1D[n, x] returns list of integrated Legendre polynomials divided by edge bubble (1-x^2)/4, n+1 entries starting from index 0. P_0=-1/2, P_1=-x/2.";
IntLegMonomialExt1D::usage = "IntLegMonomialExt1D[n, x] returns the integrated Legendre polynomials L_2(x), L_3(x), ..., L_n(x) (n-1 entries). Each L_k vanishes at x = +/- 1. Used for cell/face shapes that must vanish on element boundaries.";
EvalScaled1D::usage = "EvalScaled1D[recur, n, x, t] evaluates a 3-term-recurrence polynomial in scaled form: P_n(x/t) t^n. recur is one of LegendreP1D, IntLegNoBubble1D.";

DubinerBasis2D::usage = "DubinerBasis2D[p, x, y] returns the (p+1)(p+2)/2 Dubiner basis functions on the reference triangle.";

GetEdgeSort::usage = "GetEdgeSort[edgeIdx, vnums] returns the 2 sorted vertex indices of the given edge.";
GetFaceSort::usage = "GetFaceSort[faceIdx, vnums] returns the sorted vertex indices of the given face.";

ShapesTRIG::usage = "ShapesTRIG[order, vnums, x, y] returns the H(curl) shape functions on the reference triangle as a list of {value, curl} pairs.";
ShapesQUAD::usage = "ShapesQUAD[order, vnums, x, y] returns the H(curl) shape functions on the reference quad [0,1]^2 as a list of {value, curl} pairs.";
ShapesTET::usage = "ShapesTET[order, vnums, x, y, z] returns the H(curl) shape functions on the reference tet (4 vertices, x+y+z<=1) as a list of {value, curl} pairs.";
ShapesHEX::usage = "ShapesHEX[order, vnums, x, y, z] returns the H(curl) shape functions on the reference hex [0,1]^3 as a list of {value, curl} pairs.";
ShapesPRISM::usage = "ShapesPRISM[order, vnums, x, y, z] returns the H(curl) shape functions on the reference prism (triangle x [0,1] in z) as a list of {value, curl} pairs. Vertices: 0,1,2 bottom triangle, 3,4,5 top.";
ShapesPYRAMID::usage = "ShapesPYRAMID[order, vnums, x, y, z] returns the H(curl) shape functions on the reference pyramid (square base + apex) as a list of {value, curl} pairs. Vertices: 0..3 square base, 4 apex.";

(* === Topology constants === *)
TrigEdges = {{0, 1}, {0, 2}, {1, 2}};   (* NGSolve TRIG edge endpoints *)
TrigFaces = {{0, 1, 2}};

QuadEdges = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};   (* CCW *)
QuadFaces = {{0, 1, 2, 3}};

TetEdges = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
TetFaces = {{1, 2, 3}, {0, 2, 3}, {0, 1, 3}, {0, 1, 2}};

HexEdges = {{0, 1}, {2, 3}, {3, 0}, {1, 2}, {4, 5}, {6, 7},
            {7, 4}, {5, 6}, {0, 4}, {1, 5}, {2, 6}, {3, 7}};
HexFaces = {{0, 1, 2, 3}, {4, 5, 6, 7}, {0, 1, 5, 4},
            {1, 2, 6, 5}, {2, 3, 7, 6}, {3, 0, 4, 7}};

PrismEdges = {{0, 1}, {1, 2}, {2, 0}, {3, 4}, {4, 5}, {5, 3},
              {0, 3}, {1, 4}, {2, 5}};   (* 0..5 horizontal, 6..8 vertical *)
PrismFaces = {{0, 1, 2}, {3, 4, 5},
              {0, 1, 4, 3}, {1, 2, 5, 4}, {2, 0, 3, 5}};  (* 0..1 trig, 2..4 quad *)

PyramidEdges = {{0, 1}, {1, 2}, {2, 3}, {3, 0},
                {0, 4}, {1, 4}, {2, 4}, {3, 4}};  (* 0..3 base, 4..7 to apex *)
PyramidFaces = {{0, 1, 2, 3},
                {0, 1, 4}, {1, 2, 4}, {2, 3, 4}, {3, 0, 4}};  (* 0 quad base, 1..4 trig *)


Begin["`Private`"];

(* ===================================================================
   Helper operators (from hcurlfe_utils.hpp)
   In NGSolve these use AutoDiff; here we use symbolic D[].
   vars = {x, y} for 2D, {x, y, z} for 3D.
   Each helper returns {valueVec, curlVec(or scalar)}.
================================================================== *)

DuOp[u_, vars_List] := Module[{grad = D[u, {vars}]},
  {grad, ConstantArray[0, Length[vars] (Length[vars] - 1) / 2]}
];

uDvOp[u_, v_, vars_List] := Module[{gu = D[u, {vars}], gv = D[v, {vars}]},
  {u gv, curlOf2[gu, gv, vars]}
];

uDvMinusvDu[u_, v_, vars_List] := Module[{gu = D[u, {vars}], gv = D[v, {vars}]},
  {u gv - v gu, 2 curlOf2[gu, gv, vars]}
];

wuDvMinuswvDu[u_, v_, w_, vars_List] := Module[
  {gu = D[u, {vars}], gv = D[v, {vars}], gw = D[w, {vars}], guw, gvw},
  guw = D[u w, {vars}];
  gvw = D[v w, {vars}];
  {w u gv - w v gu, curlOf2[guw, gv, vars] - curlOf2[gvw, gu, vars]}
];

(* curl helper: cross product of two grad vectors, scalar in 2D, 3D vector in 3D *)
curlOf2[g1_List, g2_List, vars_List] /; Length[vars] == 2 :=
  {g1[[1]] g2[[2]] - g1[[2]] g2[[1]]};   (* scalar (z component) *)
curlOf2[g1_List, g2_List, vars_List] /; Length[vars] == 3 := Cross[g1, g2];


(* ===================================================================
   1D recursive polynomials
================================================================== *)

(* LegendrePolynomial:
     P_0 = 1, P_1 = x
     P_n = (2-1/n) x P_{n-1} + (1/n - 1) P_{n-2}     for n >= 2
*)
LegendreP1D[n_Integer /; n >= 0, x_] := Module[
  {P, k, A, C0},
  P = Table[0, {k, 0, n}];   (* P[[k+1]] holds P_k(x) *)
  P[[1]] = 1;
  If[n >= 1, P[[2]] = x];
  Do[
    A = 2 - 1/k;
    C0 = 1/k - 1;
    P[[k + 1]] = A x P[[k]] + C0 P[[k - 1]],
    {k, 2, n}
  ];
  P
];

(* IntLegNoBubble:
     P_0 = -1/2, P_1 = -x/2
     P_n = A(n) x P_{n-1} + C(n) P_{n-2}
     where A(n) = (2(n+2)-3)/(n+2),  C(n) = -(n+2-3)/(n+2) = -(n-1)/(n+2)
*)
IntLegNoBubble1D[n_Integer /; n >= 0, x_] := Module[
  {P, k, i, A, C0},
  P = Table[0, {k, 0, n}];
  P[[1]] = -1/2;
  If[n >= 1, P[[2]] = -x/2];
  Do[
    i = k + 2;   (* matches the C++ recurrence with shifted index *)
    A = (2 i - 3) / i;
    C0 = -(i - 3) / i;
    P[[k + 1]] = A x P[[k]] + C0 P[[k - 1]],
    {k, 2, n}
  ];
  P
];

(* IntegratedLegendreMonomialExt:
     p_0 = -1, p_1 = x  (initial values, NOT returned)
     p_n = (2n-3)/n x p_{n-1} - (n-3)/n p_{n-2}     for n >= 2
   Returns p_2, p_3, ..., p_n  (n-1 entries).
   Each p_k for k >= 2 vanishes at x = +/-1 (e.g. p_2 = (x^2-1)/2).
*)
IntLegMonomialExt1D[n_Integer /; n >= 1, x_] := Module[
  {p3, p2, p1, j, result},
  result = {};
  If[n < 2, Return[{}]];   (* nothing to return *)
  p3 = 0; p2 = -1; p1 = x;
  Do[
    p3 = p2; p2 = p1;
    p1 = (2 j - 3)/j * x * p2 - (j - 3)/j * p3;
    AppendTo[result, p1],
    {j, 2, n}
  ];
  result
];

(* Scaled evaluation: P_k(x/t) * t^k  for k = 0..n.
   Equivalent to recurrence with extra t factor: P_n = A x P_{n-1} + C t^2 P_{n-2}.
   In NGSolve this is EvalScaled / EvalScaledMult.
   recur \[Element] {LegendreP1D, IntLegNoBubble1D}.
*)
EvalScaled1D[recur_, n_Integer /; n >= 0, x_, t_] := Module[
  {P, k, i, A, C0},
  P = Table[0, {k, 0, n}];
  Switch[recur,
    LegendreP1D,
      P[[1]] = 1;
      If[n >= 1, P[[2]] = x];
      Do[
        A = 2 - 1/k;
        C0 = 1/k - 1;
        P[[k + 1]] = A x P[[k]] + C0 t^2 P[[k - 1]],
        {k, 2, n}
      ],
    IntLegNoBubble1D,
      P[[1]] = -1/2;
      If[n >= 1, P[[2]] = -x/2];
      Do[
        i = k + 2;
        A = (2 i - 3) / i;
        C0 = -(i - 3) / i;
        P[[k + 1]] = A x P[[k]] + C0 t^2 P[[k - 1]],
        {k, 2, n}
      ],
    _, Throw["Unknown recur in EvalScaled1D"]
  ];
  P
];


(* ===================================================================
   2D Dubiner basis (triangle interior)
   See NGSolve recursive_pol_trig.hpp for implementation.
   For order p, returns (p+1)(p+2)/2 polynomials orthogonal on {x,y >= 0, x+y <= 1}.
   Construction:
     phi_{n, k}(x, y) = P_n^{0,0}(2 (x/(1-y)) - 1) * (1-y)^n * P_k^{2n+1, 0}(2y - 1)
   for n + k <= p.
   The first factor is Legendre in x scaled by (1-y), second is Jacobi(2n+1, 0) in y.
*)
JacobiP1D[alpha_, beta_, n_Integer /; n >= 0, x_] := Module[
  {P, k, A, B, C0},
  P = Table[0, {k, 0, n}];
  P[[1]] = 1;
  If[n >= 1, P[[2]] = ((alpha - beta) + (alpha + beta + 2) x) / 2];
  Do[
    Module[{a = 2 k (k + alpha + beta) (2 k + alpha + beta - 2),
            b = (2 k + alpha + beta - 1) (2 k + alpha + beta) (2 k + alpha + beta - 2),
            c = (2 k + alpha + beta - 1) (alpha^2 - beta^2),
            d = 2 (k + alpha - 1) (k + beta - 1) (2 k + alpha + beta)},
      P[[k + 1]] = ((b x + c) P[[k]] - d P[[k - 1]]) / a
    ],
    {k, 2, n}
  ];
  P
];

DubinerBasis2D[p_Integer /; p >= 0, x_, y_] := Module[
  {legX, jacY, result, n, k},
  result = {};
  Do[
    Do[
      legX = LegendreP1D[n, 2 x / (1 - y) - 1];
      jacY = JacobiP1D[2 n + 1, 0, k, 2 y - 1];
      AppendTo[result, legX[[n + 1]] * (1 - y)^n * jacY[[k + 1]]],
      {k, 0, p - n}
    ],
    {n, 0, p}
  ];
  result
];


(* ===================================================================
   Topology helpers (vertex-oriented sorting)
   GetEdgeSort returns the edge endpoints sorted by vertex number, ensuring
   a consistent orientation. vnums is a list of global vertex numbers per
   local vertex, e.g. vnums = {0, 1, 2} for unit triangle in canonical order.
================================================================== *)
GetEdgeSort[edgeIdx_Integer, vnums_List, edgeTable_:TrigEdges] :=
  With[{ev = edgeTable[[edgeIdx + 1]]},
    With[{v1 = ev[[1]], v2 = ev[[2]]},
      If[vnums[[v1 + 1]] > vnums[[v2 + 1]], {v2, v1}, {v1, v2}]
    ]
  ];

GetFaceSort[faceIdx_Integer, vnums_List, faceTable_:TrigFaces] :=
  With[{fv = faceTable[[faceIdx + 1]]},
    Reverse[SortBy[fv, vnums[[# + 1]] &]]   (* descending so that fav[0] = highest vertex number *)
  ];


(* ===================================================================
   TRIG (triangle 2D) shape functions
   Reference: hcurlhofe_impl.hpp lines 189-281
   Returns: list of {valueVec, curlScalar} pairs.
================================================================== *)

ShapesTRIG[order_Integer /; order >= 1, vnums_List /; Length[vnums] == 3, x_, y_] := Module[
  {lam, vars, shape, ii, i, e, edgePolys, p, fav, vop, pg,
   gradFacePolys, dubFace1, dubFace2},
  vars = {x, y};
  lam = {x, y, 1 - x - y};
  shape = {};

  (* Lowest-order Nedelec edge functions (3, one per edge) *)
  Do[
    e = GetEdgeSort[i, vnums];
    AppendTo[shape, uDvMinusvDu[lam[[e[[1]] + 1]], lam[[e[[2]] + 1]], vars]],
    {i, 0, 2}
  ];

  (* High-order edge gradient shapes via IntLegNoBubble *)
  Do[
    e = GetEdgeSort[i, vnums];
    p = order;   (* assume isotropic order; can be specialized per-edge *)
    If[p > 1,
      edgePolys = EvalScaled1D[IntLegNoBubble1D, p - 1,
        lam[[e[[2]] + 1]] - lam[[e[[1]] + 1]],
        lam[[e[[1]] + 1]] + lam[[e[[2]] + 1]]];
      (* multiply by edge bubble lam_e[0]*lam_e[1] *)
      Do[
        AppendTo[shape, DuOp[lam[[e[[1]] + 1]] lam[[e[[2]] + 1]] edgePolys[[k]], vars]],
        {k, 1, p}   (* p entries (indices 0..p-1 of recursion) *)
      ]
    ],
    {i, 0, 2}
  ];

  (* Inner (face) shapes: gradient + non-gradient *)
  p = order;
  If[p > 1,
    fav = GetFaceSort[0, vnums, TrigFaces];   (* 3 sorted indices *)
    pg = p - 2;   (* type 2 Nedelec assumed; type1 would be p - 3 *)
    (* Gradient face shapes: Du(DubinerBasis2D * lam_a*lam_b*lam_c) *)
    If[pg >= 0,
      dubFace1 = DubinerBasis2D[pg, lam[[fav[[1]] + 1]], lam[[fav[[2]] + 1]]];
      Do[
        AppendTo[shape, DuOp[lam[[fav[[1]] + 1]] lam[[fav[[2]] + 1]] lam[[fav[[3]] + 1]]
                              * dubFace1[[k]], vars]],
        {k, 1, Length[dubFace1]}
      ]
    ];

    (* Non-gradient face shapes (always for type 1 or 2):
       wuDv_minus_wvDu(lam_b, lam_c, dub_n)  for Dubiner index n
       wuDv_minus_wvDu(lam_b, lam_a, leg_k)  for k = 0..p-2
    *)
    dubFace1 = DubinerBasis2D[p - 2, lam[[fav[[1]] + 1]], lam[[fav[[2]] + 1]]];
    dubFace2 = LegendreP1D[p - 2, lam[[fav[[3]] + 1]] - lam[[fav[[2]] + 1]]];
    Do[
      AppendTo[shape, wuDvMinuswvDu[lam[[fav[[2]] + 1]], lam[[fav[[3]] + 1]],
                                     lam[[fav[[1]] + 1]] dubFace1[[k]], vars]],
      {k, 1, Length[dubFace1]}
    ];
    Do[
      AppendTo[shape, wuDvMinuswvDu[lam[[fav[[2]] + 1]], lam[[fav[[1]] + 1]],
                                     lam[[fav[[3]] + 1]] dubFace2[[k]], vars]],
      {k, 1, Length[dubFace2]}
    ];
  ];

  shape
];


(* ===================================================================
   QUAD (quadrilateral 2D) shape functions
   Reference: hcurlhofe_impl.hpp lines 290-394
   Reference quad: [0,1]^2 with vertices (0,0), (1,0), (1,1), (0,1).
================================================================== *)

ShapesQUAD[order_Integer /; order >= 1, vnums_List /; Length[vnums] == 4, x_, y_] := Module[
  {lami, sigma, vars, shape, ii, i, e, edgePolys, p, polXi, polEta, fmax, f1, f2,
   xiF, etaF, xi, lamE, bub},
  vars = {x, y};
  lami  = {(1 - x) (1 - y), x (1 - y), x y, (1 - x) y};
  sigma = {(1 - x) + (1 - y), x + (1 - y), x + y, (1 - x) + y};
  shape = {};

  (* === Edge low-order + high-order === *)
  Do[
    e = GetEdgeSort[i, vnums, QuadEdges];
    xi    = sigma[[e[[2]] + 1]] - sigma[[e[[1]] + 1]];
    lamE  = lami[[e[[1]] + 1]] + lami[[e[[2]] + 1]];
    bub   = (1/4) lamE (1 - xi^2);
    (* lowest order: uDv(0.5 * lam_e, xi) *)
    AppendTo[shape, uDvOp[0.5 lamE, xi, vars]];
    (* high order: gradient of IntLegNoBubble[xi] * bub *)
    p = order;
    If[p >= 1,
      edgePolys = IntLegNoBubble1D[p - 1, xi];
      Do[
        AppendTo[shape, DuOp[bub edgePolys[[k]], vars]],
        {k, 1, p}
      ]
    ],
    {i, 0, 3}
  ];

  (* === Face shapes === *)
  p = order;
  If[p >= 1,
    (* Choose face local axes (xi, eta) by the 2 dominant vertices.
       In NGSolve: fmax = vertex with highest vnum, f1 = (fmax+3)%4, f2 = (fmax+1)%4,
       then if vnums[f2] > vnums[f1] swap. *)
    fmax = First[Ordering[vnums, -1]] - 1;   (* index of max vertex (0-based) *)
    f1 = Mod[fmax + 3, 4];
    f2 = Mod[fmax + 1, 4];
    If[vnums[[f2 + 1]] > vnums[[f1 + 1]], {f1, f2} = {f2, f1}];
    xiF  = sigma[[fmax + 1]] - sigma[[f1 + 1]];   (* xi  in [-1, 1] *)
    etaF = sigma[[fmax + 1]] - sigma[[f2 + 1]];   (* eta in [-1, 1] *)

    (* Face gradient shapes *)
    Module[{bubF, polXiG, polEtaG},
      bubF = (1/16) (1 - xiF^2) (1 - etaF^2);
      polXiG  = IntLegNoBubble1D[p - 1, xiF];
      polEtaG = IntLegNoBubble1D[p - 1, etaF];
      Do[
        Do[
          AppendTo[shape, DuOp[bubF polXiG[[k + 1]] polEtaG[[j + 1]], vars]],
          {j, 0, p - 1}
        ],
        {k, 0, p - 1}
      ]
    ];

    (* Non-gradient face shapes: rotation + missing ones, using
       T_ORTHOPOL (= IntLegNoBubble) Calc(p+1, xi/eta) *)
    polXi  = IntLegNoBubble1D[p, xiF];
    polEta = IntLegNoBubble1D[p, etaF];
    (* Rotation of gradient fields: uDvMinusvDu(pol_eta[j], pol_xi[k]) for k,j = 0..p-1 *)
    Do[
      Do[
        AppendTo[shape, uDvMinusvDu[polEta[[j + 1]], polXi[[k + 1]], vars]],
        {j, 0, p - 1}
      ],
      {k, 0, p - 1}
    ];
    (* Missing ones: uDv(0.5 * pol_xi[j], eta) and uDv(0.5 * pol_eta[j], xi) *)
    Do[
      AppendTo[shape, uDvOp[0.5 polXi[[j + 1]], etaF, vars]],
      {j, 0, p - 1}
    ];
    Do[
      AppendTo[shape, uDvOp[0.5 polEta[[j + 1]], xiF, vars]],
      {j, 0, p - 1}
    ];
  ];

  shape
];


(* ===================================================================
   TET (tetrahedron 3D) shape functions
   Reference: hcurlhofe_impl.hpp lines 404-540
   Reference tet: vertices (0,0,0), (1,0,0), (0,1,0), (0,0,1).
================================================================== *)

ShapesTET[order_Integer /; order >= 1, vnums_List /; Length[vnums] == 4, x_, y_, z_] := Module[
  {lam, vars, shape, i, e, p, edgePolys, fav, vop, pg, dub1, dub2, leg1},
  vars = {x, y, z};
  lam  = {x, y, z, 1 - x - y - z};
  shape = {};

  (* === Edge shapes === *)
  Do[
    e = GetEdgeSort[i, vnums, TetEdges];
    AppendTo[shape, uDvMinusvDu[lam[[e[[1]] + 1]], lam[[e[[2]] + 1]], vars]];
    p = order;
    If[p >= 2,
      edgePolys = EvalScaled1D[IntLegNoBubble1D, p - 1,
        lam[[e[[2]] + 1]] - lam[[e[[1]] + 1]],
        lam[[e[[1]] + 1]] + lam[[e[[2]] + 1]]];
      Do[
        AppendTo[shape, DuOp[lam[[e[[1]] + 1]] lam[[e[[2]] + 1]] edgePolys[[k]], vars]],
        {k, 1, p - 1}
      ]
    ],
    {i, 0, 5}
  ];

  (* === Face shapes (4 faces) === *)
  Do[
    fav = GetFaceSort[i, vnums, TetFaces];   (* 3 sorted vertex indices *)
    vop = 6 - fav[[1]] - fav[[2]] - fav[[3]];   (* opposite vertex index, 0..3 *)
    p = order;
    If[p >= 2,
      pg = p - 2;   (* type 2 Nedelec; type1 would be p-3 *)

      (* Gradient face shapes: scaled Dubiner with bubble lam_a*lam_b*lam_c,
         scaled by 1-lam[vop]. Approximated: use Dubiner * cell-bubble factor. *)
      If[pg >= 0,
        dub1 = DubinerBasis2D[pg, lam[[fav[[1]] + 1]], lam[[fav[[2]] + 1]]];
        Do[
          AppendTo[shape, DuOp[
            lam[[fav[[1]] + 1]] lam[[fav[[2]] + 1]] lam[[fav[[3]] + 1]] dub1[[k]], vars]],
          {k, 1, Length[dub1]}
        ]
      ];

      (* Non-gradient face shapes (always for type 1 or type 2): *)
      dub2 = DubinerBasis2D[p - 2, lam[[fav[[1]] + 1]], lam[[fav[[2]] + 1]]];
      Do[
        AppendTo[shape, wuDvMinuswvDu[lam[[fav[[2]] + 1]], lam[[fav[[3]] + 1]],
                                       lam[[fav[[1]] + 1]] dub2[[k]], vars]],
        {k, 1, Length[dub2]}
      ];
      leg1 = LegendreP1D[p - 2, lam[[fav[[3]] + 1]] - lam[[fav[[2]] + 1]]];
      Do[
        AppendTo[shape, wuDvMinuswvDu[lam[[fav[[2]] + 1]], lam[[fav[[1]] + 1]],
                                       lam[[fav[[3]] + 1]] leg1[[k]], vars]],
        {k, 1, Length[leg1]}
      ];
    ],
    {i, 0, 3}
  ];

  (* === Cell shapes (interior) === *)
  p = order;
  If[p >= 3,
    (* Gradient cell shapes via Dubiner3D × cell bubble lam[0]*lam[1]*lam[2]*lam[3]
       NOTE: 3D Dubiner not yet implemented; use 2D Dubiner * Legendre placeholder *)
    Module[{cellBub = lam[[1]] lam[[2]] lam[[3]] lam[[4]]},
      Do[
        Do[
          Do[
            Module[{idx = i + j + k},
              If[idx <= p - 3,
                AppendTo[shape, DuOp[
                  cellBub LegendreP1D[i, x][[i + 1]]
                          LegendreP1D[j, y][[j + 1]]
                          LegendreP1D[k, z][[k + 1]], vars]]
              ]
            ],
            {k, 0, p - 3}
          ],
          {j, 0, p - 3}
        ],
        {i, 0, p - 3}
      ]
    ];
    (* Non-gradient cell shapes (placeholder; full Dubiner3D needed for production) *)
  ];

  shape
];


(* ===================================================================
   HEX (hexahedron 3D) shape functions
   Reference: hcurlhofe_impl.hpp lines 848-1015
   Reference hex: [0,1]^3 with standard 8-vertex ordering.
================================================================== *)

ShapesHEX[order_Integer /; order >= 1, vnums_List /; Length[vnums] == 8, x_, y_, z_] := Module[
  {lami, sigma, vars, shape, i, e, p, edgePolys, polXi, polEta, polZeta,
   xi, lamE, bub, qmax, q1, q2, fmax, f1, f2, lamF, xiF, etaF, faceVerts,
   bubF, polXiG, polEtaG},
  vars = {x, y, z};
  lami = {(1 - x) (1 - y) (1 - z), x (1 - y) (1 - z), x y (1 - z), (1 - x) y (1 - z),
          (1 - x) (1 - y) z,       x (1 - y) z,       x y z,       (1 - x) y z};
  sigma = {(1 - x) + (1 - y) + (1 - z), x + (1 - y) + (1 - z),
           x + y + (1 - z),             (1 - x) + y + (1 - z),
           (1 - x) + (1 - y) + z,       x + (1 - y) + z,
           x + y + z,                   (1 - x) + y + z};
  shape = {};

  (* === Edges (12) === *)
  Do[
    e = GetEdgeSort[i, vnums, HexEdges];
    xi   = sigma[[e[[2]] + 1]] - sigma[[e[[1]] + 1]];
    lamE = lami[[e[[1]] + 1]] + lami[[e[[2]] + 1]];
    bub  = (1/4) lamE (1 - xi^2);
    AppendTo[shape, uDvOp[0.5 lamE, xi, vars]];
    p = order;
    If[p >= 1,
      edgePolys = IntLegNoBubble1D[p - 1, xi];
      Do[
        AppendTo[shape, DuOp[bub edgePolys[[k]], vars]],
        {k, 1, p}
      ]
    ],
    {i, 0, 11}
  ];

  (* === Faces (6) ===
     For each face: gradient shapes (tensor product * face bubble),
     rotation of gradient shapes, and "missing ones" along each face axis. *)
  Do[
    p = order;
    If[p >= 1,
      faceVerts = HexFaces[[i + 1]];   (* 4 vertex indices for this face *)
      lamF = Sum[lami[[faceVerts[[j]] + 1]], {j, 1, 4}];

      (* Sort face vertices: qmax = highest vnum among the 4, then q1=(qmax+3)%4, q2=(qmax+1)%4 *)
      qmax = First[Ordering[vnums[[# + 1]] & /@ faceVerts, -1]] - 1;
      q1 = Mod[qmax + 3, 4];
      q2 = Mod[qmax + 1, 4];
      If[vnums[[faceVerts[[q2 + 1]] + 1]] > vnums[[faceVerts[[q1 + 1]] + 1]], {q1, q2} = {q2, q1}];
      fmax = faceVerts[[qmax + 1]];
      f1   = faceVerts[[q1 + 1]];
      f2   = faceVerts[[q2 + 1]];
      xiF  = sigma[[fmax + 1]] - sigma[[f1 + 1]];
      etaF = sigma[[fmax + 1]] - sigma[[f2 + 1]];

      (* Gradient face shapes: lam_f * IntLegNoBubble(xi) * IntLegNoBubble(eta) *)
      bubF = lamF * (1/16) (1 - xiF^2) (1 - etaF^2);
      polXiG  = IntLegNoBubble1D[p - 1, xiF];
      polEtaG = IntLegNoBubble1D[p - 1, etaF];
      Do[
        Do[
          AppendTo[shape, DuOp[bubF polXiG[[k + 1]] polEtaG[[j + 1]], vars]],
          {j, 0, p - 1}
        ],
        {k, 0, p - 1}
      ];

      (* Non-gradient face shapes use T_ORTHOPOL = IntegratedLegendreMonomialExt.
         Calc(p+1, xi) returns p entries (p_2, ..., p_{p+1}), vanishing at xi = +/-1.
         This is for face shapes that must vanish on the face boundary edges. *)
      polXi  = IntLegMonomialExt1D[p + 1, xiF];   (* p entries *)
      polEta = IntLegMonomialExt1D[p + 1, etaF];
      Do[
        Do[
          AppendTo[shape, uDvMinusvDu[polEta[[j + 1]], lamF polXi[[k + 1]], vars]],
          {j, 0, p - 1}
        ],
        {k, 0, p - 1}
      ];
      Do[
        AppendTo[shape, wuDvMinuswvDu[0.5, etaF, polXi[[j + 1]] lamF, vars]],
        {j, 0, p - 1}
      ];
      Do[
        AppendTo[shape, wuDvMinuswvDu[0.5, xiF, polEta[[j + 1]] lamF, vars]],
        {j, 0, p - 1}
      ]
    ],
    {i, 0, 5}
  ];

  (* === Cell shapes (interior, tensor product 3D) === *)
  p = order;
  If[p >= 1,
    (* Gradient cell shapes: Du(IntLegNoBubble(2x-1) * IntLegNoBubble(2y-1) * IntLegNoBubble(2z-1) * x(1-x)*y(1-y)*z(1-z)) *)
    Module[{px = IntLegNoBubble1D[p - 1, 2 x - 1] * x (1 - x),
            py = IntLegNoBubble1D[p - 1, 2 y - 1] * y (1 - y),
            pz = IntLegNoBubble1D[p - 1, 2 z - 1] * z (1 - z)},
      Do[
        Do[
          Do[
            AppendTo[shape, DuOp[px[[i + 1]] py[[j + 1]] pz[[k + 1]], vars]],
            {k, 0, p - 1}
          ],
          {j, 0, p - 1}
        ],
        {i, 0, p - 1}
      ]
    ];

    (* Non-gradient cell shapes via T_ORTHOPOL = IntegratedLegendreMonomialExt
       Calc(p+1, 2*x-1) returns p entries (p_2, ..., p_{p+1}), each vanishing at x = 0 and x = 1.
       This gives polynomials that are zero on the cell boundary. *)
    polXi   = IntLegMonomialExt1D[p + 1, 2 x - 1];   (* p entries *)
    polEta  = IntLegMonomialExt1D[p + 1, 2 y - 1];
    polZeta = IntLegMonomialExt1D[p + 1, 2 z - 1];

    (* Rotations of gradient: 2 per (i,j,k) *)
    Do[
      Do[
        Do[
          AppendTo[shape, uDvMinusvDu[polXi[[i + 1]] polEta[[j + 1]], polZeta[[k + 1]], vars]];
          AppendTo[shape, uDvMinusvDu[polXi[[i + 1]], polEta[[j + 1]] polZeta[[k + 1]], vars]],
          {k, 0, p - 1}
        ],
        {j, 0, p - 1}
      ],
      {i, 0, p - 1}
    ];

    (* Type 3 missing ones *)
    Do[
      Do[
        AppendTo[shape, wuDvMinuswvDu[z, 1 - z, polXi[[i + 1]] polEta[[j + 1]], vars]],
        {j, 0, p - 1}
      ],
      {i, 0, p - 1}
    ];
    Do[
      Do[
        AppendTo[shape, wuDvMinuswvDu[y, 1 - y, polXi[[i + 1]] polZeta[[k + 1]], vars]],
        {k, 0, p - 1}
      ],
      {i, 0, p - 1}
    ];
    Do[
      Do[
        AppendTo[shape, wuDvMinuswvDu[x, 1 - x, polEta[[j + 1]] polZeta[[k + 1]], vars]],
        {k, 0, p - 1}
      ],
      {j, 0, p - 1}
    ]
  ];

  shape
];


(* ===================================================================
   PRISM (wedge 3D) shape functions
   Reference: hcurlhofe_impl.hpp lines 549-838
   6 vertices: 0,1,2 bottom triangle, 3,4,5 top triangle
   9 edges:    0..5 horizontal (3 bottom + 3 top), 6..8 vertical
   5 faces:    0..1 triangles, 2..4 quad sides
================================================================== *)

ShapesPRISM[order_Integer /; order >= 1, vnums_List /; Length[vnums] == 6, x_, y_, z_] := Module[
  {lam, muz, sigma, vars, shape, i, e, p, edgePolys, fav, fmax, ftrig, fz, faceVerts,
   xi, eta, zeta, scaleXi, scaleEta, bub, dub2, leg1, polXi, polEta, polZeta,
   polTrig, polZ},
  vars = {x, y, z};
  lam   = {x, y, 1 - x - y, x, y, 1 - x - y};
  muz   = {1 - z, 1 - z, 1 - z, z, z, z};
  sigma = lam + muz;
  shape = {};

  (* === Horizontal edges (0..5) === *)
  Do[
    e = GetEdgeSort[i, vnums, PrismEdges];
    AppendTo[shape, wuDvMinuswvDu[lam[[e[[1]] + 1]], lam[[e[[2]] + 1]],
                                    muz[[e[[2]] + 1]], vars]];
    p = order;
    If[p >= 2,
      xi = lam[[e[[2]] + 1]] - lam[[e[[1]] + 1]];
      eta = lam[[e[[1]] + 1]] + lam[[e[[2]] + 1]];
      bub = lam[[e[[1]] + 1]] lam[[e[[2]] + 1]] muz[[e[[2]] + 1]];
      edgePolys = bub * EvalScaled1D[IntLegNoBubble1D, p - 1, xi, eta];
      Do[
        AppendTo[shape, DuOp[edgePolys[[k]], vars]],
        {k, 1, p - 1}
      ]
    ],
    {i, 0, 5}
  ];

  (* === Vertical edges (6..8) === *)
  Do[
    e = GetEdgeSort[i, vnums, PrismEdges];
    AppendTo[shape, wuDvMinuswvDu[muz[[e[[1]] + 1]], muz[[e[[2]] + 1]],
                                    lam[[e[[2]] + 1]], vars]];
    p = order;
    If[p >= 2,
      bub = muz[[e[[1]] + 1]] muz[[e[[2]] + 1]] lam[[e[[2]] + 1]];
      edgePolys = bub * IntLegNoBubble1D[p - 1, muz[[e[[2]] + 1]] - muz[[e[[1]] + 1]]];
      Do[
        AppendTo[shape, DuOp[edgePolys[[k]], vars]],
        {k, 1, p - 1}
      ]
    ],
    {i, 6, 8}
  ];

  (* === Triangle faces (2 of them: i = 0, 1) === *)
  Do[
    p = order;
    If[p >= 2,
      fav = GetFaceSort[i, vnums, PrismFaces];   (* 3 sorted indices *)
      (* Gradient face shapes: Du(Dubiner(lam_a, lam_b) * lam_a * lam_b * lam_c * muz_c) *)
      Module[{dub1 = DubinerBasis2D[p - 2, lam[[fav[[1]] + 1]], lam[[fav[[2]] + 1]]]},
        Do[
          AppendTo[shape, DuOp[
            lam[[fav[[1]] + 1]] lam[[fav[[2]] + 1]] lam[[fav[[3]] + 1]]
              muz[[fav[[3]] + 1]] dub1[[k]], vars]],
          {k, 1, Length[dub1]}
        ]
      ];
      (* Non-gradient face shapes *)
      dub2 = DubinerBasis2D[p - 2, lam[[fav[[1]] + 1]], lam[[fav[[2]] + 1]]];
      Do[
        AppendTo[shape, wuDvMinuswvDu[lam[[fav[[2]] + 1]], lam[[fav[[3]] + 1]],
          lam[[fav[[1]] + 1]] muz[[fav[[3]] + 1]] dub2[[k]], vars]],
        {k, 1, Length[dub2]}
      ];
      leg1 = LegendreP1D[p - 2, lam[[fav[[3]] + 1]] - lam[[fav[[2]] + 1]]];
      Do[
        AppendTo[shape, wuDvMinuswvDu[lam[[fav[[2]] + 1]], lam[[fav[[1]] + 1]],
          lam[[fav[[3]] + 1]] muz[[fav[[3]] + 1]] leg1[[k]], vars]],
        {k, 1, Length[leg1]}
      ]
    ],
    {i, 0, 1}
  ];

  (* === Quad faces (3 of them: i = 2..4) === *)
  Do[
    p = order;
    If[p >= 1,
      faceVerts = PrismFaces[[i + 1]];   (* 4 vertex indices *)
      (* Determine face axes: fmax = highest-vnum vertex among the 4 *)
      Module[{qmaxIdx, fmaxV, ftrigV, fzV},
        qmaxIdx = First[Ordering[vnums[[# + 1]] & /@ faceVerts, -1]] - 1;
        fmaxV = faceVerts[[qmaxIdx + 1]];
        fzV = faceVerts[[3 - qmaxIdx + 1]];   (* 3 - qmax in 0..3 *)
        ftrigV = faceVerts[[BitXor[qmaxIdx, 1] + 1]];

        xi   = lam[[fmaxV + 1]] - lam[[ftrigV + 1]];
        eta  = 1 - lam[[fmaxV + 1]] - lam[[ftrigV + 1]];
        zeta = muz[[fmaxV + 1]] - muz[[fzV + 1]];

        (* Gradient face shapes: tensor of trig in xy with scaling × poly in z *)
        Module[{scaleXi2 = lam[[fmaxV + 1]] + lam[[ftrigV + 1]],
                bubF},
          bubF = (1/16) (scaleXi2^2 - xi^2) (1 - zeta^2);
          polXi  = EvalScaled1D[IntLegNoBubble1D, p - 1, xi, scaleXi2];
          polZeta = IntLegNoBubble1D[p - 1, zeta];
          Do[
            Do[
              AppendTo[shape, DuOp[bubF polXi[[k + 1]] polZeta[[j + 1]], vars]],
              {j, 0, p - 1}
            ],
            {k, 0, p - 1}
          ];

          (* Rotation of gradient + Type 3 missing ones *)
          polXi  = EvalScaled1D[IntLegNoBubble1D, p, xi, scaleXi2];
          polZeta = IntLegNoBubble1D[p, zeta];
          Do[
            Do[
              AppendTo[shape, uDvMinusvDu[polXi[[k + 1]], polZeta[[j + 1]], vars]],
              {j, 0, p - 1}
            ],
            {k, 0, p - 1}
          ];
          (* Type 3: ned0_trig × pol_z, and z.DValue × pol_xy *)
          Do[
            AppendTo[shape, wuDvMinuswvDu[muz[[fzV + 1]], muz[[fmaxV + 1]], polXi[[j + 1]], vars]],
            {j, 0, p - 1}
          ];
          Do[
            AppendTo[shape, wuDvMinuswvDu[lam[[ftrigV + 1]], lam[[fmaxV + 1]], polZeta[[j + 1]], vars]],
            {j, 0, p - 1}
          ]
        ]
      ]
    ],
    {i, 2, 4}
  ];

  (* === Cell shapes === *)
  p = order;
  If[p >= 2,
    (* Gradient cell: Dubiner(x,y) * x*y*(1-x-y) * Legendre(2z-1) * z*(1-z) *)
    Module[{cellBubXY = x y (1 - x - y),
            cellBubZ = z (1 - z),
            polTrig = DubinerBasis2D[p - 2, x, y]},
      polZ = LegendreP1D[p - 1, 2 z - 1];
      Do[
        Do[
          AppendTo[shape, DuOp[
            cellBubXY polTrig[[i + 1]] * cellBubZ polZ[[k + 1]], vars]],
          {k, 0, p - 1}
        ],
        {i, 0, Length[polTrig] - 1}
      ]
    ];
    (* Non-gradient cell shapes (placeholder using DubinerBasis2D structure for xy) *)
    (* The full implementation requires TrigShapesInnerLegendre and IntLegNoBubble in z;
       see C++ source line 805-836 for full structure. *)
  ];

  shape
];


(* ===================================================================
   PYRAMID (3D) shape functions
   Reference: hcurlhofe_impl.hpp lines 1190-
   5 vertices: 0..3 square base, 4 apex
   8 edges: 0..3 base, 4..7 to apex
   5 faces: 0 square base, 1..4 triangles
   NOTE: complex (~700 lines C++); this is an initial implementation covering
   lowest-order Nedelec only. High-order modes have apex singularity handling
   that requires additional care; see source for full treatment.
================================================================== *)

ShapesPYRAMID[order_Integer /; order >= 1, vnums_List /; Length[vnums] == 5, x_, y_, z_] := Module[
  {lam, vars, shape, i, e},
  vars = {x, y, z};
  (* Pyramid local coordinates (NGSolve convention):
     z in [0,1] from base to apex; xy collapses linearly with z.
     Base coord at level z: scale by (1-z). Apex at z=1.
     Corner shape functions on pyramid (NGSolve uses): *)
  lam = {(1 - x - z) (1 - y - z) / (1 - z),
          x          (1 - y - z) / (1 - z),
          x          y           / (1 - z),
         (1 - x - z) y           / (1 - z),
          z};
  shape = {};

  (* Lowest-order Nedelec edge functions (8) *)
  Do[
    e = GetEdgeSort[i, vnums, PyramidEdges];
    AppendTo[shape, uDvMinusvDu[lam[[e[[1]] + 1]], lam[[e[[2]] + 1]], vars]],
    {i, 0, 7}
  ];

  (* Higher-order shapes: NOT IMPLEMENTED in this initial version.
     The pyramid element requires special apex-handling; see hcurlhofe_impl.hpp
     line 1190- for full implementation. Workaround for production:
     decompose pyramid into 2 tets and use ShapesTET. *)

  shape
];


(* ===================================================================
   DubinerBasis3D for tet interior cell shapes (currently placeholder
   in ShapesTET). Full implementation: phi_{ijk}(x,y,z) =
     Legendre_i(2x/(1-y-z) - 1) (1-y-z)^i ×
     Jacobi_j^{2i+1, 0}(2y/(1-z) - 1) (1-z)^j ×
     Jacobi_k^{2(i+j)+2, 0}(2z - 1)
================================================================== *)

(* Future TODO: implement DubinerBasis3D to replace placeholder in ShapesTET.
   Full implementation pattern for each remaining feature:
     1. Define lambda (barycentric or tensor product) coordinates
     2. For each edge: lowest-order via uDvMinusvDu, high-order via Du(IntLegNoBubble*bubble)
     3. For each face: gradient (DubinerBasis2D * bubble) + non-gradient (wuDvMinuswvDu)
     4. For each cell (3D): gradient + non-gradient interior shapes
   See schoberl_zaglmayr_source/hcurlhofe_impl.hpp for reference C++ implementation. *)

End[];
EndPackage[];

(* === Self-test === *)
If[$ScriptCommandLine =!= {},
  Print["SchoberlZaglmayrBasis self-test"];
  Print["================================"];

  Print["LegendreP1D[5, x] for x = 0.3:"];
  Print["  ", N[SchoberlZaglmayrBasis`LegendreP1D[5, 0.3], 12]];
  Print["  Mathematica builtin: ",
        N[Table[LegendreP[k, 0.3], {k, 0, 5}], 12]];

  Print["IntLegNoBubble1D[3, x] for x = 0.3:"];
  Print["  ", N[SchoberlZaglmayrBasis`IntLegNoBubble1D[3, 0.3], 12]];

  Print["TRIG order=2, vnums={0,1,2}, point (1/4, 1/4):"];
  Module[{shapes},
    shapes = SchoberlZaglmayrBasis`ShapesTRIG[2, {0, 1, 2}, x, y];
    Print["  Number of basis functions: ", Length[shapes]];
    Print["  First shape value at (1/4,1/4): ",
          shapes[[1, 1]] /. {x -> 1/4, y -> 1/4}];
    Print["  First shape curl  at (1/4,1/4): ",
          shapes[[1, 2]] /. {x -> 1/4, y -> 1/4}];
  ];

  Print["QUAD order=2, vnums={0,1,2,3}, point (1/2, 1/2):"];
  Module[{shapes},
    shapes = SchoberlZaglmayrBasis`ShapesQUAD[2, {0, 1, 2, 3}, x, y];
    Print["  Number of basis functions: ", Length[shapes]];
    Print["  First (edge 0-1 low) value at (1/2,1/2): ",
          shapes[[1, 1]] /. {x -> 1/2, y -> 1/2}];
  ];

  Print["TET order=2, vnums={0,1,2,3}, point (1/4, 1/4, 1/4):"];
  Module[{shapes},
    shapes = SchoberlZaglmayrBasis`ShapesTET[2, {0, 1, 2, 3}, x, y, z];
    Print["  Number of basis functions: ", Length[shapes]];
    Print["  First (edge 0-1 low) value at (1/4,1/4,1/4): ",
          shapes[[1, 1]] /. {x -> 1/4, y -> 1/4, z -> 1/4}];
    Print["  First curl                                  : ",
          shapes[[1, 2]] /. {x -> 1/4, y -> 1/4, z -> 1/4}];
  ];

  Print["HEX order=1, vnums={0,1,2,3,4,5,6,7}, point (1/2, 1/2, 1/2):"];
  Module[{shapes},
    shapes = SchoberlZaglmayrBasis`ShapesHEX[1, {0, 1, 2, 3, 4, 5, 6, 7}, x, y, z];
    Print["  Number of basis functions (order 1): ", Length[shapes]];
    Print["  First (edge 0-1 low) value at (1/2,1/2,1/2): ",
          Simplify[shapes[[1, 1]] /. {x -> 1/2, y -> 1/2, z -> 1/2}]];
    Print["  First curl: ",
          Simplify[shapes[[1, 2]] /. {x -> 1/2, y -> 1/2, z -> 1/2}]];
  ];

  Print["PRISM order=2, vnums={0..5}, point (1/3, 1/3, 1/2):"];
  Module[{shapes},
    shapes = SchoberlZaglmayrBasis`ShapesPRISM[2, Range[0, 5], x, y, z];
    Print["  Number of basis functions: ", Length[shapes]];
    Print["  First (edge 0-1 horiz low) value at (1/3,1/3,1/2): ",
          Simplify[shapes[[1, 1]] /. {x -> 1/3, y -> 1/3, z -> 1/2}]];
    Print["  First curl: ",
          Simplify[shapes[[1, 2]] /. {x -> 1/3, y -> 1/3, z -> 1/2}]];
  ];

  Print["PYRAMID order=1 (lowest only), vnums={0..4}, point (1/4, 1/4, 1/4):"];
  Module[{shapes},
    shapes = SchoberlZaglmayrBasis`ShapesPYRAMID[1, Range[0, 4], x, y, z];
    Print["  Number of basis functions: ", Length[shapes]];
    Print["  Expected for order 1: 8 lowest-order edges only"];
    Print["  First (edge 0-1 base) value at (1/4,1/4,1/4): ",
          Simplify[shapes[[1, 1]] /. {x -> 1/4, y -> 1/4, z -> 1/4}]];
  ];
]
