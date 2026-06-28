(* ::Package:: *)

(* ::Title:: *)
(*RadiaBasis -- finite-element basis function reference for Radia*)


(* ::Text:: *)
(*Phase 1 (radia-mcp 0.38.0) -- Triangle + Tetrahedron, H1 Lagrange + RWG (HDiv RT0) + L2.*)
(*Mathematica is the canonical notation; SymPy / NumPy translations live in*)
(*  - tests/basis/test_basis_functions.py  (CI-verified Python ports)*)
(*  - radia.bem.efie_rwg                    (production NumPy for RWG)*)
(*  - NGSolve internal                      (production C++ for hierarchical H1)*)
(**)
(*Loading:*)
(*  Get["S:/Radia/01_GitHub/packages/radia-mcp/src/radia_mcp/mathematica/basis_functions/RadiaBasis.m"]*)
(**)
(*Reference triangle:  vertices (0,0), (1,0), (0,1).*)
(*Reference tet:       vertices (0,0,0), (1,0,0), (0,1,0), (0,0,1).*)
(*Barycentric: lam0 = 1 - xi - eta [- zeta], lam1..3 = xi/eta/zeta.*)


BeginPackage["RadiaBasis`"];

(* ----- Public symbols ----- *)
TriH1P1::usage  = "TriH1P1[xi, eta] -- Triangle H1 Lagrange P1 (3 vertex shapes).";
TriH1P2::usage  = "TriH1P2[xi, eta] -- Triangle H1 Lagrange P2 (3 vertex + 3 edge midpoints).";
TriH1P3::usage  = "TriH1P3[xi, eta] -- Triangle H1 Lagrange P3 (10 dofs).";
TriH1P4::usage  = "TriH1P4[xi, eta] -- Triangle H1 Lagrange P4 (15 dofs).";
TriH1P5::usage  = "TriH1P5[xi, eta] -- Triangle H1 Lagrange P5 (21 dofs).";
TriH1Pk::usage  = "TriH1Pk[k][xi, eta] -- Triangle H1 Lagrange P_k for arbitrary k>=1.";

TriRT0::usage   = "TriRT0[xi, eta, v, sigma] -- Triangle RWG (HDiv RT0).  v={v0,v1,v2}, sigma={s0,s1,s2}.";
TriRT0Div::usage = "TriRT0Div[v, sigma] -- divergence of Triangle RT0 (constant per triangle, returns {sigma_k / A_T}).";

TriL2P0::usage  = "TriL2P0[xi, eta] -- Triangle L2 P0 (constant).";
TriL2P1::usage  = "TriL2P1[xi, eta] -- Triangle L2 P1 (centred linear, 3 dofs).";
TriDubiner::usage = "TriDubiner[i, j, xi, eta] -- Dubiner orthogonal basis on reference triangle, indices i+j<=k.";

TetH1P1::usage  = "TetH1P1[xi, eta, zeta] -- Tetrahedron H1 Lagrange P1 (4 vertex shapes).";
TetH1P2::usage  = "TetH1P2[xi, eta, zeta] -- Tetrahedron H1 Lagrange P2 (10 dofs).";
TetH1P3::usage  = "TetH1P3[xi, eta, zeta] -- Tetrahedron H1 Lagrange P3 (20 dofs).";

TetRT0::usage   = "TetRT0[xi, eta, zeta, v, sigma] -- Tet HDiv RT0 (4 face-based basis).";
TetRT0Div::usage = "TetRT0Div[v, sigma] -- divergence of Tet RT0 (constant {sigma_k / V_T}).";

(* Verification helpers *)
TriPartitionOfUnity::usage = "TriPartitionOfUnity[basis] -- assert Sum[basis_i, i] == 1 symbolically.";
TriDOFCount::usage = "TriDOFCount[basis] -- length of basis list.";


Begin["`Private`"];

(* ----- Triangle H1 Lagrange ----- *)

TriH1P1[xi_, eta_] := Module[{lam0, lam1, lam2},
  lam0 = 1 - xi - eta; lam1 = xi; lam2 = eta;
  {lam0, lam1, lam2}
];

TriH1P2[xi_, eta_] := Module[{lam0, lam1, lam2},
  lam0 = 1 - xi - eta; lam1 = xi; lam2 = eta;
  {
    lam0 (2 lam0 - 1),     (* v0 *)
    lam1 (2 lam1 - 1),     (* v1 *)
    lam2 (2 lam2 - 1),     (* v2 *)
    4 lam0 lam1,            (* m01 *)
    4 lam1 lam2,            (* m12 *)
    4 lam2 lam0             (* m20 *)
  }
];

TriH1P3[xi_, eta_] := Module[{lam0, lam1, lam2},
  lam0 = 1 - xi - eta; lam1 = xi; lam2 = eta;
  {
    (1/2) lam0 (3 lam0 - 1) (3 lam0 - 2),
    (1/2) lam1 (3 lam1 - 1) (3 lam1 - 2),
    (1/2) lam2 (3 lam2 - 1) (3 lam2 - 2),
    (9/2) lam0 lam1 (3 lam0 - 1),
    (9/2) lam0 lam1 (3 lam1 - 1),
    (9/2) lam1 lam2 (3 lam1 - 1),
    (9/2) lam1 lam2 (3 lam2 - 1),
    (9/2) lam2 lam0 (3 lam2 - 1),
    (9/2) lam2 lam0 (3 lam0 - 1),
    27 lam0 lam1 lam2
  }
];

(* TriH1Pk: generic Lagrange P_k on triangle.  Nodal points are
   barycentric multi-indices (i,j,m) with i+j+m = k.  The basis
   function for node (i,j,m) is:
     phi_{i,j,m}(lam0,lam1,lam2) =
        L_i(k, lam0) * L_j(k, lam1) * L_m(k, lam2)
   where L_n(k, x) = Product_{r=0..n-1} ((k x - r) / (n - r))    (Lagrange on equispaced grid 0, 1/k, .., 1).      *)
LagrangeShape[k_, n_, x_] := If[n == 0, 1,
  Product[(k x - r) / (n - r), {r, 0, n - 1}]];

TriH1Pk[k_Integer?Positive] := Module[{nodes},
  Function[{xi, eta},
    Module[{lam0, lam1, lam2, idx},
      lam0 = 1 - xi - eta; lam1 = xi; lam2 = eta;
      idx = Flatten[Table[{i, j, k - i - j}, {i, 0, k}, {j, 0, k - i}], 1];
      Table[
        LagrangeShape[k, ind[[1]], lam0] *
        LagrangeShape[k, ind[[2]], lam1] *
        LagrangeShape[k, ind[[3]], lam2],
        {ind, idx}]
    ]
  ]
];

TriH1P4[xi_, eta_] := TriH1Pk[4][xi, eta];
TriH1P5[xi_, eta_] := TriH1Pk[5][xi, eta];


(* ----- Triangle HDiv RT0 (= RWG) ----- *)

TriRT0[xi_, eta_, v_List, sigma_List] := Module[
  {r, vopp, area},
  r = (1 - xi - eta) v[[1]] + xi v[[2]] + eta v[[3]];
  vopp = {v[[1]], v[[2]], v[[3]]};
  area = (1/2) Norm[Cross[v[[2]] - v[[1]], v[[3]] - v[[1]]]];
  Table[ sigma[[k]] (r - vopp[[k]]) / (2 area), {k, 1, 3}]
];

TriRT0Div[v_List, sigma_List] := Module[{area},
  area = (1/2) Norm[Cross[v[[2]] - v[[1]], v[[3]] - v[[1]]]];
  sigma / area
];


(* ----- Triangle L2 ----- *)

TriL2P0[xi_, eta_] := {1};

TriL2P1[xi_, eta_] := {1, xi - 1/3, eta - 1/3};

TriDubiner[i_Integer?NonNegative, j_Integer?NonNegative, xi_, eta_] :=
  Module[{a, b, P1, P2},
    a = If[eta == 1, 0, 2 xi / (1 - eta) - 1];
    b = 2 eta - 1;
    P1 = JacobiP[i, 0, 0, a];
    P2 = JacobiP[j, 2 i + 1, 0, b];
    P1 P2 (1 - eta)^i
];


(* ----- Tetrahedron H1 Lagrange ----- *)

TetH1P1[xi_, eta_, zeta_] := Module[{lam0, lam1, lam2, lam3},
  lam0 = 1 - xi - eta - zeta; lam1 = xi; lam2 = eta; lam3 = zeta;
  {lam0, lam1, lam2, lam3}
];

TetH1P2[xi_, eta_, zeta_] := Module[{lam0, lam1, lam2, lam3},
  lam0 = 1 - xi - eta - zeta; lam1 = xi; lam2 = eta; lam3 = zeta;
  {
    lam0 (2 lam0 - 1), lam1 (2 lam1 - 1),
    lam2 (2 lam2 - 1), lam3 (2 lam3 - 1),
    4 lam0 lam1, 4 lam0 lam2, 4 lam0 lam3,
    4 lam1 lam2, 4 lam1 lam3, 4 lam2 lam3
  }
];

(* TetH1P3 -- 20 dofs.  Node order: v0, v1, v2, v3,
   e01a, e01b, e02a, e02b, e03a, e03b, e12a, e12b, e13a, e13b, e23a, e23b,
   f012, f013, f023, f123                                              *)
TetH1P3[xi_, eta_, zeta_] := Module[{lam0, lam1, lam2, lam3, vshape},
  lam0 = 1 - xi - eta - zeta; lam1 = xi; lam2 = eta; lam3 = zeta;
  vshape[lam_] := (1/2) lam (3 lam - 1) (3 lam - 2);
  {
    vshape[lam0], vshape[lam1], vshape[lam2], vshape[lam3],
    (9/2) lam0 lam1 (3 lam0 - 1), (9/2) lam0 lam1 (3 lam1 - 1),
    (9/2) lam0 lam2 (3 lam0 - 1), (9/2) lam0 lam2 (3 lam2 - 1),
    (9/2) lam0 lam3 (3 lam0 - 1), (9/2) lam0 lam3 (3 lam3 - 1),
    (9/2) lam1 lam2 (3 lam1 - 1), (9/2) lam1 lam2 (3 lam2 - 1),
    (9/2) lam1 lam3 (3 lam1 - 1), (9/2) lam1 lam3 (3 lam3 - 1),
    (9/2) lam2 lam3 (3 lam2 - 1), (9/2) lam2 lam3 (3 lam3 - 1),
    27 lam0 lam1 lam2,
    27 lam0 lam1 lam3,
    27 lam0 lam2 lam3,
    27 lam1 lam2 lam3
  }
];


(* ----- Tetrahedron HDiv RT0 ----- *)

TetRT0[xi_, eta_, zeta_, v_List, sigma_List] := Module[
  {r, vopp, vol},
  r = (1 - xi - eta - zeta) v[[1]] + xi v[[2]] + eta v[[3]] + zeta v[[4]];
  vopp = v;   (* face k is opposite vertex k *)
  vol = (1/6) Abs[Det[
    {v[[2]] - v[[1]], v[[3]] - v[[1]], v[[4]] - v[[1]]}]];
  Table[ sigma[[k]] (r - vopp[[k]]) / (3 vol), {k, 1, 4}]
];

TetRT0Div[v_List, sigma_List] := Module[{vol},
  vol = (1/6) Abs[Det[
    {v[[2]] - v[[1]], v[[3]] - v[[1]], v[[4]] - v[[1]]}]];
  sigma / vol
];


(* ----- Verification helpers ----- *)

TriPartitionOfUnity[basis_List] := Simplify[Total[basis] - 1] === 0;
TriDOFCount[basis_List] := Length[basis];


End[];   (* `Private` *)
EndPackage[];

(* ::Section:: *)
(*Self-test (run interactively after loading)*)


(* ::Input:: *)
(* Get["S:/Radia/01_GitHub/packages/radia-mcp/src/radia_mcp/mathematica/basis_functions/RadiaBasis.m"];                     *)
(*                                                                                     *)
(* (* Partition of unity for H1 Lagrange P_k *)                                       *)
(* Table[Simplify[Total[TriH1Pk[k][xi, eta]] - 1], {k, 1, 5}]                          *)
(*   -> {0, 0, 0, 0, 0}                                                                *)
(*                                                                                     *)
(* (* DOF counts *)                                                                    *)
(* Table[Length[TriH1Pk[k][xi, eta]], {k, 1, 5}]                                       *)
(*   -> {3, 6, 10, 15, 21}                                                             *)
(*                                                                                     *)
(* (* RT0 divergence -- should equal sigma/A_T per face *)                             *)
(* v = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}};                                              *)
(* sigma = {1, 1, 1};                                                                  *)
(* TriRT0Div[v, sigma]                                                                 *)
(*   -> {2, 2, 2}    (* A_T = 1/2, sigma_k / A_T = 2 each *)                           *)
