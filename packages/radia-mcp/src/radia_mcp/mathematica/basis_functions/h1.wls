#!/usr/bin/env wolframscript
(* ::Package:: *)

(*
  h1.wls -- NGSolve-style hierarchical H1 (scalar) shape functions, clean wls.
  Zaglmayr construction (ref: README.md, ../notes_fem_hcurl.md).  Foundation:
  recursive_pol.wls (inlined here so this file runs standalone).

  Simplex elements (this file): Trig (2D), Tet (3D).
  Quad/Hex (tensor), Prism, Pyramid: separate files / sections (see plan in README).

  Per element, order p, the basis splits into:
    vertex   : the barycentric coordinate  lam_v
    edge     : scaled integrated Legendre LS_i along each edge (i = 2..p),
               vanishes on all other edges and at every vertex
    face/cell: interior bubbles  (product of the entity's lam) * complete poly
  The vertex block is a partition of unity; the higher blocks are the hierarchy.
*)

(* ---- polynomial foundation (= recursive_pol.wls) ---- *)
l[i_, x_] := LegendreP[i, x];
L[i_, x_] := L[i, x] = Simplify[Integrate[l[i - 1, z], {z, -1, x}]];
(* scaled integrated Legendre as a CLEAN homogeneous degree-n polynomial in (s,t):
   t^n L[n, s/t] = Sum_k coeff_k(L_n) s^k t^(n-k)  -- no s/t division, no 0/0 at t=0. *)
LS[n_, s_, t_] := LS[n, s, t] = Module[{u},
   Expand@Sum[Coefficient[L[n, u], u, k] s^k t^(n - k), {k, 0, n}]];

(* edge bubble on barycentric pair (la, lb): LS_i((lb-la),(la+lb)) -- on the edge
   it is the 1D integrated-Legendre bubble; zero on the other edges + all vertices. *)
edgeBubble[la_, lb_, i_] := LS[i, lb - la, la + lb];

(* ============================ TRIANGLE ============================ *)
(* barycentric: lam0 = 1-x-y, lam1 = x, lam2 = y ; vertices 0,1,2 ;
   edges {0,1},{0,2},{1,2} ; one interior face. *)
H1Trig[p_] := Module[{lam, V, E, edges, vtx, edg, intr, a, b},
  lam = {1 - x - y, x, y};
  vtx = lam;                                                   (* 3 vertex funcs *)
  edges = {{1, 2}, {1, 3}, {2, 3}};                            (* 1-indexed *)
  edg = Flatten@Table[edgeBubble[lam[[e[[1]]]], lam[[e[[2]]]], i],
                      {e, edges}, {i, 2, p}];                  (* 3(p-1) *)
  intr = Flatten@Table[lam[[1]] lam[[2]] lam[[3]] lam[[1]]^a lam[[2]]^b,
                       {a, 0, p - 3}, {b, 0, p - 3 - a}];      (* (p-1)(p-2)/2 *)
  Join[vtx, edg, intr]];

(* ============================== TET ============================== *)
(* lam0=1-x-y-z, lam1=x, lam2=y, lam3=z ; 4 vtx, 6 edges, 4 faces, 1 cell. *)
H1Tet[p_] := Module[{lam, vtx, edges, edg, faces, fac, cell, a, b, c},
  lam = {1 - x - y - z, x, y, z};
  vtx = lam;                                                   (* 4 *)
  edges = {{1,2},{1,3},{1,4},{2,3},{2,4},{3,4}};
  edg = Flatten@Table[edgeBubble[lam[[e[[1]]]], lam[[e[[2]]]], i],
                      {e, edges}, {i, 2, p}];                  (* 6(p-1) *)
  faces = {{1,2,3},{1,2,4},{1,3,4},{2,3,4}};
  fac = Flatten@Table[
     With[{l1 = lam[[f[[1]]]], l2 = lam[[f[[2]]]], l3 = lam[[f[[3]]]]},
       Table[l1 l2 l3 l1^a l2^b, {a, 0, p - 3}, {b, 0, p - 3 - a}]],
     {f, faces}];                                              (* 4*(p-1)(p-2)/2 *)
  cell = Flatten@Table[lam[[1]] lam[[2]] lam[[3]] lam[[4]] lam[[1]]^a lam[[2]]^b lam[[3]]^c,
                       {a, 0, p - 4}, {b, 0, p - 4 - a}, {c, 0, p - 4 - a - b}]; (* (p-1)(p-2)(p-3)/6 *)
  Join[vtx, edg, fac, cell]];

(* ===================== TENSOR ELEMENTS (Quad, Hex) ===================== *)
(* 1D hierarchical shapes on [0,1]: vertices {1-t, t} + bubbles L_i(2t-1)
   (integrated Legendre, vanishes at t=0,1), i=2..p.  Length p+1. *)
oneD[p_, t_] := Join[{1 - t, t}, Table[L[i, 2 t - 1], {i, 2, p}]];

(* Quad [0,1]^2: full tensor Q_p, dim (p+1)^2.  Entity split: index in {1,2} = a
   vertex factor, index >= 3 = a bubble factor (vertex/edge/face by #bubbles). *)
H1Quad[p_] := Module[{sx = oneD[p, x], sy = oneD[p, y], n = p + 1},
  Flatten@Table[sx[[i]] sy[[j]], {i, n}, {j, n}]];

(* Hex [0,1]^3: full tensor Q_p, dim (p+1)^3. *)
H1Hex[p_] := Module[{sx = oneD[p, x], sy = oneD[p, y], sz = oneD[p, z], n = p + 1},
  Flatten@Table[sx[[i]] sy[[j]] sz[[k]], {i, n}, {j, n}, {k, n}]];

(* ===== PRISM (anisotropic tensor: Trig in (x,y)  (x)  Segm in z) ===== *)
(* dim = dim(Trig_p) * (p+1) ; 6 vertices = 3 trig-vtx * 2 segm-vtx. *)
H1Prism[p_] := Module[{tr = H1Trig[p], sz = oneD[p, z]},
  Flatten@Table[tr[[i]] sz[[j]], {i, Length[tr]}, {j, Length[sz]}]];

(* ===== PYRAMID (RATIONAL; ported from NGSolve h1hofe_impl ET_PYRAMID) ===== *)
(* reference: base [0,1]^2 at z=0, apex (0,0,1).  Collapsed coords xt=x/(1-z),
   yt=y/(1-z) map the shrinking square at height z to the unit square, so the H1
   shapes are RATIONAL in (1-z) -- the hallmark of the pyramid (the FEEC-hard
   element NGSolve implements, NOT a polynomial).  Lowest order = the 5 vertex
   functions lambda3d (this block: 4 base*(1-z) collapsed-bilinear + apex z).
   Higher order (p>=2: horizontal/vertical edge + trig/quad face bubbles via the
   scaled integrated-Legendre / trig orthopols) is the next increment. *)
H1PyramidVertices[] := Module[{xt = x/(1 - z), yt = y/(1 - z), lam},
  lam = {(1 - xt) (1 - yt), xt (1 - yt), xt yt, (1 - xt) yt};
  Append[Table[(1 - z) lam[[i]], {i, 4}], z]];

(* ============================ self-test ============================ *)
Module[{ok = True, chk, rk, pts2, pts3},
  chk[name_, cond_] := (Print["  ", If[cond, "ok  ", "FAIL"], name]; ok = ok && cond);
  Print["h1.wls self-test:"];
  (* dimension formulas *)
  chk["dim Trig(p) = (p+1)(p+2)/2, p=1..5",
      And @@ Table[Length[H1Trig[p]] == (p+1)(p+2)/2, {p, 1, 5}]];
  chk["dim Tet(p)  = (p+1)(p+2)(p+3)/6, p=1..5",
      And @@ Table[Length[H1Tet[p]] == (p+1)(p+2)(p+3)/6, {p, 1, 5}]];
  (* vertex partition of unity *)
  chk["Trig vertices sum to 1", Simplify[Total[Take[H1Trig[3], 3]]] === 1];
  chk["Tet vertices sum to 1",  Simplify[Total[Take[H1Tet[3], 4]]] === 1];
  (* edge bubble vanishes at all vertices (trig, i=2..4 on edge {0,1}) *)
  chk["Trig edge bubble = 0 at all 3 vertices",
      And @@ Flatten@Table[(edgeBubble[1 - x - y, x, i] /.
            {{x -> 0, y -> 0}, {x -> 1, y -> 0}, {x -> 0, y -> 1}}) == {0, 0, 0},
          {i, 2, 4}]];
  (* linear independence: rank of values at random points == dim *)
  SeedRandom[1];
  pts2 = RandomReal[{0.05, 0.3}, {40, 2}];
  rk = MatrixRank[Table[H1Trig[4] /. {x -> pp[[1]], y -> pp[[2]]}, {pp, pts2}], Tolerance -> 1*^-9];
  chk["Trig(4) shape funcs independent (rank=15)", rk == 15];
  pts3 = RandomReal[{0.05, 0.25}, {60, 3}];
  rk = MatrixRank[Table[H1Tet[3] /. {x -> pp[[1]], y -> pp[[2]], z -> pp[[3]]}, {pp, pts3}], Tolerance -> 1*^-9];
  chk["Tet(3) shape funcs independent (rank=20)", rk == 20];
  (* tensor elements *)
  chk["dim Quad(p) = (p+1)^2, p=1..5", And @@ Table[Length[H1Quad[p]] == (p+1)^2, {p, 1, 5}]];
  chk["dim Hex(p)  = (p+1)^3, p=1..5", And @@ Table[Length[H1Hex[p]] == (p+1)^3, {p, 1, 5}]];
  chk["Quad bilinear vertices sum to 1",
      Simplify[Sum[oneD[3, x][[i]] oneD[3, y][[j]], {i, 2}, {j, 2}]] === 1];
  chk["Hex trilinear vertices sum to 1",
      Simplify[Sum[oneD[3, x][[i]] oneD[3, y][[j]] oneD[3, z][[k]], {i, 2}, {j, 2}, {k, 2}]] === 1];
  chk["1D bubbles vanish at t=0,1 (i=2..5)",
      And @@ Flatten@Table[(L[i, 2 t - 1] /. {{t -> 0}, {t -> 1}}) == {0, 0}, {i, 2, 5}]];
  SeedRandom[2];
  rk = MatrixRank[Table[H1Quad[3] /. {x -> pp[[1]], y -> pp[[2]]}, {pp, RandomReal[{0, 1}, {40, 2}]}], Tolerance -> 1*^-9];
  chk["Quad(3) independent (rank=16)", rk == 16];
  rk = MatrixRank[Table[H1Hex[2] /. {x -> pp[[1]], y -> pp[[2]], z -> pp[[3]]}, {pp, RandomReal[{0, 1}, {50, 3}]}], Tolerance -> 1*^-9];
  chk["Hex(2) independent (rank=27)", rk == 27];
  (* prism *)
  chk["dim Prism(p) = dim(Trig_p)*(p+1), p=1..4",
      And @@ Table[Length[H1Prism[p]] == Length[H1Trig[p]] (p + 1), {p, 1, 4}]];
  chk["Prism 6 vertices sum to 1",
      Simplify[Sum[H1Trig[3][[i]] oneD[3, z][[j]], {i, 3}, {j, 2}]] === 1];
  (* pyramid (rational) *)
  chk["Pyramid 5 vertices sum to 1", Simplify[Total[H1PyramidVertices[]]] === 1];
  chk["Pyramid vertices are RATIONAL (not polynomial)",
      ! And @@ (PolynomialQ[#, {x, y, z}] & /@ H1PyramidVertices[])];
  chk["Pyramid base-vertex denominator has (1-z)",
      Simplify[Denominator[Together[H1PyramidVertices[][[1]]]] /. z -> 0] =!= 0 &&
      ! FreeQ[Together[H1PyramidVertices[][[1]]], 1 - z | -1 + z]];
  rk = MatrixRank[Table[H1PyramidVertices[] /. {x -> pp[[1]] (1 - pp[[3]]), y -> pp[[2]] (1 - pp[[3]]), z -> pp[[3]]},
                       {pp, RandomReal[{0.1, 0.7}, {12, 3}]}], Tolerance -> 1*^-9];
  chk["Pyramid 5 vertices independent (rank=5)", rk == 5];
  Print[If[ok, "ALL PASS", "SOME FAILED"]];
];
