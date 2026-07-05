#!/usr/bin/env wolframscript
(* ::Package:: *)

(*
  derham.wls -- the lowest-order Whitney de Rham complex on the reference TET,
  the cohomology / loop foundation for the HDiv-VIM VIM solver.

      W0 (H1, nodal)  --grad-->  W1 (H(curl), edge)  --curl-->  W2 (H(div), face)  --div-->  W3 (L2, cell)

  The Whitney forms make this complex EXACT by construction:
    curl o grad = 0,  div o curl = 0  (d o d = 0)
  so the div-free fields ("loops", ker div) are exactly  curl(H(curl)) = curl(W1)
  -- the LOCAL (contractible) loop space, field-null by construction.  The GLOBAL
  loops = the de Rham cohomology H^1 (harmonic forms, one per topological hole;
  e.g. a C-type yoke ring has H^1 = 1) -- added on the assembled mesh as the
  tree-cotree / harmonic generators, NOT visible on a single element.

  Reference tet: lam0=1-x-y-z, lam1=x, lam2=y, lam3=z (constant gradients).
*)

lam = {1 - x - y - z, x, y, z};
grad[f_] := D[f, {{x, y, z}}];
gl = grad /@ lam;                                  (* constant gradients of barycentrics *)
crl[F_] := Curl[F, {x, y, z}];
dvg[F_] := Div[F, {x, y, z}];

(* edges (1-indexed verts) and faces of the tet *)
edges = {{1,2},{1,3},{1,4},{2,3},{2,4},{3,4}};
faces = {{2,3,4},{1,3,4},{1,2,4},{1,2,3}};

(* Whitney forms *)
W0 = lam;                                          (* H1 vertex *)
W1[{a_, b_}] := lam[[a]] gl[[b]] - lam[[b]] gl[[a]];            (* H(curl) edge *)
W2[{a_, b_, c_}] := 2 (lam[[a]] Cross[gl[[b]], gl[[c]]] +
                       lam[[b]] Cross[gl[[c]], gl[[a]]] +
                       lam[[c]] Cross[gl[[a]], gl[[b]]]);       (* H(div) face *)

(* ============================ self-test ============================ *)
Module[{ok = True, chk, e, f},
  chk[name_, cond_] := (Print["  ", If[cond, "ok  ", "FAIL"], name]; ok = ok && cond);
  Print["derham.wls self-test (Whitney de Rham complex on tet):"];

  (* d o d = 0 : curl(grad W0) = 0 *)
  chk["curl(grad W0) = 0 (4 verts)",
      And @@ Table[Simplify[crl[gl[[i]]]] === {0, 0, 0}, {i, 4}]];

  (* d o d = 0 : div(curl W1) = 0 (6 edges) -- the LOOP is div-free by construction *)
  chk["div(curl W1) = 0 (6 edges) [loops are div-free]",
      And @@ Table[Simplify[dvg[crl[W1[e]]]] === 0, {e, edges}]];

  (* curl(W1_edge) is CONSTANT = 2 grad_a x grad_b  (the edge->face incidence) *)
  chk["curl(W1[{a,b}]) = 2 grad_a x grad_b (constant)",
      And @@ Table[Simplify[crl[W1[e]] - 2 Cross[gl[[e[[1]]]], gl[[e[[2]]]]]] === {0, 0, 0}, {e, edges}]];

  (* div(W2_face) is CONSTANT (face->cell incidence; nonzero const) *)
  chk["div(W2[face]) = nonzero constant (4 faces)",
      And @@ Table[With[{d = Simplify[dvg[W2[f]]]},
                    FreeQ[d, x | y | z] && d =!= 0], {f, faces}]];

  (* W1 reproduces the de Rham pairing: W1 tangential dof = Kronecker on edges.
     The line integral of W1[e] along edge e' equals delta_{e,e'} (Whitney property). *)
  chk["W1 edge line-integral = Kronecker delta (Whitney dof)",
      With[{vc = {{0,0,0},{1,0,0},{0,1,0},{0,0,1}}},          (* tet vertex coords *)
       And @@ Flatten@Table[
         With[{a = ei[[1]], b = ei[[2]], c = ej[[1]], d = ej[[2]]},
          (* integral over edge (c,d) of W1[ei] . tangent *)
          Module[{P = vc[[c]] + t (vc[[d]] - vc[[c]]), tang = vc[[d]] - vc[[c]], val},
           val = Integrate[(W1[ei] /. {x -> P[[1]], y -> P[[2]], z -> P[[3]]}) . tang, {t, 0, 1}];
           Simplify[val] == If[ei === ej, 1, 0]]],
         {ei, edges}, {ej, edges}]]];

  Print[If[ok, "ALL PASS", "SOME FAILED"]];
];
