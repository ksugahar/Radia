#!/usr/bin/env wolframscript
(* ::Package:: *)

(*
  cohomology.wls -- the DISCRETE de Rham complex and its cohomology H^1: the GLOBAL
  "loops" the HDiv-VIM VIM solver must carry on a multiply-connected body.

  derham.wls / h1 / hcurl / hdiv give the LOCAL (per-element) FEEC basis; the div-free
  circulations there are curl(H(curl)) -- the LOCAL, contractible loops, field-null by
  construction.  But a body with a HOLE (the C-type yoke RING, a torus, a coil window)
  carries extra div-free fields that are NOT curl of anything global: the harmonic /
  cohomology forms.  Their count is the first Betti number b1 = dim H^1 (one per hole).
  These are EXACTLY the global loops  rad_hacapk.cpp::BuildLoopBasis  enumerates by
  tree-cotree (b1 = E - V + C for the 1-skeleton); each cotree edge closes one cycle.

  This file builds the signed incidence (coboundary) operators of a small simplicial
  complex,

        C^0(verts) --d0=grad--> C^1(edges) --d1=curl--> C^2(faces)

  verifies it IS a complex (d1.d0 = 0, the discrete  curl grad = 0), computes the Betti
  numbers by rank-nullity, and EXTRACTS a harmonic H^1 generator (a cycle that is not a
  gradient = the loop around the hole).  Integer incidence => EXACT rank, no tolerance.

  Betti (rank-nullity on the cochain complex):
     b0 = #V - rank(d0)                      (connected components)
     b1 = (#E - rank(d1)) - rank(d0)         (independent loops / holes)
     b2 = #F - rank(d1)                      (voids)
*)

(* ---- signed coboundary operators of a 2-complex given by a triangle list ---- *)
edgesOfTris[tris_] := Union[Flatten[(Sort /@ Subsets[#, {2}]) & /@ tris, 1]];

build2Complex[tris_] := Module[{verts, edges, vi, ei, d0, d1, d1rules},
  verts = Union[Flatten[tris]];
  edges = edgesOfTris[tris];
  vi = AssociationThread[verts -> Range[Length[verts]]];
  ei = AssociationThread[edges -> Range[Length[edges]]];
  (* d0 : grad, edge (a<b) -> -1 at a, +1 at b *)
  d0 = Normal@SparseArray[
     Join[Table[{ei[e], vi[e[[1]]]} -> -1, {e, edges}],
          Table[{ei[e], vi[e[[2]]]} -> +1, {e, edges}]],
     {Length[edges], Length[verts]}];
  (* d1 : curl, oriented triangle boundary (a->b->c->a), signed to canonical edges *)
  d1rules = Flatten@Table[
     Module[{t = tris[[fi]], cyc},
      cyc = {{t[[1]], t[[2]]}, {t[[2]], t[[3]]}, {t[[3]], t[[1]]}};
      Table[With[{p = c[[1]], q = c[[2]]},
        If[p < q, {fi, ei[{p, q}]} -> +1, {fi, ei[{q, p}]} -> -1]], {c, cyc}]],
     {fi, Length[tris]}];
  d1 = Normal@SparseArray[d1rules, {Length[tris], Length[edges]}];
  <|"verts" -> verts, "edges" -> edges, "tris" -> tris, "d0" -> d0, "d1" -> d1|>];

(* ---- 1-complex (graph): edges given directly, no faces (d1 = 0) ---- *)
build1Complex[edgeList_] := Module[{verts, edges, vi, d0},
  edges = Sort /@ edgeList;
  verts = Union[Flatten[edges]];
  vi = AssociationThread[verts -> Range[Length[verts]]];
  d0 = Normal@SparseArray[
     Join[Table[{i, vi[edges[[i, 1]]]} -> -1, {i, Length[edges]}],
          Table[{i, vi[edges[[i, 2]]]} -> +1, {i, Length[edges]}]],
     {Length[edges], Length[verts]}];
  <|"verts" -> verts, "edges" -> edges, "tris" -> {},
    "d0" -> d0, "d1" -> ConstantArray[0, {0, Length[edges]}]|>];

betti[K_] := Module[{nv = Length[K["verts"]], ne = Length[K["edges"]],
                     nf = Length[K["tris"]], r0 = MatrixRank[K["d0"]], r1},
  r1 = If[nf == 0, 0, MatrixRank[K["d1"]]];
  {nv - r0, (ne - r1) - r0, nf - r1}];     (* {b0, b1, b2} *)

(* harmonic 1-form count under an ARBITRARY edge inner product M1 (SPD): a harmonic h is
   CLOSED (d1 h = 0) and CO-CLOSED (d0^T M1 h = 0, the M1-adjoint of d0).  Hodge theorem:
   dim(harmonic) = b1 for ANY metric M1 -- the Hodge star / metric / material CANNOT change
   the topological Betti number.  This is "topology is NOT in the Hodge" (transformation
   optics moves the geometry, never the holes). *)
harmonicDim[K_, M1_] := Length[K["edges"]] -
  MatrixRank[Join[K["d1"], Transpose[K["d0"]].M1], Tolerance -> 1.*^-9];

(* harmonic H^1 generator: a cycle (d1 g = 0) that is NOT a gradient (g not in im d0) *)
isComplex[K_] := With[{nf = Length[K["tris"]]},
  nf == 0 || K["d1"].K["d0"] == ConstantArray[0, {nf, Length[K["verts"]]}]];

h1Generators[K_] := Module[{ker, cob, gens, augR, baseR},
  ker = If[Length[K["tris"]] == 0, IdentityMatrix[Length[K["edges"]]], NullSpace[K["d1"]]];
  cob = Transpose[K["d0"]];                 (* rows = im(d0) basis (coboundaries) *)
  baseR = MatrixRank[cob];
  (* greedily keep cycle vectors that raise the rank over the coboundary space *)
  gens = {};
  Do[If[MatrixRank[Join[cob, gens, {v}]] > baseR + Length[gens], AppendTo[gens, v]], {v, ker}];
  gens];

(* ============================ self-test ============================ *)
Module[{ok = True, chk, Kfill, Kcyc, Kann, K2, b, gens},
  chk[name_, cond_] := (Print["  ", If[cond, "ok  ", "FAIL"], name]; ok = ok && cond);
  Print["cohomology.wls self-test (discrete de Rham complex + Betti + H^1 generator):"];

  (* filled triangle = a disk: contractible, b = {1,0,0} *)
  Kfill = build2Complex[{{1, 2, 3}}];
  chk["d1.d0 = 0 (discrete curl grad = 0) on filled triangle", isComplex[Kfill]];
  chk["filled triangle Betti = {1,0,0} (contractible disk)", betti[Kfill] == {1, 0, 0}];

  (* cycle graph C5 (a RING, no fill) = the yoke loop: b1 = 1 *)
  Kcyc = build1Complex[{{1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 1}}];
  b = betti[Kcyc];
  chk["cycle C5 Betti b0=1, b1=1 (the global ring loop)", b[[1]] == 1 && b[[2]] == 1];
  chk["tree-cotree: b1 = E - V + C  (BuildLoopBasis formula)",
      b[[2]] == Length[Kcyc["edges"]] - Length[Kcyc["verts"]] + b[[1]]];
  gens = h1Generators[Kcyc];
  chk["C5 has exactly 1 cohomology generator (one cotree edge)", Length[gens] == 1];

  (* triangulated ANNULUS (outer 123, inner hole 456): b = {1,1,0} *)
  Kann = build2Complex[{{1, 2, 5}, {1, 5, 4}, {2, 3, 6}, {2, 6, 5}, {3, 1, 4}, {3, 4, 6}}];
  chk["annulus is a complex (d1.d0 = 0)", isComplex[Kann]];
  b = betti[Kann];
  chk["annulus Betti = {1,1,0} (one hole, faces present but hole unfilled)", b == {1, 1, 0}];
  gens = h1Generators[Kann];
  chk["annulus has exactly 1 harmonic H^1 generator (loop around hole)", Length[gens] == 1];
  chk["the generator IS a cycle:  d1 . g = 0",
      And @@ (# == 0 & /@ (Kann["d1"].gens[[1]]))];
  chk["the generator is NOT a gradient: g not in im(d0)",
      MatrixRank[Append[Transpose[Kann["d0"]], gens[[1]]]] >
        MatrixRank[Transpose[Kann["d0"]]]];

  (* two disjoint triangles: b0 = 2 (components) *)
  K2 = build2Complex[{{1, 2, 3}, {4, 5, 6}}];
  chk["two disjoint filled triangles: b0 = 2 components", betti[K2][[1]] == 2];

  (* genus-2-ish: two holes share the b1 count -- annulus + extra hole *)
  Module[{Ktwo, bb},
   (* outer 1-2-3, two inner holes 4-5-6 and 7-8-9 is fiddly to mesh by hand;
      instead: figure-eight graph (two cycles sharing a vertex) -> b1 = 2 *)
   Ktwo = build1Complex[{{1, 2}, {2, 3}, {3, 1}, {1, 4}, {4, 5}, {5, 1}}];
   bb = betti[Ktwo];
   chk["figure-eight graph: b1 = 2 (two independent loops)", bb[[2]] == 2];
   chk["figure-eight: 2 cohomology generators (2 cotree edges)",
       Length[h1Generators[Ktwo]] == 2]];

  (* ---- metric-INDEPENDENCE of Betti: the Hodge / metric CANNOT change the topology ----
     The harmonic 1-form count = b1 for ANY SPD edge inner product M1 (Hodge theorem), so a
     coordinate map / material / Hodge star moves the GEOMETRY but never the holes.  This is
     "topology is NOT in the Hodge" -- the half of transformation optics that does NOT move. *)
  Module[{Kann2, Kfig, spd},
   Kann2 = build2Complex[{{1, 2, 5}, {1, 5, 4}, {2, 3, 6}, {2, 6, 5}, {3, 1, 4}, {3, 4, 6}}]; (* annulus, b1=1 *)
   Kfig = build1Complex[{{1, 2}, {2, 3}, {3, 1}, {1, 4}, {4, 5}, {5, 1}}];                     (* figure-eight, b1=2 *)
   spd[n_, seed_] := (SeedRandom[seed];
     With[{A = RandomReal[{-1, 1}, {n, n}]}, A.Transpose[A] + n IdentityMatrix[n]]); (* random SPD *)
   chk["annulus: harmonic-1-form dim = b1 = 1 under the trivial metric M1 = I",
       harmonicDim[Kann2, IdentityMatrix[Length[Kann2["edges"]]]] == 1];
   chk["annulus: harmonic-1-form dim STILL = 1 under TWO random SPD metrics (Hodge can't change Betti)",
       harmonicDim[Kann2, spd[Length[Kann2["edges"]], 1]] == 1 &&
       harmonicDim[Kann2, spd[Length[Kann2["edges"]], 2]] == 1];
   chk["figure-eight: harmonic-1-form dim = b1 = 2 under I AND a random SPD metric (topology not in the Hodge)",
       harmonicDim[Kfig, IdentityMatrix[Length[Kfig["edges"]]]] == 2 &&
       harmonicDim[Kfig, spd[Length[Kfig["edges"]], 3]] == 2]];

  Print[If[ok, "ALL PASS", "SOME FAILED"]];
];
