#!/usr/bin/env wolframscript
(* ::Package:: *)

(*
  surface_derham.wls -- HOIBC (Higher-Order Impedance Boundary Condition) read in the
  SAME framework: TOPOLOGY (surface de Rham / Betti) + the ANALYTIC DtN.  Companion to
  weakform_hodge.wls (the Hodge/material half), hodograph.wls, canonical.wls, and
  ../basis_functions/cohomology.wls (the global loops).

  HOIBC lives on the boundary SURFACE Gamma.  Two readings, both verified here:

  [A] TOPOLOGY.  Gamma carries its own de Rham complex  H1(G) -grad_G-> H(curl,G) -curl_G-> L2(G);
      the HOIBC trial fields ARE the surface FEEC basis (RWG / Nedelec-surface = surface
      Whitney forms).  The GLOBAL impedance modes = the surface HARMONIC 1-forms, whose count
      is the surface first Betti number b1(Gamma): a SPHERE boundary (simply connected) has
      b1 = 0 (NO global loop); a cylinder/torus boundary has b1 = 1 / 2 (global loops the
      impedance must carry).  b1 is TOPOLOGICAL (metric-independent, cohomology.wls) -- the
      impedance VALUES on the loops are geometric, the COUNT is not.

  [B] ANALYTIC DtN.  The exact exterior Dirichlet-to-Neumann (Steklov) operator on a sphere of
      radius R is DIAGONAL in spherical harmonics:  from the exterior-decaying harmonic
      u_n = r^{-(n+1)} Y_n,   Lambda_n = -u_n'/u_n |_{r=R} = (n+1)/R.   Kelvin / the infinite
      element reproduce this EXACTLY (they represent r^{-(n+1)} exactly).  HOIBC is a LOCAL
      (surface-differential-operator) approximation of this nonlocal DtN: a polynomial in the
      surface Laplace-Beltrami  Z = a + b Delta_G + c Delta_G^2 + ...  Since on the sphere
      Delta_G Y_n = -n(n+1)/R^2 Y_n, an order-k HOIBC matches Lambda_n at k+1 modes and then
      DEVIATES (the "knee") -- because Lambda(n) = (n+1)/R is NOT a finite polynomial in the
      Laplace-Beltrami eigenvalue n(n+1).  Higher ORDER = more matched modes, never exact.

  So: HOIBC basis = surface de Rham (topology); HOIBC operator = a LOCAL Pade of the analytic
  surface DtN (= the surface Hodge / Steklov).  Exact = Kelvin/IE (nonlocal); HOIBC = local approx.
*)

(* ---- minimal discrete 2-complex (for surface Betti; self-contained, cf. cohomology.wls) ---- *)
edgesOfTris[tris_] := Union[Flatten[(Sort /@ Subsets[#, {2}]) & /@ tris, 1]];
build2[tris_] := Module[{verts, edges, vi, ei, d0, d1},
  verts = Union[Flatten[tris]]; edges = edgesOfTris[tris];
  vi = AssociationThread[verts -> Range[Length[verts]]];
  ei = AssociationThread[edges -> Range[Length[edges]]];
  d0 = Normal@SparseArray[Join[Table[{ei[e], vi[e[[1]]]} -> -1, {e, edges}],
                               Table[{ei[e], vi[e[[2]]]} -> +1, {e, edges}]],
     {Length[edges], Length[verts]}];
  d1 = Normal@SparseArray[Flatten@Table[
     Module[{t = tris[[f]], cyc}, cyc = {{t[[1]], t[[2]]}, {t[[2]], t[[3]]}, {t[[3]], t[[1]]}};
       Table[With[{p = c[[1]], q = c[[2]]}, If[p < q, {f, ei[{p, q}]} -> +1, {f, ei[{q, p}]} -> -1]], {c, cyc}]],
     {f, Length[tris]}], {Length[tris], Length[edges]}];
  <|"nv" -> Length[verts], "ne" -> Length[edges], "nf" -> Length[tris], "d0" -> d0, "d1" -> d1|>];
betti2[K_] := With[{r0 = MatrixRank[K["d0"]], r1 = MatrixRank[K["d1"]]},
  {K["nv"] - r0, (K["ne"] - r1) - r0, K["nf"] - r1}];     (* {b0, b1, b2} *)

(* ---- sphere-surface operators (m=0 representative Y_n = P_n(cos th)) ---- *)
Yf[n_] := LegendreP[n, Cos[th]];
lapBeltrami[f_] := (1/(R^2 Sin[th])) D[Sin[th] D[f, th], th];   (* Delta_Gamma on the sphere *)
dtnEig[n_] := Module[{u = r^(-(n + 1)) Yf[n]}, Simplify[(-D[u, r]/u) /. r -> R]];   (* exact exterior DtN *)

(* ============================ self-test ============================ *)
Module[{ok = True, chk, Ksph, Kann, hoibc1, hoibc2, a1, b1c, a2, b2c, c2, sol1, sol2},
  chk[name_, cond_] := (Print["  ", If[cond, "ok  ", "FAIL"], name]; ok = ok && cond);
  Print["surface_derham.wls self-test (HOIBC via topology + the analytic DtN):"];

  (* ---- [A] TOPOLOGY: surface de Rham + Betti = global impedance loops ---- *)
  (* surface de Rham d.d = 0: for f on Gamma, d(df)=0 (curl_Gamma grad_Gamma f = 0) *)
  chk["[A0] surface de Rham d o d = 0 (curl_Gamma grad_Gamma f = 0) on the (th,ph) chart",
      Simplify[D[D[ff[th, ph], th], ph] - D[D[ff[th, ph], ph], th]] === 0];
  (* a SPHERE boundary is simply connected: b1 = 0 => HOIBC has NO global loop (tet boundary = S^2) *)
  Ksph = build2[{{1, 2, 3}, {1, 2, 4}, {1, 3, 4}, {2, 3, 4}}];
  chk["[A1] sphere boundary (tetrahedron = S^2) Betti = {1,0,1}: b1 = 0, NO global impedance loop",
      betti2[Ksph] == {1, 0, 1}];
  (* a CYLINDER/annulus boundary is multiply connected: b1 = 1 => ONE global loop the impedance carries *)
  Kann = build2[{{1, 2, 5}, {1, 5, 4}, {2, 3, 6}, {2, 6, 5}, {3, 1, 4}, {3, 4, 6}}];
  chk["[A2] cylinder/annulus boundary Betti = {1,1,0}: b1 = 1, ONE global impedance loop",
      betti2[Kann] == {1, 1, 0}];

  (* ---- [B] ANALYTIC DtN + HOIBC as its local Laplace-Beltrami-polynomial approximation ---- *)
  (* exact exterior DtN ladder: Lambda_n = (n+1)/R, diagonal in Y_n *)
  chk["[B0] analytic exterior DtN eigenvalue Lambda_n = (n+1)/R (from u_n = r^{-(n+1)} Y_n)",
      And @@ Table[Simplify[dtnEig[n] - (n + 1)/R] === 0, {n, 0, 4}]];
  (* Laplace-Beltrami eigenvalue on the sphere: Delta_Gamma Y_n = -n(n+1)/R^2 Y_n *)
  chk["[B1] surface Laplace-Beltrami eig: Delta_Gamma Y_n = -n(n+1)/R^2 Y_n",
      And @@ Table[Simplify[lapBeltrami[Yf[n]] + (n (n + 1)/R^2) Yf[n]] === 0, {n, 1, 4}]];

  (* order-1 HOIBC  Z = a + b Delta_Gamma : eigenvalue a - b n(n+1)/R^2.  Fit to Lambda_n at n=0,1. *)
  sol1 = First@Solve[{a1 - b1c (0)/R^2 == 1/R, a1 - b1c (2)/R^2 == 2/R}, {a1, b1c}];
  hoibc1[n_] := (a1 - b1c n (n + 1)/R^2) /. sol1;
  chk["[B2] order-1 HOIBC (a + b Delta_G) MATCHES the analytic DtN at n = 0,1 (fitted there)",
      Simplify[hoibc1[0] - 0 - 1/R] === 0 && Simplify[hoibc1[1] - 2/R] === 0];
  chk["[B3] order-1 HOIBC DEVIATES at n = 2 (gives 4/R vs exact 3/R): the knee -- local =/= exact DtN",
      Simplify[hoibc1[2] - 4/R] === 0 && Simplify[hoibc1[2] - 3/R] =!= 0];

  (* order-2 HOIBC  Z = a + b Delta + c Delta^2 : matches one MORE mode (n=0,1,2), still finite *)
  sol2 = First@Solve[{a2 == 1/R, a2 - b2c (2)/R^2 + c2 (2)^2/R^4 == 2/R,
                      a2 - b2c (6)/R^2 + c2 (6)^2/R^4 == 3/R}, {a2, b2c, c2}];
  hoibc2[n_] := (a2 - b2c n (n + 1)/R^2 + c2 (n (n + 1))^2/R^4) /. sol2;
  chk["[B4] order-2 HOIBC matches n = 0,1,2 (one more mode) -- higher ORDER = more matched DtN modes",
      Simplify[hoibc2[0] - 1/R] === 0 && Simplify[hoibc2[1] - 2/R] === 0 &&
      Simplify[hoibc2[2] - 3/R] === 0];
  chk["[B5] order-2 HOIBC STILL deviates at n = 3 (finite local order can never be exact; Kelvin/IE are)",
      Simplify[hoibc2[3] - 4/R] =!= 0];

  Print[If[ok, "ALL PASS", "SOME FAILED"]];
];
