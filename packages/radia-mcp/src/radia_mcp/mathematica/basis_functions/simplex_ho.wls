#!/usr/bin/env wolframscript
(* ::Package:: *)

(*
  simplex_ho.wls -- high-order H(curl)/H(div) on the SIMPLEX (tet, triangle) as
  explicit BDM / Raviart-Thomas / Nedelec-first-kind polynomial SPACES.

  NGSolve's hdivhofe/hcurlhofe build these with a Zaglmayr integrated-Jacobi (Dubiner)
  HIERARCHICAL basis -- intricate to transcribe and basis-dependent.  For a symbolic
  VIM reference what matters is the SPACE (the 1/r field operator integrates whatever
  basis spans it), so we build the spaces in their textbook complete-polynomial form,
  which is clean and fully verifiable:

     BDM_p (H(div))  = [P_p]^d, p >= 1
     RT_k  (H(div))  = [P_k]^d  (+)  x~ * Ptilde_k        (x~ = (x,y,z), radial)
     N_k   (H(curl)) = [P_k]^d  (+)  S_k,   S_k = { p in [Ptilde_{k+1}]^d : p . x~ = 0 }

  The key study distinction is div(BDM_p) = P_(p-1), while div(RT_p) = P_p.
  BDM_p is a strict subspace of RT_p, but both have the same divergence-free
  kernel dimension.  Thus RT adds charge-carrying modes, not loop modes.  The
  de Rham maps grad(P_(k+1)) c N_k and curl(N_k) c BDM_k c RT_k are exact.
  Lowest order RT (k=0) reproduces Whitney: RT_0 = 4 tet faces / 3 trig edges
  (cf. derham.wls W2), N_0 = 6 tet / 3 trig edges (W1).  Standard BDM starts
  at p=1.

  Reference simplex (x,y,z>=0, x+y+z<=1).  These are the tetrahedral moment-path high-order
  companions to the tensor Hex/Quad spaces in hdiv.wls / hcurl.wls.
*)

gradd[f_] := D[f, {{x, y, z}}];
divv[V_]  := D[V[[1]], x] + D[V[[2]], y] + D[V[[3]], z];
curlv[V_] := Curl[V, {x, y, z}];

(* ----- monomial sets ----- *)
monos3[k_] := Flatten@Table[x^a y^b z^c, {a, 0, k}, {b, 0, k - a}, {c, 0, k - a - b}]; (* deg<=k *)
homo3[k_]  := Flatten@Table[x^a y^b z^(k - a - b), {a, 0, k}, {b, 0, k - a}];          (* deg==k *)
exps3[k_]  := Flatten[Table[{a, b, k - a - b}, {a, 0, k}, {b, 0, k - a}], 1];
polyToVec[poly_, exps_] := Module[{cr},
  cr = Association[(#[[1]] -> #[[2]]) & /@ CoefficientRules[Expand[poly], {x, y, z}]];
  Lookup[cr, exps, 0]];

(* ----- TET RT_k (Raviart-Thomas, H(div)) ----- *)
RTtet[k_] := Join[
  {#, 0, 0} & /@ monos3[k], {0, #, 0} & /@ monos3[k], {0, 0, #} & /@ monos3[k],
  ({x, y, z} #) & /@ homo3[k]];

(* ----- TET BDM_p = [P_p]^3, p>=1 (NGSolve bare HDiv on simplices) ----- *)
BDMtet[p_Integer /; p >= 1] := Join[
  {#, 0, 0} & /@ monos3[p], {0, #, 0} & /@ monos3[p], {0, 0, #} & /@ monos3[p]];

(* Uniform public names across simplex and tensor scripts. *)
HDivTetRT[k_] := RTtet[k];
HDivTetBDM[p_Integer /; p >= 1] := BDMtet[p];

(* ----- TET S_k = {p in [homo_{k+1}]^3 : p . (x,y,z) = 0} ----- *)
Stet[k_] := Module[{B, M, ns},
  B = Join[{#, 0, 0} & /@ homo3[k + 1], {0, #, 0} & /@ homo3[k + 1], {0, 0, #} & /@ homo3[k + 1]];
  M = Transpose[polyToVec[#.{x, y, z}, exps3[k + 2]] & /@ B];     (* rows = homo_{k+2}, cols = B *)
  ns = NullSpace[M];
  (#.B) & /@ ns];                                                 (* each coeff-row -> a 3-vector *)

(* ----- TET N_k (Nedelec first kind, H(curl)) ----- *)
Ntet[k_] := Join[
  {#, 0, 0} & /@ monos3[k], {0, #, 0} & /@ monos3[k], {0, 0, #} & /@ monos3[k],
  Stet[k]];

(* ----- TRIANGLE (2D, z=0 plane; vectors carry a 0 z-component) ----- *)
monos2[k_] := Flatten@Table[x^a y^(d - a), {d, 0, k}, {a, 0, d}];                      (* deg<=k *)
homo2[k_]  := Table[x^a y^(k - a), {a, 0, k}];
RTtrig[k_] := Join[
  {#, 0, 0} & /@ monos2[k], {0, #, 0} & /@ monos2[k],
  ({x, y, 0} #) & /@ homo2[k]];
BDMtrig[p_Integer /; p >= 1] := Join[
  {#, 0, 0} & /@ monos2[p], {0, #, 0} & /@ monos2[p]];
HDivTrigRT[k_] := RTtrig[k];
HDivTrigBDM[p_Integer /; p >= 1] := BDMtrig[p];
div2[V_] := D[V[[1]], x] + D[V[[2]], y];

(* dimension formulas *)
dimPk3[k_] := Binomial[k + 3, 3];
dimPk2[k_] := Binomial[k + 2, 2];
dimRTtet[k_] := 3 dimPk3[k] + Binomial[k + 2, 2];
dimBDMtet[p_] := 3 dimPk3[p];
dimRTtrig[k_] := 2 dimPk2[k] + (k + 1);
dimBDMtrig[p_] := 2 dimPk2[p];
(* dim S_k = 3 dim(homo_{k+1}) - dim(homo_{k+2}) = 3 C(k+3,2) - C(k+4,2)  (p.x~ surjective) *)
dimNtet[k_]  := 3 dimPk3[k] + (3 Binomial[k + 3, 2] - Binomial[k + 4, 2]);

(* Compact symbolic ledgers for equal-order RT/BDM studies. *)
HDivTetFamilyLedger[p_Integer /; p >= 1] := <|
  "Order" -> p,
  "BDMDoF" -> dimBDMtet[p], "RTDoF" -> dimRTtet[p],
  "BDMDivDoF" -> dimPk3[p - 1], "RTDivDoF" -> dimPk3[p],
  "SharedKerDivDoF" -> dimBDMtet[p] - dimPk3[p - 1],
  "ExtraRTChargeDoF" -> dimPk3[p] - dimPk3[p - 1]|>;

HDivTrigFamilyLedger[p_Integer /; p >= 1] := <|
  "Order" -> p,
  "BDMDoF" -> dimBDMtrig[p], "RTDoF" -> dimRTtrig[p],
  "BDMDivDoF" -> dimPk2[p - 1], "RTDivDoF" -> dimPk2[p],
  "SharedKerDivDoF" -> dimBDMtrig[p] - dimPk2[p - 1],
  "ExtraRTChargeDoF" -> dimPk2[p] - dimPk2[p - 1]|>;

(* ============================ self-test ============================ *)
Module[{ok = True, chk, pts, pts2, p, rB, rA},
  chk[name_, cond_] := (Print["  ", If[cond, "ok  ", "FAIL"], name]; ok = ok && cond);
  Print["simplex_ho.wls self-test (BDM_p / RT_p / N_p de Rham on the simplex):"];
  SeedRandom[31];
  pts = Table[{x -> RandomReal[{.05, .25}], y -> RandomReal[{.05, .25}], z -> RandomReal[{.05, .25}]}, {60}];
  SeedRandom[32];
  pts2 = Table[{x -> RandomReal[{.05, .3}], y -> RandomReal[{.05, .3}]}, {30}];

  (* ---- dimensions ---- *)
  chk["dim BDM_tet(p) = 3 C(p+3,3), p=1..3",
      And @@ Table[Length[BDMtet[p]] == dimBDMtet[p], {p, 1, 3}]];
  chk["dim RT_tet(k) = 3 C(k+3,3)+C(k+2,2), k=0..2",
      And @@ Table[Length[RTtet[k]] == dimRTtet[k], {k, 0, 2}]];
  chk["dim N_tet(k)  = 3 C(k+3,3)+(3C(k+3,2)-C(k+4,2)), k=0..2",
      And @@ Table[Length[Ntet[k]] == dimNtet[k], {k, 0, 2}]];
  chk["RT_tet(0) = 4 (faces) , N_tet(0) = 6 (edges)  [Whitney lowest order]",
      Length[RTtet[0]] == 4 && Length[Ntet[0]] == 6];

  (* ---- independence (rank = dim) ---- *)
  chk["BDM_tet(1) linearly independent (rank = dim)",
      MatrixRank[(Flatten[Table[# /. pt, {pt, pts}]] & /@ BDMtet[1]), Tolerance -> 1*^-9] == dimBDMtet[1]];
  chk["RT_tet(1) linearly independent (rank = dim)",
      MatrixRank[(Flatten[Table[# /. pt, {pt, pts}]] & /@ RTtet[1]), Tolerance -> 1*^-9] == dimRTtet[1]];
  chk["N_tet(1) linearly independent (rank = dim)",
      MatrixRank[(Flatten[Table[# /. pt, {pt, pts}]] & /@ Ntet[1]), Tolerance -> 1*^-9] == dimNtet[1]];

  (* ---- S_k constraint:  p . x~ = 0 ---- *)
  chk["S_tet(1) satisfies p . (x,y,z) = 0 (rotational, radial-free)",
      And @@ (Simplify[#.{x, y, z}] === 0 & /@ Stet[1])];

  (* ---- de Rham: grad(P_{k+1}) c N_k, curl-free ---- *)
  Module[{grads},
   grads = gradd /@ monos3[2];                    (* P_2 -> grad c N_1 *)
   chk["grad(P_{k+1}) is curl-free", And @@ (Simplify[curlv[#]] === {0, 0, 0} & /@ grads)];
   rB = MatrixRank[(Flatten[Table[# /. pt, {pt, pts}]] & /@ Ntet[1]), Tolerance -> 1*^-9];
   rA = MatrixRank[(Flatten[Table[# /. pt, {pt, pts}]] & /@ Join[Ntet[1], grads]), Tolerance -> 1*^-9];
   chk["grad(P_{k+1}) lies IN span(N_tet(k))  (de Rham: grad c H(curl))", rA == rB]];

  (* ---- de Rham: curl(N_k) c BDM_k c RT_k, div-free ---- *)
  Module[{curls},
   curls = curlv /@ Ntet[1];
   chk["curl(N_tet) is div-free", And @@ (Simplify[divv[#]] === 0 & /@ curls)];
   rB = MatrixRank[(Flatten[Table[# /. pt, {pt, pts}]] & /@ BDMtet[1]), Tolerance -> 1*^-9];
   rA = MatrixRank[(Flatten[Table[# /. pt, {pt, pts}]] & /@ Join[BDMtet[1], curls]), Tolerance -> 1*^-9];
   chk["curl(N_tet) lies IN span(BDM_tet)", rA == rB];
   rB = MatrixRank[(Flatten[Table[# /. pt, {pt, pts}]] & /@ RTtet[1]), Tolerance -> 1*^-9];
   rA = MatrixRank[(Flatten[Table[# /. pt, {pt, pts}]] & /@ Join[RTtet[1], curls]), Tolerance -> 1*^-9];
   chk["curl(N_tet) lies IN span(RT_tet)  (de Rham: curl H(curl) c H(div))", rA == rB]];

  (* ---- BDM subset, divergence spaces, and shared loop kernel ---- *)
  p = 1;
  rB = MatrixRank[(Flatten[Table[# /. pt, {pt, pts}]] & /@ RTtet[p]), Tolerance -> 1*^-9];
  rA = MatrixRank[(Flatten[Table[# /. pt, {pt, pts}]] & /@ Join[RTtet[p], BDMtet[p]]), Tolerance -> 1*^-9];
  chk["BDM_tet(p) is a strict subspace of RT_tet(p)", rA == rB && dimBDMtet[p] < dimRTtet[p]];
  chk["div(BDM_tet(p)) spans P_(p-1), p=1..2",
      And @@ Table[
       MatrixRank[Table[Simplify[divv[v]] /. pt, {v, BDMtet[p]}, {pt, pts}], Tolerance -> 1*^-9] == dimPk3[p - 1],
       {p, 1, 2}]];
  chk["div(RT_tet(p)) spans P_p, p=1..2",
      And @@ Table[
       MatrixRank[Table[Simplify[divv[v]] /. pt, {v, RTtet[p]}, {pt, pts}], Tolerance -> 1*^-9] == dimPk3[p],
       {p, 1, 2}]];
  chk["BDM_tet(p) and RT_tet(p) have the same ker(div) dimension, p=1..3",
      And @@ Table[dimBDMtet[p] - dimPk3[p - 1] == dimRTtet[p] - dimPk3[p], {p, 1, 3}]];
  chk["BDM/RT tet normal traces both span P_p on a triangular face, p=1..2",
      And @@ Flatten[Table[
        MatrixRank[Table[Simplify[-v[[3]] /. z -> 0] /. pt, {v, family[p]}, {pt, pts2}],
                   Tolerance -> 1*^-9] == dimPk2[p],
        {p, 1, 2}, {family, {BDMtet, RTtet}}]]];
  chk["tet family ledger p=1: 12 BDM vs 15 RT, shared ker(div)=11",
      With[{q = HDivTetFamilyLedger[1]},
       q["BDMDoF"] == 12 && q["RTDoF"] == 15 && q["SharedKerDivDoF"] == 11 &&
       q["ExtraRTChargeDoF"] == 3]];
  chk["tet family ledger p=6: 252 BDM vs 280 RT, shared ker(div)=196",
      With[{q = HDivTetFamilyLedger[6]},
       q["BDMDoF"] == 252 && q["RTDoF"] == 280 && q["SharedKerDivDoF"] == 196 &&
       q["ExtraRTChargeDoF"] == 28]];

  (* ---- 2D triangle BDM / RT ---- *)
  chk["dim BDM_trig(p) = 2 C(p+2,2), p=1..3",
      And @@ Table[Length[BDMtrig[p]] == dimBDMtrig[p], {p, 1, 3}]];
  chk["dim RT_trig(k) = 2 C(k+2,2)+(k+1), k=0..2",
      And @@ Table[Length[RTtrig[k]] == dimRTtrig[k], {k, 0, 2}]];
  chk["RT_trig(0) = 3 edge funcs, div const (Whitney)",
      Length[RTtrig[0]] == 3 && And @@ (FreeQ[Simplify[div2[#[[1 ;; 2]]]], x | y] & /@ RTtrig[0])];
  chk["div(BDM_trig(p)) spans P_(p-1), p=1..2",
      And @@ Table[
       MatrixRank[Table[Simplify[div2[v[[1 ;; 2]]]] /. pt, {v, BDMtrig[p]}, {pt, pts2}],
                  Tolerance -> 1*^-9] == dimPk2[p - 1], {p, 1, 2}]];
  chk["div(RT_trig(p)) spans P_p, p=1..2",
      And @@ Table[
       MatrixRank[Table[Simplify[div2[v[[1 ;; 2]]]] /. pt, {v, RTtrig[p]}, {pt, pts2}],
                  Tolerance -> 1*^-9] == dimPk2[p], {p, 1, 2}]];
  chk["BDM_trig(p) and RT_trig(p) have the same ker(div) dimension, p=1..3",
      And @@ Table[dimBDMtrig[p] - dimPk2[p - 1] == dimRTtrig[p] - dimPk2[p], {p, 1, 3}]];
  chk["BDM/RT trig normal traces both span P_p on an edge, p=1..3",
      And @@ Flatten[Table[
        MatrixRank[Table[Simplify[-v[[2]] /. y -> 0] /. {x -> xx}, {v, family[p]},
                         {xx, {.1, .25, .45, .7, .9}}], Tolerance -> 1*^-9] == p + 1,
        {p, 1, 3}, {family, {BDMtrig, RTtrig}}]]];
  chk["trig family ledger p=1: 6 BDM vs 8 RT, shared ker(div)=5",
      With[{q = HDivTrigFamilyLedger[1]},
       q["BDMDoF"] == 6 && q["RTDoF"] == 8 && q["SharedKerDivDoF"] == 5 &&
       q["ExtraRTChargeDoF"] == 2]];
  chk["trig family ledger p=6: 56 BDM vs 63 RT, shared ker(div)=35",
      With[{q = HDivTrigFamilyLedger[6]},
       q["BDMDoF"] == 56 && q["RTDoF"] == 63 && q["SharedKerDivDoF"] == 35 &&
       q["ExtraRTChargeDoF"] == 7]];

  Print[If[ok, "ALL PASS", "SOME FAILED"]];
];
