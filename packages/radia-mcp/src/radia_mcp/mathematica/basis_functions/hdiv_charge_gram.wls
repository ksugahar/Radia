#!/usr/bin/env wolframscript
(* ::Package:: *)

(*
  hdiv_charge_gram.wls -- closed-form building blocks for the production BDM1 hex charge Coulomb Gram (the
  HDiv-VIM demag operator N = B^T G B, with G_ab = INT INT c_a(x) c_b(y)/|x-y| over box hexes).

  WHY: the C++ BDM1 hex ChargeGram carries a tiny (~1e-5) NEAR-PAIR reflection-symmetry defect -- a symmetric
  cube full solve gets a spurious transverse M -- because the ref-fixed hex->6-tet split (body diagonal 0-6) +
  the Duffy near grading are reflection-asymmetric.  A symmetric-DECOMPOSITION quadrature fix costs ~16x
  (QuadBlockHex is a 24x24 sub-tet double loop).  This module is the CHEAP, EXACT alternative:

    inner  INT_B c_b(y)/|x-y| dy   = CLOSED FORM  (box Newtonian potential; corner sum of a scalar antideriv)
    outer  INT_A c_a(x) [...] dx   = SYMMETRIC tensor Gauss-Legendre over box A  (smooth integrand now)

  The symmetric outer rule commutes with the box reflections, so a symmetry-forbidden (transverse) Gram entry
  is machine-zero BY CONSTRUCTION -- no 6-tet, no Duffy, no 16x.  This derivation targets the BDM1 production
  charge image, the small Q1 = Q(x)Q(y)Q(z) trilinear space.  The tensor RT certificate in hdiv.wls has the
  same Q1 divergence target (HDivHex[1] has 21 vector polynomials spanning
  the 8-dim monomial basis {1,x,y,z,xy,xz,yz,xyz}), so the closed forms stay tractable.  Companion of the
  analytical_formulas closed-form-reference layer.

  SCOPE (verified 2026-07-05): the closed form is EXACT for DISJOINT box pairs -- the near/far blocks, i.e.
  exactly where the C++ ~1e-5 defect lives (exterior potential matches NIntegrate to ~6e-9; the transverse
  exact-symmetry entries are machine-zero).  The COINCIDENT / self block (A=B) is NOT covered by the naive
  corner-sum: the Newtonian-potential antiderivative (ArcSinh + ArcTan) picks up interior branch offsets
  (measured: closed - numeric = 3*pi/2 at the box center), so the self potential comes out wrong.  This is the
  known Waldvogel interior-branch subtlety and is left OUT OF SCOPE because the C++ self Gram is ALREADY EXACT
  (3.66e-15, decomposition-independent same-hex integral) -- the fix only needs the disjoint near blocks.  If a
  fully-closed self block is ever wanted, add the octant-dependent pi branch corrections to F3s below.
  (Fully-symbolic OUTER = the 6-fold antiderivative of 1/r; Integrate times out, so the outer stays a symmetric
  numeric rule -- which is all the exact symmetry needs.)
*)

Needs["NumericalDifferentialEquationAnalysis`"];
kern[dx_, dy_, dz_] := 1/Sqrt[dx^2 + dy^2 + dz^2];

(* ---------- the 8 three-fold antiderivatives of (prod_{k in T} u_k)/r, T subset {1,2,3}  (d^3 F_T = mono_T/r).
   T = {} -> 1/r, {2} -> v/r, {1,3} -> u w/r, {1,2,3} -> u v w/r, etc.  These span the Q1 charge space. ---------- *)
coord = {u, v, w};
mono[T_] := Times @@ coord[[T]];
subsets3 = Subsets[{1, 2, 3}];
anti = Association[];
Do[anti[T] = Integrate[mono[T]/Sqrt[u^2 + v^2 + w^2], w, u, v], {T, subsets3}];
evalF[Fexpr_, a_, b_, c_] := Fexpr /. {u -> a, v -> b, w -> c};

(* ---------- weighted box potential over box B = {{bx0,bx1},{by0,by1},{bz0,bz1}} at eval point X ----------
   A Q1 (trilinear) charge c(y) = prod_{k in S} y_k is expanded in u = X - y and integrated with the matching
   3-fold antiderivative; each box integral is the 8-corner inclusion-exclusion sum (the leading -1 is (-1)^3
   from du = -dy).  Valid for DISJOINT pairs (X outside B) -- the near/far blocks; see SCOPE re the self block. *)
cornerSum[Fexpr_, X_, box_] := Sum[
    (If[i == 1, 1, -1] If[j == 1, 1, -1] If[k == 1, 1, -1]) *
      evalF[Fexpr, X[[1]] - box[[1, If[i == 1, 2, 1]]], X[[2]] - box[[2, If[j == 1, 2, 1]]],
                   X[[3]] - box[[3, If[k == 1, 2, 1]]]],
    {i, 0, 1}, {j, 0, 1}, {k, 0, 1}];
(* INT_B (prod_{k in T} u_k)/|X-y| dy,  u = X - y : *)
boxPotMom[T_, X_, box_] := -cornerSum[anti[Sort[T]], X, box];
(* potential of the y-monomial charge c(y) = prod_{k in S} y_k :
   INT_B (prod_{k in S} y_k)/|X-y| dy = sum_{T subset S} (-1)^|T| (prod_{k in S\T} X_k) boxPotMom[T]   (y_k = X_k - u_k) *)
boxPotCharge[S_, X_, box_] := Sum[
    (-1)^Length[T] (Times @@ (X[[#]] & /@ Complement[S, T])) boxPotMom[T, X, box],
    {T, Subsets[S]}];
(* convenience wrappers (classic scalar box potential + a centered y-linear charge) *)
boxPotConst[X_, box_] := boxPotCharge[{}, X, box];                       (* INT_B 1/|X-y| dy *)
boxPotYlin[X_, box_, cy_] := boxPotCharge[{2}, X, box] - cy boxPotCharge[{}, X, box];   (* INT_B (y2-cy)/|X-y| dy *)

(* ---------- semi-closed Gram entry: G = INT_A c_a(x) * (box potential of c_b over B)(x) dx,
   via a SYMMETRIC tensor Gauss-Legendre over box A (nodes symmetric about A's center -> exact reflection sym).
   Sa, Sb are index-subsets of {1,2,3}: the Q1 monomial charges prod_{k in Sa} x_k and prod_{k in Sb} y_k. ---------- *)
gaussNodes1D[n_, {a_, b_}] := N[GaussianQuadratureWeights[n, a, b], 30];
chargeGramSemiclosed[Sa_, Sb_, A_, B_, nGauss_: 8] := Module[{gx, gy, gz},
  gx = gaussNodes1D[nGauss, A[[1]]]; gy = gaussNodes1D[nGauss, A[[2]]]; gz = gaussNodes1D[nGauss, A[[3]]];
  Sum[gx[[i, 2]] gy[[j, 2]] gz[[k, 2]] *
      (Times @@ ({gx[[i, 1]], gy[[j, 1]], gz[[k, 1]]}[[Sa]])) *
      boxPotCharge[Sb, {gx[[i, 1]], gy[[j, 1]], gz[[k, 1]]}, B],
      {i, 1, nGauss}, {j, 1, nGauss}, {k, 1, nGauss}]];

(* branch-safe symmetric Waldvogel antideriv of 1/r (exterior EXACT ~1e-9; interior needs +pi octant
   corrections -- see SCOPE).  Kept for reference / future self-block completion. *)
rr[a_, b_, c_] := Sqrt[a^2 + b^2 + c^2];
F3s[a_, b_, c_] := Module[{r = rr[a, b, c]},
   a b ArcSinh[c/Sqrt[a^2 + b^2]] + b c ArcSinh[a/Sqrt[b^2 + c^2]] + c a ArcSinh[b/Sqrt[c^2 + a^2]]
   - a^2/2 ArcTan[a r, b c] - b^2/2 ArcTan[b r, c a] - c^2/2 ArcTan[c r, a b]];

(* ---------------------------------- self-test ---------------------------------- *)
If[MemberQ[$ScriptCommandLine, "-selftest"] || True,
  Module[{ok = True, chk, dpt, A, B, lab},
   chk[name_, cond_] := (Print["  ", If[cond, "ok  ", "FAIL"], name]; ok = ok && cond);
   lab[S_] := If[S === {}, "1", StringJoin[({"x", "y", "z"}[[#]] & /@ S)]];
   Print["hdiv_charge_gram.wls self-test (full Q1 basis, disjoint blocks):"];

   (* (1) all 8 antiderivatives: numeric d^3 F_T == mono_T / r *)
   dpt = {0.4, -0.7, 1.3};
   Do[chk["d^3 F_" <> ToString[T] <> " == mono/r",
      Abs[(D[anti[T], u, v, w] /. Thread[{u, v, w} -> dpt]) - (Times @@ dpt[[T]]) (kern @@ dpt)] < 1*^-8],
    {T, subsets3}];

   (* stacked disjoint pair: A above, B below, sharing z=0 *)
   A = {{0, 1}, {0, 1}, {0, 1}}; B = {{0, 1}, {0, 1}, {-1, 0}};

   (* (2) exterior box potential closed == numeric -- validity of the closed inner integral on DISJOINT pairs *)
   With[{X = {0.3, 0.55, 1.72}},
    chk["exterior box potential closed == numeric (<1e-6)",
     Abs[Re@N@boxPotConst[X, {{0, 1}, {0, 1}, {0, 1}}] -
        NIntegrate[kern[X[[1]] - y1, X[[2]] - y2, X[[3]] - y3], {y1, 0, 1}, {y2, 0, 1}, {y3, 0, 1},
          PrecisionGoal -> 8, MaxRecursion -> 20]] < 1*^-6]];

   (* (3) closed-form Gauss-convergence + anchored value for the genuine const-const box-box Coulomb *)
   chk["disjoint G[1_A,1_B] Gauss-converged (n=8 vs 16 < 1e-6)",
      Abs[Re@N@chargeGramSemiclosed[{}, {}, A, B, 8] - Re@N@chargeGramSemiclosed[{}, {}, A, B, 16]] < 1*^-6];
   chk["disjoint G[1_A,1_B] == 0.9808850 (anchored box-box Coulomb)",
      Abs[Re@N@chargeGramSemiclosed[{}, {}, A, B, 8] - 0.9808850] < 1*^-5];

   (* (3b) FULL Q1 basis vs numeric on a SEPARATED disjoint pair -- validates the DOUBLE/TRIPLE moment charges
      (xy, xz, yz, xyz), the entries not exercised by the const/single-moment checks above.  Separated (a z-gap
      between A and Bs) so the 6-fold NIntegrate reference converges cleanly (no touching-face slowdown). *)
   With[{Bs = {{0, 1}, {0, 1}, {-2, -1}},
         gN = Function[{Sa, Sb, BB},
            NIntegrate[(Times @@ ({x1, x2, x3}[[Sa]])) (Times @@ ({y1, y2, y3}[[Sb]])) *
               kern[x1 - y1, x2 - y2, x3 - y3],
               {x1, 0, 1}, {x2, 0, 1}, {x3, 0, 1},
               {y1, BB[[1, 1]], BB[[1, 2]]}, {y2, BB[[2, 1]], BB[[2, 2]]}, {y3, BB[[3, 1]], BB[[3, 2]]},
               PrecisionGoal -> 7, MaxRecursion -> 16]]},
    Do[chk["closed==num (separated)  G[" <> lab[p[[1]]] <> "_A, " <> lab[p[[2]]] <> "_B]",
        Abs[Re@N@chargeGramSemiclosed[p[[1]], p[[2]], A, Bs, 10] - gN[p[[1]], p[[2]], Bs]] < 1*^-5],
     {p, {{{}, {1, 2}}, {{}, {1, 2, 3}}, {{1, 2}, {1, 3}}, {{3}, {2, 3}}}}]];

   (* (4) EXACT reflection symmetry: a charge odd about a box center pairs to MACHINE ZERO -- purely from the
      closed form.  These are the entries the C++ 6-tet+Duffy quadrature returns as the spurious ~1e-5. *)
   chk["EXACT-SYM  G[1_A,(y-1/2)_B] machine-zero",
      Abs[Re@N@(chargeGramSemiclosed[{}, {2}, A, B, 8] - (1/2) chargeGramSemiclosed[{}, {}, A, B, 8])] < 1*^-13];
   chk["EXACT-SYM  G[x_A,(y-1/2)_B] machine-zero",
      Abs[Re@N@(chargeGramSemiclosed[{1}, {2}, A, B, 8] - (1/2) chargeGramSemiclosed[{1}, {}, A, B, 8])] < 1*^-13];
   chk["EXACT-SYM  G[(x-1/2)_A,1_B] machine-zero",
      Abs[Re@N@(chargeGramSemiclosed[{1}, {}, A, B, 8] - (1/2) chargeGramSemiclosed[{}, {}, A, B, 8])] < 1*^-13];
   chk["EXACT-SYM  G[z_A,(y-1/2)_B] machine-zero",
      Abs[Re@N@(chargeGramSemiclosed[{3}, {2}, A, B, 8] - (1/2) chargeGramSemiclosed[{3}, {}, A, B, 8])] < 1*^-13];

   Print[If[ok, "ALL PASS (disjoint near/far blocks; self block is C++-exact, see SCOPE)", "FAILURES ABOVE"]];
  ]];
