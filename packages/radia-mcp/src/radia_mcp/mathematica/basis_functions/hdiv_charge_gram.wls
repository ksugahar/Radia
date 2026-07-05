#!/usr/bin/env wolframscript
(* ::Package:: *)

(*
  hdiv_charge_gram.wls -- closed-form building blocks for the RT1 hex charge Coulomb Gram (the HDiv-VIM demag
  operator N = B^T G B, with G_ab = INT INT c_a(x) c_b(y)/|x-y| over box hexes; c = div of the RT1 basis).

  WHY: the C++ RT1 hex ChargeGram carries a tiny (~1e-5) NEAR-PAIR reflection-symmetry defect -- a symmetric
  cube full solve gets a spurious transverse M -- because the ref-fixed hex->6-tet split (body diagonal 0-6) +
  the Duffy near grading are reflection-asymmetric.  A symmetric-DECOMPOSITION quadrature fix costs ~16x
  (QuadBlockHex is a 24x24 sub-tet double loop).  This module is the CHEAP, EXACT alternative:

    inner  INT_B c_b(y)/|x-y| dy   = CLOSED FORM  (box Newtonian potential; corner sum of a scalar antideriv)
    outer  INT_A c_a(x) [...] dx   = SYMMETRIC tensor Gauss-Legendre over box A  (smooth integrand now)

  The symmetric outer rule commutes with the box reflections, so a symmetry-forbidden (transverse) Gram entry
  is machine-zero BY CONSTRUCTION -- no 6-tet, no Duffy, no 16x.  Since the VIM is RT1-only the charges are the
  small Q1 = Q(x)Q(y)Q(z) space (see hdiv.wls HDivHex[1]), so the closed forms stay tractable.  Companion of
  the analytical_formulas closed-form-reference layer.  (Fully-symbolic OUTER = the 6-fold antiderivative of
  1/r; Integrate times out, so the outer stays a symmetric numeric rule -- which is all the exact symmetry
  needs.)
*)

(* ---------- 3-fold antiderivatives of 1/r and its first moments (d^3 F = weight/r) ---------- *)
kern[dx_, dy_, dz_] := 1/Sqrt[dx^2 + dy^2 + dz^2];
F3  = Integrate[1/Sqrt[u^2 + v^2 + w^2], w, u, v];   (* d^3 F3  = 1/r *)
Fu  = Integrate[u/Sqrt[u^2 + v^2 + w^2], w, u, v];   (* d^3 Fu  = u/r *)
Fv  = Integrate[v/Sqrt[u^2 + v^2 + w^2], w, u, v];   (* d^3 Fv  = v/r *)
Fw  = Integrate[w/Sqrt[u^2 + v^2 + w^2], w, u, v];   (* d^3 Fw  = w/r *)
SetAttributes[{scalarF}, HoldFirst];
evalF[Fexpr_, a_, b_, c_] := Fexpr /. {u -> a, v -> b, w -> c};

(* ---------- weighted box potential  INT_B  wpoly(y) / |X-y| dy , box B = {{bx0,bx1},{by0,by1},{bz0,bz1}} ----------
   For a Q1 weight wpoly(y) = prod of linear factors, expand into the constant + first moments (X-y) and use the
   matching antiderivative; the box integral is the 8-corner inclusion-exclusion sum (the -1 is the (-1)^3 from
   du = -dy).  Here we expose the const and the three (X_j - y_j) first-moment potentials -- the pieces any Q1
   charge needs. *)
cornerSum[Fexpr_, X_, box_] := Sum[
    (If[i == 1, 1, -1] If[j == 1, 1, -1] If[k == 1, 1, -1]) *
      evalF[Fexpr, X[[1]] - box[[1, If[i == 1, 2, 1]]], X[[2]] - box[[2, If[j == 1, 2, 1]]],
                   X[[3]] - box[[3, If[k == 1, 2, 1]]]],
    {i, 0, 1}, {j, 0, 1}, {k, 0, 1}];
boxPotConst[X_, box_] := -cornerSum[F3, X, box];                 (* INT_B 1/|X-y| dy *)
boxPotMomX[X_, box_]  := -cornerSum[Fu, X, box];                 (* INT_B (X1-y1)/|X-y| dy *)
boxPotMomY[X_, box_]  := -cornerSum[Fv, X, box];                 (* INT_B (X2-y2)/|X-y| dy *)
boxPotMomZ[X_, box_]  := -cornerSum[Fw, X, box];                 (* INT_B (X3-y3)/|X-y| dy *)
(* the y-linear charge (y2 - cy) potential: INT_B (y2-cy)/|X-y| = (X2-cy) boxPotConst - boxPotMomY *)
boxPotYlin[X_, box_, cy_] := (X[[2]] - cy) boxPotConst[X, box] - boxPotMomY[X, box];

(* ---------- semi-closed Gram entry: G = INT_A sigmaA(x) * (weighted box potential of B)(x) dx,
   via a SYMMETRIC tensor Gauss-Legendre over box A (nodes symmetric about A's center -> exact reflection sym) ---------- *)
Needs["NumericalDifferentialEquationAnalysis`"];
gaussNodes1D[n_, {a_, b_}] := N[GaussianQuadratureWeights[n, a, b], 30];
chargeGramSemiclosed[sigmaA_, phiB_, A_, nGauss_: 8] := Module[{gx, gy, gz},
  gx = gaussNodes1D[nGauss, A[[1]]]; gy = gaussNodes1D[nGauss, A[[2]]]; gz = gaussNodes1D[nGauss, A[[3]]];
  Sum[gx[[i, 2]] gy[[j, 2]] gz[[k, 2]] *
      sigmaA[{gx[[i, 1]], gy[[j, 1]], gz[[k, 1]]}] * phiB[{gx[[i, 1]], gy[[j, 1]], gz[[k, 1]]}],
      {i, 1, nGauss}, {j, 1, nGauss}, {k, 1, nGauss}]];

(* ---------------------------------- self-test ---------------------------------- *)
If[MemberQ[$ScriptCommandLine, "-selftest"] || True,
  Module[{ok = True, chk, dpt, A, B, Gtrans, Gconst},
   chk[name_, cond_] := (Print["  ", If[cond, "ok  ", "FAIL"], name]; ok = ok && cond);
   Print["hdiv_charge_gram.wls self-test:"];
   dpt = {0.4, -0.7, 1.3};
   chk["d^3 F3 == 1/r", Abs[(D[F3, u, v, w] /. Thread[{u, v, w} -> dpt]) - kern @@ dpt] < 1*^-9];
   chk["d^3 Fv == v/r", Abs[(D[Fv, u, v, w] /. Thread[{u, v, w} -> dpt]) - dpt[[2]] kern @@ dpt] < 1*^-9];
   (* box potential closed vs numeric, point above the unit box *)
   With[{X = {0.3, 0.55, 1.72}},
    chk["box potential closed == numeric",
     Abs[Re@N@boxPotConst[X, {{0, 1}, {0, 1}, {0, 1}}] -
        NIntegrate[kern[X[[1]] - y1, X[[2]] - y2, X[[3]] - y3], {y1, 0, 1}, {y2, 0, 1}, {y3, 0, 1}]] < 1*^-6]];
   (* the stacked 1x1x2 pair: A above, B below, share z=0 *)
   A = {{0, 1}, {0, 1}, {0, 1}}; B = {{0, 1}, {0, 1}, {-1, 0}};
   Gtrans = chargeGramSemiclosed[(1 &), boxPotYlin[#, B, 1/2] &, A, 8];   (* symmetry-forbidden -> 0 *)
   Gconst = chargeGramSemiclosed[(1 &), boxPotConst[#, B] &, A, 8];       (* genuine box-box Coulomb *)
   chk["TRANSVERSE G[const_A, ylin_B] is machine-zero (C++ quadrature ~1e-5)", Abs[Re@N@Gtrans] < 1*^-14];
   chk["nonzero  G[const_A, const_B] == numeric box-box Coulomb (~0.98088)",
       Abs[Re@N@Gconst - 0.9808850] < 1*^-5];
   Print[If[ok, "ALL PASS", "FAILURES ABOVE"]];
  ]];
