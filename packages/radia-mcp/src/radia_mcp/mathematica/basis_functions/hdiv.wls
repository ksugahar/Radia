#!/usr/bin/env wolframscript
(* ::Package:: *)

(*
  hdiv.wls -- NGSolve-style H(div) (Raviart-Thomas / BDM) shape functions, clean wls.
  THE space the HDiv-VIM solver needs: the magnetization face-normal density lives
  in H(div), and the de Rham property  div(curl(H(curl))) = 0  is exactly what makes
  the div-free magnetization circulations ("loops", ker N) field-null by construction
  -- the machine-zero loops, replacing hand-built null vectors / element engineering.

  Two complementary layers, both ported from NGSolve fem/hdivhofe_impl.hpp:

   (1) Schoeberl-Zaglmayr COVARIANT helpers (Du_Cross_Dv, ...): the exact building
       blocks NGSolve uses.  Each is a cross product of gradients => a polynomial =>
       integrable against the 1/r VIM kernel.  Their divergence identities are the
       FEEC content and are self-tested directly.

   (2) The TENSOR de Rham construction on the cube (the HEX, = VIM main target):
       H1 = W(x)W(x)W,  H(curl) = Q(x)W(x)W (+cyc),  H(div) = W(x)Q(x)Q (+cyc),
       L2 = Q(x)Q(x)Q,  with the 1D exact sequence  d: W -> Q  (W = integrated
       Legendre / H1, Q = Legendre / L2).  This spans NGSolve's RT_p hex space
       (HDiv(mesh,order=p)); a clean explicit-polynomial basis, ideal for the
       analytical VIM field operator.  The de Rham exactness is SHOWN: curl of the
       H(curl)-hex lands in H(div)-hex and is div-free to symbolic zero.

  Simplex RT0 (tet face / trig edge) = the Whitney W2 / rotated-W1 (cf. derham.wls).
  High-order simplex (integrated-Jacobi / Dubiner orthopols) = next increment.

  Reference cube/sgmt on [0,1]^d ; tet lami = {x,y,z,1-x-y-z}.
*)

(* ---- polynomial foundation (= recursive_pol.wls) ---- *)
l[i_, s_] := LegendreP[i, s];
L[i_, s_] := L[i, s] = Simplify[Integrate[l[i - 1, u], {u, -1, s}]];   (* integrated Legendre *)

(* 1D exact sequence on [0,1]:  d/dt : W(H1) -> Q(L2), surjective, ker = const. *)
oneDH1[q_, t_] := Join[{1 - t, t}, Table[L[i, 2 t - 1], {i, 2, q}]];   (* W: dim q+1, deg<=q *)
oneDL2[q_, t_] := Table[LegendreP[j, 2 t - 1], {j, 0, q - 1}];         (* Q: dim q,   deg<=q-1 *)

(* ============ Schoeberl-Zaglmayr covariant helpers (NGSolve) ============ *)
gradd[f_] := D[f, {{x, y, z}}];
divv[V_]  := D[V[[1]], x] + D[V[[2]], y] + D[V[[3]], z];
curlv[V_] := Curl[V, {x, y, z}];

(* div-free:  Cross[grad u, grad v] = curl(u grad v) *)
DuCrossDv[u_, v_]      := Cross[gradd[u], gradd[v]];
wDuCrossDv[u_, v_, w_] := w Cross[gradd[u], gradd[v]];                 (* div = grad w . (gu x gv) *)
(* RT0 cyclic face function (= Whitney W2 / 2):  div = 3 det[gu,gv,gw] = const *)
uDvDwCyclic[a_, b_, c_] := a Cross[gradd[b], gradd[c]] +
                          b Cross[gradd[c], gradd[a]] + c Cross[gradd[a], gradd[b]];
(* non-div-free "Type 2":  (u grad v - v grad u) x grad w ; div = 2 grad w . (gu x gv) *)
uDvDwMinusDuvDw[u_, v_, w_] := u Cross[gradd[v], gradd[w]] - v Cross[gradd[u], gradd[w]];
(* div-free face/cell bubble: curl of the H(curl) fn  w (u grad v - v grad u) *)
curlUDvwMinusDuvw[u_, v_, w_] := Curl[w (u gradd[v] - v gradd[u]), {x, y, z}];

(* ===================== HEX (tensor de Rham, RT_p) ====================== *)
(* x-family: W(x) (x) Q(y) (x) Q(z) e_x ; +cyclic.  dim = 3 (p+2)(p+1)^2. *)
HDivHex[p_] := Module[{Wf, Qf, xf, yf, zf},
  Wf[t_] := oneDH1[p + 1, t]; Qf[t_] := oneDL2[p + 1, t];
  xf = Flatten[Table[{Wf[x][[a]] Qf[y][[b]] Qf[z][[c]], 0, 0},
       {a, p + 2}, {b, p + 1}, {c, p + 1}], 2];
  yf = Flatten[Table[{0, Qf[x][[a]] Wf[y][[b]] Qf[z][[c]], 0},
       {a, p + 1}, {b, p + 2}, {c, p + 1}], 2];
  zf = Flatten[Table[{0, 0, Qf[x][[a]] Qf[y][[b]] Wf[z][[c]]},
       {a, p + 1}, {b, p + 1}, {c, p + 2}], 2];
  Join[xf, yf, zf]];

(* H(curl) hex, x-family only (enough to SHOW curl(H(curl)) c H(div), div-free):
   Q(x) (x) W(y) (x) W(z) e_x. *)
HCurlHexXfam[p_] := Module[{Wf, Qf},
  Wf[t_] := oneDH1[p + 1, t]; Qf[t_] := oneDL2[p + 1, t];
  Flatten[Table[{Qf[x][[a]] Wf[y][[b]] Wf[z][[c]], 0, 0},
       {a, p + 1}, {b, p + 2}, {c, p + 2}], 2]];

(* ===================== QUAD (2D tensor de Rham, RT_p) ================== *)
div2[V_]  := D[V[[1]], x] + D[V[[2]], y];
HDivQuad[p_] := Module[{Wf, Qf, xf, yf},
  Wf[t_] := oneDH1[p + 1, t]; Qf[t_] := oneDL2[p + 1, t];
  xf = Flatten[Table[{Wf[x][[a]] Qf[y][[b]], 0}, {a, p + 2}, {b, p + 1}], 1];
  yf = Flatten[Table[{0, Qf[x][[a]] Wf[y][[b]]}, {a, p + 1}, {b, p + 2}], 1];
  Join[xf, yf]];
(* 2D H1 (for the rot-grad showcase): W(x) (x) W(y). *)
H1QuadT[p_] := Module[{Wf}, Wf[t_] := oneDH1[p + 1, t];
  Flatten[Table[Wf[x][[a]] Wf[y][[b]], {a, p + 2}, {b, p + 2}]]];

(* ===================== SIMPLEX RT0 (Whitney) ========================== *)
(* TET face RT0 = uDvDwCyclic(la,lb,lc) = Whitney W2 / 2.  4 faces, div const. *)
lamTet = {x, y, z, 1 - x - y - z};
tetFaces = {{1, 2, 3}, {1, 2, 4}, {1, 3, 4}, {2, 3, 4}};
HDivTetRT0[] := uDvDwCyclic[lamTet[[#[[1]]]], lamTet[[#[[2]]]], lamTet[[#[[3]]]]] & /@ tetFaces;

(* TRIG edge RT0 = 90deg rotation of the 2D Whitney W1.  3 edges, div const. *)
lamTrig = {1 - x - y, x, y};
gl2 = D[#, {{x, y}}] & /@ lamTrig;
rot90[{u_, v_}] := {v, -u};
W1trig[{a_, b_}] := lamTrig[[a]] gl2[[b]] - lamTrig[[b]] gl2[[a]];
trigEdges = {{1, 2}, {1, 3}, {2, 3}};
HDivTrigRT0[] := rot90[W1trig[#]] & /@ trigEdges;

(* ============================ self-test ============================ *)
evalRow[f_, pts_] := Flatten[Table[f /. pt, {pt, pts}]];
valsMat[funcs_, pts_] := evalRow[#, pts] & /@ funcs;

Module[{ok = True, chk, u, v, w, idDF, idCyc, idT2, pts, rB, rA, p, dims, divs},
  chk[name_, cond_] := (Print["  ", If[cond, "ok  ", "FAIL"], name]; ok = ok && cond);
  Print["hdiv.wls self-test:"];

  (* ---- (1) SZ covariant-helper divergence identities (the FEEC content) ---- *)
  u = x^2 y + z; v = y z + x; w = x + 2 y + 3 z + x y;
  chk["Du_Cross_Dv is div-free (curl of u grad v)",
      Simplify[divv[DuCrossDv[u, v]]] === 0];
  chk["uDvDwCyclic div = 3 det[gu,gv,gw] (const for affine) ",
      With[{a = x, b = y, c = z},
       Simplify[divv[uDvDwCyclic[a, b, c]] - 3 Det[{gradd[a], gradd[b], gradd[c]}]] === 0]];
  chk["uDvDwMinusDuvDw div = 2 grad w . (gu x gv)",
      Simplify[divv[uDvDwMinusDuvDw[u, v, w]] - 2 gradd[w].Cross[gradd[u], gradd[v]]] === 0];
  chk["curlUDvwMinusDuvw is div-free (it is a curl)",
      Simplify[divv[curlUDvwMinusDuvw[u, v, w]]] === 0];

  (* ---- (2) 1D exact sequence  d: W -> Q  surjective, ker = const ---- *)
  chk["d(oneDH1[q]) spans oneDL2[q]  (q=1..5)",
      And @@ Table[
        With[{dW = D[oneDH1[q, t], t], Q = oneDL2[q, t]},
         MatrixRank[Table[Join[dW, Q] /. t -> tt, {tt, RandomReal[{0, 1}, 12]}],
                    Tolerance -> 1*^-9] == q], {q, 1, 5}]];

  (* ---- (3) HEX tensor RT_p: dimension, div in L2, RT0 normal trace ---- *)
  chk["dim HDivHex(p) = 3 (p+2)(p+1)^2, p=0..2",
      And @@ Table[Length[HDivHex[p]] == 3 (p + 2) (p + 1)^2, {p, 0, 2}]];
  chk["HDivHex(0) = 6 RT0 face funcs, each div = const",
      Length[HDivHex[0]] == 6 &&
      And @@ (FreeQ[Simplify[divv[#]], x | y | z] & /@ HDivHex[0])];
  SeedRandom[11];
  pts = Table[{x -> RandomReal[], y -> RandomReal[], z -> RandomReal[]}, {30}];
  p = 1;
  divs = Simplify[divv[#]] & /@ HDivHex[p];
  chk["div(HDivHex(1)) spans L2 = Q(x)Q(x)Q (rank (p+1)^3 = 8)",
      MatrixRank[Table[d /. pt, {d, divs}, {pt, pts}], Tolerance -> 1*^-9] == (p + 1)^3];

  (* ---- (4) de Rham SHOWCASE: curl(H(curl)-hex) c H(div)-hex AND div-free ---- *)
  Module[{curls},
   curls = curlv /@ HCurlHexXfam[1];
   chk["curl(HCurlHex x-family) is div-free (machine-zero loops)",
       And @@ (Simplify[divv[#]] === 0 & /@ curls)];
   rB = MatrixRank[valsMat[HDivHex[1], pts], Tolerance -> 1*^-9];
   rA = MatrixRank[valsMat[Join[HDivHex[1], curls], pts], Tolerance -> 1*^-9];
   chk["curl(HCurlHex) lies IN span(HDivHex)  (de Rham: curl H(curl) c H(div))",
       rA == rB]];

  (* ---- (5) QUAD 2D RT_p + 2D rot-grad showcase ---- *)
  chk["dim HDivQuad(p) = 2 (p+2)(p+1), p=0..2",
      And @@ Table[Length[HDivQuad[p]] == 2 (p + 2) (p + 1), {p, 0, 2}]];
  chk["HDivQuad(0) = 4 edge funcs, div const",
      Length[HDivQuad[0]] == 4 &&
      And @@ (FreeQ[Simplify[div2[#]], x | y] & /@ HDivQuad[0])];
  Module[{rg, pts2, ra, rb},
   (* 2D de Rham: rot(grad s) = (d_y s, -d_x s) c H(div), div-free *)
   rg = ({D[#, y], -D[#, x]} & /@ H1QuadT[1]);
   chk["rot-grad(H1Quad) is div-free (2D)", And @@ (Simplify[div2[#]] === 0 & /@ rg)];
   SeedRandom[12]; pts2 = Table[{x -> RandomReal[], y -> RandomReal[]}, {20}];
   rb = MatrixRank[(Flatten[Table[# /. pt, {pt, pts2}]] & /@ HDivQuad[1]), Tolerance -> 1*^-9];
   ra = MatrixRank[(Flatten[Table[# /. pt, {pt, pts2}]] & /@ Join[HDivQuad[1], rg]), Tolerance -> 1*^-9];
   chk["rot-grad(H1Quad) c span(HDivQuad)  (2D de Rham)", ra == rb]];

  (* ---- (6) simplex RT0 (Whitney): tet faces / trig edges ---- *)
  chk["HDivTetRT0 = 4 face funcs, each div = const (Whitney W2)",
      Length[HDivTetRT0[]] == 4 &&
      And @@ (FreeQ[Simplify[divv[#]], x | y | z] & /@ HDivTetRT0[])];
  chk["HDivTrigRT0 = 3 edge funcs, each div = const (rot Whitney W1)",
      Length[HDivTrigRT0[]] == 3 &&
      And @@ (FreeQ[Simplify[div2[#]], x | y] & /@ HDivTrigRT0[])];
  (* RT0 invariant: each face fn has constant divergence of EQUAL MAGNITUDE (= 3 |det g|).
     The SIGN is the face-orientation convention (NGSolve's fav-sort + fop pick the outward
     normal for inter-element conformity); it does not change the single-element space. *)
  chk["HDivTetRT0 four |divergences| are equal nonzero constants (RT0)",
      With[{dd = Abs[Simplify[divv[#]]] & /@ HDivTetRT0[]},
       (Equal @@ dd) && dd[[1]] =!= 0]];

  Print[If[ok, "ALL PASS", "SOME FAILED"]];
];
