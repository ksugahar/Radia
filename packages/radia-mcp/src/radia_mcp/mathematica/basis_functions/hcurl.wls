#!/usr/bin/env wolframscript
(* ::Package:: *)

(*
  hcurl.wls -- NGSolve-style H(curl) (Nedelec) shape functions, clean wls.
  H(curl) is the MIDDLE of the magnetostatic de Rham complex and the home of the
  "loops": on a contractible element

        H1  --grad-->  H(curl)  --curl-->  H(div)  --div-->  L2

  is EXACT, so  ker(curl|H(curl)) = grad(H1)  -- every curl-free edge field is a
  gradient.  The div-free magnetization circulations the HDiv-VIM solver must handle
  ("loops", ker N c H(div)) are precisely  curl(H(curl))  -- div-free BY CONSTRUCTION
  (div o curl = 0).  The part of ker(curl) that is NOT a gradient is the de Rham
  cohomology H^1 (harmonic forms, one per topological hole) -- invisible on one
  element, assembled at mesh level as the tree-cotree / harmonic generators (next
  increment).  This file gives the per-element H(curl) basis and SHOWS both exact
  maps (grad H1 c H(curl) curl-free ; curl H(curl) c H(div) div-free).

  Ported from NGSolve fem/hcurlhofe_impl.hpp.  Two layers:

   (1) Schoeberl-Zaglmayr COVARIANT helpers (Du, uDv_minus_vDu, ...): NGSolve's
       building blocks; curl identities are the FEEC content, self-tested directly.
   (2) TENSOR de Rham construction on the cube/square (Hex, Quad): the clean
       explicit-polynomial Nedelec basis spanning NGSolve's order-p H(curl).

  Simplex lowest order = Whitney W1 (tet/trig edges; cf. derham.wls).  The curl-free
  "gradient block" grad(H1 bubbles) is given explicitly (the ker-curl part).
*)

(* ---- polynomial foundation ---- *)
l[i_, s_] := LegendreP[i, s];
L[i_, s_] := L[i, s] = Simplify[Integrate[l[i - 1, u], {u, -1, s}]];
oneDH1[q_, t_] := Join[{1 - t, t}, Table[L[i, 2 t - 1], {i, 2, q}]];   (* W: H1, deg<=q *)
oneDL2[q_, t_] := Table[LegendreP[j, 2 t - 1], {j, 0, q - 1}];         (* Q: L2, deg<=q-1 *)

(* ============ Schoeberl-Zaglmayr covariant helpers (NGSolve) ============ *)
gradd[f_] := D[f, {{x, y, z}}];
curlv[V_] := Curl[V, {x, y, z}];
divv[V_]  := D[V[[1]], x] + D[V[[2]], y] + D[V[[3]], z];

Du[u_]              := gradd[u];                       (* gradient field: curl = 0 *)
uDv[u_, v_]         := u gradd[v];
uDvMinusVDu[u_, v_] := u gradd[v] - v gradd[u];        (* Whitney/Nedelec: curl = 2 gu x gv *)
wuDvMinusWvDu[u_, v_, w_] := w (u gradd[v] - v gradd[u]);

(* ===================== HEX (tensor de Rham, Nedelec_p) ================= *)
(* x-family: Q(x) (x) W(y) (x) W(z) e_x ; +cyclic.  dim = 3 (p+1)(p+2)^2. *)
HCurlHex[p_] := Module[{Wf, Qf, xf, yf, zf},
  Wf[t_] := oneDH1[p + 1, t]; Qf[t_] := oneDL2[p + 1, t];
  xf = Flatten[Table[{Qf[x][[a]] Wf[y][[b]] Wf[z][[c]], 0, 0},
       {a, p + 1}, {b, p + 2}, {c, p + 2}], 2];
  yf = Flatten[Table[{0, Wf[x][[a]] Qf[y][[b]] Wf[z][[c]], 0},
       {a, p + 2}, {b, p + 1}, {c, p + 2}], 2];
  zf = Flatten[Table[{0, 0, Wf[x][[a]] Wf[y][[b]] Qf[z][[c]]},
       {a, p + 2}, {b, p + 2}, {c, p + 1}], 2];
  Join[xf, yf, zf]];

(* H1 hex (= grad source) and H(div) hex (= curl target) for the de Rham showcase *)
H1HexT[p_] := Module[{Wf}, Wf[t_] := oneDH1[p + 1, t];
  Flatten[Table[Wf[x][[a]] Wf[y][[b]] Wf[z][[c]], {a, p + 2}, {b, p + 2}, {c, p + 2}]]];
HDivHexT[p_] := Module[{Wf, Qf, xf, yf, zf},
  Wf[t_] := oneDH1[p + 1, t]; Qf[t_] := oneDL2[p + 1, t];
  xf = Flatten[Table[{Wf[x][[a]] Qf[y][[b]] Qf[z][[c]], 0, 0}, {a, p + 2}, {b, p + 1}, {c, p + 1}], 2];
  yf = Flatten[Table[{0, Qf[x][[a]] Wf[y][[b]] Qf[z][[c]], 0}, {a, p + 1}, {b, p + 2}, {c, p + 1}], 2];
  zf = Flatten[Table[{0, 0, Qf[x][[a]] Qf[y][[b]] Wf[z][[c]]}, {a, p + 1}, {b, p + 1}, {c, p + 2}], 2];
  Join[xf, yf, zf]];

(* ===================== QUAD (2D tensor de Rham) ======================= *)
curl2[V_] := D[V[[2]], x] - D[V[[1]], y];                 (* scalar 2D curl *)
HCurlQuad[p_] := Module[{Wf, Qf, xf, yf},
  Wf[t_] := oneDH1[p + 1, t]; Qf[t_] := oneDL2[p + 1, t];
  xf = Flatten[Table[{Qf[x][[a]] Wf[y][[b]], 0}, {a, p + 1}, {b, p + 2}], 1];
  yf = Flatten[Table[{0, Wf[x][[a]] Qf[y][[b]]}, {a, p + 2}, {b, p + 1}], 1];
  Join[xf, yf]];
H1QuadT[p_] := Module[{Wf}, Wf[t_] := oneDH1[p + 1, t];
  Flatten[Table[Wf[x][[a]] Wf[y][[b]], {a, p + 2}, {b, p + 2}]]];

(* ===================== SIMPLEX lowest order (Whitney W1) =============== *)
lamTet = {x, y, z, 1 - x - y - z};
glTet = gradd /@ lamTet;
W1tet[{a_, b_}] := lamTet[[a]] glTet[[b]] - lamTet[[b]] glTet[[a]];     (* = uDvMinusVDu *)
tetEdges = {{1, 2}, {1, 3}, {1, 4}, {2, 3}, {2, 4}, {3, 4}};
HCurlTetW1[] := W1tet /@ tetEdges;                                       (* 6 edges *)

lamTrig = {1 - x - y, x, y};
gl2 = D[#, {{x, y}}] & /@ lamTrig;
W1trig[{a_, b_}] := lamTrig[[a]] gl2[[b]] - lamTrig[[b]] gl2[[a]];
trigEdges = {{1, 2}, {1, 3}, {2, 3}};
HCurlTrigW1[] := W1trig /@ trigEdges;                                    (* 3 edges *)

(* curl-free "gradient block" = grad of the H1 interior/edge bubbles (ker curl) *)
H1TetBubbles[p_] := Module[{lam, edg, edges},                            (* edge bubbles, p>=2 *)
  lam = lamTet; edges = tetEdges;
  edg = Flatten@Table[
     With[{la = lam[[e[[1]]]], lb = lam[[e[[2]]]]}, L[i, lb - la]],       (* on-edge int-Legendre *)
     {e, edges}, {i, 2, p}];
  edg];
GradBlockTet[p_] := gradd /@ H1TetBubbles[p];

(* ============================ self-test ============================ *)
Module[{ok = True, chk, u, v, w, pts, pts2, rB, rA, p},
  chk[name_, cond_] := (Print["  ", If[cond, "ok  ", "FAIL"], name]; ok = ok && cond);
  Print["hcurl.wls self-test:"];

  (* ---- (1) SZ covariant-helper curl identities (FEEC content) ---- *)
  u = x^2 y + z; v = y z + x; w = x + 2 y + 3 z + x y;
  chk["Du is curl-free (gradient field)", Simplify[curlv[Du[u]]] === {0, 0, 0}];
  chk["uDvMinusVDu curl = 2 grad u x grad v (Nedelec)",
      Simplify[curlv[uDvMinusVDu[u, v]] - 2 Cross[gradd[u], gradd[v]]] === {0, 0, 0}];

  (* ---- (2) HEX tensor Nedelec_p: dimension, 12 edges at p=0 ---- *)
  chk["dim HCurlHex(p) = 3 (p+1)(p+2)^2, p=0..2",
      And @@ Table[Length[HCurlHex[p]] == 3 (p + 1) (p + 2)^2, {p, 0, 2}]];
  chk["HCurlHex(0) has 12 funcs (cube edges)", Length[HCurlHex[0]] == 12];

  SeedRandom[21];
  pts = Table[{x -> RandomReal[], y -> RandomReal[], z -> RandomReal[]}, {40}];

  (* ---- (3) de Rham LEFT: grad(H1) c H(curl) AND curl-free ---- *)
  Module[{grads},
   grads = gradd /@ H1HexT[1];
   chk["grad(H1Hex) is curl-free", And @@ (Simplify[curlv[#]] === {0, 0, 0} & /@ grads)];
   rB = MatrixRank[(Flatten[Table[# /. pt, {pt, pts}]] & /@ HCurlHex[1]), Tolerance -> 1*^-9];
   rA = MatrixRank[(Flatten[Table[# /. pt, {pt, pts}]] & /@ Join[HCurlHex[1], grads]), Tolerance -> 1*^-9];
   chk["grad(H1Hex) lies IN span(HCurlHex)  (de Rham: grad H1 c H(curl))", rA == rB]];

  (* ---- (4) de Rham RIGHT: curl(H(curl)) c H(div) AND div-free ---- *)
  Module[{curls},
   curls = curlv /@ HCurlHex[1];
   chk["curl(HCurlHex) is div-free (loops div-free by construction)",
       And @@ (Simplify[divv[#]] === 0 & /@ curls)];
   rB = MatrixRank[(Flatten[Table[# /. pt, {pt, pts}]] & /@ HDivHexT[1]), Tolerance -> 1*^-9];
   rA = MatrixRank[(Flatten[Table[# /. pt, {pt, pts}]] & /@ Join[HDivHexT[1], curls]), Tolerance -> 1*^-9];
   chk["curl(HCurlHex) lies IN span(HDivHex)  (de Rham: curl H(curl) c H(div))", rA == rB]];

  (* ---- (5) QUAD 2D Nedelec_p ---- *)
  chk["dim HCurlQuad(p) = 2 (p+1)(p+2), p=0..2",
      And @@ Table[Length[HCurlQuad[p]] == 2 (p + 1) (p + 2), {p, 0, 2}]];
  chk["HCurlQuad(0) has 4 funcs (square edges)", Length[HCurlQuad[0]] == 4];
  Module[{grads, ra, rb},
   grads = (D[#, {{x, y}}] & /@ H1QuadT[1]);
   chk["grad(H1Quad) is 2D-curl-free", And @@ (Simplify[curl2[#]] === 0 & /@ grads)];
   SeedRandom[22]; pts2 = Table[{x -> RandomReal[], y -> RandomReal[]}, {20}];
   rb = MatrixRank[(Flatten[Table[# /. pt, {pt, pts2}]] & /@ HCurlQuad[1]), Tolerance -> 1*^-9];
   ra = MatrixRank[(Flatten[Table[# /. pt, {pt, pts2}]] & /@ Join[HCurlQuad[1], grads]), Tolerance -> 1*^-9];
   chk["grad(H1Quad) c span(HCurlQuad)  (2D de Rham)", ra == rb]];

  (* ---- (6) simplex Whitney W1 + tangential dof ---- *)
  chk["HCurlTetW1 = 6 edge funcs, each curl = const (Whitney)",
      Length[HCurlTetW1[]] == 6 &&
      And @@ (FreeQ[Simplify[curlv[#]], x | y | z] & /@ HCurlTetW1[])];
  chk["HCurlTrigW1 = 3 edge funcs", Length[HCurlTrigW1[]] == 3];
  (* tangential (edge line-integral) Kronecker dof for the trig W1 *)
  chk["HCurlTrigW1 edge line-integral = Kronecker delta (Whitney dof)",
      With[{vc = {{0, 0}, {1, 0}, {0, 1}}},
       And @@ Flatten@Table[
         With[{c = ej[[1]], d = ej[[2]]},
          Module[{P = vc[[c]] + t (vc[[d]] - vc[[c]]), tang = vc[[d]] - vc[[c]], val},
           val = Integrate[(W1trig[ei] /. {x -> P[[1]], y -> P[[2]]}).tang, {t, 0, 1}];
           Simplify[val] == If[ei === ej, 1, 0]]],
         {ei, trigEdges}, {ej, trigEdges}]]];

  (* ---- (7) gradient block = ker curl (the cohomologically-trivial loops' complement) ---- *)
  chk["grad(H1Tet edge bubbles) is curl-free (ker curl, p=3)",
      And @@ (Simplify[curlv[#]] === {0, 0, 0} & /@ GradBlockTet[3])];
  chk["grad(H1Tet edge bubbles) are independent (rank = 6(p-1), p=3)",
      MatrixRank[(Flatten[Table[# /. pt, {pt, pts}]] & /@ GradBlockTet[3]),
                 Tolerance -> 1*^-9] == 6 (3 - 1)];

  Print[If[ok, "ALL PASS", "SOME FAILED"]];
];
