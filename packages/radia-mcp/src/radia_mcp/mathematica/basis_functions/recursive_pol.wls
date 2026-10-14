#!/usr/bin/env wolframscript
(* ::Package:: *)

(*
  recursive_pol.wls -- the polynomial foundation of the NGSolve high-order
  hierarchical FEM shape functions (Zaglmayr construction).  Clean, runnable
  port of the recursions used by the H1 / H(curl) / H(div) shape functions;
  NGSolve C++ source: fem/recursive_pol.hpp.  See README.md.

  Symbols:
    l[i, x]       Legendre P_i(x)
    L[i, x]       integrated Legendre  L_i(x) = Integral_{-1}^x l_{i-1}(z) dz
                  (vanishes at x = +-1 for i >= 2; the edge/face bubble seed)
    lS[n, s, t]   scaled Legendre        t^n l[n, s/t]   (homogeneous, for
    LS[n, s, t]   scaled integrated Leg. t^n L[n, s/t]    collapsed coords)
*)

l[i_, x_] := LegendreP[i, x];

(* integrated Legendre: define symbolically AND memoise per order *)
L[i_, x_] := L[i, x] = Simplify[Integrate[l[i - 1, z], {z, -1, x}]];

(* scaled (homogeneous) forms -- used on triangles/tets where a barycentric
   pair (s, t) collapses; lS/LS are degree-n homogeneous in (s, t). *)
lS[n_, s_, t_] := lS[n, s, t] = Simplify[t^n  l[n, s/t]];
LS[n_, s_, t_] := LS[n, s, t] = Simplify[t^n  L[n, s/t]];

(* ---- self-test (run: wolframscript -file recursive_pol.wls) ---- *)
If[$ScriptCommandLine =!= {} || True,
  Module[{ok = True, chk},
   chk[name_, cond_] := (Print["  ", If[cond, "ok  ", "FAIL"], name];
                         ok = ok && cond);
   Print["recursive_pol.wls self-test:"];
   (* L_2(x) = (x^2 - 1)/2 *)
   chk["L[2,x] == (x^2-1)/2", Simplify[L[2, x] - (x^2 - 1)/2] === 0];
   (* integrated Legendre vanishes at the endpoints for i >= 2 *)
   chk["L[i,+-1] == 0, i=2..5", And @@ Table[
        Simplify[(L[i, 1]) ] === 0 && Simplify[(L[i, -1])] === 0, {i, 2, 5}]];
   (* scaled Legendre is homogeneous degree n: lS[n,a s,a t] == a^n lS[n,s,t] *)
   chk["lS homogeneous deg n", And @@ Table[
        Simplify[lS[n, a s, a t] - a^n lS[n, s, t]] === 0, {n, 0, 3}]];
   (* L recursion ties back to Legendre: d/dx L[i,x] == l[i-1,x] *)
   chk["L' == l[i-1]", And @@ Table[
        Simplify[D[L[i, x], x] - l[i - 1, x]] === 0, {i, 2, 5}]];
   Print[If[ok, "ALL PASS", "SOME FAILED"]];
  ];
];
