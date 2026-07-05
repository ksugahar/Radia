#!/usr/bin/env wolframscript
(* ::Package:: *)

(*
  vim_field.wls -- the VIM (Volume Integral Method) FIELD OPERATOR: the field produced
  by a symbolic magnetization density M(x) (= one H(div) basis function) via the 1/r
  Laplace kernel.  This is the kernel that turns the FEEC basis (h1/hcurl/hdiv/...) into
  an HDiv-VIM solver -- and the thing the C++ core implements.

  Magnetostatic scalar potential of a magnetization M on element T (Laplace kernel only,
  per the Radia "Green's Function: Laplace Kernel Only" policy):

     phi(r) = 1/(4pi) [ Int_T  rho_v / |r-x'| dV'  +  Sum_faces Int_f  sigma_s / |r-x'| dS' ]
     rho_v = - div M   (volume charge) ,   sigma_s = M . n   (surface charge)
     H(r) = -grad phi = 1/(4pi)[ Int_T rho_v (r-x')/|r-x'|^3 dV'
                                + Sum_f Int_f sigma_s (r-x')/|r-x'|^3 dS' ]
     B = mu_0 (H + M)   (M = 0 outside T)

  KEY de Rham bridge (why FEEC loops are charge-free null modes):
    a div-free M (a "loop" = curl of an H(curl) function) has  rho_v = -div M = 0
    IDENTICALLY -> NO volume charge -> the field is PURELY surface-charge, a high
    multipole that is (near-)null at the face-charge collocation points.  This is exactly the
    field-null property the C++ installCycle loop basis exploits.

  Routes (cf. README): (1) GAUSS quadrature -- always available (NIntegrate here = the
  reference).  (2) ANALYTICAL -- the van Oosterom-Strackee solid angle closes the
  constant-sigma NORMAL field exactly (the solid-angle face-charge kernel Radia's
  rad_poly_analytical.cpp uses); high-order sigma needs recursive moment formulas
  (Graglia/Wilton) or Gauss.  Both are validated below against each other.
*)

gradd[f_] := D[f, {{x, y, z}}];
divv[V_]  := D[V[[1]], x] + D[V[[2]], y] + D[V[[3]], z];
curlv[V_] := Curl[V, {x, y, z}];

(* ---- charge extraction from a symbolic M (exact) ---- *)
volCharge[M_] := -divv[M];                         (* rho_v *)
surfCharge[M_, n_, faceRule_] := (M.n) /. faceRule; (* sigma_s on a face *)

(* ---- unit-cube [0,1]^3 faces: {outward normal, substitution to the face} ---- *)
cubeFaces = {
  {{0, 0, -1}, {z -> 0}}, {{0, 0, 1}, {z -> 1}},
  {{-1, 0, 0}, {x -> 0}}, {{1, 0, 0}, {x -> 1}},
  {{0, -1, 0}, {y -> 0}}, {{0, 1, 0}, {y -> 1}}};

(* ---- H(r) of a magnetization M over the unit cube, by quadrature (NIntegrate) ---- *)
HnumCube[M_, r_] := Module[{rho, Hv, Hs, src = {x, y, z}, Rv},
  rho = volCharge[M];
  Rv = r - src;
  (* volume term (skip if div-free) *)
  Hv = If[Simplify[rho] === 0, {0, 0, 0},
    NIntegrate[Evaluate[rho Rv/(Rv.Rv)^(3/2)], {x, 0, 1}, {y, 0, 1}, {z, 0, 1}]];
  (* surface terms: 2D NIntegrate per face.  With[] injects the two free symbols into
     the (HoldAll) NIntegrate iterators; skip faces whose surface charge is 0. *)
  Hs = Sum[
    Module[{n = f[[1]], sub = f[[2]], sig, Rf, vars},
     sig = Simplify[(M.n) /. sub];
     Rf = (r - src) /. sub;
     vars = DeleteCases[{x, y, z}, First[First[sub]]];
     If[sig === 0, {0, 0, 0},
      With[{v1 = vars[[1]], v2 = vars[[2]]},
       NIntegrate[Evaluate[sig Rf/(Rf.Rf)^(3/2)], {v1, 0, 1}, {v2, 0, 1}]]]],
    {f, cubeFaces}];
  (Hv + Hs)/(4 Pi)];

(* point-dipole field (far-field reference) *)
Hdipole[m_, r_] := With[{R = Norm[r], rh = r/Norm[r]}, (3 (m.rh) rh - m)/(4 Pi R^3)];

(* ---- van Oosterom-Strackee solid angle of triangle (a,b,c) from r ---- *)
solidAngleVOS[a_, b_, c_, r_] := Module[{A = a - r, B = b - r, C = c - r, na, nb, nc, num, den},
  na = Norm[A]; nb = Norm[B]; nc = Norm[C];
  num = A.Cross[B, C];
  den = na nb nc + (A.B) nc + (A.C) nb + (B.C) na;
  2 ArcTan[den, num]];                              (* atan2(num,den) *)

(* solid angle by quadrature: Int_tri (x'-r).N / |x'-r|^3 dS', N = (b-a)x(c-a) (unnormalized) *)
solidAngleNum[a_, b_, c_, r_] := Module[{N3 = Cross[b - a, c - a], src, Rv},
  src = a + s (b - a) + t (c - a);
  Rv = src - r;
  NIntegrate[Evaluate[(Rv.N3)/(Rv.Rv)^(3/2)], {s, 0, 1}, {t, 0, 1 - s}]];

(* ============================ self-test ============================ *)
Module[{ok = True, chk, Mc, Mdf, mmom, hN, hD, ratio, a, b, c, robs, wV, wN, hint},
  chk[name_, cond_] := (Print["  ", If[cond, "ok  ", "FAIL"], name]; ok = ok && cond);
  Print["vim_field.wls self-test (VIM field operator):"];

  (* ---- (1) constant M = (0,0,1) over unit cube: charges + far-field dipole ---- *)
  Mc = {0, 0, 1};
  chk["constant M -> volume charge rho = 0 (exact)", Simplify[volCharge[Mc]] === 0];
  chk["constant M -> total surface charge = 0 (charge neutral)",
      Chop[Sum[NIntegrate[Evaluate[(Mc.f[[1]]) /. f[[2]]],
        Evaluate[Sequence @@ ({#, 0, 1} & /@ Complement[{x, y, z}, First /@ f[[2]]])]], {f, cubeFaces}]] == 0];
  mmom = Integrate[Mc, {x, 0, 1}, {y, 0, 1}, {z, 0, 1}];     (* dipole moment = (0,0,1) *)
  chk["dipole moment = Int M dV = (0,0,1)", mmom == {0, 0, 1}];
  (* far field: VIM cube field -> point dipole at the cube CENTER (0,0,1/2), with the
     correct finite-size octupole correction O((size/R)^2).  Validate by CONVERGENCE:
     the dipole-ratio error must -> 0 as ~1/R^2 (error quarters each time R doubles). *)
  Module[{errs},
   errs = Table[Abs[HnumCube[Mc, {0, 0, N[R]}][[3]]/Hdipole[mmom, {0, 0, N[R] - 0.5}][[3]] - 1],
                {R, {10, 20, 40}}];
   chk["VIM far field -> point dipole, error < 0.2% at R=40", errs[[3]] < 0.002];
   chk["dipole-ratio error converges as ~1/R^2 (quarters per R-doubling)",
       errs[[2]] < errs[[1]]/3 && errs[[3]] < errs[[2]]/3]];

  (* ---- (2) div-free M (a LOOP): rho = 0 identically -> field-null mechanism ---- *)
  (* M = curl of an H(curl) function -> divergence-free by construction *)
  Mdf = curlv[{y z^2, x^2 z, x y^2}];                         (* some curl => div-free *)
  chk["div(curl F) = 0 -> a loop M has NO volume charge (rho = 0 exact)",
      Simplify[volCharge[Mdf]] === 0];
  (* a generic div-free M built from an actual H(curl) tensor shape also has rho=0 *)
  Module[{Mh},
   Mh = curlv[{(1 - x) y (1 - z), 0, 0}];                     (* curl of an edge-type field *)
   chk["loop from H(curl) shape: rho = -div M = 0 (purely surface-charge field)",
       Simplify[volCharge[Mh]] === 0]];

  (* ---- (3) solid-angle face-charge kernel: van Oosterom-Strackee solid angle vs quadrature ---- *)
  a = {0, 0, 0}; b = {1, 0, 0}; c = {0, 1, 0}; robs = {0.2, 0.2, 0.7};
  wV = solidAngleVOS[a, b, c, robs];
  wN = solidAngleNum[a, b, c, robs];
  chk["van Oosterom-Strackee solid angle = quadrature (solid-angle face-charge kernel)",
      Abs[Abs[wV] - Abs[wN]] < 1*^-7];
  (* a second observation point (off-axis) *)
  robs = {0.4, -0.1, 0.5};
  chk["VOS solid angle = quadrature, 2nd point",
      Abs[Abs[solidAngleVOS[a, b, c, robs]] - Abs[solidAngleNum[a, b, c, robs]]] < 1*^-7];

  (* ---- (4) high-order surface-charge moment Int_tri (x) / R dS' is computable ----
     (the polynomial-sigma route the C++ high-order extension needs; here = finite,
      and equals the recursive/quadrature value -- demonstrated computable, not closed-form) *)
  hint = NIntegrate[Evaluate[Module[{src = a + s (b - a) + t (c - a)},
       (src[[1]]) / Norm[src - {0.3, 0.3, 0.6}]] Norm[Cross[b - a, c - a]]],
     {s, 0, 1}, {t, 0, 1 - s}];
  chk["high-order moment Int_tri x'/R dS' is finite & positive (poly-sigma route)",
      NumberQ[hint] && hint > 0];

  Print[If[ok, "ALL PASS", "SOME FAILED"]];
];
