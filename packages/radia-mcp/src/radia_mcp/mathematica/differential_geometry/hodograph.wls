#!/usr/bin/env wolframscript
(* ::Package:: *)

(*
  hodograph.wls -- the coordinate-transform BACKBONE of the Clebsch/hodograph
  design line, made symbolic.  Companion to weakform_hodge.wls (the Hodge /
  material-modulation half) and docs/clebsch_hodograph/HODOGRAPH_BACKBONE.md
  (the spine: "five transforms are five cells of ONE diagram, indexed by three
  axes").  Lab policy: understand the machinery through DIFFERENTIAL GEOMETRY.

  THE ONE DIAGRAM.  Magnetostatics is a closed 2-form with a (possibly nonlinear)
  Hodge-star constitutive law:
        B in Omega^2 ,  dB = 0   (de Rham: B = dA, no monopoles)   [TOPOLOGY half]
        H in Omega^1 ,  H = Star_nu B   (the material Hodge star)   [METRIC half]
  Under ANY diffeomorphism phi the two halves behave completely differently --
  THAT split is the whole story:
        dB = 0        pulls back NATURALLY (no metric data)  -> preserved by any phi
        H = Star_nu B the Hodge star pulls back with a WEIGHT -> the only thing that moves
  (the metric weight W = |det J|(J^T J)^{-1} is in weakform_hodge.wls).

  THE THREE AXES (each cell of the diagram is a choice on three independent axes):
    axis 1  which map phi          geometric Kelvin / potential-plane / field-plane / partial
    axis 2  dimension              2-D conformal weight-free; 3-D Clebsch needs helicity = 0
    axis 3  material               linear Star vs nonlinear Star_nu (saturation: stays SPD)

  Verified here, symbolically:
    [1] the split: the pullback is a cochain map (grad of pullback = chain rule) so
        dB=0 is metric-free; only Star carries the metric.
    [2] axis 1: the 2-D geometric Kelvin inversion x -> R^2 x/|x|^2 is CONFORMAL
        (J^T J = scalar*I) hence WEIGHT-FREE (W = I); the complex potential (Phi, A_z)
        is the Cauchy-Riemann conjugate (orthogonal equal-speed net).
    [3] axis 2: a Clebsch field B = grad(alpha) x grad(beta) is div-free (closed 2-form),
        A = alpha grad(beta) is its vector potential, and the helicity density A.B = 0
        pointwise -- WHILE the ABC Beltrami field (curl B = B) has A.B = |B|^2 != 0,
        the helicity obstruction (Moffatt 1969) to a global Clebsch pair in 3-D.
    [4] axis 3: in 2-D, A_z IS the Clebsch potential (B = grad(A_z) x zhat); flux lines
        are A_z level sets (B.grad(A_z) = 0), so A_z is the flux-line Hamiltonian
        (conserved), which is why the 2-D saturable problem linearises in the hodograph
        (Chaplygin) plane.  The nonlinear Star_nu stays SPD (weakform_hodge.wls) so the
        magnetic hodograph -- unlike transonic gas dynamics -- never folds (no limiting line).

  Prior art (classical; this is a working map, not an originality claim): complex-potential
  pole design (1900s); 2-D conformal weight-freeness (classical); geometric Kelvin FE open
  boundary (Wong-Ciric 1985, Freeman-Lowther 1988, Sugahara 2022; pullback form
  Nagamine-Yamaguchi-Sugahara CEFC 2026); Chaplygin/von Mises hodograph (gas dynamics);
  Clebsch + helicity (Robert 1991, Moffatt 1969).
*)

grad[f_] := D[f, {{x, y, z}}];
crl[F_]  := Curl[F, {x, y, z}];
dvg[F_]  := Div[F, {x, y, z}];
weightW[J_] := Abs[Det[J]] Inverse[Transpose[J].J];      (* = weakform_hodge.wls *)

(* ============================ self-test ============================ *)
Module[{ok = True, chk, zero, F2, J2, JtJ, alpha, beta, Bc, Ac, Babc, az, Bz2, pt},
  chk[name_, cond_] := (Print["  ", If[cond, "ok  ", "FAIL"], name]; ok = ok && cond);
  zero[m_] := Simplify[m] === ConstantArray[0, Dimensions[m]];
  Print["hodograph.wls self-test (the coordinate-transform backbone: 3 axes, symbolic):"];

  (* ---- [1] THE SPLIT: dB=0 is metric-free (pullback = cochain map, chain rule) ----
     For a 0-form f and a map phi, d(f o phi) = (J^T) (df o phi): grad of the pullback
     is the chain rule.  So d commutes with pullback -> exact/closed forms stay exact/closed
     under ANY phi, with NO metric data.  Verify on f = x^2 y + z, phi a generic linear map. *)
  chk["pullback is a cochain map: grad(f o phi) = J^T (grad f o phi)  (dB=0 metric-free)",
      Module[{phi, Jp, f, fpull, lhs, rhs},
       phi = {a11 x + a12 y + a13 z, a21 x + a22 y + a23 z, a31 x + a32 y + a33 z};
       Jp = D[phi, {{x, y, z}}];                       (* d(phi)/d(x,y,z) *)
       f = u^2 v + w;                                   (* a 0-form in target coords u,v,w *)
       fpull = f /. {u -> phi[[1]], v -> phi[[2]], w -> phi[[3]]};
       lhs = grad[fpull];
       rhs = Transpose[Jp].((D[f, {{u, v, w}}]) /. {u -> phi[[1]], v -> phi[[2]], w -> phi[[3]]});
       Simplify[lhs - rhs] === {0, 0, 0}]];

  (* ---- [2a] AXIS 1: 2-D geometric Kelvin inversion is CONFORMAL -> WEIGHT-FREE ----
     F(x,y) = R^2 (x,y)/(x^2+y^2).  J^T J = scalar * I  (conformal); hence W = I. *)
  F2 = R^2 {x, y}/(x^2 + y^2);
  J2 = D[F2, {{x, y}}];
  JtJ = Simplify[Transpose[J2].J2];
  chk["2-D Kelvin inversion is conformal: J^T J off-diagonal = 0",
      Simplify[JtJ[[1, 2]]] === 0 && Simplify[JtJ[[2, 1]]] === 0];
  chk["2-D Kelvin inversion is conformal: J^T J diagonal entries equal (= scalar*I)",
      Simplify[JtJ[[1, 1]] - JtJ[[2, 2]]] === 0];
  chk["2-D Kelvin inversion is WEIGHT-FREE: W = |det J|(J^T J)^{-1} = I  (mu' = mu0 in 2-D)",
      Simplify[weightW[J2] - IdentityMatrix[2], R > 0 && x^2 + y^2 > 0] === ConstantArray[0, {2, 2}]];

  (* ---- [2b] AXIS 1: the complex potential (Phi, A_z) is the Cauchy-Riemann conjugate net ----
     A holomorphic W = Phi + i Az has CR: Phi_x = Az_y, Phi_y = -Az_x.  Take W = z_c^3,
     z_c = x + i y: Phi = Re, Az = Im.  Then grad Phi . grad Az = 0 (orthogonal) and
     |grad Phi| = |grad Az| (equal speed) -- the hodograph net. *)
  chk["complex potential: grad(Phi) . grad(Az) = 0 (orthogonal hodograph net)",
      Module[{phiR = Re[(x + I y)^3] // ComplexExpand, azI = Im[(x + I y)^3] // ComplexExpand},
       Simplify[D[phiR, x] D[azI, x] + D[phiR, y] D[azI, y]] === 0]];
  chk["complex potential: |grad Phi| = |grad Az| (equal-speed conjugate net)",
      Module[{phiR = Re[(x + I y)^3] // ComplexExpand, azI = Im[(x + I y)^3] // ComplexExpand},
       Simplify[(D[phiR, x]^2 + D[phiR, y]^2) - (D[azI, x]^2 + D[azI, y]^2)] === 0]];

  (* ---- [3] AXIS 2: Clebsch field + the helicity obstruction (Moffatt) ----
     B = grad(alpha) x grad(beta) is an EXACT (hence closed) 2-form: div B = 0.
     A = alpha grad(beta) is a vector potential: curl A = B.  Helicity density A.B = 0
     pointwise (Clebsch is helicity-free).  ABC Beltrami: curl B = B so A=B, A.B=|B|^2 != 0. *)
  alpha = x^2 y + Sin[z]; beta = y z + x^2;
  Bc = Cross[grad[alpha], grad[beta]];
  chk["Clebsch B = grad(alpha) x grad(beta) is div-free (closed 2-form)",
      Simplify[dvg[Bc]] === 0];
  Ac = alpha grad[beta];
  chk["Clebsch A = alpha grad(beta) is a vector potential: curl A = B",
      zero[Simplify[crl[Ac] - Bc]]];
  chk["Clebsch helicity density A.B = 0 pointwise (helicity-free -> global Clebsch exists)",
      Simplify[Ac.Bc] === 0];
  (* ABC Beltrami field: curl B = B (Beltrami), helicity |B|^2 != 0 *)
  Babc = {aA Sin[z] + cC Cos[y], bB Sin[x] + aA Cos[z], cC Sin[y] + bB Cos[x]};
  chk["ABC field is Beltrami: curl B = B",
      zero[Simplify[crl[Babc] - Babc]]];
  pt = {x -> 3/10, y -> 1/2, z -> 7/10, aA -> 1, bB -> 1, cC -> 1};
  chk["ABC helicity density A.B = |B|^2 != 0 (helicity obstruction, Moffatt 1969)",
      (Simplify[Babc.Babc] /. pt) > 0];

  (* ---- [4] AXIS 3: in 2-D, A_z IS the Clebsch potential; flux lines are A_z level sets ----
     B = grad(A_z) x zhat = (A_z,y, -A_z,x, 0): div-free automatically, and B.grad(A_z) = 0
     so flux lines are A_z = const -- A_z is the flux-line Hamiltonian (conserved).  This is
     why the saturable 2-D problem linearises in the hodograph (Chaplygin) plane. *)
  az = azf[x, y];
  Bz2 = Cross[grad[az], {0, 0, 1}];
  chk["2-D B = grad(A_z) x zhat is div-free (A_z = Clebsch potential alpha, beta = z)",
      Simplify[dvg[Bz2]] === 0];
  chk["2-D flux lines are A_z level sets: B . grad(A_z) = 0 (A_z = flux-line Hamiltonian)",
      Simplify[Bz2.grad[az]] === 0];

  Print[If[ok, "ALL PASS", "SOME FAILED"]];
];
