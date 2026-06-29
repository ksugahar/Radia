#!/usr/bin/env wolframscript
(* ::Package:: *)

(*
  weakform_hodge.wls -- differential geometry made OPERATIONAL for FEM:
  the WEAK (energy) FORM as a Hodge pairing, and the HODGE STAR as the carrier
  of the MATERIAL CONSTANT -- i.e. how the material law is MODULATED under a
  coordinate map (mesh deformation / Kelvin open boundary / a hodograph chart).

  Lab policy (the reason this file exists): use differential geometry to make
  two things easy to SEE --

    (1) the magnetostatic weak form
            a(w, A) = INT_Omega  <dw, nu dA>  dOmega
        is a HODGE PAIRING.  The exterior derivative d carries only TOPOLOGY
        (it is metric-free); the inner product <.,.> -- equivalently the Hodge
        star Star -- carries the MATERIAL nu and the METRIC.  Constitutive law:
            H = nu Star(B)      (B a 2-form, H a 1-form; Star carries nu).

    (2) under a coordinate map phi (deformation / Kelvin / hodograph) ONLY the
        Hodge star moves, and it moves as an EXPLICIT WEIGHT.  The transformed
        material constant is
            nu'  =  nu |det P| (P^T P)^{-1}                          [eqA]
        where P is the pullback on the ORTHONORMAL 1-form basis.  This is the
        SAME algebraic form as the scalar Dirichlet (Laplacian) weight
            W    =     |det J| (J^T J)^{-1} ,                        [eqB]
        J = d(physical)/d(computational).  "Material modulation = the metric
        weight" is the whole point: [eqA] and [eqB] are one formula.

  [eqA] reproduces, symbolically, the Nagamine-Yamaguchi-Sugahara PULLBACK KELVIN
  law (CEFC 2026, "A Pullback-Based Formulation of Kelvin Transformation"):
    spherical  inversion (conformal)     -> nu' = (r'/R)^2 nu        (isotropic)
    cylindrical inversion (non-conformal) -> nu' = diag(1,1,(rho'/R)^4) nu (anisotropic)
  and matches W of  docs/clebsch_hodograph/demos/bidirectional_coordinate_transform_2d.py
  (symbolic_pullback_check, golden) + docs/clebsch_hodograph/HODOGRAPH_BACKBONE.md sec.1.

  d carries topology, Star carries material+metric: the pullback (kStar) is a
  COCHAIN MAP (kStar.d = d.kStar) so dB=0 survives ANY phi for free; all that
  MOVES is in Star.

  Sec.7: the metric IS the geometry, but it is NOT the curvature.  The action is a
  coordinate-invariant scalar; the metric (Hodge) carries the map AND the material;
  yet a coordinate map / Kelvin / material weight changes the metric while the space
  stays FLAT (Riemann = 0) -- a diffeomorphism cannot curve flat space.  ("Material =
  effective metric" is transformation optics; genuine curvature needs curved spacetime.)
*)

(* ================= the one formula: Hodge weight on 1-forms ================= *)
(* scalar Dirichlet (0-form energy) weight, J = coordinate Jacobian:
   INT_phys grad_x u . grad_x v dx = INT_comp (grad_xi u)^T W (grad_xi v) dxi *)
weightW[J_] := Abs[Det[J]] Inverse[Transpose[J].J];

(* material modulation factor, P = pullback on the ORTHONORMAL 1-form basis:
   nu' = nu * nuFactor[P].  SAME formula as weightW -- that is the unification. *)
nuFactor[P_] := Abs[Det[P]] Inverse[Transpose[P].P];

(* nonlinear Hodge *_nu tangent  dH/dB = nu I + nu'(|B|) (B (x) B)/|B|  (flat helper) *)
tangentNu[nu_, nup_, bx_, by_, bz_] :=
  nu IdentityMatrix[3] + nup Outer[Times, {bx, by, bz}, {bx, by, bz}]/Sqrt[bx^2 + by^2 + bz^2];

(* ---- Riemann curvature of a metric g(xs): metric is local geometry, curvature is the
   coordinate-INVARIANT 2nd-order quantity.  A coordinate map / Kelvin / material weight
   changes g but NOT the curvature (a diffeomorphism of flat space stays flat). ---- *)
christo[g_, xs_] := With[{gi = Inverse[g], n = Length[xs]},
  Table[(1/2) Sum[gi[[k, l]] (D[g[[l, j]], xs[[i]]] + D[g[[l, i]], xs[[j]]]
                              - D[g[[i, j]], xs[[l]]]), {l, n}], {k, n}, {i, n}, {j, n}]];
riem[g_, xs_] := With[{G = christo[g, xs], n = Length[xs]},   (* R^a_{bcd} *)
  Table[D[G[[a, b, d]], xs[[c]]] - D[G[[a, b, c]], xs[[d]]]
        + Sum[G[[a, c, e]] G[[e, b, d]] - G[[a, d, e]] G[[e, b, c]], {e, n}],
    {a, n}, {b, n}, {c, n}, {d, n}]];
scalarR[g_, xs_] := With[{gi = Inverse[g], Rm = riem[g, xs], n = Length[xs]},
  Sum[gi[[b, d]] Rm[[a, b, a, d]], {a, n}, {b, n}, {d, n}]];   (* g^{bd} R^a_{bad} *)

(* premetric "metric -> material" map (Plebanski / transformation optics): the material
   tensor that makes a flat-space field SEE the metric g is  chi(g) = sqrt(det g) g^{-1}.
   This is the Hodge star written as a material; ALL of {coordinate map, material, IE
   decay, PML stretch} is this one map -- "push everything onto the Hodge". *)
chiMat[g_] := Sqrt[Det[g]] Inverse[g];

(* ============================ self-test ============================ *)
Module[{ok = True, chk, zero, diagQ},
  chk[name_, cond_] := (Print["  ", If[cond, "ok  ", "FAIL"], name]; ok = ok && cond);
  zero[m_] := Simplify[m] === ConstantArray[0, Dimensions[m]];
  Print["weakform_hodge.wls self-test (weak form = Hodge pairing; Hodge star = material modulation):"];

  (* ---- (0) reference frame: no map => no modulation ---- *)
  chk["reference frame: scalar weight W(I) = I",
      zero[weightW[IdentityMatrix[3]] - IdentityMatrix[3]]];
  chk["reference frame: material factor nuFactor(I) = I  (nu' = nu)",
      zero[nuFactor[IdentityMatrix[3]] - IdentityMatrix[3]]];

  (* ---- (1) the unification [eqA]=[eqB]: material weight = (2-form pullback cofactor) energy ----
     A 2-form pulls back by the cofactor cof(P) = det(P) P^{-T} (1-forms by P, volume by det P).
     Energy invariance nu (cof(P) b).(cof(P) b)/|det P| = nu' b.b  =>  nu' = nu * cof^T cof/|det P|.
     Claim: that equals nu * |det P| (P^T P)^{-1} = nu * nuFactor[P].  (det P > 0 branch.) *)
  chk["material weight = cofactor energy:  cof(P)^T cof(P)/det P = nuFactor[P]  (the unification)",
      With[{P = DiagonalMatrix[{pa, pb, pc}]},
        With[{cof = Det[P] Inverse[Transpose[P]]},
          Simplify[Transpose[cof].cof/Det[P] - nuFactor[P], pa > 0 && pb > 0 && pc > 0]
            === ConstantArray[0, {3, 3}]]]];

  (* ---- (2) 2-D CONFORMAL map is WEIGHT-FREE  (Hodge * conformally invariant on 1-forms in 2-D) ---- *)
  chk["2-D conformal J = s R  =>  W = I  (weight-free; mu' = mu0, no Kelvin factor in 2-D)",
      With[{Jc = sC {{Cos[thC], -Sin[thC]}, {Sin[thC], Cos[thC]}}},
        Simplify[weightW[Jc] - IdentityMatrix[2], sC > 0] === ConstantArray[0, {2, 2}]]];

  (* a GENERIC 2-D map keeps a non-trivial metric weight (not the identity) *)
  chk["2-D generic J = {{a,b},{c,d}}  =>  W != I  (a real metric survives)",
      With[{Jg = {{ga, gb}, {gc, gd}}},
        ! (Simplify[weightW[Jg] - IdentityMatrix[2]] === ConstantArray[0, {2, 2}])]];

  (* ---- (3) 3-D conformal scaling does NOT cancel: W = s I  (weight survives in 3-D) ---- *)
  chk["3-D isotropic scale J = s I  =>  W = s I  (NOT weight-free: this is why 3-D Kelvin carries (r'/R)^2)",
      With[{J3 = s3 IdentityMatrix[3]},
        Simplify[weightW[J3] - s3 IdentityMatrix[3], s3 > 0] === ConstantArray[0, {3, 3}]]];

  (* ---- (4) PULLBACK KELVIN, spherical (conformal) -> isotropic nu' = (r'/R)^2 nu ----
     orthonormal 1-form basis (e_r,e_th,e_ph)=(dr, r dth, r sin th dph); inversion r->R^2/r.
     pullback on the basis:  k* e_r' = -(R/r)^2 e_r ,  k* e_th' = +(R/r)^2 e_th ,  k* e_ph' = +(R/r)^2 e_ph.
     Let sig = (R/r)^2 = (r'/R)^2 .  P = diag(-sig, sig, sig). *)
  chk["spherical Kelvin: P = diag(-sig,sig,sig)  =>  nu' = sig * I = (r'/R)^2 nu  (ISOTROPIC, eq.19)",
      With[{P = DiagonalMatrix[{-sig, sig, sig}]},
        Simplify[nuFactor[P] - sig IdentityMatrix[3], sig > 0] === ConstantArray[0, {3, 3}]]];

  (* ---- (5) PULLBACK KELVIN, cylindrical (NON-conformal) -> anisotropic nu' ----
     (e_rho,e_phi,e_z)=(drho, rho dphi, dz); inversion rho->R^2/rho, z untouched.
     k* e_rho' = -(R/rho)^2 e_rho , k* e_phi' = +(R/rho)^2 e_phi , k* e_z' = e_z.
     sig = (R/rho)^2 = (rho'/R)^2 .  P = diag(-sig, sig, 1)  =>  nu' = diag(1,1,sig^2) = diag(1,1,(rho'/R)^4). *)
  diagQ[m_, d_] := Simplify[m - DiagonalMatrix[d], sig > 0] === ConstantArray[0, {Length[d], Length[d]}];
  chk["cylindrical Kelvin: P = diag(-sig,sig,1)  =>  nu' = diag(1,1,sig^2) = diag(1,1,(rho'/R)^4) nu  (ANISOTROPIC, eq.26)",
      With[{P = DiagonalMatrix[{-sig, sig, 1}]}, diagQ[nuFactor[P], {1, 1, sig^2}]]];

  (* anisotropy IS the signature of non-conformality: the cylindrical nu' is NOT a scalar multiple of I *)
  chk["cylindrical nu' is anisotropic (not a scalar*I) -- the non-conformal signature",
      With[{P = DiagonalMatrix[{-sig, sig, 1}]},
        ! (Simplify[nuFactor[P] - (nuFactor[P][[1, 1]]) IdentityMatrix[3], sig > 0]
             === ConstantArray[0, {3, 3}])]];

  (* ---- (6) NONLINEAR Hodge  *_nu : H = nu(|B|) B.  Tangent dH/dB = nu I + nu'(|B|)(B (x) B)/|B|. ----
     Eigen-structure: longitudinal eigenvalue nu + nu'|B| (along B), transverse nu (perp to B).
     SPD (elliptic, no fold / no limiting line) iff nu >= 0 AND nu + nu'|B| >= 0  (monotone reluctivity).
     This is the differential-geometric reason magnetic saturation does NOT fold the hodograph
     (unlike transonic gas dynamics, where the analogous operator changes type at the sonic line). *)
  (* longitudinal: tangent . B = (nu + nu'|B|) B *)
  chk["nonlinear *_nu tangent: longitudinal eigenvalue = nu + nu'|B|  (tangent.B = (nu+nu'|B|) B)",
      Simplify[tangentNu[nu, nup, bx, by, bz].{bx, by, bz}
               - (nu + nup Sqrt[bx^2 + by^2 + bz^2]) {bx, by, bz}] === {0, 0, 0}];
  (* transverse: any t with t.B = 0 has tangent . t = nu t  (use t = {-by, bx, 0}) *)
  chk["nonlinear *_nu tangent: transverse eigenvalue = nu  (tangent.t = nu t, t perp B)",
      Simplify[tangentNu[nu, nup, bx, by, bz].{-by, bx, 0} - nu {-by, bx, 0}] === {0, 0, 0}];
  (* SPD <=> elliptic: both eigenvalues > 0 for monotone reluctivity (nu>0, nu'>=0).
     Implies[premise, conclusion]: monotone reluctivity (nu>0, nu'>=0, |B|>0) forces
     BOTH the transverse (nu) and longitudinal (nu + nu'|B|) eigenvalues positive. *)
  chk["nonlinear *_nu stays SPD (elliptic, no limiting line) for nu>0, nu'>=0",
      Resolve[ForAll[{nu, nup, nbmag},
                Implies[nu > 0 && nup >= 0 && nbmag > 0,
                        nu > 0 && nu + nup nbmag > 0]], Reals] === True];

  (* ---- (7) METRIC vs CURVATURE: the metric IS the geometry, but it is NOT the curvature ----
     The action INT (1/2)<dA, nu dA> is a coordinate-invariant scalar (same form in any chart);
     ALL the geometry -- the coordinate map AND the material -- lives in the metric (Hodge star).
     BUT the metric is the LOCAL geometry (lengths/angles); CURVATURE is the coordinate-INVARIANT
     2nd-order quantity derived from it.  A coordinate map / Kelvin inversion / material weight
     changes the metric, yet a diffeomorphism of flat space CANNOT create curvature: the space
     stays FLAT (Riemann = 0), the metric just looks non-trivial in the chosen chart.  (The
     deep "material = effective metric" identity is transformation optics: the field sees an
     effective geometry, the space stays flat.  Genuine curvature needs curved spacetime / GR.) *)

  (* the 2-D Kelvin inversion metric g = (R^2/|x|^2)^2 I is non-trivial but FLAT (Riemann = 0) *)
  chk["(7a) the 2-D Kelvin inversion metric (R^2/|x|^2)^2 I is FLAT (Riemann = 0): a map can't curve space",
      Simplify[Flatten[riem[(R^2/(x^2 + y^2))^2 IdentityMatrix[2], {x, y}]],
               x^2 + y^2 > 0 && R > 0] === ConstantArray[0, 2^4]];

  (* the polar metric diag(1, r^2) is non-trivial but FLAT -- non-trivial metric =/= curved *)
  chk["(7b) the polar metric diag(1, r^2) is FLAT (Riemann = 0): a non-trivial metric is NOT curvature",
      Simplify[Flatten[riem[DiagonalMatrix[{1, rr^2}], {rr, ph}]], rr > 0]
        === ConstantArray[0, 2^4]];

  (* CONTROL: the round 2-sphere genuinely IS curved -- the machinery detects it (R_scalar = 2) *)
  chk["(7c) CONTROL: the round 2-sphere diag(1, sin^2 th) is CURVED, scalar curvature = 2 (machinery detects it)",
      Simplify[scalarR[DiagonalMatrix[{1, Sin[th]^2}], {th, ph}], 0 < th < Pi] === 2];

  (* ---- (8) TRANSFORMATION OPTICS = "push everything onto the Hodge" ----
     premetric EM (Hehl-Obukhov): dF=0, dG=J are metric-free; ONLY the constitutive Hodge
     G = star F carries metric+material.  A coordinate map can therefore be traded entirely
     for a material = the metric->material map chi(g) = sqrt(det g) g^{-1}.  Our Kelvin rule
     nu' and Pendry's transformation-optics material are BOTH this one chi, on the pullback
     vs the pushforward metric -- "push everything onto the Hodge." *)
  chk["(8a) our material rule IS the premetric map:  nuFactor[P] = chi(P^T P) = sqrt(det g) g^{-1}  (= transformation optics)",
      With[{P = DiagonalMatrix[{pa, pb, pc}]},
        Simplify[nuFactor[P] - chiMat[Transpose[P].P], pa > 0 && pb > 0 && pc > 0]
          === ConstantArray[0, {3, 3}]]];
  chk["(8b) Pendry vacuum TO material  eps' = Lam Lam^T / det Lam  is the SAME chi on the pushforward metric (Lam Lam^T)^{-1}",
      With[{Lam = {{la, lb}, {lc, ld}}},
        Simplify[chiMat[Inverse[Lam.Transpose[Lam]]] - Lam.Transpose[Lam]/Det[Lam],
                 la ld - lb lc > 0] === ConstantArray[0, {2, 2}]]];

  (* ---- (9) the OTHER half of "push onto the Hodge": a coordinate map cannot CURVE space.
     The pullback metric g = J^T J of ANY diffeomorphism is FLAT (Riemann = 0) -- the map
     moves the metric (the geometry) but not the curvature.  Verified for a conformal AND a
     non-conformal map (the two bidirectional_coordinate_transform_2d.py maps): W != I, yet flat. *)
  chk["(9a) conformal map f = z + e z^2: pullback metric J^T J is FLAT (Riemann = 0)  (W=I but still: flat)",
      With[{J = D[{x + ee (x^2 - y^2), y + ee (2 x y)}, {{x, y}}]},
        Simplify[Flatten[riem[Transpose[J].J, {x, y}]]] === ConstantArray[0, 2^4]]];
  chk["(9b) NON-conformal map f = (3/2 x + 1/5 y^2, y + 1/10 x^2): pullback metric is FLAT (W != I, yet Riemann = 0)",
      With[{J = D[{3/2 x + 1/5 y^2, y + 1/10 x^2}, {{x, y}}]},
        Simplify[Flatten[riem[Transpose[J].J, {x, y}]]] === ConstantArray[0, 2^4]]];

  Print[If[ok, "ALL PASS", "SOME FAILED"]];
];
