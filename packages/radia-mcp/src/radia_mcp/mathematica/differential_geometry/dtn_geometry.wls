#!/usr/bin/env wolframscript
(* ::Package:: *)

(*
  dtn_geometry.wls -- the Dirichlet-to-Neumann (DtN / Steklov-Poincare) operator in the
  language of differential geometry.  Capstone of the differential_geometry/ suite: the DtN
  is the EXTERIOR's Hodge star condensed onto the boundary -- the "*" (geometry) object that
  Kelvin/IE compute exactly and HOIBC approximates.  Companion to weakform_hodge.wls (the
  Hodge half), surface_derham.wls (HOIBC + the surface DtN), and docs/.../PREMETRIC_PICTURE.md.

  Three differential-geometric characterisations, all verified here:

  [1] CONDENSATION (Schur complement).  The DtN is the boundary condensation of the bulk
      Hodge-Laplace form a(u,v) = INT <du, * dv>: eliminate the interior, keep the boundary.
      Operationally  Lambda = A_BB - A_BI A_II^{-1} A_IB  (a positive boundary operator).
      Continuously it is the conormal trace of the condensed exterior solution u = (R/r)^{n+1},
      and the radial DtN m(r) = -u'/u is the decaying fixed point of the RICCATI flow
      m' = m^2 - (2/r) m - n(n+1)/r^2  (the infinitesimal Schur complement of a thin shell).

  [2] PSEUDODIFFERENTIAL, order +1 = a shifted SQUARE ROOT of the surface Laplacian.  On the
      sphere  Lambda_n = (n+1)/R  while -Delta_Gamma Y_n = n(n+1)/R^2 Y_n, and EXACTLY
          Lambda = (1/2 + sqrt( R^2 (-Delta_Gamma) + 1/4 )) / R .
      The DtN is essentially sqrt(-Delta_Gamma) -- NONLOCAL.  That is why HOIBC (a local
      polynomial in Delta_Gamma) only matches finitely many modes (the knee), while Kelvin/IE
      (which solve the exterior) carry the sqrt exactly.

  [3] SELF-ADJOINT and POSITIVE.  Green's identity (the symmetry of the Hodge pairing
      <du,*dv> = <dv,*du>) gives INT_Gamma (u Lambda v - v Lambda u) = 0, so Lambda = Lambda* >= 0:
      the DtN defines the boundary H^{1/2} inner product -- it is the boundary's "metric".  Its
      spectrum is the Steklov spectrum; eigenfunctions = boundary de Rham harmonics (the two
      vector ladders n/R, (n+1)/R of surface_derham.wls [C]).

  In the premetric picture (PREMETRIC_PICTURE.md): the DtN sits in the "*" column (geometry),
  basis = boundary de Rham ("d").  It is not a new object -- it is the exterior Hodge, condensed.
*)

uExt[n_] := (R/r)^(n + 1);                                   (* exterior solution condensed to u(R)=1 *)
riccatiRHS[m_, rr_, n_] := m^2 - (2/rr) m - n (n + 1)/rr^2;  (* the radial DtN/Schur Riccati RHS *)
dtnSqrt[lap_] := (1/2 + Sqrt[lap R^2 + 1/4])/R;             (* DtN as a function of the L-B eigenvalue lap *)

(* ============================ self-test ============================ *)
Module[{ok = True, chk, aBB, aBI, aII, S},
  chk[name_, cond_] := (Print["  ", If[cond, "ok  ", "FAIL"], name]; ok = ok && cond);
  Print["dtn_geometry.wls self-test (the DtN operator in differential-geometry language):"];

  (* ---- [1] DtN = boundary CONDENSATION (Schur complement) ---- *)
  (* the condensed exterior solution's conormal trace IS the DtN:  -u'(R)/u(R) = (n+1)/R *)
  chk["[1a] DtN = conormal trace of the condensed exterior solution: -u'(R)/u(R) = (n+1)/R",
      And @@ Table[Simplify[(-D[uExt[n], r]/uExt[n] /. r -> R) - (n + 1)/R] === 0, {n, 0, 4}]];
  (* the radial DtN m=(n+1)/r is the decaying fixed point of the Riccati (the differential Schur step) *)
  chk["[1b] m = (n+1)/r solves the radial Schur/Riccati  m' = m^2 - (2/r) m - n(n+1)/r^2 (thin-shell condensation)",
      And @@ Table[Simplify[D[(n + 1)/r, r] - riccatiRHS[(n + 1)/r, r, n]] === 0, {n, 1, 4}]];
  (* the Schur-complement STRUCTURE: condensing an SPD bulk leaves a positive boundary operator *)
  aBB = {{kbb}}; aBI = {{kbi}}; aII = {{kii}};
  S = aBB - aBI.Inverse[aII].Transpose[aBI];   (* Schur complement onto the boundary *)
  chk["[1c] boundary condensation = Schur complement S = A_BB - A_BI A_II^{-1} A_IB, positive for SPD bulk",
      Simplify[S[[1, 1]] - (kbb - kbi^2/kii)] === 0 &&
      Resolve[ForAll[{kbb, kbi, kii}, Implies[kbb > 0 && kii > 0 && kbb kii - kbi^2 > 0,
                                              kbb - kbi^2/kii > 0]], Reals] === True];

  (* ---- [2] DtN = shifted SQUARE ROOT of the surface Laplacian (pseudodifferential, order +1) ---- *)
  chk["[2a] Lambda_n = (1/2 + sqrt(n(n+1)+1/4))/R = (n+1)/R  (exact shifted-sqrt of -Delta_Gamma)",
      And @@ Table[Simplify[(1/2 + Sqrt[n (n + 1) + 1/4])/R - (n + 1)/R] === 0, {n, 0, 5}]];
  chk["[2b] operator form: dtnSqrt[ n(n+1)/R^2 ] = (n+1)/R  (the DtN as a function of the L-B eigenvalue)",
      And @@ Table[Simplify[dtnSqrt[n (n + 1)/R^2] - (n + 1)/R] === 0, {n, 0, 5}]];

  (* ---- [3] SELF-ADJOINT and POSITIVE (Green's identity = the Hodge pairing symmetry) ---- *)
  (* the bulk material-Hodge pairing is symmetric a(u,v)=a(v,u) for a symmetric material tensor *)
  chk["[3a] Hodge pairing is symmetric: grad(u).nu.grad(v) = grad(v).nu.grad(u)  (=> Lambda self-adjoint)",
      With[{nuMat = {{m11, m12, m13}, {m12, m22, m23}, {m13, m23, m33}},
            gu = {ux, uy, uz}, gv = {vx, vy, vz}},
        Simplify[gu.nuMat.gv - gv.nuMat.gu] === 0]];
  (* self-adjointness on the sphere: the antisymmetric boundary form (Lambda_n - Lambda_m) INT Y_m Y_n
     vanishes for every m,n (Legendre orthogonality kills m!=n; the prefactor kills m=n) *)
  chk["[3b] self-adjoint: (Lambda_n - Lambda_m) INT_0^pi P_m P_n sin th dth = 0 for all m,n (Green + orthogonality)",
      And @@ Flatten@Table[
         Simplify[((n + 1)/R - (m + 1)/R)
                  Integrate[LegendreP[m, Cos[th]] LegendreP[n, Cos[th]] Sin[th], {th, 0, Pi}]] === 0,
         {m, 0, 3}, {n, 0, 3}]];
  (* positive: the Steklov spectrum is strictly positive *)
  chk["[3c] positive: the Steklov eigenvalue (n+1)/R > 0 (n>=0, R>0)",
      Resolve[ForAll[{n, R}, Implies[n >= 0 && R > 0, (n + 1)/R > 0]], Reals] === True];

  (* ---- [4] NONLOCAL => HOIBC is necessarily approximate (the sqrt is not a polynomial) ---- *)
  chk["[4] the DtN sqrt is NOT a polynomial in the L-B eigenvalue (3rd derivative != 0): no finite local HOIBC equals it (the knee)",
      Simplify[D[dtnSqrt[lap], {lap, 3}]] =!= 0];

  Print[If[ok, "ALL PASS", "SOME FAILED"]];
];
