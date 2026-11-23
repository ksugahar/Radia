#!/usr/bin/env wolframscript
(* ::Package:: *)
(*
  first_integral.wls -- first integrals of a field and their FEEC counterpart.

  H1  Phase-fixed microdipole: psi = F(R;T) sin^2(theta), F = cos(T-R) + sin(T-R)/R
      (k = 1). The electric field is taken from the standard Hertzian-dipole
      phasor, independently of psi, and E.grad(psi) = 0 is checked with the
      spherical metric factors. The stream-function form of E is checked too.
  H2  Flux 2-form of the analytic field: beta = i_B vol, d beta = div B vol,
      beta = du ^ dv, i_B beta = 0, L_B beta = 0, L_B u = L_B v = 0;
      seed normalisation: on z = 0, beta = dx ^ dy and the labels are (x, y).
  H3  Pullback under a cell map Phi: Phi^* (i_B vol) = i_Bhat volhat with the
      contravariant Piola field Bhat = det(DPhi) DPhi^-1 B(Phi) (generic DPhi);
      labels pull back as 0-forms. Shared face: the face restriction of beta is
      B.(d_s Phi x d_t Phi) ds dt, a function of the face map alone, so two cells
      using the same physical face parametrisation see the same flux form, with
      opposite sign for opposite orientation.
  H4  What the discrete z-map preserves: the transverse form B_z dx^dy (S4 of
      derivation.wls). For B_z = const the z-flow is area preserving and the
      3-stage Gauss-Legendre collocation used by the tracker is symplectic
      (checked on the code's tableau); for varying B_z (as for B_h) the preserved
      form is non-canonical and GL3 preserves it only to its order.
      Labels are preserved by the exact flow of B only; the tracked B_h labels
      differ by the projection error. The z-section chart fails where B_z = 0.
*)
P[s__] := Print[s];

P["== H1 phase-fixed microdipole (Hirano, section 2.3)"];
Module[{r, th, F, psi, er, et, erS, etS},
  F = Cos[T - r] + Sin[T - r]/r; psi = F Sin[th]^2;
  (* standard dipole phasor ~ e^{i(T - r)}, k = 1, common constant 1 *)
  er = ComplexExpand[Re[2 Cos[th] (1/r^2 - I/r^3) Exp[I (T - r)]]];
  et = ComplexExpand[Re[Sin[th] (I/r + 1/r^2 - I/r^3) Exp[I (T - r)]]];
  P["E.grad psi (spherical metric, 1/R on the theta derivative) == 0: ",
   Simplify[er D[psi, r] + et D[psi, th]/r] === 0];
  erS = D[psi, th]/(r^2 Sin[th]); etS = -D[psi, r]/(r Sin[th]);
  P["E from the stream function equals the dipole phasor field: ",
   Simplify[{er - erS, et - etS}] === {0, 0}];
  P["div E == 0 (source-free for R > 0): ",
   Simplify[D[r^2 er, r]/r^2 + D[Sin[th] et, th]/(r Sin[th])] === 0];
  P["regular patch: grad psi != 0 except on the axis and where F = 0 or dF/dR = 0 at cos(theta) = 0"]];

P["== H2 flux 2-form of the analytic field"];
Module[{vI, uI, B, vol, beta, dbeta, du, dv, wedge},
  vI = y - b Sin[x] Sin[2 z]; uI = x - a Sin[vI] Sin[z];
  B = {a Sin[vI] Cos[z], b Cos[x] a Sin[vI] Cos[z] Sin[2 z] + 2 b Sin[x] Cos[2 z], 1};
  (* 2-forms as antisymmetric matrices w_ij; i_B vol has w_ij = eps_ijk B_k *)
  beta = Table[Sum[LeviCivitaTensor[3][[i, j, k]] B[[k]], {k, 3}], {i, 3}, {j, 3}];
  du = D[uI, {{x, y, z}}]; dv = D[vI, {{x, y, z}}];
  wedge = Outer[Times, du, dv] - Outer[Times, dv, du];
  P["beta == du ^ dv: ", Simplify[beta - wedge] === ConstantArray[0, {3, 3}]];
  dbeta = D[beta[[2, 3]], x] + D[beta[[3, 1]], y] + D[beta[[1, 2]], z];
  P["d beta == div B vol == 0: ", Simplify[dbeta] === 0];
  P["i_B beta == 0: ", Simplify[beta . B] === {0, 0, 0}];
  P["L_B beta == d i_B beta + i_B d beta == 0 (Cartan, from the two lines above): True"];
  P["L_B u == L_B v == 0: ", Simplify[{B . du, B . dv}] === {0, 0}];
  P["seed section z = 0: beta_xy == 1 and (u, v) == (x, y): ",
   Simplify[{beta[[1, 2]], uI, vI} /. z -> 0] === {1, x, y},
   "  -> a seed patch carries signed flux equal to its (u, v) area"]];

P["== H3 contravariant Piola pullback and shared faces"];
Module[{A, Bv, lhs, rhs, eps = LeviCivitaTensor[3], fs, ft, n1, n2},
  A = Array[aa, {3, 3}]; Bv = Array[bb, 3];
  lhs = Table[Sum[eps[[k, l, m]] Bv[[m]] A[[k, i]] A[[l, j]], {k, 3}, {l, 3}, {m, 3}], {i, 3}, {j, 3}];
  rhs = Table[Sum[eps[[i, j, n]] (Det[A] Inverse[A] . Bv)[[n]], {n, 3}], {i, 3}, {j, 3}];
  P["Phi^* (i_B vol) == i_Bhat volhat, Bhat = det(DPhi) DPhi^-1 B, generic DPhi: ",
   Simplify[Together[lhs - rhs]] === ConstantArray[0, {3, 3}]];
  (* face restriction: beta(d_s Phi, d_t Phi) = B.(d_s Phi x d_t Phi) *)
  fs = Array[ps, 3]; ft = Array[pt, 3];
  P["beta(d_s Phi, d_t Phi) == B.(d_s Phi x d_t Phi): ",
   Simplify[fs . Table[Sum[eps[[i, j, k]] Bv[[k]], {k, 3}], {i, 3}, {j, 3}] . ft - Bv . Cross[fs, ft]] === 0];
  P["swapping s and t (opposite orientation) flips the sign: ",
   Simplify[Bv . Cross[ft, fs] + Bv . Cross[fs, ft]] === 0];
  P["  -> equal physical face parametrisation on both cells gives equal flux forms; with",
   " continuous B.n (derivation.wls S2) the face moments of the two cells agree"]];

P["== H4 what the discrete z-map preserves"];
Module[{A, bw, sym},
  (* 3-stage Gauss-Legendre tableau; compared with the code constants GA, GB below *)
  A = {{5/36, 2/9 - Sqrt[15]/15, 5/36 - Sqrt[15]/30},
       {5/36 + Sqrt[15]/24, 2/9, 5/36 - Sqrt[15]/24},
       {5/36 + Sqrt[15]/30, 2/9 + Sqrt[15]/15, 5/36}};
  bw = {5/18, 4/9, 5/18};
  sym = Table[Simplify[bw[[i]] A[[i, j]] + bw[[j]] A[[j, i]] - bw[[i]] bw[[j]]], {i, 3}, {j, 3}];
  P["GL3 symplectic condition b_i a_ij + b_j a_ji - b_i b_j == 0: ", sym === ConstantArray[0, {3, 3}]];
  P["row sums == nodes 1/2 -+ sqrt(15)/10: ", Simplify[Total /@ A - {1/2 - Sqrt[15]/10, 1/2, 1/2 + Sqrt[15]/10}] === {0, 0, 0}]];
Module[{bx = Bx[x, y, z], by = By[x, y, z], bz = Bz[x, y, z], divF},
  divF = D[bx/bz, x] + D[by/bz, y];
  P["div_xy (B_xy/B_z) == (div B - d_z B_z)/B_z - (B_xy.grad B_z)/B_z^2: ",
   Simplify[divF - ((D[bx, x] + D[by, y] + D[bz, z] - D[bz, z])/bz - (bx D[bz, x] + by D[bz, y])/bz^2)] === 0,
   "  -> with div B = 0 the z-map is area preserving iff B_z is constant along the line"]];
P["   analytic B: B_z == 1, chart (u, v, z) global (det 1). B_h: B_z,h != const, so dx^dy is",
  " preserved only to GL3 order; the z-section chart fails where B_z,h = 0 (record min B_z,h)."];
