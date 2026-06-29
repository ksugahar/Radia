#!/usr/bin/env wolframscript
(* ::Package:: *)

(*
  derive_quad_q2_henrotte.wls

  COEFFICIENT-FORM derivation: output only K_phi[i,j] and M_sigma_phi[i,j]
  in the monomial basis (no Vandermonde substitution). C++ runtime computes
  Vandermonde inverse and assembles the nodal/V-DOF matrices.

  Why: substituting the symbolic Vandermonde inverse into the energy creates
  expression explosion. Doing the inversion numerically at runtime is cheap
  (LU on 9x9 ~ microseconds) and keeps the symbolic part compact.

  Convention: Hessian-of-W (matches axifem_quad.py:274 reference).
  W_mag = (1/(2 mu)) Integrate(|grad A_phi|^2, 2 pi r dr dz)
        = Integrate(p_s^2/(2 pi muZ) + p_z^2/(8 pi muR s), ds dz)  for p = phi
  W_kin = (sigma/2) Integrate(A_phi^2, 2 pi r dr dz)
        = sigma/(8 pi) Integrate(phi^2/s, ds dz)
  Hessian K_ij = d^2 W_mag / dc_i dc_j = 2 * (W coefficient) for the matrix
  representation K v = lambda M v (same for M).

  GENERAL case (sa > 0):
    Basis = 9 monomials  m[i] = {1, s, s^2, z, sz, s^2 z, z^2, sz^2, s^2 z^2}
    K_phi[i,j] = Integrate(
        d m[i]/ds * d m[j]/ds / (Pi muZ)
      + d m[i]/dz * d m[j]/dz / (4 Pi muR s),
      {s, sa, sb}, {z, za, zb})            (Hessian = 2 x W coefficient)
    M_sigma_phi[i,j] = Integrate(
        m[i] * m[j] * sigma / (4 Pi s),
      {s, sa, sb}, {z, za, zb})            (Hessian = 2 x W coefficient)

  AXIS case (sa = 0):
    Restricted basis (6 monomials, all proportional to s):
      ma[i] = {s, s^2, sz, s^2 z, sz^2, s^2 z^2}
    Same integral formulas with sa = 0.  No 1/s singularity since basis ~ s.

  Output JSON contains 81 + 81 + 36 + 36 = 234 expressions in InputForm.
  C++ side parses these via a tiny Python codegen.
*)

Print["=== Q2 Henrotte K_phi / M_sigma_phi (coefficient form) ==="];
Print[""];
$startTime = AbsoluteTime[];

ClearAll[s, z, sa, sb, za, zb, muR, muZ, sigma];

(* === GENERAL: 9 monomials === *)
mons = {1, s, s^2, z, s z, s^2 z, z^2, s z^2, s^2 z^2};

Print["[GEN] Computing K_phi (9x9, symbolic integrals)..."];
$t0 = AbsoluteTime[];
KphiG = Table[
  Module[{mi = mons[[i]], mj = mons[[j]], dmis, dmjs, dmiz, dmjz, intK, intM},
    dmis = D[mi, s]; dmjs = D[mj, s];
    dmiz = D[mi, z]; dmjz = D[mj, z];
    intK = Integrate[dmis * dmjs / (Pi muZ),
                     {s, sa, sb}, {z, za, zb},
                     Assumptions -> {sb > sa > 0, zb > za}];
    intM = Integrate[dmiz * dmjz / (4 Pi muR s),
                     {s, sa, sb}, {z, za, zb},
                     Assumptions -> {sb > sa > 0, zb > za}];
    Expand[intK + intM]
  ],
  {i, 1, 9}, {j, 1, 9}];
Print["    done ", Round[AbsoluteTime[] - $t0, 0.1], " s"];

Print["[GEN] Computing M_sigma_phi (9x9)..."];
$t0 = AbsoluteTime[];
MsigPhiG = Table[
  Expand@Integrate[mons[[i]] * mons[[j]] * sigma / (4 Pi s),
                   {s, sa, sb}, {z, za, zb},
                   Assumptions -> {sb > sa > 0, zb > za, sigma > 0}],
  {i, 1, 9}, {j, 1, 9}];
Print["    done ", Round[AbsoluteTime[] - $t0, 0.1], " s"];

(* Numeric sanity for general case *)
testRulesG = {sa -> 1*10^-6, sb -> 4*10^-6, za -> 0, zb -> 1*10^-3,
              muR -> 4*Pi*10^-7, muZ -> 4*Pi*10^-7, sigma -> 5.8*10^7};
KphiGnum = N[KphiG /. testRulesG];
MsigPhiGnum = N[MsigPhiG /. testRulesG];
Print["  GEN K_phi symmetric? ", N[Norm[KphiGnum - Transpose[KphiGnum]] / Norm[KphiGnum]] < 1*10^-12];
Print["  GEN M_sigma_phi symmetric? ", N[Norm[MsigPhiGnum - Transpose[MsigPhiGnum]] / Norm[MsigPhiGnum]] < 1*10^-12];
Print["  GEN K_phi[2,2] = ", KphiG[[2, 2]]];

(* === AXIS: 6 restricted monomials === *)
monsA = {s, s^2, s z, s^2 z, s z^2, s^2 z^2};

Print[""];
Print["[AXIS] Computing K_phi_axis (6x6) at sa=0..."];
$t0 = AbsoluteTime[];
KphiA = Table[
  Module[{mi = monsA[[i]], mj = monsA[[j]], dmis, dmjs, dmiz, dmjz, intK, intM},
    dmis = D[mi, s]; dmjs = D[mj, s];
    dmiz = D[mi, z]; dmjz = D[mj, z];
    intK = Integrate[dmis * dmjs / (Pi muZ),
                     {s, 0, sb}, {z, za, zb},
                     Assumptions -> {sb > 0, zb > za}];
    (* 1/s singularity gone: dmiz, dmjz both have factor s, so dmiz*dmjz/s is polynomial *)
    intM = Integrate[dmiz * dmjz / (4 Pi muR s),
                     {s, 0, sb}, {z, za, zb},
                     Assumptions -> {sb > 0, zb > za}];
    Expand[intK + intM]
  ],
  {i, 1, 6}, {j, 1, 6}];
Print["    done ", Round[AbsoluteTime[] - $t0, 0.1], " s"];

Print["[AXIS] Computing M_sigma_phi_axis (6x6)..."];
$t0 = AbsoluteTime[];
MsigPhiA = Table[
  Expand@Integrate[monsA[[i]] * monsA[[j]] * sigma / (4 Pi s),
                   {s, 0, sb}, {z, za, zb},
                   Assumptions -> {sb > 0, zb > za, sigma > 0}],
  {i, 1, 6}, {j, 1, 6}];
Print["    done ", Round[AbsoluteTime[] - $t0, 0.1], " s"];

testRulesA = {sb -> 4*10^-6, za -> 0, zb -> 1*10^-3,
              muR -> 4*Pi*10^-7, muZ -> 4*Pi*10^-7, sigma -> 5.8*10^7};
KphiAnum = N[KphiA /. testRulesA];
MsigPhiAnum = N[MsigPhiA /. testRulesA];
Print["  AXIS K_phi symmetric? ", N[Norm[KphiAnum - Transpose[KphiAnum]] / Norm[KphiAnum]] < 1*10^-12];
Print["  AXIS M_sigma_phi symmetric? ", N[Norm[MsigPhiAnum - Transpose[MsigPhiAnum]] / Norm[MsigPhiAnum]] < 1*10^-12];

(* === Output expressions for codegen === *)
Print[""];
Print["=== Saving symbolic matrices (InputForm strings) ==="];

(* Convert each entry to InputForm string. Avoid SetDelayed pattern function. *)
KphiGstrings = Flatten[Table[ToString[InputForm[KphiG[[i, j]]]], {i, 1, 9}, {j, 1, 9}]];
MsigPhiGstrings = Flatten[Table[ToString[InputForm[MsigPhiG[[i, j]]]], {i, 1, 9}, {j, 1, 9}]];
KphiAstrings = Flatten[Table[ToString[InputForm[KphiA[[i, j]]]], {i, 1, 6}, {j, 1, 6}]];
MsigPhiAstrings = Flatten[Table[ToString[InputForm[MsigPhiA[[i, j]]]], {i, 1, 6}, {j, 1, 6}]];

Print["  K_phi_general_strings type = ", Head[KphiGstrings], ", length = ", Length[KphiGstrings]];
Print["  K_phi_general_strings[[2]] = ", KphiGstrings[[2]]];

result = <|
  "monomials_general" ->
    "{1, s, s^2, z, s*z, s^2*z, z^2, s*z^2, s^2*z^2}",
  "monomials_axis" ->
    "{s, s^2, s*z, s^2*z, s*z^2, s^2*z^2}",
  "node_order_general" ->
    "0:(sa,za) 1:(sb,za) 2:(sb,zb) 3:(sa,zb) 4:(sm,za) 5:(sb,zm) 6:(sm,zb) 7:(sa,zm) 8:(sm,zm); sm=(sa+sb)/2, zm=(za+zb)/2",
  "node_order_axis" ->
    "0:(sb,za) 1:(sb,zb) 2:(sm,za) 3:(sb,zm) 4:(sm,zb) 5:(sm,zm); sm=sb/2",
  "K_phi_general_strings" -> KphiGstrings,
  "M_sigma_phi_general_strings" -> MsigPhiGstrings,
  "K_phi_axis_strings" -> KphiAstrings,
  "M_sigma_phi_axis_strings" -> MsigPhiAstrings,
  "test_params_general" -> "sa=1e-6, sb=4e-6, za=0, zb=1e-3, mu0, sigma_Cu",
  "test_params_axis" -> "sb=4e-6, za=0, zb=1e-3, mu0, sigma_Cu"
|>;
(* Save numerical test matrices to separate JSON, plain numbers only *)
KphiGtest = N[Flatten[KphiGnum]];
MsigPhiGtest = N[Flatten[MsigPhiGnum]];
KphiAtest = N[Flatten[KphiAnum]];
MsigPhiAtest = N[Flatten[MsigPhiAnum]];
testResults = <|
  "K_phi_general_test_flat" -> KphiGtest,
  "M_sigma_phi_general_test_flat" -> MsigPhiGtest,
  "K_phi_axis_test_flat" -> KphiAtest,
  "M_sigma_phi_axis_test_flat" -> MsigPhiAtest
|>;
Export["W:/30_CauerLadderNetwork/2026_04_01_\:9577\:65b9\:5f62CLN/axifem/quad_q2_henrotte_matrices_test.json", testResults, "JSON"];
Export["W:/30_CauerLadderNetwork/2026_04_01_\:9577\:65b9\:5f62CLN/axifem/quad_q2_henrotte_matrices.json", result, "JSON"];
Print["Saved: quad_q2_henrotte_matrices.json"];
Print[""];
Print["Total elapsed: ", Round[AbsoluteTime[] - $startTime, 0.1], " s"];
