#!/usr/bin/env wolframscript
(* ::Package:: *)

(*
  bem_cln_multiport.wls

  Phase 4 Task C: multi-port admittance matrix Y_ij(s) generalisation
  of the single-port BEM-CLN framework.

  Setup: N conductors at positions r_1, ..., r_N in 3D, each with
  polarizability tensor alpha^(n)_ii(s) per principal axis i.  The
  system has 3N "ports": each conductor x 3 polarization axes
  (B_x, B_y, B_z applied to each).

  Coupled system (3N x 3N):
    c^(n)_i = alpha^(n)_i(s) * (B^applied_i + sum_{m,j != (n,i)} G^(nm)_ij c^(m)_j)
  where G^(nm)_ij is the 3D dipole-dipole tensor between conductors
  n and m:
    G^(nm)_ij = (mu_0 / (4 pi r_nm^5)) (3 r_nm,i r_nm,j - delta_ij r_nm^2)
  for r_nm = r_n - r_m (n != m), and zero on the diagonal n = m.

  Multi-port admittance matrix:
    Y_ij(s) = induced m_i per applied B_j  (in some normalised units)
  In our convention with c_n = alpha_n B^local_n, we extract:
    Y^(n,m)_ij = partial(c^(n)_i) / partial(B^applied_j_at_m)
                                                 (= matrix element)

  Also includes:
    Phase 4 B1 (bipolar 2-cyl reference): documented in this script
    as a corollary of the Phase 2.5 polarizability formula c = alpha B /
    (1 - alpha/D^2), which is the geometric series of all dipole images
    (= exact dipole-order 2-cylinder solution; multipole correction
    discussed inline).

  Run: wolframscript -file bem_cln_multiport.wls
*)

ClearAll[a, b, c, sigma, mu, dParam, basisN];

a       = 5*^-3;
b       = 2*^-3;
c       = 1*^-3;
sigma   = 5.80*^7;
mu      = 4 Pi 1*^-7;
dParam  = 2 Pi 1*^5;
basisN  = 3;

VCub    = a b c;
SCub    = 2 (a b + b c + c a);
KsibcCub = SCub Sqrt[sigma/mu];

Print["=== bem_cln_multiport.wls ==="];
Print["Phase 4 C: multi-port admittance matrix Y_ij(s) for BEM-CLN"];
Print["Cuboid: ", a 1000, " x ", b 1000, " x ", c 1000, " mm Cu, V = ", N[VCub, 6]];
Print[""];

(* === Aharoni demag factors ============================================= *)
aharoniNz[ax_, bx_, cx_] := Module[
  {A = ax/2, B = bx/2, C = cx/2, R, Rab, Rac, Rbc},
  R   = Sqrt[A^2 + B^2 + C^2]; Rab = Sqrt[A^2 + B^2];
  Rac = Sqrt[A^2 + C^2];        Rbc = Sqrt[B^2 + C^2];
  (1/Pi) (
    (B^2 - C^2)/(2 B C) Log[(R - A)/(R + A)] +
    (A^2 - C^2)/(2 A C) Log[(R - B)/(R + B)] +
    B/(2 C) Log[(Rab + A)/(Rab - A)] + A/(2 C) Log[(Rab + B)/(Rab - B)] +
    C/(2 A) Log[(Rbc - B)/(Rbc + B)] + C/(2 B) Log[(Rac - A)/(Rac + A)] +
    2 ArcTan[A B / (C R)] + (A^3 + B^3 - 2 C^3)/(3 A B C) +
    (A^2 + B^2 - 2 C^2)/(3 A B C) R +
    C/(A B) (Rac + Rbc) - (Rab^3 + Rbc^3 + Rac^3)/(3 A B C))];

Nx = N[aharoniNz[b, c, a], 30];
Ny = N[aharoniNz[a, c, b], 30];
Nz = N[aharoniNz[a, b, c], 30];

Print["Demag factors: N_x=", N[Nx, 4], ", N_y=", N[Ny, 4], ", N_z=", N[Nz, 4],
      ", sum=", N[Nx+Ny+Nz, 6]];
Print[""];

(* === Per-cuboid building blocks ======================================== *)
Mmax = 13;
modeList = Flatten[Table[{m, n, p},
  {m, 1, Mmax, 2}, {n, 1, Mmax, 2}, {p, 1, Mmax, 2}], 2];
g2MNP[m_, n_, p_]      := 512 a b c / ((m n p)^2 Pi^6);
lambdaMNP[m_, n_, p_]  := (m Pi/a)^2 + (n Pi/b)^2 + (p Pi/c)^2;
tauMNP[m_, n_, p_]     := mu sigma / lambdaMNP[m, n, p];

g2Vec  = N[Table[g2MNP[mn[[1]], mn[[2]], mn[[3]]], {mn, modeList}], 30];
tauVec = N[Table[tauMNP[mn[[1]], mn[[2]], mn[[3]]], {mn, modeList}], 30];

YcubFoster[s_] := sigma Total[g2Vec / (1 + s tauVec)];
chiFunc[s_]    := YcubFoster[s] / (sigma VCub) - 1;
alphaAxis[s_, Ni_] := VCub chiFunc[s] / (1 + Ni chiFunc[s]);

(* === 3D dipole-dipole coupling tensor =================================== *)
(* G_ij = (mu_0 / 4 pi r^5) (3 r_i r_j - delta_ij r^2) for r != 0 *)
gTensor[rvec_] := Module[{r = Norm[rvec], rij},
  rij = Outer[Times, rvec, rvec];
  (mu/(4 Pi r^5)) (3 rij - r^2 IdentityMatrix[3])];

(* === Multi-port system assembly ========================================= *)
(*
  N conductors at positions positionsList[[n]] (n=1..Nc).
  Each conductor has polarizability tensor alpha^(n) = diag(alpha_x, alpha_y, alpha_z).
  System (3 N x 3 N):
    c_n_i = alpha^(n)_i (B^app_i + sum_{m != n, j} G^(nm)_ij c^(m)_j)

  Equivalent in matrix form: (I + alpha_block . G_block) c_vec = alpha_block . B_vec
*)

buildMultiportSystem[positionsList_, s_, NfacList_] := Module[
  {Nc, alphaList, alphaBlock, Gblock, sysMat},
  Nc = Length[positionsList];
  (* Per-conductor polarizability tensors *)
  alphaList = Table[
    DiagonalMatrix[alphaAxis[s, #] & /@ NfacList],
    {Nc}];
  (* Block-diagonal alpha *)
  alphaBlock = ArrayFlatten[
    Table[If[i == j, alphaList[[i]], 0 IdentityMatrix[3]],
      {i, 1, Nc}, {j, 1, Nc}]];
  (* Coupling matrix G_block: off-diagonal blocks are dipole-dipole *)
  Gblock = ArrayFlatten[
    Table[If[i == j, 0 IdentityMatrix[3],
                     gTensor[positionsList[[i]] - positionsList[[j]]]],
      {i, 1, Nc}, {j, 1, Nc}]];
  sysMat = IdentityMatrix[3 Nc] + alphaBlock . Gblock;
  {sysMat, alphaBlock}];

(* Solve at one frequency for one applied B direction *)
solveMultiportOneB[positionsList_, s_, NfacList_, BappliedVec_] := Module[
  {sysMat, alphaBlock, rhs, cvec, Nc},
  {sysMat, alphaBlock} = buildMultiportSystem[positionsList, s, NfacList];
  Nc = Length[positionsList];
  rhs = alphaBlock . Flatten[Table[BappliedVec, {Nc}]];
  cvec = LinearSolve[sysMat, rhs];
  Partition[cvec, 3]];

(* Build the full Y_ij(s) port admittance matrix: drive each port with B=1
   in one direction, observe response c at all ports.  Result is a
   (3 Nc) x 3 matrix per applied-direction column. *)
buildYportMatrix[positionsList_, s_, NfacList_] := Module[
  {Nc, Yport, cvec, BvecX, BvecY, BvecZ},
  Nc = Length[positionsList];
  (* For each applied B direction, solve and stack the resulting c *)
  BvecX = {1, 0, 0}; BvecY = {0, 1, 0}; BvecZ = {0, 0, 1};
  Yport = Table[
    Flatten[solveMultiportOneB[positionsList, s, NfacList, Bvec]],
    {Bvec, {BvecX, BvecY, BvecZ}}];
  Transpose[Yport]];

(* === Test case: 2 cuboids at (-D/2, 0, 0) and (+D/2, 0, 0) ============== *)
Print["Multi-port admittance matrix for 2-cuboid system:"];
Print["  Conductor 1 at (-7.5, 0, 0) mm"];
Print["  Conductor 2 at (+7.5, 0, 0) mm"];
Print[""];

positions2 = {{-7.5e-3, 0, 0}, {7.5e-3, 0, 0}};
NfacList   = {Nx, Ny, Nz};

(* Probe at f = 10^6 Hz *)
With[{sProbe = I 2 Pi 1*^6},
  With[{Yp = buildYportMatrix[positions2, sProbe, NfacList]},
    Print["Y_port (6 x 3 matrix, rows = output port (cuboid 1 xyz, cuboid 2 xyz), cols = applied B axis):"];
    Print["  Y[1,1] (m_x_cyl1 per B_x_applied) = ", N[Yp[[1,1]], 4]];
    Print["  Y[2,2] (m_y_cyl1 per B_y_applied) = ", N[Yp[[2,2]], 4]];
    Print["  Y[3,3] (m_z_cyl1 per B_z_applied) = ", N[Yp[[3,3]], 4]];
    Print["  Y[4,1] (m_x_cyl2 per B_x_applied) = ", N[Yp[[4,1]], 4]];
    Print["  Y[5,2] (m_y_cyl2 per B_y_applied) = ", N[Yp[[5,2]], 4]];
    Print["  Y[6,3] (m_z_cyl2 per B_z_applied) = ", N[Yp[[6,3]], 4]];
    Print[""];
    Print["By symmetry of 2-cuboid system, off-diagonal blocks reveal coupling:"];
    Print["  Y[4,2] (m_x_cyl2 per B_y_applied, cross-coupling) = ", N[Yp[[4,2]], 4]];
    Print["  Y[5,1] (m_y_cyl2 per B_x_applied) = ", N[Yp[[5,1]], 4]];
  ]];

(* === Symmetry / reciprocity check ===================================== *)
Print[""];
Print["Reciprocity: Y_ij = Y_ji (linear passive system)"];
With[{sProbe = I 2 Pi 1*^6},
  With[{Yp = buildYportMatrix[positions2, sProbe, NfacList]},
    Print["  Y - Y^T max norm = ",
          N[Max[Abs[Yp - Transpose[Yp]]], 6]];
    Print["  (target: 0 for reciprocal system)"];
  ]];

(* === Phase 4 B1 (Bipolar 2-cyl exact reference) ========================= *)
Print[""];
Print["=== Phase 4 B1: Bipolar 2-cyl exact reference status ==="];
Print[""];
Print["The Phase 2.5 polarizability framework"];
Print["    c_n = alpha(s) B_n^applied / (1 - alpha(s)/D^2)"];
Print["is itself the EXACT dipole-order 2-cyl coupled solution"];
Print["(geometric series of all dipole images sum to this closed form)."];
Print[""];
Print["The 'next-order' correction beyond Phase 2.5 is the quadrupole"];
Print["coupling: the field from cyl 1 at cyl 2's location varies"];
Print["across cyl 2's cross-section, inducing a quadrupole moment in"];
Print["cyl 2.  Including quadrupole moments gives a 4 x 4 coupled"];
Print["system per cylinder pair, with corrections O((a/D)^2) of the"];
Print["dipole result."];
Print[""];
Print["For 5 mm cylinder at D = 15 mm separation: a/D = 1/3,"];
Print["(a/D)^2 = 0.11, so quadrupole correction is ~10 percent."];
Print["This is the residual gap of Phase 2.5 vs full bipolar exact."];
Print["Multipole-corrected reference is left as future work."];

(* === CSV export of Y matrix at several frequencies ====================== *)
Print[""];
Print["Writing multi-port Y matrix CSV..."];
csvData = {{"freq_Hz", "Y_11_re", "Y_11_im", "Y_22_re", "Y_22_im",
            "Y_33_re", "Y_33_im", "Y_44_re", "Y_44_im",
            "Y_55_re", "Y_55_im", "Y_66_re", "Y_66_im"}};
Do[
  With[{omega = 2 Pi (10.0^logf), s = I 2 Pi (10.0^logf)},
    With[{Yp = buildYportMatrix[positions2, s, NfacList]},
      AppendTo[csvData,
        {10.0^logf,
         N[Re[Yp[[1,1]]], 6], N[Im[Yp[[1,1]]], 6],
         N[Re[Yp[[2,2]]], 6], N[Im[Yp[[2,2]]], 6],
         N[Re[Yp[[3,3]]], 6], N[Im[Yp[[3,3]]], 6],
         N[Re[Yp[[4,4]]], 6], N[Im[Yp[[4,4]]], 6],
         N[Re[Yp[[5,5]]], 6], N[Im[Yp[[5,5]]], 6],
         N[Re[Yp[[6,6]]], 6], N[Im[Yp[[6,6]]], 6]}]
    ]],
  {logf, {3, 5, 7, 9, 12}}];

Export["bem_cln_multiport.csv", csvData, "CSV"];
Print["Wrote bem_cln_multiport.csv"];
Print[""];
Print["=== bem_cln_multiport.wls DONE ==="];
