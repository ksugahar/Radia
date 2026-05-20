#!/usr/bin/env wolframscript
(* ::Package:: *)

(*
  bem_cln_multiport_v2.wls

  Phase 4 final-1: clean numerical implementation of the multi-port
  admittance matrix Y_ij(s) for N coupled 3D cuboid conductors.
  Replaces bem_cln_multiport.wls (which had symbolic-evaluation
  issues).

  Each conductor has a polarizability tensor diag(alpha_x, alpha_y, alpha_z)
  derived from the demagnetization-corrected formulation
  (polarizability_tensor_3D.wls).  Inter-conductor coupling uses the
  3D dipole-dipole tensor:
     G^(nm)_ij = (mu_0 / (4 pi r^5)) (3 r_i r_j - delta_ij r^2)
  for r = r_n - r_m.

  The full system is 3N x 3N (each "port" is a (conductor, axis) pair).
  Reciprocity: Y(s) = Y(s)^T.

  Run: wolframscript -file bem_cln_multiport_v2.wls
*)

ClearAll[a, b, c, sigma, mu, dParam];

a       = 5*^-3;
b       = 2*^-3;
c       = 1*^-3;
sigma   = N[5.80*^7, 30];
mu      = N[4 Pi 1*^-7, 30];
dParam  = N[2 Pi 1*^5, 30];

VCub    = a b c;
SCub    = 2 (a b + b c + c a);
KsibcCub = SCub Sqrt[sigma/mu];

Print["=== bem_cln_multiport_v2.wls ==="];
Print["3 x 3 polarizability tensor x N cuboid multi-port admittance"];
Print[""];

(* === Aharoni demag factors (precompute to high precision) =============== *)
aharoniNz[ax_, bx_, cx_] := Module[
  {A = N[ax/2, 30], B = N[bx/2, 30], C = N[cx/2, 30], R, Rab, Rac, Rbc},
  R   = Sqrt[A^2 + B^2 + C^2]; Rab = Sqrt[A^2 + B^2];
  Rac = Sqrt[A^2 + C^2];        Rbc = Sqrt[B^2 + C^2];
  N[(1/Pi) (
    (B^2 - C^2)/(2 B C) Log[(R - A)/(R + A)] +
    (A^2 - C^2)/(2 A C) Log[(R - B)/(R + B)] +
    B/(2 C) Log[(Rab + A)/(Rab - A)] + A/(2 C) Log[(Rab + B)/(Rab - B)] +
    C/(2 A) Log[(Rbc - B)/(Rbc + B)] + C/(2 B) Log[(Rac - A)/(Rac + A)] +
    2 ArcTan[A B / (C R)] + (A^3 + B^3 - 2 C^3)/(3 A B C) +
    (A^2 + B^2 - 2 C^2)/(3 A B C) R +
    C/(A B) (Rac + Rbc) - (Rab^3 + Rbc^3 + Rac^3)/(3 A B C)), 30]];

Nx = aharoniNz[b, c, a];
Ny = aharoniNz[a, c, b];
Nz = aharoniNz[a, b, c];
NList = {Nx, Ny, Nz};

Print["Demag factors: N_x=", N[Nx, 4], ", N_y=", N[Ny, 4], ", N_z=", N[Nz, 4]];
Print[""];

(* === Per-cuboid Y_cuboid(s) (Foster sum) =============================== *)
Mmax = 13;
modeList = Flatten[Table[{m, n, p},
  {m, 1, Mmax, 2}, {n, 1, Mmax, 2}, {p, 1, Mmax, 2}], 2];
g2Vec  = N[Table[512 a b c / ((mn[[1]] mn[[2]] mn[[3]])^2 Pi^6), {mn, modeList}], 30];
tauVec = N[Table[
  mu sigma / ((mn[[1]] Pi/a)^2 + (mn[[2]] Pi/b)^2 + (mn[[3]] Pi/c)^2),
  {mn, modeList}], 30];

YcubFoster[s_?NumericQ] := N[sigma Total[g2Vec / (1 + s tauVec)], 20];
chiFunc[s_?NumericQ]    := N[YcubFoster[s] / (sigma VCub) - 1, 20];
alphaAxis[s_?NumericQ, Ni_?NumericQ] := N[VCub chiFunc[s] / (1 + Ni chiFunc[s]), 20];

(* Per-cuboid polarizability tensor at frequency s *)
alphaTensor[s_?NumericQ] := DiagonalMatrix[N[alphaAxis[s, #], 20] & /@ NList];

(* === 3D dipole-dipole coupling tensor =================================== *)
gTensor3D[rvec_] := Module[
  {r2 = N[rvec . rvec, 30], r = Norm[N[rvec, 30]], rij},
  rij = Outer[Times, N[rvec, 30], N[rvec, 30]];
  N[(mu/(4 Pi r^5)) (3 rij - r^2 IdentityMatrix[3]), 20]];

(* === Multi-port system assembly ========================================= *)
buildYportNumerical[positionsList_, s_?NumericQ] := Module[
  {Nc, alphaTensors, Iblock3N, sysMat, Yport, BvecList, ones, rhs, cvec},
  Nc = Length[positionsList];

  (* All cuboids identical: same alpha tensor everywhere *)
  alphaT = alphaTensor[s];
  alphaTensors = Table[alphaT, {Nc}];   (* per-cuboid; identical *)

  (* Coupling matrix G_block (3N x 3N): off-diagonal 3x3 dipole-dipole tensors *)
  GblockEntries = Table[
    If[i == j, ConstantArray[0, {3, 3}],
       gTensor3D[positionsList[[i]] - positionsList[[j]]]],
    {i, 1, Nc}, {j, 1, Nc}];
  Gblock = ArrayFlatten[GblockEntries];

  (* alpha_block: block-diagonal with per-cuboid 3x3 alpha *)
  alphaBlockEntries = Table[
    If[i == j, alphaTensors[[i]], ConstantArray[0, {3, 3}]],
    {i, 1, Nc}, {j, 1, Nc}];
  alphaBlock = ArrayFlatten[alphaBlockEntries];

  sysMat = IdentityMatrix[3 Nc] + alphaBlock . Gblock;

  (* For each applied B direction, solve and stack column vectors *)
  BvecList = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  ones = ConstantArray[1, Nc];
  Yport = Table[
    rhs = alphaBlock . Flatten[Table[Bvec, {Nc}]];
    cvec = LinearSolve[sysMat, rhs];
    cvec,
    {Bvec, BvecList}];

  (* Yport: 3 rows (one per applied B axis), each row a (3 Nc)-vector
     of per-port responses.  Returned as 3 x (3Nc) matrix; the i-th
     row gives the response to B applied in direction i.            *)
  Yport];

(* === Test: 2-cuboid system ============================================== *)
Print["2-cuboid Y_port matrix at f = 10^6 Hz:"];
Print[""];
positions2 = {{-7.5*^-3, 0, 0}, {7.5*^-3, 0, 0}};

With[{sProbe = I 2 Pi 1.0*^6},
  With[{Yp = buildYportNumerical[positions2, sProbe]},
    Print["Yp shape: 3 (applied B axis) x ", Length[Yp[[1]]], " (per-port response)"];
    Print[""];
    Print["Row 1: response to B_applied = (1, 0, 0):"];
    Print["  m_x_cuboid1 (Yp[1,1]) = ", N[Yp[[1, 1]], 5]];
    Print["  m_y_cuboid1 (Yp[1,2]) = ", N[Yp[[1, 2]], 5]];
    Print["  m_z_cuboid1 (Yp[1,3]) = ", N[Yp[[1, 3]], 5]];
    Print["  m_x_cuboid2 (Yp[1,4]) = ", N[Yp[[1, 4]], 5]];
    Print["  m_y_cuboid2 (Yp[1,5]) = ", N[Yp[[1, 5]], 5]];
    Print["  m_z_cuboid2 (Yp[1,6]) = ", N[Yp[[1, 6]], 5]];
    Print[""];
    Print["Row 2: response to B_applied = (0, 1, 0):"];
    Print["  m_y_cuboid1 (Yp[2,2]) = ", N[Yp[[2, 2]], 5]];
    Print["  m_y_cuboid2 (Yp[2,5]) = ", N[Yp[[2, 5]], 5]];
    Print[""];
    Print["Row 3: response to B_applied = (0, 0, 1):"];
    Print["  m_z_cuboid1 (Yp[3,3]) = ", N[Yp[[3, 3]], 5]];
    Print["  m_z_cuboid2 (Yp[3,6]) = ", N[Yp[[3, 6]], 5]];
    Print[""];
    Print["By geometry (cuboids along x), inter-cuboid coupling factor:"];
    Print["  Yp[1,4] / Yp[1,1] (x-x mutual) = ",
          N[Re[Yp[[1, 4]]/Yp[[1, 1]]], 5]];
    Print["  Yp[2,5] / Yp[2,2] (y-y mutual) = ",
          N[Re[Yp[[2, 5]]/Yp[[2, 2]]], 5]];
    Print["  Yp[3,6] / Yp[3,3] (z-z mutual) = ",
          N[Re[Yp[[3, 6]]/Yp[[3, 3]]], 5]]
  ]];

(* === Frequency sweep CSV ================================================ *)
Print[""];
Print["Writing multi-port Y CSV..."];
csvData = {{"freq_Hz",
            "Y_11_re", "Y_11_im",   (* m_x_c1 / B_x *)
            "Y_22_re", "Y_22_im",   (* m_y_c1 / B_y *)
            "Y_33_re", "Y_33_im",   (* m_z_c1 / B_z *)
            "Y_44_re", "Y_44_im",   (* m_x_c2 / B_x *)
            "Y_41_re", "Y_41_im"    (* cross: m_x_c2 from B_x_applied *)
           }};

Do[
  With[{omega = 2 Pi (10.0^logf), s = I 2 Pi (10.0^logf)},
    With[{Yp = buildYportNumerical[positions2, s]},
      AppendTo[csvData,
        {10.0^logf,
         N[Re[Yp[[1, 1]]], 6], N[Im[Yp[[1, 1]]], 6],  (* m_x_c1 from Bx *)
         N[Re[Yp[[2, 2]]], 6], N[Im[Yp[[2, 2]]], 6],  (* m_y_c1 from By *)
         N[Re[Yp[[3, 3]]], 6], N[Im[Yp[[3, 3]]], 6],  (* m_z_c1 from Bz *)
         N[Re[Yp[[1, 4]]], 6], N[Im[Yp[[1, 4]]], 6],  (* m_x_c2 from Bx *)
         N[Re[Yp[[1, 4]] - Yp[[1, 1]]], 6],            (* mutual minus self *)
         N[Im[Yp[[1, 4]] - Yp[[1, 1]]], 6]}]
    ]],
  {logf, {3, 4, 5, 6, 7, 8, 10, 12}}];

Export["bem_cln_multiport_v2.csv", csvData, "CSV"];
Print["Wrote bem_cln_multiport_v2.csv (", Length[csvData] - 1, " freq points)"];

Print[""];
Print["=== bem_cln_multiport_v2.wls DONE ==="];
