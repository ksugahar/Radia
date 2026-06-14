#!/usr/bin/env wolframscript
(* ::Package:: *)

(*
  polarizability_tensor_3D.wls

  Phase 4 A1: 3D polarizability tensor with demagnetization factors
  for a rectangular cuboid conductor.  This addresses the
  demagnetization-factor gap identified in the Phase 3 B.2 NGSolve
  cross-check (§V.H honest scope note).

  Mathematical formulation
  ------------------------
  For a cuboid Omega = [0,a] x [0,b] x [0,c] in applied uniform B in
  direction i, the induced magnetic moment is

    m_i(s) = alpha_i(s) * B_i_applied / mu_0
    alpha_i(s) = V * chi(s) / (1 + N_i * chi(s))

  where:
    V = a b c                            (cuboid volume)
    N_i = demagnetization factor along axis i (Aharoni 1998)
    N_x + N_y + N_z = 1
    chi(s) = Y_cuboid(s) / (sigma V) - 1   (intrinsic susceptibility)

  chi(s) limits:
    chi(0)        = sigma V / (sigma V) - 1 = 0        (no eddy at DC)
    chi(infinity) = 0 - 1 = -1                          (PEC)

  alpha_i(s) limits:
    alpha_i(0)        = 0
    alpha_i(infinity) = V * (-1) / (1 - N_i) = -V/(1 - N_i)

  Aharoni 1998 formula for demagnetization factor of a rectangular
  prism with sides 2a, 2b, 2c centred at origin, along z axis:
    See Aharoni, J. Appl. Phys. 83 (1998) 3432.
    N_z = (1/(pi)) * (... complicated expression in terms of a, b, c)

  For our 5x2x1 mm cuboid, we use the well-known approximate values.
  The exact Aharoni formula is implemented below.

  Run: wolframscript -file polarizability_tensor_3D.wls
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

Print["=== polarizability_tensor_3D.wls ==="];
Print["3D polarizability tensor with Aharoni demagnetization factors"];
Print["Cuboid: ", a 1000, " x ", b 1000, " x ", c 1000, " mm Cu, V = ", N[VCub, 6], " m^3"];
Print[""];

(* === Aharoni 1998 demagnetization factor formula ======================== *)
(* For a rectangular prism with sides 2A, 2B, 2C centred at origin,
   the demagnetization factor along z (the C axis) is given by
   Aharoni's closed form. Here we use sides (a, b, c) = (2A, 2B, 2C). *)

aharoniNz[ax_, bx_, cx_] := Module[
  {A = ax/2, B = bx/2, C = cx/2, R, Rab, Rac, Rbc},
  R   = Sqrt[A^2 + B^2 + C^2];
  Rab = Sqrt[A^2 + B^2];
  Rac = Sqrt[A^2 + C^2];
  Rbc = Sqrt[B^2 + C^2];
  (1/Pi) (
    (B^2 - C^2)/(2 B C) Log[(R - A)/(R + A)] +
    (A^2 - C^2)/(2 A C) Log[(R - B)/(R + B)] +
    B/(2 C) Log[(Rab + A)/(Rab - A)] +
    A/(2 C) Log[(Rab + B)/(Rab - B)] +
    C/(2 A) Log[(Rbc - B)/(Rbc + B)] +
    C/(2 B) Log[(Rac - A)/(Rac + A)] +
    2 ArcTan[A B / (C R)] +
    (A^3 + B^3 - 2 C^3)/(3 A B C) +
    (A^2 + B^2 - 2 C^2)/(3 A B C) R +
    C/(A B) (Rac + Rbc) -
    (Rab^3 + Rbc^3 + Rac^3)/(3 A B C)
  )];

(* By symmetry, N_x = aharoni(b, c, a), N_y = aharoni(a, c, b), N_z = aharoni(a, b, c) *)
Nx = aharoniNz[b, c, a];
Ny = aharoniNz[a, c, b];
Nz = aharoniNz[a, b, c];

Print["Aharoni demagnetization factors:"];
Print["  N_x (along 5 mm axis) = ", N[Nx, 6]];
Print["  N_y (along 2 mm axis) = ", N[Ny, 6]];
Print["  N_z (along 1 mm axis) = ", N[Nz, 6]];
Print["  Sum N_x+N_y+N_z       = ", N[Nx + Ny + Nz, 8], "  (target 1)"];
Print[""];

(* === Per-cuboid Foster building block (from §V.B) ====================== *)
Mmax = 13;
modeList = Flatten[Table[{m, n, p},
  {m, 1, Mmax, 2}, {n, 1, Mmax, 2}, {p, 1, Mmax, 2}], 2];

g2MNP[m_, n_, p_] := 512 a b c / ((m n p)^2 Pi^6);
lambdaMNP[m_, n_, p_] := (m Pi/a)^2 + (n Pi/b)^2 + (p Pi/c)^2;
tauMNP[m_, n_, p_]    := mu sigma / lambdaMNP[m, n, p];

g2Vec  = N[Table[g2MNP[mn[[1]], mn[[2]], mn[[3]]], {mn, modeList}], 30];
tauVec = N[Table[tauMNP[mn[[1]], mn[[2]], mn[[3]]], {mn, modeList}], 30];

YcubFoster[s_] := sigma Total[g2Vec / (1 + s tauVec)];
YcubSchur[s_]  := YcubFoster[s] + KsibcCub Sqrt[s] / (s + dParam);

(* === Intrinsic susceptibility chi(s) =================================== *)
chiFoster[s_] := YcubFoster[s] / (sigma VCub) - 1;
chiSchur[s_]  := YcubSchur[s] / (sigma VCub) - 1;

(* Verify limits *)
Print["chi(s) sanity (Foster basis):"];
Print["  chi(0+)        = ", N[Re[chiFoster[I 2 Pi 0.01]], 6], "  (target 0)"];
Print["  chi(infinity)  = ", N[Re[chiFoster[I 2 Pi 1*^14]], 6], "  (target -1)"];
Print[""];

(* === Polarizability tensor: alpha_i = V chi/(1 + N_i chi) per axis ===== *)
alphaTensorFoster[s_, Nfac_] := VCub chiFoster[s] / (1 + Nfac chiFoster[s]);
alphaTensorSchur[s_, Nfac_]  := VCub chiSchur[s]  / (1 + Nfac chiSchur[s]);

(* === PEC limits along each axis ======================================== *)
Print["alpha_i(s -> infinity) at PEC limit:"];
With[{sInf = I 2 Pi 1*^14},
  Print["  alpha_x(inf) = -V/(1-N_x) = ", N[-VCub/(1 - Nx), 8],
        "  Foster:  ", N[Re[alphaTensorFoster[sInf, Nx]], 8]];
  Print["  alpha_y(inf) = -V/(1-N_y) = ", N[-VCub/(1 - Ny), 8],
        "  Foster:  ", N[Re[alphaTensorFoster[sInf, Ny]], 8]];
  Print["  alpha_z(inf) = -V/(1-N_z) = ", N[-VCub/(1 - Nz), 8],
        "  Foster:  ", N[Re[alphaTensorFoster[sInf, Nz]], 8]];
];
Print[""];

(* === Wall-band response per axis ======================================== *)
Print["Polarizability tensor per-axis response across wall band:"];
Print[""];
Print[StringPadRight["f [Hz]", 12] <>
      StringPadRight["alpha_x", 20] <>
      StringPadRight["alpha_y", 20] <>
      "alpha_z"];
Print[StringRepeat["-", 75]];

Do[
  With[{omega = 2 Pi (10.0^logf), s = I 2 Pi (10.0^logf)},
    With[{ax = alphaTensorFoster[s, Nx],
          ay = alphaTensorFoster[s, Ny],
          az = alphaTensorFoster[s, Nz]},
      Print["f=1e", logf, "      ",
            StringPadRight[ToString[N[ax, 5]], 20],
            StringPadRight[ToString[N[ay, 5]], 20],
            N[az, 5]]
    ]],
  {logf, {2, 4, 5, 6, 7, 8, 10, 12}}];

(* === Cross-check: m_y at PEC vs NGSolve =============================== *)
Print[""];
Print["Cross-check vs NGSolve FEM (per-cuboid in B_y=1 T):"];
With[{B0 = 1.0, omega = 2 Pi 1*^7, s = I 2 Pi 1*^7},
  With[{alphay = alphaTensorFoster[s, Ny]},
    With[{mPredicted = alphay B0 / mu},
      Print["  f = 10^7 Hz, applied B_y = 1 T"];
      Print["  alpha_y(s) (Foster tensor)   = ", N[alphay, 6]];
      Print["  m_y_predicted = alpha_y * B / mu_0 = ", N[mPredicted, 6], " A m^2"];
      Print["  NGSolve B.2 (per-cuboid) gave    m_y ~ -7.17j A m^2"];
      Print["  Ratio predicted/NGSolve       = ", N[mPredicted/(-7.17 I), 4]]
    ]];
  Print["Note: NGSolve at f=10^7 is under-resolved (mesh maxh=0.5mm,"];
  Print["      skin depth=21um, c/delta=48); proper FEM cross-check"];
  Print["      requires skin-graded mesh.  See B.2 honest scope note."];
];

(* === Trace identity ====================================================== *)
Print[""];
Print["Trace identity (sum over principal axes):"];
With[{sProbe = I 2 Pi 1*^14},
  With[{trace = alphaTensorFoster[sProbe, Nx] +
                alphaTensorFoster[sProbe, Ny] +
                alphaTensorFoster[sProbe, Nz]},
    Print["  alpha_x + alpha_y + alpha_z (s=inf) = ", N[Re[trace], 6]];
    Print["  Conventional bound: -3 V = ", N[-3 VCub, 6]];
    Print["  Ratio (always <= -3V for diamag)    = ", N[Re[trace]/(-3 VCub), 6]]
  ]];

(* === CSV output ========================================================= *)
csvData = {{"f_Hz", "alpha_x_re", "alpha_x_im", "alpha_y_re", "alpha_y_im",
            "alpha_z_re", "alpha_z_im"}};
Do[
  With[{omega = 2 Pi (10.0^logf), s = I 2 Pi (10.0^logf)},
    With[{ax = alphaTensorFoster[s, Nx],
          ay = alphaTensorFoster[s, Ny],
          az = alphaTensorFoster[s, Nz]},
      AppendTo[csvData,
        {10.0^logf, N[Re[ax], 8], N[Im[ax], 8],
                    N[Re[ay], 8], N[Im[ay], 8],
                    N[Re[az], 8], N[Im[az], 8]}]
    ]],
  {logf, {2, 3, 4, 5, 6, 7, 8, 10, 12}}];

Export["polarizability_tensor_3D.csv", csvData, "CSV"];
Print[""];
Print["Wrote polarizability_tensor_3D.csv"];

Print[""];
Print["=== polarizability_tensor_3D.wls DONE ==="];
