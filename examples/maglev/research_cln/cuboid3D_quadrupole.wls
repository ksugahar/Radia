#!/usr/bin/env wolframscript
(* ::Package:: *)

(*
  cuboid3D_quadrupole.wls

  Theorem 2 (multipole-resolved SIBC asymptote) verification for the
  journal manuscript section V.E.  Quadrupole excitation on a 2D long
  cuboid (a x b Cu cross-section, c -> inf), demonstrating that the
  Schur-F K_SIBC coefficient generalises from "surface area" (dipole)
  to "boundary L^2 norm of the applied scalar potential" (multipole):

    K_SIBC^(n) = sqrt(sigma/mu) * ||u_ext||^2_{L^2(boundary Omega)}

  where u_ext is the applied scalar source in the scalar diffusion
  equation -Laplacian u + s mu sigma u = -s mu sigma u_ext.

  Dipole:     u_ext = 1                     -> ||1||^2_L^2(d-Omega) = perimeter
  Quadrupole: u_ext = (y^2 - x^2)/2         -> non-trivial boundary integral

  Eigenmode expansion (Dirichlet on [0,a]x[0,b]):
    phi_mn = (2/sqrt(ab)) sin(m pi x/a) sin(n pi y/b),  m,n = 1,3,5,...
    lambda_mn = (m pi/a)^2 + (n pi/b)^2
    tau_mn = mu sigma / lambda_mn
    Y(s) = Sum g^2_mn / (1 + s tau_mn),  g^2_mn = sigma * |<u_ext, phi_mn>|^2

  At SIBC limit: |Y(j omega)| * sqrt(omega) -> K_SIBC.

  Run: wolframscript -file cuboid3D_quadrupole.wls
*)

ClearAll[a, b, sigma, mu];

a = 5*^-3;
b = 2*^-3;
sigma = 5.80*^7;
mu = 4 Pi 1*^-7;

Print["=== cuboid3D_quadrupole.wls ==="];
Print["Geometry: 2D long cuboid, ", a 1000, " x ", b 1000, " mm Cu cross-section"];
Print["sigma = ", sigma, " S/m, mu = ", N[mu, 6], " H/m"];
Print["sqrt(sigma/mu) = ", N[Sqrt[sigma/mu], 8]];
Print[""];

(* Analytical projection helpers *)
intX2[m_, len_] := If[OddQ[m],
  (len^3/(m Pi))(1 - 4/(m Pi)^2),
  -len^3/(m Pi)];

intSin[m_, len_] := If[OddQ[m], 2 len/(m Pi), 0];

(* Quadrupole source: u_ext_quad = (y^2 - x^2)/2 *)
(* <u_ext_quad, sin(m pi x/a) sin(n pi y/b)> on rectangle [0,a]x[0,b] *)
PrawQuad[m_, n_] := (1/2)(intX2[n, b] intSin[m, a] - intX2[m, a] intSin[n, b]);

(* Dipole source (reference, must match existing slab/cuboid scripts): u_ext_dip = 1 *)
PrawDip[m_, n_] := intSin[m, a] intSin[n, b];

(* L^2-orthonormal basis projection: <u_ext, phi_mn> = (2/sqrt(ab)) * P_raw *)
gQuadMN[m_, n_] := (2/Sqrt[a b]) PrawQuad[m, n];
gDipMN[m_, n_]  := (2/Sqrt[a b]) PrawDip[m, n];

lambdaMN[m_, n_] := (m Pi/a)^2 + (n Pi/b)^2;
tauMN[m_, n_]    := mu sigma / lambdaMN[m, n];

(* Truncate Foster sum.  Quadrupole couples to ALL (m,n) except (m even,
   n even); dipole couples only to (m odd, n odd).  We include all
   (m,n) here so both excitations get full Parseval convergence. *)
Mmax = 199;
modeList = Flatten[Table[{m, n}, {m, 1, Mmax}, {n, 1, Mmax}], 1];
Nmodes = Length[modeList];

projQuad = N[Table[gQuadMN[mn[[1]], mn[[2]]], {mn, modeList}], 30];
projDip  = N[Table[gDipMN[mn[[1]], mn[[2]]],  {mn, modeList}], 30];
lamVec   = N[Table[lambdaMN[mn[[1]], mn[[2]]], {mn, modeList}], 30];
tauVec   = N[mu sigma / lamVec, 30];

Print["Foster spectrum: ", Nmodes, " modes (m, n = 1..", Mmax, ")"];
Print[""];

(* Y(s) = Sum sigma * |proj|^2 / (1 + s tau) *)
YquadFoster[s_] := sigma * Total[projQuad^2 / (1 + s tauVec)];
YdipFoster[s_]  := sigma * Total[projDip^2  / (1 + s tauVec)];

(* DC values *)
YquadDC  = N[YquadFoster[0], 20];
YdipDC   = N[YdipFoster[0], 20];

(* Analytical ||u_ext||^2_L^2(Omega) (volume L^2 norm) for cross-check *)
volL2Quad = sigma NIntegrate[((y^2 - x^2)/2)^2, {x, 0, a}, {y, 0, b}];
volL2Dip  = sigma NIntegrate[1, {x, 0, a}, {y, 0, b}];   (* = sigma * a * b *)

Print["Y(0) Foster sum (quadrupole) = ", N[YquadDC, 8]];
Print["sigma * ||u_quad||^2_L^2(Omega) = ", N[volL2Quad, 8]];
Print["Foster convergence (quad)    = ",
      N[100 YquadDC/volL2Quad, 6], " %"];
Print[""];
Print["Y(0) Foster sum (dipole)     = ", N[YdipDC, 8]];
Print["sigma * ||u_dip||^2_L^2(Omega) = sigma * ab = ", N[volL2Dip, 8]];
Print["Foster convergence (dipole)  = ",
      N[100 YdipDC/volL2Dip, 6], " %"];
Print[""];

(* Theorem 2: K_SIBC^(n) = sqrt(sigma/mu) * ||u_ext||^2_L^2(boundary Omega) *)
(* For quadrupole u_ext = (y^2 - x^2)/2 on the four faces of [0,a]x[0,b]: *)
(*   y=b face (x in [0,a]):  u_ext^2 = ((b^2 - x^2)/2)^2  *)
(*   y=0 face (x in [0,a]):  u_ext^2 = (x^2/2)^2          *)
(*   x=a face (y in [0,b]):  u_ext^2 = ((y^2 - a^2)/2)^2  *)
(*   x=0 face (y in [0,b]):  u_ext^2 = (y^2/2)^2          *)
bndryL2Quad =
  NIntegrate[((b^2 - x^2)/2)^2, {x, 0, a}] +
  NIntegrate[(x^2/2)^2,         {x, 0, a}] +
  NIntegrate[((y^2 - a^2)/2)^2, {y, 0, b}] +
  NIntegrate[(y^2/2)^2,         {y, 0, b}];

(* Dipole: ||1||^2_L^2(boundary) = perimeter = 2(a+b) *)
bndryL2Dip = 2 (a + b);

KsibcQuad = Sqrt[sigma/mu] * bndryL2Quad;
KsibcDip  = Sqrt[sigma/mu] * bndryL2Dip;

Print["Theorem 2: K_SIBC^(n) = sqrt(sigma/mu) * ||u_ext||^2_L^2(boundary Omega)"];
Print[""];
Print["Dipole boundary L^2 = perimeter = 2(a+b) = ",
      N[bndryL2Dip 1000, 6], " mm"];
Print["Predicted K_SIBC^(dipole)     = ", N[KsibcDip, 8]];
Print[""];
Print["Quadrupole boundary L^2 = ", N[bndryL2Quad, 6], " (m^5)"];
Print["Predicted K_SIBC^(quadrupole) = ", N[KsibcQuad, 6]];
Print[""];

(* Verify the SIBC asymptote: |Y(j omega)| * sqrt(omega) -> K_SIBC *)
freqList = {1*^4, 1*^5, 1*^6, 1*^7, 1*^8};

Print["Asymptote verification: |Y(j omega)| * sqrt(omega) -> K_SIBC"];
Print[""];
Print[StringPadRight["f [Hz]", 10] <>
      StringPadRight["|Y_dip|*sqrt(omega)", 22] <>
      StringPadRight["ratio_dip", 12] <>
      StringPadRight["|Y_quad|*sqrt(omega)", 22] <>
      "ratio_quad"];
Print[StringRepeat["-", 80]];

resultsTable = {{"f_Hz",
                 "Y_dip_mag", "Ydip_times_sqrt_omega", "ratio_dip",
                 "Y_quad_mag", "Yquad_times_sqrt_omega", "ratio_quad"}};

Do[
  Module[{omega, yDipMag, yQuadMag, kdipObs, kquadObs, rDip, rQuad},
    omega = 2 Pi f;
    yDipMag  = Abs[N[YdipFoster[I omega], 20]];
    yQuadMag = Abs[N[YquadFoster[I omega], 20]];
    kdipObs  = yDipMag Sqrt[omega];
    kquadObs = yQuadMag Sqrt[omega];
    rDip  = kdipObs / KsibcDip;
    rQuad = kquadObs / KsibcQuad;
    Print["f = ", N[f, 3],
          "  K_dip_obs = ", N[kdipObs, 5],
          "  r_dip = ", N[rDip, 5],
          "  K_quad_obs = ", N[kquadObs, 5],
          "  r_quad = ", N[rQuad, 5]];
    AppendTo[resultsTable,
      {N[f, 6], N[yDipMag, 8], N[kdipObs, 8], N[rDip, 8],
                N[yQuadMag, 8], N[kquadObs, 8], N[rQuad, 8]}]
  ],
  {f, freqList}];

Print[""];
Print["Summary table (Theorem 2 verification on 5x2 mm Cu, 2D long cuboid):"];
Print[""];
Print["  Mode        | K_SIBC^(n) predicted | ratio r at f=1e6 Hz"];
Print["  -----------+----------------------+--------------------"];
Module[{rDip6, rQuad6, fProbe = 1*^6, omegaProbe},
  omegaProbe = 2 Pi fProbe;
  rDip6  = Abs[N[YdipFoster[I omegaProbe], 20]]  Sqrt[omegaProbe] / KsibcDip;
  rQuad6 = Abs[N[YquadFoster[I omegaProbe], 20]] Sqrt[omegaProbe] / KsibcQuad;
  Print["  Dipole      | ", N[KsibcDip, 6], "           | ", N[rDip6, 5]];
  Print["  Quadrupole  | ", N[KsibcQuad, 6], "        | ", N[rQuad6, 5]];
];

(* Export CSV *)
Export["cuboid3D_quadrupole_KSIBC.csv", resultsTable, "CSV"];
Print[""];
Print["Wrote cuboid3D_quadrupole_KSIBC.csv"];
Print[""];
Print["=== cuboid3D_quadrupole.wls DONE ==="];
