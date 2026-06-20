#!/usr/bin/env wolframscript
(* ::Package:: *)

(*
  cuboid3D_quadrupole_3D.wls

  Phase 3 Task A: genuine 3D quadrupole verification of Theorem 3
  (multipole-resolved SIBC asymptote) on a 5x2x1 mm Cu cuboid.

  Extends cuboid3D_quadrupole.wls (which used a z-independent
  quadrupole on a 2D long cuboid) to a fully 3D quadrupole:
     u_ext = x y      (3D tesseral quadrupole, harmonic on R^3
                       since Laplacian(xy) = 0)
  on the rectangular box Omega = [0, a] x [0, b] x [0, c].
  This source has cross-axis structure and tests Theorem 3 in
  a fully 3D setting.

  Theorem 3 (V.E):
    K_SIBC^(u_ext) = sqrt(sigma/mu) * ||u_ext||^2_{L^2(d Omega)}
  with the boundary L^2 norm taken over all 6 cuboid faces.

  For u_ext = xy on [0,a] x [0,b] x [0,c]:
    Face x=0 (u_ext=0):       contribution 0
    Face x=a (u_ext=a*y):     a^2 * c * b^3/3
    Face y=0 (u_ext=0):       contribution 0
    Face y=b (u_ext=b*x):     b^2 * c * a^3/3
    Face z=0 (u_ext=xy):      a^3/3 * b^3/3
    Face z=c (u_ext=xy):      a^3/3 * b^3/3
  ||u_ext||^2_{L^2(boundary)} = c*(a^2 b^3 + a^3 b^2)/3 + 2 a^3 b^3 / 9
                            = a^2 b^2 c (a+b)/3 + 2 a^3 b^3 / 9

  Foster eigenmode sum on the 3D cuboid:
    phi_mnp(x,y,z) = (sqrt(8/(abc))) sin(m pi x/a) sin(n pi y/b) sin(p pi z/c)
    lambda_mnp = (m pi/a)^2 + (n pi/b)^2 + (p pi/c)^2
    tau_mnp    = mu sigma / lambda_mnp
    g^2_mnp    = sigma * <u_ext, phi_mnp>^2
    Y(s) = Sum_{mnp} g^2_mnp / (1 + s tau_mnp)
    Y(0) = sigma * ||u_ext||^2_{L^2(Omega)}     (Parseval)
    K_SIBC^(quad,3D) = sqrt(sigma/mu) * ||u_ext||^2_{L^2(boundary Omega)}

  For u_ext = xy, projections onto Dirichlet eigenmodes:
    <xy, sin(m pi x/a)> = intX[m] * intSin[n,b] * intSin[p,c]  (x part)
    Hmm not quite - need <xy, phi_mnp> = <x, sin_m^x> * <y, sin_n^y> * <1, sin_p^z>
    All separable.

  Run: wolframscript -file cuboid3D_quadrupole_3D.wls
*)

ClearAll[a, b, c, sigma, mu];

a = 5*^-3; b = 2*^-3; c = 1*^-3;
sigma = 5.80*^7;
mu    = 4 Pi 1*^-7;

Print["=== cuboid3D_quadrupole_3D.wls ==="];
Print["3D cuboid ", a 1000, " x ", b 1000, " x ", c 1000, " mm Cu"];
Print["3D quadrupole excitation: u_ext = x y (harmonic)"];
Print["sqrt(sigma/mu) = ", N[Sqrt[sigma/mu], 8]];
Print[""];

(* === Projection integrals =============================================== *)
(* <x, sin(m pi x/a)> on [0,a]:  intX[m, a] *)
intX[m_, len_] := If[OddQ[m],
   len^2 / (m Pi) (1 - 2/(m Pi) (-1)^((m-1)/2)),  (* odd m *)
   -len^2 / (m Pi) (-1)^(m/2)];                    (* even m, may be needed *)

(* Hmm let me recompute. Integral of x sin(m pi x/a) from 0 to a:
   Integration by parts: u = x, dv = sin(m pi x/a) dx
   = -[a x / (m pi)] cos(m pi x/a) | + (a/(m pi)) integral cos dx
   = -(a^2/(m pi)) cos(m pi) + 0 + (a/(m pi))(a/(m pi)) [sin(m pi x/a)]
   = -(a^2/(m pi)) cos(m pi) + (a^2/(m pi)^2)(sin(m pi) - 0)
   = -(a^2/(m pi)) (-1)^m + 0
   = (a^2/(m pi)) (-1)^(m+1) *)
intX[m_, len_] := (len^2 / (m Pi)) (-1)^(m+1);

(* <1, sin(m pi x/a)> = 2a/(m pi) for m odd, 0 for m even *)
intSin[m_, len_] := If[OddQ[m], 2 len / (m Pi), 0];

(* Projection: <u_ext = x*y, phi_mnp>
   = sqrt(8/(abc)) * <x, sin(m pi x/a)> * <y, sin(n pi y/b)> * <1, sin(p pi z/c)>
   = sqrt(8/(abc)) * intX[m, a] * intX[n, b] * intSin[p, c]                  *)
gMNP[m_, n_, p_] := Sqrt[8 / (a b c)] intX[m, a] intX[n, b] intSin[p, c];

lambdaMNP[m_, n_, p_] := (m Pi/a)^2 + (n Pi/b)^2 + (p Pi/c)^2;
tauMNP[m_, n_, p_]    := mu sigma / lambdaMNP[m, n, p];

(* === Mode truncation ====================================================
   For u_ext = xy:
     - m: needs ALL m (intX nonzero for all m), but symmetric in sign so all m count
     - n: same, ALL n contribute
     - p: only odd p contributes (intSin = 0 for p even)                  *)
Mmax = 99;
modeList = Flatten[
  Table[{m, n, p}, {m, 1, Mmax}, {n, 1, Mmax}, {p, 1, Mmax, 2}], 2];
Nmodes = Length[modeList];

Print["Foster modes (m,n all, p odd): ", Nmodes, " (Mmax = ", Mmax, ")"];

gVec  = N[Table[gMNP[mn[[1]], mn[[2]], mn[[3]]], {mn, modeList}], 30];
tauV  = N[Table[tauMNP[mn[[1]], mn[[2]], mn[[3]]], {mn, modeList}], 30];

(* Y(0) via Foster sum *)
YDC = sigma Total[gVec^2];

(* Analytical ||u_ext||^2_L^2(Omega) for u_ext = xy:
   = Integral of (xy)^2 dx dy dz = (a^3/3)(b^3/3)(c) = a^3 b^3 c / 9         *)
analVolL2 = sigma a^3 b^3 c / 9;

Print[""];
Print["Y(0) Foster sum               = ", N[YDC, 10]];
Print["sigma * ||xy||^2_L^2(Omega)   = ", N[analVolL2, 10]];
Print["Foster convergence (vol)      = ", N[100 YDC/analVolL2, 6], " %"];
Print[""];

(* === K_SIBC^(quad 3D xy) prediction from boundary L^2 =================== *)
(* From analytical derivation in header:
     ||xy||^2_L^2(boundary) = (a^2 b^2 c (a+b))/3 + 2 a^3 b^3 / 9          *)
bndryL2Analytical = (a^2 b^2 c (a + b))/3 + 2 a^3 b^3 / 9;

(* Cross-check via NIntegrate over all 6 faces *)
bndryL2Numerical =
  NIntegrate[(0)^2,    {y, 0, b}, {z, 0, c}] +                  (* x = 0 *)
  NIntegrate[(a y)^2,  {y, 0, b}, {z, 0, c}] +                  (* x = a *)
  NIntegrate[(0)^2,    {x, 0, a}, {z, 0, c}] +                  (* y = 0 *)
  NIntegrate[(b x)^2,  {x, 0, a}, {z, 0, c}] +                  (* y = b *)
  NIntegrate[(x y)^2,  {x, 0, a}, {y, 0, b}] +                  (* z = 0 *)
  NIntegrate[(x y)^2,  {x, 0, a}, {y, 0, b}];                   (* z = c *)

Print["||xy||^2_L^2(boundary) analytical    = ", N[bndryL2Analytical, 8]];
Print["||xy||^2_L^2(boundary) numerical     = ", N[bndryL2Numerical, 8]];
Print["match: ", N[bndryL2Analytical / bndryL2Numerical, 10]];

KsibcQuad3D = Sqrt[sigma/mu] bndryL2Analytical;
Print[""];
Print["Theorem 3 prediction:"];
Print["  K_SIBC^(quad 3D xy) = sqrt(sigma/mu) * ||xy||^2_L^2(boundary)"];
Print["                     = ", N[KsibcQuad3D, 8]];
Print[""];

(* === Foster sum verification at high omega ============================== *)
YfosterFn[s_] := sigma Total[gVec^2 / (1 + s tauV)];

Print["SIBC asymptote verification: |Y(j omega)| sqrt(omega) -> K_SIBC at high omega"];
Print[""];
Print["f [Hz]    |Y(j omega)|    |Y| sqrt(omega)    ratio r"];
Print[StringRepeat["-", 60]];

Do[
  With[{omega = 2 Pi (10.0^logf)},
    With[{s = I omega},
      With[{ymag = Abs[N[YfosterFn[s], 20]]},
        Print["f=1e", logf, "   ",
              "  |Y| = ", N[ymag, 5],
              "  K_obs = ", N[ymag Sqrt[omega], 5],
              "  r = ", N[ymag Sqrt[omega] / KsibcQuad3D, 5]]
      ]]],
  {logf, {3, 4, 5, 6, 7, 8}}];

(* === Summary table ======================================================= *)
Print[""];
Print["Summary (Theorem 3 numerical verification, 3D quadrupole xy):"];
Print[""];
Print["  Quantity                  Predicted   Foster-derived   Ratio"];
Print["  -----------------------  ----------   --------------   -----"];
Print["  ||u||^2_L^2(Omega)        ",
      N[analVolL2/sigma, 5], "        ", N[YDC/sigma, 5],
      "      ", N[YDC/analVolL2, 6]];

(* Estimate K_SIBC from Foster sum at moderate frequency *)
fProbe = 1*^7;
omegaProbe = 2 Pi fProbe;
KsibcFosterObs = Abs[N[YfosterFn[I omegaProbe], 20]] Sqrt[omegaProbe];
Print["  K_SIBC^(quad 3D)         ",
      N[KsibcQuad3D, 5], "    ", N[KsibcFosterObs, 5],
      "    ", N[KsibcFosterObs/KsibcQuad3D, 5]];

(* === Export CSV ========================================================= *)
csvData = {{"f_Hz", "Y_mag", "K_obs", "ratio_r"}};
Do[
  With[{omega = 2 Pi (10.0^logf), s = I 2 Pi (10.0^logf)},
    With[{ymag = Abs[N[YfosterFn[s], 20]]},
      AppendTo[csvData,
        {10.0^logf, N[ymag, 8], N[ymag Sqrt[omega], 8], N[ymag Sqrt[omega]/KsibcQuad3D, 8]}]
    ]],
  {logf, {2, 3, 4, 5, 6, 7, 8, 9}}];

Export["cuboid3D_quadrupole_3D.csv", csvData, "CSV"];
Print[""];
Print["Wrote cuboid3D_quadrupole_3D.csv"];

Print[""];
Print["=== cuboid3D_quadrupole_3D.wls DONE ==="];
