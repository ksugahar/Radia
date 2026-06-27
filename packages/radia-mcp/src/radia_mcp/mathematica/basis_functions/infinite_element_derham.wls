#!/usr/bin/env wolframscript
(* ::Package:: *)

(*
  infinite_element_derham.wls -- the de Rham (exact-sequence) INFINITE ELEMENT on the
  spherical EXTERIOR r >= a (low-frequency / static).  Mathematica twin of
  docs/open_boundary/INFINITE_ELEMENT_SOTA.md and the maintained
  radia.infinite_element API,
  here using Mathematica's built-in ORTHONORMAL spherical Grad/Curl/Div as the authority.

      H1 --grad--> H(curl) --curl--> H(div) --div--> L2     (exterior, mode n)

  The static infinite element decays algebraically; the de Rham version (Demkowicz & Pal,
  "An infinite element for Maxwell's equations", CMAME 164, 1998) chooses the RADIAL decay
  families so that grad / curl / div COMMUTE (the commuting diagram):
      0-form  r^-k Y                       radial powers  S0 = {n+1 .. n+P}
      1-/2-form  P = (r^-k Y,  0,     0)   (radial)        S1 = S0 + 1
                 B = (0, r^-k dY,     0)   (poloidal-tan)  S2 = S0 + 2
                 C = (0,     0, r^-k dY)   (toroidal)      S3 = S0 + 3
      3-form  r^-k Y
  i.e. the radial power shifts by +1 per form degree.  Verified EXACTLY (FullSimplify):
      grad(r^-k Y) = -k P(k+1) + 1 B(k+1)
      curl P(k)    = -1 C(k+1)
      curl B(k)    = (1-k) C(k+1)
      curl C(k)    = -n(n+1) P(k+1) + (k-1) B(k+1)        [Legendre eigenvalue]
      div  P(k)    = (2-k) r^-(k+1) Y
      div  B(k)    = -n(n+1) r^-(k+1) Y                   [Legendre eigenvalue]
      div  C(k)    = 0
  plus the complex property  curl o grad = 0,  div o curl = 0.

  The toroidal-curl and div-B closures rest on the angular-Laplace (Legendre) eigenvalue
      (1/sin) d/dth( sin dY/dth ) = -n(n+1) Y,
  checked here for general m too (so the radial families are m-independent).

  m=0 axisymmetric representative (Y = LegendreP[n, Cos[th]]).  Known construction
  (Demkowicz-Pal); implemented + self-tested, no novelty.  The 0-form tower is the scalar
  Bettess infinite element (act7_25).  NOTE: the SHIPPED de Rham open boundary is instead
  the coordinate-mapping family -- the Kelvin transformation (and the coordinate-scaling
  infinite-element domain in commercial FE) -- standard Nedelec / RT elements on the
  mapped region, de Rham inherited for free; this decay-basis exact-sequence IE is the
  academic alternative.  The HIGH-FREQUENCY (radiating) infinite element needs an
  oscillatory exp(i k r) radial basis (Astley conjugated / wave-envelope; Demkowicz-Pal
  for Maxwell) -- a separate construction, not this static one.

  Reference: Demkowicz & Pal, CMAME 164 (1998); Bettess, "Infinite Elements" (1992);
  Silvester-Hsieh (1971).  See ../README.md and the DtN-spectrum datasheet
  docs/open_boundary/DTN_SPECTRUM_COMPARISON.md.
*)

(* orthonormal spherical operators; chart {r, th, ph} is positional: radius, polar, azimuth *)
grad[f_] := Grad[f, {r, th, ph}, "Spherical"];
crl[F_]  := Curl[F, {r, th, ph}, "Spherical"];
dvg[F_]  := Div[F, {r, th, ph}, "Spherical"];

Yf[n_]      := LegendreP[n, Cos[th]];
sc0[k_, n_] := r^(-k) Yf[n];
Pf[k_, n_]  := {r^(-k) Yf[n], 0, 0};
Bf[k_, n_]  := {0, r^(-k) D[Yf[n], th], 0};
Cf[k_, n_]  := {0, 0, r^(-k) D[Yf[n], th]};

zS[s_] := TrueQ[PossibleZeroQ[FullSimplify[s, 0 < th < Pi]]];   (* scalar is zero?  *)
zV[v_] := AllTrue[v, zS];                                       (* vector is zero?  *)

(* ============================ self-test ============================ *)
Module[{ok = True, chk, P = 3, S0},
  chk[name_, cond_] := (Print["  ", If[cond, "ok  ", "FAIL"], name]; ok = ok && cond);

  Print["infinite_element_derham : exterior de Rham (exact-sequence) IE, mode n, P=3"];

  Do[
    (* (0) Legendre angular-Laplace eigenvalue (m=0) -- the closure identity *)
    chk["n=" <> ToString[n] <> ": Legendre eig (1/sin)d(sin dY) = -n(n+1) Y",
      zS[(1/Sin[th]) D[Sin[th] D[Yf[n], th], th] + n (n + 1) Yf[n]]];

    (* (1) complex property  d o d = 0 *)
    chk["n=" <> ToString[n] <> ": curl(grad u) = 0   (complex d o d = 0)",
      zV[crl[grad[sc0[n + 1, n]]]]];
    chk["n=" <> ToString[n] <> ": div(curl A) = 0   (complex d o d = 0)",
      zS[dvg[crl[Cf[n + 1, n]]]]];

    S0 = Range[n + 1, n + P];
    (* (2) grad : V0 -> V1   (radial power +1) *)
    chk["n=" <> ToString[n] <> ": grad(V0) subset V1 EXACT  (-k P + B; S0->S1)",
      AllTrue[S0, zV[grad[sc0[#, n]] - (-# Pf[# + 1, n] + Bf[# + 1, n])] &]];
    (* (3) curl : V1 -> V2 *)
    chk["n=" <> ToString[n] <> ": curl(V1) subset V2 EXACT (toroidal via Legendre -n(n+1))",
      AllTrue[S0, (zV[crl[Pf[#, n]] - (-1) Cf[# + 1, n]] &&
                   zV[crl[Bf[#, n]] - (1 - #) Cf[# + 1, n]] &&
                   zV[crl[Cf[#, n]] - (-n (n + 1) Pf[# + 1, n] + (# - 1) Bf[# + 1, n])]) &]];
    (* (4) div : V2 -> V3 *)
    chk["n=" <> ToString[n] <> ": div(V2) subset V3 EXACT  (Legendre -n(n+1) in div B)",
      AllTrue[S0, (zS[dvg[Pf[#, n]] - (2 - #) sc0[# + 1, n]] &&
                   zS[dvg[Bf[#, n]] - (-n (n + 1)) sc0[# + 1, n]] &&
                   zS[dvg[Cf[#, n]]]) &]],
    {n, 1, 3}];

  (* (5) general m: the closure eigenvalue is m-independent -> so are the radial families *)
  chk["general m: spherical-Laplace angular eig = -n(n+1) for SphericalHarmonicY (m != 0)",
    AllTrue[{{2, 1}, {2, 2}, {3, 1}, {3, -2}},
      Function[nm,
        Module[{nn = nm[[1]], mm = nm[[2]], Y},
          Y = FunctionExpand[SphericalHarmonicY[nn, mm, th, ph]];
          zS[(1/Sin[th]) D[Sin[th] D[Y, th], th] + (1/Sin[th]^2) D[Y, {ph, 2}] +
             nn (nn + 1) Y]]]]];

  Print[""];
  If[ok, Print["ALL CHECKS PASSED"], Print["SOME CHECKS FAILED"]; Exit[1]];
];
