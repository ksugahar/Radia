#!/usr/bin/env wolframscript
(*
  bem_disk_axisym_cauer_highprec.wls — high-precision variant of
  bem_disk_axisym_cauer.wls with NMOMENTS=40 and arbitrary precision
  arithmetic for the Eigensystem + Foster amplitudes.

  Track B-ii high-k diagnosis: original n_moments=20 + standard mpmath
  dps=50 produces a sign change in tau_pair at k=6 (cylinder Cu disk).
  This script doubles n_moments and uses N[..., 80] to push the high-k
  Cauer extraction further before round-off ends it.

  Output: bem_disk_axisym_cauer_highprec.json (same schema as
          bem_disk_axisym_cauer.json, just longer alphas array)
*)

mu0 = 4 Pi * 10^-7;
sigma = 5.8 * 10^7;
R = 10 * 10^-3;
t = 2 * 10^-3;
B0 = 1;

(* Same mesh as v3 *)
NR = 80;
Nz = 24;
dr = R / NR;
dz = t / Nz;
ne = NR * Nz;

NMODES = 80;          (* number of Foster modes (was 50) *)
NMOMENTS = 40;        (* number of alpha_n moments (was 20) *)
WP = 80;              (* working precision in digits *)

Print["=== BEM disk axisym for Cauer-I extraction (HIGH PRECISION) ==="];
Print["  R=", R*1000, "mm, t=", t*1000, "mm, NR=", NR, " Nz=", Nz, " ne=", ne];
Print["  modes=", NMODES, " moments=", NMOMENTS, " WP=", WP];

elements = Flatten[Table[{(i - 1/2) dr, (j - 1/2) dz - t/2},
                         {i, 1, NR}, {j, 1, Nz}], 1];
dV = Table[2 Pi elements[[i, 1]] dr dz, {i, 1, ne}];

Gaxi[r1_, z1_, r2_, z2_] := Module[{k2, K1, E1, eps = 1.0*^-15},
  k2 = 4 r1 r2 / ((r1 + r2)^2 + (z1 - z2)^2 + eps);
  K1 = EllipticK[k2];
  E1 = EllipticE[k2];
  (1/(2 Pi^2 Sqrt[r1 r2 + eps] Sqrt[k2])) * ((1 - k2/2) K1 - E1)
];

aEq = Sqrt[dr dz / Pi];
GaxiSelf[r_] := (1/(2 Pi^2 r)) * (Log[8 r / aEq] - 2);

Print["  Building K_sym (", ne, "x", ne, ") at WP=", WP, "..."];
t0 = AbsoluteTime[];

Ksym = Table[
  Module[{r1, z1, r2, z2, gval},
    {r1, z1} = elements[[i]];
    {r2, z2} = elements[[j]];
    If[i === j,
      gval = GaxiSelf[r1],
      gval = Gaxi[r1, z1, r2, z2]
    ];
    N[mu0 * gval * Sqrt[dV[[i]] * dV[[j]]], WP]
  ],
  {i, 1, ne}, {j, 1, ne}
];

Print["    K_sym built in ", AbsoluteTime[] - t0, " s"];

Print["  Solving symmetric eigenvalue problem (top ", NMODES, ", WP=", WP, ")..."];
t0 = AbsoluteTime[];
{eigvals, eigvecs} = Eigensystem[Ksym, NMODES];
Print["    Done in ", AbsoluteTime[] - t0, " s"];

posMask = Table[Re[eigvals[[k]]] > 0, {k, 1, NMODES}];
eigvalsPos = Pick[eigvals, posMask];
eigvecsPos = Pick[eigvecs, posMask];

Print["  Positive-eigenvalue modes kept: ", Length[eigvalsPos], "/", NMODES];

ord = Reverse[Ordering[Re[eigvalsPos]]];
eigvalsSorted = eigvalsPos[[ord]];
eigvecsSorted = eigvecsPos[[ord]];

tausUs = N[sigma * Re[eigvalsSorted] * 10^6, 30];

Print["  tau_1..tau_8 us = ", N[tausUs[[1;;8]], 8]];

aExtSym = Table[N[Sqrt[dV[[i]]] * B0 * elements[[i, 1]] / 2, WP], {i, 1, ne}];

Print["  Computing Foster amplitudes g_k..."];
gAmpl = Table[
  Module[{vk = eigvecsSorted[[k]]},
    Re[vk . aExtSym]
  ],
  {k, 1, Length[eigvalsSorted]}];

Print["  g_1..g_6 = ", N[gAmpl[[1;;Min[6, Length[gAmpl]]]], 6]];

tausSec = sigma * Re[eigvalsSorted];

Print["  Computing ", NMOMENTS, " alpha_n moments..."];
alphasMomLegacy = Table[
  Module[{n = nn},
    Sum[gAmpl[[k]]^2 * tausSec[[k]]^(n + 1) * (-1)^n,
        {k, 1, Length[eigvalsSorted]}]
  ],
  {nn, 0, NMOMENTS - 1}];

Print["  alpha_0..alpha_5 = ", N[alphasMomLegacy[[1;;6]], 6]];

result = <|
  "method" -> "BEM disk axisym + Foster amplitudes + alpha_n moments for Cauer-I (HIGHPREC)",
  "params" -> <|"R_mm" -> R*1000, "t_mm" -> t*1000, "sigma" -> sigma,
                "B0" -> B0, "NR" -> NR, "Nz" -> Nz, "ne" -> ne,
                "n_modes" -> Length[eigvalsSorted], "n_moments" -> NMOMENTS,
                "wp_digits" -> WP|>,
  "bem_tau_us" -> tausUs,
  "bem_eigvals_lambda" -> N[Re[eigvalsSorted], 30],
  "foster_amplitudes_g" -> N[gAmpl, 30],
  "alphas" -> N[alphasMomLegacy, 30],
  "element_r" -> Table[elements[[i, 1]], {i, 1, ne}],
  "element_z" -> Table[elements[[i, 2]], {i, 1, ne}],
  "dV" -> N[dV, 30]
|>;
Export["bem_disk_axisym_cauer_highprec.json", result, "JSON"];
Print[""];
Print["Saved: bem_disk_axisym_cauer_highprec.json (", Length[eigvalsSorted], " modes, ",
      NMOMENTS, " moments, WP=", WP, ")"];
