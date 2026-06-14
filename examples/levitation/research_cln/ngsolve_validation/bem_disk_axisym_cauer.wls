#!/usr/bin/env wolframscript
(*
  bem_disk_axisym_cauer.wls — extends bem_disk_axisym_v3_refined.wls with
  Foster amplitudes + alpha_n moments, dumps everything for Python-side
  Cauer-I CFE -> R_n, L_n extraction.

  Output: bem_disk_axisym_cauer.json
    method, params              -- same as v3
    bem_tau_us                  -- top N tau values (us)
    eigvals                     -- raw K_sym eigenvalues (positive subset)
    foster_amplitudes_g         -- g_k = <v_k, sqrt(dV) A_ext> per mode
    alphas                      -- moments alpha_0 .. alpha_{N_MOMENTS-1}
    element_r, element_z, dV    -- element geometry (for Python-side checks)

  Convention:
    J_phi(r, z) = ring current density at element centre (r_i, z_i)
    K_sym[i,j] = mu_0 G_axi(r_i, z_i; r_j, z_j) sqrt(dV_i dV_j)   (=v3 def)
    Eigenvalue lambda_k of K_sym <=> tau_k = sigma * lambda_k  (=v3 def)

    External uniform B_z drive: A_ext(r) = B_0 r / 2.
    Foster amplitude g_k = sum_i v_k[i] * sqrt(dV_i) * A_ext(r_i)
    Moments alpha_n = (-1)^n * Sum_k g_k^2 * tau_k^(n+1)
                                        (sigma^(n+1) mu_0^? ... factor folded
                                         into definition; we extract Cauer
                                         RATIOS so global factors cancel.)
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

NMODES = 50;          (* number of Foster modes to keep *)
NMOMENTS = 20;        (* number of alpha_n moments to dump *)

Print["=== BEM disk axisym for Cauer-I extraction ==="];
Print["  R=", R*1000, "mm, t=", t*1000, "mm, NR=", NR, " Nz=", Nz, " ne=", ne];
Print["  modes=", NMODES, " moments=", NMOMENTS];

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

Print["  Building K_sym (", ne, "x", ne, ")..."];
t0 = AbsoluteTime[];

Ksym = Table[
  Module[{r1, z1, r2, z2, gval},
    {r1, z1} = elements[[i]];
    {r2, z2} = elements[[j]];
    If[i === j,
      gval = GaxiSelf[r1],
      gval = Gaxi[r1, z1, r2, z2]
    ];
    mu0 * gval * Sqrt[dV[[i]] * dV[[j]]]
  ],
  {i, 1, ne}, {j, 1, ne}
];

Print["    K_sym built in ", AbsoluteTime[] - t0, " s"];

Print["  Solving symmetric eigenvalue problem (top ", NMODES, ")..."];
t0 = AbsoluteTime[];
{eigvals, eigvecs} = Eigensystem[N[Ksym], NMODES];
Print["    Done in ", AbsoluteTime[] - t0, " s"];

(* Eigensystem returns eigenvalues in DECREASING magnitude.
   eigvecs[[k]] is the k-th eigenvector (row vector). *)

(* Filter to positive eigenvalues (eddy decay modes) *)
posMask = Table[Re[eigvals[[k]]] > 0, {k, 1, NMODES}];
eigvalsPos = Pick[eigvals, posMask];
eigvecsPos = Pick[eigvecs, posMask];

Print["  Positive-eigenvalue modes kept: ", Length[eigvalsPos], "/", NMODES];

(* Sort by tau decreasing (= eigenvalue decreasing) *)
ord = Reverse[Ordering[Re[eigvalsPos]]];
eigvalsSorted = eigvalsPos[[ord]];
eigvecsSorted = eigvecsPos[[ord]];

(* tau_k = sigma * lambda_k *)
tausUs = N[sigma * Re[eigvalsSorted] * 10^6, 16];

Print["  tau_1..tau_6 us = ", N[tausUs[[1;;6]], 6]];

(* === Foster amplitudes ===
   v_k is k-th eigenvector, treated as a column. Components v_k[i] are in
   the *symmetric* basis (= J_i sqrt(dV_i)). External field driving is
   A_ext(r_i) = B_0 r_i / 2; in symmetric basis the projection coefficient
   on element i is sqrt(dV_i) * A_ext(r_i).

   g_k = (v_k . a_ext_sym) where a_ext_sym[i] = sqrt(dV_i) * (B_0 r_i / 2). *)

aExtSym = Table[Sqrt[dV[[i]]] * B0 * elements[[i, 1]] / 2, {i, 1, ne}];

Print["  Computing Foster amplitudes g_k..."];
gAmpl = Table[
  Module[{vk = eigvecsSorted[[k]]},
    Re[vk . aExtSym]
  ],
  {k, 1, Length[eigvalsSorted]}];

Print["  g_1..g_6 = ", N[gAmpl[[1;;Min[6, Length[gAmpl]]]], 6]];

(* === Moments alpha_n =================================
   alpha_n = (-1)^n * Sum_k g_k^2 * tau_k^(n+1)         (n = 0, 1, 2, ...)
   (Foster-driven impedance Taylor expansion around s = 0.)

   Note: signs and normalizations are chosen so that the Cauer-I CFE on
   c_taylor = [alpha_0, alpha_1, ...] reproduces L_n, R_n in the same
   convention as the Hiruma 3-term recursion in radia-axifemm.
*)

tausSec = sigma * Re[eigvalsSorted];   (* tau_k in seconds *)

alphasMomLegacy = Table[
  Module[{n = nn},
    Sum[gAmpl[[k]]^2 * tausSec[[k]]^(n + 1) * (-1)^n,
        {k, 1, Length[eigvalsSorted]}]
  ],
  {nn, 0, NMOMENTS - 1}];

Print["  alpha_0..alpha_5 = ", N[alphasMomLegacy[[1;;6]], 6]];

(* Save JSON =========================================== *)
result = <|
  "method" -> "BEM disk axisym + Foster amplitudes + alpha_n moments for Cauer-I",
  "params" -> <|"R_mm" -> R*1000, "t_mm" -> t*1000, "sigma" -> sigma,
                "B0" -> B0, "NR" -> NR, "Nz" -> Nz, "ne" -> ne,
                "n_modes" -> Length[eigvalsSorted], "n_moments" -> NMOMENTS|>,
  "bem_tau_us" -> tausUs,
  "bem_eigvals_lambda" -> N[Re[eigvalsSorted], 16],
  "foster_amplitudes_g" -> N[gAmpl, 16],
  "alphas" -> N[alphasMomLegacy, 16],
  "element_r" -> Table[elements[[i, 1]], {i, 1, ne}],
  "element_z" -> Table[elements[[i, 2]], {i, 1, ne}],
  "dV" -> N[dV, 16]
|>;
Export["bem_disk_axisym_cauer.json", result, "JSON"];
Print[""];
Print["Saved: bem_disk_axisym_cauer.json (", Length[eigvalsSorted], " modes, ",
      NMOMENTS, " moments)"];
