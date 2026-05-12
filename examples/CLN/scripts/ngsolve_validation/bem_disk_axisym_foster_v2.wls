#!/usr/bin/env wolframscript
(*
  Axisymmetric BEM Foster-CLN for Cu disk (v2, corrected).

  Eigenvalue problem (eddy current decay in conductor):
    σ × ∫ K_axi(P, P') × J_φ(P') × dV' = (1/τ) × J_φ(P)

  Discretization: piecewise constant J_φ on rectangular cells (NR × Nz).
  Each cell i has center (r_i, z_i) and volume dV_i = 2π r_i × dr × dz.

  Symmetric matrix form: K_sym[i,j] = μ_0 × G_axi(i,j) × √(dV_i × dV_j)
  Eigenvalue: u = (σ/τ) × K_sym × u → τ_n = σ × eig(K_sym).

  Validation: ∞-cylinder leading τ_1 ≈ μ_0 σ a²/14.68 = 497 μs (Cu, a=10mm).
*)

mu0 = 4 Pi * 10^-7;
sigma = 5.8 * 10^7;
R = 10 * 10^-3;
t = 2 * 10^-3;

NR = 30;
Nz = 6;
dr = R / NR;
dz = t / Nz;
ne = NR * Nz;

Print["=== Axisymmetric BEM Foster v2 ==="];
Print["  R=", R*1000, "mm, t=", t*1000, "mm, NR=", NR, " Nz=", Nz, " ne=", ne];

(* Element centroids *)
elements = Flatten[Table[{(i - 0.5) * dr, (j - 0.5) * dz - t/2},
                         {i, 1, NR}, {j, 1, Nz}], 1];

(* dV[i] = 2 pi r_i dr dz *)
dV = Table[2 Pi elements[[i, 1]] dr dz, {i, 1, ne}];

(* Axisymmetric kernel for A_phi from current loop:
   A_phi(r, z; r', z') = (mu_0 / pi) * sqrt(r' / r) / k * [(1 - k^2/2) K - E]
   For volume integration: G_axi(r, z; r', z') is the kernel such that
     A_phi(r, z) = mu_0 × ∫ G_axi(r, z; r', z') × J_phi(r', z') × dV'
   With dV' = 2 pi r' dr' dz', G_axi has units [m^-1] (so mu_0 × G_axi × dV'
   gives correct A_phi units).

   From Jackson: A_phi from a circular loop of current I at (r', z'):
     A_phi(r, z) = (mu_0 I / pi) × sqrt(r'/r) / k × [(1 - k^2/2) K - E]
   For volume J_phi (per volume): I_loop = J_phi × dr' × dz'
   So A_phi(r, z) = ∫ (mu_0 / pi) × sqrt(r'/r) / k × [(1 - k^2/2) K - E]
                       × J_phi(r', z') × dr' dz'
   Comparing with A_phi = mu_0 × ∫ G_axi × J × dV' = mu_0 × ∫ G_axi × J × 2π r' dr' dz':
     G_axi(r, z; r', z') = (1/(2 pi^2 r')) × sqrt(r'/r) / k × [(1 - k^2/2) K - E]
                        = (1/(2 pi^2)) × 1/sqrt(r r') / k × [(1 - k^2/2) K - E]

   This is symmetric: G_axi(r, z; r', z') = G_axi(r', z'; r, z). Good. *)

Gaxi[r1_, z1_, r2_, z2_] := Module[{k2, K1, E1, eps = 1.0*^-15},
  k2 = 4 r1 r2 / ((r1 + r2)^2 + (z1 - z2)^2 + eps);
  K1 = EllipticK[k2];
  E1 = EllipticE[k2];
  (1/(2 Pi^2 Sqrt[r1 r2 + eps] Sqrt[k2])) * ((1 - k2/2) K1 - E1)
];

(* Self-term: integrate A_phi over a small cell of dimensions dr × dz centered
   at the source. The dominant log-singular part of the kernel can be
   integrated analytically. For square cell, approximation:
     G_self ≈ (1/(2 pi^2 r)) × [ln(8 r / a_eq) - 2]
   where a_eq = sqrt(dr × dz / pi) (equivalent disk radius).
   Per volume normalization: G_self × dV_self = (1/(2 pi^2 r)) × [ln - 2] × 2π r dr dz
                                              = (dr dz / pi) × [ln(8r/a_eq) - 2] *)

aEq = Sqrt[dr dz / Pi];
GaxiSelf[r_] := (1/(2 Pi^2 r)) * (Log[8 r / aEq] - 2);

Print["  Building K_sym (symmetric kernel matrix, ", ne, "×", ne, ")..."];
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

Print["  Solving symmetric eigenvalue problem..."];
t0 = AbsoluteTime[];
{eigvals, eigvecs} = Eigensystem[N[Ksym], 8];
Print["    Done in ", AbsoluteTime[] - t0, " s"];

(* τ_n = σ × eig (largest first) *)
positiveEigs = Select[Re[eigvals], # > 0 &];
tausUs = N[sigma * positiveEigs * 10^6, 8];
tausUsSorted = Reverse[Sort[tausUs]];

Print[""];
Print["=== Foster pole spectrum (Cu disk R=", R*1000, "mm, t=", t*1000, "mm) ==="];
Do[
  If[n <= Length[tausUsSorted],
    Print["  tau_", n, " = ", NumberForm[tausUsSorted[[n]], {6, 4}], " us"]
  ],
  {n, 1, 8}
];

(* Reference: ∞-cylinder analytical (Bessel J_1 zeros) *)
j1zeros = {3.83171, 7.01559, 10.17347, 13.32369, 16.47063};
Print[""];
Print["=== Reference: ∞-cylinder PEC analytical (Bessel J_1 zeros) ==="];
Do[
  Module[{tauTrue = mu0 sigma R^2 / j1zeros[[n]]^2 * 10^6},
    Print["  tau_", n, " = ", NumberForm[tauTrue, {6, 4}], " us"]
  ],
  {n, 1, 5}
];

Print[""];
Print["=== Comparison ==="];
v22bTauB2 = {208.32, 77.90, 46.41, 33.19, 25.22, 18.23};
Print[StringJoin[
  StringPadRight["n", 4],
  StringPadRight["BEM tau us", 14],
  StringPadRight["v22b tau_B2 us", 16],
  StringPadRight["inf-cyl Bessel us", 18]
]];
Do[
  Module[{bem, v22, infcyl},
    bem = If[n <= Length[tausUsSorted], tausUsSorted[[n]], 0];
    v22 = If[n <= Length[v22bTauB2], v22bTauB2[[n]], 0];
    infcyl = If[n <= 5, mu0 sigma R^2 / j1zeros[[n]]^2 * 10^6, 0];
    Print[StringJoin[
      StringPadRight[ToString[n], 4],
      StringPadRight[ToString[NumberForm[bem, {6, 3}]], 14],
      StringPadRight[ToString[NumberForm[v22, {6, 3}]], 16],
      StringPadRight[ToString[NumberForm[infcyl, {6, 3}]], 18]
    ]]
  ],
  {n, 1, 6}
];

result = <|
  "method" -> "Axisymmetric BEM Foster v2 (Cu disk)",
  "params" -> <|"R_mm" -> R*1000, "t_mm" -> t*1000, "sigma" -> sigma,
                "NR" -> NR, "Nz" -> Nz, "ne" -> ne|>,
  "tau_us" -> tausUsSorted,
  "inf_cylinder_bessel_tau_us" -> Table[mu0 sigma R^2 / j1zeros[[n]]^2 * 10^6, {n, 1, 5}],
  "v22b_tau_B2_for_compare" -> v22bTauB2
|>;
Export["bem_disk_axisym_foster_v2.json", result, "JSON"];
Print[""];
Print["Saved: bem_disk_axisym_foster_v2.json"];
