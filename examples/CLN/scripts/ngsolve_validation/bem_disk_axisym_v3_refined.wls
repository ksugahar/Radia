#!/usr/bin/env wolframscript
(*
  Axisymmetric BEM Foster-CLN for Cu disk (v3, refined + analytical k_z).

  Improvements over v2:
    - Refined mesh (NR=80, Nz=24) for ne=1920
    - Higher precision in self-term (analytic ring-stack inductance)
    - Reference analytical: PEC disk with both k_r (Bessel) and k_z (axial)
        τ_{m,n} = μ_0 σ a² / [j_{1,m}² + (n π a/t)²]
*)

mu0 = 4 Pi * 10^-7;
sigma = 5.8 * 10^7;
R = 10 * 10^-3;
t = 2 * 10^-3;

NR = 80;
Nz = 24;
dr = R / NR;
dz = t / Nz;
ne = NR * Nz;

Print["=== Axisymmetric BEM Foster v3 (refined) ==="];
Print["  R=", R*1000, "mm, t=", t*1000, "mm, NR=", NR, " Nz=", Nz, " ne=", ne];

elements = Flatten[Table[{(i - 0.5) * dr, (j - 0.5) * dz - t/2},
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

Print["  Building K_sym (", ne, "×", ne, ", will take a while)..."];
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

Print["  Solving symmetric eigenvalue problem (top 12)..."];
t0 = AbsoluteTime[];
{eigvals, eigvecs} = Eigensystem[N[Ksym], 12];
Print["    Done in ", AbsoluteTime[] - t0, " s"];

positiveEigs = Select[Re[eigvals], # > 0 &];
tausUs = N[sigma * positiveEigs * 10^6, 8];
tausUsSorted = Reverse[Sort[tausUs]];

Print[""];
Print["=== BEM Foster pole spectrum (refined) ==="];
Do[
  If[n <= Length[tausUsSorted],
    Print["  tau_", n, " = ", NumberForm[tausUsSorted[[n]], {6, 4}], " us"]
  ],
  {n, 1, 12}
];

(* Analytical reference: PEC disk with both k_r (Bessel) and k_z (axial).
   τ_{m,n} = μ_0 σ a² / [j_{1,m}² + (n π a/t)²]
   m: radial mode (1, 2, 3, ...), n: axial mode (0, 1, 2, ...)
   For axisymmetric J_φ, n=0 means uniform in z (no z-mode).
   n>=1 with sin(n π z/t) z-dependence. *)
j1zeros = N[BesselJZero[1, Range[1, 6]]];

pecModes = Flatten[Table[
  Module[{tau, m1Sq, kzSq},
    m1Sq = j1zeros[[m]]^2;
    kzSq = (nz Pi R / t)^2;
    tau = mu0 sigma R^2 / (m1Sq + kzSq) * 10^6;
    {tau, m, nz}
  ],
  {m, 1, 6}, {nz, 0, 6}
], 1];
pecModesSorted = Reverse[SortBy[pecModes, First]];

Print[""];
Print["=== Reference: PEC disk analytical Foster poles (k_r + k_z) ==="];
Print["  τ_{m,n} = μ_0 σ a² / [j_{1,m}² + (n π a/t)²]"];
Print["  Top 12 modes sorted by τ (decreasing):"];
Do[
  Module[{rec = pecModesSorted[[k]]},
    Print["  tau_", k, " = ", NumberForm[rec[[1]], {6, 4}],
          " us  (m=", rec[[2]], ", n=", rec[[3]], ")"]
  ],
  {k, 1, 12}
];

(* Side-by-side *)
v22bTauB2 = {208.32, 77.90, 46.41, 33.19, 25.22, 18.23};
Print[""];
Print["=== Comparison ==="];
Print[StringJoin[
  StringPadRight["n", 4],
  StringPadRight["BEM v3 us", 14],
  StringPadRight["v22b us", 14],
  StringPadRight["PEC analytical us", 18],
  "PEC (m,n)"
]];
Do[
  Module[{bem, v22, pecRow, pecTau, pecMn},
    bem = If[n <= Length[tausUsSorted], tausUsSorted[[n]], 0];
    v22 = If[n <= Length[v22bTauB2], v22bTauB2[[n]], 0];
    pecRow = pecModesSorted[[n]];
    pecTau = pecRow[[1]];
    pecMn = StringJoin["(", ToString[pecRow[[2]]], ",", ToString[pecRow[[3]]], ")"];
    Print[StringJoin[
      StringPadRight[ToString[n], 4],
      StringPadRight[ToString[NumberForm[bem, {6, 3}]], 14],
      StringPadRight[ToString[NumberForm[v22, {6, 3}]], 14],
      StringPadRight[ToString[NumberForm[pecTau, {6, 3}]], 18],
      pecMn
    ]]
  ],
  {n, 1, 6}
];

(* Save JSON *)
result = <|
  "method" -> "BEM v3 refined + PEC analytical (k_r, k_z)",
  "params" -> <|"R_mm" -> R*1000, "t_mm" -> t*1000, "sigma" -> sigma,
                "NR" -> NR, "Nz" -> Nz, "ne" -> ne|>,
  "bem_tau_us" -> tausUsSorted,
  "pec_analytical_modes" -> pecModesSorted[[1;;12]],
  "v22b_tau_B2_for_compare" -> v22bTauB2
|>;
Export["bem_disk_axisym_v3_refined.json", result, "JSON"];
Print[""];
Print["Saved: bem_disk_axisym_v3_refined.json"];
