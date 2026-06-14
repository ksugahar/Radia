#!/usr/bin/env wolframscript
(*
  Axisymmetric BEM Foster-CLN for Cu disk in vacuum (2026-05-05).

  Computes Foster poles τ_n of eddy current decay in an axisymmetric
  conducting disk via:
    Generalized eigenvalue: τ × M_mass × J = σ × M_BEM × J

  where M_BEM[i, j] = mutual inductance between current ring i and j,
  using axisymmetric Green's function with complete elliptic integrals K, E.

  Output: τ_n = Foster pole spectrum (μs), to be compared with NGSolve v22b
  axisymmetric Cu disk (R=10mm, t=2mm) per-stage Hiruma L_B2 form values.
*)

mu0 = 4 Pi * 10^-7;
sigma = 5.8 * 10^7;        (* Cu *)
R = 10 * 10^-3;            (* disk radius [m] *)
t = 2 * 10^-3;             (* disk thickness [m] *)

(* Discretization *)
NR = 30;        (* radial divisions *)
Nz = 6;         (* axial divisions *)
dr = R / NR;
dz = t / Nz;
ne = NR * Nz;

Print["=== Axisymmetric BEM Foster for Cu disk ==="];
Print["  R=", R*1000, "mm, t=", t*1000, "mm, sig=", sigma];
Print["  NR=", NR, ", Nz=", Nz, ", ne=", ne];

(* Element centroids (r_i, z_i) for i = 1..ne *)
elements = Flatten[Table[{(i - 0.5) * dr, (j - 0.5) * dz - t/2},
                         {i, 1, NR}, {j, 1, Nz}], 1];

(* Mutual inductance kernel for current loops at (r1, z1) and (r2, z2):
   M = mu0 * sqrt(r1 r2) * [(2/k - k) K(k^2) - (2/k) E(k^2)]
   where k^2 = 4 r1 r2 / ((r1+r2)^2 + (z1-z2)^2)
   We compute G_BEM = M / mu0 (drop constants for eigenvalue problem since
   they appear on both sides via σ × mu0). Actually we keep mu0 for
   correct physical units of τ. *)

mutualG[r1_, z1_, r2_, z2_] := Module[{k2, K1, E1, eps = 1.0*^-12},
  k2 = 4 r1 r2 / ((r1 + r2)^2 + (z1 - z2)^2 + eps);
  K1 = EllipticK[k2];
  E1 = EllipticE[k2];
  Sqrt[r1 r2] * ((2/Sqrt[k2] - Sqrt[k2]) K1 - (2/Sqrt[k2]) E1)
];

(* Self-term approximation: for ring of radius r with cross-section dr × dz,
   self-inductance (per unit current density) is well-approximated by
   M_self ≈ mu0 r [ln(8 r / a_eq) - 7/4]
   where a_eq is equivalent cross-section radius: a_eq = sqrt(dr × dz / pi) *)

aEq = Sqrt[dr * dz / Pi];
selfM[r_] := mu0 r (Log[8 r / aEq] - 7/4);

Print["  Building BEM matrix (", ne, "×", ne, ")..."];
t0 = AbsoluteTime[];

(* M_BEM[i, j] = mu0 × G_axi(elem_i, elem_j) for i ≠ j
                 = self_M(r_i)             for i == j *)
Mbem = Table[
  Module[{r1, z1, r2, z2, val},
    {r1, z1} = elements[[i]];
    {r2, z2} = elements[[j]];
    If[i === j,
      val = selfM[r1] / (2 Pi),  (* normalize: integrand vs ring inductance *)
      val = mu0 * mutualG[r1, z1, r2, z2]
    ];
    val
  ],
  {i, 1, ne}, {j, 1, ne}
];

Print["    BEM matrix built in ", AbsoluteTime[] - t0, " s"];

(* Mass matrix: diagonal, M_mass[i, i] = 2 pi r_i dr dz (cylindrical volume) *)
Mmass = DiagonalMatrix[Table[2 Pi elements[[i, 1]] dr dz, {i, 1, ne}]];

(* Generalized eigenvalue problem:
   τ × M_mass × J = σ × M_BEM × J
   Equivalent to: σ × M_BEM × J = τ^{-1} × M_mass × J ??? Let me re-derive.

   For Foster pole decay J(t) = J_0 × exp(-t/τ):
   ∂J/∂t = -J/τ
   Faraday + Ohm: J = σ × E = -σ × ∂A/∂t (in conductor)
   A(P) = ∫ G_axi(P, P') × J(P') × dV' (vacuum vector potential from current J)
   So: J = -σ × ∂/∂t [∫ G × J × dV']
   For decay: -J/τ = -σ × ∫ G × J × dV' (with sign flip for ∂/∂t of exp(-t/τ))
   J/τ = σ × ∫ G × J × dV'

   Discretized:
     (1/τ) × M_mass × J = σ × M_BEM × J  (where M_mass = diag(2 π r dr dz))
     τ × σ × M_BEM × J = M_mass × J        (rearranged)

   Or equivalently: σ × M_BEM × J = (1/τ) × M_mass × J
   Eigenvalues e_n = 1/τ_n → τ_n = 1/e_n.
*)

Print["  Solving generalized eigenvalue problem..."];
t0 = AbsoluteTime[];
{eigvals, eigvecs} = Eigensystem[
  {N[sigma * Mbem], N[Mmass]}, 6
];
Print["    Eigensolve done in ", AbsoluteTime[] - t0, " s"];

(* Convert: e_n = 1/τ_n → τ_n = 1/e_n in seconds *)
taus = 1 / Re[eigvals];
tausUs = N[taus, 8] * 10^6;  (* μs *)
tausUsSorted = Reverse[Sort[Select[tausUs, # > 0 &]]];

Print[""];
Print["=== Foster pole spectrum (Cu disk R=", R*1000, "mm, t=", t*1000, "mm) ==="];
Print["  Top ", Length[tausUsSorted], " τ values (μs):"];
Do[
  Print["  τ_", n, " = ", tausUsSorted[[n]], " μs"],
  {n, 1, Length[tausUsSorted]}
];

(* Comparison with NGSolve v22b per-stage τ_B2 *)
v22bTauB2 = {208.32, 77.90, 46.41, 33.19, 25.22, 18.23};
Print[""];
Print["=== Comparison with NGSolve v22b axisymmetric Cu disk ==="];
Print["  Note: v22b uses Hiruma per-stage L_B2 form (canonical)"];
Print["  ", StringPadRight["n", 4], StringPadRight["BEM τ μs", 14],
        StringPadRight["v22b τ_B2 μs", 16], "ratio"];
Do[
  If[n <= Length[tausUsSorted] && n <= Length[v22bTauB2],
    Module[{ratio = tausUsSorted[[n]] / v22bTauB2[[n]]},
      Print["  ",
        StringPadRight[ToString[n], 4],
        StringPadRight[ToString[NumberForm[tausUsSorted[[n]], {6, 4}]], 14],
        StringPadRight[ToString[v22bTauB2[[n]]], 16],
        NumberForm[ratio, {5, 4}]
      ]
    ]
  ],
  {n, 1, 6}
];

(* Save JSON *)
result = <|
  "method" -> "Axisymmetric BEM Foster-CLN for Cu disk",
  "params" -> <|
    "R_mm" -> R * 1000,
    "t_mm" -> t * 1000,
    "sigma" -> sigma,
    "NR" -> NR,
    "Nz" -> Nz,
    "ne" -> ne
  |>,
  "tau_us" -> tausUsSorted,
  "v22b_tau_B2_us_for_compare" -> v22bTauB2
|>;
Export["bem_disk_axisym_foster.json", result, "JSON"];
Print[""];
Print["Saved: bem_disk_axisym_foster.json"];
