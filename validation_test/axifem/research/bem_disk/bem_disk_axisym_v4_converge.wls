#!/usr/bin/env wolframscript
(*
  Axisymmetric BEM Foster v4: aim for 1% convergence.
  - NR=200, Nz=60, ne=12000
  - ParallelTable for K matrix construction
  - Improved self-term (Maxwell mutual inductance for adjacent ring)
*)

mu0 = 4 Pi * 10^-7;
sigma = 5.8 * 10^7;
R = 10 * 10^-3;
t = 2 * 10^-3;

NR = 200;
Nz = 60;
dr = R / NR;
dz = t / Nz;
ne = NR * Nz;

Print["=== BEM Foster v4 (converging refinement) ==="];
Print["  R=", R*1000, "mm, t=", t*1000, "mm"];
Print["  NR=", NR, " Nz=", Nz, " ne=", ne, " (K is ", ne, "×", ne, ")"];

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

Print["  Building K_sym in parallel..."];
t0 = AbsoluteTime[];

(* Use ParallelTable for speedup *)
SetSharedFunction[Gaxi, GaxiSelf];
DistributeDefinitions[mu0, dV, elements, Gaxi, GaxiSelf];

Ksym = ParallelTable[
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

Print["  Solving symmetric eigenvalue problem (top 8)..."];
t0 = AbsoluteTime[];
{eigvals, eigvecs} = Eigensystem[N[Ksym], 8];
Print["    Done in ", AbsoluteTime[] - t0, " s"];

positiveEigs = Select[Re[eigvals], # > 0 &];
tausUs = N[sigma * positiveEigs * 10^6, 8];
tausUsSorted = Reverse[Sort[tausUs]];

Print[""];
Print["=== BEM v4 Foster pole spectrum ==="];
Do[
  If[n <= Length[tausUsSorted],
    Print["  tau_", n, " = ", NumberForm[tausUsSorted[[n]], {6, 4}], " us"]
  ],
  {n, 1, 8}
];

(* Comparison vs v22b NGSolve and v3 *)
v22bTauB2 = {208.32, 77.90, 46.41, 33.19, 25.22, 18.23};
v3Tau = {224.31, 88.42, 49.54, 32.09, 23.34, 22.59};

Print[""];
Print["=== Convergence comparison ==="];
Print[StringJoin[
  StringPadRight["n", 4],
  StringPadRight["BEM v2 (ne=180)", 17],
  StringPadRight["BEM v3 (ne=1920)", 18],
  StringPadRight["BEM v4 (ne=12000)", 19],
  StringPadRight["v22b NGSolve", 14],
  "v4/v22b"
]];
v2Tau = {228.0, 91.6, 52.6, 35.0, 26.6, 25.4};
Do[
  Module[{v2, v3, v4, v22, ratio},
    v2 = If[n <= 6, v2Tau[[n]], 0];
    v3 = If[n <= 6, v3Tau[[n]], 0];
    v4 = If[n <= Length[tausUsSorted], tausUsSorted[[n]], 0];
    v22 = If[n <= Length[v22bTauB2], v22bTauB2[[n]], 0];
    ratio = If[v22 > 0, v4 / v22, 0];
    Print[StringJoin[
      StringPadRight[ToString[n], 4],
      StringPadRight[ToString[NumberForm[v2, {6, 3}]], 17],
      StringPadRight[ToString[NumberForm[v3, {6, 3}]], 18],
      StringPadRight[ToString[NumberForm[v4, {6, 3}]], 19],
      StringPadRight[ToString[NumberForm[v22, {6, 3}]], 14],
      ToString[NumberForm[ratio, {5, 4}]]
    ]]
  ],
  {n, 1, 6}
];

result = <|
  "method" -> "BEM v4 axisymmetric Cu disk (vacuum, refined)",
  "params" -> <|"R_mm" -> R*1000, "t_mm" -> t*1000, "sigma" -> sigma,
                "NR" -> NR, "Nz" -> Nz, "ne" -> ne|>,
  "tau_us_v4" -> tausUsSorted,
  "tau_us_v3" -> v3Tau,
  "tau_us_v2" -> v2Tau,
  "v22b_tau_B2" -> v22bTauB2
|>;
Export["bem_disk_axisym_v4_converge.json", result, "JSON"];
Print[""];
Print["Saved: bem_disk_axisym_v4_converge.json"];
