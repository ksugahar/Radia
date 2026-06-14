#!/usr/bin/env wolframscript
(* High-precision K[0,0] computation to compare with C++ Spherical Duffy. *)

ClearAll[x, y, z, xp, yp, zp];

xi  = 2 x  - 1; eta  = 2 y  - 1; zeta  = 2 z  - 1;
xip = 2 xp - 1; etap = 2 yp - 1; zetap = 2 zp - 1;

viX = 8 (xi^2 - 1)/2 * eta * zeta;
viY = 0;
viZ = -8 xi * eta * (zeta^2 - 1)/2;

vjX = viX /. {x -> xp, y -> yp, z -> zp};
vjY = viY /. {x -> xp, y -> yp, z -> zp};
vjZ = viZ /. {x -> xp, y -> yp, z -> zp};

(* Cuboid 5x2x1 mm *)
a = 5*10^-3; b = 2*10^-3; c = 1*10^-3;
dotProduct = viX vjX a^2 + viY vjY b^2 + viZ vjZ c^2;
kernel = 1 / Sqrt[a^2 (x - xp)^2 + b^2 (y - yp)^2 + c^2 (z - zp)^2];

Print["High-precision K[0,0] for cuboid 5x2x1 mm, Block A (0,0,0):"];
Print["===================================================================="];
Print["WP=15, PG=3, AG=5, MR=6 (original .wls settings):"];
result1 = Quiet[NIntegrate[dotProduct kernel,
  {x, 0, 1}, {y, 0, 1}, {z, 0, 1}, {xp, 0, 1}, {yp, 0, 1}, {zp, 0, 1},
  Method -> {"GlobalAdaptive", "SymbolicProcessing" -> 0},
  WorkingPrecision -> 15, PrecisionGoal -> 3, AccuracyGoal -> 5,
  MaxRecursion -> 6]];
Print["  result = ", result1];

Print[];
Print["WP=20, PG=6, AG=8, MR=10 (tighter):"];
t0 = AbsoluteTime[];
result2 = Quiet[NIntegrate[dotProduct kernel,
  {x, 0, 1}, {y, 0, 1}, {z, 0, 1}, {xp, 0, 1}, {yp, 0, 1}, {zp, 0, 1},
  Method -> {"GlobalAdaptive", "SymbolicProcessing" -> 0},
  WorkingPrecision -> 20, PrecisionGoal -> 6, AccuracyGoal -> 8,
  MaxRecursion -> 10]];
Print["  result = ", result2, "  (", Round[AbsoluteTime[]-t0, 0.1], " s)"];

Print[];
Print["WP=25, PG=8, AG=10, MR=14 (very tight):"];
t0 = AbsoluteTime[];
result3 = Quiet[NIntegrate[dotProduct kernel,
  {x, 0, 1}, {y, 0, 1}, {z, 0, 1}, {xp, 0, 1}, {yp, 0, 1}, {zp, 0, 1},
  Method -> {"GlobalAdaptive", "SymbolicProcessing" -> 0},
  WorkingPrecision -> 25, PrecisionGoal -> 8, AccuracyGoal -> 10,
  MaxRecursion -> 14]];
Print["  result = ", result3, "  (", Round[AbsoluteTime[]-t0, 0.1], " s)"];

Print[];
Print["===================================================================="];
Print["C++ Spherical Duffy reports converging to ~1.243e-3"];
Print["Comparison ratio analysis (high prec vs C++):"];
Print["  WP=20:  ", N[result2 / 1.243*^-3, 6]];
Print["  WP=25:  ", N[result3 / 1.243*^-3, 6]];
