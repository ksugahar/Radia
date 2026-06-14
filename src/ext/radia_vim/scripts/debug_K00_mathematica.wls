#!/usr/bin/env wolframscript
(* Compute K[0,0] for unit cube, order 1, directly via NIntegrate.
   Reference for my C++ Spherical Duffy debug. *)

ClearAll[x, y, z, xp, yp, zp];

(* Order 1, Block A (i=j=k=0) basis function *)
xi  = 2 x  - 1; eta  = 2 y  - 1; zeta  = 2 z  - 1;
xip = 2 xp - 1; etap = 2 yp - 1; zetap = 2 zp - 1;

(* L_2(s) = (s^2-1)/2,  L_2'(s) = s,  P_1(s) = s *)
viX = 8 (xi^2 - 1)/2 * eta * zeta;       (* 8 L_2(xi) P_1(eta) L_2'(zeta) *)
viY = 0;
viZ = -8 xi * eta * (zeta^2 - 1)/2;       (* -8 L_2'(xi) P_1(eta) L_2(zeta) *)

vjX = viX /. {x -> xp, y -> yp, z -> zp};
vjY = viY /. {x -> xp, y -> yp, z -> zp};
vjZ = viZ /. {x -> xp, y -> yp, z -> zp};

(* Unit cube *)
a = 1; b = 1; c = 1;
dotProduct = viX vjX a^2 + viY vjY b^2 + viZ vjZ c^2;
kernel = 1 / Sqrt[a^2 (x - xp)^2 + b^2 (y - yp)^2 + c^2 (z - zp)^2];

Print["Computing K[0,0] for unit cube, order 1 ..."];
result = NIntegrate[dotProduct kernel,
  {x, 0, 1}, {y, 0, 1}, {z, 0, 1},
  {xp, 0, 1}, {yp, 0, 1}, {zp, 0, 1},
  Method -> {"GlobalAdaptive", "SymbolicProcessing" -> 0},
  WorkingPrecision -> 16, PrecisionGoal -> 6, AccuracyGoal -> 8,
  MaxRecursion -> 10];
Print["NIntegrate result: ", result];

(* Also compute for cuboid 5x2x1 mm, order 3 idx 0 (Block A i=j=k=0) -- same basis since same block *)
a2 = 5*10^-3; b2 = 2*10^-3; c2 = 1*10^-3;
dotProduct2 = viX vjX a2^2 + viY vjY b2^2 + viZ vjZ c2^2;
kernel2 = 1 / Sqrt[a2^2 (x - xp)^2 + b2^2 (y - yp)^2 + c2^2 (z - zp)^2];
Print[];
Print["Computing K[0,0] for cuboid 5x2x1 mm, order 1 (Block A 0,0,0) ..."];
result2 = NIntegrate[dotProduct2 kernel2,
  {x, 0, 1}, {y, 0, 1}, {z, 0, 1},
  {xp, 0, 1}, {yp, 0, 1}, {zp, 0, 1},
  Method -> {"GlobalAdaptive", "SymbolicProcessing" -> 0},
  WorkingPrecision -> 16, PrecisionGoal -> 6, AccuracyGoal -> 8,
  MaxRecursion -> 10];
Print["NIntegrate result: ", result2];
