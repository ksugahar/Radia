(* Test Aharoni extended for polynomial x^k x'^l. *)

ClearAll[x, xp, len, rho];
$MaxExtraPrecision = 50;

Print["Test 1: Aharoni base k=l=0"];
t0 = AbsoluteTime[];
F00 = Integrate[1/Sqrt[(x-xp)^2 + rho^2], {x, 0, len}, {xp, 0, len},
                Assumptions -> {len > 0, rho > 0}];
Print["  time: ", Round[AbsoluteTime[] - t0, 0.1], " s"];
Print["  F_00 = ", F00];

Print[];
Print["Test 2: k=1, l=0  (x in numerator)"];
t0 = AbsoluteTime[];
F10 = Integrate[x/Sqrt[(x-xp)^2 + rho^2], {x, 0, len}, {xp, 0, len},
                Assumptions -> {len > 0, rho > 0}];
Print["  time: ", Round[AbsoluteTime[] - t0, 0.1], " s"];
Print["  F_10 = ", F10];

Print[];
Print["Test 3: k=1, l=1  (x*xp)"];
t0 = AbsoluteTime[];
F11 = Integrate[x xp/Sqrt[(x-xp)^2 + rho^2], {x, 0, len}, {xp, 0, len},
                Assumptions -> {len > 0, rho > 0}];
Print["  time: ", Round[AbsoluteTime[] - t0, 0.1], " s"];
Print["  F_11 LeafCount: ", LeafCount[F11]];

Print[];
Print["Test 4: k=2, l=0"];
t0 = AbsoluteTime[];
F20 = Integrate[x^2/Sqrt[(x-xp)^2 + rho^2], {x, 0, len}, {xp, 0, len},
                Assumptions -> {len > 0, rho > 0}];
Print["  time: ", Round[AbsoluteTime[] - t0, 0.1], " s"];
Print["  F_20 LeafCount: ", LeafCount[F20]];

Print[];
Print["All tests done. If all complete in seconds, we can use this approach."];
