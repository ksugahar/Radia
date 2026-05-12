(* Test if Mathematica can do the inner Newton-kernel integration symbolically. *)

ClearAll[x, y, z, xp, yp, zp];
$MaxExtraPrecision = 50;

(* Test 1: simplest case - constant density on unit cube *)
Print["Test 1: int_[0,1]^3 1/|r-r'| dV' = Phi_const(x,y,z)"];
t0 = AbsoluteTime[];
phi1 = Integrate[1/Sqrt[(x-xp)^2 + (y-yp)^2 + (z-zp)^2],
                 {xp, 0, 1}, {yp, 0, 1}, {zp, 0, 1},
                 Assumptions -> {0 < x < 1, 0 < y < 1, 0 < z < 1}];
Print["  time: ", Round[AbsoluteTime[] - t0, 0.1], " s"];
Print["  result type: ", Head[phi1]];
Print["  evaluated at (0.5, 0.5, 0.5): ", N[phi1 /. {x->1/2, y->1/2, z->1/2}, 16]];

(* Test 2: linear polynomial *)
Print[];
Print["Test 2: int_[0,1]^3 yp/|r-r'| dV'"];
t0 = AbsoluteTime[];
phi2 = Integrate[yp/Sqrt[(x-xp)^2 + (y-yp)^2 + (z-zp)^2],
                 {xp, 0, 1}, {yp, 0, 1}, {zp, 0, 1},
                 Assumptions -> {0 < x < 1, 0 < y < 1, 0 < z < 1}];
Print["  time: ", Round[AbsoluteTime[] - t0, 0.1], " s"];
Print["  evaluated at (0.5, 0.5, 0.5): ", N[phi2 /. {x->1/2, y->1/2, z->1/2}, 16]];

(* Test 3: bilinear polynomial - HCurl edge 0-1 inner *)
Print[];
Print["Test 3: int_[0,1]^3 (1-yp)(1-zp)/|r-r'| dV'"];
t0 = AbsoluteTime[];
phi3 = Integrate[(1-yp)(1-zp)/Sqrt[(x-xp)^2 + (y-yp)^2 + (z-zp)^2],
                 {xp, 0, 1}, {yp, 0, 1}, {zp, 0, 1},
                 Assumptions -> {0 < x < 1, 0 < y < 1, 0 < z < 1}];
Print["  time: ", Round[AbsoluteTime[] - t0, 0.1], " s"];
Print["  result LeafCount: ", LeafCount[phi3]];
Print["  evaluated at (0.5, 0.5, 0.5): ", N[phi3 /. {x->1/2, y->1/2, z->1/2}, 16]];
