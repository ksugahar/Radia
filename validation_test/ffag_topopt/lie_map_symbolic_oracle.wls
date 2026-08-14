#!/usr/bin/env wolframscript

(* Independent symbolic oracle for the canonical Radia R/T/U/V Lie map.

   Mathematica is deliberately a derivation/golden tool only.  Production
   evaluation is C++/Python/MATLAB and does not invoke this script.
*)

ClearAll["Global`*"];

variables = {x, px, y, py, ell, delta};
epsilon = Unique["epsilon"];
h = 1/15;
beta = 4/5;
k1 = 1/3;
ks1 = 1/10;
k2 = 2/3;
ks2 = 7/30;
k3 = 4/3;
ks3 = -1;

coordinate = x + I y;
potential = h x + h^2 x^2/2 +
  ComplexExpand[Re[(k1 + I ks1) coordinate^2/2]] +
  ComplexExpand[Re[(k2 + I ks2) coordinate^3/3]] +
  ComplexExpand[Re[(k3 + I ks3) coordinate^4/4]];
massSquare = 1/beta^2 - 1;
hamiltonian = -(1 + h x) Sqrt[(1 + delta)^2 - px^2 - py^2] +
  potential + Sqrt[(1 + delta)^2 + massSquare]/beta;
scaled = hamiltonian /. Thread[variables -> epsilon variables];
series = Normal[Series[scaled, {epsilon, 0, 4}]] // Expand;
H2 = Coefficient[series, epsilon, 2];
H3 = Coefficient[series, epsilon, 3];
H4 = Coefficient[series, epsilon, 4];

tensorEntry[polynomial_, indices_List] :=
  Fold[D[#1, #2] &, polynomial, variables[[indices]]];

poisson = ConstantArray[0, {6, 6}];
poisson[[1, 2]] = 1; poisson[[2, 1]] = -1;
poisson[[3, 4]] = 1; poisson[[4, 3]] = -1;
poisson[[5, 6]] = -1; poisson[[6, 5]] = 1;
hamiltonian2Hessian = Table[D[H2, variables[[i]], variables[[j]]], {i, 6}, {j, 6}];
linearGenerator = Simplify[poisson . hamiltonian2Hessian];

checkDifferences = {
  tensorEntry[H2, {1, 1}] - (h^2 + k1),
  tensorEntry[H2, {1, 3}] - (-ks1),
  tensorEntry[H2, {1, 6}] - (-h),
  tensorEntry[H2, {6, 6}] - (1 - beta^2),
  tensorEntry[H3, {2, 2, 6}] - (-1),
  tensorEntry[H3, {1, 2, 2}] - h,
  tensorEntry[H3, {1, 1, 1}] - 2 k2,
  tensorEntry[H3, {1, 1, 3}] - (-2 ks2),
  tensorEntry[H4, {2, 2, 2, 2}] - 3,
  tensorEntry[H4, {2, 2, 4, 4}] - 1,
  tensorEntry[H4, {1, 2, 2, 6}] - (-h),
  tensorEntry[H4, {1, 1, 1, 1}] - 6 k3,
  tensorEntry[H4, {1, 1, 1, 3}] - (-6 ks3),
  tensorEntry[H4, {6, 6, 6, 6}] -
    3 beta^2 (1 - beta^2) (5 beta^2 - 1)
};
If[!And @@ (TrueQ[Simplify[#] == 0] & /@ checkDifferences),
  Print["Hamiltonian coefficient check failed"];
  Print[InputForm[Simplify[checkDifferences]]];
  Exit[1]
];

(* Non-trivial f3 self cascade plus independent f4 and f5 kicks. *)
a = 2/5;
b = -3/7;
c = 5/11;
f3sample = a q^2 p/2;
f4sample = b q^4/4;
f5sample = c q^5/5;
lie[generator_, observable_] :=
  D[observable, q] D[generator, p] - D[observable, p] D[generator, q];
flow3[observable_] := Expand[
  observable + lie[f3sample, observable] +
  lie[f3sample, lie[f3sample, observable]]/2 +
  lie[f3sample, lie[f3sample, lie[f3sample, observable]]]/6
];
qAfter4 = q;
pAfter4 = p + lie[f5sample, p] + lie[f4sample, p];
(* Coordinate-map composition is exp(f3) after exp(f4) after exp(f5).
   The pullback/substitution order is therefore the reverse of applying the
   Lie operator directly to pAfter4. *)
qMap = flow3[q] /. {q -> qAfter4, p -> pAfter4};
pMap = flow3[p] /. {q -> qAfter4, p -> pAfter4};
pMap = Total[Select[MonomialList[Expand[pMap], {q, p}],
  Total[Exponent[#, {q, p}]] <= 4 &]];

numeric[value_] := N[value, 17];
data = <|
  "schema" -> "radia.accelerator-lie-symbolic/v1",
  "coordinates" -> {"x", "px_over_p0", "y", "py_over_p0", "ell", "delta"},
  "longitudinal_pair_sign" -> -1,
  "parameters" -> <|
    "curvature" -> numeric[h], "reference_beta" -> numeric[beta],
    "normal_quadrupole" -> numeric[k1], "skew_quadrupole" -> numeric[ks1],
    "normal_sextupole" -> numeric[k2], "skew_sextupole" -> numeric[ks2],
    "normal_octupole" -> numeric[k3], "skew_octupole" -> numeric[ks3]
  |>,
  "hamiltonian_tensor_entries" -> <|
    "H2_xx" -> numeric[tensorEntry[H2, {1, 1}]],
    "H2_xy" -> numeric[tensorEntry[H2, {1, 3}]],
    "H2_x_delta" -> numeric[tensorEntry[H2, {1, 6}]],
    "H2_delta_delta" -> numeric[tensorEntry[H2, {6, 6}]],
    "H3_px_px_delta" -> numeric[tensorEntry[H3, {2, 2, 6}]],
    "H3_x_px_px" -> numeric[tensorEntry[H3, {1, 2, 2}]],
    "H3_x_x_x" -> numeric[tensorEntry[H3, {1, 1, 1}]],
    "H3_x_x_y" -> numeric[tensorEntry[H3, {1, 1, 3}]],
    "H4_px_px_px_px" -> numeric[tensorEntry[H4, {2, 2, 2, 2}]],
    "H4_px_px_py_py" -> numeric[tensorEntry[H4, {2, 2, 4, 4}]],
    "H4_x_px_px_delta" -> numeric[tensorEntry[H4, {1, 2, 2, 6}]],
    "H4_x_x_x_x" -> numeric[tensorEntry[H4, {1, 1, 1, 1}]],
    "H4_x_x_x_y" -> numeric[tensorEntry[H4, {1, 1, 1, 3}]],
    "H4_delta_delta_delta_delta" -> numeric[tensorEntry[H4, {6, 6, 6, 6}]]
  |>,
  "linear_generator_entries" -> <|
    "A_x_px" -> numeric[linearGenerator[[1, 2]]],
    "A_px_x" -> numeric[linearGenerator[[2, 1]]],
    "A_px_y" -> numeric[linearGenerator[[2, 3]]],
    "A_px_delta" -> numeric[linearGenerator[[2, 6]]],
    "A_y_py" -> numeric[linearGenerator[[3, 4]]],
    "A_py_x" -> numeric[linearGenerator[[4, 1]]],
    "A_py_y" -> numeric[linearGenerator[[4, 3]]],
    "A_ell_x" -> numeric[linearGenerator[[5, 1]]],
    "A_ell_delta" -> numeric[linearGenerator[[5, 6]]]
  |>,
  "lie_sample" -> <|
    "f3_a" -> numeric[a], "f4_b" -> numeric[b], "f5_c" -> numeric[c],
    "q_map_coefficients" -> <|
      "q" -> numeric[Coefficient[qMap, q]],
      "q2" -> numeric[Coefficient[qMap, q, 2]],
      "q3" -> numeric[Coefficient[qMap, q, 3]],
      "q4" -> numeric[Coefficient[qMap, q, 4]]
    |>,
    "p_map_coefficients" -> <|
      "p" -> numeric[Coefficient[pMap, p] /. q -> 0],
      "q_p" -> numeric[Coefficient[Coefficient[pMap, p], q]],
      "q2_p" -> numeric[Coefficient[Coefficient[pMap, p], q, 2]],
      "q3_p" -> numeric[Coefficient[Coefficient[pMap, p], q, 3]],
      "q3" -> numeric[Coefficient[pMap /. p -> 0, q, 3]],
      "q4" -> numeric[Coefficient[pMap /. p -> 0, q, 4]]
    |>
  |>
|>;

output = FileNameJoin[{DirectoryName[$InputFileName], "lie_map_symbolic_reference.json"}];
Export[output, data, "RawJSON"];
Print["wrote " <> output];
