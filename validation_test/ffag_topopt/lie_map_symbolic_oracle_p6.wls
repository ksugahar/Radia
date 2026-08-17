#!/usr/bin/env wolframscript

(* Independent symbolic oracle for the ARBITRARY-JET (p=6 HCurl A-map) Lie map.

   Extends lie_map_symbolic_oracle.wls (body multipoles, fourth order) to a
   general median-plane-symmetric transverse jet through total degree six:

       a_y(x,y) = Sum cy[i,j] x^i y^j   (1 <= i+j <= 6, j ODD),
       a_s(x,y) = Sum cs[i,j] x^i y^j   (1 <= i+j <= 6, j EVEN),

   inserted into the exact curvilinear parent Hamiltonian

       H = -(1+hh*x) Sqrt[(1+delta)^2 - px^2 - (py - a_y)^2]
           - a_s + Sqrt[(1+delta)^2 + m^2]/beta.

   hh is the Hamiltonian metric curvature: the production tracker uses the
   metric 1 - h_orbit*x, i.e. hh = -h_orbit, and feeds the COVARIANT
   longitudinal component (metric * physical A_s) as a_s.  The design-orbit
   condition H1 = 0 fixes cs[1,0] = -hh (and cs[0,1] is absent by parity).

   Deliverables written to lie_map_symbolic_reference_p6.json:
     1. exact H2..H6 tensor entries for the committed general jet
        (H6 is the p=6 novelty: degree-6 a_s, degree-5 a_y cross, a_y^2);
     2. the linear map R = MatrixExp[L A] (A = J.Hess H2), 25-digit;
     3. variational-jet map contractions R.u, T[u,u], U[u,u,u], V[u^4],
        W[u^5] for committed directions, from the exact flow of the
        degree<=6-truncated Hamiltonian around the fixed design orbit z=0
        (W is the fifth-order tensor a p=6 jet supports);
     4. a 1-DOF fifth-order Dragt--Finn sample: the composed map of
        KNOWN generators f3..f6 through phase-space degree 5, so a future
        dragt_finn_factorize_fifth_order must recover the committed f6.

   Self-tests: H1 vanishes symbolically; the multipole special case
   reproduces the committed fourth-order oracle's closed-form entries;
   R is symplectic; MatrixExp agrees with the variational jet's delta1.

   Mathematica is deliberately a derivation/golden tool only.  Production
   evaluation is C++/Python/MATLAB and does not invoke this script.
*)

ClearAll["Global`*"];

variables = {x, px, y, py, ell, delta};
epsilon = Unique["epsilon"];
hh = 1/15;
beta = 4/5;
segmentLength = 1/10;

(* ---- committed general parity-respecting jet ------------------------- *)
csCoefficient[i_, j_] := (-1)^(i + j) (3 i + 2 j + 2)/(13 + 7 i + 5 j);
cyCoefficient[i_, j_] := (-1)^(i + j + 1) (2 i + 3 j + 1)/(11 + 5 i + 7 j);

jetPairs[maxDegree_, parity_] := Flatten[
  Table[
    If[1 <= i + j <= maxDegree && Mod[j, 2] == parity, {i, j}, Nothing],
    {i, 0, maxDegree}, {j, 0, maxDegree}
  ], 1];

asPairs = jetPairs[6, 0];   (* a_s: even in y *)
ayPairs = jetPairs[6, 1];   (* a_y: odd in y *)

asJet = Sum[
  If[pair === {1, 0}, -hh, csCoefficient @@ pair] *
    x^pair[[1]] y^pair[[2]],
  {pair, asPairs}];
ayJet = Sum[cyCoefficient @@ pair x^pair[[1]] y^pair[[2]], {pair, ayPairs}];

massSquare = 1/beta^2 - 1;
hamiltonian = -(1 + hh x) Sqrt[(1 + delta)^2 - px^2 - (py - ayJet)^2] -
  asJet + Sqrt[(1 + delta)^2 + massSquare]/beta;

scaled = hamiltonian /. Thread[variables -> epsilon variables];
series = Normal[Series[scaled, {epsilon, 0, 6}]] // Expand;
H1 = Coefficient[series, epsilon, 1];
H2 = Coefficient[series, epsilon, 2];
H3 = Coefficient[series, epsilon, 3];
H4 = Coefficient[series, epsilon, 4];
H5 = Coefficient[series, epsilon, 5];
H6 = Coefficient[series, epsilon, 6];

If[!TrueQ[Simplify[H1] == 0],
  Print["design-orbit condition failed: H1 = ", InputForm[Simplify[H1]]];
  Exit[1]
];

tensorEntry[polynomial_, indices_List] :=
  Fold[D[#1, #2] &, polynomial, variables[[indices]]];

(* ---- self-test: the multipole special case reproduces the committed
       fourth-order oracle's closed forms ------------------------------- *)
Module[{k1 = 1/3, ks1 = 1/10, k2 = 2/3, ks2 = 7/30, k3 = 4/3, ks3 = -1,
        k4 = 5/4, ks4 = -2/5, potential, hamMultipole, scaledM, seriesM,
        H2m, H3m, H4m, H5m, checks},
  potential = hh x + hh^2 x^2/2 +
    ComplexExpand[Re[(k1 + I ks1) (x + I y)^2/2]] +
    ComplexExpand[Re[(k2 + I ks2) (x + I y)^3/3]] +
    ComplexExpand[Re[(k3 + I ks3) (x + I y)^4/4]] +
    ComplexExpand[Re[(k4 + I ks4) (x + I y)^5/5]];
  (* general engine with a_y = 0 and a_s = -potential *)
  hamMultipole = -(1 + hh x) Sqrt[(1 + delta)^2 - px^2 - py^2] +
    potential + Sqrt[(1 + delta)^2 + massSquare]/beta;
  scaledM = hamMultipole /. Thread[variables -> epsilon variables];
  seriesM = Normal[Series[scaledM, {epsilon, 0, 5}]] // Expand;
  {H2m, H3m, H4m, H5m} =
    Coefficient[seriesM, epsilon, #] & /@ {2, 3, 4, 5};
  checks = {
    tensorEntry[H2m, {1, 1}] - (hh^2 + k1),
    tensorEntry[H2m, {1, 3}] - (-ks1),
    tensorEntry[H2m, {1, 6}] - (-hh),
    tensorEntry[H2m, {6, 6}] - (1 - beta^2),
    tensorEntry[H3m, {2, 2, 6}] - (-1),
    tensorEntry[H3m, {1, 2, 2}] - hh,
    tensorEntry[H3m, {1, 1, 1}] - 2 k2,
    tensorEntry[H4m, {2, 2, 2, 2}] - 3,
    tensorEntry[H4m, {1, 1, 1, 1}] - 6 k3,
    tensorEntry[H5m, {1, 1, 1, 1, 1}] - 24 k4,
    tensorEntry[H5m, {6, 6, 6, 6, 6}] -
      15 beta^4 (1 - beta^2) (3 - 7 beta^2)
  };
  If[!And @@ (TrueQ[Simplify[#] == 0] & /@ checks),
    Print["multipole special-case check failed"];
    Print[InputForm[Simplify[checks]]];
    Exit[1]
  ]
];

(* ---- linear map and symplecticity ------------------------------------ *)
poisson = ConstantArray[0, {6, 6}];
poisson[[1, 2]] = 1; poisson[[2, 1]] = -1;
poisson[[3, 4]] = 1; poisson[[4, 3]] = -1;
poisson[[5, 6]] = -1; poisson[[6, 5]] = 1;
hessian = Table[D[H2, variables[[i]], variables[[j]]], {i, 6}, {j, 6}];
linearGenerator = poisson . hessian;
linearMap = MatrixExp[N[segmentLength linearGenerator, 30]];
symplecticResidual = Max[Abs[
  Transpose[linearMap] . poisson . linearMap - poisson]];
If[symplecticResidual > 10^-24,
  Print["MatrixExp linear map lost symplecticity: ",
    N[symplecticResidual, 3]];
  Exit[1]
];

(* ---- exact flow variational jet around z = 0 ------------------------- *)
truncated = H2 + H3 + H4 + H5 + H6;
vectorField = poisson . Table[D[truncated, v], {v, variables}];

(* z(s; e) = Sum_k e^k w_k(s)/k!, matched order-by-order in e.  z=0 is a
   fixed trajectory (H1=0), so w_0 = 0 and the jet is exact for the
   truncated Hamiltonian. *)
jetOrder = 5;
wSymbols = Table[w[k, i], {k, 1, jetOrder}, {i, 1, 6}];
zSubstitution = Thread[variables ->
  Table[Sum[epsilon^k wSymbols[[k, i]]/k!, {k, 1, jetOrder}], {i, 1, 6}]];
fieldSeries = Normal[Series[vectorField /. zSubstitution,
  {epsilon, 0, jetOrder}]];
jetRHS = Table[
  k! Coefficient[fieldSeries[[i]], epsilon, k],
  {k, 1, jetOrder}, {i, 1, 6}];

solveJet[u_List] := Module[{functions, rules, equations, initial, solution},
  functions = Table[wf[k, i], {k, 1, jetOrder}, {i, 1, 6}];
  rules = Flatten[Table[
    wSymbols[[k, i]] -> wf[k, i][s], {k, 1, jetOrder}, {i, 1, 6}]];
  equations = Flatten[Table[
    wf[k, i]'[s] == (jetRHS[[k, i]] /. rules),
    {k, 1, jetOrder}, {i, 1, 6}]];
  initial = Flatten[Table[
    wf[k, i][0] == If[k == 1, u[[i]], 0],
    {k, 1, jetOrder}, {i, 1, 6}]];
  solution = NDSolve[Join[equations, initial],
    Flatten[functions], {s, 0, segmentLength},
    WorkingPrecision -> 30, AccuracyGoal -> 25, PrecisionGoal -> 25,
    MaxSteps -> 100000][[1]];
  Table[wf[k, i][segmentLength] /. solution, {k, 1, jetOrder}, {i, 1, 6}]
];

(* committed probe directions *)
probeU = {1, 1/2, -1/3, 1/4, 1/5, -1/6};
probeV = {-1/2, 1/3, 1/4, -1/5, 1/6, 1/7};
jetU = solveJet[probeU];
jetV = solveJet[probeV];

(* delta1(L) must equal R.u *)
linearDefectU = Max[Abs[jetU[[1]] - linearMap . probeU]];
linearDefectV = Max[Abs[jetV[[1]] - linearMap . probeV]];
If[Max[linearDefectU, linearDefectV] > 10^-20,
  Print["variational delta1 disagrees with MatrixExp: ",
    N[{linearDefectU, linearDefectV}, 3]];
  Exit[1]
];

(* ---- 1-DOF fifth-order Dragt--Finn factorization sample -------------- *)
(* Known generators; the composed map through degree 5 is the golden a
   future dragt_finn_factorize_fifth_order must invert exactly. *)
f3sample = (2/5) q^2 p/2;
f4sample = (-3/7) q^4/4 + (5/9) q^2 p^2/2;
f5sample = (5/11) q^5/5 + (-4/13) q^3 p^2/3;
f6sample = (7/13) q^6/6 + (3/17) q^4 p^2/4;
lie[generator_, observable_] :=
  D[observable, q] D[generator, p] - D[observable, p] D[generator, q];
truncate5[expr_] := Total[Select[MonomialList[Expand[expr], {q, p}],
  Total[Exponent[#, {q, p}]] <= 5 &]];
lieFlow[generator_, iterates_][observable_] := truncate5[
  Sum[Nest[lie[generator, #] &, observable, n]/n!, {n, 0, iterates}]];
(* M = exp(:f3:) after exp(:f4:) after exp(:f5:) after exp(:f6:).
   Coordinate substitution composes in the reverse order. *)
flow3 = lieFlow[f3sample, 4];
flow4 = lieFlow[f4sample, 2];
flow5 = lieFlow[f5sample, 1];
flow6 = lieFlow[f6sample, 1];
qStage6 = flow6[q]; pStage6 = flow6[p];
qStage5 = truncate5[flow5[q] /. {q -> qStage6, p -> pStage6}];
pStage5 = truncate5[flow5[p] /. {q -> qStage6, p -> pStage6}];
qStage4 = truncate5[flow4[q] /. {q -> qStage5, p -> pStage5}];
pStage4 = truncate5[flow4[p] /. {q -> qStage5, p -> pStage5}];
qMap = truncate5[flow3[q] /. {q -> qStage4, p -> pStage4}];
pMap = truncate5[flow3[p] /. {q -> qStage4, p -> pStage4}];

mapCoefficients[expr_] := Association[Flatten[Table[
  If[i + j >= 1 && i + j <= 5,
    "q" <> ToString[i] <> "p" <> ToString[j] ->
      N[Coefficient[Coefficient[expr, q, i], p, j], 17],
    Nothing],
  {i, 0, 5}, {j, 0, 5}]]];

(* ---- JSON ------------------------------------------------------------ *)
(* Machine doubles only: arbitrary-precision tinies such as 0.``-29 export
   as the invalid JSON token "0.e-29", and every downstream consumer parses
   float64 anyway. *)
numeric[value_] := With[{machine = N[value]},
  If[machine == 0., 0, machine]];
jetAssociation[pairs_, valueFunction_] := Association[
  ("c" <> ToString[#[[1]]] <> "_" <> ToString[#[[2]]] ->
    numeric[valueFunction @@ #]) & /@ pairs];

data = <|
  "schema" -> "radia.accelerator-lie-symbolic-p6/v1",
  "coordinates" -> {"x", "px_over_p0", "y", "py_over_p0", "ell", "delta"},
  "longitudinal_pair_sign" -> -1,
  "conventions" -> <|
    "metric" -> "1 + hh*x with hh = -h_orbit (tracker metric 1 - h*x)",
    "a_s" -> "covariant longitudinal component; H1=0 fixes cs_1_0 = -hh",
    "parity" -> "a_s even in y, a_y odd in y (median-plane symmetric)",
    "map" -> "z_out = R z + T[z,z]/2 + U/6 + V/24 + W/120 + O(z^6)"
  |>,
  "parameters" -> <|
    "hh" -> numeric[hh], "reference_beta" -> numeric[beta],
    "segment_length" -> numeric[segmentLength]
  |>,
  "as_jet_coefficients" -> jetAssociation[asPairs,
    Function[{i, j}, If[{i, j} === {1, 0}, -hh, csCoefficient[i, j]]]],
  "ay_jet_coefficients" -> jetAssociation[ayPairs, cyCoefficient],
  "hamiltonian_tensor_entries" -> <|
    "H2_xx" -> numeric[tensorEntry[H2, {1, 1}]],
    "H2_xy" -> numeric[tensorEntry[H2, {1, 3}]],
    "H2_yy" -> numeric[tensorEntry[H2, {3, 3}]],
    "H2_y_py" -> numeric[tensorEntry[H2, {3, 4}]],
    "H2_x_delta" -> numeric[tensorEntry[H2, {1, 6}]],
    "H2_delta_delta" -> numeric[tensorEntry[H2, {6, 6}]],
    "H3_x_x_x" -> numeric[tensorEntry[H3, {1, 1, 1}]],
    "H3_x_y_y" -> numeric[tensorEntry[H3, {1, 3, 3}]],
    "H3_x_y_py" -> numeric[tensorEntry[H3, {1, 3, 4}]],
    "H3_y_y_py" -> numeric[tensorEntry[H3, {3, 3, 4}]],
    "H3_x_px_px" -> numeric[tensorEntry[H3, {1, 2, 2}]],
    "H3_px_px_delta" -> numeric[tensorEntry[H3, {2, 2, 6}]],
    "H4_x_x_x_x" -> numeric[tensorEntry[H4, {1, 1, 1, 1}]],
    "H4_x_x_y_py" -> numeric[tensorEntry[H4, {1, 1, 3, 4}]],
    "H4_y_y_py_py" -> numeric[tensorEntry[H4, {3, 3, 4, 4}]],
    "H4_px_px_py_py" -> numeric[tensorEntry[H4, {2, 2, 4, 4}]],
    "H5_x_x_x_x_x" -> numeric[tensorEntry[H5, {1, 1, 1, 1, 1}]],
    "H5_x_x_x_y_py" -> numeric[tensorEntry[H5, {1, 1, 1, 3, 4}]],
    "H5_x_px_px_px_px" -> numeric[tensorEntry[H5, {1, 2, 2, 2, 2}]],
    "H6_x_x_x_x_x_x" -> numeric[tensorEntry[H6, {1, 1, 1, 1, 1, 1}]],
    "H6_x_x_x_x_y_y" -> numeric[tensorEntry[H6, {1, 1, 1, 1, 3, 3}]],
    "H6_x_x_x_x_y_py" -> numeric[tensorEntry[H6, {1, 1, 1, 1, 3, 4}]],
    "H6_y_y_y_y_py_py" -> numeric[tensorEntry[H6, {3, 3, 3, 3, 4, 4}]],
    "H6_px_px_px_px_px_px" -> numeric[tensorEntry[H6, {2, 2, 2, 2, 2, 2}]],
    "H6_delta6" -> numeric[tensorEntry[H6, {6, 6, 6, 6, 6, 6}]]
  |>,
  "linear_map" -> <|
    "matrix" -> Map[numeric, linearMap, {2}],
    "symplectic_residual" -> numeric[symplecticResidual]
  |>,
  "variational_map_contractions" -> <|
    "probe_u" -> numeric /@ probeU,
    "probe_v" -> numeric /@ probeV,
    "Ru" -> (numeric /@ jetU[[1]]),
    "Tuu" -> (numeric /@ jetU[[2]]),
    "Uuuu" -> (numeric /@ jetU[[3]]),
    "Vuuuu" -> (numeric /@ jetU[[4]]),
    "Wuuuuu" -> (numeric /@ jetU[[5]]),
    "Rv" -> (numeric /@ jetV[[1]]),
    "Tvv" -> (numeric /@ jetV[[2]]),
    "Uvvv" -> (numeric /@ jetV[[3]]),
    "Vvvvv" -> (numeric /@ jetV[[4]]),
    "Wvvvvv" -> (numeric /@ jetV[[5]])
  |>,
  "dragt_finn_fifth_order_sample" -> <|
    "f3" -> <|"q2p" -> numeric[2/5 / 2]|>,
    "f4" -> <|"q4" -> numeric[-3/7 / 4], "q2p2" -> numeric[5/9 / 2]|>,
    "f5" -> <|"q5" -> numeric[5/11 / 5], "q3p2" -> numeric[-4/13 / 3]|>,
    "f6" -> <|"q6" -> numeric[7/13 / 6], "q4p2" -> numeric[3/17 / 4]|>,
    "composition" ->
      "exp(:f3:) after exp(:f4:) after exp(:f5:) after exp(:f6:)",
    "q_map_coefficients" -> mapCoefficients[qMap],
    "p_map_coefficients" -> mapCoefficients[pMap]
  |>
|>;

output = FileNameJoin[{DirectoryName[$InputFileName],
  "lie_map_symbolic_reference_p6.json"}];
Export[output, data, "RawJSON"];
Print["wrote " <> output];
