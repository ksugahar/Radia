#!/usr/bin/env wolframscript
(* ::Package:: *)

(*
  derive_quad_q1_henrotte.wls

  Derive closed-form 4x4 element matrices for the FEMM/Henrotte axisymmetric
  Q1 element on an axis-aligned rectangle in (r, z).

  Basis (Q1 in s = r^2 reference coords):
      phi(r, z) = a + b*s + c*z + d*s*z      where s = r^2

  Nodes: (s_a, z_a), (s_b, z_a), (s_b, z_b), (s_a, z_b)  with s_a = r_a^2, s_b = r_b^2

  Outputs:
    - Mz: stiffness from B_z component                  (b, d coefficients)
    - Mr: stiffness from B_r component (with log terms) (c, d coefficients)
    - M_sigma: sigma mass matrix in V coords            (full 4x4)

  Energy form:
    W_z = (1/(2 mu_z)) Integrate( B_z^2, 2 pi r dr dz )    with B_z = (b + d z)/pi
    W_r = (1/(2 mu_r)) Integrate( B_r^2, 2 pi r dr dz )    with B_r = -(c + d s)/(2 pi sqrt(s))

  Useful change of variables: 2 pi r dr dz = pi ds dz
*)

Print["=== Deriving FEMM/Henrotte Q1 axis-aligned rectangle element matrices ==="];
Print[""];

(* Symbols *)
ClearAll[a, b, c, d, sa, sb, za, zb, ra, rb, muR, muZ, sigma];

(* Phi interpolation: linear in (s, z), bilinear including s*z term *)
phi[s_, z_, a_, b_, c_, d_] := a + b s + c z + d s z;

(* Nodal values *)
phi0 = phi[sa, za, a, b, c, d];  (* node 0: (sa, za) *)
phi1 = phi[sb, za, a, b, c, d];  (* node 1: (sb, za) *)
phi2 = phi[sb, zb, a, b, c, d];  (* node 2: (sb, zb) *)
phi3 = phi[sa, zb, a, b, c, d];  (* node 3: (sa, zb) *)

(* Solve coefficients (a,b,c,d) in terms of nodal values phi0..phi3 *)
sol = Solve[{
    p0 == phi0, p1 == phi1, p2 == phi2, p3 == phi3
  }, {a, b, c, d}][[1]];

aSol = a /. sol;
bSol = b /. sol;
cSol = c /. sol;
dSol = d /. sol;

Print["a = ", FullSimplify[aSol]];
Print["b = ", FullSimplify[bSol]];
Print["c = ", FullSimplify[cSol]];
Print["d = ", FullSimplify[dSol]];
Print[""];

(* Field components *)
(* B_z = (1/(2 pi r)) d phi/dr = (1/(2 pi r)) * 2 r * d phi/ds = d phi/ds / pi *)
(* For phi = a + b s + c z + d s z, d phi/ds = b + d z, so B_z = (b + d z)/pi *)
(* B_r = -(1/(2 pi r)) d phi/dz = -(c + d s)/(2 pi sqrt(s)) *)

Bz[s_, z_] := (b + d z)/Pi;
Br[s_, z_] := -(c + d s)/(2 Pi Sqrt[s]);

(* Energy: W = (1/2)(B_z^2/mu_z + B_r^2/mu_r) integrated over 2 pi r dr dz = pi ds dz *)
(* Integration limits: s in [sa, sb], z in [za, zb] *)

(* W_z integrand: (B_z^2 / mu_z) * pi ds dz / 2 *)
(*    = (1/(2 mu_z))((b + d z)/pi)^2 * pi ds dz *)
(*    = (b + d z)^2 / (2 pi mu_z) ds dz *)

Wz = Integrate[
  (b + d z)^2 / (2 Pi muZ),
  {s, sa, sb}, {z, za, zb}
];
Wz = FullSimplify[Wz];
Print["W_z (energy from B_z), in (b, d) coefficients:"];
Print[Wz];
Print[""];

(* W_r integrand: (B_r^2 / mu_r) * pi ds dz / 2 *)
(*    = (1/(2 mu_r)) * (c + d s)^2 / (4 pi^2 s) * pi ds dz *)
(*    = (c + d s)^2 / (8 pi mu_r s) ds dz *)

Wr = Integrate[
  (c + d s)^2 / (8 Pi muR s),
  {s, sa, sb}, {z, za, zb}
];
Wr = FullSimplify[Wr];
Print["W_r (energy from B_r), in (c, d) coefficients:"];
Print[Wr];
Print[""];

(* Total energy in coefficients *)
W = Wz + Wr;

(* Substitute coefficients back to (phi0, phi1, phi2, phi3) *)
W_phi = W /. {a -> aSol, b -> bSol, c -> cSol, d -> dSol};
W_phi = FullSimplify[W_phi];

(* Differentiate to get element matrix in phi coords *)
phiVars = {p0, p1, p2, p3};
M_phi = Table[D[D[W_phi, phiVars[[i]]], phiVars[[j]]], {i, 1, 4}, {j, 1, 4}];
M_phi = FullSimplify[M_phi];

Print["Element matrix M_phi (in phi nodal value coords):"];
Print[MatrixForm[M_phi]];
Print[""];

(* Convert to V coords: phi_j = 2 pi r_j V_j  =>  V_j = phi_j / (2 pi r_j)  *)
(* M_V = T . M_phi . T  where  T_jj = 2 pi r_j *)
T = DiagonalMatrix[{2 Pi ra, 2 Pi rb, 2 Pi rb, 2 Pi ra}];
M_V = T . M_phi . T;
M_V = FullSimplify[M_V];

Print["Element matrix M_V (in V = A_phi nodal coords):"];
Print[MatrixForm[M_V]];
Print[""];

(* Substitute sa = ra^2, sb = rb^2 *)
M_V_r = M_V /. {sa -> ra^2, sb -> rb^2};
M_V_r = FullSimplify[M_V_r, Assumptions -> {ra > 0, rb > ra, zb > za}];
Print["Element matrix M_V (with sa=ra^2, sb=rb^2):"];
Print[MatrixForm[M_V_r]];
Print[""];

(* Sigma mass matrix in V coords *)
(* W_sigma = (sigma/2) Integrate( A_phi^2, 2 pi r dr dz ) *)
(*         = (sigma/2) Integrate( phi^2 / (2 pi r)^2, 2 pi r dr dz ) *)
(*         = (sigma/(8 pi^2)) Integrate( phi^2 / r * 2 dr dz ) *)
(*         = (sigma/(4 pi^2)) Integrate( phi^2 / r dr dz ) *)
(*  In s coords: dr = ds/(2 sqrt(s)), so dr/r = ds/(2 s) *)
(*         = (sigma/(8 pi^2)) Integrate( phi^2 / s ds dz ) *)
(*  Wait: phi^2 / r * dr * dz = phi^2 * (1/r) * (ds/(2 r)) * dz = phi^2 * ds*dz/(2 r^2) = phi^2 * ds*dz/(2 s) *)
(*  So W_sigma = (sigma/(4 pi^2)) Integrate[phi^2 / (2 s) ds dz] = (sigma/(8 pi^2)) Integrate[phi^2/s ds dz] *)

Print["Computing sigma mass matrix..."];
phiSym[s_, z_] := a + b s + c z + d s z;
WsigCoef = Integrate[
  sigma * phiSym[s, z]^2 / (8 Pi^2 s),
  {s, sa, sb}, {z, za, zb}
];
WsigCoef = FullSimplify[WsigCoef];
WsigPhi = WsigCoef /. {a -> aSol, b -> bSol, c -> cSol, d -> dSol};
WsigPhi = FullSimplify[WsigPhi];

MsigPhi = Table[D[D[WsigPhi, phiVars[[i]]], phiVars[[j]]], {i, 1, 4}, {j, 1, 4}];
MsigPhi = FullSimplify[MsigPhi];

MsigV = T . MsigPhi . T;
MsigV = FullSimplify[MsigV];
MsigV_r = MsigV /. {sa -> ra^2, sb -> rb^2};
MsigV_r = FullSimplify[MsigV_r, Assumptions -> {ra > 0, rb > ra, zb > za, sigma > 0, muR > 0, muZ > 0}];

Print["Sigma mass matrix M_sigma_V (rb > ra > 0):"];
Print[MatrixForm[MsigV_r]];
Print[""];

(* Numerical sanity check: pick concrete values, evaluate, verify symmetry/PSD *)
testRules = {ra -> 1*^-3, rb -> 2*^-3, za -> 0, zb -> 1*^-3,
             muR -> 4 Pi*^-7, muZ -> 4 Pi*^-7, sigma -> 5.8*^7};

K_test = N[M_V_r /. testRules];
M_test = N[MsigV_r /. testRules];

Print["Numerical test (ra=1mm, rb=2mm, za=0, zb=1mm, mu=mu0, sigma=Cu):"];
Print["K_test = ", MatrixForm[Chop[K_test]]];
Print["M_test = ", MatrixForm[Chop[M_test]]];
Print["K symmetric? ", Norm[K_test - Transpose[K_test]] < 1*^-15];
Print["M symmetric? ", Norm[M_test - Transpose[M_test]] < 1*^-15];
Print["K eigenvalues: ", Eigenvalues[K_test]];
Print["M eigenvalues: ", Eigenvalues[M_test]];

(* Export results: structured association *)
result = <|
  "M_V_general" -> ToString[M_V_r, InputForm],
  "M_sigma_V_general" -> ToString[MsigV_r, InputForm],
  "K_test_numerical" -> Chop[K_test],
  "M_test_numerical" -> Chop[M_test],
  "test_params" -> testRules
|>;
Export["W:/30_CauerLadderNetwork/2026_04_01_\:9577\:65b9\:5f62CLN/axifem/quad_q1_henrotte_matrices.json", result, "JSON"];
Print[""];
Print["Saved: quad_q1_henrotte_matrices.json"];
