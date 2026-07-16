#!/usr/bin/env wolframscript
(* ::Package:: *)

(*
  evrs_derham_bridge.wls -- symbolic FEEC bridge for the EVRS / EID eddy-current
  reduction idea.

  The research construction starts from a high-order HCurl parent space, keeps a
  small Eddy-Visible Response Space (EVRS), and uses T as the HCurl unknown:

      T in H(curl),      J = curl T.

  This file verifies the de Rham facts that make that reduction structurally
  safe:

      H1 --grad--> HCurl --curl--> HDiv --div--> L2,
      curl(grad phi) = 0,
      div(curl T) = 0,

  The second half of this script makes the corresponding DISCRETE T-method
  algebra explicit:

      phi_coeff --G--> T_coeff --C--> J_coeff --D--> rho_coeff,
      C G = 0,   D C = 0,
      T_coeff = Q z,   J_coeff = C Q z,
      R_T = C^T M_R C,   L_T = C^T M_L C.

  Therefore any EVRS projection Q preserves solenoidal current, gauge-null
  gradients, symmetry, and the passive Galerkin quadratic-form structure.

  For the SIBC/surface-Omega branch,

      K = n x grad_Gamma Omega,       div_Gamma K = 0

  on a flat reference face.  The checks are symbolic and intentionally small:
  they do not test numerical accuracy of a solver; they lock the topology and
  calculus identities behind the EVRS / eddy-bubble terminology.
*)

lam = {1 - x - y - z, x, y, z};
grad[f_] := D[f, {{x, y, z}}];
crl[F_] := Curl[F, {x, y, z}];
dvg[F_] := Div[F, {x, y, z}];
gl = grad /@ lam;

edges = {{1, 2}, {1, 3}, {1, 4}, {2, 3}, {2, 4}, {3, 4}};
faces = {{2, 3, 4}, {1, 3, 4}, {1, 2, 4}, {1, 2, 3}};

W1[{a_, b_}] := lam[[a]] gl[[b]] - lam[[b]] gl[[a]];
W2[{a_, b_, c_}] := 2 (lam[[a]] Cross[gl[[b]], gl[[c]]] +
                       lam[[b]] Cross[gl[[c]], gl[[a]]] +
                       lam[[c]] Cross[gl[[a]], gl[[b]]]);

evalRow[f_, pts_] := Flatten[Table[f /. pt, {pt, pts}]];
valsMat[funcs_, pts_] := evalRow[#, pts] & /@ funcs;
zeroMatrix[m_, n_] := ConstantArray[0, {m, n}];
basisPts = {
  {x -> 0, y -> 0, z -> 0},
  {x -> 1, y -> 0, z -> 0},
  {x -> 0, y -> 1, z -> 0},
  {x -> 0, y -> 0, z -> 1},
  {x -> 1/5, y -> 1/7, z -> 1/11},
  {x -> 1/3, y -> 1/5, z -> 1/7}
};
symMatrix[n_, prefix_] := Table[
  If[i <= j,
    Symbol[prefix <> ToString[i] <> ToString[j]],
    Symbol[prefix <> ToString[j] <> ToString[i]]
  ],
  {i, n}, {j, n}
];

solveInBasis[vec_, basis_] := Module[{cc, eqs, sol},
  cc = Simplify[LeastSquares[Transpose[valsMat[basis, basisPts]], evalRow[vec, basisPts]]];
  If[
    Simplify[vec - Total[Table[cc[[i]] basis[[i]], {i, Length[basis]}]]] =!= {0, 0, 0},
    Print["  FAIL basis solve residual for ", vec];
    Exit[1]
  ];
  cc
];

(* Flat surface Gamma = {z=0}, n=e_z.  Surface-Omega current is n x grad_Gamma Omega. *)
surfaceK[omega_] := {-D[omega, y], D[omega, x], 0};
surfaceDiv[K_] := D[K[[1]], x] + D[K[[2]], y];

Module[
  {ok = True, chk, w1, w2, coeff, T, J, phiCoeff, phi, omega,
   pts, curls, rHDiv, rWithCurls, Gmat, Cmat, Dmat, qMat, zVec, Tcoeff,
   Jcoeff, mr, ml, rT, lT, rEVRS, lEVRS, pMat, bT, bEVRS, quadR, quadREvrs},

  chk[name_, cond_] := (
    Print["  ", If[cond, "ok  ", "FAIL"], name];
    ok = ok && cond
  );

  Print["evrs_derham_bridge.wls self-test (EVRS/EID de Rham bridge):"];

  w1 = W1 /@ edges;
  w2 = W2 /@ faces;

  (* Gauge / gradient directions do not create current.  This is why HCurl gauge
     cleanup and nograds/tree-cotree choices are compatible with the EVRS branch. *)
  phiCoeff = Array[a, 4];
  phi = phiCoeff.lam;
  chk["curl(grad phi) = 0 for arbitrary H1 scalar potential",
      Simplify[crl[grad[phi]]] === {0, 0, 0}];

  (* EVRS is an arbitrary linear subspace of HCurl.  Therefore each retained
     coordinate T_i still maps to a solenoidal HDiv-compatible current J_i. *)
  coeff = Array[q, Length[w1]];
  T = Total[Table[coeff[[i]] w1[[i]], {i, Length[w1]}]];
  J = Simplify[crl[T]];
  chk["div(curl T_EVRS) = 0 for arbitrary retained HCurl combination",
      Simplify[dvg[J]] === 0];

  (* The image of curl(W1) lands in the Whitney W2 / RT0 HDiv space. *)
  curls = Simplify[crl /@ w1];
  pts = {
    {x -> 0, y -> 0, z -> 0},
    {x -> 1, y -> 0, z -> 0},
    {x -> 0, y -> 1, z -> 0},
    {x -> 0, y -> 0, z -> 1},
    {x -> 1/5, y -> 1/7, z -> 1/11},
    {x -> 1/3, y -> 1/5, z -> 1/7}
  };
  rHDiv = MatrixRank[valsMat[w2, pts]];
  rWithCurls = MatrixRank[valsMat[Join[w2, curls], pts]];
  chk["curl(Whitney HCurl) lies in Whitney HDiv / RT0 span",
      rWithCurls == rHDiv];

  chk["each EVRS current basis vector is divergence-free before any VIM projection",
      And @@ (Simplify[dvg[#]] === 0 & /@ curls)];

  (* ----- Discrete T-method matrices on the Whitney tet. -----
     G maps H1 nodal coefficients to HCurl edge coefficients.
     C maps HCurl edge coefficients of T to HDiv face coefficients of J=curl T.
     D maps HDiv face coefficients to cell L2 divergence.
  *)
  Gmat = Transpose[Table[solveInBasis[grad[lam[[i]]], w1], {i, Length[lam]}]];
  Cmat = Transpose[Table[solveInBasis[curls[[i]], w2], {i, Length[w1]}]];
  Dmat = {Simplify[dvg /@ w2]};

  chk["discrete de Rham chain uses variable T: C.G = 0 (curl grad phi = 0)",
      Simplify[Cmat.Gmat] === zeroMatrix[Length[w2], Length[lam]]];
  chk["discrete de Rham chain: D.C = 0 (div curl T = 0)",
      Simplify[Dmat.Cmat] === zeroMatrix[1, Length[w1]]];

  (* Arbitrary EVRS projection.  Q is not assumed special: it can be any retained
     response-visible subspace, so the identity proves EID elimination cannot
     break charge conservation. *)
  qMat = Array[q, {Length[w1], 3}];
  zVec = Array[z, 3];
  Tcoeff = qMat.zVec;
  Jcoeff = Simplify[Cmat.Tcoeff];
  chk["EVRS T = Q z remains solenoidal: D.C.Q.z = 0 for arbitrary Q",
      Simplify[Dmat.Jcoeff] === {0}];

  (* T-method curl-only matrices.  mr/ml are generic symmetric current-space
     matrices standing in for resistance, inductance, BEM/VIM, or CLN moments.
     The proofs are algebraic, not numerical. *)
  mr = symMatrix[Length[w2], "r"];
  ml = symMatrix[Length[w2], "l"];
  rT = Simplify[Transpose[Cmat].mr.Cmat];
  lT = Simplify[Transpose[Cmat].ml.Cmat];
  chk["R_T = C^T M_R C and L_T = C^T M_L C are symmetric for symmetric current matrices",
      Simplify[rT - Transpose[rT]] === zeroMatrix[Length[w1], Length[w1]] &&
      Simplify[lT - Transpose[lT]] === zeroMatrix[Length[w1], Length[w1]]];
  chk["curl-only T matrices annihilate gradient gauge directions",
      Simplify[rT.Gmat] === zeroMatrix[Length[w1], Length[lam]] &&
      Simplify[lT.Gmat] === zeroMatrix[Length[w1], Length[lam]]];

  rEVRS = Simplify[Transpose[qMat].rT.qMat];
  lEVRS = Simplify[Transpose[qMat].lT.qMat];
  quadR = Simplify[Tcoeff.rT.Tcoeff - Jcoeff.mr.Jcoeff];
  quadREvrs = Simplify[zVec.rEVRS.zVec - Jcoeff.mr.Jcoeff];
  chk["EVRS Galerkin resistance identity: (Q z)^T R_T (Q z) = J^T M_R J",
      quadR === 0 && quadREvrs === 0];
  chk["EVRS projected matrices stay symmetric",
      Simplify[rEVRS - Transpose[rEVRS]] === zeroMatrix[3, 3] &&
      Simplify[lEVRS - Transpose[lEVRS]] === zeroMatrix[3, 3]];

  (* Port forcing written in current space.  A gradient-only T has zero current,
     so it sees no port RHS.  This is the algebraic port version of the gauge
     statement and the correct shape for ReducedPortAdmittance. *)
  pMat = Array[p, {Length[w2], 2}];
  bT = Simplify[Transpose[Cmat].pMat];
  bEVRS = Simplify[Transpose[qMat].bT];
  chk["current-space port RHS is gauge-invariant: G^T C^T P = 0",
      Simplify[Transpose[Gmat].bT] === zeroMatrix[Length[lam], 2]];
  chk["EVRS port RHS has reduced size nEVRS x nPorts",
      Dimensions[bEVRS] === {3, 2}];

  (* Surface-Omega / SIBC branch: the rotated surface gradient is tangential and
     surface-solenoidal on a flat face.  This is the surface de Rham analogue of
     J=curl T. *)
  omega = b0 + b1 x + b2 y + b3 x^2 + b4 x y + b5 y^2;
  chk["surface-Omega K = n x grad_Gamma Omega is tangential",
      Simplify[surfaceK[omega][[3]]] === 0];
  chk["div_Gamma(n x grad_Gamma Omega) = 0 for arbitrary quadratic Omega",
      Simplify[surfaceDiv[surfaceK[omega]]] === 0];

  Print[If[ok, "ALL PASS", "SOME FAILED"]];
  If[! ok, Exit[1]];
];
