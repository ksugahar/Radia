(* ::Package:: *)

(* hex_cuboid_eigenvalue.wls
   (b) Build curl-curl + mass matrices for cuboid 5x2x1 mm using FAST monomial
   integration via CoefficientRules. Solve generalized eigenvalue problem,
   compare to Phase 7-11 reference 19.77 us.

   Key optimization: replace symbolic Integrate (sec/entry) with fast monomial
   integration (ms/entry) via CoefficientRules.

   Date: 2026-05-02
*)

Get["w:\\30_CauerLadderNetwork\\2026_04_01_\:9577\:65b9\:5f62CLN\\ngsolve_validation\\schoberl_zaglmayr_basis.wls"];

ClearAll[x, y, z];
$MaxExtraPrecision = 100;

(* === Fast monomial integration over [0,1]^3 === *)
fastIntegrate3D[poly_] := Module[{p = Expand[poly], rules},
  If[p === 0, Return[0]];
  rules = CoefficientRules[p, {x, y, z}];
  Total[
    Function[r, r[[2]] / ((r[[1, 1]] + 1) (r[[1, 2]] + 1) (r[[1, 3]] + 1))]
      /@ rules
  ]
];

(* === Build component-wise matrices fast === *)
componentMatrixSym[shapes_, slot_, comp_] := Module[{n = Length[shapes], M, exprs},
  exprs = Table[Expand[shapes[[i, slot, comp]]], {i, 1, n}];
  M = Table[
    If[j < i, 0, fastIntegrate3D[exprs[[i]] exprs[[j]]]],
    {i, 1, n}, {j, 1, n}
  ];
  M + Transpose[M] - DiagonalMatrix[Diagonal[M]]
];

ord = 1;
Print["Building HEX H(curl) reference component matrices via fast integration"];
Print["  order = ", ord];
shapes = SchoberlZaglmayrBasis`ShapesHEX[ord, Range[0, 7], x, y, z];
n = Length[shapes];
Print["  N basis = ", n];

t0 = AbsoluteTime[];
MrefX = componentMatrixSym[shapes, 1, 1];   (* slot 1 = value, comp 1 = x *)
Print["  MrefX done in ", Round[AbsoluteTime[] - t0, 0.1], " s"];
t0 = AbsoluteTime[];
MrefY = componentMatrixSym[shapes, 1, 2];
Print["  MrefY done in ", Round[AbsoluteTime[] - t0, 0.1], " s"];
t0 = AbsoluteTime[];
MrefZ = componentMatrixSym[shapes, 1, 3];
Print["  MrefZ done in ", Round[AbsoluteTime[] - t0, 0.1], " s"];

t0 = AbsoluteTime[];
SrefX = componentMatrixSym[shapes, 2, 1];   (* slot 2 = curl, comp 1 = curl_x *)
SrefY = componentMatrixSym[shapes, 2, 2];
SrefZ = componentMatrixSym[shapes, 2, 3];
Print["  Sref XYZ done in ", Round[AbsoluteTime[] - t0, 0.1], " s"];

(* === Physical cuboid 5 x 2 x 1 mm === *)
sigma = SetPrecision[58 10^6, 30];
mu0   = 4 Pi SetPrecision[10^-7, 30];
{aBox, bBox, cBox} = SetPrecision[{5 10^-3, 2 10^-3, 1 10^-3}, 30];

(* Covariant H(curl) Piola transform for axis-aligned scaling J = diag(a,b,c):
     phi_phys = J^{-T} phi_ref = (phi_x/a, phi_y/b, phi_z/c)
     curl_phys = J/det(J) curl_ref = (curl_x/(bc), curl_y/(ac), curl_z/(ab))

   Mass matrix entry:
     M_phys[i,j] = (bc/a) Mref_x + (ac/b) Mref_y + (ab/c) Mref_z

   Stiffness matrix entry:
     S_phys[i,j] = (a/(bc)) Sref_x + (b/(ac)) Sref_y + (c/(ab)) Sref_z
*)
abc = aBox bBox cBox;
Mphys = N[(bBox cBox / aBox) MrefX + (aBox cBox / bBox) MrefY + (aBox bBox / cBox) MrefZ, 30];
Sphys = N[(aBox / (bBox cBox)) SrefX + (bBox / (aBox cBox)) SrefY + (cBox / (aBox bBox)) SrefZ, 30];

Print[];
Print["Physical cuboid: ", N[aBox 10^3, 4], " x ", N[bBox 10^3, 4], " x ", N[cBox 10^3, 4], " mm"];

(* === Generalized eigenvalue solve === *)
Print[];
Print["Solving generalized eigenvalue problem K v = lambda M v ..."];
t0 = AbsoluteTime[];
{eigvals, eigvecs} = Eigensystem[{Sphys, Mphys}];
Print["  Solved in ", Round[AbsoluteTime[] - t0, 0.1], " s"];

(* Filter physical eigenvalues:
   - Pure gradient modes have curl = 0, giving lambda = 0
   - Real positive lambda are eddy-current eigenvalues
*)
realEigs = Re[eigvals];
nonzeroEigs = Select[realEigs, # > 10^-6 &];
sortedEigs = Sort[nonzeroEigs];

Print[];
Print["Number of nonzero eigenvalues: ", Length[nonzeroEigs], " / ", Length[eigvals]];
Print["First 10 nonzero eigenvalues lambda (ascending):"];
Do[Print["  lambda_", i, " = ", N[sortedEigs[[i]], 12]], {i, 1, Min[10, Length[sortedEigs]]}];

Print[];
Print["Time constants tau_n = 1/(sigma mu0 lambda_n) [us]:"];
Do[
  Module[{tau = sigma mu0 / sortedEigs[[i]]},
    Print["  tau_", i, " = ", N[10^6 tau, 12], " us"]
  ],
  {i, 1, Min[10, Length[sortedEigs]]}
];

Print[];
Print["============================================"];
Print["Phase 7-11 spectral basis: tau_1 = 19.77 us"];
Print["This run leading tau:       ",
      N[10^6 sigma mu0 / sortedEigs[[1]], 8], " us"];
Print["============================================"];
