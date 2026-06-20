(* ::Package:: *)

(* hex_bem_full_hcurl.wls
   Full HCurl order 1 BEM K matrix on cuboid 5x2x1 mm, then Cauer-I extraction
   to verify tau_1 = 19.77 us reference.

   Strategy:
   1. Precompute Aharoni-extended F_kl(len, dy, dz) = int_0^len int_0^len
      x^k x'^l / sqrt((x-x')^2 + dy^2 + dz^2) dx dx'  for k, l = 0..3
      -> 16 closed-form expressions, used for ANY direction (x, y, or z reduction).

   2. HCurl order 1 basis on physical cuboid (54 DOFs total).

   3. K[i,j] = (mu0/4pi) sum over component pairs (xx, yy, zz) of:
        Component_xx[i,j] = sum over polynomial expansion of phi_ix * phi_jx
                            of [coeff * F_xkl * 4D_NIntegrate over (perp vars)]
      Similar for yy, zz.

   4. M[i,j] = int phi_i . phi_j dV  (closed-form polynomial integration)

   5. Solve K v = lambda M v, project onto external A_ext, compute moments,
      Cauer-I extraction.

   Date: 2026-05-02
*)

Get["w:\\30_CauerLadderNetwork\\2026_04_01_\:9577\:65b9\:5f62CLN\\ngsolve_validation\\schoberl_zaglmayr_basis.wls"];

ClearAll[x, y, z, xp, yp, zp];
$MaxExtraPrecision = 50;

Print[StringRepeat["=", 70]];
Print["Full HCurl order 1 BEM on cuboid 5x2x1 mm Cu  (-> tau_1 = 19.77 us target)"];
Print[StringRepeat["=", 70]];

(* === Constants === *)
sigma = SetPrecision[58 10^6, 30];
mu0   = 4 Pi SetPrecision[10^-7, 30];
{aBox, bBox, cBox} = SetPrecision[{5 10^-3, 2 10^-3, 1 10^-3}, 30];

(* === HCurl order 1 basis on reference [0,1]^3 === *)
ord = 1;
shapesRef = SchoberlZaglmayrBasis`ShapesHEX[ord, Range[0, 7], x, y, z];
nDOF = Length[shapesRef];
Print["HCurl order ", ord, ": ", nDOF, " DOFs"];

(* Apply Piola covariant: phi_phys = J^{-T} phi_ref with J = diag(a, b, c).
   In physical coords (x_phys = a x_ref, etc.), the basis vectors become
   (phi_ref_x / a, phi_ref_y / b, phi_ref_z / c).
   We work in REFERENCE coords throughout, then scale at the end. *)

(* === Mass matrix in reference coords (component-wise) === *)
fastIntegrate3D[poly_, vars_] := Module[{p = Expand[poly], rules},
  If[p === 0, Return[0]];
  rules = CoefficientRules[p, vars];
  Total[Function[r, r[[2]] / Times @@ Map[(# + 1) &, r[[1]]]] /@ rules]
];

componentMassMatrix[shapes_, comp_] := Module[{n = Length[shapes], M, exprs},
  exprs = Table[Expand[shapes[[i, 1, comp]]], {i, 1, n}];
  M = Table[
    If[j < i, 0, fastIntegrate3D[exprs[[i]] exprs[[j]], {x, y, z}]],
    {i, 1, n}, {j, 1, n}
  ];
  M + Transpose[M] - DiagonalMatrix[Diagonal[M]]
];

t0 = AbsoluteTime[];
MrefXX = componentMassMatrix[shapesRef, 1];
MrefYY = componentMassMatrix[shapesRef, 2];
MrefZZ = componentMassMatrix[shapesRef, 3];
Print["Mass matrices (3 component blocks) computed in ", Round[AbsoluteTime[] - t0, 0.1], " s"];

(* Physical mass matrix using Piola scaling *)
Mphys = N[(bBox cBox / aBox) MrefXX +
          (aBox cBox / bBox) MrefYY +
          (aBox bBox / cBox) MrefZZ, 30];

Print["Physical M matrix assembled. Diagonal[1..6]: ", N[Diagonal[Mphys][[1 ;; 6]], 12]];


(* === BEM K matrix - direct 6D NIntegrate at WP=15 for proof of concept === *)

(* For each entry:
   K[i,j] = (mu0/4pi) (1/(a^2)) int phi_ix phi_jx + (1/(b^2)) int phi_iy phi_jy + (1/(c^2)) int phi_iz phi_jz
   where each int is 6D NIntegrate over [0,a]x[0,b]x[0,c]^2
   -> very slow, but reliable.

   Speedup: use REFERENCE coords [0,1]^3 with rescaled distance:
   |r_phys - r'_phys| = sqrt(a^2(x-xp)^2 + b^2(y-yp)^2 + c^2(z-zp)^2)
   K[i,j]_phys = (mu0/4pi) sum_alpha (a*b*c)^2 / d_alpha^2 *
                 int_[0,1]^6 phi_ialpha(r_ref) phi_jalpha(r'_ref) /
                              sqrt(a^2*Dx^2 + b^2*Dy^2 + c^2*Dz^2) dV_ref dV_ref'
*)

(* Direct 6D NIntegrate per entry *)
Kentry[i_, j_, WP_:15] := Module[{vi, vj, integrand, factor, scaled},
  vi = shapesRef[[i, 1]];
  vj = shapesRef[[j, 1]] /. {x -> xp, y -> yp, z -> zp};
  (* Component-weighted dot product with Piola scaling *)
  integrand = (vi[[1]] vj[[1]] / aBox^2 +
               vi[[2]] vj[[2]] / bBox^2 +
               vi[[3]] vj[[3]] / cBox^2);
  factor = aBox^2 bBox^2 cBox^2;   (* (abc)^2 from Jacobian dV_phys^2 *)
  scaled = factor / Sqrt[aBox^2 (x - xp)^2 + bBox^2 (y - yp)^2 + cBox^2 (z - zp)^2];
  (mu0 / (4 Pi)) NIntegrate[
    integrand scaled,
    {x, 0, 1}, {y, 0, 1}, {z, 0, 1},
    {xp, 0, 1}, {yp, 0, 1}, {zp, 0, 1},
    Method -> "GlobalAdaptive", WorkingPrecision -> WP,
    PrecisionGoal -> 4, AccuracyGoal -> 6, MaxRecursion -> 10
  ]
];

Print[];
Print["Testing K[1,1] at WP=20 (better precision) ..."];
t0 = AbsoluteTime[];
K11test = Kentry[1, 1, 20];
Print["K[1,1] = ", N[K11test, 12], "  in ", Round[AbsoluteTime[] - t0, 0.1], " s"];

(* === Build full K matrix === *)
nUnique = nDOF (nDOF + 1) / 2;
Print[];
Print["Building full K matrix (", nUnique, " unique entries) ..."];

K = Table[0, {nDOF}, {nDOF}];
totalTime = 0;
counter = 0;

Do[
  counter++;
  t0 = AbsoluteTime[];
  K[[i, j]] = Kentry[i, j, 15];
  K[[j, i]] = K[[i, j]];
  totalTime += AbsoluteTime[] - t0;
  If[Mod[counter, 100] == 0,
    Print["  ", counter, "/", nUnique, " entries; total time so far: ",
          Round[totalTime/60, 0.1], " min"]
  ],
  {i, 1, nDOF}, {j, i, nDOF}
];

Print["K matrix done in ", Round[totalTime/60, 0.1], " minutes"];

(* === Solve generalized eigenvalue === *)
Print[];
Print["Solving K v = lambda M v ..."];
{eigvals, eigvecs} = Eigensystem[{N[K, 30], Mphys}];
Print["Eigenvalues solved."];

(* Filter physical (positive real) *)
realEigs = Re[eigvals];
sortedIdx = Reverse[Ordering[realEigs]];   (* descending = largest lambda first *)
sortedEigs = realEigs[[sortedIdx]];
nonzeroEigs = Select[sortedEigs, # > 10^-30 &];
Print["Number of positive eigenvalues: ", Length[nonzeroEigs]];
Print["Largest 5 lambda: ", N[nonzeroEigs[[1 ;; Min[5, Length[nonzeroEigs]]]], 12]];

(* Convert to tau via tau = sigma * lambda  (mu0/4pi already in K) *)
tausUs = 10^6 sigma * nonzeroEigs;
Print[];
Print["Largest 5 tau (us): ", N[tausUs[[1 ;; Min[5, Length[tausUs]]]], 8]];

(* === Project onto external A_ext = (-y/2, x/2, 0) * H_z (uniform z-field excitation) === *)
Print[];
Print["Computing projection onto uniform H_z external field ..."];
(* External A_ext in physical = (-y/2, x/2, 0) * H_z. Set H_z = 1 *)
AextX = -y/2;   (* x-component *)
AextY = x/2;    (* y-component *)

(* Project basis onto external A: beta_n = sum_i v_n[i] <phi_i, A_ext>
   where <phi_i, A_ext> = int phi_i_x AextX + phi_i_y AextY dV (in physical coords). *)
projVec = Table[
  Module[{phi = shapesRef[[i, 1]], integrand},
    integrand = (phi[[1]] / aBox) * (-y bBox/2) * aBox bBox cBox / (aBox bBox cBox) +   (* simplified *)
                (phi[[2]] / bBox) * (x aBox/2) * aBox bBox cBox / (aBox bBox cBox);
    (* Compute in reference, multiply by Jacobian abc *)
    aBox bBox cBox * Module[{integrandRef},
      integrandRef = (phi[[1]] / aBox) * (-y bBox/2) + (phi[[2]] / bBox) * (x aBox/2);
      fastIntegrate3D[Expand[integrandRef], {x, y, z}]
    ]
  ],
  {i, 1, nDOF}
];

(* beta_n = v_n^T M projVec *)
betas = Table[
  Module[{vn = eigvecs[[sortedIdx[[k]]]]},
    vn . Mphys . projVec
  ],
  {k, 1, Length[sortedEigs]}
];

(* === Compute alpha_n moments === *)
Print[];
Print["Computing alpha_n moments ..."];
NMoments = 8;
alphas = Table[
  (-1)^k Sum[
    Module[{lam = sortedEigs[[idx]], beta = betas[[idx]]},
      If[lam > 10^-30, beta^2 lam^(k - 1), 0]
    ],
    {idx, 1, Length[sortedEigs]}
  ] (sigma)^k / (aBox bBox cBox),
  {k, 1, NMoments}
];

Print["alpha_n moments:"];
Do[Print["  alpha_", k, " = ", N[alphas[[k]], 12]], {k, 1, NMoments}];

(* tau_1 = alpha_2 / (-alpha_1) = Pade[0,1] *)
Print[];
tau1Pade = N[10^6 alphas[[2]] / (-alphas[[1]]), 12];
Print["tau_1 (Pade[0,1]) = ", tau1Pade, " us"];
Print["Phase 7-11 reference: tau_1 = 16.98 us (Pade[0,1]), 19.77 us (Cauer K=8)"];
