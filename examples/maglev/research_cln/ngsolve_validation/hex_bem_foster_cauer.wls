(* ::Package:: *)

(* hex_bem_foster_cauer.wls
   Improvement of HCurl BEM result toward Phase 7-11 reference 19.77 us.

   Pipeline (same as bem_foster_spectral_v6_hybrid.py):
   1. K, M matrices at WP=20 (HCurl order 1, 54 DOFs)
   2. Generalized eigenvalue problem -> (lambda_n, phi_n)
   3. M-orthonormalize eigenvectors
   4. Project external A_ext = (-y/2, x/2, 0) onto each mode -> g_n
   5. Compute alpha_n moments: alpha_n = (-1)^n (mu0 sigma)^n / V Σ g_k² λ_k^(n-1)
      where lambda_k_phys = lambda_k_mine / mu0  (since K has mu0/(4pi) factor)
   6. Cauer-I continued fraction extraction -> tau_1 leading pole

   Date: 2026-05-03
*)

Get["w:\\30_CauerLadderNetwork\\2026_04_01_\:9577\:65b9\:5f62CLN\\ngsolve_validation\\schoberl_zaglmayr_basis.wls"];

ClearAll[x, y, z, xp, yp, zp];
$MaxExtraPrecision = 50;

Print[StringRepeat["=", 70]];
Print["HCurl order 1 BEM -> Foster projection + Cauer-I (target 19.77 us)"];
Print[StringRepeat["=", 70]];

(* === Constants === *)
sigma = SetPrecision[58 10^6, 30];
mu0   = 4 Pi SetPrecision[10^-7, 30];
{aBox, bBox, cBox} = SetPrecision[{5 10^-3, 2 10^-3, 1 10^-3}, 30];
V = aBox bBox cBox;

(* === HCurl order 1 basis === *)
ord = 1;
shapesRef = SchoberlZaglmayrBasis`ShapesHEX[ord, Range[0, 7], x, y, z];
nDOF = Length[shapesRef];
Print["HCurl order ", ord, ": ", nDOF, " DOFs"];

(* === Mass matrices === *)
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
Mphys = N[(bBox cBox / aBox) MrefXX +
          (aBox cBox / bBox) MrefYY +
          (aBox bBox / cBox) MrefZZ, 30];
Print["M matrix in ", Round[AbsoluteTime[] - t0, 0.1], " s"];

(* === K matrix at WP=20 === *)
Kentry[i_, j_, WP_] := Module[{vi, vj, integrand, factor, scaled},
  vi = shapesRef[[i, 1]];
  vj = shapesRef[[j, 1]] /. {x -> xp, y -> yp, z -> zp};
  integrand = (vi[[1]] vj[[1]] / aBox^2 +
               vi[[2]] vj[[2]] / bBox^2 +
               vi[[3]] vj[[3]] / cBox^2);
  factor = aBox^2 bBox^2 cBox^2;
  scaled = factor / Sqrt[aBox^2 (x - xp)^2 + bBox^2 (y - yp)^2 + cBox^2 (z - zp)^2];
  (mu0 / (4 Pi)) NIntegrate[
    integrand scaled,
    {x, 0, 1}, {y, 0, 1}, {z, 0, 1},
    {xp, 0, 1}, {yp, 0, 1}, {zp, 0, 1},
    Method -> "GlobalAdaptive", WorkingPrecision -> WP,
    PrecisionGoal -> 4, AccuracyGoal -> 6, MaxRecursion -> 8
  ]
];

Print[];
Print["Building K matrix (WP=20) ..."];
nUnique = nDOF (nDOF + 1) / 2;
K = Table[0, {nDOF}, {nDOF}];
totalTime = 0;
counter = 0;
Do[
  counter++;
  t0 = AbsoluteTime[];
  K[[i, j]] = Kentry[i, j, 15];
  K[[j, i]] = K[[i, j]];
  totalTime += AbsoluteTime[] - t0;
  If[Mod[counter, 200] == 0,
    Print["  ", counter, "/", nUnique, " entries; ",
          Round[totalTime/60, 0.1], " min so far"]
  ],
  {i, 1, nDOF}, {j, i, nDOF}
];
Print["K matrix done in ", Round[totalTime/60, 0.1], " min"];
Kmat = N[K, 30];

(* === Solve generalized eigen === *)
Print[];
Print["Solving Eigensystem[{K, M}] ..."];
{eigvals, eigvecs} = Eigensystem[{Kmat, Mphys}];
Print["Eigensystem solved, ", Length[eigvals], " modes"];

(* M-orthonormalize *)
eigvecsNorm = Table[
  Module[{v = eigvecs[[k]]}, v / Sqrt[Re[v . Mphys . v]]],
  {k, 1, Length[eigvals]}
];

(* === Project A_ext = (-y/2, x/2, 0) (uniform B_z = 1 field) === *)
(* In physical: phi_phys_alpha = phi_ref_alpha / d_alpha
   <phi_i, A_ext>_phys = int phi_i_x_phys (-y/2) + phi_i_y_phys (x/2) dV_phys
                      = int (phi_ref_x/a)(-b y_ref/2) + (phi_ref_y/b)(a x_ref/2) (abc) dV_ref
                      = -bc/2 int phi_ref_x y_ref dV_ref + ac/2 int phi_ref_y x_ref dV_ref *)

bVec = Table[
  -bBox cBox/2 * fastIntegrate3D[Expand[shapesRef[[i, 1, 1]] y], {x, y, z}] +
   aBox cBox/2 * fastIntegrate3D[Expand[shapesRef[[i, 1, 2]] x], {x, y, z}],
  {i, 1, nDOF}
];

(* g_k = ⟨A_ext, phi_k⟩_M-normalized = v_k_normalized . bVec *)
gVec = Table[Re[eigvecsNorm[[k]] . bVec], {k, 1, Length[eigvals]}];

(* === alpha_n moments (Phase 7-11 formula) ===
   lambda_phys = lambda_mine / mu0  (since K has mu0/(4pi) inside)
   alpha_n = (-1)^n (mu0 sigma)^n / V * Σ_k g_k^2 lambda_phys^(n-1) *)

lambdaPhys = Re[eigvals] / mu0;
mu0sigma = mu0 sigma;

NMoments = 8;
alphas = Table[
  (-1)^nM Sum[gVec[[k]]^2 lambdaPhys[[k]]^(nM - 1), {k, 1, Length[eigvals]}] mu0sigma^nM / V,
  {nM, 1, NMoments}
];

Print[];
Print["alpha_n moments:"];
Do[Print["  alpha_", k, " = ", N[alphas[[k]], 12]], {k, 1, NMoments}];

Print[];
tau1Pade = N[10^6 alphas[[2]] / (-alphas[[1]]), 12];
Print["tau_1 (Pade[0,1] = alpha_2/-alpha_1) = ", tau1Pade, " us"];
Print["Phase 7-11 reference Pade[0,1]: ~ 16.98 us"];

(* === Cauer-I continued fraction extraction === *)
Print[];
Print["Cauer-I continued fraction extraction ..."];

(* Taylor coefficients of chi(s): c_n = alpha_n  (where chi(s) = sum alpha_n s^n)
   But alpha_n in Phase 7-11 includes (-1)^n already, so c_n = alpha_n *)

invertTaylor[c_] := Module[{n = Length[c], e, k, s},
  e = Table[0, {n}];
  e[[1]] = 1/c[[1]];
  Do[
    s = Sum[c[[j+1]] e[[k - j + 1]], {j, 1, k}];
    e[[k + 1]] = -s / c[[1]],
    {k, 1, n - 1}
  ];
  e
];

cauerCF[cTaylor_, maxStages_] := Module[{c = cTaylor, pList = {}, e, k},
  Do[
    If[Length[c] < 2 || Abs[c[[1]]] < 10^-200, Break[]];
    e = invertTaylor[c];
    AppendTo[pList, e[[1]]];
    c = Drop[e, 1],
    {k, 1, maxStages}
  ];
  pList
];

cauerTruncatePoles[pList_, K_] := Module[{N0, D0, k},
  If[K < 1, Return[{}]];
  N0 = {1};
  D0 = {pList[[K]]};
  Do[
    Module[{sN = Prepend[N0, 0],
            DnewPk = pList[[k]] * D0,
            maxLen, Dnew},
      maxLen = Max[Length[DnewPk], Length[sN]];
      Dnew = Table[0, {maxLen}];
      Do[Dnew[[i]] += DnewPk[[i]], {i, 1, Length[DnewPk]}];
      Do[Dnew[[i]] += sN[[i]], {i, 1, Length[sN]}];
      N0 = D0; D0 = Dnew
    ],
    {k, K - 1, 1, -1}
  ];
  {N0, D0}
];

(* Use alpha_n as Taylor coefficients of chi(s) *)
(* Need alpha_0 too. From Phase 7-11: M_0 = chi(0) = sum a_k tau_k (some moment).
   For us: alpha_0 = c_0 = sum_k g_k^2 lambda_phys_k^(-1) ... wait *)
(* Actually for Phase 7-11: alpha_n for n >= 1, with alpha_0 = M_0 separately.
   M_0 = chi(s=0) = sum a_k tau_k.
   For our pipeline: alpha_n includes lambda^(n-1), so for n=0 we'd have lambda^(-1).
   But we don't compute alpha_0 directly. From cauer_depth_test.py:
     M0 = mp.mpf("3.058972587463514598008059342935820") (hardcoded from Phase 7 P=8 WP=30)
   For HCurl, we need to compute the corresponding M0.
   Formula: M_0 = (mu0 sigma)^0 / V * sum_k g_k^2 lambda_phys^(-1) doesn't work since alpha_n
            assumes (-1)^n factor and integer powers >= 1.
   Actually: chi(s) = sum a_k tau_k / (1 + s tau_k) = sum a_k tau_k - sum a_k tau_k^2 s + ...
   So coefficient of s^0: M_0 = sum a_k tau_k (with proper a_k).
   For Phase 7-11 normalization: chi has units of seconds (?), unclear.
*)

(* Try Cauer with alpha_n directly as Taylor coefficients of chi *)
(* This may need an alpha_0 prepended *)
(* For now, use the alpha_n's and see if Cauer extraction gives sensible poles *)

(* Construct Taylor: c_n = alpha_(n+1) starting from n=0. So chi(s) ~ alpha_1 + alpha_2 s + alpha_3 s^2 + ... *)
(* Wait that's not right either. Let me re-derive. *)
(* chi(s) := sum a_k tau_k / (1 + s tau_k) *)
(* Expand: = sum a_k tau_k (1 - s tau_k + s^2 tau_k^2 - ...) *)
(* = sum_{n>=0} (-s)^n sum_k a_k tau_k^(n+1) *)
(* So c_n = (-1)^n M_n where M_n = sum a_k tau_k^(n+1) *)
(* alpha_n in Phase 7-11 = (-1)^n (mu0 sigma)^n / V sum_k g_k^2 lambda_phys^(n-1) *)
(* For tau_k = mu0 sigma lambda_k, alpha_n = (-1)^n / V sum g_k^2 (mu0 sigma)^n lambda^(n-1) = (-1)^n / V (mu0 sigma) sum g_k^2 tau_k^(n-1) (mu0 sigma)^(n-1) / (mu0 sigma)^(n-1) ... messy *)

(* Pragmatic: just use alpha_n as Taylor coefficients of an effective chi-function and run Cauer.
   The Cauer poles will be the τ_n in some normalization. *)

cTaylor = Prepend[Table[alphas[[k]], {k, 1, NMoments}], 0];
(* Alpha-prepend with placeholder; recompute properly: *)

(* Following cauer_depth_test.py exactly:
   c_taylor = [M_0] + [(-1)^n M_n for n >= 1]
   where M_0 from Phase 7 = 3.058...
   But M_n = (-1)^n alpha_n  (since alpha_n = (-1)^n M_n in their convention)
   So c_n = M_n for n >= 1. But c_0 = M_0 separately.
   Hmm: c_taylor in cauer_depth_test = [M_0, M_1, M_2, ..., M_8]
        and M_n = (-1)^n alpha_n for n >= 1.

   For us: M_0 should be  (mu0 sigma) / V * sum g^2 * lambda_phys^(-1) ?
   Or M_0 = sum a_k tau_k (Foster expansion at s=0).
*)

(* Recompute M_0 *)
M0 = Sum[gVec[[k]]^2 / lambdaPhys[[k]], {k, 1, Length[eigvals]}] mu0 sigma / V;
M0 = Re[M0];
Print["M_0 = ", N[M0, 12]];

cTaylorFull = Prepend[Table[(-1)^n alphas[[n]], {n, 1, NMoments}], M0];
Print["Taylor coefficients (M_0..M_8):"];
Do[Print["  c_", k - 1, " = ", N[cTaylorFull[[k]], 12]], {k, 1, Length[cTaylorFull]}];

(* Cauer extraction *)
pList = cauerCF[cTaylorFull, NMoments];
Print[];
Print["Cauer p_k stages: ", Length[pList]];
Do[Print["  p_", k, " = ", N[pList[[k]], 12]], {k, 1, Length[pList]}];

(* Cauer truncation -> poles *)
findRealNegRoots[coeffs_] := Module[{cHi, roots},
  cHi = Reverse[coeffs];
  roots = Quiet[NSolve[Sum[cHi[[k+1]] t^(Length[cHi]-1-k), {k, 0, Length[cHi]-1}] == 0, t,
                       WorkingPrecision -> 25]];
  Select[t /. roots, (Im[N[#]] == 0) && (Re[N[#]] < 0) &]
];

Print[];
Print["Cauer truncation -> recovered tau_k [us]:"];
Do[
  Module[{ND, Nc, Dc, poles, taus},
    ND = cauerTruncatePoles[pList, K];
    {Nc, Dc} = ND;
    poles = findRealNegRoots[Dc];
    If[Length[poles] > 0,
      taus = Sort[-1/Re[poles], Greater];
      Print["  K=", K, ": ", Length[poles], " poles, taus = ",
            N[10^6 # & /@ taus, 10]],
      Print["  K=", K, ": no negative real poles"]
    ]
  ],
  {K, 2, Length[pList]}
];

Print[];
Print["============================================"];
Print["Phase 7-11 reference (P=4, K=8): tau_1 = 19.71 us"];
Print["Phase 7-11 reference (P=6): tau_1 = 19.77 us"];
Print["============================================"];
