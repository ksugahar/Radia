(* sphere_analytical_cln.wls
   Mathematica analytical computation of Cauer ladder for sphere with
   uniform B applied.

   Setup: Cu sphere radius a, σ uniform, uniform B_0 = 1 T applied.

   Analytical results:
     1) R_0 = 15 / (2π σ B_0² a⁵)  [from ∫σ|A_imposed|² dV with A = (B_0/2) r sin(θ)]
     2) Stoll Foster spectrum: τ_n = μ_0 σ a² / (n² π²), n = 1, 2, ...
     3) Foster amplitudes: from χ(s) = -1 + Σ (6/(n²π²))/(1 + sτ_n)
     4) Hankel-Padé Cauer ladder from Stoll Foster moments
*)

(* === Physical parameters === *)
mu0 = 4 Pi 10^-7;
sigma = 58 10^6;        (* Cu *)
B0 = 1;
a = 10/1000;            (* 10 mm *)

Print["=" <> StringRepeat["=", 70]];
Print[" Sphere analytical CLN: a = ", a 1000, " mm, σ = ", sigma, " S/m"];
Print["=" <> StringRepeat["=", 70]];
Print[];

(* === R_0 from Joule loss formula === *)
(* ∫σ|A_imposed|² dV = ∫_{sphere} σ × (B_0 r sin(θ)/2)² dV
                     = σ × B_0²/4 × ∫_0^a r⁴ dr × ∫_0^π sin³θ dθ × ∫_0^2π dφ
                     = σ × B_0²/4 × a⁵/5 × 4/3 × 2π
                     = (2π σ B_0² a⁵) / 15
*)
intsigmaAsq = (2 Pi sigma B0^2 a^5)/15;
R0canonical = 1/intsigmaAsq;

Print["1. R_0 from Joule loss formula (1T 印加, ω² で割った損失逆数):"];
Print["   R_0 = 15 / (2π σ B_0² a⁵)"];
Print["   Analytical R_0 = ", N[R0canonical, 10], " (dimensionless units)"];
Print[];

(* === Stoll Foster spectrum === *)
tauStoll[n_] := mu0 sigma a^2 / (n^2 Pi^2);

Print["2. Stoll Foster pole spectrum:"];
Print["   τ_n = μ_0 σ a² / (n² π²)"];
Print[Table[
  Row[{"   τ_", n, " = ", N[tauStoll[n] 10^6, 6], " μs"}],
  {n, 1, 6}] // TableForm];
Print[];

(* === Foster amplitudes === *)
(* χ(s) = -1 + Σ a_n / (1 + s τ_n) where a_n = 6/(n² π²) *)
fostAmplitude[n_] := 6/(n^2 Pi^2);

Print["3. Stoll Foster amplitudes (χ(s) representation):"];
Print["   a_n = 6 / (n² π²)"];
Print[Table[
  Row[{"   a_", n, " = ", N[fostAmplitude[n], 8]}],
  {n, 1, 6}] // TableForm];
Print["   Σ a_n = ", N[Sum[fostAmplitude[n], {n, 1, Infinity}], 10],
       " (= 1, satisfies χ(0) = 0)"];
Print[];

(* === Compute α_n moments for Cauer extraction === *)
(* α_n = Σ_k a_k × τ_k^(n+1) × (-1)^n *)
(* Note: in our convention, χ(s) = -1 + Σ a_k/(1+sτ_k), and we want Cauer-I
         expansion of (1 - χ(s)) or similar to get the Cauer ladder.
         The α_n moments are derivatives at s=0. *)

(* For numerical Hankel-Padé, sum truncated Foster series: *)
nModes = 100;     (* modes to include in Foster sum *)
nMoments = 20;    (* moments to compute *)

alphas = Table[
  (-1)^k Sum[fostAmplitude[n] tauStoll[n]^(k + 1), {n, 1, nModes}],
  {k, 0, nMoments - 1}];

Print["4. α_n moments (from truncated Stoll Foster, n_modes=", nModes, "):"];
Print[Table[
  Row[{"   α_", k, " = ", ScientificForm[N[alphas[[k+1]], 10]]}],
  {k, 0, 5}] // TableForm];
Print[];

(* === Cauer-I extraction via continued fraction inversion === *)
(* From Taylor coefficients α_n of f(s) = α_0 + α_1 s + α_2 s² + ...:
   Continued fraction expansion at s=0 gives Cauer-I rungs (R_0, L_1, R_2, ...).
   Algorithm:
     p_0 = 1/α_0  (= "first rung")
     residue = 1/(f(s)) - p_0 (with leading term subtracted)
     ... iterate
*)

invertTaylor[c_List] := Module[{n, e, k, j},
  n = Length[c];
  e = Table[0, {n}];
  e[[1]] = 1/c[[1]];
  Do[
    e[[k + 1]] = -Sum[c[[j + 1]] e[[k - j + 1]], {j, 1, k}]/c[[1]],
    {k, 1, n - 1}];
  e
];

cauerExtract[cTaylor_List, maxStages_] := Module[
  {c, pList, e, stage},
  c = cTaylor;
  pList = {};
  Do[
    If[Length[c] < 2 || Abs[c[[1]]] < 10^-300, Return[pList]];
    e = invertTaylor[c];
    pList = Append[pList, e[[1]]];
    c = Drop[e, 1];
    , {stage, 1, maxStages}];
  pList
];

(* Apply to alphas *)
pCauer = cauerExtract[alphas, nMoments - 1];
Print["5. Cauer-I ladder extraction from Stoll α_n moments:"];
Print["   Number of p-coefficients: ", Length[pCauer]];
Print["   Number of (R, L) pairs:   ", Floor[Length[pCauer]/2]];
Print[];

(* Convert p-coefficients to R, L pairs *)
(* p[2k] = R_{2k}, p[2k+1] = 1/L_{2k+1} *)
nPairs = Floor[Length[pCauer]/2];
Print["6. Cauer-I ladder (R_{2k}, L_{2k+1}, τ_pair = L/R):"];
Print[
  Row[{StringPadLeft["k", 4], " | ", StringPadLeft["R_2k", 18], " | ",
       StringPadLeft["L_2k+1", 18], " | ", StringPadLeft["τ_pair [μs]", 18]}]];
Do[
  Module[{R2k, Lin, L, tauUs},
    R2k = pCauer[[2 k + 1]];
    Lin = pCauer[[2 k + 2]];
    L = If[Abs[Lin] > 10^-300, 1/Lin, Indeterminate];
    tauUs = If[NumericQ[L] && R2k != 0, L / R2k 10^6, Indeterminate];
    Print[
      Row[{StringPadLeft[ToString[k], 4], " | ",
           StringPadLeft[ScientificForm[N[R2k, 8]], 18], " | ",
           StringPadLeft[ScientificForm[N[L, 8]], 18], " | ",
           StringPadLeft[ScientificForm[N[tauUs, 8]], 18]}]];
  ], {k, 0, nPairs - 1}];

Print[];
Print["==== Comparison: τ_pair vs Stoll Foster τ_n ===="];
Print["Note: stage k τ_pair is Cauer rung k, NOT Foster pole k+1."];
Print["For low k, Cauer τ ≈ leading Foster τ + correction."];

(* === Comparison reference === *)
Print[];
Print["=" <> StringRepeat["=", 70]];
Print["Reference: Hiruma+Kelvin (FEM axisym ORDER=3, mesh.Curve(3)):"];
Print["   stage 0: 728.85 μs (vs Stoll τ_1 = 738.48, gap -1.3%)"];
Print["   stage 1: 171.51 μs (Cauer rung k=1 — different from Stoll τ_2)"];
Print["Reference: Kameari+Kelvin (FEM axisym):"];
Print["   stage 0: 694.14 μs (vs Stoll τ_1 = 738.48, gap -6.0%)"];
Print["   stage 1: 154.60 μs (Cauer rung k=1)"];
