(* sphere_analytical_cln_highstage.wls
   High-precision (60-digit) Cauer-I extraction for sphere from Stoll
   Bessel analytical Foster spectrum. Goes to ~30 stages.

   Uniform B_0 = 1 T applied to Cu sphere (a=10mm, σ=58e6).

   Stoll closed form (= analytical BEM-Foster for the sphere):
     τ_n = μ_0 σ a² / (n² π²),     n = 1, 2, ...
     a_n = 6 / (n² π²)             (Foster amplitudes)

   Joule loss (1T applied):
     R_0 = 15 / (2π σ B_0² a⁵)

   Cauer-I extraction: Hankel-Padé continued-fraction inversion of the
   Stoll Taylor-series moments
     α_k = (-1)^k Σ_n a_n × τ_n^(k+1).
*)

prec = 240;      (* working precision in digits *)
nModes = 12000;  (* Stoll modes summed for each moment *)
nMoments = 200;  (* Taylor moments computed *)
maxStages = 60;  (* Cauer (R,L) pairs requested *)

(* === Physical parameters === *)
mu0 = 4 Pi 10^-7;
sigma = SetPrecision[58 10^6, prec];
B0 = SetPrecision[1, prec];
a = SetPrecision[10/1000, prec];

Print["================================================================"];
Print[" Sphere analytical CLN (high-stage)"];
Print["   a = ", N[a 1000, 6], " mm,  σ = ", N[sigma, 6], " S/m,  B_0 = ", N[B0, 6], " T"];
Print["   working precision = ", prec, " digits, modes = ", nModes,
      ", moments = ", nMoments, ", stages = ", maxStages];
Print["================================================================"];

R0 = 15/(2 Pi sigma B0^2 a^5);
Print[];
Print["1. R_0 from Joule loss (1T 印加, ω² で割った損失逆数):"];
Print["   R_0 = 15 / (2π σ B_0² a⁵) = ", N[R0, 12]];

(* === Stoll spectrum (Foster) === *)
tauStoll[n_] := SetPrecision[mu0 sigma a^2/(n^2 Pi^2), prec];
fostAmplitude[n_] := SetPrecision[6/(n^2 Pi^2), prec];

Print[];
Print["2. Stoll Foster pole spectrum (closed form, = analytical BEM-Foster):"];
Print[Table[
  Row[{"   τ_", n, " = ", N[tauStoll[n] 10^6, 12], " μs"}],
  {n, 1, 10}] // TableForm];

(* === Taylor moments of f(s) = Σ a_n τ_n / (1 + sτ_n) === *)
(* d^k/ds^k at s=0:  α_k = (-1)^k Σ a_n τ_n^(k+1) (no factorial here) *)

alphas = Table[
  (-1)^k Sum[fostAmplitude[n] tauStoll[n]^(k + 1), {n, 1, nModes}],
  {k, 0, nMoments - 1}];

Print[];
Print["3. First α_k (high precision):"];
Do[
  Print["   α_", k, " = ", ScientificForm[N[alphas[[k+1]], 12]]],
  {k, 0, 5}];

(* === Hankel-Padé / continued-fraction inversion === *)
invertTaylor[c_List] := Module[{n, e, k, j},
  n = Length[c];
  e = Table[0, {n}];
  e[[1]] = 1/c[[1]];
  Do[
    e[[k + 1]] = -Sum[c[[j + 1]] e[[k - j + 1]], {j, 1, k}]/c[[1]],
    {k, 1, n - 1}];
  e
];

cauerExtract[cTaylor_List, maxS_] := Module[
  {c, pList, e, stage},
  c = cTaylor;
  pList = {};
  Do[
    If[Length[c] < 2 || Abs[c[[1]]] < 10^-(prec - 5),
       Print["   [stop at stage ", stage - 1, ": numerical floor]"];
       Return[pList]];
    e = invertTaylor[c];
    pList = Append[pList, e[[1]]];
    c = Drop[e, 1];
    , {stage, 1, maxS}];
  pList
];

pCauer = cauerExtract[alphas, 2 maxStages];

nPairs = Floor[Length[pCauer]/2];
Print[];
Print["4. Cauer-I ladder (R_{2k}, L_{2k+1}, τ_pair = L/R):"];
Print[];
Print["    k |       R_2k        |       L_2k+1      |  τ_pair [μs]"];
Print["   ---+-------------------+-------------------+---------------"];
Do[
  Module[{R2k, Lin, L, tauUs},
    R2k = pCauer[[2 k + 1]];
    Lin = pCauer[[2 k + 2]];
    L = If[Abs[Lin] > 10^-(prec - 10), 1/Lin, Indeterminate];
    tauUs = If[NumericQ[L] && R2k =!= 0, L/R2k 10^6, Indeterminate];
    Print[
      Row[{StringPadLeft[ToString[k], 5], " | ",
           StringPadLeft[ToString[ScientificForm[N[R2k, 10], NumberFormat -> (Row[{#1, "e", #3}] &)]], 17], " | ",
           StringPadLeft[ToString[ScientificForm[N[L, 10], NumberFormat -> (Row[{#1, "e", #3}] &)]], 17], " | ",
           StringPadLeft[ToString[N[tauUs, 10]], 14]}]];
  ], {k, 0, nPairs - 1}];

(* === Identify the largest reliable stage (no Indeterminate / Infinity) === *)
isFinite[x_] := NumericQ[x] && Abs[x] < 10^200 && Abs[x] > 10^-200;
nReliable = nPairs;
Do[
  If[Length[pCauer] >= 2 k + 2,
    If[!isFinite[pCauer[[2 k + 1]]] || !isFinite[pCauer[[2 k + 2]]],
      nReliable = k; Break[]]],
  {k, 0, nPairs - 1}];

Print[];
Print["Reliable Cauer stages with prec=", prec, " : k = 0..", nReliable - 1];

(* === Save high-precision results to JSON (reliable stages only) === *)
results = <|
  "params" -> <|"a_m" -> N[a, 20], "sigma_Spm" -> N[sigma, 20],
                "B0_T" -> N[B0, 20], "mu0" -> N[mu0, 20]|>,
  "R0" -> N[R0, 20],
  "n_reliable" -> nReliable,
  "stoll_tau_us" -> Table[N[tauStoll[n] 10^6, 16], {n, 1, 60}],
  "stoll_a_n"   -> Table[N[fostAmplitude[n], 16], {n, 1, 60}],
  "alpha_taylor" -> Table[N[alphas[[k+1]], 16], {k, 0, nMoments - 1}],
  "cauer_R" -> Table[N[pCauer[[2 k + 1]], 16], {k, 0, nReliable - 1}],
  "cauer_L_inv" -> Table[N[pCauer[[2 k + 2]], 16], {k, 0, nReliable - 1}],
  "cauer_L" -> Table[N[1/pCauer[[2 k + 2]], 16], {k, 0, nReliable - 1}],
  "cauer_tau_us" -> Table[N[(1/pCauer[[2 k + 2]])/pCauer[[2 k + 1]] 10^6, 16],
                          {k, 0, nReliable - 1}]
|>;

SetDirectory[NotebookDirectory[] /. _NotebookDirectory -> DirectoryName[$InputFileName]];
outPath = FileNameJoin[{Directory[], "sphere_analytical_cln_highstage.json"}];
Print[];
Print["Writing to: ", outPath];
Export[outPath, results, "JSON"];
Print["File exists after Export: ", FileExistsQ[outPath]];
If[FileExistsQ[outPath], Print["Size: ", FileByteCount[outPath], " bytes"]];

Print[];
Print["================================================================"];
Print[" Saved: sphere_analytical_cln_highstage.json"];
Print["================================================================"];
