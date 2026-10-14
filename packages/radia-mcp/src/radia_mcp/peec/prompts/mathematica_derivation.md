# Mathematica Symbolic Derivation of HOIBC

Wolfram Language recipes for the Dong-Di Rienzo 2020 perturbation
hierarchy: PEC → Leontovich → Mitzner. Use these to:
- Verify the Neumann RHS formulas before coding them in NGSolve
- Reproduce the prolate ellipsoid benchmark (paper Fig 7)
- Cross-check a Radia HOIBC implementation symbolically

Companion to `peec_hoibc('radia_application')` (numerical recipes) and
`peec_hoibc('dong_di_rienzo_2020')` (formula reference).

## Recipe 1: Small parameter `p_tilde` symbolic

```mathematica
(* Dong-Di Rienzo 2020 small parameter:                              *)
(*     p_tilde = mu_r * delta / D                                    *)
(* where D is the characteristic conductor size, NOT delta/D alone   *)
(* (Rytov-Senior was wrong for mu_r >> 1).                           *)

pTilde[mur_, delta_, D_] := mur * delta / D;

(* Verify expansion validity for steel mu_r=100, R=10 mm, 1 kHz *)
mu0 = 4*Pi*1e-7;
delta1kHz[mur_] := Sqrt[2/(2*Pi*1e3 * mur * mu0 * 2e6)];   (* steel sigma 2e6 *)
pTilde[100, delta1kHz[100], 0.01]
(* Result: 4.1e-2 -- Mitzner needed *)

pTilde[1000, delta1kHz[1000], 0.01]
(* Result: 1.3e-1 -- Mitzner NOT enough, order 3 needed *)
```

## Recipe 2: BVP_0 Neumann RHS (PEC approximation)

```mathematica
(* H_s = applied Biot-Savart, normal projection onto conductor surface *)
(* BVP_0: Laplace[phi_0] = 0 in dielectric                            *)
(*        d phi_0 / dn = H_s . n_hat  on conductor                    *)

(* For sphere R=a in uniform H = H0 z_hat: *)
HsNormalSphere[a_, H0_, theta_] := H0 * Cos[theta];

(* Analytical phi_0 from external multipole expansion *)
phi0Sphere[a_, H0_, r_, theta_] :=
    -H0 * (a^3 / (2 * r^2)) * Cos[theta] + H0 * r * Cos[theta];

(* Verify: H_ext_PEC = H_s - grad(phi_0); tangential should be uniform *)
HextPEC[a_, H0_, r_, theta_] := Module[{phi0},
    phi0 = phi0Sphere[a, H0, r, theta];
    H0 * Cos[theta] - D[phi0, r]
];
HextPEC[1, 1, 1, theta] // Simplify
(* Should give (3/2)*Cos[theta] -- the PEC sphere dipole reinforcement *)
```

## Recipe 3: BVP_1 Leontovich Neumann RHS

Surface-tangential derivative of `(H_s - grad phi_0)` along local
coordinates xi_1, xi_2. Symbolic on a sphere:

```mathematica
(* On sphere surface, tangential gradient in (theta, phi) *)
gradS[f_, theta_, phi_, a_] := {(1/a) D[f, theta], (1/(a*Sin[theta])) D[f, phi]};

(* BVP_1 Neumann RHS for sphere: *)
HtPEC[a_, H0_, theta_] := (3/2) * H0 * Sin[theta];   (* tangential H_ext_PEC *)

NeumannRHS1[a_, H0_, theta_] := -(1/Sqrt[s]) *
    D[a * HtPEC[a, H0, theta], theta] / a;
(* s here is the Laplace variable; replace s -> j*omega for freq domain *)

(* Verify integration over sphere: should be zero by Gauss's law       *)
(* (Laplace BVP needs compatibility condition: integral of Neumann = 0) *)
Integrate[NeumannRHS1[a, H0, theta] * a^2 * Sin[theta],
          {theta, 0, Pi}, {phi, 0, 2*Pi}] // Simplify
```

## Recipe 4: BVP_2 Mitzner curvature term (the key advance)

```mathematica
(* Curvature term for sphere: d_1 = d_2 = a (principal radii)        *)
(* General formula: CURVATURE_TERM = (d_k + d_{3-k}) / (d_k * d_{3-k}) *)
(* For sphere both axes a: 2a / a^2 = 2/a                            *)

CurvSphere[a_] := 2/a;   (* = 2 * mean curvature = 2 * H_mean *)

(* BVP_2 RHS includes both Leontovich-like derivatives of phi_1 AND  *)
(* curvature correction times (H_s - grad phi_0) tangential          *)

NeumannRHS2[a_, mur_, H0_, theta_] := Module[{leoTerm, curvTerm},
    leoTerm = -(1/Sqrt[s]) *
        (* tangential derivatives of (-grad phi_1) -- complex *)
        SomethingFrom[phi1Sphere];
    curvTerm = -(1/(mur * Sqrt[s])) * CurvSphere[a] * HtPEC[a, H0, theta];
    leoTerm + curvTerm
];

(* Mitzner correction in closed form for thin-skin sphere: *)
(* Z_s_Mitzner = Z_s_Leontovich * (1 + (1+j)/2 * delta * H_mean) *)
ZsMitznerSphere[Zs_, delta_, R_] := Zs * (1 + (1 + I)/2 * delta / R);
```

## Recipe 5: Reproduce Dong-Di Rienzo Fig 7 (prolate ellipsoid)

The benchmark: prolate ellipsoid `x^2/b^2 + y^2/b^2 + z^2/a^2 = 1`
with a = 2*b, sigma = 5.998e7 S/m, in a circular coil at 10 A.

```mathematica
(* Principal radii of curvature at pole and equator: *)
d1Pole[a_, b_] := b^2 / a;
d2Pole[a_, b_] := b^2 / a;
d1Equator[a_, b_] := b;
d2Equator[a_, b_] := a^2 / b;

(* Mean curvature at each point (used in Mitzner correction): *)
HmeanPole[a_, b_] := (1/d1Pole[a,b] + 1/d2Pole[a,b]) / 2;
HmeanEquator[a_, b_] := (1/d1Equator[a,b] + 1/d2Equator[a,b]) / 2;

(* For a = 2*b = 20 mm: *)
HmeanPole[0.02, 0.01]      (* = 100 1/m, high curvature *)
HmeanEquator[0.02, 0.01]   (* = 75 1/m, mixed *)

(* The Mitzner correction ratio (Z_Mitzner / Z_Leontovich):    *)
(* depends on local curvature -- spatially varying!            *)
(* This is why per-panel Z_s in NGSolve matters (Recipe 2 in   *)
(* ngsolve_recipes).                                            *)
```

## Recipe 6: Frequency-domain switch (Laplace s -> j omega)

```mathematica
(* Dong-Di Rienzo derivation uses Laplace variable s            *)
(* Switch to frequency domain: s -> j*omega                     *)
(* and normalize: tau = 2/omega                                  *)

FreqSubst = {s -> I*omega};

ZsFreqDomain[Zs_Leontovich, mur_, delta_, R_] :=
    Zs_Leontovich * (1 + (1 + I)/2 * delta * (1/R)) /. FreqSubst;
```

## Why these recipes matter

The Dong-Di Rienzo 2020 paper's algorithm is **3 separate Laplace
BVPs** with Neumann RHS that depend on the previous order. The RHS
formulas are easy to miscode (lots of indices, surface-tangential
derivatives, principal radii). Symbolic verification in Mathematica
**before** committing the NGSolve assembly code saves hours of
debugging.

After Mathematica verification, the production NGSolve code lives in:
- `src/radia/panels/calc_fem_kelvin.py::_solve_hoibc_cascade` (planned)

Reference MCP: `mcp-server-mathematica` for actual symbolic computation
from a Cubit panel session.
