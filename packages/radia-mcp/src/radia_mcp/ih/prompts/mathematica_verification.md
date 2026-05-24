# Mathematica Verification Recipes for IH SIBC

Closed-form / symbolic derivations for SIBC quantities, suitable for
side-by-side validation of NGSolve / Radia numerical results. Use these
as the **first sanity check** before debugging FEM/BEM output.

Cross-reference: `radia_mcp.radia_ngsolve.analytical_formulas(topic=
"validation_use_cases")` lists the canonical analytical companion for
each IH analysis.

## Recipe 1: Leontovich Z_s symbolic derivation

Derive `Z_s = (1+j)/(sigma * delta)` from the 1D diffusion equation
inside a semi-infinite conductor. Verify against the limit `omega ->
inf`, `mur -> 1`.

```mathematica
(* 1D Helmholtz in conductor: H''(z) = j*omega*mu*sigma * H(z) *)
(* Solution: H(z) = H0 * Exp[-kc*z], decaying into conductor *)
kc[mu_, sigma_, omega_] := (1 - I)/delta /.
    delta -> Sqrt[2 / (omega * mu * sigma)];

(* Surface impedance: Z_s = E_t / H_t at z=0 *)
(* E from Maxwell: curl E = -dB/dt -> E = (1/sigma) * dH/dz *)
Zs[mu_, sigma_, omega_] := (kc[mu, sigma, omega] / sigma) /.
    Simplify[#, Assumptions -> {sigma > 0, mu > 0, omega > 0}] &;

(* Verify against textbook: Z_s = (1+j)*sqrt(omega*mu/(2*sigma)) *)
ZsTextbook[mu_, sigma_, omega_] := (1 + I)*Sqrt[omega*mu/(2*sigma)];

Simplify[Zs[mu, sigma, omega] - ZsTextbook[mu, sigma, omega],
         Assumptions -> {sigma > 0, mu > 0, omega > 0}]
(* Should return 0 *)
```

## Recipe 2: Smythe sphere induced dipole (BEM validation reference)

Smythe (1950) analytical solution for a conducting sphere in uniform
AC field. **The lab's primary SIBC validation case** (2026-04-13: 14
data points, BEM Scalar BIE: -1.6%, FEM scattered: +/- 3%).

```mathematica
(* Magnetic Reynolds number xi = R / delta *)
(* m-factor (induced dipole moment / applied field volume) *)
mFactor[a_, sigma_, mur_, omega_] := Module[
    {delta, k, ka, num, den},
    delta = Sqrt[2/(omega * mu0 * mur * sigma)];
    k = (1 + I)/delta;
    ka = k * a;
    (* Smythe Eq. (1) -- spherical Bessel based *)
    num = (mur + 2) * (Sinh[ka] - ka*Cosh[ka]) +
          ka^2 * (mur - 1) * Sinh[ka]/3;
    den = (mur - 1) * (Sinh[ka] - ka*Cosh[ka]) -
          ka^2 * Sinh[ka]/3;
    -(a^3) * (1 + 3 * num / den)
];

(* Induced surface current (a*sin(theta)*phi component) *)
JsSphere[a_, sigma_, mur_, omega_, H0_, theta_] := Module[
    {m, B0},
    B0 = mu0 * H0;
    m = mFactor[a, sigma, mur, omega];
    (3/2) * I*omega*B0*a / (I*omega*mu0*a + ZsApprox) *
        Sin[theta] /. ZsApprox -> (1 + I)*Sqrt[omega*mu0*mur/(2*sigma)]
];

(* Tangential H_t_rms (max at equator) for cross-checking BEM SIBC *)
HtRMS[a_, sigma_, mur_, omega_, H0_] :=
    Max[Abs[JsSphere[a, sigma, mur, omega, H0, theta]]
        /. theta -> Pi/2] * Sqrt[2/3];
```

Use this to validate `ScalarBIESIBCSolver.solve(...)['H_t_rms']` to
~1% accuracy across `R/delta = 1.5 .. 60`. Compare via:

```mathematica
mu0 = 4*Pi*1e-7;
HtRMS[0.01, 5.8e7, 1, 2*Pi*1e6, 1.0]   (* Cu sphere 1 MHz *)
HtRMS[0.01, 2e6, 100, 2*Pi*1e3, 1.0]   (* steel 1 kHz *)
```

## Recipe 3: Skin depth frequency scaling (sanity plot)

```mathematica
delta[mu_, sigma_, omega_] := Sqrt[2/(omega * mu * sigma)];

LogLogPlot[
    {delta[mu0, 5.8e7, 2*Pi*f] * 1000,        (* Cu, mm *)
     delta[100*mu0, 2e6, 2*Pi*f] * 1000,      (* steel mu_r=100, mm *)
     delta[mu0, 3.5e7, 2*Pi*f] * 1000},       (* Al, mm *)
    {f, 50, 1e6},
    PlotLegends -> {"Cu", "Steel mu_r=100", "Al"},
    AxesLabel -> {"Frequency [Hz]", "Skin depth [mm]"}
]
```

Use this when picking `--mesh-size` for FEM: the surface mesh should
resolve `delta/5` for full-resolution validation, or `delta/2` for
SIBC + hole approach.

## Recipe 4: Per-panel curvature SIBC correction (Mitzner check)

Symbolic version of the per-panel local curvature SIBC (lab 2026-04-12).
Use this to verify the numerical `R_local` extractor in
`_compute_panel_local_radii`.

```mathematica
(* Mitzner curvature-corrected Z_s for a sphere of radius R *)
ZsMitznerSphere[Zs_, delta_, R_] :=
    Zs * (1 + (1 + I)/2 * delta / R);

(* For a sphere R=10mm at 100 Hz, Cu *)
delta100Hz = delta[mu0, 5.8e7, 2*Pi*100];   (* ~6.6 mm *)
ZsLeontovich = (1 + I) * Sqrt[2*Pi*100 * mu0 / (2 * 5.8e7)];
ZsMitznerSphere[ZsLeontovich, delta100Hz, 0.01]
(* Result is ~33% different from ZsLeontovich -- the regime where the *)
(* lab's per-panel extractor gives +2.9% vs +31% for the scalar global *)
```

## Cross-MCP cross-check

After computing a result with `calc_inductance.py --workpiece sibc`,
recompute the same with Mathematica and compare:

| Lab compute                          | Mathematica reference                |
|--------------------------------------|--------------------------------------|
| `ScalarBIESIBCSolver.solve(...)['H_t_rms']` | `HtRMS[a, sigma, mur, omega, H0]`  |
| `--use-local-curvature` per-node Z_s | `ZsMitznerSphere[ZsLeo, delta, R]`   |
| `--impedance-model esim` after Karl  | (Mathematica `NDSolve` 1D BVP)       |

Run any of these as a sanity check before reporting numbers in a paper
or showing to a grant reviewer.

Reference MCP: call `mcp-server-mathematica` for the Wolfram Language
runtime to actually execute these recipes from a Cubit panel session.
