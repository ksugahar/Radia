"""Litz wire + transmission line knowledge.

Distilled from W:/.../39_部品材料/02_リッツ線/ (44 files) and
W:/.../39_部品材料/03_伝送線路/ (22 files).

This knowledge layer covers the analytical, semi-analytical, and
numerical machinery used in the lab to predict AC resistance,
inductance, internal impedance, and radiation of multi-strand
conductors (Litz wire) and multiconductor transmission lines (MTL).

The primary cross-references are:
  - radia_mcp.peec.carstensen_ac_copper_loss
  - radia_mcp.motor.hollaus_eddy
  - radia_mcp.ih (workpiece SIBC + ESIM)
  - radia_mcp.pcb (compensation network selection)

All formulas are written in cp932-safe ASCII (no Unicode math).
"""


# Authoritative topic enum for the dispatcher tool (wired into
# `<short>_topics()` via common.register_topics_tool).
TOPICS: dict[str, str] = {
    "overview":               "Skin depth + AC loss landscape, when to pick Litz vs hollow vs foil.",
    "litz_construction":      "Strand / 1st-bundle / 2nd-bundle / cable hierarchy, packing factor, twist pitch.",
    "skin_proximity_basics":  "Single-wire skin Bessel, proximity orthogonality, internal vs external proximity.",
    "dowell_formula":         "1D Dowell (1966) foil/multi-layer formula; rectangular conductor AC R / DC R.",
    "ferreira_strands":       "2D Bessel-function model for round strands (Ferreira 1992 lineage, M2/M3/M4).",
    "umetani_multi_level":    "Multi-level twisted Litz, measurable-parameter analytical model (Umetani 2021).",
    "homogenization":         "FEM homogenization (Igarashi complex-permeability, Ollendorff, foil approximation).",
    "multiconductor_tl":      "MTL theory: telegrapher equations, P/L coefficients, common/normal/antenna modes.",
    "fem_thin_wire":          "FEM and MoM for thin-wire / rectangular TL impedance; nonuniform grid for skin.",
    "abc_pml":                "Asymptotic / PML-like absorbing boundary conditions for open-region MTL FEM.",
    "tlm_method":             "Transmission-line modelling method (Johns) for full Maxwell on circuit network.",
    "rosskopf_coupled":       "Coupled FEM+PEEC numeric approach for inductive components with Litz (Rosskopf 2018).",
    "ac_loss_measurement":    "Impedance analyzer / network analyzer measurement protocol for Litz AC R.",
    "all":                    "Concatenation of all sections above.",
}


OVERVIEW = r"""
# Litz wire and transmission line analysis -- overview

This subpackage groups two physically related topics:

1. **Litz wire** -- hierarchically twisted multi-strand conductor used to
   suppress AC resistance increase in the 1 kHz - 1 MHz band
   (induction heating, WPT, motor windings, high-frequency
   transformers).
2. **Multiconductor transmission lines (MTL)** -- distributed-parameter
   model R, L, G, C per unit length for parallel conductors at MHz - GHz
   (PCB traces, cable bundles, microstrips, busbars).

The two topics share the same underlying physics: skin effect,
proximity effect, internal impedance from cross-section eddy currents,
and external inductance / capacitance from the surrounding field.
Many of the analytical building blocks (Bessel-function single-wire
internal impedance, Belevitch proximity matrix, Dowell layered-foil
formula) appear in both.

## Skin depth (the master length scale)

The plane-wave skin depth in a homogeneous conductor

```
delta = sqrt(2 / (omega * mu * sigma))
      = 1 / sqrt(pi * f * mu * sigma)
```

is the depth at which the current density has decayed to 1/e of the
surface value, for a half-space conductor.

| Frequency band | Application       | Cu skin depth (mm) |
|----------------|-------------------|--------------------|
| 50 Hz          | Power grid        | 9.3                |
| 1 kHz          | Audio, induction  | 2.1                |
| 10 kHz         | IH (steel coils)  | 0.66               |
| 85 kHz         | EV WPT (SAE J2954)| 0.23               |
| 1 MHz          | WPT, high-freq IH | 0.066              |
| 6.78 MHz       | A4WP / Rezence    | 0.025              |
| 13.56 MHz      | RFID, WPT         | 0.018              |
| 100 MHz        | RF                | 0.0066             |

For *magnetic* conductors (mu_r > 1) skin depth contracts by
`sqrt(mu_r)` -- e.g. for sigma_steel ~ 4e6 S/m and mu_r ~ 100, the
delta at 50 Hz collapses from ~12 mm (in non-magnetic Cu-equivalent)
to ~0.6 mm.  This is the central reason why iron / steel workpieces
have such high AC loss in the kHz IH band.

## When to pick which mitigation

| Approach              | Use band      | Pro                                 | Con                                    |
|-----------------------|---------------|--------------------------------------|-----------------------------------------|
| **Litz wire**         | 1 kHz - 1 MHz | many fine strands, Dowell formula    | terminations, cost, hierarchical loss   |
| **Hollow conductor**  | > 100 kHz     | water-cooled, large current capacity | only outer wall conducts                |
| **Foil winding**      | broadband     | 1D AC distribution, planar           | end effects, hard to interleave         |
| **Magnetic-plated**   | NMI sensors   | mu-plating reduces flux concentration| extra process, frequency-limited       |
| **Solid + design margin** | DC, very low f | simplest                           | unacceptable above ~1 kHz for high current |

For the lab's WPT and induction heating work, Litz is the default at
20 kHz - 200 kHz, with hollow-conductor used above ~100 kHz when the
current per turn exceeds ~50 A_rms.

## Analytical / numerical method hierarchy (Meng & Biela 2020)

Loss prediction methods are typically grouped as:

1. **Analytical** (microseconds per evaluation):
   - 1D Dowell (1966) -- foil / multi-layer (M1 model)
   - 2D Ferreira (1992) -- round wires (M2 model, M3 Bartoli,
     M4 Tourkhani extensions)
   - Multi-level Sullivan & Zhang / Umetani 2021
2. **Semi-analytical** (seconds per evaluation):
   - FEM for the macroscopic H-field + Carstensen / Ferreira loss
     post-processing per strand position
   - Igarashi (2017) homogenized complex permeability
3. **Numerical / coupled** (minutes - hours per evaluation):
   - PEEC for inter-strand currents (Rosskopf 2018, ETH coupled approach)
   - Full 3D FEM with every strand meshed (Earth Simulator,
     Kawase 2016)

Each level trades accuracy for speed; the lab convention is to
*always* start with the analytical Dowell / Ferreira number, then
escalate only if the design margin is < 30 % or if a published
measurement disagrees by > 50 %.

## Cross-references

- `radia_mcp.peec.carstensen_ac_copper_loss` -- post-process FE H-field
  with Carstensen / Ferreira per-strand loss formula
- `radia_mcp.motor.hollaus_eddy` -- MSFEM for laminated cores
  (analogous eddy-current PDE in laminations)
- `radia_mcp.ih` -- workpiece SIBC + ESIM
- `radia_mcp.pcb.coil_compensation` -- Litz wire in WPT coils
"""


LITZ_CONSTRUCTION = r"""
# Litz wire construction and geometric parameters

A Litz wire is *not* a simple bundle of N parallel strands -- it is a
hierarchically twisted structure designed so that each strand
statistically occupies every radial position of every bundle level
once per pitch.  This is what makes the bundle-level (Type 1)
external proximity loss small.

## Hierarchy

Standard nomenclature (Sullivan / Umetani 2021):

```
strand               (level 0, single isolated copper wire,
                      diameter d_s, electrically isolated by varnish)
   |
   v  N_1 strands twisted together with pitch p_1
1st bundle           (level 1, "operating bundle")
   |
   v  N_2 first bundles twisted together with pitch p_2
2nd bundle           (level 2)
   |
   v  ... up to typically 5 levels
Litz wire / cable    (top level)
```

Typical industrial Litz spec (Tohoku, Pack, Rubadue):
- d_s = 0.05 - 0.20 mm (40 AWG - 32 AWG)
- N_1 = 5 - 20 strands per 1st bundle
- N_2 = 5 strands per 2nd bundle
- Total strand count 50 - 5000

## Packing factor

The packing factor of a bundle is the ratio of total copper
cross-section to bundle cross-section:

```
beta = N * pi * r_s^2 / (pi * r_Li^2) = N * (d_s / d_Li)^2
```

Realistic values:
- Loose Litz:  beta = 0.3 - 0.4
- Standard:    beta = 0.5 - 0.6
- Compact / served: beta = 0.7 (rare, hard to twist)

A high beta lowers DC resistance per unit cross-section but increases
internal proximity loss (PP_int scales with the bundle radius r_Li
and the strand current density).

## Twist pitch

The pitch p_k at level k is the axial distance over which one full
twist (360 deg of rotation) occurs.  The strand inclination angle is:

```
alpha_k = atan(2 * pi * r_Lk / p_k)
```

where r_Lk is the bundle radius at level k.  For p_k >> r_Lk the
inclination is small (< 5 deg) and the 2D cross-section approximation
(A1.1 in Meng & Biela's taxonomy) holds.  Modern Litz typically has
inclination ~ 1 - 3 deg per level.

The longitudinal eddy current induced by axial flux (the Type 1 / Type
2 bundle-level loss) is suppressed exactly when every strand spends
equal time at every radial position over a length of one pitch -- i.e.
when the cable is "perfectly twisted" (A1.1).

## Strand insulation

Strand isolation is required to prevent inter-strand eddy currents
that would otherwise short-circuit the twisting benefit.  Common
coatings:

| Insulation | Operating temp | Notes                                |
|------------|----------------|--------------------------------------|
| Polyurethane (PUR) | 130 deg C | most common, solderable below 380 C  |
| Polyamide-imide (PAI) | 200 deg C | high-temp Litz                     |
| Polyimide  | 240 deg C      | thinnest film, highest reliability   |

The insulation thickness (typically 5 - 10 % of d_s) directly limits
the packing factor.

## Magnetic plating (lab specialty)

A separate research line places a thin (~ 1 mu m) layer of ferromagnetic
material (Fe-Ni, Co, NiCo) on the *outside* of the copper strand.
This redirects the magnetic flux *around* the conductor (mu_plating *
delta_plating ~ external flux), reducing the surface flux that drives
proximity eddy currents.

References:
- `02_リッツ線/02_proximity_effect/Extending the Operating Distance of
  Inductive Proximity Sensor Using Magnetoplated Wire.pdf`
- 渦電流を考慮したマクロ磁気特性細線の解析 (lab Japanese refs)

## Design knobs (the Stage-1 enumeration for any Litz panel)

The "what is changeable" list for a Litz selection panel:

```
1. d_s        strand diameter (m)
2. N_total    total strand count
3. levels     number of twisting levels (1 - 5)
4. beta_k     packing factor at each level
5. p_k        pitch at each level (m)
6. sigma      strand material conductivity (S/m)
7. mu_r       relative permeability (1 for Cu, > 1 if magnetic-plated)
8. r_Li_max   outer Litz bundle radius (m), often constrained
9. T_op       operating temperature -> sigma(T) correction
```

For a fixed (d_s, beta, levels), the design optimum on N_total has a
well-known plateau:
- Below the optimum, AC R ~ 1/N (DC limit dominates)
- Above the optimum, internal proximity loss dominates and AC R *grows*
  with further N

Sullivan (1999, 2014) and Umetani (2021) provide closed-form
optimum-N formulas; Meng & Biela (2020) show that all four 1D / 2D
models agree on the optimum location to within ~ 10 % for beta < 0.5.
"""


SKIN_PROXIMITY_BASICS = r"""
# Skin and proximity effects -- the underlying physics

## Single round wire skin effect (Belevitch / Kelvin)

For an isolated round wire of radius a carrying a uniform-current
contribution I_skin and immersed in vacuum (no external field), the
quasi-static axial electric field satisfies a Bessel equation.  The
classical result (Kelvin 1889, Belevitch 1977 review) is:

```
Z_skin / R_dc = (k * a) / 2 * J_0(k * a) / J_1(k * a)
```

where `k = (1 - j) / delta` is the complex wave number in the
conductor and `R_dc = 1 / (sigma * pi * a^2)` is the DC resistance per
unit length.

In terms of the normalised radius `gamma_s = a / delta * sqrt(2)`
(used as `gamma_s = d_s / (sqrt(2) * delta)` in Meng & Biela), this
yields the Kelvin form using the real-valued Kelvin functions ber, bei:

```
F_R,S(f) = gamma_s / 2 *
   (ber(gamma_s) * bei'(gamma_s) - bei(gamma_s) * ber'(gamma_s))
   / (ber'(gamma_s)^2 + bei'(gamma_s)^2)
```

and `R_AC = R_dc * F_R,S(f)`.

Asymptotic limits:
- Low frequency (gamma_s << 1):
  `F_R,S ~ 1 + (gamma_s^4) / 48 = 1 + (1/48) * (a/delta)^4`
  (the leading correction is O(f^2))
- High frequency (gamma_s >> 1):
  `F_R,S ~ (1/4) * gamma_s = a / (2 * delta)`
  (corresponds to current confined in skin of thickness delta:
  R_AC = 1 / (sigma * 2 * pi * a * delta))

## Proximity effect (single wire in uniform external field)

When the same isolated round wire is immersed in a *uniform*
time-harmonic external field H_e perpendicular to the wire axis
(no axial current), an eddy-current dipole flows around the wire
generating a loss per unit length:

```
P_p = G_R,S(f) * H_e^2
```

where

```
G_R,S(f) = - 2 * pi * gamma_s / sigma *
   (ber2(gamma_s) * ber'(gamma_s) + bei2(gamma_s) * bei'(gamma_s))
   / (ber2(gamma_s)^2 + bei2(gamma_s)^2)
```

(ber2, bei2 are the second-order Kelvin functions; the sign is
conventional, the product is positive in the dissipative regime.)
This is the "external proximity" loss per wire.

Asymptotic:
- Low frequency: `G_R,S ~ (pi * mu_0^2 * sigma * omega^2 * a^4) / 8`
  (scales with f^2 and a^4)
- High frequency: saturates at `2 * pi * a / (sigma * delta)`

## Orthogonality (A1.3 in Meng & Biela 2020)

Because the skin-effect current density inside a wire is an *even*
function of the polar angle phi (it depends only on r), and the
proximity eddy current is *odd* in phi (it has cos(phi) dependence),
the integral of their product over the cross-section is zero.
Therefore total loss is the sum:

```
P_total = P_skin + P_proximity_internal + P_proximity_external
```

Each term can be computed independently with its own analytical or
numerical formula, and superposed.  This is the foundation of the
Ferreira and Bartoli Litz models.

## Internal vs external proximity (Ferreira 1992)

Inside a *bundle* of strands, each strand "sees" the magnetic field
generated by:
- All *other* strands of the same bundle (internal proximity, H_in)
- The currents in *other* turns / bundles of the device (external
  proximity, H_e)

For a Litz bundle of N strands twisted with packing factor beta and
carrying a *total* current I (so each strand carries I/N), the
internal field follows the Ampere-law radial profile:

```
H_in(r) = (N * (I/N)) * r / (2 * pi * r_Li^2) = I * r / (2 * pi * r_Li^2)
                                                 for r < r_Li
```

The loss per strand from H_in integrated over the bundle cross-section
(Ferreira 1992, eqn 17 of Meng & Biela survey):

```
P_p,int,M2 = N_s * G_R,S(f) * I^2 / (8 * pi^2 * r_Li^2)
```

where N_s is the number of strands per bundle (here N).  This term
scales as I^2 / r_Li^2 -- a *smaller* bundle radius gives *higher*
internal proximity loss.

The external proximity field H_e is the layer-position-dependent
contribution from other turns; for a transformer with N_l turns per
layer in a winding window of breadth b, the n-th layer has:

```
H_e,n = (2*n - 1) / 2 * N_l * I / b
```

so the *outermost* layer (n closest to the core boundary) carries
zero external field, while the innermost (n = N_layer) sees the full
N_l * I / b magnitude squared.

## Why Litz works

Perfect twisting (A1.1) forces every strand to occupy every position
in the bundle for an equal fraction of the cable length.  The
*differential* longitudinal voltage from the linked external flux
therefore averages to zero over one pitch -> no bundle-level loop
current, only the strand-level eddy currents survive.  The
bundle-level "Type 1" loss (which would scale as (r_Li / delta)^4
* f^2 for solid wire) is thereby eliminated.

The remaining "Type 2" losses (strand internal proximity from H_in,
strand skin effect, strand external proximity from H_e) all scale at
the *strand* skin depth gamma_s = d_s / (sqrt(2) * delta), which is
small when d_s << delta.  This is the design rule: pick d_s such that
d_s / delta < 2 at the highest harmonic of interest.
"""


DOWELL_FORMULA = r"""
# Dowell formula (1966) -- 1D layered foil winding

Dowell's classic paper [Dowell 1966, IEEE Proc.] gives a closed-form
AC-to-DC resistance ratio for a *layered foil* winding under the
assumption of a 1D magnetic field varying only across the layer
stack.

## Original formulation

For a foil winding with m layers, each of thickness h, carrying
sinusoidal current at angular frequency omega:

```
R_AC / R_DC = X + (2/3) * (m^2 - 1) * Y

X = xi * (sinh(2*xi) + sin(2*xi)) / (cosh(2*xi) - cos(2*xi))
Y = xi * (sinh(xi) - sin(xi)) / (cosh(xi) + cos(xi))

xi = h / delta
```

Interpretation:
- `X * R_DC` is the *skin-effect* contribution per layer
- `(2/3) * (m^2 - 1) * Y * R_DC` is the *external proximity* (MMF
  building up across layer stack)

At low frequency (xi << 1):
- `X -> 1`, `Y -> xi^4 / 3`
- `R_AC / R_DC -> 1 + (2/9) * (m^2 - 1) * xi^4`

At high frequency (xi >> 1):
- `X -> xi`, `Y -> xi`
- `R_AC / R_DC -> xi * (1 + (2/3) * (m^2 - 1)) = xi * (2*m^2 + 1) / 3`

The proximity term dominates dramatically at high f and high m.

## Adaptation to round Litz (M1 model, Wojda 2010, 2012)

For round Litz wire, Wojda re-derived the Dowell formula by replacing
the foil thickness with the strand-equivalent dimension:

```
Delta_s = (pi/4)^(3/4) * (d_s / delta) * sqrt(d_s / p_s)
```

where p_s is the strand-to-strand pitch (= d_s / sqrt(beta) for a
hexagonal close-packed bundle, with beta the packing factor).  The
single-layer Dowell loss factor is:

```
F_F,S(f) = Delta_s / 2 *
   (sinh(Delta_s) + sin(Delta_s)) / (cosh(Delta_s) - cos(Delta_s))
G_F,S(f) = (b * Delta_s) / (h_s * sigma_eff) *
   (sinh(Delta_s) - sin(Delta_s)) / (cosh(Delta_s) + cos(Delta_s))
```

with `sigma_eff = sigma * eta` (eta = sqrt(pi)/2 * d_s/p_s is the
porosity factor that homogenizes the strand array into an effective
foil of conductivity sigma_eff and thickness h_s).

This is the "M1" model in the Meng & Biela 2020 survey.

## Range of validity

The Dowell formula is *exact* for an infinitely-wide foil with
sinusoidal MMF; for round Litz it is *approximate*:
- Accurate to ~ 10 % when d_s < delta and beta < 0.5
- Increasingly inaccurate when d_s > delta (the foil approximation
  cannot capture intra-strand current density variation)
- M1 systematically *under-predicts* internal proximity loss
  vs. Ferreira M2 (Bessel) by ~ 5 - 20 % for typical beta = 0.5

For careful design work, prefer Ferreira M2 / Bartoli M3 / Tourkhani
M4 (see `ferreira_strands`).  Use Dowell M1 only for hand calculations
or sanity checks.

## Lab implementation

Dowell is implemented for *rectangular* solid conductors (foil
winding) in:

```
src/core/rad_peec_surface_impedance.cpp
```

as the high-frequency SIBC for PEEC + FEM coupling.  It is invoked
by `calc_inductance.py --coil-solver peec` when a rectangular
cross-section is detected.

## Reference

- Dowell, P.L., "Effects of eddy currents in transformer windings",
  Proc. IEE, 113(8):1387-1394, 1966.
- Wojda, R.P. & Kazimierczuk, M.K., "Winding resistance of litz-wire
  and multi-strand inductors", IET Power Electron. 5(2):257-268, 2012.
- Meng, Q. & Biela, J., "Survey and Comparison of 1D/2D analytical
  Models of HF Losses in Litz Wire", EPE-ECCE Europe 2020,
  DOI 10.23919/EPE20ECCEEurope43536.2020.9215651.
"""


FERREIRA_STRANDS = r"""
# Ferreira 2D Bessel model for round strands (M2/M3/M4)

Ferreira [1990, 1992, 1994] derived the exact 2D solution for a round
conductor in an external uniform magnetic field, then summed
contributions strand-by-strand.  Three lineage extensions exist
(Meng & Biela 2020 notation):

- **M2 Ferreira**: original Bessel-based model
- **M3 Bartoli** [Bartoli 1995]: + internal porosity factor eta_in
- **M4 Tourkhani** [Tourkhani 2008]: + extended external field model

## M2 -- Ferreira (1992)

Two loss types per layer (n = layer index from inside to outside,
N_l = turns per layer, b = winding window breadth):

Skin effect (same per strand, independent of position):
```
P_skin,M2 = R_DC * F_R,S(f) * I^2 / 2
```

External proximity:
```
P_p,ext,M2 = N_s * G_R,S(f) * H_e,M2^2
H_e,M2     = (2*n - 1) / 2 * N_l * I / b
```

Internal proximity:
```
P_p,int,M2 = N_s * G_R,S(f) * I^2 / (8 * pi^2 * r_Li^2)
```

The Kelvin-function factors F_R,S and G_R,S are defined in
`skin_proximity_basics`.

## M3 -- Bartoli (1995)

Identical to M2 for skin loss.  Adds an *internal porosity factor*:

```
eta_in = sqrt(pi) * d_s / (2 * t_s)
```

where t_s is the strand-to-strand pitch.  The internal proximity
formula becomes:

```
P_p,int,M3 = N_s * eta_in^2 * G_R,S(f) * I^2 / (8 * pi^2 * r_Li^2)
            = N_s * G_R,S(f) * I^2 / (8 * pi^2 * r_b^2)
```

(using the relation `N_s * t_s^2 = beta * pi * r_b^2`).  In effect M3
replaces the geometric bundle radius r_Li with the "true" bundle
radius r_b that excludes the inter-strand insulation gap.

The external proximity uses the foil-equivalent breadth t_Li:

```
H_e,M3 = (2*n - 1) / 2 * N_l * I / (N_l * t_Li)
       = (2*n - 1) * I / (2 * t_Li)
```

This emphasises that, in the foil-equivalent geometry, the layer
breadth is the *strand* extent, not the *winding window* breadth.

## M4 -- Tourkhani (2008)

Same skin and internal loss as M2 / M3.  External H_e is computed via
a refined formula accounting for the parabolic field profile inside
the layer:

```
H_e,M4 = N_l * I / b_l * sqrt((2*n-1)^2 / 4 + 1/16)^(1/2)
```

This is a small correction that becomes important only for layer
thickness h_l comparable to the layer height d_Li (very dense
winding).

## Which model when?

Meng & Biela's benchmark study (2020 EPE-ECCE) compared M1-M4 against
2D FEM for a reference Litz parameter set:
- N_s = 60, d_s = 100 mu m, p_s = 5 mm, r_Li = 0.7 mm, m = 4 layers
- Frequency band 10 kHz - 5 MHz

Findings (paraphrased):

| Model | Error vs 2D FEM | Best when             | Worst when            |
|-------|-----------------|-----------------------|-----------------------|
| M1 Wojda/Dowell | 5 - 25 % | d_s << delta, beta < 0.4 | dense packing, gamma_s > 2 |
| M2 Ferreira     | 2 - 10 % | most cases             | very high m (> 10)    |
| M3 Bartoli      | 2 - 8 %  | high beta (> 0.5)      | very low f            |
| M4 Tourkhani    | 2 - 7 %  | dense layer (h_l/d_Li > 0.8) | sparse packing |

Recommendation: **M3 Bartoli is the production choice** for the
lab's WPT and IH design work; M2 Ferreira is the textbook reference
for tutorials and verification.

## Reference

- Ferreira, J.A., "Improved analytical modeling of conductive
  losses in magnetic components", IEEE Trans. Power Electron.
  9(1):127-131, 1994.
- Bartoli, M., Reatti, A., Kazimierczuk, M.K., "Modeling iron-powder
  inductors at high frequencies", IEEE IAS Annu. Mtg., 1994.
- Tourkhani, F. & Viarouge, P., "Accurate analytical model of winding
  losses in round Litz wire windings", IEEE Trans. Magn. 37(1):538,
  2001.
- Sullivan, C.R., "Optimal choice for number of strands in a Litz-wire
  transformer winding", IEEE Trans. Power Electron. 14(2):283-291,
  1999.
"""


UMETANI_MULTI_LEVEL = r"""
# Multi-level twisted Litz -- measurable parameter model (Umetani 2021)

Most 1D / 2D analytical models (M1-M4) assume a *single-level* twist:
N strands twisted together in one go.  Industrial Litz wires with
high strand count (> 100) are nearly always *multi-level*: strands ->
1st bundle -> 2nd bundle -> ... -> cable.

Sullivan & Zhang [Sullivan 2014] extended the analytical theory to
multi-level Litz, but the formula requires *packing factors at every
level* (beta_1, beta_2, ...), which are not directly measurable.

Umetani et al. [Umetani 2021, IEEE TIA 57(3):2407] re-derived the
multi-level Litz copper loss model using *only parameters measurable
with basic testing instruments*:

- d_s          strand diameter (caliper)
- d_Lk         bundle outer diameter at each level (caliper)
- p_k          pitch at each level (visual count over a known length)
- N_k          number of sub-bundles at level k (visual count)
- DC R         DC resistance (4-wire ohm-meter)
- sigma        from material spec sheet

## Hierarchical formulation

The Litz wire has L levels of twisting.  Define:

```
r_s    = strand radius (level 0)
r_Lk   = bundle radius at level k
n_k    = number of sub-bundles twisted at level k
N_total = product of all n_k (= total strand count)
```

The copper loss per unit length of the *total Litz wire* is built up
inductively:

1. **Strand level (level 0)**: each strand sees its own skin effect
   (F_R,S applied to I/N_total) and its own internal field H_strand
   from the surrounding strands of the same 1st bundle.

2. **1st bundle level**: each 1st bundle of N_1 strands has an
   effective bundle current I_bundle = N_1 * (I/N_total) and an
   effective bundle radius r_L1.  This bundle is exposed to the
   field H_bundle from neighbouring 1st bundles within the 2nd
   bundle.

3. ... continue up the hierarchy ...

L. **Cable level**: the entire cable carries I and is exposed to the
   external field H_e from the device (e.g. core MMF).

## Specifics of the Umetani simplification

The key contribution is the *approximation* that:
- The lowest-level twisting (strand -> 1st bundle) has *many* strands
  (typically > 5), so the bundle behaves like an isotropic Litz bundle
  (M3 Bartoli applies).
- The *higher* levels of twisting have *few* sub-bundles (typically
  <= 5), so the discrete number of sub-bundles cannot be averaged out
  -- each sub-bundle's contribution must be summed explicitly.

This leads to a hybrid model:
- M3 Bartoli at level 1
- Explicit summation over the n_k sub-bundles at levels 2, 3, ...

The model also includes the *twist pitch* effect (strand inclination
alpha_k) on the effective skin / proximity parameters; this raises
the AC R by typically 2 - 8 % vs the perfect-twist limit.

## Conditions for validity

Umetani lists 5 assumptions:
1. Strand surface is non-conducting (isolated)
2. Twisting is in 2 or more levels
3. Lowest-level twist > 5 strands per bundle
4. Higher levels twist <= 5 sub-bundles per bundle
5. Litz length >> max twist pitch (so each level's twist completes
   many rotations)

## Lab applicability

The Umetani model is the recommended production model for the lab's
high-frequency transformer design (resonant DC-DC, MFTM, EV WPT
ferrite-core coils) when the Litz has 3+ levels.  For 1-2 level Litz,
M3 Bartoli is adequate and faster.

## Reference

- Umetani, K., Kawahara, S., Acero, J., Sarnago, H., Lucia, O., Hiraki,
  E., "Analytical Formulation of Copper Loss of Litz Wire With
  Multiple Levels of Twisting Using Measurable Parameters", IEEE Trans.
  Industry Applications, 57(3):2407-2420, 2021,
  DOI 10.1109/TIA.2021.3063993.
- Sullivan, C.R. & Zhang, R.Y., "Simplified design method for Litz
  wire", IEEE APEC 2014, pp. 2667-2674.
"""


HOMOGENIZATION = r"""
# Homogenization-based FEM for Litz wire / multi-turn coils

Native FEM of every strand in a Litz wire is intractable (a 1000-strand
bundle with 10 turns = 10^4 strands, each meshed with ~ 10 elements
across delta -> ~ 10^6 elements per cross-section).  Homogenization
replaces the bundle with an *effective* anisotropic material.

## Igarashi (2017) complex-permeability homogenization

The key insight [Igarashi 2017, IEEE TMAG 53(1):7400107]: a single
round conductor of radius a, carrying a uniform current and immersed
in a uniform time-harmonic external B_0, behaves *macroscopically*
like a magnetised material with complex permeability:

```
mu_eff(omega) = mu_0 * [1 + 2 * (Q(k*a) - mu_r) / (mu_r * Q(k*a) + 2)]
```

where `k = (1-j)/delta`, mu_r is the strand material's relative
permeability (= 1 for Cu), and Q is a Bessel-function ratio:

```
Q(z) = 1 - 2 * J_1(z) / (z * J_0(z))
```

Asymptotic limits:
- Low frequency: `mu_eff -> mu_0 * mu_r` (no eddy effect)
- High frequency: `mu_eff -> mu_0 * (1 - 2 / (mu_r * (k*a)^2))` ->
  perfect diamagnet (mu_eff -> 0)

The imaginary part of mu_eff is the proximity loss.

For a multi-turn coil with packing factor beta, the *homogenized*
permeability over the coil cross-section follows the **Ollendorff
formula** (homogenization of a periodic array of cylinders):

```
mu_homog = mu_air * [1 + 2 * beta * (mu_eff - mu_air) /
                       (mu_air + mu_eff - beta * (mu_eff - mu_air))]
```

The coil is then meshed as a uniform diamagnetic block with mu_homog
and the (frequency-independent) DC resistance R_DC + (frequency
correction from skin Z_skin).

## Verification

Igarashi's paper compares the homogenized model against:
1. Native FEM with every strand meshed
2. Impedance analyzer measurement

Agreement is within 2 % up to f * delta < d_s (the strand skin depth
condition).  Above that, the Ollendorff homogenization starts to fail
because intra-strand current density variation matters.

## Lab variants

The lab has implemented multiple variants for different applications:

- 均質化法を用いた磁性メッキ巻線の有限要素解析 -- magnetic-plated
  Litz, with mu_r > 1 in the Igarashi formula
- 不整合メッシュを用いた有限要素解析のための均質化法 -- non-conforming
  mesh for the homogenized coil region
- 異方性媒質を含む不整合有限要素メッシュ均質化法 -- anisotropic
  permeability tensor (mu_par along strand axis, mu_perp transverse)
- Voxel-based FEM with homogenization for 3D Litz characterisation

References (all in W:/.../02_リッツ線/99_misc/01_均質化_homogenization/):
- Eddy Current Analysis of Litz Wire Using Homogenization-Based FEM
- Fast 3-D Analysis of Eddy Current in Litz Wire Using Integral Equation
- 均質化法を用いた磁性メッキ巻線の有限要素解析

## Application to Radia's PEEC + FEM coupling

In the lab's PEEC + FEM workflow (calc_inductance.py --coil-solver peec
--vol workpiece.vol), the coil itself is represented by *filaments*
(no mesh), so coil-side homogenization is unnecessary.  However, when
a *coil-meshed* approach is needed (e.g. for tightly-coupled WPT pads
where the leakage flux through the coil matters), the Igarashi
homogenization should be applied via a Materials block in the .vol
file with the complex permeability evaluated at the operating
frequency.

A planned implementation route: extend `MatLin(mu_r)` to support
complex `mu_r` and use it inside the FEM solve.

## Reference

- Igarashi, H., "Semi-Analytical Approach for Finite-Element Analysis
  of Multi-Turn Coil Considering Skin and Proximity Effects",
  IEEE Trans. Magn. 53(1):7400107, 2017.
- Moreau, O., Popiel, L., Pages, J.L., "Proximity losses computation
  with a 2D complex permeability modelling", IEEE Trans. Magn.
  34(5):3616-3619, 1998.
- Ollendorff, F., "Magnetostatik der Massekerne", Arch. Elektrotech.
  25:436-447, 1931.
"""


MULTICONDUCTOR_TL = r"""
# Multiconductor transmission line theory

A multiconductor transmission line (MTL) is a system of N+1 parallel
conductors (one reference / ground + N "signal" lines) carrying
quasi-TEM waves.  The per-unit-length R, L, G, C matrices describe
the system.

## Telegrapher equations (matrix form)

For N signal conductors:

```
d V(z, omega) / d z = - (R(omega) + j*omega*L(omega)) * I(z, omega)
d I(z, omega) / d z = - (G(omega) + j*omega*C(omega)) * V(z, omega)
```

V, I are N x 1 phasor vectors; R, L, G, C are N x N symmetric matrices.

In the time domain (for transient analysis):

```
d V/d z = - R * I - L * d I/d t
d I/d z = - G * V - C * d V/d t
```

## P / L coefficient formulation (Toki & Sato 2011)

Toki & Sato [J. Phys. Soc. Japan 81(1):014202, 2012] re-derived MTL
theory directly from the retarded potentials of Maxwell's equations
in Lorenz gauge.  They obtained:

- **Coefficients of potential** P_ij: how charge q_j on conductor j
  contributes to the scalar potential V at conductor i
- **Coefficients of inductance** L_ij: how current I_j contributes to
  the vector potential A at conductor i

Both P and L are derived from the same Neumann-type double-integral
formula:

```
P_ij = (1 / (4 * pi * epsilon_0)) * (1 / l) * int int [1 / |r - r'|] dS dS'
L_ij = (mu_0 / (4 * pi))         * (1 / l) * int int [1 / |r - r'|] dl dl'
```

over the conductor surfaces (P) or centerlines (L).  Their ratio is:

```
L_ij * P_ij = (1 / c^2) -> sqrt(P_ij / L_ij) = c (speed of light)
```

This means TEM waves on an MTL system propagate at c (vacuum case)
*regardless of the conductor cross-section*, as long as the medium
is vacuum.

The capacitance matrix C is the inverse of the matrix-of-potential
coefficients:

```
C = epsilon_0 * P^{-1}
```

For dielectric-filled regions, replace epsilon_0 with eps(r).

## Three-conductor system and propagation modes

Toki & Sato emphasize that a *three-conductor* MTL (two signal + one
ground) is the minimum that exhibits the full mode structure:

| Mode    | Definition                  | Physical meaning              |
|---------|-----------------------------|-------------------------------|
| Normal  | I_1 = -I_2, I_g = 0         | Differential, "useful" signal |
| Common  | I_1 = I_2, I_g = -2*I_1     | Common-mode, drives ground    |
| Antenna | I_1 = I_2 = I_g (or pattern with EM radiation) | Radiates to space  |

In a symmetric arrangement (lines 1 and 2 placed symmetrically about
the ground line), the *normal* mode decouples from the common and
antenna modes -- no mode conversion at lossless discontinuities.
This is the design principle for differential signalling (LVDS, USB,
HDMI, Ethernet) on PCBs.

In an *asymmetric* arrangement, mode mixing occurs and a portion of
the normal-mode signal converts to common mode at every
discontinuity -- the main source of EMI and crosstalk in real PCBs.

## Tilde potentials and radiation

Toki & Sato split the retarded scalar / vector potentials into:
- A *singular* part -> gives rise to P and L coefficients (TEM mode)
- A *regular* (tilde) part `V_tilde`, `A_tilde` -> represents
  electromagnetic radiation away from the line

For an isolated single-conductor line, the tilde potentials are the
source of the antenna radiation -- this directly couples MTL theory
to antenna theory and the EMF method (radiation resistance).

This formulation is mathematically clean and is the reference
framework for the lab's PCB EMI calculations (transmission line
trace ↔ near-field probe coupling).

## Parameter extraction methods

For a given cross-section geometry, R, L, G, C can be computed by:

1. **Analytical** (microstrip, coaxial): closed-form per textbook
   (Pozar, Paul "Analysis of Multiconductor Transmission Lines").
2. **2D quasi-static FEM**: separate solves
   - Electrostatic Laplace with Dirichlet V on each conductor -> C
   - Magnetostatic Laplace with surface currents -> L_ext
3. **2D MoM / integral equation** (Matsuki & Matsushima 2012):
   non-uniform skin-depth-adapted grid in the conductor cross-section
   -> internal impedance including skin and proximity
4. **PEEC** (lab method): 3D extension via partial inductance / capacitance
5. **BEM SIBC** (Huangfu & Di Rienzo 2020): boundary-element formulation
   with high-order surface impedance for frequency-dependent loss

## Reference

- Paul, C.R., "Analysis of Multiconductor Transmission Lines",
  2nd ed., Wiley-IEEE, 2007.
- Toki, H. & Sato, K., "Multiconductor Transmission-Line Theory with
  Electromagnetic Radiation", J. Phys. Soc. Japan 81(1):014202, 2012.
- Toki, H. & Sato, K., "Three Conductor Transmission Line Theory and
  Origin of Electromagnetic Radiation and Noise" (companion paper).
"""


FEM_THIN_WIRE = r"""
# FEM and MoM for thin-wire / rectangular TL impedance

For transmission lines with rectangular conductors (PCB traces, IC
interconnect, busbars) or arbitrary cross-section bundles, the
per-unit-length impedance must be computed by 2D quasi-static FEM or
MoM that resolves the *skin layer* on every conductor surface.

## Matsuki & Matsushima (2012) non-uniform grid MoM

For a single rectangular conductor of width w and thickness t in
vacuum, the axial current density satisfies a 2D Helmholtz equation:

```
[d^2/dx^2 + d^2/dy^2] A_z = - mu_0 * J_z   (inside conductor)
[d^2/dx^2 + d^2/dy^2 + beta_0^2] A_z = 0   (outside, beta_0 = omega/c)
```

with the constitutive relation `J_z = sigma * (-j*omega*A_z - d V/d z)`
inside.  In integral-equation form (Matsushima):

```
A_z(r) = (mu_0 / (4j)) * int_S Jz(r') H_0^(2)(beta_0 |r - r'|) dS'
```

The Method of Moments (MoM) discretizes S with constant basis
functions.  For accuracy at high frequency (delta << t), Matsuki &
Matsushima propose a *non-uniform grid* designed around the skin
effect:

```
grid spacing h_n = h_min * (h_max / h_min)^(n / N_grid)
                   for n = 0, ..., N_grid - 1
```

with `h_min ~ delta / 3` near the surface and `h_max ~ w / 8` in the
bulk.  This geometric grading captures the exponential decay of J_z
near the surface without wasting nodes in the bulk.

Results (Matsuki & Matsushima 2012, PIER M Vol. 23):
- Convergence to < 1 % error with ~ 50 unknowns per side at f/delta ~ 100
- Conformal-mapping analytical solution recovered to 0.1 % in HF limit

## Extension to multiconductor + lossy substrate (Matsuki 2012b)

For VLSI interconnect, the conductors sit above a *lossy semi-infinite
substrate* (silicon, sigma_s, eps_s).  Matsuki & Matsushima [2012, PIER B
Vol. 43] extend the formulation by adding a *Sommerfeld-integral
correction* to the free-space Hankel-function kernel:

```
g_substrate(r, r') = g_free(r, r') + g_correction(r, r'; sigma_s, eps_s)
```

The correction is evaluated analytically in the quasi-static regime
using Struve and Neumann special functions, avoiding the costly
numerical Sommerfeld integral.

Energy-conservation check:
```
R_total = R_internal + R_external (substrate eddy current loss)
```
with relative error < 1 % across f = 1 MHz - 100 GHz.

## Higher-order asymptotic boundary condition (Khebir, Kouki, Mittra 1990)

For open-region FEM with multilayered dielectric MTL, the artificial
outer boundary must absorb the outward-radiating wave.  Khebir et al.
[IEEE TMTT 38(10):1433, 1990] derived a higher-order asymptotic
boundary condition (HABC) from the general solution of Laplace's
equation in polar coordinates:

```
B_m u = (sum_{n_i} (d/d_rho + n_i / rho)) u  (m-th order operator)
```

The third-order operator is:

```
B_3 u = u_rho + ((n+1)/rho) * u_rho + ((n+1)/rho^2) * u
```

Applied at a circular artificial boundary, B_3 yields significantly
more accurate FE solutions than a perfectly-conducting shield at the
same location, while preserving the sparsity of the FE matrix.

This ABC is the lab's recommended choice for open-region 2D MTL FEM
in NGSolve.

## FDTD MTL with SIBC (Huangfu & Di Rienzo 2020)

For *time-domain* MTL simulation with frequency-dependent losses, the
classical Paul FDTD method assumes well-separated conductors and only
the skin effect.  Huangfu & Di Rienzo [IEEE TMAG 56(1):7502304, 2020]
extend this with:

1. A *partial* internal impedance matrix Z_i,p(s) computed via BEM
   with a high-order SIBC (Yuferev, Ida 2018) capturing both skin
   *and* proximity effects.
2. Time-domain convolution `Z_i(t) * I(z, t)` via inverse Laplace
   transform of Z_i(s), using vector-fitting rational approximation
   to make the convolution O(N) per timestep.
3. External inductance L_e from the standard low-frequency formulation
   (frequency-independent).

The resulting FDTD scheme handles transient signals (e.g. PEEC + power
electronics) with full skin + proximity dispersion.

## Lab applicability

For the lab's PCB EMI and busbar work:
- Use Matsuki MoM (or PEEC) for *per-unit-length* impedance extraction
  at the cross-section level
- Use ngsolve.bem with shielded HCurl + Compact HX preconditioner for
  full 3D MTL simulation
- Use Huangfu FDTD when transient response (e.g. switching surge) is
  needed

The Trefftz-like method of Giussani et al. [COMPUMAG 2019] for
parallel magnetic + conductive wires (submarine power cable armor) is
an alternative when the wires are *round* and the field is
*transverse*: no mesh required, basis functions are exact Bessel
solutions inside each wire and decaying multipoles outside.

## Reference

- Matsuki, M. & Matsushima, A., "Improved numerical method for
  computing internal impedance of a rectangular conductor and
  discussions of its high frequency behavior", PIER M Vol. 23, pp.
  139-152, 2012.
- Matsuki, M. & Matsushima, A., "Efficient impedance computation for
  multiconductor transmission lines of rectangular cross section",
  PIER B Vol. 43, pp. 373-391, 2012.
- Khebir, A., Kouki, A.B., Mittra, R., "Higher Order Asymptotic
  Boundary Condition for the Finite Element Modeling of Two-
  Dimensional Transmission Line Structures", IEEE TMTT 38(10):1433,
  1990.
- Huangfu, Y., Di Rienzo, L., Wang, S., "FDTD Formulation Based on
  High-Order Surface Impedance Boundary Conditions for
  Frequency-Dependent Lossy Multi-Conductor Transmission Lines",
  IEEE Trans. Magn. 56(1):7502304, 2020.
- Giussani, L., Di Rienzo, L., Bechis, M., de Falco, C., "MoM
  Computation of Losses of an Array of Parallel Magnetic and
  Conductive Wires in a Transverse Magneto-quasi-static Field",
  COMPUMAG 2019, paper 1902020.
"""


ABC_PML = r"""
# Absorbing boundary conditions for open-region MTL

Open-region 2D MTL FEM requires terminating the artificial outer
boundary with an absorbing boundary condition (ABC) so that the
solution mimics an unbounded domain.  Three approaches are commonly
used in the lab.

## Asymptotic boundary conditions (Khebir et al. 1990)

For the *quasi-TEM* mode (no propagating wave, only decaying near
field), the FE potential at the outer boundary satisfies an
asymptotic series in inverse powers of rho:

```
u(rho, phi) = sum_{n=0}^{infty} c_n(phi) / rho^n
```

Khebir, Kouki, Mittra [TMTT 1990] derived second-order and
third-order operators that *annihilate* the first few terms of this
series exactly:

```
B_2 u = u_rho + (n+1)/rho * u     (annihilates 1/rho^0)
B_3 u = u_rho_rho + (2n+3)/rho * u_rho + (n+1)(n+2)/rho^2 * u
```

Applied at the FE truncation boundary, these yield:
- B_2 absorbs the dominant 1/rho monopole tail
- B_3 absorbs monopole + dipole (much more accurate for
  multi-conductor structures with non-zero net charge)

Crucially these ABCs are *local* (no integral kernel) and *preserve
the sparsity* of the FE matrix.  Their accuracy is comparable to a
perfectly-matched layer (PML) at a fraction of the implementation
cost.

## PML-like ABC for transmission line termination (Lauer et al. 1999)

For *time-domain* FDTD or FE simulation of MTL absorbers (used to
match the load end of a line and suppress reflection back into the
simulation domain), a Perfectly-Matched-Layer-like construction is
appropriate.

Lauer, Wien, Waldow, Wolff [EMC Zurich 1999] propose a 1D equivalent
circuit (EC) model for layered-media MTL absorbers (parallel-plate
TL with layered fill) that can be analytically optimized for minimum
reflection at the design frequency:

```
Type A: air-filled
  L_0 = mu_0 * h / w,  C_0 = epsilon_0 * w / h,  Y_c0 = sqrt(C_0/L_0)
Type B: partial dielectric (substrate eps_rd)
  L_B = L_0,  C_B = ((1+eps_rd)/2) * C_0,  Y_cB = Y_c0 * sqrt((1+eps_rd)/2)
```

The discrete conductivity profile sigma(z) of the lossy backing
layers is then chosen to satisfy:

```
Y(z) = j*omega*C(z) + sigma(z)*w/h
ABS(Gamma) -> minimum at f_0
```

For a 2-cell microstrip TL absorber on FR-4 (eps_r = 12.9), the
reduction in S_11 is from -20 dB (uniform profile) to better than
-60 dB (optimized profile) -- a 3-cell PML at the line termination
acts equivalent to a 30-cell uniform absorber.

This is a useful pattern for any FE simulation of *trace -> antenna*
or *trace -> radiator* coupling where the conductor leaves the
simulation domain.

## Reference

- Khebir, A., Kouki, A.B., Mittra, R., "Higher Order Asymptotic
  Boundary Condition for the Finite Element Modeling of Two-
  Dimensional Transmission Line Structures", IEEE Trans. MTT
  38(10):1433-1438, 1990.
- Lauer, A., Wien, A., Waldow, P., Wolff, I., "Enhanced PML-like
  ABCs for layered media transmission line termination", IEEE EMC
  Zurich 1999 Proceedings, pp. ...
- Berenger, J.P., "A perfectly matched layer for the absorption of
  electromagnetic waves", J. Comput. Phys. 114:185-200, 1994.
"""


TLM_METHOD = r"""
# Transmission-line modelling (TLM) method

The Transmission-Line Modelling method (also called Transmission-Line
Matrix method, TLM) is a *numerical technique* that solves Maxwell's
equations on a mesh of *interconnected two-wire transmission lines*.
It exploits the field-circuit equivalence of Maxwell's equations and
distributed transmission-line equations.

## History and motivation

- Schwinger (1940s) introduced the equivalent circuit method for
  waveguide discontinuities.
- Johns (1971) generalized to time-domain scattering problems on a
  3D mesh of TLs, founding modern TLM.
- Sadiku & Cagba [IEEE TCAS 37(8):991, 1990] published a pedagogical
  introduction making TLM accessible to engineers without specialised
  numerical methods background.

## Principle

Each *node* of the TLM mesh is a junction of N (typically 4 in 2D, 12
in 3D) transmission lines whose characteristic admittance Y_c is
chosen to satisfy the local Maxwell discretization.  Voltage pulses
propagate along each TL with delay equal to one timestep `Delta_t`;
at each node, an *S-matrix scattering* event redistributes incoming
pulses to outgoing pulses according to local material properties
(epsilon, mu, sigma).

The basic mesh recursion:

```
V_out[node, t+Delta_t] = S(node) * V_in[node, t]
V_in[neighbour, t+Delta_t] = V_out[node, t+Delta_t]  (propagation)
```

The S-matrix at each node is derived from current continuity (KCL)
and the local circuit element (R, L, G, C).

## Telegrapher equation recovery

For a single 2-conductor TL of length Delta_z with R, L, G, C per
unit length:

```
dV/dz = - R * I - L * dI/dt
dI/dz = - G * V - C * dV/dt
```

After eliminating I:

```
d^2 V/dz^2 = LC * d^2 V/dt^2 + (RC + GL) * dV/dt + RG * V
```

This is the *damped wave equation*, and its lossy-dielectric Maxwell
analogue:

```
d^2 E/dx^2 = mu * epsilon * d^2 E/dt^2 + mu * sigma * dE/dt
```

is recovered by identifying R <-> 0, L <-> mu, G <-> sigma, C <->
epsilon.  Any local material can be modelled by adjusting (R, L, G, C)
of the equivalent TL stubs at each node.

## When to use TLM vs FDTD vs FE

| Feature                | TLM                | FDTD               | FE (NGSolve)       |
|------------------------|--------------------|--------------------|--------------------|
| Time-domain native     | yes                | yes                | no (requires Newmark) |
| Frequency-domain       | post-FFT           | post-FFT           | direct             |
| Arbitrary materials    | trivial per-node   | trivial per-cell   | requires assembly  |
| Unstructured mesh      | hard               | hard               | trivial            |
| Stability              | unconditional      | conditional (CFL)  | unconditional (implicit) |
| Programming effort     | low                | low                | moderate           |
| 3D large-scale         | competitive        | competitive        | best (with HACApK) |
| Lab use                | rare (legacy code) | EMC scans          | production         |

The lab does not use TLM directly, but the method is worth knowing
for two reasons:
1. It is the conceptual underpinning of PEEC (PEEC = TLM with
   lumped-element nodes obtained from partial inductance / capacitance
   integrals, no distributed wave propagation).
2. Vendor EMC tools (CST Microwave Studio, IMST EMPIRE) use TLM as
   their core 3D EM solver; understanding TLM helps in interpreting
   their results.

## Reference

- Johns, P.B., Beurle, R.L., "Numerical solution of 2-dimensional
  scattering problems using a transmission-line matrix", Proc. IEE
  118(9):1203-1208, 1971.
- Sadiku, M.N.O. & Cagba, L., "A Simple Introduction to the
  Transmission-Line Modeling", IEEE Trans. Circuits Syst.
  37(8):991-999, 1990.
- Hoefer, W.J.R., "The transmission-line matrix method -- theory and
  applications", IEEE Trans. MTT 33(10):882-893, 1985.
"""


ROSSKOPF_COUPLED = r"""
# Coupled FEM + PEEC for inductive components with Litz (Rosskopf 2018)

Rosskopf's PhD thesis [Erlangen-Nuernberg, 2018] developed an
end-to-end coupled numerical approach for frequency-dependent
inductor / transformer loss prediction with hierarchically-twisted
Litz wire.  The thesis is the most comprehensive published reference
on coupled FEM + PEEC for Litz.

## The coupled approach (high-level)

The simulation pipeline is split into two sequential stages:

### Stage 1 -- System-level FEM (macro-scale)

The complete inductive component (winding, core, air gaps, surrounding
parts) is meshed *coarsely* with the Litz wire replaced by a
*homogenised* equivalent conductor having uniform current density.
A 3D FE simulation of the system gives the magnetic field
distribution H(r) throughout the device, accounting for:

- Geometrical effects (core saturation curve, air-gap fringing)
- Material magnetic properties (mu_r, mu_complex)
- Interaction with surrounding shielding / mounting hardware

### Stage 2 -- Strand-level PEEC (micro-scale)

Along the length of the conductor, a *high number* (~ 50 - 500) of
*2D cross-section cuts* are extracted from the FE solution.  At each
cut, the local *external* magnetic field H_e(z) at the conductor
position is recorded.

For each cross-section cut, a PEEC simulation of the Litz wire
*structure* (every strand, every twisting level) is run with H_e(z)
as the imposed external field.  PEEC returns:
- Per-strand current distribution at this z
- Local AC resistance R_AC(z) and power dissipation

Integrating R_AC(z) over the conductor length gives the *total*
frequency-dependent winding loss P(f).

## Implementation specifics

Rosskopf uses:
- FEM solver: COMSOL or ANSYS Maxwell 3D for Stage 1
- PEEC solver: custom in-house code (Erlangen Fraunhofer IISB) for
  Stage 2, with the partial inductance / capacitance matrices
  precomputed once per Litz cross-section template

Key optimisations:
- The PEEC matrices depend only on the Litz *cross-section geometry*
  (strand positions, packing, twist pitch).  For a given Litz spec,
  these are computed once and reused for every cross-section cut.
- The Stage-2 PEEC solves at every cut share the same matrix; only
  the RHS (from H_e(z)) changes.  -> precomputed LU.

## Validation

The thesis validates against experiments on:
- Single-turn air coils with various Litz specs (50 kHz - 1 MHz)
- Multi-turn air-core inductors
- Transformer windings on ferrite cores

Agreement up to 1 MHz is within 5 % for air-core devices.  For
ferrite-core devices, additional core losses degrade the agreement
to 10 - 20 % (core loss model is the limiting factor, not the Litz
model).

## Lab applicability

The Rosskopf approach is the *gold standard* for the lab's
high-frequency transformer design.  The 2-stage decomposition matches
the Radia 4-Layer Architecture perfectly:

- Stage 1 -> Radia FEM (calc_fem_coilmesh.py with homogenized Litz)
- Stage 2 -> Radia PEEC (calc_inductance.py --coil-solver peec)

The lab has not yet implemented the *coupling* (Stage 2 driven by
Stage 1 H_e), but it is the natural next step for the WPT / IH coil
optimisation.

## Reference

- Rosskopf, A.H., "Calculation of frequency dependent power losses
  in inductive systems with litz wire conductors by a coupled
  numeric approach", PhD thesis, Friedrich-Alexander-Universitaet
  Erlangen-Nuernberg, 2018.
  PDF at W:/.../02_リッツ線/01_AC_loss/.
"""


AC_LOSS_MEASUREMENT = r"""
# AC loss measurement protocol for Litz / multiconductor

Validating an AC loss prediction (analytical, semi-analytical, or
coupled-numerical) against measurement requires careful experimental
methodology, because the AC resistance is typically tens to hundreds
of times smaller than the inductive reactance at the operating
frequency.

## Equipment

For the lab's typical 1 kHz - 1 MHz band:

| Instrument            | Use                                  | Lab unit               |
|-----------------------|--------------------------------------|------------------------|
| LCR meter             | 1 kHz - 1 MHz impedance              | HP 4284A, Keysight E4980A |
| Impedance analyzer    | 100 Hz - 100 MHz Z spectrum          | HP 4194A, Keysight E4990A |
| Network analyzer + ZA | up to 1 GHz with ZA fixture          | Keysight E5061B + 16047 |
| 4-wire ohm-meter      | DC R reference                       | Keithley 2750, Yokogawa GP10 |

## Sample preparation

- **Terminations**: solder Litz ends to *large* copper pads (cm-scale
  patches) so termination resistance << device R.
- **DC R reference**: measure with 4-wire Kelvin connection, log the
  exact sample temperature (Cu rho/d T is +0.39 %/K -- a 10 K
  drift = 4 % R error).
- **Coil geometry**: for a coil sample, fix the geometry mechanically
  (no flexing during the scan); flexing changes the air-gap and
  shifts the inductive reactance, which projects into the R via
  the LCR meter's R - X separation.

## Measurement procedure

1. Calibrate the analyzer with the recommended OPEN / SHORT / LOAD
   standards at the test cable end where the sample connects.
2. *Zero* the test fixture (e.g. with the 16047A or 16089E Kelvin
   clip), then connect the sample.
3. Sweep frequency in log spacing (typically 10 points/decade) from
   100 Hz to the highest frequency of interest.
4. Record R_meas(f), L_meas(f), and the corresponding Q = omega*L/R.
   If Q > 100, the LCR meter's R measurement may be unreliable -- use
   a Network Analyzer S-parameter measurement instead.
5. Subtract R_termination + R_lead from R_meas to obtain R_winding(f).

## AC resistance ratio

Compute F_R(f) = R_AC(f) / R_DC and plot on a log-log scale vs f.
A typical Litz curve has:

- F_R = 1 below the strand corner frequency f_strand = ~ delta_strand^2
  * sigma * mu_0 / 4
- F_R rises slowly (~ 1 + (f/f_strand)^2) up to the bundle corner
  frequency f_bundle
- F_R rises faster (~ (f/f_strand)^4 or steeper) above f_bundle, due
  to internal proximity in the cross-section

## Common measurement pitfalls

| Pitfall                             | Symptom                                |
|-------------------------------------|----------------------------------------|
| Q > 100 in LCR mode                 | R noise floor dominates                |
| Sample temperature drift            | R_DC shifts ~ 0.4 %/K                  |
| Inductive lead loop                 | extra L, X *and* eddy R in lead loops  |
| Termination resistance              | overestimates R, especially at low f   |
| Magnetic environment (steel desk)   | extra core loss, raises R              |
| Coil flexing during scan            | L wobble couples to R via Q            |

The lab convention: report R_AC(f) with explicit uncertainty bounds
+/- delta_R from a *repeated* measurement series (typically 5 sweeps,
1-sigma standard deviation reported).

## Cross-check

For a published Litz spec (manufacturer data sheet typically gives
F_R at 100 kHz, 500 kHz, 1 MHz):
- Lab measured F_R should agree with vendor value within 10 %
- If discrepancy > 20 %, suspect (in order):
  1. Termination resistance not subtracted
  2. Sample flexing or different geometry than vendor's reference
  3. Vendor data sheet uses different ambient (room temp, no enclosure)
  4. Measurement frequency or technique mismatch

After cross-check, the measured F_R(f) is the ground truth for
validating any analytical / numerical model.

## Reference

- IEEE Std 389-1990, "IEEE Recommended Practice for Testing
  Electronics Transformers and Inductors".
- Keysight Application Note 1369-1, "Impedance Measurement Handbook,
  6th Edition", 2020.
- IEC 60404-6, "Magnetic materials -- Part 6: Methods of measurement
  of the magnetic properties of magnetically soft metallic and
  powder materials at frequencies in the range 20 Hz to 100 kHz by
  the use of ring specimens".
"""


def get_knowledge(topic: str = "overview") -> str:
    """Dispatch Litz / transmission line topics.

    Topics:
        overview              - Skin depth + AC loss summary (DEFAULT)
        litz_construction     - Hierarchical Litz structure / packing / twist
        skin_proximity_basics - Single-wire Bessel + Belevitch orthogonality
        dowell_formula        - 1D Dowell (1966) + Wojda M1 round adaptation
        ferreira_strands      - 2D Bessel M2/M3/M4 models for round strands
        umetani_multi_level   - Measurable-parameter multi-level Litz (2021)
        homogenization        - Igarashi complex-permeability + Ollendorff
        multiconductor_tl     - MTL theory: telegrapher, P/L, modes (Toki/Sato)
        fem_thin_wire         - MoM / FEM for rectangular TL impedance
        abc_pml               - Asymptotic / PML-like absorbing BCs
        tlm_method            - Transmission-line modelling (Johns / Sadiku)
        rosskopf_coupled      - Coupled FEM + PEEC for inductive components
        ac_loss_measurement   - Lab measurement protocol for AC R validation
        all                   - Everything (concatenated)
    """
    topic = topic.lower().strip()

    # Aliases / legacy names
    if topic in ("overview", "intro"):
        return OVERVIEW
    if topic in ("litz_construction", "litz", "litz_wire", "wire",
                 "construction", "structure"):
        return LITZ_CONSTRUCTION
    if topic in ("skin_proximity_basics", "skin", "proximity",
                 "skin_effect", "proximity_effect", "basics",
                 "belevitch", "kelvin"):
        return SKIN_PROXIMITY_BASICS
    if topic in ("dowell_formula", "dowell", "m1", "wojda", "foil"):
        return DOWELL_FORMULA
    if topic in ("ferreira_strands", "ferreira", "m2", "m3", "m4",
                 "bartoli", "tourkhani"):
        return FERREIRA_STRANDS
    if topic in ("umetani_multi_level", "umetani", "multi_level",
                 "multilevel", "sullivan"):
        return UMETANI_MULTI_LEVEL
    if topic in ("homogenization", "igarashi", "ollendorff",
                 "complex_permeability"):
        return HOMOGENIZATION
    if topic in ("multiconductor_tl", "mtl", "transmission_line",
                 "transmission", "telegrapher", "toki", "sato"):
        return MULTICONDUCTOR_TL
    if topic in ("fem_thin_wire", "fem", "mom", "moments",
                 "thin_wire", "matsuki", "matsushima"):
        return FEM_THIN_WIRE
    if topic in ("abc_pml", "abc", "pml", "khebir", "boundary"):
        return ABC_PML
    if topic in ("tlm_method", "tlm", "johns", "sadiku"):
        return TLM_METHOD
    if topic in ("rosskopf_coupled", "rosskopf", "coupled",
                 "fem_peec_coupled"):
        return ROSSKOPF_COUPLED
    if topic in ("ac_loss_measurement", "measurement", "lcr",
                 "impedance_analyzer", "validation"):
        return AC_LOSS_MEASUREMENT
    if topic == "all":
        return "\n\n".join([
            OVERVIEW,
            LITZ_CONSTRUCTION,
            SKIN_PROXIMITY_BASICS,
            DOWELL_FORMULA,
            FERREIRA_STRANDS,
            UMETANI_MULTI_LEVEL,
            HOMOGENIZATION,
            MULTICONDUCTOR_TL,
            FEM_THIN_WIRE,
            ABC_PML,
            TLM_METHOD,
            ROSSKOPF_COUPLED,
            AC_LOSS_MEASUREMENT,
        ])

    available = ", ".join(sorted(TOPICS.keys()))
    return f"Unknown topic '{topic}'. Available: {available}."
