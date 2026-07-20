"""Carstensen AC copper-loss knowledge for radia_mcp.peec.

Reference: Christian Carstensen, "Eddy Currents in Windings of
Switched Reluctance Machines", PhD dissertation, RWTH Aachen ISEA
(De Doncker primary advisor, Hameyer second reviewer), Nov 2007.
Chapter 3 (Fundamentals of Analytic Eddy Current Loss Calculation).

Lab library path:
  public-safe curated corpus
    Eddy Currents in Windings of Switched Reluctance Machines.pdf
  (187 pages, full PhD monograph)

Companion knowledge: `radia_mcp.motor.henrotte_lineage_knowledge.
carstensen_2007` for the broader research-arc context.

Coverage:
  motivation        — Why analytic loss matters for slot windings
  skin_round        — Skin effect in round conductor (Bessel formula)
  proximity_round   — Proximity effect from external H field
  dowell_layered    — Layered slot-winding formula (Dowell)
  combined_skin_prox — Total per-layer AC loss factor
  field_mapping     — How to map FE solution to per-layer H_slot
  python_recipe     — Wiring it into radia post-processing
  validation_cases  — Reference test cases for verification
  motor_application — End-to-end recipe for motor stator winding
"""

from __future__ import annotations


MOTIVATION = """\
## Why analytic AC copper loss matters for slot windings

A typical motor stator has 24-72 slots, each carrying 50-200 turns of
round copper wire (typical winding factor ~50-70%).  At any
appreciable inverter switching frequency:

- **Wire diameter** is typically 0.5-2 mm
- **Skin depth in Cu at 10 kHz**: delta = sqrt(2 / (omega mu sigma))
  = sqrt(2 / (62832 * 4pi*1e-7 * 5.8e7)) = 0.66 mm
- So d/delta is around 1-3 at the carrier frequency, AND each
  individual conductor sits in the high-tangential-H field of the
  slot environment

**Naive FE handling** would mesh every individual strand: with 100
turns × 24 slots × 8 elements per cross-section = 19,200 elements
across all stator coils ALONE.  Add through-thickness elements
for the 1-2 mm diameter (3-5 elements) and you're at >100,000
elements just for the copper.  Combined with the stator + rotor +
air-gap mesh, this kills the FE solve.

**Analytic AC loss formulae** solve this elegantly:

1. Mesh the slot as a SINGLE conductor region (no individual strands)
2. Use FE to compute the **tangential H field** in the slot
3. Apply analytical skin + proximity factors PER LAYER of conductors
4. Sum to get total R_AC per phase
5. Cross-check with measurement, no FE strand-meshing needed

This is the same trade-off as **Hollaus MSFEM** for laminated iron:
analytical homogenization replaces explicit fine-scale resolution.

### Where Carstensen sits in the lab's PEEC stack

| Method | Use case | radia_mcp module |
|--------|----------|------------------|
| PEEC with `nwinc/nhinc` strand subdivision | <= 10 strands, full mutual coupling | `peec.multi_filament` |
| Analytical skin+proximity (Carstensen-Dowell) | 100+ strands per slot, stranded coils | **this module** |
| SIBC (Bessel for round, Dowell for rect) | Solid conductors only, no strands | `peec.sibc` |
| Full 3D FE with strand resolution | Research only, infeasible for production | n/a |
"""

SKIN_ROUND = """\
## Skin effect in a round conductor (Carstensen Ch 3.1)

Reference: Carstensen 2007 PhD, Sec 3.1, pp. 32-39.  The classical
result for a round conductor of radius R carrying total current I at
angular frequency omega.

### Bessel-function solution

In a single round wire with axial current and no external field, the
internal current density satisfies a 1D radial diffusion equation
whose solution is a Bessel function of complex argument.

Define the dimensionless variable:

  xi = sqrt(omega * mu * sigma) * R = R / delta * sqrt(2)

(where delta = sqrt(2 / (omega mu sigma)) is the skin depth).

The Kelvin / Bessel-related per-unit-length AC resistance ratio is:

  R_AC / R_DC = (xi/2) * (ber(xi) * bei'(xi) - bei(xi) * ber'(xi))
                       / (ber'(xi)^2 + bei'(xi)^2)

where `ber, bei` are Kelvin functions (real and imaginary parts of
`J_0(j^(3/2) xi)`).

> **radia-ngsolve FE cross-check.** This exact ber/bei ratio is the regression
> reference for `radia_mcp.radia_ngsolve.solve.solve_planar_eddy` (complex A_z
> with a current/voltage-driven NumberSpace), matched to **0.07 %** at q=4
> (validation_test/radia_mcp/test_planar_eddy.py).  For multi-strand / litz bundles where the
> closed-form single-wire ratio breaks down, `solve_planar_eddy_multi` adds the
> **proximity** effect (series/parallel conductors sharing one A field) -- the
> +20 % Rac rise of two adjacent wires is in validation_test/radia_mcp/test_planar_proximity.py.
> Use the analytic ratio for quick estimates; the FE solvers for real geometries.

In `scipy` this evaluates to:
```python
from scipy.special import jv, yv
import numpy as np

def round_wire_R_AC_ratio(R, omega, sigma, mu_r=1.0):
    mu = 4e-7 * np.pi * mu_r
    delta = np.sqrt(2.0 / (omega * mu * sigma))
    xi = R / delta * np.sqrt(2.0)
    # Use scipy.special.iv for modified Bessel (CLAUDE.md NOTE: NOT jv!)
    # Actually Kelvin functions are jv with complex argument:
    # J_0( j^(3/2) xi ) = ber(xi) + j bei(xi)
    # Easier: use the asymptotic form for large xi:
    if xi < 0.5:
        # Low-freq limit: R_AC/R_DC = 1 + (1/12) (R/delta)^4
        return 1.0 + xi**4 / 48.0
    elif xi > 5:
        # High-freq limit: R_AC/R_DC = R/(2 delta) + 1/4 + 3/(16*delta/R)
        return R / (2 * delta) + 0.25
    else:
        # Exact via mpmath Kelvin functions:
        from mpmath import ber, bei, berprime, beiprime
        b = ber(xi); bi = bei(xi); bp = berprime(xi); bip = beiprime(xi)
        return float(xi/2 * (b * bip - bi * bp) / (bp**2 + bip**2))
```

### Asymptotic forms

| Regime | R/delta | R_AC / R_DC |
|--------|---------|-------------|
| Low frequency | < 0.5 | 1 + (R/delta)^4 / 48 |
| Transition | 0.5 - 5 | exact Kelvin (see above) |
| High frequency | > 5 | R / (2 delta) + 1/4 |

### Physical interpretation

- Low f: current uniform across cross-section, R = R_DC
- High f: current concentrated in surface layer of thickness delta,
  effective cross-section = 2 pi R delta, so R_AC = R_DC * R/(2 delta)
- Transition: smooth interpolation via Kelvin functions

For Cu at 50 Hz, R = 1 mm wire: delta = 9.3 mm >> R, so R_AC = R_DC.
For Cu at 10 kHz, R = 1 mm: delta = 0.66 mm < R, so R_AC > R_DC.
"""

PROXIMITY_ROUND = """\
## Proximity effect in a round conductor (Carstensen Ch 3.2)

Reference: Carstensen 2007 PhD, Sec 3.2.  When a round conductor is
placed in an EXTERNAL tangential H field H_ext (e.g. from adjacent
turns in a slot), eddy currents are induced even at ZERO net current
through the wire.  This is the "proximity effect".

### Proximity loss per wire

For a single round conductor of radius R in an external uniform H
field of amplitude H_ext, the time-averaged power dissipation is:

  P_prox = (R_DC * I_skin_unit^2) * G_prox(xi)

where `I_skin_unit = pi R^2 * (omega * mu_0 * H_ext)` is a reference
current and `G_prox(xi)` is the proximity factor:

  G_prox(xi) = -(2 xi / R_DC) * (ber_2(xi) ber'(xi) + bei_2(xi) bei'(xi))
                              / (ber(xi)^2 + bei(xi)^2)

(with ber_2, bei_2 the 2nd-order Kelvin functions).

In the limit xi << 1 (low f):
  G_prox(xi) ~ pi R^2 * (xi^4 / 192) * (sigma / R_DC)
             = (1/96) * (R/delta)^4 * R_DC * (sigma * pi R^2)^2

Equivalent simpler form: P_prox = (pi^2 / 96) * R^4 * sigma * (omega * H_ext)^2 / R_DC

### Combined with skin effect

For a wire carrying its own current I AND sitting in external H_ext,
the TOTAL loss is the SUM of skin and proximity contributions
(linear superposition is exact for the integral-equation problem):

  P_total = P_skin(I) + P_prox(H_ext)
          = R_AC(omega) * I^2 + G_prox(xi) * H_ext^2 / sigma_eq

### Why H_ext matters for slot windings

In a stator slot with N_layers of conductors stacked between slot
walls, each layer experiences a DIFFERENT H field:

- Bottom-of-slot layer (closest to air-gap): H_ext is small
  (slot opening field)
- Middle layers: H_ext increases as you sum contributions from
  layers below
- Top-of-slot layer (deepest in slot): H_ext is LARGEST (sum of
  all conductor currents in layers below)

The TOP layer experiences ~N_layers times the H_ext of the bottom
layer, so its proximity loss is N_layers^2 times larger.  This is
the **dominant** loss mechanism in tall slot windings.
"""

DOWELL_LAYERED = """\
## Dowell layered-slot formula (Carstensen Ch 3.2.2)

Reference: Carstensen 2007 PhD, Sec 3.2.2, building on Dowell 1966.

Dowell's classical formula for a slot containing m **layers** of
rectangular conductors of height h, with one-turn-per-layer (each
layer carries the full phase current I_phase):

  R_AC / R_DC = phi(xi)
              + (2/3) * (m^2 - 1) * psi(xi)

where xi = h / delta * sqrt(2) (one-direction skin depth in the
conductor) and:

  phi(xi) = xi * (sinh(2 xi) + sin(2 xi)) / (cosh(2 xi) - cos(2 xi))
  psi(xi) = xi * (sinh(xi)   - sin(xi))   / (cosh(xi)   + cos(xi))

### Interpretation

- phi(xi) = pure skin-effect factor (single layer carrying I_phase)
- (2/3) * (m^2 - 1) * psi(xi) = proximity factor that grows
  QUADRATICALLY with the number of layers m

For low frequency (xi -> 0):
  phi(xi) -> 1
  psi(xi) -> xi^4 / 3
  R_AC / R_DC -> 1 + (2/9) * (m^2 - 1) * xi^4

For high frequency (xi -> infinity):
  phi(xi) -> xi
  psi(xi) -> xi
  R_AC / R_DC -> xi + (2/3) * (m^2 - 1) * xi
              = (1 + (2/3)(m^2 - 1)) * xi
              = ((2 m^2 + 1) / 3) * xi  → very large

So for m = 10 layers at xi = 2, R_AC / R_DC ≈ (201/3) * 2 ≈ 134.
This is the "winding-loss disaster" that Carstensen's PhD addresses.

### Generalization to multi-turn-per-layer windings

If each of the m layers contains N_perlayer parallel strands of
height h_strand, replace h with h_strand in the xi definition.
The m^2 scaling becomes (m * N_perlayer)^2 for the total ampere-
turns per layer.

### Carstensen's specific contribution

Beyond Dowell, Carstensen extends to:
- **Trapezoidal slots** (typical motor geometry, not rectangular)
- **Non-uniform conductor distribution** (filling factor < 1)
- **Round wires** in addition to rectangular strands (with an
  effective h calculation)
- **Arbitrary current waveforms** (Fourier-decomposed PWM)
"""

COMBINED_SKIN_PROX = """\
## Combined skin + proximity per-layer formula

Reference: Carstensen Sec 3.2.3, extending Dowell.

For layer k (counted 1 to m from the slot bottom), the relevant H
field magnitudes inside that layer vary as the **stair-step**
function across the layer.  Carstensen's per-layer breakdown:

  R_AC,k / R_DC = phi(xi) + ((k-1)^2 + (k-1) * k * (2/3)) * (h/d_per_layer)^2 * psi(xi)

Summing over all m layers and normalizing by m:

  R_AC,total / R_DC,total
    = phi(xi) + (2/3) * (m^2 - 1) * psi(xi)

This recovers the Dowell formula.  The PER-LAYER breakdown is
useful for thermal analysis (hotspot prediction at the top layer)
and for optimization (consider air-cooling channels mid-slot).

### Python implementation

```python
import numpy as np

def carstensen_per_layer_R_AC_factor(layer_k, m_total, h, omega, sigma, mu_r=1.0):
    '''
    Carstensen per-layer AC resistance factor for slot winding.

    Args:
        layer_k : int, layer index (1-based) from slot bottom
        m_total : int, total number of layers in the slot
        h       : conductor height [m] (along the field direction)
        omega   : angular frequency [rad/s]
        sigma   : conductor conductivity [S/m]
        mu_r    : conductor relative permeability (1.0 for Cu/Al)

    Returns:
        R_AC,k / R_DC,k  (dimensionless)
    '''
    mu = 4e-7 * np.pi * mu_r
    delta = np.sqrt(2.0 / (omega * mu * sigma))
    xi = h / delta * np.sqrt(2.0)

    if xi < 1e-3:
        # Low-freq limit
        phi = 1.0
        psi = xi**4 / 3.0
    else:
        sh = np.sinh(2 * xi); ch = np.cosh(2 * xi)
        sn = np.sin(2 * xi);  cn = np.cos(2 * xi)
        sh1 = np.sinh(xi);    ch1 = np.cosh(xi)
        sn1 = np.sin(xi);     cn1 = np.cos(xi)
        phi = xi * (sh + sn) / (ch - cn)
        psi = xi * (sh1 - sn1) / (ch1 + cn1)

    # Carstensen per-layer:  phi + [(k-1)^2 + 2/3 (k-1)*k] * psi
    layer_factor = (layer_k - 1)**2 + (2.0 / 3.0) * (layer_k - 1) * layer_k
    return phi + layer_factor * psi


def carstensen_total_R_AC_ratio(m_layers, h, omega, sigma, mu_r=1.0):
    '''Dowell formula for whole-slot R_AC / R_DC ratio.'''
    return sum(carstensen_per_layer_R_AC_factor(k, m_layers, h, omega, sigma, mu_r)
               for k in range(1, m_layers + 1)) / m_layers
```

### Sanity checks

For m=1 (single layer): R_AC/R_DC = phi(xi) (pure skin).  At low f:
phi ~ 1, so R_AC ~ R_DC.  At high f: phi ~ xi = h/delta * sqrt(2),
so R_AC ~ R_DC * h / delta.

For m=10 layers at xi=1: phi(1) ~ 1.04, psi(1) ~ 0.027
  R_AC/R_DC = 1.04 + (2/3)(99)(0.027) = 1.04 + 1.78 = 2.82
"""

FIELD_MAPPING = """\
## Mapping FE solution to per-layer H_slot

The Carstensen formula needs `H` at each conductor layer.  In a 2D
motor FE, this comes from the standard A_z formulation:

```python
# After NGSolve solve:
from ngsolve import grad, Integrate

# B = (dA/dy, -dA/dx) in 2D
gA = grad(A_gridfunction)
B_in_slot = (gA[1], -gA[0])     # B vector inside slot
H_in_slot_mag = sqrt(B_in_slot[0]**2 + B_in_slot[1]**2) / mu_0

# Average H along each layer line:
for k in range(1, m_layers + 1):
    # Define layer-k as the horizontal line y = slot_bottom + (k-0.5)*d_layer
    layer_y = slot_bottom + (k - 0.5) * d_layer
    # Mesh selector: cells whose centroid is near y = layer_y
    # Use NGSolve's BBND (1D in 2D) line integral:
    H_avg_k = (Integrate(H_in_slot_mag, mesh,
                         definedon=mesh.BBoundaries(f"slot_layer_{k}"))
               / slot_width)
```

For this to work, the mesh must have **1D sidesets** at each layer
position (BBoundaries in NGSolve).  In the toy mesh
`build_test_motor_mesh.py`, the stator slot is a single region;
adding per-layer sidesets requires re-meshing with the slot
subdivided into m horizontal strips.

### Practical workaround for the first cut

When per-layer sidesets are not available, use **slot-averaged H**:

  H_avg = (1 / slot_volume) * integral_slot(|H|)

Then approximate per-layer as constant H_avg.  This gives a
slightly conservative estimate (over-counts proximity at the top
layer, under-counts at the bottom) but stays within 20% of the
proper per-layer result for typical motor slot geometries.

### Connection to PEEC nwinc/nhinc

The PEEC framework already has `nwinc` (number of width subdivisions)
and `nhinc` (number of height subdivisions) for splitting a single
filament into multi-strand bundles.  Setting nhinc = m_layers and
using the resulting per-strand currents in a PEEC self-consistent
solve gives the SAME answer as the Carstensen analytic formula in
the limit of fine subdivision — but at higher cost.

Carstensen is preferred when:
- m_layers >= 10 (PEEC strand subdivision becomes expensive)
- The slot geometry is captured by 2D FE anyway (re-use H field)
- Production-speed analysis is needed (post-process, no extra solve)
"""

PYTHON_RECIPE = """\
## Python recipe — drop-in post-processor

The complete analytical AC loss calculation, suitable as a
post-processor after `calc_motor_transient.py` or
`calc_fem_kelvin.py`:

```python
'''Carstensen AC copper loss for a motor stator winding.

Inputs:
    H_slot_per_layer  (numpy array, len=m_layers) -- |H| averaged
                      over each conductor layer in the slot, from FE
    I_phase           -- RMS phase current [A]
    R_DC              -- DC phase resistance [Ohm]
    h                 -- conductor height in field direction [m]
    sigma             -- conductor conductivity [S/m]
    freqs             -- list of frequencies to evaluate [Hz]

Output:
    R_AC_per_freq     -- AC resistance at each frequency [Ohm]
    P_AC_per_layer    -- power dissipation per layer [W]
'''

import numpy as np
from scipy.special import jv  # Kelvin functions via Bessel of complex arg

MU_0 = 4e-7 * np.pi

def kelvin_factors(xi):
    '''Return (phi, psi) for Carstensen's Dowell formula.'''
    if xi < 1e-3:
        return 1.0, xi**4 / 3.0
    sh2 = np.sinh(2 * xi); ch2 = np.cosh(2 * xi)
    sn2 = np.sin(2 * xi);  cn2 = np.cos(2 * xi)
    sh1 = np.sinh(xi);     ch1 = np.cosh(xi)
    sn1 = np.sin(xi);      cn1 = np.cos(xi)
    phi = xi * (sh2 + sn2) / (ch2 - cn2)
    psi = xi * (sh1 - sn1) / (ch1 + cn1)
    return phi, psi


def carstensen_R_AC(H_slot_per_layer, I_phase, R_DC, h, sigma, freqs):
    '''Whole-slot AC resistance vs frequency.'''
    R_AC_list = []
    for f in freqs:
        omega = 2 * np.pi * f
        delta = np.sqrt(2.0 / (omega * MU_0 * sigma))
        xi = h / delta * np.sqrt(2.0)
        phi, psi = kelvin_factors(xi)
        # Per-layer R_AC contribution
        m = len(H_slot_per_layer)
        R_AC_layers = []
        for k, H_k in enumerate(H_slot_per_layer, start=1):
            # H ratio: H_k / I_phase gives "geometric H factor"
            H_ratio = H_k / max(I_phase, 1e-30)
            # Carstensen per-layer formula in terms of H field:
            R_AC_layer = R_DC * (phi + ((k - 1)**2 + (2.0/3.0)*(k - 1)*k)
                                  * psi * H_ratio**2)
            R_AC_layers.append(R_AC_layer)
        R_AC_total = sum(R_AC_layers) / m
        R_AC_list.append(R_AC_total)
    return np.array(R_AC_list)


def carstensen_P_loss(I_phase_rms, H_slot_per_layer, R_DC, h, sigma, f):
    '''Power dissipation per layer [W].'''
    omega = 2 * np.pi * f
    delta = np.sqrt(2.0 / (omega * MU_0 * sigma))
    xi = h / delta * np.sqrt(2.0)
    phi, psi = kelvin_factors(xi)
    m = len(H_slot_per_layer)
    P_per_layer = []
    for k, H_k in enumerate(H_slot_per_layer, start=1):
        H_ratio = H_k / max(I_phase_rms, 1e-30)
        R_AC_k = R_DC * (phi + ((k-1)**2 + (2.0/3.0)*(k-1)*k)
                          * psi * H_ratio**2)
        P_per_layer.append(R_AC_k * I_phase_rms**2)
    return P_per_layer
```

### Where to plumb this

Two natural places:

1. **Post-processor** in `calc_motor_transient.py`: after the main
   FE solve, extract H per layer and append Carstensen AC loss
   to the JSON output.
2. **Standalone CLI** `calc_carstensen_loss.py`: takes a JSON
   from a prior FE run + winding parameters, outputs R_AC, P_AC.

Path 1 is more user-friendly but couples lamination + transient +
winding logic.  Path 2 keeps concerns separated.

For Stage-2 cleanliness, prefer Path 2.
"""

VALIDATION_CASES = """\
## Validation cases (Carstensen Ch 3.4)

The PhD provides four reference cases for verifying the analytical
formula against measurement or full FE:

### Case 1: Single round wire, no external field (pure skin)

- Geometry: R = 0.5 mm Cu wire, free in air
- Frequency: 1-100 kHz sweep
- Expected: R_AC/R_DC follows Bessel-function curve (skin_round
  topic), no proximity contribution.
- Cross-check: scipy's iv (modified Bessel) directly + measurement
  (network analyzer).

### Case 2: 2 parallel wires, opposing currents (proximity dominant)

- Geometry: two R = 0.5 mm wires, 2 mm apart, opposite-direction
  currents (so H field between them is strong)
- Verifies the H-field external proximity formula.
- Expected: 5-15× R_DC at 10 kHz.

### Case 3: Slot winding, m = 5 layers

- Rectangular slot, 5 layers of 8 strands each, single-phase AC
- Verifies Dowell formula for layered geometry
- Expected: ~2-4× R_DC at 5 kHz, scaling like (1 + (2/3)*(m^2 - 1)*psi)

### Case 4: Full SRM stator winding under PWM excitation

- Real geometry, 5-pole SRM with bundled stranded coils
- PWM at 20 kHz carrier
- Compare:
  - Analytical Carstensen (this module)
  - Full 3D FE with strand resolution (reference)
  - Measurement (calorimetric loss)
- Expected agreement: ±10% on total loss

### Implementation strategy for lab

Build a test fixture `tests/peec/test_carstensen_dowell.py` that
verifies Cases 1, 2, 3 against closed-form / Wolfram-Alpha
reference values.  Case 4 needs a real measurement which we don't
have in lab, but can be reproduced from Carstensen's PhD Table 3.x
(reference data published in the thesis).
"""

MOTOR_APPLICATION = """\
## End-to-end recipe for motor stator AC copper loss

Combining the pieces:

### Stage-1: Get the slot H field from FE

Run `calc_motor_transient.py` (Lange-Henrotte-Hameyer) or
`calc_fem_kelvin.py` for the motor.  The output A_z field gives
B = (dA/dy, -dA/dx) everywhere.

In the **stator slot region** specifically, compute |B|/mu_0 = |H|
averaged per layer.

### Stage-2: Compute per-layer R_AC + P_AC via Carstensen

```python
from radia.panels.calc_carstensen_loss import carstensen_R_AC, carstensen_P_loss

# Winding parameters (per slot, per phase)
m_layers = 8                 # 8 layers of conductors
h_wire = 1.5e-3              # 1.5 mm conductor height per layer
R_DC = 0.040                 # 40 mOhm DC phase resistance
sigma_Cu = 5.8e7
I_phase_RMS = 20             # 20 A RMS

# H per layer (from FE)
H_layers = [120, 180, 240, 290, 340, 380, 420, 450]   # A/m, by layer index 1-8

# AC loss per phase
P_layers = carstensen_P_loss(I_phase_RMS, H_layers, R_DC, h_wire, sigma_Cu, freq=10000)
P_total_per_slot = sum(P_layers)
P_total_motor = P_total_per_slot * n_slots * n_phases  # 24 slots × 3 phases
```

### Stage-3: Add to JSON output

The total motor loss breakdown becomes:

```json
{
    "iron_losses_W": 1023.0,     // from Hollaus EM (calc_motor_lamination.py)
    "PM_eddy_losses_W": 45.2,    // future, from solid-PM FE
    "copper_DC_losses_W": 1200,  // I_phase^2 * R_DC * 3
    "copper_AC_extra_losses_W":  // from Carstensen
       340.5,
    "copper_AC_per_layer_W": [12, 18, 28, 38, 50, 62, 72, 60],
    "total_motor_loss_W": 2608.7
}
```

### Heat-source coupling for thermal analysis

The PER-LAYER AC loss breakdown is **critical** for thermal
analysis — the top layer (closest to slot opening) is hottest and
limits motor current rating.  Pass `P_layers` directly as a
piecewise source term to a thermal FE (NGSolve heat) solve.

The lab has standard heat solver infrastructure in
`calc_heat.py` / `calc_heat_axisym.py`; coupling Carstensen
P_layers in is a 1-day integration project.

### Loss-vs-frequency table

For controller design, sweep frequencies (PWM harmonics):

```python
freqs = np.logspace(np.log10(50), np.log10(50000), 50)
R_AC_sweep = carstensen_R_AC(H_layers, I_phase_RMS, R_DC, h_wire, sigma_Cu, freqs)
```

This R_AC(f) table feeds directly into the controller-side
"frequency-aware loss model" that the Lange-HH transient framework
needs for accurate iron-loss + copper-loss prediction at PWM
harmonics.
"""

SECTIONS = {
    "motivation": MOTIVATION,
    "skin_round": SKIN_ROUND,
    "proximity_round": PROXIMITY_ROUND,
    "dowell_layered": DOWELL_LAYERED,
    "combined_skin_prox": COMBINED_SKIN_PROX,
    "field_mapping": FIELD_MAPPING,
    "python_recipe": PYTHON_RECIPE,
    "validation_cases": VALIDATION_CASES,
    "motor_application": MOTOR_APPLICATION,
}


def get_carstensen_knowledge(topic: str = "motivation") -> str:
    """Return Carstensen analytical AC copper-loss knowledge.

    Args:
        topic: section name; "all" returns everything.
    """
    t = topic.strip().lower()
    if t == "all":
        return "\n\n---\n\n".join(SECTIONS[k] for k in SECTIONS)
    if t not in SECTIONS:
        valid = ", ".join(sorted(SECTIONS))
        return f"Unknown topic {t!r}. Valid topics: {valid}\n"
    return SECTIONS[t]
