"""Electromagnetic Non-Destructive Testing (NDT) knowledge.

This module is the substantive knowledge layer for the
`radia_mcp.ndt` subpackage. It complements the searchable
`bibliography_index_knowledge.py` catalog by providing distilled,
topic-organized prose drawn from the configured literature library at
`public-safe curated corpus`
(approximately 145 PDFs across 12 subfolders).

The catalog covers:
  - 01_eddy_current_NDT  (14 PDFs) -- the core technique
  - 02_hysteresis_models (25 PDFs) -- vector hysteresis for steel NDE
  - 03_integral_equation_BEM (17 PDFs) -- BEM for eddy current crack
  - 04_FEM_formulations  (22 PDFs) -- A, A-phi, T-Omega for eddy current
  - 05_open_boundary_PML (24 PDFs) -- PML, ABC for unbounded NDT geometry
  - 06_multiply_connected (7 PDFs) -- Topology (cut surfaces) for h
  - 07_multi_physics_coupled (7 PDFs) -- Darwin, magneto-mechanical
  - 08_scattering_antenna (23 PDFs) -- MoM, FDTD for microwave NDE
  - 00_metal_detection (5 PDFs) -- landmine clearance, EMI spectroscopy
  - 10_machines_related (4 PDFs) -- inductor, MFL on rotating machinery
  - 11_textbook_references (1 PDF) -- Itoh `Computer Techniques for EM`

Wired into the dispatcher tool `ndt(topic=...)` in server.py via
`register_topics_tool`.
"""

# Single source of truth for the dispatcher's accepted topic enum.
# Wired into the `ndt_topics()` MCP tool via common.register_topics_tool.
TOPICS: dict[str, str] = {
    "overview": "Electromagnetic NDT landscape: ECT, MFL, microwave, "
                "EMI; physical regimes (DC -- GHz)",
    "eddy_current_basics":
        "Skin depth, impedance plane, lift-off, induction principle",
    "probe_types":
        "Pancake / bobbin / split-D / array / infinity-coil; absolute "
        "vs differential; reflection vs transmission",
    "defect_models":
        "Volumetric flaw, surface-breaking crack, sub-surface defect; "
        "back-side defect detection",
    "signal_processing":
        "Impedance plane analysis, multi-frequency mixing, calibration "
        "standards (EDM notch), POD/POS",
    "fem_formulation":
        "A-V vs A-A* vs T-Omega; edge elements; gauge fixing for eddy "
        "current FEM",
    "integral_equation":
        "BEM / VIE for cracks; line integral with loop electric + "
        "surface magnetic currents; reduced-order formulations",
    "multiply_connected":
        "Cut surfaces, tree-cotree, h formulation in regions with "
        "holes (tubes, plates with rivet holes)",
    "open_boundary":
        "PML, asymptotic BC, infinite elements, strategic dual image, "
        "FEM-BEM coupling for unbounded NDT",
    "hysteresis_in_ndt":
        "Vector Play / Stop / Jiles-Atherton; magneto-elastic coupling; "
        "MFL of plastically deformed steel",
    "mfl_magnetic_flux_leakage":
        "DC + Hall sensor array for ferromagnetic pipe / rail / plate; "
        "yoke magnetization; PIG inspection",
    "microwave_emi_landmine":
        "Metal detector (DC -- 100 kHz), GPR (100 MHz -- 3 GHz), EMI "
        "spectroscopy for UXO / landmine clearance",
    "ml_for_ndt":
        "Inverse problem (defect parameter estimation), POD, deep "
        "learning for crack characterization",
    "all":
        "Concatenated dump of all topics above",
}


# ---------------------------------------------------------------------------
# OVERVIEW
# ---------------------------------------------------------------------------

OVERVIEW = r"""
# Electromagnetic Non-Destructive Testing (NDT / NDE)

Non-destructive Testing (NDT, also called Non-Destructive Evaluation /
NDE) measures material properties or detects flaws without damaging
the part. The electromagnetic family is the largest in industrial use,
covering everything from gold detection in metal detectors to in-pipe
rail/pipeline inspection at speeds of 10 m/s.

This subpackage covers the electromagnetic-only methods. X-ray,
ultrasonic, and dye-penetrant NDT are out of scope.

## Method taxonomy

| Method | Frequency band | Penetration | Best for |
|--------|---------------|-------------|----------|
| **Metal detector** | DC -- 100 kHz | 0.5 -- 1 m air | Buried metal, landmines |
| **Eddy current testing (ECT)** | 1 kHz -- 10 MHz | one skin depth | Surface / near-surface flaws |
| **Magnetic flux leakage (MFL)** | DC magnetization + transit | full thickness | Pipes, rails, storage tanks |
| **Microwave / GHz NDE** | 1 -- 100 GHz | wavelength-scale | Composites, dielectrics |
| **GPR (ground penetrating radar)** | 100 MHz -- 3 GHz | 1 -- 30 m soil | Buried utilities, UXO |
| **Electromagnetic Induction Spectroscopy** | 100 Hz -- 100 kHz | 0.5 -- 5 m soil | UXO classification |
| **NMR MOUSE** | 10 -- 20 MHz | mm-scale | Polymer crosslink, moisture |

## Physical regimes

The unifying physics is Maxwell's equations in the quasi-static (MQS /
Darwin) regime for ECT/MFL/metal detection, and full-wave for
microwave/GPR. Skin effect dominates the eddy current methods:

  delta = sqrt(2 / (mu omega sigma))

For carbon steel (mu_r ~ 100, sigma ~ 5e6 S/m) at 10 kHz the skin
depth is ~0.7 mm. For aluminium (mu_r = 1, sigma = 3.5e7 S/m) at
10 kHz it is ~2.7 mm. This sets the depth of inspection -- a single
ECT scan can only see one skin depth into the wall, and that is why
multi-frequency mixing and pulsed-eddy-current are used to extend the
inspection depth.

## Applications by industry

| Industry | Method | Inspected |
|----------|--------|-----------|
| **Aviation** | ECT array | fastener holes, fatigue cracks under rivets |
| **Power generation** | bobbin ECT | heat-exchanger tubes (Inconel, SS, Cu-Ni) |
| **Rail** | MFL + ECT | wheels, axles, rail web cracks |
| **Pipeline** | MFL PIG, RFT | corrosion pitting, ID/OD distinction |
| **Petrochemical** | ECT, MFL | storage tank floor, vessel walls |
| **Auto / Medical** | impedance ECT | bearing race, surgical implants |
| **Defense / Humanitarian** | GPR + metal detector + EMI | mines, UXO |

## Lab connections (Sugahara Lab, Kindai University)

The configured literature corpus holds approximately 145 NDT papers, dominated by:

- Hokkaido / Igarashi school -- vector potential ECT (Igarashi 1998,
  thin magnetic conductive layer IBC).
- AEM Society (Japan Society for Applied Electromagnetics) -- vector
  hysteresis JSAEM rotational-saturation benchmarks.
- 10 Mitsubishi-NIFS heliotron papers (related stellarator MFL).
- Korean MFLT (Kang 2008) hall-sensor array for ferromagnetic crack.
- Bowler (Surrey) -- half-space Green's function for ECT coils.

## Cross-references

- `radia_mcp.ih`  -- induction heating shares the eddy current /
  workpiece physics (SIBC, ESIM)
- `radia_mcp.magnetic_materials`  -- BH curve, vector hysteresis
- `radia_mcp.fem.potential_formulations.t_omega`  -- ECT formulation
- `radia_mcp.bem` -- BEM for unbounded geometry (open NDT)
- `radia_mcp.radia_ngsolve` -- Kelvin transformation, PML, AND validated
  eddy-current FEM solvers usable directly as ECT forward models:
  `solve.solve_planar_eddy` (complex A_z, current/voltage-driven probe coil;
  validated vs round-wire Kelvin Rac to 0.07%), `solve.solve_planar_eddy_multi`
  (multi-conductor skin + PROXIMITY -- array probes / litz windings),
  `solve.solve_eddy_current_harmonic_APhi` (3-D A-Phi).  The FEMM Kelvin corpus
  (`literature_semantic_search(collection="femm_kelvin")`) bundles
  planar_coil_sigma crack-scan examples (50-5000 Hz, crack 0-2 mm) -- a
  ready-made ECT lift-off / crack-depth benchmark set.
"""


# ---------------------------------------------------------------------------
# EDDY CURRENT BASICS
# ---------------------------------------------------------------------------

EDDY_CURRENT_BASICS = r"""
# Eddy Current Testing (ECT) -- Fundamentals

Eddy Current Testing detects discontinuities by measuring the change
in coil impedance when the alternating magnetic field of a probe is
disturbed by a flaw in a conducting workpiece.

## Faraday + Ampere loop

1. Probe coil carries AC current I(t) at frequency f (1 kHz -- 10 MHz)
2. Coil creates time-varying H field above the workpiece surface
3. By Faraday's law dB/dt induces an EMF that drives circulating
   "eddy" currents J_e in the conductor
4. By Lenz's law, J_e produces a secondary B that opposes the primary
5. Secondary B couples back into the probe coil, changing its
   complex impedance Z = R + j omega L

A defect perturbs J_e (forces detour around / through it). This
changes the secondary B, which appears as a change in Z. ECT measures
this Delta Z.

## Skin depth

For a homogeneous half-space at angular frequency omega:

  delta = sqrt(2 / (mu_0 mu_r sigma omega))

| Material         | sigma (S/m) | mu_r | delta @ 10 kHz |
|------------------|-------------|------|----------------|
| Copper           | 5.8e7       | 1    | 0.66 mm        |
| Aluminium 7075   | 1.9e7       | 1    | 1.16 mm        |
| Stainless 316L   | 1.4e6       | 1    | 4.26 mm        |
| Inconel 600      | 1.0e6       | 1    | 5.03 mm        |
| Carbon steel A36 | 7.0e6       | 100  | 0.19 mm        |
| Ferritic SS 430  | 1.7e6       | 600  | 0.16 mm        |

Note that ferromagnetic steels have a 10x to 30x smaller delta than
their non-magnetic counterparts at the same conductivity, because of
the mu_r factor. This is why ferromagnetic carbon steel tubes are
almost always inspected by MFL, not ECT -- ECT cannot see deeper than
~0.2 mm.

## The impedance plane

The standard ECT calibration diagram plots normalized impedance:

  X = (omega L) / (omega L_0)      (normalized inductive reactance)
  R = (R - R_0) / (omega L_0)      (normalized resistance, after
                                    subtracting coil DC resistance)

Where L_0 is the air-only inductance of the bare probe. The locus
traced as a probe moves from air to workpiece, or as the workpiece
properties vary, forms a characteristic "S" curve.

Key parameters that move the operating point on the impedance plane:

- **Lift-off**: probe distance above surface; large effect, mostly
  on |Z|, moderate phase
- **Conductivity sigma**: moves along a "conductivity curve" (mostly
  vertical for high sigma, mostly horizontal for low sigma)
- **Permeability mu_r**: dominates for ferromagnetic materials,
  shifts the entire S-curve
- **Flaw**: changes (sigma, mu_r) locally; appears as a small
  perpendicular excursion off the conductivity curve

Mastering the impedance plane (Forster 1954, Libby 1971) lets the
operator separate lift-off (a nuisance) from real flaws by phase
discrimination -- they appear at different angles.

## Analytical solutions

For circular probe over half-space / plate / tube, analytical
solutions exist via Hankel-transform integrals:

- **Dodd & Deeds (1968)** -- closed-form integral for circular coil
  above a layered half-space; foundation of ECT modeling.
- **Bowler (Surrey, 1986)** -- "Eddy current calculations using
  half-space Green's functions" -- scalar+dyadic Green's function for
  arbitrary time-harmonic source over imperfect conductor.
- **Burke (1990s)** -- analytical solutions for crack as
  boundary-condition modification to Dodd-Deeds.

Even with these, numerical methods (FEM, BEM, VIE) are needed when:
- the geometry is non-axisymmetric (split-D, racetrack, butterfly)
- the workpiece is non-planar (tube, fastener hole, fillet)
- the crack is finite (not a sheet) or in a junction

## Lab reference

public-safe curated corpus

- "Analytical Solutions to Eddy-Current Probe-Coil Problems" (the
  closed-form reference manual)
- "Numerical Analysis of Eddy Current Testing by Integral Equation
  Method" (foundational IE for ECT)
- Bowler "Eddy current calculations using half-space Green's
  functions" (1986)
"""


# ---------------------------------------------------------------------------
# PROBE TYPES
# ---------------------------------------------------------------------------

PROBE_TYPES = r"""
# ECT Probe Types

Probe geometry is the single largest factor in ECT sensitivity --
each probe topology is matched to a specific defect type, scan
direction, and workpiece geometry. Choosing the wrong probe is the
most common error in field ECT.

## Major families

### Pancake (flat spiral) coil

Single planar spiral, axis perpendicular to surface.
- Use: surface scan of flat plates, wing skins
- Sensitivity: all in-plane orientations (omni)
- Liability: very lift-off sensitive
- Field shape: roughly axial below coil; eddy currents are concentric

### Bobbin (encircling) coil

Coaxial coil that fits inside a tube; axis is the tube axis.
- Use: heat-exchanger tube inside-diameter scan
- Sensitivity: circumferential cracks (eddy currents are
  circumferential, so a circumferential crack interrupts them)
- Blind to: axial cracks (parallel to current flow, no disturbance)
- Configurations: absolute (single coil), differential (two
  coils in series-opposition, drift-rejecting)

### Pencil / point probe

Small (1 -- 2 mm) pancake on a hand-held pen.
- Use: hand scan of weld toes, small fillets, fastener heads
- Sensitivity: localized; resolves cracks ~0.5 mm long
- Liability: slow (manual), demands operator skill

### Split-D (orthogonal) coil

Two D-shaped half-coils placed back-to-back.
- Use: directional sensitivity to cracks perpendicular to the split
- Mode: drive one half, measure the other (transmit/receive)
- Advantage: directional crack discrimination
- Form: rectangular D pair is common in aerospace fastener-hole ECT

### Plus-point (Zetec)

Two perpendicular D coils crossing at the center.
- Use: weld inspection in heat-exchangers
- Sensitivity: both circumferential and axial cracks
- Cancellation: nominal field at center is zero, so only flaws
  appear; reduces lift-off noise

### Array / multiplexed coils

Linear array of N coils (4 -- 64 typical) covering a wide aperture.
- Use: rapid scanning of large areas (tank floors, plates)
- Modes: multiplexed transmit (one at a time), or driver-receiver
  pairs (one drives, neighbours pick up)
- Speed: 100 mm/s scan with one pass instead of 64 passes

### Infinity-coil sensor (8-shape)

Two coils wound in a figure-8 (infinity) pattern.
- Lab reference: "infinity-coil eddy current sensor for back-side
  defect detection" (Japanese paper in the lab library)
- Use: surface-mount coil for back-side (far-side) flaw detection
- Principle: figure-8 winding cancels lift-off noise but is
  asymmetric to a back-side defect under one lobe

### Hall-effect probe (DC + Hall array)

Permanent magnet or DC coil magnetizes the part, Hall sensor reads
flux leakage. Used for MFL (see mfl_magnetic_flux_leakage topic).

## Absolute vs differential

- **Absolute coil**: one coil; reads true impedance value; sensitive
  to absolute conductivity / thickness / lift-off (good for sorting
  / sizing)
- **Differential pair**: two coils in series-opposition; reads
  difference -- gradient of defect signal as probe moves -- which
  rejects DC drift and slow lift-off changes (good for cracks but
  blind to long parallel flaws like wall thinning)

## Transmit-receive (TR) configuration

Two-coil topology where one drives and the other senses, sharing
the same flux path. Decouples drive impedance from receive
sensitivity. Used in:
- Cracker probe (aerospace fastener)
- Plus-point weld probe (heat-exchanger)
- Array probes (multiplexed driver pairs)

## Probe design tradeoffs

| Parameter | Effect |
|-----------|--------|
| Larger coil | Deeper penetration, lower spatial resolution |
| Smaller coil | Finer resolution, surface-only sensitivity |
| Higher freq | Smaller skin depth, more sensitivity to surface flaws |
| Lower freq | More penetration, less surface noise |
| More turns | Higher impedance, better SNR, slower frequency response |

## Lab reference

public-safe curated corpus
- "An Eddy Current Testing Platform System for Pipe Defect Inspection
  Based on an Optimized Eddy Current Technique Probe Design" --
  modern parametric coil optimization
- Japanese infinity-coil paper for back-side defect detection
"""


# ---------------------------------------------------------------------------
# DEFECT MODELS
# ---------------------------------------------------------------------------

DEFECT_MODELS = r"""
# Defect Models in Electromagnetic NDT

The mathematical model of a flaw determines what an ECT / MFL signal
looks like, what code the inverse problem solves, and what the
simulation forecast tells the operator. Three families dominate.

## Family 1: Sheet-like (zero-volume) crack

The crack is modeled as an infinitesimally thin sheet across which
the tangential E and normal J vanish (i.e. it is an electric-current
barrier). The classical analytical work is Burke (1985, 1988) for
surface-breaking, semi-elliptical cracks treated as additional BC
on top of Dodd-Deeds analytical solution.

- Pros: linear in the crack-opening displacement; fast forward
  solution
- Cons: assumes infinite thinness; cannot model the crack mouth's
  finite gap, oxide products, or partial closure
- Use: validation of array probes against EDM notch standards

## Family 2: Volumetric flaw (3D blob)

The crack is meshed as a volume with locally modified material
properties:
- For ECT: sigma_flaw < sigma_bulk (typically sigma_flaw = 0 for an
  open void, or sigma_flaw = sigma_oxide for tight cracks)
- For MFL: mu_r_flaw = 1 (the missing iron leaves an air-permeability
  region inside an iron host)

- Pros: physically realistic for corrosion pitting, wall thinning,
  metallurgical inclusions
- Cons: requires meshing the inside of the defect, which is small
  (~50 micrometres for a tight crack); adaptive mesh refinement
  needed

## Family 3: Boundary-condition-modified surface

The defect is the boundary itself -- e.g. a pit on the OD of a tube
is modeled by deforming the OD surface inward. Combined with SIBC
(surface impedance BC), no inside-the-wall mesh is needed.

- Pros: keeps the mesh in air only; very fast for parametric studies
- Cons: cannot capture sub-surface defects; loses information once
  the surface is meshed coarser than the defect

## Surface-breaking vs sub-surface

The detectable depth depends on the skin depth delta:

| Depth | ECT signature |
|-------|---------------|
| Surface-breaking | Strong, sharp; large phase response |
| < 1 delta | Diminished; phase ~ 45 deg |
| 1 -- 2 delta | Detectable with low-freq excitation |
| > 2 delta | Below ECT sensitivity; use MFL or RFT |

This is why multi-frequency ECT mixes a high-frequency channel
(surface noise rejection) with a low-frequency channel (deeper
penetration).

## Back-side defect detection

A back-side defect (defect on the far side of a single-skin
workpiece relative to the probe) is the hardest ECT problem:
- The signal must traverse the wall twice (probe -> defect -> probe)
- It is attenuated by exp(-2 d / delta) where d is wall thickness
- For carbon steel at 10 kHz, d = 1 mm gives ~3000x attenuation

Solutions:
1. **Infinity-coil sensor** (lab reference) -- figure-8 winding
   cancels lift-off but preserves back-side defect signal
2. **Remote-field eddy current (RFT)** -- transmit and receive coils
   axially separated by 2 -- 3 tube diameters; receive coil sees
   only the through-wall field. Equal sensitivity to ID and OD
   defects.
3. **Pulsed eddy current (PEC)** -- step-function drive,
   time-domain reception; later time samples = deeper material.

## Sizing -- the inverse problem

Forward: given crack (length, depth, opening), compute Delta Z
Inverse: given Delta Z signal, recover crack (length, depth, opening)

Two approaches:
- **Calibration curves**: empirical, from EDM-notch standards of
  known dimensions. Industry standard for safety-critical use.
- **Model-based inversion**: solve a constrained nonlinear least
  squares with the forward model in the loop. Used for unique
  geometries not covered by standards.

## Lab reference

public-safe curated corpus
- "Eddy-Current Analysis of a Conductor With a Crack by Line
  Integral Equations With Loop Electric and Magnetic Currents as
  Unknowns" -- BEM for cracks (loop-current basis)
- "Numerical Analysis of Eddy Current Testing by Integral Equation
  Method" -- VIE forward model for ECT
- Korean MFL hall-sensor paper -- crack evaluation on ferromagnetic
  material via Hall-array MFL (Kang 2008)
"""


# ---------------------------------------------------------------------------
# SIGNAL PROCESSING
# ---------------------------------------------------------------------------

SIGNAL_PROCESSING = r"""
# Signal Processing in Eddy Current Testing

Raw ECT signal is the complex impedance Z(t) = R(t) + j omega L(t)
as a function of probe position. The signal-processing chain converts
this into a yes/no flaw decision and (optionally) a sizing estimate.

## Stage 1: Calibration

Before scanning, the probe is run over a reference standard with
known artificial flaws:
- **EDM notches**: rectangular cuts of known length, depth, width
  (typically 0.5 mm wide, 10 mm long, 10%/20%/40% through-wall depth)
- **ASME reference holes**: through-wall drilled holes
- **Tube standards**: drilled flat-bottom holes, OD/ID grooves

Calibration normalizes the impedance and rotates the impedance plane
so that lift-off lies along the horizontal axis and crack signals
lie along the vertical (phase-rotation to 0/90 degrees).

## Stage 2: Filtering

- **High-pass**: removes slow lift-off drift, geometric variations,
  ambient temperature drift
- **Low-pass**: removes electrical noise above the defect bandwidth
- **Median filter**: removes salt-and-pepper interconnect noise

## Stage 3: Phase rotation and gain

After filtering, the trace is plotted on the impedance plane (X-Y
oscilloscope mode) or as horizontal/vertical channels vs scan
position. Phase rotation aligns the lift-off direction with the X
axis so that:
- X channel: lift-off, geometry, edge effects (rejected)
- Y channel: crack signal (accepted)

This is the foundation of phase discrimination.

## Stage 4: Multi-frequency mixing

Common pattern: scan at two frequencies simultaneously (e.g. 100 kHz
+ 400 kHz), then linearly combine them to suppress a known nuisance
signature.

  S_mix = w1 * S_f1 + w2 * S_f2

Where w1, w2 are chosen so that the nuisance signal (e.g. tube
support plate, tube-sheet expansion, dents) cancels and the defect
signal survives. Standard in nuclear steam-generator tube
inspection.

## Stage 5: Detection threshold

For each segment of the scan, declare "flaw" if |Z_filtered| exceeds
an amplitude threshold AND falls within a phase window matching the
calibration response.

Setting the threshold is the tradeoff between:
- **POD (Probability of Detection)** -- want high
- **POFA (Probability of False Alarm)** -- want low

The receiver operating characteristic (ROC) curve plots POD vs POFA;
operating point chosen by industry code (e.g. ASME XI for nuclear).

## Stage 6: Sizing

If a flaw is detected, estimate its dimensions via:
- **Voltage normalization**: peak Y-channel voltage proportional to
  defect cross-section (rough)
- **Phase-amplitude curves**: depth from phase angle, length from
  width of the signal lobe (more refined, from calibration tables)
- **Model-based inversion**: solve the forward problem in a loop to
  match the observed signature (research / unique inspections)

## Probability of Detection (POD)

POD is typically reported as a function of defect size
(POD vs depth, POD vs length). Industry codes require:
- 90% POD at 90% confidence for the smallest detectable defect
- Reporting threshold at, e.g., 40% through-wall for safety-critical
  tubes

POD studies use the MIL-HDBK-1823A (Mil Handbook 1823A) methodology:
- A set of N known specimens with carefully measured flaw sizes
- Multiple operator runs (round-robin)
- Logistic regression to fit POD(size; a, b) curve

## Modern techniques

- **Wavelet denoising**: improves SNR for tube support plate signals
- **PCA / ICA**: separate flaw signal from coherent background
- **Deep learning**: CNN classifier on impedance plane "images"
  (lab reference: 06_非破壊検査 includes a Japanese deep-learning
  textbook for the lab)

## Lab reference

- public-safe curated corpus -- platform paper
  with optimized probe + signal processing
- public-safe curated corpus -- general
  deep-learning textbook cross-referenced for ML-NDT
"""


# ---------------------------------------------------------------------------
# FEM FORMULATION
# ---------------------------------------------------------------------------

FEM_FORMULATION = r"""
# FEM Formulations for Eddy Current NDT

Eddy current problems in NDT are 3D and time-harmonic (most ECT) or
transient (pulsed ECT, MFL transit). The choice of unknown field
(potential formulation) sets the algorithmic structure of the FEM
code.

## Quick comparison

| Formulation | Unknowns | Air region | Eddy region | Constraints |
|-------------|----------|------------|-------------|-------------|
| A-V | A everywhere + V in cond | OK | OK | Coulomb gauge |
| A* | A everywhere | OK | OK | Coulomb gauge (penalty) |
| T-Omega | T in cond + Omega everywhere | scalar | vector | Multiply-connected cut |
| h | H field directly | hard | natural | Multiply-connected cut |
| h-phi | H in cond + phi in air | natural | natural | Multiply-connected cut |

## A-V formulation

The classical workhorse. In the eddy-current region:

  curl(nu curl A) + sigma (dA/dt + grad V) = J_source     (Ampere)
                          div(sigma (dA/dt + grad V)) = 0  (charge conservation)

In air:

  curl(nu_0 curl A) = J_source

The Coulomb gauge `div A = 0` is enforced explicitly via the
penalty method (Bossavit 1988) or via the choice of A in a
divergence-free space (edge elements; see Three Dimensional Eddy
Current Calculation Using Edge Elements paper).

- Pros: simplest physical interpretation; mature solver chain
- Cons: dense ill-conditioned matrix; needs many DOFs in air
- Lab paper: "Numerical Investigation of Penalty Factor for Coulomb
  Gauge of Magnetic Vector Potential for 3D Fields" (Lin, Heyun, Tu,
  Sun, IEEE Trans. Magn. 1999)

## A-A* (Total / reduced)

The total vector potential A is split into source A_s (from coils,
analytical) + reduced A_r (from induced currents in workpiece, FEM).

  A = A_s + A_r
  A_s = (mu_0 / 4 pi) integral [J_source / r] dV (Biot-Savart)

Only A_r needs FEM. Avoids meshing the coil interior; coil is a
source-point list.

- Lab paper: "Total and Reduced Magnetic Vector Potentials and
  Electrical Scalar Potential for Eddy Current Calculation"

## T-Omega (current vector + scalar potential)

Replaces A with the current vector potential T (J = curl T) in the
conductor and the magnetic scalar potential Omega in air
(H = -grad Omega + T).

- Pros: only the conductor needs vector DOFs; air has scalar DOFs
  (1/3 the unknowns of A-V); naturally handles unbounded air
  (paired with BEM)
- Cons: needs a cut surface for multiply-connected conductors
  (tubes, plates with holes)
- Lab paper: "Three-dimensional nonlinear electromagnetic field
  computations, using scalar potentials"

## h formulation

H is the unknown. Natural for SC magnets, but for eddy current it
needs:
- HCurl elements (Nedelec edge elements)
- Multiply-connected cut surfaces if conductor has holes
- Lab paper: "Formulation of the Eddy Current Problem in Multiply
  Connected Regions in Terms of h" (key 6-page paper)

## Edge elements (Nedelec / Whitney)

The breakthrough that made vector formulations practical (Bossavit
1988, Nedelec 1980). DOFs are the line integrals of A or H along
edges, automatically:
- Tangentially continuous across faces (right for vector fields)
- Normally discontinuous (allows jumps in D or B)
- Naturally divergence-free in the lowest-order space (kills
  spurious modes that plagued nodal vector elements)

Lab papers:
- "Three Dimensional Eddy Current Calculation Using Edge Elements
  for Magnetic Vector Potential"
- "Hierarchal vector basis functions of arbitrary order for
  triangular and tetrahedral finite elements"
- "High Order Finite Element Methods for Electromagnetic Field
  Computation"

## Tree-cotree

For h or T-Omega with multiply connected geometry, the tree-cotree
decomposition (Albanese-Rubinacci 1988) provides the discrete cut
surface algorithmically:
- Build a spanning tree of the mesh edges
- Cotree edges are independent DOFs (loop currents)
- Tree edges are dependent (loop sums must satisfy Kirchhoff)
- Lab paper: "Tree-cotree decompositions for first-order complete
  tangential vector finite elements"
- Lab paper: "Application of Tree-Cotree Splitting to the Time-Domain
  Finite-Element Analysis of Electromagnetic Problems"

## Gauge fixing summary

| Gauge | Method | Reference |
|-------|--------|-----------|
| Coulomb (A-V) | Penalty +int div A * v dV | Bossavit, Lin-Tu |
| Coulomb (A) | Two-MVP regularization | Plagne-Albanese (lab paper) |
| None (h) | Tree-cotree cut | Albanese-Rubinacci |
| None (T-Omega) | Cut surface | Brown 1992 |

## Hybrid FE-BE for unbounded

Pair FEM (in conductor + immediate neighbourhood) with BEM (on the
outer boundary) so that:
- FEM handles nonlinear materials, complex shapes
- BEM handles the unbounded air region with one matrix
- Coupling via consistent boundary integrals

Lab paper: "Novel Boundary Element Formulation in Hybrid FE-BE
Methods for Electromagnetic Field Computations"
"""


# ---------------------------------------------------------------------------
# INTEGRAL EQUATION
# ---------------------------------------------------------------------------

INTEGRAL_EQUATION = r"""
# Integral Equation Methods for Eddy Current NDT

For NDT problems with unbounded air and small flaws, integral
equation (IE) methods are competitive with FEM:

- Only the conductor surface OR volume is meshed (not the air)
- Outgoing-wave / decay BC is automatically built in
- Multiply-connected geometry needs explicit cuts (same as h-FEM)

## Three classes

### 1. Surface IE (BEM proper)

DOFs on the conductor surface. Two unknowns: electric current J_s and
magnetic current M_s on the surface.

EFIE (Electric Field IE):
  n x E_inc = n x integral [j omega mu G(r,r') J_s + ...] dS'

MFIE (Magnetic Field IE):
  n x H_inc = J_s/2 + n x PV integral [grad G x J_s] dS'

Lab papers:
- "Three-dimensional eddy current analysis by the boundary element
  method using vector potential"
- "Eddy Current Analysis by BEM Formulated by Loop Electric and
  Surface Magnetic Currents Utilizing Edge Boundary Conditions"
- "Eddy Current Analysis by Integral Equation Method Utilizing Loop
  Electric and Surface Magnetic Currents as Unknowns" -- lab style
  loop+surface basis (compatible with HDiv-VIM philosophy)
- "Vector Potential Boundary Element Method for Three Dimensional
  Magnetostatic"

### 2. Volume IE (VIE)

DOFs are the equivalent currents (J for conductors, M for magnetic
materials) on the volume mesh.
- Lab paper: "A new vector element in the volume integral equation
  method for nonlinear magnetostatics"
- Lab paper: "A Locally Corrected Nystrom Formulation for the
  Magnetostatic Volume Integral Equation"
- Lab paper: "Integral Analysis of a Magnetic Field for an Arbitrary
  Geometry Coil Width Rectangular Cross Section"

VIE matches the Radia HDiv-VIM / face-charge volume routes
charge) philosophy -- see `radia_mcp.radia_ngsolve`.

### 3. Line IE (LIE)

For cracks modeled as 1D segments embedded in a half-space.
- Lab paper: "Eddy-Current Analysis of a Conductor With a Crack by
  Line Integral Equations With Loop Electric and Magnetic Currents
  as Unknowns"

Particularly efficient for fastener-hole ECT where the only
3D-relevant feature is the bolt-hole crack itself; the rest is
analytical (Dodd-Deeds half-space).

## Low-frequency breakdown

The classical EFIE has a notorious singularity at omega -> 0
(electrostatic limit): the matrix becomes ill-conditioned because the
charge contribution dominates over the current contribution by a
factor 1 / omega.

Modern fixes:
- **Loop-Tree / Loop-Star basis**: split RWG into solenoidal (loops)
  + irrotational (tree edges), preconditioning each block separately
- **Quasi-Helmholtz projectors**: orthogonal projection onto loop
  and star subspaces; lab paper:
  "Handling the Low-frequency Breakdown of the PMCHWT Integral
  Equation with the quasi-Helmholtz Projectors"
- **Taylor-Orthogonal basis**: modifies the basis to remove the
  singularity; lab paper:
  "Stable Discretization of the Electric-Magnetic Field Integral
  Equation With the Taylor-Orthogonal Basis Functions"

For eddy-current NDT (kHz -- MHz range) the low-freq breakdown is
NOT severe enough to require these fixes -- straightforward EFIE
works -- but they matter for DC magnetostatic IE.

## Multiply connected breakdown

The static MFIE has a null space for non-simply-connected surfaces
(tubes, plates with holes). Lab papers:
- "Nullspace of the Static MFIE Operator and its Effect on the
  Numerical Solution of BIEs"
- "Low Frequency Scaling of the Mixed MFIE for Scatterers with a
  Non-Simply Connected Surface"

Solution: add a cut surface inside the hole (Kotiuga 1991 -- see
"Clebsch_potentials_and_the_visualization_of_three-dimensional_solenoidal_vector_fields"
in the bibliography catalog).

## Acceleration

Dense BEM/VIE matrices are O(N^2). Acceleration:
- **FMM** (Fast Multipole) -- now retired in Radia (see CLAUDE.md);
  doesn't help compact NDT geometries
- **ACA / H-matrix** -- preferred (used in Radia via HACApK)
- **Locally Corrected Nystrom** -- high-order quadrature with
  singular-corner correction (lab paper above)

## Lab reference

public-safe curated corpus (17 PDFs)
public-safe curated corpus (BEM-flavoured papers)
- Dense IE solver: "Solution of Dense Systems of Linear Equations
  Arising from Integral-Equation Formulations" (bibliography catalog)
"""


# ---------------------------------------------------------------------------
# MULTIPLY CONNECTED
# ---------------------------------------------------------------------------

MULTIPLY_CONNECTED = r"""
# Multiply Connected Conductors in NDT

A multiply connected conductor is one in which a closed curve in the
conductor cannot be shrunk to a point without leaving the conductor.
Real NDT examples:

- **Tube**: doughnut topology -- a loop going around the tube hole
  is not contractible
- **Plate with rivet holes**: every hole gives a non-contractible
  loop
- **Punched coil bracket**: similar
- **Multi-aperture EM shield**: each window is a non-contractible
  loop

In all these the scalar magnetic potential Omega is multi-valued
(jumps by the enclosed current upon traversing the non-contractible
loop). To make Omega single-valued for FEM, the conductor must be
**cut** with a surface that closes the topology.

## Why this matters for NDT

ECT of multi-hole plates and bobbin-probe ECT of tubes both fail
silently in h, Omega, and T-Omega formulations if cuts are not
introduced. The simulation appears to converge but the integral of H
around the hole gives zero (instead of the enclosed current),
producing 50% magnitude errors that look plausible.

Detection: a hand-check that the integral of H along a
non-contractible loop equals the encircled net current. If 0, the cut
is missing.

## Cut surface

A cut is a 2-manifold S inside the conductor whose boundary lies on
the conductor surface, such that:
  conductor \ S is simply connected.

After cutting, Omega is well-defined; A double-valued condition is
imposed on the cut faces: Omega^+ - Omega^- = I_encl (the encircled
current).

## Algorithmic construction (tree-cotree)

Tree-cotree algorithm (Albanese-Rubinacci 1988; Bossavit-Mayergoyz
1992):
1. Build a spanning tree of the dual graph (edges = faces, nodes =
   elements)
2. The cotree (faces not in the tree) corresponds to independent
   loops
3. For each genus-g hole, g cotree faces form a cut surface

Lab papers:
- "Tree-cotree decompositions for first-order complete tangential
  vector finite elements" -- the algorithmic recipe
- "Application of Tree-Cotree Splitting to the Time-Domain
  Finite-Element Analysis of Electromagnetic Problems"

## Topology terminology

- **Genus g**: number of holes (tube = genus 1, two-hole plate =
  genus 2)
- **First Betti number b1**: number of independent non-contractible
  loops (= genus for orientable 2-manifolds)
- **Cut surface**: 2-chain whose boundary lies on the conductor
  boundary and whose removal makes the conductor simply connected

Number of cuts needed = b1 = number of independent loops.

## h formulation cuts

Lab paper: "A Novel h Approach for Solving Eddy-Current Problems in
Multiply Connected Regions" -- the dedicated reference; describes
explicit cut construction for h-field FEM.

Lab paper: "An approach to solution of three-dimensional voltage
driven and multiply connected eddy current problems" -- handles the
external-circuit constraint (voltage source on a loop).

## T-Omega cuts

Lab paper: "Fast Solution of 3-D Eddy-Current Problems in Multiply
Connected Domains by a, v-phi and t-phi Formulations With
Multigrid" -- comparison of A-V vs T-phi (= T-Omega) with cut and
multigrid acceleration.

## Vector potential T cuts

Lab paper: "Treatment of multiply connected regions in two-component
electric vector potentials formulations" -- T with two scalar gauges
+ cuts.

## Circuit-field coupling

In ECT of a coil array, each coil is a non-contractible loop driven
by an external voltage source. The IE formulation needs:
- A loop electric basis (Hertz loops -- divergence-free volume
  basis)
- A separate flux-linkage variable per loop (Faraday voltage)
- Lab paper: "Circuits_fields coupling and multiply connected
  domains in integral formulations"

## Lab reference

public-safe curated corpus (7 PDFs, all explicitly
on this topic)
"""


# ---------------------------------------------------------------------------
# OPEN BOUNDARY
# ---------------------------------------------------------------------------

OPEN_BOUNDARY = r"""
# Open Boundary Treatment for NDT FEM

ECT and MFL almost always have unbounded air -- the probe and
workpiece are in free space (or in an inspection chamber far larger
than the probe). FEM intrinsically discretizes a bounded domain;
truncation introduces reflection artefacts. The four families of
solution are:

| Method | Use | Cost | Accuracy |
|--------|-----|------|----------|
| **Box ABC** | Quick test | Cheap | Poor near box |
| **PML** | Wave problems | Cheap | Excellent for wave |
| **Infinite element** | Static / quasi-static | Cheap | Moderate |
| **Asymptotic BC** | Static / quasi-static | Cheap | Good (high order) |
| **Kelvin transform** | Static | Cheap | Excellent |
| **FEM-BEM** | All | Expensive | Excellent |

## PML for ECT / MFL

The Perfectly Matched Layer (PML; Berenger 1994) was developed for
wave problems where reflected energy is the issue. For ECT/MFL in
the quasi-static regime:

- **PML for FDTD** (Berenger 1994) -- the original; lab paper:
  "Perfectly Matched Layer Termination for Finite-Element Meshes"
- **PML for FEM frequency-domain** -- lab paper:
  "An Application of the Perfectly Matched Layer (PML) Concept to
  the Finite Element Method Frequency Domain Analysis of Scattering
  Problems"
- **Anisotropic PML** -- lab paper:
  "An Anisotropic Perfectly Matched Layer-Absorbing Medium for the
  Truncation of FDTD Lattices"
- **Static-field PML** -- the topic is non-trivial because PML was
  designed for outgoing waves and "static" has no wave; lab papers:
  - "Perfectly Matched Layer for Static and Quasistatic Fields A
    False Method" (warning paper -- the topic is subtle)
  - "Perfectly Matched Layers in Static Fields"
  - "Optimization of Perfectly Matched Layer for Laplace's Equation"
- **Eddy-current PML** -- lab paper:
  "Perfectly Matched Anisotropic Layer for the Numerical Analysis
  of Unbounded Eddy-Current Problems"

## Asymptotic Boundary Condition (ABC)

Khebir-Ramahi-Mittra (1990) for FEM, applied per-order:
- 1st-order: Sommerfeld-like
- 2nd-order: corrects for surface waves
- Higher: residual decays as 1/r^(2n+1)

Lab papers:
- "Asymptotic Boundary Conditions for Finite Element Analysis of 2D
  and 3D Electrical Field Problems"
- "Asymptotic Boundary Conditions for the Finite Element Modeling of
  Axisymmetric Electrical Field Problems"
- "Improvised Absorbing Boundary Conditions for Three-dimensional
  Electromagnetic Finite Elements"
- "Improvised Asymptotic Boundary Conditions for Electrostatic
  Finite Elements"

## Strategic dual image method (Saito-Hayano)

Lab paper: "The strategic dual image method for the open boundary
electromagnetic field problem" by Y. Saito, K. Takahashi, S. Hayano
(Hosei University). A clever image-based truncation that drops the
inverse-distance error by placing virtual image sources at a chosen
distance.

## Kelvin transformation

Map exterior of a sphere onto the interior of a "kelvin" sphere via
inversion through a reference sphere:
  r -> R^2 / r

The infinite exterior becomes a finite ball; the field decays to zero
at infinity is mapped to zero at the kelvin center. Standard in
Radia (`radia_mcp.radia_ngsolve.kelvin_*`) and for accelerator magnets.

Lab paper:
- "An open boundary Technique for axisymmetric and three dimensional
  magnetic and electric field problems"

## Asymptotic + infinite element review

Lab paper: "A Review of Finite Element Open Boundary Techniques for
Static and Quasi-Static Electromagnetic Field Problems" -- the
authoritative survey, covers all methods listed here.

Lab paper: "A case study on open boundary techniques for
electromagnetic field problems with translational symmetry" --
practical comparison.

## FEM-BEM coupling (DBCI etc)

Pair FEM in the interior with BEM on the truncating surface; the BEM
provides the exact decay condition. Lab papers:
- "Comparing FEM-BEM and FEM-DBCI for open-boundary electrostatic
  field problems"
- "Improving the Accuracy of the Integral Equation in the Hybrid
  FEM-DBCI Method for Open Boundary Electrostatic Problems"

DBCI = Difference-Based Coupled Integral; a variant that reduces the
BEM cost by writing only the field difference across the FEM-BEM
boundary.

## Trefftz methods

Lab paper: "Trefftz Approximations A New Framework for Nonreflecting
Boundary Conditions" -- frontier method using exact-solution basis
functions in the truncated region.

## Lab reference

public-safe curated corpus (24 PDFs)
"""


# ---------------------------------------------------------------------------
# HYSTERESIS IN NDT
# ---------------------------------------------------------------------------

HYSTERESIS_IN_NDT = r"""
# Hysteresis in NDT

Ferromagnetic NDT (carbon steel pipes, rails, ship hulls) is a
nonlinear problem because:
- B-H curve is hysteretic (~1 -- 2 T saturation, ~50 -- 500 A/m
  coercivity for low-carbon steel)
- The magnetization history matters: virgin curve vs major loop vs
  minor loop give different mu_r
- Plastic deformation, residual stress, and decarburization change
  the BH curve locally -- this is the basis of *magnetic Barkhausen
  noise* and *MFL stress sensing*

Linear MFL/ECT simulations underestimate signal because they ignore
the magnetization-saturation-yield path of the iron near the flaw.

## Scalar vs vector hysteresis

- **Scalar** (Preisach, Jiles-Atherton) -- single direction; OK for
  toroid cores, narrow-gauge wire, simple sheet
- **Vector** -- rotational excitation produces a rotating M that
  lags rotating H; NDT needs vector hysteresis because a moving
  probe rotates the local H direction

Sugahara Lab specialty: vector hysteresis. The 25 PDFs in
02_hysteresis_models/ are dominated by:
- Matsuo / Shimoji vector Stop and Play models (JSAEM benchmark)
- Egger energy-based vector hysteresis (variational, convex)
- Hantila polarization (also in Radia, see `radia_mcp.electromagnet`)

## Vector Play / Stop models

Both are isotropic-vector extensions of the scalar Play/Stop
operators. They satisfy:
- **Rotational saturation**: at high H_max the locus of M traces a
  circle (proper saturation)
- **Mirror symmetry**: M(-H) = -M(H)

Lab papers:
- "Isotropic vector hysteresis represented by superposition of stop
  hysteron models"
- "Correction to Isotropic vector hysteresis represented by
  superposition of stop hysteron models" -- correction sheet
- "Rotational Saturation Properties of Isotropic Vector Hysteresis
  Models Using Vectorized Stop and Play Hysterons"
- "Comparison of Rotational Hysteretic Properties of Isotropic
  Vector Stop Models"

## Anisotropic vector hysteresis

Real steel sheet is rolled and is therefore (i) anisotropic (rolling
direction differs from transverse), (ii) graded by depth (skin
deformation from rolling). Anisotropy is added via:
- A weighted-orientation distribution function over isotropic
  models -- lab paper: "Anisotropic Vector Hysteresis Model Using an
  Isotropic Vector Play Model"
- US Patent 10830836 (lab IP): "Systems and Methods for Anisotropic
  Vector Hysteresis Analysis"

## Energy-based vector hysteresis (Egger, Henrotte, Steentjes)

Constitutive law is derived from a convex thermodynamic energy
functional U(M); evolution is M_dot = -dU/dM (with friction).
- Lab papers:
  - "An energy-based vector hysteresis model for ferromagnetic
    materials" (Henrotte)
  - "2006 An energy based vector hysteresis model"
  - "On the Energy-Based Variational Model for Vector M"
  - "A Vector Play Model for Finite-Element Eddy-Current Analysis
    Using the Newton-Raphson Method"
  - "Test-case transformer for the energy-based vector hysteresis
    model"

This is the model line directly used by Radia's
`MatPlayHysteresis` / `MatEnergyHysteresis` (see `CLAUDE.md`,
Hysteresis Materials section).

## Inverse Jiles-Atherton

Jiles-Atherton (1986) is scalar; vector extensions exist but are
non-convex. The "Inverse" form lets you solve B-input (FEM-friendly)
instead of H-input.
- Lab paper: "Inverse Jiles-Atherton Vector Hysteresis Model"

## Magneto-elastic coupling

Mechanical stress changes the BH curve (Villari effect for tension,
Joule effect for compression). MFL stress measurement uses this:
- Lab paper: "A Magneto-elastic Vector-Play Model Including
  Piezomagnetic Behavior"

## After-effect (magnetic aging)

Slow time-domain recovery of magnetization after a field change,
caused by domain-wall thermal activation.
- Lab papers:
  - "Aftereffect and Accommodation"
  - "Magnetic After Effect of Cold Rolled Iron" (1 and 2)
  - "Thermal fluctuation after-effect in soft magnetic materials"

Relevant for low-noise NDT where the probe stops and waits -- the
DC signal drifts as the iron relaxes.

## US patents in the lab library

- 10345398 "Systems and Methods for Vector Hysteresis Analysis"
- 10830836 "Systems and Methods for Anisotropic Vector Hysteresis
  Analysis"

These cover the lab's hysteresis identification + simulation IP.

## Lab reference

public-safe curated corpus (25 PDFs)
- Also `radia_mcp.magnetic_materials` for non-NDT hysteresis context
- Hysteresis models in Radia: `MatPlayHysteresis` / `MatEnergyHysteresis`
"""


# ---------------------------------------------------------------------------
# MFL (MAGNETIC FLUX LEAKAGE)
# ---------------------------------------------------------------------------

MFL_MAGNETIC_FLUX_LEAKAGE = r"""
# Magnetic Flux Leakage (MFL)

MFL is the dominant electromagnetic NDT method for ferromagnetic
infrastructure -- pipelines, storage-tank floors, train rails, ship
hulls, and wire ropes. For carbon steel, ECT is limited to ~0.2 mm
depth (skin effect), but MFL works on the full wall thickness because
it uses DC magnetization.

## Principle

1. Magnetize the steel to near saturation (B ~ 1.5 -- 2 T) using
   a permanent magnet yoke or DC electromagnet
2. Most of the flux flows through the steel
3. At a defect (corrosion pit, crack, wall thinning) the local
   permeability drops to mu_r = 1; the flux is forced to detour
4. Part of the flux escapes ("leaks") into the air above the surface
5. A sensor placed above the surface measures B_normal, B_axial, or
   B_transverse

The amplitude and shape of the leakage signal are correlated with
defect size, depth, and orientation.

## Sensor types

| Sensor | B range | Use |
|--------|---------|-----|
| **Hall effect** | mT -- T | DC + AC, low-cost, planar arrays |
| **Inductive coil** | wide | Velocity-dependent (PIG must move) |
| **GMR / TMR** | nT -- mT | Very high sensitivity, modern arrays |
| **Fluxgate** | nT -- mT | Lab-grade |

Hall arrays dominate industrial use. Lab paper:
"Nondestructive testing and crack evaluation of ferromagnetic
material by using the linearly integrated hall sensor array" (Kang,
J. Mech. Sci. Tech. 2008) describes a 64-element Hall array for
ferromagnetic plate crack inspection.

## Magnetization yoke design

- **Permanent-magnet yoke** -- NdFeB blocks + iron pole pieces;
  always-on, simple, but the magnetic field is fixed (cannot vary
  excitation level)
- **DC electromagnet yoke** -- coil-and-iron-yoke; tunable
  magnetization, heavier
- **Saddle yoke** -- two C-cores on the OD of a pipe, used in
  in-line PIG inspection
- **Inverse-magnetized yoke** -- two opposite-polarity magnets
  produce a vertical-B above the surface; common for pipeline PIG

The yoke must produce B ~ 1.5 -- 2 T in the wall to drive deep
defects above the leakage threshold. Below ~1 T the signal is
unreliable.

## Pipeline PIG (Pipeline Inspection Gauge)

A "smart pig" is a tool that travels through a pipeline carrying
the magnetizing yoke and sensor array. Standard PIG sizes:
- 6" (150 mm) -- 48" (1.2 m) diameter
- 4 -- 6 m/s transit speed (limited by sensor sampling and yoke
  saturation lag)
- 100 -- 1000 km in one run

PIG MFL is the industry's primary defense against pipeline
corrosion. False positives have been reduced by:
- Multi-sensor (axial + transverse + radial) MFL
- Sensor fusion with the same-PIG ultrasonic data
- ML defect classification (pit vs groove vs gouge)

## ID vs OD distinction

A challenge: an internal pit and an external pit can produce similar
leakage signals (the field leaks through the same wall opening).
Standard distinguishers:

- **Saturation level**: heavily-saturated yoke; ID and OD signals
  differ in phase / shape
- **Two-sensor stack**: near-wall sensor + far-wall sensor; OD pit
  saturates near sensor more
- **High-resolution velocity**: as the PIG passes, ID and OD signals
  have different rise/fall times

## Wire-rope MFL (mine elevator, crane cables)

Wire rope is steel-strand cable used in mining elevators, port
cranes, and aerial trams. MFL inspectors clamp around the cable,
saturate it longitudinally with a permanent-magnet head, and detect
broken individual strands by Hall signals at the cable surface.

Standard: International Standard ISO 4309 for crane wire rope
inspection.

## Rail MFL

Train rails are inspected by track-mounted MFL at typical scan
speeds 10 -- 50 km/h. Detects:
- Transverse cracks in the rail head (cause of derailment if
  unchecked)
- Squat / shelling defects (rolling contact fatigue)
- Internal cracks at the web-foot transition

## Tank floor MFL

Above-ground storage tank floors are scanned with hand-held or
crawler-mounted MFL arrays. Targets corrosion pitting from the
soil-side (bottom) of the floor plate, which is invisible from above.

## Comparison with ECT

| Feature | MFL | ECT |
|---------|-----|-----|
| Depth | Full wall (10 -- 30 mm typical) | One skin depth |
| Speed | 5 -- 50 km/h (rail), 4 -- 6 m/s (PIG) | 100 mm/s (manual) |
| Power | High (DC magnet) | Low (mA coil) |
| Ferro target | Yes -- exclusive | No -- not for carbon steel |
| Crack signal | Leakage amplitude | Impedance change |
| Sizing | Voltage / pattern fitting | Phase / amplitude |

## Lab reference

- public-safe curated corpus Hall-array paper
- public-safe curated corpus -- the hysteresis
  underpinning MFL nonlinearity
- `radia_mcp.electromagnet` -- DC yoke design via Hantila method
"""


# ---------------------------------------------------------------------------
# MICROWAVE / EMI / LANDMINE
# ---------------------------------------------------------------------------

MICROWAVE_EMI_LANDMINE = r"""
# Microwave, EMI Spectroscopy, and Metal Detection (Humanitarian)

The low-frequency / high-power end of electromagnetic NDT covers
metal detection, humanitarian demining, ground-penetrating radar
(GPR), and high-frequency microwave NDE for dielectrics.

## Metal detector (DC -- 100 kHz)

Operating principle: a search coil at the head of a handheld or
walk-through unit drives an AC field that induces eddy currents in
any nearby metal object. The disturbed field perturbs the coil
impedance, alerting the user.

Variants:
- **VLF (very-low-frequency)** ~ 3 -- 30 kHz; hobbyist gold detectors
- **PI (pulsed induction)** -- step-function drive, transient
  reception; immune to mineralized soil, depth ~ 1.5 m for car-key
  targets
- **MD-BFO (beat-frequency-oscillator)** -- legacy, declining use
- **Walk-through frame** -- 1 -- 5 kHz, multi-coil array

## Buried target ID by EMI spectroscopy

For humanitarian demining (anti-personnel mine clearance):
- A metal-detector pulse alone gives only yes/no; it cannot
  distinguish a shell casing (target) from a soda-can pull-tab
  (clutter)
- **Electromagnetic Induction Spectroscopy (EMI-S)** measures the
  frequency response of the buried object impedance from 100 Hz to
  100 kHz. Different objects have different relaxation spectra
  (poles of the impedance transfer function) that act as a
  signature.
- Lab paper: "Electromagnetic Induction Spectroscopy" -- the
  fundamental paper
- Lab paper: "Electromagnetic Induction Spectroscopy for Clearing
  Landmines" -- humanitarian application
- Lab paper: "Method and Procedure in Metal Detection and
  Distinctive Discrimination of Size, Shape, and Location of Hidden
  Single and Multiple Buried Objects"

A typical signature comprises 2 -- 4 magnetic poles with relaxation
times spanning 1 microsecond -- 1 ms. The pole pattern is invariant to
position, orientation (largely), and depth -- it depends only on
target geometry and material. The forward model is a low-frequency
multipole expansion; the inverse problem is a sparse-spectral fit.

## Ground-penetrating radar (GPR)

GPR uses a UWB pulse (centre frequency 100 MHz -- 3 GHz, bandwidth
~ 100%) to image subsurface structure. Targets:
- Utility location (buried gas / water / power lines)
- Concrete inspection (rebar mapping, void detection)
- Archaeology (buried walls, graves)
- UXO + landmine (large metal + dielectric contrast)

Penetration vs resolution:
- 100 MHz: 30 m in dry sand, 15 cm resolution
- 1 GHz: 1 -- 2 m in dry sand, 5 cm resolution
- 3 GHz: 0.3 m, 1 cm resolution

GPR signal processing:
- Migration (Kirchhoff or F-K) to focus hyperbolic reflections
- 3D imaging from multi-channel arrays
- Polarimetric (HH, HV, VH, VV) for target orientation

## Microwave NDE for dielectrics

For composites (CFRP, GFRP), polymers, and ceramics where eddy
currents are absent. Methods:
- **Reflection coefficient** at sample surface (1 -- 100 GHz)
- **Free-space transmission** (sample between two horn antennas)
- **Open-ended waveguide / coaxial probe** (point contact)
- **Cavity perturbation** (Q-change inside resonator)

Detects:
- Disbonds and delaminations (impedance discontinuity)
- Moisture (dielectric loss tan delta peaks at GHz)
- Cure state of resins (epsilon shift during polymerization)
- Fiber breakage in CFRP (polarization rotation)

## NMR-MOUSE for surface chemistry

Lab paper: "The NMR MOUSE, a Mobile Universal Surface Explorer"
(Eidmann, Savelsberg, Blumler, Blumich, J. Magn. Reson. 1996)

A portable single-sided NMR scanner that uses the strong B0 gradient
in the stray field of a permanent magnet to slice-select a thin
surface layer of the sample (~10 -- 100 micrometres deep). Inducts
a 10 -- 20 MHz RF field via a small surface coil and measures the
proton (or rubber, oil) free induction decay. Used for:

- Polymer crosslink density (T2 relaxation)
- Moisture content in concrete, wood
- Aging of rubber tires, conveyor belts
- Oil saturation in rock cores

Pure surface-NMR is non-eddy-current NDT, but is included in the
lab library as an adjacent technique.

## Lab reference

- public-safe curated corpus (5 PDFs -- metal detection +
  landmine EMI spectroscopy)
- public-safe curated corpus (23 PDFs --
  microwave scattering, MoM, FDTD)
- "Effects of building materials and structures on radiowave
  propagation above about 100 MHz" -- field-strength prediction in
  buildings (00_金属探知機/ folder)
- "Electromagnetic analysis of horizontal wire in two-layered soil"
  -- GPR / pipe-locator antenna analysis
"""


# ---------------------------------------------------------------------------
# ML FOR NDT
# ---------------------------------------------------------------------------

ML_FOR_NDT = r"""
# Machine Learning for NDT

Inverse problems in NDT (recover defect parameters from a signal) are
nonlinear and ill-posed. Classical methods rely on calibration tables
or model-based inversion with manual feature engineering. Modern
deep learning replaces both with end-to-end learned regression /
classification.

## Three ML problem classes

### Class 1: Detection (binary classification)

Decide flaw / no-flaw from a scan segment.
- Input: 1D impedance trace, 2D C-scan image, 3D volume scan
- Output: binary label
- Common: CNN on impedance-plane scatter plot images; XGBoost on
  hand-engineered features (peak, RMS, zero-crossings)

### Class 2: Sizing (regression)

Predict crack length, depth, width.
- Input: same as above
- Output: continuous vector (length, depth, opening)
- Common: CNN -> dense regression head; physics-informed loss
  (forward model evaluated for training-time consistency)

### Class 3: Reconstruction (image -> image)

Recover 2D / 3D defect shape from scan data.
- Input: 2D B(x,y) image from MFL or 1D impedance trace from ECT
- Output: 2D / 3D voxel grid of defect (binary or distance map)
- Common: U-Net / DenseNet; GAN-based super-resolution; diffusion
  models (recent)

## Training data sources

ML for NDT lives or dies by training data:

1. **Real specimens**: a fixed set of EDM-notch standards and
   destructive-tested service-aged samples. Hundreds to thousands of
   labelled scans.
2. **Synthetic FEM/BEM**: forward-model thousands of geometries
   under a parametric defect distribution. Cheap, but suffers from
   "sim-to-real" gap.
3. **Hybrid (transfer learning)**: pretrain on synthetic, fine-tune
   on real -- the dominant production paradigm.

## Domain-specific architectures

- **DeepImage** (Lemos 2018) -- CNN regression on ECT impedance
  images for crack length
- **MFL-Net** (Liu 2019) -- 1D CNN on hall-array MFL traces
- **Physics-informed neural networks (PINN)** -- enforce Maxwell's
  equations as a regularizer on the network

## Sugahara Lab connection

Lab library includes:
- "ディープラーニング、ビッグデータ、機械学習" (Asakawa, Shinyo
  Sha) -- general Japanese deep-learning textbook
- "General theory of higher-order decomposition of exponential
  operator" (Suzuki 1992) -- foundational paper on operator
  splitting; cited in physics-informed ML

Cross-references:
- `radia_mcp.bayesian_opt` -- BBO for ECT probe design
- `radia_mcp.gnn` -- graph NN for mesh-based forward models
- `radia_mcp.evolutionary` -- GA for inverse problem
- `radia_mcp.data_assimilation` -- Kalman filter for transient ECT

## Open challenges

1. **Sim-to-real**: forward simulators omit probe-cable artefacts,
   coil winding randomness, surface roughness; real-world deployment
   needs domain adaptation
2. **Uncertainty quantification**: industry codes demand
   probability-of-detection curves; deep nets must produce
   calibrated confidence
3. **Adversarial robustness**: a small probe-cable kink should not
   change the network's defect call
4. **Explainability**: regulatory bodies (FAA, NRC, ASME) need to
   audit the network's reasoning

## Lab reference

- public-safe curated corpus
- General textbooks in public-safe curated corpus
"""


# ---------------------------------------------------------------------------
# DISPATCHER
# ---------------------------------------------------------------------------

_SECTIONS: dict[str, str] = {
    "overview": OVERVIEW,
    "eddy_current_basics": EDDY_CURRENT_BASICS,
    "probe_types": PROBE_TYPES,
    "defect_models": DEFECT_MODELS,
    "signal_processing": SIGNAL_PROCESSING,
    "fem_formulation": FEM_FORMULATION,
    "integral_equation": INTEGRAL_EQUATION,
    "multiply_connected": MULTIPLY_CONNECTED,
    "open_boundary": OPEN_BOUNDARY,
    "hysteresis_in_ndt": HYSTERESIS_IN_NDT,
    "mfl_magnetic_flux_leakage": MFL_MAGNETIC_FLUX_LEAKAGE,
    "microwave_emi_landmine": MICROWAVE_EMI_LANDMINE,
    "ml_for_ndt": ML_FOR_NDT,
}


_ALIASES: dict[str, str] = {
    # convenient short forms
    "intro": "overview",
    "ect": "eddy_current_basics",
    "eddy_current": "eddy_current_basics",
    "skin_depth": "eddy_current_basics",
    "probes": "probe_types",
    "coil": "probe_types",
    "defects": "defect_models",
    "crack": "defect_models",
    "back_side": "defect_models",
    "signal": "signal_processing",
    "impedance_plane": "signal_processing",
    "fem": "fem_formulation",
    "a_v": "fem_formulation",
    "t_omega": "fem_formulation",
    "edge_elements": "fem_formulation",
    "bem": "integral_equation",
    "ie": "integral_equation",
    "vie": "integral_equation",
    "multiconnected": "multiply_connected",
    "cut": "multiply_connected",
    "tree_cotree": "multiply_connected",
    "pml": "open_boundary",
    "kelvin": "open_boundary",
    "abc": "open_boundary",
    "hysteresis": "hysteresis_in_ndt",
    "play_model": "hysteresis_in_ndt",
    "stop_model": "hysteresis_in_ndt",
    "jsaem": "hysteresis_in_ndt",
    "mfl": "mfl_magnetic_flux_leakage",
    "pig": "mfl_magnetic_flux_leakage",
    "rail": "mfl_magnetic_flux_leakage",
    "pipeline": "mfl_magnetic_flux_leakage",
    "metal_detector": "microwave_emi_landmine",
    "gpr": "microwave_emi_landmine",
    "landmine": "microwave_emi_landmine",
    "emi": "microwave_emi_landmine",
    "nmr_mouse": "microwave_emi_landmine",
    "ml": "ml_for_ndt",
    "deep_learning": "ml_for_ndt",
    "inverse": "ml_for_ndt",
}


def get_knowledge(topic: str = "overview") -> str:
    """Dispatch NDT knowledge by topic.

    Args:
        topic: One of TOPICS keys, or a shorthand alias.

    Returns:
        Markdown-formatted knowledge block, or an error message
        listing valid topics if the given topic is unknown.
    """
    if topic is None:
        topic = "overview"
    key = topic.lower().strip()
    if key in _ALIASES:
        key = _ALIASES[key]
    if key == "all":
        ordered = [
            "overview",
            "eddy_current_basics",
            "probe_types",
            "defect_models",
            "signal_processing",
            "fem_formulation",
            "integral_equation",
            "multiply_connected",
            "open_boundary",
            "hysteresis_in_ndt",
            "mfl_magnetic_flux_leakage",
            "microwave_emi_landmine",
            "ml_for_ndt",
        ]
        return "\n\n".join(_SECTIONS[k] for k in ordered)
    if key in _SECTIONS:
        return _SECTIONS[key]
    return (
        f"Unknown topic '{topic}'. Valid topics: "
        + ", ".join(sorted(_SECTIONS.keys()))
        + ", all"
    )
