"""Fusion magnet knowledge.

Sourced from peer-reviewed papers and design reports archived at
W:/03_文献・論文/00_電磁界解析/99_アプリケーション/08_核融合/.  Each
topic constant ends with a "References" block listing the specific
PDFs / journal articles used to populate it, so future contributors
(and the AI assistant) can re-open the source instead of
re-discovering numbers from scratch.
"""

# Single source of truth for the dispatcher's accepted topic enum.
# Wired into the `fusion_reactor_topics()` MCP tool via
# common.register_topics_tool -- see server.py.
TOPICS: dict[str, str] = {
    "overview":         "Confinement landscape (tokamak / stellarator / FRC / "
                        "inertial), operating regime (5-13 T, DC + disruption "
                        "transients), Radia/NGSolve role in fusion EM analysis",
    "tokamak":          "ITER TF/PF/CS coil system, Nb3Sn/NbTi cable-in-conduit, "
                        "JT-60SA error field correction, disruption forces, "
                        "halo currents, vacuum vessel eddies",
    "stellarator":      "Modular vs helical (LHD, W7-X), 3D twisted toroidal "
                        "confinement, quasi-symmetry classes (QA / QI / QH), "
                        "stage-1 / stage-2 optimization workflow",
    "iter_coil_system": "ITER TF/PF/CS/CC coil parameters (18 TF, 6 PF, 6 CS, "
                        "18 CC), Nb3Sn cable-in-conduit conductor (CICC), "
                        "operating B-field, force levels, joint resistance",
    "w7x_modular":      "Wendelstein 7-X engineering: 50 non-planar + 20 planar "
                        "NbTi coils, narrow support elements (NSE), lateral "
                        "support elements (LSE), 4 K cryoplant, 18-year build",
    "lhd_helical":      "Large Helical Device (NIFS, Japan): continuous helical "
                        "winding pole pair, NbTi forced-flow CICC, Mitsubishi "
                        "fabrication lineage, operation since 1998",
    "sc_cable":         "Cable-in-conduit conductor (CICC) for fusion magnets: "
                        "NbTi (4.5 K, B < ~8 T) vs Nb3Sn (4.5 K, B up to 13 T) "
                        "vs HTSC current leads, supercritical He cooling",
    "stellarator_coil_design": "Coil design codes (NESCOIL/NESVD/REGCOIL/COILOPT/"
                        "FOCUS), Biot-Savart current potential, Hessian "
                        "sensitivity analysis, error field tolerance",
    "error_field":      "Manufacturing/assembly tolerances -> non-axisymmetric "
                        "B perturbations -> locked modes / disruptions; "
                        "Monte-Carlo error-field prediction; EFCC design "
                        "(JT-60SA, ITER); m/n=2/1 mode control",
    "transient_eddy":   "Kameari thin-conductor current-function FEM (1981); "
                        "Cariddi-type T-Omega/A-V integral formulations; "
                        "QR-recursive H-matrix compression for ITER-scale "
                        "passive structures; disruption-induced eddy forces",
    "rmp_elm_control":  "Resonant Magnetic Perturbation (RMP) for ELM "
                        "suppression: in-vessel coils, stochastic edge layer, "
                        "n=3 perturbation, separatrix island chains",
    "all":              "Concatenated dump of every topic above",
}


OVERVIEW = r"""
# Fusion reactor magnet systems

This subpackage collects working knowledge -- physics regime, design
numbers, formulations -- that Radia/NGSolve users hit when they
analyze fusion-class magnets.  The aim is not to replace the IAEA
"Encyclopaedia of Fusion" but to provide the *minimum quantitative
context* a magnet engineer or applied mathematician needs before
opening a panel.

## Confinement principles

| Type                      | Plasma confinement                 | Magnet system                          |
|---------------------------|-------------------------------------|-----------------------------------------|
| **Tokamak**               | Toroidal B + induced poloidal      | TF + PF + CS coils                     |
| **Stellarator**           | Fully 3D twisted toroidal B        | Modular non-planar / helical coils     |
| **Spherical tokamak (ST)**| Low aspect ratio toroidal          | TF center-post + PF                    |
| **FRC** (Field-Reversed)  | Compact toroid, theta-pinch        | Theta-pinch + mirror coils             |
| **Mirror** (open)         | Magnetic mirror trap (open)        | End mirror coils + central solenoid    |
| **Inertial Confinement**  | None (laser/X-ray ablation)        | None (drivers are lasers/Z-pinch)      |

Magnetic confinement (tokamak + stellarator) is the dominant magnet
engineering challenge of the field and the focus of this subpackage.
Inertial confinement is out of scope.

## Operating regime

- **B field at plasma**: 2.5 T (W7-X) to 5.3 T (ITER) to ~12 T
  (compact tokamak proposals SPARC / CFETR).  Conductor peak fields
  are higher (ITER TF peak 11.8 T, CS peak 13 T).
- **Frequency**: DC for steady state.  Slow ramp (~10 s plasma current
  ramp-up, ~tens of ms scenario evolution).
- **Disruption transient**: millisecond plasma current drop induces
  large eddy currents in the vacuum vessel + cryostat + first wall,
  generating MN-class forces on the structure (ITER vacuum vessel:
  some hundreds of tonnes net force during major disruption).
- **AC ripple**: TF ripple from discrete coils (~18 in ITER, 50 in
  W7-X non-planar) -- spatial harmonic, not temporal.

## Radia / NGSolve / Cariddi role in fusion EM analysis

| Analysis                         | Tool family                                                |
|----------------------------------|------------------------------------------------------------|
| Coil shape optimization          | NESCOIL, REGCOIL, COILOPT++, FOCUS (Princeton)             |
| Plasma equilibrium / MHD         | VMEC, NIMROD, M3D-C1                                       |
| 3D coil B-field (filament)       | Biot-Savart -- can be done with `radia_mcp.peec` filaments |
| Iron yoke / passive ferritic     | `radia_mcp.electromagnet` Hantila path                      |
| Vacuum vessel eddy current       | Cariddi (CREATE), `radia_mcp.fem.time_domain_axisym`        |
| Error-field forecasting (Monte)  | Same as eddy current, with Monte-Carlo coil perturbations  |
| H-matrix acceleration            | `radia_mcp.bem.h_matrix` (ACA+) for dense BEM/integral     |

Radia's MMM/MSC is **not** the natural first tool for fusion; the
field is dominated by Cariddi-type integral eddy current codes
(Rubinacci, Villone et al., CREATE Consortium).  Radia is useful
when (a) the geometry has soft iron (TBM modules, magnetic shields,
NB ferritic boxes) and (b) the user wants natural open-boundary
field evaluation without PML.

## Cross-references

- `radia_mcp.accelerator` -- adjacent superconducting magnet research
  (LHC dipole, undulators); shares Nb3Sn / NbTi cable physics
- `radia_mcp.electromagnet` -- DC magnet design (Hantila polarization
  for nonlinear iron, e.g. JT-60SA NB ferritic shielding)
- `radia_mcp.fem.potential_formulations` -- HCurl A and H formulations
  for SC eddy current
- `radia_mcp.bem.h_matrix` -- ACA / H-matrix for ITER-scale BEM
- `radia_mcp.matrix_solvers.preconditioners.amg` -- AMG for huge
  ITER-scale linear systems

## Bibliography (entry-level)

- Wesson, J. "Tokamaks" (Oxford, 4th ed. 2011) -- standard textbook
- Stacey, W.M. "Fusion: An Introduction to the Physics and Technology
  of Magnetic Confinement Fusion" (Wiley, 2nd ed. 2010)
- Iiyoshi, A. et al. "Overview of the Large Helical Device project"
  Nucl. Fusion 39 (1999) 1245
- Bosch, H.-S. et al. "Construction of Wendelstein 7-X -- Engineering
  a Steady-State Stellarator" IEEE Trans. Plasma Sci. 38 (2010) 265
"""


TOKAMAK = r"""
# Tokamak (ITER, JT-60SA and similar)

The tokamak is the axisymmetric magnetic confinement device with
the largest installed base and largest current investment (ITER,
EUR ~20 B construction).  Plasma current is driven inductively by
the central solenoid (transformer-like primary), generating the
poloidal field that combines with the externally-applied toroidal
field to twist field lines and confine charged particles.

## Coil system

| Coil               | Function                                | Field at plasma / peak | Conductor (ITER) |
|--------------------|------------------------------------------|------------------------|------------------|
| **TF** (Toroidal)  | Main confining B; sets sigma_axisymmetry| 5.3 T / 11.8 T peak    | Nb3Sn CICC       |
| **PF** (Poloidal)  | Plasma position / vertical stability    | varies (~kT level)     | NbTi CICC        |
| **CS** (Central Solenoid) | Inductive Ip drive (flux swing)  | up to 13 T peak        | Nb3Sn CICC       |
| **Correction (CC)**| Error field correction                  | ~10 mT order           | NbTi CICC        |
| **In-vessel ELM control** | RMP for edge localized modes     | ~mT at plasma edge     | water-cooled Cu  |
| **Test Blanket Module (TBM)** | Tritium breeding qualification | local ferritic | -                |

For details on ITER TF/CS conductor design see `iter_coil_system`.
For correction coil design and error-field budget see `error_field`.

## ITER specifics (Cadarache, France)

35-nation fusion demonstrator (EU/JP/US/RU/CN/KR/IN); construction
~2010-2030s, first plasma originally 2025, now revised to ~2034
following the 2024 baseline revision.

| Parameter                | Value                          |
|--------------------------|--------------------------------|
| Major radius R           | 6.2 m                          |
| Minor radius a           | 2.0 m                          |
| Plasma volume            | 840 m^3                        |
| Plasma current Ip        | 15 MA (nominal)                |
| Toroidal B on axis       | 5.3 T                          |
| Fusion power Pfus        | 500 MW (Q = 10 target)         |
| Burn duration            | ~400 s nominal, 1000 s long pulse |
| Total mass               | ~23 000 t                      |
| TF coil count            | 18                             |
| PF coil count            | 6 (PF1 - PF6)                  |
| CS module count          | 6 (stacked, independently powered) |
| Correction coil sets     | 18 (6 top + 6 side + 6 bottom) |

(Source: ITER Organization Public Web, JADA/QST "ITER計画" pamphlet
W:/.../01_ITER/Iter計画.pdf)

## JT-60SA (Naka, Japan)

EU-JP Broader Approach satellite tokamak, supports ITER scientifically
and demonstrates steady-state advanced scenarios:

| Parameter                | Value                          |
|--------------------------|--------------------------------|
| Major radius R           | 2.96 m                         |
| Plasma current Ip        | 5.5 MA (LN scenario)           |
| Toroidal B on axis       | 2.25 T                         |
| TF coil count            | 18 (NbTi)                      |
| PF coil count            | 10 (4 CS + 6 EF)               |
| In-vessel EFCC           | 18 coils, 3 rows x 6 (Cu, 30/45 kA-turns) |

(Reference: Matsunaga et al. "In-vessel coils for magnetic error
field correction in JT-60SA" Fusion Eng. Des. 98-99 (2015) 1113.)

JT-60SA is also the testbed for the W7-X / ITER HTSC current lead
technology (Fietz et al. IEEE TAS 2009).

## Disruption physics and EM consequences

When plasma confinement collapses (millisecond scale, MHD instability
or impurity-driven radiative collapse):

1. **Thermal quench** (TQ, ~1 ms): plasma energy dumps to first wall
   via radiation + conduction.
2. **Current quench** (CQ, ~10-100 ms in ITER): Ip drops, dIp/dt
   induces eddy currents in:
   - vacuum vessel (304L SS, ~60 mm thick double wall)
   - cryostat (304L SS shell)
   - blanket back-plate
   - TF coil casings (eddy heats casing)
3. **Halo currents**: poloidal current that flows partly through
   the plasma and partly through the wall, generating large local
   forces (Halo Toroidal Peaking Factor ~2-3).
4. **Runaway electrons**: relativistic electron beam carrying tens
   of percent of pre-disruption Ip; localized damage on first wall.

ITER disruption-induced force budget:
- Vacuum vessel: net vertical force ~150-200 t during VDE (Vertical
  Displacement Event)
- TF coil casing eddy heating ~ kJ-scale per disruption
- ITER design: ~3000 disruptions over machine lifetime allowed

Lab cross-reference: transient eddy currents in passive 3D
structures are computed with `radia_mcp.fem.time_domain_axisym.
time_domain_fem` (axisymmetric) or, for full 3D with periodicity,
Cariddi-style QR-recursive H-matrix (`transient_eddy` topic).

## TF ripple

With finite number N_TF of toroidal field coils, the toroidal B is
not perfectly axisymmetric: there is a residual N_TF-periodic ripple.

ITER (N_TF = 18): peak-to-peak ripple at outer plasma edge
  ~0.5% (acceptable for fast-ion confinement after ferritic insert
  panels are installed at LFS).
JT-60SA (N_TF = 18): ~1% peak-to-peak.

TF ripple is mitigated by ferritic inserts in the vacuum vessel
shadow.  This is one of the few places Radia's MMM (nonlinear iron)
is genuinely useful for tokamak design.

## References

- ITER Organization Public Web (https://www.iter.org)
- "ITER計画" pamphlet, QST / National Institutes for Quantum and
  Radiological Science (Naka), 2017
  W:/.../01_ITER/Iter計画.pdf
- Matsunaga, G. et al. "In-vessel coils for magnetic error field
  correction in JT-60SA" Fusion Eng. Des. 98-99 (2015) 1113-1117.
  doi:10.1016/j.fusengdes.2015.06.024
  W:/.../research/InVesselCoils_JT60SA_ErrorField.pdf
- Hender, T.C. et al. "Chapter 3: MHD stability, operational limits
  and disruptions" Nucl. Fusion 47 (2007) S128 -- canonical error
  field threshold reference (TMEI metric)
- 2017 Ohm 12月号 "ITER特集" (Japanese review issue)
  W:/.../01_ITER/2017_Ohm_12月号.pdf
- 原型炉_誤差磁場_補正コイル_日立 (Hitachi DEMO error field study)
  W:/.../research/原型炉_誤差磁場_補正コイル_日立.pdf

## Cross-references

- `radia_mcp.electromagnet` -- Hantila polarization for nonlinear
  ferritic inserts (TF ripple mitigation, NB ferritic shielding)
- `radia_mcp.team_benchmark.eddy_current.problem_24` -- transient
  nonlinear eddy current benchmark (analogue of disruption eddies)
- `radia_mcp.bem.h_matrix` -- H-matrix compression for vacuum vessel
  integral formulation
- `radia_mcp.fem.time_domain_axisym` -- axisymmetric transient FEM
  (good first pass for VDE / centered disruption)
"""


STELLARATOR = r"""
# Stellarator -- twisted toroidal confinement

Stellarators achieve confinement **without** a net plasma current
(unlike tokamak).  The 3D rotational transform is generated entirely
by stationary external coils.  This buys steady-state operation and
intrinsic immunity to disruptions, at the cost of dramatically more
complex coil geometry.

## Coil topology classes

| Class                     | Examples              | Coil layout                                  |
|---------------------------|------------------------|---------------------------------------------|
| **Classical stellarator** | Wendelstein 7-AS, NCSX| External helical windings + TF              |
| **Heliotron / Torsatron** | LHD (NIFS Japan)      | Continuous helical pole pair, no TF needed  |
| **Heliac**                | TJ-II                 | Circular conductor with helical center      |
| **Modular**               | W7-X, HSX, NCSX (cancelled) | Discrete non-planar coils                 |
| **Quasi-axisymmetric (QA)**| CFQS, NCSX            | Modular coils designed for B(theta) symmetric in Boozer coords |
| **Quasi-helically sym (QH)**| HSX                  | Modular coils, helical B symmetry           |
| **Quasi-isodynamic (QI)** | W7-X (approximate)    | Modular, optimized for low neoclassical    |

The distinction matters because **the symmetry class determines
neoclassical transport**: any continuous symmetry of B in magnetic
coordinates implies a conserved canonical momentum and hence
intrinsically small radial drift of trapped particles.  Tokamak is
the trivial case (axisymmetry); stellarators must engineer one of
QA / QH / QI by coil shaping.

(Source: Landreman, Simons Summer School 2020; reviewed in
Zhu et al. "Towards simpler coils for optimized stellarators"
W:/.../research/Simpler_Coils_Stellarators.pdf)

## Stage-1 / Stage-2 optimization workflow

Modern stellarator design splits into two coupled optimization
stages, iterated until consistent:

**Stage 1 -- plasma equilibrium optimization** (free-boundary VMEC):
- 0D parameters (R, a, B0, beta, volume)
- Neoclassical transport (target: quasi-symmetric B-spectrum)
- Ideal MHD stability (Mercier, ballooning)
- Energetic particle confinement (alpha losses < ~5%)
- Bootstrap current minimization (for QI configurations)
- Divertor topology (island divertor for W7-X-type)

**Stage 2 -- coil optimization** (NESCOIL / REGCOIL / FOCUS):
- Buildable coils that reproduce Stage-1 plasma boundary B.n = 0
- Coil-to-coil spacing for diagnostic ports
- Coil-to-plasma distance for blankets (in a reactor)
- Maximum curvature, minimum length (cost / fabricability proxies)
- Linking number = number of toroidal field-periods

Coil complexity is largely determined by the Stage-1 physics
requirements -- you cannot wish away non-planar coils if you want
high beta + good confinement.  Engineering design (cable, casing,
support) comes after Stage-2 converges.

(Source: Zhu et al. Simpler_Coils_Stellarators.pdf, slide 5;
Garabedian PNAS 2007.)

## LHD (Large Helical Device, NIFS Toki, Japan)

The largest heliotron-class device in operation.

| Parameter                  | Value                          |
|----------------------------|--------------------------------|
| Major radius R             | 3.9 m                          |
| Minor radius a             | 0.5-0.65 m (adjustable)        |
| Toroidal field on axis     | 2.5-3.0 T                      |
| Helical pole count         | l = 2 (two continuous helical windings) |
| Toroidal field period      | m = 10                         |
| Helical conductor          | NbTi CICC, ~13 kA              |
| Three pairs of PF rings    | Inner Vertical (IV), Inner Shaping (IS), Outer Vertical (OV) |
| First plasma               | 1998-03-31                     |

LHD's signature design choice is the **continuous helical winding**:
the helical pole pair was wound directly onto the cryostat support
structure rather than assembled from discrete modules.  This avoided
the W7-X assembly tolerance problem but locked the geometry --
LHD cannot be "re-tuned" to a different rotational transform
profile without replacing the helical pole, only by adjusting the
ratio of helical and vertical field currents.

LHD demonstrated:
- 1 hour 1 MW steady-state plasma (2005-2008 campaigns)
- ne = 10^21 m^-3 super-dense core (NIFS world record class)
- Deuterium operation since 2017

Mitsubishi was the prime contractor for LHD's helical pole and PF
coil assembly.  See `10_三菱_ヘリカル型/` folder for the design
genealogy (LHD plus the "Mitsubishi heliotron" lineage of smaller
experimental devices).

## W7-X (Wendelstein 7-X, IPP Greifswald)

The largest and most fully optimized stellarator built to date.
Construction took 18 years (1996-2014); first plasma 2015-12-10.

| Parameter                  | Value                          |
|----------------------------|--------------------------------|
| Average major radius       | 5.5 m                          |
| Average minor radius       | 0.53 m                         |
| Magnetic axis              | Helical (3D), not planar       |
| Field on axis (nominal)    | 2.5 T (designed for 3 T)       |
| Rotational transform iota  | 0.72 < iota/(2 pi) < 1.25     |
| Field-period count         | 5 (pentagon symmetry)          |
| Non-planar coils (NPC)     | 50 (10 modules x 5 different types) |
| Planar coils (PC)          | 20 (2 different types)         |
| Coil conductor             | NbTi CICC in Al-alloy conduit  |
| Cryoplant capacity         | ~7 kW @ 4.5 K                  |
| ECRH heating               | 10 MW @ 140 GHz, CW (10 gyrotrons) |
| Divertor                   | Island divertor (m=5 island chain) |
| Pulse length goal          | 30 min @ 10 MW ECRH            |

W7-X is **quasi-isodynamic** (QI), an approximate symmetry that
minimizes bootstrap current and neoclassical transport
simultaneously.  Its non-planar coils were the engineering bottleneck:
108-turn NbTi cable wound around a 3D template, embedded in a cast
steel casing.  Coil position tolerance was 1 mm / 5 m, requiring
the narrow support elements described in `w7x_modular`.

(Source: Bosch et al. IEEE Trans. Plasma Sci. 38 (2010) 265,
W:/.../11_三菱_モジュラー型/Construction of Wendelstein 7-X.pdf)

## Other / next-generation stellarators

- **HSX** (Wisconsin) -- quasi-helically symmetric, 4 field periods,
  R=1.2 m, B=1 T.  First QH demonstrator.
- **CFQS** (Chengdu, China / NIFS collab) -- Chinese First Quasi-
  axisymmetric Stellarator; modular coils under fabrication 2020s.
  (Liu Nucl. Fusion 2020)
- **Helicity / Type One / Renaissance Fusion** -- private-sector
  stellarator startups using simpler coils + HTS.  See
  `企業による核融合研究の最近の動向` (Japanese review of fusion
  startups, W:/.../10_三菱_ヘリカル型/).

## References

- Alexander, R. and Garabedian, P. "Choice of coils for a fusion
  reactor" PNAS 104 (2007) 12250-12252.
  W:/.../research/Choice_of_coils_fusion.pdf
- Zhu, C. et al. "Towards simpler coils for optimized stellarators"
  Talk slides (Princeton PPPL).
  W:/.../research/Simpler_Coils_Stellarators.pdf
- Zhu, C. et al. "New method to design stellarator coils without
  the winding surface" Nucl. Fusion 58 (2018) 016008 -- FOCUS code
- Bosch, H.-S. et al. "Construction of Wendelstein 7-X -- Engineering
  a Steady-State Stellarator" IEEE Trans. Plasma Sci. 38 (2010) 265
  W:/.../11_三菱_モジュラー型/Construction of Wendelstein 7-X.pdf
- LHDプラズマのモデリングと理論解析 (LHD plasma modelling and
  theoretical analysis), NIFS
  W:/.../10_三菱_ヘリカル型/LHDプラズマのモデリングと理論解析.pdf
- ヘリカル型核融合炉設計 (Helical fusion reactor design)
  W:/.../10_三菱_ヘリカル型/ヘリカル型核融合炉設計.pdf
- モジュラーコイルを用いたプラズマの閉じ込め (Plasma confinement
  using modular coils)
  W:/.../11_三菱_モジュラー型/モジュラーコイルを用いたプラズマの閉じ込め.pdf

## Cross-references

- `stellarator_coil_design` -- the math of going from plasma B
  boundary to deliverable coil filaments
- `error_field` -- coil position tolerance budgets
- `w7x_modular` -- W7-X engineering deep-dive (NSE, LSE, cryoplant)
- `lhd_helical` -- LHD design and Mitsubishi heliotron lineage
"""


ITER_COIL_SYSTEM = r"""
# ITER coil system: numbers and conductors

This topic pins the actual specifications a magnet engineer needs
when sizing an analysis against the ITER reference design.  All
numbers from the ITER Organization Public Web and the JADA/QST
"ITER計画" pamphlet (2017).

## TF (Toroidal Field) coils

| Parameter                       | Value                          |
|---------------------------------|--------------------------------|
| Coil count                      | 18                             |
| Conductor                       | Nb3Sn cable-in-conduit (CICC)  |
| Operating temperature           | 4.5 K (supercritical He, 5 bar)|
| Operating current               | 68 kA                          |
| Stored energy (full system)     | 41 GJ                          |
| Peak field on conductor         | 11.8 T                         |
| Field on plasma axis            | 5.3 T (Bt0)                    |
| Coil overall height             | ~16.5 m                        |
| Coil overall width              | ~9 m                           |
| Single coil mass                | ~310 t (winding + casing)      |
| Casing material                 | 316LN-IG forging               |
| Casing mass                     | ~200 t per coil                |
| Winding pack mass               | ~110 t per coil                |
| Out-of-plane support            | Inner Intercoil Structure (IIS)|
| Total TF conductor length       | ~83 km                         |

Japan (JADA / QST) supplies 25% of the TF conductor and 100% of
the 18 TF casings (cast/forged 316LN-IG, machined at Mitsubishi
Heavy Industries Futami Works).  EU supplies the remaining 75%
TF conductor + 9 of 18 winding packs; Japan supplies the other 9
TF winding packs (47% of "lower winding integration").

## CS (Central Solenoid)

| Parameter                       | Value                          |
|---------------------------------|--------------------------------|
| Module count                    | 6 (stacked, independent power) |
| Conductor                       | Nb3Sn CICC                     |
| Operating current               | 40 kA                          |
| Peak field on conductor         | 13.0 T                         |
| Total flux swing                | ~280 V-s                       |
| Module height                   | ~2.1 m each                    |
| Total stack height              | ~12.5 m                        |
| Outer diameter                  | ~4.2 m                         |
| Total CS mass                   | ~1000 t                        |
| Manufacture                     | US (General Atomics)           |
| Conductor                       | JADA 100% supply               |

The CS is the "transformer primary" that inductively drives plasma
current.  Operating peak field 13 T at 4.5 K is at the edge of
Nb3Sn's capability and required a multi-year qualification program
(2012-2020) including the CSIC (CS Insert Coil) test at the Naka
Test Facility, Japan.

## PF (Poloidal Field) coils

| Parameter                       | Value                          |
|---------------------------------|--------------------------------|
| Coil count                      | 6 (PF1 through PF6, ring coils)|
| Conductor                       | NbTi CICC                      |
| Operating temperature           | 4.5 K                          |
| Operating current               | ~45 kA (varies by coil)        |
| Peak field on conductor         | ~6 T                           |
| Largest coil diameter (PF3, PF4)| ~24 m                          |
| Largest coil mass               | ~400 t                         |
| Manufacture                     | EU + Russia (China for PF6)    |

PF coils control plasma position, vertical stability, and shape
(elongation, triangularity, X-point position).  Because they
operate at lower peak field than TF/CS, NbTi is sufficient and
much cheaper than Nb3Sn.

## CC (Correction Coils)

| Parameter                       | Value                          |
|---------------------------------|--------------------------------|
| Set count                       | 18 (6 top + 6 side + 6 bottom) |
| Conductor                       | NbTi CICC                      |
| Operating current               | up to 10 kA                    |
| Purpose                         | Cancel m/n=2/1 and m/n=3/2 error fields |
| Field at plasma                 | ~10 mT order (small correction)|

The CC system is the active half of the ITER error-field control
loop.  Tolerance to error field is ~5x10^-5 * B_T at the q=2
surface for a 2/1 mode -- see `error_field` topic for the analysis
methodology.

## Cable-in-conduit conductor (CICC) physics

ITER all-CICC magnet system uses the same conductor architecture
adapted to NbTi (PF, CC) or Nb3Sn (TF, CS):

- Cable: 1000-1500 superconducting strands twisted in 5-6
  cabling stages (sub-cable / petal / final cable)
- Wrapped in stainless-steel jacket forming a square or round
  conduit
- Supercritical He flows through the void fraction (~30%) for
  cooling
- Joints between segments: < 1 nano-ohm resistance, cooled by He flow
- Strain on Nb3Sn: limited to ~0.3% by jacket reinforcement
  (Nb3Sn is brittle and strain-sensitive; Jc drops with strain)

See `sc_cable` topic for the NbTi vs Nb3Sn vs HTS trade-off.

## Plasma-facing components (relevant to EM analysis)

- **Vacuum vessel**: 316L SS double wall, 60 mm typical thickness,
  ~5200 t.  Dominant passive eddy current path during disruption.
- **Blanket modules**: 440 modules, each ~4.5 t, water-cooled
  CuCrZr / Be.  Eddy currents in the blanket back-plate are
  significant.
- **Divertor cassettes**: 54 cassettes at the bottom of vacuum
  vessel, W mono-block targets.
- **Cryostat**: 316L SS, single wall, ~30 m diameter x 30 m tall,
  ~3800 t.  Slow eddy currents during ramp-up.

## References

- ITER Organization Public Web (https://www.iter.org/mach)
- "ITER計画" pamphlet, JADA / QST (2017)
  W:/.../01_ITER/Iter計画.pdf
- Mitchell, N. et al. "The ITER Magnet System" IEEE Trans. Appl.
  Supercond. 18 (2008) 435 -- canonical reference for ITER magnets
- Devred, A. et al. "Status of ITER conductor development and
  production" IEEE Trans. Appl. Supercond. 22 (2012) 4804909
- ITER Research Plan within the Staged Approach, ITR-18-003 (2018)

## Cross-references

- `tokamak` -- broader context
- `sc_cable` -- CICC architecture in detail
- `error_field` -- where the CC current targets come from
- `transient_eddy` -- disruption-induced VV eddy current
"""


W7X_MODULAR = r"""
# Wendelstein 7-X engineering

W7-X is the engineering tour-de-force of magnetic-confinement fusion:
50 individually-fabricated non-planar superconducting coils held to
sub-mm tolerance over a 16 m pentagon, cooled to 4 K, and run for
30 min steady state.  Built by IPP Greifswald + a German-Italian-
British industrial consortium 1996-2014, first plasma 2015-12-10.

## Magnet system topology

Pentagon symmetry: 5 identical modules, each module = 2 flip-symmetric
half-modules.  Per half-module:

- 5 non-planar coils (NPC), each from one of 5 design types
- 2 planar coils (PC), each from one of 2 design types

Total: 50 NPC (5 types x 10 of each) + 20 PC (2 types x 10 of each)
= 70 superconducting coils.

The 5 NPC types per half-module are the irreducible 3D coil shapes
required to generate the W7-X target equilibrium.  The 2 PC types
add adjustable axisymmetric flux to tune the magnetic configuration
(rotational transform, iota; magnetic axis radius).

## Conductor and winding pack

- Conductor: NbTi cable in Al-alloy conduit
- 108 windings per NPC winding pack (6 cable lengths joined in series)
- Joint resistance: < 1 nOhm (low-resistance solder-type joints,
  also serve as He inlet/outlet ports)
- Cable cooled by supercritical He through the void fraction
- Winding pack embedded in cast 316LN-IG steel casing for force
  reaction
- Each NPC tested individually at CEA Saclay at 4 K, full current,
  before integration (49 of 50 tested as of Bosch et al. 2010)

NPC supplier: Babcock Noell (Germany) + Ansaldo Superconduttori
(Italy) consortium.
PC supplier: Tesla Engineering Ltd (UK).

## Coil support structure

The mechanical challenge of W7-X is the force network between the
50 NPC.  Each coil experiences:
- Radial Lorentz force up to 4.4 MN (pulling coil away from plasma)
- Bending moment up to 450 kN-m (3D coil shape -> twisted force)
- Inter-coil forces up to 1.5 MN compression on the NSE

Three classes of structural element:

### Central Support Ring (CSR)

Pentagon-shaped welded ring of 5 modules, each from 2 half-modules.
Cast 316LN-IG steel extensions weld onto the ring; coils bolt onto
the extensions via long, slender Inconel bolts and sleeves that
allow slight flange opening to reduce peak stress.  Each NPC bolts
to the CSR in 2 places.

### Lateral Support Elements (LSE)

Bolted Inconel "bridges" on the **outboard** side of the torus,
between neighbouring coils within a module.  Rigid connection;
the engineering challenge was control of weld shrinkage and
distortion to comply with assembly tolerances.

### Narrow Support Elements (NSE)

**Sliding** elements on the **inboard** side, where coils approach
each other to within a few cm.  Take up to 1.5 MN compression
while permitting up to 5 mm relative sliding and 1 degree tilt
during magnet energization.  300 NSE in total (3-7 per coil pair,
depending on coil pair type).

NSE design (Heinemann et al. 2006, W:/.../11_三菱_モジュラー型/
Design of Narrow Support Elements...pdf):
- Central pad: Al-bronze alloy (DIN 2.0966 / 2.0978), 60 mm or
  73 mm diameter depending on load class
- Sliding surface roughness Rz = 0.7 microm (mirror polish)
- Lubrication: 6 microm MoS2 PVD coating + burnished MoS2 powder
  on counter-side
- Dust cover: Al-bronze, MoS2-coated
- Wave spring: Inconel 718, presses dust cover against counter-side
- Initial gaps 0-4.5 mm at 4 K zero current; close as coils energize
- Dry-air purge during assembly to protect MoS2 from humidity

Validated up to 4300 cycles at >= 1.0 MN compression in cryo-vacuum,
friction coefficient mu < 0.08 in the "fixed pad" design.  The
"fixed pad" design was preferred over "loose pad" because plastic
deformation of the spherical pad is welcome -- it self-adjusts
to compensate for assembly errors.

## Cryoplant

7 kW equivalent cooling capacity at 4.5 K, scalable to 10 kW for
several hours by expanding LHe from a 10000 L tank.

Cooling loops:
- Loop 1: CICC supply (supercritical He, 4 K, 5 bar)
- Circuit 3: CSR + coil housings (4 K, 3.5 bar; 80 parallel
  sub-loops with throttle valves)
- Circuit 5: Thermal shield + HTSC current leads (50-60 K, 14-18 bar
  inlet, throttled to 2-3 bar)

Heat exchangers downstream of circulators reduce supply T to 3.4 K
minimum via three sub-atmospheric He baths (3.3 / 3.8 / 4.4 K).

## HTSC current leads

14 leads, each 20 kA, bridging cryostat (4 K) to room-temperature
busbars.  HTSC (BSCCO-2223 or YBCO) insert reduces conduction heat
load between 65 K and 4 K compared to all-Cu leads.

Developed in collaboration with Karlsruhe Institute of Technology
(KIT), with the special "bottom-up" geometry (cold end on top,
warm end at bottom) because the W7-X power supplies sit directly
below the cryostat.

(Fietz et al. "High temperature superconductor current leads for
Wendelstein 7-X and JT-60SA" IEEE Trans. Appl. Supercond. 19
(2009) 2202)

## ECRH heating system

10 MW at 140 GHz continuous wave (resonant with the 2.5 T design
field, second harmonic X-mode = X2).  Ten 1-MW gyrotrons.

For half-field operation (1.25 T -> 70 GHz second harmonic = 0.7 T
fundamental, or 1.25 T fundamental = 35 GHz which is too low),
the gyrotrons are tunable to 103.6 GHz output, fundamental O-mode
at lower B.

Quasi-optical transmission via single-beam-waveguide (SBWG) + multi-
beam-waveguide (MBWG): 5 single-beam mirrors per gyrotron for
polarization conditioning, then combined into a 7-channel MBWG
transmitting 40 m to in-vessel launchers.  Total losses ~3% over
40 m -- one of the cleanest high-power quasi-optical systems built.

(Erckmann et al. Fusion Sci. Technol. 52 (2007) 291)

## Divertor

W7-X uses an **island divertor**: the plasma edge naturally has an
m=5 island chain (matches the 5-fold symmetry).  10 divertor units
(5 top, 5 bottom), each 12 target modules.

- High heat flux targets: CuCrZr cooling + CFC tiles, bonded by
  Plansee SE technology; survive 10 MW/m^2 + 10000 full-power cycles
- Baffles: actively cooled CuCrZr + clamped graphite tiles, 0.5 MW/m^2
- 10 cryopumps in divertor chamber (based on ASDEX Upgrade design)
- 10 in-vessel control coils behind baffle for island position tuning

## Construction timeline

- 1996: Project authorization
- 1998-2014: Component manufacture (15+ years)
- 2005-2014: Assembly in Greifswald (assembly process used 5
  sequential rigs)
- Total assembly effort: ~1 000 000 man-hours up to March 2014
  (Bosch et al. Nucl. Fusion 57 (2017) 116015), dominated by
  tolerance management for coil positioning
- 2015-12-10: First helium plasma
- 2018: First fully-cooled divertor operation
- 2021-present: Long-pulse operation toward 30 min steady-state

## Lessons for future stellarators

The W7-X cost growth (~1.06 B Euro per Wikipedia) and 18-year
construction time motivated the **simpler coils** research program
at PPPL / UW-Madison / etc.: redesign Stage-1 to allow easier
Stage-2 coils, accepting some physics performance compromise.
This is one of the founding motivations for the FOCUS code
(`stellarator_coil_design` topic).

Permanent magnets (NdFeB) as part of the shaping field have been
proposed as a way to entirely avoid the most twisted modular coils
(Zhu et al. arXiv 2020).  Toroidal field is still produced by
simple (often circular) TF coils -- the Ampere's-Law argument
shows you cannot generate the full B with magnets alone.

## References

- Bosch, H.-S. et al. "Construction of Wendelstein 7-X -- Engineering
  a Steady-State Stellarator" IEEE Trans. Plasma Sci. 38 (2010) 265
  W:/.../11_三菱_モジュラー型/Construction of Wendelstein 7-X.pdf
- Heinemann, B. et al. "Design of Narrow Support Elements for Non
  Planar Coils of Wendelstein 7-X" Proc. 21st IEEE SOFE (2006)
  W:/.../11_三菱_モジュラー型/Design of Narrow Support Elements...pdf
- Bosch, H.-S. et al. "Final integration, commissioning and start
  of the Wendelstein 7-X stellarator operation" Nucl. Fusion 57
  (2017) 116015
- Alternative conceptual design of a magnet support structure for
  plasma fusion (W:/.../11_三菱_モジュラー型/Alternative.../...)
- First operational phase of the superconducting magnet system of
  W7-X (W:/.../11_三菱_モジュラー型/First operational phase...)
- Engineering Lessons Learned in the Assembly, Commissioning,
  Initial Operation and in the Further Upgrading of W7-X
  (W:/.../11_三菱_モジュラー型/Engineering Lessons Learned...)
"""


LHD_HELICAL = r"""
# Large Helical Device (LHD) -- NIFS Toki, Japan

LHD is the world's largest heliotron-class superconducting
stellarator, in operation since 1998-03-31 at the National Institute
for Fusion Science (NIFS), Toki, Gifu prefecture.  Built by a
Mitsubishi-led consortium under the Japanese national large-device
program announced in the late 1980s.

## Heliotron principle vs modular stellarator

The heliotron is a stellarator subclass distinguished by **continuous
helical poloidal-field windings** instead of discrete modular coils.
The two-helical-pole (l=2) heliotron, of which LHD is the largest
example, generates both the rotational transform and a significant
fraction of the toroidal field from the same two helical conductors
wound in opposite senses around the torus.

Advantages relative to modular stellarator (e.g. W7-X):
- No coil-coil tolerance problem: the helical winding is one
  continuous piece, eliminating the W7-X NSE / LSE complication
- No assembly tolerance bottleneck (W7-X spent ~10^6 man-hours
  on coil positioning)
- Larger plasma access area (no modular coil "fingers" in the way)

Disadvantages:
- Cannot be re-optimized for a different rotational transform
  profile by changing coil currents alone -- the geometry is fixed
  at fabrication
- Larger neoclassical transport than QA/QH/QI modular optimization
  (LHD is not fully optimized for low neoclassical)
- High build cost concentrated in one massive monolithic winding

## LHD parameters

| Parameter                  | Value                          |
|----------------------------|--------------------------------|
| Major radius R             | 3.9 m (adjustable 3.6-4.1 m)  |
| Average minor radius       | 0.5-0.65 m                     |
| Toroidal field B0          | 2.5 T (operational), 3.0 T (design) |
| Helical pole              | l = 2 (two helical windings)   |
| Toroidal period            | m = 10 (10-fold symmetry)      |
| Helical conductor          | NbTi forced-flow CICC          |
| Helical operating current  | ~13 kA                         |
| Helical winding mass       | ~270 t (total)                 |
| PF ring coil pairs         | 3 (IV / IS / OV: inner-vertical, inner-shaping, outer-vertical) |
| PF conductor               | NbTi forced-flow CICC          |
| Cryostat outer diameter    | ~13.5 m                        |
| Cryostat height            | ~9 m                           |
| Magnet stored energy       | ~0.9 GJ                        |
| Cryoplant capacity         | 5.65 kW @ 4.4 K                |
| First plasma               | 1998-03-31                     |
| Deuterium operation        | 2017-03 onward                 |

## Key operational milestones

- 1998-03-31: First plasma (hydrogen)
- 2005-2008: 1 hour x 1 MW steady-state operation
- 2007: Te = 10 keV achieved with NBI + ECH
- 2017: Deuterium operation begins (improved confinement, alpha
  particle physics)
- 2019-2020: ne = 1.2 x 10^21 m^-3 super-dense core mode
  (gas-puff fueled, ICRF + NBI heated) -- a NIFS world-record
  class result for steady-state high density

(Source: Iiyoshi et al. "Overview of the Large Helical Device
project" Nucl. Fusion 39 (1999) 1245; "LHDプラズマのモデリングと
理論解析" W:/.../10_三菱_ヘリカル型/)

## Mitsubishi heliotron lineage

The Japanese heliotron program traces back to Heliotron-D (1969),
Heliotron-E (1980, Kyoto University), and CHS (Compact Helical
System, 1989, NIFS).  Mitsubishi Heavy Industries and Mitsubishi
Atomic Power Industries (now MHI Nuclear Engineering, MNE) were
prime contractors for the large devices in this lineage:

- Heliotron-E: 0.4 m minor radius, 2 T field, conventional copper
- CHS: 0.2 m minor radius, 0.9 T field, conventional copper
- LHD: 0.5+ m minor radius, 2.5+ T field, superconducting

The same Mitsubishi engineering team built the NIFS Compact Helical
System (CHS) before LHD, and contributed Cariddi-style transient
eddy current methods (Kameari "Transient Eddy Current Analysis on
Thin Conductors with Arbitrary Connections and Shapes" J. Comput.
Phys. 42 (1981) 124 -- written while Kameari was at Mitsubishi
Atomic Power Industries, Omiya).  See `transient_eddy` topic for
the Kameari method.

## References (Mitsubishi heliotron lineage folder)

W:/.../10_三菱_ヘリカル型/:
- ヘリカル型核融合炉設計 (Helical fusion reactor design)
- LHDプラズマのモデリングと理論解析 (LHD modelling)
- 核融合科学研究所（仮称）設立と大型ヘリカル装置計画 (NIFS
  establishment and LHD plan, vintage 1980s document)
- 大型ヘリカル装置計画から (From the LHD project)
- 核融合における技術革新(1) -- 閉じ込め装置本体 (Innovation in
  fusion confinement devices)
- 核融合科学研究所 setup and LHD plan
- 企業による核融合研究の最近の動向 (Recent trends in private-
  sector fusion R&D in Japan)
- 平成14年度 研究・企画情報センター作業会報告 (NIFS 2002
  work-group report on fusion DEMO)
- コイル軸のシフトが大きいヘリカルヘリアック装置 (Helical-Heliac
  with large coil axis shift)
- ヘリオトロンプラズマにおける抵抗性交換型不安定性と磁気島との
  相互作用に関する研究 (Resistive interchange + magnetic island
  interactions in heliotron plasmas)

## Cross-references

- `stellarator` -- broader stellarator taxonomy and Stage-1/2
  optimization
- `w7x_modular` -- modular alternative (contrasting design choices)
- `transient_eddy` -- Kameari thin-conductor FEM (Mitsubishi Atomic
  Power Industries provenance)
- `radia_mcp.matrix_solvers` -- iterative solvers for large CICC
  network problems
"""


SC_CABLE = r"""
# Superconducting cable for fusion magnets

All large-scale fusion magnets in operation or under construction
use the **cable-in-conduit conductor (CICC)** architecture: a
multi-stage twisted cable of superconducting strands flowing inside
a structural conduit, with supercritical helium cooling through
the void fraction.

## Superconductor materials

### NbTi (low-temperature, "workhorse")

- T_c0 = 9.2 K (zero field), B_c2(4.2 K) ~ 14 T
- Practical operating limit: ~8 T at 4.5 K (with margin)
- Ductile, drawable: standard wire metallurgy works
- Cheap, mature, used in MRI/NMR/LHC dipoles since 1980s
- ITER: PF coils, correction coils
- W7-X: all 70 coils (NPC + PC)
- LHD: helical and PF coils
- JT-60SA: TF, CS, EF coils (lower peak field than ITER)

### Nb3Sn (intermediate, "high-field workhorse")

- T_c0 = 18 K, B_c2(4.2 K) ~ 23 T
- Practical operating limit: ~13 T at 4.5 K
- A15 intermetallic, **brittle and strain-sensitive**:
  Jc(epsilon) drops by ~50% at epsilon = +/- 0.5%
- "Wind and react": wind ductile Nb + Sn precursor, then react at
  650 C to form Nb3Sn -- cannot bend the finished conductor
- Internal-tin or bronze-route filaments
- ITER: TF, CS coils (where peak field > 8 T forces Nb3Sn)
- LHC HL-LHC inner triplet: Nb3Sn quadrupoles (11-12 T)

### Nb3Al (research)

- T_c0 = 18.5 K, B_c2(4.2 K) ~ 26 T
- Less strain-sensitive than Nb3Sn (better for racetrack and
  D-shape coils)
- Manufacturing not yet at industrial scale -- JAEA / NIMS
  research program (mass cabling demonstrated 2010s)
- Not yet used in operational fusion magnet

### HTSC (REBCO, BSCCO)

- T_c ~ 90-110 K
- Practical use in fusion magnets so far: **current leads only**
  (W7-X 20 kA leads, JT-60SA 25 kA leads -- Fietz et al.
  IEEE Trans. Appl. Supercond. 19 (2009) 2202)
- HTS reduces conduction heat load from room-T busbar into the
  4 K helium bath by ~factor 3 vs all-Cu lead
- Future fusion magnets (SPARC, ARC, ST40 follow-ons) target
  REBCO tape stacks at 20 K, 20 T -- compact tokamak revolution
  enabled by commercial REBCO availability (2020-)

## CICC architecture

Typical ITER-style Nb3Sn TF CICC:
1. Single Nb3Sn strand: ~0.8 mm diameter, 1000-3000 filaments of
   Nb3Sn embedded in Cu matrix (Cu:non-Cu ratio ~1.4)
2. Sub-cable: 3 strands twisted (twist pitch ~50 mm)
3. Petal: 3 sub-cables x 4 (= 36 strands) twisted
4. Final cable: 6 petals (= 1080 strands total) around a central
   spiral He cooling channel
5. Wrapped in Inconel or 316L SS jacket forming a square or round
   conduit
6. Void fraction ~30% for He flow

W7-X NPC CICC (NbTi):
- NbTi strand ~0.6 mm dia
- Final cable in Al-alloy conduit (lighter than steel)
- 108 windings per NPC

LHD helical CICC (NbTi):
- Conventional rope-lay NbTi in 316L SS conduit
- Pool-boiling He bath cooling for the lower-current-density PF
  coils; forced-flow CICC for the high-density helical winding

## Joint resistance

CICC segments must be joined every 100-200 m (limited by cabling
machine capacity, not physics).  Joint resistance budget:
- ITER target: < 5 nOhm per joint
- W7-X achieved: < 1 nOhm
- Heat dissipation: Q = I^2 R = (68 kA)^2 x 1 nOhm = 4.6 W per
  joint -- manageable for ~hundreds of joints with He flow

Joints are typically **shaking-hands soldered** lap joints, with
strands of both segments interleaved and Sn-Ag soldered.  Modern
joints (post-2010) use compact "stainless-soldered" overlap with
intermediate Cu sole-plate.

## Quench protection

If a small region of the cable transitions to normal state
(quench), Joule heating melts the conductor unless current is
dumped fast.  Quench protection:
- Voltage-tap-based quench detection (typical V_th ~100 mV)
- Quench detection triggers fast-discharge circuit
- Stored magnetic energy dumped into external resistor bank
  (typical time constant tau = L/R ~ 5-15 s for ITER, ~30 s for
  W7-X)
- Internal temperature rise during dump: target T_max < 150 K
  (limited by jacket thermal capacity)
- "Hot-spot temperature" Tmax = the analytical quench limit, set
  by integral J^2 dt = MIITs (Mega-Ampere^2 - milliseconds)

## Cross-references

- `iter_coil_system` -- where CICC parameters come from for ITER
- `w7x_modular` -- W7-X NbTi CICC and HTSC current leads
- `tokamak` -- broader context
- `radia_mcp.peec` -- when modeling the AC loss of CICC at ramp
  (filament-level Biot-Savart for self / mutual inductance)
"""


STELLARATOR_COIL_DESIGN = r"""
# Stellarator coil design: math and codes

The stellarator coil design problem is, mathematically:

Given a desired plasma boundary on which the magnetic field
B . n = 0 (so plasma equilibrium is supported), find a system of
external coils (filaments) that reproduce this normal-component-
zero condition to acceptable accuracy.

This is an **ill-posed inverse problem**: a smooth target B on the
plasma surface can in principle be reproduced by infinitely many
coil systems, but most are not buildable.  Regularization toward
"smooth, sparse, low-curvature" coil sets is what makes the
problem tractable.

## Mathematical formulation

Plasma boundary parameterized in Fourier coordinates:
  r + i z = sum_{m,n} b_{mn} exp(-i m u + i n v) * exp(i u)

Coils live on a "winding surface" or "control surface" S enclosing
the plasma:
  r + i z = sum_{m,n} c_{mn} exp(-i m u + i n v) * exp(i u)

The external B from a current distribution kappa on S follows
the surface-Biot-Savart law:

  B(r) = (mu_0 / 4 pi) curl integral_S kappa(u,v) cross N_c / R dS

with N_c the normal of S and R = |r - r'|.  The level curves of
kappa become coil filaments.

Expand kappa as Fourier series:
  kappa(u, v) = kappa_0 v + sum_{m,n} kappa_{mn} sin(m u + n v)

with kappa_0 the (given) net poloidal current and kappa_{mn} the
unknown harmonics.  Imposing B . N_b = 0 on the plasma boundary
gives a linear system in kappa_{mn} -- but the matrix is poorly
conditioned because high (m, n) coil harmonics decay exponentially
between S and the plasma surface, so the inverse amplifies noise.

Regularization: restrict to a small number of low-order coefficients,
both in the plasma boundary description (b_{mn}) and in kappa.
This is the "limited Fourier" approach of Alexander & Garabedian
(PNAS 2007, W:/.../research/Choice_of_coils_fusion.pdf).

## Design code lineage

| Code         | Year | Approach                                    |
|--------------|------|---------------------------------------------|
| **NESCOIL**  | 1987 | Merkel: surface current potential, fixed winding surface |
| **NESVD**    | 2003 | SVD regularization of NESCOIL                |
| **REGCOIL**  | 2015 | Tikhonov regularization, current sheet      |
| **ONSET**    | 1990s| Nonlinear optimization on top of NESCOIL    |
| **COILOPT**  | 1990s| Filament-level optimization, fixed surface  |
| **COILOPT++**| 2014 | C++ rewrite, smoother optimization          |
| **FOCUS**    | 2018 | Free 3D filament representation, analytical 1st/2nd derivatives, gradient-based + Newton |

FOCUS (Zhu, Hudson, Gates et al.) is the modern reference:
- Filaments described as 3D Fourier curves (no fixed winding surface)
- Analytical gradients and Hessians of the objective with respect
  to filament Fourier coefficients
- Faster than finite-difference, by 10x+ at full coil resolution
- Multiple optimization algorithms (gradient descent, NCG, LM,
  modified Newton) selectable
- 36 users worldwide as of 2020 (PPPL distribution list)

(Source: Zhu et al. "Towards simpler coils for optimized
stellarators" slide deck, W:/.../research/Simpler_Coils_Stellarators.pdf;
Zhu et al. Nucl. Fusion 58 (2018) 016008.)

## Objectives in FOCUS / REGCOIL

Typical multi-objective cost:
- Normal field error |B . n|^2 on plasma boundary
- Toroidal magnetic flux mismatch (kappa_0 sets this)
- Resonant perturbation amplitudes B_mn(at resonant surface)
- Quasi-symmetry preserving (Boozer B-spectrum)
- Coil length penalty (cost / fabricability)
- Coil curvature penalty (radius-of-curvature lower bound)
- Coil-surface separation (assembly access)
- Coil-coil separation (diagnostic port space)
- (Stage 1 + 2 coupling: full free-boundary VMEC + STELLOPT for
  joint optimization)

## Hessian sensitivity analysis

After convergence to a coil set, compute the Hessian H of the
objective with respect to coil-position parameters.  Eigenvectors
of H identify the most error-prone deformation modes; their
eigenvalues give the cost penalty per unit deformation.

Application (Zhu et al. Nucl. Fusion 59 (2019) 126007):
- Apply principal eigenvector deformation to a CFQS coil set at
  amplitude xi = 0.01
- Observe whether magnetic islands at the n/m = 4/11 rational
  surface enlarge or shrink
- Result: "red parts" (high-sensitivity coil regions) need tight
  tolerance; "blue parts" can be relaxed -- replaces the W7-X
  practice of uniform 1/1000 R tolerance (Bosch et al. NF 2017)

This is the **adaptive tolerance** approach: tighten where the
physics is sensitive, loosen where it is not, saving fabrication
cost and assembly time.

## Permanent-magnet stellarators

Once Stage-1 says "I need this 3D shaping field", you can ask:
can permanent magnets generate part of it, sparing some modular
coils?

Answer (Helander et al. NF 2020, Zhu et al. arXiv 2020):
- Toroidal field must come from coils (Ampere's-law argument:
  loop integral of B inside the plasma is nonzero, but
  magnetization integral is zero in the plasma region)
- Shaping field (the 3D part) **can** come from NdFeB permanent
  magnets arranged in Halbach-like patterns
- NdFeB B_r up to 1.55 T (Matsuura JMM 2006)
- NdFeB H_ci up to 5 T when cooled (so they can sit inside a few-T
  background without demagnetizing)
- Kumada IEEE 2004: 5.16 T peak field demonstrated in a permanent-
  magnet quadrupole (laboratory dipole magnet for accelerators)
- Stellarator application: PM4Stell project (PPPL + collaborators)

This is the most active research direction in modular stellarator
simplification, and the closest tie-in for `radia_mcp` users:
NdFeB magnetization can be modeled as fixed M = Br / mu_0 in Radia
without solving anything (see `radia_mcp.electromagnet`).

## radia's stream-function solver vs these codes (2026-06-21)

radia ships the SAME object as NESCOIL/REGCOIL/FOCUS -- the winding-surface
CURRENT POTENTIAL psi (surface current K = n x grad psi) -- via
`radia.stream_function` + `calc_streamfunction.py` (knowledge:
`radia_mcp.streamfunction`).  Honest positioning:

**DESIGN -- AT PARITY.**  docs/stream_function/fusion.md and the docs-local
demo_regcoil_fusion.py (+_advanced) helpers reproduce the REGCOIL workflow:
current potential, the Tikhonov
L-curve (B.n residual vs peak |grad psi|), the net-current secular term
Psi = psi + (G/2pi) zeta + (I/2pi) theta (the 2 first-cohomology generators,
b1 = 2 - chi, gmsh-free), VMEC-shaped boundaries (li383 / NCSX-like, B.n ~4e-8),
coil Lorentz force/stress, and the FOCUS standoff + winding-SHAPE levers
(kappa=2 plasma: conformal winding cuts coil complexity ~34%).  B.n ~2e-9 on
producible targets; no fusion-specific solver code (design-matrix rows are just
the plasma-normal Biot-Savart of the winding surface).

**SCALE / METHOD -- ACA + ridge-TSVD.**  radia's design solver is the dense
Biot-Savart design matrix COMPRESSED by ACA (an H-matrix) and inverted by a
ridge (Tikhonov) TRUNCATED-SVD pseudo-inverse (`radia.stream_function.aca_tsvd`;
knowledge topic `aca_tsvd`).  The ridge IS REGCOIL's lambda (same regularisation
role), but the ACA compression + the FE-direct surface make it a SCALABLE
regularised least-squares on an ARBITRARY meshed winding surface -- not REGCOIL's
dense Fourier least-squares on a parameterised torus.  Goldens (a)/(b):
test_regcoil_parity_deliverable_golden.py, test_regcoil_iron_differentiator_golden.py.

**DELIVERABLE -- radia is a strict SUPERSET (design-to-manufacture, not
design-only).**  NESCOIL/REGCOIL/FOCUS stop at psi (a field and its contours).
radia continues: contours -> SINGLE-STROKE wire (grad-psi winding orientation,
so even l>=2 saddle shims chain without the common current cancelling) ->
SHEET-METAL distort (bend the discrete wire against its own Biot-Savart) ->
STEP CAD (OCC WriteStep) -> PEEC.  The PEEC step is a full circuit-extraction
SOLVER (`radia.peec_*`: L, R, C, M + SPICE netlist, MMM coupling), not just an
inductance number.  So from ONE run radia emits a windable conductor + CAD + a
circuit model that the design codes do not.

**PHYSICS -- radia is a SUPERSET (iron).**  NESCOIL/REGCOIL/FOCUS are ALL
free-space (vacuum Biot-Savart) and cannot handle magnetic material.  radia's
material-aware kernel (Kelvin-FEM DtN transfer M = M_free + M_react, the
`--iron-vol` path) designs coils WITH iron yoke/shield/core (free-space design
misses ~77% with iron; material-aware matches ~1e-4) -- opening the domains
REGCOIL structurally cannot enter (and, for fusion, iron flux returns /
structures).

**HONEST NUANCE (do not over-claim).**  The single-stroke win is decisive for
SINGLE-CONDUCTOR coils (MRI gradient/shim, induction heating -- same current-
potential math).  Stellarator MODULAR coils are intentionally SEPARATE coils
(equal-Delta-Psi contours), so "one wire" is not their manufacturing step;
there radia still adds STEP CAD + PEEC L beyond REGCOIL's psi, and its distort
is analogous to FOCUS filament-shape optimization (an addition vs NESCOIL/
REGCOIL, not unique vs FOCUS).  Precise claim: "REGCOIL is a special case
(vacuum, toroidal, design-only) of radia's SF; design at parity, deliverable +
iron a superset; a complete win for single-conductor coils."

**EARN THE CLAIM BY MEASUREMENT (Repository-First, not paper-chasing):**
  (a) vacuum-parity golden -- same target: B.n parity with REGCOIL AND the same
      run emits a *.step + PEEC L (locks "design equal, deliverable beyond").
  (b) iron differentiator on a standard problem -- REGCOIL misses, radia matches
      FEM/measurement (productionise the demo_ee/ff bridge).
  (c) manufacturing end-to-end -- target -> psi -> single-stroke -> distort ->
      STEP -> PEEC, wound-wire field verified.
Open gaps (honest): no SIMSOPT/STELLOPT Stage-1+2 integration; performance
head-to-head at community mode counts unmeasured; full winding-surface
optimization loop not built.

## Cross-references

- `stellarator` -- overview, Stage-1 / Stage-2 framework
- `error_field` -- linking Hessian sensitivity to error-field
  tolerance budget
- `w7x_modular` -- the engineering reality after Stage-2
- `radia_mcp.electromagnet` -- PM-stellarator coil modeling
  (Radia ObjHexahedron with fixed [Mx, My, Mz])
- `radia_mcp.streamfunction` -- the SF engine that IS the current
  potential (topics: fusion, material_aware, low_turn, single_stroke)
"""


ERROR_FIELD = r"""
# Error field control in tokamaks and stellarators

Plasma equilibrium and stability are remarkably sensitive to small
non-axisymmetric perturbations of the magnetic field.  In tokamaks,
the threshold for **mode locking** (a tearing mode locking to a
quasi-static perturbation, often leading to disruption) is at the
~10^-4 relative-amplitude level for the m/n = 2/1 component.
Stellarators are similarly sensitive to high-m resonant harmonics
that can destroy outer flux surfaces.

## Sources of error field

| Source                                         | Magnitude (ITER scale)           |
|------------------------------------------------|-----------------------------------|
| Coil manufacturing tolerances (winding pack)   | ~0.1-0.5 mT (Monte Carlo)         |
| Coil assembly tolerances (rigid displacement)  | ~0.1 mT                           |
| Magnetic materials (NB shield, TBM ferritic)   | ~0.05 mT                          |
| Iron yokes around NB injector (compensation)   | ~0.05 mT (with compensation coils)|
| Bus-bar stray fields                           | bifilar routing reduces to ~0.01 mT|
| Geomagnetism                                   | 0.05 mT (constant background)     |
| Cryostat / building steel                      | < 0.01 mT                         |

Total expected ITER error field (TMEI): O(0.5 mT) at q=2 surface
in absence of correction.  Target after CC correction:
< 0.05 mT (5 x 10^-5 * B_t = 5 x 10^-5 * 5.3 T).

## Tolerance metric: TMEI (Three Mode Error Index)

The dominant tearing-mode-locking modes are m/n = 1/1, 2/1, 3/1.
The relevant scalar metric is

  BTMEI = sqrt(0.2 * |B11|^2 + |B21|^2 + 0.8 * |B31|^2)

where Bmn is the normal component of the m, n harmonic on the q=2
magnetic surface (Hender et al. NF 47 (2007) S128).

Threshold for mode locking and induced disruption: BTMEI / B_t
~ 1 x 10^-4 in tokamaks similar to JET (and used as the working
estimate for JT-60SA and ITER, where direct experimental data is
not yet available).

## Monte-Carlo error-field prediction

Manufacturing tolerances are specified as bounds (e.g. +/- 2 mm
position tolerance at high-field side of TFC, +/- 4 mm at LFS).
Translating tolerances to mode amplitudes requires Monte-Carlo
sampling:

1. For each coil, sample a random deformation within tolerance
   (typical: 8-node spline function for winding-pack shape errors,
   6-degree-of-freedom rigid body displacement for assembly errors)
2. Compute the perturbation B-field from the deformed coil current
   filaments using Biot-Savart
3. Decompose B at the q=2 surface into (m, n) Fourier harmonics
4. Repeat O(1000-250 000) times to build histograms of TMEI

JT-60SA EFCC design used 5000 trials for the TMEI distribution and
250 000 trials for the correction-current optimization (Matsunaga
et al. Fusion Eng. Des. 98-99 (2015) 1113).

Outcome on JT-60SA (scenario-LN, Ip = 5.5 MA, B_t R = 6.66 T-m):
- TMEI from TF coils alone: 0.28 mT (95%-cumulative)
- TMEI from PF coils alone: 0.24 mT
- TMEI combined: 0.34 mT
- Target after EFCC: < 0.1 mT (margin below 0.2 mT mode-lock
  threshold)

EFCC design solution:
- 18 in-vessel coils in 3 toroidal rows of 6 (upper, middle, lower)
- Each coil: 35 turns of copper conductor (water-cooled, not
  superconducting -- they sit inside the cryostat at the first
  wall back side)
- Current: 30 kA-turn upper/lower, 45 kA-turn middle
- Optimization: least-squares minimization of m = 1, 2, 3 normal-
  component-Bn at q = 2 surface, with anti-series connections
  between toroidally opposite coils

## Error field for stellarators

Stellarators have higher inherent error-field sensitivity than
tokamaks because their magnetic field is **entirely** produced by
the coil system (no plasma current to wash out small perturbations).

W7-X assembly tolerance: 1 mm / 5 m for coil-coil positioning,
3 mm / 5 m for global coil set position.  These were achieved
through the NSE / LSE structure described in `w7x_modular`.

Kisslinger & Andreeva (Fusion Eng. Des. 2005) computed the
Fourier components in B_mn (in units of 10^-4 / B_00) generated
by 1 cm rigid shifts and 1 deg inclination angles applied to a
single W7-X module.  Useful as a sensitivity matrix for
prioritization.

CFQS (Chengdu) used FOCUS Hessian analysis to identify which
coil segments need tight tolerance and which can be relaxed --
"adaptive tolerance" instead of uniform tolerance.  See
`stellarator_coil_design` topic.

## Mitigation hierarchy

1. **Passive**: design coils that are inherently tolerant
   (use sensitivity Hessian to pick the right Stage-2 minimum)
2. **Assembly QC**: laser-tracking measurements of as-built coil
   positions; up-front error-field prediction from as-built data
3. **Compensation magnetic shields** (ferritic) on the LFS (TF
   ripple mitigation) and around NB injectors
4. **Active error-field correction coils** (CC in ITER, EFCC in
   JT-60SA) -- these can also serve as RMP coils for ELM control,
   see `rmp_elm_control`
5. **Real-time mode-locking detection** + active feedback

## References

- Hender, T.C. et al. "Chapter 3: MHD stability, operational limits
  and disruptions" Nucl. Fusion 47 (2007) S128 -- TMEI metric
- Matsunaga, G. et al. "In-vessel coils for magnetic error field
  correction in JT-60SA" Fusion Eng. Des. 98-99 (2015) 1113-1117
  W:/.../research/InVesselCoils_JT60SA_ErrorField.pdf
- Knaster, J. et al. "Assessment of the JT-60SA error field and
  correction" Fusion Eng. Des. 86 (2011) 1053
- Park, J.-K. et al. "Computation of three-dimensional tokamak and
  spherical torus equilibria" Nucl. Fusion 48 (2008) 045006
- Zhu, C. et al. "Hessian matrix approach for determining error
  field sensitivity to coil deviations" Plasma Phys. Control.
  Fusion 60 (2018) 054016
- 原型炉_誤差磁場_補正コイル_日立 (Hitachi DEMO error field study)
  W:/.../research/原型炉_誤差磁場_補正コイル_日立.pdf

## Cross-references

- `iter_coil_system` -- where the CC coils sit in ITER
- `tokamak` -- broader context, including disruption physics that
  is triggered by mode locking
- `stellarator_coil_design` -- Hessian sensitivity in coil
  optimization
- `rmp_elm_control` -- the EFCC also doubles as the ELM-control
  RMP coil
"""


TRANSIENT_EDDY = r"""
# Transient eddy current analysis in fusion devices

Fusion magnets and passive structures (vacuum vessel, cryostat,
blanket, divertor cassettes, casings) experience large transient
magnetic fields during plasma scenarios (current ramp-up, scenario
evolution) and especially during **disruptions** (millisecond
plasma current collapse).  Resulting eddy currents drive mechanical
forces and Joule heating that must be quantified.

## Physics regime

- Time scales: ~ms (disruption CQ) to ~minutes (slow scenario ramp)
- Length scales: ITER vacuum vessel is ~20 m diameter
- Penetration: typical SS wall thickness 30-60 mm, skin depth at
  100 Hz (disruption fundamental) ~30 mm -> partial penetration
- Quasistatic regime: wavelength c / f >> device size for all
  relevant frequencies (Laplace kernel suffices, no full-wave)

## Numerical methods family

### Volume integral formulations (Cariddi family)

The dominant method for fusion-scale 3D passive structures
(Rubinacci, Villone et al., CREATE consortium / Univ. Naples / Univ.
Cassino).

T-formulation (electric vector potential):

  J = curl T,  with T = sum_i I_i N_i

Edge-element shape functions N_i for T with tree-cotree gauge
ensure J is solenoidal (J . n = 0 on the boundary minus the
electrodes).  Galerkin discretization gives:

  L dI/dt + R I + F^T phi = -dVs/dt
  F I = I_h
  I(0) = I_0

where
- L is the (dense) inductance matrix (Biot-Savart kernel)
- R is the (sparse) resistance matrix (local Ohm's law)
- F couples electrode current constraints (sparse, surface only)
- phi is the unknown electrode potential
- Vs is the externally-applied source potential

Implicit time stepping gives a dense system to solve each step,
prohibitive for N > ~10000 DOF.

(Source: Cau, Chiariello, Rubinacci, Scalera, Tamburrino, Ventre,
Villone "A Fast Matrix Compression Method for Large Scale Numerical
Modelling of Rotationally Symmetric 3D Passive Structures in
Fusion Devices" Energies 15 (2022) 3214.
W:/.../research/FastMatrixCompression_Fusion_3D.pdf)

### Thin-conductor current-function formulation (Kameari 1981)

When the conductor is thin compared to the skin depth, treat it as
a 2D current sheet with a scalar **current function** V:

  J = grad V cross n

Energy integrals (magnetic energy, mutual energy with external
field, Joule loss) are integrated over the conductor surface
discretized into triangular FE.

  M dV/dt + R V = E_ext(t)

with M the dense surface inductance matrix and R the sparse
resistance matrix.  Same eigen-analysis structure as Cariddi but
on a 2D surface mesh.

This is the **Kameari method** (Mitsubishi Atomic Power Industries,
Omiya, 1980; published J. Comput. Phys. 42 (1981) 124).  It was
adopted as the standard for INTOR-J primary shield analysis and
has been extended by Bossavit (1980), Blum et al. (1985), and
others.  It is still the right method for the first-pass
eddy-current analysis of thin shells in modern devices.

(Source: Kameari "Transient Eddy Current Analysis on Thin
Conductors with Arbitrary Connections and Shapes" J. Comput. Phys.
42 (1981) 124-140.
W:/.../research/Transient_Eddy_Thin_Conductors.pdf)

### Boundary conditions
- Conductors are insulated from each other except at junctions
- At a connected junction line C with multiple subconductors,
  sum of currents flowing in = 0 (Kirchhoff's first law,
  generalized to current functions)
- Symmetry: V(T x) = +/- V(x) depending on the parity of the
  external field under T (rotation / inversion / reflection)

## H-matrix / fast methods for ITER-scale

For a fully-3D ITER vacuum vessel + blanket back-plate model with
N ~ 10^5-10^6 DOF, the dense L matrix is ~10^12 entries -- 8 TB at
double precision, infeasible.

Three families of acceleration:

1. **FFT-based** (toeplitz exploitation): works for periodic grids,
   poor fit for fusion-vessel arbitrary geometry
2. **Fast multipole methods (FMM)**: O(N) but constants are large
   for the surface integrals; competitive but not dominant
3. **QR-recursive H-matrix** (CREATE consortium): O(N log N)
   memory and O(N log^2 N) MatVec, dominant in fusion practice

QR-recursive algorithm (Cau et al. 2022):
1. Partition mesh elements into hierarchical boxes by spatial
   bisection
2. Classify box-box pairs as near-field (geometrically adjacent)
   or far-field (geometrically separated by a few box widths)
3. Near-field interactions: store exactly (sparse pattern)
4. Far-field interactions: compute QR factorization of the L_bc
   block with truncation tolerance -> low-rank representation
5. MatVec L * x = sum of near-field exact + sum of far-field
   low-rank products

Extensions in Cau 2022 paper:
- DOF-based partition (instead of element-based) to handle the
  electrode-boundary case cleanly
- Symmetric / periodic boundary condition embedded in the matrix F
  (force J at one sector edge = J at opposite sector edge)
- "Small box" remedy: merge under-populated leaf boxes to keep
  the recursive QR efficient
- Applied to ITER full 3D vacuum vessel + port plug structures

## Application: ITER VDE (Vertical Displacement Event)

Vertical displacement of the plasma column during a disruption
generates large axisymmetric eddy currents in the vacuum vessel,
poloidal in direction.  Net vertical force on VV:

- ITER design VDE: 150-200 t (1.5-2 MN) net vertical force
- Direction: opposes the plasma motion (Lenz)
- Duration: ~100-300 ms (CQ + VDE drift time)
- Structural constraint: VV + cryostat support must withstand
  this force without yielding

Halo current contribution: a fraction of plasma current flows in
the open field-line "halo" region, contacting the wall and
flowing in the SS structure.  Halo Toroidal Peaking Factor
(HTPF) ~2-3 -> local current ~3x toroidal average -> local force
hot spots.

## Tools available in `radia_mcp`

- `radia_mcp.fem.time_domain_axisym.time_domain_fem` -- axisymmetric
  transient FEM (good for VDE and centered disruptions where
  symmetry is preserved)
- `radia_mcp.bem.h_matrix` -- H-matrix BEM accelerator (analog of
  the Cau et al. QR-recursive method for surface integral formulations)
- `radia_mcp.matrix_solvers.preconditioners.amg` -- AMG
  preconditioner for the dense / quasi-dense system
- `radia_mcp.team_benchmark.eddy_current.problem_24` -- TEAM 24
  transient eddy current benchmark (analogue of disruption eddies
  in a simple geometry)

## References

- Kameari, A. "Transient Eddy Current Analysis on Thin Conductors
  with Arbitrary Connections and Shapes" J. Comput. Phys. 42
  (1981) 124-140
  W:/.../research/Transient_Eddy_Thin_Conductors.pdf
- Bossavit, A. "Computational Electromagnetism" Academic Press
  (1998) -- canonical text, chapter on the T-Omega formulation
- Albanese, R. and Rubinacci, G. "Integral formulation for 3D
  eddy-current computation using edge elements" IEE Proc. A 135
  (1988) 457 -- Cariddi foundational paper
- Cau, F. et al. "A Fast Matrix Compression Method for Large Scale
  Numerical Modelling of Rotationally Symmetric 3D Passive
  Structures in Fusion Devices" Energies 15 (2022) 3214
  W:/.../research/FastMatrixCompression_Fusion_3D.pdf
- Hammerschmid, D. et al. "Transient eddy current and force
  analysis for ITER" Fusion Eng. Des. 86 (2011) 1808 -- ITER
  VDE benchmark

## Cross-references

- `tokamak` -- disruption physics that drives the transient
- `iter_coil_system` -- the passive structures (VV, blanket,
  cryostat) that carry the eddy currents
- `radia_mcp.bem.h_matrix` -- the H-matrix toolkit
- `radia_mcp.fem.time_domain_axisym` -- axisymmetric transient FEM
- `radia_mcp.matrix_solvers` -- linear-system tooling
"""


RMP_ELM_CONTROL = r"""
# Resonant Magnetic Perturbation (RMP) for ELM control

Edge Localized Modes (ELMs) are quasi-periodic relaxation events of
the H-mode pedestal: pressure builds at the edge until peeling-
ballooning instability triggers, ejecting ~5-10% of the pedestal
energy in <1 ms.  For ITER and DEMO, even Type-I ELMs would deliver
unacceptable transient heat loads (~1 GW/m^2 peak) to the divertor
targets -- ELM suppression or mitigation is mandatory.

## Why RMPs work

A small non-axisymmetric perturbation imposed on the H-mode edge
(with toroidal mode number n = 1-4 typically) creates magnetic
islands at the edge resonant rational surfaces (q = m/n).  If
the islands overlap (Chirikov criterion), the edge magnetic
structure becomes **stochastic** -- field lines wander rather than
closing into flux surfaces.

Stochastic field lines enhance parallel transport, reducing the
pedestal pressure gradient below the peeling-ballooning stability
threshold.  ELMs are either:
- **Suppressed** (vanish entirely)
- **Mitigated** (smaller and more frequent, lower per-ELM
  energy deposition on divertor)

First demonstration: DIII-D 2003 (Evans et al. PRL 92 (2004) 235003,
Nucl. Fusion 45 (2005) 595).  Now standard tool on every major
tokamak.

## RMP coil hardware

RMPs are produced by in-vessel coils sitting behind the first wall:

| Device   | RMP coil layout                                |
|----------|-------------------------------------------------|
| DIII-D   | I-coils (12 inboard) + C-coils (6 outboard)    |
| ASDEX-U  | B-coils (8 upper + 8 lower)                    |
| KSTAR    | IVCC (3 toroidal rows x 4 coils each)          |
| JET      | EFCC (4 coils external to vessel)              |
| ITER     | ELM coils (27 in-vessel + 18 ex-vessel CC)     |
| JT-60SA  | EFCC (3 rows x 6 coils, same as error field)   |

The JT-60SA philosophy is dual-use: the same in-vessel coils that
correct error fields (`error_field` topic) also generate RMPs for
ELM control.  This is enabled by independently powering the 18
coils and choosing currents that produce either a low-n correction
field (m/n = 1/1, 2/1, 3/1) or a higher-n RMP (n = 3 typically).

## n = 3 ELM control on JT-60SA

(Source: Matsunaga et al. Fusion Eng. Des. 98-99 (2015) 1113-1117,
W:/.../research/InVesselCoils_JT60SA_ErrorField.pdf)

JT-60SA EFCC for ELM control:
- n = 3 toroidal mode (3 sets of co-phased coils around the torus)
- "Even parity": upper and lower in phase, middle out of phase
- EFCC current 5 kA: stochastic region at psi_n ~ 0.97 (very edge)
- EFCC current 10 kA: stochastic region extends to psi_n ~ 0.94
  (covers the pedestal region)
- Poincare plots show magnetic islands at rational surfaces; with
  increasing current, neighbouring islands overlap and the edge
  becomes fully stochastic

10 kA EFCC current is sufficient to produce a stochastic layer
that covers the typical pedestal extent.  In-coil current density
is well below the 30 kA / 45 kA maximum used for static error
field correction.

## RMP design trade-offs

- **n number**: higher n localizes the perturbation more strongly
  to the edge (good for ELM control, less core impact) but requires
  more coils to resolve the toroidal harmonic
- **Up-down phasing**: even parity (current direction same top
  and bottom) generates symmetric stochastization; odd parity is
  used for rotation control / plasma response studies
- **Coil location**: in-vessel coils (closest to plasma, lowest
  current required, but radiation-hardened design needed) vs
  ex-vessel coils (cheaper, but need much higher current to
  produce the same edge perturbation due to shielding by the wall)
- **Vacuum approximation vs plasma response**: simple ray-tracing
  / Poincare in vacuum field gives the geometric island structure;
  plasma response codes (MARS-F, M3D-C1, JOREK) compute amplification
  / screening by the plasma itself, often producing very different
  edge perturbation than vacuum

## Plasma response amplification (vacuum approximation limit)

In the simple vacuum-field calculation used by JT-60SA EFCC design:
- Field-line tracing applies static n = 3 RMP on top of 2D
  equilibrium
- 3D perturbation B is superposed linearly (vacuum approximation)
- Result: geometric island chain at n = 3 resonant surfaces
- Caveat: real plasma may amplify the perturbation by O(10) at
  some surfaces (resistive ballooning resonance) or screen it
  entirely at others (rotation-induced screening)

Full plasma-response analysis is out of scope for `radia_mcp` --
it requires linear MHD codes (MARS-F, IPEC) or nonlinear codes
(M3D-C1, JOREK) coupled to the vacuum perturbation field that
`radia_mcp` and similar tools provide as a boundary input.

## ELM control alternatives

- **Pellet pacing**: high-frequency pellet injection triggers
  small ELMs at controlled timing
- **Vertical kicks**: small vertical position oscillations trigger
  controlled ELMs
- **Lithium injection**: lower edge density gradient (DIII-D)
- **Negative triangularity**: avoid H-mode entirely (TCV, DIII-D)
- **Snowflake / super-X divertor**: increase wetted area to
  reduce per-ELM peak heat flux

The community consensus for ITER is **RMP suppression** as primary
strategy, with pellet pacing as backup.

## References

- Evans, T.E. et al. "Suppression of large edge-localized modes in
  high-confinement DIII-D plasmas with a stochastic magnetic
  boundary" Phys. Rev. Lett. 92 (2004) 235003
- Evans, T.E. et al. "Suppression of large edge localized modes
  with edge resonant magnetic fields in high confinement DIII-D
  plasmas" Nucl. Fusion 45 (2005) 595
- Shimada, M. et al. "Chapter 1: Overview and summary" Nucl. Fusion
  47 (2007) S1 -- ITER physics basis chapter on ELMs
- Matsunaga, G. et al. "In-vessel coils for magnetic error field
  correction in JT-60SA" Fusion Eng. Des. 98-99 (2015) 1113-1117
  W:/.../research/InVesselCoils_JT60SA_ErrorField.pdf -- dual-use
  EFCC for error correction + RMP ELM control
- Suzuki, Y. et al. (IAEA-CN-221 TH/P 7-37, 2014) -- RMP-induced
  field structure analysis for JT-60SA

## Cross-references

- `error_field` -- the same coil hardware does both jobs
- `tokamak` -- the disruption-and-ELM physics context
- `iter_coil_system` -- ITER ELM coil hardware
"""


def get_knowledge(topic: str = "overview") -> str:
    """Dispatch fusion topics.

    Topics:
        overview                  - Fusion confinement landscape (DEFAULT)
        tokamak                   - Tokamak (ITER, JT-60SA) coil system
        stellarator               - Stellarator overview + Stage-1/2
        iter_coil_system          - ITER TF/PF/CS/CC numbers, CICC
        w7x_modular               - W7-X engineering: NSE, LSE, cryoplant
        lhd_helical               - LHD heliotron, Mitsubishi lineage
        sc_cable                  - NbTi / Nb3Sn / HTSC cable physics
        stellarator_coil_design   - NESCOIL/REGCOIL/FOCUS, Hessian
        error_field               - Tolerance budgets, TMEI, EFCC design
        transient_eddy            - Kameari, Cariddi, H-matrix
        rmp_elm_control           - Resonant magnetic perturbation
        all                       - Everything
    """
    topic = topic.lower().strip()
    if topic in ("overview", "intro"):
        return OVERVIEW
    if topic in ("tokamak", "iter"):
        return TOKAMAK
    if topic in ("stellarator", "helical", "stellarator_overview"):
        return STELLARATOR
    if topic in ("iter_coil_system", "iter_coils", "iter_tf", "iter_cs", "iter_pf"):
        return ITER_COIL_SYSTEM
    if topic in ("w7x_modular", "w7x", "wendelstein", "wendelstein7x", "w7-x"):
        return W7X_MODULAR
    if topic in ("lhd_helical", "lhd", "heliotron", "mitsubishi"):
        return LHD_HELICAL
    if topic in ("sc_cable", "cicc", "nbti", "nb3sn", "htsc", "superconductor"):
        return SC_CABLE
    if topic in ("stellarator_coil_design", "coil_design", "focus", "nescoil", "regcoil"):
        return STELLARATOR_COIL_DESIGN
    if topic in ("error_field", "tmei", "efcc", "correction"):
        return ERROR_FIELD
    if topic in ("transient_eddy", "eddy", "disruption", "kameari", "cariddi"):
        return TRANSIENT_EDDY
    if topic in ("rmp_elm_control", "rmp", "elm", "elm_control"):
        return RMP_ELM_CONTROL
    if topic == "all":
        return "\n\n".join([
            OVERVIEW,
            TOKAMAK,
            STELLARATOR,
            ITER_COIL_SYSTEM,
            W7X_MODULAR,
            LHD_HELICAL,
            SC_CABLE,
            STELLARATOR_COIL_DESIGN,
            ERROR_FIELD,
            TRANSIENT_EDDY,
            RMP_ELM_CONTROL,
        ])
    return (
        f"Unknown topic '{topic}'. Available: "
        "overview, tokamak, stellarator, iter_coil_system, w7x_modular, "
        "lhd_helical, sc_cable, stellarator_coil_design, error_field, "
        "transient_eddy, rmp_elm_control, all."
    )
