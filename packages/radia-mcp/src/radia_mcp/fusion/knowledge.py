"""Fusion magnet knowledge."""

OVERVIEW = r"""
# Fusion reactor magnet systems

## Confinement principles

| Type | Plasma confinement | Magnet system |
|------|---------------------|---------------|
| **Tokamak** | Toroidal + induced poloidal | TF + PF + CS coils |
| **Stellarator** | 3D twisted toroidal | Modular / helical coils |
| **FRC** (Field-Reversed Config) | Compact toroid | Theta-pinch coil |
| **Inertial** | None (laser ablation) | none (drivers are lasers) |

## Operating regime

- **B field**: 5-13 T (LHC: 8.3 T, ITER TF: 11.8 T, SCMaglev EDS: 5 T)
- **Frequency**: DC for steady-state, slow ramp during startup
- **Eddy currents**: induced in passive structures during disruption
  (millisecond transient → forces / heating)

## Cross-references

- `radia_mcp.accelerator` — adjacent superconducting magnet research
- `radia_mcp.electromagnet` — DC magnet design (Hantila polarization)
- `radia_mcp.fem.potential_formulations.h_formulation` — HCurl H formulation for SC
"""


TOKAMAK = r"""
# Tokamak (ITER and similar)

## Coil system

| Coil | Function | Field |
|------|----------|-------|
| **TF** (Toroidal Field) | Main confining B | 5-12 T at plasma |
| **PF** (Poloidal Field) | Plasma position / shape control | various |
| **CS** (Central Solenoid) | Induce plasma current via flux swing | up to 13 T |
| **Correction coils** | Error field correction | ~1 T |
| **In-vessel coils** | Edge current control (ELM) | ~1 T |

## ITER specifics

ITER is the 35-nation fusion demonstrator (Cadarache, France):
- 18 TF coils (Nb3Sn superconducting)
- 6 PF coils
- 6 CS modules
- Plasma volume: 840 m³
- First plasma: target 2034 (delayed from earlier dates)

## Disruption forces

When plasma disrupts (loss of confinement, millisecond scale):
- Massive eddy currents induced in vacuum vessel + cryostat
- Forces up to 60 MN per coil
- Lab cross-link: transient eddy current FEM
  → `radia_mcp.fem.time_domain_axisym.time_domain_fem`

## References (W:/.../08_核融合/01_ITER/)

- 3 ITER papers (NP / IF coils, error field correction, etc.)
- 原型炉における超伝導コイル由来の誤差磁場及び補正コイル必要電流値の評価 (日立)

## Cross-references

- `radia_mcp.electromagnet` — Hantila for NL iron
- `radia_mcp.team_benchmark.eddy_current.problem_24` — transient nonlinear eddy
"""


STELLARATOR = r"""
# Stellarator — twisted toroidal confinement

Stellarators achieve confinement **without** a plasma current (unlike
tokamak). Magnetic field is fully 3D-twisted by stationary coils.

## Coil types

| Type | Examples |
|------|----------|
| **Helical** | LHD (Japan, Mitsubishi heliotron) |
| **Modular** | W7-X (Wendelstein 7-X, Germany) |
| **Heliotron** | Mitsubishi heliotron papers (lab connection) |

## LHD (Large Helical Device, NIFS Japan)

- Built by Mitsubishi (lab connection)
- Operating since 1998
- Continuous helical winding around torus

## W7-X (Wendelstein 7-X, IPP Greifswald)

- 50 non-planar modular coils (3D shape, no symmetry)
- Optimized for neoclassical confinement
- World's largest stellarator

## Mitsubishi research lineage (lab connection)

W:/.../99_アプリケーション/08_核融合/:
- 10_三菱_ヘリカル型_文献/ (16 files, 461 MB) — heliotron studies
- 11_三菱_モジュラー型_文献/ (20 files, 71 MB) — modular coil studies

These reflect Mitsubishi-NIFS collaboration on coil design + EM
analysis.

## Specific papers

- Choice of coils for a fusion reactor (lab core ref)
- Towards simpler coils for optimized stellarators
- A Fast Matrix Compression Method for Large Scale Numerical
  Modelling of Rotationally Symmetric 3D Passive Structures
- Transient Eddy Current Analysis on Thin Conductors with Arbitrary Connections
- 核融合研究50年の進展を振り返って (50-year review)

## Cross-references

- `radia_mcp.bem.h_matrix` — ACA / H-matrix for fusion-scale BEM
- `radia_mcp.matrix_solvers.preconditioners.amg` — AMG for huge systems
- `radia_mcp.accelerator` — adjacent SC magnet expertise
"""


def get_knowledge(topic: str = "overview") -> str:
    """Dispatch fusion topics.

    Topics:
        overview      - Fusion confinement landscape (DEFAULT)
        tokamak       - Tokamak (ITER) coil system
        stellarator   - Stellarator (LHD, W7-X, Mitsubishi heliotron lineage)
        all           - Everything
    """
    topic = topic.lower().strip()
    if topic in ("overview", "intro"):
        return OVERVIEW
    if topic in ("tokamak", "iter"):
        return TOKAMAK
    if topic in ("stellarator", "helical", "lhd", "w7x", "mitsubishi"):
        return STELLARATOR
    if topic == "all":
        return "\n\n".join([OVERVIEW, TOKAMAK, STELLARATOR])
    return f"Unknown topic '{topic}'. Available: overview, tokamak, stellarator, all."
