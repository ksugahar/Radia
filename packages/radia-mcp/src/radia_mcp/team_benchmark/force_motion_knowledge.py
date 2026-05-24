"""TEAM force/motion problems: 17, 20, 23, 28, 33b — lab core."""

OVERVIEW = r"""
# TEAM force / motion / levitation problems

| # | Name | Force class | Lab use |
|---|------|------------|---------|
| 17 | Jumping Ring | Lorentz transient | demo / education |
| **20** | 3-D Static Force ★ | Maxwell stress, magnetostatic | lab production |
| **23** | Forces in Permanent Magnets ★ | PM force | lab production |
| 28 | Electrodynamic Levitation | levitation, force balance | maglev research |
| **33b** | Electric Local Force ★ | local force formulations | per-node verification |

★ = lab core force validation problems.
"""


PROBLEM_20 = r"""
# TEAM Problem 20 — 3-D Static Force Problem (★ lab core)

## Geometry

Two ferromagnetic blocks (or block + plate) with a permanent magnet or
coil between. Compute the static force on each block.

```
Block 1: mu_r = 1000, dimensions 50x40x20 mm
Block 2: mu_r = 1000, dimensions 50x40x20 mm
Gap: 5 mm (varies in parametric study)
Source: PM (Nd magnetization) or DC coil
```

## What it tests

- ★ Maxwell stress tensor implementation
- Force on iron in static B field
- Comparison of force methods (Maxwell, Coulomb, Kameari, VWP)
- Lab-implemented in `calc_em_force.py`

## Published reference

Multi-lab benchmark with measured forces at multiple gap sizes.

## Lab usage

★ THE canonical force benchmark for `radia_mcp.differential_forms.forces`.
All 7 implemented methods (Maxwell, Coulomb nodal, Kameari, Lorentz,
Arkkio, Henrotte, VWP) are validated here.

| Method | Lab implementation | Validates |
|--------|---------------------|-----------|
| Maxwell stress | `calc_em_force.py --method maxwell` | global force |
| Coulomb nodal | `calc_em_force.py --method coulomb_nodal` | per-node |
| Kameari (1993) | `calc_em_force.py --method kameari` | local distribution |
| Henrotte (2004) | `calc_em_force.py --method henrotte` | EM force density |

## Reference

[LOCAL] `05_TEAM_benchmark/problem20/problem20_3-D Static Force Problem.pdf`
[LOCAL] `05_TEAM_benchmark/problem20/problem20_tech_3ax.pdf` (technical spec)

## Cross-reference

- `radia_mcp.differential_forms.forces` — all 7 force methods
- `radia_mcp.electromagnet` — Hantila solver applied to problem 20
- `memory/feedback_em_force_team20_validation.md` (lab memory)
"""


PROBLEM_23 = r"""
# TEAM Problem 23 — Forces in Permanent Magnets (★ lab core)

## Geometry

Two permanent magnets with various pole arrangements; compute the
mutual force as a function of separation and orientation.

```
PM1: NdFeB, 30x30x10 mm, Br = 1.4 T
PM2: same, separation parametric
Pole arrangements: attractive, repulsive, parallel, perpendicular
```

## What it tests

- ★ PM force computation (no current source, just M)
- Treatment of constant-M material in force integral
- Verification of force on M-fixed body

## Lab usage

★ Canonical benchmark for PM force in `calc_em_force.py`.
Used to validate:
- `permanent_magnet_force.py` lab implementation
- Radia's `ObjHexahedron(verts, [Mx,My,Mz])` PM source
- Per-element force density export

## Reference

[LOCAL] `05_TEAM_benchmark/problem23/problem23_Forces in Permanent Magnets.pdf`

## Cross-reference

- `radia_mcp.differential_forms.forces.permanent_magnet_force` ★
- `radia_mcp.magnetic_materials.permanent_magnet` — PM datasheets
"""


PROBLEM_33B = r"""
# TEAM Problem 33b — Electric Local Force (★ lab per-node validation)

## Geometry

C-shaped electromagnet pulling a movable plunger. Iron + coil + plunger
with measured forces at multiple positions.

```
C-yoke: nonlinear steel (BH provided)
Coil: known turns + current
Plunger: same steel, movable along x
Position scan: x = 0..50 mm
```

## What it tests

- ★ Local force formulations (force per node / per element)
- Comparison between MULTIPLE published numerical methods
- "Experimental Validation of Electric Local Force Formulations"
  (V3 version is the most cited)
- Verification of PER-NODE Kameari force decomposition

## Lab usage

★ Used for `radia_mcp.differential_forms.forces.local_methods`:
- Per-node Kameari (Kameari 1993)
- Iino-Okamoto two-stage error correction
- Eggshell method
- Sensitivity-method VWP (Pile 2018)

These methods provide DISTRIBUTION of force, not just global integral.

## Reference

[LOCAL] `05_TEAM_benchmark/problem33b/problem-33b_V3_Experimental Validation of Electric Local Force Formulations.pdf`

## Cross-reference

- `radia_mcp.differential_forms.forces.local_methods` ★
- `radia_mcp.differential_forms.forces.iino_okamoto`
- Lab memory: `tests/panels/test_calc_em_force_team33b_golden.py`
"""


PROBLEM_28 = r"""
# TEAM Problem 28 — Electrodynamic Levitation Device

## Geometry

Conducting plate above an axisymmetric coil set with AC excitation.
The plate floats due to eddy current Lorentz force.

```
Plate: aluminum, 50x50x5 mm
Coil: 3-axis axisymmetric (multiple windings)
Excitation: AC at 50 Hz, 10 A rms
Equilibrium height: ~5 mm
```

## What it tests

- ★ Levitation force balance (gravity vs eddy Lorentz)
- 3D harmonic eddy current with motion (force depends on height)
- Comparison with measured levitation height

## Lab usage

- Cross-link to `radia_mcp.maglev` (future MCP) for maglev physics
- Cross-link to `radia_mcp.ih` for eddy current in similar config

## Reference

[LOCAL] `05_TEAM_benchmark/problem28/` (11 files, 15 MB)
- problem28 description
- numerical solutions, measured data

## Cross-references

- `radia_mcp.differential_forms.forces.lorentz` — Lorentz force
- `radia_mcp.ih.sibc` — analogous eddy + force physics
"""


def get_force_motion_knowledge(topic: str = "overview") -> str:
    """Dispatch TEAM force / motion problem topics.

    Topics:
        overview     - All force/motion problems summary (DEFAULT)
        problem_20   - ★ 3-D Static Force Problem (lab core)
        problem_23   - ★ Forces in Permanent Magnets (lab core)
        problem_33b  - ★ Electric Local Force (per-node validation)
        problem_28   - Electrodynamic Levitation
        all          - Everything
    """
    topic = topic.lower().strip().replace("-", "_")
    if topic in ("overview", "all_force"):
        return OVERVIEW
    if topic in ("problem_20", "20", "problem20", "static_force"):
        return PROBLEM_20
    if topic in ("problem_23", "23", "problem23", "pm_force"):
        return PROBLEM_23
    if topic in ("problem_33b", "33b", "problem33b", "local_force"):
        return PROBLEM_33B
    if topic in ("problem_28", "28", "problem28", "levitation"):
        return PROBLEM_28
    if topic == "all":
        return "\n\n".join([OVERVIEW, PROBLEM_20, PROBLEM_23, PROBLEM_33B,
                            PROBLEM_28])
    return (f"Unknown topic '{topic}'. Available: overview, problem_20, "
            "problem_23, problem_33b, problem_28, all.")
