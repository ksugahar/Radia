"""TEAM eddy current problems: 1a, 1b, 2, 3, 4, 5, 7, 9, 21, 24."""

OVERVIEW = r"""
# TEAM eddy current problems summary

| # | Name | Eddy class | Lab use |
|---|------|------------|---------|
| 1a, 1b | FELIX Long Cylinder | 3D transient | early benchmark |
| 2 | Infinitely Long Cylinder | harmonic, analytic | quick validation |
| 3 | Bath Plate with 2 Holes | harmonic, multiply-connected | A-V Biro-Preis test |
| 4 | FELIX Brick | 3D transient | step response |
| 5 | Bath Cube | 3D + open BC | Kelvin transform test |
| **7** | Asymmetrical Conductor w/ Hole ★ | 3D harmonic | ★ most cited |
| 9 | Velocity Effects | moving conductor | v × B term |
| 21 | Family (Stray Field) | 漂遊磁界 | TEAM 21 series |
| 24 | Nonlinear Time-Transient Test Rig | transient + nonlinear | rig hardware ref |

★ = THE standard 3D harmonic eddy current benchmark.
"""


PROBLEM_07 = r"""
# TEAM Problem 7 — Asymmetrical Conductor with a Hole (★ canonical 3D eddy)

## Geometry

A thick aluminum plate with an OFF-CENTER circular hole, excited by
a sinusoidal coil current above the plate.

```
Plate: 280 x 140 x 19 mm, sigma = 35.5 MS/m (aluminum)
Hole: r=18 mm, center off-axis by 40 mm
Coil: 100 x 100 mm above plate, 2742 turns, 30 A * sin(2π * 50 Hz)
Air gap: 12 mm above plate, 22 mm below
```

## What it tests

- ★ 3D harmonic eddy current FEM/BEM
- Complex geometry (off-center hole)
- Skin depth resolution near hole edge
- 50 Hz response (relevant to power frequency)

## Reference solutions

Multiple published curves of Eddy current density on plate top surface
along measurement lines. Experimental measurements + numerical results
from Compumag/Cefc proceedings.

## Lab usage

- ★ Standard benchmark for any new eddy current implementation
- Used to validate `calc_fem_kelvin.py` FEM-SIBC + Kelvin
- Used to validate ngsolve.bem single-trace (Weggler) at low freq
- 50 Hz → quasi-static regime, perfect for lab MQS focus

## Reference

[LOCAL] `05_TEAM_benchmark/problem07/problem7_Asymmetrical Conductor with a Hole.pdf` (63 MB - huge file)

## Cross-references

- `radia_mcp.ih.sibc` — SIBC formulation tested against problem 7
- `radia_mcp.bem.surface_ie` — ngbem against problem 7
- `radia_mcp.fem.potential_formulations.a_omega` — A-V formulation tested
"""


PROBLEM_09 = r"""
# TEAM Problem 9 — Velocity Effects and Low Level Fields

## Geometry

Conducting bar in axial motion within a coil's magnetic field. Tests
the v × B term coupling.

```
Bar: cylindrical conductor, σ varying
Coil: external, excites B_0
Bar velocity: imposed (parametric)
```

## What it tests

- Moving conductor eddy current
- v × B contribution to E
- Lorentz force on bar
- Low-amplitude regime (linear scaling test)

## Lab usage

- Validates `radia_mcp.differential_forms.forces.moving_conductor`
- Validates motor air-gap analysis (analogous physics)

## Reference

[LOCAL] `05_TEAM_benchmark/problem09/problem9_Velocity Effects and Low Level Fields.pdf`

## Cross-reference

- `radia_mcp.motor` — rotor motion in motor analysis
- `radia_mcp.differential_forms.forces` — Lorentz force theory
"""


PROBLEM_24 = r"""
# TEAM Problem 24 — Nonlinear Time-Transient Rotational Test Rig

## Geometry

Actual hardware rig at IGTE Graz: rotational machine for transient
nonlinear iron testing.

```
Rotor: solid steel cylinder
Stator: 3-phase winding
Operating: time-stepped transient with rotor rotation
Iron: nonlinear BH curve provided
```

## What it tests

- ★ Time-transient FEM with nonlinear iron
- Rotor-stator coupling (relative motion)
- Newton-Raphson convergence over many time steps
- Comparison with rig measurements

## Lab usage

- Validates `radia.calc_motor_transient`
- Validates Hollaus MSFEM with TD coupling
- Cross-link to `radia_mcp.motor.calc_motor_transient`

## Reference

[LOCAL] `05_TEAM_benchmark/problem24/problem24_Nonlinear Time-Transient Rotational Test Rig.pdf`

## Cross-references

- `radia_mcp.motor` — motor transient solver
- `radia_mcp.fem.time_domain_axisym.time_domain_fem` — TD-FEM theory
"""


def get_eddy_current_knowledge(topic: str = "overview") -> str:
    """Dispatch TEAM eddy current problem topics.

    Topics:
        overview     - All eddy current problems summary (DEFAULT)
        problem_07   - ★ Asymmetrical Conductor with Hole (canonical 3D)
        problem_09   - Velocity Effects (moving conductor)
        problem_24   - Nonlinear Time-Transient Rig
        all          - Everything
    """
    topic = topic.lower().strip().replace("-", "_")
    if topic in ("overview", "all_eddy"):
        return OVERVIEW
    if topic in ("problem_07", "07", "problem7", "asymmetrical"):
        return PROBLEM_07
    if topic in ("problem_09", "09", "problem9", "velocity"):
        return PROBLEM_09
    if topic in ("problem_24", "24", "problem24", "nonlinear_transient"):
        return PROBLEM_24
    if topic == "all":
        return "\n\n".join([OVERVIEW, PROBLEM_07, PROBLEM_09, PROBLEM_24])
    return (f"Unknown topic '{topic}'. Available: overview, problem_07, "
            "problem_09, problem_24, all.")
