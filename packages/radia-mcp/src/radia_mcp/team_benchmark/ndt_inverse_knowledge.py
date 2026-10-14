"""TEAM NDT + inverse / optimization problems: 8, 15, 22, 25, 27, 31."""

NDT_PROBLEMS = r"""
# TEAM NDT (Non-Destructive Testing) problems

| # | Name | NDT class |
|---|------|-----------|
| 8 | Coil Above a Crack | classic ECT |
| 15 | Rectangular Slot in Thick Plate | canonical slot |
| 27 | NDT Deep Flaws (EDDY CURRENT NDT) | depth-dependent |

These problems define standardized NDT scenarios for validating eddy
current testing solvers.

## Problem 8: Coil Above a Crack

```
Plate: 12.22 mm thick, sigma = 30.6 MS/m
Crack: 12.6 mm long x 0.28 mm wide x 5 mm deep
Coil: 3 mm above plate, 900 turns @ 900 Hz
```

Test: induced voltage in coil as it scans over crack. Reference data
from EM-NDE workshops.

## Problem 15: Rectangular Slot in Thick Plate

```
Plate: large, sigma = 30.6 MS/m
Slot: rectangular notch
Coil: scanning above slot
```

## Problem 27: EDDY CURRENT NDT AND DEEP FLAWS

Multi-layered plate with cracks at increasing depths. Tests:
- Skin depth penetration
- Crack-depth resolution
- Probe coil design

## Lab usage

★ Cross-link to `radia_mcp.ndt` (already implemented):
- `eddy_current_NDT` topic in ndt subpackage
- Used as benchmark for coil scan + signal processing

## References

[LOCAL] `05_TEAM_benchmark/problem08/problem8_Coil Above a Crack.pdf`
[LOCAL] `05_TEAM_benchmark/problem15/problem15_Rectangular Slot in a Thick Plate.pdf`
[LOCAL] `05_TEAM_benchmark/problem27/problem27_EDDY CURRENT NDT AND DEEP FLAWS.pdf`

## Cross-reference

- `radia_mcp.ndt.eddy_current_NDT` — production NDT solver
- `radia_mcp.bem.surface_ie` — BEM for NDT (lighter than FEM)
"""


INVERSE_OPT = r"""
# TEAM inverse / optimization problems

| # | Name | Type |
|---|------|------|
| 22 | SMES Optimization | shape opt + cost function |
| 25 | Die Press Optimization | topology opt |
| 31 | Low-Freq Currents Reconstruction | inverse problem |

## Problem 22: SMES Optimization

```
SMES coil (Superconducting Magnetic Energy Storage)
Design variables: coil cross-section, position
Cost: stored energy / stray field tradeoff
Constraints: superconductor critical current, mechanical stress
```

Classic shape optimization benchmark; published reference Pareto-optimal
designs.

## Problem 25: Die Press Optimization

Topology optimization of die-press pole shape to minimize iron mass
while maintaining required field.

## Problem 31: Low-Freq Currents Reconstruction

Inverse problem: given measured boundary B values, reconstruct internal
current distribution. Classic ill-posed inverse problem; tests
regularization techniques.

## Lab usage

★ Cross-link to `radia_mcp.topology_optimization`:
- Problem 22 → parametric shape opt benchmark
- Problem 25 → SIMP / ON-OFF / level set methods
- Problem 31 → inverse design / sensitivity analysis

## References

[LOCAL] `05_TEAM_benchmark/problem22/` (SMES, 2 files inc 'Family' variant)
[LOCAL] `05_TEAM_benchmark/problem25/problem25_Optimization of Die Press Model.pdf`
[LOCAL] `05_TEAM_benchmark/problem31/problem31_Reconstruction of Low Frequency Currents.pdf`

## Cross-reference

- `radia_mcp.topology_optimization` — topology opt MCP
- `radia_mcp.electromagnet` — accelerator magnet shape opt
"""


def get_ndt_inverse_knowledge(topic: str = "ndt_problems") -> str:
    """Dispatch TEAM NDT / inverse problem topics.

    Topics:
        ndt_problems   - Problems 8, 15, 27 (NDT) (DEFAULT)
        inverse_opt    - Problems 22, 25, 31 (inverse / optimization)
        all            - Both
    """
    topic = topic.lower().strip()
    if topic in ("ndt_problems", "ndt", "ect", "crack"):
        return NDT_PROBLEMS
    if topic in ("inverse_opt", "inverse", "optimization", "topology"):
        return INVERSE_OPT
    if topic == "all":
        return "\n\n".join([NDT_PROBLEMS, INVERSE_OPT])
    return f"Unknown topic '{topic}'. Available: ndt_problems, inverse_opt, all."
