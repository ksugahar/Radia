"""TEAM magnetostatic problems: 6, 13."""

PROBLEM_06 = r"""
# TEAM Problem 6 — Sphere in Uniform Magnetic Field

## Geometry

A magnetic sphere (mu_r) in an externally applied uniform B_0.
The classic textbook problem with closed-form solution.

```
Sphere: r = 100 mm, mu_r = 1000 (linear)
B_ext = 1 T uniform (in z direction)
Air region: 500 mm box around sphere
```

## Analytic solution

Internal field (uniform):
```
B_in = B_ext * 3 mu_r / (mu_r + 2)    along z
```

For mu_r = 1000: B_in ≈ 2.997 T (3x boost vs B_ext)

External field: dipole superposed on B_ext.

## What it tests

- ★ Basic FEM/BEM/HDiv-VIM implementation correctness
- Open boundary handling (Kelvin transform vs PML vs ABC)
- Material discontinuity at sphere surface
- Higher-order convergence (analytic → exact ref)

## Lab usage

- Used for HDiv-VIM / FEM-Kelvin validation
- Used for FEM-Kelvin open-boundary validation
- Sub-millimeter B accuracy at center → < 0.01% error target

## Reference

[LOCAL] `05_TEAM_benchmark/problem06/problem6_Sphere in Uniform Magnetic Field.pdf` (10 MB)
[REF] L.R. Turner, IEEE TMAG 1985-1988 (TEAM proceedings)
"""


PROBLEM_13 = r"""
# TEAM Problem 13 — 3-D Non-Linear Magnetostatic Model (★ canonical nonlinear)

## Geometry

A C-core electromagnet with nonlinear steel plates and coil winding.

```
Iron yoke: C-shape, 250x250x100 mm, nonlinear BH curve
Air gap: 50 mm
Coil: 12 turns, 100 A current
Vertical plate: parallel to coil axis (key feature)
```

## What it tests

- ★ Nonlinear magnetostatic FEM (Picard / Newton-Raphson)
- 3D BH curve interpolation
- Multiply-connected iron yoke (cohomology)
- Convergence of nonlinear iteration

## Published reference solutions

| Solver | Method | Flux density along ref line |
|--------|--------|----------------------------|
| Various FEM (1988+) | A formulation | published in Compumag |
| Volume IE (Nakata) | published volume-integral reference | published |
| BEM (Onuki) | surface integral | published |
| Hybrid FEM-BEM | combined | published |

## Lab usage

- ★ Validates `radia_mcp.electromagnet` (Hantila polarization solver)
- Validates HDiv-VIM / reduced-FEM nonlinear magnetic-material paths
- Used for nonlinear Picard convergence study (lab benchmarks)

## Reference

[LOCAL] `05_TEAM_benchmark/problem13/problem13_3-D Non-Linear Magnetostatic Model.pdf` (15 MB)
[REF] T. Nakata et al., IEEE TMAG 24 (1988) — original specification

## Cross-references

- `radia_mcp.electromagnet` — Hantila polarization method
- `radia_mcp.radia_ngsolve.hdiv_vim` — HDiv-VIM magnetic material path
- `radia_mcp.matrix_solvers.preconditioners.shifted` — convergence for nonlinear
"""


def get_magnetostatic_knowledge(topic: str = "problem_06") -> str:
    """Dispatch TEAM magnetostatic problem topics.

    Topics:
        problem_06  - Sphere in Uniform Field (analytic ref, DEFAULT)
        problem_13  - ★ 3-D Non-Linear Magnetostatic Model (lab core)
        all         - Both
    """
    topic = topic.lower().strip().replace("-", "_")
    if topic in ("problem_06", "06", "problem6", "sphere"):
        return PROBLEM_06
    if topic in ("problem_13", "13", "problem13", "nonlinear"):
        return PROBLEM_13
    if topic == "all":
        return "\n\n".join([PROBLEM_06, PROBLEM_13])
    return f"Unknown topic '{topic}'. Available: problem_06, problem_13, all."
