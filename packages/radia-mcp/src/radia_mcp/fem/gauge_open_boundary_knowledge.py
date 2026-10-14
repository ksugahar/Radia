"""Gauging + open boundary techniques (cross-links matrix_solvers and bem)."""

GAUGE = r"""
# Gauging for HCurl A formulations

A is defined only up to gradient: A' = A + ∇φ gives same B = curl A.
This freedom appears as null space in the discrete operator.

## Gauge choices

| Gauge | Constraint | When |
|-------|-----------|------|
| **Tree-cotree (Albanese-Rubinacci 1988)** | A=0 on spanning tree edges | small/medium meshes |
| **Coulomb (∇·A = 0)** | penalty (∇·u, ∇·v) | simple but adds positive operator |
| **A-V (Biro-Preis 1989)** | scalar V supplements A | multiply-connected conductor |
| **Ungauged + AMS** | live with kernel, use Hiptmair-Xu AMS | ★ lab default |
| **Two-Step Maxwell (Ostrowski-Hiptmair 2021)** | generalized tree-cotree, freq-stable | full Maxwell DC to MHz |

## Lab default: ungauged + CompactAMS

[CODE] `radia.sparsesolv_ngsolve.CompactAMSPreconditioner`

The Hiptmair-Xu AMS preconditioner has a **gradient subspace correction**:
it applies a scalar AMG solve on H1 to handle the gradient null space.
Result: no gauge is needed at the system level — the preconditioner
absorbs the kernel.

| Tradeoff | gauged | ungauged + AMS |
|----------|--------|----------------|
| System size | smaller (no extra scalar) | larger by H1 |
| Setup | gauge logic | none |
| Iteration count | similar | ~10% more |
| Code complexity | medium | low ★ |

## When to use which

| Problem | Recommended |
|---------|-------------|
| Static magnetostatic, simple | ungauged + CompactAMS ★ |
| Eddy current, multiply connected conductors | A-V Biro-Preis (V handles loop ambig) |
| Need DC-to-MHz frequency stability | Two-Step Maxwell (Ostrowski-Hiptmair) |
| Legacy code / textbook problem | tree-cotree |

## References (this folder)

[LOCAL] `10_FEM_定式化/11_ゲージ/` (11 files, 136 MB)

## Cross-references

- `radia_mcp.matrix_solvers.em_specific` — solver-side view (gauge from
  preconditioner perspective)
- `radia_mcp.bem.low_freq` — analogous gauging in BEM (Loop-Star,
  Calderón)
"""


KELVIN_TRANSFORM = r"""
# Kelvin transformation — open boundary for FEM (lab default)

Convert unbounded exterior R^3 \ V_compact to a bounded "Kelvin ball"
of radius R_K via the inversion map r' = R_K^2 / |r|.

## Why use Kelvin

| Alternative | Issue |
|-------------|-------|
| Truncated outer box + Dirichlet | accuracy ~1/R_far → expensive mesh |
| PML (Perfectly Matched Layer) | tuning required, fails near zero freq |
| Asymptotic BC | accuracy ~ 1/R_far^n → still needs large mesh |
| BEM coupling (FEM-BEM hybrid) | implementation complexity |
| **Kelvin transform** | exact for static; clean FEM matrix; no parameter |

## How it works

The Laplace equation Δφ = 0 in R^3 transforms to:
    Δ' φ' = 0   in Kelvin ball

with appropriate Jacobian weighting. Combine the interior FEM problem
in V_compact with the Kelvin FEM problem in V_kelvin via Robin matching
on the interface sphere.

## Implementation in radia/NGSolve

```python
# calc_fem_kelvin.py (radia panel)
# Solve inside V_iron + air + V_kelvin (compact ball)
# Apply Kelvin transform on interface → open BC
# All in pure FEM (NGSolve)
```

The lab uses this for:
- IH workpiece (radia.calc_fem_kelvin)
- Accelerator magnets with open boundary
- PEEC + FEM hybrid validation

## Identification of Kelvin pair

For the Kelvin BC, the FE space must identify the interior interface
with the Kelvin-mapped exterior interface. This is done at the
C++ level in `cubit_plugin/Mark.cpp` (or NGSolve `Periodic` on
labelled boundaries with identification map).

## Verify-first

Per CLAUDE.md "Verify-First Policy": always check the FES setup
BEFORE running a 5-minute physics solve. The minimum verify trio:

```python
print(mesh.GetMaterials(), mesh.GetBoundaries())
# Kelvin BC must show slaved DOFs > 0
fb = H1(mesh, order=p)
fp = Periodic(fb)  # or Kelvin-identified version
slaved = sum(fb.FreeDofs()) - sum(fp.FreeDofs())
# Functional test:
gfu = GridFunction(fp); gfu.vec[:] = 0
gfu.Set(1.0, definedon=mesh.Boundaries("kelvin_int"))
r = Integrate(gfu*gfu, mesh, definedon=mesh.Boundaries("kelvin_ext")) \
  / Integrate(gfu*gfu, mesh, definedon=mesh.Boundaries("kelvin_int"))
# r must equal 1.0
```

## Cross-references

- `radia_mcp.radia_ngsolve.kelvin_transformation` — concrete code recipe
- `radia_mcp.ih` — Kelvin in IH workpiece SIBC
- `radia_mcp.mathematica.recipes.KELVIN_TRANSFORM` — symbolic derivation
"""


ABC_PML = r"""
# Absorbing Boundary Conditions (ABC) and PML — open boundary alternatives

## Asymptotic BC (ABC)

Truncate the domain at radius R_far and enforce a Robin BC that
approximates the radiation condition:

    ∂u/∂n + α u = 0    on outer boundary

For Laplace problems: α = 1/r_far (1st order) or higher derivatives.

| Order | Accuracy | When |
|-------|----------|------|
| 1st (Sommerfeld) | O(1/R_far) | quick prototyping |
| 2nd (BGT) | O(1/R_far^2) | better; standard |
| 3rd+ (Engquist-Majda) | higher | rarely needed |

## PML (Perfectly Matched Layer)

Berenger 1994. Add an absorbing layer with complex coordinate stretching:

    x → x + i σ(x)/ω    inside PML

Theoretically reflectionless for ALL frequencies and incident angles
(continuum). Discretized PML has reflections, especially at low freq
(ω→0 → 1/ω blowup in coordinate stretch).

| PML variant | Use |
|-------------|-----|
| Split-field (Berenger original) | classic FDTD |
| Stretched-coordinate (Chew-Weedon 1994) | FEM, BEM, FDTD |
| UPML (Sacks-Kingsland-Lee-Lee 1995) | physically realizable medium |
| CFS-PML (complex frequency-shifted) | low-freq stable |

## When NOT to use PML

✗ **Low frequency / static** — PML fails as ω→0 (use Kelvin or BEM instead)
✗ Iterative solvers may struggle with PML's complex matrix
✗ Geometry close to PML → reflections

## Lab use

- Lab is **Laplace kernel only** (MQS/Darwin) per CLAUDE.md
- PML is for full-wave (Helmholtz), NOT for lab production
- Kelvin transform is the production open-boundary technique

## References (this folder)

[LOCAL] `10_FEM_定式化/26_開領域解法/` (71 files, 425 MB) — comprehensive
        open boundary collection including ABC, PML, Kelvin, FEM-BEM hybrid

## Cross-reference

- `radia_mcp.bem.fem_bem_hybrid` — when to use FEM-BEM coupling instead
"""


def get_gauge_open_boundary_knowledge(topic: str = "gauge") -> str:
    """Dispatch gauge / open-boundary topics.

    Topics:
        gauge              - Gauge choices for HCurl A (DEFAULT)
        kelvin_transform   - Kelvin transformation (★ lab default for open BC)
        abc_pml            - Asymptotic BC and PML alternatives
        all                - Everything
    """
    topic = topic.lower().strip()
    if topic in ("gauge", "gauging"):
        return GAUGE
    if topic in ("kelvin_transform", "kelvin", "kelvin_transformation"):
        return KELVIN_TRANSFORM
    if topic in ("abc_pml", "abc", "pml", "open_boundary"):
        return ABC_PML
    if topic == "all":
        return "\n\n".join([GAUGE, KELVIN_TRANSFORM, ABC_PML])
    return (f"Unknown topic '{topic}'. Available: gauge, kelvin_transform, "
            "abc_pml, all.")
