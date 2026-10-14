"""TEAM special problems: 32 (vector hysteresis), 30b/34 (motor), 18/19/29 (HF)."""

HYSTERESIS_32 = r"""
# TEAM Problem 32 — Vector Hysteresis (★ canonical hysteresis test)

## Geometry

T-shaped magnetic core with applied rotating B field. Vector hysteresis
(non-aligned H and B trajectories) must be captured.

```
T-core: laminated electrical steel
Excitation: ROTATING applied B field (B_x, B_y vary periodically)
Measurement: H(t) trajectory vs B(t) trajectory → vector lissajous
```

## What it tests

- ★ Vector hysteresis modeling (NOT scalar hysteresis)
- Power loss prediction under rotational excitation
- Comparison: Energy-Based (Henrotte) vs Vector Play vs Vector Stop
- Comparison with measured rotational loss

## Lab usage

★ THE benchmark for `radia_mcp.magnetic_materials.hysteresis_models`:
- **MatEnergyHysteresis** (★ lab CORE method, B-input Stop+Energy)
- MatPlayHysteresis (alternative)
- Compare with measured vector loci

## Reference

[LOCAL] `05_TEAM_benchmark/problem32/problem32_A Test-Case for Validation of Magnetic Field Analysis with Vector Hysteresis.pdf`

## Cross-references

- `radia_mcp.magnetic_materials.hysteresis_models.lab_core` ★
- `radia_mcp.magnetic_materials.hysteresis_models.energy_based`
- `radia_mcp.magnetic_materials.hysteresis_models.vector`
"""


MOTORS_30B_34 = r"""
# TEAM Problems 30b, 34 — Motor analytic references

## Problem 30b: Single and Three-Phase Induction Motor

```
3-phase induction motor with analytical reference solution
(from Boldea's classic IM textbook)
Used to validate motor FE solvers against published analytic
```

## Problem 34: Spherical Induction Motor with Anisotropic Rotor

```
Novel spherical IM geometry with anisotropic μ_r tensor in rotor
Tests anisotropic material handling + 3D rotor field
```

## Lab usage

- Cross-link to `radia_mcp.motor` for motor design
- Validate Hollaus MSFEM with motor geometry
- Frljić anisotropic μ_eff testing (problem 34)

## References

[LOCAL] `05_TEAM_benchmark/problem30b/problem30b_Analytic Analysis of Single and Three Phase Induction Motors.pdf`
[LOCAL] `05_TEAM_benchmark/problem34/problem34_Spherical induction motor with anisotropic rotor.pdf`

## Cross-references

- `radia_mcp.motor` — motor analysis MCP
- `radia_mcp.motor.hollaus_eddy` — Frljić anisotropic
"""


HIGH_FREQ_18_19_29 = r"""
# TEAM Problems 18, 19, 29 — High-frequency cavity / RF (out of lab scope)

| # | Name | Physics |
|---|------|---------|
| 18 | Waveguide Loaded Cavity | full Maxwell, resonant modes |
| 19 | Microwave Field in Loaded Cavity | scattering + resonance |
| 29 | Whole Body Cavity Resonator | full-wave biology |

## Why out of lab scope

Per CLAUDE.md: lab is **Laplace kernel only** (MQS / Darwin). High-freq
Helmholtz / full-wave is NOT a radia focus area.

## Reference for users who need HF

If you need to validate a full-wave Maxwell solver:
- Use ngsolve.bem with Helmholtz kernel
- These problems are standard benchmarks for HF tools (HFSS, FEKO, CST)

## References

[LOCAL] `05_TEAM_benchmark/problem18/problem18_Waveguide Loaded Cavity.pdf`
[LOCAL] `05_TEAM_benchmark/problem19/problem19_Microwave Field in a Loaded Cavity.pdf`
[LOCAL] `05_TEAM_benchmark/problem29/problem29_Whole body cavity resonator.pdf`
"""


REFERENCES = r"""
# TEAM Workshop key reference documents

Located in `_references/`:

| File | What |
|------|------|
| Proceeding of IGTE Symposium 1990.pdf | First TEAM compilation |
| Solution of TEAM workshop problem No. 7 by FEM.pdf | Problem 7 FEM solution |
| An Integral Formulation for the Computation of 3-D Eddy Current Using Facet Elements.pdf | Method paper |
| Basic equations for deformation under external efforts.pdf | Force-coupled solid |

## Where to find more TEAM info

- **TEAM Workshop home**: https://www.compumag.org/wp/team/
- Active problems: published in Compumag / Cefc / IEEE TMAG proceedings
- New problems added periodically (most recent: problem 36+)

## Lab integration

The lab uses TEAM problems as:
1. **Validation fixtures** for solver evidence (`validation_test/panels/`)
2. **Citation** in research papers ("validated against TEAM problem N")
3. **Inter-comparison** when comparing methods (Maxwell vs Coulomb force)
"""


def get_special_knowledge(topic: str = "hysteresis_32") -> str:
    """Dispatch TEAM special problem topics.

    Topics:
        hysteresis_32     - ★ Problem 32 vector hysteresis (lab core)
        motors_30b_34     - Motor analytic references
        high_freq         - Problems 18, 19, 29 (HF, out of lab scope)
        references        - Key reference documents in _references/
        all               - Everything
    """
    topic = topic.lower().strip().replace("-", "_")
    if topic in ("hysteresis_32", "32", "problem_32", "vector_hysteresis"):
        return HYSTERESIS_32
    if topic in ("motors_30b_34", "motors", "30b", "34"):
        return MOTORS_30B_34
    if topic in ("high_freq", "hf", "cavity"):
        return HIGH_FREQ_18_19_29
    if topic in ("references", "refs"):
        return REFERENCES
    if topic == "all":
        return "\n\n".join([HYSTERESIS_32, MOTORS_30B_34, HIGH_FREQ_18_19_29,
                            REFERENCES])
    return (f"Unknown topic '{topic}'. Available: hysteresis_32, motors_30b_34, "
            "high_freq, references, all.")
