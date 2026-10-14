"""TEAM Workshop catalog — all 30 problems organized by physics class."""

OVERVIEW = r"""
# TEAM Workshop benchmark catalog

The TEAM (Testing Electromagnetic Analysis Methods) Workshop is the
international standard validation suite for EM-FEM/BEM solvers.
Started 1985; problems continually added/updated.

## All TEAM problems (local copies in 05_TEAM_benchmark/)

| # | Name | Physics | Year | Folder |
|---|------|---------|------|--------|
| 1a | FELIX Long Cylinder | 3D eddy current | 1985 | problem01a/ |
| 1b | FELIX Long Cylinder | 3D eddy current | 1985 | problem01b/ |
| 2 | Infinitely Long Cylinder | sinusoidal eddy | 1985 | problem02/ |
| 3 | Bath Plate with 2 Holes | 3D eddy + holes | 1985 | problem03/ |
| 4 | FELIX Brick | 3D eddy transient | 1985 | problem04/ |
| 5 | Bath Cube | 3D eddy + open BC | 1985 | problem05/ |
| 6 | Sphere in Uniform Field | analytic ref | 1985 | problem06/ |
| 7 | Asymmetrical Conductor with Hole | 3D eddy, holes ★ | 1986 | problem07/ |
| 8 | Coil Above a Crack | NDT eddy current | 1986 | problem08/ |
| 9 | Velocity Effects | moving conductor | 1986 | problem09/ |
| 13 | 3-D Non-Linear Magnetostatic | nonlinear yoke ★ | 1988 | problem13/ |
| 15 | Rectangular Slot in Thick Plate | NDT (crack) | 1989 | problem15/ |
| 17 | Jumping Ring | Lorentz force transient | 1989 | problem17/ |
| 18 | Waveguide Loaded Cavity | high-freq | 1990 | problem18/ |
| 19 | Microwave Field in Loaded Cavity | high-freq | 1990 | problem19/ |
| 20 | 3-D Static Force Problem | force computation ★ | 1990 | problem20/ |
| 21 | Family (Stray Field) | 漂遊磁界 | 1992 | problem21/ |
| 22 | SMES Optimization | optimization | 1995 | problem22/ |
| 23 | Forces in Permanent Magnets | PM force ★ | 1996 | problem23/ |
| 24 | Nonlinear Time-Transient Test Rig | nonlinear + transient | 1996 | problem24/ |
| 25 | Die Press Optimization | optimization | 1996 | problem25/ |
| 27 | NDT Deep Flaws | NDT (eddy NDT) | 2000 | problem27/ |
| 28 | Electrodynamic Levitation | levitation force | 2000 | problem28/ |
| 29 | Whole Body Cavity Resonator | high-freq biology | 2001 | problem29/ |
| 30b | Single & Three-Phase Induction Motor | analytic IM ref | 2001 | problem30b/ |
| 31 | Low Freq Currents Reconstruction | inverse problem | 2002 | problem31/ |
| 32 | Vector Hysteresis | vector hysteresis ★ | 2002 | problem32/ |
| 33b | Electric Local Force | local force formulation ★ | 2007 | problem33b/ |
| 34 | Anisotropic Spherical IM | anisotropic motor | 2008 | problem34/ |
| 36 | (recent) | recent | 2014+ | problem36/ |

★ = high lab relevance (force / nonlinear / vector hysteresis)

## Physics class breakdown

| Class | Problems | Lab MCP relevance |
|-------|----------|-------------------|
| **Magnetostatic** (linear, nonlinear) | 6, 13, 20, 23 | `radia_mcp.electromagnet`, `radia_mcp.radia_ngsolve.hdiv_vim` |
| **Eddy current** (3D, transient, harmonic) | 1a, 1b, 2, 3, 4, 5, 7, 9, 21, 24 | `radia_mcp.ih`, `radia_mcp.peec`, `radia_mcp.fem` |
| **Force / motion / levitation** | 17, 20, 23, 28, 33b | `radia_mcp.differential_forms.forces` |
| **NDT (cracks, flaws)** | 8, 15, 27 | `radia_mcp.ndt` |
| **Inverse / optimization** | 22, 25, 31 | `radia_mcp.topology_optimization` |
| **Hysteresis** | 32 | `radia_mcp.magnetic_materials.hysteresis_models` |
| **High-freq cavity / RF** | 18, 19, 29 | (out of lab scope; ngsolve.bem) |
| **Motor analytic ref** | 30b, 34 | `radia_mcp.motor` |
"""


VALIDATION_MAP = r"""
# Validation map: technique → which TEAM problem to use

When you implement a new EM solver feature, the corresponding TEAM
problem provides the reference solution.

## Magnetostatic solvers

| Feature | Validate with | Why |
|---------|---------------|-----|
| Linear magnetostatic (Phi) | problem06 (Sphere) | analytic comparison |
| Nonlinear yoke | problem13 (3-D Non-Linear Magnetostatic) | published reference |
| HDiv-VIM / reduced FEM magnetic material | problem13 | nonlinear reference |
| HDiv-VIM sphere / open boundary | problem06 | exact ref |

## Eddy current solvers

| Feature | Validate with | Why |
|---------|---------------|-----|
| 3D harmonic eddy current | problem07 (Asymmetrical w/ Hole) | most cited; complex geometry |
| 3D transient eddy current | problem04 (FELIX Brick) | step input response |
| Cylindrical harmonic | problem02 (Infinitely Long Cylinder) | analytic |
| Multiply connected eddy | problem03 (Bath Plate w/ 2 Holes) | topology test |
| Moving conductor | problem09 (Velocity Effects) | v × B coupling |

## Force methods

| Feature | Validate with | Why |
|---------|---------------|-----|
| Maxwell stress | problem20 (3-D Static Force) ★ | lab-implemented; clean ref |
| Permanent magnet force | problem23 (Forces in PMs) ★ | exact PM source |
| Local force formulations | problem33b (Experimental Validation Local Force) ★ | per-node verification |
| Coulomb / nodal force | problem20 or problem33b | comparison set |
| Eggshell + virtual material | problem20 + problem33b | local vs global |

## NDT

| Feature | Validate with | Why |
|---------|---------------|-----|
| Coil above crack | problem08 | classic NDT |
| Conductor with slot | problem15 | rect slot canonical |
| Deep flaws | problem27 | depth-dependent ECT |

## Hysteresis / nonlinear

| Feature | Validate with | Why |
|---------|---------------|-----|
| Vector hysteresis | problem32 ★ | THE canonical test |
| Nonlinear iron (BH curve) | problem13, problem24 | static + transient |
| Time-periodic hysteresis | problem24 | transient nonlinear |

## Inverse / optimization

| Feature | Validate with | Why |
|---------|---------------|-----|
| Topology optimization | problem25 (Die Press) | shape opt benchmark |
| SMES coil opt | problem22 | parametric opt |
| Current reconstruction | problem31 | inverse problem |

## Cross-reference

- `radia_mcp.electromagnet` — uses problem23 for PM force validation
- `radia_mcp.differential_forms.forces` — uses problem20 + problem33b
- `radia_mcp.magnetic_materials.hysteresis_models` — uses problem32
- `radia_mcp.topology_optimization` — uses problem22, problem25
"""


def get_catalog_knowledge(topic: str = "overview") -> str:
    """Dispatch TEAM catalog topics.

    Topics:
        overview         - All 30 problems by physics class (DEFAULT)
        validation_map   - Which TEAM problem validates technique X
        all              - Everything
    """
    topic = topic.lower().strip()
    if topic in ("overview", "all_problems", "catalog"):
        return OVERVIEW
    if topic in ("validation_map", "validation", "which_problem"):
        return VALIDATION_MAP
    if topic == "all":
        return "\n\n".join([OVERVIEW, VALIDATION_MAP])
    return (f"Unknown topic '{topic}'. Available: overview, validation_map, all.")
