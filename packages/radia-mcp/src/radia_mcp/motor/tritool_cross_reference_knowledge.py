"""Tri-tool cross-reference: FEMM / JMAG / radia-ngsolve (相互学習).

Cross-learning knowledge that ties the Sugahara-lab's three motor-FEA tools
together so each strengthens the others:

- **FEMM** (open-source 2D, Compumag standard) — ground-truth yardstick,
  pyFEMM headless cross-check, rich .fem model assets.
- **JMAG** (JSOL commercial, JMAG-Designer) — highest-capability: vector-Play
  hysteresis, native motion+slide, multi-slice skew + cage circuit, coupled
  IH, compiled user-subroutine DLL hooks, motor topology optimization.
- **radia-ngsolve** (in-house, NGSolve) — order-2 H1 accuracy, open Python
  extensibility, EM<->thermal auto-coupling via the lab IH solver.

Two dedicated knowledge MCP servers now hold each tool's distilled practice:

| Server | Location | Entry point |
|--------|----------|-------------|
| mcp-server-femm | lab-internal | `femm_*` tools |
| mcp-server-jmag | lab-internal | `jmag_*` tools |

This module is the radia-ngsolve-side mirror of their shared `cross_ref.json`:
it tells radia_mcp.motor where each capability maps across the three tools and
which strengthening steps have the highest payoff for radia-motor.

Built 2026-06 from: groups.io/g/femm member browsing (synthesis, no verbatim
bulk ingest) + pyFEMM 0.1.3 cross-checks + the lab's internal commercial-tool
corpus (kept lab-private).
"""

from __future__ import annotations


OVERVIEW = """\
## Tri-tool cross-learning — FEMM / JMAG / radia-ngsolve

The lab runs three motor-FEA tools. The point of cross-learning (相互学習) is
to import each tool's strengths into the others and use all three to raise the
verification bar:

- Import **FEMM's field-tested rules** (cogging period = 360/LCM(slots,poles),
  >=4 air-gap mesh layers, weighted-stress-tensor torque) and **JMAG's advanced
  capabilities** (vector-Play hysteresis, multi-slice skew + cage circuit,
  coupled IH) into radia-ngsolve.
- Use **radia-ngsolve's high-order accuracy** (order-2 H1) and **EM<->thermal
  auto-coupling** to complement and, where possible, exceed the others.
- Triangulate every result across the three tools (e.g. transverse-magnetized
  cylinder: analytic 0.2513 T / radia -0.05%; a stored independent 2D reference
  lands at -0.60%).

### Two new dedicated knowledge servers

- **mcp-server-femm** (lab-internal): FEMM motor-FEA practice, the pyFEMM API
  + ground-truth cross-check harness, .fem/.ans formats, pitfalls. Query with
  `femm_motor_topic`, `femm_pyfemm_api`, `femm_cross_reference`,
  `femm_radia_motor_roadmap`.
- **mcp-server-jmag** (lab-internal): JMAG study-type taxonomy, the Python
  scripting API, IM multi-slice skew, vector-hysteresis, user-subroutine
  hooks, lab know-how. Query with `jmag_study_types`, `jmag_scripting_api`,
  `jmag_cross_reference`.

(The FEMM/JMAG servers are lab-internal — commercial/license-constrained
content stays out of this published package; query them for tool specifics.)

Both servers carry the SAME `cross_ref.json`; this module is its radia-side
view. See also `motor_femm_transient` (Lange-Henrotte-Hameyer 2009) which is
the FEMM transient method already captured here.
"""

CAPABILITY_MATRIX = """\
## Capability matrix (who does what, who is strongest)

| Capability | FEMM | JMAG | radia-ngsolve | Strongest |
|---|---|---|---|---|
| 2D magnetostatic (A_z) | planar, 1st-order tri | Static2D | solve_planar_magnetostatic (order-2 H1) | **radia** (0.2512 vs stored ref 0.2498) |
| Nonlinear B-H | BH table + N-R | BhTable + N-R | _nonlinear (Picard ν(\\|B\\|)) | tie |
| Eddy / time-harmonic | complex A_z | freq-response + FEMコンダクタ | solve_planar_eddy (jωσ) | tie |
| Torque / force | mo_blockintegral(22/18/19/11) | nodal-force / Maxwell-stress | eggshell_torque_2d | tie |
| Rotation / motion | none (re-mesh per angle, script sweep) | Transient2D + slide band + 6-DOF motion | rotor-angle sweep (magdir rotate) | **JMAG** |
| IM multi-slice skew + cage | slice method + sinc factor | integrated N-slice + cage circuit | stack 2D solves; external cage circuit | **JMAG** |
| IM circuit-param ID (Meeker) | static L(ωs) -> M/Ll/Rr fit | quasi-steady IM (fixed slip) | solve_planar_eddy slip sweep + anti-periodic 1/4 | FEMM=authority / radia reproduced |
| Ld/Lq (frozen-permeability) | manual freeze + re-solve | turnkey dq | freeze converged ν -> dq perturb -> inductance_2d | JMAG turnkey / radia reachable |
| Core loss (rotating, B harmonics) | time-harmonic only; .ans hack | iron-loss + current-control harmonics | per-element FFT -> Bertotti/iGSE | **radia possible** (replace .ans hack) |
| Demag (knee) | manual post-check | demag study + variable-magnet subroutine | per-element op-point (M-dir) + knee flag | JMAG turnkey / radia can exceed FEMM |
| Hysteresis (Play / Vector-Play) | none | Play/Vector-Play vector hysteresis | none (would need vector-Play operator) | **JMAG-only** |
| EM<->thermal coupling | none | coupled IH (mag<->heat) | solve_thermal + loss heat_source | JMAG turnkey / radia capable |
| Custom physics extension | Lua only | C/C++ user-subroutine DLL | Python CoefficientFunction/integrator (no DLL) | **radia** (most open) / JMAG (compiled) |
| Mesh | internal Triangle (coarse auto) | internal + Nastran .bdf import | netgen (high-order, spline3 arcs) | radia/JMAG |
| Topology optimization | none | turnkey motor topology-opt | NGSolve research code | **JMAG** |
"""

RADIA_CAN_EXCEED = """\
## Where radia-ngsolve can match or exceed FEMM / JMAG

1. **2D magnetostatic accuracy** — order-2 H1 beats 1st-order triangles
   (transverse-magnet cross-check: radia -0.05% vs a stored 2D ref at -0.60%).
2. **Rotating-machine core loss** — per-element FFT over a rotor sweep ->
   Bertotti/Steinmetz/iGSE (the lab magnetic-materials/electromagnet loss
   models) replaces FEMM's .ans hand-edit hack; competitive with JMAG.
3. **Per-element demag knee check** — evaluate the operating point along the
   magnetization direction (B·m̂, H·m̂) and flag irreversible-demag risk;
   FEMM is manual, JMAG needs a subroutine.
4. **EM<->thermal auto-coupling** — feed loss density into solve_thermal
   (convection / radiation BC) in one package; FEMM cannot couple.
5. **Open custom physics** — write custom material / source / loss operators
   directly in Python (no DLL build), consistent with the lab's
   "embrace real-solver modification" policy.
"""

JMAG_ONLY = """\
## JMAG-only capabilities (genuine gaps for the open tools)

These are where JMAG is the differentiator; radia/FEMM cannot easily match:

1. **Play / Vector-Play vector hysteresis** with nested minor-loop families.
2. **Turnkey coupled IH / EM-thermal-structural chains.**
3. **Compiled user-subroutine hooks** spanning magnetization, hysteresis,
   loss, source, circuit, and force.
4. **Native motion + slide air-gap band + N-layer multi-slice skew** with an
   integrated cage circuit.
5. **Built-in motor topology optimization.**

For each, the radia path is "research it" rather than "wire it" — e.g. a
vector-Play hysteron operator as a Python CoefficientFunction is a future
study, not a near-term capability.
"""

FEMM_ROLE = """\
## FEMM's role — the open-source 2D yardstick

FEMM is the open-source, Compumag-standard 2D reference. Even though radia's
order-2 H1 can be more accurate, FEMM's value is as a **widely-validated
implementation** for cross-checking, plus its rich **.fem model asset base**.

- Cross-check headless via **pyFEMM** (`import femm` -> local FEMM install,
  `openfemm(1)`). Repo tests guard with skipif(no-femm); validated results
  live in the lab-internal mcp-server-femm.
- **.fem strategy**: a one-way .fem importer (geometry+materials+BC -> netgen
  re-mesh) is the high-value option to reuse Compumag/FEMM assets. The mesh is
  a non-issue because .fem carries no mesh (FEMM meshes at solve time). A .ans
  parser is redundant (pyFEMM covers it). Scope as import-for-reuse, NOT a
  drop-in replacement (lab policy: no .fem drop-in / no GUI clone). See
  mcp-server-femm `femm_fem_parser_strategy`.
"""

ROADMAP = """\
## radia-motor strengthening roadmap (from cross-learning)

Prioritized capabilities to build into radia-ngsolve, derived from where the
open/commercial tools are strong and radia is reachable. All build on
**already-committed** functions to stay clean of the entangled solve.py WIP
(the cogging and IM tests at 2989bc15 / d33e7ac5 followed this discipline).

1. **Per-element demag knee post-processor** — on top of the existing nonlinear
   PM solve (`test_nonlinear_magnet`), evaluate (B·m̂, H·m̂) per element and
   flag points below the knee. Surpasses FEMM's manual check; mirrors JMAG's
   variable-magnet subroutine. Pure post-processing, no solver change.
2. **Rotor-sweep core loss** — per-element FFT of B over a rotor sweep, then
   Bertotti/iGSE via the lab loss models. periodic_h1 sector (= one pole pitch)
   + air-gap Bn-FFT. Replaces FEMM's .ans hack; competes with JMAG iron-loss.
3. **Ld/Lq frozen-permeability dq** — freeze converged ν(\\|B\\|) as a fixed CF,
   dq current perturbation, `inductance_2d`. Prerequisite for MTPA / T-N.
4. **Full IM 1/4-model** — anti-periodic (`periodic_h1(antiperiodic=True)`,
   validated in test_motor_sector_antiperiodic) + 3-phase distributed winding
   + torque-slip curve / L(ωs) fit -> M, Ll, Rr. Completes the Meeker method;
   the coupling vs slip was already reproduced (d33e7ac5).
5. **EM<->thermal coupling** — feed solve_planar_eddy / core-loss density into
   solve_thermal (convection/radiation BC) for a one-package loss->temperature
   loop. FEMM cannot; JMAG turnkey; radia capable via the lab IH solver.

### Status
- Committed (clean): cogging/reluctance torque (2989bc15), IM rotor-coupling vs
  slip (d33e7ac5).
- Uncommitted handoff (user owns git placement): periodic_h1 in solve.py +
  test_planar_periodic + test_motor_sector_antiperiodic + femm_parity periodic
  rows.
- Next: item 1 (demag knee) is the lowest-risk highest-clarity increment.
"""

SECTIONS = {
    "overview": OVERVIEW,
    "capability_matrix": CAPABILITY_MATRIX,
    "radia_can_exceed": RADIA_CAN_EXCEED,
    "jmag_only": JMAG_ONLY,
    "femm_role": FEMM_ROLE,
    "roadmap": ROADMAP,
}


def get_tritool_cross_reference(topic: str = "overview") -> str:
    """Return tri-tool (FEMM/JMAG/radia-ngsolve) cross-reference knowledge.

    Args:
        topic: section name; "all" returns everything.
    """
    t = topic.strip().lower()
    if t == "all":
        return "\n\n---\n\n".join(SECTIONS[k] for k in SECTIONS)
    if t not in SECTIONS:
        valid = ", ".join(sorted(SECTIONS))
        return f"Unknown topic {t!r}. Valid topics: {valid}\n"
    return SECTIONS[t]
