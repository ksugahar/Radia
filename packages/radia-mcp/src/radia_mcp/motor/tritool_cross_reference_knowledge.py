"""Motor-FEA cross-reference: radia-ngsolve against the open 2D yardstick.

Maps radia-ngsolve motor capabilities against FEMM (the open-source 2D
reference whose parity is the project's stated goal) and lists the
strengthening steps with the highest payoff for radia-motor.
"""

from __future__ import annotations


OVERVIEW = """\
## Motor-FEA cross-reference — radia-ngsolve and FEMM

FEMM is the widely used open-source 2D magnetics code and the yardstick for
radia-ngsolve's 2D motor capability. Cross-referencing serves two purposes:

- Adopt **published, field-tested rules** (cogging period =
  360/LCM(slots,poles), >=4 air-gap mesh layers, weighted-stress-tensor torque).
- Use **radia-ngsolve's high-order accuracy** (order-2 H1) and **EM<->thermal
  coupling** to exceed the 2D yardstick where possible.

Every result is triangulated against an analytic reference where one exists
(e.g. transverse-magnetized cylinder: analytic 0.2513 T, radia -0.05 %).
See also `motor_transient_circuit` (Lange-Henrotte-Hameyer 2009).
"""

CAPABILITY_MATRIX = """\
## Capability matrix

| Capability | FEMM | radia-ngsolve |
|---|---|---|
| 2D magnetostatic (A_z) | planar, 1st-order tri | solve_planar_magnetostatic (order-2 H1) |
| Nonlinear B-H | BH table + N-R | monotone PCHIP B-H law |
| Eddy / time-harmonic | complex A_z | solve_planar_eddy (j omega sigma) |
| Torque / force | block integrals | eggshell_torque_2d |
| Rotation / motion | re-mesh per angle, script sweep | rotor-angle sweep (magdir rotate) |
| IM multi-slice skew + cage | slice method + sinc factor | stacked 2D solves; external cage circuit |
| IM circuit-parameter ID | static L(omega_s) -> M/Ll/Rr fit | solve_planar_eddy slip sweep + anti-periodic 1/4 |
| Ld/Lq (frozen permeability) | manual freeze + re-solve | freeze converged nu -> dq perturb -> inductance_2d |
| Core loss (rotating, B harmonics) | time-harmonic only | per-element FFT -> Bertotti/iGSE |
| Demag (knee) | manual post-check | per-element op-point (M-dir) + knee flag |
| EM<->thermal coupling | none | solve_thermal + loss heat source |
| Custom physics extension | Lua only | Python CoefficientFunction/integrator |
| Mesh | internal Triangle | Netgen (high-order, curved arcs) |
"""

RADIA_CAN_EXCEED = """\
## Where radia-ngsolve can match or exceed the 2D yardstick

1. **2D magnetostatic accuracy** — order-2 H1 beats 1st-order triangles
   (transverse-magnet check: radia -0.05 % against the analytic value).
2. **Rotating-machine core loss** — per-element FFT over a rotor sweep ->
   Bertotti/Steinmetz/iGSE loss models.
3. **Per-element demag knee check** — evaluate the operating point along the
   magnetization direction (B·m̂, H·m̂) and flag irreversible-demag risk.
4. **EM<->thermal coupling** — feed loss density into solve_thermal
   (convection / radiation BC) in one package.
5. **Open custom physics** — custom material / source / loss operators in
   Python, no compiled plug-in build.
"""

OPEN_GAPS = """\
## Open gaps for radia-ngsolve

Research items rather than wiring tasks:

1. Play / Vector-Play vector hysteresis with nested minor-loop families.
2. Native motion with a sliding air-gap band and N-layer multi-slice skew with
   an integrated cage circuit.
3. Motor topology optimization as a turnkey workflow.
"""

FEMM_ROLE = """\
## FEMM's role — the open-source 2D yardstick

FEMM is the open-source 2D reference. Even though radia's order-2 H1 can be
more accurate, FEMM's value is as a **widely-validated implementation** for
cross-checking, plus its rich **.fem model asset base**.

- Cross-check headless via **pyFEMM** (`import femm` -> local FEMM install,
  `openfemm(1)`). Repo tests guard with skipif(no-femm).
- **.fem strategy**: a one-way .fem importer (geometry+materials+BC -> Netgen
  re-mesh) is the high-value option to reuse existing assets. The mesh is a
  non-issue because .fem carries no mesh. Scope it as import-for-reuse, not a
  drop-in replacement or GUI clone.
"""

ROADMAP = """\
## radia-motor strengthening roadmap

1. **Per-element demag knee post-processor** — on top of the nonlinear PM
   solve, evaluate (B·m̂, H·m̂) per element and flag points below the knee.
   Pure post-processing, no solver change.
2. **Rotor-sweep core loss** — per-element FFT of B over a rotor sweep, then
   Bertotti/iGSE. periodic_h1 sector (= one pole pitch) + air-gap Bn-FFT.
3. **Ld/Lq frozen-permeability dq** — freeze converged nu(|B|) as a fixed CF,
   dq current perturbation, `inductance_2d`. Prerequisite for MTPA / T-N.
4. **Full IM 1/4-model** — anti-periodic (`periodic_h1(antiperiodic=True)`)
   + 3-phase distributed winding + torque-slip curve / L(omega_s) fit -> M, Ll,
   Rr.
5. **EM<->thermal coupling** — feed solve_planar_eddy / core-loss density into
   solve_thermal for a one-package loss->temperature loop.
"""

SECTIONS = {
    "overview": OVERVIEW,
    "capability_matrix": CAPABILITY_MATRIX,
    "radia_can_exceed": RADIA_CAN_EXCEED,
    "open_gaps": OPEN_GAPS,
    "femm_role": FEMM_ROLE,
    "roadmap": ROADMAP,
}


def get_tritool_cross_reference(topic: str = "overview") -> str:
    """Return motor-FEA cross-reference knowledge.

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
