"""ELF/MAGIC public motor corpus bridge for radia_mcp.motor.

This module is deliberately public-safe: it records only the public MCP
surface and validation labels of ELF-mcp-server, not solver outputs,
private paths, product-run numeric references, or provenance.
"""

from __future__ import annotations


OVERVIEW = """\
## ELF/MAGIC motor bridge for radia-motor

The ELF/MAGIC public MCP package now has a strong motor authoring corpus:

- 652 public motor input cases across 37 `application/motor/...` families.
- 18 motor families expose an explicit observable contract
  (`FLUM`, `HBRM`, `HBCU`, `OHM2` where appropriate).
- 19 additional motor families are proxy-energy validated.
- All 37 motor families are currently `ngsolve_proxy_energy`; there are no
  motor-specific `gold_numeric_invariant` families yet.

So the right reading is:

- **Breadth is strong**: SPM/BLDC, IPM, IM, SRM, SynRM, AFPM, wound-field,
  stepper, linear PM, reluctance, and hysteresis motor patterns are present.
- **Depth is the next training target**: radia-motor should provide stronger
  independent anchors for back-EMF, cogging torque, Ld/Lq, MTPA,
  induction-machine slip loss, and reluctance torque.

This bridge does not execute ELF/MAGIC. It gives radia-motor a shared public
map so that prompt routing, validation planning, and example selection stay
aligned with the ELF public corpus.
"""


COVERAGE_MATRIX = """\
## ELF motor corpus -> radia-motor target matrix

| ELF/MAGIC motor archetype | Public corpus status | radia-motor target |
|---|---:|---|
| SPM / BLDC / PM pickup | broad coverage; 18 observable-contract families include key PM decks | `back_emf`, `cogging_torque`, `ld_lq`, `mtpa` |
| IPM | hairpin, interior-PM, fractional, and static-torque families | `ld_lq`, `mtpa`, `field_weakening`, `demag_margin` |
| Induction machine | cage, rotor-bar, fractional-sector families | `induction_machine`, `deep_bar`, `airgap_eddy_machine` |
| SRM / SR motor | switched-reluctance, 6/4, 8/6, 12/8, 12/16 variants | `reluctance_torque`, `saturating_inductance`, angle-current maps |
| SynRM / reluctance motor | flux-barrier, static-torque, reluctance families | `synchronous_power_angle`, `mtpa`, `cross_saturation` |
| AFPM / axial flux | axial-flux and AFPM linearized-airgap families | `build123d_pmsm_field`, `airgap_machine_rotation` |
| Wound-field / stepper / linear / hysteresis | specialist 10-case families | `power_angle`, `skew_factor`, `hysteresis_motor_loss` |

Use the ELF corpus as runnable input-deck authoring patterns. Use radia-motor
and radia-ngsolve as independent physics anchors that can turn those patterns
from proxy-validated to stronger numeric-invariant stories.
"""


INSUFFICIENCY_AUDIT = """\
## What is still insufficient?

The public motor set is not short of examples; it is short of **gold motor
invariants**. That is a good problem to have.

Current gap:

- `gold_numeric_invariant` motor families: 0.
- Strongest existing motor label: `silver_observable_contract`.
- Existing validation level: `ngsolve_proxy_energy` for all 37 motor families.

Meaning:

- It is fair to say the MCP can route users to many motor input patterns.
- It is not yet fair to claim full numerical agreement for motor torque,
  back-EMF, Ld/Lq, slip torque, or loss without an additional private/local
  solve and public-safe reduction.

Best next upgrades:

1. SPM/IPM: flux linkage -> back-EMF constant, cogging order, Ld/Lq saliency,
   MTPA current-angle checks.
2. IM: slip-frequency rotor eddy loss, Thevenin torque-slip curve,
   deep-bar trend checks.
3. SRM/SynRM: angle-current torque maps, nonlinear inductance, reluctance
   torque sign and periodicity.
4. Hysteresis motor: separate material hysteresis-loop validation from motor
   geometry routing; keep raw material/product references out of public text.
"""


ROUTING_PLAYBOOK = """\
## Prompt routing playbook

When a user asks for a motor simulation:

1. On the ELF side, call `elf_motor_readiness()` to inspect breadth and gaps.
2. Call `elf_sample_decks_route(goal)` and prefer the first family if the
   requested subtype is explicit (`SPM`, `IPM`, `IM`, `SRM`, `SynRM`, etc.).
3. Call `elf_local_simulation_handoff(goal)` when the next step is a
   user-local ELF/MAGIC run.
4. On the radia side, choose the independent anchor:
   - `ngsolve_usage("back_emf")` for PM flux-linkage / EMF.
   - `ngsolve_usage("cogging")` for PM/reluctance torque periodicity.
   - `ngsolve_usage("ld_lq")`, `ngsolve_usage("mtpa")`, and
     `ngsolve_usage("field_weakening")` for IPM/SPM control-layer checks.
   - `ngsolve_usage("induction_machine")` for IM torque-slip behavior.
   - `motor_em_force_recipe("method_choice")` for force/torque extraction.

The two servers should be used as a pair:

- ELF/MAGIC MCP: authoring decks, syntax, product input patterns, public
  validated examples, and local-runner handoff.
- radia-motor / radia-ngsolve: independent physics checks, dq theory,
  field-solver recipes, and public-safe validation reductions.
"""


RADIA_STRENGTHENING_QUEUE = """\
## radia-motor strengthening queue from ELF motor gaps

Near-term additions that would make both servers smarter:

1. Add a compact `motor_validation_router` in radia-motor that maps
   `SPM/IPM/IM/SRM/SynRM` prompt terms to the right radia-ngsolve usage
   snippets and tests.
2. Add one public-safe numeric invariant per motor archetype:
   - SPM: sinusoidal flux-linkage -> back-EMF constant.
   - IPM: `Ld != Lq`, MTPA angle improves torque-per-amp.
   - IM: rotor loss rises with slip; torque-slip peak is finite.
   - SRM/SynRM: torque periodicity and reluctance torque sign.
   - Hysteresis motor: loop-area loss positivity and frequency/material
     separation.
3. Feed those reduced invariants back to ELF-mcp-server as quality-label
   upgrade notes, without publishing solver logs or product-run raw values.

That is the answer to "are we training radia-motor at the same time?":
not automatically from an ELF release, but yes once this bridge and the
corresponding radia tests are updated in the radia-mcp repo.
"""


JMAG_COVERAGE_REALITY = """\
## Can radia-motor cover every JMAG motor capability?

No. radia-motor can cover many **physics anchors**, but it is not a full
JMAG replacement.

Covered or reachable in radia-motor / radia-ngsolve:

- 2D magnetostatic and magnetodynamic A-formulation patterns.
- Air-gap torque, cogging/reluctance torque, and sector/periodic reductions.
- PM flux linkage, back-EMF constants, Ld/Lq, MTPA, and field-weakening
  control-layer checks.
- Induction-machine equivalent-circuit and slip-frequency eddy-current
  anchors.
- SRM/SynRM reluctance-torque and nonlinear inductance trend anchors.
- Research-grade lamination, hysteresis, thermal, and force extraction
  components when the model is kept explicit and testable.

Not fully covered as a turnkey production motor workflow:

- Native production moving-band motion with tightly coupled drive circuits.
- Integrated multi-slice skew plus cage/end-ring circuit workflows.
- Vector-Play / advanced vector hysteresis as a production material model.
- Complete 3D end-region machine workflows with industrial pre/post tooling.
- Built-in motor topology optimization workflows with product-grade templates.
- Turnkey demagnetization, loss, thermal, and structural chains as one GUI
  study type.

The right claim is therefore:

- **radia-motor can be an excellent independent verification and research
  backend for motor physics.**
- **radia-motor should not claim turnkey production-motor coverage.**
"""


AGE_VS_FIELD_STRATEGY = """\
## AGE vs 2D field quick check / ELF-like evaluator

For 2D rotating machines, the main radia-motor path should be NGSolve AGE:

- A-formulation in 2D.
- Sliding / air-gap element style coupling.
- Periodic or anti-periodic sector reduction.
- Direct access to B, A, flux linkage, torque stress, energy, eddy loss,
  and dq perturbation quantities.

That is the natural backbone for SPM/IPM/IM/SRM/SynRM examples because it can
handle conductors, nonlinear reluctivity, and controlled weak forms.

A 2D Radia-style field-kernel / BEM-like evaluator would still be valuable, especially
as an ELF/MAGIC bridge:

- Fast PM and coil flux-linkage sweeps.
- Boundary/source-panel intuition close to ELF/MAGIC-style input decks.
- Lightweight co-energy and pickup-flux regression anchors.
- Good prompt-time authoring checks before a heavier AGE solve.

But a simple 2D field quick check would not automatically cover everything:

- Eddy currents need conductor dynamics or an impedance/operator extension.
- Nonlinear iron needs iteration and a robust material law.
- Hysteresis needs an internal state model.
- IM cage/end-ring behavior needs circuit coupling.
- Motion/skew/end effects still need AGE/multi-slice/3D methods.

So the practical strategy is hybrid:

1. Use **NGSolve AGE** as the authoritative 2D motor solve path.
2. Add a small **2D magnetic-circuit/BEM-like public validation backend** for ELF-like
   PM/coil/reluctance anchors, if we want fast deck-level feedback.
3. Cross-check reduced quantities only: FLUM-like flux linkage, co-energy,
   torque sign/periodicity, back-EMF constants, Ld/Lq saliency, and slip-loss
   trends.
"""


LINEAR_MOTOR_DUAL_LANE = """\
## Linear motor dual-lane replay contract

Linear PM motor examples should train radia-motor as a **dual-lane** target,
not as a one-off deck replay:

1. Use `NGSolve+AGE` for the field/air-gap finite-element lane:
   `linear_thrust`, `linear_pm_flux`, and field-energy consistency.
2. Use `HDiv-MMM + HCurl eddy-bubble` for the independent mixed-system lane:
   `linear_pm_flux`, eddy-current response, and thrust sign/trend checks.
3. Keep the source-tool run as `product_local_reference`; publish only the
   scrubbed lesson, not raw local logs or benchmark rows.
4. A local direct-solver replay may confirm the source deck is runnable, but
   radia-motor learning is accepted only after both radia lanes have verification
   commands and timing metadata attached.

This is the same rule as rotary motors, but linear motors make the force axis
explicit: call it thrust/force rather than torque unless the reduced case is a
rotary surrogate.
"""


ROTARY_MOTOR_FAMILY_SWEEP = """\
## Rotary motor family-sweep replay contract

Rotary motor source-tool examples should train radia-motor as a **family
sweep**, not as one generic motor bucket. A useful 30-case replay covers at
least these five families:

| Family | AGE lane focus | HDiv-MMM + HCurl eddy-bubble focus |
|---|---|---|
| SPM / SPMSM | PM flux linkage, back-EMF, cogging/torque periodicity | pickup flux, demag/source-field trend |
| IPM / hairpin | `Ld/Lq`, MTPA, field-weakening, demag margin | PM plus saliency flux-linkage trend |
| Induction / rotor bar | slip-frequency eddy response and torque-slip trend | source-field and reduced-response trend |
| SRM | reluctance torque sign, angle-current map, saturation | coenergy and force-or-torque trend |
| SynRM | saliency torque, cross-saturation, power-angle checks | coenergy and reduced reluctance trend |
| BLDC / outer-rotor BLDC | PM flux, cogging order, winding/slot polarity | PM source-field and pickup-flux trend |
| Fractional SPMSM | sector periodicity, winding factor, harmonic aliasing | reduced flux-linkage periodicity |
| AFPM / linearized axial flux | unfolded air-gap flux and thrust/torque trend | source-panel flux and demag trend |

The public MCP learning rule is the same as other product-local slots:

1. Record only `product_local_reference` as the source class.
2. Accept radia-motor learning only when `ngsolve_age` and
   `hdiv_mmm_hcurl_eddy_bubble` are both represented by verification commands
   on the same geometry/material/excitation identity.
3. Keep raw product rows, logs, paths, and benchmark numbers in the private
   source lane.
4. Promote only the family coverage, observable names, tolerances, and gate
   commands to public radia-mcp knowledge.
"""


RUN_ARTIFACT_CONTRACT = """\
## Source-tool run artifact contract

For local product runs used as private cross-validation seeds, keep the
artifact roles explicit:

| Suffix | Role | Public-safe use |
|---|---|---|
| `.mai` | analysis/control input deck | route and classify the motor family |
| `.meg` | compiled geometry/mesh input | check that the source-native run is solver-ready |
| `.mei` | mesh-script input, when present | never treat it as the solver result |
| `.mao` | primary execution log | parse version, status, timings, BMAX, FLUM, force/torque rows |
| `.mag` | field/result file | optional field post-processing and probes |
| `.mat`, `.mac`, `.mas` | auxiliary solver outputs | keep private unless reduced to public-safe invariants |

The important correction is that `.mao` is the primary run-log artifact.
Do not count `.mei` as a run result; it is part of the mesh-generation input
route.  When a private source-tool slot is promoted into radia-motor learning,
the public artifact should contain only reduced quantities and generic lessons:
which observable was checked, which radia lane was used, whether AGE and the
HDiv-MMM/HCurl eddy-bubble system both ran, tolerances, version/date/timing
metadata, and a
scrubbed summary.
"""


SECTIONS = {
    "overview": OVERVIEW,
    "coverage_matrix": COVERAGE_MATRIX,
    "insufficiency_audit": INSUFFICIENCY_AUDIT,
    "routing_playbook": ROUTING_PLAYBOOK,
    "radia_strengthening_queue": RADIA_STRENGTHENING_QUEUE,
    "jmag_coverage_reality": JMAG_COVERAGE_REALITY,
    "age_vs_field_strategy": AGE_VS_FIELD_STRATEGY,
    "linear_motor_dual_lane": LINEAR_MOTOR_DUAL_LANE,
    "rotary_motor_family_sweep": ROTARY_MOTOR_FAMILY_SWEEP,
    "run_artifact_contract": RUN_ARTIFACT_CONTRACT,
}


def get_elf_magic_bridge(topic: str = "overview") -> str:
    """Return public-safe ELF/MAGIC motor bridge knowledge for radia-motor."""
    t = topic.strip().lower()
    if t == "all":
        return "\n\n---\n\n".join(SECTIONS[key] for key in SECTIONS)
    if t not in SECTIONS:
        valid = ", ".join(sorted(SECTIONS))
        return f"Unknown topic {t!r}. Valid topics: {valid}\n"
    return SECTIONS[t]
