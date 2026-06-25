"""Public motor deck corpus bridge for radia_mcp.motor.

This module is deliberately public-safe: it records only the public MCP
surface and validation labels of an external deck corpus, not solver outputs,
private paths, product-run numeric references, or provenance.
"""

from __future__ import annotations


OVERVIEW = """\
## Public motor deck bridge for radia-motor

The external public deck package now has a strong motor authoring corpus:

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

This bridge does not execute any external solver. It gives radia-motor a
shared public map so that prompt routing, validation planning, and example
selection stay aligned with the public deck corpus.
"""


COVERAGE_MATRIX = """\
## Public motor deck corpus -> radia-motor target matrix

| Motor deck archetype | Public corpus status | radia-motor target |
|---|---:|---|
| SPM / BLDC / PM pickup | broad coverage; 18 observable-contract families include key PM decks | `back_emf`, `cogging_torque`, `ld_lq`, `mtpa` |
| IPM | hairpin, interior-PM, fractional, and static-torque families | `ld_lq`, `mtpa`, `field_weakening`, `demag_margin` |
| Induction machine | cage, rotor-bar, fractional-sector families | `induction_machine`, `deep_bar`, `airgap_eddy_machine` |
| SRM / SR motor | switched-reluctance, 6/4, 8/6, 12/8, 12/16 variants | `reluctance_torque`, `saturating_inductance`, angle-current maps |
| SynRM / reluctance motor | flux-barrier, static-torque, reluctance families | `synchronous_power_angle`, `mtpa`, `cross_saturation` |
| AFPM / axial flux | axial-flux and AFPM linearized-airgap families | `build123d_pmsm_field`, `airgap_machine_rotation` |
| Wound-field / stepper / linear / hysteresis | specialist 10-case families | `power_angle`, `skew_factor`, `hysteresis_motor_loss` |

Use the deck corpus as runnable input-deck authoring patterns. Use radia-motor
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

1. On the deck-corpus side, inspect motor readiness to understand breadth
   and gaps.
2. Route the prompt to public sample decks and prefer the first family if the
   requested subtype is explicit (`SPM`, `IPM`, `IM`, `SRM`, `SynRM`, etc.).
3. Use a local-simulation handoff only when the next step is a user-local
   product run.
4. On the radia side, choose the independent anchor:
   - `ngsolve_usage("back_emf")` for PM flux-linkage / EMF.
   - `ngsolve_usage("cogging")` for PM/reluctance torque periodicity.
   - `ngsolve_usage("ld_lq")`, `ngsolve_usage("mtpa")`, and
     `ngsolve_usage("field_weakening")` for IPM/SPM control-layer checks.
   - `ngsolve_usage("induction_machine")` for IM torque-slip behavior.
   - `motor_em_force_recipe("method_choice")` for force/torque extraction.

The two servers should be used as a pair:

- External deck MCP: authoring decks, syntax, product input patterns, public
  validated examples, and local-runner handoff.
- radia-motor / radia-ngsolve: independent physics checks, dq theory,
  field-solver recipes, and public-safe validation reductions.
"""


RADIA_STRENGTHENING_QUEUE = """\
## radia-motor strengthening queue from public deck gaps

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
3. Feed those reduced invariants back to the deck-corpus server as quality-label
   upgrade notes, without publishing solver logs or product-run raw values.

That is the answer to "are we training radia-motor at the same time?":
not automatically from an external deck release, but yes once this bridge and the
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


AGE_VS_MMM_STRATEGY = """\
## AGE vs 2D MMM / deck-like evaluator

For 2D rotating machines, the main radia-motor path should be NGSolve AGE:

- A-formulation in 2D.
- Sliding / air-gap element style coupling.
- Periodic or anti-periodic sector reduction.
- Direct access to B, A, flux linkage, torque stress, energy, eddy loss,
  and dq perturbation quantities.

That is the natural backbone for SPM/IPM/IM/SRM/SynRM examples because it can
handle conductors, nonlinear reluctivity, and controlled weak forms.

A 2D Radia-style MMM / BEM-like evaluator would still be valuable, especially
as a public deck bridge:

- Fast PM and coil flux-linkage sweeps.
- Boundary/source-panel intuition close to public deck-style inputs.
- Lightweight co-energy and pickup-flux regression anchors.
- Good prompt-time authoring checks before a heavier AGE solve.

But a simple 2D MMM would not automatically cover everything:

- Eddy currents need conductor dynamics or an impedance/operator extension.
- Nonlinear iron needs iteration and a robust material law.
- Hysteresis needs an internal state model.
- IM cage/end-ring behavior needs circuit coupling.
- Motion/skew/end effects still need AGE/multi-slice/3D methods.

So the practical strategy is hybrid:

1. Use **NGSolve AGE** as the authoritative 2D motor solve path.
2. Add a small **2D MMM/BEM-like public validation backend** for deck-like
   PM/coil/reluctance anchors, if we want fast deck-level feedback.
3. Cross-check reduced quantities only: FLUM-like flux linkage, co-energy,
   torque sign/periodicity, back-EMF constants, Ld/Lq saliency, and slip-loss
   trends.
"""


SECTIONS = {
    "overview": OVERVIEW,
    "coverage_matrix": COVERAGE_MATRIX,
    "insufficiency_audit": INSUFFICIENCY_AUDIT,
    "routing_playbook": ROUTING_PLAYBOOK,
    "radia_strengthening_queue": RADIA_STRENGTHENING_QUEUE,
    "jmag_coverage_reality": JMAG_COVERAGE_REALITY,
    "age_vs_mmm_strategy": AGE_VS_MMM_STRATEGY,
}


def get_deck_bridge(topic: str = "overview") -> str:
    """Return public-safe motor deck bridge knowledge for radia-motor."""
    t = topic.strip().lower()
    if t == "all":
        return "\n\n---\n\n".join(SECTIONS[key] for key in SECTIONS)
    if t not in SECTIONS:
        valid = ", ".join(sorted(SECTIONS))
        return f"Unknown topic {t!r}. Valid topics: {valid}\n"
    return SECTIONS[t]
