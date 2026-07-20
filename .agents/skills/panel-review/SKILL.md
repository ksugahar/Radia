---
name: panel-review
description: Review Radia application interfaces for Simulink-block wiring, DesignSpec-to-calc parity, explicit-trigger execution, durable result artifacts, and the temporary IH notebook comparison. Use simulink-app-health as the production gate.
---

# Application Interface Review

Review the masked block, `radia.simulink.application` runner, DesignSpec,
headless `calc_*.py`, and tests as one contract. Prioritize silent parameter
drops, per-step Python calls, missing failure artifacts, and mismatched result
keys. MEX is optional and must not be assumed production-ready.

Run `simulink-app-health`. For IH, also run `ipynb-gui-health` and confirm both
interfaces delegate to the same headless computation.
