# Optimization examples

Optional optimization examples for CAE-style objectives.  These examples may
use external optimizers such as Optuna, but radia-mcp itself does not vendor or
serve Optuna; install the official public `optuna` / `optuna-mcp` packages
separately when needed.

| Example | Shows | Capabilities used |
|---|---|---|
| [`validation_optuna_waveguide_slab.py`](validation_optuna_waveguide_slab.py) | Validation-class Optuna study for a waveguide dielectric slab: enqueue analytic half-wave candidates, minimize `|S11|`, record feasible trials | `optuna`, `waveguide_dielectric_slab_sparams`, `best_feasible_record`, `constraint_violation` |

```powershell
python validation_optuna_waveguide_slab.py
```

The public example contains only open, analytic waveguide formulas and optional
plain Optuna usage. Study/trial/dashboard MCP operation belongs to the official
external `optuna/optuna-mcp` server.
