---
name: panel-cli-diff
description: Check Radia application-block settings against each DesignSpec and calc CLI to catch silently dropped parameters, unknown flags, and defaults that bypass the Simulink configuration. Use with simulink-app-health after changing a block, DesignSpec, or calc_TOPIC.py parser.
---

# Application-to-CLI Diff

Trace every public setting through this chain:

```text
masked Simulink block
  -> versioned application config JSON
  -> DesignSpec
  -> DesignSpec.build_command()
  -> calc_<topic>.py argparse
```

Report settings dropped between layers, emitted flags absent from argparse,
argparse knobs unreachable through the DesignSpec, conflicting defaults, and
implicit solver switches.

Run:

```powershell
python tools/audit_new_panel_contract.py
python -m pytest tests/test_simulink_application.py -q
```

Then run `simulink-app-health`. No application keeps a notebook production
adapter; docs notebooks are checked separately by `ipynb-gui-health`.
