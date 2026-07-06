---
name: panel-cli-diff
description: Retired desktop-panel CLI diff notes. Current Radia panels are notebook workbenches; use ipynb-gui-health, tools/audit_new_panel_contract.py, and validation_test/panels/test_notebook_workbench.py instead.
---

# panel-cli-diff (retired desktop-panel note)

> **CURRENT POLICY (2026-07-06)**: production Radia analysis panels are
> `src/radia/panels/notebooks/radia_*.ipynb` workbenches backed by
> `DesignSpec` and `CommandWorkbench`. Do not add new `src/radia/radia_*.py`
> desktop panels. Use `python tools/audit_new_panel_contract.py` and
> `python -m pytest validation_test/panels/test_notebook_workbench.py -q`
> for the active panel contract.

Radia's 4-layer architecture (see AGENTS.md § Cubit Panel Architecture)
treats each `calc_*.py` as a pure headless CLI, and each
`radia_*.py::ModePanel` assembles an argv list for that CLI.  Three
classes of silent bug accumulate at that boundary:

1. **Silent drop** — panel has widget `X` whose value never appears in
   any `_build_*_command()`.  User edits the field; nothing reaches
   the calc.
2. **Silent default** — calc has CLI flag `--foo` but no `_build_*_command()`
   emits it.  Calc uses argparse default unconditionally.
3. **argparse reject** — panel emits `--bar` that no calc script
   declares.  Subprocess exits 2 at run time.

This skill finds all three statically, per method in each panel.

## When to run

- AFTER editing `src/radia/radia_*.py`
- AFTER editing any `src/radia/panels/calc_*.py` argparse block
- BEFORE shipping a panel (add to `deploy` L3 layer)

## Output (example)

```
PANEL  radia_ih.py::IHPanel
  method=PEEC+BEM  -> calc_peec_bem.py
    panel->calc OK         --peec-step --peec-nwinc --peec-nhinc
                           --frequency --current --coil-sigma
                           --vol --wp-label --sigma --half-thickness
                           --mu-r --impedance-model --peec-solver
                           --h1-order
    unused widgets:        coil_mu_r (no matching CLI)
    silent defaults:       (none)
  method=FEM A-V   -> calc_fem_coilmesh.py
    panel->calc OK         --vol --frequency --current --coil-sigma
                           --coil-mu-r --sigma --mu-r --half-thickness
                           --fes-order --solver --sibc-bnd --source-bnd
                           --sink-bnd --coil-mat --impedance-model
    silent defaults:       --esim-max-iter --esim-tol   (OK, ESIM-only)
    REJECT (calc rejects):  (none)
  ...
```

Fail the check if any `REJECT` row appears, or any `silent drop /
default` row is not explicitly waived in the panel source via a
`# CLI-DIFF: ignore <flag> -- <why>` comment.

## How it works

### 1. Extract calc_*.py argparse surface

```python
import ast
def cli_args(path: str) -> set[str]:
    tree = ast.parse(open(path).read())
    flags = set()
    for node in ast.walk(tree):
        if (isinstance(node, ast.Call) and
            isinstance(node.func, ast.Attribute) and
            node.func.attr == "add_argument" and
            node.args and isinstance(node.args[0], ast.Constant) and
            isinstance(node.args[0].value, str) and
            node.args[0].value.startswith("--")):
            flags.add(node.args[0].value)
    return flags
```

Note: `add_material_args(parser, ...)` and similar helpers add flags
indirectly.  The scanner follows `from em_material import add_material_args`
and inlines those too — the helpers live in a small whitelist in the
skill (currently `em_material.add_material_args`).

### 2. Extract panel `_build_*_command()` emission set

Parse each `_build_*_command` method and collect every string literal
(or `str(self.val(key))`) passed positionally after a `"--flag"`
element in the `cmd` list.  Conservative: if a flag is only added in
an `if` branch, it is reported as "conditional" and not counted as
always-emitted.

### 3. Extract widget inventory

Every `self.add_line / add_combo / add_spin / add_browse / add_check`
call registers a widget key.  A widget is "used" if its key appears
inside any `self.val("KEY")` (or `build_command` body) that is
reachable for at least one method.

### 4. Diff

- `calc_flags \ panel_emits` → silent-default (may be OK; waived if
  matches a `# CLI-DIFF: ignore` comment)
- `panel_emits \ calc_flags` → REJECT (always a bug)
- `widgets_defined \ widgets_read` → orphan widget

## The checker script

`check_panel_cli.py` (committed to
`tests/panels/check_panel_cli.py`) implements this.  It emits
machine-readable JSON and a human summary.  Exit code 0 = clean,
non-zero = issues found.

```bash
python tests/panels/check_panel_cli.py  # defaults to all radia_*.py
python tests/panels/check_panel_cli.py --panel radia_ih.py  # just IH
python tests/panels/check_panel_cli.py --strict  # fail on any silent-default too
```

## Integration with deploy / notebook workbench tests

- `deploy` skill L3 "Launcher widget matrix" gets an extra line:
  `tests/panels/check_panel_cli.py --strict` must pass before ship.
- `validation_test/panels/test_notebook_workbench.py` checks DesignSpec /
  Workbench wiring and the no-PySide notebook contract; `panel-cli-diff`
  remains useful only where a generated panel/CLI compatibility check still
  exists.

## Waiving a silent-default

If you intentionally let a CLI flag fall back to its argparse default
(common for advanced / debug flags not worth a widget), annotate the
panel source:

```python
# CLI-DIFF: ignore --reg --shift-eps --nthreads -- advanced solver knobs
# CLI-DIFF: ignore --output --msh-output -- auto-generated by panel
```

The checker reads these comments per panel file; flags listed are
silently-defaulted on purpose.
