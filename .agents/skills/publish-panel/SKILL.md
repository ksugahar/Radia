---
name: publish-panel
description: Release checklist for a new Radia GUI panel mode (e.g. adding a Method option to IHPanel, or a new radia_*.py panel entirely). Enforces the 6 requirements -- UI readable on real desktop Qt, launcher surfaces the mode, sample jou + STEP/vol exists, MCP docs updated, QA tests pass, deploy verified -- before the mode is considered shipped. Use AFTER implementing the Python and BEFORE saying "done".
---

# publish-panel

Do NOT claim a new panel mode is "shipped" until ALL six items below
pass.  The PEEC-inductance regression of 2026-04-20 (section headers
vertically clipped on real desktop, no launcher path, no sample,
MCP doc missing one alias) was caused by skipping four of these steps.

## 1. UI readable on REAL desktop Qt (not offscreen)

Offscreen Qt uses a different font metric path.  Rows that render
correctly offscreen can be vertically clipped on the real desktop.
Always grab a screenshot via a real Qt platform and read it:

```powershell
# From repo root
$env:QT_QPA_PLATFORM = $null  # force native platform, not offscreen
python -c "
import sys
sys.path.insert(0, 'src/radia'); sys.path.insert(0, 'src/radia/panels')
from PySide6.QtWidgets import QApplication, QLabel
app = QApplication([])
import radia_ih, time
for tag, method in [('ind', radia_ih.METHOD_PEEC_IND),
                    ('bem', radia_ih.METHOD_PEEC_BEM),
                    ('fem', radia_ih.METHOD_FEM_FULL)]:
    w = radia_ih.IHWindow('')
    w._panel._method_combo.setCurrentText(method)
    w.show()
    for _ in range(10): app.processEvents()
    time.sleep(0.3); app.processEvents()
    w.grab().save(f'temp/panel_real_{tag}.png')
    # Also measure every section header height
    for lbl in w.findChildren(QLabel):
        if lbl.text().startswith('<b>'):
            h = lbl.size().height()
            fm_h = lbl.fontMetrics().boundingRect('Mg').height()
            ok = 'OK' if h >= fm_h + 10 else 'CLIPPED'
            print(f'  {tag} {lbl.text()!r:<35} h={h} need>={fm_h+10} [{ok}]')
    w.close()
"
```

Then open each `temp/panel_real_*.png` with the `Read` tool and verify:
- Every section header (`<b>Title</b>`) is fully visible (not clipped)
- No section header overlaps with the next row's content
- Run / Stop / Open GMSH buttons are inside the window
- Combo box current value is fully visible
- No orphan section headers (title with no content rows below)

**Common bug**: `margin-top: Npx` in QLabel stylesheet is painted INSIDE
the widget rect but NOT added to `sizeHint()`.  The QFormLayout then
gives the row a too-small height and the bold text is clipped.  Use
`setFixedHeight(fontMetrics + 10)` instead.

## 2. Launcher (.ccl) can reach the new mode

```
Cubit [Solve menu] -> [C++ Qt5 .ccl launcher dialog]
                       -> [IHWindow / EMWindow / ...]
                          -> [Method combo picks the new mode]
```

The `.ccl` launcher (`src/cubit_plugin/RadiaComp.cpp`) gates panel
launch on "Analysis" type + required labels + `.vol` path.  A new
mode is **user-reachable** only when:

- [ ] The `.ccl` Analysis combo includes the mode (or a parent mode
      that routes to it) -- AND the label / .vol requirements match
      what the mode actually needs (e.g. PEEC-inductance does NOT
      need a .vol; the launcher must allow STEP-only entry).
- [ ] If the new mode is a sub-option selected inside the IHWindow
      Method combo, the parent Analysis launches IHWindow cleanly
      and the Method combo is visible immediately on first launch.

**Every launcher widget must be mode-appropriate** (2026-04-21 miss:
"Mesh order: 2" stayed visible for PEEC-inductance which has no
mesh).  For STEP-only modes (`NEEDS_VOL = False`):

| Widget           | Show when needsVol=True | Show when needsVol=False |
|------------------|-------------------------|--------------------------|
| Analysis combo   | always                  | always                   |
| Mesh order combo | yes (.vol export order) | **NO** (no mesh export)  |
| Labels group     | yes (block/sideset req) | **NO** (no labels used)  |
| Output .vol path | yes                     | **NO** (no .vol)         |
| OK / Cancel      | always                  | always                   |

Every new mode-specific widget in the launcher **must** have a
`setVisible(ms.needsVol)` call inside the `updateLabels` lambda.
The deploy skill's L3 static audit (`7d2`) parses
`RadiaComp.cpp` and warns on any `QWidget* *Row` / `*Group` that
lacks such a call.

**Verification**: launch Cubit manually, click through the launcher,
confirm the new mode is selectable without confusing error messages
about "missing .vol" when .vol is not needed.  Take a screenshot of
the launcher dialog in each mode and verify only the relevant
widgets are visible.

## 3. Sample in `panels/samples/`

Every panel mode needs a working minimal example that a new user can
run in < 5 minutes:

- [ ] `panels/samples/{panel}_{mode}.jou`      -- Cubit script that
      produces the input file(s) the mode needs
- [ ] `panels/samples/{panel}_{mode}_coil.step` (or `.vol`, etc.) --
      pre-built output of the .jou so the user can skip Cubit
- [ ] Header comment in the .jou lists: required labels / what the
      method does / expected output values / runtime estimate

**Naming convention**: `ih_bem_sample.jou`, `ih_fem_kelvin_sample.jou`,
`ih_peec_inductance.jou`, `em_sample.jou`, `pcb_sample.jou`.

**Verification**: `pip install radia[cubit]` -> `cubit-plugin-install`
-> Cubit [Solve] -> pick the new mode -> Browse to the sample ->
Run.  Expected output matches the header comment.

## 4. MCP documentation

Both mcp-server-cubit and mcp-server-radia-ngsolve must gain a topic
(or extend an existing one) that documents:

- [ ] What the new mode is for (physics / use case)
- [ ] Required inputs (STEP / .vol / labels)
- [ ] Expected outputs (result schema)
- [ ] Caveats / limits (e.g. thin-skin only, 2D axisym only)
- [ ] Topic aliases so `cubit_docs(topic="helix")` and
      `ngsolve_usage(topic="peec-inductance")` etc. all resolve

**Verification**:
```bash
python -c "
import sys
sys.path.insert(0, 'packages/radia-mcp/src')
sys.path.insert(0, 'src/radia')
from radia_mcp.cubit.cubit_scripting_knowledge import get_cubit_documentation
from mcp_server.radia_ngsolve.ngsolve_knowledge import get_ngsolve_documentation
# Expect matching content for your topic:
print(get_cubit_documentation('YOUR_ALIAS')[:200])
print(get_ngsolve_documentation('YOUR_ALIAS')[:200])
"
```

## 5. `tests/panels/panel_qa.py` registers the new panel mode

- [ ] Add `(tag, WindowClass, combo_attr_or_key, value)` entry to
      `get_panel_registry()` in `tests/panels/panel_qa.py`
- [ ] `pytest tests/panels/test_panel_qa.py -v` passes on all 7 checks
- [ ] `temp/panel_{tag}.png` generated and looks right

If the mode's UI differs structurally (e.g. hides a whole section),
the section-header collapse + `_set_row_visible` logic must keep
no-orphan-section-headers passing.

## 6. Deploy verified on 100号機 (golden-range, not just "exit 0")

SMB `cp` is **forbidden** (see deploy skill "SMB cp は禁止").  Deploy
must go through the wheel-first pip install path so what 100号機
receives is what PyPI would give a public user.

- [ ] Build wheel on LAB (`Build_Wheel.ps1 -DryRun`) and confirm
      `dist/radia-X.Y.Z-cp312-cp312-win_amd64.whl` exists
- [ ] On 100号機: force-close all Cubit (see deploy Stage 2(a)) →
      `pip install --no-deps --force-reinstall W:\...\dist\*.whl`
- [ ] `pytest tests/panels/test_panel_qa.py -v` on 100号機 PASS
- [ ] **Register the new mode in the deploy skill's Panel Mode Matrix**
      (section "Panel mode matrix" + Panel Mode matrix table in
      `deploy/SKILL.md`): add a row with
      `panel | calc script | sample input | golden output range`
- [ ] **Add a golden-range assertion** in the Sample E2E runner for
      the new mode.  The range must distinguish the correct physics
      from a silently-wrong fallback (e.g. L off by 10x from a
      wrong cross-section area).  `status: "ok"` from the subprocess
      is NOT evidence; a numeric range check IS.
- [ ] Run the full Sample E2E runner on 100号機 and confirm the new
      row PASSes its golden-range check
- [ ] Cubit [Solve] on 100号機 can find the new mode in the launcher
      (requires .ccl rebuild + deploy if launcher changed) — manual
      click-through, NOT just `Test-Path` on .ccl bytes

Refer to `.Codex/skills/deploy/SKILL.md` "Panel mode matrix" and
"Sample E2E runner" sections for the golden-range test pattern.

## Anti-checklist: things that are NOT enough

These each *feel* like "done" but don't count:

- ❌ Offscreen Qt shows the panel correctly
      (real desktop renders differently -- see #1)
- ❌ Python headless test passes
      (UI can still be unreadable even if `_widgets` are present)
- ❌ `calc_*.py --help` works
      (users never call it this way)
- ❌ MCP knows the name but not the aliases users will try
      (check `boolean / glue / fuse / trampoline` style aliases)
- ❌ It worked on LAB but was not deployed to 100号機
      (other lab members can't see it)
- ❌ Deploy to 100号機 succeeded but nobody tried Cubit launcher
      (launcher blocks entry to the mode)
- ❌ **`calc_*.py` exits 0 with `status: "ok"` JSON but the numeric
      output was not range-checked** (2026-04-21 Kubota 4.78 nH
      regression: L was 100x wrong, panel "worked", user got a
      useless number and lost trust in the panel)
- ❌ **Sample validated only on a toy geometry** (single-loop torus)
      but shipped without a multi-turn / realistic-complexity test
      (the failure mode is almost always in the complex case)

## See also

- `.Codex/skills/deploy/SKILL.md` Step 8 -- panel visual-render QA
  (run as part of every deploy, catches partial regressions)
- `.Codex/skills/panel-qt-test/SKILL.md` -- headless panel test layer
- `.Codex/skills/verify-deploy/SKILL.md` -- signature-level verify
  that 100号機 Python actually loads your edits
- mcp-server-radia-ngsolve `panel_gui_pitfalls` -- common pitfalls
  catalogue (combo_state, mode_switch, gmsh_viz, subprocess_args,
  sample_jou, learn_edition_cap, ...)
