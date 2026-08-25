"""
Radia install / deploy policy and recipes — 2-tier configuration
(LAB + 100号機 editable / mdx + hibino PyPI), reversible migration steps,
and the non-obvious gotchas that cause silent breakage.

Read this when:
* Setting up a new lab machine.
* Migrating an existing machine between editable / PyPI install.
* Diagnosing "import works but pip says wrong version", "DLL load failed",
  or "panel subprocess dies with ModuleNotFoundError: PySide6" on a
  freshly-deployed machine.

The MCP server exposes this via install_deploy(topic=...). Topics:
two_tier, lab_editable, hyaku_editable, mdx_pypi, hibino_pypi, gui_extra,
editable_to_pypi_migration, pypi_to_editable_migration,
metadata_sync, pyd_dll_bootstrap, cubit_plugin_layers,
first_install_panel_registration, common_failure_modes
"""

INSTALL_DEPLOY = """\
# Radia install / deploy policy (2026-05-02 simplified to 2-tier)

This document is the canonical record of how the Radia deployment
machines (LAB, 100号機, mdx, hibino) are installed and kept in sync. It
supersedes:
* the 2026-04-24 policy (LAB+100号機 = NAS editable, mdx = PyPI),
* the 2026-05-01 3-tier policy (LAB editable, mdx editable, 100号機 PyPI).

The current 2-tier policy is:
* LAB + 100号機 are NAS editable installs.
* mdx + hibino are PyPI wheel consumers.

The user set the 100号機 editable / mdx PyPI split explicitly, and on
2026-06-25 added hibino as another PyPI consumer reached by `ssh hibino`.
The release gate is QUAD: four machines must agree at Phase 9.

Topics: two_tier, lab_editable, hyaku_editable, mdx_pypi, hibino_pypi, gui_extra,
        editable_to_pypi_migration, pypi_to_editable_migration,
        metadata_sync, pyd_dll_bootstrap, cubit_plugin_layers,
        first_install_panel_registration, common_failure_modes

============================================================
## two_tier — current 2-tier configuration (2026-05-02)
============================================================

| Machine  | Install                                  | Source                                            |
|----------|------------------------------------------|---------------------------------------------------|
| LAB      | editable (`pip install -e`)              | NAS `public-safe curated corpus`                        |
| 100号機  | editable (`pip install -e`)              | managed internal share |
| mdx      | PyPI wheel consumer                      | `C:\\Program Files\\Python312\\Lib\\site-packages` |
| hibino   | PyPI wheel consumer via `ssh hibino`     | `C:\\Program Files\\Python312\\Lib\\site-packages` |

Roles:

* **LAB** — dev iteration; Build.ps1 outputs are immediately live via
  the editable install.
* **100号機** — the lab's shared 21-user box. It reads the NAS source
  editably so LAB-side fixes can be reflected by metadata refresh plus
  plugin deploy, not by waiting for PyPI.
* **mdx** — PyPI verification point, isolated from NAS, used as a wheel
  consumer in `release_quad.py phase8e` and the Phase 9 drift probe.
* **hibino** — PyPI verification / deployment point, reached by
  `ssh hibino`. It installs `radia`, `cubit-mesh-export`, and
  `radia-mcp` from PyPI and is included in `release_quad.py phase9`.
  Cubit is optional on hibino: if Coreform Cubit 2025.12+ is absent,
  `release_quad.py phase8 --target hibino` skips the plugin/smoke part
  after the PyPI install. Its bare `python` command is a Windows Store
  alias, so release probes and deploys use `py -3.12`.

The 2026-05-01 3-tier policy (mdx editable from a local clone) was
retired 2026-05-02 because mdx editable required (a) a base64-over-ssh
push of `.pyd` artifacts (`tools/push_pyds_to_mdx.py`) since mdx has
no `gh` CLI, (b) manual deletion of legacy site-packages shadows, and
(c) `cubit-plugin-install --all-users` had to run anyway to deploy
the Cubit plugin — net zero gain over PyPI install.

============================================================
## lab_editable — LAB editable install
============================================================

LAB has the canonical NAS-mapped working tree at `public-safe curated corpus`.
All editing happens here. Build.ps1 compiles C++ artifacts directly
into the local clone (`src/radia/_radia_pybind.pyd`,
`packages/cubit-mesh-export/.../cubit_mesh_curver.pyd`, etc.) and the
editable install picks them up immediately.

After every release, refresh the LAB metadata:

```powershell
pip install -e . --no-deps --no-cache-dir
pip install -e packages/cubit-mesh-export --no-deps --no-cache-dir
pip install -e packages/radia-mcp --no-deps --no-cache-dir
```

This is the **metadata_sync** step (see topic). Don't skip it after a
release — `importlib.metadata.version("radia")` will lag.

LAB ↔ 100号機 NAS-share is preserved (filesystem same), and 100号機 is
also editable from that share. mdx and hibino are wheel consumers and
do not read the NAS source during the release gate.

============================================================
## hyaku_editable — 100号機 NAS editable install
============================================================

100号機 (the lab's shared workstation) receives Radia from the managed
internal share via editable install. The SSH host is configured outside
the public package.  Use the workstation's mapped drive form
configured on that machine for the editable path, not a self-referential
UNC path back to the same workstation: loading `_radia_pybind.pyd` directly
from that UNC form can access-violate on Windows even though the same binary
imports from a local path or PyPI wheel.

After every release or source-side deploy, `release_quad.py phase8
--target 100` runs the equivalent of:

```powershell
cat << 'PS' | ssh <shared-workstation-host> 'pwsh -ExecutionPolicy Bypass -Command -'
# Stop locks
Get-Process -ErrorAction SilentlyContinue | Where-Object {
    $_.Name -like 'mcp-server*' -or
    $_.ProcessName -in 'coreform_cubit','cubit'
} | ForEach-Object { Stop-Process -Id $_.Id -Force }
Start-Sleep -Seconds 2

pip install -e "<target-mapped-repo>" `
    -e "<target-mapped-repo>\\packages\\cubit-mesh-export" `
    -e "<target-mapped-repo>\\packages\\radia-mcp" `
    --no-deps --no-cache-dir

cubit-plugin-install --all-users
cubit-plugin-install --verify-only --all-users
cubit-smoke-test
PS
```

The `--all-users` flag writes the Cubit plugin to
`C:\\Program Files\\Coreform Cubit 2025.12\\bin\\plugins\\` and registers
the Radia toolbar startup shim for all existing profiles plus the Default
profile (admin privilege required).  `--verify-only --all-users` checks
both binary hashes and panel startup registration.

============================================================
## mdx_pypi — mdx PyPI install
============================================================

mdx is a PyPI verification point and intentionally does not read NAS
source during the release gate. It installs `radia` and
`cubit-mesh-export`; `radia-mcp` is intentionally omitted on mdx and
appears as `N/A` in the Phase 9 drift table.

```powershell
# Stop locks (python / mcp-server / Cubit hold .pyd open)
Get-Process -ErrorAction SilentlyContinue | Where-Object {
    $_.Name -like 'mcp-server*' -or
    $_.ProcessName -in 'coreform_cubit','cubit'
} | ForEach-Object { Stop-Process -Id $_.Id -Force }
Start-Sleep -Seconds 2

pip install --upgrade --no-cache-dir \\
    'radia[cubit]==<X.Y.Z>' \\
    'cubit-mesh-export==<X.Y.Z>'

cubit-plugin-install --all-users
cubit-plugin-install --verify-only --all-users
```

mdx HAS Cubit installed at
`C:\\Program Files\\Coreform Cubit 2025.12` (Coreform Cubit Learn
Edition).  `cubit-plugin-install --verify-only --all-users` validates
plugin binary sha256 hashes and Radia panel startup registration.
Headless smoke test: the `cubit-smoke-test` CLI (shipped with
cubit-mesh-export >= 0.5.3) runs `coreform_cubit -batch -nographics`
and verifies `export netgen` writes a valid .vol.

The retired editable workflow at `C:\\Radia\\01_GitHub` may stay for
fixture reads (`tests/panels/golden/*.step` etc.) but Python no
longer imports from it.  `tools/push_pyds_to_mdx.py` is preserved
for branch-test scenarios only (testing a pre-release `.pyd` without
cutting a tag); it's NOT part of the standard release path anymore.

============================================================
## hibino_pypi — hibino PyPI install
============================================================

hibino is the second PyPI consumer machine in the release-quad gate. Use
the SSH alias `hibino`; do not pass Japanese paths or long PowerShell
commands as SSH arguments. Use `py -3.12` for Python on hibino; the
bare `python` command is the Windows Store alias and will only print
`Python`.  `C:\\Python312\\Scripts` may be absent from PATH, so
`release_quad.py` derives the Scripts directory from `py -3.12 -c
"import sys; print(sys.executable)"` before calling
`cubit-plugin-install` or `cubit-smoke-test`.  If Coreform Cubit
2025.12+ is not installed, hibino remains a package/MCP verification
target and the Cubit plugin/smoke step is skipped.

```powershell
cat << 'PS' | ssh hibino 'pwsh -ExecutionPolicy Bypass -Command -'
Get-Process -ErrorAction SilentlyContinue | Where-Object {
    $_.Name -like 'mcp-server*' -or
    $_.ProcessName -in 'coreform_cubit','cubit'
} | ForEach-Object { Stop-Process -Id $_.Id -Force }
Start-Sleep -Seconds 2

py -3.12 -m pip install --upgrade --force-reinstall --no-deps --no-cache-dir `
    'radia[cubit]==<X.Y.Z>' `
    'radia-mcp==<X.Y.Z>' `
    'cubit-mesh-export==<X.Y.Z>'

cubit-plugin-install --all-users
cubit-plugin-install --verify-only --all-users
cubit-smoke-test
PS
```

`tools/release_quad.py phase8 --target hibino` runs this recipe through
the common PyPI deploy helper. `tools/release_quad.py phase9` includes
hibino in the version/hash drift table.

============================================================
## gui_extra — retired after Simulink application migration
============================================================

The old `gui` extra was used by standalone PySide6 desktop panels.  It is no
longer part of the production deploy recipe. The canonical human surface is
the Radia Simulink library; normal Radia Python on LAB / 100号機 / mdx / hibino
should not install PySide6. Notebook workbenches are retired for all applications.

Production install:

```powershell
pip install --upgrade 'radia[cubit]==<X.Y.Z>' 'radia-mcp==<X.Y.Z>' 'cubit-mesh-export==<X.Y.Z>'
cubit-plugin-install --all-users
cubit-plugin-install --verify-only --all-users
```

Application interface health gates:

```powershell
python -m pytest tests/test_simulink_application.py -q
python -m pytest tests/test_application_interface_manifest.py -q
```

Boundary check:

```powershell
python -c "import importlib.util; print('PySide6', 'FOUND' if importlib.util.find_spec('PySide6') else 'MISS')"
```

`MISS` is the expected result in normal Radia Python.  Do not uninstall or
delete Coreform Cubit's bundled PySide6 under its own installation directory;
Cubit owns that embedded runtime.

============================================================
## editable_to_pypi_migration — switch a machine from editable to PyPI
============================================================

Use this when moving a machine from editable to PyPI, e.g. mdx or
hibino. 100号機 is intentionally editable in the current release-quad
policy.

```powershell
# Stop locks (mcp-server, Cubit, jupyter, python)
Get-Process -ErrorAction SilentlyContinue | Where-Object {
    $_.Name -in 'python','pythonw','jupyter','jupyter-lab','ipykernel' -or
    $_.Name -like 'mcp-server*' -or
    $_.ProcessName -in 'coreform_cubit','cubit'
} | ForEach-Object { Stop-Process -Id $_.Id -Force -ErrorAction SilentlyContinue }
Start-Sleep -Seconds 3

# Uninstall editable (frees the .pth file, keeps clone intact)
pip uninstall -y radia radia-mcp cubit-mesh-export

# Delete legacy single-file shadow if present (cubit_mesh_export.py)
$shadow = "C:\\Program Files\\Python312\\Lib\\site-packages\\cubit_mesh_export.py"
$shadowPyc = "C:\\Program Files\\Python312\\Lib\\site-packages\\__pycache__\\cubit_mesh_export.cpython-312.pyc"
if (Test-Path $shadow) { Remove-Item $shadow -Force }
if (Test-Path $shadowPyc) { Remove-Item $shadowPyc -Force }

# Install from PyPI. The Simulink/headless route intentionally avoids PySide6 in
# normal Radia Python; see the gui_extra topic.
pip install --no-cache-dir \\
    'radia[cubit]==<X.Y.Z>' \\
    'radia-mcp==<X.Y.Z>' \\
    'cubit-mesh-export==<X.Y.Z>'

# Cubit plugin: regular-file deploy from PyPI wheel
cubit-plugin-install --all-users
cubit-plugin-install --verify-only --all-users
# sha256 plus Radia toolbar startup registration must match
```

Sanity:

```python
import importlib.metadata as md
import radia
print(radia.__file__)
# Expected: C:\\Program Files\\Python312\\Lib\\site-packages\\radia\\__init__.py
# (NOT a NAS or local clone path)

import json
info = json.loads(md.distribution("radia").read_text("direct_url.json") or "{}")
assert not info.get("dir_info", {}).get("editable"), "still editable!"
```

The 2026-05-01 100号機 migration and 2026-05-02 mdx migration both
completed cleanly with this recipe (5 sha256 [OK]'s on
cubit-plugin-install --verify-only).

============================================================
## pypi_to_editable_migration — switch a machine from PyPI to editable
============================================================

Reverse direction. Useful for putting a development branch on a
specific machine for testing.  Note: this is NOT the standard
workflow as of 2026-05-02 (only LAB is editable in production).

```powershell
# Stop locks
Get-Process -ErrorAction SilentlyContinue | Where-Object {
    $_.Name -in 'python','pythonw','jupyter','jupyter-lab','ipykernel' -or
    $_.Name -like 'mcp-server*'
} | ForEach-Object { Stop-Process -Id $_.Id -Force -ErrorAction SilentlyContinue }
Start-Sleep -Seconds 3

# Make sure local clone exists
$dest = 'C:\\Radia\\01_GitHub'
if (-not (Test-Path $dest)) {
    git clone https://github.com/ksugahar/Radia.git $dest
}
Set-Location $dest
git fetch origin --tags
git checkout v<X.Y.Z>

# Uninstall PyPI version + install editable
pip uninstall -y radia radia-mcp cubit-mesh-export
pip install --no-deps --no-cache-dir -e .
pip install --no-deps --no-cache-dir -e packages/cubit-mesh-export
pip install --no-deps --no-cache-dir -e packages/radia-mcp
```

If the target is mdx (no NAS, no MSVC), then run
`tools/push_pyds_to_mdx.py` from LAB to drop the C++ binaries.

============================================================
## metadata_sync — pip metadata vs `radia.__version__`
============================================================

Editable installs have a separate metadata record (the
`<package>-<version>.dist-info` directory). Without explicit refresh,
`importlib.metadata.version("radia")` and `pip list` can lag behind
`radia.__version__` (which is read from the local clone's
`src/radia/__init__.py`). After every release on an editable
machine (in 2-tier this is LAB only):

```powershell
pip install -e . --no-deps --no-cache-dir
pip install -e packages/cubit-mesh-export --no-deps --no-cache-dir
pip install -e packages/radia-mcp --no-deps --no-cache-dir
```

`tools/release_quad.py done` runs a Phase 9 cross-machine probe
that compares `__version__` and `pip list` outputs. A DRIFT row in
that report usually means metadata_sync was skipped on LAB (the
only editable machine left).

For mdx / hibino (PyPI install), the metadata is synced
automatically by `pip install --upgrade radia==<X.Y.Z>` because
PyPI installs are not editable.
mdx intentionally lacks `radia-mcp`; hibino installs it.

============================================================
## pyd_dll_bootstrap — `cubit_mesh_curver.pyd` requires `import radia` first
============================================================

`cubit_mesh_curver.pyd` (the cubit-mesh-export C++ extension) depends
on Netgen DLLs (`nglib.dll`, `ngcore.dll`). On import these DLLs must
be findable via the Windows DLL search path.

`radia/__init__.py` calls `os.add_dll_directory(...)` for the
NGSolve / Netgen / MKL DLL directories. So:

```python
import radia                    # First — sets DLL paths
from cubit_mesh_export import cubit_mesh_curver   # Now this loads
```

Without the `import radia` first, you get
`ImportError: DLL load failed while importing cubit_mesh_curver: ...`.

This issue applied historically to mdx editable; on a PyPI install
(mdx and hibino now) the wheel installer registers DLL search
paths via its own hooks, so the bootstrap is automatic.

============================================================
## cubit_plugin_layers — Cubit plugin lives in TWO independent places
============================================================

The Cubit plugin (`cubit_mesh_export.ccm`,
`cubit_mesh_curver.cp312-win_amd64.pyd`, `nglib.dll`, `ngcore.dll`) is
deployed to `C:\\Program Files\\Coreform Cubit 2025.12\\bin\\plugins\\`
by `cubit-plugin-install`.  The old Qt5 `.ccl` path was retired in
radia 4.80.0; any remaining `cubit_mesh_export.ccl` or `radia_cubit.ccl`
is stale and `cubit-plugin-install --verify-only` reports it as a failure.
This is INDEPENDENT of the Python editable install: Cubit reads its
plugin directory directly, not via Python's import system.

Two layers:

| Layer | Location | Used by |
|-------|----------|---------|
| Python import | `<package>/cubit_mesh_export/` | `from cubit_mesh_export import ...` |
| Cubit plugin  | `Coreform Cubit 2025.12\\bin\\plugins\\` | Cubit GUI / `coreform_cubit -batch` |

Implications:

* On mdx / hibino (PyPI install), `pip install --upgrade
  cubit-mesh-export==X.Y.Z` updates the Python side. You MUST then
  run `cubit-plugin-install --all-users` to update the Cubit side and
  panel startup registration.
* On LAB (editable install), Build.ps1 writes the binaries to BOTH
  the local clone AND the Cubit\\bin\\plugins location (via the
  Build.ps1 post-build step). No separate cubit-plugin-install
  needed in routine dev.

============================================================
## first_install_panel_registration — why first install used to miss panels
============================================================

The Cubit plugin and the Radia Solve menu are separate layers:

1. Cubit loads `cubit_mesh_export.ccm` from
   `Coreform Cubit 2025.12\\bin\\plugins\\`.
2. Cubit replays each user's `.cubit` startup file.
3. The startup file plays a generated Radia shim:
   * current-user install:
     `%LOCALAPPDATA%\\Radia\\Cubit\\radia_startup.py`
   * `--all-users` install:
     `%ProgramData%\\Radia\\Cubit\\radia_startup.py`
4. The shim loads `site-packages\\radia\\panels\\register_toolbar.py`,
   which creates the PySide6 Solve menu entries.

Historical failure: the installer wrote `startup.py` inside the Python
package tree and `--verify-only` checked only plugin hashes.  A deploy
could therefore be "green" while the user's `.cubit` file still pointed
at an old or missing startup script.  Current policy:

```powershell
cubit-plugin-install --all-users
cubit-plugin-install --verify-only --all-users
```

`--verify-only --all-users` fails if any targeted `.cubit` startup block,
generated startup shim, `register_toolbar.py` path, or Cubit.ini
`plugin\\paths` entry is missing.  Do not report a first install as done
until this passes.

============================================================
## common_failure_modes — symptoms and fixes
============================================================

| Symptom                                    | Root cause                                                                | Fix |
|-------------------------------------------|---------------------------------------------------------------------------|-----|
| Simulink application block does not produce `result.json` | config / DesignSpec / Python environment / calc argv drift | Run `pytest tests/test_simulink_application.py -q` and MATLAB Simulink tests; inspect `run.log`; do not install PySide |
| `radia.__version__` says X.Y.Z but `pip list` says A.B.C | metadata_sync skipped on LAB (editable)                                | `pip install -e <path> --no-deps --no-cache-dir` |
| `ImportError: cannot import name 'check' from 'cubit_mesh_export'` | legacy `cubit_mesh_export.py` shadow at site-packages                | Delete `Lib/site-packages/cubit_mesh_export.py` and `__pycache__/cubit_mesh_export.cpython-312.pyc`, then re-install |
| `AttributeError: module 'cubit_mesh_export' has no attribute '__version__'` | same legacy shadow                                                     | same |
| `ImportError: DLL load failed while importing cubit_mesh_curver` | `import radia` not done first (DLL bootstrap) — only on legacy editable | `import radia` before `from cubit_mesh_export import cubit_mesh_curver`; on PyPI installs this is automatic |
| `cubit-smoke-test` fails with "Learn Edition restriction"     | mdx (or 100号機) has Cubit Learn Edition; harmless ERROR line on export | Ignore — `export netgen` writes the .vol successfully despite the message; the smoke test tolerates this |
| Phase 9 cross-machine drift on `cubit_mesh_export.ccm` / `cubit_mesh_curver.pyd` | LAB has been doing local Cubit-plugin rebuilds; mdx / 100号機 have the published PyPI binary which differs slightly (PE timestamps, embedded paths) | Acceptable drift if no `src/cubit_plugin/` source change since cubit-mesh-export tag |
| Cubit Solve menu missing after first install | `.cubit` startup block, generated Radia startup shim, or Cubit.ini `plugin\\paths` was not registered | Re-run `cubit-plugin-install --all-users`; then require `cubit-plugin-install --verify-only --all-users` to pass |
| `cubit-plugin-install --verify-only` reports `[SIZE]`, `[HASH]`, `[STALE]`, or panel-registration failure | binary deploy partial / Cubit was running during last install / stale Qt5 `.ccl` / startup registration missing | Stop Cubit, re-run `cubit-plugin-install --all-users` |
| panel windows die in headless test with `QFontDatabase: Cannot find font directory`     | benign warning on `QT_QPA_PLATFORM=offscreen`                          | Ignore — not a failure; for production rendering Qt finds system fonts |

============================================================
## Tooling reference
============================================================

* `.claude/skills/deploy/SKILL.md` — deployment background for 100号機,
  mdx, and hibino.
* `tools/release_quad.py` — Phase 8a-8e deploy
  steps after a release.
* `tools/push_pyds_to_mdx.py` — LAB → mdx C++ artifact pusher (with
  lock-killer prelude).  RETIRED 2026-05-02 from standard release
  flow; preserved for branch-test scenarios only.
* `tools/release_quad.py done` — Phase 9 cross-machine drift
  probe (the Definition Of Done gate).
"""


_TOPICS = (
    "two_tier",
    "lab_editable",
    "hyaku_editable",
    "mdx_pypi",
    "hibino_pypi",
    "gui_extra",
    "editable_to_pypi_migration",
    "pypi_to_editable_migration",
    "metadata_sync",
    "pyd_dll_bootstrap",
    "cubit_plugin_layers",
    "first_install_panel_registration",
    "common_failure_modes",
)

# Backward-compat alias map (old topic name → new topic name).  Lets
# callers that learned the old API continue to function during the
# 2-tier transition.
_TOPIC_ALIASES = {
    "three_tier": "two_tier",
    "hyaku_pypi": "hyaku_editable",
    "mdx_editable": "mdx_pypi",
    "hibino": "hibino_pypi",
}


def get_install_deploy_documentation(topic: str = "") -> str:
    """Return the install/deploy knowledge.

    Args:
        topic: Empty for the full document, or one of the entries in
               ``_TOPICS`` above for a single section.  The historical
               topic names ``three_tier`` and ``mdx_editable`` map to
               the current ``two_tier`` and ``mdx_pypi`` for callers
               that learned the older API.
    """
    if not topic:
        return INSTALL_DEPLOY

    topic = _TOPIC_ALIASES.get(topic, topic)

    if topic not in _TOPICS:
        return (f"Unknown topic: {topic!r}. Available topics:\n"
                f"  {', '.join(_TOPICS)}\n\n"
                f"Pass empty string for the full document.")

    headers = []
    pos = 0
    while True:
        next_pos = INSTALL_DEPLOY.find("\n## ", pos)
        if next_pos < 0:
            break
        line_end = INSTALL_DEPLOY.find("\n", next_pos + 1)
        line = INSTALL_DEPLOY[next_pos + 1:line_end]
        for t in _TOPICS:
            if line.startswith(f"## {t} "):
                headers.append((t, next_pos + 1))
                break
        pos = next_pos + 1

    req_starts = [off for kw, off in headers if kw == topic]
    if not req_starts:
        return f"Topic {topic!r} declared but not found in document."
    section_start = INSTALL_DEPLOY.rfind("## ===", 0, req_starts[0])
    if section_start < 0:
        section_start = req_starts[0]

    section_end = len(INSTALL_DEPLOY)
    last_req = req_starts[-1]
    for kw, off in headers:
        if kw != topic and off > last_req:
            delim = INSTALL_DEPLOY.rfind("## ===", 0, off)
            section_end = delim if delim > 0 else off
            break

    return INSTALL_DEPLOY[section_start:section_end].rstrip() + "\n"
