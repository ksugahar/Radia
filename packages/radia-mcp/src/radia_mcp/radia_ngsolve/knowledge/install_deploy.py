"""
Radia install / deploy policy and recipes — 3-tier configuration
(LAB / mdx / 100号機), reversible migration steps, and the non-obvious
gotchas that cause silent breakage.

Read this when:
* Setting up a new lab machine.
* Migrating an existing machine between editable / PyPI install.
* Diagnosing "import works but pip says wrong version" or
  "DLL load failed" on a freshly-deployed machine.

The MCP server exposes this via install_deploy(topic=...). Topics:
three_tier, lab_editable, mdx_editable, hyaku_pypi,
editable_to_pypi_migration, pypi_to_editable_migration,
metadata_sync, pyd_dll_bootstrap, cubit_plugin_layers,
common_failure_modes
"""

INSTALL_DEPLOY = """\
# Radia install / deploy policy (2026-05-01 reconfigured, 2026-05-02
# completed)

This document is the canonical record of how the three Radia
machines (LAB, mdx, 100号機) are installed and kept in sync. It
supersedes any earlier description that says "100号機 = editable,
mdx = PyPI" — that was the old policy.

Topics: three_tier, lab_editable, mdx_editable, hyaku_pypi,
        editable_to_pypi_migration, pypi_to_editable_migration,
        metadata_sync, pyd_dll_bootstrap, cubit_plugin_layers,
        common_failure_modes

============================================================
## three_tier — current 3-tier configuration (2026-05-01)
============================================================

| Machine  | Install                                   | Source              |
|----------|-------------------------------------------|---------------------|
| LAB      | editable (`pip install -e`)               | NAS `S:\\Radia\\01_GitHub` |
| mdx      | editable (`pip install -e`)               | local clone `C:\\Radia\\01_GitHub` |
| 100号機  | PyPI (`pip install radia==X.Y.Z`)         | `C:\\Program Files\\Python312\\Lib\\site-packages` (machine-wide, 21 users) |

Roles:

* **LAB** — dev iteration; Build.ps1 outputs are immediately live.
* **mdx** — LAB-independent dev iteration + headless verification +
  Cubit headless smoke test (cubit-smoke-test). mdx HAS Cubit
  installed at `C:\\Program Files\\Coreform Cubit 2025.3` and the
  Radia plugin (`cubit-plugin-install` deploys to the bin\\plugins
  directory).
* **100号機** — end-to-end PyPI wheel + Cubit plugin verification on
  21-user shared lab box. Released as PyPI by `cubit-plugin-install
  --all-users`.

The OLD configuration (LAB+100=editable from NAS, mdx=PyPI) was
swapped on 2026-05-01. The reason: 100号機 is the production
deployment target where 21 lab users actually run Radia, so it
should test what they see (PyPI). mdx is for headless dev iteration
where editable is more convenient.

============================================================
## lab_editable — LAB editable install
============================================================

LAB has the canonical NAS-mapped working tree at `S:\\Radia\\01_GitHub`.
All editing happens here. Build.ps1 compiles C++ artifacts directly
into the local clone (`src/radia/_radia_pybind.pyd`,
`packages/cubit-mesh-export/.../radia_cubit_mesh.pyd`, etc.) and the
editable install picks them up immediately.

After every release, refresh the LAB metadata:

```powershell
pip install -e . --no-deps --no-cache-dir
pip install -e packages/cubit-mesh-export --no-deps --no-cache-dir
pip install -e packages/radia-mcp --no-deps --no-cache-dir
```

This is the **metadata_sync** step (see topic). Don't skip it after a
release — `importlib.metadata.version("radia")` will lag.

LAB ↔ 100号機 NAS-share is preserved (filesystem same), but 100号機's
Python install path is now independent (PyPI), so LAB edits do NOT
auto-propagate to 100号機 anymore.

============================================================
## mdx_editable — mdx editable install (full recipe)
============================================================

mdx has its own git clone at `C:\\Radia\\01_GitHub`. Editable install
for all 3 packages. C++ binaries come from LAB via base64-over-ssh
push (mdx has no MSVC + Intel MKL + NGSolve build environment).

**Critical gotchas** (each was a 30-min debug):

1. `cubit_mesh_export.py` legacy single-file shadow at
   `C:\\Program Files\\Python312\\Lib\\site-packages\\cubit_mesh_export.py`
   silently overrides the editable package. `pip uninstall` does NOT
   touch it. MUST be deleted manually before `pip install -e`.

2. mdx has **no `gh` CLI**. `download_binaries.sh` reports success
   but writes nothing. Use `tools/push_pyds_to_mdx.py` from LAB
   instead.

3. `radia_cubit_mesh.pyd` lives in TWO paths on the local clone:
   `src/radia/` AND
   `packages/cubit-mesh-export/src/cubit_mesh_export/`. Both must be
   populated; the `tools/push_pyds_to_mdx.py` mirrors automatically.

4. Three packages, not two: `radia`, `cubit-mesh-export`,
   `radia-mcp`. Don't forget radia-mcp.

5. `python.exe` / `pythonw.exe` / `jupyter*` / `mcp-server-*` lock
   `.pyd` files when they have imported radia. The push tool's
   `stop_locks_on_mdx()` kills them all before pushing.

Full reference: `.claude/skills/deploy/SKILL.md` Stage 2 mdx
(Steps A/B/C/D + 7-row failure-mode table).

============================================================
## hyaku_pypi — 100号機 PyPI install
============================================================

100号機 (the lab's shared 21-user box) receives Radia exclusively via
PyPI install. All three packages go to machine-wide site-packages.

After every release, the admin runs:

```powershell
# Stop locks
Get-Process -ErrorAction SilentlyContinue | Where-Object {
    $_.Name -like 'mcp-server*' -or
    $_.ProcessName -in 'coreform_cubit','cubit'
} | ForEach-Object { Stop-Process -Id $_.Id -Force }
Start-Sleep -Seconds 2

# Upgrade from PyPI (NOT NAS source)
pip install --upgrade --no-deps --no-cache-dir \\
    radia==<X.Y.Z> radia-mcp==<X.Y.Z> cubit-mesh-export==<X.Y.Z>

# Deploy Cubit plugin from PyPI wheel for all 21 users
cubit-plugin-install --all-users
```

The `--all-users` flag writes the Cubit plugin to
`C:\\Program Files\\Coreform Cubit 2025.3\\bin\\plugins\\` (admin
privilege required) so all users see the same version.

============================================================
## editable_to_pypi_migration — switch a machine from editable to PyPI
============================================================

Use this when, e.g., 100号機 was on NAS-editable and you're moving it
to PyPI (the 2026-05-01 migration). Procedure:

```powershell
# Stop locks (mcp-server, Cubit, jupyter)
Get-Process -ErrorAction SilentlyContinue | Where-Object {
    $_.Name -like 'mcp-server*' -or
    $_.ProcessName -in 'coreform_cubit','cubit','jupyter','jupyter-lab'
} | ForEach-Object { Stop-Process -Id $_.Id -Force }
Start-Sleep -Seconds 2

# Uninstall editable (frees the .pth file, keeps clone intact)
pip uninstall -y radia radia-mcp cubit-mesh-export

# Install from PyPI
pip install --no-cache-dir --no-deps \\
    radia==<X.Y.Z> radia-mcp==<X.Y.Z> cubit-mesh-export==<X.Y.Z>

# Cubit plugin: regular-file deploy from PyPI wheel
cubit-plugin-install --all-users
cubit-plugin-install --verify-only   # all sha256 must match
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

The 2026-05-01 100号機 migration completed cleanly with this
recipe (5 sha256 [OK]'s on cubit-plugin-install --verify-only).

============================================================
## pypi_to_editable_migration — switch a machine from PyPI to editable
============================================================

Use this when, e.g., mdx had a PyPI install and you want it editable
(the 2026-05-01 mdx migration was easier — mdx already had editable).

```powershell
# Stop locks
Get-Process -ErrorAction SilentlyContinue | Where-Object {
    $_.Name -in 'python','pythonw','jupyter','jupyter-lab','ipykernel' -or
    $_.Name -like 'mcp-server*'
} | ForEach-Object { Stop-Process -Id $_.Id -Force -ErrorAction SilentlyContinue }
Start-Sleep -Seconds 3

# Delete legacy single-file shadow if present (see mdx_editable topic)
$shadow = "C:\\Program Files\\Python312\\Lib\\site-packages\\cubit_mesh_export.py"
$shadowPyc = "C:\\Program Files\\Python312\\Lib\\site-packages\\__pycache__\\cubit_mesh_export.cpython-312.pyc"
if (Test-Path $shadow) { Remove-Item $shadow -Force }
if (Test-Path $shadowPyc) { Remove-Item $shadowPyc -Force }

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

Then run `tools/push_pyds_to_mdx.py` from LAB to drop the C++ binaries.

============================================================
## metadata_sync — pip metadata vs `radia.__version__`
============================================================

Editable installs have a separate metadata record (the
`<package>-<version>.dist-info` directory). Without explicit refresh,
`importlib.metadata.version("radia")` and `pip list` can lag behind
`radia.__version__` (which is read from the local clone's
`src/radia/__init__.py`). After every release on an editable
machine:

```powershell
# Re-running pip install -e regenerates the .dist-info / .pth file
pip install -e . --no-deps --no-cache-dir
pip install -e packages/cubit-mesh-export --no-deps --no-cache-dir
pip install -e packages/radia-mcp --no-deps --no-cache-dir
```

`tools/release_triple.py done` runs a Phase 9 cross-machine probe
that compares `__version__` and `pip list` outputs. A DRIFT row in
that report usually means metadata_sync was skipped on one machine.

For 100号機 (PyPI install), the metadata is synced automatically by
`pip install --upgrade radia==<X.Y.Z>` because PyPI installs are
not editable.

============================================================
## pyd_dll_bootstrap — `radia_cubit_mesh.pyd` requires `import radia` first
============================================================

`radia_cubit_mesh.pyd` (the cubit-mesh-export C++ extension) depends
on Netgen DLLs (`nglib.dll`, `ngcore.dll`). On import these DLLs must
be findable via the Windows DLL search path.

`radia/__init__.py` calls `os.add_dll_directory(...)` for the
NGSolve / Netgen / MKL DLL directories. So:

```python
import radia                    # First — sets DLL paths
from cubit_mesh_export import radia_cubit_mesh   # Now this loads
```

Without the `import radia` first, you get
`ImportError: DLL load failed while importing radia_cubit_mesh: ...`.

Workarounds (none of them are encouraged; just `import radia` first):

* Manually call `os.add_dll_directory(r"C:\\Program Files\\Python312\\Lib\\site-packages\\netgen\\lib")` before the import.
* Add the netgen lib dir to PATH at the OS level.

This issue applies to mdx editable specifically; on a PyPI install
(100号機) the wheel installer handles DLL search paths via its own
hooks.

============================================================
## cubit_plugin_layers — Cubit plugin lives in TWO independent places
============================================================

The Cubit plugin (`radia_cubit.ccm`, `radia_cubit.ccl`,
`radia_cubit_mesh.cp312-win_amd64.pyd`, `nglib.dll`, `ngcore.dll`) is
deployed to `C:\\Program Files\\Coreform Cubit 2025.3\\bin\\` by
`cubit-plugin-install`. This is INDEPENDENT of the Python editable
install: Cubit reads its plugin directory directly, not via Python's
import system.

Two layers:

| Layer | Location | Used by |
|-------|----------|---------|
| Python import | `<package>/cubit_mesh_export/` | `from cubit_mesh_export import ...` |
| Cubit plugin  | `Coreform Cubit 2025.3\\bin\\plugins\\` (+ ccl in `bin\\`) | Cubit GUI / `coreform_cubit -batch` |

Implications:

* On 100号機 (PyPI install), `pip install --upgrade
  cubit-mesh-export==X.Y.Z` updates the Python side. You MUST then
  run `cubit-plugin-install --all-users` to update the Cubit side.
* On mdx (editable install), `pip install -e ...` updates Python.
  `tools/push_pyds_to_mdx.py` updates the local clone's `.pyd`.
  `cubit-plugin-install` updates the Cubit side. All three are
  required after a release that bumps the C++ plugin.
* On LAB (editable install), Build.ps1 writes the binaries to BOTH
  the local clone AND the Cubit\\bin\\plugins location (via the
  Build.ps1 post-build step). No separate cubit-plugin-install
  needed in routine dev.

============================================================
## common_failure_modes — symptoms and fixes
============================================================

| Symptom                                    | Root cause                                                                | Fix |
|-------------------------------------------|---------------------------------------------------------------------------|-----|
| `radia.__version__` says X.Y.Z but `pip list` says A.B.C | metadata_sync skipped on editable machine                              | `pip install -e <path> --no-deps --no-cache-dir` |
| `ImportError: cannot import name 'check' from 'cubit_mesh_export'` | legacy `cubit_mesh_export.py` shadow at site-packages                | Delete `Lib/site-packages/cubit_mesh_export.py` and `__pycache__/cubit_mesh_export.cpython-312.pyc`, then `pip install -e packages/cubit-mesh-export --no-deps --no-cache-dir` |
| `AttributeError: module 'cubit_mesh_export' has no attribute '__version__'` | same legacy shadow                                                     | same |
| `ImportError: DLL load failed while importing radia_cubit_mesh` | `import radia` not done first (DLL bootstrap)                          | `import radia` before `from cubit_mesh_export import radia_cubit_mesh` |
| `WriteAllBytes: WinError 32 process cannot access the file` (during `tools/push_pyds_to_mdx.py`) | `python.exe` / `mcp-server-*` holds the .pyd open                    | Tool's `stop_locks_on_mdx()` kills them; if it still fails, kill manually + retry |
| `cubit-smoke-test` fails on 100号機 with "Learn Edition restriction"     | 100号機 has Cubit Learn Edition; harmless ERROR line on radia_export   | Ignore — `radia_export netgen` writes the .vol successfully despite the message; the smoke test tolerates this |
| `Phase 9 cross-machine drift on cme/radia_cubit.ccl` | LAB has been doing local Cubit-plugin rebuilds; mdx / 100号機 have the published PyPI binary which differs slightly (PE timestamps, embedded paths) | Acceptable drift if no `src/cubit_plugin/` source change since cubit-mesh-export tag |
| `cubit-plugin-install --verify-only` reports `[SIZE]` or `[HASH]` mismatch | binary deploy partial / Cubit was running during last install         | Stop Cubit, re-run `cubit-plugin-install --all-users` |

============================================================
## Tooling reference
============================================================

* `.claude/skills/deploy/SKILL.md` — Stage 2 mdx editable
  (Steps A/B/C/D), Stage 3 100号機 PyPI.
* `.claude/skills/release-triple/SKILL.md` — Phase 8a-8e deploy
  steps after a release.
* `tools/push_pyds_to_mdx.py` — LAB → mdx C++ artifact pusher (with
  lock-killer prelude).
* `tools/release_triple.py done` — Phase 9 cross-machine drift
  probe (the Definition Of Done gate).
"""


_TOPICS = (
    "three_tier",
    "lab_editable",
    "mdx_editable",
    "hyaku_pypi",
    "editable_to_pypi_migration",
    "pypi_to_editable_migration",
    "metadata_sync",
    "pyd_dll_bootstrap",
    "cubit_plugin_layers",
    "common_failure_modes",
)


def get_install_deploy_documentation(topic: str = "") -> str:
    """Return the install/deploy knowledge.

    Args:
        topic: Empty for the full document, or one of the entries in
               ``_TOPICS`` above for a single section.
    """
    if not topic:
        return INSTALL_DEPLOY

    if topic not in _TOPICS:
        return (f"Unknown topic: {topic!r}. Available topics:\n"
                f"  {', '.join(_TOPICS)}\n\n"
                f"Pass empty string for the full document.")

    # Find ``## <topic> ``-prefixed headers and slice between the
    # previous ``## ===`` delimiter and the next distinct-topic header.
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
