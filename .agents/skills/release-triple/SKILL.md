---
name: release-triple
description: Simultaneously bump, tag, and PyPI-publish all three Radia packages (radia, cubit-mesh-export, radia-mcp) in one coordinated release. Use when the user asks for a release, tag push, PyPI upload, or wants to cut a new version. Handles version bumps across four files, one commit for all three, three separate tags, and push (which triggers auto-publish via GitHub Actions → PyPI OIDC).
---

# release-triple

## Canonical entry point: `tools/release_triple.py`

The full Phase 0 → 9 flow is automated. Use the orchestrator instead
of running ad-hoc commands:

```bash
python tools/release_triple.py preflight   # read-only state report
python tools/release_triple.py phase0      # mandatory clean rebuild
python tools/release_triple.py phase8      # deploy LAB + 100号機, verify, smoke
python tools/release_triple.py phase8e     # mdx PyPI upgrade (after PyPI live)
python tools/release_triple.py phase9      # cross-machine consistency probe
python tools/release_triple.py all         # phase8 -> phase8e -> phase9
python tools/release_triple.py done        # DEFINITION OF DONE: preflight + phase9
```

**Definition of Done**: a release is NOT done until
`python tools/release_triple.py done` exits 0. Do not say
"release complete" before this command prints
`DEFINITION OF DONE met`. The 2026-04-14 incidents all involved
a release being declared done while at least one of the gates was
quietly failing.

Each subcommand enforces its own preconditions:
- `phase8` refuses to run if `phase0` was not done since the last
  `src/cubit_plugin/` change (mtime check).
- `phase8e` refuses to run if PyPI hasn't propagated to the repo's
  current version (it polls `pip index versions`).
- `phase9` exits non-zero if any of the 11 fields drifts across LAB /
  100号機 / mdx.

The phases below (still narrated for human review) are now the
INTERNAL CONTRACT of the script. Read them when the script flags
something — they explain the WHY of each gate.

One atomic release covering all three packages in this monorepo:

| Package | Version file(s) | Tag prefix | PyPI name |
|---------|-----------------|------------|-----------|
| **radia** | `pyproject.toml` + `src/radia/__init__.py` | `v` | `radia` |
| **cubit-mesh-export** | `packages/cubit-mesh-export/pyproject.toml` | `cubit-mesh-export-v` | `cubit-mesh-export` |
| **radia-mcp** | `packages/radia-mcp/pyproject.toml` | `radia-mcp-v` | `radia-mcp` |

All three go out together so downstream users (100号機, mdx, external)
see a consistent set. radia-mcp depends on `radia`; cubit-mesh-export
ships binaries that radia's panels expect — mismatched versions are how
the 2026-04-14 sideset bug happened.

## Phase 0: mandatory clean rebuild of the Cubit plugin

**POLICY (2026-04-14)**: every triple release starts with a clean
rebuild of the Cubit C++ plugin. The plugin is no longer under heavy
development, so the rebuild cost is small compared to the risk of
shipping stale `.ccl/.ccm/.pyd` in the `cubit-mesh-export` wheel (which
is how the 2026-04-14 "source/sink not read" incident happened).

**Fast path (2–5 min): rebuild Cubit plugin targets only.** The full
`Build.ps1 -Rebuild` rebuilds Radia's _radia_pybind.pyd too (~10 min,
46 object files) and the Radia core rarely needs rebuilding for a
Cubit-plugin-only change. Skip ahead to the targeted build:

```bash
rm -rf S:/Radia/01_GitHub/src/cubit_plugin/build-pyd \
       S:/Radia/01_GitHub/src/cubit_plugin/build-ccm

cmd.exe /c 'call "C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\VC\Auxiliary\Build\vcvars64.bat" && \
  set "CUBIT_DIR=C:\Program Files\Coreform Cubit 2025.3\cmake" && \
  set "NETGEN_DIR=C:\Program Files\Python312\Lib\site-packages\netgen" && \
  set "CMAKE=C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe" && \
  set "SRC=S:\Radia\01_GitHub\src\cubit_plugin" && \
  mkdir "%SRC%\build-pyd" & cd /d "%SRC%\build-pyd" && \
  "%CMAKE%" -G Ninja -DCMAKE_BUILD_TYPE=Release -DCubit_DIR="%CUBIT_DIR%" -DNETGEN_DIR="%NETGEN_DIR%" "%SRC%" && \
  "%CMAKE%" --build . --config Release --target radia_cubit_mesh -j && \
  mkdir "%SRC%\build-ccm" & cd /d "%SRC%\build-ccm" && \
  "%CMAKE%" -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=cl -DCMAKE_CXX_COMPILER=cl -DCubit_DIR="%CUBIT_DIR%" -DNETGEN_DIR="%NETGEN_DIR%" "%SRC%" && \
  "%CMAKE%" --build . --config Release --target radia_cubit_ccm -j && \
  "%CMAKE%" --build . --config Release --target radia_cubit_ccl -j'

# Then manually propagate the fresh binaries (Build.ps1 does this
# automatically; the targeted build above skips it).
cp S:/Radia/01_GitHub/src/cubit_plugin/build-pyd/radia_cubit_mesh.cp312-win_amd64.pyd \
   S:/Radia/01_GitHub/src/radia/radia_cubit_mesh.pyd
cp S:/Radia/01_GitHub/src/cubit_plugin/build-pyd/radia_cubit_mesh.cp312-win_amd64.pyd \
   S:/Radia/01_GitHub/packages/cubit-mesh-export/src/cubit_mesh_export/radia_cubit_mesh.pyd
cp S:/Radia/01_GitHub/src/cubit_plugin/build-ccm/radia_cubit.ccm S:/Radia/01_GitHub/src/radia/
cp S:/Radia/01_GitHub/src/cubit_plugin/build-ccm/radia_cubit.ccm S:/Radia/01_GitHub/packages/cubit-mesh-export/src/cubit_mesh_export/
cp S:/Radia/01_GitHub/src/cubit_plugin/build-ccm/radia_cubit.ccl S:/Radia/01_GitHub/src/radia/
cp S:/Radia/01_GitHub/src/cubit_plugin/build-ccm/radia_cubit.ccl S:/Radia/01_GitHub/packages/cubit-mesh-export/src/cubit_mesh_export/
```

**Full rebuild fallback (10 min, rarely needed):**

```bash
powershell.exe -ExecutionPolicy Bypass -File "S:/Radia/01_GitHub/Build.ps1" -Rebuild
```

`Build.ps1 -Rebuild` deletes the CMake build dirs
(`src/cubit_plugin/build-pyd`, `build-ccm`, `build-msvc`) and rebuilds
from scratch, then propagates the artifacts to BOTH source-of-truth
directories per the Binary Propagation Policy:

- `src/radia/radia_cubit.ccl / .ccm / radia_cubit_mesh.pyd`
- `packages/cubit-mesh-export/src/cubit_mesh_export/...`

Verify before continuing:

```bash
# All three artifacts should have mtime within the last few minutes.
for f in \
  src/radia/radia_cubit.ccl \
  src/radia/radia_cubit.ccm \
  src/radia/radia_cubit_mesh.pyd \
  packages/cubit-mesh-export/src/cubit_mesh_export/radia_cubit.ccl \
  packages/cubit-mesh-export/src/cubit_mesh_export/radia_cubit.ccm \
  packages/cubit-mesh-export/src/cubit_mesh_export/radia_cubit_mesh.pyd; do
  stat -c "%y %n" "$f" 2>/dev/null || dir "$f"
done
```

If any artifact is older than the latest `src/cubit_plugin/*.cpp`
edit, rebuild did not propagate — investigate before proceeding.
Skip this phase only if no file under `src/cubit_plugin/` has changed
since the last release tag (verify with `git log --since=<last-tag>
-- src/cubit_plugin/`).

## Phase 1: decide version bumps

Check unreleased commits per package:

```bash
git log --oneline v<LAST>..HEAD --not -- packages/cubit-mesh-export/ packages/radia-mcp/  # radia
git log --oneline cubit-mesh-export-v<LAST>..HEAD -- packages/cubit-mesh-export/ src/cubit_plugin/ src/radia/radia_cubit.c*
git log --oneline radia-mcp-v<LAST>..HEAD -- packages/radia-mcp/ src/radia/mcp_server/
```

Bump rules (SemVer-ish, be honest about scope):

- **patch** (`X.Y.Z+1`): only bug fixes, no API change.
- **minor** (`X.Y+1.0`): new features, accumulated work, backward compatible.
- **major** (`X+1.0.0`): breaking API changes. Rare; confirm with user.

If the user says "release" without specifying, look at the log and
**propose** minor vs patch before bumping. Example past call (2026-04-14):
all three went minor because each had 5-20 commits since the last tag.

## Phase 2: bump the four version lines

```bash
# radia
# pyproject.toml:        version = "X.Y.Z"
# src/radia/__init__.py: __version__ = "X.Y.Z"

# cubit-mesh-export
# packages/cubit-mesh-export/pyproject.toml: version = "X.Y.Z"

# radia-mcp
# packages/radia-mcp/pyproject.toml: version = "X.Y.Z"
```

Keep radia's `pyproject.toml` and `__init__.py` in lockstep — they get
read by different tools and mismatch causes cryptic wheel-install bugs.
For radia-mcp the same lockstep applies:
`packages/radia-mcp/pyproject.toml` AND
`packages/radia-mcp/src/radia_mcp/__init__.py` BOTH have `version` /
`__version__` strings; the `--selftest` matrix and the
`test_radia_mcp_version_consistency` test compare them.

Use Edit tool on the specific `name = "<pkg>"` + `version = ...` pair
to avoid editing the wrong `version` line (there are several `version`
occurrences per `pyproject.toml`).

## Phase 2.5: pre-flight CI validation (MANDATORY 2026-05-03)

**POLICY (added 2026-05-03 after v4.27.0/v4.27.1 both round-tripped
through tag CI failures)**: run the same gates locally that CI runs,
BEFORE the first `git push`.  Each round-trip via tag CI burns 5-10
min and a version number; the local checks below take ~30 s and
catch every CI failure mode we hit on the v4.27.x retries.

### Gate 1: regenerate radia-mcp TOOLS.md and verify it matches

The CI test `test_radia_mcp_tools_doc_up_to_date` checks the shipped
`packages/radia-mcp/docs/TOOLS.md` against a re-generated version.
Whenever **any** subpackage gains/loses a tool, **changes name**,
or **is removed entirely** (e.g. the 0.38 -> 0.39 elf extraction),
TOOLS.md goes stale and the gate fails.

ALSO: `gen_tools_doc.py` itself has a hard-coded `SERVERS` list of
`(subpkg_name, console_script, blurb)` tuples.  If a subpackage was
removed without also editing this list, regeneration writes back an
"(import failed: ModuleNotFoundError)" stub — the test then fails
even after a regeneration.  This was the v4.27.1 -> v4.27.2 trip.

```bash
# 1. Regenerate
python packages/radia-mcp/scripts/gen_tools_doc.py

# 2. If TOOLS.md changed -- something was stale -- stage and verify:
git diff --stat packages/radia-mcp/docs/TOOLS.md
python -m pytest tests/mcp_server/test_tools_doc.py -xvs

# 3. If the regen produced "(import failed:" lines, edit
# packages/radia-mcp/scripts/gen_tools_doc.py SERVERS list to remove
# the dead subpackage entry, regen, and verify again.
```

### Gate 2: collect-only sweep with the CI exclusion set

CI imports every test module to collect them, even when --ignore'd
ones are skipped.  A stale import in a NON-ignored test module fails
collection and the whole CI run fails before a single test runs.
This was the v4.27.0 -> v4.27.1 trip (a leftover
`tests/mcp_server/test_build123d_gmsh_elf_responses.py` imported
`radia_mcp.elf` which had been deleted).

The CI ignore-set is in `.github/workflows/build-test.yml` "Run basic
tests" block.  Run the same sweep locally:

```bash
python -m pytest \
  --ignore=tests/cubit \
  --ignore=tests/panels \
  --ignore=tests/test_far_field_accuracy.py \
  --ignore=tests/test_rad_ngsolve_function.py \
  --ignore=tests/test_tetrahedral_solver.py \
  --ignore=tests/test_batch_evaluation.py \
  --ignore=tests/test_curlA_equals_B.py \
  --ignore=tests/test_mesh_import.py \
  --ignore=tests/test_scalar_bie_sibc.py \
  --collect-only -q
```

Expected: ``N tests collected in <1 minute>``, **NO** ``error`` lines,
**NO** ``Interrupted: N errors during collection``.  If any test
module fails to import: `grep -r "from <removed_module>" tests/` to
find the stale reference, fix or delete the test, and re-sweep.

### Gate 3: actually RUN the affected tests

```bash
# If radia-mcp / mcp_server tests changed:
python -m pytest tests/mcp_server/ -x

# If radia BEM / panel code changed:
python -m pytest tests/bem/ tests/panels/ -x  # tests/panels needs Cubit -- run on LAB
```

The test names that catch each historical regression class:

| Last seen | Test that catches it |
|-----------|----------------------|
| v4.27.0 (elf import) | `pytest tests/mcp_server/ --collect-only` |
| v4.27.1 (TOOLS.md stale) | `pytest tests/mcp_server/test_tools_doc.py` |
| v4.25.0 (panel CLI flag) | `panel-cli-diff` skill |
| v4.x.x (sideset regression) | `cubit-smoke-test` (in cubit-mesh-export) |
| Multi-machine drift | `tools/release_triple.py phase9` |

### Gate 4: workflow-runner cross-check the .yml --selftest matrix

If `--selftest` matrix in `.github/workflows/build-test.yml` was
edited (subpackage removed/added), grep that the listed servers all
import:

```bash
python -c "
import importlib
servers = ['radia_mcp.radia_ngsolve.server',
           'radia_mcp.cubit.server',
           'radia_mcp.build123d.server',
           'radia_mcp.gmsh.server',
           'radia_mcp.interop.server']
for s in servers:
    try: importlib.import_module(s); print(f'OK   {s}')
    except Exception as e: print(f'FAIL {s}: {e}')"
```

Every server must print `OK`.  Subpackages no longer in the
codebase must NOT be in the .yml.

**Only proceed to Phase 3 once all four gates above are green.**
The 2026-05-03 v4.27.x round-trip cost burned three patch numbers
and ~30 min of CI on each cycle; running the four gates locally
takes ~1 minute.

## Phase 3: stage exactly the release files

DO NOT use `git add -A` or `git add .`. The user typically has unstaged
WIP in other files, and untracked files should never enter a release
commit. Stage this precise set:

```bash
git add \
  pyproject.toml \
  src/radia/__init__.py \
  packages/cubit-mesh-export/pyproject.toml \
  packages/radia-mcp/pyproject.toml \
  <any fix files included in this release>
```

Then `git status --short` and eyeball the staged set. Unstaged WIP
(with `M ` vs ` M` distinction) should NOT be in the staged group.

Note: `.Codex/` is gitignored in this repo — skill changes are local
only and do not go into the release commit.

## Phase 4: commit with one composite message

Use HEREDOC to preserve formatting (AGENTS.md policy). Title names all
three tags; body explains the motivation plus per-package highlights.

```bash
git commit -m "$(cat <<'EOF'
Release v<X.Y.Z> / cubit-mesh-export-v<A.B.C> / radia-mcp-v<D.E.F>

<One-paragraph motivation: what drove this release, often a specific
incident or a feature crossing the "ship it" line. Lead with the WHY.>

- radia <X.Y.Z>: <highlights>
- cubit-mesh-export <A.B.C>: <highlights>
- radia-mcp <D.E.F>: <highlights>

Co-Authored-By: Codex Opus 4.6 (1M context) <noreply@anthropic.com>
EOF
)"
```

Verify the commit landed cleanly: `git log --oneline -2`.

## Phase 5: three annotated tags

```bash
git tag -a v<X.Y.Z>                  -m "radia v<X.Y.Z>: <one-line summary>"
git tag -a cubit-mesh-export-v<A.B.C> -m "cubit-mesh-export v<A.B.C>: <one-line summary>"
git tag -a radia-mcp-v<D.E.F>        -m "radia-mcp v<D.E.F>: <one-line summary>"

git tag --contains HEAD  # sanity: should list all three
```

## Phase 6: push main + all three tags

```bash
git push origin main
git push origin v<X.Y.Z> cubit-mesh-export-v<A.B.C> radia-mcp-v<D.E.F>
```

A git pre-push hook uploads .pyd/.ccl/.ccm binaries to a dedicated
`binaries` GitHub Release — this is expected output
("Uploading N binary files to GitHub Releases"), not an error.

## Phase 7: PyPI propagation is automatic (but verify)

Three GitHub Actions workflows trigger in parallel off the tag pushes:

| Tag pattern | Workflow | Target |
|-------------|----------|--------|
| `v*` | `release.yml` | PyPI `radia` |
| `cubit-mesh-export-v*` | `release-cubit-mesh-export.yml` | PyPI `cubit-mesh-export` |
| `radia-mcp-v*` | `release-radia-mcp.yml` | PyPI `radia-mcp` |

All three workflows wait for the `CI` workflow on the same ref to
succeed (`on: workflow_run`). Each uses `pypa/gh-action-pypi-publish`
with OIDC (`environment: pypi`) — no stored API tokens.

Monitor:

```bash
gh run list --limit 8
gh run watch <run-id>  # for specific run
gh release view v<X.Y.Z> --json url,assets 2>&1 || echo "not yet"

# PyPI propagation (usually <60 s after workflow success)
python -m pip index versions radia 2>&1 | head -3
python -m pip index versions cubit-mesh-export 2>&1 | head -3
python -m pip index versions radia-mcp 2>&1 | head -3
```

## Phase 8: downstream redeploy (MANDATORY on 100号機)

The user policy (feedback_deploy_responsibility.md) requires Codex to
deploy to 100号機 for **all users** after every release — do not leave
it to the lab users to `pip install` themselves.

### Phase 8a (REQUIRED FIRST): stop every Cubit process on the target

**POLICY (2026-04-14)**: Cubit holds `.ccm` / `.ccl` / `.pyd` files
open while running, and `shutil.copy` silently fails (OSError, swallowed
by `_clean_old_plugins`) leaving the user with a half-updated plugin.
The 2026-04-14 deploy hit this: the .ccm copy failed and the .ccl was
partially updated, forcing a retry after the user closed Cubit. Do NOT
run `cubit-plugin-install` while Cubit is up.

Check and stop Cubit on 100号機 (shared lab machine, any user may
have it open):

```bash
cat << 'PS' | ssh 192.168.11.100 'pwsh -ExecutionPolicy Bypass -Command -'
$procs = Get-Process -ErrorAction SilentlyContinue |
  Where-Object { $_.Name -like "*ubit*" -or $_.ProcessName -like "*oreform*" }
if ($procs) {
  $procs | Select-Object Id,ProcessName,StartTime,
    @{Name='User';Expression={(Get-CimInstance Win32_Process -Filter "ProcessId=$($_.Id)").GetOwner().User}} |
    Format-Table -AutoSize
  Write-Host "ACTION REQUIRED: Cubit is running — ask the owner(s) to save + close before proceeding."
} else {
  Write-Host "OK: no Cubit process running, safe to deploy."
}
# Also verify no file handles on the .ccm (belt-and-suspenders)
$ccm = "C:/Program Files/Coreform Cubit 2025.3/bin/plugins/radia_cubit.ccm"
if (Test-Path $ccm) {
  try { $fs = [System.IO.File]::Open($ccm,'Open','ReadWrite','None'); $fs.Close(); "  .ccm not locked" }
  catch { "  .ccm IS LOCKED: $($_.Exception.Message)" }
}
PS
```

If Cubit is running on **either LAB or 100号機**, force-terminate it.
The deploy-failure policy applies to every machine Codex deploys to
(LAB and 100号機 currently; mdx is verify-only). Waiting for any user
to close Cubit defeats the "Codex does deploy" policy and risks
leaving the same `.ccm` lock that caused the 2026-04-14 incident.
Cubit is re-launchable; users on these boxes save proactively.

LAB (local):
```bash
pwsh -Command "Get-Process -ErrorAction SilentlyContinue | Where-Object {
  \$_.ProcessName -eq 'coreform_cubit' -or \$_.ProcessName -eq 'cubit'
} | ForEach-Object { Write-Host \"Killing \$(\$_.ProcessName) pid=\$(\$_.Id)\"; Stop-Process -Id \$_.Id -Force }; Start-Sleep -Seconds 2"
```

100号機 (over SSH):
```bash
cat << 'PS' | ssh 192.168.11.100 'pwsh -ExecutionPolicy Bypass -Command -'
$procs = Get-Process -ErrorAction SilentlyContinue | Where-Object {
  $_.ProcessName -eq "coreform_cubit" -or $_.ProcessName -eq "cubit"
}
foreach ($p in $procs) {
  Write-Host "Force-terminating $($p.ProcessName) pid=$($p.Id) (started $($p.StartTime))"
  Stop-Process -Id $p.Id -Force
}
Start-Sleep -Seconds 2
PS
```

Then re-run the preflight (`cubit-plugin-install --check-only`) on the
target machine and proceed to 8b. mdx is exempt only because we do NOT
run cubit-plugin-install there (mdx is for PyPI-install verification,
not Cubit-bound deploys).

### Phase 8b (2026-05-02 reconfigured): deploy strategy after PyPI is live

**POLICY UPDATE (2026-05-02)**: install layout reconfigured from
"LAB = NAS editable, mdx = local clone editable, 100号機 = PyPI" to
**"LAB = NAS editable, 100号機 + mdx = PyPI"** (2-tier).  100号機 and
mdx now follow IDENTICAL deploy steps — both are PyPI-from-canonical
verification points with Cubit plugin deployed via `cubit-plugin-install
--all-users`.

**LAB** (editable, refresh metadata after release):
```bash
# Stop any mcp-server-*.exe to release Scripts/mcp-server-*.exe locks.
pwsh -Command "Get-Process -ErrorAction SilentlyContinue | Where-Object { \$_.Name -like 'mcp-server*' } | ForEach-Object { Stop-Process -Id \$_.Id -Force }; Start-Sleep -Seconds 2"

$repo = "S:/Radia/01_GitHub"   # NAS-mapped drive on LAB
# Re-editable to refresh metadata (radia-mcp __version__ etc.) after release.
pip install -e $repo --no-deps --no-cache-dir
pip install -e $repo/packages/cubit-mesh-export --no-deps --no-cache-dir
pip install -e $repo/packages/radia-mcp --no-deps --no-cache-dir
# No cubit-plugin-install on LAB normally; LAB Build.ps1 already populates
# src/radia/radia_cubit.{ccl,ccm} and packages/cubit-mesh-export/.../radia_cubit_mesh.pyd.
```

**100号機 / mdx** (PyPI install — IDENTICAL canonical deploy as of
2026-05-02).  These are the two end-to-end PyPI wheel verification
points, both Cubit-equipped, both running the standalone PySide6
panels.  Wait until `pip index versions radia / radia-mcp /
cubit-mesh-export` confirms the new versions are live on PyPI before
running these:

100号機 (lab production, 21+ users):
```bash
cat << 'PS' | ssh 192.168.11.100 'pwsh -ExecutionPolicy Bypass -Command -'
# Stop any mcp-server-*.exe and Cubit before pip install.
Get-Process -ErrorAction SilentlyContinue | Where-Object { $_.Name -like 'mcp-server*' -or $_.ProcessName -eq 'coreform_cubit' } |
  ForEach-Object { Stop-Process -Id $_.Id -Force }
Start-Sleep -Seconds 2

# Upgrade from PyPI (NOT NAS source).  Pin the version explicitly to
# avoid race with PyPI propagation.  [cubit,gui] extras pull
# cubit-mesh-export AND PySide6 (the standalone panel needs PySide6).
pip install --upgrade --no-cache-dir 'radia[cubit,gui]==<X.Y.Z>' 'radia-mcp==<X.Y.Z>' 'cubit-mesh-export==<X.Y.Z>'
# Deploy Cubit plugin from PyPI wheel as regular files for all 21 lab users.
cubit-plugin-install --all-users
PS
```

mdx (cross-machine consistency verification, no end-user load):
```bash
cat << 'PS' | ssh mdx 'pwsh -ExecutionPolicy Bypass -Command -'
# Stop any mcp-server-*.exe / Cubit / python.exe holding .pyd locks.
Get-Process -ErrorAction SilentlyContinue | Where-Object {
  $_.Name -like 'mcp-server*' -or $_.ProcessName -eq 'coreform_cubit' -or $_.ProcessName -eq 'cubit'
} | ForEach-Object { Stop-Process -Id $_.Id -Force }
Start-Sleep -Seconds 2

# Identical to 100号機: PyPI install + Cubit plugin deploy.
pip install --upgrade --no-cache-dir 'radia[cubit,gui]==<X.Y.Z>' 'radia-mcp==<X.Y.Z>' 'cubit-mesh-export==<X.Y.Z>'
cubit-plugin-install --all-users
PS
```

Note: `radia[cubit]` ALONE does NOT pull PySide6 — it only pulls
`cubit-mesh-export`.  Without `[gui]`, the standalone Radia-NGSolve
panel windows (radia_ih, radia_em, radia_pcb, radia_heat) cannot
launch.  Always use `radia[cubit,gui]` for production deploy on
Cubit-equipped machines.

Mandatory pip flags:

- `--no-cache-dir`: pip's wheel cache (`%LOCALAPPDATA%/pip/Cache`) will
  happily return a STALE wheel even when versions change. Always
  include this for any pip install in the release flow.

**Stale editable state migration**: if a machine still has the old
"editable from NAS" or "editable from local clone" install, the
`pip install --upgrade` from PyPI may NOT replace the editable. To
migrate: `pip uninstall -y radia radia-mcp cubit-mesh-export` first,
then run the PyPI install above. (The 2026-05-01 100号機 migration
and the 2026-05-02 mdx migration both confirmed this is the safe
sequence.)

### Phase 8c: verify deployment state matches package source

```bash
# LAB
cubit-plugin-install --verify-only

# 100号機
cat << 'PS' | ssh 192.168.11.100 'pwsh -ExecutionPolicy Bypass -Command -'
cubit-plugin-install --verify-only
PS

# mdx (added 2026-05-02 — mdx now also runs Cubit plugin via PyPI)
cat << 'PS' | ssh mdx 'pwsh -ExecutionPolicy Bypass -Command -'
cubit-plugin-install --verify-only
PS
```

`--verify-only` (cubit-mesh-export >= 0.5.4) re-checks AFTER the install
that every expected file is present in Cubit's plugin/bin dir, that
every deployed file's sha256 matches the package source, and that no
stale `radia_cubit.ccl` lingers in `bin/plugins/` (the historical
shadow location). Output looks like:

```
Verification:
  [OK] radia_cubit.ccm  (677888 bytes, sha256 match)
  [OK] radia_cubit.ccl  (280064 bytes, sha256 match)
  [OK] radia_cubit_mesh.cp312-win_amd64.pyd  (208384 bytes, sha256 match)
  [OK] nglib.dll        (...)
  [OK] ngcore.dll       (...)
  compat: radia X.Y.Z <-> cubit-mesh-export A.B.C compatible
  [OK] every expected binary present and matches package source.
```

Any `[MISS] / [SIZE] / [HASH] / [STALE]` line is a deploy failure —
re-run cubit-plugin-install (after closing Cubit) and re-verify.

### Phase 8d: verify with cubit-smoke-test — INCLUDES launcher GUI test

cubit-smoke-test (>= 0.5.8) now drives `radia_export verify_launcher`
inside the Cubit -batch session, AND parses `VERIFY_LAUNCHER:` lines
from stdout. The headless launcher probe asserts exactly the same
state the Qt dialog computes:

- `get_model_labels()` returns every block + sideset name. Detects
  the 2026-04-14 `get_entity_name("sideset",...)` regression in
  seconds.
- `default_vol_path(jouPath)` is sane (catches "/" → "/.vol" bug).
- Every label in `--expect` and `--expect-materials` is reachable
  through the same code path the dialog runs (shared
  `LauncherLogic.hpp`).

Mandatory: a deploy is NOT done until cubit-smoke-test prints
`[OK] round-trip healthy` AND `verify_launcher: ... status=OK`. The
2026-04-14 incidents (sideset MISSING + .vol path "/.vol") would
each have been caught at this gate.

```bash
cat << 'PS' | ssh 192.168.11.100 'pwsh -ExecutionPolicy Bypass -Command -'
cubit-smoke-test
PS
cubit-smoke-test   # also on LAB
```

`cubit-smoke-test` (shipped with cubit-mesh-export >= 0.5.3) runs Cubit
in `-batch -nographics` mode on the canonical `ih_bem_sample.jou`,
exports a .vol via `radia_export netgen`, and independently asserts:

- **bcnames** (sideset labels) contains `source`, `sink`, `sibc`
- **materials** (block names) contains `coil`, `workpiece`, `air`

These are two separate sections of the Netgen .vol text format — do not
try to find a block name in bcnames or vice versa (this was the 0.5.2 →
0.5.3 false-FAIL).

Exit 0 = the deploy is genuinely healthy (the 2026-04-14 incident would
have been caught in seconds by this). Non-zero = STOP, investigate, do
NOT call the release done.

Known quirk: Cubit's `-batch -nographics` mode sometimes segfaults
during mesh teardown AFTER `radia_export netgen` has already written
the .vol. `cubit-smoke-test` tolerates a non-zero Cubit exit code as
long as the .vol is present and valid — it prints `[WARN] Cubit exited
0xC0000005 ...` and continues.

Other sanity checks that `cubit-plugin-install` now does automatically
(v0.5.2 hardening) and that you should still visually confirm:

- Every expected file printed as `[OK]` (no `[FAIL]`, no silent skips).
- Preflight block reports `[OK] no Cubit processes, no locked plugin files`.
- Compat line shows `radia X.Y.Z <-> cubit-mesh-export A.B.C compatible`.
- `C:/radia_panel_log.txt` at next Cubit launch has no `PLUGIN_STALE:`
  line.

See also `deploy` skill for the full 2-stage story (LAB editable +
100号機/mdx PyPI).

## Phase 8e (2026-05-02 reconfigured): mdx PyPI verification

**OLD (2026-05-01)**: mdx ran editable from a local clone at
`C:\Radia\01_GitHub`.
**NEW (2026-05-02)**: mdx is a **PyPI install** target identical to
100号機 — same `pip install --upgrade --no-cache-dir 'radia[cubit,gui]
==<X.Y.Z>' ...` + `cubit-plugin-install --all-users` recipe (already
covered in Phase 8b).

mdx now also runs `cubit-plugin-install --verify-only` (Phase 8c) and
the standalone PySide6 panel smoke (Phase 8d-equivalent).  This phase
collapsed into Phase 8b/8c/8d when the 2-tier policy was adopted; it
exists here for backward reference only.

If mdx still has the previous editable install (local clone at
`C:\Radia\01_GitHub`), the migration recipe is:

```bash
cat << 'PS' | ssh mdx 'pwsh -ExecutionPolicy Bypass -Command -'
# Stop locks first.
Get-Process -ErrorAction SilentlyContinue | Where-Object {
  $_.Name -like 'mcp-server*' -or $_.ProcessName -eq 'coreform_cubit' -or $_.ProcessName -eq 'cubit'
} | ForEach-Object { Stop-Process -Id $_.Id -Force }
Start-Sleep -Seconds 2

# Tear down editable install before PyPI takeover.
pip uninstall -y radia radia-mcp cubit-mesh-export

# Identical to 100号機.
pip install --no-cache-dir 'radia[cubit,gui]==<X.Y.Z>' 'radia-mcp==<X.Y.Z>' 'cubit-mesh-export==<X.Y.Z>'
cubit-plugin-install --all-users
cubit-plugin-install --verify-only
PS
```

The local clone at `C:\Radia\01_GitHub` may stay (useful for fixture
reads, e.g. `tests/panels/golden/*.step`) but Python no longer reads
from it.

## Phase 9: cross-machine consistency verification (MANDATORY)

The three machines now use only TWO install paths (2026-05-02 reconfigured):

| Machine | Install source | Why |
|---------|----------------|-----|
| LAB | NAS repo (`pip install -e S:/Radia/01_GitHub`) | dev iteration, gets Phase 0 build artifacts immediately |
| 100号機 | PyPI (`pip install --upgrade 'radia[cubit,gui]==X.Y.Z' ...`) | shared lab box, end-to-end PyPI wheel + Cubit plugin verification |
| mdx | PyPI (`pip install --upgrade 'radia[cubit,gui]==X.Y.Z' ...`) | second PyPI verification point, isolated from NAS (cross-machine consistency probe) |

The point of Phase 9 is to PROVE the three install paths converge to
the same runtime code. Drift here means the NAS build and the PyPI
wheel diverged — a Phase 2 / Phase 3 / CI bug that must be fixed
before declaring the release done.

Compared per file:

- same `radia.__version__` / `cubit_mesh_export.__version__`
- same compat constants
- same SHA256 on the canonical Python files (`register_toolbar.py`,
  `calc_inductance.py`, `calc_heating_bem.py`, `bem_inductance.py`)
- same SHA256 on the bundled Cubit plugin (`radia_cubit.ccm`,
  `radia_cubit.ccl`) — proves the NAS-built binary and the PyPI-shipped
  binary are bit-identical

Why each file matters:

| File | Why drift would matter |
|------|------------------------|
| `register_toolbar.py` | startup hook + freshness/compat checks |
| `calc_inductance.py`  | IH BEM solver entry point |
| `calc_heating_bem.py` | BEM-SIBC solver entry point |
| `bem_inductance.py`   | core EFIE saddle-point routine |
| `radia_cubit.ccm`     | Cubit menu / APREPRO commands |
| `radia_cubit.ccl`     | Cubit Qt component (launch_radia_ngsolve) |

The script below queries all three machines in one pass and prints a
matrix; any **DRIFT** row means the environments are not equivalent.

```bash
cat << 'PYALL' > /tmp/cross_machine_verify.py
"""Print radia.__version__ + content-hash of canonical files.

Files are hashed AFTER newline normalization (CRLF -> LF) so the LAB
NAS install (Git's LF) and the mdx PyPI install (Windows CI runner's
git autocrlf -> CRLF in wheel) report identical hashes. Binary
artifacts (.ccm/.ccl) hash byte-for-byte.
"""
import hashlib, os
import importlib.metadata as md


def sha256_text_normalized(path):
    """Hash the file with all CRLF/CR collapsed to LF — equivalent
    text content reports the same hash regardless of line endings."""
    h = hashlib.sha256()
    try:
        with open(path, "rb") as f:
            data = f.read()
        data = data.replace(b"\\r\\n", b"\\n").replace(b"\\r", b"\\n")
        h.update(data)
        return h.hexdigest()[:12]
    except Exception:
        return "MISSING"


def sha256_bytes(path):
    h = hashlib.sha256()
    try:
        with open(path, "rb") as f:
            for c in iter(lambda: f.read(65536), b""):
                h.update(c)
        return h.hexdigest()[:12]
    except Exception:
        return "MISSING"


def ver(name):
    try:
        return md.version(name)
    except Exception:
        return "MISSING"


import radia, cubit_mesh_export
rad_dir = os.path.dirname(radia.__file__)
cme_dir = os.path.dirname(cubit_mesh_export.__file__)
print(f"VER radia              = {radia.__version__}")
print(f"VER cubit-mesh-export  = {cubit_mesh_export.__version__}")
print(f"VER radia-mcp          = {ver('radia-mcp')}")
print(f"COMPAT cme  -> radia   = [{cubit_mesh_export.COMPAT_RADIA_MIN}, {cubit_mesh_export.COMPAT_RADIA_MAX}]")
print(f"COMPAT rad  -> cme     = [{radia.COMPAT_CUBIT_MESH_EXPORT_MIN}, {radia.COMPAT_CUBIT_MESH_EXPORT_MAX}]")
for rel in ["panels/register_toolbar.py", "panels/calc_inductance.py",
            "panels/calc_heating_bem.py", "bem_inductance.py"]:
    print(f"SHA radia/{rel:35s} = {sha256_text_normalized(os.path.join(rad_dir, rel))}")
for rel in ["radia_cubit.ccm", "radia_cubit.ccl"]:
    print(f"SHA cme/{rel:35s} = {sha256_bytes(os.path.join(cme_dir, rel))}")
PYALL

# Run the same probe on each machine, capture into separate files.
python /tmp/cross_machine_verify.py > /tmp/probe_lab.txt
cat /tmp/cross_machine_verify.py | ssh 192.168.11.100 'python -' > /tmp/probe_100.txt
cat /tmp/cross_machine_verify.py | ssh mdx              'python -' > /tmp/probe_mdx.txt

# Side-by-side comparison. Trim whitespace before equality test —
# pasted columns can have trailing spaces on Windows.
paste -d'|' /tmp/probe_lab.txt /tmp/probe_100.txt /tmp/probe_mdx.txt | awk -F'|' '
BEGIN { printf "%-44s | %-22s | %-22s | %-22s\n", "field", "LAB", "100号機", "mdx"
        print "----------------------------------------------------------------------------------------------------------------" }
{ a=$1; b=$2; c=$3
  key=a; sub(/ *=.*/,"",key); sub(/^[A-Z]+ /,"",key)
  va=a; sub(/.*= */,"",va); gsub(/[ \\r\\n]+$/,"",va)
  vb=b; sub(/.*= */,"",vb); gsub(/[ \\r\\n]+$/,"",vb)
  vc=c; sub(/.*= */,"",vc); gsub(/[ \\r\\n]+$/,"",vc)
  mark = (va==vb && vb==vc) ? "OK" : "DRIFT"
  printf "%-44s | %-22s | %-22s | %-22s  [%s]\n", key, va, vb, vc, mark
}'
```

Expected: every row prints `[OK]`. Any `[DRIFT]` row indicates one
machine is out of sync — investigate and re-deploy that machine.

Acceptable drifts (document if seen):

- LAB / 100号機 may have a NEWER version than mdx if a release was
  cut but PyPI propagation is in flight. Re-run after PyPI updates.
- mdx will lack PySide6-only modules (no GUI install path); attempting
  to import `radia.radia_ih` raises ModuleNotFoundError. The probe
  script intentionally avoids GUI-only files.

## Pitfalls to avoid

- **No `git add -A`.** User has WIP; release commits must be surgical.
- **Don't forget `__init__.py`.** `pyproject.toml` alone is not enough
  for radia — the `__version__` string is imported at runtime and
  anyone checking it will see the old version.
- **Don't amend past the commit.** If a hook or CI fails, fix forward
  with a new commit, new tag (e.g. `v4.5.1`). Rewriting a pushed tag
  is forbidden unless the user explicitly authorizes it.
- **Tag order matters.** Push main BEFORE the tags. Pushing tags to a
  commit that isn't on the remote yet breaks CI's `workflow_run`
  dispatch on head_branch detection.
- **PyPI is immutable.** Once a version is published you cannot reuse
  the number. If CI fails after push, bump the patch and try again.
  Never reuse a version.
- **Don't force-push main.** Even if the release commit has a typo.
- **`pip` cache returns stale wheels.** When pip-installing from a
  source dir that was previously installed at a lower version, pip's
  wheel cache happily returns the OLD wheel. The version number in
  pyproject.toml gets ignored. Always pass `--no-cache-dir` for
  NAS-source `pip install` (Phase 8b commands already include it).
  2026-04-14: hit this with radia-mcp 0.2.2 sticking around after
  bumping to 0.3.0.
- **`mcp-server-*.exe` locks block radia-mcp install.** Each MCP server
  console script under `Python312/Scripts/` is held open while the
  server runs (Codex itself uses these). `pip install --force-reinstall
  radia-mcp` will fail with `WinError 32: process cannot access the
  file`. Stop them first:
  ```bash
  pwsh -Command "Get-Process | Where-Object { \$_.Name -like 'mcp-server*' } | Stop-Process -Force"
  ```
  Phase 8b on LAB and 100号機 both include this kill step. Do it
  before any radia-mcp install.
