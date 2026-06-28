---
name: deploy
description: Deploy Radia + Cubit plugin to remote machines (100号機, mdx) via SSH
---

# Deploy (2-Stage: LAB editable, 100号機+mdx PyPI)

**Philosophy**: only the developer (LAB) gets the editable convenience.
Every other machine receives the **published PyPI wheel** so RECORD /
entry-points / package-data / cubit-plugin-install regular-file deploy
are exercised the same way external users see them.  SMB `cp` shortcuts
remain **forbidden** -- they hide bugs that only surface under a real
`pip install`.

```
Stage 1  LAB                : pip install -e .                  editable, fastest feedback
                              ↓ the DEVELOPER workflow
Stage 2  100号機 / mdx      : pip install 'radia[cubit]'        PyPI, multi-user distribution
                              ↓ the PUBLIC PyPI test (real-world install)
```

**POLICY (2026-05-02, user-set, REPLACES 2026-05-01 3-tier)**: 100号機
と mdx は **両方 PyPI install** に統一。3-tier (LAB editable / mdx
editable / 100号機 PyPI) を 2-tier (LAB editable / 100号機 + mdx PyPI)
へ簡素化。
- **LAB**: `pip install -e .` (NAS source).  developer iteration.
- **100号機**: `pip install 'radia[cubit]==X.Y.Z'` + `cubit-plugin-install
  --all-users`.  21 ラボメンバの本番。
- **mdx**: 100号機 と完全に同じ recipe.  PyPI cross-machine consistency
  probe (release-qud Phase 9).  C++ 変更を mdx に届けるには `v*
  tag → CD/PyPI` リリースが必要 (LAB の Build.ps1 出力は届かない).

**Cubit plugin (100号機 + mdx)**: 両方とも PyPI wheel から regular-file
deploy (`cubit-plugin-install --all-users`).  symlink 体制は使わない.

**Notebook policy (2026-06-28)**: production deploy uses `radia[cubit]` and
does not add the old GUI extra.  The canonical panel surface is the Jupyter
notebook workbench (`ipynb-gui-health`), and normal Radia Python should not
acquire PySide6.  Coreform Cubit's bundled PySide6 remains protected because
Cubit owns that embedded runtime.

**Retired** (2026-05-02): mdx editable + `tools/push_pyds_to_mdx.py`
base64-over-ssh push は通常運用から除外.  `push_pyds_to_mdx.py` は
branch-test 用に残置 (PyPI tag を切らずに pre-release `.pyd` を mdx
で試す場合のみ使用).

Each stage uses pip install (editable or from-PyPI).  **Never** copy
files to site-packages by hand or via SMB -- if pip install doesn't see
it, the real distribution channel won't either.

3 packages in monorepo (ksugahar/Radia), independent PyPI release:

| Package | Install (Stage 3) | Purpose |
|---------|---------|---------|
| **radia** | `pip install radia` | C++ core + Python (FEM/BEM/PEEC, panels) |
| **cubit-mesh-export** | `pip install cubit-mesh-export` | High-order curved mesh export from Cubit |
| **radia-mcp** | `pip install radia-mcp` | MCP servers (radia, ngsolve, cubit, gmsh) |

## Check Matrix

Each check layer catches a specific class of failure.  Upper layers
depend on lower layers -- skipping L0/L1 means higher-layer failures
become mysterious.  All three stages run the full L0-L4 ladder.

| Layer | Category | Check | Catches |
|-------|----------|-------|---------|
| **L0** | Syntax | `python -m compileall src/radia/` | NameError / SyntaxError from recent edits |
| **L0** | Install intent | `os.path.samefile(src/radia/__init__.py, radia.__file__)` | Editable drift to non-editable copy |
| **L0** | Wheel manifest | `pyproject.toml` package-data lists every `*.jou/*.step/*.png/*.sol/*.vol/*.msh` in src tree | Samples/assets missing from wheel (2026-04-12 calc_mesh_eval.py incident) |
| **L1** | File layout | `import radia; os.path.isfile(...sample.jou)` | Samples / panel_registry.json missing after install |
| **L1** | Hash integrity | SHA-256 (CRLF normalised) src vs installed | Content drift (NOT size — LF↔CRLF false-alarms) |
| **L2** | Import smoke | import each `radia.<app>_notebook` workbench | ImportError / circular import in notebook workbench modules |
| **L2** | Subprocess smoke | `python panels/calc_*.py --help` for each | argparse rejects flags, missing helper imports |
| **L2** | Signature probe | `inspect.getsource(f)` contains expected token | Edit silently reverted by linter / merge |
| **L3** | Notebook QA | `pytest validation_test/panels/test_notebook_workbench.py` | Notebook contract drift (DesignSpec, result artifact, no-PySide rule) |
| **L3** | Cubit launcher | .ccl binary contains new regex strings (`grep -a`) | Old .ccl deployed; new TITLE/NEEDS_VOL not parsed |
| **L3** | Launcher widget matrix | Static source check: every mode-specific launcher widget has `setVisible(ms.needsVol)` | "Mesh order: 2" shown for STEP-only mode (2026-04-21 orderRow miss) |
| **L4** | End-to-end | Cubit -batch .jou → radia_export netgen → .vol | Full pipeline regresses (plugin missing, solver crash) |
| **L4** | Notebook smoke | Run the notebook workbench/headless calc path with real inputs → JSON | notebook workbench + calc script actually solve |
| **L4** | Panel mode matrix | For EACH active notebook mode, run its canonical sample E2E on 100号機 | New mode shipped without validation (2026-04-21 PEEC-inductance 4.78 nH regression) |
| **L4** | Golden-range numeric | Assert `L_coil_nH` / `P_total_W` / ... falls in a known range (NOT just exit 0) | Silently-wrong fallback path (e.g. wrong cross-section area → L 100× off) |
| **L4** | Multi-user ACL sanity | On 100号機 scan `C:\Users\*\AppData\Local\Coreform\Cubit\Coreform\licenses\renewals` and any `C:\temp\radia_*` file; owner MUST match the enclosing user.  Admin-owned files under a user AppData are the 2026-04-21 false-path (`debug-remote-user` SKILL) and block that user's Cubit entirely. | Broken-owner user-cache files from past Admin deploys that shipped but only affect non-admin users |

**Stage 1 (LAB)** runs L0-L3 on LAB only.
**Stage 2 (100号機 / mdx)** runs L0-L4 on each machine after `pip install
'radia[cubit]==X.Y.Z' radia-mcp==X.Y.Z cubit-mesh-export==X.Y.Z` from
PyPI + `cubit-plugin-install --all-users`. Both machines run the IDENTICAL
recipe.

Layer failures cascade: a broken L0 makes L1-L4 results meaningless.
Run the layers in order and stop on the first FAIL.

## Multi-User Hygiene (100号機, MANDATORY L4)

100号機 is a **shared multi-user Windows box** with ~20 local lab
accounts (keiko / kubota / yano / ...).  Any deploy step that writes
while running as Administrator **must not leave files in another
user's profile or in an ACL-restricted shared location**.  Two bugs
have now been observed:

1. **2026-04-21 false-path (skill `debug-remote-user`)**: admin copied
   `renewals` into each user's `AppData\Local\Coreform\Cubit\Coreform\
   licenses\` with owner `BUILTIN\Administrators`.  Non-admin users
   could not read their own cache → RLM -102 every Cubit launch →
   Cubit silently failed to start → "radia-ngsolve panel crashes"
   report.  Propagated to **8 of 9 non-admin users** before detection.

2. **2026-04-24**: `panel_log.py` hard-coded
   `PANEL_LOG_PATH = "C:/radia_panel_log.txt"`.  The root `C:\` grants
   `Users:ReadAndExecute` only, so the first user to run Cubit (Admin)
   owned the file and every non-admin subsequent write silently
   failed.  `init_panel_log(truncate=True)` rotation attempts (rename
   + remove) also failed silently.  Non-admin `C:\` C++ logs
   (`C:\radia_ccl_log.txt`, `C:\compact_netgen_debug.log`) had the
   same class of bug — fopen returned NULL for every non-admin Cubit
   load, so the diagnostic logs were Admin-only by accident.

### POLICY — shared-file rules on a multi-user Windows box

**Never do this**:
- Admin-side `cp` / Copy-Item / `New-Item` into another user's
  `%LOCALAPPDATA%` or `%USERPROFILE%`.  Owner will be
  `BUILTIN\Administrators` and UAC will deny the user read/write.
- Hardcode a shared log / cache path at `C:\` root.  Root ACL is
  `Users:ReadAndExecute` on standard installs.
- Hardcode a shared log in `C:\temp\foo.txt`.  `C:\temp` has
  `Users:CreateFiles` but **once Admin creates the file, the file
  inherits Admin owner**, so non-admin append will fail.

**Do this instead**:
- **Per-user filename** under `C:\temp\`:
  `C:\temp\radia_<thing>_<user>.txt`.  Each user creates their own,
  owns it, writes it.  `C:\temp` policy is AGENTS.md Temp Directory.
- **Per-user `%LOCALAPPDATA%`**: each user's Cubit / Python writes
  into their own `LOCALAPPDATA` in their own context.  Admin scripts
  **never** reach into another user's LOCALAPPDATA — use a logon-
  trigger Scheduled Task so the work runs as the target user.
- **Machine-wide under `C:\ProgramData\`** (rarely): only when the
  content is genuinely read-only reference data.  Explicit ACL
  grant `Users:ReadAndExecute`; **never** `Users:Modify` on a shared
  location.

### L4 sanity sweep — before and after every Stage 2 deploy

The same sweep script detects both bug classes above.  Run as
Administrator on 100号機.  **deploy is not "done" until this is
clean** (no BROKEN rows and no non-`C:\temp\*_<user>.*` hardcoded
paths found by grep).

```powershell
# (1) Broken-owner cache sweep.  Any AppData file whose owner is
#     Administrators / SYSTEM / another user is broken.
cat << 'PS' | ssh 100 'pwsh -ExecutionPolicy Bypass -Command -'
$candidates = @(
  'AppData\Local\Coreform\Cubit\Coreform\licenses\renewals'
)
$bad = @()
Get-ChildItem 'C:\Users' -Directory -EA SilentlyContinue |
  Where-Object { $_.Name -notin @('All Users','Default','Default User','Public') } |
  ForEach-Object {
    foreach ($rel in $candidates) {
      $p = Join-Path $_.FullName $rel
      if (-not (Test-Path $p)) { continue }
      $acl = Get-Acl $p
      $idx = $acl.Owner.LastIndexOf('\')
      $short = if ($idx -ge 0) { $acl.Owner.Substring($idx + 1) } else { $acl.Owner }
      if ($short -ine $_.Name) {
        $bad += [PSCustomObject]@{User=$_.Name;Path=$p;Owner=$acl.Owner}
      }
    }
  }
if ($bad.Count -gt 0) {
  Write-Output "BROKEN-OWNER files detected:"
  $bad | Format-Table -AutoSize
  Write-Output "Fix: takeown /F <path> /A ; icacls <path> /grant Administrators:F ; Remove-Item <path>"
  Write-Output "Then have each affected user log OFF and ON to rebuild with correct owner."
  exit 1
}
Write-Output "OK: all per-user caches have correct owner."
PS
```

```powershell
# (2) Hardcoded shared-root path sweep.  Any log / cache file
#     opened at C:\ root (not C:\temp\) is the 2026-04-24 bug class.
ssh 100 'grep -rnE ''"C:\\\\[A-Za-z_][A-Za-z0-9_]*\.(log|txt|json|cfg|ini|cache|tmp)"'' W:\00_CAE\Radia\01_GitHub\src' | grep -v 'C:\\\\temp\\\\'
# expected: no output.  Any hit must be converted to
# C:\temp\<name>_<user>.<ext> (per-user filename in the lab temp dir).
```

If either sweep finds a hit, **fail the deploy loudly** — this has
silently corrupted a large fraction of users twice in 4 days.

## Binary Propagation Policy

**POLICY**: 全バイナリ (.pyd, .ccm, .ccl) は **`src/radia/` と
`packages/cubit-mesh-export/src/cubit_mesh_export/` の 2 箇所** に
存在する必要がある (editable install + wheel + cubit-plugin-install
の全てがここから読む)。

### Source-of-truth ルール

| Phase                      | Source                          | Destination(s)                                       |
|----------------------------|--------------------------------|------------------------------------------------------|
| C++ build (Cubit plugin)   | `cmake --build` 出力 dir       | `src/radia/` + `packages/cubit-mesh-export/...`      |
| C++ build (radia core)     | `build-msvc/Release/*.pyd`     | `src/radia/`                                         |
| Python source              | `src/radia/**/*.py`            | (no copy — read in place via editable install)       |
| Wheel install (100号機)    | `pip install -e` snapshot     | `site-packages/radia/`                                |
| Plugin install (Cubit)     | `cubit-mesh-export` package    | `Coreform Cubit/bin/{plugins/}` via `cubit-plugin-install` |

### Build flow (CRITICAL)

```
1. Build.ps1                                      # MSVC + ninja
     -> build-msvc/Release/_radia_pybind.pyd
     -> src/cubit_plugin/build-pyd/radia_cubit_mesh.pyd
     -> src/cubit_plugin/build-ccm/radia_cubit.ccm
     -> src/cubit_plugin/build-ccm/radia_cubit.ccl   (Qt5 GUI)
2. Build.ps1 が src/radia/ + packages/cubit-mesh-export/ の両方にコピー
3. pip install -e .  +  pip install -e packages/cubit-mesh-export
4. cubit-plugin-install     # cubit-mesh-export -> Cubit/bin/
```

### よくある罠

1. **手動 cmake は別 dir** に出力する: `cmake --build build-cubit-plugin`
   は `build-cubit-plugin/Release/` を作るが Build.ps1 は
   `src/cubit_plugin/build-ccm/` を見ている → コピー先不整合。
2. **radia core build がコケると ccm/ccl もコピーされない** (Build.ps1
   は順番実行)。回避: 修正後の C++ 部分だけ手動で cmake → 出力を
   `src/cubit_plugin/build-ccm/` と `src/radia/` の両方に手動コピー。
3. **`pip install --force-reinstall --no-deps W:\...01_GitHub`** は
   source dir を **その時点の snapshot** で配信する。source dir 内の
   ファイルが古いと古いものが配信される。先に `src/radia/` を
   最新化してから pip install すること。

### SMB `cp` は禁止 (2026-04-21 policy)

「pip install を迂回して src/radia/*.py を site-packages に直接
 cp -r する」ショートカットは **絶対に使わない**。

**Why**:
- PyPI 配布で end-user が受け取る形態と違う → 実地で動くか判らない
- pip の RECORD / metadata / entry_points が壊れる (旧ファイルが残る)
- editable install が静かに non-editable に上書きされる事故 (PEEC-inductance 2026-04-21)
- wheel に入り損ねる package-data に気付けない

**The rule** (2026-05-02 updated):
- **Stage 1 (LAB)**: `pip install -e .` のみ
- **Stage 2 (100号機 / mdx)**: `pip install 'radia[cubit]==X.Y.Z'
  radia-mcp==X.Y.Z cubit-mesh-export==X.Y.Z` (PyPI) + `cubit-plugin-install
  --all-users` のみ.  両機 完全に同じ recipe.

ショートカットで動いても *distribution* を検証したことにはならない。
Stage 2 を 2 マシンで二重化することで PyPI 互換性 (wheel manifest /
package-data / entry-points / cubit-plugin-install regular-file
deploy) を release-qud Phase 9 cross-machine consistency probe で
ピン留めする。

## radia-mcp Policy

**LAB**: editable (`pip install -e packages/radia-mcp`).  ローカル
レポジトリの `packages/radia-mcp/src/radia_mcp/*/server.py` を直接参照。
**git push で即時反映**ではない (LAB の editable は src tree 直読み)。

**100号機 / mdx**: `pip install radia-mcp==X.Y.Z` (PyPI).  ローカル
レポジトリ不可。

## cub5_to_vol.py Policy

cub5_to_vol.py (Path B: Python 実装) は **Radia リポジトリには配布しない**。
netgen fork (ksugahar/netgen) 側で管理・サポートする。
Radia の配布 wheel には bridge.py, cubit_netgen_bridge.py, cub5_to_vol.py は含まれない。

## Deployment Strategy per Machine (2026-05-02 simplified)

| Component | LAB (ローカル) | 100号機 / mdx (PyPI 配布試験) |
|-----------|--------------|--------------------------|
| **radia** | `pip install -e .` | `pip install 'radia[cubit]==X.Y.Z'` |
| **cubit-mesh-export** | `pip install -e packages/cubit-mesh-export` | (`radia[cubit]` 経由で auto-install) |
| **radia-mcp** | `pip install -e packages/radia-mcp` | `pip install radia-mcp==X.Y.Z` |
| **PySide6** | 通常 Radia Python には入れない | Cubit 同梱の PySide6 は保護対象 (削除禁止) |
| **Cubit plugin** | (LAB Build.ps1 出力 = src/radia/ に直接) | `cubit-plugin-install --all-users` (regular-file copy from PyPI wheel) |

### Notebook migration: `[gui]` is NOT part of production deploy

Since the panel operating surfaces moved to Jupyter notebook workbenches,
production installs should not add the old GUI extra just to get PySide6.
Use `pip install radia[cubit]` for the Cubit/plugin path and verify notebooks
via `ipynb-gui-health`.

Do **not** uninstall or delete Coreform Cubit's bundled PySide6 under
`C:/Program Files/Coreform Cubit*/bin/python3/lib/site-packages`; Cubit owns
that embedded runtime.  The no-PySide target applies to the normal Radia Python
environment on LAB / 100号機 / mdx / hibino, not to Cubit's private Python.

### 100号機 / mdx PyPI install の詳細

100号機:
```
ssh 100 "pip install --upgrade --no-cache-dir 'radia[cubit]==X.Y.Z' radia-mcp==X.Y.Z cubit-mesh-export==X.Y.Z"
ssh 100 'cubit-plugin-install --all-users'
ssh 100 'cubit-plugin-install --verify-only'
```

mdx (identical recipe):
```
ssh mdx "pip install --upgrade --no-cache-dir 'radia[cubit]==X.Y.Z' radia-mcp==X.Y.Z cubit-mesh-export==X.Y.Z"
ssh mdx 'cubit-plugin-install --all-users'
ssh mdx 'cubit-plugin-install --verify-only'
```

新機能をラボメンバに届けるには **PyPI release が必要**:
1. LAB で fix → 動作確認 → 必要なら別 dev-clone / mdx PyPI で多重確認
2. `pyproject.toml` の version 上げ → git tag `vX.Y.Z`
3. `git push --tags` → CD が PyPI publish
4. PyPI 反映を `pip index versions radia` で確認 (~5 分)
5. 100号機 で `pip install --upgrade 'radia[cubit]==X.Y.Z' ...`
6. mdx も同じ recipe (cross-machine consistency probe)
7. kubota 等が 100号機 GUI から動作確認 (Step 8 panel QA)

## Remote Access

| Target | SSH alias | Default shell | Python | Build chain |
|--------|-----------|---------------|--------|-------------|
| LAB (this machine) | `ssh lab` | bash + pwsh | `python` | MSVC + Intel oneAPI + NGSolve dev (full) |
| mdx (163.220.176.228) | `ssh mdx` | **pwsh 7.6** | `python` | binary sync from LAB (no compile) |
| 100号機 (192.168.11.100) | `ssh 100` | **pwsh 7.6** | `python` | none — PyPI wheel only |

- All machines use **pwsh** (PowerShell 7). NOT PowerShell 5.
- SSH config: `~/.ssh/config` (aliases `100`, `mdx`, `lab`).
- File transfer to mdx: scp may be blocked by mdx host policy; use the
  base64 `[IO.File]::WriteAllBytes(...)` SSH-stdin pattern instead.
- File transfer to 100号機: SMB (`//192.168.11.100/c$/...`) or W: drive
  (LAB ↔ 100号機 share the W:\ NAS mount).
- `{}` in f-strings conflicts with pwsh ScriptBlock.  Use `heredoc | ssh
  host 'pwsh -Command -'` and prefer an `Out-File` to disk + `python
  C:/temp/<script>.py` to keep Python source out of the pwsh parser.

## Arguments

| Argument | Action |
|----------|--------|
| (empty / `all`) | Stage 1 + 2 |
| `local` | Stage 1 only |
| `100` or `100号機` | Stage 2 only |
| `release` | git push + tag |
| `mdx` or `228` | Stage 3 only |

---

## Stage 1: ローカル ビルド・検証

**POLICY**: 以下の手順は**この順序で不可分に実行**する。途中でスキップや順序変更しない。

```bash
cd S:/Radia/01_GitHub

# 1. Kill ALL first (Cubit が .ccl をロック、Python が .pyd をロック)
taskkill //F //IM coreform_cubit.exe
taskkill //F //IM python.exe
# WARNING: SSH-based Stop-Process can only see / kill processes in
# the SSH session, NOT the user's interactive desktop session.  On
# 100号機, a user with Cubit open in their RDP session will NOT be
# killed by this SSH call -- the silent failure looks like success.
# Use WMI for cross-session DETECTION and ask the user to manually
# close Cubit if it is running.
ssh 100 'pwsh -Command "
  $p = Get-CimInstance -ClassName Win32_Process -Filter \"Name = ''coreform_cubit.exe''\"
  if ($p) {
    Write-Host \"WARN: Cubit is running on 100号機 (PID=$($p.ProcessId)).\"
    Write-Host \"      SSH cannot kill cross-session processes.\"
    Write-Host \"      Please close Cubit MANUALLY before continuing the deploy.\"
  } else {
    Write-Host \"Cubit not running on 100号機 -- OK to proceed.\"
  }
"'

# 2. Build (C++ core + ccm + ccl + pyd → src/radia/ AND packages/cubit-mesh-export/)
pwsh -ExecutionPolicy Bypass -File Build.ps1 -Verbose

# 2b. Build the wheel (the artifact that Stage 2 / 3 install)
#     This is the SINGLE SOURCE OF TRUTH for downstream deploy.
#     Same wheel goes to 100号機 (via W:) AND to PyPI (via CI OIDC).
pwsh -ExecutionPolicy Bypass -File Build_Wheel.ps1 -DryRun
# -> dist/radia-X.Y.Z-cp312-cp312-win_amd64.whl
# Also build cubit-mesh-export wheel:
cd packages/cubit-mesh-export
pwsh -ExecutionPolicy Bypass -Command "pip wheel . --no-deps --no-build-isolation --no-cache-dir -w dist"
cd ../..
# -> packages/cubit-mesh-export/dist/cubit_mesh_export-X.Y.Z-*.whl

# 3. Editable install (pip が packages/ 内のバイナリを参照可能にする)
pip install -e . --no-deps
pip install -e packages/cubit-mesh-export --no-deps
# radia-mcp: NOT installed (direct repo reference)

# 4. Deploy plugin (packages/ 内の最新バイナリを Cubit にコピー)
cubit-plugin-install

# 5. Verify (basic import + version)
python -c "import radia; print(radia.__version__)"
python -c "import cubit_mesh_export; print(cubit_mesh_export.__version__)"

# 5b. CRITICAL: verify radia was installed in EDITABLE mode.
#     `pip install .` (non-editable) silently copies src/radia/ into
#     site-packages at install time -- subsequent `src/radia/` edits
#     become INVISIBLE to any code that does `import radia` (including
#     the Cubit .ccl launcher which calls `importlib.util.find_spec`).
#
#     Symptom: deploy looks green, but Cubit Solve menu doesn't show
#     new panel modes you just added (2026-04-21 PEEC-inductance
#     regression). Or panels/calc_*.py adds a relative path to
#     sys.path that ONLY resolves under the source-tree layout, then
#     fails silently in wheel-install layout (2026-05-01 kubota
#     PEEC-BEM `assemble_dense` regression -- bem_sibc_solver lived in
#     examples/ and was unreachable from site-packages/radia/panels/).
#
# HARD GATE (do NOT skip): two independent checks must both pass.
#   (i) `direct_url.json` exists and has `dir_info.editable == true`.
#       This is the canonical PEP 660 marker.  `pip install .` writes
#       `dir_info: {}` (empty) -- non-editable.  `pip install -e .`
#       writes `dir_info: {"editable": true}` -- editable.
#  (ii) `os.path.samefile(src_py, radia.__file__)` is True (defends
#       against an editable .pth that points at a stale clone).
# If EITHER check fails: `pip uninstall -y radia` then
# `pip install -e . --no-deps`.  Do NOT cp files around -- that
# violates Stage 1's editable contract.
python -c "
import radia, os, sys, json
import importlib.metadata as md

src_py = os.path.join(r'S:\Radia\01_GitHub\src\radia', '__init__.py')
sp_py  = radia.__file__
print(f'src: {src_py}')
print(f'sp : {sp_py}')

# (i) PEP 660 marker
try:
    raw = md.distribution('radia').read_text('direct_url.json')
    info = json.loads(raw or '{}')
    editable_flag = bool(info.get('dir_info', {}).get('editable'))
except Exception as e:
    print(f'ERROR: cannot read direct_url.json: {e}')
    sys.exit(1)
print(f'direct_url.json editable = {editable_flag}')
if not editable_flag:
    print('ERROR: radia was installed non-editable (pip install . instead of pip install -e .).')
    print('       Fix: pip uninstall -y radia ; pip install -e . --no-deps')
    sys.exit(1)

# (ii) samefile cross-check
try:
    ok = os.path.samefile(src_py, sp_py)
except OSError:
    ok = False
if not ok:
    print('ERROR: radia.__file__ does not resolve to src/radia/__init__.py.')
    print('       The editable .pth is pointing at a stale clone.')
    print('       Fix: pip uninstall -y radia ; pip install -e . --no-deps')
    sys.exit(1)

print('EDITABLE install OK (PEP 660 marker + samefile both pass)')
sys.exit(0)
"

# 6. Verify samples / GUI assets are reachable from the panel
python -c "
import radia, os
panels = os.path.join(os.path.dirname(radia.__file__), 'panels')
samples = os.path.join(panels, 'samples')
assert os.path.isfile(os.path.join(samples, 'ih_bem_sample.jou')), 'ih_bem_sample.jou missing'
assert os.path.isfile(os.path.join(samples, 'ih_fem_sample.jou')), 'ih_fem_sample.jou missing'
assert os.path.isfile(os.path.join(samples, 'ih_fem_kelvin_sample.jou')), 'ih_fem_kelvin_sample.jou missing'
assert os.path.isfile(os.path.join(samples, 'em_sample.jou')), 'em_sample.jou missing'
assert os.path.isfile(os.path.join(samples, 'pcb_sample.jou')), 'pcb_sample.jou missing'
assert os.path.isfile(os.path.join(panels, 'panel_registry.json')), 'panel_registry.json missing'
print('OK: panels/samples + panel_registry.json reachable from radia package')
"

# 7. Run verify-deploy: confirm Python actually loads the freshly-edited
#    src/radia/**/*.py files (mtime + signature checks). This catches
#    stale .pyc, shadow site-packages installs, and UNC vs S: drive
#    confusion. See `.Codex/skills/verify-deploy/SKILL.md`.
python -c "
import inspect, os
import importlib

# Phase 1: install state
import radia
print(f'  radia.__file__ = {radia.__file__}')

# Phase 2: mtime sanity
for name in ['radia.bem_coupled_solver', 'radia.notebook_workbench',
             'radia.ih_notebook', 'radia.em_notebook', 'radia.bem_inductance']:
    m = importlib.import_module(name)
    st = os.stat(m.__file__)
    print(f'  {name:<32} mtime={int(st.st_mtime)} size={st.st_size}')

# Phase 4: signature check (BEM coupled per-DOF f_back)
from radia.bem_coupled_solver import (
    CoupledBEMSolver, assemble_back_reaction_RHS)
src = inspect.getsource(assemble_back_reaction_RHS)
assert 'CoefficientFunction' in src and 'LinearForm' in src, \
    'BEM coupled per-DOF f_back not loaded'
print('  OK: BEM coupled per-DOF f_back is active')

# Phase 6: notebook workbench wiring
from radia.ih_notebook import IHWorkbench
from radia.em_notebook import EMWorkbench
from radia.pcb_notebook import PCBWorkbench
from radia.motor_notebook import MotorWorkbench
from radia.streamfunction_notebook import StreamFunctionWorkbench
for cls in (IHWorkbench, EMWorkbench, PCBWorkbench, MotorWorkbench, StreamFunctionWorkbench):
    assert hasattr(cls, 'build_command') or hasattr(cls, 'run_local'), cls.__name__
print('  OK: notebook workbench modules import')
"
```

```bash
# 7b. L0 — Syntax check all Python files (catches NameError / typo
#     in any radia_*.py / panels/*.py before runtime)
python -m compileall -q src/radia 2>&1 | tail -5
# Exit non-zero + stderr if any file fails.

# 7c. L0 — Wheel manifest audit (WARN; FAIL runs in Wheel Sanity Check)
#     Scope: panels/samples/ only.  The broader src tree has legacy
#     orphan .vol/.msh files (historical mesh debug outputs) that we
#     clean up at wheel-build time via exclude-package-data.
#     This narrow check catches: user adds a new .jou / .step /
#     .sol /.vol to panels/samples/ but forgets to extend the
#     pyproject.toml glob.  The full, repo-wide audit runs in the
#     release-gate "Wheel Sanity Check" section below.
python -c "
import fnmatch, sys, tomllib
from pathlib import Path
root = Path(r'S:\Radia\01_GitHub')
src = root / 'src' / 'radia'
samples = src / 'panels' / 'samples'
cfg = tomllib.loads((root / 'pyproject.toml').read_text(encoding='utf-8'))
pd = cfg.get('tool', {}).get('setuptools', {}).get('package-data', {}).get('radia', [])
excl = cfg.get('tool', {}).get('setuptools', {}).get('exclude-package-data', {}).get('radia', [])

USER_SAMPLES = ('.jou', '.step', '.stp', '.sol', '.py', '.json')
failed = []
for f in samples.rglob('*'):
    if not f.is_file(): continue
    if '__pycache__' in f.parts: continue
    if f.suffix not in USER_SAMPLES: continue
    rel = f.relative_to(src).as_posix()
    if any(fnmatch.fnmatch(rel, pat) for pat in pd): continue
    if any(fnmatch.fnmatch(rel, pat) for pat in excl): continue
    failed.append(rel)
if failed:
    print(f'ERROR: {len(failed)} panel sample file(s) would NOT ship in wheel:')
    for rel in failed[:15]:
        print(f'  {rel}')
    print('  Fix: add a glob to pyproject.toml [tool.setuptools.package-data].radia')
    sys.exit(1)
print(f'OK: all panels/samples/ user-visible files matched by package-data')
"

# 7d. L2 — Import smoke for every radia_*.py (catches ImportError
#     in new panel modules before the user hits "Solve")
python -c "
import importlib, sys
from pathlib import Path
failures = []
for f in sorted((Path('src/radia')).glob('radia_*.py')):
    name = 'radia.' + f.stem
    try:
        importlib.import_module(name)
        print(f'  {name:<40} OK')
    except Exception as e:
        print(f'  {name:<40} FAIL -- {e}')
        failures.append(name)
if failures:
    sys.exit(1)
"

# 7d2. L3 — Launcher widget visibility contract (static source audit)
#      Every QWidget row in the Radia-NGSolve launcher dialog must
#      either be ALWAYS VISIBLE (mode-independent) or gated by
#      `setVisible(ms.needsVol)` / `setVisible(ms.someFlag)`.  A row
#      with no setVisible call is treated as always-visible and may
#      show nonsense for STEP-only modes (e.g. "Mesh order: 2" when
#      no mesh is exported — 2026-04-21 orderRow miss).
python -c "
import re, sys
ccl = open(r'S:\Radia\01_GitHub\src\cubit_plugin\RadiaComp.cpp',
           encoding='utf-8').read()
# Find the Radia-NGSolve launcher function body (between 'Radia-NGSolve'
# QDialog and the end of its updateLabels lambda scope).
m = re.search(r'dlg\.setWindowTitle\(\"Radia-NGSolve\"\)(.*?)dlg\.exec\(\)',
              ccl, re.DOTALL)
if not m:
    print('WARN: could not locate Radia-NGSolve dialog body')
    sys.exit(0)
body = m.group(1)
# Enumerate QWidget* row containers declared in this scope
rows = re.findall(r'QWidget\s*\*\s*(\w+Row|\w+Group)\s*=', body)
# Each must have a setVisible(...) inside updateLabels
updateLabels = re.search(r'updateLabels\s*=\s*\[[^\]]*\]\(\)\s*\{(.*?)\};',
                          body, re.DOTALL)
gated = set(re.findall(r'(\w+)->setVisible', updateLabels.group(1))) if updateLabels else set()
missing = [r for r in rows if r not in gated]
print(f'Launcher row widgets: {rows}')
print(f'Rows gated by updateLabels: {sorted(gated)}')
if missing:
    print(f'WARN: these rows are never hidden (always visible):')
    for r in missing:
        print(f'  {r}')
    print('  If they should hide for STEP-only modes, add setVisible(ms.needsVol).')
else:
    print('OK: every launcher row is explicitly gated by updateLabels')
"

# 7e. L2 — Subprocess argparse smoke: run --help on every calc_*.py.
#     If a calc script raises at import time (missing helper, wrong
#     path), --help errors out instead of printing the flags.
pwsh -Command '
  $panels = "S:\Radia\01_GitHub\src\radia\panels"
  $fails = @()
  foreach ($f in Get-ChildItem "$panels\calc_*.py" -File) {
    $out = & python $f.FullName --help 2>&1
    if ($LASTEXITCODE -ne 0) {
      Write-Host "  $($f.Name) FAIL ($LASTEXITCODE)"
      Write-Host "    $($out | Select-Object -First 3)"
      $fails += $f.Name
    } else {
      Write-Host "  $($f.Name) OK"
    }
  }
  if ($fails.Count -gt 0) { exit 1 }
'

# 7f. L1 — Static probe for non-package imports in panels/calc_*.py.
#     `--help` runs the module top-level only; lazy imports inside
#     solve_*() bodies are NOT exercised.  The 2026-05-01 kubota
#     PEEC-BEM `assemble_dense` regression is the canonical case:
#     calc_peec_bem.py did `from bem_sibc_solver import ...` (bare name)
#     after sys.path.insert(BEM_REF) where BEM_REF was a relative
#     path that resolves under source-tree but NOT under
#     site-packages/radia/panels/ (the wheel-install layout).
#     Subprocess --help did not import bem_sibc_solver (lazy, inside
#     the solver function).  A path-resolution bug shipped silently.
#
#     RULE: every `from <X> import ...` in src/radia/panels/calc_*.py
#     where X is bare (no dot) MUST be one of:
#       (a) a stdlib module
#       (b) a sibling .py in the same panels/ dir
#       (c) a top-level dependency listed in pyproject.toml
#     ANY other bare-name import is rejected -- those rely on a
#     sys.path.insert() of a relative dir, which is the kubota
#     class.  Use `from radia.<name> import ...` instead.
python -c "
import ast, sys, tomllib
from pathlib import Path

repo = Path(r'S:\Radia\01_GitHub')
src_radia = repo / 'src' / 'radia'
panels = src_radia / 'panels'

# panels/calc_common.py adds src/radia/ to sys.path at module load.
# So bare 'from em_material import ...' resolves to radia/em_material.py.
# Whitelist: stdlib + sibling in panels/ + sibling in src/radia/ + dep.
panel_siblings = {p.stem for p in panels.glob('*.py') if p.stem != '__init__'}
radia_siblings = {p.stem for p in src_radia.glob('*.py') if p.stem != '__init__'}

# Stdlib (Python 3.12)
stdlib = set(sys.stdlib_module_names) | {'__future__'}

# pyproject.toml top-level deps
cfg = tomllib.loads((repo / 'pyproject.toml').read_text(encoding='utf-8'))
deps = []
for dep in cfg.get('project', {}).get('dependencies', []):
    deps.append(dep.split()[0].split('>')[0].split('=')[0].split('<')[0].strip())
for grp in cfg.get('project', {}).get('optional-dependencies', {}).values():
    for dep in grp:
        deps.append(dep.split()[0].split('>')[0].split('=')[0].split('<')[0].strip())
deps_top = {d.replace('-', '_').lower() for d in deps if d}
# Common indirect imports we know are present. PySide6 is intentionally absent
# from normal Radia Python after notebook panel migration.
deps_top |= {'mcp', 'cubit_mesh_export', 'pyvista', 'matplotlib',
             'scipy', 'pytest'}
deps_top_lower = {x.lower() for x in deps_top}

allowed = stdlib | panel_siblings | radia_siblings | deps_top_lower
fails = []
for f in sorted(panels.glob('calc_*.py')):
    src = f.read_text(encoding='utf-8')
    tree = ast.parse(src, filename=str(f))
    for node in ast.walk(tree):
        if not isinstance(node, ast.ImportFrom):
            continue
        mod = node.module or ''
        # `from X import` -- only top-level X matters here
        top = mod.split('.')[0] if mod else ''
        if not top:
            continue  # relative import 'from . import ...'
        if top.lower() in allowed:
            continue
        if top.startswith('_'):
            continue
        # Resolve path: bare top-level (no dot) AND not whitelisted
        if '.' not in mod:
            fails.append(f'{f.name}:{node.lineno}: from {mod} import ...')

if fails:
    print('ERROR: panels/calc_*.py has non-package bare imports.')
    print('  These work only in source-tree layout; they break in wheel install.')
    print('  Use `from radia.<name> import ...` or move the helper into src/radia/.')
    for line in fails:
        print(f'  {line}')
    sys.exit(1)
print(f'OK: panels/calc_*.py imports all resolve via stdlib / sibling / declared dep')
"

# 7g. L0 — Version consistency: pyproject.toml <-> __init__.py for ALL 3 packages.
#     The 2026-05-02 cubit-mesh-export 0.7.3 release failed at the CI
#     wheel-verify step because pyproject.toml was bumped to 0.7.3 but
#     __init__.py was left at 0.7.2 -- the CI caught it (the publish
#     step did NOT ship the broken wheel), but the developer had to
#     bump to 0.7.4 + retry the whole tag-CI-publish loop (~30 min).
#
#     Worse: the same scan found radia-mcp 0.36.6 had ALREADY shipped to
#     PyPI with a stale __init__.__version__ = "0.34.2" because the
#     release-radia-mcp.yml CI had no version-consistency check at all
#     -- so `python -c "import radia_mcp; print(radia_mcp.__version__)"`
#     on the published wheel returned "0.34.2" while pip shows 0.36.6.
#     This silent mismatch persisted for 1 release cycle before being
#     caught by an audit, and required a 0.36.7 patch release.
#
#     LESSON: don't rely on the per-package release CI to catch this.
#     Run the scan locally in Stage 1 BEFORE pushing the version-bump
#     commit -- catches the bug 0 minutes after introduction, no
#     wasted PyPI tags.
#
#     Scan rule: each package's pyproject.toml [project].version MUST
#     match its top-level __init__.py `__version__ = "..."` literal.
#     Failure here aborts the deploy.
python -c "
import re, sys, tomllib
from pathlib import Path
repo = Path(r'S:\Radia\01_GitHub')
packages = [
    ('radia',
     repo / 'pyproject.toml',
     repo / 'src' / 'radia' / '__init__.py'),
    ('cubit-mesh-export',
     repo / 'packages' / 'cubit-mesh-export' / 'pyproject.toml',
     repo / 'packages' / 'cubit-mesh-export' / 'src'
          / 'cubit_mesh_export' / '__init__.py'),
    ('radia-mcp',
     repo / 'packages' / 'radia-mcp' / 'pyproject.toml',
     repo / 'packages' / 'radia-mcp' / 'src' / 'radia_mcp' / '__init__.py'),
]
fails = []
for name, pyproject_path, init_path in packages:
    if not pyproject_path.is_file():
        fails.append(f'{name}: pyproject.toml missing at {pyproject_path}')
        continue
    cfg = tomllib.loads(pyproject_path.read_text(encoding='utf-8'))
    pyproject_ver = cfg.get('project', {}).get('version')
    if not pyproject_ver:
        fails.append(f'{name}: no [project].version in pyproject.toml')
        continue
    if not init_path.is_file():
        fails.append(f'{name}: __init__.py missing at {init_path}')
        continue
    init_text = init_path.read_text(encoding='utf-8')
    m = re.search(r'^__version__\s*=\s*[\"\\\']([^\"\\\']+)', init_text, re.M)
    if not m:
        fails.append(f'{name}: no __version__ literal in __init__.py')
        continue
    init_ver = m.group(1)
    if pyproject_ver != init_ver:
        fails.append(
            f'{name}: pyproject.toml=\"{pyproject_ver}\" '
            f'vs __init__.py=\"{init_ver}\" '
            f'-- bump both, not just one')
    else:
        print(f'  OK {name:25s} {pyproject_ver}')
if fails:
    print()
    print('ERROR: pyproject.toml <-> __init__.py version mismatch:')
    for f in fails:
        print(f'  {f}')
    print()
    print('Fix: bump BOTH the pyproject.toml [project].version AND the')
    print('     __init__.__version__ literal in lockstep, then re-run.')
    sys.exit(1)
print('OK: all 3 packages have consistent pyproject.toml / __init__.py versions')
"
```

**順序を守る理由**:
- 1 の前にビルドすると .pyd/.ccl がロックされてコピー失敗
- 3 の前に 4 を実行すると古いバイナリがデプロイされる
- 6 の verify を省略すると wheel install で `.jou` 不在に気付かない (panel が起動しても sample 一覧が空になる)
- 7b-g の L0/L1/L2 チェックは **コードレビューより先**に走る (static → run-time の順)
- 個別実行は禁止。必ず 1→2→3→4→5→5b→6→7→7b→7c→7d→7e→7f→7g→8 の順序

### Step 8 — Notebook Panel QA (MANDATORY when notebook panel code changed)

**POLICY**: `src/radia/*_notebook.py`, `src/radia/*_design.py`,
`src/radia/notebook_workbench.py`, `src/radia/panels/notebooks/*.ipynb`, or
`panel_notebook_manifest.json` を編集したら、**`ipynb-gui-health`** を正の
panel gate として使う。Jupyter notebook workbench が canonical surface であり、
normal Radia Python に PySide6 を追加して旧 desktop panel を検証対象に戻して
はいけない。

**自動チェック**:

```powershell
python -m pytest validation_test/panels/test_notebook_workbench.py -q
```

このテストは次を同時に確認する:
- active notebook が PySide6/PyQt に backslide していない
- `DesignSpec(...)` cell が初期値の source of truth で、JSON は run artifact
- `CommandWorkbench.run_local()` が `radia_result.v2` の `result.json` を残す
- `calc_*.py` への argv と notebook manifest が一致している

旧PySide panel QA は legacy adapter を直接修正した時の補助チェックに限る。
通常のリリース判定は notebook workbench と
`cubit-plugin-install --verify-only` / `cubit-smoke-test` で行う。

### Panel Samples Quality Gate (Stage 1, MANDATORY)

**POLICY** (2026-04-21): `src/radia/panels/samples/` は wheel に入って
end-user の Browse 初期値になる。**徹底的にデバッグした高品質なもの
だけ** を置くこと。新しい sample を `panels/samples/` に追加する場合、
Stage 1 で以下を通すまで deploy しない:

| Step | Action |
|------|--------|
| 1 | Cubit で `play <sample>.jou` が最後まで通る |
| 2 | 必要なら `radia_export netgen/.../gmsh` が正しい `.vol` / `.msh` を吐く |
| 3 | 対応するパネルの Run で期待値が出る (LAB 上で subprocess 実行) |
| 4 | `tests/panels/test_<feature>_golden.py` が PASS (数値を lock) |
| 5 | `pyproject.toml` の `package-data` に glob が入っている (wheel manifest audit 7c で確認) |

**WIP / 研究用 sample** は `examples/` か `tests/**/fixtures/` に置く。
`panels/samples/` に置くのは **release action** であって dev convenience
ではない。

**Why**: end-user は「壊れた sample」と「user error」を区別できない。
試して失敗 → 「このパネル使えない」という判断になり、**壊れた sample 1 つで
パネル全体の信頼が失われる**。PEEC-inductance 2026-04-21 で
`3turncoil.stp` が walker-hang 状態のまま samples に入っていて Kubota が
初回に遭遇した教訓。

### Wheel に含まれる GUI アセット (samples + registry)

Build.ps1 はバイナリだけでなく **GUI アセット** も wheel に詰める。
これらは pyproject.toml の `[tool.setuptools.package-data]` で
明示的に列挙する必要がある (setuptools は `*.py` 以外を自動収集しない):

```toml
[tool.setuptools.package-data]
radia = [
    "*.pyd", "*.ccm", "*.ccl", "*.png", "*.ico",
    "panels/samples/*.jou",      # IH/EM/PCB の例題ジャーナル
    "panels/panel_registry.json", # パネル定義 (radia_*.py が参照)
    "CONVENTIONS.md",
]
```

新しい例題 .jou を追加したら、上記リストを更新せずに wheel を
作っても **wheel には入らない**。`pip install` 後の Cubit パネルで
sample が選べなくなるので必ず確認すること。

確認コマンド (wheel build 後):
```bash
python -c "
import zipfile
z = zipfile.ZipFile('dist/radia-X.Y.Z-cp312-cp312-win_amd64.whl')
for n in sorted(z.namelist()):
    if 'samples' in n or 'panel_registry' in n:
        print(' ', n)
"
```

## Stage 2: 100号機 + mdx (PyPI install)

**POLICY (2026-05-02)**: The 3-tier model (LAB editable, mdx editable,
100号機 PyPI) collapsed to 2-tier: **LAB editable + 100号機/mdx PyPI**.
mdx editable is retired (2026-05-01 → 2026-05-02 lifespan).  Both
100号機 and mdx now run the IDENTICAL PyPI install recipe + Cubit
plugin regular-file deploy.  Verified on mdx 2026-05-02:
`cubit-plugin-install --verify-only` 5/5 [OK] sha256 match, Phase
C-heavy fixture L_coil = 138.159 nH.

### Stage 2 — 100号機 / mdx PyPI install (COMPLETE RECIPE 2026-05-02)

Run AFTER `pip index versions radia` confirms the new wheel is live.
Both 100号機 and mdx use exactly this recipe.  Substitute the SSH
target (`192.168.11.100` for 100号機, `mdx` for mdx).  Do not add the old
GUI extra for PySide6; notebook workbenches are the canonical panel surface.

```bash
TARGET=192.168.11.100   # or 'mdx'
VERSION=4.24.0          # version being released

cat << PS | ssh "$TARGET" 'pwsh -ExecutionPolicy Bypass -Command -'
\$ErrorActionPreference = 'Stop'

# Step 1: stop locks (python.exe / mcp-server-* / Cubit hold .pyd open).
Get-Process -ErrorAction SilentlyContinue | Where-Object {
  \$_.Name -like 'mcp-server*' -or \$_.ProcessName -eq 'coreform_cubit' -or \$_.ProcessName -eq 'cubit'
} | ForEach-Object { Stop-Process -Id \$_.Id -Force }
Start-Sleep -Seconds 2

# Step 2: tear down any prior install state.  uninstall first so
# editable-from-NAS or editable-from-clone leftovers do not shadow the
# PyPI wheel.
pip uninstall -y radia radia-mcp cubit-mesh-export 2>&1 | Select-String "Successfully|already|not installed"

# Step 3: PyPI install.  The cubit extra pulls cubit-mesh-export.
# Notebook panel QA is handled without installing PySide6 into normal Python.
pip install --no-cache-dir 'radia[cubit]==$VERSION' 'radia-mcp==$VERSION' 'cubit-mesh-export==$VERSION'

# Step 4: deploy Cubit plugin as regular files (NOT symlinks).
cubit-plugin-install --all-users
cubit-plugin-install --verify-only
PS
```

The `--verify-only` output should show:
```
[OK] radia_cubit.ccm    (NNNNNN bytes, sha256 match)
[OK] radia_cubit.ccl    (NNNNNN bytes, sha256 match)
[OK] radia_cubit_mesh.cp312-win_amd64.pyd  (NNNNNN bytes, sha256 match)
[OK] nglib.dll          (NNNNNN bytes, sha256 match)
[OK] ngcore.dll         (NNNNNN bytes, sha256 match)
compat: radia X.Y.Z <-> cubit-mesh-export A.B.C compatible
[OK] every expected binary present and matches package source.
```

End-to-end smoke (run on the same SSH target after install):

```bash
ssh "$TARGET" 'cubit-smoke-test'      # Cubit batch + radia_export netgen + .vol parse
ssh "$TARGET" 'python -c "
from radia.ih_notebook import IHWorkbench
from radia.em_notebook import EMWorkbench
from radia.pcb_notebook import PCBWorkbench
from radia.motor_notebook import MotorWorkbench
from radia.streamfunction_notebook import StreamFunctionWorkbench
for cls in (IHWorkbench, EMWorkbench, PCBWorkbench, MotorWorkbench, StreamFunctionWorkbench):
    print(cls.__name__, \"ok\")
"'
```

Both probes are mandatory before declaring the deploy complete.  On a source
checkout, also run `python -m pytest validation_test/panels/test_notebook_workbench.py -q`.
See
also `release-qud` Phase 8c (verify) and Phase 8d (smoke).

### Stage 2 — LEGACY (mdx editable, 2026-05-01 only) -- RETIRED 2026-05-02

**Critical points** -- if any of these is missed, the editable install
silently breaks (the 2026-05-01 incident took 3 hours to diagnose):

1. **mdx has NO MSVC, NO Intel MKL, NO `gh` CLI.**  C++ binaries must
   come from LAB via **base64-over-ssh push**, NOT from
   `download_binaries.sh` (gh-CLI-dependent).
2. **`cubit_mesh_export` legacy shadow** at
   `C:\Program Files\Python312\Lib\site-packages\cubit_mesh_export.py`
   (single-file legacy module) MUST be deleted before editable install
   -- pip uninstall does NOT touch it because it was placed outside
   pip metadata, and the shadow takes precedence over the editable
   package's `__init__.py`.
3. **Three packages, not two**: `radia`, `cubit-mesh-export`, and
   `radia-mcp`.  The earlier recipe missed `radia-mcp`.
4. **`radia_cubit_mesh.pyd` lives in TWO places** on the local clone --
   `src/radia/` AND `packages/cubit-mesh-export/src/cubit_mesh_export/`
   -- and BOTH must be populated.  LAB `Build.ps1` writes both; on
   mdx (no build) we push both via base64-over-ssh.
5. After the recipe, run **the full verification probe** which checks
   editable=True for all 3 packages, the C++ kernel imports, and the
   `radia_cubit_mesh` extension imports.

#### Step A — first-time setup OR full recovery on mdx

Run this from LAB.  Idempotent: re-runnable any time mdx drifts.

```powershell
# A1. Fix mdx-side shadows + install editable.
cat << 'PS' | ssh mdx 'pwsh -ExecutionPolicy Bypass -Command -'
# Stop locks
Get-Process -ErrorAction SilentlyContinue | Where-Object { $_.Name -like 'mcp-server*' } | ForEach-Object { Stop-Process -Id $_.Id -Force }
Start-Sleep -Seconds 2

# Delete legacy single-file shadow at site-packages (NOT touched by pip uninstall).
$shadow = "C:\Program Files\Python312\Lib\site-packages\cubit_mesh_export.py"
$shadowPyc = "C:\Program Files\Python312\Lib\site-packages\__pycache__\cubit_mesh_export.cpython-312.pyc"
if (Test-Path $shadow) { Remove-Item $shadow -Force; Write-Host "removed legacy shadow" }
if (Test-Path $shadowPyc) { Remove-Item $shadowPyc -Force; Write-Host "removed shadow .pyc" }

# Ensure local clone exists
$dest = 'C:\Radia\01_GitHub'
if (-not (Test-Path $dest)) {
    New-Item -ItemType Directory 'C:\Radia' -Force | Out-Null
    git clone https://github.com/ksugahar/Radia.git $dest
}
Set-Location $dest

# Uninstall any stale (editable or PyPI)
pip uninstall -y radia cubit-mesh-export radia-mcp 2>&1 | Select-String "Successfully|already|not installed"

# Install editable: ALL THREE packages
pip install --no-deps --no-cache-dir -e . 2>&1 | Select-String "Successfully"
pip install --no-deps --no-cache-dir -e packages/cubit-mesh-export 2>&1 | Select-String "Successfully"
pip install --no-deps --no-cache-dir -e packages/radia-mcp 2>&1 | Select-String "Successfully"
PS
```

#### Step B — push C++ binaries from LAB to mdx

Run this from LAB after every `Build.ps1` whose output mdx should see.

```bash
python tools/push_pyds_to_mdx.py
```

The tool:

* Stops any `python.exe` / `pythonw.exe` / `jupyter*` / `mcp-server-*`
  on mdx first (otherwise WriteAllBytes fails with `WinError 32: process
  cannot access the file because it is being [used by another process]`
  -- Python keeps the loaded `.pyd` open).
* Pushes every `.pyd` to its expected location(s).  Note that
  `radia_cubit_mesh.pyd` is mirrored to TWO paths -- `src/radia/` AND
  `packages/cubit-mesh-export/src/cubit_mesh_export/`.  The `src/radia/`
  copy is for `download_binaries.sh` parity (in case mdx ever gains the
  `gh` CLI); the `packages/.../cubit_mesh_export/` copy is what
  `cubit-mesh-export` editable install expects at import time.

Tool source: [`tools/push_pyds_to_mdx.py`](../../../tools/push_pyds_to_mdx.py).
The tool is destructive (terminates in-flight Python sessions on mdx);
mdx is a verification target only so this is acceptable.

#### Step C — day-to-day pull (after LAB git push)

```bash
ssh mdx 'git -C C:\Radia\01_GitHub fetch origin && git -C C:\Radia\01_GitHub checkout main && git -C C:\Radia\01_GitHub pull --ff-only'
# Then re-push C++ pyds (Step B) if Build.ps1 ran on LAB.
```

#### Step D — full editable verification

```bash
cat << 'PY' | ssh mdx 'python -'
import importlib.metadata as md
import os, json, math, sys

# 1. radia editable
import radia
expected_src = r"C:\Radia\01_GitHub\src\radia\__init__.py"
assert os.path.samefile(radia.__file__, expected_src), f"radia not from local clone: {radia.__file__}"
info = json.loads(md.distribution("radia").read_text("direct_url.json") or "{}")
assert info.get("dir_info", {}).get("editable"), "radia not editable"
print(f"[OK] radia v{radia.__version__} editable from {radia.__file__}")

# 2. _radia_pybind C++ kernel imports
from radia._radia_pybind import _average_demag_tensor, _AssembleSLDL_Galerkin
import numpy as np
A = _average_demag_tensor(np.array([0.,0.,0.]), np.array([1.,1.,1.]),
                          np.array([0.,0.,0.]), np.array([1.,1.,1.]))
assert abs(A.trace() + 1.0) < 1e-13, f"self-cube trace wrong: {A.trace()}"
print(f"[OK] _radia_pybind C++ kernel works (self-cube trace={A.trace():.4f})")

# 3. radia-mcp editable
print(f"[OK] radia-mcp v{md.version('radia-mcp')}")

# 4. cubit-mesh-export editable + C++ extension
import cubit_mesh_export
expected_cme = r"C:\Radia\01_GitHub\packages\cubit-mesh-export\src\cubit_mesh_export\__init__.py"
assert os.path.samefile(cubit_mesh_export.__file__, expected_cme), \
    f"cubit_mesh_export NOT from local clone (legacy shadow?): {cubit_mesh_export.__file__}"
from cubit_mesh_export import check, radia_cubit_mesh
print(f"[OK] cubit-mesh-export v{cubit_mesh_export.__version__} editable + radia_cubit_mesh.pyd loads")

# 5. v4.23.0 IH formula (regression check on the C++ binding)
from radia.analytical_formulas import cylinder_axial_eddy_loss
P = cylinder_axial_eddy_loss(1000.0, 0.02, 2*math.pi*10000, 5.8e7)
expected_P = 1.6121e+00
assert abs(P - expected_P) / expected_P < 0.01, f"IH formula drift: {P} vs {expected_P}"
print(f"[OK] IH cylinder Cu 20mm 10kHz = {P:.4e} W/m (expected ~1.6121)")

print("=== mdx editable verification ALL OK ===")
PY
```

#### Common failure modes (from 2026-05-01/02 incidents)

| Symptom | Root cause | Fix |
|---------|-----------|-----|
| `ImportError: cannot import name 'check' from 'cubit_mesh_export'` | legacy `cubit_mesh_export.py` shadow at site-packages | Step A removes it |
| `AttributeError: module 'cubit_mesh_export' has no attribute '__version__'` | same legacy shadow | Step A removes it |
| `ImportError: DLL load failed: ... radia_cubit_mesh.pyd` | `.pyd` only in `src/radia/` not in `packages/.../cubit_mesh_export/` (or vice versa) | Step B pushes BOTH |
| `radia.__version__` works but `from radia._radia_pybind import ...` fails | `_radia_pybind.pyd` missing or stale | Step B re-pushes |
| `from radia._radia_pybind import _AssembleSLDL_Galerkin` ImportError | mdx pyd predates a recent C++ change on LAB | Step B re-pushes after `Build.ps1` |
| `download_binaries.sh` reports success but no .pyds appear | mdx has no `gh` CLI; the script silently exits with no error | Don't use it on mdx -- use Step B instead |

---

### 100号機 の drive layout (legacy reference — Stage 2 was 100号機 until 2026-05-01)

100号機 **の中**で見たときの drive 状況:

| Drive  | Exists | Notes                                                |
|--------|--------|------------------------------------------------------|
| `S:\`  | NO     | LAB-only. 100号機からは見えない. **使うな**           |
| `W:\`  | YES    | 100号機ローカルドライブ. repo は `W:\00_CAE\Radia\01_GitHub` |
| `C:\`  | YES    | OS / Python / Coreform Cubit                         |

LAB から 100号機の C: にアクセスするときは SMB:
`//192.168.11.100/c$/...`

```bash
ssh 100 'Test-Path S:\'              # -> False (S: は LAB-only)
ssh 100 'Test-Path W:\00_CAE\...'    # -> True
```

**LAB 上で** 100号機の repo を直接見たい場合: SMB share `//192.168.11.100/work/00_CAE/Radia/01_GitHub` (= 100号機の `W:\00_CAE\Radia\01_GitHub`).

### Stage 2 — 100号機: editable + Cubit plugin symlinks (verify only, no install per deploy)

**POLICY (2026-04-24, user-set)**: 100号機 は **修正をすぐに反映する
editable** にしておく。一回 first-time editable install + Cubit
plugin symlink を貼ったあとは、LAB の編集 (Python) と LAB の
`Build.ps1` 出力 (.ccm / .ccl / .pyd) が即 100号機 で live になる。
deploy ごとに `pip install` も `cubit-plugin-install` も走らせない。
**かなり開発が進んだときだけ** PyPI に上げて、その配布形態の
検証は **mdx (Stage 3)** で行う。100号機 は PyPI 配布の試験対象では
ない。

This is the **all-users multi-user test** — the editable install
lives in the system-wide `C:\Program Files\Python312\Lib\site-packages\`
(direct_url.json points at `W:\…\src\radia`) which every local user
on 100号機 reads.  SMB `cp` shortcuts are still forbidden (see the
"SMB `cp` は禁止" section above).

#### (a) Pre-flight: FORCE-CLOSE all Cubit sessions

**POLICY**: Deploy takes priority.  Notify each Cubit user with a
30-second warning popup, then **force-kill** any remaining
`coreform_cubit.exe` via `Stop-Process -Id <pid> -Force` (cross-
session kill-by-PID works with admin rights — it's `-Name` that's
session-scoped).  Do NOT wait for user cooperation — `.ccl` must be
replaceable or the deploy is worthless.

User data loss risk: users lose unsaved Cubit state.  This is
acceptable because the skill emits the warning popup first, and
the lab agreement is that deploy is a critical operation.

```powershell
cat << 'PS' | ssh 100 'pwsh -Command -'
$p = @(Get-CimInstance Win32_Process -Filter "Name='coreform_cubit.exe'" -EA SilentlyContinue)
if ($p.Count -eq 0) {
  Write-Host "Cubit not running -- safe to deploy"
  exit 0
}
Write-Host "Cubit is running; force-close procedure:"
$sessions = @{}
foreach ($proc in $p) {
  $owner = ($proc.GetOwner()).User
  $sess  = $proc.SessionId
  $sessions[$sess] = $owner
  Write-Host ("  PID=$($proc.ProcessId)  User=$owner  Session=$sess")
}
# 1) warn every affected user with a popup (30-s countdown)
foreach ($kv in $sessions.GetEnumerator()) {
  try {
    msg $kv.Key /TIME:30 "[Radia deploy] Cubit will be force-closed in 30 seconds. SAVE YOUR WORK NOW."
    Write-Host "  notified session $($kv.Key) ($($kv.Value))"
  } catch {
    Write-Host "  notify FAILED for session $($kv.Key) -- msg may be disabled"
  }
}
# 2) wait 30 s for voluntary close
Start-Sleep -Seconds 30
# 3) force-kill anything still running (kill-by-PID works cross-session)
$still = @(Get-CimInstance Win32_Process -Filter "Name='coreform_cubit.exe'" -EA SilentlyContinue)
foreach ($proc in $still) {
  $pid_ = $proc.ProcessId
  try {
    Stop-Process -Id $pid_ -Force -EA Stop
    Write-Host "  force-killed PID $pid_"
  } catch {
    # Fallback via WMI Terminate (works when Stop-Process is blocked)
    $proc.Terminate() | Out-Null
    Write-Host "  force-killed PID $pid_ (WMI Terminate fallback)"
  }
}
# 4) verify all gone
Start-Sleep -Seconds 2
$left = @(Get-CimInstance Win32_Process -Filter "Name='coreform_cubit.exe'" -EA SilentlyContinue)
if ($left.Count -gt 0) {
  Write-Host "ERROR: $($left.Count) Cubit process(es) survived force-kill"
  exit 1
}
Write-Host "OK: all Cubit sessions closed, proceeding with deploy"
PS
```

If the user wants to avoid losing work, they see the 30-second
popup and Ctrl+S / File → Save in time.  After deploy, re-notify:

```powershell
ssh 100 'pwsh -Command "
  query session | Select-String \"Active\" | ForEach-Object {
    \$sid = ((\$_.Line -split \"\s+\") | Where-Object { \$_ -match \"^\d+\$\" })[0]
    if (\$sid) { msg \$sid /TIME:30 \"[Radia deploy] Complete. Safe to reopen Cubit.\" }
  }
"'
```

#### (b) Editable install + Cubit plugin symlinks (no pip install / no cubit-plugin-install)

**POLICY (2026-04-24)**: 100号機 is **editable + symlink**.  No wheel
install, no `cubit-plugin-install`.  LAB `Build.ps1` output (in
`src/radia/` and `packages/cubit-mesh-export/...`) becomes effective
on 100号機 immediately because:

- `pip install -e W:\00_CAE\Radia\01_GitHub` (radia) and
  `pip install -e W:\00_CAE\Radia\01_GitHub\packages\cubit-mesh-export`
  (cubit-mesh-export) are **one-time** setup in machine-wide
  `C:\Program Files\Python312\Lib\site-packages\`.  `import radia`
  resolves directly to `W:\…\src\radia\__init__.py`.
- `<Cubit>\bin\radia_cubit.ccl` /
  `<Cubit>\bin\plugins\radia_cubit.ccm` /
  `<Cubit>\bin\plugins\radia_cubit_mesh.cp312-win_amd64.pyd`
  are **SymbolicLinks** to the corresponding files on W:.  Cubit
  loads the live LAB build through them.

**DO NOT** run `cubit-plugin-install` on 100号機 — the installer copies
regular files over the symlinks and the editable-deploy chain breaks.
PyPI-style regular-file deploy is exercised on **mdx** (Stage 3), not
100号機.

**DO NOT** run `pip install` (wheel or source) on 100号機 unless one of
the following first-time / repair conditions applies:

| Condition | Action |
|-----------|--------|
| Brand-new 100号機 install | Run "First-time editable install" below. |
| `pyproject.toml` deps changed (NGSolve / MKL bump) | `pip install -e W:\…\01_GitHub --upgrade` |
| Editable drifted (samefile check fails in Step (c)) | `pip uninstall -y radia cubit-mesh-export; pip install -e W:\…` |
| `<Cubit>\bin\…\radia_cubit*` is a regular file (not SymbolicLink) | Re-create the 3 symlinks (see below). |

Otherwise the only Stage 2 work per deploy is **verification** that
the editable install + symlinks are still healthy.

##### First-time editable install (one-time setup, admin only)

```powershell
cat << 'PS' | ssh 100 'pwsh -ExecutionPolicy Bypass -Command -'
# 1. Editable install for ALL local users (machine-wide site-packages)
pip install --no-deps -e "W:\00_CAE\Radia\01_GitHub"
pip install --no-deps -e "W:\00_CAE\Radia\01_GitHub\packages\cubit-mesh-export"

# 2. Cubit plugin symlinks (admin token required for mklink)
$cubit = "C:\Program Files\Coreform Cubit 2025.3"
$repo  = "W:\00_CAE\Radia\01_GitHub"
# Remove any pre-existing regular file before mklink
Remove-Item -Force -EA SilentlyContinue `
    "$cubit\bin\radia_cubit.ccl", `
    "$cubit\bin\plugins\radia_cubit.ccm", `
    "$cubit\bin\plugins\radia_cubit_mesh.cp312-win_amd64.pyd"
# mklink (use cmd /c — pwsh New-Item -ItemType SymbolicLink to UNC paths is fragile)
cmd /c mklink "$cubit\bin\radia_cubit.ccl"                           "$repo\src\radia\radia_cubit.ccl"
cmd /c mklink "$cubit\bin\plugins\radia_cubit.ccm"                   "$repo\src\radia\radia_cubit.ccm"
cmd /c mklink "$cubit\bin\plugins\radia_cubit_mesh.cp312-win_amd64.pyd" `
              "$repo\packages\cubit-mesh-export\src\cubit_mesh_export\radia_cubit_mesh.pyd"
PS
```

ACL note: `W:\00_CAE\Radia\01_GitHub` already has
`NT AUTHORITY\Authenticated Users Modify`, so all 21 local users can
read the editable install + symlink targets.  Only Admin can do
`mklink`, so symlink maintenance always runs as Administrator via SSH.

##### Per-deploy verification (the normal path)

For an ordinary LAB build → 100号機 propagation, **no install or
plugin-deploy command runs**.  Just verify the editable + symlink
chain is intact, then go to Step (c) and beyond.  The
non-destructive preflight is `cubit-plugin-install --check-only`
(introduced in cubit-mesh-export v0.5.2 — does NOT touch any files):

```bash
cat << 'PS' | ssh 100 'pwsh -Command -'
"=== Cubit plugin: must be SymbolicLink ==="
$paths = @(
  "C:\Program Files\Coreform Cubit 2025.3\bin\radia_cubit.ccl",
  "C:\Program Files\Coreform Cubit 2025.3\bin\plugins\radia_cubit.ccm",
  "C:\Program Files\Coreform Cubit 2025.3\bin\plugins\radia_cubit_mesh.cp312-win_amd64.pyd"
)
$bad = @()
foreach ($p in $paths) {
  $f = Get-Item $p -EA SilentlyContinue
  if (-not $f) { $bad += "MISSING: $p"; continue }
  if ($f.LinkType -ne "SymbolicLink") { $bad += "REGULAR FILE (expected SymbolicLink): $p"; continue }
  $real = Get-Item $f.Target -EA SilentlyContinue
  if (-not $real) { $bad += "BROKEN LINK: $p -> $($f.Target)"; continue }
  Write-Host ("  OK  {0}" -f $f.Name)
  Write-Host ("      -> {0} ({1} bytes, {2})" -f $f.Target, $real.Length, $real.LastWriteTime)
}
if ($bad) { Write-Host "FAIL:"; $bad | ForEach-Object { Write-Host "  $_" }; exit 1 }

"=== Editable install: import radia must resolve to W:\..\src\radia ==="
python -c "import radia, os; assert os.path.samefile(radia.__file__, r'W:\00_CAE\Radia\01_GitHub\src\radia\__init__.py'); print('  OK editable: ' + radia.__file__)"
python -c "import cubit_mesh_export, os; assert os.path.samefile(cubit_mesh_export.__file__, r'W:\00_CAE\Radia\01_GitHub\packages\cubit-mesh-export\src\cubit_mesh_export\__init__.py'); print('  OK editable: ' + cubit_mesh_export.__file__)"

"=== cubit-plugin-install --check-only (NON-DESTRUCTIVE preflight) ==="
cubit-plugin-install --check-only
PS
```

Expected output: every symlink resolves, both `samefile()` asserts
pass, preflight prints `[OK] preflight clean`.  Any FAIL means the
editable + symlink chain has drifted — repair per the table above.

#### (c) Verify: editable-install samefile check (NOT hash, NOT size)

Stage 2 is **editable**, so `radia.__file__` and the source file in
`W:\…\src\radia\__init__.py` MUST resolve to the same inode (samefile
is True).  This is stronger than a hash compare — there is no second
copy of the file at all.

If samefile fails, the editable install was silently overwritten by a
non-editable `pip install` (the 2026-04-21 PEEC-inductance regression
class, and the 2026-05-01 kubota PEEC-BEM `assemble_dense` regression
where the wheel install put `panels/calc_peec_bem.py` into site-
packages with a `BEM_REF` relative path that resolved to a non-
existent location, masking the real version of `bem_sibc_solver`).
Repair: `pip uninstall -y radia cubit-mesh-export` then re-run the
first-time editable install in Step (b).

The HARD GATE is two-pronged: (i) PEP 660 marker
`direct_url.json::dir_info.editable == true`, AND (ii) `samefile()`
between import and source.  Either failing means a `pip install`
(non-`-e`) silently shipped a frozen copy.

```bash
cat << 'PY' | ssh 100 'python -'
import os, sys, json
import radia, cubit_mesh_export
import importlib.metadata as md
src_root = r"W:\00_CAE\Radia\01_GitHub\src\radia"
cme_root = r"W:\00_CAE\Radia\01_GitHub\packages\cubit-mesh-export\src\cubit_mesh_export"

print(f"radia v{radia.__version__}  loc={radia.__file__}")
print(f"cubit_mesh_export v{cubit_mesh_export.__version__}  loc={cubit_mesh_export.__file__}")

fails = []
def check_editable(pkg_name):
    """PEP 660 marker check: direct_url.json must say editable=True."""
    try:
        raw = md.distribution(pkg_name).read_text('direct_url.json')
        info = json.loads(raw or '{}')
        flag = bool(info.get('dir_info', {}).get('editable'))
    except Exception as e:
        print(f"  {pkg_name:<48} [NO direct_url.json: {e}]")
        fails.append(pkg_name + " editable marker")
        return
    tag = "OK" if flag else "DRIFT (non-editable)"
    print(f"  {pkg_name:<48} editable={flag} [{tag}]")
    if not flag:
        fails.append(pkg_name + " editable marker")

def check(import_path, src_path, label):
    try:
        ok = os.path.samefile(import_path, src_path)
    except OSError:
        ok = False
    tag = "OK" if ok else "DRIFT"
    print(f"  {label:<48} [{tag}]")
    if not ok:
        print(f"    import: {import_path}")
        print(f"    src   : {src_path}")
        fails.append(label)

check_editable("radia")
check_editable("cubit-mesh-export")
check(radia.__file__, os.path.join(src_root, "__init__.py"), "radia.__init__ is editable from W:")
check(cubit_mesh_export.__file__, os.path.join(cme_root, "__init__.py"), "cubit_mesh_export.__init__ is editable from W:")

# Spot-check a few panel files: if these aren't samefile, something
# is shadowing the editable install (rogue site-packages copy, conda env, etc.)
for rel in ["radia_ih.py", "radia_gui_base.py", "radia_peec_inductance.py",
            "panels/calc_peec_inductance.py", "panels/panel_registry.json"]:
    src = os.path.join(src_root, rel.replace("/", os.sep))
    # Compute the would-be import path: editable means same path as src
    check(src, src, f"src present: {rel}") if not os.path.isfile(src) else None
    if os.path.isfile(src):
        print(f"  src present: {rel:<40} ({os.path.getsize(src)} bytes)")

sys.exit(0 if not fails else 1)
PY
```

Any DRIFT line means `pip install -e` was silently overwritten or
shadowed.  Do NOT proceed — fix per the table in Step (b).

#### (d) End-to-end verify (same as before, check notebook QA works)

```bash
cat << 'PS' | ssh 100 'pwsh -Command -'
$env:QT_QPA_PLATFORM = "offscreen"
Set-Location 'W:\00_CAE\Radia\01_GitHub'
python -m pytest validation_test/panels/test_notebook_workbench.py -q
PS
```

Expected: notebook workbench contract PASS on 100号機.

#### (e) Post-deploy verification: "really deployed" HARD GATE

**Worst failure mode**: deploy looks green but isn't actually in
effect.  Old `.ccl` in a zombie Cubit process, wrong Python on
PATH, stale `.pyc` cache, pip install silently skipped.  These
produce silent wrong state that propagates into the next bug hunt.

Run these checks AFTER `pip install` + `cubit-plugin-install`.
Every check FAILs LOUDly — no warnings-treated-as-OK.

```powershell
cat << 'PS' | ssh 100 'pwsh -Command -'
$fail = @()

# --- Check 1: no zombie Cubit older than the .ccl ---
# Any Cubit process running NOW must have started AFTER the .ccl's
# LastWriteTime.  If one predates, it's running the old plugin.
$ccl = Get-Item "C:\Program Files\Coreform Cubit 2025.3\bin\radia_cubit.ccl"
$cubs = @(Get-CimInstance Win32_Process -Filter "Name='coreform_cubit.exe'" -EA SilentlyContinue)
foreach ($c in $cubs) {
  $t = ([Management.ManagementDateTimeConverter]::ToDateTime($c.CreationDate))
  if ($t -lt $ccl.LastWriteTime) {
    $fail += "ZOMBIE: PID $($c.ProcessId) started $t (before .ccl update at $($ccl.LastWriteTime))"
  }
}

# --- Check 2: .ccl on disk contains NEW strings ---
# grep for strings known to be in the latest .ccl build.  If a string
# is missing, either the build was skipped or Cubit/bin has stale .ccl.
$cclText = Get-Content -LiteralPath $ccl.FullName -Raw -Encoding Byte | ForEach-Object {
  [System.Text.Encoding]::ASCII.GetString($_)
}
# (ScriptBlock-safe: PowerShell joins byte output into one string)
$cclBin = [System.IO.File]::ReadAllBytes($ccl.FullName)
$cclStr = [System.Text.Encoding]::ASCII.GetString($cclBin)
$expected_strings = @("NEEDS_VOL", "STEP-only", "radia_*.py")
foreach ($s in $expected_strings) {
  if ($cclStr -notmatch [regex]::Escape($s)) {
    $fail += "CCL MISSING STRING: '$s' not found in $($ccl.FullName)"
  }
}

# --- Check 3: import radia resolves to W:\..\src\radia (editable) ---
# editable install → samefile against the W: source.  If samefile is
# False, the editable was overwritten by a non-editable `pip install .`
# (the 2026-04-21 PEEC-inductance class of regression).
$check3 = python -c "import os, radia; print('OK' if os.path.samefile(radia.__file__, r'W:\00_CAE\Radia\01_GitHub\src\radia\__init__.py') else 'FAIL: ' + radia.__file__)" 2>$null
if ($check3 -notmatch "^OK") { $fail += "WRONG RADIA: $check3" }
$check3b = python -c "import os, cubit_mesh_export; print('OK' if os.path.samefile(cubit_mesh_export.__file__, r'W:\00_CAE\Radia\01_GitHub\packages\cubit-mesh-export\src\cubit_mesh_export\__init__.py') else 'FAIL: ' + cubit_mesh_export.__file__)" 2>$null
if ($check3b -notmatch "^OK") { $fail += "WRONG CME: $check3b" }

# --- Check 4: every radia_*.py present in W:\..\src\radia ---
# editable install means W: IS the install location, so just confirm
# the source files exist (no separate site-packages copy to enumerate).
$src_on_nas = "W:\00_CAE\Radia\01_GitHub\src\radia"
$missing = @()
foreach ($f in Get-ChildItem "$src_on_nas\radia_*.py" -File -EA SilentlyContinue) {
  if (-not (Test-Path -LiteralPath $f.FullName)) {
    $missing += $f.Name
  }
}
if ($missing) { $fail += "MISSING FROM W:\src\radia: $($missing -join ', ')" }

# --- Check 5: clear .pyc cache so next import re-compiles ---
# Editable install drops __pycache__ next to the .py files on W:.
$pycache = Join-Path $src_on_nas "__pycache__"
if (Test-Path $pycache) {
  Remove-Item -Recurse -Force $pycache -EA SilentlyContinue
  Write-Host "  cleared $pycache (forces fresh bytecode on next import)"
}
# Also clear panels/__pycache__ (editable picks up .py from W:\..\panels too)
$pycache2 = Join-Path $src_on_nas "panels\__pycache__"
if (Test-Path $pycache2) {
  Remove-Item -Recurse -Force $pycache2 -EA SilentlyContinue
  Write-Host "  cleared $pycache2"
}

# --- Check 6: NEEDS_VOL on new modules is actually enforced ---
# For each radia_*.py with NEEDS_VOL = False, confirm Python's module
# exposes it -- catches "file on disk but module stale in memory".
$peec_ok = python -c "import radia.radia_peec_inductance as m; assert m.NEEDS_VOL is False; print('OK')" 2>&1
if ($peec_ok -notmatch "OK") {
  $fail += "PEEC-inductance module did not expose NEEDS_VOL=False: $peec_ok"
}

# --- Check 7: Cubit plugin is still SymbolicLink (cubit-plugin-install
#     was NOT accidentally run, which would replace symlinks with
#     regular-file copies and break editable-deploy of the C++ side) ---
foreach ($p in @(
    "C:\Program Files\Coreform Cubit 2025.3\bin\radia_cubit.ccl",
    "C:\Program Files\Coreform Cubit 2025.3\bin\plugins\radia_cubit.ccm",
    "C:\Program Files\Coreform Cubit 2025.3\bin\plugins\radia_cubit_mesh.cp312-win_amd64.pyd")) {
  $f = Get-Item $p -EA SilentlyContinue
  if (-not $f) { $fail += "PLUGIN MISSING: $p"; continue }
  if ($f.LinkType -ne "SymbolicLink") { $fail += "REGULAR FILE (must be SymbolicLink): $p" }
}

# --- Report ---
if ($fail.Count -gt 0) {
  Write-Host ""
  Write-Host "=== DEPLOY VERIFICATION FAILED ==="
  $fail | ForEach-Object { Write-Host "  $_" }
  Write-Host ""
  Write-Host "DEPLOY IS NOT IN EFFECT.  Do NOT consider this deploy shipped."
  exit 1
}
Write-Host "OK: post-deploy verification passed -- deploy is in effect."
PS
```

If Check 1 fails, force-kill the zombie (step `(a)` was insufficient;
rerun it).  If Check 2 fails, the `.ccl` on W: is stale — rerun
`Build.ps1` on LAB; no copy/redeploy step is needed because the file
on `<Cubit>\bin\` is a symlink to W:.  Check 3 catches a shadow
Python install on PATH or an editable that was silently overwritten
by a non-editable `pip install .` — repair via `pip uninstall -y
radia cubit-mesh-export` followed by the first-time editable install
in Step (b).  Check 4-6 catch missing source files / stale bytecode
/ module-state-out-of-sync.  Check 7 catches an accidental
`cubit-plugin-install` run that replaced symlinks with regular-file
copies — repair by re-creating the 3 symlinks per Step (b).

Do NOT declare the deploy complete until this step exits 0.

### End-to-end functional test on 100号機

**POLICY** (2026-04-21): **Every radia_*.py panel mode MUST run its
canonical sample end-to-end AND the output numeric value must pass a
golden-range check** before the deploy is considered shipped.  A
subprocess that exits 0 is NOT enough — silently-wrong results
(e.g. L off by 100x from a fallback path that accepted a bogus
cross-section area) pass the `status: "ok"` check but the user sees
garbage.

This closes the 2026-04-21 Kubota regression where my deploy
pre-check all returned green (imports, hashes, QA, .ccl strings,
NEEDS_VOL) but calc_peec_inductance had never been exercised on
any real input file.  The user ran the mode and got an
unhandleable native crash on his first STEP.

#### Panel mode matrix

| radia_*.py | calc_*.py | sample input | golden output range |
|-----------|-----------|--------------|---------------------|
| `radia_ih` (BEM) | `calc_heating_bem.py` | `ih_bem_sample.jou` → .vol | `P_total_W ≈ 11.55 ± 1` |
| `radia_ih` (FEM-SIBC) | `calc_fem_kelvin.py` | `ih_fem_kelvin_sample.jou` → .vol | L ≈ 84 nH ± 5 (Cu 50 kHz) |
| `radia_peec_inductance` | `calc_peec_inductance.py` | `ih_peec_inductance_coil.step` (simple torus) | L ≈ 75 – 90 nH (Cu 50 kHz, n_peri=16) |
| `radia_peec_inductance` (multi-turn) | `calc_peec_inductance.py` | `3turnCoil.jou` (explicit centerline) | L ≈ 400 – 500 nH (Cu 150 kHz, n_peri=16) |
| `radia_em` | `calc_em.py` | `em_sample.jou` → .vol | (add range here when first run) |
| `radia_pcb` | `calc_pcb.py` | `pcb_sample.jou` → .vol | (add range here when first run) |

**When a new notebook panel mode is added**, the deploy skill MUST be updated
to add a new matrix row BEFORE declaring the feature shipped.
`ipynb-gui-health` and the relevant `validation_test/` golden lane enforce
this.

Every row must run on 100号機 (Stage 2) AND on mdx (Stage 3).
Stage 1 (LAB) runs at least the rows for modes edited in this diff.

#### Sample E2E runner (scripted)

```bash
cat << 'PS' | ssh 100 'pwsh -Command -'
# Environment — editable install: $sp resolves to W:\..\src\radia
$sp = python -c "import os, radia; print(os.path.dirname(radia.__file__))"
$cubit = "C:\Program Files\Coreform Cubit 2025.3\bin\coreform_cubit.exe"
$tmp = "C:\temp"
if (-not (Test-Path $tmp)) { New-Item -Path $tmp -ItemType Directory | Out-Null }
$fails = @()

function Assert-InRange($value, $lo, $hi, $label) {
  if ($value -ge $lo -and $value -le $hi) {
    Write-Host ("  OK     {0}: {1:F3} (range [{2}, {3}])" -f $label, $value, $lo, $hi)
  } else {
    Write-Host ("  FAIL   {0}: {1:F3} NOT in [{2}, {3}]" -f $label, $value, $lo, $hi) -ForegroundColor Red
    $script:fails += "${label}: ${value} outside [${lo}, ${hi}]"
  }
}

# --- Row 1: PEEC-inductance, simple torus ---
Write-Host "`n=== PEEC-inductance: single-loop torus ==="
$step = Join-Path $sp 'panels\samples\ih_peec_inductance_coil.step'
if (-not (Test-Path $step)) {
  # Regenerate from .jou if the .step was not shipped (user regenerates on first play).
  $jou = Join-Path $sp 'panels\samples\ih_peec_inductance.jou'
  Copy-Item $jou $tmp\ih_peec_inductance.jou
  Push-Location $tmp
  & $cubit -batch -nographics -nojournal ih_peec_inductance.jou | Out-Null
  Pop-Location
  $step = Join-Path $tmp 'ih_peec_inductance_coil.step'
}
$out = python (Join-Path $sp 'panels\calc_peec_inductance.py') `
  --peec-step $step --peec-n-peri 16 --frequency 50000 `
  --current 1.0 --coil-sigma 5.8e7 2>&1
$json = ($out | Select-String '^\{' | Select-Object -Last 1).ToString()
$L_nH = ($json | ConvertFrom-Json).L_coil_nH
Assert-InRange $L_nH 75 90 "PEEC torus L_nH (50 kHz Cu)"

# --- Row 2: PEEC-inductance, multi-turn explicit-centerline .jou ---
Write-Host "`n=== PEEC-inductance: 3-turn explicit centerline ==="
$jou = Join-Path $sp 'panels\samples\3turnCoil.jou'
$out = python (Join-Path $sp 'panels\calc_peec_inductance.py') `
  --peec-step $jou --peec-n-peri 16 --frequency 150000 `
  --current 1.0 --coil-sigma 5.8e7 2>&1
$json = ($out | Select-String '^\{' | Select-Object -Last 1).ToString()
$L_nH = ($json | ConvertFrom-Json).L_coil_nH
Assert-InRange $L_nH 400 500 "PEEC 3-turn L_nH (150 kHz Cu)"

# --- Row 3: IH BEM sample (legacy) ---
Write-Host "`n=== IH BEM: gapped-torus sample ==="
$jou = Join-Path $sp 'panels\samples\ih_bem_sample.jou'
$wrapper = Join-Path $tmp 'test_ih.jou'
"playback `"$jou`"`r`nradia_export netgen `"$tmp/ih_test.vol`" order 1 overwrite`r`n" |
    Set-Content -Path $wrapper -Encoding ASCII
& $cubit -batch -nographics -nojournal $wrapper | Out-Null
$out = python (Join-Path $sp 'panels\calc_heating_bem.py') `
    --vol $tmp\ih_test.vol --wp-label workpiece `
    --coil-radius 0.110 --coil-current 100 --gap-deg 5 --frequency 50000 `
    --material steel --half-thickness 0.050 --esim-geometry cylinder `
    --impedance-model linear --mu-r 100 2>&1
$json = ($out | Select-String '^\{' | Select-Object -Last 1).ToString()
$P = ($json | ConvertFrom-Json).P_total_W
Assert-InRange $P 10 13 "IH BEM P_total_W"

# --- Row N: add new panel modes here ---

# --- Summary ---
if ($fails.Count -gt 0) {
  Write-Host "`n=== E2E PANEL MATRIX FAILED ===" -ForegroundColor Red
  $fails | ForEach-Object { Write-Host "  $_" -ForegroundColor Red }
  Write-Host "`nDEPLOY IS NOT IN EFFECT.  Panel output values outside golden range."
  exit 1
}
Write-Host "`nOK: all panel-mode E2E rows passed golden-range checks."
PS
```

**Golden-range principle** (2026-04-21):
- `status: "ok"` from a subprocess is NOT sufficient evidence the
  feature works.  L = 4.78 nH vs L = 426 nH both exit 0.
- Every mode asserts the numeric OUTPUT falls in a range that
  distinguishes correct physics from a silently-wrong fallback path.
- Ranges are WIDE (±10-20 %) to tolerate mesh / algorithm noise but
  tight enough that order-of-magnitude regressions (like the
  longest-edge section() bug) FAIL loudly.
- Updating a range requires a line in the commit message explaining
  what physics / algorithm change justifies the new range.

**Release (CI 通過後):**
```bash
git push origin main
gh run list --limit 3   # wait for CI success

# cubit-mesh-export (if changed)
git tag cubit-mesh-export-vX.Y.Z
git push origin cubit-mesh-export-vX.Y.Z

# radia (if changed)
git tag vX.Y.Z
git push origin vX.Y.Z

# radia-mcp (if changed)
git tag radia-mcp-vX.Y.Z
git push origin radia-mcp-vX.Y.Z
```

## Stage 3: 100号機 (PyPI install, multi-user lab seats) {#stage-3-100go-pypi}

**POLICY FLIP 2026-05-01**: This stage WAS "mdx PyPI" until 2026-05-01.
The roles were swapped: **100号機 now hosts the PyPI installation**
(`pip install radia==X.Y.Z`) so that kubota and other lab members exercise
the public distribution channel as their daily workflow, while **mdx** is
the editable iteration sandbox (Stage 2).  Cubit plugin on 100号機 now
goes through `cubit-plugin-install` (regular-file deploy from PyPI wheel),
NOT the symlink chain.

```powershell
# One-time on each release: install latest published version on 100号機
ssh 100 'C:/Progra~1/Python312/python.exe -m pip install --upgrade radia cubit-mesh-export radia-mcp'
ssh 100 'cubit-plugin-install'
```

The full L0–L4 ladder (see Check Matrix) runs against this PyPI install.
Any failure here BLOCKS the release for the lab — kubota / keiko / yano
all run this `pip install`.  The legacy mdx-PyPI instructions below are
kept as a historical reference only (mdx is now editable, see Stage 2).

### Legacy mdx PyPI instructions (HISTORICAL — pre-2026-05-01)

```bash
# (No longer used — mdx switched to editable.)
ssh mdx 'pip install --upgrade radia cubit-mesh-export radia-mcp'

# Deploy Cubit plugin + panels
ssh mdx 'cubit-plugin-install'

# Verify
cat << 'PY' | ssh mdx 'python -'
import os, radia
from cubit_mesh_export import __version__ as cm_ver
print(f"radia v{radia.__version__}")
print(f"cubit-mesh-export v{cm_ver}")
for name, path in [("ccm", r"C:\Program Files\Coreform Cubit 2025.3\bin\plugins\radia_cubit.ccm"),
                   ("ccl", r"C:\Program Files\Coreform Cubit 2025.3\bin\radia_cubit.ccl")]:
    sz = os.path.getsize(path) if os.path.isfile(path) else 0
    print(f"  {name}: {sz:,} bytes")
PY
```

**Pre-PyPI deploy** (wheel pipe, PyPI 公開前にテストしたい場合):
```bash
# cubit-mesh-export
cat packages/cubit-mesh-export/dist/cubit_mesh_export-*.whl | ssh mdx 'python -c "import sys; open(r\"C:\tmp\cm.whl\",\"wb\").write(sys.stdin.buffer.read())"'
ssh mdx 'pip install --force-reinstall C:\tmp\cm.whl'

# radia
cat dist/radia-*.whl | ssh mdx 'python -c "import sys; open(r\"C:\tmp\radia.whl\",\"wb\").write(sys.stdin.buffer.read())"'
ssh mdx 'pip install --force-reinstall --no-deps C:\tmp\radia.whl'

# Deploy plugin + panels
ssh mdx 'cubit-plugin-install'
```

---

## Wheel Sanity Check (run before any PyPI release)

`pip wheel` / `Build_Wheel.ps1` の出力を検証する。
orphan ファイル混入を防ぐ。

### radia wheel

```bash
cd S:/Radia/01_GitHub
rm -rf src/radia.egg-info
pip wheel . --no-deps --no-build-isolation --no-cache-dir -w C:/tmp/wheel_check

python -c "
import zipfile, glob, sys
whl = glob.glob('C:/tmp/wheel_check/*.whl')[0]
z = zipfile.ZipFile(whl)
names = z.namelist()
print(f'Wheel: {whl.split(\"/\")[-1]} ({len(names)} files)')

# 1. Must contain
must = ['radia/_radia_pybind.pyd', 'radia/radia_cubit.ccm', 'radia/radia_cubit.ccl',
        'radia/panels/panel_registry.json', 'radia/panels/calc_mesh_eval.py',
        'radia/panels/samples/ih_bem_sample.jou', 'radia/panels/samples/ih_bem_sample_B.sol',
        'radia/gmsh_post_export.py', 'radia/radia_gui_base.py']
ok = True
for m in must:
    found = any(m in n for n in names)
    if not found: print(f'  MISSING: {m}'); ok = False
if ok: print('  All required files present')

# 2. Must NOT contain (orphans)
bad = [n for n in names if ('/_' in n and n != 'radia/__init__.py'
       and not '/mcp_server/' in n  # __init__.py in mcp subpackages OK
       and not '/__pycache__/' in n)
       or '_pconv' in n]
# Filter out legitimate __init__.py
bad = [b for b in bad if not b.endswith('__init__.py')]
if bad:
    print(f'  ORPHAN FILES ({len(bad)}):')
    for b in bad: print(f'    {b}')
    ok = False
else:
    print('  No orphan files')

# 3. .jou count check (only public samples)
jous = [n for n in names if n.endswith('.jou')]
underscore_jous = [j for j in jous if '/_' in j]
if underscore_jous:
    print(f'  BAD: {len(underscore_jous)} debug _*.jou in wheel')
    ok = False
else:
    print(f'  .jou count: {len(jous)} (all public)')

sys.exit(0 if ok else 1)
"
```

### cubit-mesh-export wheel (MANDATORY: .pyd + .ccm + .ccl)

The Cubit plugin needs ALL THREE binaries.  PyPI 0.6.0 shipped without
`radia_cubit_mesh.pyd` because both CI workflows had a "log only" path
for missing pyd (2026-04-24 incident: mdx left in half-deployed state).
This local check is the third guard (CI build + CI release + local).

```bash
cd S:/Radia/01_GitHub/packages/cubit-mesh-export
rm -rf dist build src/cubit_mesh_export.egg-info
pwsh -Command "pip wheel . --no-deps --no-build-isolation --no-cache-dir -w dist"

python -c "
import zipfile, glob, sys, re
whl = glob.glob('dist/*.whl')[0]
print(f'Wheel: {whl}')

# Wheel filename must be platform-tagged (cp312-cp312-win_amd64), never
# py3-none-any — the latter implies pure-Python wheel and would silently
# accept a missing .pyd.  Same constraint as the CI build-test.yml retag.
fname = whl.replace('\\\\', '/').rsplit('/', 1)[-1]
m = re.match(r'cubit_mesh_export-([^-]+)-cp312-cp312-win_amd64\\.whl', fname)
if not m:
    print(f'  ERROR: wheel filename not cp312/win_amd64 platform-tagged: {fname}')
    sys.exit(1)
print(f'  OK platform tag: cp312/cp312/win_amd64 (version {m.group(1)})')

z = zipfile.ZipFile(whl)
names = z.namelist()

# All three plugin binaries are MANDATORY.  cubit-plugin-install does
# NOT recover from a missing one — it would either skip silently or
# error mid-deploy, leaving the Cubit/bin tree half-updated.
required = [
    'cubit_mesh_export/radia_cubit.ccm',
    'cubit_mesh_export/radia_cubit.ccl',
    'cubit_mesh_export/radia_cubit_mesh.pyd',
]
ok = True
for r in required:
    found = any(n.endswith(r) or n == r for n in names)
    if not found:
        print(f'  MISSING (refuse to publish): {r}'); ok = False
    else:
        info = z.getinfo(next(n for n in names if n.endswith(r) or n == r))
        print(f'  OK {r}: {info.file_size:>8} bytes')

# Console-scripts entry points must include the 3 commands
ep = z.read(next(n for n in names if n.endswith('entry_points.txt'))).decode()
for cmd in ('cubit-plugin-install', 'cubit-smoke-test', 'check-vol'):
    if cmd not in ep:
        print(f'  MISSING entry_point: {cmd}'); ok = False
    else:
        print(f'  OK entry_point: {cmd}')

sys.exit(0 if ok else 1)
"
```

**OK なら**: `gh run list` + tag push で PyPI 公開
**NG なら**: orphan を source tree / site-packages から削除し、egg-info を再生成

### Why this triple guard exists

The 2026-04-24 mdx audit found that PyPI cubit-mesh-export 0.6.0 was
missing `radia_cubit_mesh.pyd` despite being tagged `cp312-cp312-win_amd64`.
Three independent layers had a "log warning only" path that let a
broken wheel slip through:

1. **`build-test.yml`** copy step warned but did not exit on missing src .pyd
2. **`release-cubit-mesh-export.yml`** verify step logged "No .pyd in wheel (compact-netgen build path, expected in CI)" but published anyway
3. **`deploy/SKILL.md` Wheel Sanity Check** validated the radia wheel only — never inspected the cubit-mesh-export wheel

All three were tightened to fail hard on missing .pyd.  The "compact
netgen build path doesn't ship pyd" assumption was wrong — every
release path must bundle the pyd, no exceptions (AGENTS.md "No
Fallbacks — Fail Fast, Fail Loud").

---

## Notes

- **一般ユーザー**: `pip install radia[cubit]` + `cubit-plugin-install`
- **MCP (LAB/100号機)**: pip install 不要。Codex 設定で直接 `server.py` を指定。git push で即反映。
- **MCP (mdx/一般)**: `pip install radia-mcp`
- **cub5_to_vol.py**: Radia には配布しない。netgen fork (ksugahar/netgen) で管理。
- **Version bump**: `pyproject.toml` AND `__init__.py` (all 3 packages independently)
- **Qt5 SDK**: `src/cubit_plugin/qt5_sdk/` (repo-local, gitignored)
- **netgen fork DLLs**: vol-viewer (netgen GUI) にのみ必要。次回 PyPI netgen リリースで解消
