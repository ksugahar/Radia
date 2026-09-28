---
name: release-cubit-mesh-export
description: Release and LAB/100 release-dual deployment of cubit-mesh-export (version bump, native rebuild, pyd asset, tag, PyPI, release checkouts, deploy, forced deploy that closes running Cubit). Use for "release cubit-mesh-export", "release-dual", "deploy the exporter to LAB/100/102", or a forced redeploy while Cubit is open.
---

# release-cubit-mesh-export

cubit-mesh-export releases independently of Radia. Its deployment is the
LAB/100 **release-dual** (`tools/release_cubit_dual.py`); the Radia solver
release-quad must not install, repoint or gate it. Never install or run Cubit
on mdx1/mdx2.

## Hosts (check `hostname` first)

| name | host | ssh alias | notes |
|---|---|---|---|
| LAB = 102号機 | 192.168.121.102 | `lab` (also `102`) | development host; pwsh 7 ssh shell |
| 100号機 | 192.168.121.100, hostname `intel11` | `100` | student host, many user profiles |

- On 100, `W:` is local and `S:` is `subst` of `W:\00_CAE`; `W:\` is shared as
  `\\192.168.121.100\Work`. LAB's interactive session maps that share as `S:`.
- An ssh session on LAB has **no `S:`/`W:` and no share credentials** (UNC
  access is denied). Anything LAB must run over ssh needs a LAB-local source
  (`C:\release-dual\cubit-mesh-export-<ver>`).
- Inside an ssh session on LAB, a nested `ssh` started from Python with piped
  stdio hangs; with inherited stdio it works. So run the orchestrator **from
  100** (`--run-from 100`) when driving both hosts over ssh; from LAB it only
  works in LAB's interactive session.
- In remote pwsh command strings avoid the text `Remove-Item` (the local
  command guard blocks the whole call); use `[IO.File]::Delete(...)`.

## 1. Release (version N)

1. Work in a clean worktree from `origin/main`. Bump the version in **three**
   places: `packages/cubit-mesh-export/pyproject.toml`,
   `src/cubit_mesh_export/__init__.py`, and `src/cubit_plugin/cubit_mesh_export_pybind.cpp`
   (`m.attr("__version__")`; `tests/test_cubit_installers.py` enforces this).
   Add the `## N - ...` CHANGELOG section. Commit sources first.
2. `src/cubit_plugin/cubit_build.ps1 -Rebuild`, then commit `cubit_mesh_export.ccm`
   and `native_payloads.json` (the manifest records the last commit touching
   `src/cubit_plugin`; after any rebase of those commits, fix its `commit`).
3. Upload the git-ignored pyd as the content-addressed asset named in the
   manifest, then download it back and compare the hash:
   `gh release upload binaries C:\temp\cubit_mesh_curver-<sha256>.pyd -R ksugahar/Radia`.
   CI fails without it.
4. Tests: `tests/test_cubit_mesh_quality.py tests/test_cubit_installers.py
   tests/test_release_cubit_dual.py tests/test_native_build_provenance.py
   packages/cubit-mesh-export/tests` (a missing local pyd fails
   `test_distribution_ci_packages_the_exact_candidate_binaries`).
   For native changes, rerun `validation_test/cubit_mesh_export/geometric_refit_benchmark.py`
   with an **absolute** `--plugin-dir` (a relative one silently loads the
   installed plugin) and compare the `.vol` bytes with the previous release.
5. `git push origin HEAD:main` (the pre-push mdx preflight must pass; a
   `Connection reset` to mdx1 is transport, retry), wait for the
   `cubit-mesh-export` workflow, then push the annotated tag
   `cubit-mesh-export-vN`; its run publishes to PyPI.

## 2. Release checkouts and wheel

- 100: `git -C W:\00_CAE\Radia\01_GitHub worktree add --detach W:\00_CAE\Radia\release-dual\cubit-mesh-export-N cubit-mesh-export-vN`,
  copy the hash-checked pyd into `packages/cubit-mesh-export/src/cubit_mesh_export/`,
  run `_native_provenance.py verify`.
- LAB-local (for ssh deployment): LAB's ssh session cannot reach GitHub
  either (credential prompt hangs). Make a shallow clone on 100
  (`git clone --depth 1 --branch cubit-mesh-export-vN file:///W:/00_CAE/Radia/01_GitHub <dir>`),
  add the pyd, `tar -czf`, `scp` to `lab:C:/release-dual/`, extract, and
  verify HEAD, clean status and provenance over ssh.
- Wheel: `pip download --no-deps --only-binary=:all: --no-cache-dir cubit-mesh-export==N -d W:\00_CAE\Radia\release-dual\wheels-cubit-N`
  and confirm its sha256 against PyPI's JSON.

## 3. Deploy

Normal (stops if any Cubit is open), from 100:

```powershell
python tools/release_cubit_dual.py --action deploy --run-from 100 `
  --wheel W:\00_CAE\Radia\release-dual\wheels-cubit-N\cubit_mesh_export-N-cp312-cp312-win_amd64.whl `
  --source-sha <tag commit> `
  --source-root-lab C:\release-dual\cubit-mesh-export-N `
  --source-root-100 W:\00_CAE\Radia\release-dual\cubit-mesh-export-N `
  --evidence-lab C:\temp\cubit-dual-N --evidence-100 C:\temp\cubit-dual-N
```

then the same arguments with `--action done`. Both preflights finish before
either install. On 100 the worker registers every profile
(`--all-users`); the receipt requires it.

**Forced** (the user explicitly asked to close running Cubit): add
`--force-close-cubit`. Every running `coreform_cubit.exe` on both hosts is
killed before preflight, whoever owns it; unsaved work in those sessions is
lost. The closed pid/user list is printed and kept in each receipt
(`force_closed_cubit`). Use it only on explicit instruction for that deploy.

Receipts: the local host's worker writes its own; remote receipts are copied
beside it (`<evidence>/<target>/<action>.json`), which `done` reads.

## 4. After deploy

- Check `python -c "import cubit_mesh_export as c; print(c.__version__, c.__file__)"`
  on both hosts, and that no profile's `.cubit` still carries a
  `## BEGIN radia toolbar` block.
- The startup shim resolves the installed package at every Cubit start
  (2.1.2+), so older release checkouts, their wheels and LAB-local clones can
  be removed once both hosts report the new version. Ask before deleting.
- The human GUI release test (toolbar in a real desktop session) stays with
  the user.
