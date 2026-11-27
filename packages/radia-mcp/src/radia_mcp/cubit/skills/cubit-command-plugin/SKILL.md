---
name: cubit-command-plugin
description: |
  Build, load and ship a native Coreform Cubit C++ command plugin: a DLL
  that adds new command verbs (for example `export <format> "file"`) to
  Cubit's command parser, so they run from the command line, a journal
  file or `cubit.cmd(...)`. Use when a Cubit feature needs in-process C++
  (mesh writers that walk the whole mesh, geometry callbacks into ACIS,
  journal-replayable verbs), when a plugin fails to load or crashes
  Cubit, or when a plugin must be rebuilt for a new Cubit release. Not
  for GUI-only buttons or dialogs: use the Python/PySide6 toolbar path
  (`cubit_toolbar_guide`) for those.
---

# cubit-command-plugin

Distilled from the radia-mcp Netgen `export netgen` native plugin
(`cubit_mesh_export.dll`, 2026) and verified against the Coreform Cubit
2026.8 SDK. The template under `templates/` builds and passes `smoke.ps1`
on Cubit 2026.8 with MSVC 14.44.35207, so start from it instead of from
memory.

Get the files with the `cubit_scaffold_command_plugin` tool, or copy
`templates/` by hand.

## 1. Decide first: does this need C++?

Use a C++ command plugin only when at least one holds:

- The verb must be journal-replayable and parsed like a built-in
  (`play x.jou`, `cubit.cmd('export foo "a.txt"')`).
- The work must run inside Cubit's process, for example walking the full
  mesh through `MeshExportInterface` or projecting nodes back onto ACIS
  geometry once per node.

Everything else belongs outside the plugin:

- GUI menus, buttons and dialogs: in-process Python/PySide6 toolbar.
- Format conversion, solver adapters and checks: ordinary Python,
  JSON in and JSON out, outside Cubit.

The plugin binds to the Cubit, ACIS and compiler ABI. Every Cubit
release can break it, so keep its surface as small as possible.

## 2. Anatomy (four pure virtuals plus one macro)

```
MyPlugin : CubitCommandInterface      get_keys(), create_command(key)
  CUBIT_PLUGIN(MyPlugin)              exports cubit_plugin_instance() + marker
HelloCommand : CubitCommand           get_syntax(), get_syntax_help(),
                                      get_help(), execute(CubitCommandData&)
```

- `CUBIT_PLUGIN(Class)` must appear in exactly one .cpp. It exports
  `cubit_plugin_instance()` and a string `cubit_plugin_verification;debug=no`.
  The loader scans the file for that string. Without the macro, the plugin
  is skipped silently. A Debug build carries `debug=yes` instead. Release
  Cubit does not support it, and mixing CRTs crashes Cubit.
- `get_keys()` returns internal keys only. `create_command(key)` returns
  `new` commands, which Cubit owns and deletes. Return `nullptr` for an
  unknown key.
- A command object persists between calls. Keep per-call state local to
  `execute()` or reset it there.
- `execute()` returns `false` on failure. Print the reason with
  `PRINT_ERROR` (from `CubitMessage.hpp`) first.
- Plugin verbs are parsed before Cubit's legacy handlers, so a verb that
  matches a built-in command masks it. Extending `export` with a new
  keyword (`export nodecount "f"`) is fine. Do not reuse a keyword that
  Cubit already exports, such as `stl`, `abaqus` or `exodus`.

## 3. Syntax strings (modified BNF)

| Form | Meaning |
|---|---|
| `<string:label='filename',help='<filename>'>` | required quoted string |
| `<value:label='order',help='<1 or 2>'>` | int or double |
| `[overwrite]` | optional keyword; neighbouring options in any order |
| `{dimension 2\|dimension 3}` / `{a\|b:label='x'}` | choice |
| `(...)`, `...` | group, repeat |
| also | `toggle`, `location`, `axis`, `plane`, `bounding_box` |

Read values in `execute()`:

```cpp
std::string f; data.get_string("filename", f);
int order = 1; data.get_value("order", order);    // false and unchanged if absent
bool ow = data.find_keyword("overwrite");
```

The complete grammar is documented in `<Cubit>/include/cubit/CubitCommandInterface.hpp`
(the `parserSyntax` page). The official examples are in `<Cubit>/sdk_examples/`
(CommandPlugin, CGMPlugin, GUIPlugin, CubitLib).

## 4. Reading the mesh

```cpp
std::unique_ptr<CubitBaseInterface, void(*)(CubitBaseInterface*)> raw(
    CubitInterface::get_interface("MeshExport"),
    [](CubitBaseInterface* p){ if (p) CubitInterface::release_interface(p); });
auto* iface = dynamic_cast<MeshExportInterface*>(raw.get());
if (!iface) { PRINT_ERROR("MeshExport interface is not available
"); return false; }
iface->set_use_sequential_ids(false);            // keep Cubit's ids (void)
if (!iface->initialize_export()) { PRINT_ERROR("initialize_export failed
"); return false; }
```

- Always release the interface. The RAII holder covers every early return.
- `get_coords(start, n, x, y, z, ids)` and
  `get_block_elements(start, n, block, types, handles)` return batches.
  Loop until they return 0.
- `BlockHandle` is a pointer (`ElementOutputGroup*`), not an id. Map it
  with `iface->id_from_handle(handle)`. The SDK only promises a list of
  blocks. With mesh that is not assigned to any block, the
  Netgen export plugin observed an extra group that is not a user
  block. Filter the handles against
  `CubitInterface::parse_cubit_list("block", "all")`.
- Check the bool returns of `initialize_export()` and
  `get_block_size()`. Otherwise a failure is reported as a successful
  undercount.
- For sidesets, use `CubitInterface::parse_cubit_list("sideset","all")`,
  `get_sideset_surfaces(id)` and `get_connectivity("tri", id)`. Free-mesh
  faces from Sculpt or Exodus have no parent surface, so read them
  directly.
- Cubit passes file names as UTF-8 `std::string`. On Japanese Windows,
  `std::filesystem::path(std::string)` decodes them as cp932, so the file
  fails to open. Build the path through UTF-16 with `utf8_path.hpp`, and
  reject invalid UTF-8 instead of falling back to the code page.
  `smoke.ps1` writes into a `日本語` directory to prove this.
- Handle `overwrite` yourself. Without it, open with an exclusive create
  (`_wfopen_s(..., L"wbx")`) rather than "check exists, then truncate".
  Check the write result and the result of `fclose`.

## 5. Build (Windows)

`templates/build.ps1` does all of this and fails on anything else:

- Release only. Cubit is Release, and Debug/Release CRT mixing crashes it.
- Use a pinned MSVC toolset (default 14.44.35207), selected through
  `vcvars64.bat -vcvars_ver=` and checked against `cl.exe`. There is no
  silent fallback.
- Use the CMake and Ninja that ship with VS2022. Pass
  `-DCubit_DIR=<Cubit>/cmake`, then call `find_package(Cubit REQUIRED CONFIG)`
  and link `cubiti cubit_util`. Read the release from `<Cubit>/BUILD_INFO`
  (`short=`), never from the directory name.
- Define `_DISABLE_CONSTEXPR_MUTEX_CONSTRUCTOR` for the code built into
  the plugin. Cubit loads its own bundled `MSVCP140.dll`. Microsoft
  supports a runtime that is at least as new as the toolset only; when
  the bundled runtime is older, a constexpr-constructed `std::mutex` from
  MSVC 14.40 or later can crash in `lock()`. The Netgen export plugin hit this
  with Cubit 2025.12, which bundles 14.31. Cubit 2026.8 bundles 14.42,
  so with toolset 14.44 this define is a targeted workaround, not a
  guaranteed crash. Check the bundled runtime's version
  (`<Cubit>/bin/MSVCP140.dll`) when you change the toolset.
- Use `/Zi` and `/DEBUG:FULL` in Release, plus
  `/PDBALTPATH:<name>.pdb`, so crash dumps can be symbolicated without
  leaking the build path.
- If the plugin needs geometry, link `cubit_geom`. In 2026.8 this is a
  split target: `find_package(CubitCGM)` plus an alias to
  `cubit_cgmcore`. Delay-load it (`/DELAYLOAD:cubit_geom.dll` plus
  `delayimp`) so the plugin loads before the geometry engine.
- Large vendored sources (for example a Netgen subset) need `/bigobj`.
  Use `/utf-8` for non-ASCII comments.
- "Could NOT find CubitGeom" during configure is harmless for a plugin
  that does not use geometry.

## 6. Loading

- Cubit loads plugins only at startup, from `CUBIT_PLUGIN_DIR` or from
  the install's plugin directory. In 2026.8 that directory is
  `<install>/plugins`, beside `bin/` (`bin/cubitx.exe` exists). Older
  releases use `bin/plugins`.
- The 2026.8 SDK example builds a normal MODULE (`.dll`), and the
  radia-mcp installer (`python -m radia_mcp.cubit.install`) deploys the
  Netgen export plugin as `.dll` on 2026.8.
  Older releases also loaded `.ccm` files and took `-commandplugindir`
  (Coreform forum, "Example plugin not loading", 2023). Whether 2026.8
  still accepts `.ccm` has not been tested, so deploy `.dll`. A package
  can keep `.ccm` as its canonical name and copy it as `.dll` when
  deploying.
- Headless runs: 2026.8 uses `cubitx.exe -nographics -batch -nojournal
  -noinitfile driver.jou` with `CUBIT_PLUGIN_DIR=<staged dir>`. Older
  releases use `coreform_cubit.com ... -commandplugindir <dir>`.
  - Close stdin. Under an MCP stdio server, stdin is the JSON-RPC pipe,
    and Cubit waits on it.
  - `-noinitfile` keeps the user's init file out of the run.
- Python started with `cubit.init([...])` (the in-process daemon) does
  not load C++ command plugins. Run plugin verbs through the console
  launcher. cubit-mcp's `session.run_headless_journal(commands,
  command_plugin_directory=dir)` does the staging for you.
- Test a rebuilt plugin from a staged copy, never by overwriting the
  deployed DLL. Another user's Cubit may have it open, and Windows locks
  loaded DLLs.
- Rename hazard: when a plugin file is renamed, delete the old file.
  Otherwise Cubit loads both and registers the verbs twice.
- Observed symptom: with a stale machine-wide `CUBIT_PLUGIN_DIR`, Cubit
  exited with code 2 and reported `-commandplugindir` and its value as
  files it could not open. That is one cause of exit code 2, not the
  only one. Read the log, fix the cause, and never accept exit code 2 as
  success. `smoke.ps1` requires exit code 0.

## 7. Verify

1. `build.ps1`: checks the marker (`debug=no`) inside the DLL and prints
   its SHA256.
2. `smoke.ps1` (Cubit 2026.8 only): stages the DLL, meshes a brick, and
   runs both verbs. It passes only when the exit code is 0, both exact
   log lines appear, and the output file has exactly the expected counts.
   The whole run, including the output pipes, is bounded by the timeout.
   On failure it keeps the work directory. Expected output:
   `SMOKE OK exit=0 hello=2 MYPLUGIN_EXPORT nodes=125 blocks=1
   elements=64` (brick x 1, size 0.25).
3. For a real writer, also read its output back with an independent
   reader (for example `python -m radia.vol_check` for Netgen `.vol`). A file that
   exists is not proof that the export is correct.

## 8. Ship

- Build only from a clean commit that has been pushed. Record the
  commit, the toolset, the Cubit `BUILD_INFO` and the SHA256 of each
  binary next to the payload.
- The installer copies the payload, compares SHA256 values after
  copying, deletes files with retired names, and offers a
  `--verify-only` mode.
- Support exactly one Cubit release. Moving to a new release replaces
  the supported version (rebuild and rerun the regression); it does not
  add a second one.
- Debug logging: write to a per-user file (for example
  `C:\temp\<plugin>_%USERNAME%.log`). A shared file at `C:\` is not
  writable for non-admin users.

## Failure table

| Symptom | Cause |
|---|---|
| Verb unknown, nothing in log | `CUBIT_PLUGIN` missing, wrong directory, not a `.dll` on 2026.8, Cubit not restarted, or the run used `cubit.init` Python |
| Plugin missing or crashing from a Debug build | marker says `debug=yes`; Debug CRT against Release Cubit |
| Crash in `std::mutex::lock` | toolset newer than Cubit's bundled `MSVCP140.dll` (>= 14.40 vs older runtime) without `_DISABLE_CONSTEXPR_MUTEX_CONSTRUCTOR` |
| Crash on load after a Cubit update | ABI drift: rebuild against the new SDK |
| Verbs registered twice | an old renamed plugin file is still in the plugin directory |
| Built-in command stops working | the plugin syntax masks it; rename the keyword |
| Headless run hangs | stdin left open under an MCP or pipe parent |
| Exit code 2 with `-commandplugindir` "files" | observed with a stale machine-wide `CUBIT_PLUGIN_DIR`; check the log for other causes |
| Open fails only on Japanese paths | path built from a UTF-8 `std::string` without the UTF-16 conversion |
