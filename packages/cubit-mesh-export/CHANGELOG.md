# Changelog

All notable changes to `cubit-mesh-export` — the high-order curved
mesh export package for Coreform Cubit (Netgen / GMSH / Nastran /
VTK / MEG / FEMEEM writers + Python bindings for consistency checks).

## 0.6.0 — Japanese / Unicode path support

Released 2026-04-22.

All 6 mesh exporters now correctly write `.vol` / `.msh` / `.bdf` /
`.vtk` / `.meg` / FEMEEM `in.dat` to paths containing non-ASCII
characters (Japanese, Korean, Greek etc.) on any Windows codepage.

Before: `radia_export netgen "C:/temp/日本語/coil.vol"` raised
  `"No mapping for the Unicode character exists in the target
   multi-byte code page."` and wrote no file.
After: same command writes the file (22,521 bytes on the reference
  sphere test).

### Implementation

New `src/cubit_plugin/utf8_path.hpp` provides one helper:

```cpp
std::filesystem::path u8_string_to_path(const std::string &s);
```

On Windows: `MultiByteToWideChar(CP_UTF8)` → `std::wstring` →
`std::filesystem::path(wstring)`.  Side-steps `std::string` →
`std::filesystem::path` implicit conversion (which uses the
system codepage = cp932 on Japanese Windows) and the previously-
used `CP_ACP` narrow-API pattern.

Applied to all 6 exporters:

* `ExportNetgenCommand.cpp` — `.vol` via `ng_mesh->Save()` + `.vol.json`
* `ExportGmshCommand.cpp` — `.msh` v4.1
* `ExportNastranCommand.cpp` — `.bdf`
* `ExportVtkCommand.cpp` — `.vtk`
* `ExportMegCommand.cpp` — `.meg` (FEMEEM / MAGIC)
* `ExportFemeemCommand.cpp` — `in.dat` + `node.dat` + ... (4 files)

### Smoke test harness

`smoke_test.py` driver wrapper now terminates with explicit
`exit 0` after `radia_export netgen`.  On slower boxes (100号機)
Cubit's headless teardown access-violates (exit code 0xC0000005)
before the mesh DB destructor flushes the `.vol` writer; the
explicit exit forces shutdown through the normal exit handler so
the file is closed first.  Removes a flake that had the smoke
test failing ~50 % on 100号機 while passing on LAB.

### Compatibility

- Requires `radia >= 4.7.0` (matching plugin binary bundle version).
- C++ plugin rebuilt; `radia_cubit.ccl` / `radia_cubit.ccm` /
  `radia_cubit_mesh.pyd` refreshed in the bundle.

## 0.5.x and earlier

See git log.  Key points: 0.5.0 split from `radia` as an independent
package; 0.5.6 moved compact_netgen snapshot in-tree; subsequent
patch releases iterated on NetgenCurver projection robustness.
