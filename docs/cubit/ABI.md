# Cubit Native ABI Boundary

Policy adopted 2026-11-07. Scope: `cubit-mesh-export` and every consumer of
Cubit output.

## Why the native plugin exists

High-order curving needs a geometry trampoline. `BuildCurvedElements` calls a
`CallbackGeometry` projection back into Cubit's ACIS kernel once per
high-order node, so that step must run in-process in C++. It is therefore
bound to the Cubit, ACIS, Netgen and compiler ABI, and a Cubit update can make
the `.ccm` fail to load or crash. Supporting Netgen/NGSolve curved meshes makes
this cost unavoidable.

## Rules

1. Keep the ABI-bound native surface minimal: mesh export and the curving
   trampoline only. Netgen is statically linked (`compact_netgen`), so no
   runtime `nglib.dll` mismatch is possible.
2. Everything else stays out of the native plugin. CAD-capture consumers and
   adapters for solvers that mesh from CAD themselves run outside Cubit as
   JSON-in/JSON-out Python bound to the capture SHA256. They need no curving
   trampoline and carry no Cubit ABI; the target mesher places its own
   high-order nodes.
3. In-Cubit Python is limited to short capture scripts through the public
   `cubit` module. A capture records the Cubit version and the API calls it
   used, so an API change after an update is detected rather than silent.
4. After every Cubit update, rebuild the plugin and rerun a curved-surface
   regression that records high-order node deviation from the CAD surface.
5. For CAD routes, verify the target's reconstructed radii and surfaces against
   the capture samples instead of relying on the trampoline.

## Not yet implemented (target)

- A small `extern "C"` function table at the native boundary.
- A load-time check of Cubit, ACIS and Netgen build identity that fails with a
  clear error instead of crashing. Until it exists, a mismatch can still crash
  Cubit (see "Interrupt Detected" in the export troubleshooting notes).
