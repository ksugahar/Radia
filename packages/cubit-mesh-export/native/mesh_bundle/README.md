# Experimental neutral native mesh acquisition

This optional C++ command acquires a linear volume mesh without importing
Python, a Python extension, or a downstream application into the CAD process.
It uses the licensed SDK's native interfaces. The source is separate from the
existing signed export payload; building it does not replace that payload.
It is not part of the currently released plugin installation.

Configure this directory with CMake, the current supported MSVC toolchain and
`Cubit_DIR` pointing to the installed SDK's `cmake` directory. Build the
`cubit_mesh_bundle` target. Use only a process-local plugin directory for
qualification. Installation, signing and distribution are separate operations.

Run a fresh batch/nographics process, not an attached GUI:

```text
export neutral_mesh "new-directory" units "mm"
```

Units are mandatory (`m` or `mm`); coordinates are not rescaled. The initial
profile supports explicit TETRA4/TET4 and HEX8 blocks, complete and unique volume
membership, and geometry-backed outward skin. Unsupported orders, missing or
ambiguous ownership, duplicate cells and unavailable orientation fail loudly.
Both geometry-associated and directly assigned blocks are expanded by the
native entity-list parser. Block IDs and names are preserved.

The output consists of `native-capture.json` and a same-mesh native Exodus
`alternate-mesh.e`. This is an acquisition, not an analysis-ready bundle.
An external Python process independently decodes Exodus, compares typed ordered
connectivity and coordinates, checks skin and Jacobians, writes and reparses a
placeholder BDF, and validates the versioned JSON schema before sealing the six
bundle members. `analysis_ready` and downstream verification remain false.
No downstream program is started, connected, configured or solved.

An explicit metadata-only profile preserves named exterior sidesets and
nodesets, including their source IDs, faces and nodes:

```text
export neutral_mesh "new-label-directory" units "mm" labels_only
```

The external `label_catalog` path checks those names and memberships against the
native Exodus file. It transfers no physical boundary values. The no-set sealed
bundle profile rejects nonempty sets; it never silently drops them. A consumer
must explicitly support the catalog or a future qualified set profile. Internal
or curve sidesets remain unsupported in this initial 3D profile.

The external validator has a direct script interface for a dedicated runtime:

```text
python mesh_bundle.py --capture native-capture.json --receipt execution-receipt.json --schema mesh-bundle.schema.json --output new-bundle
python mesh_bundle.py --capture native-capture.json --receipt execution-receipt.json --labels-only --output new-label-catalog.json
```

It requires NumPy/SciPy, pyNastran and jsonschema in that external runtime, not in
Cubit. Current independent qualification uses pyNastran 1.4.1 and jsonschema
4.25.1. Classic NetCDF Exodus is supported; unsupported encodings fail rather
than switching to a different mesh. Receipts are evidence, not authentication:
hashes and producer flags do not replace independent execution review.

The mesh contract schema contains no solver-specific settings. Materials in the
BDF are placeholders and must never be treated as an accepted physical model.
No interoperability with any particular downstream product is asserted.
