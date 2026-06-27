# build123d -> Netgen -> GMSH validation corpus

This directory holds the executable validation-class checks promoted from the
retired build123d/Netgen/GMSH examples topic.

The reusable implementation is `radia_mcp.build123d.pipeline`; these scripts
exercise geometry contracts, measurement rows, named region preservation, and
force/moment reference reductions. They are heavier than ordinary pytest unit
tests, so they live under `validation_test/` and write JSON sidecars next to the
scripts.

Run examples:

```powershell
python validation_test/build123d_netgen_gmsh_flow/validation_build123d_bbox_clearance_audit.py
python validation_test/build123d_netgen_gmsh_flow/validation_coaxial_region_stack.py --quick
python validation_test/build123d_netgen_gmsh_flow/validation_halbach_region_sweep.py --quick
python validation_test/build123d_netgen_gmsh_flow/validation_racetrack_plate_air_region.py --quick
```

The old docs-layer source archive has been retired. Use this directory for the
validation corpus and `radia_mcp.build123d.pipeline` for the reusable API.
