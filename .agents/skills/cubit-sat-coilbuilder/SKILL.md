---
name: cubit-sat-coilbuilder
description: "Audit ACIS SAT electromagnet CAD in Coreform Cubit and construct a verified CoilBuilder source from explicitly declared coil centerlines, turns, and currents. Use for SAT-to-CoilBuilder migrations; do not use it to infer windings blindly from solid bodies."
---

# Cubit SAT to CoilBuilder

## Purpose

Use this skill to migrate an electromagnet CAD assembly delivered as ACIS SAT
into a reproducible Radia `CoilBuilder` and `CoilBuilderHDivSource` model.
It turns the CAD into an auditable source definition without silently changing
the physical winding contract.

SAT represents solids.  It does **not** uniquely specify a conductor
centerline, number of turns, series connection, current direction, or intended
symmetry.  Treat a SAT coil body as evidence for geometry and a candidate for a
source region, never as sufficient information to generate a current path.

## Required Inputs

Before creating a source, collect and check all of the following:

- the original `.SAT` file, retained read-only;
- the CAD length unit and the conversion to metres;
- the Cubit version and imported model artifact;
- a mapping from CAD/Cubit coil volumes to named physical windings;
- an explicit, closed centerline for every electrical circuit;
- turns, signed current, conductor cross-section, and intended series/parallel
  connection; and
- independent acceptance probes: JMAG reference values where available, plus
  the physical locations and components to compare.

Use [references/designspec.md](references/designspec.md) as the source-model
contract.  A missing item is a blocking error, not a reason to guess.

## Workflow

### 1. Audit the SAT before import

Run the static audit without editing the CAD:

```powershell
python .agents/skills/cubit-sat-coilbuilder/scripts/audit_sat.py `
  "W:\ffag\model.SAT" `
  --output "C:\temp\model.sat-audit.json"
```

The report records the SHA-256, ACIS header, declared and parsed body counts,
and any SAT name attributes.  Its output is inventory evidence only; it is not
a material or winding classification.

For the legacy FFAG source, audit the chosen revision explicitly.  The current
latest FFAG SAT revision reports 18 SAT bodies.  This is an import starting
point, not proof that all 18 are physical iron or coils.

### 2. Import and map in Cubit

Import the SAT into the supported Cubit release.  Save the import journal and
export a `cubit_volume_inventory.json` containing, for every volume, the Cubit
volume id, SAT/CAD name, bounding box, and assigned role (`iron`, `coil`,
`void`, or `discard`).  Inspect all coil candidates visually.

Do not promote legacy JMAG external-air volumes to an HDiv FEM air domain.
For open-region HDiv calculations, retain the iron and current geometry and use
the Kelvin exterior treatment.  Keep an air gap only where it is a physical
local gap needing mesh resolution.

### 3. Declare, then construct, the coil source

Create an explicit `radia.sat-coilbuilder-design/v1` DesignSpec.  Build each
path in metres with `CoilBuilder`, and build the common source with
`CoilBuilderHDivSource.from_coilbuilders(...)`.

```python
from radia.accelerator_magnet_topopt import CoilBuilderHDivSource
from radia.coil_builder import CoilBuilder

coil = (
    CoilBuilder(current=design_current_a)
    .set_start(start_m, orientation=orientation)
    .set_cross_section(width_m, height_m)
    # Add only DesignSpec-declared straight and arc pieces here.
)
source = CoilBuilderHDivSource.from_coilbuilders(
    [coil], n_arc=200, closure_tolerance=1.0e-9
)
```

Use this `source` for the HDiv incident-field/RHS route.  Use
`source.to_radia_object()` only after its closure check passes when native
Radia tracking needs the equivalent finite-filament object.  Do not separately
discretize a solid-current coil for one path and a filament coil for the other.

### 4. Validate in increasing order of cost

1. Verify every path is continuous and closed at the declared tolerance.
2. Check the direct analytic source field at the DesignSpec probes.
3. Check the HDiv RHS / reconstructed field against the same direct source.
4. Compare `source.to_radia_object()` fields and native tracking with the same
   finite-segment source at the probes.
5. Compare the accepted full HDiv model with JMAG at declared probe locations
   and report components, norms, mesh order, and current convention.
6. For FFAG work, certify the periodic sector pairing and then compare the
   closed orbit by independent native tracking and direct-field tracking.

Store numerical results as validation JSON and publish the final, executed
result-bearing notebook.  A visual CAD overlay, an imported SAT, or a matching
coil bounding box is never a field-accuracy acceptance result.

## Prohibitions

- Do not infer turns, polarity, connectivity, or centerlines from a solid body.
- Do not overwrite the supplied SAT or hide its unit conversion.
- Do not use a legacy air box merely because it exists in a JMAG assembly.
- Do not accept a field comparison that silently falls back to another source
  model.
- Do not claim JMAG agreement until the probe locations, component convention,
  and both current conventions are recorded.
