---
name: kelvin-periodic-contract-verify
description: Verify persisted Kelvin periodic metadata in an exported VOL without changing the mesh, guessing a translation or repairing point identifications.
---

# Verify the producer's periodic contract

Use this for a mesh whose exporter was explicitly asked to write the named
`kt` group. The producer owns matching and persistence; a consumer only reads.
Do not call `IdentifyPeriodicBoundaries`, `AddPointIdentification`, `Curve`,
or a post-hoc helper to make a missing or inconsistent contract pass.

1. Bind the candidate file hash, producer identity, explicit master/slave BC
   patch map, translation, coordinate unit and approved geometric matching
   parameter. Obtain the expected pair set from an independent reference.
   An oracle table tied to a different mesh or numbering is unusable.
2. Load the file read-only in the qualified Netgen/NGSolve runtime. Require
   `GetNrIdentifications() == 1`, every identification triple's number to be
   `1`, type `IdentificationType.PERIODIC`, and the exact expected directed
   point-pair set. Check the actual persisted `identificationnames` section
   contains exactly one name, `kt`; do not invent a Python name getter.
3. Distinguish namespaces: a surface element references a FaceDescriptor;
   that descriptor's BC property gives the boundary-condition number. Neither
   is a CAD surface ID. Verify the producer's explicit CAD/FD/BC mapping.
4. Verify point numbering/connectivity and curved coordinates/Jacobians using
   the frozen evaluation rule and exact reference hashes. Equal element counts
   do not establish geometry identity. Do not curve or remesh after loading.
5. Report structural verification independently from field accuracy/adoption.
   Missing names, extra groups, incomplete or ambiguous pairs and geometry
   differences fail. Preserve the evidence and fix the producer; never widen
   a matching distance or silently reconstruct identifications in the consumer.

The native `export netgen` command requires all periodic fields explicitly.
Ordinary exports without that request do not imply periodic solver readiness.
Native builds and licensed exports use their separately admitted resource frame.
