---
block: Nonlinear HDiv-MMM Reactor
library: radia_simulink_library
referenceBlock: radia_simulink_library/Reduced Models/Nonlinear HDiv-MMM Reactor
categories:
  - signal-processing
metadataQuality: high
policyStatus: approved
source: extracted-mask-description
---

# Nonlinear HDiv-MMM Reactor

## Summary

Native nonlinear HDiv-MMM reactor. Input is winding current; outputs are terminal voltage, flux linkage, differential.... From radia_simulink_library.

## Identity

- Library: radia_simulink_library
- ReferenceBlock: radia_simulink_library/Reduced Models/Nonlinear HDiv-MMM Reactor
- MaskType: Radia Nonlinear HDiv-MMM Reactor
- BlockType: M-S-Function

## Use When

- user needs native nonlinear hdiv-mmm reactor. input is winding current; outputs are terminal voltage, flux linkage, d...
- The user asks for a validated nonlinear hdiv-mmm reactor.

## Avoid When

- The user explicitly asks to construct logic from primitive blocks.
- The required behavior is outside the documented scope of this block.

## Inputs / Outputs

Unknown from extracted metadata.

## Notes

Prefer this block over constructing equivalent logic from primitives when the intent matches and the project policy allows library reuse.
