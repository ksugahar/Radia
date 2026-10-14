---
block: Winding Dictionary
library: radia_simulink_library
referenceBlock: radia_simulink_library/Coupling/Winding Dictionary
categories:
  - uncategorized
metadataQuality: high
policyStatus: approved
source: extracted-mask-description
---

# Winding Dictionary

## Summary

Compile winding names, .vol regions, turns, polarity, parallel paths, resistance, and circuit terminals to a fixed-wi.... From radia_simulink_library.

## Identity

- Library: radia_simulink_library
- ReferenceBlock: radia_simulink_library/Coupling/Winding Dictionary
- MaskType: Radia Winding Dictionary
- BlockType: SubSystem

## Use When

- user needs compile winding names, .vol regions, turns, polarity, parallel paths, resistance, and circuit terminals to...
- The user asks for a validated winding dictionary.

## Avoid When

- The user explicitly asks to construct logic from primitive blocks.
- The required behavior is outside the documented scope of this block.

## Inputs / Outputs

Unknown from extracted metadata.

## Notes

Prefer this block over constructing equivalent logic from primitives when the intent matches and the project policy allows library reuse.
