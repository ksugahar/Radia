---
block: Material Dictionary
library: radia_simulink_library
referenceBlock: radia_simulink_library/Material Models/Material Dictionary
categories:
  - uncategorized
metadataQuality: high
policyStatus: approved
source: extracted-mask-description
---

# Material Dictionary

## Summary

Compile MATLAB material and region dictionaries against a Netgen .vol mesh. The output is a fixed-width numeric Bus f.... From radia_simulink_library.

## Identity

- Library: radia_simulink_library
- ReferenceBlock: radia_simulink_library/Material Models/Material Dictionary
- MaskType: Radia Material Dictionary
- BlockType: SubSystem

## Use When

- user needs compile matlab material and region dictionaries against a netgen .vol mesh. the output is a fixed-width nu...
- The user asks for a validated material dictionary.

## Avoid When

- The user explicitly asks to construct logic from primitive blocks.
- The required behavior is outside the documented scope of this block.

## Inputs / Outputs

Unknown from extracted metadata.

## Notes

Prefer this block over constructing equivalent logic from primitives when the intent matches and the project policy allows library reuse.
