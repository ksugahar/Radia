# Signal Processing Blocks

Use these blocks for signal processing blocks.

## Recommended Blocks

### Motor Angle Family

- Block: [[blocks/Motor_Angle_Family]]
- Library: radia_simulink_library
- ReferenceBlock: radia_simulink_library/Reduced Models/Motor Angle Family
- Description: Periodic motor ROM with native MEX interpolation, state update, and torque evaluation.
- Use when: user needs periodic motor rom with native mex interpolation, state update, and torque evaluation..
- Metadata quality: medium

### Nonlinear HDiv-MMM Reactor

- Block: [[blocks/Nonlinear_HDiv-MMM_Reactor]]
- Library: radia_simulink_library
- ReferenceBlock: radia_simulink_library/Reduced Models/Nonlinear HDiv-MMM Reactor
- Description: Native nonlinear HDiv-MMM reactor. Input is winding current; outputs are terminal voltage, flux linkage, differential...
- Use when: user needs native nonlinear hdiv-mmm reactor. input is winding current; outputs are terminal voltage, flux linkage, d...
- Avoid when: user asks only for a primitive nonlinear hdiv-mmm reactor experiment.
- Metadata quality: high

## Low-Confidence Blocks

### Field Stats

- Library: radia_simulink_library
- ReferenceBlock: radia_simulink_library/Utilities/Field Stats
- Metadata quality: low
- Guidance: Available from customer library. Select when intent matches.

## Related Categories

- [[uncategorized]]
- [[plant-models]]
- [[power]]
- [[control]]
