# Uncategorized Blocks

Use these blocks for uncategorized blocks.

## Recommended Blocks

### Field Study Configuration

- Block: [[blocks/Field_Study_Configuration]]
- Library: radia_simulink_library
- ReferenceBlock: radia_simulink_library/Coupling/Field Study Configuration
- Description: Compile electrostatic, multi-conductor force, current-flow, steady/transient-heat, or linear/nonlinear harmonic-eddy ...
- Use when: user needs compile electrostatic, multi-conductor force, current-flow, steady/transient-heat, or linear/nonlinear har...
- Avoid when: user asks only for a primitive field study configuration experiment.
- Metadata quality: high

### Winding Dictionary

- Block: [[blocks/Winding_Dictionary]]
- Library: radia_simulink_library
- ReferenceBlock: radia_simulink_library/Coupling/Winding Dictionary
- Description: Compile winding names, .vol regions, turns, polarity, parallel paths, resistance, and circuit terminals to a fixed-wi...
- Use when: user needs compile winding names, .vol regions, turns, polarity, parallel paths, resistance, and circuit terminals to...
- Avoid when: user asks only for a primitive winding dictionary experiment.
- Metadata quality: high

### Material Dictionary

- Block: [[blocks/Material_Dictionary]]
- Library: radia_simulink_library
- ReferenceBlock: radia_simulink_library/Material Models/Material Dictionary
- Description: Compile MATLAB material and region dictionaries against a Netgen .vol mesh. The output is a fixed-width numeric Bus f...
- Use when: user needs compile matlab material and region dictionaries against a netgen .vol mesh. the output is a fixed-width nu...
- Avoid when: user asks only for a primitive material dictionary experiment.
- Metadata quality: high

### Electromagnet

- Block: [[blocks/Electromagnet]]
- Library: radia_simulink_library
- ReferenceBlock: radia_simulink_library/Applications/Electromagnet
- Description: Electromagnet batch analysis. Create settings with radia.simulink.writeApplicationConfig; a rising trigger executes t...
- Use when: user needs electromagnet batch analysis. create settings with radia.simulink.writeapplicationconfig; a rising trigger...
- Metadata quality: medium

### Field Study

- Block: [[blocks/Field_Study]]
- Library: radia_simulink_library
- ReferenceBlock: radia_simulink_library/Applications/Field Study
- Description: Field Study batch analysis. Create settings with radia.simulink.writeApplicationConfig; a rising trigger executes the...
- Use when: user needs field study batch analysis. create settings with radia.simulink.writeapplicationconfig; a rising trigger e...
- Metadata quality: medium

### PCB PEEC

- Block: [[blocks/PCB_PEEC]]
- Library: radia_simulink_library
- ReferenceBlock: radia_simulink_library/Applications/PCB PEEC
- Description: PCB PEEC batch analysis. Create settings with radia.simulink.writeApplicationConfig; a rising trigger executes the va...
- Use when: user needs pcb peec batch analysis. create settings with radia.simulink.writeapplicationconfig; a rising trigger exec...
- Metadata quality: medium

### Stream Function

- Block: [[blocks/Stream_Function]]
- Library: radia_simulink_library
- ReferenceBlock: radia_simulink_library/Applications/Stream Function
- Description: Stream Function batch analysis. Create settings with radia.simulink.writeApplicationConfig; a rising trigger executes...
- Use when: user needs stream function batch analysis. create settings with radia.simulink.writeapplicationconfig; a rising trigg...
- Metadata quality: medium

### LTspice Circuit

- Block: [[blocks/LTspice_Circuit]]
- Library: radia_simulink_library
- ReferenceBlock: radia_simulink_library/LTspice/LTspice Circuit
- Use when: user asks for ltspice circuit.
- Metadata quality: medium

### Temperature-Dependent BH

- Block: [[blocks/Temperature-Dependent_BH]]
- Library: radia_simulink_library
- ReferenceBlock: radia_simulink_library/Material Models/Temperature-Dependent BH
- Description: Temperature-dependent BH law. Select formula or LUT in radia_bh_config.
- Use when: user needs temperature-dependent bh law. select formula or lut in radia_bh_config..
- Metadata quality: medium

### Adjoint Topology Optimization

- Block: [[blocks/Adjoint_Topology_Optimization]]
- Library: radia_simulink_library
- ReferenceBlock: radia_simulink_library/Optimization/Adjoint Topology Optimization
- Description: Checked adjoint topology optimization through the MMA/SQP runner.
- Use when: user needs checked adjoint topology optimization through the mma/sqp runner..
- Metadata quality: medium

### Optuna Optimization

- Block: [[blocks/Optuna_Optimization]]
- Library: radia_simulink_library
- ReferenceBlock: radia_simulink_library/Optimization/Optuna Optimization
- Description: Incremental MATLAB Optuna study with Simulink-native telemetry.
- Use when: user needs incremental matlab optuna study with simulink-native telemetry..
- Metadata quality: medium

### Sheet Metal Optimization

- Block: [[blocks/Sheet_Metal_Optimization]]
- Library: radia_simulink_library
- ReferenceBlock: radia_simulink_library/Optimization/Sheet Metal Optimization
- Description: Native NGSolve/Cubit sheet-metal optimization through radia.optuna.SheetMetalRunner.
- Use when: user needs native ngsolve/cubit sheet-metal optimization through radia.optuna.sheetmetalrunner..
- Metadata quality: medium

## Low-Confidence Blocks

### Electromagnet Topology Optimization

- Library: radia_simulink_library
- ReferenceBlock: radia_simulink_library/Applications/Electromagnet Topology Optimization
- Metadata quality: low
- Guidance: Available from customer library. Select when intent matches.

### Stream Function Optimization

- Library: radia_simulink_library
- ReferenceBlock: radia_simulink_library/Applications/Stream Function Optimization
- Metadata quality: low
- Guidance: Available from customer library. Select when intent matches.

### Optuna Monitor

- Library: radia_simulink_library
- ReferenceBlock: radia_simulink_library/Optimization/Optuna Monitor
- Metadata quality: low
- Guidance: Available from customer library. Select when intent matches.

## Related Categories

- [[plant-models]]
- [[signal-processing]]
- [[power]]
