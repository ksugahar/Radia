# Customer Library Reuse Index

This project has declared reusable Simulink libraries. Prefer these blocks when they match the modeling intent.

## Policy

The active policy mode is defined in `.satk/block-policy.json`.

- Always use customer library blocks when available. 
- Do NOT make domain-level judgments about library relevance.
- Never fall back to built-in primitives if the same block exists in a declared library.
- Only use built-in blocks when NO equivalent exists in any declared customer library after your search.
- Do not invent customer block names.
- If uncertain, inspect the relevant category page or ask the user.
- CRITICAL: Before using ANY block, search this index and the category pages for that specific block type first

## Libraries

- radia_simulink_library: Radia-NGSolve solver, circuit, and Simulink integration blocks.

## Commonly Used Blocks

Showing 9 of 22 blocks — see category pages for the full catalog.

- [[blocks/Hysteretic_LTspice_Plant]] — Hysteretic LTspice Plant from radia_simulink_library
- [[blocks/Motor]] — Motor batch analysis. Create settings with radia.simulink.writeApplicationConfig; a rising trigger executes the valid... from radia_simulink_library
- [[blocks/Nonlinear_HDiv-MMM_Reactor]] — Native nonlinear HDiv-MMM reactor. Input is winding current; outputs are terminal voltage, flux linkage, differential... from radia_simulink_library
- [[blocks/Motor_Angle_Family]] — Periodic motor ROM with native MEX interpolation, state update, and torque evaluation. from radia_simulink_library
- [[blocks/Field_Study_Configuration]] — Compile electrostatic, multi-conductor force, current-flow, steady/transient-heat, or linear/nonlinear harmonic-eddy ... from radia_simulink_library
- [[blocks/Winding_Dictionary]] — Compile winding names, .vol regions, turns, polarity, parallel paths, resistance, and circuit terminals to a fixed-wi... from radia_simulink_library
- [[blocks/Material_Dictionary]] — Compile MATLAB material and region dictionaries against a Netgen .vol mesh. The output is a fixed-width numeric Bus f... from radia_simulink_library
- [[blocks/Electromagnet]] — Electromagnet batch analysis. Create settings with radia.simulink.writeApplicationConfig; a rising trigger executes t... from radia_simulink_library
- [[blocks/Field_Study]] — Field Study batch analysis. Create settings with radia.simulink.writeApplicationConfig; a rising trigger executes the... from radia_simulink_library

## Categories

- [[uncategorized]] (15 blocks) — blocks with insufficient metadata for confident categorization
- [[plant-models]] (3 blocks) — physical dynamics, thermal, mechanical models
- [[signal-processing]] (3 blocks) — filters, scaling, interpolation, signal conditioning
- [[power]] (1 blocks) — inverters, converters, motors, power stage components
