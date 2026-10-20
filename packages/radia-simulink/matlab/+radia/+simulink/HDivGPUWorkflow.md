# Shared HDiv GPU workflow block

`radia_hdiv_gpu_sfun` wraps `radia.gpu.PreparedHDivNonlinear`. It owns one prepared GPU context for a simulation and releases it at termination, including error termination. It does not implement a separate constitutive or Newton algorithm.

The first input is a Boolean rising-edge trigger. The second is a finite field scale. Outputs are FE coefficients, mean magnetization in A/m, original-equation relative residual, Newton iteration count and an int32 status. Status 0 means idle; 2 means a completed solve. Outputs are held and become visible at the next discrete sample. SimState restoration is unsupported.

Supply the explicit GPU backend, checked FE setup, order, degree-of-freedom count, applied field, BH table, supported reference preconditioner, positive sample time and output parent directory. Each simulation creates a separate child directory containing `run.log` and atomically replaced `result.json`. The normalized host setup retains NGSolve ownership of FE maps and quadrature. Unsupported options fail rather than selecting another backend.

This block is a field-solve adapter. It does not replace a coil voltage/flux/differential-inductance reactor model. A production SLX, distribution qualification and release review remain separate gates. The CPU oracle used to certify validation results is not a GPU-user prerequisite.

`radia_hdiv_gpu_workflow.slx` is a standalone sample using a synthetic hex order-2 checked normalized FE payload in its model workspace. The two Boolean rising edges request the same applied field, so the second request demonstrates prepared-state reuse. The model has no Python call during simulation. After configuring the qualified native payload and its DLL dependencies, open the model from the installed MATLAB resource directory and simulate it. Its declared relative output parent is `radia_hdiv_gpu_runs` in the current writable working directory; each invocation creates a distinct child. Change `hdivConfig.RunDirectory` or use `Simulink.SimulationInput` to choose another destination.

The model was authored and checked with the official Simulink Agentic Toolkit. Workspace data use ordinary numeric structures rather than serialized time-series interpolation function handles. Its saved settings contain no machine-specific paths. This sample's validation does not certify mixed pyramid meshes or the separate coil reactor example.
